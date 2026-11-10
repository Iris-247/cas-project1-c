/*
 * cluster.c -- Part I: grid, NFW dark matter, BCG, hydrostatic equilibrium,
 * temperature profile and turbulent pressure support.
 */
#include "cluster.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void *xcalloc(size_t n, size_t size)
{
    void *p = calloc(n, size);
    if (!p) {
        fprintf(stderr, "error: out of memory (%zu x %zu bytes)\n", n, size);
        exit(EXIT_FAILURE);
    }
    return p;
}

int cluster_init(Cluster *c, int n, double rmax_kpc,
                 double rVir_kpc, double r500_kpc)
{
    memset(c, 0, sizeof *c);

    if (n < 8) {
        fprintf(stderr, "error: jmax must be at least 8\n");
        return -1;
    }
    c->g.n    = n;
    c->g.rmin = 0.0;
    c->g.rmax = rmax_kpc * KPC_CM;
    c->g.dr   = (c->g.rmax - c->g.rmin) / (double)(n - 1);

    c->rs       = 435.7 * KPC_CM;
    c->rhoNFW0  = 7.35e-26;
    c->MBCG     = 1.0e12 * MSUN_G;
    c->aBCG     = 12.0 / (1.0 + sqrt(2.0)) * KPC_CM;
    c->rVir     = rVir_kpc * KPC_CM;
    c->r500     = r500_kpc * KPC_CM;
    c->mu       = MU_GAS;
    c->Tmg      = 8.9e7;
    c->vt       = 500.0e5;            /* 500 km/s */
    c->fbTarget = 0.16;
    c->rho0     = 0.0;
    c->maxTurbIter = 100;
    c->jVir     = (int)floor(c->rVir / c->g.dr + 0.5);
    if (c->jVir > n - 1) c->jVir = n - 1;
    if (c->jVir < 1)     c->jVir = 1;

    c->g.r   = xcalloc((size_t)n, sizeof(double));
    c->g.rc  = xcalloc((size_t)(n - 1), sizeof(double));
    c->g.vol = xcalloc((size_t)(n - 1), sizeof(double));

    c->rhoNFW  = xcalloc((size_t)(n - 1), sizeof(double));
    c->mDM     = xcalloc((size_t)n, sizeof(double));
    c->mDMana  = xcalloc((size_t)n, sizeof(double));
    c->mBCG    = xcalloc((size_t)n, sizeof(double));
    c->mTot    = xcalloc((size_t)n, sizeof(double));
    c->T       = xcalloc((size_t)(n - 1), sizeof(double));
    c->rho     = xcalloc((size_t)(n - 1), sizeof(double));
    c->rhoAna  = xcalloc((size_t)(n - 1), sizeof(double));
    c->pnt     = xcalloc((size_t)(n - 1), sizeof(double));
    c->scratch = xcalloc((size_t)(n - 1), sizeof(double));

    return 0;
}

void cluster_free(Cluster *c)
{
    free(c->g.r);
    free(c->g.rc);
    free(c->g.vol);
    free(c->rhoNFW);
    free(c->mDM);
    free(c->mDMana);
    free(c->mBCG);
    free(c->mTot);
    free(c->T);
    free(c->rho);
    free(c->rhoAna);
    free(c->pnt);
    free(c->scratch);
    memset(c, 0, sizeof *c);
}

void build_grid(Cluster *c)
{
    const int n = c->g.n;
    for (int j = 0; j < n; ++j)
        c->g.r[j] = c->g.rmin + (double)j * c->g.dr;
    for (int j = 0; j < n - 1; ++j) {
        c->g.rc[j]  = 0.5 * (c->g.r[j] + c->g.r[j + 1]);
        c->g.vol[j] = (4.0 * PI_C / 3.0) *
                      (pow(c->g.r[j + 1], 3.0) - pow(c->g.r[j], 3.0));
    }
}

double nfw_mass_analytic(const Cluster *c, double r)
{
    if (r <= 0.0) return 0.0;
    const double x = r / c->rs;
    return 4.0 * PI_C * c->rhoNFW0 * pow(c->rs, 3.0) *
           (log1p(x) - x / (1.0 + x));
}

void build_dark_matter(Cluster *c)
{
    const int n = c->g.n;

    for (int j = 0; j < n - 1; ++j) {
        const double x = c->g.rc[j] / c->rs;
        c->rhoNFW[j] = c->rhoNFW0 / (x * (1.0 + x) * (1.0 + x));
    }

    /* Shell 0 straddles r=0 where NFW diverges.  Use the analytic mass at
     * r[1] for that first shell, then sum the remaining shells numerically.
     * This removes the O(dr) bias of the innermost shell. */
    c->mDM[0] = 0.0;
    c->mDM[1] = nfw_mass_analytic(c, c->g.r[1]);
    for (int j = 1; j < n - 1; ++j)
        c->mDM[j + 1] = c->mDM[j] + c->rhoNFW[j] * c->g.vol[j];

    c->mDMana[0] = 0.0;
    for (int j = 1; j < n; ++j)
        c->mDMana[j] = nfw_mass_analytic(c, c->g.r[j]);
}

void build_bcg(Cluster *c)
{
    const int n = c->g.n;
    const double a = c->aBCG;
    for (int j = 0; j < n; ++j) {
        const double r = c->g.r[j];
        c->mBCG[j] = c->MBCG * r * r / ((r + a) * (r + a));
        c->mTot[j] = c->mDM[j] + c->mBCG[j];
    }
}

double temperature_at(const Cluster *c, double r)
{
    const double x = r / c->r500;
    const double y = pow(x / 0.045, 1.9);
    return c->Tmg * 1.35 * (y + 0.45) / (y + 1.0) /
           pow(1.0 + pow(x / 0.6, 2.0), 0.45);
}

void build_temperature(Cluster *c)
{
    for (int j = 0; j < c->g.n - 1; ++j)
        c->T[j] = temperature_at(c, c->g.rc[j]);
}

/* Density at the first shell centre, given the central density rho0.
 * For the isothermal model the analytic NFW solution gives the exact offset
 * between r=0 and rc[0]; for the other models a half-cell hydrostatic
 * estimate is used (the correction is < 1 per cent because rc[0] ~ 0.3 kpc). */
static double first_shell_rho(const Cluster *c, ModelKind kind, double rho0)
{
    const double rc = c->g.rc[0];
    double dlnrho;

    if (kind == MODEL_ISO) {
        const double A = 4.0 * PI_C * G_CGS * c->rhoNFW0 * c->rs * c->rs *
                         c->mu * MP_CGS / (KB_CGS * c->Tmg);
        const double x = rc / c->rs;
        const double f = 1.0 - log1p(x) / x;
        return rho0 * exp(-A * f);
    }

    const double M = nfw_mass_analytic(c, rc) +
                     c->MBCG * rc * rc / ((rc + c->aBCG) * (rc + c->aBCG));
    const double T0 = (kind == MODEL_ISO_BCG) ? c->Tmg : c->T[0];
    dlnrho = G_CGS * M * c->mu * MP_CGS / (2.0 * KB_CGS * T0 * rc);
    return rho0 * exp(-dlnrho);
}

/* One outward sweep of
 *   dP/dr = -rho G M(<r)/r^2,   P = rho kT/(mu mp) [+ (1/3) rho vt^2]
 * discretised with P_{j+1} = P_j - dr * G M_{j+1}/r_{j+1}^2 * (rho_j+rho_{j+1})/2.
 * With a constant vt the turbulent term is linear in rho, so this is a
 * positive, second-order-in-space scheme.  For the iterative variant prev[]
 * holds the previous iterate used for P_nt. */
static void integrate_pressure(const Cluster *c, const double *mass,
                               const double *Tarr, double Tconst,
                               double shell0, int use_pnt,
                               const double *prev, double *rho, double *pnt)
{
    const int n = c->g.n;
    const double dr = c->g.dr;
    const double v2 = c->vt * c->vt;

    rho[0] = shell0;
    if (use_pnt && pnt) pnt[0] = (1.0 / 3.0) * shell0 * v2;

    for (int j = 0; j < n - 2; ++j) {
        const double Tj  = Tarr ? Tarr[j]     : Tconst;
        const double Tj1 = Tarr ? Tarr[j + 1] : Tconst;
        const double Cth_j  = KB_CGS * Tj  / (c->mu * MP_CGS);
        const double Cth_j1 = KB_CGS * Tj1 / (c->mu * MP_CGS);
        const double pntj   = use_pnt ? (1.0 / 3.0) * prev[j] * v2 : 0.0;
        const double pntj1  = use_pnt ? (1.0 / 3.0) * prev[j + 1] * v2 : 0.0;
        const double pj     = rho[j] * Cth_j + pntj;
        const double A      = dr * G_CGS * mass[j + 1] /
                              (c->g.r[j + 1] * c->g.r[j + 1]);
        rho[j + 1] = (pj - 0.5 * A * rho[j] - pntj1) /
                     (Cth_j1 + 0.5 * A);
        if (use_pnt && pnt)
            pnt[j + 1] = (1.0 / 3.0) * rho[j + 1] * v2;
    }
}

void integrate_model(Cluster *c, ModelKind kind, double rho0,
                     double *rho, double *pnt, int *iters)
{
    const int n = c->g.n;
    const double *mass = (kind == MODEL_ISO) ? c->mDM : c->mTot;
    const double *Tarr = (kind == MODEL_NONISO || kind == MODEL_TURB)
                         ? c->T : NULL;
    const double shell0 = first_shell_rho(c, kind, rho0);

    *iters = 0;

    if (kind != MODEL_TURB) {
        integrate_pressure(c, mass, Tarr, c->Tmg, shell0, 0, NULL, rho, NULL);
        if (pnt)
            for (int j = 0; j < n - 1; ++j) pnt[j] = 0.0;
        return;
    }

    /* Fixed point: start from the non-turbulent solution, then repeatedly
     * use P_nt = (1/3) rho_prev vt^2.  For constant vt the map is linear in
     * rho and converges to the self-consistent solution. */
    integrate_pressure(c, mass, Tarr, c->Tmg, shell0, 0, NULL,
                       c->scratch, NULL);

    for (int it = 1; it <= c->maxTurbIter; ++it) {
        integrate_pressure(c, mass, Tarr, c->Tmg, shell0, 1,
                           c->scratch, rho, pnt);
        double rel = 0.0;
        for (int j = 0; j < n - 1; ++j) {
            const double d = fabs(rho[j] - c->scratch[j]);
            const double s = fabs(rho[j]) > 1.0e-300 ? fabs(rho[j]) : 1.0e-300;
            if (d / s > rel) rel = d / s;
        }
        memcpy(c->scratch, rho, (size_t)(n - 1) * sizeof(double));
        *iters = it;
        if (rel < 1.0e-13) break;
    }
}

void compute_gas_mass(const Cluster *c, const double *rho, double *mgas)
{
    mgas[0] = 0.0;
    for (int j = 0; j < c->g.n - 1; ++j)
        mgas[j + 1] = mgas[j] + rho[j] * c->g.vol[j];
}

double fb_at_rvir(const Cluster *c, const double *rho)
{
    double mgas = 0.0;
    for (int j = 0; j < c->jVir; ++j)
        mgas += rho[j] * c->g.vol[j];

    const double md = c->mDM[c->jVir];
    const double mb = c->mBCG[c->jVir];
    return (mgas + mb) / (md + mb + mgas);
}

double solve_rho0(Cluster *c, ModelKind kind)
{
    const int n = c->g.n;
    const int nc = n - 1;
    double *tmp  = xcalloc((size_t)nc, sizeof(double));
    double *tpnt = xcalloc((size_t)nc, sizeof(double));
    double lo = 1.0e-28, hi = 1.0e-23;
    int iters = 0;

    for (int it = 0; it < 200; ++it) {
        const double mid = sqrt(lo * hi);
        integrate_model(c, kind, mid, tmp, tpnt, &iters);
        const double f = fb_at_rvir(c, tmp);
        if (f < c->fbTarget) lo = mid;
        else                 hi = mid;
        if (hi / lo < 1.0 + 1.0e-13) break;
    }

    free(tmp);
    free(tpnt);
    return sqrt(lo * hi);
}

double hse_residual(const Cluster *c, const double *rho, const double *pnt,
                    const double *mass, const double *Tarr, double *mean)
{
    const int n = c->g.n;
    double maxres = 0.0, sum = 0.0;
    int count = 0;

    /* Independent check: centred pressure difference at the shell centres,
     * while the integration itself uses fluxes at the boundaries. */
    for (int j = 1; j < n - 2; ++j) {
        if (c->g.rc[j] < 10.0 * KPC_CM) continue;
        const double Tjm1 = Tarr ? Tarr[j-1] : c->Tmg;
        const double Tjp1 = Tarr ? Tarr[j+1] : c->Tmg;
        const double pjm1 = rho[j-1] * KB_CGS * Tjm1 / (c->mu * MP_CGS) +
                            (pnt ? pnt[j-1] : 0.0);
        const double pjp1 = rho[j+1] * KB_CGS * Tjp1 / (c->mu * MP_CGS) +
                            (pnt ? pnt[j+1] : 0.0);
        const double dpdr = (pjp1 - pjm1) /
                            (c->g.rc[j+1] - c->g.rc[j-1]);
        const double Mc = 0.5 * (mass[j] + mass[j + 1]);
        const double grav = rho[j] * G_CGS * Mc /
                            (c->g.rc[j] * c->g.rc[j]);
        const double res = fabs(dpdr + grav) / fabs(grav);
        if (res > maxres) maxres = res;
        sum += res;
        ++count;
    }

    if (mean) *mean = count ? sum / (double)count : 0.0;
    return maxres;
}
