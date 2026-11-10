/*
 * diffusion.c -- Part II: explicit finite-volume diffusion of iron on the
 * shell grid, with the SNIa + stellar-wind source term.
 *
 *   d rho_Fe/dt = 1/1.4 * 1/r^2 d/dr ( r^2 D rho dZ_Fe/dr ) + S_Fe(r,t)
 *   Z_Fe = 1.4 rho_Fe / rho
 *
 * The scheme is the one used in the lecture: fluxes at cell boundaries,
 * zero-flux (outflow) boundary conditions, and dt = C dr^2/(2 D), C = 0.4.
 */
#include "diffusion.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define FE_DT_SAFETY   0.40
#define FE_TNOW_YR     13.7e9
#define FE_TSTART_YR    8.7e9

void fe_default_params(FeParams *p)
{
    p->kappa0      = 1.32e29;      /* 0.11 * 260 km/s * 15 kpc */
    p->vary        = 0;
    p->rscala      = 30.0 * KPC_CM;
    p->coreFactor  = 0.6;
    p->snu         = 0.15;
    p->alphaEvolve = 0;
}

double fe_observed_zfe(const Cluster *c, double r, int excess)
{
    (void)c;
    const double x  = r / (80.0 * KPC_CM);
    const double x3 = x * x * x;
    double z = 0.30 * 1.4 * ZSUN * (2.2 + x3) / (1.0 + x3);

    if (excess) {
        z -= 0.40 * ZSUN;
        if (z < 0.0) z = 0.0;
    }
    return z;
}

/* Total source coefficient per unit stellar density:
 * alpha_star * Z_star/1.4 + alpha_SN * M_Fe,Ia/M_Ia, at t=t_now (s^-1). */
static double fe_alpha_total(double snu)
{
    const double astar = 4.7e-20;                    /* s^-1 */
    const double asn   = 4.436e-20 * (snu / 7.5);    /* s^-1 */
    return astar * ZSUN / 1.4 + asn * 0.74 / 1.4;
}

double fe_source_total_theory(const Cluster *c, double snu, double t)
{
    const double rmaxm = c->g.r[c->g.n - 1];
    const double a = c->aBCG;
    const double mbcg_tot = c->MBCG * rmaxm * rmaxm /
                            ((rmaxm + a) * (rmaxm + a));
    return fe_alpha_total(snu) * mbcg_tot * t / MSUN_G;
}

static double fe_source_coef(const FeParams *p, double t)
{
    double ts = 1.0;
    if (p->alphaEvolve) {
        ts = (FE_TSTART_YR * YEAR_S + t) / (FE_TNOW_YR * YEAR_S);
    }
    const double astar = 4.7e-20 * pow(ts, -1.26);
    const double asn   = 4.436e-20 * (p->snu / 7.5) * pow(ts, -1.10);
    return astar * ZSUN / 1.4 + asn * 0.74 / 1.4;
}

static double kappa_at(const FeParams *p, double r)
{
    if (!p->vary) return p->kappa0;
    const double q = r / p->rscala;
    return p->kappa0 * (1.0 - p->coreFactor * exp(-q * q));
}

static void fe_record(const Cluster *c, const double *rhoFe,
                      FeSnapshot *s, double t)
{
    double m100 = 0.0, mtot = 0.0;
    for (int j = 0; j < c->g.n - 1; ++j) {
        const double m = rhoFe[j] * c->g.vol[j];
        mtot += m;
        if (c->g.rc[j] <= 100.0 * KPC_CM) m100 += m;
    }
    s->t      = t;
    s->mfe100 = m100 / MSUN_G;
    s->mfetot = mtot / MSUN_G;
}

typedef struct {
    const Cluster  *c;
    const FeParams *par;
    const double   *rhoGas;
    double *zfe, *rhoFe, *D, *flux, *rstar;
    int     useSource, useDiffusion;
    double  dr;
} FeWork;

static int fe_step(FeWork *w, double t, double dt)
{
    const Cluster *c = w->c;
    const int n  = c->g.n;
    const int nc = n - 1;
    const double dr = c->g.dr;

    if (w->useDiffusion) {
        for (int b = 0; b < n; ++b) {
            if (b == 0 || b == n - 1) {
                w->flux[b] = 0.0;              /* zero-flux boundaries */
            } else {
                const double grad = (w->zfe[b] - w->zfe[b - 1]) / dr;
                const double rb   = 0.5 * (w->rhoGas[b - 1] + w->rhoGas[b]);
                const double Db   = 0.5 * (w->D[b - 1] + w->D[b]);
                const double r    = c->g.r[b];
                w->flux[b] = r * r * Db * rb * grad;
            }
        }
    }

    const double coef = w->useSource ? fe_source_coef(w->par, t) : 0.0;

    for (int j = 0; j < nc; ++j) {
        const double r0 = c->g.r[j], r1 = c->g.r[j + 1];
        const double dV = (r1 * r1 * r1 - r0 * r0 * r0) / 3.0;
        double delta = 0.0;

        if (w->useDiffusion)
            delta += (dt / 1.4) * (w->flux[j + 1] - w->flux[j]) / dV;
        if (w->useSource)
            delta += dt * coef * w->rstar[j];

        w->rhoFe[j] += delta;
        if (!isfinite(w->rhoFe[j]) || w->rhoFe[j] < 0.0) {
            fprintf(stderr,
                    "error: negative/non-finite rho_Fe at r=%.3f kpc "
                    "(rho_Fe=%.6e, step dt=%.3e s)\n",
                    c->g.rc[j] / KPC_CM, w->rhoFe[j], dt);
            return -1;
        }
        w->zfe[j] = 1.4 * w->rhoFe[j] / w->rhoGas[j];
    }
    return 0;
}

int fe_run(const Cluster *c, const double *rhoGas,
           int useSource, int useDiffusion, FeInit init,
           const FeParams *par,
           const double *snapTimes, int nsnap,
           FeSnapshot *snaps, const char *path)
{
    const int nc = c->g.n - 1;
    FeWork w;
    double *zsave = NULL;
    double Dmax = 0.0;
    int rc = 0;

    if (nsnap < 1) return -1;

    memset(&w, 0, sizeof w);
    w.c = c;
    w.par = par;
    w.rhoGas = rhoGas;
    w.useSource = useSource;
    w.useDiffusion = useDiffusion;
    w.dr = c->g.dr;

    w.zfe   = calloc((size_t)nc, sizeof(double));
    w.rhoFe = calloc((size_t)nc, sizeof(double));
    w.D     = calloc((size_t)nc, sizeof(double));
    w.flux  = calloc((size_t)c->g.n, sizeof(double));
    w.rstar = calloc((size_t)nc, sizeof(double));
    zsave   = calloc((size_t)nsnap * (size_t)nc, sizeof(double));
    if (!w.zfe || !w.rhoFe || !w.D || !w.flux || !w.rstar || !zsave) {
        fprintf(stderr, "error: out of memory in fe_run\n");
        rc = -1;
        goto done;
    }

    for (int j = 0; j < nc; ++j) {
        const double z = (init == FE_INIT_OBSERVED)
                         ? fe_observed_zfe(c, c->g.rc[j], 1)
                         : 0.0;
        w.zfe[j]   = z;
        w.rhoFe[j] = rhoGas[j] * z / 1.4;
        w.D[j]     = kappa_at(par, c->g.rc[j]);
        if (w.D[j] > Dmax) Dmax = w.D[j];

        /* Shell-averaged stellar density: (M_*(r_{j+1}) - M_*(r_j))/V_j.
         * This makes the volume integral of the source exactly M_* (<rmax),
         * which matters because the Hernquist profile diverges as 1/r. */
        w.rstar[j] = (c->mBCG[j + 1] - c->mBCG[j]) / c->g.vol[j];
    }

    double dtStable = 1.0e7 * YEAR_S;
    if (useDiffusion) {
        if (!(Dmax > 0.0)) {
            fprintf(stderr, "error: diffusion coefficient must be positive\n");
            rc = -1;
            goto done;
        }
        dtStable = FE_DT_SAFETY * c->g.dr * c->g.dr / (2.0 * Dmax);
    }

    printf("  Fe run: source=%s diffusion=%s init=%s  (dt=%.3e s, "
           "%.1f steps/Gyr)\n",
           useSource ? "yes" : "no", useDiffusion ? "yes" : "no",
           init == FE_INIT_OBSERVED ? "observed" : "zero",
           dtStable, GYR_S / dtStable);

    double t = 0.0;
    int is = 0;
    while (is < nsnap) {
        const double target = snapTimes[is];

        /* record an initial (t=0) snapshot if requested */
        if (t >= target - 1.0e-12 * fabs(target)) {
            t = target;
            fe_record(c, w.rhoFe, &snaps[is], t);
            for (int j = 0; j < nc; ++j)
                zsave[(size_t)is * (size_t)nc + (size_t)j] =
                    w.zfe[j] / ZSUN;
            ++is;
            continue;
        }

        double dt = dtStable;
        int hit = 0;
        if (t + dt >= target) {
            dt = target - t;
            hit = 1;
        }
        if (fe_step(&w, t, dt) != 0) { rc = -1; goto done; }
        t = hit ? target : t + dt;
        if (hit) {
            fe_record(c, w.rhoFe, &snaps[is], t);
            for (int j = 0; j < nc; ++j)
                zsave[(size_t)is * (size_t)nc + (size_t)j] =
                    w.zfe[j] / ZSUN;
            ++is;
        }
    }

    if (path) {
        FILE *f = fopen(path, "w");
        if (!f) {
            fprintf(stderr, "error: cannot open %s\n", path);
            rc = -1;
            goto done;
        }
        fprintf(f, "# Part II: iron abundance Z_Fe/Z_sun, "
                   "source=%d diffusion=%d init=%s\n",
                useSource, useDiffusion,
                init == FE_INIT_OBSERVED ? "observed" : "zero");
        fprintf(f, "# r_kpc");
        for (int k = 0; k < nsnap; ++k)
            fprintf(f, " Z_t%.3fGyr", snaps[k].t / GYR_S);
        fprintf(f, "\n");
        for (int j = 0; j < nc; ++j) {
            fprintf(f, "%.8e", c->g.rc[j] / KPC_CM);
            for (int k = 0; k < nsnap; ++k)
                fprintf(f, " %.8e",
                        zsave[(size_t)k * (size_t)nc + (size_t)j]);
            fprintf(f, "\n");
        }
        fclose(f);
    }

    printf("    %8s %16s %16s\n", "t[Gyr]", "MFe(<100kpc)", "MFe(tot)");
    for (int k = 0; k < nsnap; ++k)
        printf("    %8.3f %16.6e %16.6e\n", snaps[k].t / GYR_S,
               snaps[k].mfe100, snaps[k].mfetot);

    if (useDiffusion && !useSource && nsnap > 1) {
        const double drift = fabs(snaps[nsnap - 1].mfetot - snaps[0].mfetot) /
                             snaps[0].mfetot;
        printf("    Fe-mass conservation over the whole grid: %.3e\n", drift);
    }

done:
    free(w.zfe);
    free(w.rhoFe);
    free(w.D);
    free(w.flux);
    free(w.rstar);
    free(zsave);
    return rc;
}
