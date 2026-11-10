/*
 * main.c -- driver for the C implementation of the Computational
 * Astrophysics 1st project:
 *
 *   Part I : hydrostatic equilibrium of the hot gas in a galaxy cluster
 *            (NFW dark matter + Hernquist BCG + isothermal / T(r) gas)
 *   Part II: turbulent diffusion of iron in the ICM, with the SNIa and
 *            stellar-wind source terms.
 *
 * Build:  make          Run:  ./project1 --check
 */
#include "cluster.h"
#include "diffusion.h"

#include <errno.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _WIN32
#include <direct.h>
#define MAKE_DIR(p) _mkdir(p)
#else
#include <sys/stat.h>
#define MAKE_DIR(p) mkdir((p), 0777)
#endif

typedef struct {
    int    jmax;
    double rmax_kpc, rVir_kpc, r500_kpc;
    double fb;
    double rho0_override;
    double vt_kms;
    double snu;
    double kappa;
    int    kappa_vary;
    int    alpha_evolve;
    int    fe;
    int    check;
    int    conv;
    int    turb_iters;
    const char *outdir;
} Options;

static void usage(const char *prog)
{
    printf("Usage: %s [options]\n", prog);
    printf("  -n, --jmax N        number of grid boundaries (default 5000)\n");
    printf("  -R, --rmax KPC      outer radius (default 2800 kpc)\n");
    printf("      --rvir KPC      virial radius for fb (default 2800 kpc)\n");
    printf("      --r500 KPC      R500 for T(r) (default 1400 kpc)\n");
    printf("      --fb X          target baryon fraction (default 0.16)\n");
    printf("      --rho0 X        fixed central gas density [g/cm^3]\n");
    printf("      --vt KM/S       turbulent velocity (default 500)\n");
    printf("      --turb-iters N  max fixed-point iterations (default 100)\n");
    printf("      --snu X         SNIa rate in SNu (default 0.15)\n");
    printf("      --kappa X       diffusion coefficient [cm^2/s]\n");
    printf("      --kappa-vary    use kappa(r), reduced in the core\n");
    printf("      --alpha-evolve  use the time-dependent alpha(t) source\n");
    printf("      --no-fe         skip Part II\n");
    printf("      --check         run the self-tests and exit\n");
    printf("      --conv          grid-convergence table and exit\n");
    printf("  -o, --out DIR       output directory (default output)\n");
    printf("  -h, --help          this message\n");
}

static int need_value(int argc, char **argv, int *i, const char *opt)
{
    (void)argv;
    if (*i + 1 >= argc) {
        fprintf(stderr, "error: %s needs a value\n", opt);
        exit(EXIT_FAILURE);
    }
    return ++(*i);
}

static Options parse_args(int argc, char **argv)
{
    Options o;
    o.jmax = 5000;
    o.rmax_kpc = 2800.0;
    o.rVir_kpc = 2800.0;
    o.r500_kpc = 1400.0;
    o.fb = 0.16;
    o.rho0_override = -1.0;
    o.vt_kms = 500.0;
    o.snu = 0.15;
    o.kappa = 1.32e29;
    o.kappa_vary = 0;
    o.alpha_evolve = 0;
    o.fe = 1;
    o.check = 0;
    o.conv = 0;
    o.turb_iters = 100;
    o.outdir = "output";

    for (int i = 1; i < argc; ++i) {
        const char *a = argv[i];
        if (!strcmp(a, "-h") || !strcmp(a, "--help")) {
            usage(argv[0]);
            exit(EXIT_SUCCESS);
        } else if (!strcmp(a, "-n") || !strcmp(a, "--jmax")) {
            o.jmax = atoi(argv[need_value(argc, argv, &i, a)]);
        } else if (!strcmp(a, "-R") || !strcmp(a, "--rmax")) {
            o.rmax_kpc = atof(argv[need_value(argc, argv, &i, a)]);
        } else if (!strcmp(a, "--rvir")) {
            o.rVir_kpc = atof(argv[need_value(argc, argv, &i, a)]);
        } else if (!strcmp(a, "--r500")) {
            o.r500_kpc = atof(argv[need_value(argc, argv, &i, a)]);
        } else if (!strcmp(a, "--fb")) {
            o.fb = atof(argv[need_value(argc, argv, &i, a)]);
        } else if (!strcmp(a, "--rho0")) {
            o.rho0_override = atof(argv[need_value(argc, argv, &i, a)]);
        } else if (!strcmp(a, "--vt")) {
            o.vt_kms = atof(argv[need_value(argc, argv, &i, a)]);
        } else if (!strcmp(a, "--turb-iters")) {
            o.turb_iters = atoi(argv[need_value(argc, argv, &i, a)]);
        } else if (!strcmp(a, "--snu")) {
            o.snu = atof(argv[need_value(argc, argv, &i, a)]);
        } else if (!strcmp(a, "--kappa")) {
            o.kappa = atof(argv[need_value(argc, argv, &i, a)]);
        } else if (!strcmp(a, "--kappa-vary")) {
            o.kappa_vary = 1;
        } else if (!strcmp(a, "--alpha-evolve")) {
            o.alpha_evolve = 1;
        } else if (!strcmp(a, "--no-fe")) {
            o.fe = 0;
        } else if (!strcmp(a, "--check")) {
            o.check = 1;
            o.fe = 1;
        } else if (!strcmp(a, "--conv")) {
            o.conv = 1;
        } else if (!strcmp(a, "-o") || !strcmp(a, "--out")) {
            o.outdir = argv[need_value(argc, argv, &i, a)];
        } else {
            fprintf(stderr, "error: unknown option %s (try --help)\n", a);
            exit(EXIT_FAILURE);
        }
    }

    if (o.jmax < 8) {
        fprintf(stderr, "error: --jmax must be at least 8\n");
        exit(EXIT_FAILURE);
    }
    if (!(o.fb > 0.0 && o.fb < 1.0)) {
        fprintf(stderr, "error: --fb must be in (0,1)\n");
        exit(EXIT_FAILURE);
    }
    if (!(o.rmax_kpc > 0.0) || !(o.rVir_kpc > 0.0) ||
        o.rVir_kpc > o.rmax_kpc) {
        fprintf(stderr, "error: need 0 < rvir <= rmax\n");
        exit(EXIT_FAILURE);
    }
    if (!(o.r500_kpc > 0.0)) {
        fprintf(stderr, "error: --r500 must be positive\n");
        exit(EXIT_FAILURE);
    }
    if (o.snu < 0.0) {
        fprintf(stderr, "error: --snu must be non-negative\n");
        exit(EXIT_FAILURE);
    }
    if (!(o.kappa > 0.0)) {
        fprintf(stderr, "error: --kappa must be positive\n");
        exit(EXIT_FAILURE);
    }
    if (o.vt_kms < 0.0) {
        fprintf(stderr, "error: --vt must be non-negative\n");
        exit(EXIT_FAILURE);
    }
    if (o.turb_iters < 1) {
        fprintf(stderr, "error: --turb-iters must be at least 1\n");
        exit(EXIT_FAILURE);
    }
    return o;
}

static int make_dir(const char *path)
{
    if (MAKE_DIR(path) == 0) return 0;
    if (errno == EEXIST) return 0;
    fprintf(stderr, "error: cannot create directory %s\n", path);
    return -1;
}

static FILE *dat_open(const char *dir, const char *name)
{
    char path[2048];
    snprintf(path, sizeof path, "%s/%s", dir, name);
    FILE *f = fopen(path, "w");
    if (!f) fprintf(stderr, "error: cannot open %s\n", path);
    return f;
}

static int check(const char *name, int ok, const char *fmt, ...)
{
    va_list ap;
    printf("  [%s] %-36s ", ok ? "PASS" : "FAIL", name);
    va_start(ap, fmt);
    vprintf(fmt, ap);
    va_end(ap);
    printf("\n");
    return ok ? 0 : 1;
}

static int run_convergence(const Options *o)
{
    const int ns[] = { 1250, 2500, 5000, 10000, 20000 };

    printf("# Grid convergence of the isothermal model (DM only)\n");
    printf("# %8s %12s %18s\n", "jmax", "dr [kpc]", "rel.err @100 kpc");

    for (size_t k = 0; k < sizeof ns / sizeof ns[0]; ++k) {
        Cluster c;
        if (cluster_init(&c, ns[k], o->rmax_kpc, o->rVir_kpc,
                         o->r500_kpc) != 0)
            return EXIT_FAILURE;
        c.fbTarget = o->fb;
        build_grid(&c);
        build_dark_matter(&c);
        build_bcg(&c);
        build_temperature(&c);

        const double rho0 = solve_rho0(&c, MODEL_ISO);
        const int nc = c.g.n - 1;
        double *rho = calloc((size_t)nc, sizeof(double));
        int it = 0;
        integrate_model(&c, MODEL_ISO, rho0, rho, NULL, &it);

        const double A = 4.0 * PI_C * G_CGS * c.rhoNFW0 * c.rs * c.rs *
                         c.mu * MP_CGS / (KB_CGS * c.Tmg);
        int jb = 0;
        double best = 1.0e30;
        for (int j = 0; j < nc; ++j) {
            const double d = fabs(c.g.rc[j] - 100.0 * KPC_CM);
            if (d < best) { best = d; jb = j; }
        }
        const double x = c.g.rc[jb] / c.rs;
        const double f = 1.0 - log1p(x) / x;
        const double rana = rho0 * exp(-A * f);
        const double rel = fabs(rho[jb] - rana) / rana;

        printf("  %8d %12.4f %18.3e\n", ns[k],
               c.g.dr / KPC_CM, rel);

        free(rho);
        cluster_free(&c);
    }
    return EXIT_SUCCESS;
}

int main(int argc, char **argv)
{
    const Options o = parse_args(argc, argv);

    if (o.conv) return run_convergence(&o);
    if (make_dir(o.outdir) != 0) return EXIT_FAILURE;

    Cluster c;
    if (cluster_init(&c, o.jmax, o.rmax_kpc, o.rVir_kpc,
                     o.r500_kpc) != 0)
        return EXIT_FAILURE;
    c.fbTarget = o.fb;
    c.vt = o.vt_kms * 1.0e5;
    c.maxTurbIter = o.turb_iters;

    build_grid(&c);
    build_dark_matter(&c);
    build_bcg(&c);
    build_temperature(&c);

    const int nc = c.g.n - 1;
    const int n  = c.g.n;
    double *rhoIso    = calloc((size_t)nc, sizeof(double));
    double *rhoBCG    = calloc((size_t)nc, sizeof(double));
    double *rhoNonIso = calloc((size_t)nc, sizeof(double));
    double *rhoTurb   = calloc((size_t)nc, sizeof(double));
    double *pntTurb   = calloc((size_t)nc, sizeof(double));
    double *tmp       = calloc((size_t)nc, sizeof(double));
    double *tpnt      = calloc((size_t)nc, sizeof(double));
    if (!rhoIso || !rhoBCG || !rhoNonIso || !rhoTurb || !pntTurb ||
        !tmp || !tpnt) {
        fprintf(stderr, "error: out of memory\n");
        return EXIT_FAILURE;
    }

    printf("============================================================\n");
    printf(" Computational Astrophysics 2026-27 -- 1st project (C port)\n");
    printf("============================================================\n");
    printf("grid: jmax=%d  rmin=0  rmax=%.1f kpc  dr=%.4f kpc\n",
           n, c.g.rmax / KPC_CM, c.g.dr / KPC_CM);
    printf("NFW : rho0=%.3e g/cm3  rs=%.1f kpc  rvir=%.1f kpc\n",
           c.rhoNFW0, c.rs / KPC_CM, c.rVir / KPC_CM);
    printf("BCG : M*=%.2e Msun  a=%.3f kpc (r_half=12 kpc)\n",
           c.MBCG / MSUN_G, c.aBCG / KPC_CM);
    printf("gas : mu=%.2f  Tmg=%.2e K  vt=%.0f km/s\n\n",
           c.mu, c.Tmg, c.vt / 1.0e5);

    /* ---------------- Part I ---------------- */

    const char *mode = o.rho0_override > 0.0 ? "fixed" : "solved for fb=0.16";

    const double rho0_iso = (o.rho0_override > 0.0)
                            ? o.rho0_override : solve_rho0(&c, MODEL_ISO);
    int it = 0;
    integrate_model(&c, MODEL_ISO, rho0_iso, rhoIso, NULL, &it);
    const double fbIso = fb_at_rvir(&c, rhoIso);
    printf("[1] isothermal, DM only      : rho0=%.6e g/cm3 (%s), fb=%.6f\n",
           rho0_iso, mode, fbIso);

    const double Aiso = 4.0 * PI_C * G_CGS * c.rhoNFW0 * c.rs * c.rs *
                        c.mu * MP_CGS / (KB_CGS * c.Tmg);
    for (int j = 0; j < nc; ++j) {
        const double x = c.g.rc[j] / c.rs;
        const double f = 1.0 - log1p(x) / x;
        c.rhoAna[j] = rho0_iso * exp(-Aiso * f);
    }

    const double rho0_bcg = (o.rho0_override > 0.0)
                            ? o.rho0_override : solve_rho0(&c, MODEL_ISO_BCG);
    integrate_model(&c, MODEL_ISO_BCG, rho0_bcg, rhoBCG, NULL, &it);
    const double fbBCG = fb_at_rvir(&c, rhoBCG);
    printf("[2] isothermal, DM + BCG     : rho0=%.6e g/cm3, fb=%.6f\n",
           rho0_bcg, fbBCG);

    const double rho0_noniso = (o.rho0_override > 0.0)
                               ? o.rho0_override : solve_rho0(&c, MODEL_NONISO);
    integrate_model(&c, MODEL_NONISO, rho0_noniso, rhoNonIso, NULL, &it);
    const double fbNonIso = fb_at_rvir(&c, rhoNonIso);
    printf("[3] T(r), DM + BCG           : rho0=%.6e g/cm3, fb=%.6f\n",
           rho0_noniso, fbNonIso);

    const double rho0_turb = (o.rho0_override > 0.0)
                             ? o.rho0_override : solve_rho0(&c, MODEL_TURB);
    int itTurb = 0;
    integrate_model(&c, MODEL_TURB, rho0_turb, rhoTurb, pntTurb, &itTurb);
    const double fbTurb = fb_at_rvir(&c, rhoTurb);
    printf("[4] T(r) + P_nt, DM + BCG    : rho0=%.6e g/cm3, fb=%.6f "
           "(%d iterations)\n\n", rho0_turb, fbTurb, itTurb);

    /* ---------------- diagnostics ---------------- */

    const double mass_err = fabs(c.mDM[n - 1] - c.mDMana[n - 1]) /
                            c.mDMana[n - 1];

    double ana_max = 0.0;
    int j100 = 0;
    double best = 1.0e30;
    for (int j = 0; j < nc; ++j) {
        const double r = c.g.rc[j];
        const double d = fabs(r - 100.0 * KPC_CM);
        if (d < best) { best = d; j100 = j; }
        if (r > 10.0 * KPC_CM && r < 0.95 * c.g.rmax) {
            const double rel = fabs(rhoIso[j] - c.rhoAna[j]) / c.rhoAna[j];
            if (rel > ana_max) ana_max = rel;
        }
    }
    const double ana100 = fabs(rhoIso[j100] - c.rhoAna[j100]) / c.rhoAna[j100];

    double mean;
    const double hseIso  = hse_residual(&c, rhoIso, NULL, c.mDM, NULL, &mean);
    const double hseBCG  = hse_residual(&c, rhoBCG, NULL, c.mTot, NULL, &mean);
    const double hseNon  = hse_residual(&c, rhoNonIso, NULL, c.mTot, c.T, &mean);
    const double hseTurb = hse_residual(&c, rhoTurb, pntTurb, c.mTot, c.T,
                                        &mean);

    printf("Checks / diagnostics\n");
    printf("  M_DM(%.2f Mpc): num=%.6e  analytic=%.6e  rel.err=%.2e Msun\n",
           c.rVir / KPC_CM, c.mDM[n - 1] / MSUN_G,
           c.mDMana[n - 1] / MSUN_G, mass_err);
    printf("  isothermal rho vs analytic: max rel.err (r>10 kpc)=%.2e, "
           "at %.1f kpc=%.2e\n", ana_max, c.g.rc[j100] / KPC_CM, ana100);
    printf("  HSE residual: iso=%.2e  BCG=%.2e  T(r)=%.2e  turb=%.2e\n",
           hseIso, hseBCG, hseNon, hseTurb);

    printf("\nTurbulent fixed-point convergence (P_nt from previous iterate)\n");
    double rel3 = 0.0;
    for (int k = 1; k <= 5; ++k) {
        c.maxTurbIter = k;
        int itk = 0;
        integrate_model(&c, MODEL_TURB, rho0_turb, tmp, tpnt, &itk);
        double rel = 0.0;
        for (int j = 0; j < nc; ++j) {
            const double d = fabs(tmp[j] - rhoTurb[j]) / rhoTurb[j];
            if (d > rel) rel = d;
        }
        printf("  %d iteration(s): max rel. change vs converged = %.3e\n",
               k, rel);
        if (k == 3) rel3 = rel;
    }
    c.maxTurbIter = o.turb_iters;

    const double turbCenter = rhoTurb[0] / rhoNonIso[0];
    const double turbOuter  = rhoTurb[nc - 1] / rhoNonIso[nc - 1];
    printf("  turbulence effect: rho(centre) x%.4f, rho(outer) x%.4f\n\n",
           turbCenter, turbOuter);

    /* ---------------- write Part I tables ---------------- */

    FILE *f = dat_open(o.outdir, "grid.dat");
    if (!f) return EXIT_FAILURE;
    fprintf(f, "# r_left[kpc] r_center[kpc]\n");
    for (int j = 0; j < nc; ++j)
        fprintf(f, "%.8e %.8e\n", c.g.r[j] / KPC_CM, c.g.rc[j] / KPC_CM);
    fclose(f);

    f = dat_open(o.outdir, "mass.dat");
    if (!f) return EXIT_FAILURE;
    fprintf(f, "# r[kpc] M_DM_num[Msun] M_DM_analytic[Msun]\n");
    for (int j = 0; j < n; ++j)
        fprintf(f, "%.8e %.8e %.8e\n", c.g.r[j] / KPC_CM,
                c.mDM[j] / MSUN_G, c.mDMana[j] / MSUN_G);
    fclose(f);

    f = dat_open(o.outdir, "gas.dat");
    if (!f) return EXIT_FAILURE;
    fprintf(f, "# r_center[kpc] rho_iso rho_analytic rho_NFW [g/cm3]\n");
    for (int j = 0; j < nc; ++j)
        fprintf(f, "%.8e %.8e %.8e %.8e\n", c.g.rc[j] / KPC_CM,
                rhoIso[j], c.rhoAna[j], c.rhoNFW[j]);
    fclose(f);

    f = dat_open(o.outdir, "gas_BCG.dat");
    if (!f) return EXIT_FAILURE;
    fprintf(f, "# r_center[kpc] rho_iso_BCG rho_NFW [g/cm3]\n");
    for (int j = 0; j < nc; ++j)
        fprintf(f, "%.8e %.8e %.8e\n", c.g.rc[j] / KPC_CM,
                rhoBCG[j], c.rhoNFW[j]);
    fclose(f);

    f = dat_open(o.outdir, "gas_BCG_non-iso.dat");
    if (!f) return EXIT_FAILURE;
    fprintf(f, "# r_center[kpc] rho_noniso rho_NFW [g/cm3]\n");
    for (int j = 0; j < nc; ++j)
        fprintf(f, "%.8e %.8e %.8e\n", c.g.rc[j] / KPC_CM,
                rhoNonIso[j], c.rhoNFW[j]);
    fclose(f);

    f = dat_open(o.outdir, "temperature.dat");
    if (!f) return EXIT_FAILURE;
    fprintf(f, "# r_center[kpc] T[K] T/Tmg\n");
    for (int j = 0; j < nc; ++j)
        fprintf(f, "%.8e %.8e %.8e\n", c.g.rc[j] / KPC_CM,
                c.T[j], c.T[j] / c.Tmg);
    fclose(f);

    f = dat_open(o.outdir, "gas_BCG_non-iso-turb.dat");
    if (!f) return EXIT_FAILURE;
    fprintf(f, "# r_center[kpc] rho_turb rho_noniso P_nt P_nt/P_th [cgs]\n");
    for (int j = 0; j < nc; ++j) {
        const double pth = rhoTurb[j] * KB_CGS * c.T[j] / (c.mu * MP_CGS);
        fprintf(f, "%.8e %.8e %.8e %.8e %.8e\n", c.g.rc[j] / KPC_CM,
                rhoTurb[j], rhoNonIso[j], pntTurb[j], pntTurb[j] / pth);
    }
    fclose(f);

    /* ---------------- Part II ---------------- */

    FeSnapshot sA[4], sB[2], sC[2];
    memset(sA, 0, sizeof sA);
    memset(sB, 0, sizeof sB);
    memset(sC, 0, sizeof sC);

    const double tA[4] = { 0.0, 1.0 * GYR_S, 2.0 * GYR_S, 5.0 * GYR_S };
    const double tB[2] = { 1.0 * GYR_S, 5.0 * GYR_S };
    double drift = 0.0, srcRel = 0.0;

    if (o.fe) {
        FeParams fp;
        fe_default_params(&fp);
        fp.kappa0 = o.kappa;
        fp.snu = o.snu;
        fp.vary = o.kappa_vary;
        fp.alphaEvolve = o.alpha_evolve;

        char path[512];

        f = dat_open(o.outdir, "fe_observed.dat");
        if (!f) return EXIT_FAILURE;
        fprintf(f, "# r_center[kpc] Z_obs/Zsun Z_exc/Zsun\n");
        for (int j = 0; j < nc; ++j)
            fprintf(f, "%.8e %.8e %.8e\n", c.g.rc[j] / KPC_CM,
                    fe_observed_zfe(&c, c.g.rc[j], 0) / ZSUN,
                    fe_observed_zfe(&c, c.g.rc[j], 1) / ZSUN);
        fclose(f);

        printf("Part II: iron diffusion (D=%.3e cm2/s, SNu=%.2f)\n",
               fp.kappa0, fp.snu);

        snprintf(path, sizeof path, "%s/fe_diffusion.dat", o.outdir);
        if (fe_run(&c, rhoNonIso, 0, 1, FE_INIT_OBSERVED, &fp,
                   tA, 4, sA, path) != 0)
            return EXIT_FAILURE;
        if (sA[0].mfetot > 0.0)
            drift = fabs(sA[3].mfetot - sA[0].mfetot) / sA[0].mfetot;

        printf("\n");
        snprintf(path, sizeof path, "%s/fe_source.dat", o.outdir);
        if (fe_run(&c, rhoNonIso, 1, 0, FE_INIT_ZERO, &fp,
                   tB, 2, sB, path) != 0)
            return EXIT_FAILURE;
        const double theory = fe_source_total_theory(&c, fp.snu, 5.0 * GYR_S);
        if (theory > 0.0)
            srcRel = fabs(sB[1].mfetot - theory) / theory;

        printf("\n");
        snprintf(path, sizeof path, "%s/fe_both.dat", o.outdir);
        if (fe_run(&c, rhoNonIso, 1, 1, FE_INIT_ZERO, &fp,
                   tB, 2, sC, path) != 0)
            return EXIT_FAILURE;

        f = dat_open(o.outdir, "fe_masses.dat");
        if (!f) return EXIT_FAILURE;
        fprintf(f, "# t[Gyr] experiment M_Fe(<100kpc)[Msun] M_Fe(tot)[Msun]\n");
        for (int k = 0; k < 4; ++k)
            fprintf(f, "%.6f diffusion_only %.8e %.8e\n",
                    sA[k].t / GYR_S, sA[k].mfe100, sA[k].mfetot);
        for (int k = 0; k < 2; ++k)
            fprintf(f, "%.6f source_only %.8e %.8e\n",
                    sB[k].t / GYR_S, sB[k].mfe100, sB[k].mfetot);
        for (int k = 0; k < 2; ++k)
            fprintf(f, "%.6f source_diffusion %.8e %.8e\n",
                    sC[k].t / GYR_S, sC[k].mfe100, sC[k].mfetot);
        fclose(f);
    }

    /* ---------------- checks ---------------- */

    int fails = 0;
    if (o.check) {
        printf("\nSelf-tests\n");
        fails += check("DM mass vs analytic", mass_err < 1.0e-4,
                       "max rel.err=%.2e", mass_err);
        fails += check("isothermal vs analytic", ana_max < 5.0e-3,
                       "max rel.err=%.2e", ana_max);
        if (o.rho0_override > 0.0) {
            printf("  [SKIP] fb targets: --rho0 was given\n");
        } else {
            fails += check("fb(iso)=0.16", fabs(fbIso - 0.16) < 1.0e-8,
                           "fb=%.10f", fbIso);
            fails += check("fb(BCG)=0.16", fabs(fbBCG - 0.16) < 1.0e-8,
                           "fb=%.10f", fbBCG);
            fails += check("fb(noniso)=0.16", fabs(fbNonIso - 0.16) < 1.0e-8,
                           "fb=%.10f", fbNonIso);
            fails += check("fb(turb)=0.16", fabs(fbTurb - 0.16) < 1.0e-8,
                           "fb=%.10f", fbTurb);
        }
        fails += check("HSE residual, T(r)", hseNon < 5.0e-2,
                       "max=%.2e", hseNon);
        fails += check("HSE residual, turbulence", hseTurb < 5.0e-2,
                       "max=%.2e", hseTurb);
        fails += check("turbulence: 3 iterations enough", rel3 < 1.0e-2,
                       "rel.change=%.2e", rel3);
        fails += check("turbulence flattens rho", turbCenter < 1.0 &&
                       turbOuter > 1.0, "centre x%.4f, outer x%.4f",
                       turbCenter, turbOuter);
        if (o.rho0_override <= 0.0)
            fails += check("BCG raises rho0", rho0_bcg > rho0_iso,
                           "%.3e -> %.3e", rho0_iso, rho0_bcg);
        if (o.fe) {
            if (o.alpha_evolve)
                printf("  [SKIP] Fe source total vs analytic: alpha(t) is on\n");
            else
                fails += check("Fe source total vs analytic", srcRel < 1.0e-5,
                               "rel.err=%.2e", srcRel);
            fails += check("Fe mass conserved (diffusion)", drift < 1.0e-10,
                           "drift=%.2e", drift);
            fails += check("Fe initial M(<100kpc) > 0",
                           sA[0].mfe100 > 0.0, "%.3e Msun", sA[0].mfe100);
            fails += check("Fe source+diffusion < source only",
                           sC[1].mfe100 < sB[1].mfe100,
                           "%.3e < %.3e", sC[1].mfe100, sB[1].mfe100);
        }
    }

    /* ---------------- summary ---------------- */

    f = dat_open(o.outdir, "summary.txt");
    if (!f) return EXIT_FAILURE;
    fprintf(f, "Project 1 (C) -- summary\n");
    fprintf(f, "grid: jmax=%d rmax=%.2f kpc dr=%.5f kpc\n",
            n, c.g.rmax / KPC_CM, c.g.dr / KPC_CM);
    fprintf(f, "rho0_iso      = %.10e g/cm3  fb=%.8f\n", rho0_iso, fbIso);
    fprintf(f, "rho0_BCG      = %.10e g/cm3  fb=%.8f\n", rho0_bcg, fbBCG);
    fprintf(f, "rho0_noniso   = %.10e g/cm3  fb=%.8f\n", rho0_noniso, fbNonIso);
    fprintf(f, "rho0_turb     = %.10e g/cm3  fb=%.8f\n", rho0_turb, fbTurb);
    fprintf(f, "M_DM(rmax)    = %.6e Msun (analytic %.6e, rel.err %.2e)\n",
            c.mDM[n - 1] / MSUN_G, c.mDMana[n - 1] / MSUN_G, mass_err);
    fprintf(f, "isothermal analytic max rel.err (r>10 kpc) = %.3e\n", ana_max);
    fprintf(f, "HSE residuals: iso=%.2e BCG=%.2e T(r)=%.2e turb=%.2e\n",
            hseIso, hseBCG, hseNon, hseTurb);
    fprintf(f, "turbulence: centre x%.4f outer x%.4f, %d iterations\n",
            turbCenter, turbOuter, itTurb);
    if (o.fe) {
        fprintf(f, "Fe diffusion only: M(<100kpc) %.3e -> %.3e -> %.3e Msun "
                   "(t=1,2,5 Gyr)\n", sA[1].mfe100, sA[2].mfe100,
                sA[3].mfe100);
        fprintf(f, "Fe source only   : M(<100kpc) %.3e -> %.3e Msun "
                   "(t=1,5 Gyr), total %.3e -> %.3e\n",
                sB[0].mfe100, sB[1].mfe100, sB[0].mfetot, sB[1].mfetot);
        fprintf(f, "Fe source+diff   : M(<100kpc) %.3e -> %.3e Msun "
                   "(t=1,5 Gyr), total %.3e -> %.3e\n",
                sC[0].mfe100, sC[1].mfe100, sC[0].mfetot, sC[1].mfetot);
        fprintf(f, "Fe mass conservation (diffusion only) = %.3e\n", drift);
    }
    fclose(f);

    free(rhoIso);
    free(rhoBCG);
    free(rhoNonIso);
    free(rhoTurb);
    free(pntTurb);
    free(tmp);
    free(tpnt);
    cluster_free(&c);

    if (o.check) {
        printf("\nSelf-tests: %s\n", fails ? "FAILED" : "all passed");
        return fails ? EXIT_FAILURE : EXIT_SUCCESS;
    }
    printf("\nWrote tables to %s/\n", o.outdir);
    return EXIT_SUCCESS;
}
