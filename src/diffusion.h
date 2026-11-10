/*
 * diffusion.h -- Part II: turbulent diffusion of iron in the ICM.
 */
#ifndef DIFFUSION_H
#define DIFFUSION_H

#include "cluster.h"

typedef struct {
    double kappa0;      /* diffusion coefficient at r -> infinity, cm^2/s  */
    int    vary;        /* use kappa(r) = kappa0 (1 - coreFactor exp(-(r/rs)^2)) */
    double rscala;      /* cm                                             */
    double coreFactor;  /* 0.6 in the lecture's example                   */
    double snu;         /* SNIa rate in SNu (0.15 fiducial)               */
    int    alphaEvolve; /* use alpha(t) instead of alpha(t_now)           */
} FeParams;

typedef struct {
    double t;           /* s                                              */
    double mfe100;      /* Fe mass inside 100 kpc, Msun                   */
    double mfetot;      /* Fe mass over the whole grid, Msun              */
} FeSnapshot;

typedef enum { FE_INIT_ZERO = 0, FE_INIT_OBSERVED = 1 } FeInit;

void fe_default_params(FeParams *p);

/* Evolve rho_Fe on the same shell grid as Part I.
 *   useSource    : add S_Fe(r,t)
 *   useDiffusion : add the divergence term
 *   init         : Z_Fe = 0 or the observed Perseus profile minus Z_out
 * snapTimes[] is in seconds.  The profile table goes to path (may be NULL).
 * Returns 0 on success, -1 on a numerical failure. */
int fe_run(const Cluster *c, const double *rhoGas,
           int useSource, int useDiffusion, FeInit init,
           const FeParams *par,
           const double *snapTimes, int nsnap,
           FeSnapshot *snaps, const char *path);

/* Observed Perseus fit (with the lecture's 1.4 factor); if excess != 0 the
 * outer plateau Z_out = 0.4 Z_sun is subtracted and the result clipped at 0.
 * Returns Z_Fe (not divided by Z_sun). */
double fe_observed_zfe(const Cluster *c, double r, int excess);

/* Analytic expectation for the source-only run: alpha_tot * M_BCG(<rmax) * t. */
double fe_source_total_theory(const Cluster *c, double snu, double t);

#endif /* DIFFUSION_H */
