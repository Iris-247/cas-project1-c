/*
 * cluster.h -- Part I of the "hydrostatic equilibrium in a galaxy cluster"
 * project, implemented in C.
 *
 * All quantities are in cgs unless a name says otherwise (the helper
 * conversion factors KPC_CM / MSUN_G are used at the boundaries).
 */
#ifndef CLUSTER_H
#define CLUSTER_H

#include <stddef.h>

#define KPC_CM   3.084e21
#define YEAR_S   3.15576e7
#define GYR_S    (1.0e9 * YEAR_S)
#define G_CGS    6.67430e-8
#define KB_CGS   1.380649e-16
#define MP_CGS   1.6726219e-24
#define MSUN_G   1.9885e33
#define MU_GAS   0.60
#define ZSUN     1.8e-3
#define PI_C     3.14159265358979323846

/* Radial grid: n boundaries, n-1 shell centres. */
typedef struct {
    int     n;          /* number of grid boundaries (jmax in the Fortran)   */
    double  rmin;       /* cm                                                */
    double  rmax;       /* cm                                                */
    double  dr;         /* cm (uniform)                                      */
    double *r;          /* n   boundary radii                                */
    double *rc;         /* n-1 shell-centre radii                            */
    double *vol;        /* n-1 shell volumes 4pi/3 (r[j+1]^3 - r[j]^3)       */
} Grid;

typedef enum {
    MODEL_ISO = 0,      /* isothermal, DM only                                */
    MODEL_ISO_BCG,      /* isothermal, DM + BCG                               */
    MODEL_NONISO,       /* T(r), DM + BCG                                     */
    MODEL_TURB          /* T(r) + non-thermal pressure, DM + BCG              */
} ModelKind;

typedef struct {
    Grid g;

    /* model parameters (cgs) */
    double rs;          /* NFW scale radius                                   */
    double rhoNFW0;     /* NFW characteristic density                         */
    double MBCG;        /* total BCG stellar mass                             */
    double aBCG;        /* Hernquist scale length (12 kpc half-mass -> a)     */
    double rVir;        /* virial radius                                      */
    double r500;        /* R500 = Rvir/2                                      */
    double mu;          /* mean molecular weight                              */
    double Tmg;         /* isothermal / mean temperature                       */
    double vt;          /* turbulent velocity                                 */
    double fbTarget;    /* target baryon fraction (0.16)                      */
    double rho0;        /* solved central gas density                         */
    int    maxTurbIter; /* cap for the turbulent fixed-point iteration        */

    /* profiles; "n" arrays are defined on boundaries, "n-1" at centres */
    double *rhoNFW;     /* n-1  NFW density at shell centres                   */
    double *mDM;        /* n    numerically integrated DM mass                 */
    double *mDMana;     /* n    analytic NFW mass                              */
    double *mBCG;       /* n    Hernquist BCG mass                             */
    double *mTot;       /* n    mDM + mBCG                                     */
    double *T;          /* n-1  temperature profile                            */
    double *rho;        /* n-1  scratch / current gas density                  */
    double *rhoAna;     /* n-1  analytic isothermal density                    */
    double *pnt;        /* n-1  turbulent pressure                             */
    double *scratch;    /* n-1  work array for the turbulent iteration         */

    int jVir;           /* grid boundary nearest to rVir                       */
} Cluster;

int  cluster_init(Cluster *c, int n, double rmax_kpc,
                  double rVir_kpc, double r500_kpc);
void cluster_free(Cluster *c);

void build_grid(Cluster *c);
void build_dark_matter(Cluster *c);
void build_bcg(Cluster *c);
void build_temperature(Cluster *c);

double nfw_mass_analytic(const Cluster *c, double r);
double temperature_at(const Cluster *c, double r);

/* Integrate the (non-isothermal) hydrostatic equation outward.  rho0 is the
 * central density; rho[] receives the shell-centre densities.  For MODEL_TURB
 * pnt[] receives P_nt and *iters the number of fixed-point iterations used. */
void integrate_model(Cluster *c, ModelKind kind, double rho0,
                     double *rho, double *pnt, int *iters);

void   compute_gas_mass(const Cluster *c, const double *rho, double *mgas);
double fb_at_rvir(const Cluster *c, const double *rho);
double solve_rho0(Cluster *c, ModelKind kind);

/* max |dP/dr + rho G M/r^2| / |rho G M/r^2| over interior shells;
 * mean is optional.  pnt and Tarr may be NULL, mass is mDM or mTot. */
double hse_residual(const Cluster *c, const double *rho, const double *pnt,
                    const double *mass, const double *Tarr, double *mean);

#endif /* CLUSTER_H */
