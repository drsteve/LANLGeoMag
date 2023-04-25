#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <Lgm/Lgm_DynamicMemory.h>
#include <Lgm/Lgm_ElapsedTime.h>
#include <Lgm/Lgm_CTrans.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_statistics_double.h> 

#define  ZeroNegativeEigenvalues 1
#define  ReverseNegativeEigenvalues 0

#define VERBOSE 0

typedef struct GprInfo {


    int     ForceSymmetry;


    double  sigma_f, l;
    double  log_sigma_f, log_l;

    double  sigma_b, sigma_nu;
    double  log_sigma_b, log_sigma_nu;
    double  c;

    int         n;
    gsl_vector *x;
    gsl_vector *y;
    gsl_vector *sigma_n_vec;
    gsl_vector *sigma_n_2_vec;
    double      sigma_n;
    double      sigma_n_2;

    gsl_matrix *K;
    gsl_matrix *K_dl;
    gsl_matrix *K_dsig;
    gsl_matrix *K_dsigb;
    gsl_matrix *K_dsignu;
    gsl_matrix *K_dc;

    gsl_matrix *K_inv;
    gsl_vector *r;
    gsl_matrix *A;
    gsl_matrix *C;

    //double logp, dlogp_dsig, dlogp_dl;
    double logp;
    double dlogp_dlog_sig;
    double dlogp_dlog_l;
    double dlogp_dlog_sigb;
    double dlogp_dlog_signu;
    double dlogp_dc;


    /*
     * These hold the final mean and credibility intervals
     */
    int         n_star;
    gsl_vector *x_star;
    gsl_vector *y_hat;
    gsl_vector *y_cred_lo;
    gsl_vector *y_cred_hi;

} GprInfo;

GprInfo *InitGprInfo( int, int, int );
void FreeGprInfo( GprInfo * );
void ComputeEvidenceAndDerivs( double, double, double, double, double, GprInfo * );
void GPR( GprInfo * );

