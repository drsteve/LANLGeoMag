#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "DynamicMemory.h"
#include "Lgm/Lgm_ElapsedTime.h"
#include "Lgm/Lgm_CTrans.h"
#include "Lgm/Lgm_FluxToPsd.h"



#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_statistics_double.h> 

#define  ZeroNegativeEigenvalues 1
#define  ReverseNegativeEigenvalues 0

typedef struct GprInfo {

    double  sigma_f, l;
    double  log_sigma_f, log_l;

    gsl_vector *x;
    gsl_vector *y;
    gsl_vector *sigma_n_vec;
    gsl_vector *sigma_n_2_vec;
    double      sigma_n;
    double      sigma_n_2;

    gsl_matrix *K;
    gsl_matrix *K_dl;
    gsl_matrix *K_dsig;
    gsl_matrix *K_inv;
    gsl_vector *r;
    gsl_matrix *A;
    gsl_matrix *C;

    //double logp, dlogp_dsig, dlogp_dl;
    double logp, dlogp_dlog_sig, dlogp_dlog_l;


    /*
     * These hold the final mean and credibility intervals
     */
    gsl_vector *y_hat;
    gsl_vector *y_cred_lo;
    gsl_vector *y_cred_hi;

} GprInfo;








unsigned char Ctab_Red[] = { 0, 76, 78, 79, 80, 81, 83, 84, 85, 86, 88, 89, 90, 92, 93, 94, 95, 97, 98, 99, 100, 102, 103, 104, 106, 107, 108, 109, 111, 112, 113, 114, 116, 117, 118, 119, 121, 122, 123, 125, 126, 127, 128, 130, 131, 132, 133, 135, 136, 137, 139, 140, 141, 142, 144, 145, 146, 147, 149, 150, 151, 152, 154, 155, 156, 158, 159, 160, 161, 163, 164, 165, 166, 168, 169, 170, 172, 173, 174, 175, 177, 178, 179, 180, 182, 183, 184, 186, 187, 188, 189, 191, 192, 193, 194, 196, 197, 198, 199, 201, 202, 203, 205, 206, 207, 208, 210, 211, 212, 213, 215, 216, 217, 219, 220, 221, 222, 224, 225, 226, 227, 229, 230, 231, 232, 234, 235, 236, 238, 238, 238, 238, 238, 238, 238, 239, 239, 239, 239, 239, 239, 239, 240, 240, 240, 240, 240, 240, 240, 241, 241, 241, 241, 241, 241, 241, 242, 242, 242, 242, 242, 242, 242, 243, 243, 243, 243, 243, 243, 243, 243, 244, 244, 244, 244, 244, 244, 244, 245, 245, 245, 245, 245, 245, 245, 246, 246, 246, 246, 246, 246, 246, 247, 247, 247, 247, 247, 247, 247, 248, 248, 248, 248, 248, 248, 248, 249, 249, 249, 249, 249, 249, 249, 250, 250, 250, 250, 250, 250, 250, 251, 251, 251, 251, 251, 251, 251, 252, 252, 252, 252, 252, 252, 252, 253, 253, 253, 252, 251, 250, 250, 249, 248, 247, 246, 245, 245, 244, 243, 242, 241, 240, 240, 239, 238, 237 };
unsigned char Ctab_Grn[] = { 0, 23, 23, 24, 24, 25, 25, 26, 26, 27, 27, 28, 29, 29, 30, 30, 31, 32, 32, 33, 33, 34, 35, 35, 36, 37, 37, 38, 38, 39, 40, 40, 41, 42, 43, 43, 44, 45, 45, 46, 47, 47, 48, 49, 50, 50, 51, 52, 53, 53, 54, 55, 56, 56, 57, 58, 59, 60, 60, 61, 62, 63, 64, 64, 65, 66, 67, 68, 69, 69, 70, 71, 72, 73, 74, 75, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 124, 125, 126, 127, 128, 129, 130, 132, 133, 134, 135, 136, 137, 139, 140, 141, 142, 143, 144, 146, 147, 148, 149, 150, 151, 153, 154, 155, 156, 157, 158, 160, 161, 162, 163, 164, 165, 167, 168, 169, 170, 171, 172, 174, 175, 176, 177, 178, 179, 181, 182, 183, 184, 185, 186, 188, 189, 190, 191, 192, 193, 194, 196, 197, 198, 199, 200, 201, 203, 204, 205, 206, 207, 208, 209, 211, 212, 213, 214, 215, 216, 217, 219, 220, 221, 222, 223, 224, 225, 227, 228, 229, 230, 231, 232, 233, 235, 236, 237, 238, 239, 240, 241, 242, 244, 245, 246, 247, 248, 249, 250, 251, 253, 253, 253, 253, 253, 254, 254, 254, 254, 254, 254, 254, 255, 255, 255, 255, 255, 255, 255, 255 };
unsigned char Ctab_Blu[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 10, 10, 10, 10, 10, 11, 11, 11, 11, 11, 12, 12, 12, 12, 12, 13, 13, 13, 13, 13, 14, 14, 14, 14, 14, 15, 15, 15, 15, 15, 16, 16, 16, 16, 16, 17, 17, 17, 17, 17, 18, 18, 18, 18, 19, 19, 19, 19, 19, 20, 20, 20, 20, 20, 21, 21, 21, 21, 21, 22, 22, 22, 22, 23, 23, 23, 23, 23, 24, 24, 24, 24, 24, 25, 25, 25 };
/**
 *   Routine to write out a GIF image
 */
void DumpGif3( char *FilenameBase, int W, int H, double **Image ){

    double  Val, Val2, Min, Max, dVal;
    int     w, h;
    char    Filename[1024];
    unsigned char *uImage, uVal;
    FILE    *fp_gif, *fp_info;

    int     LogScale;

    LogScale = TRUE;
    LogScale = FALSE;


    // Determine Min/Max values...
    Min =  9e99;
    Max = -9e99;
    for ( w=0; w<W; w++ ){
        for ( h=0; h<H; h++ ) {

            Val = Image[h][w];
            if ( LogScale ) {
                Val2 = (Val > 0.0) ? log10( Val ) : -9e99;
            } else {
                Val2 = Image[h][w];
            }
            if (Val2 > Max) Max = Val2;
            if (Val2 < Min) Min = Val2;
            //if ((Val2 < Min)&&(Val > 0.0)) Min = Val2;

        }
    }

    printf("Min, Max = %g %g\n", Min, Max);
//Min = 1.9;
//Max = 3.5;

/*
Min = -3.0;
Max = 1.0;
*/

    sprintf( Filename, "%s.info", FilenameBase);
    fp_info = fopen( Filename, "w" );
    fprintf( fp_info, "Min: %g\n", Min );
    fprintf( fp_info, "Max: %g\n", Max );
    fclose( fp_info );



    uImage = (unsigned char *)calloc( W*H, sizeof(unsigned char) );

    for ( w=0; w<W; w++ ){
        for ( h=0; h<H; h++ ) {

            if ( LogScale ) {
                Val = Image[h][w] > 0.0 ? log10( Image[h][w] ) : -1e31;
            } else {
                Val = Image[h][w];
            }
            if ( Val < Min ) {
                uVal = 0;
            } else {
                dVal = (Val - Min)/(Max-Min)*255.0;
                uVal = (dVal > 255.0) ? 255 : (unsigned char)dVal;
            }

            *(uImage + W*(h) + W-1-w) = uVal;

        }
    }

    sprintf( Filename, "%s.gif", FilenameBase);
    fp_gif = fopen(Filename, "w");
    WriteGIF(fp_gif, (byte *)uImage, 0, W, H, Ctab_Red, Ctab_Grn, Ctab_Blu, 256, 0, "");
    fclose(fp_gif);

    free( uImage );



    // dump a colorbar image
if (0==1){
    W = 10; H = 256;
    uImage = (unsigned char *)calloc( W*H, sizeof(unsigned char) );
    for ( w=0; w<W; w++ ){
        for ( h=0; h<H; h++ ) {
            *(uImage + W*(H-1-h) + w) = h;
        }
    }
    sprintf( Filename, "%s_Bar.gif", FilenameBase);
    fp_gif = fopen(Filename, "w");
    WriteGIF(fp_gif, (byte *)uImage, 0, W, H, Ctab_Red, Ctab_Grn, Ctab_Blu, 256, 0, "");
    fclose(fp_gif);
    free( uImage );
}

}


double kernel_function( double x, double y, double sigma_f, double l ) {

    double  diff, Norm2;

    /*
     * Compute the norm for x-y. norm = sqrt( sum( (x[i]-y[i])^2 ) ^1/2
     * since we end up squaring the norm later dont need sqrt.
     */
    diff = x - y;
    Norm2 = diff*diff;
    // Norm = sqrt(Norm2);

    return( sigma_f * exp( -Norm2 / (2.0*l*l) ) );

}

double kernel_function_dk_dsig( double x, double y, double sigma_f, double l ) {

    double  diff, Norm2;

    /*
     * Compute the norm for x-y. norm = sqrt( sum( (x[i]-y[i])^2 ) ^1/2
     * since we end up squaring the norm later dont need sqrt.
     */
    diff = x - y;
    Norm2 = diff*diff;
    // Norm = sqrt(Norm2);

    return( exp( -Norm2 / (2.0*l*l) ) );

}

double kernel_function_dk_dl( double x, double y, double sigma_f, double l ) {

    double  diff, Norm2;
    double  l2;

    /*
     * Compute the norm for x-y. norm = sqrt( sum( (x[i]-y[i])^2 ) ^1/2
     * since we end up squaring the norm later dont need sqrt.
     */
    diff = x - y;
    Norm2 = diff*diff;
    l2 = l*l;
    // Norm = sqrt(Norm2);

    return( Norm2*sigma_f/(l*l2) * exp( -Norm2 / (2.0*l2) ) );

}


/*
 * This is a composite kernel function k = k(x, x') + k(x',x)
 * which should enforce symmetry
 */
double symm_kernel_function( double x, double y, double sigma_f, double l ) {
    return( kernel_function( x, y, sigma_f, l ) + kernel_function( -x, y, sigma_f, l ) );
}
double symm_kernel_function_dk_dsig( double x, double y, double sigma_f, double l ) {
    return( kernel_function_dk_dsig( x, y, sigma_f, l ) + kernel_function_dk_dsig( -x, y, sigma_f, l ) );
}
double symm_kernel_function_dk_dl( double x, double y, double sigma_f, double l ) {
    return( kernel_function_dk_dl( x, y, sigma_f, l ) + kernel_function_dk_dl( -x, y, sigma_f, l ) );
}




double InitGprInfo( int n, GprInfo *Info ) {


    Info->x      = gsl_vector_calloc( n );
    Info->y      = gsl_vector_calloc( n );
    Info->sigma_n_vec = gsl_vector_calloc( n );
    Info->sigma_n_2_vec = gsl_vector_calloc( n );

    Info->K      = gsl_matrix_calloc( n, n );
    Info->K_dl   = gsl_matrix_calloc( n, n );
    Info->K_dsig = gsl_matrix_calloc( n, n );
    Info->K_inv  = gsl_matrix_calloc( n, n );
    Info->r      = gsl_vector_calloc( n );
    Info->A      = gsl_matrix_calloc( n, n );
    Info->C      = gsl_matrix_calloc( n, n );

}

/*
 * This is written up in Gaussian Processes for Machine Learning Carl Edward
 * Rasmussen and Christopher K. I. Williams MIT Press, 2006. ISBN-10
 * 0-262-18253-X, ISBN-13 978-0-262-18253-9. 
 * 
 * See equation 5.8. Also see equation 2.30 and the appnedix for some math
 * machinery.
 *
 * From eqn 5.8, the so-called "marginal liklehood" or "evidence" is given by:
 *
 *  log ( p(y|X,Theta) ) = -1/2 y^T K_y^-1 y - 1/2 log |K_y| - n/2 log(2 pi)
 *
 *  The matrix K_y = K_f + sigma_n^2 I
 *  and the quantity Theta refers to the hyperparameters.
 *
 *  To facilitate faster optimization, this routine computes the evidence and
 *  its derivative wrt to the 2 hyperparameter values.
 *  From eqn 5.9, we have:
 *
 *      \partial log( p(y|X,Theta) )\partial\theta_j 
 *                 = 1/2 Tr( (alpha alpha^T - K^-1)\partial K\partial\theta_j )
 *
 *  where alpha = K^-1 y
 *
 *   We need to provide a better way of including the uncertainties. Here we just assume sigma_n
 *
 */
//double ComputeEvidenceAndDerivs( gsl_vector *x, gsl_vector *y, double sigma_n, double sigma_f, double l ) {
void ComputeEvidenceAndDerivs( double log_sigma_f, double log_l, GprInfo *Info ) {


    int     n, i, j, Method;
    double  log_K_det, Trace, E1, E2, E3, val, xi, xj;
    double  sig, sigma_f, l;
    gsl_vector *x, *y;

    x = Info->x;
    y = Info->y;

    
    sigma_f = exp( log_sigma_f );
    l       = exp( log_l );

    /*
     * Compute the following matrix quantities: K, K_dl, K_dsig.
     * The last two are the partial derivs: \partial K\partial\theta_j
     */
    n = x->size;

//    gsl_matrix *K      = gsl_matrix_calloc( n, n );
//    gsl_matrix *K_dl   = gsl_matrix_calloc( n, n );
//    gsl_matrix *K_dsig = gsl_matrix_calloc( n, n );
    for ( i=0; i<n; i++ ) {
        xi = gsl_vector_get( x, i );
        for ( j=0; j<n; j++ ) {
            xj = gsl_vector_get( x, j);
            if ( i==j ) {
                sig = gsl_vector_get( Info->sigma_n_2_vec, i );
                gsl_matrix_set( Info->K, i, j, symm_kernel_function( xi, xj, sigma_f, l ) + sig );
            } else {
                gsl_matrix_set( Info->K, i, j, symm_kernel_function( xi, xj, sigma_f, l ) );
            }
            gsl_matrix_set( Info->K_dl,   i, j,   symm_kernel_function_dk_dl( xi, xj, sigma_f, l ) );
            gsl_matrix_set( Info->K_dsig, i, j, symm_kernel_function_dk_dsig( xi, xj, sigma_f, l ) );
        }
    }
   

    /*
     *  Compute det(K) and K^-1
     *  From a Cholesky decomp, det(K) = sum of squared diagonals.
     */
    Method = ReverseNegativeEigenvalues;
Method = ZeroNegativeEigenvalues;
//    gsl_matrix *K_inv = gsl_matrix_calloc( n, n );
    CholeskyLikeDecomp( Info->K, Info->K_inv, Method );
    for ( log_K_det = 0.0, i=0; i<n; i++ ) {
        // K_inv is not the inverse (yet) 
        // -- its the L in the Cholesky decomp (K = L L^T)
        val = gsl_matrix_get( Info->K_inv, i, i ); 
        log_K_det += log(val);
    }
    //log_K_det *= 2.0;
    gsl_linalg_cholesky_invert( Info->K_inv );


    /*
     *  Compute -1/2 y^T K_inv y. First do r = K_inv y, then do -1/2 y^T r.
     */
//    gsl_vector *r = gsl_vector_calloc( n );
    gsl_blas_dgemv( CblasNoTrans, 1.0, Info->K_inv, Info->y, 0.0, Info->r );
    gsl_blas_ddot( y, Info->r, &E1 ); 
    E1 *= -0.5;
    
    /*
     *  Compute -1/2 log|K|
     */
    //E2 = -0.5*log(K_det);
    E2 = -log_K_det;


    /*
     *  Compute -n/2 log(2 pi)
     */
    E3 = (double)n * -0.5 * log(2.0*M_PI);

    
    /*
     *  Compute  log( p(y|X,Theta) ) (log evidence)
     */
    Info->logp = E1 + E2 + E3;
    printf("sigma_f, l = %g %g E1, E2, E3 = %g %g %g\n", sigma_f, l, E1, E2, E3);
Info->logp *= -1.0;;



    /*
     *  Compute the alpha alpha^T - K_inv matrix. Note that
     *  K_inv y was computed as r above.
     *  dger computes a x y^T + A, so for convenience, lets
     *  compute -1.0 r r^T + K_inv (which is the negative of what we want)
     */
//    gsl_matrix *A = gsl_matrix_calloc( n, n );
    gsl_matrix_memcpy( Info->A, Info->K_inv ); // initialize A with K_inv
    gsl_blas_dger( -1.0, Info->r, Info->r, Info->A ); // this is now -1*( alpha alpha^T - K_inv )
    


    /*
     * Now compute the derivs.
     * dgemm computes aAB + bC 
     * set a=-1 to revrse sign to fix sign change above
     */
//    gsl_matrix *C = gsl_matrix_calloc( n, n );

    gsl_blas_dgemm( CblasNoTrans, CblasNoTrans, -1.0, Info->A, Info->K_dl, 0.0, Info->C );
    for ( Trace = 0.0, i=0; i<n; i++ ) Trace += gsl_matrix_get( Info->C, i, i );
    Info->dlogp_dlog_l = l * 0.5*Trace;
    //Info->dlogp_dlog_l = 0.5*Trace;
Info->dlogp_dlog_l *= -1.0;

    
    gsl_blas_dgemm( CblasNoTrans, CblasNoTrans, -1.0, Info->A, Info->K_dsig, 0.0, Info->C );
    for ( Trace = 0.0, i=0; i<n; i++ ) Trace += gsl_matrix_get( Info->C, i, i );
    Info->dlogp_dlog_sig = sigma_f * 0.5*Trace;
    //Info->dlogp_dlog_sig = 0.5*Trace;
Info->dlogp_dlog_sig *= -1.0;
    
    

    
    /*
     * We are done. The 3 quantities we are after are; logp, dlogp_dl, dlogo_dsig.
     */
    //printf("l, sigma_f = %g %g logp, dlogp_dl, dlogp_dsig = %g %g %g\n", l, sigma_f, Info->logp, Info->dlogp_dl, Info->dlogp_dsig);

    return( Info->logp );


}


/*
 * The x's here are the hyper-parameters not the x's of the Gpr function
 */
void HyperParam_FdF( gsl_vector *x, void *Data, double *F, gsl_vector *dF ) {

    double log_sigma_f, log_l;

    /*
     * Typecast the Data pointer to our GprInfo Structure type.
     */
    GprInfo *g = (GprInfo *)Data;

    /*
     * Extract the hyper parameters for this call
     */
    log_sigma_f = gsl_vector_get( x, 0 );
    log_l       = gsl_vector_get( x, 1 );

//    if ( (fabs(l) < 1e-6) || (fabs(sigma_f) < 1e-6) ) {
//        *F = GSL_NAN;
//        gsl_vector_set( dF, 0, GSL_NAN );
//        gsl_vector_set( dF, 1, GSL_NAN );
//    }

    ComputeEvidenceAndDerivs( log_sigma_f, log_l, g );

    *F = g->logp;
    gsl_vector_set( dF, 0, g->dlogp_dlog_sig );
    gsl_vector_set( dF, 1, g->dlogp_dlog_l );
    //printf("HyperParam_FdF: log_sigma_f, log_l = %g %g, F, dF = %g %g %g\n", log_sigma_f, log_l, g->logp, g->dlogp_dlog_sig, g->dlogp_dlog_l );


}
double HyperParam_F( gsl_vector *x, void *Data ) {

    double log_sigma_f, log_l;

    /*
     * Typecast the Data pointer to our GprInfo Structure type.
     */
    GprInfo *g = (GprInfo *)Data;

    /*
     * Extract the hyper parameters for this call
     */
    log_sigma_f = gsl_vector_get( x, 0 );
    log_l       = gsl_vector_get( x, 1 );

//    if ( (fabs(l) < 1e-6) || (fabs(sigma_f) < 1e-6) ) {
//        return( GSL_NAN );
//    }

    ComputeEvidenceAndDerivs( log_sigma_f, log_l, g );
    //printf("HyperParam_F: log_sigma_f, log_l = %g %g, F = %g\n", log_sigma_f, log_l, g->logp );


    return( g->logp );

}
void HyperParam_dF( gsl_vector *x, void *Data, gsl_vector *dF ) {

    double log_sigma_f, log_l;

    /*
     * Typecast the Data pointer to our GprInfo Structure type.
     */
    GprInfo *g = (GprInfo *)Data;

    /*
     * Extract the hyper parameters for this call
     */
    log_sigma_f = gsl_vector_get( x, 0 );
    log_l       = gsl_vector_get( x, 1 );

//    if ( (fabs(l) < 1e-6) || (fabs(sigma_f) < 1e-6) ) {
//        gsl_vector_set( dF, 0, GSL_NAN );
//        gsl_vector_set( dF, 1, GSL_NAN );
//    }

    ComputeEvidenceAndDerivs( log_sigma_f, log_l, g );

    gsl_vector_set( dF, 0, g->dlogp_dlog_sig );
    gsl_vector_set( dF, 1, g->dlogp_dlog_l );
    //printf("HyperParam_dF: log_sigma_f, log_l = %g %g, dF = %g %g\n", log_sigma_f, log_l, g->dlogp_dlog_sig, g->dlogp_dlog_l );


}


void compute_cov_matrices( gsl_vector *x, gsl_vector *x_star,   double sigma_f, double l, gsl_matrix *K, gsl_matrix *K_star2, gsl_matrix *K_star ) {

    int i, j, n, n_star;
    double xx;

    n      = x->size;
    n_star = x_star->size;

    for ( i=0; i<n; i++ ) {
        xx = gsl_vector_get( x, i );
        for ( j=0; j<n; j++ ) {
            gsl_matrix_set( K, i, j, symm_kernel_function( xx, gsl_vector_get( x, j), sigma_f, l ) );
        }
    }

    for ( i=0; i<n_star; i++ ) {
        xx = gsl_vector_get( x_star, i );
        for ( j=0; j<n_star; j++ ) {
            gsl_matrix_set( K_star2, i, j, symm_kernel_function( xx, gsl_vector_get( x_star, j), sigma_f, l ) );
        }
    }

    for ( i=0; i<n_star; i++ ) {
        xx = gsl_vector_get( x_star, i );
        for ( j=0; j<n; j++ ) {
            gsl_matrix_set( K_star, i, j, symm_kernel_function( xx, gsl_vector_get( x, j ), sigma_f, l ) );
        }
    }

}





double func( double x ) {
    double s, f;
//    return( sin(4.0*M_PI*x) + sin(7.0*M_PI*x) );

    s = sin( M_PI*x );
    f = pow( s, 6.7 )+0.1;

    return( f );
}



/*
 * For large matrices, the normal Cholesky decomposition seems to fail.
 *
 * This routine solves this as follows: 
 * 
 *  1. Compute eignevalue decomposition A = U Lambda U^T
 * 
 *  2. Then set all negative eignevalues to zero (or small +ve number?)
 * 
 *  3. Then consider A^prime = U Lambda^prime U^T
 * 
 *  4. Note that this can be rewritten as A^prime = (U sqrt(Lambda^prime))
 *     (U sqrt(Lambda^prime))^T
 * 
 *  5. Also note that if you do a QR decomp on (U sqrt(Lambda^prime))^T ( to
 *     get = Q R ), then A^prime = (Q R)^T (Q R) = R^T R because Q is orthogonal.
 *     So Aprime = R^T R. 
 *
 *  6. Use the L = R^T matrix to get a Cholesky-like decomp.  One additional
 *     problem (potentially) is that in normal Cholesky Decomp, the diagonal
 *     elements are typical forced to be positive. If we find a negative 
 *     diagonal in our L, we just multiply the whole column by -1. This should
 *     leave L L^T unchanged.
 */
int CholeskyLikeDecomp( gsl_matrix *A_orig, gsl_matrix *L, int Method ) {

    int     n, i, j;
    double  e, val, val2;

    n = A_orig->size1;
    if ( n != A_orig->size2 ){
        printf("A matrix is not square.\n");
        return(0);
    }

    if ( n < 2 ){
        printf("A matrix size is <= 1.\n");
        return(0);
    }

    /*
     *  Perform eignevalue decomposition (Compute Lambda and U )
     */
    gsl_vector *eval              = gsl_vector_calloc( n );
    gsl_vector *eval_prime        = gsl_vector_calloc( n );
    gsl_matrix *A                 = gsl_matrix_calloc( n, n ); gsl_matrix_memcpy( A, A_orig );
    gsl_matrix *U                 = gsl_matrix_calloc( n, n );
    gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc( n );
    gsl_eigen_symmv( A, eval, U, w ); gsl_eigen_symmv_free( w );


    // set all negative eignevalues to zero (or small positive number?)
// If we get a negative eignevalue, then first try to switch its sign and also
// the sign of the eigenvector.

    for (i=0; i<n; i++){
        e = gsl_vector_get( eval, i );
        if ( Method == ZeroNegativeEigenvalues ) {
            if ( e <= 0.0 ){
                gsl_vector_set( eval_prime, i, 0.0 );
                //gsl_vector_set( eval_prime, i, 1e-15 );
            } else {
                gsl_vector_set( eval_prime, i, e );
            }
        } else if ( Method == ReverseNegativeEigenvalues ) {
            if ( e < 0.0 ){
                gsl_vector_set( eval_prime, i, -e );
                for (j=0; j<n; j++){
                    val = gsl_matrix_get( U, i, j );
                    gsl_matrix_set( U, i, j, -val );
                }
            } else {
                gsl_vector_set( eval_prime, i, e );
            }
        }
    }

    // Compute C = U sqrt(Lambda^prime)
    gsl_matrix *C    = gsl_matrix_calloc( n, n );
    //gsl_matrix *Z    = gsl_matrix_calloc( n, n );
    for (i=0; i<n; i++){ // for each row
        for (j=0; j<n; j++){ // and column
            e = gsl_vector_get( eval_prime, j );
            gsl_matrix_set( C, i, j, gsl_matrix_get( U, i, j ) * sqrt(e) );
        }
    }
    gsl_matrix_transpose( C );

    // We can check to see that C^T C recovers A
    //gsl_blas_dgemm( CblasTrans, CblasNoTrans, 1.0, C, C, 0.0, Z );
    


    /*
     * Now do a QR decomp on C
     */
    gsl_vector *tau = gsl_vector_calloc( n );
    gsl_matrix *Q   = gsl_matrix_calloc( n, n );
    gsl_matrix *R   = gsl_matrix_calloc( n, n );
    gsl_linalg_QR_decomp( C, tau );
    gsl_linalg_QR_unpack( C, tau, Q, R);

    // We can check to see that R^T R recovers A  and that Q Q^T is the identity matrix.
    //gsl_blas_dgemm( CblasTrans, CblasNoTrans, 1.0, R, R, 0.0, Z );
    //gsl_blas_dgemm( CblasTrans, CblasNoTrans, 1.0, Q, Q, 0.0, Z );


    /*
     *  Compute lower triangular matrix from R (L = R^T)
     */
    gsl_matrix_transpose( R );
    gsl_matrix_memcpy( L, R );




    /*
     * In normal Cholesky decompition, the diagonal entries should be positive
     * (they dont all have to be mathematically though.) It appears that the
     * Cholesky-like decomp we do here doesnt care about them being positive. To
     * force them to be positive, I think you can do the following: If you find a
     * negative diagonal, then changethe sign of the whole column.
     */

    for( j=0; j<n; j++ ){

        val = gsl_matrix_get( L, j, j );

        if ( val < 0.0 ) {

            // set diagonal element positive.
            gsl_matrix_set( L, j, j, fabs(val) );

            // reverse sign of other elements in the column 
            for( i=j+1; i<n; i++ ){ 
                val2 = gsl_matrix_get( L, i, j );
                gsl_matrix_set( L, i, j, -1.0*val2 );
            }

        }

    }



    return(1);


}



int main(){


    gsl_vector  *f;
    double      dx, val;
    int         i, j, ii, d, n;
    double      eps, sig;
    FILE        *fp_out, *fp_in;
    int         Method;
    gsl_rng     *rng = gsl_rng_alloc( gsl_rng_taus );

    GprInfo *Info = calloc( 1, sizeof(*Info) );

    Method = ZeroNegativeEigenvalues;

double xin[3000], yin[3000], dyin[3000];
double ymax;
n = 0; ymax = -9e99;
fp_in = fopen( "Test_PAD.dat", "r");
while ( fscanf( fp_in, "%lf %lf %lf", &xin[n], &yin[n], &dyin[n] ) != EOF ){
    if (yin[n] > ymax) ymax = yin[n];
    ++n;
}
fclose(fp_in);
// make pitch angle symmetric across 90.
//for (i=0; i<n; i++ ) {
//    xin[2*n-1-i] = 180.0 - xin[i];
//    yin[2*n-1-i] = yin[i];
//    dyin[2*n-1-i] = dyin[i];
//}
//n *=2;
//int N = n;


//if (0==1){

// add mirror image (for 180-360).
//for (i=0; i<N; i++ ) {
//    xin[2*n-1-i] = 360.0 - xin[i];
//    yin[2*n-1-i] = yin[i];
//    dyin[2*n-1-i] = dyin[i];
//}
//n += N;

// add mirror image (for 360-540).
//for (i=0; i<N; i++ ) {
//    xin[n+i] = 360.0 + xin[i];
//    yin[n+i] = yin[i];
//    dyin[n+i] = dyin[i];
//}
//n += N;
//}


InitGprInfo( n, Info );
Info->sigma_n = 0.1;
Info->sigma_n_2 = Info->sigma_n*Info->sigma_n;
for (i=0; i<n; i++ ) {
    gsl_vector_set( Info->x, i, xin[i]-90.0 ); 
    gsl_vector_set( Info->y, i, yin[i]/ymax ); 
    sig = dyin[i]/ymax;
    gsl_vector_set( Info->sigma_n_vec, i, sig ); 
    gsl_vector_set( Info->sigma_n_2_vec, i, sig*sig ); 
}




    if (0==1){
    d = 1;
    n = 100;



    InitGprInfo( n, Info );

    // Set Error std-dev
    Info->sigma_n = 0.4;
    Info->sigma_n_2 = Info->sigma_n*Info->sigma_n;

    //dx = 1.0/(double)(n-1);
    dx = 0.8/(double)(n-1);
    f  = gsl_vector_calloc( n );
    for ( i=0; i<n; i++ ) {
        val = 0.1+ i*dx;
        gsl_vector_set( Info->x, i, val );
        gsl_vector_set( f, i, func(val) );
    }

    fp_out = fopen( "Fofx.dat", "w" );
    for (i=0; i<n; i++ ){
        fprintf( fp_out, "%g %g\n", gsl_vector_get( Info->x, i ), gsl_vector_get( f, i ) );
    }
    fclose( fp_out );

    /*
     * Compute Gaussian errors 
     */
    gsl_vector *epsilon = gsl_vector_calloc( n );
    for (i=0; i<n; i++ ) {
        eps = gsl_ran_gaussian( rng, Info->sigma_n );
        gsl_vector_set( epsilon, i, eps );
        gsl_vector_set( Info->y, i, gsl_vector_get( f, i ) + eps ); 
    }
    gsl_vector_free( f );
    gsl_vector_free( epsilon );
}

    fp_out = fopen( "SampleData.dat", "w" );
    for (i=0; i<n; i++ ){
        fprintf( fp_out, "%g %g\n", gsl_vector_get( Info->x, i ), gsl_vector_get( Info->y, i ) );
    }
    fclose( fp_out );







    /*
     * Determine optimal hyper-parameteres: sigma_f, l.
     */
    gsl_multimin_function_fdf HyperParam_Func;
    HyperParam_Func.n      = 2;
    HyperParam_Func.f      = &HyperParam_F;
    HyperParam_Func.df     = &HyperParam_dF;
    HyperParam_Func.fdf    = &HyperParam_FdF;
    HyperParam_Func.params = (void *)Info;

    // starting point
    gsl_vector *hyper = gsl_vector_calloc( 2 );
    gsl_vector_set( hyper, 0, log(2.0) ); // starting point for log_sigma_f
    gsl_vector_set( hyper, 1, log(0.4) ); // starting point for log_l
gsl_vector_set( hyper, 1, log(20.) ); // starting point for log_l
    //gsl_vector_set( hyper, 0, (2.0) ); // starting point for log_sigma_f
    //gsl_vector_set( hyper, 1, (0.4) ); // starting point for log_l

    //const gsl_multimin_fdfminimizer_type *T = gsl_multimin_fdfminimizer_conjugate_fr;
    //const gsl_multimin_fdfminimizer_type *T = gsl_multimin_fdfminimizer_conjugate_pr;
    //const gsl_multimin_fdfminimizer_type *T = gsl_multimin_fdfminimizer_vector_bfgs;
    const gsl_multimin_fdfminimizer_type *T = gsl_multimin_fdfminimizer_vector_bfgs2;
    gsl_multimin_fdfminimizer            *s = gsl_multimin_fdfminimizer_alloc( T, 2 );
    gsl_multimin_fdfminimizer_set( s, &HyperParam_Func, hyper, 0.01, 1e-2 );



    size_t iter = 0;
    int status;
    do {
        iter++;
        status = gsl_multimin_fdfminimizer_iterate( s );
        if (status) break;
        status = gsl_multimin_test_gradient( s->gradient, 1e-2 );
        if (status == GSL_SUCCESS) printf( "Minimum found at:\n" );
        printf( "%5d %.5f %.5f %10.5f\n", iter, gsl_vector_get( s->x, 0 ), gsl_vector_get( s->x, 1 ), s->f );
    } while ( (status == GSL_CONTINUE) && (iter < 100) );

    gsl_vector *best_hyper;
    best_hyper = gsl_multimin_fdfminimizer_x( s ); // best_hyer doesnt need to be free'd
//    Info->log_sigma_f = fabs( gsl_vector_get( best_hyper, 0 ) );
//    Info->log_l       = fabs( gsl_vector_get( best_hyper, 1 ) );
Info->log_sigma_f = gsl_vector_get( best_hyper, 0 );
Info->log_l       = gsl_vector_get( best_hyper, 1 );
    Info->sigma_f = exp( Info->log_sigma_f );
    Info->l       = exp( Info->log_l );
// Info->sigma_f = ( Info->log_sigma_f );
// Info->l       = ( Info->log_l );
//Info->sigma_f = .20;
//Info->l       = 28.0;
//Info->l = 5.0;
printf("sigma_f, l = %g %g\n", Info->sigma_f, Info->l);
    gsl_multimin_fdfminimizer_free( s );
    gsl_vector_free( hyper );







    int        n_star;
    double     **My_K, **My_K_star2, **My_K_star;

    n_star  = 100;
//    Info->l = 0.100259;
//    Info->sigma_f = 2.22554;

    gsl_vector *x_star = gsl_vector_calloc( n_star );
    //LGM_ARRAY_2D(       My_K, n,           n, double );
    //LGM_ARRAY_2D( My_K_star2, n_star, n_star, double );
    //LGM_ARRAY_2D(  My_K_star, n_star,      n, double );

//dx = 1.0/(double)(n_star-1); for ( i=0; i<n_star; i++ ) { gsl_vector_set( x_star, i, i*dx ); }
dx = 180.0/(double)(n_star-1); for ( i=0; i<n_star; i++ ) { gsl_vector_set( x_star, i, i*dx-90.0 ); }

    gsl_matrix *K_star2 = gsl_matrix_calloc( n_star, n_star );
    gsl_matrix *K = gsl_matrix_calloc( n, n );
    gsl_matrix *K_star = gsl_matrix_calloc( n_star, n );
    compute_cov_matrices( Info->x, x_star, Info->sigma_f, Info->l, K, K_star2, K_star );

    //DumpGif3( "K_star2", n_star, n_star, My_K_star2 );
    //DumpGif3( "K_star", n_star, n, My_K_star );
    //DumpGif3( "K", n, n, My_K );


    /*
     * Do a Cholesky-like decomposition
     */
    gsl_matrix *L = gsl_matrix_calloc( n_star, n_star );
    CholeskyLikeDecomp( K_star2, L, Method );


    /*
     *  Assume means are zero
     */
    gsl_vector *mu = gsl_vector_calloc( n_star );
    for (i=0; i<n_star; i++) gsl_vector_set( mu, i, 0.0 );

    gsl_vector *z_star = gsl_vector_calloc( n_star );
    //rng = gsl_rng_alloc( gsl_rng_taus );


    /*
     * Visualize the Cholesky decomp of the K_star2 matrix
     */
    //double **MyL;
    //LGM_ARRAY_2D( MyL, n_star, n_star, double );
    //for (i=0; i<n_star; i++) {
    //    for (j=0; j<n_star; j++) {
    //        MyL[i][n_star-1-j] = gsl_matrix_get( L, i, j );
    //    }
    //}
    //DumpGif3( "L", n_star, n_star, MyL );


    /*
     * Draw some sample priors
     */
    fp_out = fopen( "Zofx_star.dat", "w" );
    for ( i=0; i<5; i++ ) {
        /*
         * Compute Gaussian errors 
         */
        gsl_ran_multivariate_gaussian( rng, mu, L, z_star );

        for (ii=0; ii<n_star; ii++ ){
            fprintf( fp_out, "%g %g\n", gsl_vector_get( x_star, ii ), gsl_vector_get( z_star, ii ) );
        }
    }
    fclose( fp_out );

//exit(0);

    /*
     * Now, lets compute the posterior.
     *
     *        f_star|X,y,X_star ~ N( fbar_star, cov(f_star) )
     *
     *  where,
     *
     *      fbar_star = K(x_star, x) ( K(x,x) + sigma_n^2 I )^-1 y
     *
     *      cov(f_star) = K(x_star, x_star) - K(x_star, x) ( K(x,x) + sigma_n^2 I )^-1 K(x,x_star)
     *
     *  So the first thing we need to do is find the inverse;
     *
     *      A = ( K(x,x) + sigma_n^2 I )^-1
     *
     *
     */
    gsl_matrix *A  = gsl_matrix_calloc( n, n );
    gsl_matrix *L1 = gsl_matrix_calloc( n, n );

    gsl_matrix_memcpy( A, K ); // copy over the whole K matrix
    for (i=0; i<n; i++) {
        // modify the diagonal elements
        gsl_matrix_set( A, i, i, gsl_matrix_get( K, i, i ) + gsl_vector_get( Info->sigma_n_2_vec, i ) );
    }
    CholeskyLikeDecomp( A, L1, Method );
    gsl_linalg_cholesky_invert( L1 );

    // Get f_bar_star = K(x_star, x) ( K(x,x) + sigma_n^2 I )^-1 y
    gsl_matrix *R          = gsl_matrix_calloc( n_star, n ); // zero matrix
    gsl_vector *f_bar_star = gsl_vector_calloc( n_star ); 
    gsl_blas_dgemm( CblasNoTrans, CblasNoTrans, 1.0, K_star, L1, 0.0, R );
    gsl_blas_dgemv( CblasNoTrans, 1.0, R, Info->y, 0.0, f_bar_star );


    // Get cov_f_star = K_star2 -   K_star ( K + sigma_n^2 I )^-1 K_star^T
    gsl_matrix *R2         = gsl_matrix_calloc( n, n_star ); // zero matrix
    gsl_matrix *cov_f_star = gsl_matrix_calloc( n_star, n_star ); 
    gsl_matrix_memcpy( cov_f_star, K_star2 );
    gsl_blas_dgemm( CblasNoTrans, CblasTrans, 1.0, L1, K_star, 0.0, R2 );
    gsl_blas_dgemm( CblasNoTrans, CblasNoTrans, -1.0, K_star, R2, 1.0, cov_f_star );

    // Decompose the cov_f_star 
    gsl_matrix *cov_f_star_L = gsl_matrix_calloc( n_star, n_star );
    CholeskyLikeDecomp( cov_f_star, cov_f_star_L, Method );



    /*
     * Allocate space for the final avg and credibility intervals.
     * And for workspace variables to compute avg and std-dev.
     */
    Info->y_hat     = gsl_vector_calloc( n_star );
    Info->y_cred_lo = gsl_vector_calloc( n_star );
    Info->y_cred_hi = gsl_vector_calloc( n_star );

    int n_samples = 1000;
    double **y_dist, yavg, ystd;
    LGM_ARRAY_2D( y_dist, n_star, n_samples, double );

    /*
     * Draw some sample priors
     */
    fp_out = fopen( "Zofx_star.dat", "w" );
    for ( i=0; i<n_samples; i++ ) {
        /*
         * Compute Gaussian errors 
         */
        gsl_ran_multivariate_gaussian( rng, f_bar_star, cov_f_star_L, z_star );

        for (ii=0; ii<n_star; ii++ ){
            y_dist[ii][i] = gsl_vector_get( z_star, ii );
            fprintf( fp_out, "%g %g\n", gsl_vector_get( x_star, ii ), gsl_vector_get( z_star, ii ) );
        }
    }
    fclose( fp_out );

    for (ii=0; ii<n_star; ii++ ){
        yavg = gsl_stats_mean( y_dist[ii], 1, n_samples );
        ystd = gsl_stats_sd(   y_dist[ii], 1, n_samples );
        gsl_vector_set( Info->y_hat, ii, yavg );
        gsl_vector_set( Info->y_cred_lo, ii, yavg - 2.0*ystd );
        gsl_vector_set( Info->y_cred_hi, ii, yavg + 2.0*ystd );
    }

    LGM_ARRAY_2D_FREE( y_dist );
    

    fp_out = fopen( "y_hat.dat", "w" );
    for (ii=0; ii<n_star; ii++ ) fprintf( fp_out, "%g %g\n", gsl_vector_get( x_star, ii ), gsl_vector_get( Info->y_hat, ii ) );
    fclose( fp_out );

    fp_out = fopen( "y_cred_lo.dat", "w" );
    for (ii=0; ii<n_star; ii++ ) fprintf( fp_out, "%g %g\n", gsl_vector_get( x_star, ii ), gsl_vector_get( Info->y_cred_lo, ii ) );
    fclose( fp_out );
    
    fp_out = fopen( "y_cred_hi.dat", "w" );
    for (ii=0; ii<n_star; ii++ ) fprintf( fp_out, "%g %g\n", gsl_vector_get( x_star, ii ), gsl_vector_get( Info->y_cred_hi, ii ) );
    fclose( fp_out );



    


    /*
     * Cleanup
     */
    gsl_vector_free( z_star );
    gsl_vector_free( mu );
    gsl_vector_free( f_bar_star );

    gsl_matrix_free( K );
    gsl_matrix_free( K_star );
    gsl_matrix_free( K_star2 );
    gsl_matrix_free( L );
    gsl_matrix_free( A );
    gsl_matrix_free( L1 );
    gsl_matrix_free( R );
    gsl_matrix_free( R2 );
    gsl_matrix_free( cov_f_star );
    gsl_matrix_free( cov_f_star_L );

    gsl_rng_free( rng );

    //LGM_ARRAY_2D_FREE( My_K );
    //LGM_ARRAY_2D_FREE( My_K_star2 );
    //LGM_ARRAY_2D_FREE( My_K_star );

    return(0);

}
