#include "Lgm/GPR.h"

double kernel_function( double x, double y, double sigma_f, double l, double sigma_b, double sigma_nu, double c ) {

    double  diff, Norm2;
    double  k_rbf, k_lin;

    /*
     * Compute the norm for x-y. norm = sqrt( sum( (x[i]-y[i])^2 ) ^1/2
     * since we end up squaring the norm later dont need sqrt.
     */
    diff = x - y;
    Norm2 = diff*diff;
    // Norm = sqrt(Norm2);

    k_rbf = sigma_f * exp( -Norm2 / (2.0*l*l) );
//    k_lin = sigma_b*sigma_b + sigma_nu*sigma_nu*(x-c)*(y-c);

    return( k_rbf + k_lin );

}

/*
 * Derivs for k_lin
 */
double kernel_function_dk_dsigb( double x, double y, double sigma_b, double sigma_nu, double c ) {
    return( 2.0*sigma_b );
}

double kernel_function_dk_dsignu( double x, double y, double sigma_b, double sigma_nu, double c ) {
    return( 2.0*sigma_nu*(x-c)*(y-c) );
}

double kernel_function_dk_dc( double x, double y, double sigma_b, double sigma_nu, double c ) {
    return( -sigma_nu*sigma_nu*( x + y - 2.0*c ) );
}




/*
 * Derivs for k_rbf
 */
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
double symm_kernel_function( double x, double y, double sigma_f, double l, double sigma_b, double sigma_nu, double c ) {
    return( kernel_function( x, y, sigma_f, l, sigma_b, sigma_nu, c ) + kernel_function( -x, y, sigma_f, l, sigma_b, sigma_nu, c ) );
}
double symm_kernel_function_dk_dsig( double x, double y, double sigma_f, double l ) {
    return( kernel_function_dk_dsig( x, y, sigma_f, l ) + kernel_function_dk_dsig( -x, y, sigma_f, l ) );
}
double symm_kernel_function_dk_dl( double x, double y, double sigma_f, double l ) {
    return( kernel_function_dk_dl( x, y, sigma_f, l ) + kernel_function_dk_dl( -x, y, sigma_f, l ) );
}
double symm_kernel_function_dk_dsigb( double x, double y, double sigma_b, double sigma_nu, double c ) {
    return( kernel_function_dk_dsigb( x, y, sigma_b, sigma_nu, c ) + kernel_function_dk_dsigb( -x, y, sigma_b, sigma_nu, c ) );
}
double symm_kernel_function_dk_dsignu( double x, double y, double sigma_b, double sigma_nu, double c ) {
    return( kernel_function_dk_dsignu( x, y, sigma_b, sigma_nu, c ) + kernel_function_dk_dsignu( -x, y, sigma_b, sigma_nu, c ) );
}
double symm_kernel_function_dk_dc( double x, double y, double sigma_b, double sigma_nu, double c ) {
    return( kernel_function_dk_dc( x, y, sigma_b, sigma_nu, c ) + kernel_function_dk_dc( -x, y, sigma_b, sigma_nu, c ) );
}





GprInfo *InitGprInfo( int n, int n_star, int ForceSymmetry ) {

    GprInfo *Info = calloc( 1, sizeof(*Info) );

    Info->ForceSymmetry = ( ForceSymmetry > 0 ) ? 1 : 0;


    Info->n      = n;
    Info->x      = gsl_vector_calloc( n );
    Info->y      = gsl_vector_calloc( n );
    Info->sigma_n_vec = gsl_vector_calloc( n );
    Info->sigma_n_2_vec = gsl_vector_calloc( n );

    Info->K        = gsl_matrix_calloc( n, n );
    Info->K_dl     = gsl_matrix_calloc( n, n );
    Info->K_dsig   = gsl_matrix_calloc( n, n );
    Info->K_dsigb  = gsl_matrix_calloc( n, n );
    Info->K_dsignu = gsl_matrix_calloc( n, n );
    Info->K_dc     = gsl_matrix_calloc( n, n );
    Info->K_inv    = gsl_matrix_calloc( n, n );
    Info->r        = gsl_vector_calloc( n );
    Info->A        = gsl_matrix_calloc( n, n );
    Info->C        = gsl_matrix_calloc( n, n );


    Info->n_star    = n_star;
    Info->x_star    = gsl_vector_calloc( n_star );
    Info->y_hat     = gsl_vector_calloc( n_star );
    Info->y_cred_lo = gsl_vector_calloc( n_star );
    Info->y_cred_hi = gsl_vector_calloc( n_star );

    return( Info );

}


void FreeGprInfo( GprInfo *Info ) {


    gsl_vector_free( Info->x );
    gsl_vector_free( Info->y );
    gsl_vector_free( Info->sigma_n_vec );
    gsl_vector_free( Info->sigma_n_2_vec );

    gsl_matrix_free( Info->K      );
    gsl_matrix_free( Info->K_dl   );
    gsl_matrix_free( Info->K_dsig );
    gsl_matrix_free( Info->K_dsigb );
    gsl_matrix_free( Info->K_dsignu );
    gsl_matrix_free( Info->K_dc );
    gsl_matrix_free( Info->K_inv  );
    gsl_vector_free( Info->r      );
    gsl_matrix_free( Info->A      );
    gsl_matrix_free( Info->C      );

    gsl_vector_free( Info->x_star    );
    gsl_vector_free( Info->y_hat     );
    gsl_vector_free( Info->y_cred_lo );
    gsl_vector_free( Info->y_cred_hi );

    free( Info );

    return;

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

    gsl_vector_free( eval       );
    gsl_vector_free( eval_prime );
    gsl_matrix_free( A          );
    gsl_matrix_free( U          );
    gsl_matrix_free( C          );
    gsl_vector_free( tau        );
    gsl_matrix_free( Q          );
    gsl_matrix_free( R          );


    return(1);


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
 *  (5 when we add in the 3 for linear kernel)
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
void ComputeEvidenceAndDerivs( double log_sigma_f, double log_l, double log_sigma_b, double log_sigma_nu, double c, GprInfo *Info ) {


    int     n, i, j, Method;
    double  log_K_det, Trace, E1, E2, E3, val, xi, xj;
    double  sig, sigma_f, l, sigma_b, sigma_nu;
    gsl_vector *x, *y;

    x = Info->x;
    y = Info->y;

    
    sigma_f  = exp( log_sigma_f );
    l        = exp( log_l );
    sigma_b  = exp( log_sigma_b );
    sigma_nu = exp( log_sigma_nu );

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
                if ( Info->ForceSymmetry ) {
                    gsl_matrix_set( Info->K, i, j, symm_kernel_function( xi, xj, sigma_f, l, sigma_b, sigma_nu, c ) + sig );
                } else {
                    gsl_matrix_set( Info->K, i, j, kernel_function( xi, xj, sigma_f, l, sigma_b, sigma_nu, c ) + sig );
                }
            } else {
                if ( Info->ForceSymmetry ) {
                    gsl_matrix_set( Info->K, i, j, symm_kernel_function( xi, xj, sigma_f, l, sigma_b, sigma_nu, c ) );
                } else {
                    gsl_matrix_set( Info->K, i, j, kernel_function( xi, xj, sigma_f, l, sigma_b, sigma_nu, c ) );
                }
            }
            if ( Info->ForceSymmetry ) {
                gsl_matrix_set( Info->K_dl,     i, j, symm_kernel_function_dk_dl( xi, xj, sigma_f, l ) );
//printf("symm_kernel_function_dk_dl( %g, %g, %g, %g ) = %g\n", xi, xj, sigma_f, l, symm_kernel_function_dk_dl( xi, xj, sigma_f, l ));
                gsl_matrix_set( Info->K_dsig,   i, j, symm_kernel_function_dk_dsig( xi, xj, sigma_f, l ) );
//printf("symm_kernel_function_dk_dsig( %g, %g, %g, %g ) = %g\n",  xi, xj, sigma_f, l, symm_kernel_function_dk_dsig( xi, xj, sigma_f, l ));
                gsl_matrix_set( Info->K_dsigb,  i, j, symm_kernel_function_dk_dsigb( xi, xj, sigma_b, sigma_nu, c ) );
                gsl_matrix_set( Info->K_dsignu, i, j, symm_kernel_function_dk_dsignu( xi, xj, sigma_b, sigma_nu, c ) );
                gsl_matrix_set( Info->K_dc,     i, j, symm_kernel_function_dk_dc( xi, xj, sigma_b, sigma_nu, c ) );
            } else {
                gsl_matrix_set( Info->K_dl,     i, j, kernel_function_dk_dl( xi, xj, sigma_f, l ) );
                gsl_matrix_set( Info->K_dsig,   i, j, kernel_function_dk_dsig( xi, xj, sigma_f, l ) );
                gsl_matrix_set( Info->K_dsigb,  i, j, kernel_function_dk_dsigb( xi, xj, sigma_b, sigma_nu, c ) );
                gsl_matrix_set( Info->K_dsignu, i, j, kernel_function_dk_dsignu( xi, xj, sigma_b, sigma_nu, c ) );
                gsl_matrix_set( Info->K_dc,     i, j, kernel_function_dk_dc( xi, xj, sigma_b, sigma_nu, c ) );
            }
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
        //printf("val = %g\n", val);
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
    if (VERBOSE > 0 ) printf("sigma_f, l = %g %g E1, E2, E3 = %g %g %g\n", sigma_f, l, E1, E2, E3);
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
    



    gsl_blas_dgemm( CblasNoTrans, CblasNoTrans, -1.0, Info->A, Info->K_dsigb, 0.0, Info->C );
    for ( Trace = 0.0, i=0; i<n; i++ ) Trace += gsl_matrix_get( Info->C, i, i );
    Info->dlogp_dlog_sigb = sigma_b * 0.5*Trace;
    Info->dlogp_dlog_sigb *= -1.0;
    
    gsl_blas_dgemm( CblasNoTrans, CblasNoTrans, -1.0, Info->A, Info->K_dsignu, 0.0, Info->C );
    for ( Trace = 0.0, i=0; i<n; i++ ) Trace += gsl_matrix_get( Info->C, i, i );
    Info->dlogp_dlog_signu = sigma_nu * 0.5*Trace;
    Info->dlogp_dlog_signu *= -1.0;

    gsl_blas_dgemm( CblasNoTrans, CblasNoTrans, -1.0, Info->A, Info->K_dc, 0.0, Info->C );
    for ( Trace = 0.0, i=0; i<n; i++ ) Trace += gsl_matrix_get( Info->C, i, i );
    Info->dlogp_dc = 0.5*Trace;
    Info->dlogp_dc *= -1.0;

    
    /*
     * We are done. The 3 quantities we are after are; logp, dlogp_dl, dlogo_dsig.
     */
/*
printf("l, sigma_f, sigma_b, sigma_nu, c  = %g %g %g %g %g logp, dlogp_dl, dlogp_dsig, dlogp_dlog_sigb, dlogp_dsig, dlogp_dc = %g %g %g %g %g %g\n", 
l, sigma_f, sigma_b, sigma_nu, c, 
Info->logp, 
Info->dlogp_dlog_l, Info->dlogp_dlog_sig, 
Info->dlogp_dlog_sigb, Info->dlogp_dlog_signu, Info->dlogp_dc );
*/

    return;


}


/*
 * The x's here are the hyper-parameters not the x's of the Gpr function
 */
void HyperParam_FdF( const gsl_vector *x, void *Data, double *F, gsl_vector *dF ) {

    double log_sigma_f, log_l;
    double log_sigma_b, log_sigma_nu, c;

    /*
     * Typecast the Data pointer to our GprInfo Structure type.
     */
    GprInfo *g = (GprInfo *)Data;

    /*
     * Extract the hyper parameters for this call
     */
    log_sigma_f  = gsl_vector_get( x, 0 );
log_l = log( 20.0 );
//    log_l        = gsl_vector_get( x, 1 );
//    log_sigma_b  = gsl_vector_get( x, 2 );
//    log_sigma_nu = gsl_vector_get( x, 3 );
//    c            = gsl_vector_get( x, 4 );

//    if ( (fabs(l) < 1e-6) || (fabs(sigma_f) < 1e-6) ) {
//        *F = GSL_NAN;
//        gsl_vector_set( dF, 0, GSL_NAN );
//        gsl_vector_set( dF, 1, GSL_NAN );
//    }

    ComputeEvidenceAndDerivs( log_sigma_f, log_l, log_sigma_b, log_sigma_nu, c,    g );

    *F = g->logp;
    gsl_vector_set( dF, 0, g->dlogp_dlog_sig );
//    gsl_vector_set( dF, 1, g->dlogp_dlog_l );
//    gsl_vector_set( dF, 2, g->dlogp_dlog_sigb );
//    gsl_vector_set( dF, 3, g->dlogp_dlog_signu );
//    gsl_vector_set( dF, 4, g->dlogp_dc );
    //printf("HyperParam_FdF: log_sigma_f, log_l = %g %g, F, dF = %g %g %g\n", log_sigma_f, log_l, g->logp, g->dlogp_dlog_sig, g->dlogp_dlog_l );


}
double HyperParam_F( const gsl_vector *x, void *Data ) {

    double log_sigma_f, log_l, log_sigma_b, log_sigma_nu, c;

    /*
     * Typecast the Data pointer to our GprInfo Structure type.
     */
    GprInfo *g = (GprInfo *)Data;

    /*
     * Extract the hyper parameters for this call
     */
    log_sigma_f  = gsl_vector_get( x, 0 );
//    log_l        = gsl_vector_get( x, 1 );
//    log_sigma_b  = gsl_vector_get( x, 2 );
//    log_sigma_nu = gsl_vector_get( x, 3 );
//    c            = gsl_vector_get( x, 4 );

//    if ( (fabs(l) < 1e-6) || (fabs(sigma_f) < 1e-6) ) {
//        return( GSL_NAN );
//    }

log_l = log( 20.0 );
    ComputeEvidenceAndDerivs( log_sigma_f, log_l, log_sigma_b, log_sigma_nu, c,   g );
    //printf("HyperParam_F: log_sigma_f, log_l = %g %g, F = %g\n", log_sigma_f, log_l, g->logp );


    return( g->logp );

}
void HyperParam_dF( const gsl_vector *x, void *Data, gsl_vector *dF ) {

    double log_sigma_f, log_l, log_sigma_b, log_sigma_nu, c;

    /*
     * Typecast the Data pointer to our GprInfo Structure type.
     */
    GprInfo *g = (GprInfo *)Data;

    /*
     * Extract the hyper parameters for this call
     */
    log_sigma_f  = gsl_vector_get( x, 0 );
//    log_l        = gsl_vector_get( x, 1 );
//    log_sigma_b  = gsl_vector_get( x, 2 );
//    log_sigma_nu = gsl_vector_get( x, 3 );
//    c            = gsl_vector_get( x, 4 );

//    if ( (fabs(l) < 1e-6) || (fabs(sigma_f) < 1e-6) ) {
//        gsl_vector_set( dF, 0, GSL_NAN );
//        gsl_vector_set( dF, 1, GSL_NAN );
//    }

log_l = log( 20.0 );
    ComputeEvidenceAndDerivs( log_sigma_f, log_l, log_sigma_b, log_sigma_nu, c,   g );

    gsl_vector_set( dF, 0, g->dlogp_dlog_sig );
//    gsl_vector_set( dF, 1, g->dlogp_dlog_l );
//    gsl_vector_set( dF, 2, g->dlogp_dlog_sigb );
//    gsl_vector_set( dF, 3, g->dlogp_dlog_signu );
//    gsl_vector_set( dF, 4, g->dlogp_dc );
    //printf("HyperParam_dF: log_sigma_f, log_l = %g %g, dF = %g %g\n", log_sigma_f, log_l, g->dlogp_dlog_sig, g->dlogp_dlog_l );


}

//void compute_cov_matrices( gsl_vector *x, gsl_vector *x_star,   double sigma_f, double l, gsl_matrix *K, gsl_matrix *K_star2, gsl_matrix *K_star ) {
void compute_cov_matrices( double sigma_f, double l, double sigma_b, double sigma_nu, double c, gsl_matrix *K, gsl_matrix *K_star2, gsl_matrix *K_star, GprInfo *Info ) {

    int i, j, n, n_star;
    double xx;

    n      = Info->x->size;
    n_star = Info->x_star->size;

    for ( i=0; i<n; i++ ) {
        xx = gsl_vector_get( Info->x, i );
        for ( j=0; j<n; j++ ) {
            if ( Info->ForceSymmetry ) {
                gsl_matrix_set( K, i, j, symm_kernel_function( xx, gsl_vector_get( Info->x, j), sigma_f, l, sigma_b, sigma_nu, c ) );
            } else {
                gsl_matrix_set( K, i, j,      kernel_function( xx, gsl_vector_get( Info->x, j), sigma_f, l, sigma_b, sigma_nu, c ) );
            }
        }
    }

    for ( i=0; i<n_star; i++ ) {
        xx = gsl_vector_get( Info->x_star, i );
        for ( j=0; j<n_star; j++ ) {
            if ( Info->ForceSymmetry ) {
                gsl_matrix_set( K_star2, i, j, symm_kernel_function( xx, gsl_vector_get( Info->x_star, j), sigma_f, l, sigma_b, sigma_nu, c ) );
            } else {
                gsl_matrix_set( K_star2, i, j,      kernel_function( xx, gsl_vector_get( Info->x_star, j), sigma_f, l, sigma_b, sigma_nu, c ) );
            }
        }
    }

    for ( i=0; i<n_star; i++ ) {
        xx = gsl_vector_get( Info->x_star, i );
        for ( j=0; j<n; j++ ) {
            if ( Info->ForceSymmetry ) {
                gsl_matrix_set( K_star, i, j, symm_kernel_function( xx, gsl_vector_get( Info->x, j ), sigma_f, l, sigma_b, sigma_nu, c ) );
            } else {
                gsl_matrix_set( K_star, i, j,      kernel_function( xx, gsl_vector_get( Info->x, j ), sigma_f, l, sigma_b, sigma_nu, c ) );
            }
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






void GPR( GprInfo *Info ){


    int         i, ii, n, n_star;
    FILE        *fp_out;
    int         Method;
    gsl_rng     *rng = gsl_rng_alloc( gsl_rng_taus );


    Method = ZeroNegativeEigenvalues;

    n      = Info->n;
    n_star = Info->n_star;


    /*
     * Determine optimal hyper-parameteres: sigma_f, l.
     */
    gsl_multimin_function_fdf HyperParam_Func;
    HyperParam_Func.n      = 1;
    HyperParam_Func.f      = &HyperParam_F;
    HyperParam_Func.df     = &HyperParam_dF;
    HyperParam_Func.fdf    = &HyperParam_FdF;
    HyperParam_Func.params = (void *)Info;

    // starting point
    gsl_vector *hyper = gsl_vector_calloc( 1 );
    gsl_vector_set( hyper, 0, log(0.05) ); // starting point for log_sigma_f
//    gsl_vector_set( hyper, 1, log(20.0) ); // starting point for log_l
//    gsl_vector_set( hyper, 2, log(0.1) ); // starting point for log_sigb
//    gsl_vector_set( hyper, 3, log(0.1) ); // starting point for log_signu
//    gsl_vector_set( hyper, 4, 0.0 );      // starting point for c

    //const gsl_multimin_fdfminimizer_type *T = gsl_multimin_fdfminimizer_conjugate_fr;
    //const gsl_multimin_fdfminimizer_type *T = gsl_multimin_fdfminimizer_conjugate_pr;
    //const gsl_multimin_fdfminimizer_type *T = gsl_multimin_fdfminimizer_vector_bfgs;
    const gsl_multimin_fdfminimizer_type *T = gsl_multimin_fdfminimizer_vector_bfgs2;
    gsl_multimin_fdfminimizer            *s = gsl_multimin_fdfminimizer_alloc( T, 1 );
    gsl_multimin_fdfminimizer_set( s, &HyperParam_Func, hyper, 0.01, 1e-2 );


if (0==1){
    size_t  iter = 0;
    int     status;
    do {
        iter++;
        status = gsl_multimin_fdfminimizer_iterate( s ); if (status) break;
        status = gsl_multimin_test_gradient( s->gradient, 1e-2 ); if ((VERBOSE>0)&&(status == GSL_SUCCESS)) printf( "Minimum found at:\n" );
        if (VERBOSE>0) printf( "%5d %.5f %.5f %10.5f\n", (int)iter, gsl_vector_get( s->x, 0 ), gsl_vector_get( s->x, 1 ), s->f );
    } while ( (status == GSL_CONTINUE) && (iter < 100) );
}

    gsl_vector *best_hyper;
    best_hyper = gsl_multimin_fdfminimizer_x( s ); // best_hyer doesnt need to be free'd
    Info->log_sigma_f  = gsl_vector_get( best_hyper, 0 );
//    if ( Info->log_sigma_f > 20.0 ) Info->log_sigma_f = 20.0;
    
//    Info->log_l        = gsl_vector_get( best_hyper, 1 );
 //   Info->log_sigma_b  = gsl_vector_get( best_hyper, 2 );
//    Info->log_sigma_nu = gsl_vector_get( best_hyper, 3 );
//    Info->c            = gsl_vector_get( best_hyper, 4 );

    Info->sigma_f  = exp( Info->log_sigma_f );
    Info->l        = exp( Info->log_l );
//Info->sigma_f  = 4.0;
//Info->l        = 0.050;
//Info->sigma_f  = 0.01;
Info->l        = 20.0;
if ( Info->sigma_f > 1.0 ) Info->sigma_f = 1.0;
Info->sigma_f = 0.1;
//Info->l = 20.0;
//Info->l = 45.0;
//Info->l = 60.0;
// Info->sigma_f = 1e6;
//Info->l = 10.0;
//Info->l        = 10.0;
//    Info->sigma_b  = exp( Info->log_sigma_b );
//    Info->sigma_nu = exp( Info->log_sigma_nu );
//Info->sigma_f = 0.5;
//if (Info->l < 0.25) Info->l = 0.25;
 Info->sigma_b  = 0.0;
 Info->sigma_nu  = 0.0;
 Info->c  = 0.0;
    if (VERBOSE>0) printf("sigma_f, l, sigma_b, sigma_nu, c  = %g %g %g %g %g\n", Info->sigma_f, Info->l, Info->sigma_b, Info->sigma_nu, Info->c );
//Info->sigma_f = 10.0;
//Info->l = 30.0;
    gsl_multimin_fdfminimizer_free( s );
    gsl_vector_free( hyper );


    /*
     * This probably should be a user input -- user will want to specify points here.
     */
//    int        n_star;
//    n_star  = 100;
//    gsl_vector *x_star = gsl_vector_calloc( n_star );
//    dx = 180.0/(double)(n_star-1); 
//    for ( i=0; i<n_star; i++ ) { gsl_vector_set( x_star, i, i*dx-90.0 ); }



    gsl_matrix *K_star2 = gsl_matrix_calloc( n_star, n_star );
    gsl_matrix *K = gsl_matrix_calloc( n, n );
    gsl_matrix *K_star = gsl_matrix_calloc( n_star, n );
    compute_cov_matrices( Info->sigma_f, Info->l, Info->sigma_b, Info->sigma_nu, Info->c, K, K_star2, K_star, Info );


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


    /*
     * Draw some sample priors -- this isnt really needed except for visualization...
     */
    if (0==1){ // Could return these to user upon request...
        fp_out = fopen( "Zofx_star.dat", "w" );
        for ( i=0; i<5; i++ ) {
            /*
             * Compute Gaussian errors 
             */
            gsl_ran_multivariate_gaussian( rng, mu, L, z_star );

            for (ii=0; ii<n_star; ii++ ){
                fprintf( fp_out, "%g %g\n", gsl_vector_get( Info->x_star, ii ), gsl_vector_get( z_star, ii ) );
            }
        }
        fclose( fp_out );
    }


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
     * Allocate space for workspace variables to compute avg and std-dev.
     */

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
            fprintf( fp_out, "%g %g\n", gsl_vector_get( Info->x_star, ii ), gsl_vector_get( z_star, ii ) );
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

    return;

}


