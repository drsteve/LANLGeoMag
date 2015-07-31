/*! \file Lgm_DFI_RBF.c
 *  \brief Routines to perform Scalar of data using Radial Basis Functions.
 *
 *
 *
 */

#include "Lgm/Lgm_RBF.h"

//#define LGM_RBF_SOLVER  LGM_CHOLESKY_DECOMP
#define LGM_RBF_SOLVER  LGM_PLU_DECOMP
//#define LGM_RBF_SOLVER  LGM_SVD


/** Given \f$\vec{v} = (x, y, z)\f$ and \f$\vec{v}_0 = (x_0, y_0, z_0)\f$ this
 *  routine computes the the RBF, \f$\psi(x-x_0, y-y_0, z-z_0)\f$.
 *
 *  \param[in]          v  -   input target vector
 *  \param[in]         v0  -   input reference vector
 *  \param[in]        eps  -   smoothing factor in scalar RBF
 *  \param[in,out]    Psi  -   double RBF.
 *  \param[in]        rbf  -   pointer to structure containing info for RBF interpolation. 
 *
 *  \return  void
 *
 *  \author  M. G. Henderson
 *  \date    July 15, 2015
 *
 *
 */
void    Lgm_Vec_RBF_Psi( Lgm_Vector *v, Lgm_Vector *v0, double *Psi, Lgm_Vec_RBF_Info *rbf  ) {

    double  x, y, z, x2, y2, z2, r2, f, g;
    double  e;

    e = rbf->eps;
    x = v->x - v0->x;
    y = v->y - v0->y;
    z = v->z - v0->z;
    x2 = x*x; y2 = y*y; z2 = z*z;
    r2 = x2 + y2 + z2;

    if ( rbf->RadialBasisFunction == LGM_RBF_GAUSSIAN ) {

        *Psi  = exp( -e*r2 );

    } else if ( rbf->RadialBasisFunction == LGM_RBF_MULTIQUADRIC ) {

        f   = 1.0 + e*r2;
        *Psi = sqrt( f );


    } else if ( rbf->RadialBasisFunction == LGM_RBF_INV_MULTIQUADRIC ) {

        f   = 1.0 + e*r2;
        *Psi = sqrt( f );

    } else {
        printf("Lgm_Vec_RBF_Pdi: Unknown value for rbf->RadialBasisFunction. (Got %d)\n", rbf->RadialBasisFunction );
    }

    return;

}


void    Lgm_Vec_RBF_Derivs( Lgm_Vector *v, Lgm_Vector *v0, double *dPdx, double *dPdy, double *dPdz, Lgm_Vec_RBF_Info *rbf  ) {

    double  psi, x, y, z, x2, y2, z2, r2, f, g;
    double  e;

    e = rbf->eps;
    x = v->x - v0->x;
    y = v->y - v0->y;
    z = v->z - v0->z;
    x2 = x*x; y2 = y*y; z2 = z*z;
    r2 = x2 + y2 + z2;

    if ( rbf->RadialBasisFunction == LGM_RBF_GAUSSIAN ) {

        psi  = exp( -e*r2 );
        g = -2.0*e*psi;
        *dPdx = g*x;
        *dPdy = g*y;
        *dPdz = g*z;

    } else if ( rbf->RadialBasisFunction == LGM_RBF_MULTIQUADRIC ) {

        f   = 1.0 + e*r2;
        psi = sqrt( f );
        g = e/psi;
        *dPdx = g*x;
        *dPdy = g*y;
        *dPdz = g*z;


    } else if ( rbf->RadialBasisFunction == LGM_RBF_INV_MULTIQUADRIC ) {

        f   = 1.0 + e*r2;
        psi = sqrt( f );
        g = -e/(f*psi);
        *dPdx = g*x;
        *dPdy = g*y;
        *dPdz = g*z;

    } else {
        printf("Lgm_Vec_RBF_Derivs: Unknown value for rbf->RadialBasisFunction. (Got %d)\n", rbf->RadialBasisFunction );
    }

    return;

}



/** From a vector-field dataset, compute the vector-valued weighting factors,
 *  \f$\vec{c}_j\f$. Info is returned in the rbf structure.
 *
 *
 *  \param[in]                       v   -   pointer to an array of position vectors.
 *  \param[in]                       B   -   pointer to array of corresponding field vectors.
 *  \param[in]                       n   -   number of (v, B) pairs defined.
 *  \param[in]                     eps   -   smoothing factor in scalar RBF.
 *  \param[in]      RadialBasisFunction  -   RBF to use. Can be LGM_RBF_GAUSSIAN, LGM_RBF_MULTIQUADRIC
 *
 *  \return  pointer to structure containing info for RBF interpolation. User
 *           is responsible for freeing with Lgm_Vec_RBF_Free().
 *
 *  \author  M. G. Henderson
 *  date    July 15, 2015
 *
 *
 */
Lgm_Vec_RBF_Info *Lgm_Vec_RBF_Init( unsigned long int *I_data, Lgm_Vector *v, Lgm_Vector *B, int n, double eps, int RadialBasisFunction ) {

    int              i, j, ii, jj, p, q, s;
    double           *d, **a, Psi, val;
    gsl_matrix       *A, *V;
    gsl_vector       *Dx, *Dy, *Dz, *cx, *cy, *cz, *S, *Work;
    Lgm_Vec_RBF_Info *rbf;


    /*
     * Do all three components separately, but at the sane time.
     */
    A = gsl_matrix_calloc( n, n );
    cx = gsl_vector_alloc( n );
    cy = gsl_vector_alloc( n );
    cz = gsl_vector_alloc( n );
    Dx = gsl_vector_calloc( n );
    Dy = gsl_vector_calloc( n );
    Dz = gsl_vector_calloc( n );


    /*
     * Save info needed to do an evaluation.
     */
    rbf = ( Lgm_Vec_RBF_Info *)calloc( 1, sizeof(*rbf) );
    rbf->RadialBasisFunction = RadialBasisFunction;
    rbf->eps = eps;
    rbf->n   = n;
    LGM_ARRAY_1D( rbf->LookUpKey, n, unsigned long int);
    LGM_ARRAY_1D( rbf->v, n, Lgm_Vector);
    LGM_ARRAY_1D( rbf->cx, n, double);
    LGM_ARRAY_1D( rbf->cy, n, double);
    LGM_ARRAY_1D( rbf->cz, n, double);
    for ( i=0; i<n; i++ ) {
        rbf->LookUpKey[i] = I_data[i];
        rbf->v[i] = v[i];
    }
    // This subtraction doesntm seem to work out very well...?
//    rbf->Bx0 = B[0].x;
//    rbf->By0 = B[0].y;
//    rbf->Bz0 = B[0].z;

double Bbkg;
for ( Bbkg = 0.0, i=0; i<n; i++ ) Bbkg += B[i].x; rbf->Bx0 = Bbkg/(double)n;
for ( Bbkg = 0.0, i=0; i<n; i++ ) Bbkg += B[i].y; rbf->By0 = Bbkg/(double)n;
for ( Bbkg = 0.0, i=0; i<n; i++ ) Bbkg += B[i].z; rbf->Bz0 = Bbkg/(double)n;
//    rbf->Bx0 = 0.0;
//    rbf->By0 = 0.0;
//    rbf->Bz0 = 0.0;
    
    /*
     * Fill D arrays. (Subtract off the field at the nearest point v[0] -- See
     * McNally [2011].) We add this field back on later.
     */
    for (i=0; i<n; i++){
        gsl_vector_set( Dx, i, B[i].x - rbf->Bx0 );
        gsl_vector_set( Dy, i, B[i].y - rbf->By0 );
        gsl_vector_set( Dz, i, B[i].z - rbf->Bz0 );
    }


    /*
     *                                                                         [  row0  ]
     * Fill A matrix (same for all 3 conponents). In C, order is A[row][col] = [  row1  ]
     *                                                                         [  row2  ]
     */
    for ( i=0; i<n; i++ ) { // locate start row for subarray
        for ( j=i; j<n; j++ ) { // locate start column for subarray
            // Get PSi( v_i - v_j )
            Lgm_Vec_RBF_Psi( &v[i], &v[j], &Psi, rbf );
            gsl_matrix_set( A, i, j, Psi );
            if ( i!=j) gsl_matrix_set( A, j, i, Psi );
        }
    }

    /*
    for (i=0; i<n; i++ ) {
        printf("v%02d = %8g %8g %8g   B%02d = %8g %8g %8g\n", i, v[i].x, v[i].y, v[i].z, i, B[i].x, B[i].y, B[i].z );
    }
    for (i=0; i<n; i++){
        for (j=0; j<n; j++){
            printf("%8g ", gsl_matrix_get(A, i, j ) );
        }
        printf("\n");
    }
    */


    /*
     * Now we need to solve the system of equation;
     *
     *      d = ac
     *
     *  for c. (One for each comp.)
     *
     *  First create gsl_vector and gsl_matrix views of the d and A arrays.
     *  Then compute Cholesky decomposition of the a array. Then solve the
     *  system to get c.
     *
     */
    if ( LGM_RBF_SOLVER == LGM_CHOLESKY_DECOMP ){
        gsl_linalg_cholesky_decomp( A );
        gsl_linalg_cholesky_solve( A, Dx, cx ); // x-comp
        gsl_linalg_cholesky_solve( A, Dy, cy ); // y-comp
        gsl_linalg_cholesky_solve( A, Dz, cz ); // z-comp
    } else if ( LGM_RBF_SOLVER == LGM_PLU_DECOMP ){
        gsl_permutation *P = gsl_permutation_alloc( n );
        gsl_linalg_LU_decomp( A, P, &s );
        gsl_linalg_LU_solve( A, P, Dx, cx ); // x-comp
        gsl_linalg_LU_solve( A, P, Dy, cy ); // y-comp
        gsl_linalg_LU_solve( A, P, Dz, cz ); // z-comp
        gsl_permutation_free( P );
    } else if ( LGM_RBF_SOLVER == LGM_SVD ){
        V    = gsl_matrix_calloc( n, n );
        S    = gsl_vector_alloc( n );
        Work = gsl_vector_alloc( n );
        gsl_linalg_SV_decomp( A, V, S, Work );
        gsl_linalg_SV_solve( A, V, S, Dx, cx ); //x-comp
        gsl_linalg_SV_solve( A, V, S, Dy, cy ); //y-comp
        gsl_linalg_SV_solve( A, V, S, Dz, cz ); //z-comp
        gsl_vector_free( Work );
        gsl_vector_free( S );
        gsl_matrix_free( V );
    }

    for (i=0; i<n; i++){
        rbf->cx[i] = gsl_vector_get( cx, i );
        rbf->cy[i] = gsl_vector_get( cy, i );
        rbf->cz[i] = gsl_vector_get( cz, i );
    }

//exit(0);

    
    gsl_vector_free( Dx );
    gsl_vector_free( Dy );
    gsl_vector_free( Dz );
    gsl_vector_free( cx );
    gsl_vector_free( cy );
    gsl_vector_free( cz );
    gsl_matrix_free( A );

    return( rbf );

}

/** Free a previously allocated Lgm_DFI_RBF_Info structure.
 *
 *
 *  \param[in]     rbf   -   pointer to structure containing info for RBF interpolation. 
 *
 *  \return  void
 *
 *  \author  M. G. Henderson
 *  \date    January 24, 2012
 *
 *
 */
void    Lgm_Vec_RBF_Free( Lgm_Vec_RBF_Info *rbf ) {
    LGM_ARRAY_1D_FREE( rbf->LookUpKey );
    LGM_ARRAY_1D_FREE( rbf->v );
    LGM_ARRAY_1D_FREE( rbf->cx );
    LGM_ARRAY_1D_FREE( rbf->cy );
    LGM_ARRAY_1D_FREE( rbf->cz );
    free( rbf );
    return;
}







/**  Compute Divergence Free interpolant at the specified position vector. The
 *   weights given in \f$\vec{c}\f$ must have been pre-computed with
 *   Lgm_DFI_RBF_Init().
 *
 *
 *  \param[in]        v   -   position vector to compute B at.
 *  \param[out]       B   -   interpolated value of B at v.
 *  \param[out]     rbf   -   pointer to initialized Lgm_DFI_RBF_Info structure.
 *
 *  \return  void
 *
 *  \author  M. G. Henderson
 *  \date    January 24, 2012
 *
 *
 */
void    Lgm_Vec_RBF_Eval( Lgm_Vector *v, Lgm_Vector *B, Lgm_Vec_RBF_Info *rbf ) {

    int         j;
    double      psi;


    B->x = B->y = B->z = 0.0;
    for ( j=0; j<rbf->n; j++ ){
        Lgm_Vec_RBF_Psi( v, &rbf->v[j], &psi, rbf );
        B->x += psi * rbf->cx[j];
        B->y += psi * rbf->cy[j];
        B->z += psi * rbf->cz[j];
    }

    // Add the subtracted "background" back in (that we subtracted prior to the interp).
    B->x += rbf->Bx0;
    B->y += rbf->By0;
    B->z += rbf->Bz0;

    return;
    
}





/**  Compute Vector interpolant at the specified position vector. The
 *   weights given in \f$\vec{cx}, \vec{cy}, \vec{cz}\f$ must have been pre-computed with
 *   Lgm_Vec_RBF_Init().
 *
 *
 *  \param[in]         v   -   position vector to compute B at.
 *  \param[out]     dBdx   -   interpolated value of dB/dx at v.
 *  \param[out]     dBdy   -   interpolated value of dB/dy at v.
 *  \param[out]     dBdz   -   interpolated value of dB/dz at v.
 *  \param[out]      rbf   -   pointer to initialized Lgm_Vec_RBF_Info structure.
 *
 *  \return  void
 *
 *  \author  M. G. Henderson
 *  \date    July 15, 2015
 *
 *
 */
void    Lgm_Vec_RBF_Derivs_Eval( Lgm_Vector *v, Lgm_Vector *dBdx, Lgm_Vector *dBdy, Lgm_Vector *dBdz, Lgm_Vec_RBF_Info *rbf ) {

    int         j;
    double      dPdx, dPdy, dPdz;


    dBdx->x = dBdx->y = dBdx->z = 0.0;
    dBdy->x = dBdy->y = dBdy->z = 0.0;
    dBdz->x = dBdz->y = dBdz->z = 0.0;
    for ( j=0; j<rbf->n; j++ ){

        Lgm_Vec_RBF_Derivs( v, &rbf->v[j], &dPdx, &dPdy, &dPdz, rbf );

        dBdx->x += dPdx * rbf->cx[j]; 
        dBdx->y += dPdx * rbf->cy[j]; 
        dBdx->z += dPdx * rbf->cz[j]; 

        dBdy->x += dPdy * rbf->cx[j]; 
        dBdy->y += dPdy * rbf->cy[j]; 
        dBdy->z += dPdy * rbf->cz[j]; 

        dBdz->x += dPdz * rbf->cx[j];
        dBdz->y += dPdz * rbf->cy[j];
        dBdz->z += dPdz * rbf->cz[j];

    }


    return;
    
}


