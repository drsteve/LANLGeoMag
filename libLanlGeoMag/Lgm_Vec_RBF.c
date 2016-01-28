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
 *  \param[in]          e  -   input epsilon factor
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
void    Lgm_Vec_RBF_Psi( Lgm_Vector *v, Lgm_Vector *v0, double e, double *Psi, int RadialBasisFunction  ) {

    double  r, omr, omr2, omr4;
    double  x, y, z, x2, y2, z2, r2, f, g;

    x = v->x - v0->x;
    y = v->y - v0->y;
    z = v->z - v0->z;
    x2 = x*x; y2 = y*y; z2 = z*z;
    r2 = x2 + y2 + z2;

    if ( RadialBasisFunction == LGM_RBF_GAUSSIAN ) {

        *Psi  = exp( -e*r2 );

    } else if ( RadialBasisFunction == LGM_RBF_MULTIQUADRIC ) {

        f   = 1.0 + e*r2;
        *Psi = sqrt( f );


    } else if ( RadialBasisFunction == LGM_RBF_INV_MULTIQUADRIC ) {

        f   = 1.0 + e*r2;
        *Psi = sqrt( f );

    } else if ( RadialBasisFunction == LGM_RBF_WENDLAND30 ) {
        
        r = sqrt(e*r2);
        omr = 1.0 - r; if ( omr < 0.0 ) omr = 0.0;
        *Psi = omr*omr;

    } else if ( RadialBasisFunction == LGM_RBF_WENDLAND31 ) {
        
        r = sqrt(e*r2);
        omr = 1.0 - r; if ( omr < 0.0 ) omr = 0.0;
        omr2 = omr*omr;
        *Psi = omr2*omr2*(4.0*r+1.0);

    } else if ( RadialBasisFunction == LGM_RBF_WENDLAND32 ) {
        
        r = sqrt(e*r2);
        omr = 1.0 - r; if ( omr < 0.0 ) omr = 0.0;
        omr2 = omr*omr;
        *Psi = omr2*omr2*omr2*(35.0*r2 + 18.0*r + 3.0);

    } else if ( RadialBasisFunction == LGM_RBF_WENDLAND33 ) {
        
        r = sqrt(e*r2);
        omr = 1.0 - r; if ( omr < 0.0 ) omr = 0.0;
        omr2 = omr*omr;
        omr4 = omr2*omr2;
        *Psi = omr4*omr4*(32.0*r2*r + 25.0*r2 + 8.0*r + 1.0);


    } else {
        printf("Lgm_Vec_RBF_Psi: Unknown value for RadialBasisFunction. (Got %d)\n", RadialBasisFunction );
    }

    return;

}


// elliptical version
void    Lgm_Vec_RBF_Psi2( Lgm_Vector *v, Lgm_Vector *v0, double ex, double ey, double ez, double *Psi, int RadialBasisFunction  ) {

    double  r, omr, omr2, omr4;
    double  x, y, z, x2, y2, z2, r2, f, g;

    x = v->x - v0->x;
    y = v->y - v0->y;
    z = v->z - v0->z;
    x2 = x*x; y2 = y*y; z2 = z*z;
    r2 = ex*x2 + ey*y2 + ez*z2;

    if ( RadialBasisFunction == LGM_RBF_GAUSSIAN ) {

        *Psi  = exp( -r2 );

    } else if ( RadialBasisFunction == LGM_RBF_MULTIQUADRIC ) {

        f   = 1.0 + r2;
        *Psi = sqrt( f );


    } else if ( RadialBasisFunction == LGM_RBF_INV_MULTIQUADRIC ) {

        f   = 1.0 + r2;
        *Psi = sqrt( f );


    } else if ( RadialBasisFunction == LGM_RBF_WENDLAND30 ) {
        
        r = sqrt(r2);
        omr = 1.0 - r; if ( omr < 0.0 ) omr = 0.0;
        *Psi = omr*omr;

    } else if ( RadialBasisFunction == LGM_RBF_WENDLAND31 ) {

        
        r = sqrt(r2);
        omr = 1.0 - r; if ( omr < 0.0 ) omr = 0.0;
        omr2 = omr*omr;
        *Psi = omr2*omr2*(4.0*r+1.0);

    } else if ( RadialBasisFunction == LGM_RBF_WENDLAND32 ) {
        
        r = sqrt(r2);
        omr = 1.0 - r; if ( omr < 0.0 ) omr = 0.0;
        omr2 = omr*omr;
        *Psi = omr2*omr2*omr2*(35.0*r2 + 18.0*r + 3.0);

    } else if ( RadialBasisFunction == LGM_RBF_WENDLAND33 ) {
        
        r = sqrt(r2);
        omr = 1.0 - r; if ( omr < 0.0 ) omr = 0.0;
        omr2 = omr*omr;
        omr4 = omr2*omr2;
        *Psi = omr4*omr4*(32.0*r2*r + 25.0*r2 + 8.0*r + 1.0);

    } else {
        printf("Lgm_Vec_RBF_Psi2: Unknown value for RadialBasisFunction. (Got %d)\n", RadialBasisFunction );
    }

    return;

}


void    Lgm_Vec_RBF_Derivs( Lgm_Vector *v, Lgm_Vector *v0, double e, double *dPdx, double *dPdy, double *dPdz, int RadialBasisFunction  ) {

    double  psi, x, y, z, x2, y2, z2, r2, f, g;

    x = v->x - v0->x;
    y = v->y - v0->y;
    z = v->z - v0->z;
    x2 = x*x; y2 = y*y; z2 = z*z;
    r2 = x2 + y2 + z2;

    if ( RadialBasisFunction == LGM_RBF_GAUSSIAN ) {

        psi  = exp( -e*r2 );
        g = -2.0*e*psi;
        *dPdx = g*x;
        *dPdy = g*y;
        *dPdz = g*z;

    } else if ( RadialBasisFunction == LGM_RBF_MULTIQUADRIC ) {

        f   = 1.0 + e*r2;
        psi = sqrt( f );
        g = e/psi;
        *dPdx = g*x;
        *dPdy = g*y;
        *dPdz = g*z;


    } else if ( RadialBasisFunction == LGM_RBF_INV_MULTIQUADRIC ) {

        f   = 1.0 + e*r2;
        psi = sqrt( f );
        g = -e/(f*psi);
        *dPdx = g*x;
        *dPdy = g*y;
        *dPdz = g*z;

    } else if ( RadialBasisFunction == LGM_RBF_WENDLAND30 ) {

        *dPdx = 0.0;
        *dPdy = 0.0;
        *dPdz = 0.0;

    } else if ( RadialBasisFunction == LGM_RBF_WENDLAND31 ) {

        *dPdx = 0.0;
        *dPdy = 0.0;
        *dPdz = 0.0;

    } else if ( RadialBasisFunction == LGM_RBF_WENDLAND32 ) {

        *dPdx = 0.0;
        *dPdy = 0.0;
        *dPdz = 0.0;

    } else if ( RadialBasisFunction == LGM_RBF_WENDLAND33 ) {

        *dPdx = 0.0;
        *dPdy = 0.0;
        *dPdz = 0.0;

    } else {
        printf("Lgm_Vec_RBF_Derivs: Unknown value for RadialBasisFunction. (Got %d)\n", RadialBasisFunction );
    }

    return;

}



// elliptical version.
void    Lgm_Vec_RBF_Derivs2( Lgm_Vector *v, Lgm_Vector *v0, double ex, double ey, double ez, double *dPdx, double *dPdy, double *dPdz, int RadialBasisFunction  ) {

    double  psi, x, y, z, x2, y2, z2, r2, f, g, omr, r;

    x = v->x - v0->x;
    y = v->y - v0->y;
    z = v->z - v0->z;
    x2 = x*x; y2 = y*y; z2 = z*z;
    r2 = ex*x2 + ey*y2 + ez*z2;

    if ( RadialBasisFunction == LGM_RBF_GAUSSIAN ) {

        psi  = exp( -r2 );
        g = -2.0*psi;
        *dPdx = ex*x*g;
        *dPdy = ey*y*g;
        *dPdz = ez*z*g;

    } else if ( RadialBasisFunction == LGM_RBF_MULTIQUADRIC ) {

        f   = 1.0 + r2;
        psi = sqrt( f );
        *dPdx = ex*x/psi;
        *dPdy = ey*y/psi;
        *dPdz = ez*z/psi;


    } else if ( RadialBasisFunction == LGM_RBF_INV_MULTIQUADRIC ) {

        f   = 1.0 + r2;
        psi = 1.0/sqrt( f );
        g = 1.0/(f*psi);
        *dPdx = -ex*x*g;
        *dPdy = -ey*y*g;
        *dPdz = -ez*z*g;

    } else if ( RadialBasisFunction == LGM_RBF_WENDLAND30 ) {

        r = sqrt(r2);
        omr = 1.0 - r; if ( omr < 0.0 ) omr = 0.0;
        psi = omr*omr;

        *dPdx = -2.0*ex*x*omr/r;
        *dPdy = -2.0*ey*y*omr/r;
        *dPdz = -2.0*ez*z*omr/r;

    } else if ( RadialBasisFunction == LGM_RBF_WENDLAND31 ) {

        *dPdx = 0.0;
        *dPdy = 0.0;
        *dPdz = 0.0;

    } else if ( RadialBasisFunction == LGM_RBF_WENDLAND32 ) {

        *dPdx = 0.0;
        *dPdy = 0.0;
        *dPdz = 0.0;

    } else if ( RadialBasisFunction == LGM_RBF_WENDLAND33 ) {

        *dPdx = 0.0;
        *dPdy = 0.0;
        *dPdz = 0.0;

    } else {
        printf("Lgm_Vec_RBF_Derivs2: Unknown value for RadialBasisFunction. (Got %d)\n", RadialBasisFunction );
    }

    return;

}


/*
 * D = 0 just gives a constant term
 * D = 1 gives a tri-linear polynomial
 * D = 2 gives a tri-quadratic polynomial
 *
 * M is returned as the number of terms computed.
 */
void PolyTerms( int D, int *M, double x, double y, double z, double *q ) {

    int     i, j, k;
    double  qx[3], qy[3], qz[3];
    /*
     * Set polynomial terms
     */
    qx[0] = 1.0; qx[1] = x; qx[2] = x*x;
    qy[0] = 1.0; qy[1] = y; qy[2] = y*y;
    qz[0] = 1.0; qz[1] = z; qz[2] = z*z;
    for ( *M=0, i=0; i<=D; i++ ) {
        for ( j=0; j<=D; j++ ) {
            for ( k=0; k<=D; k++ ) {
                q[*M] = qx[i]*qy[j]*qz[k];
                ++(*M);
            }
        }
    }

}



/** From a vector-field dataset, compute the vector-valued weighting factors,
 *  \f$\vec{c}_j\f$. Info is returned in the rbf structure.
 *
 *
 *  \param[in]                       v   -   pointer to an array of position vectors.
 *  \param[in]                       B   -   pointer to array of corresponding field vectors.
 *  \param[in]                     eps   -   smoothing factors in scalar RBF.
 *  \param[in]                       n   -   number of (v, B) pairs defined.
 *  \param[in]                   DoPoly  -   Flag to use simultaneous linear polynomial fit as well.
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
Lgm_Vec_RBF_Info *Lgm_Vec_RBF_Init( unsigned long int *I_data, Lgm_Vector *v, Lgm_Vector *B, double *eps_x, double *eps_y, double *eps_z, int n, int DoPoly, int RadialBasisFunction ) {

    unsigned long int Bytes;
    int              i, j, k, ii, jj, p, q, s, M;
    double           *d, **a, Psi, val, Qk[27];
    gsl_matrix       *A, *V, *AA, *A_new, *AA_new, *P_new;
    gsl_vector       *Dx, *Dy, *Dz, *cx, *cy, *cz, *S, *Work, *resid;
    gsl_vector       *Dx_new, *Dy_new, *Dz_new, *cx_new, *cy_new, *cz_new, *resid_new;
    gsl_permutation  *P;
    Lgm_Vec_RBF_Info *rbf;



    /*
     * Do all three components separately, but at the same time.
     */
    A  = gsl_matrix_calloc( n, n );
    AA = gsl_matrix_calloc( n, n );
    resid = gsl_vector_alloc( n );

    cx = gsl_vector_alloc( n );
    cy = gsl_vector_alloc( n );
    cz = gsl_vector_alloc( n );

    Dx = gsl_vector_calloc( n );
    Dy = gsl_vector_calloc( n );
    Dz = gsl_vector_calloc( n );

    if ( DoPoly ) {

        M = 27;
        M = 8;
        P_new  = gsl_matrix_calloc( M, n );

    } else {

        M = 0;

    }

    A_new  = gsl_matrix_calloc( n+M, n+M );
    AA_new = gsl_matrix_calloc( n+M, n+M );
    resid_new = gsl_vector_alloc( n+M );

    cx_new = gsl_vector_alloc( n+M );
    cy_new = gsl_vector_alloc( n+M );
    cz_new = gsl_vector_alloc( n+M );

    Dx_new = gsl_vector_calloc( n+M );
    Dy_new = gsl_vector_calloc( n+M );
    Dz_new = gsl_vector_calloc( n+M );


    /*
     * Save info needed to do an evaluation.
     */
    Bytes = sizeof(*rbf);
    rbf = ( Lgm_Vec_RBF_Info *)calloc( 1, Bytes );


    rbf->RadialBasisFunction = RadialBasisFunction;
    rbf->n      = n;
    rbf->DoPoly = DoPoly;

    LGM_ARRAY_1D( rbf->LookUpKey, n, unsigned long int);    Bytes += n*sizeof( unsigned long int );
    LGM_ARRAY_1D( rbf->v,     n, Lgm_Vector);               Bytes += n*sizeof( Lgm_Vector );
    LGM_ARRAY_1D( rbf->eps_x, n, double);                   Bytes += n*sizeof( double );
    LGM_ARRAY_1D( rbf->eps_y, n, double);                   Bytes += n*sizeof( double ); 
    LGM_ARRAY_1D( rbf->eps_z, n, double);                   Bytes += n*sizeof( double );
    LGM_ARRAY_1D( rbf->cx,    n, double);                   Bytes += n*sizeof( double );
    LGM_ARRAY_1D( rbf->cy,    n, double);                   Bytes += n*sizeof( double );
    LGM_ARRAY_1D( rbf->cz,    n, double);                   Bytes += n*sizeof( double );
    LGM_ARRAY_1D( rbf->cx_new,  n, double);                 Bytes += n*sizeof( double );
    LGM_ARRAY_1D( rbf->cy_new,  n, double);                 Bytes += n*sizeof( double );
    LGM_ARRAY_1D( rbf->cz_new,  n, double);                 Bytes += n*sizeof( double );

    if ( DoPoly ) {
        LGM_ARRAY_1D( rbf->gx_new,  M, double);                 Bytes += M*sizeof( double );
        LGM_ARRAY_1D( rbf->gy_new,  M, double);                 Bytes += M*sizeof( double );
        LGM_ARRAY_1D( rbf->gz_new,  M, double);                 Bytes += M*sizeof( double );
    }

    for ( i=0; i<n; i++ ) {
        rbf->LookUpKey[i] = I_data[i];
        rbf->v[i]         = v[i];
        rbf->eps_x[i]       = eps_x[i];
        rbf->eps_y[i]       = eps_y[i];
        rbf->eps_z[i]       = eps_z[i];
    }


    rbf->size  = (double)Bytes/1.0e6;


    // This subtraction doesntm seem to work out very well...?
//    rbf->Bx0 = B[n-1].x;
//    rbf->By0 = B[n-1].y;
//    rbf->Bz0 = B[n-1].z;
//    rbf->Bx0 = 0.0;
//    rbf->By0 = 0.0;
//    rbf->Bz0 = 0.0;

//double Bbkg;
//for ( Bbkg = 0.0, i=0; i<n; i++ ) Bbkg += B[i].x; rbf->Bx0 = Bbkg/(double)n;
//for ( Bbkg = 0.0, i=0; i<n; i++ ) Bbkg += B[i].y; rbf->By0 = Bbkg/(double)n;
//for ( Bbkg = 0.0, i=0; i<n; i++ ) Bbkg += B[i].z; rbf->Bz0 = Bbkg/(double)n;
//printf("Bavg = %g %g %g\n", rbf->Bx0, rbf->By0, rbf->Bz0 );
    rbf->Bx0 = 0.0;
    rbf->By0 = 0.0;
    rbf->Bz0 = 0.0;

//double min;
//for ( min=9e99, i=0; i<n; i++ ) { if ( B[i].x < min )  min = B[i].x; } rbf->Bx0 = min;
//for ( min=9e99, i=0; i<n; i++ ) { if ( B[i].y < min )  min = B[i].y; } rbf->By0 = min;
//for ( min=9e99, i=0; i<n; i++ ) { if ( B[i].z < min )  min = B[i].z; } rbf->Bz0 = min;


    
    /*
     * Fill D arrays. (Subtract off the field at the nearest point v[0] -- See
     * McNally [2011].) We add this field back on later.
     */
    for (i=0; i<n; i++){
        gsl_vector_set( Dx, i, B[i].x - rbf->Bx0 );
        gsl_vector_set( Dy, i, B[i].y - rbf->By0 );
        gsl_vector_set( Dz, i, B[i].z - rbf->Bz0 );
    }



    ii = 0;

    if ( DoPoly ) {
        for (i=0; i<M; i++){
            gsl_vector_set( Dx_new, ii, 0.0 );
            gsl_vector_set( Dy_new, ii, 0.0 );
            gsl_vector_set( Dz_new, ii, 0.0 );
            ++ii;
        }
    }

    for (i=0; i<n; i++){
        gsl_vector_set( Dx_new, ii, B[i].x - rbf->Bx0 );
        gsl_vector_set( Dy_new, ii, B[i].y - rbf->By0 );
        gsl_vector_set( Dz_new, ii, B[i].z - rbf->Bz0 );
        ++ii;
    }


    /*
     *                                                                         [  row0  ]
     * Fill A matrix (same for all 3 conponents). In C, order is A[row][col] = [  row1  ]
     *                                                                         [  row2  ]
     *
     *  NOTE: If the RBF shape parameter, epsilon is different for each RBF,
     *  the A matrix wont be symmetric.
     *
     */
    for ( i=0; i<n; i++ ) { // locate start row for subarray
        for ( j=0; j<n; j++ ) { // locate start column for subarray
            // Get PSi( v_i - v_j )
            Lgm_Vec_RBF_Psi2( &v[i], &v[j], eps_x[j], eps_y[j], eps_z[j], &Psi, rbf->RadialBasisFunction );
            gsl_matrix_set( A, i, j, Psi );
        }
    }




    if ( DoPoly ) {

        // P matrix
        for ( j=0; j<n; j++ ) {
            PolyTerms( 1, &M, v[j].x, v[j].y, v[j].z, Qk );
            for ( i=0; i<M; i++ ) { 
                gsl_matrix_set( P_new, i, j, Qk[i] );
                //gsl_matrix_set( P_new, j, i, Qk[i] );
            }
        }

        // Zero Block
        for ( i=0; i<M; i++ ) { 
            for ( j=0; j<M; j++ ) {
                gsl_matrix_set( A_new, i, j, 0.0 );
            }
        }

        // P Block
        for ( i=0; i<M; i++ ) { 
            for ( j=0; j<n; j++ ) {
                gsl_matrix_set( A_new, i, j+M, gsl_matrix_get(P_new, i, j) );
            }
        }

        // P^T Block
        for ( i=0; i<n; i++ ) { 
            for ( j=0; j<M; j++ ) {
                gsl_matrix_set( A_new, i+M, j, gsl_matrix_get(P_new, j, i) );
            }
        }

        gsl_matrix_free( P_new );

    }

    // Psi Block
    for ( i=0; i<n; i++ ) { 
        for ( j=0; j<n; j++ ) {
            gsl_matrix_set( A_new, i+M, j+M, gsl_matrix_get(A, i, j) );
        }
    }

    

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




//        printf("\nBefore LU solver\n");
//        for (i=0; i<n; i++ ) { printf("v%02d = %8g %8g %8g   B%02d = %8g %8g %8g\n", i, v[i].x, v[i].y, v[i].z, i, B[i].x, B[i].y, B[i].z ); }
//        printf("A=\n"); for (i=0; i<n; i++){ for (j=0; j<n; j++){ printf("%8g ", gsl_matrix_get(A, i, j ) ); } printf("\n"); } printf("\n");
//        printf("Dx= "); for (i=0; i<n; i++){ printf("%8g ", gsl_vector_get(Dx, i ) ); } printf("\n");

////        gsl_matrix_memcpy( AA, A );
////        P = gsl_permutation_alloc( n );
////        gsl_linalg_LU_decomp( A, P, &s );
////        gsl_linalg_LU_solve( A, P, Dx, cx ); // x-comp
//////        printf("cx= "); for (i=0; i<n; i++){ printf("%8g ", gsl_vector_get(cx, i ) ); } printf("\n");

////        for (ii=0; ii<5; ii++){
////            gsl_linalg_LU_refine(  AA, A, P, Dx, cx, resid );
//////            printf("\ncx= "); for (i=0; i<n; i++){ printf("%8g ", gsl_vector_get(cx, i ) ); } printf("\n");
//////            printf("rx= "); for (i=0; i<n; i++){ printf("%8g ", gsl_vector_get(rx, i ) ); } printf("\n");
////        }
////
////        gsl_linalg_LU_solve( A, P, Dy, cy ); // y-comp
////        for (ii=0; ii<5; ii++) gsl_linalg_LU_refine(  AA, A, P, Dy, cy, resid );
////
////        gsl_linalg_LU_solve( A, P, Dz, cz ); // z-comp
////        for (ii=0; ii<5; ii++) gsl_linalg_LU_refine(  AA, A, P, Dy, cy, resid );


//        printf("\nA*cx = ");
//double sum;
//        for (i=0; i<n; i++){ 
//            sum = 0.0;
//            for (j=0; j<n; j++){ 
//                sum += gsl_matrix_get(AA, i, j )*gsl_vector_get(cx, j );
//            }
//            printf("%8g ", sum ); 
//        } 
//        printf("\n");
//        printf("After LU solver\n\n");

////        gsl_permutation_free( P );




if (1==1){
        gsl_matrix_memcpy( AA_new, A_new );
//printf("A_new=\n"); for (i=0; i<n+M; i++){ for (j=0; j<n+M; j++){ printf("%8g ", gsl_matrix_get(A_new, i, j ) ); } printf("\n"); } printf("\n");
        P = gsl_permutation_alloc( n+M );
        gsl_linalg_LU_decomp( A_new, P, &s );
        gsl_linalg_LU_solve( A_new, P, Dx_new, cx_new ); // x-comp
        for (ii=0; ii<1; ii++) gsl_linalg_LU_refine(  AA_new, A_new, P, Dx_new, cx_new, resid_new );
//printf("cx_new = ");
//for(ii=0; ii<n+M; ii++){
//    printf(" %g", gsl_vector_get(cx_new, ii ) );
//}
//printf("\n");

        gsl_linalg_LU_solve( A_new, P, Dy_new, cy_new ); // y-comp
        for (ii=0; ii<1; ii++) gsl_linalg_LU_refine(  AA_new, A_new, P, Dy_new, cy_new, resid_new );

        gsl_linalg_LU_solve( A_new, P, Dz_new, cz_new ); // z-comp
        for (ii=0; ii<1; ii++) gsl_linalg_LU_refine(  AA_new, A_new, P, Dz_new, cz_new, resid_new );
        gsl_permutation_free( P );

        gsl_matrix_free( AA );
        gsl_matrix_free( AA_new );
}




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


    for (i=0; i<n; i++){
        rbf->cx_new[i] = gsl_vector_get( cx_new, i+M );
        rbf->cy_new[i] = gsl_vector_get( cy_new, i+M );
        rbf->cz_new[i] = gsl_vector_get( cz_new, i+M );
    }
    for (i=0; i<M; i++){
        rbf->gx_new[i] = gsl_vector_get( cx_new, i );
        rbf->gy_new[i] = gsl_vector_get( cy_new, i );
        rbf->gz_new[i] = gsl_vector_get( cz_new, i );
    }
//exit(0);

    
    gsl_vector_free( Dx );
    gsl_vector_free( Dy );
    gsl_vector_free( Dz );
    gsl_vector_free( resid );
    gsl_vector_free( cx );
    gsl_vector_free( cy );
    gsl_vector_free( cz );
    gsl_matrix_free( A );

    gsl_vector_free( Dx_new );
    gsl_vector_free( Dy_new );
    gsl_vector_free( Dz_new );
    gsl_vector_free( resid_new );
    gsl_vector_free( cx_new );
    gsl_vector_free( cy_new );
    gsl_vector_free( cz_new );
    gsl_matrix_free( A_new );

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
    LGM_ARRAY_1D_FREE( rbf->v   );
    LGM_ARRAY_1D_FREE( rbf->eps_x );
    LGM_ARRAY_1D_FREE( rbf->eps_y );
    LGM_ARRAY_1D_FREE( rbf->eps_z );
    LGM_ARRAY_1D_FREE( rbf->cx  );
    LGM_ARRAY_1D_FREE( rbf->cy  );
    LGM_ARRAY_1D_FREE( rbf->cz  );

    LGM_ARRAY_1D_FREE( rbf->cx_new  );
    LGM_ARRAY_1D_FREE( rbf->cy_new  );
    LGM_ARRAY_1D_FREE( rbf->cz_new  );

    if ( rbf->DoPoly ) {
        LGM_ARRAY_1D_FREE( rbf->gx_new  );
        LGM_ARRAY_1D_FREE( rbf->gy_new  );
        LGM_ARRAY_1D_FREE( rbf->gz_new  );
    }

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

    int         j, M;
    double      psi, Qk[27];


    B->x = B->y = B->z = 0.0;
/*
    for ( j=0; j<rbf->n; j++ ){
        Lgm_Vec_RBF_Psi( v, &rbf->v[j], rbf->eps[j], &psi, rbf->RadialBasisFunction );
        B->x += psi * rbf->cx[j];
        B->y += psi * rbf->cy[j];
        B->z += psi * rbf->cz[j];
    }
*/

    for ( j=0; j<rbf->n; j++ ){
        Lgm_Vec_RBF_Psi2( v, &rbf->v[j], rbf->eps_x[j], rbf->eps_y[j], rbf->eps_z[j], &psi, rbf->RadialBasisFunction );
        B->x += psi * rbf->cx_new[j];
        B->y += psi * rbf->cy_new[j];
        B->z += psi * rbf->cz_new[j];
    }

    if ( rbf->DoPoly ) {
        PolyTerms( 1, &M, v->x, v->y, v->z, Qk );
        for ( j=0; j<M; j++ ){
            B->x += rbf->gx_new[j]*Qk[j];
            B->y += rbf->gy_new[j]*Qk[j];
            B->z += rbf->gz_new[j]*Qk[j];
        }
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
/*
    for ( j=0; j<rbf->n; j++ ){

        Lgm_Vec_RBF_Derivs( v, &rbf->v[j], rbf->eps[j], &dPdx, &dPdy, &dPdz, rbf->RadialBasisFunction );

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
*/

    for ( j=0; j<rbf->n; j++ ){

        Lgm_Vec_RBF_Derivs2( v, &rbf->v[j], rbf->eps_x[j], rbf->eps_y[j], rbf->eps_z[j], &dPdx, &dPdy, &dPdz, rbf->RadialBasisFunction );

        dBdx->x += dPdx * rbf->cx_new[j]; 
        dBdx->y += dPdx * rbf->cy_new[j]; 
        dBdx->z += dPdx * rbf->cz_new[j]; 

        dBdy->x += dPdy * rbf->cx_new[j]; 
        dBdy->y += dPdy * rbf->cy_new[j]; 
        dBdy->z += dPdy * rbf->cz_new[j]; 

        dBdz->x += dPdz * rbf->cx_new[j];
        dBdz->y += dPdz * rbf->cy_new[j];
        dBdz->z += dPdz * rbf->cz_new[j];

    }

    if ( rbf->DoPoly ) {
        double x, y, z;
        x = v->x; y = v->y; z = v->z;
        //Poly term dB/dx
        dBdx->x += rbf->gx_new[4] + rbf->gx_new[5]*z + rbf->gx_new[6]*y + rbf->gx_new[7]*y*z;
        dBdx->y += rbf->gy_new[4] + rbf->gy_new[5]*z + rbf->gy_new[6]*y + rbf->gy_new[7]*y*z;
        dBdx->z += rbf->gz_new[4] + rbf->gz_new[5]*z + rbf->gz_new[6]*y + rbf->gz_new[7]*y*z;
        
        //Poly term dB/dy
        dBdy->x += rbf->gx_new[2] + rbf->gx_new[3]*z + rbf->gx_new[6]*x + rbf->gx_new[7]*x*z;
        dBdy->y += rbf->gy_new[2] + rbf->gy_new[3]*z + rbf->gy_new[6]*x + rbf->gy_new[7]*x*z;
        dBdy->z += rbf->gz_new[2] + rbf->gz_new[3]*z + rbf->gz_new[6]*x + rbf->gz_new[7]*x*z;
        
        //Poly term dB/dz
        dBdz->x += rbf->gx_new[1] + rbf->gx_new[3]*y + rbf->gx_new[5]*x + rbf->gx_new[7]*x*y;
        dBdz->y += rbf->gy_new[1] + rbf->gy_new[3]*y + rbf->gy_new[5]*x + rbf->gy_new[7]*x*y;
        dBdz->z += rbf->gz_new[1] + rbf->gz_new[3]*y + rbf->gz_new[5]*x + rbf->gz_new[7]*x*y;
    }


    return;
    
}


