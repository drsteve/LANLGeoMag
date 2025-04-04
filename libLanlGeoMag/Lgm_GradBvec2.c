#include <math.h>
#include <stdio.h>
#include "Lgm/Lgm_MagModelInfo.h"

/**
 *  \brief
 *      Compute the gradient of B-vec at a given point.
 *
 *  \details
 *      Similar to Lgm_GradBvec(), except this routine also computes time
 *      derivatives of \f$ \vec{B} \f$.  Computes \f$ {\nabla} \vec{B}, which
 *      is a tensor.\f$. To make the indices consistent with the Lgm_GradBvec()
 *      version, the time derivs are put in at the highest index.
 *
 *      Each element of the m array holds a Lgm_MagModelInfo structure that is
 *      initialized with the value of the desired discrete time array.  E.g.,
 *      if you wanted to run for a whole day at dt = 10s, then you need an
 *      array that is 8640 elements, and preinitialize the Lgm_MagModelInfo
 *      structure, for that index with that time. The derivatives are done with
 *      finite differences by using the mInfo structures that range over the
 *      appropriate indices. E.g., if j = 362 and you choose a six-point
 *      scheme, then time derivs will be based on indices 362-3 -> 362+3. Some
 *      things to make sure of that may or may not be checked: (1) Make sure
 *      the j-N > 0 and j+N < n; (2) Make sure the delta-Ts are the same and
 *      are the same as given in the dt argument. 
 *    
 *
 *
 *
 *      \param[in]      j           The current time index into the array of times that is implied by the m array.
 *      \param[in]      u0          Position (in GSM) to use.
 *      \param[out]     GradBvec    The computed gradient of Bvec at position u0, as well as the dB/dt vector
 *      \param[in]      DerivScheme Derivative scheme to use (can be one of LGM_DERIV_SIX_POINT, LGM_DERIV_FOUR_POINT, or LGM_DERIV_TWO_POINT).
 *      \param[in]      h           The delta (in Re) to use for grid spacing in the derivative scheme.
 *      \param[in]      n           The number of elements in the m array.
 *      \param[in]      dt          The Delta-time value used in the m array.
 *      \param[in,out]  m           An array of properly initialized and configured Lgm_MagModelInfo structures. 
 *
 *
 *      \author         Mike Henderson
 *      \date           2011
 *
 */
void Lgm_GradBvec2( int j, Lgm_Vector *u0, double GradBvec[4][4], int DerivScheme, double h, int n, double dt, Lgm_MagModelInfo **m ) {

    double      H, s, d, dd[7], D;
    double      ft[7], fx[7], fy[7], fz[7];
    double      gt[7], gx[7], gy[7], gz[7];
    double      qt[7], qx[7], qy[7], qz[7];
    int         i, N;
    double      dbx_dt, dbx_dx, dbx_dy, dbx_dz;
    double      dby_dt, dby_dx, dby_dy, dby_dz;
    double      dbz_dt, dbz_dx, dbz_dy, dbz_dz;
    Lgm_Vector  u, Bvec;


    /*
     * Select the derivative scheme to use.
     *      f_0^(1) = 1/(60h)  ( f_3 - 9f_2 + 45f_1  - 45f_-1 + 9f_-2 - f_-3 )
     * See page 450 of CRC standard Math tables 28th edition.
     */
    switch ( DerivScheme ) {
        case LGM_DERIV_SIX_POINT:
            N = 3;
            break;
        case LGM_DERIV_FOUR_POINT:
            N = 2;
            break;
        case LGM_DERIV_TWO_POINT:
            N = 1;
            break;
    }

    if ( ( j-N < 0 ) || ( j+N > n-1 ) ){
        printf( "Lgm_GradBvec2: Error, requested time index not in range. j = %d, N = %d, n = %d, j-N = %d\n", j, N, n, j-N);
        exit(0);
    }

    for (i=-N; i<=N; ++i) dd[i] = m[j+i]->c->TT.Time;
    for (i=-N; i<N; ++i) {
        D = (dd[i+1] - dd[i])*3600;
        if ( fabs( D-dt ) > 1e-12 ) {
            printf( "Lgm_GradBvec2: Error, times are inconsistent. D, dt = %.14lf %.14lf fabs(D-dt) = %g\n", D, dt, fabs(D-dt));
            exit(0);
        }
    }





    /*
     *  Compute all first order derivatives of magnetic field components.
     *  I.e., using subscript notation for derivs, we need these quantities:
     *     bx_t, bx_x, bx_y, bx_z
     *     by_t, by_x, by_y, by_z
     *     bz_t, bz_x, bz_y, bz_z
     */
    if (m[j]->VerbosityLevel > 0) printf("\t\tLgm_GradB: Computing GradB with DerivScheme = %d,  h = %g", DerivScheme, h);

    /*
     *  Bx Components
     */
    for (i=-N; i<=N; ++i){
        if ( N != 0 ) {
            m[j+i]->Bfield( u0, &Bvec, m[j+i] );     // variation in Time
            ft[i+N] = Bvec.x;
        }
    }
    for (i=-N; i<=N; ++i){
        if ( N != 0 ) {
            u = *u0; H = (double)i*h; u.x += H; // variation in X-COMP
            m[j]->Bfield( &u, &Bvec, m[j] );
            fx[i+N] = Bvec.x;
        }
    }
    for (i=-N; i<=N; ++i){
        if ( N != 0 ) {
            u = *u0; H = (double)i*h; u.y += H; // variation in Y-COMP
            m[j]->Bfield( &u, &Bvec, m[j] );
            fy[i+N] = Bvec.x;
        }
    }
    for (i=-N; i<=N; ++i){
        if ( N != 0 ) {
            u = *u0; H = (double)i*h; u.z += H; // variation in Z-COMP
            m[j]->Bfield( &u, &Bvec, m[j] );
            fz[i+N] = Bvec.x;
        }
    }

    /*
     *  By Components
     */
    for (i=-N; i<=N; ++i){
        if ( N != 0 ) {
            m[j+i]->Bfield( u0, &Bvec, m[j+i] );     // Variation in TIME
            gt[i+N] = Bvec.y;
        }
    }
    for (i=-N; i<=N; ++i){
        if ( N != 0 ) {
            u = *u0; H = (double)i*h; u.x += H; // variation in X-COMP
            m[j]->Bfield( &u, &Bvec, m[j] );
            gx[i+N] = Bvec.y;
        }
    }
    for (i=-N; i<=N; ++i){
        if ( N != 0 ) {
            u = *u0; H = (double)i*h; u.y += H; // variation in Y-COMP
            m[j]->Bfield( &u, &Bvec, m[j] );
            gy[i+N] = Bvec.y;
        }
    }
    for (i=-N; i<=N; ++i){
        if ( N != 0 ) {
            u = *u0; H = (double)i*h; u.z += H; // variation in Z-COMP
            m[j]->Bfield( &u, &Bvec, m[j] );
            gz[i+N] = Bvec.y;
        }
    }


    /*
     *  Bz Components
     */
    for (i=-N; i<=N; ++i){
        if ( N != 0 ) {
            m[j+i]->Bfield( u0, &Bvec, m[j+i] );     // variation in Time
            qt[i+N] = Bvec.z;
        }
    }
    for (i=-N; i<=N; ++i){
        if ( N != 0 ) {
            u = *u0; H = (double)i*h; u.x += H; // variation in X-COMP
            m[j]->Bfield( &u, &Bvec, m[j] );
            qx[i+N] = Bvec.z;
        }
    }
    for (i=-N; i<=N; ++i){
        if ( N != 0 ) {
            u = *u0; H = (double)i*h; u.y += H; // variation in Y-COMP
            m[j]->Bfield( &u, &Bvec, m[j] );
            qy[i+N] = Bvec.z;
        }
    }
    for (i=-N; i<=N; ++i){
        if ( N != 0 ) {
            u = *u0; H = (double)i*h; u.z += H; // variation in Z-COMP
            m[j]->Bfield( &u, &Bvec, m[j] );
            qz[i+N] = Bvec.z;
        }
    }
    

    if (DerivScheme == LGM_DERIV_SIX_POINT){

        d = 1.0/(60.0*dt);
        s = 1.0/(60.0*h);


        dbx_dt = (ft[6] - 9.0*ft[5] + 45.0*ft[4] - 45.0*ft[2] + 9.0*ft[1] - ft[0])*d;
        dbx_dx = (fx[6] - 9.0*fx[5] + 45.0*fx[4] - 45.0*fx[2] + 9.0*fx[1] - fx[0])*s;
        dbx_dy = (fy[6] - 9.0*fy[5] + 45.0*fy[4] - 45.0*fy[2] + 9.0*fy[1] - fy[0])*s;
        dbx_dz = (fz[6] - 9.0*fz[5] + 45.0*fz[4] - 45.0*fz[2] + 9.0*fz[1] - fz[0])*s;

        dby_dt = (gt[6] - 9.0*gt[5] + 45.0*gt[4] - 45.0*gt[2] + 9.0*gt[1] - gt[0])*d;
        dby_dx = (gx[6] - 9.0*gx[5] + 45.0*gx[4] - 45.0*gx[2] + 9.0*gx[1] - gx[0])*s;
        dby_dy = (gy[6] - 9.0*gy[5] + 45.0*gy[4] - 45.0*gy[2] + 9.0*gy[1] - gy[0])*s;
        dby_dz = (gz[6] - 9.0*gz[5] + 45.0*gz[4] - 45.0*gz[2] + 9.0*gz[1] - gz[0])*s;

        dbz_dt = (qt[6] - 9.0*qt[5] + 45.0*qt[4] - 45.0*qt[2] + 9.0*qt[1] - qt[0])*d;
        dbz_dx = (qx[6] - 9.0*qx[5] + 45.0*qx[4] - 45.0*qx[2] + 9.0*qx[1] - qx[0])*s;
        dbz_dy = (qy[6] - 9.0*qy[5] + 45.0*qy[4] - 45.0*qy[2] + 9.0*qy[1] - qy[0])*s;
        dbz_dz = (qz[6] - 9.0*qz[5] + 45.0*qz[4] - 45.0*qz[2] + 9.0*qz[1] - qz[0])*s;

    } else if (DerivScheme == LGM_DERIV_FOUR_POINT){

        d = 1.0/(12.0*dt);
        s = 1.0/(12.0*h);

        dbx_dt = (-ft[4] + 8.0*ft[3] - 8.0*ft[1] + ft[0])*d;
        dbx_dx = (-fx[4] + 8.0*fx[3] - 8.0*fx[1] + fx[0])*s;
        dbx_dy = (-fy[4] + 8.0*fy[3] - 8.0*fy[1] + fy[0])*s;
        dbx_dz = (-fz[4] + 8.0*fz[3] - 8.0*fz[1] + fz[0])*s;

        dby_dt = (-gt[4] + 8.0*gt[3] - 8.0*gt[1] + gt[0])*d;
        dby_dx = (-gx[4] + 8.0*gx[3] - 8.0*gx[1] + gx[0])*s;
        dby_dy = (-gy[4] + 8.0*gy[3] - 8.0*gy[1] + gy[0])*s;
        dby_dz = (-gz[4] + 8.0*gz[3] - 8.0*gz[1] + gz[0])*s;

        dbz_dt = (-qt[4] + 8.0*qt[3] - 8.0*qt[1] + qt[0])*d;
        dbz_dx = (-qx[4] + 8.0*qx[3] - 8.0*qx[1] + qx[0])*s;
        dbz_dy = (-qy[4] + 8.0*qy[3] - 8.0*qy[1] + qy[0])*s;
        dbz_dz = (-qz[4] + 8.0*qz[3] - 8.0*qz[1] + qz[0])*s;

    } else if (DerivScheme == LGM_DERIV_TWO_POINT){

        d = 1.0/(2.0*dt);
        s = 1.0/(2.0*h);

        dbx_dt = (ft[2] - ft[0])*d;
        dbx_dx = (fx[2] - fx[0])*s;
        dbx_dy = (fy[2] - fy[0])*s;
        dbx_dz = (fz[2] - fz[0])*s;

        dby_dt = (gt[2] - gt[0])*d;
        dby_dx = (gx[2] - gx[0])*s;
        dby_dy = (gy[2] - gy[0])*s;
        dby_dz = (gz[2] - gz[0])*s;

        dbz_dt = (qt[2] - qt[0])*d;
        dbz_dx = (qx[2] - qx[0])*s;
        dbz_dy = (qy[2] - qy[0])*s;
        dbz_dz = (qz[2] - qz[0])*s;

    }
    if (m[j]->VerbosityLevel > 0) {
        printf("             = (%8g %8g %8g)\n", dbx_dx, dbx_dy, dbx_dz );
        printf("   Grad Bvec = (%8g %8g %8g)\n", dby_dx, dby_dy, dby_dz );
        printf("             = (%8g %8g %8g)\n", dbz_dx, dbz_dy, dbz_dz );
    }

    GradBvec[0][0] =    1.0;  GradBvec[0][1] = 0.0;     GradBvec[0][2] = 0.0;     GradBvec[0][3] = 0.0; 
    GradBvec[1][0] = dbx_dt;  GradBvec[1][1] = dbx_dx;  GradBvec[1][2] = dbx_dy;  GradBvec[1][3] = dbx_dz;
    GradBvec[2][0] = dby_dt;  GradBvec[2][1] = dby_dx;  GradBvec[2][2] = dby_dy;  GradBvec[2][3] = dby_dz;
    GradBvec[3][0] = dbz_dt;  GradBvec[3][1] = dbz_dx;  GradBvec[3][2] = dbz_dy;  GradBvec[3][3] = dbz_dz;




    return;

}
