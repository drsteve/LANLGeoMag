#include <stdio.h>
#include <math.h>
#include "Lgm/Lgm_MagModelInfo.h"

/**
 *  \brief
 *      Compute the gradient of B-vec at a given point.
 *
 *  \details
 *      Computes \f$ \nabla B \f$.
 *
 *
 *      \param[in]      u0          Position (in GSM) to use.
 *      \param[out]     GradBvec    The computed gradient of Bvec at position u0 -- note: this is a tensor
 *      \param[in]      DerivScheme Derivative scheme to use (can be one of LGM_DERIV_SIX_POINT, LGM_DERIV_FOUR_POINT, or LGM_DERIV_TWO_POINT).
 *      \param[in]      h           The delta (in Re) to use for grid spacing in the derivative scheme.
 *      \param[in,out]  m           A properly initialized and configured Lgm_MagModelInfo structure.
 *
 *
 *      \author         Mike Henderson
 *      \date           2023
 *
 */
void Lgm_GradBvec( Lgm_Vector *u0, double GradBvec[3][3], int DerivScheme, double h, Lgm_MagModelInfo *m ) {

    double      H, s;
    double      fx[7], fy[7], fz[7];
    double      gx[7], gy[7], gz[7];
    double      qx[7], qy[7], qz[7];
    int         i, N;
    double      dbx_dx, dbx_dy, dbx_dz;
    double      dby_dx, dby_dy, dby_dz;
    double      dbz_dx, dbz_dy, dbz_dz;
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


    /*
     *  Compute all first order derivatives of magnetic field components.
     *  I.e., using subscript notation for derivs, we need these quantities:
     *     bx_t, bx_x, bx_y, bx_z
     *     by_t, by_x, by_y, by_z
     *     bz_t, bz_x, bz_y, bz_z
     */
    if (m->VerbosityLevel > 0) printf("\t\tLgm_GradB: Computing GradB with DerivScheme = %d,  h = %g", DerivScheme, h);

    /*
     *  Bx Components
     */
    for (i=-N; i<=N; ++i){
        if ( N != 0 ) {
            u = *u0; H = (double)i*h; u.x += H; // variation in X-COMP
            m->Bfield( &u, &Bvec, m );
            fx[i+N] = Bvec.x;
        }
    }
    for (i=-N; i<=N; ++i){
        if ( N != 0 ) {
            u = *u0; H = (double)i*h; u.y += H; // variation in Y-COMP
            m->Bfield( &u, &Bvec, m );
            fy[i+N] = Bvec.x;
        }
    }
    for (i=-N; i<=N; ++i){
        if ( N != 0 ) {
            u = *u0; H = (double)i*h; u.z += H; // variation in Z-COMP
            m->Bfield( &u, &Bvec, m );
            fz[i+N] = Bvec.x;
        }
    }

    /*
     *  By Components
     */
    for (i=-N; i<=N; ++i){
        if ( N != 0 ) {
            u = *u0; H = (double)i*h; u.x += H; // variation in X-COMP
            m->Bfield( &u, &Bvec, m );
            gx[i+N] = Bvec.y;
        }
    }
    for (i=-N; i<=N; ++i){
        if ( N != 0 ) {
            u = *u0; H = (double)i*h; u.y += H; // variation in Y-COMP
            m->Bfield( &u, &Bvec, m );
            gy[i+N] = Bvec.y;
        }
    }
    for (i=-N; i<=N; ++i){
        if ( N != 0 ) {
            u = *u0; H = (double)i*h; u.z += H; // variation in Z-COMP
            m->Bfield( &u, &Bvec, m );
            gz[i+N] = Bvec.y;
        }
    }


    /*
     *  Bz Components
     */
    for (i=-N; i<=N; ++i){
        if ( N != 0 ) {
            u = *u0; H = (double)i*h; u.x += H; // variation in X-COMP
            m->Bfield( &u, &Bvec, m );
            qx[i+N] = Bvec.z;
        }
    }
    for (i=-N; i<=N; ++i){
        if ( N != 0 ) {
            u = *u0; H = (double)i*h; u.y += H; // variation in Y-COMP
            m->Bfield( &u, &Bvec, m );
            qy[i+N] = Bvec.z;
        }
    }printf("             = (%8g %8g %8g)\n", dbx_dx, dbx_dy, dbx_dz );
        printf("   Grad Bvec = (%8g %8g %8g)\n", dby_dx, dby_dy, dby_dz );
        printf("             = (%8g %8g %8g)\n", dbz_dx, dbz_dy, dbz_dz );
    for (i=-N; i<=N; ++i){
        if ( N != 0 ) {
            u = *u0; H = (double)i*h; u.z += H; // variation in Z-COMP
            m->Bfield( &u, &Bvec, m );
            qz[i+N] = Bvec.z;
        }
    }
    

    if (DerivScheme == LGM_DERIV_SIX_POINT){

        s = 1.0/(60.0*h);

        dbx_dx = (fx[6] - 9.0*fx[5] + 45.0*fx[4] - 45.0*fx[2] + 9.0*fx[1] - fx[0])*s;
        dbx_dy = (fy[6] - 9.0*fy[5] + 45.0*fy[4] - 45.0*fy[2] + 9.0*fy[1] - fy[0])*s;
        dbx_dz = (fz[6] - 9.0*fz[5] + 45.0*fz[4] - 45.0*fz[2] + 9.0*fz[1] - fz[0])*s;

        dby_dx = (gx[6] - 9.0*gx[5] + 45.0*gx[4] - 45.0*gx[2] + 9.0*gx[1] - gx[0])*s;
        dby_dy = (gy[6] - 9.0*gy[5] + 45.0*gy[4] - 45.0*gy[2] + 9.0*gy[1] - gy[0])*s;
        dby_dz = (gz[6] - 9.0*gz[5] + 45.0*gz[4] - 45.0*gz[2] + 9.0*gz[1] - gz[0])*s;

        dbz_dx = (qx[6] - 9.0*qx[5] + 45.0*qx[4] - 45.0*qx[2] + 9.0*qx[1] - qx[0])*s;
        dbz_dy = (qy[6] - 9.0*qy[5] + 45.0*qy[4] - 45.0*qy[2] + 9.0*qy[1] - qy[0])*s;
        dbz_dz = (qz[6] - 9.0*qz[5] + 45.0*qz[4] - 45.0*qz[2] + 9.0*qz[1] - qz[0])*s;

    } else if (DerivScheme == LGM_DERIV_FOUR_POINT){

        s = 1.0/(12.0*h);

        dbx_dx = (-fx[4] + 8.0*fx[3] - 8.0*fx[1] + fx[0])/(12.0*h);
        dbx_dy = (-fy[4] + 8.0*fy[3] - 8.0*fy[1] + fy[0])/(12.0*h);
        dbx_dz = (-fz[4] + 8.0*fz[3] - 8.0*fz[1] + fz[0])/(12.0*h);

        dby_dx = (-gx[4] + 8.0*gx[3] - 8.0*gx[1] + gx[0])*s;
        dby_dy = (-gy[4] + 8.0*gy[3] - 8.0*gy[1] + gy[0])*s;
        dby_dz = (-gz[4] + 8.0*gz[3] - 8.0*gz[1] + gz[0])*s;

        dbz_dx = (-qx[4] + 8.0*qx[3] - 8.0*qx[1] + qx[0])*s;
        dbz_dy = (-qy[4] + 8.0*qy[3] - 8.0*qy[1] + qy[0])*s;
        dbz_dz = (-qz[4] + 8.0*qz[3] - 8.0*qz[1] + qz[0])*s;

    } else if (DerivScheme == LGM_DERIV_TWO_POINT){

        s = 1.0/(2.0*h);
        dbx_dx = (fx[2] - fx[0])*s;
        dbx_dy = (fy[2] - fy[0])*s;
        dbx_dz = (fz[2] - fz[0])*s;

        dby_dx = (gx[2] - gx[0])*s;
        dby_dy = (gy[2] - gy[0])*s;
        dby_dz = (gz[2] - gz[0])*s;

        dbz_dx = (qx[2] - qx[0])*s;
        dbz_dy = (qy[2] - qy[0])*s;
        dbz_dz = (qz[2] - qz[0])*s;

    }
    if (m->VerbosityLevel > 0) {
        printf("             = (%8g %8g %8g)\n", dbx_dx, dbx_dy, dbx_dz );
        printf("   Grad Bvec = (%8g %8g %8g)\n", dby_dx, dby_dy, dby_dz );
        printf("             = (%8g %8g %8g)\n", dbz_dx, dbz_dy, dbz_dz );
    }

    GradBvec[0][0] = dbx_dx;  GradBvec[0][1] = dbx_dy;  GradBvec[0][2] = dbx_dz;
    GradBvec[1][0] = dby_dx;  GradBvec[1][1] = dby_dy;  GradBvec[1][2] = dby_dz;
    GradBvec[2][0] = dbz_dx;  GradBvec[2][1] = dbz_dy;  GradBvec[2][2] = dbz_dz;

    return;

}
