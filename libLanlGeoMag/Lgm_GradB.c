#include <stdio.h>
#include <math.h>
#include "Lgm/Lgm_MagModelInfo.h"


/**
 *  \brief
 *      Compute the gradient of |B-vec| at a given point.
 *
 *  \details
 *      Computes \f$ \nabla B \f$.
 *
 *
 *      \param[in]      u0          Position (in GSM) to use.
 *      \param[out]     GradB       The computed gradient of B at position u0.
 *      \param[in]      DerivScheme Derivative scheme to use (can be one of LGM_DERIV_SIX_POINT, LGM_DERIV_FOUR_POINT, or LGM_DERIV_TWO_POINT).
 *      \param[in]      h           The delta (in Re) to use for grid spacing in the derivative scheme.
 *      \param[in,out]  m           A properly initialized and configured Lgm_MagModelInfo structure.
 *
 *
 *      \author         Mike Henderson
 *      \date           2011
 *
 */
void Lgm_GradB( Lgm_Vector *u0, Lgm_Vector *GradB, int DerivScheme, double h, Lgm_MagModelInfo *m ) {

    double      B, H, fx[7], fy[7], fz[7];
    int         i, N;
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

    if (m->VerbosityLevel > 0) printf("\t\tLgm_GradB: Computing GradB with DerivScheme = %d,  h = %g", DerivScheme, h);
    for (i=-N; i<=N; ++i){
        if ( N != 0 ) {
            u = *u0; H = (double)i*h; u.x += H;
            m->Bfield( &u, &Bvec, m );
            B = Lgm_Magnitude( &Bvec );
            fx[i+N] = B;
        }
    }
    for (i=-N; i<=N; ++i){
        if ( N != 0 ) {
            u = *u0; H = (double)i*h; u.y += H;
            m->Bfield( &u, &Bvec, m );
            B = Lgm_Magnitude( &Bvec );
            fy[i+N] = B;
        }
    }
    for (i=-N; i<=N; ++i){
        if ( N != 0 ) {
            u = *u0; H = (double)i*h; u.z += H;
            m->Bfield( &u, &Bvec, m );
            B = Lgm_Magnitude( &Bvec );
            fz[i+N] = B;
        }
    }


    if (DerivScheme == LGM_DERIV_SIX_POINT){
        GradB->x = (fx[6] - 9.0*fx[5] + 45.0*fx[4] - 45.0*fx[2] + 9.0*fx[1] - fx[0])/(60.0*h);
        GradB->y = (fy[6] - 9.0*fy[5] + 45.0*fy[4] - 45.0*fy[2] + 9.0*fy[1] - fy[0])/(60.0*h);
        GradB->z = (fz[6] - 9.0*fz[5] + 45.0*fz[4] - 45.0*fz[2] + 9.0*fz[1] - fz[0])/(60.0*h);
    } else if (DerivScheme == LGM_DERIV_FOUR_POINT){
        GradB->x = (-fx[4] + 8.0*fx[3] - 8.0*fx[1] + fx[0])/(12.0*h);
        GradB->y = (-fy[4] + 8.0*fy[3] - 8.0*fy[1] + fy[0])/(12.0*h);
        GradB->z = (-fz[4] + 8.0*fz[3] - 8.0*fz[1] + fz[0])/(12.0*h);
    } else if (DerivScheme == LGM_DERIV_TWO_POINT){
        GradB->x = (fx[2] - fx[0])/(2.0*h);
        GradB->y = (fy[2] - fy[0])/(2.0*h);
        GradB->z = (fz[2] - fz[0])/(2.0*h);
    }
    if (m->VerbosityLevel > 0) printf("   GradB = (%g %g %g)\n", GradB->x, GradB->y, GradB->z );

    return;

}


/**
 *  \brief
 *      Compute the d|B|/dcomp at a given point.
 *
 *  \details
 *      Computes \f$ \partial B/\partial comp \f$.
 *
 *
 *      \param[in]      u0          Position (in GSM) to use.
 *      \param[in]      comp        seclects component (0=x, 1=y, 2=z)
 *      \param[out]     dBdx        The computed partial deriviative dB/dx at position u0.
 *      \param[in]      DerivScheme Derivative scheme to use (can be one of LGM_DERIV_SIX_POINT, LGM_DERIV_FOUR_POINT, or LGM_DERIV_TWO_POINT).
 *      \param[in]      h           The delta (in Re) to use for grid spacing in the derivative scheme.
 *      \param[in,out]  m           A properly initialized and configured Lgm_MagModelInfo structure.
 *
 *
 *      \author         Mike Henderson
 *      \date           2023
 *
 */
void Lgm_dBdcomp( Lgm_Vector *u0, int comp, double *dBdcomp, int DerivScheme, double h, Lgm_MagModelInfo *m ) {

    double      B, H, fc[7];
    int         i, N;
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

    if (m->VerbosityLevel > 0) printf("\t\tLgm_dBdcomp: Computing dBdcomp for comp=%d with DerivScheme = %d,  h = %g", comp, DerivScheme, h);
    for (i=-N; i<=N; ++i){
        if ( N != 0 ) {
            u = *u0; H = (double)i*h; 
            if      ( comp == 0 ) { u.x += H; }
            else if ( comp == 1 ) { u.y += H; }
            else if ( comp == 2 ) { u.z += H; }
            else { printf("\t\tLgm_dBdcomp: Invalid component. comp = %d (must be 0, 1, or 2).\n", comp); exit(0); }
            m->Bfield( &u, &Bvec, m );
            B = Lgm_Magnitude( &Bvec );
            fc[i+N] = B;
        }
    }

    if (DerivScheme == LGM_DERIV_SIX_POINT){
        *dBdcomp = (fc[6] - 9.0*fc[5] + 45.0*fc[4] - 45.0*fc[2] + 9.0*fc[1] - fc[0])/(60.0*h);
    } else if (DerivScheme == LGM_DERIV_FOUR_POINT){
        *dBdcomp = (-fc[4] + 8.0*fc[3] - 8.0*fc[1] + fc[0])/(12.0*h);
    } else if (DerivScheme == LGM_DERIV_TWO_POINT){
        *dBdcomp = (fc[2] - fc[0])/(2.0*h);
    }
    if (m->VerbosityLevel > 0) printf("   dBdcomp = %g (component = %d)\n", *dBdcomp, comp);

    return;

}
void Lgm_dBdx( Lgm_Vector *u0, double *dBdx, int DerivScheme, double h, Lgm_MagModelInfo *m ) { Lgm_dBdcomp( u0, 0, dBdx, DerivScheme, h, m ); }
void Lgm_dBdy( Lgm_Vector *u0, double *dBdy, int DerivScheme, double h, Lgm_MagModelInfo *m ) { Lgm_dBdcomp( u0, 1, dBdy, DerivScheme, h, m ); }
void Lgm_dBdz( Lgm_Vector *u0, double *dBdz, int DerivScheme, double h, Lgm_MagModelInfo *m ) { Lgm_dBdcomp( u0, 2, dBdz, DerivScheme, h, m ); }


/**
 *  \brief
 *      Compute the gradient of the diverence of b-hat at a given point.
 *
 *  \details
 *      Computes \f$ \nabla \nabla\cdot \hat{b} \f$.
 *
 *
 *      \param[in]      u0          Position (in GSM) to use.
 *      \param[out]     Grad_Divb   The computed gradient of Divb at position u0.
 *      \param[in]      DerivScheme Derivative scheme to use (can be one of LGM_DERIV_SIX_POINT, LGM_DERIV_FOUR_POINT, or LGM_DERIV_TWO_POINT).
 *      \param[in]      h           The delta (in Re) to use for grid spacing in the derivative scheme.
 *      \param[in,out]  m           A properly initialized and configured Lgm_MagModelInfo structure.
 *
 *
 *      \author         Mike Henderson
 *      \date           2023
 *
 */
void Lgm_Grad_Divb( Lgm_Vector *u0, Lgm_Vector *Grad_Divb, int DerivScheme, double h, Lgm_MagModelInfo *m ) {

    double      Divb, H, fx[7], fy[7], fz[7];
    int         i, N;
    Lgm_Vector  u;



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

    if (m->VerbosityLevel > 0) printf("\t\tLgm_GradB: Computing GradB with DerivScheme = %d,  h = %g", DerivScheme, h);
    for (i=-N; i<=N; ++i){
        if ( N != 0 ) {
            u = *u0; H = (double)i*h; u.x += H;
            Lgm_Divb( &u, &Divb, DerivScheme, h, m );
            fx[i+N] = Divb;
        }
    }
    for (i=-N; i<=N; ++i){
        if ( N != 0 ) {
            u = *u0; H = (double)i*h; u.y += H;
            Lgm_Divb( &u, &Divb, DerivScheme, h, m );
            fy[i+N] = Divb;
        }
    }
    for (i=-N; i<=N; ++i){
        if ( N != 0 ) {
            u = *u0; H = (double)i*h; u.z += H;
            Lgm_Divb( &u, &Divb, DerivScheme, h, m );
            fz[i+N] = Divb;
        }
    }


    if (DerivScheme == LGM_DERIV_SIX_POINT){
        Grad_Divb->x = (fx[6] - 9.0*fx[5] + 45.0*fx[4] - 45.0*fx[2] + 9.0*fx[1] - fx[0])/(60.0*h);
        Grad_Divb->y = (fy[6] - 9.0*fy[5] + 45.0*fy[4] - 45.0*fy[2] + 9.0*fy[1] - fy[0])/(60.0*h);
        Grad_Divb->z = (fz[6] - 9.0*fz[5] + 45.0*fz[4] - 45.0*fz[2] + 9.0*fz[1] - fz[0])/(60.0*h);
    } else if (DerivScheme == LGM_DERIV_FOUR_POINT){
        Grad_Divb->x = (-fx[4] + 8.0*fx[3] - 8.0*fx[1] + fx[0])/(12.0*h);
        Grad_Divb->y = (-fy[4] + 8.0*fy[3] - 8.0*fy[1] + fy[0])/(12.0*h);
        Grad_Divb->z = (-fz[4] + 8.0*fz[3] - 8.0*fz[1] + fz[0])/(12.0*h);
    } else if (DerivScheme == LGM_DERIV_TWO_POINT){
        Grad_Divb->x = (fx[2] - fx[0])/(2.0*h);
        Grad_Divb->y = (fy[2] - fy[0])/(2.0*h);
        Grad_Divb->z = (fz[2] - fz[0])/(2.0*h);
    }
    if (m->VerbosityLevel > 0) printf("   Grad_Divb = (%g %g %g)\n", Grad_Divb->x, Grad_Divb->y, Grad_Divb->z );

    return;

}


/**
 *  \brief
 *      Compute the Laplacian of |B| at a given point.
 *
 *  \details
 *      Computes \f$ \nabla^2 |B| \f$.
 *
 *
 *      \param[in]      u0          Position (in GSM) to use.
 *      \param[out]     LaplacianB  The computed gradient of Divb at position u0.
 *      \param[in]      DerivScheme Derivative scheme to use (can be one of LGM_DERIV_SIX_POINT, LGM_DERIV_FOUR_POINT, or LGM_DERIV_TWO_POINT).
 *      \param[in]      h           The delta (in Re) to use for grid spacing in the derivative scheme.
 *      \param[in,out]  m           A properly initialized and configured Lgm_MagModelInfo structure.
 *
 *
 *      \author         Mike Henderson
 *      \date           2023
 *
 */
void Lgm_LaplacianB( Lgm_Vector *u0, double *LaplacianB, int DerivScheme, double h, Lgm_MagModelInfo *m ) {

    double      H, d2Bdx2, d2Bdy2, d2Bdz2;
    double      dBdx, dBdy, dBdz;
    double      fx[7], fy[7], fz[7];
    double      gx[7], gy[7], gz[7];
    double      hx[7], hy[7], hz[7];
    int         i, N;
    Lgm_Vector  u;



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

    if (m->VerbosityLevel > 0) printf("\t\tLgm_GradB: Computing GradB with DerivScheme = %d,  h = %g", DerivScheme, h);
    for (i=-N; i<=N; ++i){
        if ( N != 0 ) {
            u = *u0; H = (double)i*h; u.x += H;
            Lgm_dBdx( &u, &dBdx, DerivScheme, h, m );
            fx[i+N] = dBdx;
        }
    }
    for (i=-N; i<=N; ++i){
        if ( N != 0 ) {
            u = *u0; H = (double)i*h; u.y += H;
            Lgm_dBdy( &u, &dBdy, DerivScheme, h, m );
            gy[i+N] = dBdy;
        }
    }
    for (i=-N; i<=N; ++i){
        if ( N != 0 ) {
            u = *u0; H = (double)i*h; u.z += H;
            Lgm_dBdz( &u, &dBdz, DerivScheme, h, m );
            hz[i+N] = dBdz;
        }
    }


    if (DerivScheme == LGM_DERIV_SIX_POINT){
        d2Bdx2 = (fx[6] - 9.0*fx[5] + 45.0*fx[4] - 45.0*fx[2] + 9.0*fx[1] - fx[0])/(60.0*h);
        d2Bdy2 = (gy[6] - 9.0*gy[5] + 45.0*gy[4] - 45.0*gy[2] + 9.0*gy[1] - gy[0])/(60.0*h);
        d2Bdz2 = (hz[6] - 9.0*hz[5] + 45.0*hz[4] - 45.0*hz[2] + 9.0*hz[1] - hz[0])/(60.0*h);
    } else if (DerivScheme == LGM_DERIV_FOUR_POINT){
        d2Bdx2 = (-fx[4] + 8.0*fx[3] - 8.0*fx[1] + fx[0])/(12.0*h);
        d2Bdy2 = (-gy[4] + 8.0*gy[3] - 8.0*gy[1] + gy[0])/(12.0*h);
        d2Bdz2 = (-hz[4] + 8.0*hz[3] - 8.0*hz[1] + hz[0])/(12.0*h);
    } else if (DerivScheme == LGM_DERIV_TWO_POINT){
        d2Bdx2 = (hz[2] - hz[0])/(2.0*h);
        d2Bdy2 = (hz[2] - hz[0])/(2.0*h);
        d2Bdz2 = (hz[2] - hz[0])/(2.0*h);
    }
    *LaplacianB = d2Bdx2 + d2Bdy2 + d2Bdz2;
    if (m->VerbosityLevel > 0) printf("   LaplacianB = %g\n", *LaplacianB );

    return;

}


/**
 *  \brief
 *      Compute the gradient of B and it's parallel and perpendicular components at a given point.
 *
 *  \details
 *      Computes \f$ \nabla B \f$, \f$ (\nabla B)_\parallel \f$, and \f$ (\nabla B)_\perp \f$.
 *
 *
 *      \param[in]      u0          Position (in GSM) to use.
 *      \param[out]     GradB       The computed curl of B at position u0.
 *      \param[out]     GradB_para  The computed parallel component of curl of B at position u0.
 *      \param[out]     GradB_perp  The computed perpendicular component of curl of B at position u0.
 *      \param[in]      DerivScheme Derivative scheme to use (can be one of LGM_DERIV_SIX_POINT, LGM_DERIV_FOUR_POINT, or LGM_DERIV_TWO_POINT).
 *      \param[in]      h           The delta (in Re) to use for grid spacing in the derivative scheme.
 *      \param[in,out]  m           A properly initialized and configured Lgm_MagModelInfo structure.
 *
 *
 *      \author         Mike Henderson
 *      \date           2011
 *
 */
void Lgm_GradB2( Lgm_Vector *u0, Lgm_Vector *GradB, Lgm_Vector *GradB_para, Lgm_Vector *GradB_perp, int DerivScheme, double h, Lgm_MagModelInfo *m ) {

    double      g;
    Lgm_Vector  Bvec;


    Lgm_GradB( u0, GradB, DerivScheme, h, m );


    m->Bfield( u0, &Bvec, m );
    Lgm_NormalizeVector( &Bvec );

    // Compute parallel component of GradB
    g = Lgm_DotProduct( &Bvec, GradB );
    *GradB_para = *GradB;
    Lgm_ScaleVector( GradB_para, g );


    // Compute perpendicular component of GradB (GradB_perp = GradB - GradB_para)
    Lgm_VecSub( GradB_perp, GradB, GradB_para );


    return;

}


void Lgm_B_Cross_GradB_Over_B( Lgm_Vector *u0, Lgm_Vector *A, int DerivScheme, double h, Lgm_MagModelInfo *m ) {

    double      B;
    Lgm_Vector  GradB, Bvec;

    m->Bfield( u0, &Bvec, m );
    B = Lgm_Magnitude( &Bvec );

    Lgm_GradB( u0, &GradB, DerivScheme, h, m );
    Lgm_CrossProduct( &Bvec, &GradB, A );
    Lgm_ScaleVector( A, 1.0/B );
    
}


/**
 *  \brief
 *      Compute the curl of B-vec at a given point.
 *
 *  \details
 *      Computes \f$ \nabla\times B \f$.
 *
 *
 *      \param[in]      u0          Position (in GSM) to use.
 *      \param[out]     CurlB       The computed curl of B at position u0.
 *      \param[in]      DerivScheme Derivative scheme to use (can be one of LGM_DERIV_SIX_POINT, LGM_DERIV_FOUR_POINT, or LGM_DERIV_TWO_POINT).
 *      \param[in]      h           The delta (in Re) to use for grid spacing in the derivative scheme.
 *      \param[in,out]  m           A properly initialized and configured Lgm_MagModelInfo structure.
 *
 *
 *      \author         Mike Henderson
 *      \date           2011
 *
 */
void Lgm_CurlB( Lgm_Vector *u0, Lgm_Vector *CurlB, int DerivScheme, double h, Lgm_MagModelInfo *m ) {

    double      B, H;
    double      fxy[7], fxz[7], fyx[7], fyz[7], fzx[7], fzy[7];
    double      dBxdy, dBxdz, dBydx, dBydz, dBzdx, dBzdy;
    int         i, N;
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

    if (m->VerbosityLevel > 0) printf("\t\tLgm_CurlB: Computing CurlB with DerivScheme = %d,  h = %g", DerivScheme, h);
    // dBx/dy and dBx/dz
    for (i=-N; i<=N; ++i){
        if ( N != 0 ) {
            u = *u0; H = (double)i*h; u.y += H;
            m->Bfield( &u, &Bvec, m );
            fxy[i+N] = Bvec.x;
        }
    }
    for (i=-N; i<=N; ++i){
        if ( N != 0 ) {
            u = *u0; H = (double)i*h; u.z += H;
            m->Bfield( &u, &Bvec, m );
            fxz[i+N] = Bvec.x;
        }
    }

    // dBy/dx and dBy/dz
    for (i=-N; i<=N; ++i){
        if ( N != 0 ) {
            u = *u0; H = (double)i*h; u.x += H;
            m->Bfield( &u, &Bvec, m );
            fyx[i+N] = Bvec.y;
        }
    }
    for (i=-N; i<=N; ++i){
        if ( N != 0 ) {
            u = *u0; H = (double)i*h; u.z += H;
            m->Bfield( &u, &Bvec, m );
            fyz[i+N] = Bvec.y;
        }
    }

    // dBz/dx and dBz/dy
    for (i=-N; i<=N; ++i){
        if ( N != 0 ) {
            u = *u0; H = (double)i*h; u.x += H;
            m->Bfield( &u, &Bvec, m );
            fzx[i+N] = Bvec.z;
        }
    }
    for (i=-N; i<=N; ++i){
        if ( N != 0 ) {
            u = *u0; H = (double)i*h; u.y += H;
            m->Bfield( &u, &Bvec, m );
            fzy[i+N] = Bvec.z;
        }
    }

    if (DerivScheme == LGM_DERIV_SIX_POINT){
        dBxdy = (fxy[6] - 9.0*fxy[5] + 45.0*fxy[4] - 45.0*fxy[2] + 9.0*fxy[1] - fxy[0])/(60.0*h);
        dBxdz = (fxz[6] - 9.0*fxz[5] + 45.0*fxz[4] - 45.0*fxz[2] + 9.0*fxz[1] - fxz[0])/(60.0*h);

        dBydx = (fyx[6] - 9.0*fyx[5] + 45.0*fyx[4] - 45.0*fyx[2] + 9.0*fyx[1] - fyx[0])/(60.0*h);
        dBydz = (fyz[6] - 9.0*fyz[5] + 45.0*fyz[4] - 45.0*fyz[2] + 9.0*fyz[1] - fyz[0])/(60.0*h);

        dBzdx = (fzx[6] - 9.0*fzx[5] + 45.0*fzx[4] - 45.0*fzx[2] + 9.0*fzx[1] - fzx[0])/(60.0*h);
        dBzdy = (fzy[6] - 9.0*fzy[5] + 45.0*fzy[4] - 45.0*fzy[2] + 9.0*fzy[1] - fzy[0])/(60.0*h);
    } else if (DerivScheme == LGM_DERIV_FOUR_POINT){
        dBxdy = (-fxy[4] + 8.0*fxy[3] - 8.0*fxy[1] + fxy[0])/(12.0*h);
        dBxdz = (-fxz[4] + 8.0*fxz[3] - 8.0*fxz[1] + fxz[0])/(12.0*h);

        dBydx = (-fyx[4] + 8.0*fyx[3] - 8.0*fyx[1] + fyx[0])/(12.0*h);
        dBydz = (-fyz[4] + 8.0*fyz[3] - 8.0*fyz[1] + fyz[0])/(12.0*h);

        dBzdx = (-fzx[4] + 8.0*fzx[3] - 8.0*fzx[1] + fzx[0])/(12.0*h);
        dBzdy = (-fzy[4] + 8.0*fzy[3] - 8.0*fzy[1] + fzy[0])/(12.0*h);
    } else if (DerivScheme == LGM_DERIV_TWO_POINT){
        dBxdy = (fxy[2] - fxy[0])/(2.0*h);
        dBxdz = (fxz[2] - fxz[0])/(2.0*h);

        dBydx = (fyx[2] - fyx[0])/(2.0*h);
        dBydz = (fyz[2] - fyz[0])/(2.0*h);

        dBzdx = (fzx[2] - fzx[0])/(2.0*h);
        dBzdy = (fzy[2] - fzy[0])/(2.0*h);
    }
    CurlB->x = dBzdy - dBydz;
    CurlB->y = dBxdz - dBzdx;
    CurlB->z = dBydx - dBxdy;
//printf("NRM dBdx = %g %g %g\n", -999.0, dBydx, dBzdx );
//printf("NRM dBdy = %g %g %g\n", dBxdy, -999.0, dBzdy );
//printf("NRM dBdz = %g %g %g\n", dBxdz, dBydz, -999.0 );
    if (m->VerbosityLevel > 0) printf("   CurlB = (%g %g %g)\n", CurlB->x, CurlB->y, CurlB->z );

    return;

}



/**
 *  \brief
 *      Compute the curl of b-hat at a given point.
 *
 *  \details
 *      Computes \f$ \nabla\times b \f$.
 *
 *
 *      \param[in]      u0          Position (in GSM) to use.
 *      \param[out]     CurlB       The computed curl of B at position u0.
 *      \param[in]      DerivScheme Derivative scheme to use (can be one of LGM_DERIV_SIX_POINT, LGM_DERIV_FOUR_POINT, or LGM_DERIV_TWO_POINT).
 *      \param[in]      h           The delta (in Re) to use for grid spacing in the derivative scheme.
 *      \param[in,out]  m           A properly initialized and configured Lgm_MagModelInfo structure.
 *
 *
 *      \author         Mike Henderson
 *      \date           2011
 *
 */
void Lgm_Curlb( Lgm_Vector *u0, Lgm_Vector *Curlb, int DerivScheme, double h, Lgm_MagModelInfo *m ) {

    double      Bmag, H;
    double      fxy[7], fxz[7], fyx[7], fyz[7], fzx[7], fzy[7];
    double      dbxdy, dbxdz, dbydx, dbydz, dbzdx, dbzdy;
    int         i, N;
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

    if (m->VerbosityLevel > 0) printf("\t\tLgm_Curlb: Computing Curlb with DerivScheme = %d,  h = %g\n", DerivScheme, h);
    // dbx/dy and dbx/dz
    for (i=-N; i<=N; ++i){
        if ( N != 0 ) {
            u = *u0; H = (double)i*h; u.y += H;
            m->Bfield( &u, &Bvec, m );
            Bmag = Lgm_Magnitude( &Bvec );
            fxy[i+N] = Bvec.x/Bmag;
        }
    }
    for (i=-N; i<=N; ++i){
        if ( N != 0 ) {
            u = *u0; H = (double)i*h; u.z += H;
            m->Bfield( &u, &Bvec, m );
            Bmag = Lgm_Magnitude( &Bvec );
            fxz[i+N] = Bvec.x/Bmag;
        }
    }

    // dby/dx and dby/dz
    for (i=-N; i<=N; ++i){
        if ( N != 0 ) {
            u = *u0; H = (double)i*h; u.x += H;
            m->Bfield( &u, &Bvec, m );
            Bmag = Lgm_Magnitude( &Bvec );
            fyx[i+N] = Bvec.y/Bmag;
        }
    }
    for (i=-N; i<=N; ++i){
        if ( N != 0 ) {
            u = *u0; H = (double)i*h; u.z += H;
            m->Bfield( &u, &Bvec, m );
            Bmag = Lgm_Magnitude( &Bvec );
            fyz[i+N] = Bvec.y/Bmag;
        }
    }

    // dbz/dx and dbz/dy
    for (i=-N; i<=N; ++i){
        if ( N != 0 ) {
            u = *u0; H = (double)i*h; u.x += H;
            m->Bfield( &u, &Bvec, m );
            Bmag = Lgm_Magnitude( &Bvec );
            fzx[i+N] = Bvec.z/Bmag;
        }
    }
    for (i=-N; i<=N; ++i){
        if ( N != 0 ) {
            u = *u0; H = (double)i*h; u.y += H;
            m->Bfield( &u, &Bvec, m );
            Bmag = Lgm_Magnitude( &Bvec );
            fzy[i+N] = Bvec.z/Bmag;
        }
    }

    if (DerivScheme == LGM_DERIV_SIX_POINT){
        dbxdy = (fxy[6] - 9.0*fxy[5] + 45.0*fxy[4] - 45.0*fxy[2] + 9.0*fxy[1] - fxy[0])/(60.0*h);
        dbxdz = (fxz[6] - 9.0*fxz[5] + 45.0*fxz[4] - 45.0*fxz[2] + 9.0*fxz[1] - fxz[0])/(60.0*h);

        dbydx = (fyx[6] - 9.0*fyx[5] + 45.0*fyx[4] - 45.0*fyx[2] + 9.0*fyx[1] - fyx[0])/(60.0*h);
        dbydz = (fyz[6] - 9.0*fyz[5] + 45.0*fyz[4] - 45.0*fyz[2] + 9.0*fyz[1] - fyz[0])/(60.0*h);

        dbzdx = (fzx[6] - 9.0*fzx[5] + 45.0*fzx[4] - 45.0*fzx[2] + 9.0*fzx[1] - fzx[0])/(60.0*h);
        dbzdy = (fzy[6] - 9.0*fzy[5] + 45.0*fzy[4] - 45.0*fzy[2] + 9.0*fzy[1] - fzy[0])/(60.0*h);
    } else if (DerivScheme == LGM_DERIV_FOUR_POINT){
        dbxdy = (-fxy[4] + 8.0*fxy[3] - 8.0*fxy[1] + fxy[0])/(12.0*h);
        dbxdz = (-fxz[4] + 8.0*fxz[3] - 8.0*fxz[1] + fxz[0])/(12.0*h);

        dbydx = (-fyx[4] + 8.0*fyx[3] - 8.0*fyx[1] + fyx[0])/(12.0*h);
        dbydz = (-fyz[4] + 8.0*fyz[3] - 8.0*fyz[1] + fyz[0])/(12.0*h);

        dbzdx = (-fzx[4] + 8.0*fzx[3] - 8.0*fzx[1] + fzx[0])/(12.0*h);
        dbzdy = (-fzy[4] + 8.0*fzy[3] - 8.0*fzy[1] + fzy[0])/(12.0*h);
    } else if (DerivScheme == LGM_DERIV_TWO_POINT){
        dbxdy = (fxy[2] - fxy[0])/(2.0*h);
        dbxdz = (fxz[2] - fxz[0])/(2.0*h);

        dbydx = (fyx[2] - fyx[0])/(2.0*h);
        dbydz = (fyz[2] - fyz[0])/(2.0*h);

        dbzdx = (fzx[2] - fzx[0])/(2.0*h);
        dbzdy = (fzy[2] - fzy[0])/(2.0*h);
    }
    Curlb->x = dbzdy - dbydz;
    Curlb->y = dbxdz - dbzdx;
    Curlb->z = dbydx - dbxdy;
//printf("NRM dbdx = %g %g %g\n", -999.0, dbydx, dbzdx );
//printf("NRM dbdy = %g %g %g\n", dbxdy, -999.0, dbzdy );
//printf("NRM dbdz = %g %g %g\n", dbxdz, dbydz, -999.0 );
    if (m->VerbosityLevel > 0) printf("   Curlb = (%g %g %g)\n", Curlb->x, Curlb->y, Curlb->z );

    return;

}


/**
 *  \brief
 *      Compute the curl of B and it's parallel and perpendicular components at a given point.
 *
 *  \details
 *      Computes \f$ \nabla\times B \f$, \f$ (\nabla\times B)_\parallel \f$, and \f$ (\nabla\times B)_\perp \f$.
 *
 *
 *      \param[in]      u0          Position (in GSM) to use.
 *      \param[out]     CurlB       The computed curl of B at position u0.
 *      \param[out]     CurlB_para  The computed parallel component of curl of B at position u0.
 *      \param[out]     CurlB_perp  The computed perpendicular component of curl of B at position u0.
 *      \param[in]      DerivScheme Derivative scheme to use (can be one of LGM_DERIV_SIX_POINT, LGM_DERIV_FOUR_POINT, or LGM_DERIV_TWO_POINT).
 *      \param[in]      h           The delta (in Re) to use for grid spacing in the derivative scheme.
 *      \param[in,out]  m           A properly initialized and configured Lgm_MagModelInfo structure.
 *
 *
 *      \author         Mike Henderson
 *      \date           2011
 *
 */
void Lgm_CurlB2( Lgm_Vector *u0, Lgm_Vector *CurlB, Lgm_Vector *CurlB_para, Lgm_Vector *CurlB_perp, int DerivScheme, double h, Lgm_MagModelInfo *m ) {

    double      g;
    Lgm_Vector  Bvec;


    Lgm_CurlB( u0, CurlB, DerivScheme, h, m );


    m->Bfield( u0, &Bvec, m );
    Lgm_NormalizeVector( &Bvec );

    // Compute parallel component of CurlB
    g = Lgm_DotProduct( &Bvec, CurlB );
    *CurlB_para = *CurlB;
    Lgm_ScaleVector( CurlB_para, g );


    // Compute perpendicular component of CurlB (CurlB_perp = CurlB - CurlB_para)
    Lgm_VecSub( CurlB_perp, CurlB, CurlB_para );


    return;

}





/**
 *  \brief
 *      Compute the divergence of B at a given point.
 *
 *  \details
 *      Computes \f$ \nabla\cdot B \f$. By Maxwell's equations it should be
 *      zero. This routine is more for verification/validation purposes. Note
 *      that the dipole fields all return veru small numbers (essentially
 *      zero). T87 and T89 do also, T04s and OP77 dont seem to though...
 *
 *
 *      \param[in]      u0          Position (in GSM) to use.
 *      \param[out]     DivB        The computed divergence of B at position u0.
 *      \param[in]      DerivScheme Derivative scheme to use (can be one of LGM_DERIV_SIX_POINT, LGM_DERIV_FOUR_POINT, or LGM_DERIV_TWO_POINT).
 *      \param[in]      h           The delta (in Re) to use for grid spacing in the derivative scheme.
 *      \param[in,out]  m           A properly initialized and configured Lgm_MagModelInfo structure.
 *
 *
 *      \author         Mike Henderson
 *      \date           2011
 *
 */
void Lgm_DivB( Lgm_Vector *u0, double *DivB, int DerivScheme, double h, Lgm_MagModelInfo *m ) {

    double      B, H;
    double      fxx[7], fyy[7], fzz[7];
    double      dBxdx, dBydy, dBzdz;
    int         i, N;
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

    if (m->VerbosityLevel > 0) printf("\t\tLgm_DivB: Computing DivB with DerivScheme = %d,  h = %g", DerivScheme, h);
    // dBx/dx
    for (i=-N; i<=N; ++i){
        if ( N != 0 ) {
            u = *u0; H = (double)i*h; u.x += H;
            m->Bfield( &u, &Bvec, m );
            fxx[i+N] = Bvec.x;
        }
    }

    // dBy/dy
    for (i=-N; i<=N; ++i){
        if ( N != 0 ) {
            u = *u0; H = (double)i*h; u.y += H;
            m->Bfield( &u, &Bvec, m );
            fyy[i+N] = Bvec.y;
        }
    }

    // dBz/dz
    for (i=-N; i<=N; ++i){
        if ( N != 0 ) {
            u = *u0; H = (double)i*h; u.z += H;
            m->Bfield( &u, &Bvec, m );
            fzz[i+N] = Bvec.z;
        }
    }

    if (DerivScheme == LGM_DERIV_SIX_POINT){
        dBxdx = (fxx[6] - 9.0*fxx[5] + 45.0*fxx[4] - 45.0*fxx[2] + 9.0*fxx[1] - fxx[0])/(60.0*h);
        dBydy = (fyy[6] - 9.0*fyy[5] + 45.0*fyy[4] - 45.0*fyy[2] + 9.0*fyy[1] - fyy[0])/(60.0*h);
        dBzdz = (fzz[6] - 9.0*fzz[5] + 45.0*fzz[4] - 45.0*fzz[2] + 9.0*fzz[1] - fzz[0])/(60.0*h);
    } else if (DerivScheme == LGM_DERIV_FOUR_POINT){
        dBxdx = (-fxx[4] + 8.0*fxx[3] - 8.0*fxx[1] + fxx[0])/(12.0*h);
        dBydy = (-fyy[4] + 8.0*fyy[3] - 8.0*fyy[1] + fyy[0])/(12.0*h);
        dBzdz = (-fzz[4] + 8.0*fzz[3] - 8.0*fzz[1] + fzz[0])/(12.0*h);
    } else if (DerivScheme == LGM_DERIV_TWO_POINT){
        dBxdx = (fxx[2] - fxx[0])/(2.0*h);
        dBydy = (fyy[2] - fyy[0])/(2.0*h);
        dBzdz = (fzz[2] - fzz[0])/(2.0*h);
    }
    *DivB = dBxdx + dBydy + dBzdz;
    if (m->VerbosityLevel > 0) printf("   dBxdx, dBydy, dBzdz, DivB = %g %g %g %g\n", dBxdx, dBydy, dBzdz, *DivB );

    return;

}

/**
 *  \brief
 *      Compute the divergence of b-hat at a given point.
 *
 *  \details
 *      Computes \f$ \nabla\cdot \hat{b} \f$. Note that while by Maxwell's
 *      equations, \f$ \nabla\cdot \vec{B} \f$ should be zero, the divergence of
 *      the b-field unit vector is not. This quantity is needed by some
 *      calculations (e.g. calculating higher order terms in the first
 *      adiabatic invariant ala Burby et al. 2013.)
 *
 *
 *      \param[in]      u0          Position (in GSM) to use.
 *      \param[out]     Divb        The computed divergence of b-hat at position u0.
 *      \param[in]      DerivScheme Derivative scheme to use (can be one of LGM_DERIV_SIX_POINT, LGM_DERIV_FOUR_POINT, or LGM_DERIV_TWO_POINT).
 *      \param[in]      h           The delta (in Re) to use for grid spacing in the derivative scheme.
 *      \param[in,out]  m           A properly initialized and configured Lgm_MagModelInfo structure.
 *
 *
 *      \author         Mike Henderson
 *      \date           2023
 *
 */
void Lgm_Divb( Lgm_Vector *u0, double *Divb, int DerivScheme, double h, Lgm_MagModelInfo *m ) {

    double      B, H;
    double      fxx[7], fyy[7], fzz[7];
    double      dbxdx, dbydy, dbzdz;
    int         i, N;
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

    if (m->VerbosityLevel > 0) printf("\t\tLgm_DivB: Computing DivB with DerivScheme = %d,  h = %g", DerivScheme, h);
    // dBx/dx
    for (i=-N; i<=N; ++i){
        if ( N != 0 ) {
            u = *u0; H = (double)i*h; u.x += H;
            m->Bfield( &u, &Bvec, m );
            Lgm_NormalizeVector( &Bvec );
            fxx[i+N] = Bvec.x;
        }
    }

    // dBy/dy
    for (i=-N; i<=N; ++i){
        if ( N != 0 ) {
            u = *u0; H = (double)i*h; u.y += H;
            m->Bfield( &u, &Bvec, m );
            Lgm_NormalizeVector( &Bvec );
            fyy[i+N] = Bvec.y;
        }
    }

    // dBz/dz
    for (i=-N; i<=N; ++i){
        if ( N != 0 ) {
            u = *u0; H = (double)i*h; u.z += H;
            m->Bfield( &u, &Bvec, m );
            Lgm_NormalizeVector( &Bvec );
            fzz[i+N] = Bvec.z;
        }
    }

    if (DerivScheme == LGM_DERIV_SIX_POINT){
        dbxdx = (fxx[6] - 9.0*fxx[5] + 45.0*fxx[4] - 45.0*fxx[2] + 9.0*fxx[1] - fxx[0])/(60.0*h);
        dbydy = (fyy[6] - 9.0*fyy[5] + 45.0*fyy[4] - 45.0*fyy[2] + 9.0*fyy[1] - fyy[0])/(60.0*h);
        dbzdz = (fzz[6] - 9.0*fzz[5] + 45.0*fzz[4] - 45.0*fzz[2] + 9.0*fzz[1] - fzz[0])/(60.0*h);
    } else if (DerivScheme == LGM_DERIV_FOUR_POINT){
        dbxdx = (-fxx[4] + 8.0*fxx[3] - 8.0*fxx[1] + fxx[0])/(12.0*h);
        dbydy = (-fyy[4] + 8.0*fyy[3] - 8.0*fyy[1] + fyy[0])/(12.0*h);
        dbzdz = (-fzz[4] + 8.0*fzz[3] - 8.0*fzz[1] + fzz[0])/(12.0*h);
    } else if (DerivScheme == LGM_DERIV_TWO_POINT){
        dbxdx = (fxx[2] - fxx[0])/(2.0*h);
        dbydy = (fyy[2] - fyy[0])/(2.0*h);
        dbzdz = (fzz[2] - fzz[0])/(2.0*h);
    }
    *Divb = dbxdx + dbydy + dbzdz;
    if (m->VerbosityLevel > 0) printf("   dbxdx, dbydy, dbzdz, Divb = %g %g %g %g\n", dbxdx, dbydy, dbzdz, *Divb );

    return;

}



/**
 *  \brief
 *      Compute the instantaneous gradient and curvature drift velocity of the guiding center.
 *
 *
 *
 *  \details
 *      The gradient and curvature drift velocity of the guiding center is given by;
 *          \f[
 *              V_s = {mv^2\over qB^2}\left(\left(1-{B\over
 *              2B_m}\right){\vec{B}\times\nabla B\over B} + \left(1-{B\over
 *              B_m}\right){(\nabla\times \vec{B})_perp} \right)
 *          \f]
 *      (e.g. see Sukhtina, M.A., "On the calculation of the magnetic drift
 *      velocity of particles with arbitrary pitch angles, Planet Space Sci.
 *      41, 327--331, 1993.) With \f$\eta = T/E_0\f$, where \f$T\f$ is the
 *      particle's kinetic energy, and \f$E_0\f$ its rest energy, this can be
 *      rewritten as,
 *          \f[
 *              V_s = {\eta (2+\eta)\over (1+\eta)}{E_0\over qB^2}\left(\left(1-{B\over
 *              2B_m}\right){\vec{B}\times\nabla B\over B} + \left(1-{B\over
 *              B_m}\right){(\nabla\times \vec{B})_perp} \right)
 *          \f]
 *
 *
 *      \param[in]      u0          Position (in GSM) to use. [Re]
 *      \param[out]     Vel         The computed drift velocity in GSM coords. [km/s]
 *      \param[in]      q           Charge of the particle. [C]
 *      \param[in]      T           Kinetic energy of the particle. [MeV]
 *      \param[in]      E0          Rest energy of the particle. [MeV]
 *      \param[in]      Bm          Magnitude of B-field at mirror point for particle. [nT]
 *      \param[in]      DerivScheme Derivative scheme to use (can be one of LGM_DERIV_SIX_POINT, LGM_DERIV_FOUR_POINT, or LGM_DERIV_TWO_POINT).
 *      \param[in]      h           The delta (in Re) to use for grid spacing in the derivative scheme.
 *      \param[in,out]  m           A properly initialized and configured Lgm_MagModelInfo structure.
 *
 *
 *      \author         Mike Henderson
 *      \date           2011
 *
 */
int Lgm_GradAndCurvDriftVel( Lgm_Vector *u0, Lgm_Vector *Vel, Lgm_MagModelInfo *m ) {

    double      B, BoverBm, g, eta, h, Beta, Beta2, Gamma;
    double      q, T, E0, Bm;
    int         DerivScheme;
    Lgm_Vector  CurlB, CurlB_para, CurlB_perp, GradB, Bvec, Q, R, S, W, Z;
    Lgm_Vector  B_Cross_GradB;

    q           = m->Lgm_VelStep_q;
    T           = m->Lgm_VelStep_T;
    E0          = m->Lgm_VelStep_E0;
    Bm          = m->Lgm_VelStep_Bm;
    h           = m->Lgm_VelStep_h;
    DerivScheme = m->Lgm_VelStep_DerivScheme;

//printf("Lgm_GradAndCurvDriftVel: u0 = %g %g %g\n", u0->x, u0->y, u0->z);
    m->Bfield( u0, &Bvec, m );
    B = Lgm_Magnitude( &Bvec );
    BoverBm = B/Bm;
//printf("BoverBm = %g\n", BoverBm);


    /*
     * Compute B cross GradB  [nT/Re]
     */
    Lgm_GradB( u0, &GradB, DerivScheme, h, m );
    Lgm_CrossProduct( &Bvec, &GradB, &B_Cross_GradB );


    /*
     * Compute (curl B)_perp [nT/Re]
     */
    Lgm_CurlB2( u0, &CurlB, &CurlB_para, &CurlB_perp, DerivScheme, h, m );


    /*
     * Compute instantaneous gradient/curv drift Velocity.
     */
    Q = B_Cross_GradB;
    g = (1.0-0.5*BoverBm)/B;
//printf("g1 = %g\n", g);
    Lgm_ScaleVector( &Q, g ); // [nT/Re]

    R = CurlB_perp;
    g = 1.0-BoverBm;
//printf("g2 = %g\n", g);
    Lgm_ScaleVector( &R, g ); // [nT/Re]

    Lgm_VecAdd( &S, &Q, &R ); // [nT/Re]
    eta = T/E0; // kinetic energy over rest energy
    Gamma = eta + 1.0;
    Beta2 = 1.0 - 1.0/(Gamma*Gamma);
    Beta  = sqrt( Beta2 );
    g = Gamma*Beta2*E0*1e6*Joules_Per_eV/(q*B*B); // [ N m / (A s nT^2) ]  
    // gS would have units of [ N m / (A s nT Re) ], so convert nT and Re to SI
    // to get a final answer in m/s. Then convert to km/s.
    g *= 1e9/(Re*1e6);        
    g /= Re;        
//Re/s
//printf("g3 = %g\n", g);
    Lgm_ScaleVector( &S, g ); // [ N m / (A s T m) ] = [ N m / (A s (N/(A m) m) ] = [ m/s ]
// now its Re/s


    // Add in parallel velocity (Note that v = c*Beta)
    if (BoverBm >= 1.0) {
        g = 0.0;
    } else {
        g = -LGM_c/(1000.0*Re)*Beta*sqrt(1.0-BoverBm)/B;
    }
//printf("T = %g, Beta = %g g = %g\n", T, Beta, g);
    W = Bvec;
    Lgm_ScaleVector( &W, g );
//printf("W = %g %g %g\n", W.x, W.y, W.z);

    Lgm_VecAdd( &Z, &S, &W );
    

    *Vel = Z;


    return(1);

}


/**
 *  \brief
 *      Compute the field radius of curvature vector.
 *
 *  \details
 *      The gradient and curvature drift velocity.
 *
 *
 *      \param[in]      u0          Position (in GSM) to use. [Re]
 *      \param[out]     Rc          The radius of curvature vector. [Re]
 *      \param[in]      DerivScheme Derivative scheme to use (can be one of LGM_DERIV_SIX_POINT, LGM_DERIV_FOUR_POINT, or LGM_DERIV_TWO_POINT).
 *      \param[in,out]  m           A properly initialized and configured Lgm_MagModelInfo structure.
 *
 *
 *      \author         Mike Henderson
 *      \date           2023
 *
 */
void Lgm_Curvature( Lgm_Vector *q, Lgm_Vector *dbds, int DerivScheme, double h, Lgm_MagModelInfo *m ){

    double      H, Htry, Hdid, Hnext, sgn, s, fs[7], gs[7], hs[7];
    int         i, N, reset;
    Lgm_Vector  u, u_scale, b;


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
        default:
            N = 3;
    }

    u_scale.x = u_scale.y = u_scale.z = 1.0;
    if (m->VerbosityLevel > 0) printf("\t\tLgm_Curvature: Computing Curvature of B with DerivScheme = %d,  h = %g", DerivScheme, h);
    for (i=-N; i<=N; ++i){
        if ( i != 0 ) {
            u = *q; H = (double)i*h;
            reset = TRUE;
            if ( H < 0.0 ) { Htry = -H; sgn = -1.0; } else { Htry = H; sgn = 1.0; }
            Lgm_MagStep( &u, &u_scale, Htry, &Hdid, &Hnext, sgn, &s, &reset, m->Bfield, m );
            m->Bfield( &u, &b, m );
            Lgm_NormalizeVector( &b );
            fs[i+N] = b.x; //Re^-1
            gs[i+N] = b.y; //Re^-1
            hs[i+N] = b.z; //Re^-1
        }
    }


    if (DerivScheme == LGM_DERIV_SIX_POINT){
        dbds->x = (fs[6] - 9.0*fs[5] + 45.0*fs[4] - 45.0*fs[2] + 9.0*fs[1] - fs[0])/(60.0*h); // Re^-1
        dbds->y = (gs[6] - 9.0*gs[5] + 45.0*gs[4] - 45.0*gs[2] + 9.0*gs[1] - gs[0])/(60.0*h); // Re^-1
        dbds->z = (hs[6] - 9.0*hs[5] + 45.0*hs[4] - 45.0*hs[2] + 9.0*hs[1] - hs[0])/(60.0*h); // Re^-1
    } else if (DerivScheme == LGM_DERIV_FOUR_POINT){
        dbds->x = (-fs[4] + 8.0*fs[3] - 8.0*fs[1] + fs[0])/(12.0*h); // Re^-1
        dbds->y = (-gs[4] + 8.0*gs[3] - 8.0*gs[1] + gs[0])/(12.0*h); // Re^-1
        dbds->z = (-hs[4] + 8.0*hs[3] - 8.0*hs[1] + hs[0])/(12.0*h); // Re^-1
    } else if (DerivScheme == LGM_DERIV_TWO_POINT){
        dbds->x = (fs[2] - fs[0])/(2.0*h); // Re^-1
        dbds->y = (gs[2] - gs[0])/(2.0*h); // Re^-1
        dbds->z = (hs[2] - hs[0])/(2.0*h); // Re^-1
    }
//    if (m->VerbosityLevel > 0) printf("   GradB = (%g %g %g)\n", GradB->x, GradB->y, GradB->z );
    

    return;

}

void Lgm_Curvature2( Lgm_Vector *q, Lgm_Vector *dbds, int DerivScheme, double h, Lgm_MagModelInfo *m ){

    double      H, Htry, Hdid, Hnext, sgn, s, fs[7], gs[7], hs[7];
    double      fxx[7], fxy[7], fxz[7];
    double      fyx[7], fyy[7], fyz[7];
    double      fzx[7], fzy[7], fzz[7];
    double      dbx_dx, dby_dx, dbz_dx, dbx_dy, dby_dy, dbz_dy, dbx_dz, dby_dz, dbz_dz;
    int         i, N, reset;
    Lgm_Vector  u0, u, u_scale, b;


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
        default:
            N = 3;
    }

    if (m->VerbosityLevel > 0) printf("\t\tLgm_Curvature: Computing Curvature of B with DerivScheme = %d,  h = %g", DerivScheme, h);

    // x derivs
    for (i=-N; i<=N; ++i){
        if ( i != 0 ) {
            u = *q; H = (double)i*h; u.x += H;
            m->Bfield( &u, &b, m );
            Lgm_NormalizeVector( &b );
            fxx[i+N] = b.x; //Re^-1
            fyx[i+N] = b.y; //Re^-1
            fzx[i+N] = b.z; //Re^-1
        }
    }
    // y derivs
    for (i=-N; i<=N; ++i){
        if ( i != 0 ) {
            u = *q; H = (double)i*h; u.y += H;
            m->Bfield( &u, &b, m );
            Lgm_NormalizeVector( &b );
            fxy[i+N] = b.x; //Re^-1
            fyy[i+N] = b.y; //Re^-1
            fzy[i+N] = b.z; //Re^-1
        }
    }
    // z derivs
    for (i=-N; i<=N; ++i){
        if ( i != 0 ) {
            u = *q; H = (double)i*h; u.z += H;
            m->Bfield( &u, &b, m );
            Lgm_NormalizeVector( &b );
            fxz[i+N] = b.x; //Re^-1
            fyz[i+N] = b.y; //Re^-1
            fzz[i+N] = b.z; //Re^-1
        }
    }


    if (DerivScheme == LGM_DERIV_SIX_POINT){
        dbx_dx = (fxx[6] - 9.0*fxx[5] + 45.0*fxx[4] - 45.0*fxx[2] + 9.0*fxx[1] - fxx[0])/(60.0*h); // Re^-1
        dby_dx = (fyx[6] - 9.0*fyx[5] + 45.0*fyx[4] - 45.0*fyx[2] + 9.0*fyx[1] - fyx[0])/(60.0*h); // Re^-1
        dbz_dx = (fzx[6] - 9.0*fzx[5] + 45.0*fzx[4] - 45.0*fzx[2] + 9.0*fzx[1] - fzx[0])/(60.0*h); // Re^-1
        dbx_dy = (fxy[6] - 9.0*fxy[5] + 45.0*fxy[4] - 45.0*fxy[2] + 9.0*fxy[1] - fxy[0])/(60.0*h); // Re^-1
        dby_dy = (fyy[6] - 9.0*fyy[5] + 45.0*fyy[4] - 45.0*fyy[2] + 9.0*fyy[1] - fyy[0])/(60.0*h); // Re^-1
        dbz_dy = (fzy[6] - 9.0*fzy[5] + 45.0*fzy[4] - 45.0*fzy[2] + 9.0*fzy[1] - fzy[0])/(60.0*h); // Re^-1
        dbx_dz = (fxz[6] - 9.0*fxz[5] + 45.0*fxz[4] - 45.0*fxz[2] + 9.0*fxz[1] - fxz[0])/(60.0*h); // Re^-1
        dby_dz = (fyz[6] - 9.0*fyz[5] + 45.0*fyz[4] - 45.0*fyz[2] + 9.0*fyz[1] - fyz[0])/(60.0*h); // Re^-1
        dbz_dz = (fzz[6] - 9.0*fzz[5] + 45.0*fzz[4] - 45.0*fzz[2] + 9.0*fzz[1] - fzz[0])/(60.0*h); // Re^-1
    } else if (DerivScheme == LGM_DERIV_FOUR_POINT){
        dbx_dx = (-fxx[4] + 8.0*fxx[3] - 8.0*fxx[1] + fxx[0])/(12.0*h); // Re^-1
        dby_dx = (-fyx[4] + 8.0*fyx[3] - 8.0*fyx[1] + fyx[0])/(12.0*h); // Re^-1
        dbz_dx = (-fzx[4] + 8.0*fzx[3] - 8.0*fzx[1] + fzx[0])/(12.0*h); // Re^-1
        dbx_dy = (-fxy[4] + 8.0*fxy[3] - 8.0*fxy[1] + fxy[0])/(12.0*h); // Re^-1
        dby_dy = (-fyy[4] + 8.0*fyy[3] - 8.0*fyy[1] + fyy[0])/(12.0*h); // Re^-1
        dbz_dy = (-fzy[4] + 8.0*fzy[3] - 8.0*fzy[1] + fzy[0])/(12.0*h); // Re^-1
        dbx_dz = (-fxz[4] + 8.0*fxz[3] - 8.0*fxz[1] + fxz[0])/(12.0*h); // Re^-1
        dby_dz = (-fyz[4] + 8.0*fyz[3] - 8.0*fyz[1] + fyz[0])/(12.0*h); // Re^-1
        dbz_dz = (-fzz[4] + 8.0*fzz[3] - 8.0*fzz[1] + fzz[0])/(12.0*h); // Re^-1
    } else if (DerivScheme == LGM_DERIV_TWO_POINT){
        dbx_dx = (fxx[2] - fxx[0])/(2.0*h); // Re^-1
        dby_dx = (fyx[2] - fyx[0])/(2.0*h); // Re^-1
        dbz_dx = (fzx[2] - fzx[0])/(2.0*h); // Re^-1
        dbx_dy = (fxy[2] - fxy[0])/(2.0*h); // Re^-1
        dby_dy = (fyy[2] - fyy[0])/(2.0*h); // Re^-1
        dbz_dy = (fzy[2] - fzy[0])/(2.0*h); // Re^-1
        dbx_dz = (fxz[2] - fxz[0])/(2.0*h); // Re^-1
        dby_dz = (fyz[2] - fyz[0])/(2.0*h); // Re^-1
        dbz_dz = (fzz[2] - fzz[0])/(2.0*h); // Re^-1
    }
//    if (m->VerbosityLevel > 0) printf("   GradB = (%g %g %g)\n", GradB->x, GradB->y, GradB->z );
    
    m->Bfield( q, &b, m );
    Lgm_NormalizeVector( &b );
    dbds->x = b.x * dbx_dx + b.y * dbx_dy + b.z * dbx_dz;
    dbds->y = b.x * dby_dx + b.y * dby_dy + b.z * dby_dz;
    dbds->z = b.x * dbz_dx + b.y * dbz_dy + b.z * dbz_dz;

    return;

}






/**
 *  \brief
 *      Compute the instantaneous gradient and curvature drift velocity of the guiding center.
 *
 *  \details
 *      The gradient and curvature drift velocity of the guiding center is given by;
 *
 *
 *      \param[in]      u0          Position (in GSM) to use. [Re]
 *      \param[in]      v           Total velocity (in GSM) to use. [m/s]
 *      \param[out]     v_GradB     The computed GradB drift velocity in GSM coords. [m/s]
 *      \param[out]     v_CurvB     The computed Curvature drift velocity in GSM coords. [m/s]
 *      \param[out]     v_CurvB     The computed Curvature drift velocity in GSM coords. [m/s]
 *      \param[in,out]  m           A properly initialized and configured Lgm_MagModelInfo structure.
 *
 *
 *      \author         Mike Henderson
 *      \date           2023
 *
 */
int Lgm_Vdrift_GradB_Curv( Lgm_Vector *u0, Lgm_Vector *v_full, double Mass, double Charge, Lgm_Vector *v_GradB, Lgm_Vector *v_CurvB, Lgm_Vector *Rc_vec, Lgm_MagModelInfo *m ) {

    double      Gamma;
    double      q, T, E0, Bm;
    Lgm_Vector  Kappa, B_vec, GradB;
    double      K, Rc, Rc2, B, B2, B3, u, u2, v, v2, w2, g;


    m->Bfield( u0, &B_vec, m );
    B_vec.x *= 1e-9; B_vec.y *= 1e-9; B_vec.z *= 1e-9;  // T
    B = Lgm_Magnitude( &B_vec );
    B2 = B*B; B3 = B2*B;

    v = Lgm_Magnitude( v_full ); // speed in m/s
    v2 = v*v;
    u = Lgm_DotProduct( &B_vec, v_full)/B; // u_parallel in m/s
    u2 = u*u;
    w2 = v2 - u2; // u_perp^2 in m^2/s^2
    Gamma = 1.0/sqrt( 1.0 - v2/(LGM_c*LGM_c) );

Lgm_Vector dbds;
    Lgm_Curvature( u0, &dbds, LGM_DERIV_SIX_POINT, 1e-4, m);
    g = 1.0/(Re*1e3); dbds.x *= g; dbds.y *= g; dbds.z *= g;
//    K = Lgm_Magnitude( &Kappa );
 //   Rc = 1.0/K;
//    Rc2 = Rc*Rc;
//printf("Lgm_Vdrift_GradB_Curv: Rc = %g\n", Rc);
//    *Rc_vec = Kappa;
    
//    Lgm_NormalizeVector( Rc_vec );
//    Lgm_ScaleVector( Rc_vec, -Rc );

//Lgm_Vector Rc_vec_over_Rc2;
//Rc_vec_over_Rc2.x = -dbds.x;
//Rc_vec_over_Rc2.y = -dbds.y;
//Rc_vec_over_Rc2.z = -dbds.z;
//    Lgm_CrossProduct( &Rc_vec_over_Rc2, &B_vec, v_CurvB );
Lgm_CrossProduct( &B_vec, &dbds, v_CurvB );

    Lgm_ScaleVector( v_CurvB, Gamma*Mass*u2/(Charge*B2) );





    Lgm_GradB( u0, &GradB, LGM_DERIV_SIX_POINT, 1e-4, m );
    g = 1e-9/(Re*1e3); GradB.x *= g; GradB.y *= g; GradB.z *= g;  // T/m
    Lgm_CrossProduct( &B_vec, &GradB, v_GradB );
    Lgm_ScaleVector( v_GradB, Gamma*Mass*w2/(2.0*Charge*B3) );

    Lgm_ScaleVector( Rc_vec, Rc );

//printf("u0 = %g %g %g\n", u0->x, u0->y, u0->z );
//printf("v_GradB = %g %g %g\n", v_GradB->x, v_GradB->y, v_GradB->z );
//printf("v_CurvB = %g %g %g\n", v_CurvB->x, v_CurvB->y, v_CurvB->z );

    return(1);

}
