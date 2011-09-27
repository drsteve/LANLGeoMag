#include <stdio.h>
#include <math.h>
#include "Lgm/Lgm_MagModelInfo.h"


/**
 *  \brief
 *      Compute the gradient of B at a given point.
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
        u = *u0; H = (double)i*h; u.x += H;
        m->Bfield( &u, &Bvec, m );
        B = Lgm_Magnitude( &Bvec );
        fx[i+N] = B;
    }
    for (i=-N; i<=N; ++i){
        u = *u0; H = (double)i*h; u.y += H;
        m->Bfield( &u, &Bvec, m );
        B = Lgm_Magnitude( &Bvec );
        fy[i+N] = B;
    }
    for (i=-N; i<=N; ++i){
        u = *u0; H = (double)i*h; u.z += H;
        m->Bfield( &u, &Bvec, m );
        B = Lgm_Magnitude( &Bvec );
        fz[i+N] = B;
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
 *      Compute the curl of B at a given point.
 *
 *  \details
 *      Computes \f$ \nabla\times B \f$.
 *
 *
 *      \param[in]      u0          Position (in GSM) to use.
 *      \param[out]     GradB       The computed curl of B at position u0.
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
        u = *u0; H = (double)i*h; u.y += H;
        m->Bfield( &u, &Bvec, m );
        fxy[i+N] = Bvec.x;
    }
    for (i=-N; i<=N; ++i){
        u = *u0; H = (double)i*h; u.z += H;
        m->Bfield( &u, &Bvec, m );
        fxz[i+N] = Bvec.x;
    }

    // dBy/dx and dBy/dz
    for (i=-N; i<=N; ++i){
        u = *u0; H = (double)i*h; u.x += H;
        m->Bfield( &u, &Bvec, m );
        fyx[i+N] = Bvec.y;
    }
    for (i=-N; i<=N; ++i){
        u = *u0; H = (double)i*h; u.z += H;
        m->Bfield( &u, &Bvec, m );
        fyz[i+N] = Bvec.y;
    }

    // dBz/dx and dBz/dy
    for (i=-N; i<=N; ++i){
        u = *u0; H = (double)i*h; u.x += H;
        m->Bfield( &u, &Bvec, m );
        fzx[i+N] = Bvec.z;
    }
    for (i=-N; i<=N; ++i){
        u = *u0; H = (double)i*h; u.y += H;
        m->Bfield( &u, &Bvec, m );
        fzy[i+N] = Bvec.z;
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
    if (m->VerbosityLevel > 0) printf("   CurlB = (%g %g %g)\n", CurlB->x, CurlB->y, CurlB->z );

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
        u = *u0; H = (double)i*h; u.x += H;
        m->Bfield( &u, &Bvec, m );
        fxx[i+N] = Bvec.x;
    }

    // dBy/dy
    for (i=-N; i<=N; ++i){
        u = *u0; H = (double)i*h; u.y += H;
        m->Bfield( &u, &Bvec, m );
        fyy[i+N] = Bvec.y;
    }

    // dBz/dz
    for (i=-N; i<=N; ++i){
        u = *u0; H = (double)i*h; u.z += H;
        m->Bfield( &u, &Bvec, m );
        fzz[i+N] = Bvec.z;
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
 *      Compute the instantaneous gradient and curvature drift velocity of the guiding center.
 *
 *
 *
 *  \details
 *      The gradient and curvature drift velocity of the guiding center is given by;
 *          \f[
 *              V_s = {mv^2\overqB^2}\left{\left(1-{B\over
 *              2B_m}\right){\vec{B}\times\nabla B\over B} + \left(1-{B\over
 *              B_m}\right){(\nabla\times \vec{B})_perp} \right}
 *          \f]
 *      (e.g. see Sukhtina, M.A., "One the calculation of the magnetic srift
 *      velocity of particles with arbitrary pitch angles, Planet Space Sci.
 *      41, 327--331, 1993.) With \f$\eta = T/E_0\f$, where \f$T\f$ is the
 *      particle's kinetic energy, and \f$E_0\f$ its rest energy, this can be
 *      rewritten as,
 *          \f[
 *              V_s = {\eta (2+\eta)\over (1+\eta)}E_0\overqB^2}\left{\left(1-{B\over
 *              2B_m}\right){\vec{B}\times\nabla B\over B} + \left(1-{B\over
 *              B_m}\right){(\nabla\times \vec{B})_perp} \right}
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
void Lgm_GradAndCurvDriftVel( Lgm_Vector *u0, Lgm_Vector *Vel, double q, double T, double E0, double Bm, int DerivScheme, double h, Lgm_MagModelInfo *m ) {

    double      B, BoverBm, g, eta;
    Lgm_Vector  CurlB, CurlB_para, CurlB_perp, GradB, Bvec, Q, R, S;
    Lgm_Vector  B_Cross_GradB;

    m->Bfield( u0, &Bvec, m );
    B = Lgm_Magnitude( &Bvec );
    BoverBm = B/Bm;


    /*
     * Compute (B cross GradB)/B 
     */
    Lgm_GradB( u0, &GradB, DerivScheme, h, m );
    Lgm_CrossProduct( &Bvec, &GradB, &B_Cross_GradB );


    /*
     * Compute (curl B)_perp
     */
    Lgm_CurlB2( u0, &CurlB, &CurlB_para, &CurlB_perp, DerivScheme, h, m );


    /*
     * Compute instantaneous gradient/curv drift Velocity.
     */
    Q = B_Cross_GradB;
    g = (1.0-0.5*BoverBm)/B;
    Lgm_ScaleVector( &Q, g );

    R = CurlB_perp;
    g = 1.0-BoverBm;
    Lgm_ScaleVector( &R, g );

    Lgm_VecAdd( &S, &Q, &R );
    eta = T/E0; // kinetic energy over rest energy
    g = eta*(2.0+eta)/(1.0+eta)*E0/(q*B*B)*LGM_e*1e9/Re;
    Lgm_ScaleVector( &S, g );

}
