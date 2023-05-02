#include "Lgm/Lgm_MagModelInfo.h"                                                                                                                                                   

/*
 * Get the a, c, and b unit vectors, given q, phi and mInfo. q is the position
 * of the particle, phi is the gyrophase (defined in terms of the two
 * coordinate systems constructed here.)
 *
 * Set up the coordinate system that described the particle position
 *
 * We need to define a coordinate system that is perpendicular to the
 * B-field. See Figure 1 of Weyssow and Balescu, (1986), Hamiltonian theory
 * of guiding center motion revisited, J. Plasma Physics, 35, 449-471.
 *
 *      q: position of the particle
 *      b: magnetic field direction at q
 *
 *       r: outward radial unit vector
 *      e1: westward unit vector that is perp to both r and b. I.e. r x b = e1
 *      e2: completes RHS system e1 and b. I.e. b x e1 = e2
 *
 *
 *  NOTE: In general we would need to take care of the case where the radius is
 *  +/- b-hat...
 *
 *
 */
void Lgm_Set_Particle_Frames( Lgm_Vector *q, double Mass, double RestEnergy, double Energy, double Phi, double Delta, Lgm_Vector *a, Lgm_Vector *c, Lgm_Vector *b, Lgm_Vector *v, Lgm_MagModelInfo *mInfo ) {

    Lgm_Vector  r, e1, e2; 
    //double      m;
    double      SinPhi, CosPhi, vel, w, u, gamma;
    //double      m;

    r.x = q->x; r.y = q->y; r.z = q->z;
    Lgm_NormalizeVector( &r );  //r-hat
    
    mInfo->Bfield( q, b, mInfo ); 
    Lgm_NormalizeVector( b );  //b-hat

    Lgm_CrossProduct( &r, b, &e1 ); 
    Lgm_NormalizeVector( &e1 ); // e1-hat

    Lgm_CrossProduct( b, &e1, &e2 ); 
    Lgm_NormalizeVector( &e2 ); // e2-hat


    /*
     * Now construct a moving local frame. It is rotated by angle phi
     * around b-hat. I.e.;
     *
     *      c-hat = -sin(phi) e1-hat - cos(phi) e2_hat
     *      a-hat =  cos(phi) e1-hat - sin(phi) e2_hat
     *
     */
    CosPhi = cos(Phi*RadPerDeg); SinPhi = sin(Phi*RadPerDeg); // these are computed over and over again needlessly.

    // c-hat
    a->x = -SinPhi*e1.x - CosPhi*e2.x;
    a->y = -SinPhi*e1.y - CosPhi*e2.y;
    a->z = -SinPhi*e1.z - CosPhi*e2.z;

    // a-hat
    c->x = CosPhi*e1.x - SinPhi*e2.x;
    c->y = CosPhi*e1.y - SinPhi*e2.y;
    c->z = CosPhi*e1.z - SinPhi*e2.z;

    //m = Mass; // kg
    gamma = 1.0 + Energy/RestEnergy;
    vel = sqrt( 1.0 - 1.0/(gamma*gamma) ) * LGM_c; // m/s
    w = vel*sin(Delta*RadPerDeg); // v_perp, m/s
    u = vel*cos(Delta*RadPerDeg); // v_par, m/s

    // Construct velocity vector
    v->x = u*b->x + w*c->x; // m/s
    v->y = u*b->y + w*c->y; // m/s
    v->z = u*b->z + w*c->z; // m/s


    return;

}


/**
 *  \brief
 *      Compute the gradient of Curl_b at a given point. Note - this is a tensor (and b is b-hat).
 *
 *  \details
 *      Computes \f$ \nabla(\nabla\crossb)\f$.
 *
 *
 *      \param[in]      u0          Position (in GSM) to use.
 *      \param[out]     Grad_Cb     The computed gradient of CurlB at position u0.
 *      \param[in]      DerivScheme Derivative scheme to use (can be one of LGM_DERIV_SIX_POINT, LGM_DERIV_FOUR_POINT, or LGM_DERIV_TWO_POINT).
 *      \param[in]      h           The delta (in Re) to use for grid spacing in the derivative scheme.
 *      \param[in,out]  m           A properly initialized and configured Lgm_MagModelInfo structure.
 *
 *
 *      \author         Mike Henderson
 *      \date           2023
 *
 */
void Lgm_Grad_Cb( Lgm_Vector *u0, double Grad_Cb[3][3], int DerivScheme, double h, Lgm_MagModelInfo *m ) {

    double      H, fx[7][3], fy[7][3], fz[7][3];
    int         i, N;
    Lgm_Vector  u, Curlb;



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

    if (m->VerbosityLevel > 0) printf("\t\tLgm_GradB: Computing GradB with DerivScheme = %d,  h = %g", DerivScheme, h);
    for (i=-N; i<=N; ++i){
        if ( N != 0 ) {
            u = *u0; H = (double)i*h; u.x += H;
            Lgm_Curlb( &u, &Curlb, LGM_DERIV_SIX_POINT, 1e-4, m );
            fx[i+N][0] = Curlb.x; //Re^-1
            fx[i+N][1] = Curlb.y; //Re^-1
            fx[i+N][2] = Curlb.z; //Re^-1
        }
    }
    for (i=-N; i<=N; ++i){
        if ( N != 0 ) {
            u = *u0; H = (double)i*h; u.y += H;
            Lgm_Curlb( &u, &Curlb, LGM_DERIV_SIX_POINT, 1e-4, m );
            fy[i+N][0] = Curlb.x; // Re^-1
            fy[i+N][1] = Curlb.y; // Re^-1
            fy[i+N][2] = Curlb.z; // Re^-1
        }
    }
    for (i=-N; i<=N; ++i){
        if ( N != 0 ) {
            u = *u0; H = (double)i*h; u.z += H;
            Lgm_Curlb( &u, &Curlb, LGM_DERIV_SIX_POINT, 1e-4, m );
            fz[i+N][0] = Curlb.x; // Re^-1
            fz[i+N][1] = Curlb.y; // Re^-1
            fz[i+N][2] = Curlb.z; // Re^-1
        }
    }


    if (DerivScheme == LGM_DERIV_SIX_POINT){
        for (i=0; i<3; i++ ){
//THESE ARE IN  REVERSE ORDER I THNK? CHECK
            Grad_Cb[0][i] = (fx[6][i] - 9.0*fx[5][i] + 45.0*fx[4][i] - 45.0*fx[2][i] + 9.0*fx[1][i] - fx[0][i])/(60.0*h); // Re^-2
            Grad_Cb[1][i] = (fy[6][i] - 9.0*fy[5][i] + 45.0*fy[4][i] - 45.0*fy[2][i] + 9.0*fy[1][i] - fy[0][i])/(60.0*h); // Re^-2
            Grad_Cb[2][i] = (fz[6][i] - 9.0*fz[5][i] + 45.0*fz[4][i] - 45.0*fz[2][i] + 9.0*fz[1][i] - fz[0][i])/(60.0*h); // Re^-2
        }
    } else if (DerivScheme == LGM_DERIV_FOUR_POINT){
        for (i=0; i<3; i++ ){
            Grad_Cb[0][i] = (-fx[4][i] + 8.0*fx[3][i] - 8.0*fx[1][i] + fx[0][i])/(12.0*h); // Re^-2
            Grad_Cb[1][i] = (-fy[4][i] + 8.0*fy[3][i] - 8.0*fy[1][i] + fy[0][i])/(12.0*h); // Re^-2
            Grad_Cb[2][i] = (-fz[4][i] + 8.0*fz[3][i] - 8.0*fz[1][i] + fz[0][i])/(12.0*h); // Re^-2
        }
    } else if (DerivScheme == LGM_DERIV_TWO_POINT){
        for (i=0; i<3; i++ ){
            Grad_Cb[0][i] = (fx[2][i] - fx[0][i])/(2.0*h); // Re^-2
            Grad_Cb[1][i] = (fy[2][i] - fy[0][i])/(2.0*h); // Re^-2
            Grad_Cb[2][i] = (fz[2][i] - fz[0][i])/(2.0*h); // Re^-2
        }
    }
//    if (m->VerbosityLevel > 0) printf("   GradB = (%g %g %g)\n", GradB->x, GradB->y, GradB->z );

    return;

}


/**
 *  \brief
 *      Compute the gradient of Gradient{ b-hat } at a given point. Note - this is a tensor 
 *
 *  \details
 *      Computes \f$ \nabla(\hat{b})\f$.
 *
 *
 *      \param[in]      u0          Position (in GSM) to use.
 *      \param[out]     Grad_b      The computed gradient of b at position u0.
 *      \param[in]      DerivScheme Derivative scheme to use (can be one of LGM_DERIV_SIX_POINT, LGM_DERIV_FOUR_POINT, or LGM_DERIV_TWO_POINT).
 *      \param[in]      h           The delta (in Re) to use for grid spacing in the derivative scheme.
 *      \param[in,out]  m           A properly initialized and configured Lgm_MagModelInfo structure.
 *
 *
 *      \author         Mike Henderson
 *      \date           2023
 *
 */
void Lgm_Grad_b( Lgm_Vector *u0, double Grad_b[3][3], int DerivScheme, double h, Lgm_MagModelInfo *m ) {

    double      H, fx[7][3], fy[7][3], fz[7][3];
    int         i, N;
    Lgm_Vector  u, b;



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

    if (m->VerbosityLevel > 0) printf("\t\tLgm_GradB: Computing GradB with DerivScheme = %d,  h = %g", DerivScheme, h);
    for (i=-N; i<=N; ++i){
        if ( N != 0 ) {
            u = *u0; H = (double)i*h; u.x += H;
            m->Bfield( &u, &b, m );
            Lgm_NormalizeVector( &b );
            fx[i+N][0] = b.x; 
            fx[i+N][1] = b.y; 
            fx[i+N][2] = b.z; 
        }
    }
    for (i=-N; i<=N; ++i){
        if ( N != 0 ) {
            u = *u0; H = (double)i*h; u.y += H;
            m->Bfield( &u, &b, m );
            Lgm_NormalizeVector( &b );
            fy[i+N][0] = b.x; 
            fy[i+N][1] = b.y; 
            fy[i+N][2] = b.z; 
        }
    }
    for (i=-N; i<=N; ++i){
        if ( N != 0 ) {
            u = *u0; H = (double)i*h; u.z += H;
            m->Bfield( &u, &b, m );
            Lgm_NormalizeVector( &b );
            fz[i+N][0] = b.x; 
            fz[i+N][1] = b.y; 
            fz[i+N][2] = b.z; 
        }
    }


    if (DerivScheme == LGM_DERIV_SIX_POINT){
        for (i=0; i<3; i++ ){
            Grad_b[0][i] = (fx[6][i] - 9.0*fx[5][i] + 45.0*fx[4][i] - 45.0*fx[2][i] + 9.0*fx[1][i] - fx[0][i])/(60.0*h); // Re^-1
            Grad_b[1][i] = (fy[6][i] - 9.0*fy[5][i] + 45.0*fy[4][i] - 45.0*fy[2][i] + 9.0*fy[1][i] - fy[0][i])/(60.0*h); // Re^-1
            Grad_b[2][i] = (fz[6][i] - 9.0*fz[5][i] + 45.0*fz[4][i] - 45.0*fz[2][i] + 9.0*fz[1][i] - fz[0][i])/(60.0*h); // Re^-1
        }
    } else if (DerivScheme == LGM_DERIV_FOUR_POINT){
        for (i=0; i<3; i++ ){
            Grad_b[0][i] = (-fx[4][i] + 8.0*fx[3][i] - 8.0*fx[1][i] + fx[0][i])/(12.0*h); // Re^-1
            Grad_b[1][i] = (-fy[4][i] + 8.0*fy[3][i] - 8.0*fy[1][i] + fy[0][i])/(12.0*h); // Re^-1
            Grad_b[2][i] = (-fz[4][i] + 8.0*fz[3][i] - 8.0*fz[1][i] + fz[0][i])/(12.0*h); // Re^-1
        }
    } else if (DerivScheme == LGM_DERIV_TWO_POINT){
        for (i=0; i<3; i++ ){
            Grad_b[0][i] = (fx[2][i] - fx[0][i])/(2.0*h); // Re^-1
            Grad_b[1][i] = (fy[2][i] - fy[0][i])/(2.0*h); // Re^-1
            Grad_b[2][i] = (fz[2][i] - fz[0][i])/(2.0*h); // Re^-1
        }
    }
//    if (m->VerbosityLevel > 0) printf("   GradB = (%g %g %g)\n", GradB->x, GradB->y, GradB->z );

    return;

}

/**
 *  \brief
 *      Compute the gradient of Gradient{ GradB } at a given point. Note - this is a tensor 
 *
 *  \details
 *      Computes \f$ \nabla( \nabla B )\f$.
 *
 *
 *      \param[in]      u0          Position (in GSM) to use.
 *      \param[out]     Grad_GradB  The computed gradient of GradB at position u0.
 *      \param[in]      DerivScheme Derivative scheme to use (can be one of LGM_DERIV_SIX_POINT, LGM_DERIV_FOUR_POINT, or LGM_DERIV_TWO_POINT).
 *      \param[in]      h           The delta (in Re) to use for grid spacing in the derivative scheme.
 *      \param[in,out]  m           A properly initialized and configured Lgm_MagModelInfo structure.
 *
 *
 *      \author         Mike Henderson
 *      \date           2023
 *
 */
void Lgm_Grad_GradB( Lgm_Vector *u0, double Grad_GradB[3][3], int DerivScheme, double h, Lgm_MagModelInfo *m ) {

    double      H, fx[7][3], fy[7][3], fz[7][3];
    int         i, N;
    Lgm_Vector  u, GradB;



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

    if (m->VerbosityLevel > 0) printf("\t\tLgm_GradB: Computing GradB with DerivScheme = %d,  h = %g", DerivScheme, h);
    for (i=-N; i<=N; ++i){
        if ( N != 0 ) {
            u = *u0; H = (double)i*h; u.x += H;
            Lgm_GradB( &u, &GradB, LGM_DERIV_SIX_POINT, 1e-4, m );
            fx[i+N][0] = GradB.x; 
            fx[i+N][1] = GradB.y; 
            fx[i+N][2] = GradB.z; 
        }
    }
    for (i=-N; i<=N; ++i){
        if ( N != 0 ) {
            u = *u0; H = (double)i*h; u.y += H;
            Lgm_GradB( &u, &GradB, LGM_DERIV_SIX_POINT, 1e-4, m );
            fy[i+N][0] = GradB.x; 
            fy[i+N][1] = GradB.y; 
            fy[i+N][2] = GradB.z; 
        }
    }
    for (i=-N; i<=N; ++i){
        if ( N != 0 ) {
            u = *u0; H = (double)i*h; u.z += H;
            Lgm_GradB( &u, &GradB, LGM_DERIV_SIX_POINT, 1e-4, m );
            fz[i+N][0] = GradB.x; 
            fz[i+N][1] = GradB.y; 
            fz[i+N][2] = GradB.z; 
        }
    }


    if (DerivScheme == LGM_DERIV_SIX_POINT){
        for (i=0; i<3; i++ ){
            Grad_GradB[0][i] = (fx[6][i] - 9.0*fx[5][i] + 45.0*fx[4][i] - 45.0*fx[2][i] + 9.0*fx[1][i] - fx[0][i])/(60.0*h); // Re^-2
            Grad_GradB[1][i] = (fy[6][i] - 9.0*fy[5][i] + 45.0*fy[4][i] - 45.0*fy[2][i] + 9.0*fy[1][i] - fy[0][i])/(60.0*h); // Re^-2
            Grad_GradB[2][i] = (fz[6][i] - 9.0*fz[5][i] + 45.0*fz[4][i] - 45.0*fz[2][i] + 9.0*fz[1][i] - fz[0][i])/(60.0*h); // Re^-2
        }
    } else if (DerivScheme == LGM_DERIV_FOUR_POINT){
        for (i=0; i<3; i++ ){
            Grad_GradB[0][i] = (-fx[4][i] + 8.0*fx[3][i] - 8.0*fx[1][i] + fx[0][i])/(12.0*h); // Re^-2
            Grad_GradB[1][i] = (-fy[4][i] + 8.0*fy[3][i] - 8.0*fy[1][i] + fy[0][i])/(12.0*h); // Re^-2
            Grad_GradB[2][i] = (-fz[4][i] + 8.0*fz[3][i] - 8.0*fz[1][i] + fz[0][i])/(12.0*h); // Re^-2
        }
    } else if (DerivScheme == LGM_DERIV_TWO_POINT){
        for (i=0; i<3; i++ ){
            Grad_GradB[0][i] = (fx[2][i] - fx[0][i])/(2.0*h); // Re^-2
            Grad_GradB[1][i] = (fy[2][i] - fy[0][i])/(2.0*h); // Re^-2
            Grad_GradB[2][i] = (fz[2][i] - fz[0][i])/(2.0*h); // Re^-2
        }
    }
//    if (m->VerbosityLevel > 0) printf("   GradB = (%g %g %g)\n", GradB->x, GradB->y, GradB->z );

    return;

}




/**
 *  \brief
 *      Compute the d{b dot v}/dcomp at a given point.
 *
 *  \details
 *      Computes \f$ \partial (\hat{b}\cdot\vec{v})/\partial comp \f$.
 *
 *
 *      \param[in]      u0          Position (in GSM) to use.
 *      \param[in]      v           Velocity vector (in GSM) to use.
 *      \param[in]      comp        seclects component (0=x, 1=y, 2=z)
 *      \param[out]     b dot v     The computed partial deriviative { d dot v}/dcomp at position u0.
 *      \param[in]      DerivScheme Derivative scheme to use (can be one of LGM_DERIV_SIX_POINT, LGM_DERIV_FOUR_POINT, or LGM_DERIV_TWO_POINT).
 *      \param[in]      h           The delta (in Re) to use for grid spacing in the derivative scheme.
 *      \param[in,out]  m           A properly initialized and configured Lgm_MagModelInfo structure.
 *
 *
 *      \author         Mike Henderson
 *      \date           2023
 *
 */
void Lgm_dbdotvdcomp( Lgm_Vector *u0, Lgm_Vector *v, int comp, double *dbdotvdcomp, int DerivScheme, double h, Lgm_MagModelInfo *m ) {

    double      bdotv, H, fc[7];
    int         i, N=3;
    Lgm_Vector  u, b;

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

    if (m->VerbosityLevel > 0) printf("\t\tLgm_dbdotvdcomp: Computing dbdotvdcomp for comp=%d with DerivScheme = %d,  h = %g", comp, DerivScheme, h);
    for (i=-N; i<=N; ++i){
        if ( N != 0 ) {
            u = *u0; H = (double)i*h; 
            if      ( comp == 0 ) { u.x += H; }
            else if ( comp == 1 ) { u.y += H; }
            else if ( comp == 2 ) { u.z += H; }
            else { printf("\t\tLgm_dbdotvdcomp: Invalid component. comp = %d (must be 0, 1, or 2).\n", comp); exit(0);}
            m->Bfield( &u, &b, m );
            Lgm_NormalizeVector( &b );  //b-hat
            bdotv = Lgm_DotProduct( &b, v );

            fc[i+N] = bdotv;
        }
    }

    if (DerivScheme == LGM_DERIV_SIX_POINT){
        *dbdotvdcomp = (fc[6] - 9.0*fc[5] + 45.0*fc[4] - 45.0*fc[2] + 9.0*fc[1] - fc[0])/(60.0*h);
    } else if (DerivScheme == LGM_DERIV_FOUR_POINT){
        *dbdotvdcomp = (-fc[4] + 8.0*fc[3] - 8.0*fc[1] + fc[0])/(12.0*h);
    } else if (DerivScheme == LGM_DERIV_TWO_POINT){
        *dbdotvdcomp = (fc[2] - fc[0])/(2.0*h);
    }
    if (m->VerbosityLevel > 0) printf("   dbdotvdcomp = %g (component = %d)\n", *dbdotvdcomp, comp);

    return;

}
void Lgm_dbdotvdx( Lgm_Vector *u0, Lgm_Vector *v, double *dbdotvdx, int DerivScheme, double h, Lgm_MagModelInfo *m ) { Lgm_dbdotvdcomp( u0, v, 0, dbdotvdx, DerivScheme, h, m ); }
void Lgm_dbdotvdy( Lgm_Vector *u0, Lgm_Vector *v, double *dbdotvdy, int DerivScheme, double h, Lgm_MagModelInfo *m ) { Lgm_dbdotvdcomp( u0, v, 1, dbdotvdy, DerivScheme, h, m ); }
void Lgm_dbdotvdz( Lgm_Vector *u0, Lgm_Vector *v, double *dbdotvdz, int DerivScheme, double h, Lgm_MagModelInfo *m ) { Lgm_dbdotvdcomp( u0, v, 2, dbdotvdz, DerivScheme, h, m ); }
/**
 *  \brief
 *      Compute the Laplacian of b dot v at a given point.
 *
 *  \details
 *      Computes \f$ \nabla^2 (\hat{b}\cdot\vec{v}) \f$.
 *
 *
 *      \param[in]      u0             Position (in GSM) to use.
 *      \param[in]      v              Velocity vector (in GSM) to use.
 *      \param[out]     Laplacian_bdotv  The computed Laplacian of bdotv at position u0.
 *      \param[in]      DerivScheme    Derivative scheme to use (can be one of LGM_DERIV_SIX_POINT, LGM_DERIV_FOUR_POINT, or LGM_DERIV_TWO_POINT).
 *      \param[in]      h              The delta (in Re) to use for grid spacing in the derivative scheme.
 *      \param[in,out]  m              A properly initialized and configured Lgm_MagModelInfo structure.
 *
 *
 *      \author         Mike Henderson
 *      \date           2023
 *
 */
void Lgm_Laplacian_bdotv( Lgm_Vector *u0, Lgm_Vector *v, double *Laplacian_bdotv, int DerivScheme, double h, Lgm_MagModelInfo *m ) {

    double      H, d2bdotvdx2, d2bdotvdy2, d2bdotvdz2;
    double      dbdotvdx, dbdotvdy, dbdotvdz;
    double      fx[7], gy[7], hz[7];
    int         i, N=3;
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
        default:
            N = 3;
    }
   

    dbdotvdx = dbdotvdy = dbdotvdz = 0.0;
    if (m->VerbosityLevel > 0) printf("\t\tLgm_GradB: Computing GradB with DerivScheme = %d,  h = %g", DerivScheme, h);
    for (i=-N; i<=N; ++i){
        if ( N != 0 ) {
            u = *u0; H = (double)i*h; u.x += H;
            Lgm_dbdotvdx( &u, v, &dbdotvdx, DerivScheme, h, m );
            fx[i+N] = dbdotvdx;
        }
    }
    for (i=-N; i<=N; ++i){
        if ( N != 0 ) {
            u = *u0; H = (double)i*h; u.y += H;
            Lgm_dbdotvdy( &u, v, &dbdotvdy, DerivScheme, h, m );
            gy[i+N] = dbdotvdy;
        }
    }
    for (i=-N; i<=N; ++i){
        if ( N != 0 ) {
            u = *u0; H = (double)i*h; u.z += H;
            Lgm_dbdotvdz( &u, v, &dbdotvdz, DerivScheme, h, m );
            hz[i+N] = dbdotvdz;
        }
    }


    d2bdotvdx2 = d2bdotvdy2 = d2bdotvdz2 = 0.0;
    if (DerivScheme == LGM_DERIV_SIX_POINT){
        d2bdotvdx2 = (fx[6] - 9.0*fx[5] + 45.0*fx[4] - 45.0*fx[2] + 9.0*fx[1] - fx[0])/(60.0*h);
        d2bdotvdy2 = (gy[6] - 9.0*gy[5] + 45.0*gy[4] - 45.0*gy[2] + 9.0*gy[1] - gy[0])/(60.0*h);
        d2bdotvdz2 = (hz[6] - 9.0*hz[5] + 45.0*hz[4] - 45.0*hz[2] + 9.0*hz[1] - hz[0])/(60.0*h);
    } else if (DerivScheme == LGM_DERIV_FOUR_POINT){
        d2bdotvdx2 = (-fx[4] + 8.0*fx[3] - 8.0*fx[1] + fx[0])/(12.0*h);
        d2bdotvdy2 = (-gy[4] + 8.0*gy[3] - 8.0*gy[1] + gy[0])/(12.0*h);
        d2bdotvdz2 = (-hz[4] + 8.0*hz[3] - 8.0*hz[1] + hz[0])/(12.0*h);
    } else if (DerivScheme == LGM_DERIV_TWO_POINT){
        d2bdotvdx2 = (hz[2] - hz[0])/(2.0*h);
        d2bdotvdy2 = (hz[2] - hz[0])/(2.0*h);
        d2bdotvdz2 = (hz[2] - hz[0])/(2.0*h);
    }
    *Laplacian_bdotv = d2bdotvdx2 + d2bdotvdy2 + d2bdotvdz2;
    if (m->VerbosityLevel > 0) printf("   Laplacian_bdotv = %g\n", *Laplacian_bdotv );

    return;

}


/**
 *  \brief
 *      Compute the gradient of (b-hat dot v) at a given point. Note: This
 *      routine is designed to implement the Burby 1st invariant calculation
 *      and in those formulea, the v is to be taken as a constant with respect
 *      to derivatives. This could also be implemented using existing routines
 *      via some vector/tensor identities, but here I am just sticking literally
 *      to the Burby terms.
 *
 *  \details
 *      Computes \f$ \nabla( \hat{B}\cdot\vec{v} ) \f$.
 *
 *
 *      \param[in]      u0          Position (in GSM) to use.
 *      \param[in]      v           Velocity vector.
 *      \param[out]     Grad_bdotv  The computed gradient of bdotv at position u0.
 *      \param[in]      DerivScheme Derivative scheme to use (can be one of LGM_DERIV_SIX_POINT, LGM_DERIV_FOUR_POINT, or LGM_DERIV_TWO_POINT).
 *      \param[in]      h           The delta (in Re) to use for grid spacing in the derivative scheme.
 *      \param[in,out]  m           A properly initialized and configured Lgm_MagModelInfo structure.
 *
 *
 *      \author         Mike Henderson
 *      \date           2023
 *
 */
void Lgm_Grad_bdotv( Lgm_Vector *u0, Lgm_Vector *v, Lgm_Vector *Grad_bdotv, int DerivScheme, double h, Lgm_MagModelInfo *m ) {

    double      H, fx[7], fy[7], fz[7];
    int         i, N=3;
    Lgm_Vector  u, b;



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

    if (m->VerbosityLevel > 0) printf("\t\tLgm_GradB: Computing GradB with DerivScheme = %d,  h = %g", DerivScheme, h);
    for (i=-N; i<=N; ++i){
        if ( N != 0 ) {
            u = *u0; H = (double)i*h; u.x += H;
            m->Bfield( &u, &b, m );
            Lgm_NormalizeVector( &b );
            fx[i+N] = Lgm_DotProduct( &b, v );
        }
    }
    for (i=-N; i<=N; ++i){
        if ( N != 0 ) {
            u = *u0; H = (double)i*h; u.y += H;
            m->Bfield( &u, &b, m );
            Lgm_NormalizeVector( &b );
            fy[i+N] = Lgm_DotProduct( &b, v );
        }
    }
    for (i=-N; i<=N; ++i){
        if ( N != 0 ) {
            u = *u0; H = (double)i*h; u.z += H;
            m->Bfield( &u, &b, m );
            Lgm_NormalizeVector( &b );
            fz[i+N] = Lgm_DotProduct( &b, v );
        }
    }


    if (DerivScheme == LGM_DERIV_SIX_POINT){
        Grad_bdotv->x = (fx[6] - 9.0*fx[5] + 45.0*fx[4] - 45.0*fx[2] + 9.0*fx[1] - fx[0])/(60.0*h);
        Grad_bdotv->y = (fy[6] - 9.0*fy[5] + 45.0*fy[4] - 45.0*fy[2] + 9.0*fy[1] - fy[0])/(60.0*h);
        Grad_bdotv->z = (fz[6] - 9.0*fz[5] + 45.0*fz[4] - 45.0*fz[2] + 9.0*fz[1] - fz[0])/(60.0*h);
    } else if (DerivScheme == LGM_DERIV_FOUR_POINT){
        Grad_bdotv->x = (-fx[4] + 8.0*fx[3] - 8.0*fx[1] + fx[0])/(12.0*h);
        Grad_bdotv->y = (-fy[4] + 8.0*fy[3] - 8.0*fy[1] + fy[0])/(12.0*h);
        Grad_bdotv->z = (-fz[4] + 8.0*fz[3] - 8.0*fz[1] + fz[0])/(12.0*h);
    } else if (DerivScheme == LGM_DERIV_TWO_POINT){
        Grad_bdotv->x = (fx[2] - fx[0])/(2.0*h);
        Grad_bdotv->y = (fy[2] - fy[0])/(2.0*h);
        Grad_bdotv->z = (fz[2] - fz[0])/(2.0*h);
    }
    if (m->VerbosityLevel > 0) printf("   Grad_bdotv = (%g %g %g)\n", Grad_bdotv->x, Grad_bdotv->y, Grad_bdotv->z );

    return;

}

/**
 *  \brief
 *      Compute the gradient of Gradient{ GradB } at a given point. Note - this is a tensor 
 *
 *  \details
 *      Computes \f$ \nabla( \nabla B )\f$.
 *
 *
 *      \param[in]      u0               Position (in GSM) to use.
 *      \param[out]     Grad_Grad_bdotv  The computed gradient of GradB at position u0.
 *      \param[in]      DerivScheme      Derivative scheme to use (can be one of LGM_DERIV_SIX_POINT, LGM_DERIV_FOUR_POINT, or LGM_DERIV_TWO_POINT).
 *      \param[in]      h                The delta (in Re) to use for grid spacing in the derivative scheme.
 *      \param[in,out]  m                A properly initialized and configured Lgm_MagModelInfo structure.
 *
 *
 *      \author         Mike Henderson
 *      \date           2023
 *
 */
void Lgm_Grad_Grad_bdotv( Lgm_Vector *u0, Lgm_Vector *v, double Grad_Grad_bdotv[3][3], int DerivScheme, double h, Lgm_MagModelInfo *m ) {

    double      H, fx[7][3], fy[7][3], fz[7][3];
    int         i, N;
    Lgm_Vector  u, Grad_bdotv;



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

    if (m->VerbosityLevel > 0) printf("\t\tLgm_GradB: Computing GradB with DerivScheme = %d,  h = %g", DerivScheme, h);
    for (i=-N; i<=N; ++i){
        if ( N != 0 ) {
            u = *u0; H = (double)i*h; u.x += H;
            Lgm_Grad_bdotv( &u, v, &Grad_bdotv, LGM_DERIV_SIX_POINT, 1e-4, m );
            fx[i+N][0] = Grad_bdotv.x; 
            fx[i+N][1] = Grad_bdotv.y; 
            fx[i+N][2] = Grad_bdotv.z; 
        }
    }
    for (i=-N; i<=N; ++i){
        if ( N != 0 ) {
            u = *u0; H = (double)i*h; u.y += H;
            Lgm_Grad_bdotv( &u, v, &Grad_bdotv, LGM_DERIV_SIX_POINT, 1e-4, m );
            fy[i+N][0] = Grad_bdotv.x; 
            fy[i+N][1] = Grad_bdotv.y; 
            fy[i+N][2] = Grad_bdotv.z; 
        }
    }
    for (i=-N; i<=N; ++i){
        if ( N != 0 ) {
            u = *u0; H = (double)i*h; u.z += H;
            Lgm_Grad_bdotv( &u, v, &Grad_bdotv, LGM_DERIV_SIX_POINT, 1e-4, m );
            fz[i+N][0] = Grad_bdotv.x; 
            fz[i+N][1] = Grad_bdotv.y; 
            fz[i+N][2] = Grad_bdotv.z; 
        }
    }


    if (DerivScheme == LGM_DERIV_SIX_POINT){
        for (i=0; i<3; i++ ){
            Grad_Grad_bdotv[0][i] = (fx[6][i] - 9.0*fx[5][i] + 45.0*fx[4][i] - 45.0*fx[2][i] + 9.0*fx[1][i] - fx[0][i])/(60.0*h); // Re^-2
            Grad_Grad_bdotv[1][i] = (fy[6][i] - 9.0*fy[5][i] + 45.0*fy[4][i] - 45.0*fy[2][i] + 9.0*fy[1][i] - fy[0][i])/(60.0*h); // Re^-2
            Grad_Grad_bdotv[2][i] = (fz[6][i] - 9.0*fz[5][i] + 45.0*fz[4][i] - 45.0*fz[2][i] + 9.0*fz[1][i] - fz[0][i])/(60.0*h); // Re^-2
        }
    } else if (DerivScheme == LGM_DERIV_FOUR_POINT){
        for (i=0; i<3; i++ ){
            Grad_Grad_bdotv[0][i] = (-fx[4][i] + 8.0*fx[3][i] - 8.0*fx[1][i] + fx[0][i])/(12.0*h); // Re^-2
            Grad_Grad_bdotv[1][i] = (-fy[4][i] + 8.0*fy[3][i] - 8.0*fy[1][i] + fy[0][i])/(12.0*h); // Re^-2
            Grad_Grad_bdotv[2][i] = (-fz[4][i] + 8.0*fz[3][i] - 8.0*fz[1][i] + fz[0][i])/(12.0*h); // Re^-2
        }
    } else if (DerivScheme == LGM_DERIV_TWO_POINT){
        for (i=0; i<3; i++ ){
            Grad_Grad_bdotv[0][i] = (fx[2][i] - fx[0][i])/(2.0*h); // Re^-2
            Grad_Grad_bdotv[1][i] = (fy[2][i] - fy[0][i])/(2.0*h); // Re^-2
            Grad_Grad_bdotv[2][i] = (fz[2][i] - fz[0][i])/(2.0*h); // Re^-2
        }
    }
//    if (m->VerbosityLevel > 0) printf("   GradB = (%g %g %g)\n", GradB->x, GradB->y, GradB->z );

    return;

}



void Lgm_Dyad( Lgm_Vector *a, Lgm_Vector *b, double r[3][3] ) {

    r[0][0] = a->x*b->x; r[0][1] = a->x*b->y; r[0][2] = a->x*b->z;
    r[1][0] = a->y*b->x; r[1][1] = a->y*b->y; r[1][2] = a->y*b->z;
    r[2][0] = a->z*b->x; r[2][1] = a->z*b->y; r[2][2] = a->z*b->z;

    return;
}

void Lgm_VectorDotTensor( Lgm_Vector *v, double T[3][3], Lgm_Vector *Result ) {
    Result->x = v->x*T[0][0]  +  v->y*T[1][0]  +  v->z*T[2][0];
    Result->y = v->x*T[0][1]  +  v->y*T[1][1]  +  v->z*T[2][1];
    Result->z = v->x*T[0][2]  +  v->y*T[1][2]  +  v->z*T[2][2];
}

void Lgm_TensorDotVector( double T[3][3], Lgm_Vector *v, Lgm_Vector *Result ) {
    Result->x = T[0][0]*v->x  +  T[0][1]*v->y  +  T[0][2]*v->z;
    Result->y = T[1][0]*v->x  +  T[1][1]*v->y  +  T[1][2]*v->z;
    Result->z = T[2][0]*v->x  +  T[2][1]*v->y  +  T[2][2]*v->z;
}
void Lgm_TensorDotTensor( double A[3][3], double B[3][3], double R[3][3] ) {
    int i, j, k;
    for ( i=0; i<3; i++ ) {
        for ( j=0; j<3; j++ ) {
            R[i][j] = 0.0;
            for (k=0; k<3; k++) R[i][j] += A[i][k]*B[k][j];

        }
    }
}
double Lgm_TraceTensor( double T[3][3] ) {
    return( T[0][0] + T[1][1] + T[2][2] );
}
void Lgm_TensorTranspose( double T[3][3], double T_transpose[3][3] ) {
    int i, j;
    for ( i=0; i<3; i++ ) {
        for ( j=0; j<3; j++ ) {
            T_transpose[i][j] = T[j][i];
        }
    }
}
double Lgm_DoubleDot( double A[3][3], double B[3][3] ) {
    int i, j;
    double sum;
    for (sum=0.0, i=0; i<3; i++ ) {
        for ( j=0; j<3; j++ ) {
            sum += A[i][j]*B[i][j];
        }
    }
    return( sum );
}



double Lgm_Burby( Lgm_Vector *q, Lgm_Vector *v, double Gamma, double Mass, double Charge, double *mu0_out, double *mu1_out, double *mu2_out, Lgm_MagModelInfo *mInfo ){


    Lgm_Vector  b;              // b-hat(q)
    Lgm_Vector  GradB;          // Grad{ |B-vec| }(q)
    double      B;              // |B-vec|((q)
    Lgm_Vector  L, N, K, H;     // Burby's "lembda, eta, Kappa" vectors
    double      Grad_b[3][3];   // tensor
    Lgm_Vector  v_cross_b, b_cross_K, GradB_cross_b;
    double      G, M, f1, f2, f3, f4, m, e, eps, eps2, mu, mu0, mu1, mu2;
    int         i, j;
    double      B2, B3, v2, u, u2, u3, u4, w2, g, g2;

    e = Charge; //Coulombs (or s·A)
    //c = LGM_c;  // Speed of light  m/s
    m = Mass;   // Mass of particle  kg

    mInfo->Bfield( q, &b, mInfo );  // Get B-vec, nT
    B = Lgm_NormalizeVector(&b);    // Normalize vector and return the magntiude as B (in nT)
    B *= 1e-9; // T
    B2 = B*B;
    B3 = B2*B;

    Lgm_GradB( q, &GradB, LGM_DERIV_SIX_POINT, 1e-5, mInfo ); // Compute gradient of B-field at position q
    g = 1e-9/(Re*1e3);
    GradB.x *= g; // T/m
    GradB.y *= g; // T/m
    GradB.z *= g; // T/m

    // Compute Grad{ ln|B| } = Grad{B}/B
    H = GradB; Lgm_ScaleVector( &H, 1.0/B );
    
    // Compute G (Burby's "gamma" scalar. Note that Grad{ln|B|} = Grad{|B|}/B )
    G = Lgm_DotProduct( v, &GradB )/B; // 1/s

    // Compute M (Burby's "mu" scalar)
    Lgm_CrossProduct( v, &b, &v_cross_b );
    M = 0.5*Lgm_DotProduct( &v_cross_b, &v_cross_b )/B; // m^2/s^2 / T

    // Compute u (Burby's "v_par" scalar)
    u = Lgm_DotProduct( &b, v );
    u2 = u*u; // m^2/s^2
    u3 = u2*u;
    u4 = u2*u2;
    v2 = Lgm_DotProduct( v, v );
    w2 = v2 - u2;
    eps  = Gamma*m/e; // Gamma is given to us.
    eps2 = eps*eps;

    // Compute L (Burby's "lambda" vector)
    g = 1.0/(Re*1e3);
    Lgm_Grad_b( q, Grad_b, LGM_DERIV_SIX_POINT, 1e-5, mInfo );
    for (i=0;i<3;i++){ for (j=0;j<3;j++){ Grad_b[i][j] *= g; } } // 1/m
    Lgm_VectorDotTensor( v, Grad_b, &L );  // This is a vector dotted with a tensor = vector

    // Compute N (Burby's "eta" vector)
    Lgm_TensorDotVector( Grad_b, v, &N );  // This is a tensor dotted vector = vector

    // Compute K (Burby's "kappa" vector = b dot Grad(b-hat}
    Lgm_VectorDotTensor( &b, Grad_b, &K );  // This is a vector dotted with a tensor = vector



    /*
     *  mu0 
     */
    mu0 = Gamma*m*w2/(2.0*B); // J/T

    
    // Compute terms in brackets in Burby's eqn 29
    f1 =  0.25*u  *  Lgm_DotProduct( &L, &v_cross_b );
    f2 = -0.75*u  * Lgm_DotProduct( &v_cross_b, &N );
    Lgm_CrossProduct( &b, &K, &b_cross_K );
    f3 = -1.25*u2 * Lgm_DotProduct( &b_cross_K, v );
    Lgm_CrossProduct( &GradB, &b, &GradB_cross_b );
    f4 = M * Lgm_DotProduct( &GradB_cross_b, v );
    //mu1 = 0.5*Gamma*m*(f1 + f2 + f3 + f4)/(B*B);
    mu1 = Gamma*m*(f1 + f2 + f3 + f4)/(B*B); // seems to be off by exactly 2


    //printf("BURBY: mu1 = %g\n", mu1);
//    printf("BURBY: mu1 = %.10g MeV/G\n", eps*mu1 * 6.242e12 / 1.0e4);




    /* 
     * mu2
     */
    mu2  = -71.0*M* Lgm_DotProduct(&N, &N)/(128.0*B2);
    mu2 += -5.0*u*M * Lgm_DotProduct(&L, &H)/(3.0*B2);
    mu2 += 10.0*u*M*Lgm_DotProduct(&H, &N)/(3.0*B2);
    mu2 += -715.0*u*M*Lgm_DotProduct(&K, &N)/(192.0*B2);
    double Kdotv;
    Kdotv = Lgm_DotProduct(&K, v);
    mu2 += 71.0*M*Kdotv*Kdotv/(128.0*B2);
    mu2 += -3.0*M*G*G/(2.0*B2);
    mu2 += 71.0*u*M*Lgm_DotProduct(&L, &K)/(64.0*B2);
    double LdotL;
    LdotL = Lgm_DotProduct( &L, &L );
    mu2 += -7.0*M*LdotL/(128.0*B2);
    double LdotN;
    LdotN = Lgm_DotProduct( &L, &N );
    mu2 += 25.0*M*LdotN/(64.0*B2);
    double Ldotv;
    Ldotv = Lgm_DotProduct( &L,  v );
    mu2 += 5.0*u*Kdotv*Ldotv/(12.0*B3);
    mu2 += -2.0*u*G*Ldotv/(3.0*B3);
    mu2 += Ldotv*Ldotv/(8.0*B3);
    double Divb;
    Lgm_Divb( q, &Divb, LGM_DERIV_SIX_POINT, 1e-5, mInfo ); Divb /= (Re*1e3);
    mu2 += 217.0*u*M*Kdotv*Divb/(32.0*B2);
    mu2 += -5.0*u*M*G*Divb/(3.0*B2);
    mu2 += 23.0*M*Ldotv*Divb/(32.0*B2);

    double Laplacian_bdotv;
    Lgm_Laplacian_bdotv( q, v, &Laplacian_bdotv, LGM_DERIV_SIX_POINT, 1e-5, mInfo );
    g = 1.0/(Re*1e3); g2 = g*g; Laplacian_bdotv *= g2;
    mu2 += -5.0*u*M*Laplacian_bdotv/(6.0*B2);

    Lgm_Vector Grad_Divb;
    Lgm_Grad_Divb( q, &Grad_Divb, LGM_DERIV_SIX_POINT, 1e-5, mInfo );
    g = 1.0/(Re*1e3); g2 = g*g; Grad_Divb.x *= g2; Grad_Divb.y *= g2; Grad_Divb.z *= g2;
    mu2 += 5.0*u*M*Lgm_DotProduct( v, &Grad_Divb )/(6.0*B2);

    double vv[3][3];
    double Grad_GradB[3][3];
    Lgm_Dyad( v, v, vv );
    Lgm_Grad_GradB( q, Grad_GradB, LGM_DERIV_SIX_POINT, 1e-5, mInfo );
    g = 1.0/(Re*1e3); g2 = 1e-9*g*g;
    for (i=0;i<3;i++){ for (j=0;j<3;j++){ Grad_GradB[i][j] *= g2; } } // 1/m
    mu2 += M/(2.0*B3) * Lgm_DoubleDot( vv, Grad_GradB );

    double bb[3][3];
    double Grad_Grad_bdotv[3][3];
    Lgm_Dyad( &b, &b, bb );
    Lgm_Grad_Grad_bdotv( q, v, Grad_Grad_bdotv, LGM_DERIV_SIX_POINT, 1e-5, mInfo );
    g = 1.0/(Re*1e3); g2 = g*g;
    for (i=0;i<3;i++){ for (j=0;j<3;j++){ Grad_Grad_bdotv[i][j] *= g2; } } // 1/m
    mu2 += 5.0*u*M/(6.0*B2) * Lgm_DoubleDot( bb, Grad_Grad_bdotv );
    mu2 += u/(6.0*B3) * Lgm_DoubleDot( vv, Grad_Grad_bdotv );
    mu2 += -25.0*u2*M/(12.0*B2) * Lgm_DotProduct( &b, &Grad_Divb );

    double NdotN;
    NdotN = Lgm_DotProduct( &N, &N );
    mu2 += -5.0*u2*NdotN/(8.0*B3);
    mu2 += 20.0*u2*M*Lgm_DotProduct( &K, &H )/(3.0*B2);
    double KdotK = Lgm_DotProduct( &K, &K );
    mu2 += -3013.0*u2*M*KdotK/(384.0*B2);
    mu2 += 5.0*u2*Kdotv*Kdotv/(6.0*B3);
    mu2 += -5.0*u2*Kdotv*G/(6.0*B3);
    mu2 += -29.0*u2*LdotL/(24.0*B3);
    mu2 += 5.0*u2*LdotN/(2.0*B3);
    mu2 += -5.0*u2*Ldotv*Divb/(12.0*B3);
    mu2 += -25.0*u2*M*Divb*Divb/(24.0*B2);

    double Grad_b_T[3][3];
    double Gradb_dot_GradbT[3][3];
    double Trace_Gradb_dot_GradbT;
    Lgm_TensorTranspose( Grad_b, Grad_b_T );
    Lgm_TensorDotTensor( Grad_b, Grad_b_T, Gradb_dot_GradbT );
    Trace_Gradb_dot_GradbT = Lgm_TraceTensor( Gradb_dot_GradbT );
    //printf("Trace_Gradb_dot_GradbT = %g\n", Trace_Gradb_dot_GradbT);
    mu2 += 55.0*u2*M/(24.0*B2)*Trace_Gradb_dot_GradbT;

    double Gradb_dot_Gradb[3][3];
    double Trace_Gradb_dot_Gradb;
    Lgm_TensorDotTensor( Grad_b, Grad_b, Gradb_dot_Gradb );
    Trace_Gradb_dot_Gradb = Lgm_TraceTensor( Gradb_dot_Gradb );
    //printf("Trace_Gradb_dot_Gradb = %g\n", Trace_Gradb_dot_Gradb);
    mu2 += -35.0*u2*M/(8.0*B2)*Trace_Gradb_dot_Gradb;

    // get bv dyad
    double bv[3][3];
    Lgm_Dyad( &b, v, bv );
    mu2 += 5.0*u2/(12.0*B3) * Lgm_DoubleDot( bv, Grad_Grad_bdotv );

    double LdotK;
    LdotK = Lgm_DotProduct( &L, &K );
    mu2 += 5.0*u3*LdotK/(3.0*B3);
    mu2 += 5.0*u3*Kdotv*Divb/(12.0*B3);
    mu2 += 5.0*u3/(12.0*B3) * Lgm_DoubleDot( bb, Grad_Grad_bdotv );

    mu2 += 25.0*u4*KdotK/(24.0*B3);
    mu2 += -5.0*M*M/(4.0*B)*Lgm_DotProduct( &b, &Grad_Divb );
    double LaplacianB;
    Lgm_LaplacianB( q, &LaplacianB, LGM_DERIV_SIX_POINT, 1e-5, mInfo );
    g = 1.0/(Re*1e3); g2 = 1e-9*g*g; LaplacianB *= g2;
    mu2 += -5.0*M*M*LaplacianB/(4.0*B2);
    mu2 += 15.0*M*M* Lgm_DotProduct( &H, &H )/(4.0*B);
    mu2 += -5.0*M*M*Lgm_DotProduct( &K, &H )/(4.0*B);
    mu2 += -11.0*M*M*KdotK/(64.0*B);
    mu2 += -133.0*M*M*Divb*Divb/(32.0*B);
    mu2 += 11.0*M*M/(64.0*B)*Trace_Gradb_dot_GradbT;
    mu2 += -5.0*M*M/(64.0*B)*Trace_Gradb_dot_Gradb;
    
    //mu2 *= 0.5*Gamma*m;
    mu2 *= Gamma*m;

//    printf("BURBY: mu2 = %g MeV/G\n", eps2*mu2 * 6.242e12 / 1.0e4);

//    printf("BURBY: mu  = %g MeV/G\n", (mu0 + eps*mu1 + eps2*mu2) * 6.242e12 / 1.0e4);


    mu = (mu0 + eps*mu1 + eps2*mu2) * 6.242e12 / 1.0e4;

    *mu0_out = mu0*6.242e12/1.0e4;
    *mu1_out = eps*mu1*6.242e12/1.0e4;
    *mu2_out = eps2*mu2*6.242e12/1.0e4;

    return( mu );
    
}

