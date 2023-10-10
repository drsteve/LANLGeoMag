#include "Lgm_ParticleInfo.h"


/*
 * 6D system for Full Orbit equations of motion.
 */
int f_full( double t, const double y[], double dydt[], void *params ) {

    double  g, gamma;
    Lgm_Vector E, B, uxB, P, u;
    Lgm_ParticleInfo *p;
    p = (Lgm_ParticleInfo *)params;

    /*
     *  Electric field
     */
    E.x = E.y = E.z = 0.0;



    /*
     *  The position vector (in Re)
     */
    P.x = y[0]/1000.0/Re; P.y = y[1]/1000.0/Re; P.z = y[2]/1000.0/Re;

    /*
     *  The velocity vector (m/s)
     */
    u.x = y[3]; u.y = y[4]; u.z = y[5];
    p->gamma = sqrt( 1.0 + (u.x*u.x + u.y*u.y + u.z*u.z)/p->c2 );
    

    /*
     *  Compute uxB
     */
    p->mInfo->Bfield( &P, &B, p->mInfo );  // B in nT
    B.x *= 1e-9; B.y *= 1e-9; B.z *= 1e-9; // B in T
    Lgm_CrossProduct( &u, &B, &uxB );

    /*
     * Now set up the system of equations to pass to solver
     */
    g = p->q/p->m/p->gamma;
    dydt[0] = u.x/p->gamma;
    dydt[1] = u.y/p->gamma;
    dydt[2] = u.z/p->gamma;
    dydt[3] = g * ( E.x + uxB.x );
    dydt[4] = g * ( E.y + uxB.y );
    dydt[5] = g * ( E.z + uxB.z );

    return( GSL_SUCCESS );

}

/*
 * 6D system for Full Orbit equations of motion. Jacobian
 *    The Jacobian is stored as J(i,j) = dfdy[ i*6 + j ]
 */
int J_full( double t, const double y[], double *dfdy, double dydt[], void *params ) {

    double  h, g, gamma;
    Lgm_Vector E, B, uxB, P, u;
    Lgm_ParticleInfo *p;
    p = (Lgm_ParticleInfo *)params;

    /*
     *  Electric field
     */
    E.x = E.y = E.z = 0.0;



    /*
     *  The position vector (in Re)
     */
    P.x = y[0]/1000.0/Re; P.y = y[1]/1000.0/Re; P.z = y[2]/1000.0/Re;

    /*
     *  The velocity vector (m/s)
     */
    u.x = y[3]; u.y = y[4]; u.z = y[5];
    gamma = sqrt( 1.0 + (u.x*u.x + u.y*u.y + u.z*u.z)/p->c2 );
    

    /*
     *  Compute uxB
     */
    p->mInfo->Bfield( &P, &B, p->mInfo );  // B in nT
    B.x *= 1e-9; B.y *= 1e-9; B.z *= 1e-9; // B in T
    Lgm_CrossProduct( &u, &B, &uxB );

    /*
     * Now set up the system of equations to pass to solver
     */
    g = p->q/p->m/gamma;
    dydt[0] = u.x/gamma;
    dydt[1] = u.y/gamma;
    dydt[2] = u.z/gamma;
    dydt[3] = g * ( E.x + uxB.x );
    dydt[4] = g * ( E.y + uxB.y );
    dydt[5] = g * ( E.z + uxB.z );


    h = 1.0/gamma;
    dfdy[ 0*6 + 0 ] = h;   dfdy[ 0*6 + 1 ] = 0.0; dfdy[ 0*6 + 2 ] = 0.0; dfdy[ 0*6 + 3 ] = 0.0;     dfdy[ 0*6 + 4 ] = 0.0;      dfdy[ 0*6 + 5 ] = 0.0;
    dfdy[ 1*6 + 0 ] = 0.0; dfdy[ 1*6 + 1 ] = h;   dfdy[ 1*6 + 2 ] = 0.0; dfdy[ 1*6 + 3 ] = 0.0;     dfdy[ 1*6 + 4 ] = 0.0;      dfdy[ 1*6 + 5 ] = 0.0;
    dfdy[ 2*6 + 0 ] = 0.0; dfdy[ 2*6 + 1 ] = 0.0; dfdy[ 2*6 + 2 ] = h;   dfdy[ 2*6 + 3 ] = 0.0;     dfdy[ 2*6 + 4 ] = 0.0;      dfdy[ 2*6 + 5 ] = 0.0;
    dfdy[ 3*6 + 0 ] = 0.0; dfdy[ 3*6 + 1 ] = 0.0; dfdy[ 3*6 + 2 ] = 0.0; dfdy[ 3*6 + 3 ] = 0.0;     dfdy[ 3*6 + 4 ] =  g*B.z;   dfdy[ 3*6 + 5 ] = -g*B.y;
    dfdy[ 4*6 + 0 ] = 0.0; dfdy[ 4*6 + 1 ] = 0.0; dfdy[ 4*6 + 2 ] = 0.0; dfdy[ 4*6 + 3 ] = -g*B.z;  dfdy[ 4*6 + 4 ] = 0.0;      dfdy[ 4*6 + 5 ] =  g*B.x;
    dfdy[ 5*6 + 0 ] = 0.0; dfdy[ 5*6 + 1 ] = 0.0; dfdy[ 5*6 + 2 ] = 0.0; dfdy[ 5*6 + 3 ] =  g*B.y;  dfdy[ 5*6 + 4 ] = -g*B.x;   dfdy[ 5*6 + 5 ] = 0.0;




    return( GSL_SUCCESS );

}





/*
 * 4D system for Guiding Center equations of motion.
 */
int f_gc( double t, const double y[], double dydt[], void *params ) {

    double      g, h;       // temporary vars
    double      b_dot_CurlE;
    Lgm_Vector  P, Q;

    Lgm_ParticleInfo *p;
    p = (Lgm_ParticleInfo *)params;


    p->p_par = y[3];


    /*
     * Determine B, E, Grad_B, Curl_b, dbdt for given time and position.
     * Lets try pure static dipole with no E-field to start.
     */
    p->X.x = y[0]/1000.0/Re; p->X.y = y[1]/1000.0/Re; p->X.z = y[2]/1000.0/Re; // convert position from m -> Re
    //p->mInfo->RBF_CompGradAndCurl = TRUE;
    p->mInfo->Bfield( &p->X, &p->B, p->mInfo );                                // B in nT
    //p->mInfo->RBF_CompGradAndCurl = FALSE;
    p->b = p->B; p->Bmag = Lgm_NormalizeVector( &p->b );                       // b is unit vector Bmag in nT
    //p->Grad_B = p->mInfo->RBF_Grad_B;
    //p->Curl_b = p->mInfo->RBF_Curl_B;

    Lgm_GradB( &p->X, &p->Grad_B, LGM_DERIV_SIX_POINT, 1e-3, p->mInfo );       // Grad_B in nT/Re
    g = 1e9*Re*1000.0; p->Grad_B.x /= g; p->Grad_B.y /= g; p->Grad_B.z /= g;   // Grad_B now in T/m

    Lgm_Curlb( &p->X, &p->Curl_b, LGM_DERIV_SIX_POINT, 1e-3, p->mInfo );       // Curl_b in 1/Re
    g = Re*1000.0; p->Curl_b.x /= g; p->Curl_b.y /= g; p->Curl_b.z /= g;       // Curl_b in 1/m

    //p->Curl_E.x /= g; p->Curl_E.y /= g; p->Curl_E.z /= g;                    // Curl_E in V/m^2
p->Curl_E.x = p->Curl_E.y = p->Curl_E.z = 0.0;
    p->Bmag *= 1e-9;                                                           // Bmag now in T
    p->B.x *= 1e-9; p->B.y *= 1e-9; p->B.z *= 1e-9;                            // B now in T

    /*
     * Compute dbhat/dt = (-CurlE + bhat(bhat dot CurlE) ) /Bmag
     */
p->Curl_E.x = p->Curl_E.y = p->Curl_E.z = 0.0;
b_dot_CurlE = Lgm_DotProduct( &p->b, &p->Curl_E );
Q.x = p->b.x * b_dot_CurlE; Q.y = p->b.y * b_dot_CurlE; Q.z = p->b.z * b_dot_CurlE;
Lgm_VecSub( &P, &Q, &p->Curl_E );
p->dbdt.x = P.x/p->Bmag; p->dbdt.y = P.y/p->Bmag; p->dbdt.z = P.z/p->Bmag; // dbdt in s^-1
//printf("dbdt = %g %g %g\n", p->dbdt.x, p->dbdt.y, p->dbdt.z);
    
    //p->E.x = p->mInfo->RBF_E.x;                                                 // Ex in V/m
    //p->E.y = p->mInfo->RBF_E.y;                                                 // Ey in V/m
    //p->E.z = p->mInfo->RBF_E.z;                                                 // Ez in V/m
p->E.x = 0.0; p->E.y = 0.0; p->E.z = 0.0;                                 // E in V/m
    
    
    /*
     * Compute The relativistic factor gamma and its gradient.  See Eq (6) of
     * Tao et al. [2007] for gamma. Its gradient is (see text immediately after
     * Eq (17) of Tao et al. [2007]);
     *        Grad gamma = (mu*Grad B)/(gamma m c^2) 
     */
    g = p->p_par/p->mc;
    p->gamma = sqrt( 1.0 + 2.0*p->mu*p->Bmag/p->mc2 + g*g );

//double hh = 2.0*p->mu*p->Bmag/p->mc2 + g*g;
//double KineticEnergy_n = p->m*hh*p->c2/(p->gamma + 1.0)*6.241509e+15;
//printf("mu = %g KineticEnergy_n=%g\n", p->mu, KineticEnergy_n);
    g = p->mu/(p->gamma*p->mc2);
    p->Grad_gamma.x = g*p->Grad_B.x;
    p->Grad_gamma.y = g*p->Grad_B.y;
    p->Grad_gamma.z = g*p->Grad_B.z;
//printf("Grad_gamma = %g %g %g\n", p->Grad_gamma.x, p->Grad_gamma.y, p->Grad_gamma.z);

    /*
     * Compute "Effective B Field", B*.  Eq (15) of Tao et al. [2007].  (Bs =
     * Bstar).
     */
    //g = p->c*p->p_par/p->q;
    g = p->p_par/p->q;
    p->Bs.x = p->B.x + g*p->Curl_b.x;
    p->Bs.y = p->B.y + g*p->Curl_b.y;
    p->Bs.z = p->B.z + g*p->Curl_b.z;
//printf("B = %g %g %g\n", p->B.x, p->B.y, p->B.z);


    /*
     * Compute "Effective E Field", E*.  Eq (16) of Tao et al. [2007].  (Es =
     * Estar).
     * Note that in Cary and Brizard, the formula for E* = E - 1/q (mc^2Grad_Gamma - p_par dbdt) (see Eqn 6.27)
     * I.e. the sign of dbdt term is opposite.
     */
    g = p->p_par/p->q;
    h = p->mc2/p->q;
    //p->Es.x = p->E.x - g*p->dbdt.x - h*p->Grad_gamma.x;
    //p->Es.y = p->E.y - g*p->dbdt.y - h*p->Grad_gamma.y;
    //p->Es.z = p->E.z - g*p->dbdt.z - h*p->Grad_gamma.z;
    p->Es.x = p->E.x + g*p->dbdt.x - h*p->Grad_gamma.x; //Cary and Brizard, [2009]
    p->Es.y = p->E.y + g*p->dbdt.y - h*p->Grad_gamma.y; //Cary and Brizard, [2009]
    p->Es.z = p->E.z + g*p->dbdt.z - h*p->Grad_gamma.z; //Cary and Brizard, [2009]

//printf("CurlE = %g %g %g   g*p->dbdt =  %g %g %g    h*p->Grad_gamma = %g %g %g\n", p->mInfo->RBF_Curl_E.x, p->mInfo->RBF_Curl_E.y, p->mInfo->RBF_Curl_E.z, g*p->dbdt.x, g*p->dbdt.y, g*p->dbdt.z, h*p->Grad_gamma.x, h*p->Grad_gamma.y, h*p->Grad_gamma.z );

    /*
     * Compute B*_par.  Eq (17) of Tao et al. [2007].
     */
    p->Bs_par = Lgm_DotProduct( &p->Bs, &p->b );


    /*
     * Compute E* x b (needed in Eqn (13) of Tao et al. [2007].
     */
    Lgm_CrossProduct( &p->Es, &p->b, &p->Es_cross_b );



    /****************************************************************************
     *  Relativistic Guiding Center equations.  Eqs (13-14) of Tao et al.       *
     *  [2007].  Xdot = velocity of guiding center, v_par_dot is rate of change *
     *  of parallel velocity.                                                   *
     *                                                                          *
     ****************************************************************************/
    g = p->p_par/(p->gamma*p->m);
    //h = p->c;
    h = 1.0; // this makes the formula consistent with the MKSA system
    p->X_dot.x = ( g*p->Bs.x + h*p->Es_cross_b.x ) / p->Bs_par;
    p->X_dot.y = ( g*p->Bs.y + h*p->Es_cross_b.y ) / p->Bs_par;
    p->X_dot.z = ( g*p->Bs.z + h*p->Es_cross_b.z ) / p->Bs_par;
    
    p->p_par_dot = p->q * Lgm_DotProduct( &p->Es, &p->Bs ) / p->Bs_par;
    //p->v_par_dot = p->q * Lgm_DotProduct( &p->Es, &p->Bs ) / p->Bs_par / p->m;


    /*
     * Now set up the system of equations to pass to solver
     */
    dydt[0] = p->X_dot.x;
    dydt[1] = p->X_dot.y;
    dydt[2] = p->X_dot.z;
    dydt[3] = p->p_par_dot;
//printf("Derivs: %g %g %g %g\n", dydt[0], dydt[1], dydt[2], dydt[3] );

    return( GSL_SUCCESS );

}


int Lgm_TraceParticle_GCA( Lgm_Vector *ParticlePosition, long int *nParticlePosition  ) {

    Lgm_ParticleInfo    p;
    double v_par, v_per, v_tot, p_par;
    long int Date;
    double   UTC;
    double Ra, MLT, Phi, Lat;
    Lgm_Vector u_sm, u0;
    int Flag;
    Lgm_Vector v1, v2, v3, Bvec, un;
    double AlphaEq, g, Blocal, Beq, Tb_approx, gamma, gamma_n, KineticEnergy_n;


    // Setup particle properties
    p.KE         = 10000.0/1000.0;    // Initial Kinetic Energy MeV
    AlphaEq      = 85.0; g = sin( AlphaEq*RadPerDeg );
    p.q          = 1.60217657e-19;  // Particle charge,  C
    p.m          = 9.10938291e-31;  // electron mass,  kg
    //p.m          = 1.67262178e-27;  // proton mass,  kg
    p.c          = 2.99792458e8;    // Speed of light, m/s
    p.c2         = p.c * p.c;       // c^2
    p.mc         = p.m * p.c;
    p.mc2        = p.m * p.c * p.c;



    // Date -- really should be taken from elsewhere..
    Date = 20020418;
    UTC  = 5.0 + 30.0/60.0 + 0.0/3600.0;

    p.mInfo = Lgm_InitMagInfo( );
    Lgm_Set_Coord_Transforms( Date, UTC, p.mInfo->c );
    p.mInfo->Kp = 2;
    Lgm_MagModelInfo_Set_MagModel( LGM_IGRF, LGM_EXTMODEL_T89c, p.mInfo );
//    Lgm_MagModelInfo_Set_MagModel( LGM_CDIP, LGM_EXTMODEL_NULL, p.mInfo );

    // Set initial position - make user selectable
    Ra  = 6.60000;
    MLT = -1.5;
    Phi = (MLT*15.0 - 180.0)*RadPerDeg; Lat = 0.0*RadPerDeg;
    u_sm.x = Ra*cos( Phi )*cos(Lat); u_sm.y = Ra*sin( Phi )*cos(Lat); u_sm.z = Ra*sin(Lat);
    Lgm_Convert_Coords( &u_sm, &u0, SM_TO_GSM, p.mInfo->c );



    Flag = Lgm_Trace( &u0, &v1, &v2, &v3, 120.0, 1e-7, 1e-7, p.mInfo );
    p.mInfo->Bfield( &u0, &Bvec, p.mInfo ); Blocal = Lgm_Magnitude( &Bvec );
    p.mInfo->Bfield( &v3, &Bvec, p.mInfo ); Beq    = Lgm_Magnitude( &Bvec );
    //printf("Blocal, Beq = %g %g\n", Blocal, Beq);


    p.B0         = Blocal;       // Initial |B| in nT.
    p.Alpha      = asin( sqrt( Blocal/Beq * g*g ) )*DegPerRad;
//v_tot        =  sqrt( Lgm_v2overc2( p.KE, LGM_Ep0 ) ) * p.c;
v_tot        =  sqrt( Lgm_v2overc2( p.KE, LGM_Ee0 ) ) * p.c;
    v_par        = v_tot * cos( p.Alpha*RadPerDeg );
    v_per        = v_tot * sin( p.Alpha*RadPerDeg );
    gamma        = 1.0/sqrt( 1.0 - v_tot*v_tot/p.c2 );

g = gamma*v_tot;
gamma_n = sqrt( 1.0 + g*g/p.c2 );
KineticEnergy_n = p.m*g*g/(gamma_n + 1.0)*6.241509e+15;
printf("KineticEnergy_n=%g\n", KineticEnergy_n);

//p.mu         = Lgm_Ek_to_Mu( p.KE, p.Alpha, p.B0, LGM_Ep0 ); // mu in MeV/G
p.mu         = Lgm_Ek_to_Mu( p.KE, p.Alpha, p.B0, LGM_Ee0 ); // mu in MeV/G
//printf("E = %g MeV, Alpha = %g Degrees, B = %g nT,  LGM_Ep0 = %g MeV  Mu = %g MeV/G\n", p.KE, p.Alpha, p.B0, LGM_Ep0, p.mu );
printf("E = %g MeV, Alpha = %g Degrees, B = %g nT,  LGM_Ee0 = %g MeV  Mu = %g MeV/G\n", p.KE, p.Alpha, p.B0, LGM_Ee0, p.mu );
    p.mu         *= 1.60217657e-9;  // mu in J/T


    p_par = gamma*p.m*v_par;



    int nSteps;
    double  y[4];
    nSteps = 50000;
    *nParticlePosition = 0;
    ParticlePosition[0] = u0; ++(*nParticlePosition);
    y[0] = u0.x*1000.0*Re;
    y[1] = u0.y*1000.0*Re;
    y[2] = u0.z*1000.0*Re;
    y[3] = p_par;



    // Setup ODE solver in GSL
    void *params;
    params = (void *)&p;
    const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rk8pd;
    //const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rk4;

    gsl_odeiv2_step    *s = gsl_odeiv2_step_alloc( T, 4 );
    gsl_odeiv2_control *c = gsl_odeiv2_control_standard_new( 1e-6, 1e-6, 1.0, 1.0 );
    gsl_odeiv2_evolve  *e = gsl_odeiv2_evolve_alloc( 4 );

    gsl_odeiv2_system sys = { f_gc, NULL, 4, params };

FILE *fp;
fp = fopen( "TestParticle.txt", "w");



    //See Hess or Anderson, 1966 (Annual Review of Nuc. Sci.)
//Tb_approx = 0.085* Ra/sqrt( Lgm_v2overc2( p.KE, LGM_Ep0 ) ) * (1.3 - 0.56*sin(p.Alpha*RadPerDeg) );
Tb_approx = 0.085* Ra/sqrt( Lgm_v2overc2( p.KE, LGM_Ee0 ) ) * (1.3 - 0.56*sin(p.Alpha*RadPerDeg) );
    printf("Approximate Bounce Period in Dipole: %g seconds\n", Tb_approx );


    g = gamma*v_tot;
    gamma_n = sqrt( 1.0 + g*g/p.c2 );
    KineticEnergy_n = p.m*g*g/(gamma_n + 1.0)*6.241509e+15;
    p.mInfo->Bfield( &ParticlePosition[0], &Bvec, p.mInfo ); Blocal = Lgm_Magnitude( &Bvec );
    fprintf( fp, "%g %g %g %g %g %g\n", 0.0, ParticlePosition[0].x, ParticlePosition[0].y, ParticlePosition[0].z, Blocal, KineticEnergy_n  );
printf( "START GCA: %g %g %g %g %g %g\n", 0.0, ParticlePosition[0].x, ParticlePosition[0].y, ParticlePosition[0].z, Blocal, KineticEnergy_n  );

    
    double h = 1e-6, t, ti, dT;
    long int i;
    int status;
    t  = 0.0;
    dT = 0.01;
dT = 0.001;
//    ti = 3600.0;
    for ( i = 1; i <= nSteps; i++ ) {

        ti = i*dT;
        while ( t < ti ) {

            if ( h > (ti-t) ) h = ti-t;

            status = gsl_odeiv2_evolve_apply( e, c, s, &sys, &t, ti, &h, y );

//printf( "h = %g   ti = %g  %g %g %g %g %g\n", h, ti, t, y[0]/1000.0/Re, y[1]/1000.0/Re, y[2]/1000.0/Re, y[3] );

            if ( status != GSL_SUCCESS ) {
                h = 1e-6;
                printf ("error, return value=%d\n", status);
                //exit(0);
                break;
            } else {
                dT = 0.001;
            }
        }

        g = gamma*v_tot;
        gamma_n = sqrt( 1.0 + g*g/p.c2 );
        KineticEnergy_n = p.m*g*g/(gamma_n + 1.0)*6.241509e+15;

        double xxx, yyy, zzz;
        xxx = y[0]/1000.0/Re;
        yyy = y[1]/1000.0/Re;
        zzz = y[2]/1000.0/Re;

        if ( ( xxx > 15.0) || ( xxx < -200.0) || ( fabs( yyy ) > 30.0 ) || ( fabs( zzz ) > 30.0 ) ) {
            printf("Particle lost: ParticlePosition[%d] = %g %g %g\n", i, ParticlePosition[i].x, ParticlePosition[i].y, ParticlePosition[i].z );
            break;
        }

        if ((i%1 == 0)&&( *nParticlePosition<50000)){
            ParticlePosition[*nParticlePosition].x = xxx;
            ParticlePosition[*nParticlePosition].y = yyy;
            ParticlePosition[*nParticlePosition].z = zzz;

            p.mInfo->Bfield( &ParticlePosition[*nParticlePosition], &Bvec, p.mInfo ); Blocal = Lgm_Magnitude( &Bvec );
            fprintf( fp, "%g %g %g %g %g %g\n", t, xxx, yyy, zzz, Blocal, KineticEnergy_n  );
            fflush( fp );

            ++(*nParticlePosition);
        }


    }

    //gsl_odeiv2_driver_free( d );
    gsl_odeiv2_evolve_free( e );
    gsl_odeiv2_control_free( c );
    gsl_odeiv2_step_free( s );

fclose(fp);


}



/*
 * Implements the relativistic "Boris push" to update the velocity.
 *  Inputs: dt, unm1o2, E, B, p
 *  Outputs: unp1o2, gamma_np1o2, gamma_n
 *  Note that gamma_n can be used to compute the "tiome centered kinetic energy". (See Birdsall and Langdon page 357.)
 */
int UpdateVelocityBoris( Lgm_Vector *unm1o2, Lgm_Vector *unp1o2, double *gamma_np1o2, double *gamma_n, Lgm_Vector *E, Lgm_Vector *B, double dt, Lgm_ParticleInfo *p ) {

    Lgm_Vector  um, up, uu, w;
    Lgm_Vector  s, t, e, bhat;
    double      g, Bmag;

    /*
     * Compute the "half acceleration"
     */
    g = 0.5*dt*p->q/p->m;
    e.x = g*E->x; e.y = g*E->y; e.z = g*E->z;
    

    /*
     * Compute u- by adding the half-accel, e to u_n-1/2
     */
    um.x = unm1o2->x + e.x; um.y = unm1o2->y + e.y; um.z = unm1o2->z + e.z;
    

    /*
     * Compute gamma_n = sqrt( 1 + (um/c)^2 )
     */
    g = um.x*um.x + um.y*um.y + um.z*um.z; // |um|^2
    *gamma_n = sqrt( 1.0 + g/p->c2 );


    /*
     * Compute the vectors needed to implement the rotation from v- to v+
     *     The angle of rotation is theta where tan(theta/2) where theta = q dt B/(m gamma-)
     */
    // normal approximation
    g = p->q*dt/(2.0* (*gamma_n) * p->m);
    t.x = g*B->x; t.y = g*B->y; t.z = g*B->z;

    // This is for Boris push with gyrophase correction
//    bhat = *B; Bmag = Lgm_NormalizeVector( &bhat );
//    g = tan(p->q*dt*Bmag/(2.0*(*gamma_n) * p->m));
//    t.x = g*bhat.x; t.y = g*bhat.y; t.z = g*bhat.z;


    g = 2.0/(1.0 + t.x*t.x + t.y*t.y + t.z*t.z);
    s.x = g*t.x; s.y = g*t.y; s.z = g*t.z;
    

    /*
     * Do Boris rotation
     */
    Lgm_CrossProduct( &um, &t, &w );   // w = um x t
    uu.x = um.x + w.x; uu.y = um.y + w.y; uu.z = um.z + w.z;  // u^prime = u- + u- x t

    Lgm_CrossProduct( &uu, &s, &w );   // w = u^prime x s
    up.x = um.x + w.x; up.y = um.y + w.y; up.z = um.z + w.z;  // u+ = u- + u^prime x s


    /*
     * Add final half acceleration to get u_n+1/2
     */
    unp1o2->x = up.x + e.x; unp1o2->y = up.y + e.y; unp1o2->z = up.z + e.z;


    /*
     *  Finally compute and save the gamma_n+1/2 value (needed for position advancing)
     */
    g = Lgm_Magnitude( unp1o2 )/p->c;
    *gamma_np1o2 = sqrt( 1.0 + g*g );
    //printf("unp1o2 = %g %g %g,    p->c = %g    g, gamma_np1o2 = %g %g\n", unp1o2->x, unp1o2->y, unp1o2->z, p->c, g, *gamma_np1o2);
    

}



/*
 * Simple euler scheme to update velocity.
 * Needed to get Boris leapfrog scheme going..
 */
int UpdateVelocityEuler( Lgm_Vector *u0, Lgm_Vector *u1, Lgm_Vector *E, Lgm_Vector *B, double dt, Lgm_ParticleInfo *p ) {


    Lgm_Vector  G, u_cross_B;
    double      g, gamma0;

    /*
     *  Compute u x B
     */
    Lgm_CrossProduct( u0, B, &u_cross_B );

    /*
     *  Compute F/m
     */
    g = Lgm_Magnitude( u0 );
    gamma0 = sqrt( 1.0 + g*g/p->c2 );
    g = p->q/p->m;
    G.x = g * (E->x + u_cross_B.x/gamma0); G.y = g * (E->y + u_cross_B.y/gamma0); G.z = g * (E->z + u_cross_B.z/gamma0);

    /*
     * Update u0 -> u1
     */
    u1->x = u0->x + dt * G.x; u1->y = u0->y + dt*G.y; u1->z = u0->z + dt*G.z;

}

int UpdatePositionBoris( Lgm_Vector *wn, Lgm_Vector *wnp1, Lgm_Vector *unp1o2, double gamma_np1o2, double dt, Lgm_ParticleInfo *p ) {

    wnp1->x = wn->x + unp1o2->x * dt / gamma_np1o2;
    wnp1->y = wn->y + unp1o2->y * dt / gamma_np1o2;
    wnp1->z = wn->z + unp1o2->z * dt / gamma_np1o2;
//printf("gamma_np1o2 = %g\n", gamma_np1o2);

}



int Lgm_TraceParticle_Boris( Lgm_Vector *ParticlePosition, long int *nParticlePosition  ) {

    Lgm_ParticleInfo    p;
    double v_tot, v_par, v_per;
    long int Date;
    double   UTC;
    double Ra, MLT, Phi, Lat;
    Lgm_Vector  Rc_hat_pqb, Rc_hat_gsm, v_per_hat_pqb, v_per_hat_gsm;
    Lgm_Vector w_sm, w0, w_off, ww;
    int Flag;
    Lgm_Vector u0, u1;
    Lgm_Vector v1, v2, v3, Bvec, bhat, B, E;
    double AlphaEq, g, Blocal, Beq, Tb_approx, t, dt, Rc, phi;

    Lgm_Vector  unm1o2, unp1o2, wn, wnp1, un;
    double      gamma_np1o2, gamma, gamma_n, gamma_n2, KineticEnergy_n;
    int         i;


    // Setup particle properties
    p.KE         = 5000.0/1000.0;    // Initial Kinetic Energy MeV
    AlphaEq      = 85.0; g = sin( AlphaEq*RadPerDeg );
    p.q          = 1.60217657e-19;  // Particle charge,  C
    //p.m          = 9.10938291e-31;  // electron mass,  kg
    p.m          = 1.67262178e-27;  // proton mass,  kg
    p.c          = 2.99792458e8;    // Speed of light, m/s
    p.c2         = p.c * p.c;       // c^2
    p.mc         = p.m * p.c;
    p.mc2        = p.m * p.c * p.c;



    // Date -- really should be taken from elsewhere..
    Date = 20020418;
    UTC  = 5.0 + 30.0/60.0 + 0.0/3600.0;

    p.mInfo = Lgm_InitMagInfo( );
    Lgm_Set_Coord_Transforms( Date, UTC, p.mInfo->c );
    p.mInfo->Kp = 2;
    Lgm_MagModelInfo_Set_MagModel( LGM_IGRF, LGM_EXTMODEL_T89c, p.mInfo );
//    Lgm_MagModelInfo_Set_MagModel( LGM_CDIP, LGM_EXTMODEL_NULL, p.mInfo );

    // Set initial position - make user selectable
    Ra  = 6.60000;
    MLT = -1.5;
    Phi = (MLT*15.0 - 180.0)*RadPerDeg; Lat = 0.0*RadPerDeg;
    w_sm.x = Ra*cos( Phi )*cos(Lat); w_sm.y = Ra*sin( Phi )*cos(Lat); w_sm.z = Ra*sin(Lat);
    Lgm_Convert_Coords( &w_sm, &w0, SM_TO_GSM, p.mInfo->c );



    Flag = Lgm_Trace( &w0, &v1, &v2, &v3, 120.0, 1e-7, 1e-7, p.mInfo );
    p.mInfo->Bfield( &v3, &Bvec, p.mInfo ); Beq = Lgm_Magnitude( &Bvec );
    p.mInfo->Bfield( &w0, &Bvec, p.mInfo ); B = Bvec; bhat = Bvec; Blocal = Lgm_NormalizeVector( &bhat );
    //printf("Blocal, Beq = %g %g\n", Blocal, Beq);
    B.x *= 1e-9; B.y *= 1e-9; B.z *= 1e-9;


    p.B0         = Blocal;       // Initial |B| in nT.
    p.Alpha      = asin( sqrt( Blocal/Beq * g*g ) )*DegPerRad;
    g = sin( p.Alpha * RadPerDeg );
    double B_mirror = Blocal/(g*g);
//v_tot        = sqrt( Lgm_v2overc2( p.KE, LGM_Ep0 ) ) * p.c;
v_tot        = sqrt( Lgm_v2overc2( p.KE, LGM_Ee0 ) ) * p.c;
    v_par        = v_tot * cos( p.Alpha*RadPerDeg );
    v_per        = v_tot * sin( p.Alpha*RadPerDeg );

//p.mu         = Lgm_Ek_to_Mu( p.KE, p.Alpha, p.B0, LGM_Ep0 ); // mu in MeV/G
p.mu         = Lgm_Ek_to_Mu( p.KE, p.Alpha, p.B0, LGM_Ee0 ); // mu in MeV/G
//printf("E = %g MeV, Alpha = %g Degrees, B = %g nT,  LGM_Ep0 = %g MeV  Mu = %g MeV/G     B_mirror = %g\n", p.KE, p.Alpha, p.B0, LGM_Ep0, p.mu, B_mirror );
printf("E = %g MeV, Alpha = %g Degrees, B = %g nT,  LGM_Ee0 = %g MeV  Mu = %g MeV/G     B_mirror = %g\n", p.KE, p.Alpha, p.B0, LGM_Ee0, p.mu, B_mirror );
    p.mu         *= 1.60217657e-9;  // mu in J/T



    FILE *fp;
    fp = fopen( "TestParticle2.txt", "w"); //Boris


    double  Tc;
    //See Hess or Anderson, 1966 (Annual Review of Nuc. Sci.)
    Tc = 2.0*M_PI * p.m/fabs(p.q * p.B0*1e-9);
    printf("Approximate Cycoltron Period in Dipole: %g seconds\n", Tc );

    //printf("w0 = %g %g %g\n", w0.x, w0.y, w0.z);






// FROM HERE
// We are shifting from GCA to Full orbit positions
// Initial GCA position is w0
    /*
     * Our original pointm is the GC point. Here, we shift to full orbit position.
     */
    Lgm_Set_GSM_TO_PQB( &w0, p.mInfo ); // This will set up the PQB coord system for this location

    /*
     *  Need to pick a phase angle in the perp plane to get where to offset particle
     */
    phi = 90.0; Rc_hat_pqb.x = cos(phi*RadPerDeg); Rc_hat_pqb.y = sin(phi*RadPerDeg); Rc_hat_pqb.z = 0.0;
    Lgm_MatTimesVec( p.mInfo->Apqb_to_gsm, &Rc_hat_pqb, &Rc_hat_gsm ); // dont need to shift system for velocity
    //printf("Rc_hat_pqb = %g %g %g    Rc_hat_gsm = %g %g %g\n", Rc_hat_pqb.x, Rc_hat_pqb.y, Rc_hat_pqb.z,   Rc_hat_gsm.x, Rc_hat_gsm.y, Rc_hat_gsm.z);


    v_per_hat_pqb.x = -sin(phi*RadPerDeg); v_per_hat_pqb.y = cos(phi*RadPerDeg); v_per_hat_pqb.z = 0.0;;
    Lgm_MatTimesVec( p.mInfo->Apqb_to_gsm, &v_per_hat_pqb, &v_per_hat_gsm ); // dont need to shift system for velocity
    //printf("v_per_hat_pqb = %g %g %g    v_per_hat_gsm = %g %g %g\n", v_per_hat_pqb.x, v_per_hat_pqb.y, v_per_hat_pqb.z,   v_per_hat_gsm.x, v_per_hat_gsm.y, v_per_hat_gsm.z);


    /*
     * We need to offset the starting position of the particle based on its gyroradius.
     */
    gamma = 1.0/sqrt( 1.0 - v_tot*v_tot/p.c2);
    Rc = fabs(gamma*p.m*v_per/(p.q*Blocal*1e-9)); // particle gyro-radius in meters
    //printf("Blocal = %g   Rc = %g\n", Blocal, Rc);


    w_off.x = w0.x*1000.0*Re + Rc * Rc_hat_gsm.x; // m
    w_off.y = w0.y*1000.0*Re + Rc * Rc_hat_gsm.y; // m
    w_off.z = w0.z*1000.0*Re + Rc * Rc_hat_gsm.z; // m
    //printf("w_off = %g %g %g\n", w_off.x/1000.0/Re, w_off.y/1000.0/Re, w_off.z/1000.0/Re);
    

    /*
     * Recompute local b-field related quantities.
     * 
     * We are shifting where the particle is (from the GC to its full orbit helix).
     * A problem with this is that the b-hast is slightly different at this ne location.
     * If we just use the previous b-hat and v-perp directions, thinks will be a bit off.
     * (You have to be careful here or KE wont be conserved.)
     */
    ww.x = w_off.x/1000.0/Re; ww.y = w_off.y/1000.0/Re; ww.z = w_off.z/1000.0/Re;
    
    // New (slightly different) b-hat.
    p.mInfo->Bfield( &ww, &bhat, p.mInfo ); Blocal = Lgm_NormalizeVector( &bhat ); 

    // New (slightly different) v_per_hat_pqb and v_per_hat_gsm vector.
    Lgm_Set_GSM_TO_PQB( &ww, p.mInfo ); 
    Lgm_MatTimesVec( p.mInfo->Apqb_to_gsm, &v_per_hat_pqb, &v_per_hat_gsm ); 

    // Recompute local pitch angle (because |Blocal| is slightly different).
    p.Alpha = asin( sqrt( Blocal/B_mirror ) )*DegPerRad;
    printf("Blocal, p.Alpha = %g, %g\n", Blocal, p.Alpha);

    // Reset u0 = w_off
    w0 = w_off;
// TO HERE
// We are shifting from GCA to Full orbit positions
















    // partition of v_tot between v_par and v_per.
    v_par        = v_tot * cos( p.Alpha*RadPerDeg );
    v_per        = v_tot * sin( p.Alpha*RadPerDeg );


    // initial velocity vector u = gamma*v
    //printf("v_par, v_per = %g %g    Bvec = %g %g %g    v_per_hat_gsm = %g %g %g\n", v_par, v_per, Bvec.x, Bvec.y, Bvec.z, v_per_hat_gsm.x, v_per_hat_gsm.y, v_per_hat_gsm.z);
    u0.x = gamma*(v_par * bhat.x + v_per * v_per_hat_gsm.x);
    u0.y = gamma*(v_par * bhat.y + v_per * v_per_hat_gsm.y);
    u0.z = gamma*(v_par * bhat.z + v_per * v_per_hat_gsm.z);
    //printf("bhat dot v_per_hat_gsm  = %g    Lgm_Magnitude( &v_per_hat_gsm ) = %g\n", Lgm_DotProduct( &bhat, &v_per_hat_gsm ), Lgm_Magnitude( &v_per_hat_gsm ) );


Lgm_Vector vvhat;
vvhat = u0;
Lgm_NormalizeVector( &vvhat );
double NewPA;
NewPA = acos( Lgm_DotProduct( &bhat, &vvhat ) ) * DegPerRad;


printf("Lgm_Magnitude( &bhat ) = %g\n", Lgm_Magnitude( &bhat ));
printf("Lgm_Magnitude( &vvhat ) = %g\n", Lgm_Magnitude( &vvhat ));
printf("bhat = %g %g %g\n", bhat.x, bhat.y, bhat.z);
printf("vvhat = %g %g %g\n", vvhat.x, vvhat.y, vvhat.z);
printf("Pitch ANGLE:  %g\n",    NewPA);
//exit(0);

    /*
     *  Backup the velocity by 1/2 time step. And then do Boris stepping
     */
    wn = w0;
    E.x = E.y = E.z = 0.0;
    t = 0.0; dt = 0.00001;
    UpdateVelocityEuler( &u0, &unm1o2, &E, &B, -dt/2.0, &p );
//printf("u0 = %g %g %g\n", u0.x, u0.y, u0.z);
//printf("unm1o2 = %g %g %g\n", unm1o2.x, unm1o2.y, unm1o2.z);
    *nParticlePosition = 0;
    ParticlePosition[0].x = wn.x/1000.0/Re;
    ParticlePosition[0].y = wn.y/1000.0/Re;
    ParticlePosition[0].z = wn.z/1000.0/Re;
    ++(*nParticlePosition);

g = Lgm_Magnitude( &u0 );
gamma_n2 = 1.0 + g*g/p.c2;
gamma_n = sqrt( gamma_n2 );
KineticEnergy_n = p.m*g*g/(gamma_n + 1.0)*6.241509e+15;
printf("|u0| = %.12lf\n", g);
printf("|u0|^2 = %.12lf\n", g*g);
printf("c^2 = %.12lf\n", p.c2);
printf("m = %g\n", p.m);
printf("gamma_n = %.12lf\n", gamma_n);
printf("KineticEnergy_n = %.12lf\n", KineticEnergy_n);
p.mInfo->Bfield( &ParticlePosition[0], &Bvec, p.mInfo ); Blocal = Lgm_Magnitude( &Bvec );
fprintf( fp, "%g %g %g %g %g %g\n", 0.0, ParticlePosition[0].x, ParticlePosition[0].y, ParticlePosition[0].z, Blocal, KineticEnergy_n ); fflush( fp );
printf( "START Boris: %g %g %g %g %g %g\n", 0.0, ParticlePosition[0].x, ParticlePosition[0].y, ParticlePosition[0].z, Blocal, KineticEnergy_n  );


    for (i=1; i<1320000; i++){

        
        ww.x = wn.x/1000.0/Re; ww.y = wn.y/1000.0/Re; ww.z = wn.z/1000.0/Re; // Re
        p.mInfo->Bfield( &ww, &B, p.mInfo ); // B in nT
        B.x *= 1e-9; B.y *= 1e-9; B.z *= 1e-9; // B in T

        // Update veolicty
        UpdateVelocityBoris( &unm1o2, &unp1o2, &gamma_np1o2, &gamma_n, &E, &B, dt, &p ); 

        // update position
//printf("gamma_n, gamma_np1o2 = %g %g\n", gamma_n, gamma_np1o2);
        UpdatePositionBoris( &wn, &wnp1, &unp1o2, gamma_np1o2, dt, &p );
        t += dt;

        KineticEnergy_n = p.m*g*g/(gamma_n + 1.0)*6.241509e+15;

        if ((i%1 == 0)&&( *nParticlePosition<50000)){
            ParticlePosition[*nParticlePosition].x = wnp1.x/1000.0/Re;
            ParticlePosition[*nParticlePosition].y = wnp1.y/1000.0/Re;
            ParticlePosition[*nParticlePosition].z = wnp1.z/1000.0/Re;

            p.mInfo->Bfield( &ParticlePosition[*nParticlePosition], &Bvec, p.mInfo ); Blocal = Lgm_Magnitude( &Bvec );
            fprintf( fp, "%g %g %g %g %g %g\n", t, ParticlePosition[*nParticlePosition].x, ParticlePosition[*nParticlePosition].y, ParticlePosition[*nParticlePosition].z, Blocal, KineticEnergy_n );
            fflush( fp );



            ++(*nParticlePosition);
        }
        
        unm1o2 = unp1o2; // set previous vel. to current vel.
        wn     = wnp1;   // set previous pos. to current pos.

    }

}




int Lgm_TraceParticle_Full( Lgm_Vector *ParticlePosition, long int *nParticlePosition,
Lgm_Vector *ParticlePosition2, long int *nParticlePosition2  ) {

    Lgm_ParticleInfo    p;
    double v_tot, v_par, v_per;
    long int Date;
    double   UTC;
    double Ra, MLT, Phi, Lat;
    Lgm_Vector u_sm, u0;
    int Flag;
    Lgm_Vector v1, v2, v3, Bvec, bhat, B, E;
    Lgm_Vector  Rc_hat_pqb, Rc_hat_gsm, v_per_hat_pqb, v_per_hat_gsm;
    Lgm_Vector w_sm, w0, w_off, ww, un;
    double AlphaEq, g, Blocal, Beq, Tb_approx, t, dt, Rc, phi, gamma;

    double      gamma_n, KineticEnergy_n;


    // Setup particle properties
    p.KE         = 5000.0/1000.0;    // Initial Kinetic Energy MeV
    AlphaEq      = 85.0; g = sin( AlphaEq*RadPerDeg );
    p.q          = 1.60217657e-19;  // Particle charge,  C
    //p.m          = 9.10938291e-31;  // electron mass,  kg
    p.m          = 1.67262178e-27;  // proton mass,  kg
    p.c          = 2.99792458e8;    // Speed of light, m/s
    p.c2         = p.c * p.c;       // c^2
    p.mc         = p.m * p.c;
    p.mc2        = p.m * p.c * p.c;



    // Date -- really should be taken from elsewhere..
    Date = 20020418;
    UTC  = 5.0 + 30.0/60.0 + 0.0/3600.0;

    p.mInfo = Lgm_InitMagInfo( );
    Lgm_Set_Coord_Transforms( Date, UTC, p.mInfo->c );
    p.mInfo->Kp = 2;
    Lgm_MagModelInfo_Set_MagModel( LGM_IGRF, LGM_EXTMODEL_T89c, p.mInfo );
//    Lgm_MagModelInfo_Set_MagModel( LGM_CDIP, LGM_EXTMODEL_NULL, p.mInfo );

    // Set initial position - make user selectable
    Ra  = 6.60000;
    MLT = -1.5;
    Phi = (MLT*15.0 - 180.0)*RadPerDeg; Lat = 0.0*RadPerDeg;
    u_sm.x = Ra*cos( Phi )*cos(Lat); u_sm.y = Ra*sin( Phi )*cos(Lat); u_sm.z = Ra*sin(Lat);
    Lgm_Convert_Coords( &u_sm, &w0, SM_TO_GSM, p.mInfo->c );



    Flag = Lgm_Trace( &w0, &v1, &v2, &v3, 120.0, 1e-7, 1e-7, p.mInfo );
    p.mInfo->Bfield( &v3, &Bvec, p.mInfo ); Beq = Lgm_Magnitude( &Bvec );
    p.mInfo->Bfield( &w0, &Bvec, p.mInfo ); B = Bvec; bhat = Bvec; Blocal = Lgm_NormalizeVector( &bhat );
    //printf("Blocal, Beq = %g %g\n", Blocal, Beq);
    B.x *= 1e-9; B.y *= 1e-9; B.z *= 1e-9;


    p.B0         = Blocal;       // Initial |B| in nT.
    p.Alpha      = asin( sqrt( Blocal/Beq * g*g ) )*DegPerRad;
    g = sin( p.Alpha * RadPerDeg );
    double B_mirror = Blocal/(g*g);
//v_tot        = sqrt( Lgm_v2overc2( p.KE, LGM_Ep0 ) ) * p.c;
v_tot        = sqrt( Lgm_v2overc2( p.KE, LGM_Ee0 ) ) * p.c;
v_par        = v_tot * cos( p.Alpha*RadPerDeg );
v_per        = v_tot * sin( p.Alpha*RadPerDeg );

//p.mu         = Lgm_Ek_to_Mu( p.KE, p.Alpha, p.B0, LGM_Ep0 ); // mu in MeV/G
p.mu         = Lgm_Ek_to_Mu( p.KE, p.Alpha, p.B0, LGM_Ee0 ); // mu in MeV/G
//printf("E = %g MeV, Alpha = %g Degrees, B = %g nT,  LGM_Ep0 = %g MeV  Mu = %g MeV/G     B_mirror = %g\n", p.KE, p.Alpha, p.B0, LGM_Ep0, p.mu, B_mirror );
printf("E = %g MeV, Alpha = %g Degrees, B = %g nT,  LGM_Ee0 = %g MeV  Mu = %g MeV/G     B_mirror = %g\n", p.KE, p.Alpha, p.B0, LGM_Ee0, p.mu, B_mirror );
    p.mu         *= 1.60217657e-9;  // mu in J/T


gamma_n = 1.0/sqrt(1.0-v_tot*v_tot/p.c2);
KineticEnergy_n = (gamma_n - 1.0)*p.m*p.c2*6.241509e+15;
printf("gamma_n, KE = %lf,   %lf\n", gamma_n, KineticEnergy_n);

double uuu;
uuu = v_tot*gamma_n;
KineticEnergy_n = p.m*uuu*uuu/(1.0+gamma_n)*6.241509e+15;
printf("uuu = %g   gamma_n = %g   KE = %lf\n", uuu, gamma_n, KineticEnergy_n);



//exit(0);




    double  Tc;
    //printf("w0 = %g %g %g\n", w0.x, w0.y, w0.z);
    Lgm_Set_GSM_TO_PQB( &w0, p.mInfo ); // This will set up the PQB coord system for this location

    /*
     *  Need to pick a phase angle in the perp plane to get where to offset particle
     */
    phi = 90.0; Rc_hat_pqb.x = cos(phi*RadPerDeg); Rc_hat_pqb.y = sin(phi*RadPerDeg); Rc_hat_pqb.z = 0.0;
    Lgm_MatTimesVec( p.mInfo->Apqb_to_gsm, &Rc_hat_pqb, &Rc_hat_gsm ); // dont need to shift system for velocity
    //printf("Rc_hat_pqb = %g %g %g    Rc_hat_gsm = %g %g %g\n", Rc_hat_pqb.x, Rc_hat_pqb.y, Rc_hat_pqb.z,   Rc_hat_gsm.x, Rc_hat_gsm.y, Rc_hat_gsm.z);


    v_per_hat_pqb.x = -sin(phi*RadPerDeg); v_per_hat_pqb.y = cos(phi*RadPerDeg); v_per_hat_pqb.z = 0.0;;
    Lgm_MatTimesVec( p.mInfo->Apqb_to_gsm, &v_per_hat_pqb, &v_per_hat_gsm ); // dont need to shift system for velocity
    //printf("v_per_hat_pqb = %g %g %g    v_per_hat_gsm = %g %g %g\n", v_per_hat_pqb.x, v_per_hat_pqb.y, v_per_hat_pqb.z,   v_per_hat_gsm.x, v_per_hat_gsm.y, v_per_hat_gsm.z);

printf("bhat dot v_per_hat_gsm = %g\n", Lgm_DotProduct( &bhat, &v_per_hat_gsm) );


    //See Hess or Anderson, 1966 (Annual Review of Nuc. Sci.)
    Tc = 2.0*M_PI * p.m/fabs(p.q * p.B0*1e-9);
    //printf("Approximate Cycoltron Period in Dipole: %g seconds\n", Tc );

    /*
     * We need to offset the starting position of the particle based on its gyroradius.
     */
printf("A. gamma = %lf\n", gamma_n);
    gamma = 1.0/sqrt( 1.0 - v_tot*v_tot/p.c2);
printf("B. gamma = %lf\n", gamma);
    Rc = fabs(gamma*p.m*v_per/(p.q*Blocal*1e-9)); // particle gyro-radius in meters
    //printf("Blocal = %g v_per = %g   p.q = %g   Rc = %g\n", Blocal, v_per, p.q, Rc);


    w_off.x = w0.x*1000.0*Re + Rc * Rc_hat_gsm.x; // m
    w_off.y = w0.y*1000.0*Re + Rc * Rc_hat_gsm.y; // m
    w_off.z = w0.z*1000.0*Re + Rc * Rc_hat_gsm.z; // m
    //printf("w_off = %g %g %g\n", w_off.x/1000.0/Re, w_off.y/1000.0/Re, w_off.z/1000.0/Re);
    






    /*
     * Recompute local b-field related quantities.
     * 
     * We are shifting where the particle is (from the GC to its full orbit helix).
     * A problem with this is that the b-hast is slightly different at this ne location.
     * If we just use the previous b-hat and v-perp directions, thinks will be a bit off.
     * (You have to be careful here or KE wont be conserved.)
     */
    ww.x = w_off.x/1000.0/Re; ww.y = w_off.y/1000.0/Re; ww.z = w_off.z/1000.0/Re;
    
    // New (slightly different) b-hat.
    p.mInfo->Bfield( &ww, &bhat, p.mInfo ); Blocal = Lgm_NormalizeVector( &bhat ); 

    // New (slightly different) v_per_hat_pqb and v_per_hat_gsm vector.
    Lgm_Set_GSM_TO_PQB( &ww, p.mInfo ); 
    Lgm_MatTimesVec( p.mInfo->Apqb_to_gsm, &v_per_hat_pqb, &v_per_hat_gsm ); 

    // Recompute local pitch angle (because |Blocal| is slightly different).
    p.Alpha = asin( sqrt( Blocal/B_mirror ) )*DegPerRad;
    printf("Blocal, p.Alpha = %g, %g\n", Blocal, p.Alpha);

    // Recompute (slightly different partitionin of v_tot between v_par and v_per.
    v_par        = v_tot * cos( p.Alpha*RadPerDeg );
    v_per        = v_tot * sin( p.Alpha*RadPerDeg );
    //printf("A. v_tot = %g\n", v_tot ) ;
    //printf("B. v_tot = %g\n", sqrt( v_par*v_par + v_per*v_per ) );

    

    



    // initial velocity vector u = gamma*v
    //printf("v_par, v_per = %g %g    Bvec = %g %g %g    v_per_hat_gsm = %g %g %g\n", v_par, v_per, Bvec.x, Bvec.y, Bvec.z, v_per_hat_gsm.x, v_per_hat_gsm.y, v_per_hat_gsm.z);
    u0.x = gamma*(v_par * bhat.x + v_per * v_per_hat_gsm.x);
    u0.y = gamma*(v_par * bhat.y + v_per * v_per_hat_gsm.y);
    u0.z = gamma*(v_par * bhat.z + v_per * v_per_hat_gsm.z);
    //printf("bhat dot v_per_hat_gsm  = %g    Lgm_Magnitude( &v_per_hat_gsm ) = %g\n", Lgm_DotProduct( &bhat, &v_per_hat_gsm ), Lgm_Magnitude( &v_per_hat_gsm ) );



printf("C. v_tot, |u0|/gamma = %g %g\n", v_tot, Lgm_Magnitude( &u0 )/gamma );
Lgm_Vector uuu0;
uuu0.x = v_par * bhat.x + v_per * v_per_hat_gsm.x;
uuu0.y = v_par * bhat.y + v_per * v_per_hat_gsm.y;
uuu0.z = v_par * bhat.z + v_per * v_per_hat_gsm.z;
printf("D. v_tot, |uuu0| = %g %g\n", v_tot, Lgm_Magnitude( &uuu0 ) );
printf("|bhat| = %g\n", Lgm_Magnitude( &bhat ) );
printf("|v_per_hat_gsm| = %g\n", Lgm_Magnitude( &v_per_hat_gsm ) );
printf("E. v_tot = %g\n", sqrt( v_par*v_par +v_per*v_per ) );





    double  y[6];
    *nParticlePosition = 0;
    ParticlePosition[0].x = w_off.x/1000.0/Re; 
    ParticlePosition[0].y = w_off.y/1000.0/Re; 
    ParticlePosition[0].z = w_off.z/1000.0/Re; 
    ++(*nParticlePosition);
    y[0] = w_off.x;
    y[1] = w_off.y;
    y[2] = w_off.z;
    y[3] = u0.x;
    y[4] = u0.y;
    y[5] = u0.z;


    InitParticleIntegrators( &p );




    g = Lgm_Magnitude( &u0 );
    gamma_n = sqrt( 1.0 + g*g/p.c2 );
    printf("g = %g   gamma, gamma_n = %lf %lf\n", g,  gamma, gamma_n);
    KineticEnergy_n = p.m*g*g/(gamma_n + 1.0)*6.241509e+15;
    printf("g = %g , gamma_n = %lf , KineticEnergy_n = %g\n", g, gamma_n, KineticEnergy_n);





    // save initial position
    //p.mInfo->Bfield( &ParticlePosition[0], &Bvec, p.mInfo ); Blocal = Lgm_Magnitude( &Bvec );
    //fprintf( fp, "%g %g %g %g %g %g\n", 0.0, ParticlePosition[0].x, ParticlePosition[0].y, ParticlePosition[0].z, Blocal, KineticEnergy_n  );


    
    int nSteps;
    double ti, dT;
    int i, status;
    nSteps = 132000;
    t  = 0.0;
    dT = 0.00001;
    for ( i = 1; i <= nSteps; i++ ) {
        ti = i*dT;
        status = gsl_odeiv2_driver_apply( p.Driver_Full, &t, ti, y );
        if ( status != GSL_SUCCESS ) {
            printf ("error, return value=%d\n", status);
            break;
        } 
        
       

        double xxx, yyy, zzz;
        xxx = y[0]/1000.0/Re;
        yyy = y[1]/1000.0/Re;
        zzz = y[2]/1000.0/Re;
        if ( (xxx > 15.0) || (xxx < -200.0) || ( fabs( yyy ) > 30.0 ) || ( fabs( zzz ) > 30.0 ) ) {
            printf("Particle lost: ParticlePosition[%d] = %g %g %g\n", i, xxx, yyy, zzz );
            break;
        }

        un.x = y[3]; un.y = y[4]; un.z = y[5];
        g = Lgm_Magnitude( &un );
        gamma_n = sqrt( 1.0 + g*g/p.c2 );
        KineticEnergy_n = p.m*g*g/(gamma_n + 1.0)*6.241509e+15;



        if ((i%1 == 0)&&( *nParticlePosition<50000)){

            ParticlePosition[ *nParticlePosition ].x = xxx;
            ParticlePosition[ *nParticlePosition ].y = yyy;
            ParticlePosition[ *nParticlePosition ].z = zzz;

            p.mInfo->Bfield( &ParticlePosition[ *nParticlePosition ], &Bvec, p.mInfo ); Blocal = Lgm_Magnitude( &Bvec );
            //fprintf( fp, "%g %g %g %g %g %g\n", t, xxx, yyy, zzz, Blocal, KineticEnergy_n  );
            //fflush( fp );

            ++(*nParticlePosition);
        }

    }


}



/*
 * This sets up GSL ODE drivers for the Full Orbit and GCA 
 * equations.
 */

int InitParticleIntegrators( Lgm_ParticleInfo *p ) {

    void *params; params  = (void *)p;

    /*
     * Setup ODE solver for 6D Full orbit equations in GSL
     */
    p->Sys_Full.function  = f_full;
    p->Sys_Full.jacobian  = J_full;
    p->Sys_Full.dimension = 6;
    p->Sys_Full.params    = params;
    p->Driver_Full        = gsl_odeiv2_driver_alloc_y_new( &(p->Sys_Full), gsl_odeiv2_step_rk8pd, 1e-6, 1e-6, 0.0 );


    /*
     * Setup ODE solver for 4D Guiding Center Approximation (GCA) equations in GSL
     */
    p->Sys_GCA.function  = f_gc;
    p->Sys_GCA.jacobian  = NULL; // dont use methods that require jacobian (not easy to compute)
    p->Sys_GCA.dimension = 4;
    p->Sys_GCA.params    = params;
    p->Driver_GCA        = gsl_odeiv2_driver_alloc_y_new( &(p->Sys_GCA), gsl_odeiv2_step_rk8pd, 1e-6, 1e-6, 0.0 );


}





/*
To implement domain decomp on the time slices, we will need to push a particle
from ti -> ti+1, where these could be fairly large steps.  I.e., they would be as
large as you can fit into the memory of a single core.  If you have very high
spatial resolution, then the full files are pretty large.  Further if we have
high time resolution, then we cant fit very many time slices in memory at once.
Thus we need to decompose the time over the full run.  To do this, we could
allocate banks of processors to handle only certain time slices...  Then we can
use the appropriate processors to push over the associated times...  For
example, if we had NxM array of processors, we could reserve columms for
feeding particles in at the "top" The rows would be time intervals. I.e. row 1
would be populated with processors that are all loaded with times slices to go
from t0->t1. Etc...  When a particle finishes in row 1 it moves on the to next
row unmtil done....
*/

/*
 *  Push a particle by a timestep dt (in seconds).
 *
 *  The particle characteristics are all contained in the Lgm_ParticleInfo
 *  structure.
 *  
 *  The positions and velocities are also contained in the Lgm_ParticleInfo
 *  structure.
 */
int Lgm_PushParticle( double dT, Lgm_ParticleInfo *p ) {

    double  y[6];
    int     status;    // error status from gsl routines
    double  t  = 0.0;


    if ( p->CurrentIntegrator == FULL_ORBIT ) {

        y[0] = p->x.x; y[1] = p->x.y; y[2] = p->x.z; // Full orbit particle position
        y[3] = p->u.x; y[4] = p->u.y; y[5] = p->u.z; // full orbit particle gamma*v
        status = gsl_odeiv2_driver_apply( p->Driver_Full, &t, dT, y );
        if ( status != GSL_SUCCESS ) { printf ("error, return value=%d\n", status); } 
        p->x.x = y[0]; p->x.y = y[1]; p->x.z = y[2];
        p->u.x = y[3]; p->u.y = y[4]; p->u.z = y[5]; 

    } else if ( p->CurrentIntegrator == GCA ) {

        y[0] = p->X.x; y[1] = p->X.y; y[2] = p->X.z; // GCA position
        y[3] = p->p_par;                             // GCA parallel momentum
        status = gsl_odeiv2_driver_apply( p->Driver_GCA, &t, dT, y );
        p->X.x = y[0]; p->X.y = y[1]; p->X.z = y[2]; // GCA position
        p->p_par = y[3];
    
    } else {
        printf( "Integrator not defined\n" );
        exit(0);
    }



}


/*
 *
 *      KineticEnergy in keV
 *
 *
 *
 */
int Lgm_InitParticle( int Species, double KineticEnergy, double Charge, double AlphaEq, Lgm_Vector *x0, Lgm_ParticleInfo *p ){

    double v_tot, v_par, v_per;
    long int Date;
    double   UTC;
    double Ra, MLT, Phi, Lat;
    Lgm_Vector u_sm, u0, w0;
    int Flag;
    Lgm_Vector v1, v2, v3, Bvec, bhat, B, E;
    Lgm_Vector  Rc_hat_pqb, Rc_hat_gsm, v_per_hat_pqb, v_per_hat_gsm;
    Lgm_Vector w_sm, w_off, ww, un;
    double SinAlpha, CosAlpha, g, Blocal, Beq, Tb_approx, t, dt, Rc, phi, gamma;

    double      gamma_n, KineticEnergy_n;


    p->KE      = KineticEnergy/1000.0;   // Initial Kinetic Energy, MeV
    p->AlphaEq = AlphaEq; g = sin( AlphaEq*RadPerDeg ); // Initial equatorial pitch angle
    
    if ( Species == ELECTRON ) {
        p->m = 9.10938291e-31;           // electron mass,  kg
        p->q = Charge*1.60217657e-19;    // Particle charge,  C
    } else if ( Species == PROTON ) {
        p->m = 1.67262178e-27;           // proton mass,  kg
        p->q = Charge*1.60217657e-19;    // Particle charge,  C
    }
        

    p->c          = 2.99792458e8;    // Speed of light, m/s
    p->c2         = p->c * p->c;       // c^2
    p->mc         = p->m * p->c;       // mc
    p->mc2        = p->m * p->c * p->c; // mc^2


// INIT_MISC
p->mInfo = Lgm_InitMagInfo( );
p->mInfo->Kp = 2;
Lgm_MagModelInfo_Set_MagModel( LGM_IGRF, LGM_EXTMODEL_T89c, p->mInfo );
//Lgm_MagModelInfo_Set_MagModel( LGM_CDIP, LGM_EXTMODEL_NULL, p->mInfo );

// UPDATE_TIME
Date = 20020418;
UTC  = 5.0 + 30.0/60.0 + 0.0/3600.0;
Lgm_Set_Coord_Transforms( Date, UTC, p->mInfo->c );

u_sm = *x0;
printf("x0 = %g %g %g\n", x0->x, x0->y, x0->z);
Lgm_Convert_Coords( &u_sm, &u0, SM_TO_GSM, p->mInfo->c );
Ra  = 6.60000;
MLT = -1.5;
Phi = (MLT*15.0 - 180.0)*RadPerDeg; Lat = 0.0*RadPerDeg;
u_sm.x = Ra*cos( Phi )*cos(Lat); u_sm.y = Ra*sin( Phi )*cos(Lat); u_sm.z = Ra*sin(Lat);
Lgm_Convert_Coords( &u_sm, &u0, SM_TO_GSM, p->mInfo->c );

w0 = u0;
Flag = Lgm_Trace( &w0, &v1, &v2, &v3, 120.0, 1e-7, 1e-7, p->mInfo );
p->mInfo->Bfield( &v3, &Bvec, p->mInfo ); Beq = Lgm_Magnitude( &Bvec );
p->mInfo->Bfield( &w0, &Bvec, p->mInfo ); B   = Bvec; bhat = Bvec; Blocal = Lgm_NormalizeVector( &bhat );
//printf("Blocal, Beq = %g %g\n", Blocal, Beq);
B.x *= 1e-9; B.y *= 1e-9; B.z *= 1e-9; 


    // GCA position
    p->X.x = u0.x*1000.0*Re;
    p->X.y = u0.y*1000.0*Re;
    p->X.z = u0.z*1000.0*Re;

    // GCA p_par
// this v_par is not right.
// We should have two tracks ...k
    v_tot    = sqrt( Lgm_v2overc2( p->KE, LGM_Ep0 ) ) * p->c;
    SinAlpha = sqrt( Blocal/Beq * g*g );
    CosAlpha = sqrt( 1.0 - SinAlpha*SinAlpha);
    v_par    = v_tot * CosAlpha;
    v_per    = v_tot * SinAlpha;
    gamma    = 1.0/sqrt(1.0-v_tot*v_tot/p->c2);
    p->p_par = gamma*p->m*v_par;









p->B0         = Blocal;       // Initial |B| in nT.
p->Alpha      = asin( sqrt( Blocal/Beq * g*g ) )*DegPerRad;
g = sin( p->Alpha * RadPerDeg );
double B_mirror = Blocal/(g*g);
v_par        = v_tot * cos( p->Alpha*RadPerDeg );
v_per        = v_tot * sin( p->Alpha*RadPerDeg );

p->mu         = Lgm_Ek_to_Mu( p->KE, p->Alpha, p->B0, LGM_Ep0 ); // mu in MeV/G
printf("E = %g MeV, Alpha = %g Degrees, B = %g nT,  LGM_Ee0 = %g MeV  Mu = %g MeV/G     B_mirror = %g\n", p->KE, p->Alpha, p->B0, LGM_Ep0, p->mu, B_mirror );
p->mu         *= 1.60217657e-9;  // mu in J/T


gamma_n = 1.0/sqrt(1.0-v_tot*v_tot/p->c2);
KineticEnergy_n = (gamma_n - 1.0)*p->m*p->c2*6.241509e+15;
printf("gamma_n, KE = %lf,   %lf\n", gamma_n, KineticEnergy_n);

double uuu;
uuu = v_tot*gamma_n;
KineticEnergy_n = p->m*uuu*uuu/(1.0+gamma_n)*6.241509e+15;
printf("uuu = %g   gamma_n = %g   KE = %lf\n", uuu, gamma_n, KineticEnergy_n);







    //printf("w0 = %g %g %g\n", w0.x, w0.y, w0.z);
    Lgm_Set_GSM_TO_PQB( &w0, p->mInfo ); // This will set up the PQB coord system for this location

    /*
     *  Need to pick a phase angle in the perp plane to get where to offset particle
     */
    phi = 90.0; Rc_hat_pqb.x = cos(phi*RadPerDeg); Rc_hat_pqb.y = sin(phi*RadPerDeg); Rc_hat_pqb.z = 0.0;
    Lgm_MatTimesVec( p->mInfo->Apqb_to_gsm, &Rc_hat_pqb, &Rc_hat_gsm ); // dont need to shift system for velocity
    //printf("Rc_hat_pqb = %g %g %g    Rc_hat_gsm = %g %g %g\n", Rc_hat_pqb.x, Rc_hat_pqb.y, Rc_hat_pqb.z,   Rc_hat_gsm.x, Rc_hat_gsm.y, Rc_hat_gsm.z);


    v_per_hat_pqb.x = -sin(phi*RadPerDeg); v_per_hat_pqb.y = cos(phi*RadPerDeg); v_per_hat_pqb.z = 0.0;;
    Lgm_MatTimesVec( p->mInfo->Apqb_to_gsm, &v_per_hat_pqb, &v_per_hat_gsm ); // dont need to shift system for velocity
    //printf("v_per_hat_pqb = %g %g %g    v_per_hat_gsm = %g %g %g\n", v_per_hat_pqb.x, v_per_hat_pqb.y, v_per_hat_pqb.z,   v_per_hat_gsm.x, v_per_hat_gsm.y, v_per_hat_gsm.z);



    /*
     * We need to offset the starting position of the particle based on its gyroradius.
     */
    gamma = 1.0/sqrt( 1.0 - v_tot*v_tot/p->c2);
    Rc = fabs(gamma*p->m*v_per/(p->q*Blocal*1e-9)); // particle gyro-radius in meters


    w_off.x = w0.x*1000.0*Re + Rc * Rc_hat_gsm.x; // m
    w_off.y = w0.y*1000.0*Re + Rc * Rc_hat_gsm.y; // m
    w_off.z = w0.z*1000.0*Re + Rc * Rc_hat_gsm.z; // m
    

    /*
     * Recompute local b-field related quantities.
     * 
     * We are shifting where the particle is (from the GC to its full orbit helix).
     * A problem with this is that the b-hast is slightly different at this ne location.
     * If we just use the previous b-hat and v-perp directions, thinks will be a bit off.
     * (You have to be careful here or KE wont be conserved.)
     */
    ww.x = w_off.x/1000.0/Re; ww.y = w_off.y/1000.0/Re; ww.z = w_off.z/1000.0/Re;
    
    // New (slightly different) b-hat.
    p->mInfo->Bfield( &ww, &bhat, p->mInfo ); Blocal = Lgm_NormalizeVector( &bhat ); 

    // New (slightly different) v_per_hat_pqb and v_per_hat_gsm vector.
    Lgm_Set_GSM_TO_PQB( &ww, p->mInfo ); 
    Lgm_MatTimesVec( p->mInfo->Apqb_to_gsm, &v_per_hat_pqb, &v_per_hat_gsm ); 

    // Recompute local pitch angle (because |Blocal| is slightly different).
    
    p->Alpha = asin( sqrt( Blocal/B_mirror ) )*DegPerRad;
    printf("Blocal, p->Alpha = %g, %g\n", Blocal, p->Alpha);

    // Recompute (slightly different partitionin of v_tot between v_par and v_per.
    v_par        = v_tot * cos( p->Alpha*RadPerDeg );
    v_per        = v_tot * sin( p->Alpha*RadPerDeg );


    // initial Full Orbit position
    p->x.x = w_off.x; p->x.y = w_off.y; p->x.z = w_off.z;

    // initial Full Orbit velocity vector u = gamma*v
    p->u.x = gamma*(v_par * bhat.x + v_per * v_per_hat_gsm.x);
    p->u.y = gamma*(v_par * bhat.y + v_per * v_per_hat_gsm.y);
    p->u.z = gamma*(v_par * bhat.z + v_per * v_per_hat_gsm.z);





}



