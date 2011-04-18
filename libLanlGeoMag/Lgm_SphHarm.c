/*
 *     Lgm_SphHarm.c - Copyright (c) 2011 Michael G. Henderson <mghenderson@lanl.gov>
 *
 *     Permission to use, copy, modify, distribute, and sell this software and its
 *     documentation for any purpose is hereby granted without fee, provided that
 *     the above copyright notice appear in all copies and that both that
 *     copyright notice and this permission notice appear in supporting
 *     documentation.  No representations are made about the suitability of this
 *     software for any purpose.  It is provided "as is" without express or
 *     implied warranty.
 *
 *
 *
 */

#include "Lgm/Lgm_CTrans.h"
#include "Lgm/Lgm_SphHarm.h"
#define TINY 1.0e-25
#include <omp.h>

double Lgm_fac1( int n ){                                                                                                                                                                            
    int     i;                                                                                                                                                                                 
    double  f;                                                                                                                                                                                 
    if (n<=0)return(1.0);                                                                                                                                                                            
    for (f=1.0, i=2; i<=n; i++) f *= i;                                                                                                                                                              
    return(f);                                                                                                                                                                                     
}                                                                                                                                                                                                  
                                                                                                                                                                                                   
double Lgm_fac2( int n ){                                                                                                                                                                            
    int     i;                                                                                                                                                                                 
    double  f;                                                                                                                                                                                 
    if (n<=0)return(1.0);                                                                                                                                                                            
    for (f=1.0, i=1; i<=n; i += 2) f *= i;                                                                                                                                                           
    return(f);                                                                                                                                                                                     
}                                                                                                                                                                                                  
                                                                                                                                                                                                   
double Lgm_kdelta( int i, int j ) {                                                                                                                                                                        
    if ( i==j ) return( 1.0 );                                                                                                                                                                       
    else return( 0.0 );                                                                                                                                                                              
}                           

void  Lgm_GaussToSchmidtSemiNorm( int n, int m, double c_gauss, double *c_schmidt, double *Snm ) {
    double f1 = Lgm_fac1(n-m);
    *Snm = Lgm_fac2( 2*n-1 )/f1 * sqrt( (2.0-Lgm_kdelta(m,0))*f1/Lgm_fac1(n+m) );
    *c_schmidt = c_gauss/(*Snm);
    return;
}
void  Lgm_SchmidtSemiNormToGauss( int n, int m, double c_schmidt, double *c_gauss, double *Snm ) {
    double f1 = Lgm_fac1(n-m);
    *Snm = Lgm_fac2( 2*n-1 )/f1 * sqrt( (2.0-Lgm_kdelta(m,0))*f1/Lgm_fac1(n+m) );
    *c_gauss = c_schmidt * *Snm;
    return;
}



void Lgm_JensenCain1960( Lgm_Vector *v, Lgm_Vector *B, Lgm_CTrans *c ) {

    double    x[7], y[7], z[7], t[7], e;
    int       i, n;
    Lgm_Vector    w;

    /*
     *  Since there is a singularity at the pole, we should
     *  check to see how close we are to it. If we are too
     *  close, calculate B by interpolating to it from points
     *  away from the pole.
     *
     */
    if ( fabs( v->y*RadPerDeg ) < 1e-4 ) {
        for (n=0, i=-3; i<=3; ++i){
            if (i != 0) {
                t[n] = 1e-4*(double)i;
                //t[n] = 5.0*RadPerDeg*(double)i;
                w.x = v->x; w.y = t[n]; w.z = v->z;
                _Lgm_JensenCain1960( &w, B, c );
                x[n] = B->x; y[n] = B->y; z[n] = B->z;
                ++n;
            }
        }
        Lgm_PolFunInt( t, x, 6, v->y, &(B->x), &e);
        Lgm_PolFunInt( t, y, 6, v->y, &(B->y), &e);
        Lgm_PolFunInt( t, z, 6, v->y, &(B->z), &e);

    } else {

        _Lgm_JensenCain1960( v, B, c );

    }

}



/*
 */
void    _Lgm_JensenCain1960( Lgm_Vector *v, Lgm_Vector *B, Lgm_CTrans *c ) {

    double          r, Theta, Phi, B_r, B_theta, B_phi;
    double          st, ct, sp, cp, t;
    double          val, val2, Cmp[13], Smp[13], f2[13], rinv, Pnn[13], dPnn[13];
    int             N;
    register double P_n_m, P_nm1_m, P_nm1_mm1, P_nm2_m;
    register double dP_n_m, dP_nm1_m, dP_nm1_mm1, dP_nm2_m;
    register int    n, m;
    
    /*
     *  Jensen Cain goes up to N=6.
     */
    N = SphHarm_MaxN[LGM_JENSEN_CAIN_1960];
    Lgm_InitSphHarm( LGM_JENSEN_CAIN_1960, c->Lgm_IGRF_g, c->Lgm_IGRF_h, N, c->Lgm_IGRF_FirstCall, c);




    /*
     *  (r, theta, phi) are the spherical input coords.
     *
     *         theta is geographic co-latitude in degrees.
     *         phi is geographic longitude in degrees.
     *
     *          st = sin( Theta )       ct = cos( Theta )
     *
     *          sp = sin( Phi )         cp = cos( Phi )
     *
     *          r in units of Re
     */

    /*
     *  Compute cos(theta), sin(theta), cos(phi), sin(phi)
     *  Precompute the Cos(m*phi) and Sin(m*phi) using recurance relations
     */
    st = sin( v->y ); ct = cos( v->y);
    sp = sin( v->z ); cp = cos( v->z );
    Lgm_InitTrigmp( cp, sp, Cmp, Smp, N );



    /*
     *  Precompute f2_n = (1/r)^(n+2) efficiently
     */
    rinv = 1.0/v->x;
    f2[0] = t = rinv*rinv; // f2[0] is never used(?)

    for (n=1; n<=N; ++n){
        t *= rinv;
        f2[n] = t;
    }



    if ( c->Lgm_IGRF_FirstCall ) {
        Lgm_InitK( c->Lgm_IGRF_K, N );
        Lgm_InitS( c->Lgm_IGRF_S, N );
    }



    /*
     *  Use the recursions:
     *
     *      P0,0 = 1
     *      Pn,n = sin(theta)Pn-1,n-1
     *      Pn,m = cos(theta)Pn-1,m - Kn,m Pn-2,m
     *      Kn,m = ((n-1)^2 - m^2)/((2n-1)(2n-3)) n >1
     *      Kn,m = 0  n=1
     *
     *      dP0,0 = 0
     *      dPn,n = sin(theta)dPn-1,n-1 + cos(theta)*Pn-1,n-1
     *      dPn,m = cos(theta)dPn-1,m - sin(theta)Pn-1,m - Kn,m dPn-2,m
     *
     *
     *  So for each m, we start with a Pnn and work up.
     *
     *         n=0   n=1   n=2   n=3         n=N
     * 
     * m=0     P00   P10   P20   P30   ...   PN0
     * m=1           P11   P21   P31   ...   PN1
     * m=2                 P22   P32   ...   PN2
     * m=3                       P33   ...   PN3
     *                                        .
     *                                        .
     *                                        .
     * m=N                                   PNN
     * 
     *
     * This can be parallelized if we precompute all of the P_n,n's first.
     * Then all of the m's can be done independantly.
     *
     *
     */

    // precompute P_n_n's and dP_n_n's
    Pnn[0] = P_nm1_mm1 = 1.0;
    dPnn[0] = dP_nm1_mm1 = 0.0;
    for ( m=1; m<= N; ++m ) {
        Pnn[m] = st*P_nm1_mm1;
        dPnn[m] = st*dP_nm1_mm1 + ct*P_nm1_mm1;
        P_nm1_mm1  = Pnn[m]; dP_nm1_mm1 = dPnn[m];
    }


    // initialize sums
    B_r = B_theta = B_phi = 0.0;

    for ( m=0; m<= N; ++m ) {

        P_n_m    = Pnn[m];
        dP_n_m   = dPnn[m];
        P_nm1_m  = P_n_m;  P_nm2_m   = 0.0;
        dP_nm1_m = dP_n_m; dP_nm2_m  = 0.0;

        for ( n=m; n<= N; ++n ) {

            if ( n != m ) {
                P_n_m = ct*P_nm1_m - c->Lgm_IGRF_K[n][m]*P_nm2_m;
                dP_n_m = ct*dP_nm1_m - st*P_nm1_m - c->Lgm_IGRF_K[n][m]*dP_nm2_m;
                P_nm2_m  = P_nm1_m;  P_nm1_m  = P_n_m;
                dP_nm2_m = dP_nm1_m; dP_nm1_m = dP_n_m;
            }
                    
            if ( n > 0 ){
                val = c->Lgm_IGRF_g[n][m]*Cmp[m] + c->Lgm_IGRF_h[n][m]*Smp[m];
                val2 = c->Lgm_IGRF_S[n][m]*f2[n];

                B_r     += (val2 * (double)(n+1)*val*P_n_m);
                B_theta += (val2 * val*dP_n_m);
                B_phi   += (val2 * (double)m*(-c->Lgm_IGRF_g[n][m]*Smp[m] + c->Lgm_IGRF_h[n][m]*Cmp[m])*P_n_m);
            }
        
        }
    }



    B->x = B_r;
    B->y = -B_theta;
    B->z = -B_phi/st;

    c->Lgm_IGRF_FirstCall = FALSE;

}


void Lgm_InitSphHarm( int Model, double g[13][13], double h[13][13], int N, int Flag, Lgm_CTrans *c ){

    double          Year;
    double          g0, g1, h0, h1, gs, hs, y0, y1;
    double          H0, H02, Lx, Ly, Lz, E;
    int             j, j0, j1, n, m;


    Year = c->UTC.fYear;
    if ( Year < 1.0 ) {
        printf("Year not set?!   ( Year = %g )\n", Year);
        exit(-1);
    }


    /*
     *  Set Model based purely on the coeffs.
     */
//    if ( (fabs(Year - c->Lgm_IGRF_OldYear) > 0.0) || Flag ) {


        for (n=0; n<=N; ++n){
            for (m=0; m<=n; ++m){
                g[n][m] = SphHarm_g[Model][n][m]; 
                h[n][m] = SphHarm_h[Model][n][m];
            }
        }


        /*
         *   Compute the various IGRF dependent things like position of CD 
         */
        H02 = g[1][0]*g[1][0] + g[1][1]*g[1][1] + h[1][1]*h[1][1];
        H0  = sqrt(H02);

        /*
         *  Compute dipole moments.
         */
        c->M_cd = H0;
        c->M_cd_McIllwain = 0.311653e5;;

        c->CD_gcolat = M_PI - acos(g[1][0]/H0);
        c->CD_glon   = atan(h[1][1]/g[1][1]);


        /*
         *   Compute the Eccentric dipole offset vector Chapman and Bartels,
         *   [1962] is a good reference on this. But, at least the version of
         *   the text that I have appears to have alot of typos (e.g. the
         *   indices on one of the g coeffs are reversed and the L values are
         *   mixed up.) A better reference is Akasofu and Chapman [1972].
         *   Plus, I'm sure Jacobs' books are good on this too.
         */
        Lx = -g[1][1]*g[2][0] + (g[1][1]*g[2][2] + h[1][1]*h[2][2] + g[1][0]*g[2][1])*M_SQRT_3;
        Ly = -h[1][1]*g[2][0] + (g[1][1]*h[2][2] - h[1][1]*g[2][2] + g[1][0]*h[2][1])*M_SQRT_3;
        Lz = 2.0*g[1][0]*g[2][0] + (g[1][1]*g[2][1] + h[1][1]*h[2][1])*M_SQRT_3;
        E  = (Lx*g[1][1] + Ly*h[1][1] + Lz*g[1][0])/(4.0*H02);

        c->ED_x0 = (Lx-g[1][1]*E)/(3.0*H02); // in units of Re
        c->ED_y0 = (Ly-h[1][1]*E)/(3.0*H02); // in units of Re
        c->ED_z0 = (Lz-g[1][0]*E)/(3.0*H02); // in units of Re

//    } 

    c->Lgm_IGRF_OldYear = Year;


}


/*
 *   $Id: Lgm_SphHarm.c 139 2011-01-27 21:01:34Z mgh $
 */
