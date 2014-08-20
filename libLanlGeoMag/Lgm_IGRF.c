/*
 *     IGFC.c - Copyright (c) 1999-2006 Michael G. Henderson <mghenderson@lanl.gov>
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
 *
 *    This is my own implementation. It isnt *radically* optimized yet.
 *    But it does use recurrence relations for lots of stuff. 
 *    (Could probably be made faster by fine tuning and could also use
 *    Clenshaw's recurrence formula to do some of the sums...?)
 *
 *
 *    Be careful when testing this code against other IGRF codes. Many of the
 *    available versions assume that you are inputing "geodetic coordinates", not
 *    geocentric coordinates. So in some other codes, your Theta and Phi are assummed
 *    to be geodetic and then they are converted to geocentric, then B is calculated
 *    and then the B is converted back to geodetic. This is NOT what we want here.
 *    Here, I assume the input coords are in geocentric and the output field is also.
 *
 *    The Earth model is the WGS84 spheroid. (Equatorial radius a=6378.137 km,
 *    Polar Radius b=6356.752 km. These are also defined in LgmCtrans.h as WGS84_A
 *    and WGS84_B)
 *
 *    I think this is now reentrant and thread-safe. But we should do some more
 *    checking....
 *
 */

#include "Lgm/Lgm_CTrans.h"
#include "Lgm/Lgm_IGRF.h"
#define TINY 1.0e-25
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#if USE_OPENMP
#include <omp.h>
#endif




void Lgm_IGRF( Lgm_Vector *vin, Lgm_Vector *B, Lgm_CTrans *c ) {

    double      x[7], y[7], z[7], t[7], e;
    int         i, n;
    Lgm_Vector  v, w;


    /*
     * For various reasons, LGM adopted a definition of an Earth radius as the
     * WGS84_A value (the equatorial radius). However, in the mathematical
     * formulation of the IGRF model, the value of "a, the magnetic reference
     * spherical radius" in the scalar potential equation is taken as 6371.2
     * km.  Therefore, to get the proper r value into the IGRF equations, we
     * need to rescale the input r value. See Finlay et al., International
     * Geomagnetic Reference Field: the eleventh generation, Geophys. J. Int.
     * (2010) 183, 1216–1230, doi: 10.1111/j.1365-246X.2010.04804.x 
     */
    v = *vin;
    v.x *= Re;     // convert r back to km
    v.x /= 6371.2; // convert r to have units of IGRF_A









    /*
     *  Since there is a singularity at the pole, we should
     *  check to see how close we are to it. If we are too
     *  close, calculate B by interpolating to it from points
     *  away from the pole.
     *
     */
    if ( fabs( v.y*RadPerDeg ) < 1e-4 ) {
        for (n=0, i=-3; i<=3; ++i){
            if (i != 0) {
                t[n] = 1e-4*(double)i;
                //t[n] = 5.0*RadPerDeg*(double)i;
                w.x = v.x; w.y = t[n]; w.z = v.z;
                _Lgm_IGRF4( &w, B, c );
                x[n] = B->x; y[n] = B->y; z[n] = B->z;
                ++n;
            }
        }
        Lgm_PolFunInt( t, x, 6, v.y, &(B->x), &e);
        Lgm_PolFunInt( t, y, 6, v.y, &(B->y), &e);
        Lgm_PolFunInt( t, z, 6, v.y, &(B->z), &e);

    } else {

        _Lgm_IGRF4( &v, B, c );

    }

}



void    _Lgm_IGRF( Lgm_Vector *v, Lgm_Vector *B, Lgm_CTrans *c ) {

    double          r, Theta, Phi, B_r, B_theta, B_phi;
    double          st, ct, sp, cp, t;
    double          P[14][14], dP[14][14], Cmp[14], Smp[14], f1[14][14], f2[14], rinv;
    int             N;
    register double nsum, msum;
    register int    n, m;
    



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
    r     = v->x;
    Theta = v->y;
    Phi   = v->z;
    st = sin( Theta ); ct = cos( Theta );
    sp = sin( Phi );   cp = cos( Phi );




    /*
     *  IGRF 2005 goes up to N=13, but 10 ought to be good enough for our purposes?
     */
    N = 13;
    Lgm_InitIGRF(c->Lgm_IGRF_g, c->Lgm_IGRF_h, N, c->Lgm_IGRF_FirstCall, c);



    /*
     *  Compute all of the Legendre-related values.
     */
    Lgm_InitPnm( ct, st, c->Lgm_IGRF_R, P, dP, N, c );
    Lgm_InitTrigmp( cp, sp, Cmp, Smp, N );




    /*
     * precompute some quantities
     */
    {
#if USE_OPENMP
        //#pragma omp parallel private(m)
        //#pragma omp for
#endif
        for (n=1; n<=N; ++n){
            for (m=0; m<=n; ++m) f1[n][m] = c->Lgm_IGRF_g[n][m]*Cmp[m] + c->Lgm_IGRF_h[n][m]*Smp[m];
        }
    }

    /*
     *  Precompute f2_n = (1/r)^(n+2)
     */
    {
        rinv = 1.0/r;
        t = rinv*rinv;
        f2[0] = t;  // f2[0] is never used, could comment out?

        for (n=1; n<=N; ++n){
            t *= rinv;
            f2[n] = t;
        }
    }
 


    /*
     *   Compute B_r
     */
    {
        for (nsum=0.0, n=1; n<=N; ++n){
            for (msum=0.0, m=0; m<=n; ++m) msum += f1[n][m]*P[n][m];
            nsum += (n+1)*f2[n]*msum;
        }
    }
    B_r = nsum;




    /*
     *   Compute B_theta
     */
    for (nsum=0.0, n=1; n<=N; ++n){
        for (msum=0.0, m=0; m<=n; ++m) msum += f1[n][m]*dP[n][m];
        nsum += f2[n]*msum;
    }
    B_theta = -nsum;



    /*
     *   Compute B_phi
     */
    for (nsum=0.0, n=1; n<=N; ++n){
        for (msum=0.0, m=0; m<=n; ++m) msum += (c->Lgm_IGRF_g[n][m]*Smp[m] - c->Lgm_IGRF_h[n][m]*Cmp[m])*m*P[n][m];
        nsum += f2[n]*msum;
    }
    B_phi = nsum/st;




    B->x = B_r;
    B->y = B_theta;
    B->z = B_phi;

    c->Lgm_IGRF_FirstCall = FALSE;

}



/*
 *  Using a different recurrence relation
 */
void    _Lgm_IGRF2( Lgm_Vector *v, Lgm_Vector *B, Lgm_CTrans *c ) {

    double          r, Theta, Phi, B_r, B_theta, B_phi;
    double          st, ct, sp, cp, t;
    double          val, val2, Cmp[14], Smp[14], f2[14], rinv;
    int             N;
    register double P_n_m, P_nm1_m, P_nm1_mm1, P_nm2_m;
    register double dP_n_m, dP_nm1_m, dP_nm1_mm1, dP_nm2_m;
    register int    n, m;
    
    /*
     *  IGRF 2005 goes up to N=13, but 10 ought to be good enough for our purposes?
     */
    N = 13;
    Lgm_InitIGRF( c->Lgm_IGRF_g, c->Lgm_IGRF_h, N, c->Lgm_IGRF_FirstCall, c);



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
    r     = v->x;
    Theta = v->y;
    Phi   = v->z;
    st = sin( Theta ); ct = cos( Theta );
    sp = sin( Phi );   cp = cos( Phi );


    if ( c->Lgm_IGRF_FirstCall ) {
        Lgm_InitK( c->Lgm_IGRF_K, N );
        Lgm_InitS( c->Lgm_IGRF_S, N );
    }



    /*
     *  Precompute the Cos(m*phi) and Sin(m*phi) using recurance relations
     */
    Lgm_InitTrigmp( cp, sp, Cmp, Smp, N );


    /*
     *  Precompute f2_n = (1/r)^(n+2) efficiently
     */
    {
        rinv = 1.0/r;
        t = rinv*rinv;
        f2[0] = t;  // f2[0] is never used, could comment out?

        for (n=1; n<=N; ++n){
            t *= rinv;
            f2[n] = t;
        }
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
     */


    B_r = B_theta = B_phi = 0.0;
    P_nm1_mm1  = 1.0; P_nm1_m    = 1.0; // Initially these are P_0,0
    dP_nm1_mm1 = 0.0; dP_nm1_m   = 0.0; // Initially these are dP_0,0

    for ( m=0; m<= N; ++m ) {
        for ( n=1; n<= N; ++n ) {
            if ( n >= m ) {

                if ( n == m ) {
                    P_n_m = st*P_nm1_mm1;
                    dP_n_m = st*dP_nm1_mm1 + ct*P_nm1_mm1;
                    P_nm1_mm1  = P_n_m;  P_nm1_m  = P_n_m;  P_nm2_m  = 0.0;
                    dP_nm1_mm1 = dP_n_m; dP_nm1_m = dP_n_m;  dP_nm2_m = 0.0;
                } else if ( n == 1 ) {
                    P_n_m = ct*P_nm1_m; 
                    dP_n_m = ct*dP_nm1_m - st*P_nm1_m;
                    P_nm2_m  = P_nm1_m;  P_nm1_m  = P_n_m;
                    dP_nm2_m = dP_nm1_m; dP_nm1_m = dP_n_m;
                } else {
                    P_n_m = ct*P_nm1_m - c->Lgm_IGRF_K[n][m]*P_nm2_m;
                    dP_n_m = ct*dP_nm1_m - st*P_nm1_m - c->Lgm_IGRF_K[n][m]*dP_nm2_m;
                    P_nm2_m  = P_nm1_m;  P_nm1_m  = P_n_m;
                    dP_nm2_m = dP_nm1_m; dP_nm1_m = dP_n_m;
                }
                    
                val = c->Lgm_IGRF_g[n][m]*Cmp[m] + c->Lgm_IGRF_h[n][m]*Smp[m];
                val2 = c->Lgm_IGRF_S[n][m]*f2[n];
                B_r     += val2 * (double)(n+1)*val*P_n_m;
                B_theta += val2 * val*dP_n_m;
                B_phi   += val2 * (double)m*(-c->Lgm_IGRF_g[n][m]*Smp[m] + c->Lgm_IGRF_h[n][m]*Cmp[m])*P_n_m;
            
            }
        }
    }

    B->x = B_r;
    B->y = -B_theta;
    B->z = -B_phi/st;

    c->Lgm_IGRF_FirstCall = FALSE;

}



/*
 *  Same recurrence relation as _Lgm_IGRF2, but paralellized
 */
void    _Lgm_IGRF3( Lgm_Vector *v, Lgm_Vector *B, Lgm_CTrans *c ) {

    double          r, Theta, Phi, B_r, B_theta, B_phi;
    double          st, ct, sp, cp, t;
    double          val, val2, Cmp[14], Smp[14], f2[14], rinv, Pnn[14], dPnn[14], g, h;
    int             N, nmin;
    //register double P_n_m, P_nm1_m, P_nm1_mm1, P_nm2_m;
    //register double dP_n_m, dP_nm1_m, dP_nm1_mm1, dP_nm2_m;
    //register int    n, m;
    double P_n_m, P_nm1_m, P_nm1_mm1, P_nm2_m;
    double dP_n_m, dP_nm1_m, dP_nm1_mm1, dP_nm2_m;
    int    n, m;
    
    /*
     *  IGRF 2005 goes up to N=13, but 10 ought to be good enough for our purposes?
     */
    N = 13;
    Lgm_InitIGRF( c->Lgm_IGRF_g, c->Lgm_IGRF_h, N, c->Lgm_IGRF_FirstCall, c);



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
    t = rinv*rinv;
    f2[0] = t;  // f2[0] is never used, could comment out?

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



//    { // BEGIN PARALLEL
        // this granularity here is too small. Parallel version actually runs slower. But serially, its faster than 
        // previous versions.
#if USE_OPENMP
        //#pragma omp parallel firstprivate(m,n,P_nm1_m,P_nm2_m,dP_nm1_m,dP_nm2_m,P_n_m,dP_n_m,val,val2) 
        //#pragma omp for schedule(static,1) reduction(+:B_r,B_theta,B_phi)
#endif
        for ( m=0; m<= N; ++m ) {

            P_n_m    = Pnn[m];
            dP_n_m   = dPnn[m];
            P_nm1_m  = P_n_m;  P_nm2_m   = 0.0;
            dP_nm1_m = dP_n_m; dP_nm2_m  = 0.0;

            for ( n=m; n<= N; ++n ) {

                if ( n != m ) {
                    if ( n == 1 ) {
                        P_n_m = ct*P_nm1_m; 
                        dP_n_m = ct*dP_nm1_m - st*P_nm1_m;
                        P_nm2_m  = P_nm1_m;  P_nm1_m  = P_n_m;
                        dP_nm2_m = dP_nm1_m; dP_nm1_m = dP_n_m;
                    } else {
                        P_n_m = ct*P_nm1_m - c->Lgm_IGRF_K[n][m]*P_nm2_m;
                        dP_n_m = ct*dP_nm1_m - st*P_nm1_m - c->Lgm_IGRF_K[n][m]*dP_nm2_m;
                        P_nm2_m  = P_nm1_m;  P_nm1_m  = P_n_m;
                        dP_nm2_m = dP_nm1_m; dP_nm1_m = dP_n_m;
                    }
                }
                        
                if ( n > 0 ){
                    val = c->Lgm_IGRF_g[n][m]*Cmp[m] + c->Lgm_IGRF_h[n][m]*Smp[m];
                    val2 = c->Lgm_IGRF_S[n][m]*f2[n];

                    //B_r     = B_r     + (val2 * (double)(n+1)*val*P_n_m);
                    //B_theta = B_theta + (val2 * val*dP_n_m);
                    //B_phi   = B_phi   + (val2 * (double)m*(-c->Lgm_IGRF_g[n][m]*Smp[m] + c->Lgm_IGRF_h[n][m]*Cmp[m])*P_n_m);

                    B_r     += (val2 * (double)(n+1)*val*P_n_m);
                    B_theta += (val2 * val*dP_n_m);
                    B_phi   += (val2 * (double)m*(-c->Lgm_IGRF_g[n][m]*Smp[m] + c->Lgm_IGRF_h[n][m]*Cmp[m])*P_n_m);
                }
            
            }
        }

//    } // END PARALLEL

    B->x = B_r;
    B->y = -B_theta;
    B->z = -B_phi/st;

    c->Lgm_IGRF_FirstCall = FALSE;

}



/*
 *  Same recurrence relation as _Lgm_IGRF3, but attempt at making it more efficient
 */
void    _Lgm_IGRF4( Lgm_Vector *v, Lgm_Vector *B, Lgm_CTrans *c ) {

    double          r, Theta, Phi, B_r, B_theta, B_phi;
    double          st, ct, sp, cp, t, gnm, hnm, Knm;
    double          val, val2, val3, Cmp_m, Smp_m, Cmp[14], Smp[14], f2[14], rinv, Pnn[14], dPnn[14];
    int             N;
    register double P_n_m, P_nm1_m, P_nm1_mm1, P_nm2_m;
    register double dP_n_m, dP_nm1_m, dP_nm1_mm1, dP_nm2_m;
    register int    n, m;
    
    /*
     *  IGRF 2005 goes up to N=13, but 10 ought to be good enough for our purposes?
     */
    N = 13;
    Lgm_InitIGRF( c->Lgm_IGRF_g, c->Lgm_IGRF_h, N, c->Lgm_IGRF_FirstCall, c);



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

            gnm = c->Lgm_IGRF_g[n][m];
            hnm = c->Lgm_IGRF_h[n][m];

            if ( n != m ) {
                Knm = c->Lgm_IGRF_K[n][m];
                P_n_m = ct*P_nm1_m - Knm*P_nm2_m;
                dP_n_m = ct*dP_nm1_m - st*P_nm1_m - Knm*dP_nm2_m;
                P_nm2_m  = P_nm1_m;  P_nm1_m  = P_n_m;
                dP_nm2_m = dP_nm1_m; dP_nm1_m = dP_n_m;
            }
                    
            if ( n > 0 ){
                Cmp_m = Cmp[m]; Smp_m = Smp[m];
                val = gnm*Cmp_m + hnm*Smp_m;
                val2 = c->Lgm_IGRF_S[n][m]*f2[n];
                val3 = val2 * val;

                B_r     += (val3*(n+1)*P_n_m);
                B_theta += (val3*dP_n_m);
                B_phi   += (val2 * m*(-gnm*Smp_m + hnm*Cmp_m)*P_n_m);
            }
        
        }
    }



    B->x = B_r;
    B->y = -B_theta;
    B->z = -B_phi/st;

    c->Lgm_IGRF_FirstCall = FALSE;

}





void Lgm_InitPnm( double ct, double st, double R[14][14], double P[14][14], double dP[14][14], int N, Lgm_CTrans *c ) {

    double         Pmm, Pmp1m, Pnm, Pnm1m, Pnm2m, a, b, f, x, x2;
    int            NmM, Nm1, Mp1;
//    static double  TwoNm1_Over_NmM[13][13], NpMm1_Over_NmM[13][13];
//    static int     FirstTimeThrough=TRUE;
    register int   n, m, i;

    x = ct;

    /*
     *  Initialize
     */
/*
we dont really need this if we set everything below
*/
    if ( c->Lgm_IGRF_FirstCall ) {
        for (n=0; n<=N; ++n){
            for (m=0; m<=N; ++m){
                P[n][m] = 0.0;
                dP[n][m] = 0.0;
                R[n][m] = 0.0;
            }
        }
    }

    /*
     *  Now calculate the Schmidt normalization factors ("Ratios").
     *  There's almost certainly a better way to do this, but
     *  we really only ever need to do this once (they are constants).
     *
     *  One solution would be to precompute them an stuff them into 
     *  a static variable.
     */
    if ( c->Lgm_IGRF_FirstCall ) {
        for (n=0; n<=N; ++n){
            for (m=0; m<=n; ++m){

                NmM = n-m;
                Nm1 = n-1;

                if (m==0) {
                    R[n][m] = 1.0;
                } else {
                    R[n][m] =  sqrt( 2.0*Lgm_Factorial(NmM)/Lgm_Factorial(n+m) );
                }

                if (n>m+1){
                    c->Lgm_IGRF_TwoNm1_Over_NmM[n][m] = (double)(Nm1+n)/(double)NmM;
                    c->Lgm_IGRF_NpMm1_Over_NmM[n][m]  = (double)(Nm1+m)/(double)NmM;
                }

            }
        }
    } 


    x2 = x*x;
    b = 1.0 - x2;
    a = sqrt( b );


    /*
     *  Compute all of the Pnm
     */
    for (m=0; m<=N; ++m){

        /*
         *  For a given m, we need to first compute Pm,m
         *  Set it to 1 first. If m == 0, then leave it at
         *  1. If m>0, then compute Pmm
         */
        Pmm = 1.0;
        if (m>0) {
            f = 1.0;
            // are we computing the following too many times? I.e., could we be more efficient?
            for (i=1; i<=m; ++i){
                Pmm *= f*a; f += 2.0;
            }
        }
        P[m][m] = Pmm*R[m][m];


        if (m == N) {

            /*  In this case Pnmp1 should be zero -- we are done with
             *  this value of m. Move one to next value of m...
             */


        } else {

            Mp1 = m+1;
            Pmp1m = x*(m+Mp1)*Pmm;
            P[Mp1][m] = Pmp1m*R[Mp1][m];
            if (N == (Mp1)) {

                /*  There arent any more values to compute for 
                 *  this value of m. Move one to next value of m...
                 */


            } else {

                Pnm1m = Pmp1m;
                Pnm2m = Pmm;
                for (n=m+2; n<=N; ++n){

                    // Pnm  = ( x*(2*n-1)*Pnm1m - (n+m-1)*Pnm2m )/ (n-m);
                    Pnm     =  x*c->Lgm_IGRF_TwoNm1_Over_NmM[n][m]*Pnm1m - c->Lgm_IGRF_NpMm1_Over_NmM[n][m]*Pnm2m;
                    P[n][m] = Pnm*R[n][m];

                    Pnm2m = Pnm1m;
                    Pnm1m = Pnm;

                }

            }

        }

    }


    Lgm_InitdPnm( P, dP, N, c );


}



void Lgm_InitdPnm( double P[14][14], double dP[14][14], int N, Lgm_CTrans *c ) {

    int             n, m;
//    static double   SqrtNM1[13][13], SqrtNM2[13][13];
//    static int      FirstTimeThrough=TRUE;


//c->Lgm_IGRF_FirstCall = TRUE;
    if ( c->Lgm_IGRF_FirstCall ) Lgm_InitSqrtFuncs( c->Lgm_IGRF_SqrtNM1, c->Lgm_IGRF_SqrtNM2, N );

    for (n=0; n<=N; ++n){
        for (m=n; m>=0; --m){
    
            if (m==0) {

                //dP[n][0] = -sqrt(0.5*n*(n+1)) * P[n][1];
                dP[n][0] = c->Lgm_IGRF_SqrtNM1[n][m] * P[n][1];

            } else {

                //q = ( (m-1) == 0 ) ? 2 : 1;

                if (m==n){
                    //dP[n][m] = 0.5* sqrt((double)(q*(n+m)*(n-m+1))) * P[n][m-1];
                    dP[n][m] = c->Lgm_IGRF_SqrtNM1[n][m] * P[n][m-1];
                } else {
                    //dP[n][m] = 0.5*( sqrt((double)(q*(n+m)*(n-m+1))) * P[n][m-1] - sqrt((double)((n+m+1)*(n-m))) * P[n][m+1] );
                    dP[n][m] = c->Lgm_IGRF_SqrtNM1[n][m] * P[n][m-1] - c->Lgm_IGRF_SqrtNM2[n][m]* P[n][m+1];
                }

            }

        }
    }


}



void Lgm_InitSqrtFuncs( double SqrtNM1[14][14], double SqrtNM2[14][14], int N ) {

    int    n, m, q, NmM, NpM;

    for (n=0; n<=N; ++n){
        for (m=n; m>=0; --m){

            if (m==0) {

                SqrtNM1[n][m] = -sqrt(0.5*n*(n+1));

            } else {

                q = ( (m-1) == 0 ) ? 2 : 1;

                NpM = n+m;
                NmM = n-m;
                SqrtNM1[n][m] = 0.5*sqrt((double)(q*NpM*(NmM+1)));
                SqrtNM2[n][m] = 0.5*sqrt((double)((NpM+1)*NmM));

            }

        }
    }

}





double Lgm_Factorial( int k ) {

    double  f;
    int     i;

    if ( k == 0 ) return( 1.0 );

    for (f=1.0, i=2; i<=k; ++i) f *= (double)i;

    return( f );

}




/*
 *  Compute all of the sin(m phi) and cos(m phi) values
 *  efficiently using recurrences. Doesnt require any costly
 *  trig functions.
 */
void Lgm_InitTrigmp( double cp, double sp, double *Cmp, double *Smp, int N ) {

    double    b;
    double    Cm, Cmm1, Cmm2, Sm, Smm1, Smm2;
    int        m;

    b = 2.0*cp;
    Cmp[0] = 1.0;  Cmp[1] = cp;
    Cmm2 = Cmp[0]; Cmm1 = Cmp[1];
    Smp[0] = 0.0;  Smp[1] = sp;
    Smm2 = Smp[0]; Smm1 = Smp[1];
    for (m=2; m<=(N+1); ++m){
        Cm = b*Cmm1 - Cmm2;
        Cmp[m] = Cm;
        Cmm2 = Cmm1;
        Cmm1 = Cm;

        Sm = b*Smm1 - Smm2;
        Smp[m] = Sm;
        Smm2 = Smm1;
        Smm1 = Sm;
    }
    

}
void Lgm_InitS( double S[14][14], int N ) {

    register int n, m;
    int          nm1, nm1_2;

    for ( n=0; n<=N; n++ ) {
        for ( m=0; m<=N; m++ ) {
            S[n][m] = 0.0;
        }
    }

    
    S[0][0] = 1.0;

    for ( m=0; m<=N; m++ ) {
        for ( n=1; n<=N; n++ ) {

            if ( m<=n){
                if ( m==0 ) {
                    S[n][0] = S[n-1][0]*(double)(2*n-1)/(double)n;
                } else if ( m==1) {
                    S[n][m] = S[n][m-1]*sqrt( (double)(2*(n-m+1))/((double)(n+m)) );
                } else {
                    S[n][m] = S[n][m-1]*sqrt( (double)(n-m+1)/((double)(n+m)) );
                }
            }

//printf("S[%d][%d] = %g\n", n, m, S[n][m] );
        }
    }
    
    

}

void Lgm_InitK( double K[14][14], int N ) {

    register int n, m;
    int          nm1, nm1_2;


    for ( n=1; n<=N; n++ ) {
        for ( m=0; m<=N; m++ ) {
            if ( n == 1 ) {
                K[n][m] = 0.0;
            } else {
                nm1 = n-1; nm1_2 = nm1*nm1;
                K[n][m] = (double)((nm1_2 - m*m))/((double)((2*n-1)*(2*n-3)));
            }
//printf("K[%d][%d] = %g\n", n, m, K[n][m] );
        }
    }
    

}


void Lgm_InitIGRF( double g[14][14], double h[14][14], int N, int Flag, Lgm_CTrans *c ){

    double          Year;
    double          g0, g1, h0, h1, gs, hs, y0, y1;
    double          H0, H02, Lx, Ly, Lz, E;
    int             j, j0, j1, n, m;


    /* Get Year from Lgm_CTrans structure */
    Year = c->UTC.fYear;
    if ( Year < 1.0 ) {
        printf("Year not set?!   ( Year = %g )\n", Year);
        exit(-1);
    }


    /*
     *  Set IGRF Model based on the current epoch.
     */
    if ( (fabs(Year - c->Lgm_IGRF_OldYear) > 0.0) || Flag ) {


        if ((Year >= IGRF_epoch[0])&&(Year <= IGRF_epoch[IGRF_nModels-1])) {

            /*
             *   Interpolation.  
             *   Use linear for now. 
             *   Rational or Polynomial function extrapolation would be better.
             *
             */
            j = IGRF_nModels-1;
            while ( Year < IGRF_epoch[j] ){ --j; }
            j0 = j;
            j1 = j+1;

        } else {


            /*
             * extrapolation 
             */
            if (Year > IGRF_epoch[IGRF_nModels-1]){
                j0 = IGRF_nModels-2; j1 = IGRF_nModels-1;
            } else {
                j0 = 0; j1 = 1;
            }

        }

        if ( j1 >= IGRF_nModels ){
            j1 = IGRF_nModels-1;
            j0 = j1-1;
        }

        for (n=0; n<=N; ++n){
            for (m=0; m<=n; ++m){
                g0 = IGRF_g[j0][n][m]; g1 = IGRF_g[j1][n][m];
                h0 = IGRF_h[j0][n][m]; h1 = IGRF_h[j1][n][m];

                y0 = IGRF_epoch[j0];   y1 = IGRF_epoch[j1];
                gs = (g1-g0)/(y1-y0);  hs = (h1-h0)/(y1-y0);
            
                g[n][m] = gs*(Year - y0) + g0;
                h[n][m] = hs*(Year - y0) + h0;
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
        c->M_cd_McIllwain = 31165.3;
        c->M_cd_2010      = 29950.1686985232;

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

    } 

    c->Lgm_IGRF_OldYear = Year;


}




void Lgm_PolFunInt( double *xa, double *ya, int n, double x, double *y, double *dy ) {

    int     i, m, ns=0;
    double  den, dif, dift, ho, hp, w;
    double  *c, *d;

    c = (double *)calloc(n, sizeof(double));
    d = (double *)calloc(n, sizeof(double));

    dif = fabs(x - xa[0]);
    for (i=0; i<n; ++i) {
        if ( (dift = fabs(x-xa[i])) < dif) {
            ns  = i;
            dif = dift;
        }
        c[i] = ya[i];
        d[i] = ya[i];
    }

    *y = ya[ns--];

    for (m=1; m<n; ++m) {
        for (i=0; i<(n-m); ++i) {
            ho = xa[i]-x;
            hp = xa[i+m]-x;
            w  = c[i+1]-d[i];
            if ( (den=ho-hp) == 0.0) printf("Error in routine polint");
            den  = w/den;
            d[i] = hp*den;
            c[i] = ho*den;
        }
        *y += (*dy = (2*ns < (n-m-1) ? c[ns+1] : d[ns--]));
    }

    free(d);
    free(c);

}


void Lgm_RatFunInt( double *xa, double *ya, int n, double x, double *y, double *dy ) {

    int     i, m, ns=0;
    double  w, t, hh, h, dd;
    double  *c, *d;

    c = (double *)calloc(n, sizeof(double));
    d = (double *)calloc(n, sizeof(double));

    hh = fabs( x - xa[0] );
    for (i=0; i<n; i++) {
        h = fabs( x - xa[i] );
        if ( h == 0.0 ){
            *y  = ya[i];
            *dy = 0.0;
            free(c);
            free(d);
            return;
        } else if ( h< hh ) {
            ns = i;
            hh = h;
        }
        c[i] = ya[i];
        d[i] = ya[i] + TINY;
    }

    *y = ya[ns--];

    for (m=1; m<n; ++m) {
        for (i=0; i<(n-m); ++i) {
            w  = c[i+1] - d[i];
			h  = xa[i+m] - x;
			t  = (xa[i] - x)*d[i]/h;
			dd = t - c[i+1];
            if ( dd == 0.0 ) {
                printf("Error in routine Lgm_RatFunInt");
                *y  = ya[i];
                *dy = 0.0;
                free(c);
                free(d);
                return;
            }
            dd = w/dd;
            d[i]=c[i+1]*dd;
			c[i]=t*dd;
		}
		*y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
    }

    free(d);
    free(c);
    return;

}
