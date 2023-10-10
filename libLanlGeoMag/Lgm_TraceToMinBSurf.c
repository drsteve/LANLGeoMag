/*! \file Lgm_TraceToMinBSurf.c
 *
 *  \brief Routines to find the Minimum-B surface along a field line.
 *
 *
 *
 *  \author M.G. Henderson
 *  \date   1999
 *
 *
 *
 */


/* trace_test, Copyright (c) 1999 Michael G. Henderson <mghenderson@lanl.gov>
 *
 *
 *
 * Permission to use, copy, modify, distribute, and sell this software and its
 * documentation for any purpose is hereby granted without fee, provided that
 * the above copyright notice appear in all copies and that both that
 * copyright notice and this permission notice appear in supporting
 * documentation.  No representations are made about the suitability of this
 * software for any purpose.  It is provided "as is" without express or
 * implied warranty.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include "Lgm/Lgm_MagModelInfo.h"

/*
 *  TODO:
 *  TODO:
 *  TODO:
 *  TODO: In Lgm_Trace(), we now do a pre-trace with one mof the TraceLine()
 *  routiners in order to get an approximation to the global Bmin. Then
 *  Lgm_TraceToMinBSurf() is called which will converge on the nearest min --
 *  i.e. the global min. 
 *
 *  The problem with this is that running the Lgm_TraceToMinBSurf() routine
 *  standalone for such FLs will likely give different results (not always, but
 *  sometimes.)
 *
 *  The solution is to take that "pre-trace" stuff out of Lgm_Trace() and put
 *  it here so calls to Lgm_Trace() and Lgm_TraceToMinBSurf() yield the same
 *  results no matter what input position you start at.
 */


int Lgm_TraceToMinBSurf( Lgm_Vector *u, Lgm_Vector *v, double Htry, double tol, Lgm_MagModelInfo *Info ) {

    Lgm_Vector	u_scale;
    double	    Htry_max, Hdid, Hnext, Hmin, Hmax, s=0.0, sgn, r2;
    double	    Sa, Sb, Sc, d1, d2;
    double	    Ba, Bb, Bc, B, B2, R;
    Lgm_Vector	Btmp;
    Lgm_Vector	Pa, Pb, Pc, P, P2;
    int		    done, reset=TRUE;
    double s2 = 0.0;



    Hmax = 0.5;
//    Hmin = 0.00001;
    Hmin = 1e-7;
    Hmin = 1e-10;
    u_scale.x = u_scale.y = u_scale.z = 1.0;


    /*
     *  Bracket the minimum. We want to find two points along
     *  the field line such that the location of |B|_min is
     *  gauranteed to lie between them. To do this, we need to find
     *  three points; Pa, Pb, and Pc such that;
     *
     *		     Pc > Pb > Pa
     *
     *		and  |B( Pa )|   >   |B( Pb )|
     *		and  |B( Pc )|   >   |B( Pb )|
     *
     *
     */

    done = FALSE; 
    Sa = Sb = Sc = 0.0;
    Ba = Bb = Bc = 0.0;

    /*
     *  Set the start point, Pa and |B(Pa)|;
     *  I.e., Pa, Ba, Sa
     */
    Pa   = *u;
    R    = Lgm_Magnitude( &Pa );
    Info->Bfield( &Pa, &Btmp, Info );
    Ba   = Lgm_Magnitude( &Btmp );
    Sa   = 0.0;

//printf("P, B (initial) = %g %g %g %g\n", Pa.x, Pa.y, Pa.z, Ba); 

    /*
     *  Get an initial Htry that is safe -- i.e. start off slowly
     *  We dont really know where we are, so be conservative on the first try.
     *  If if ( Lgm_MagStep() gives back an Hnext thats higher, we'll crank Htry up then...
     */
    Htry = 0.9*(R-1.0);        // This computes Htry as 90% of the distance to the Earth's surface (could be small if we are already close!)
    if (Htry > 0.01) Htry = 0.01; // If its bigger than 0.01 reset it to 0.01 -- to be safe.





    /*
     *  Determine which way to go....
     *  Set Pb and |B(Pb)|.
     *  Try a step up the field line. If that doesnt give us
     *  a |B| less than Ba, then try down the field line.
     *  Use a fairly small stepsize to start with.
     *  This determined which direction we should move in (via value of sgn)
     */
    P = Pa;
    //reset = TRUE;
    if ( Lgm_MagStep( &P, &u_scale, Htry, &Hdid, &Hnext, -1.0, &s, &reset, Info->Bfield, Info ) < 0 ) return(-1);
    Info->Bfield( &P, &Btmp, Info );
    B = Lgm_Magnitude( &Btmp );
//printf("NEG: P, B  = %g %g %g %g\n", P.x, P.y, P.z, B); 

    if ( B < Ba ) {

	    Pb  = P;
	    Bb  = B;
	    Sb  = Hdid;
	    sgn = -1.0;     // We should move in direction opposite the field direction.

    } else {

        P2 = Pa; //reset = TRUE;
        if ( Lgm_MagStep( &P2, &u_scale, Htry, &Hdid, &Hnext, 1.0, &s2, &reset, Info->Bfield, Info ) < 0 ) return(-1);
        Info->Bfield( &P2, &Btmp, Info );
        B2 = Lgm_Magnitude( &Btmp );
//printf("POS: P2, B2  = %g %g %g %g\n", P2.x, P2.y, P2.z, B2); 

	    if ( B2 < Ba ) {
	        Pb  = P2;
            Bb  = B2;
            Sb  = Hdid;
	        sgn = 1.0;      // We should move in direction with the field direction.
	    } else {
	        /*
	         *  We must have already bracketed the min.
	         */
	        Pb   = Pa;  Bb = Ba; Sb = Sa;
	        Pa   = P;   Ba = B;  Sa = -s;
	        Pc   = P2;  Bc = B2; Sc = s2;
	        sgn  = 1.0;
	        done = TRUE;
	    }

    }

//printf("B = %g %g %g   done = %d\n", Ba, Bb, Bc, done);








    /*
     *  Keep stepping along the field line until we have a bracket.
     */
    while (!done) {

	    P = Pb;
        if ( Lgm_MagStep( &P, &u_scale, Htry, &Hdid, &Hnext, sgn, &s, &reset, Info->Bfield, Info ) < 0 ) return(-1);
        Info->Bfield( &P, &Btmp, Info );
        B = Lgm_Magnitude( &Btmp );

	    if ( B < Bb ) {
	        Pa = Pb; Ba = Bb; Sa = Sb;
	        Pb = P;  Bb = B;  Sb = Sa + Hdid;
            if (   (P.x > Info->OpenLimit_xMax) || (P.x < Info->OpenLimit_xMin) || (P.y > Info->OpenLimit_yMax) || (P.y < Info->OpenLimit_yMin)
                    || (P.z > Info->OpenLimit_zMax) || (P.z < Info->OpenLimit_zMin) || ( s > 1000.0 ) ) {
		        /*
		         *  Open FL!
		         */
		        v->x = v->y = v->z = 0.0;
		        return(0);
	        } else {
                Htry = Hnext;   // adaptively reset Htry

                /*
                 *  Dont attempt steps bigger than 0.25*distance to Earths surface
                 *  Also respect Hmin and Hmax.
                 */
                R = Lgm_Magnitude( &P );
                Htry_max = 0.25*fabs(R-1.0);
                if (Htry > Htry_max)  Htry = Htry_max;
                if      (Htry < Hmin) Htry = Hmin;
                else if (Htry > Hmax) Htry = Hmax;
	        }
	    } else {
	        Pc = P;  Bc = B; Sc = Sb + Hdid;
	        done = TRUE;
	    }


    }
	
    if ( (Bb > Bc) || (Bb > Ba) ) {
        // no bracket
        return(0);
    }
//printf("Sa, Sb, Sc = %lf %lf %lf\n", Sa, Sb, Sc);
//printf("Ba, Bb, Bc = %lf %lf %lf\n", Ba, Bb, Bc);




    /*
     *  We have a bracket. Now go in for the kill.
     *  Use golden-section search to converge toward minimum.
     *  (Sa, Sb, Sc) are the distances of the triple points along
     *  the FL.
     */
if (1==1){
    done = FALSE;
//reset=TRUE;
    while (!done) {

	    d1 = Sb - Sa;
	    d2 = Sc - Sb;
	    if ( (Sc-Sa) < tol ) {

	        done = TRUE;

	    } else if ( d1 > d2 ) {

	        //P = Pa; Htry = 0.381966011*d1;
	        P = Pa; Htry = 0.5*d1;
//printf("A. Sa, Sb, Sc = %g %g %g   d1, d2 = %g %g Htry = %g tol = %g Sc-Sa = %g\n", Sa, Sb, Sc, d1, d2, Htry, tol, Sc-Sa);
            if ( Lgm_MagStep( &P, &u_scale, Htry, &Hdid, &Hnext, sgn, &s, &reset, Info->Bfield, Info ) < 0 ) return(-1);
            Info->Bfield( &P, &Btmp, Info );
            B = Lgm_Magnitude( &Btmp );
//printf("A. B = %g\n", B);

	        if ( B < Bb ) {
                Pc = Pb; Bc = Bb;  Sc = Sb;
		        Pb = P; Bb = B; Sb = Sa + Hdid;
	        } else {
		        Pa = P; Ba = B; Sa += Hdid;
	        }
	
	    } else {

	        //P = Pb; Htry = 0.381966011*d2;
	        P = Pb; Htry = 0.5*d2;
//printf("B. Sa, Sb, Sc = %g %g %g   d1, d2 = %g %g Htry = %g tol = %g Sc-Sa = %g\n", Sa, Sb, Sc, d1, d2, Htry, tol, Sc-Sa);
            if ( Lgm_MagStep( &P, &u_scale, Htry, &Hdid, &Hnext, sgn, &s, &reset, Info->Bfield, Info ) < 0 ) return(-1);
            Info->Bfield( &P, &Btmp, Info );
            B = Lgm_Magnitude( &Btmp );
//printf("B. P = %g %g %g B = %g\n", P.x, P.y, P.z, B);

	        if ( B < Bb ) {
                Pa = Pb; Ba = Bb;  Sa = Sb;
		        Pb = P; Bb = B; Sb += Hdid;
	        } else {
		        Pc = P; Bc = B; Sc = Sb + Hdid;
	        }

	    }

    }
}

//reset = TRUE;

    /*
     * Try Brent's method
THIS DOESNT SEEM TO GIVE RESULTS THAT ARE AS GOOD. FIND OUT WHY....
     */
if (0==1){
//printf("Sa, Sb, Sc = %g %g %g  Ba, Bb, Bc = %g %g %g   tol = %g\n", Sa, Sb, Sc, Ba, Bb, Bc, tol);
    double      Smin, Bmin;
    Lgm_Vector  Pmin;
    BrentFuncInfoP    f;

    f.u_scale = u_scale;
    f.Htry    = Htry;
    f.sgn     = sgn;
    f.reset   = reset;
    f.Info    = Info;
    Lgm_BrentP( Sa, Sb, Sc, Bb, Pa, Pb, Pc, &f, tol, &Smin, &Bmin, &Pmin );
    Bb = Bmin;
    Sb = Smin;
    Pb = Pmin;
//printf("Sa, Sb, Sc = %g %g %g  Ba, Bb, Bc = %.15lf %.15lf %.15lf   tol = %g\n", Sa, Sb, Sc, Ba, Bb, Bc, tol);
}





    /*
     *  Take location of the Min-B surface to be Pb.
     */
    *v = Pb;

    //Info->Trace_s = Sb*sgn;
    Info->Trace_s = -Sb*sgn;

    if ( Info->VerbosityLevel > 2 ) printf("TraceToMinBSurf(): Number of Bfield evaluations = %ld\n", Info->Lgm_nMagEvals );

    return( 1 );

}


/*
 * This is a very similar version, but uses knowledge of where the northern
 * and/or southern footpoints are. Used by Lgm_Trace -- may not be convenient
 * otherwise.
 *
 *      v1    -- is the southern footpoint if it exists.
 *      v2    -- is the northern footpoint if it exists.
 *      Flag1 -- true if v1 valid, flase if not.
 *      Flag2 -- true if v2 valid, flase if not.
 *      S_FL  -- FL length. Either this is distance along FL between v1 and v2, or
 *               it is distance from v1 or v2 to the boundary (where FLs are consered
 *               open if they cross.)
 *      v     -- will be returned as the position of the Bmin point.
 *
 */
int Lgm_TraceToMinBSurf2( Lgm_Vector *v1, Lgm_Vector *v2, int Flag1, int Flag2, double S_FL, Lgm_Vector *v3, double *s3, double Htry, double tol, Lgm_MagModelInfo *Info ) {

    Lgm_Vector	u_scale;
    double	    Htry_max, Hdid, Hnext, Hmin, Hmax, s=0.0, sgn, r2;
    double	    Sa, Sb, Sc, d1, d2;
    double	    Ba, Bb, Bc, B, B2, R;
    Lgm_Vector	Btmp;
    Lgm_Vector	Pa, Pb, Pc, P, P2;
    int		    done, reset=TRUE;
    double      s2 = 0.0;



    Hmax = 0.5;
Hmax = 10.0;
//    Hmin = 0.00001;
    Hmin = 1e-7;
    Hmin = 1e-10;
    u_scale.x = u_scale.y = u_scale.z = 1.0;



   /*
    *    Pre-trace the FL to find approximate location of the true global Bmin. This is necessary for FLs
    *    that have multiple minima on the -- probably wasteful otherwise.
    *    We could also add code to detect number of minima here.
    *    
    *    1. Copy the Info structure so we dont mess anything up.
    *    2. Start from the southern footpoint (v1) and trace up to north (we know the distance.)
    *    3. save the various info as well as how far we are in s.
    */
    double tSS, tBmin, s_approx;
    Lgm_Vector u_approx;
    int        ii, iBmin;
    Lgm_MagModelInfo *Info2 = Lgm_CopyMagInfo( Info );

    //printf("Info2->s[%d] = %g\n", iBmin, Info2->s[iBmin]);
    tSS = S_FL;
    
    if ( Flag1 ) {
        // v1 (southern footpoint) is valid.
        //Lgm_TraceLine3( v1, tSS, 500, 1.0, 1e-7, FALSE, Info2 );
        Lgm_TraceLine3( v1, tSS, 500, 1.0, 1e-1, FALSE, Info2 );
    } else if ( Flag2 ) {
        // v2 (northern footpoint) is valid.
        //Lgm_TraceLine3( v2, tSS, 500, 1.0, 1e-7, FALSE, Info2 );
        Lgm_TraceLine3( v2, tSS, 500, 1.0, 1e-1, FALSE, Info2 );
    } else {
        // neither v1 nort v2 is valid, why did we get called?
        Lgm_FreeMagInfo( Info2 ); // free Info2 before we go.
        return(-1); 
    }
    tBmin = 9e99;
    iBmin = -1;
    for (ii=0; ii<Info2->nPnts; ii++){
        if (Info2->Bmag[ii] <= tBmin) {
            tBmin = Info2->Bmag[ii];
            iBmin = ii;
        }
    }
    if ( (iBmin >= 0) && (iBmin < Info2->nPnts) ){
        u_approx.x = Info2->Px[iBmin];
        u_approx.y = Info2->Py[iBmin];
        u_approx.z = Info2->Pz[iBmin];
        s_approx   = Info2->s[iBmin];  // distancev along FL from v1 to approximate global Bmin.
    } else {
        u_approx = *v1;
        s_approx = 0.0;
    }
    //printf("u_approx = %g %g %g\n", u_approx.x, u_approx.y, u_approx.z);
    //printf("Info2->s[%d] = %g\n", iBmin, Info2->s[iBmin]);
    //printf("Info2->Bmag[%d] = %g\n", iBmin-1, Info2->Bmag[iBmin-1]);
    //printf("Info2->Bmag[%d] = %g\n", iBmin, Info2->Bmag[iBmin]);
    //printf("Info2->Bmag[%d] = %g\n", iBmin+1, Info2->Bmag[iBmin+1]);
    Lgm_FreeMagInfo( Info2 );

    
    /*
     *  Bracket the minimum. We want to find two points along
     *  the field line such that the location of |B|_min is
     *  gauranteed to lie between them. To do this, we need to find
     *  three points; Pa, Pb, and Pc such that;
     *
     *		     Pc > Pb > Pa
     *
     *		and  |B( Pa )|   >   |B( Pb )|
     *		and  |B( Pc )|   >   |B( Pb )|
     *
     *
     */

    done = FALSE; 
    Sa = Sb = Sc = 0.0;
    Ba = Bb = Bc = 0.0;

    /*
     *  Set the start point, Pa and |B(Pa)|;
     *  I.e., Pa, Ba, Sa
     */
    //Pa   = *u;
    Pa   = u_approx;  // estimate for location of Bmin
    R    = Lgm_Magnitude( &Pa );
    Info->Bfield( &Pa, &Btmp, Info );
    Ba   = Lgm_Magnitude( &Btmp ); 
    Sa   = 0.0;


    /*
     *  Get an initial Htry that is safe -- i.e. start off slowly
     *  We dont really know where we are, so be conservative on the first try.
     *  If if ( Lgm_MagStep() gives back an Hnext thats higher, we'll crank Htry up then...
     */
    Htry = 0.9*(R-1.0);        // This computes Htry as 90% of the distance to the Earth's surface (could be small if we are already close!)
    if (Htry > 0.01) Htry = 0.01; // If its bigger than 0.01 reset it to 0.01 -- to be safe.


    /*
     *  Determine which way to go....
     *  Set Pb and |B(Pb)|.
     *  Try a step up the field line. If that doesnt give us
     *  a |B| less than Ba, then try down the field line.
     *  Use a fairly small stepsize to start with.
     *  This determined which direction we should move in (via value of sgn)
     */
    P = Pa;
    //reset = TRUE;
    if ( Lgm_MagStep( &P, &u_scale, Htry, &Hdid, &Hnext, -1.0, &s, &reset, Info->Bfield, Info ) < 0 ) return(-1);
    Info->Bfield( &P, &Btmp, Info );
    B = Lgm_Magnitude( &Btmp );

    if ( B < Ba ) {

	    Pb  = P;
	    Bb  = B;
	    Sb  = Hdid;
	    sgn = -1.0;     // We should move in direction opposite the field direction.

    } else {

        P2 = Pa; //reset = TRUE;
        if ( Lgm_MagStep( &P2, &u_scale, Htry, &Hdid, &Hnext, 1.0, &s2, &reset, Info->Bfield, Info ) < 0 ) return(-1);
        Info->Bfield( &P2, &Btmp, Info );
        B2 = Lgm_Magnitude( &Btmp );

	    if ( B2 < Ba ) {
	        Pb  = P2;
            Bb  = B2;
            Sb  = Hdid;
	        sgn = 1.0;      // We should move in direction with the field direction.
	    } else {
	        /*
	         *  We must have already bracketed the min.
	         */
	        Pb   = Pa;  Bb = Ba; Sb = Sa;
	        Pa   = P;   Ba = B;  Sa = -s;
	        Pc   = P2;  Bc = B2; Sc = s2;
	        sgn  = 1.0;
	        done = TRUE;
	    }

    }



    /*
     *  Keep stepping along the field line until we have a bracket.
     */
    while (!done) {

	    P = Pb;
        if ( Lgm_MagStep( &P, &u_scale, Htry, &Hdid, &Hnext, sgn, &s, &reset, Info->Bfield, Info ) < 0 ) return(-1);
        Info->Bfield( &P, &Btmp, Info );
        B = Lgm_Magnitude( &Btmp );

	    if ( B < Bb ) {
	        Pa = Pb; Ba = Bb; Sa = Sb;
	        Pb = P;  Bb = B;  Sb = Sa + Hdid;
            if (   (P.x > Info->OpenLimit_xMax) || (P.x < Info->OpenLimit_xMin) || (P.y > Info->OpenLimit_yMax) || (P.y < Info->OpenLimit_yMin)
                    || (P.z > Info->OpenLimit_zMax) || (P.z < Info->OpenLimit_zMin) || ( s > 1000.0 ) ) {
		        /*
		         *  Open FL!
		         */
		        //v3->x = v3->y = v3->z = 0.0;
                // return bext estimate so far
                *v3 = P; 
                Info->Bmin = B;
                
		        return(0);
	        } else {
                Htry = Hnext;   // adaptively reset Htry

                /*
                 *  Dont attempt steps bigger than 0.25*distance to Earths surface
                 *  Also respect Hmin and Hmax.
                 */
                R = Lgm_Magnitude( &P );
                Htry_max = 0.25*fabs(R-1.0);
                if (Htry > Htry_max)  Htry = Htry_max;
                if      (Htry < Hmin) Htry = Hmin;
                else if (Htry > Hmax) Htry = Hmax;
	        }
	    } else {
	        Pc = P;  Bc = B; Sc = Sb + Hdid;
	        done = TRUE;
	    }


    }
	
    if ( (Bb > Bc) || (Bb > Ba) ) {
        // no bracket
        return(0);
    }


    /*
     *  We have a bracket. Now go in for the kill.
     *  Use golden-section search to converge toward minimum.
     *  (Sa, Sb, Sc) are the distances of the triple points along
     *  the FL.
     */
if (1==1){
    done = FALSE;
    //reset=TRUE;
    while (!done) {

	    d1 = Sb - Sa;
	    d2 = Sc - Sb;
	    if ( (Sc-Sa) < tol ) {

	        done = TRUE;

	    } else if ( d1 > d2 ) {

	        //P = Pa; Htry = 0.381966011*d1;
	        P = Pa; Htry = 0.5*d1;
            //printf("A. Sa, Sb, Sc = %g %g %g   d1, d2 = %g %g Htry = %g tol = %g Sc-Sa = %g\n", Sa, Sb, Sc, d1, d2, Htry, tol, Sc-Sa);
            if ( Lgm_MagStep( &P, &u_scale, Htry, &Hdid, &Hnext, sgn, &s, &reset, Info->Bfield, Info ) < 0 ) return(-1);
            Info->Bfield( &P, &Btmp, Info );
            B = Lgm_Magnitude( &Btmp );
            //printf("A. B = %g\n", B);

	        if ( B < Bb ) {
                Pc = Pb; Bc = Bb;  Sc = Sb;
		        Pb = P; Bb = B; Sb = Sa + Hdid;
	        } else {
		        Pa = P; Ba = B; Sa += Hdid;
	        }
	
	    } else {

	        //P = Pb; Htry = 0.381966011*d2;
	        P = Pb; Htry = 0.5*d2;
            //printf("B. Sa, Sb, Sc = %g %g %g   d1, d2 = %g %g Htry = %g tol = %g Sc-Sa = %g\n", Sa, Sb, Sc, d1, d2, Htry, tol, Sc-Sa);
            if ( Lgm_MagStep( &P, &u_scale, Htry, &Hdid, &Hnext, sgn, &s, &reset, Info->Bfield, Info ) < 0 ) return(-1);
            Info->Bfield( &P, &Btmp, Info );
            B = Lgm_Magnitude( &Btmp );
            //printf("B. P = %g %g %g B = %g\n", P.x, P.y, P.z, B);

	        if ( B < Bb ) {
                Pa = Pb; Ba = Bb;  Sa = Sb;
		        Pb = P; Bb = B; Sb += Hdid;
	        } else {
		        Pc = P; Bc = B; Sc = Sb + Hdid;
	        }

	    }

    }
}

    //reset = TRUE;

    /*
     * Try Brent's method
THIS DOESNT SEEM TO GIVE RESULTS THAT ARE AS GOOD. FIND OUT WHY....
     */
if (0==1){
//printf("Sa, Sb, Sc = %g %g %g  Ba, Bb, Bc = %g %g %g   tol = %g\n", Sa, Sb, Sc, Ba, Bb, Bc, tol);
    double      Smin, Bmin;
    Lgm_Vector  Pmin;
    BrentFuncInfoP    f;

    f.u_scale = u_scale;
    f.Htry    = Htry;
    f.sgn     = sgn;
    f.reset   = reset;
    f.Info    = Info;
    Lgm_BrentP( Sa, Sb, Sc, Bb, Pa, Pb, Pc, &f, tol, &Smin, &Bmin, &Pmin );
    Bb = Bmin;
    Sb = Smin;
    Pb = Pmin;
//printf("Sa, Sb, Sc = %g %g %g  Ba, Bb, Bc = %.15lf %.15lf %.15lf   tol = %g\n", Sa, Sb, Sc, Ba, Bb, Bc, tol);
}





    /*
     *  Take location of the Min-B surface to be Pb.
     */
    *v3 = Pb;

    //Info->Trace_s = Sb*sgn;
    //Info->Trace_s = -Sb*sgn;
    *s3 = s_approx - Sb*sgn;
    //printf("s_approx, s3 = %g %g\n", s_approx, *s3);

    if ( Info->VerbosityLevel > 2 ) printf("TraceToMinBSurf(): Number of Bfield evaluations = %ld\n", Info->Lgm_nMagEvals );

    return( 1 );

}


