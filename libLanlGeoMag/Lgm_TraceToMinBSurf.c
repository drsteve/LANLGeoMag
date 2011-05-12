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

int Lgm_TraceToMinBSurf( Lgm_Vector *u, Lgm_Vector *v, double Htry, double tol, Lgm_MagModelInfo *Info ) {

    Lgm_Vector	u_scale;
    double	    Htry_max, Hdid, Hnext, Hmin, Hmax, s=0.0, sgn, r2;
    double	    Sa, Sb, Sc, s2, d1, d2;
    double	    Ba, Bb, Bc, B, B2, R;
    Lgm_Vector	Btmp;
    Lgm_Vector	Pa, Pb, Pc, P, P2;
    int		    done, reset=TRUE;



    Hmax = 0.5;
//    Hmin = 0.00001;
    Hmin = 1e-7;
    u_scale.x =  100.0;  u_scale.y = 100.0; u_scale.z = 100.0;


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

    done = FALSE; Sa = Sb = Sc = 0.0;

    /*
     *  Set the start point, Pa and |B(Pa)|;
     *  I.e., Pa, Ba, Sa
     */
    Pa   = *u;
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
    if ( Lgm_MagStep( &P, &u_scale, Htry, &Hdid, &Hnext, 1.0e-7, -1.0, &s, &reset, Info->Bfield, Info ) < 0 ) return(-1);
    Info->Bfield( &P, &Btmp, Info );
    B = Lgm_Magnitude( &Btmp );

    if ( B < Ba ) {

	    Pb  = P;
	    Bb  = B;
	    Sb  = Hdid;
	    sgn = -1.0;     // We should move in direction opposite the field direction.

    } else {

        P2 = Pa; //reset = TRUE;
        if ( Lgm_MagStep( &P2, &u_scale, Htry, &Hdid, &Hnext, 1.0e-7, 1.0, &s2, &reset, Info->Bfield, Info ) < 0 ) return(-1);
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
	        Pb   = P;  Bb   = B;  Sb   = s;
	        Pc   = P2; Bc   = B2; Sc   = s2;
	        sgn  = 1.0;
	        done = TRUE;
	    }

    }









    /*
     *  Keep stepping along the field line until we have a bracket.
     */
    while (!done) {

	    P = Pb;
        if ( Lgm_MagStep( &P, &u_scale, Htry, &Hdid, &Hnext, 1.0e-7, sgn, &s, &reset, Info->Bfield, Info ) < 0 ) return(-1);
        Info->Bfield( &P, &Btmp, Info );
        B = Lgm_Magnitude( &Btmp );

	    if ( B < Bb ) {
	        Pa = Pb; Ba = Bb; Sa = Sb;
	        Pb = P;  Bb = B;  Sb = Sa + Hdid;
            if (   (P.x > Info->OpenLimit_xMax) || (P.x < Info->OpenLimit_xMin) || (P.y > Info->OpenLimit_yMax) || (P.y < Info->OpenLimit_yMin)
                    || (P.z > Info->OpenLimit_zMax) || (P.z < Info->OpenLimit_zMin) ) {
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
	





    /*
     *  We have a bracket. Now go in for the kill.
     *  Use golden-section search to converge toward minimum.
     *  (Sa, Sb, Sc) are the distances of the triple points along
     *  the FL.
     */
if (0==1){
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
            if ( Lgm_MagStep( &P, &u_scale, Htry, &Hdid, &Hnext, 1.0e-7, sgn, &s, &reset, Info->Bfield, Info ) < 0 ) return(-1);
            Info->Bfield( &P, &Btmp, Info );
            B = Lgm_Magnitude( &Btmp );

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
            if ( Lgm_MagStep( &P, &u_scale, Htry, &Hdid, &Hnext, 1.0e-7, sgn, &s, &reset, Info->Bfield, Info ) < 0 ) return(-1);
            Info->Bfield( &P, &Btmp, Info );
            B = Lgm_Magnitude( &Btmp );

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
     */
if (1==1){
//printf("Sa, Sb, Sc = %g %g %g  Ba, Bb, Bc = %g %g %g   tol = %g\n", Sa, Sb, Sc, Ba, Bb, Bc, tol);
    double      Smin, Bmin;
    Lgm_Vector  Pmin;
    FuncInfo    f;

    f.u_scale = u_scale;
    f.Htry    = Htry;
    f.sgn     = sgn;
    f.reset   = reset;
    f.Info    = Info;
    Lgm_Brent( Sa, Sb, Sc, Bb, Pa, Pb, Pc, &f, tol, &Smin, &Bmin, &Pmin );
    Bb = Bmin;
    Sb = Smin;
    Pb = Pmin;
//printf("Sa, Sb, Sc = %g %g %g  Ba, Bb, Bc = %g %g %g   tol = %g\n", Sa, Sb, Sc, Ba, Bb, Bc, tol);
}





    /*
     *  Take location of the Min-B surface to be Pb.
     */
    *v = Pb;

    Info->Trace_s = Sb*sgn;

    if ( Info->VerbosityLevel > 2 ) printf("TraceToMinBSurf(): Number of Bfield evaluations = %d\n", Info->Lgm_nMagEvals );

    return( 1 );

}


