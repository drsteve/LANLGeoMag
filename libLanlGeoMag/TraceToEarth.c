/* TraceToEarth, Copyright (c) 1999 Michael G. Henderson <mghenderson@lanl.gov>
 *
 *    - Attempts to trace FL down to (some height above) the Earth.
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



int Lgm_TraceToEarth( Lgm_Vector *u, Lgm_Vector *v, double H0, double sgn, double tol, Lgm_MagModelInfo *Info ) {

    Lgm_Vector	u_scale;
    double	    Htry, Htry_max, Hdid, Hnext, Hmin, Hmax, s;
    double	    Sa=0.0, Sc=0.0, d;
    double	    Rtarget, R, Fa, Fb, Fc, F;
    double	    Ra, Rb, Rc;
    Lgm_Vector	Pa, Pc, P;
    int		    done, reset;
    Pa.x = Pa.y = Pa.z = 0.0;
    Pc.x = Pc.y = Pc.z = 0.0;
    P.x = P.y = P.z = 0.0;


    if (Info->SavePoints) fprintf(Info->fp, "%f \t%f\t %f\t 3\n", u->x, u->y, u->z);

    /*
     *  H0 is in km above Earth's surface. Convert to geocentric
     *  radius.
     */
    Rtarget = H0/Re + 1.0;



    Hmax = Info->Hmax;
    Hmin = 0.001;
    u_scale.x =  100.0;  u_scale.y = 100.0; u_scale.z = 100.0;
    R = Ra = Rb = Rc = 0.0;
    F = Fa = Fb = Fc = 0.0;


    /*
     *   Add code to check to see if the start point is already below the target height.
     *   If it is, we can trace up until we are beyond it or some such thing...
     */
    done  = FALSE;
    reset = TRUE;
    P = *u;
    R = Lgm_Magnitude( &P );
    F = R - Rtarget;
    if (F <= 0.0 ){
        /* 
         *  We are already at or below target height.
         *  Trace until we are not.
         */
        Htry  = 0.1;
        while ( !done ) {
            Lgm_MagStep( &P, &u_scale, Htry, &Hdid, &Hnext, 1.0e-7, sgn, &s, &reset, Info->Bfield, Info );
            R = Lgm_Magnitude( &P );
            F = R - Rtarget;
            if ( F > 0.0 ){
                done = TRUE;
            }
        }

    }
  







    /*
     *  Bracket the zero first. 
     *  We want to stop in the ionosphere which we represent as a sphere at a
     *  certain height above the Earth (we assume the Earth is a sphere not an
     *  oblate spheroid). We need to find 2 points along the field line such
     *  that the intersection of the FL with the sphere is guaranteed to lie
     *  between them.  To do this, we need to find three points; Pa and Pc such
     *  that;
     *  
     *		     Pc > Pa
     *
     *		and  R( Pa ) - Rtarget   >   0.0
     *		and  R( Pc ) - Rtarget   <   0.0
     *
     *
     * 
     *  Set the start point, Pa and Fa. (Fa is Euclidean distance between sphere and
     *  point Pa.)
     *
     */
    Pa   = P;
    Ra   = Lgm_Magnitude( &Pa );
    Sa   = 0.0;


    /*
     *  Get an initial Htry that is safe -- i.e. start off slowly
     *  We dont really know where we are, so be conservative on the first try.
     *  If Lgm_MagStep() gives back an Hnext thats higher, we'll crank Htry up then...
     */
    Htry = 0.9*(Ra-1.0);	    // This computes Htry as 90% of the distance to the Earth's surface (could be small if we are already close!)
    if (Htry > 0.1) Htry = 0.1; // If its bigger than 0.1 reset it to 0.1 -- to be safe.
    


    /*
     *  Keep stepping along the field line until we drop below the target
     *  radius, Rtarget.  (Or bail if its open). This completes the endpoints of the
     *  bracket pair. We get these points first, because there are field lines
     *  that have more than one local minimum in fabs(R-Rtarget). So find Pa and Pc
     *  first and then complete the triple later. CHECK on this. I changed this
     *  from a minima search to a simple root finder. Is it OK still?
     */
    P     = Pa;
    done  = FALSE;
    reset = TRUE;
    while ( !done ) {

        Lgm_MagStep( &P, &u_scale, Htry, &Hdid, &Hnext, 1.0e-7, sgn, &s, &reset, Info->Bfield, Info );
        R = Lgm_Magnitude( &P );
	    F =  R - Rtarget;
	    if ((F > 0.0) && (Info->SavePoints)) fprintf(Info->fp, "%f \t%f\t %f\t 2\n", P.x, P.y, P.z);

            if (   (P.x > Info->OpenLimit_xMax) || (P.x < Info->OpenLimit_xMin) || (P.y > Info->OpenLimit_yMax) || (P.y < Info->OpenLimit_yMin)
                || (P.z > Info->OpenLimit_zMax) || (P.z < Info->OpenLimit_zMin) ) {
	        /*
	         *  Open FL!
	         */
	        v->x = v->y = v->z = 0.0;
	        return(0);
	    } else if ( F < 0.00 ) {
	        done = TRUE;
	        Pc = P;
	        Rc = R;
	        Fc = F;
	        Sc += Hdid;
	    } else {
	        Pa = P;
	        Fa = F;
	        Ra = R;
	        Sa = 0.0;
	    }

        Htry = Hnext; // adaptively reset Htry

	    /*
	     *  Go no farther than some small distance below 
	     *  the target Radius. Also respect Hmin and Hmax.
	     */
	    Htry_max = 0.9*fabs(R-1.0);
	    if (Htry > Htry_max)  Htry = Htry_max;
	    if      (Htry < Hmin) Htry = Hmin;
	    else if (Htry > Hmax) Htry = Hmax;


    }




    /*
     *  We have a bracket. Now go in for the kill.
     *
     *  (Note: If all we wanted to know is whether or not the line hits the
     *  Earth, we could stop here: it must hit the Earth or we wouldnt
     *  have a minimum bracketed.)
     *
     *  Use golden-section search to converge toward minimum.
     *  (Sa, Sb, Sc) are the distances of the triple points along
     *  the FL. For convenience, we maintain Sa = 0.0, so that
     *  Sb and Sc are always the distances relative to Pa. (And
     *  so Sc is always the width of the bracketed interval. So when
     *  Sc gets small enough we bail out and take Pb as the min).
     *
     */
    done  = FALSE;
    reset = TRUE;
    while (!done) {

        d = Sc - Sa;
	
        if ( fabs(F) < tol ) {
            done = TRUE;
        } else {

            P = Pa; Htry = 0.5*d;
            Lgm_MagStep( &P, &u_scale, Htry, &Hdid, &Hnext, 1.0e-7, sgn, &s, &reset, Info->Bfield, Info );
            R = Lgm_Magnitude( &P );
            F =  R - Rtarget;
            if ( F >= 0.0 ) {
                Pa = P; Fa = F; Sa += Hdid;
            } else {
                Pc = P; Fc = F; Sc = Sa + Hdid;
            }

        }
    }



    /*
     *  Take average as the final answer.
     */
    v->x = 0.5*(Pa.x + Pc.x); v->y = 0.5*(Pa.y + Pc.y); v->z = 0.5*(Pa.z + Pc.z);
    if (Info->SavePoints) fprintf(Info->fp, "%f \t%f\t %f\t 2\n", v->x, v->y, v->z);

    return( 1 );

}

/*
 *   $Id$
 */
