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
#include "Lgm/Lgm_WGS84.h"




int Lgm_TraceToEarth( Lgm_Vector *u, Lgm_Vector *v, double TargetHeight, double sgn, double tol, Lgm_MagModelInfo *Info ) {

    Lgm_Vector	u_scale;
    double	    Htry, Htry_max, Hdid, Hnext, Hmin, Hmax, s;
    double	    Sa=0.0, Sc=0.0, d;
    double	    Rinitial, Fa, Fb, Fc, F;
    double      Height, StartHeight;
    double	    Height_a, Height_b, Height_c, HeightPlus, HeightMinus, direction;
    Lgm_Vector	Pa, Pc, P, w;
    int		    done, reset, AboveTargetHeight;

    reset = TRUE;

    Info->Trace_s = 0.0;
    Sa = Sc = 0.0;

    /*
     * Determine our initial geocentric radius in km. (u is assumed to be in
     * units of Re where we define Re to be WGS84_A.)
     */
    Rinitial = WGS84_A*Lgm_Magnitude( u ); // km


    /*
     * The Earth is a spheroid. It is more flattened at the poles than at the
     * equator.  Check to see if we are initially at a height where we need to
     * worry about it. For WGS84 Earth shape model, the equatorial radius is
     * WGS84_A and polar radius is WGS84_B (WGS84_B is about 20km smaller than
     * WGS84_A).
     */
    if ( Rinitial < WGS84_B ) {

        // must be inside the Earth, which is no good -- bail with
        // LGM_INSIDE_EARTH error code
        return( LGM_INSIDE_EARTH );

    } else {

        // We are at least at or above WGS84_B. We could still be in trouble in
        // terms of being inside the Earth, so we have to check more
        // thouroughly now.  Determine Geodetic Height
        Lgm_Convert_Coords( u, &w, GSM_TO_WGS84, Info->c );
        Lgm_WGS84_to_GeodHeight( &w, &Height );
        //Height = WGS84_A*(Lgm_Magnitude( &w )-1.0);

        if ( Height < 0.0 )  {

            // inside the Earth, which is no good -- bail with error
            return( LGM_INSIDE_EARTH );

        }

    }


    /*
     *  If we get here, the initial point must be above the surface of the
     *  (spheroidal) Earth.
     *
     *  Now check to see if we are currently above or below the target height.
     */
    AboveTargetHeight = ( Height > TargetHeight ) ? TRUE : FALSE;



    /*
     * Initialize some things
     */
    Pa.x = Pa.y = Pa.z = 0.0;
    Pc.x = Pc.y = Pc.z = 0.0;
    P = *u;
    Hmax = Info->Hmax;
    Hmin = 0.001;
    u_scale.x =  100.0;  u_scale.y = 100.0; u_scale.z = 100.0;
    Height = Height_a = Height_b = Height_c = 0.0;
    F = Fa = Fb = Fc = 0.0;



    /*
     * Save initial point if we need to 
     */
    if (Info->SavePoints) fprintf(Info->fp, "%f \t%f\t %f\t 3\n", u->x, u->y, u->z);



    /*
     * If we are above the target height, a potential problem will occur if the
     * field line is open in the direction of tracing.  Fortunately, that
     * situation is easy to detect and is taken care of later...
     *
     * Here we need to test for and deal with the opposite problem. If we are
     * already below the target height, it is possible that the field line
     * never gets high enough to reach the target height. This is more tricky
     * to deal with because we can run the risk of getting below the surface of
     * the Earth.  We need to try to get back above the target height to test
     * for this.  If we cannot get above we need to bail out of the routine
     * altogether.
     */
    StartHeight = Height;
    if ( !AboveTargetHeight ){


        /*
         *  Determine which direction along the field line will take us higher.
         */
        Htry = 0.01;

        // sgn = +1
        P = *u; Lgm_MagStep( &P, &u_scale, Htry, &Hdid, &Hnext, 1.0e-7, 1.0, &s, &reset, Info->Bfield, Info );
        Lgm_Convert_Coords( &P, &w, GSM_TO_WGS84, Info->c );
        Lgm_WGS84_to_GeodHeight( &w, &HeightPlus );
        //HeightPlus = WGS84_A*(Lgm_Magnitude( &w )-1.0);

        // sgn = -1
        P = *u; Lgm_MagStep( &P, &u_scale, Htry, &Hdid, &Hnext, 1.0e-7, -1.0, &s, &reset, Info->Bfield, Info );
        Lgm_Convert_Coords( &P, &w, GSM_TO_WGS84, Info->c );
        Lgm_WGS84_to_GeodHeight( &w, &HeightMinus );
        //HeightMinus = WGS84_A*(Lgm_Magnitude( &w )-1.0);

        direction = ( HeightPlus > HeightMinus ) ? 1.0 : -1.0;



        /* 
         *  We are already at or below target height. Trace until we are not.
         */
        done  = FALSE;
        reset = TRUE;
        while ( !done ) {
            Htry = fabs(0.9*(TargetHeight - Height));	    // This computes Htry as 90% of the distance to the TargetHeight
            if (Htry > 0.1) Htry = 0.1; // If its bigger than 0.1 reset it to 0.1 -- to be safe.
            Lgm_MagStep( &P, &u_scale, Htry, &Hdid, &Hnext, 1.0e-7, direction, &s, &reset, Info->Bfield, Info );
            Sa += Hdid;
            Lgm_Convert_Coords( &P, &w, GSM_TO_WGS84, Info->c );
            Lgm_WGS84_to_GeodHeight( &w, &Height );
            F = Height - TargetHeight;
            if ( F > 0.0 ){
                done = TRUE;
            } else if ( Height < StartHeight ) {
                // We are going back down again -- Target Height unreachable? -- Bail out
                return( LGM_TARGET_HEIGHT_UNREACHABLE );
            }
        }

    }
  







    /*
     *  Bracket the zero first. 
     *  We want to stop in the ionosphere at a certain height above the Earth.
     *  We need to find 2 points along the field line such that the
     *  intersection of the FL with the spheroid is guaranteed to lie between
     *  them. To do this, we need to find two points; Pa and Pc such that;
     *  
     *		and  Height( Pa ) - TargetHeight   >   0.0
     *		and  Height( Pc ) - TargetHeight   <   0.0
     *
     *  I.e., F = Height - TargetHeight has opposite signs
     * 
     *  Set the start point, Pa and Fa. (Fa is difference between Height_a and
     *  TargetHeight.)
     *
     */
    Pa   = P;
    Lgm_Convert_Coords( &Pa, &w, GSM_TO_WGS84, Info->c );
    Lgm_WGS84_to_GeodHeight( &w, &Height_a );
    //Height_a = WGS84_A*(Lgm_Magnitude( &w )-1.0);

//printf( "Lgm_TraceToEarth: Pa = <%g, %g, %g> Re\n", Pa.x, Pa.y, Pa.z );
//printf( "Lgm_TraceToEarth: Height_a = %g\n", Height_a );

    /*
     *  Get an initial Htry that is safe -- i.e. start off slowly
     *  We dont really know where we are, so be conservative on the first try.
     *  If Lgm_MagStep() gives back an Hnext thats higher, we'll crank Htry up then...
     */
    Htry = 0.9*Height_a;	    // This computes Htry as 90% of the distance to the Earth's surface (could be small if we are already close!)
    if (Htry > 0.1) Htry = 0.1; // If its bigger than 0.1 reset it to 0.1 -- to be safe.
    


    /*
     *  Keep stepping along the field line until we drop below the target
     *  height, TargetHeight.  (Or bail if its open). This completes the endpoints of the
     *  bracket pair. We get these points first, because there are field lines
     *  that have more than one local minimum in fabs(Height-TargetHeight). So find Pa and Pc
     *  first and then complete the triple later. CHECK on this. I changed this
     *  from a minima search to a simple root finder. Is it OK still?
     */
    P     = Pa;
    done  = FALSE;
//    reset = TRUE;
    while ( !done ) {

        Lgm_MagStep( &P, &u_scale, Htry, &Hdid, &Hnext, 1.0e-7, sgn, &s, &reset, Info->Bfield, Info );
        Lgm_Convert_Coords( &P, &w, GSM_TO_WGS84, Info->c );
        Lgm_WGS84_to_GeodHeight( &w, &Height );
	    F =  Height - TargetHeight;
	    if ((F > 0.0) && (Info->SavePoints)) fprintf(Info->fp, "%f \t%f\t %f\t 2\n", P.x, P.y, P.z);

            if (   (P.x > Info->OpenLimit_xMax) || (P.x < Info->OpenLimit_xMin) || (P.y > Info->OpenLimit_yMax) || (P.y < Info->OpenLimit_yMin)
                || (P.z > Info->OpenLimit_zMax) || (P.z < Info->OpenLimit_zMin) ) {
	        /*
	         *  Open FL!
	         */
	        v->x = v->y = v->z = 0.0;
	        return(0);
	    } else if ( F < 0.0 ) {
	        done = TRUE;
	        Pc = P;
	        Height_c = Height;
	        Fc = F;
	        Sc = Sa + Hdid;
	    } else {
	        Pa = P;
	        Fa = F;
	        Height_a = Height;
	        Sa += Hdid;
	    }

        Htry = Hnext; // adaptively reset Htry

	    /*
	     *  Go no farther than some small distance below 
	     *  the target Height. Also respect Hmin and Hmax.
	     */
	    Htry_max = 0.9*Height;
	    if (Htry > Htry_max)  Htry = Htry_max;
	    if      (Htry < Hmin) Htry = Hmin;
	    else if (Htry > Hmax) Htry = Hmax;


    }




    /*
     *  We have a bracket. Now go in for the kill.
     *
     *  (Note: If all we wanted to know is whether or not the line hits the
     *  TargetHeight, we could stop here: it must or we wouldnt have a minimum
     *  bracketed.)
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
//printf( "Lgm_TraceToEarth: About to go in for kill \n");
    while (!done) {

        d = Sc - Sa;
	
        if ( fabs(F) < tol ) {
            done = TRUE;
        } else {

            P = Pa; Htry = 0.5*d;
            Lgm_MagStep( &P, &u_scale, Htry, &Hdid, &Hnext, 1.0e-7, sgn, &s, &reset, Info->Bfield, Info );
            Lgm_Convert_Coords( &P, &w, GSM_TO_WGS84, Info->c );
            Lgm_WGS84_to_GeodHeight( &w, &Height );
            //Height = WGS84_A*(Lgm_Magnitude( &w )-1.0);
	        F =  Height - TargetHeight;
//printf( "Lgm_TraceToEarth: Inside kill loop. Htry, F = %g %g %g\n", Htry, tol, F);
            if ( F >= 0.0 ) {
                Pa = P; Fa = F; Sa += Hdid;
            } else {
                Pc = P; Fc = F; Sc = Sa + Hdid;
            }

        }
    }
    Info->Trace_s = Sa;



    /*
     *  Take average as the final answer.
     */
    v->x = 0.5*(Pa.x + Pc.x); v->y = 0.5*(Pa.y + Pc.y); v->z = 0.5*(Pa.z + Pc.z);
    if (Info->SavePoints) fprintf(Info->fp, "%f \t%f\t %f\t 2\n", v->x, v->y, v->z);

    return( 1 );

}
