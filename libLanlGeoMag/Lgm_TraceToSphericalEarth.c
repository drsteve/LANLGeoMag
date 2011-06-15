/* TraceToSphericalEarth, Copyright (c) 2010 Michael G. Henderson <mghenderson@lanl.gov>
 *
 *    - Attempts to trace FL down to (some height above) the Earth.
 *
 *     This is a version of Lgm_TraceToEarth, but in this version the Earth is
 *     approximated by a sphere of radius WGS84_A (i.e. the equatorial radius
 *     in the WGS84 model). This is useful for computing magnetic flux of a
 *     drift shell since a spherical surface is much easier to compute B dot dA
 *     on. The magnetic flux doesnt care what surface we use, so using a
 *     spherical earth for that isnt an approximation...
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


double seFunc( Lgm_Vector *P, double TargetHeight, Lgm_MagModelInfo *Info ){

    Lgm_Vector  w;
    double      Height, F;

    Height = WGS84_A*(Lgm_Magnitude( P )-1.0);
    F =  Height - TargetHeight;

    return( F );

}



int Lgm_TraceToSphericalEarth( Lgm_Vector *u, Lgm_Vector *v, double TargetHeight, double sgn, double tol, Lgm_MagModelInfo *Info ) {

    Lgm_Vector	u_scale;
    double	    Htry, Htry_max, Hdid, Hnext, Hmin, Hmax, s;
    double	    Sa=0.0, Sc=0.0, d;
    double	    Rinitial, Fa, Fb, Fc, F;
    double      Height, StartHeight;
    double	    Height_a, Height_b, Height_c, HeightPlus, HeightMinus, direction;
    Lgm_Vector	Pa, Pc, P;
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
     *
     * Note that although we are trying to trace to a spherical Earth in this
     * routine, we still dont want to drop below the surface of the real Earth
     * (IGRF doesnt like that so much).
     *
     */
    if ( Rinitial < WGS84_B ) {

        // must be inside the Earth, which is no good -- bail with
        // LGM_INSIDE_EARTH error code
        return( LGM_INSIDE_EARTH );

    } else {

        // We are at least at or above WGS84_B. We could still be in trouble in
        // terms of being inside the Earth, so we have to check more
        // thouroughly now.  Determine Geodetic Height
        Lgm_WGS84_to_GeodHeight( u, &Height );

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
     *  (reset Height to be distance in km above the spherical approx of the Earth.)
     */
    Height = Rinitial - WGS84_A; // distance above spherical Earth
    AboveTargetHeight = ( Height > TargetHeight ) ? TRUE : FALSE;



    /*
     * Initialize some things
     */
    Pa.x = Pa.y = Pa.z = 0.0;
    Pc.x = Pc.y = Pc.z = 0.0;
    P = *u;
    Hmax = Info->Hmax;
    Hmin = 0.001;
    Hmin = 1e-8;
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
        Htry = 0.001;

        // sgn = +1
        P = *u; if ( Lgm_MagStep( &P, &u_scale, Htry, &Hdid, &Hnext, 1.0e-7, 1.0, &s, &reset, Info->Bfield, Info ) < 0 ) return(-1);
        HeightPlus = WGS84_A*(Lgm_Magnitude( &P )-1.0);

        // sgn = -1
        P = *u; if ( Lgm_MagStep( &P, &u_scale, Htry, &Hdid, &Hnext, 1.0e-7, -1.0, &s, &reset, Info->Bfield, Info ) < 0 ) return(-1);
        HeightMinus = WGS84_A*(Lgm_Magnitude( &P )-1.0);

        direction = ( HeightPlus > HeightMinus ) ? 1.0 : -1.0;



        /*
         *  We are already at or below target height. Trace until we are not.
         */
        done  = FALSE;
        //reset = TRUE;
        while ( !done ) {
            Htry = fabs(0.9*(TargetHeight - Height));	    // This computes Htry as 90% of the distance to the TargetHeight
            if (Htry > 0.1) Htry = 0.1; // If its bigger than 0.1 reset it to 0.1 -- to be safe.
            if ( Lgm_MagStep( &P, &u_scale, Htry, &Hdid, &Hnext, 1.0e-7, direction, &s, &reset, Info->Bfield, Info ) < 0 ) return(-1);
            Sa += Hdid;
            Info->Trace_s += Hdid;
            Height = WGS84_A*(Lgm_Magnitude( &P )-1.0);
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
     *  We want to stop in the ionosphere at a certain height above the Earth
     *  (assumed to be a spehere of radius WGS84_A in this routine).  We need
     *  to find 2 points along the field line such that the intersection of the
     *  FL with the sphere is guaranteed to lie between them. To do this, we
     *  need to find two points; Pa and Pc such that;
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
    Height_a = WGS84_A*(Lgm_Magnitude( &Pa )-1.0);
    Fa   = Height_a - TargetHeight;

    /*
     *  Get an initial Htry that is safe -- i.e. start off slowly
     *  We dont really know where we are, so be conservative on the first try.
     *  If if ( Lgm_MagStep() gives back an Hnext thats higher, we'll crank Htry up then...
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
    //reset = TRUE;
    while ( !done ) {

        if ( Lgm_MagStep( &P, &u_scale, Htry, &Hdid, &Hnext, 1.0e-7, sgn, &s, &reset, Info->Bfield, Info ) < 0 ) return(-1);
        Height = WGS84_A*(Lgm_Magnitude( &P )-1.0);
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
	        Fc = F;
	        Height_c = Height;
	        //Sc += Hdid;
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

    if ( ((Fc > 0.0) && (Fa > 0.0)) || ((Fc < 0.0) && (Fa < 0.0)) ) {
        // No bracket
        return(0);
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
    //reset = TRUE;
if (0==1){
    done  = FALSE;
    while (!done) {

        d = Sc - Sa;
	
        //if ( fabs(F) < tol ) {
        if ( fabs(d) < tol ) {
            done = TRUE;
        } else {

            P = Pa; //Htry = LGM_1M_1O_GOLD*d; // LGM_1M_1O_GOLD is 0.381966...
//            if ( Htry > 2.0*fabs(F) ) {
                //Htry = LGM_1M_1O_GOLD*d; // LGM_1M_1O_GOLD is 0.381966...
                Htry = 0.5*fabs(d); // LGM_1M_1O_GOLD is 0.381966...
//            }
            if ( Lgm_MagStep( &P, &u_scale, Htry, &Hdid, &Hnext, tol, sgn, &s, &reset, Info->Bfield, Info ) < 0 ) return(-1);
            Height = WGS84_A*(Lgm_Magnitude( &P )-1.0);
	        F =  Height - TargetHeight;
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
    //printf( "Lgm_TraceToSphericalEarth: v = %g %g %g\n", v->x, v->y, v->z);

}


if (1==1){
    /*
     * Try Brent's method
     */
//printf("Sa, Sb = %g %g  Fa, Fb = %g %g   tol = %g\n", Sa, Sb, Fa, Fb, tol);
    double      Sz, Fz;
    Lgm_Vector  Pz;
    BrentFuncInfoP    f;

    f.u_scale = u_scale;
    f.Htry    = Htry;
    f.sgn     = sgn;
    f.reset   = reset;
    f.Info    = Info;
    f.func    = &seFunc;
    f.Val     = TargetHeight;
    Lgm_zBrentP( Sa, Sc, Fa, Fc, Pa, Pc, &f, tol, &Sz, &Fz, &Pz );
    Fc = Fz;
    Sc = Sz;
    Pc = Pz;
//printf("Sa, Sb = %g %g  Fa, Fb = %g %g   tol = %g\n", Sa, Sb, Fa, Fb, tol);

    v->x = Pz.x; v->y = Pz.y; v->z = Pz.z;
    Info->Trace_s = Sz;

}




    if (Info->SavePoints) fprintf(Info->fp, "%f \t%f\t %f\t 2\n", v->x, v->y, v->z);

    if ( Info->VerbosityLevel > 2 ) printf("Lgm_TraceToSphericalEarth(): Number of Bfield evaluations = %d\n", Info->Lgm_nMagEvals );

    return( 1 );

}
