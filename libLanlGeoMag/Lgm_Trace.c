/* 
 * Trace, Copyright (c) 1999 Michael G. Henderson <mghenderson@lanl.gov>
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


/*
 *  Attempt to trace to Earth in both directions.
 *  If the field line is closed, then also trace to equatorial plane.
 *  
 *      v1 is northern footprint (at 120 km)
 *      v2 is southern footprint (at 120 km)
 *      v3 is equatorial crossing point  (Minimum B)
 *  
 *  
 *  Return values. Function returns:
 *  
 *  
 *      LGM_OPEN_IMF     (==0) if field line is open at both ends (i.e. an IMF FL)
 *      LGM_CLOSED       (==1) if field line is closed (v1, v2, v3 all valid)
 *      LGM_OPEN_N_LOBE  (==2) if field line is northern lobe FL (v1 valid)
 *      LGM_OPEN_S_LOBE  (==3) if field line is southern lobe FL (v2 valid)
 *
 *      LGM_INSIDE_EARTH (==-1) initial point is inside Earth (no points valid)
 *
 *
 *  Changed June 13, 2008
 *      Added TOL1 and TOL2 args. They wer hard set at 0.01 and 1e-7
 *  Changed October 14, 2008
 *      Added Height argument.  (Height is height in km above Earth to stop at)
 *  Original version quite old (early 90's?)
 *  
 */
int Lgm_Trace( Lgm_Vector *u, Lgm_Vector *v1, Lgm_Vector *v2, Lgm_Vector *v3, double Height, double TOL1, double TOL2, Lgm_MagModelInfo *Info ) {

    int		    i, reset, flag1, flag2, InitiallyBelowTargetHeight, done;
    double	    sgn=1.0, R, Rtarget, Rinitial, Rplus;
    Lgm_Vector	w, Bvec;
    double      h, h2_inv, F[7], s, Hdid, Hnext, Htry;
    double      GeodLat, GeodLong, GeodHeight;
    Lgm_Vector  u_scale, P;



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

    } else if ( Rinitial < WGS84_A ) {

        // we are between WGS84_A and WGS84_B ( i.e. WGS84_B < Rinitial <
        // WGS84_A).  we could still be in trouble in terms of being inside the
        // Earth, so we have to check more thouroughly now. Determine Geodetic
        // Lat, Lon and Height
        Lgm_WGS84_to_GEOD( u, &GeodLat, &GeodLong, &GeodHeight );

        if ( GeodHeight < 0.0 )  {

            // inside the Earth, which is no good -- bail with error
            return( LGM_INSIDE_EARTH );

        }

    }


    // determine the desired target height in km
    Rtarget  = WGS84_A + Height;


    /*
     * Note, that if the initial point is already below the target height it
     * must be attached to the Earth at least at that point.  So keep track of
     * status of initial point. Also in this case, we will have trouble trying
     * to trace in one of the directions.  A solution in this case is to try to
     * trace in the proper direction to get use above the target height -- then
     * attempt the trace in both directions.
     */

    // determine if we are initially below the target height
    InitiallyBelowTargetHeight = ( Rinitial <= Rtarget ) ? TRUE : FALSE; 

    // determine which direction will take us above the target R
    if (InitiallyBelowTargetHeight){
        Htry = 0.01; // step finely
        reset = TRUE;
        P = *u;
        Lgm_MagStep( &P, &u_scale, Htry, &Hdid, &Hnext, 1.0e-7, 1.0, &s, &reset, Info->Bfield, Info );
        Rplus = Lgm_Magnitude( &P );
        if (Rplus > Rinitial ) {
            sgn =  1.0;
        } else {
            sgn = -1.0;
        }

        done = FALSE;
        Htry = 0.01; // step finely
        reset = TRUE;
        P = *u;
        while( !done ) {
            Lgm_MagStep( &P, &u_scale, Htry, &Hdid, &Hnext, 1.0e-7, sgn, &s, &reset, Info->Bfield, Info );
            R = Lgm_Magnitude( &P );
            if (R > Rtarget) {
                done = TRUE;
            } else if ( R < Rinitial ) {
                done = TRUE;
                return( -1 );
            }
            
        }

    }
    







    /*
     * Try to get to the Earth by tracing along the field in either direction.
     */
    sgn = -1.0;
    Info->Hmax = 10.0;
Info->Hmax = 0.50;
Info->Hmax = 0.10;
    flag2 = Lgm_TraceToEarth(  u, v2, Height, -sgn, TOL1, Info );
    flag1 = Lgm_TraceToEarth(  u, v1, Height,  sgn, TOL1, Info );

    

    Info->Smin = -9e99;
    Info->Bmin = -9e99;


    if ( flag1 && flag2 ) {

	    /*
	     *  Closed FL -- attempt to trace to Eq Plane.
	     */
        Lgm_TraceToMinBSurf( v1, v3, TOL1, TOL2, Info );
        Info->Pmin = *v3;
        Info->Smin = Info->Trace_s;     // save location of Bmin. NOTE:  Smin is measured from the southern footpoint.
        Info->Bfield( v3, &Bvec, Info );
        Info->Bvecmin = Bvec;
        Info->Bmin = Lgm_Magnitude( &Bvec );

        Info->Ellipsoid_Footprint_Pn = *v2;
        Info->Bfield( v2, &Bvec, Info );
        Info->Ellipsoid_Footprint_Bvecn = Bvec;
        Info->Ellipsoid_Footprint_Bn    = Lgm_Magnitude( &Bvec );

        Info->Ellipsoid_Footprint_Ps = *v1;
        Info->Bfield( v1, &Bvec, Info );
        Info->Ellipsoid_Footprint_Bvecs = Bvec;
        Info->Ellipsoid_Footprint_Bs    = Lgm_Magnitude( &Bvec );

    } else if ( flag1 ) {

        Info->Ellipsoid_Footprint_Ps = *v1;
        Info->Bfield( v1, &Bvec, Info );
        Info->Ellipsoid_Footprint_Bvecs = Bvec;
        Info->Ellipsoid_Footprint_Bs    = Lgm_Magnitude( &Bvec );

	    Lgm_Convert_Coords( v1, &w, GSM_TO_SM, Info->c );
	    return( (w.z > 0.0) ? LGM_OPEN_N_LOBE : LGM_OPEN_S_LOBE );


    } else if ( flag2 ) {

        Info->Ellipsoid_Footprint_Pn = *v2;
        Info->Bfield( v2, &Bvec, Info );
        Info->Ellipsoid_Footprint_Bvecn = Bvec;
        Info->Ellipsoid_Footprint_Bn    = Lgm_Magnitude( &Bvec );

	    Lgm_Convert_Coords( v2, &w, GSM_TO_SM, Info->c );
	    return( (w.z > 0.0) ? LGM_OPEN_N_LOBE : LGM_OPEN_S_LOBE );

    } else if ( InitiallyBelowTargetHeight ) {

	    Lgm_Convert_Coords( u, &w, GSM_TO_SM, Info->c );
	    return( (w.z > 0.0) ? LGM_OPEN_N_LOBE : LGM_OPEN_S_LOBE );

    } else {

	    /*
	     *  IMF
 	     */
	    return( LGM_OPEN_IMF );

    }




    /*
     *  Finally, if we need to, lets determine the curvature of |B| at Bmin.
     * 
     *  Why would you need this? Well, the short answer is that we dont always
     *  need it. But if we ever want to compute the Sb integral for nearly
     *  equatorially mirroring particles and we do it with an interpolated
     *  array of B(s), the resulting integrals can be very inaccurate. On the
     *  other hand, for nearly equatorially mirroring particles, Sb can be
     *  obtained simply as an expansion that involves both Bmin and d^2B/ds^2.
     * 
     *  Thus, compute d^2B/ds^2 if we need it later.
     *  Use a 7-point formula;
     *
     *        2    (2)      1  (                                                          )
     *      h    F      =  --- (  2F   - 27F  +  270F  - 490F + 270F    - 27F    +  2F    )
     *             0       180 (    3       2        1       0      -1       -2        -3 )
     *      
     *                            F[0]    F[1]     F[2]    F[3]   F[4]      F[5]    F[6]   
     */


    h      = 0.01;     // Re
    h2_inv = 10000.0;  // Re^(-2)
    F[3] = Info->Bmin;
    u_scale.x =  100.0;  u_scale.y = 100.0; u_scale.z = 100.0; reset = TRUE;
    for (i=-3; i<=3; i++){
        if (i != 0) {
            s = i*h;
            P = Info->Pmin;
            Lgm_MagStep( &P, &u_scale, s, &Hdid, &Hnext, 1.0e-7, -1.0, &s, &reset, Info->Bfield, Info );
            Info->Bfield( &P, &Bvec, Info );
            F[i+3] = Lgm_Magnitude( &Bvec );
            //printf("F[%d] = %g Hdid = %g\n", i, F[i+3], Hdid);
        }
    }
    Info->d2B_ds2 = h2_inv/180.0 * ( 2.0*F[0] - 27.0*F[1] + 270.0*F[2] - 490.0*F[3] + 270.0*F[4] - 27.0*F[5] + 2.0*F[6] );
    //printf("d2B_ds2 = %g\n", d2B_ds2);

    Info->Sb0 = M_PI*M_SQRT2*sqrt( Info->Bmin/Info->d2B_ds2 );
    //printf("Sb0 = %g\n", Info->Sb0);
    



    return( LGM_CLOSED );



}


/*
 *   $Id: Lgm_Trace.c 143 2011-02-09 22:09:13Z mgh $
 */
