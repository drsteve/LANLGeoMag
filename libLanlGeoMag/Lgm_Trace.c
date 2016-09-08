/*! \file Lgm_Trace.c
 *
 *  \brief Convenience routine to find north and south footpoints and Bmin location in magnetopshere.
 *
 *
 *
 *  \author M.G. Henderson
 *  \date   1999
 *
 *
 *
 */


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


/**
 *  \brief
 *      This routine attempts to trace to the Earth (from the input position) in both
 *      directions along the field line.  If the field line is closed (meaning
 *      that we hit the Earth in both directions), then the routine will also
 *      trace to the minimum-B point.
 *  
 *  \detail
 *      This is a convenience routine that calls the routines Lgm_TraceToEarth()
 *      and Lgm_TraceToMinBSurf() (and performs some aditional calculations also).
 *      Note that the routine Lgm_TraceToEarth() traces to a height above the WGS84
 *      ellipsoid. Thus the target height here ( the input parameter Height ) is
 *      the geodetic height above the WGS84 ellipsoid. If you want to trace to a
 *      spherical representation of the Earth, use Lgm_TraceToSphericalEarth()
 *      instead.
 *  
 *      \param[in]       u     Input position vector in GSM coordinates.
 *      \param[out]     v1     Southern footpoint (where field line crosses the given geodetic height in the south) in GSM coordinates.
 *      \param[out]     v2     Northern footpoint (where field line crosses the given geodetic height in the north) in GSM coordinates.
 *      \param[out]     v3     Minimum B point along the field line in GSM coordinates.
 *      \param[in]      Height The altitude (in km) above the WGS84 ellispoid (i.e. geodetic height) used to define what we mean by the footpoint altitude.
 *      \param[in]      TOL1   Tolerance for tracing to the Earth.
 *      \param[in]      TOL2   Tolerance for tracing to the Min-B point.
 *      \param[in,out]  Info   Properly initialized/configured Lgm_MagModelInfo structure.
 *  
 *      \return 
 *          - LGM_OPEN_IMF     ( = 0 ) if field line is open at both ends (i.e. an IMF FL). In this case, nonem of the v1, v2, v3 values are valid.
 *          - LGM_CLOSED       ( = 1 ) if field line is closed. In this case all of the v1, v2, v3 are valid.
 *          - LGM_OPEN_N_LOBE  ( = 2 ) if field line is northern lobe FL (v2 valid).
 *          - LGM_OPEN_S_LOBE  ( = 3 ) if field line is southern lobe FL (v1 valid).
 *          - LGM_INSIDE_EARTH ( = -1 ) initial point is inside Earth (no points valid).
 *          - LGM_TARGET_HEIGHT_UNREACHABLE ( = -2 ) field line never got above the target height (no points valid).
 *          - LGM_BAD_TRACE    ( = -3 ) Lgm_MagStep() was unable to make a non-zero step (B-field zero?).
 *
 *
 *  Upon return, the following elements of the Info structure will be set or changed;
 *      - Hmax (maximum step size used for tracing. Not really useful to the user.)
 *      - Pmin. Position of Bmin along the field line. In Re relative to GSM coord system. (Same as v3.)
 *      - Smin. Distance from southern footpoint to Pmin along the FL. (in Re).
 *      - Stotal. Distance from southern footpoint to northern footpoint along the FL. (in Re).
 *      - Bvecmin. The GSM components of the magnetic field vector at Pmin. Units are nT.
 *      - Bmin. The magnitude of Bvecmin. Units of nT.
 *      - Ellipsoid_Footprint_Pn. GSM position of northern footpoint. (Same as v2.)
 *      - Ellipsoid_Footprint_Bvecn. GSM components of the magnetic field vector at Ellipsoid_Footprint_Pn. Units of nT.
 *      - Ellipsoid_Footprint_Bn. Magnitude of Ellipsoid_Footprint_Bvecn. Units of nT.
 *      - Ellipsoid_Footprint_Ps. GSM position of southern footpoint. (Same as v1.)
 *      - Ellipsoid_Footprint_Bvecs. GSM components of the magnetic field vector at Ellipsoid_Footprint_Ps. Units of nT.
 *      - Ellipsoid_Footprint_Bs. Magnitude of Ellipsoid_Footprint_Bvecs. Units of nT.
 *      - d2B_ds2. Second derivative of B wrt s at Pmin. Units of nT/Re/Re.
 *      - Sb0. Approximation for the Sb integral for nearly equatorially mirroring particles. Units of Re. See eqn 2.13b in Roederer.
 *
 *
 *  \author         M. Henderson
 *  \date           1999-2011
 *
 *
 *
 *  Changes:
 *      - April 1, 2011
 *          Added doxygen documenation.
 *      - June 13, 2008
 *          Added TOL1 and TOL2 args. They were hard set at 0.01 and 1e-7
 *      - October 14, 2008
 *          Added Height argument.  (Height is height in km above Earth to stop at)
 *
 *  Original version quite old (early 90's?)
 *  
 */
int Lgm_Trace( Lgm_Vector *u, Lgm_Vector *v1, Lgm_Vector *v2, Lgm_Vector *v3, double Height, double TOL1, double TOL2, Lgm_MagModelInfo *Info ) {

    int		    i, reset, flag1, flag2, InitiallyBelowTargetHeight, done;
    double	    sgn=1.0, R, Rtarget, Rinitial, Rplus, H, Hinitial;
    Lgm_Vector	w, Bvec;
    double      h, h_inv, h2_inv, F[7], Px[7], Py[7], Pz[7], s, Hdid, Hnext, Htry;
    double      d2Px_ds2, d2Py_ds2, d2Pz_ds2;
    double      GeodLat, GeodLong, GeodHeight;
    Lgm_Vector  u_scale, P, gpp;



    /*
     * Determine our initial geocentric radius in km. (u is assumed to be in
     * units of Re where we define Re to be WGS84_A.)
     */
    Rinitial = WGS84_A*Lgm_Magnitude( u ); // km

//Info->Bfield( u, &Bvec, Info );
//printf("u = %g %g %g B = %.15lf\n", u->x, u->y, u->z, Lgm_Magnitude(&Bvec));

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
    //InitiallyBelowTargetHeight = ( Rinitial <= Rtarget ) ? TRUE : FALSE; 

    Lgm_WGS84_to_GeodHeight( u, &Hinitial );
    InitiallyBelowTargetHeight = ( Hinitial <= Height ) ? TRUE : FALSE;
    

    // determine which direction will take us above the target R
    if (InitiallyBelowTargetHeight){
        Htry = 0.01; // step finely
        reset = TRUE;
        P = *u;
        if ( Lgm_MagStep( &P, &u_scale, Htry, &Hdid, &Hnext, 1.0, &s, &reset, Info->Bfield, Info ) < 0 ) return(LGM_BAD_TRACE);
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
            if ( Lgm_MagStep( &P, &u_scale, Htry, &Hdid, &Hnext, sgn, &s, &reset, Info->Bfield, Info ) < 0) return(LGM_BAD_TRACE);
            Lgm_WGS84_to_GeodHeight( &P, &H );
            if (H > Height) {
                done = TRUE;
            } else if ( H <  Hinitial ) {
                done = TRUE;
                return( -1 );
            }
            
        }

    }
    


    Info->Smin   = LGM_FILL_VALUE;
    Info->Snorth = LGM_FILL_VALUE;
    Info->Ssouth = LGM_FILL_VALUE;
    Info->Bmin   = LGM_FILL_VALUE;





    /*
     * Try to get to the Earth by tracing along the field in either direction.
     */
    sgn = -1.0;
    Info->Hmax = 10.0;
Info->Hmax = 0.50;
Info->Hmax = 0.10;

    
    flag2 = Lgm_TraceToEarth(  u, v2, Height, -sgn, TOL1, Info );
    Info->Snorth = Info->Trace_s;     // save distance from u to northern footpoint location.
    Info->v2_final = *v2;


    flag1 = Lgm_TraceToEarth(  u, v1, Height,  sgn, TOL1, Info );
    Info->Ssouth = Info->Trace_s;     // save distance from u to southern footpoint location.
    Info->v1_final = *v1;

    Info->Stotal = LGM_FILL_VALUE;
    Info->Smin   = LGM_FILL_VALUE;
    

//printf("HERE ************************   flag1 = %d flag2 = %d\n", flag1, flag2);


    if ( flag1 && flag2 ) {

	    /*
	     *  Closed FL -- attempt to trace to Eq Plane.
	     */
        //Lgm_TraceToMinBSurf( v1, v3, TOL1, TOL2, Info );
        //Lgm_TraceToMinBSurf( v1, v3, 0.1, TOL2, Info );
        Lgm_TraceToMinBSurf( u, v3, 0.1, TOL2, Info );
        Info->v3_final = *v3;
        Info->Pmin = *v3;
        //Info->Smin = Info->Trace_s;     // save location of Bmin. NOTE:  Smin is measured from the southern footpoint.
        Info->Bfield( v3, &Bvec, Info );
        Info->Bvecmin = Bvec;
        Info->Bmin = Lgm_Magnitude( &Bvec );
        //printf("Bmin = %.15lf\n",  Info->Bmin );

        /*
         * Various FL arc lengths...
         *  Snorth - distance from S/C to northern footpoint (set above)
         *  Ssouth - distance from S/C to southern footpoint (set above)
         *  Stotal - Total FL length ( Snorth + Ssouth ) (only compute for closed FL)
         *  Smin   - distance from southern footpoint to S/C (Ssouth - what we
         *           got from Lgm_TraceToMinBSurf() because we started at S/C)
         */
        Info->Stotal = Info->Snorth + Info->Ssouth; // Total FL length
        Info->Smin   = Info->Ssouth - Info->Trace_s;  // length from south foot to S/C
        Info->Trace_s = Info->Stotal;
//printf("Info->Ssouth, Info->Trace_s = %g %g\n", Info->Ssouth, Info->Trace_s);
        

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

        /*
         * It is not a closed FL, but it may still have a min-B
         * Try to find it.
         */
        if ( Lgm_TraceToMinBSurf( v1, v3, 0.1, TOL2, Info ) ) {
            Info->v3_final = *v3;
            Info->Pmin = *v3;
            Info->Bfield( v3, &Bvec, Info );
            Info->Bvecmin = Bvec;
            Info->Bmin = Lgm_Magnitude( &Bvec );
        } else {
            Info->v3_final = *v3;
            v3->x = v3->y = v3->z = -1e31;
        }

	    return( (w.z > 0.0) ? LGM_OPEN_N_LOBE : LGM_OPEN_S_LOBE );


    } else if ( flag2 ) {

        Info->Ellipsoid_Footprint_Pn = *v2;
        Info->Bfield( v2, &Bvec, Info );
        Info->Ellipsoid_Footprint_Bvecn = Bvec;
        Info->Ellipsoid_Footprint_Bn    = Lgm_Magnitude( &Bvec );

	    Lgm_Convert_Coords( v2, &w, GSM_TO_SM, Info->c );

        /*
         * It is not a closed FL, but it may stilkl have a min-B
         * Try to find it.
         */
        if ( Lgm_TraceToMinBSurf( v2, v3, 0.1, TOL2, Info ) ) {
            Info->v3_final = *v3;
            Info->Pmin = *v3;
            Info->Bfield( v3, &Bvec, Info );
            Info->Bvecmin = Bvec;
            Info->Bmin = Lgm_Magnitude( &Bvec );
        } else {
            Info->v3_final = *v3;
            v3->x = v3->y = v3->z = -1e31;
        }

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
     *  Finally, if we need to, lets determine the d^2|B|/ds^2 at Bmin.
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

    if ( Info->ComputeSb0 ) {

        h      = 0.01;     // Re
        h_inv  = 1000.0;   // Re^(-1)
        h2_inv = 10000.0;  // Re^(-2)
        F[3]  = Info->Bmin;
        Px[3] = Info->Pmin.x;
        Py[3] = Info->Pmin.y;
        Pz[3] = Info->Pmin.z;
        u_scale.x =  u_scale.y = u_scale.z = 1.0; reset = TRUE;
        for (i=-3; i<=3; i++){
            if (i != 0) {
                s = i*h;
                P = Info->Pmin;
                if ( Lgm_MagStep( &P, &u_scale, s, &Hdid, &Hnext, -1.0, &s, &reset, Info->Bfield, Info ) < 0 ) return(LGM_BAD_TRACE);
                Info->Bfield( &P, &Bvec, Info );
                F[i+3] = Lgm_Magnitude( &Bvec );
                //printf("F[%d] = %g Hdid = %g\n", i, F[i+3], Hdid);

                Px[i+3] = P.x;
                Py[i+3] = P.y;
                Pz[i+3] = P.z;


            }
        }
        Info->d2B_ds2 = h2_inv/180.0 * ( 2.0*F[0] - 27.0*F[1] + 270.0*F[2] - 490.0*F[3] + 270.0*F[4] - 27.0*F[5] + 2.0*F[6] );
        //printf("d2B_ds2 = %g\n", d2B_ds2);

        Info->Sb0 = M_PI*M_SQRT2*sqrt( Info->Bmin/Info->d2B_ds2 );
        //printf("Sb0 = %g\n", Info->Sb0);


        // 2nd deriv. page 450 CRC standard math tables.
        d2Px_ds2 = h2_inv/180.0 * ( 2.0*Px[0] - 27.0*Px[1] + 270.0*Px[2] - 490.0*Px[3] + 270.0*Px[4] - 27.0*Px[5] + 2.0*Px[6] );
        d2Py_ds2 = h2_inv/180.0 * ( 2.0*Py[0] - 27.0*Py[1] + 270.0*Py[2] - 490.0*Py[3] + 270.0*Py[4] - 27.0*Py[5] + 2.0*Py[6] );
        d2Pz_ds2 = h2_inv/180.0 * ( 2.0*Pz[0] - 27.0*Pz[1] + 270.0*Pz[2] - 490.0*Pz[3] + 270.0*Pz[4] - 27.0*Pz[5] + 2.0*Pz[6] );

        /*
         * Construct g" = (d2Px_ds2, d2Py_ds2, d2Pz_ds2).
         *
         *  Then K = |g"|  and R = 1/K
         */
        gpp.x = d2Px_ds2; gpp.y = d2Py_ds2; gpp.z = d2Pz_ds2;
        Info->Kappa = Lgm_Magnitude( &gpp );
        Info->RofC  = 1.0/Info->Kappa;
    
    }




    return( LGM_CLOSED );



}
