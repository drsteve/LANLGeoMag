/* TraceLine, Copyright (c) 2007 Michael G. Henderson <mghenderson@lanl.gov>
 *
 *  This assumes a spherical Earth. Use Lgm_TraceToEarth() if you need to trace
 *  to a height aboveb the ellipsoid.
 *
 *    - This routine is intended to trace a field line and save the points in
 *    an array. The idea is that the user can use the adaptive tracing routines
 *    to identify critical points on a FL, (and/or if its open, etc.). Here we
 *    just assume that the user wants to trace out the line and save the points
 *    for later use (for example, to plot the FL or to be used to compute
 *    integrals like I, Sb, etc., along the FL).
 *
 * Some Cautions:
 *
 *    (1) We need to be careful when using the arrays to interpolate though!
 *    The strategy here is to trace with equal intervals in s, and then to tack
 *    on the last point of the trace. Thus, the spacing between the last point
 *    and the previous to last one will have an odd spacing. You need to keep
 *    this in mind when doing interpolation.
 *
 *
 *    (2) Because of the equal spacing, various critical points wont be
 *    captured very precisely.  For example, the location and value of Bmin
 *    will in general fall between grid points. In order to ensure that this
 *    point is representable by a linear interpolation, we need to add it into
 *    the interpolation scheme when we interpolate over the interval that
 *    contains the point.
 *
 * Permission to use, copy, modify, distribute, and sell this software and its
 * documentation for any purpose is hereby granted without fee, provided that
 * the above copyright notice appear in all copies and that both that copyright
 * notice and this permission notice appear in supporting documentation.  No
 * representations are made about the suitability of this software for any
 * purpose.  It is provided "as is" without express or implied warranty.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include "Lgm/Lgm_MagModelInfo.h"


/*
 *  Set the interpolation type to use. (Note, the searches in gsl are not as
 *  optimal as they could be. Probably losing some efficiency there. A
 *  taylor-made interpolator may do better?)
 */
//#define GSL_INTERP  gsl_interp_linear
//#define GSL_INTERP  gsl_interp_akima
#define GSL_INTERP  gsl_interp_cspline




int Lgm_TraceLine( Lgm_Vector *u, Lgm_Vector *v, double H0, double sgn, double tol, int AddBminPoint, Lgm_MagModelInfo *Info ) {

    Lgm_Vector	u_scale;
    double	    Htry, Hdid, Hnext, Hmin, Hmax, s, ss;
    double	    Sa=0.0, Sc=0.0, d;
    double	    R0, R, Fa, Fb, Fc, F;
    double	    Ra, Rb, Rc;
    Lgm_Vector	Pa, Pc, P, Bvec, Bcdip;
    int		    done, reset, n, SavePnt;




    Pa.x = Pa.y = Pa.z = 0.0;
    Pc.x = Pc.y = Pc.z = 0.0;
    P.x  = P.y  = P.z  = 0.0;


    /*
     *  H0 is in km above Earth's surface (assumed to be spherical here).
     *  Convert to geocentric radius.
     */
    R0 = H0/Re + 1.0;



    Htry = Info->Hmax;  // we want to step with constant increments.
    Hmin = 0.0001;      // This may be necessary to find the endpoint.
    Hmax = Info->Hmax;  // Dont use step bigger than this.
    u_scale.x =  10.0;  u_scale.y = 1.0; u_scale.z = 10.0;
    R = Ra = Rb = Rc = 0.0;
    F = Fa = Fb = Fc = 0.0;

    /*
     *  Save first point
     */
    n = 0; Info->nPnts = n;
    ss = 0.0;
    Info->Bfield( u, &Bvec, Info );
    Info->s[n]    = ss;                         // save arc length
    Info->Px[n]   = u->x;                       // save 3D position vector.
    Info->Py[n]   = u->y;                       //
    Info->Pz[n]   = u->z;                       //
    Info->Bvec[n] = Bvec;                       // save 3D B-field vector.
    Info->Bmag[n] = Lgm_Magnitude( &Bvec );     // save field strength (and increment counter)
    Lgm_B_cdip( u, &Bcdip, Info );
    Info->BminusBcdip[n] = Info->Bmag[n] - Lgm_Magnitude( &Bcdip );     // save field strength (and increment counter)
    ++n;
    if (n > LGM_MAX_INTERP_PNTS){
	    print("Warning: n > LGM_MAX_INTERP_PNTS (%d)\n", LGM_MAX_INTERP_PNTS);
    }



    /*
     *   Add code to check to see if the start point is already below the target height.
     *   If it is, we can trace up until we are beyond it or some such thing...
     */
    done  = FALSE;
    reset = TRUE;
    P = *u;
    R = Lgm_Magnitude( &P );
    F = R - R0;
    if (F < 0.0 ){
        // we are already below target height.
        // Trace until we are not.
        Htry  = 0.1;
        while ( !done ) {
            if ( Lgm_MagStep( &P, &u_scale, Htry, &Hdid, &Hnext, 1.0e-7, sgn, &s, &reset, Info->Bfield, Info ) < 0 ) { printf("BAILING 1\n"); return(-1);}
            ss += Hdid;  // Note that we should trap conditions where Hdid != Htry since this will be a problem...
            R = Lgm_Magnitude( &P );
            F = R - R0;
            if ( F > 0.0 ){
                done = TRUE;
            }
        }

        /*
         * Save this (new) point only if its different from the first one)
         */
        if ( ss > Info->s[n-1] ) {
            Info->Bfield( &P, &Bvec, Info );
            Info->s[n]    = ss;                         // save arc length
            Info->Px[n]   = P.x;                        // save 3D position vector.
            Info->Py[n]   = P.y;                        //
            Info->Pz[n]   = P.z;                        //
            Info->Bvec[n] = Bvec;                       // save 3D B-field vector.
            Info->Bmag[n] = Lgm_Magnitude( &Bvec );     // save field strength (and increment counter)
            Lgm_B_cdip( &P, &Bcdip, Info ); 
            Info->BminusBcdip[n] = Info->Bmag[n] - Lgm_Magnitude( &Bcdip );     // save field strength (and increment counter)
            ++n;
        }

    }






    Pa   = P;
    Ra   = Lgm_Magnitude( &Pa );
    Sa   = 0.0;

    /*
     *  Keep stepping along the field line until we drop below the target
     *  radius, R0.  (Or bail if its open). This completes the endpoints of the
     *  bracket pair. We get these points first, because there are field lines
     *  that have more than one local minimum in fabs(R-R0). So find Pa and Pc
     *  first and then complete the triple later. CHECK on this. I changed this
     *  from a minima search to a simple root finder. Is it OK still?
     */
    P       = Pa;
    done    = FALSE;
    reset   = TRUE;
    SavePnt = TRUE;
    while ( !done ) {

        if ( Lgm_MagStep( &P, &u_scale, Htry, &Hdid, &Hnext, 1.0e-7, sgn, &s, &reset, Info->Bfield, Info ) < 0 ) { printf("BAILING 2\n"); return(-1);}


        R = Lgm_Magnitude( &P );
	    F =  R - R0;

	    if (   (P.x > Info->OpenLimit_xMax) || (P.x < Info->OpenLimit_xMin) || (P.y > Info->OpenLimit_yMax) || (P.y < Info->OpenLimit_yMin) 
	        || (P.z > Info->OpenLimit_zMax) || (P.z < Info->OpenLimit_zMin) ) {
	        /*
	         *  Open FL!
	         */
	        v->x = v->y = v->z = 0.0;
	        return(0);
	    } else if ( (F < 0.0)&&(ss>0.01) ){
	        done = TRUE;
	        Pc = P;
	        Rc = R;
	        Fc = F;
	        Sc = Sa + Hdid;
	    } else {
	        Pa = P;
	        Fa = F;
	        Ra = R;
	        Sa = 0.0;
            ss += Hdid;  
            /*
             * Compute field strength and save the results in the array.
             */
            if ( SavePnt && (ss > Info->s[n-1]) ) {
                Info->Bfield( &P, &Bvec, Info );
                Info->s[n]    = ss;                         // save arc length
                Info->Px[n]   = P.x;                        // save 3D position vector.
                Info->Py[n]   = P.y;                        //
                Info->Pz[n]   = P.z;                        //
                Info->Bvec[n] = Bvec;                        // save 3D B-field vector.
                Info->Bmag[n] = Lgm_Magnitude( &Bvec );     // save field strength (and increment counter)
                Lgm_B_cdip( &P, &Bcdip, Info ); 
                Info->BminusBcdip[n] = Info->Bmag[n] - Lgm_Magnitude( &Bcdip );     // save field strength (and increment counter)
                ++n;
            }
	    }



	    /*
         *  We have been stepping along at a constant Htry. But when we get
         *  close to the end we do not want to drop below the surface of the
         *  Earth -- some models crash and burn if you do that. To make sure we
         *  bracket the minimum, we need to reduce Htry if F is less than our
         *  current altitude above the Earth.  Note that F may well be < 0, in
         *  which case we already have a bracket and we will bail out of this
         *  loop anyway.
         *
         *  If we have to reduce Htry, then we also want to stop saving the
         *  points.  The final point will be the actual end of the FL. (I.e. we
         *  dont want all of these smaller steps cluttering up our array at the
         *  end.)
	     */
        if ( Htry > (R-1.0) ) {
            Htry = 0.5*(R-1.0);
            SavePnt = FALSE;
         } else {
            SavePnt = TRUE;
	     }

    }

    /*
     *  We have a bracket. Now go in for the kill. This just finds where mthe
     *  FL intersects the traget altitude. And this point will become the final
     *  value in our saved arrays.
     *
     *  Use golden-section search to converge toward minimum.  (Sa, Sb, Sc) are
     *  the distances of the triple points along the FL. For convenience, we
     *  maintain Sa = 0.0, so that Sb and Sc are always the distances relative
     *  to Pa. (And so Sc is always the width of the bracketed interval. So
     *  when Sc gets small enough we bail out and take Pb as the min).
     *
     */
    done  = FALSE;
    reset = TRUE;
    while (!done) {

	    d = Sc - Sa;
	    if ( fabs(d) < tol ) {
	        done = TRUE;
	    } else {

            P = Pa; Htry = 0.5*d;
            if ( Lgm_MagStep( &P, &u_scale, Htry, &Hdid, &Hnext, 1.0e-7, sgn, &s, &reset, Info->Bfield, Info ) < 0 ) { printf("BAILING 3\n");return(-1);}
            R = Lgm_Magnitude( &P );
            F =  R - R0;
            if ( F >= 0.0 ) {
                Pa = P; Fa = F; Sa += Hdid;
            } else {
                Pc = P; Fc = F; Sc = Sa + Hdid;
            }
	    }
    }


    /*
     *  Compute a final value based on what we have. One possibility is to
     *  average Pa and Pc.  This is what we used for a long time, but it may
     *  not do well when it comes to interpolating near the endpoints. Another
     *  way is the take the point we know is beyond the mirror point altitude.
     *  That would (probably?) ensure the real mirror point value is contained
     *  in our interp array.
     */
    //v->x = 0.5*(Pa.x + Pc.x); v->y = 0.5*(Pa.y + Pc.y); v->z = 0.5*(Pa.z + Pc.z);
    *v = Pc; ss += Sc;


    /*
     *  Save final point.
     */
    if ( ss > Info->s[n-1] ) {
        Info->Bfield( v, &Bvec, Info );
        Info->s[n]    = ss;                         // save arc length
        Info->Px[n]   = v->x;                       // save 3D position vector.
        Info->Py[n]   = v->y;                       //
        Info->Pz[n]   = v->z;                       //
        Info->Bvec[n] = Bvec;                       // save 3D B-field vector.
        Info->Bmag[n] = Lgm_Magnitude( &Bvec );     // save field strength (and increment counter)
        Lgm_B_cdip( v, &Bcdip, Info ); 
        Info->BminusBcdip[n] = Info->Bmag[n] - Lgm_Magnitude( &Bcdip );     // save field strength (and increment counter)
        ++n;
    }
    Info->nPnts     = n;                       // set total number of points in the array.
    Info->ds        = Info->Hmax;              // spacing in s for the array -- will help to know this
                                               // when trying to interpolate.







    /*
     *  Add the Smin, Bmin point. Only do this if AddBminPoint is TRUE
     *  This will only make sense if these values are legitimate for this FL.
     *  Perhaps it would be better to force user to do this elesewhere.
     */
    if ( AddBminPoint ) {
        printf("1) ADDING NEW POINT\n");
        // MUST ADD Bcdip for this too!
        AddNewPoint( Info->Smin, Info->Bmin, &Info->Pmin, Info );
    }
    //printf("1) F >>>>>>>>> Info->nPnts = %d <<<<<<<<<<\n", Info->nPnts);



    return( 1 );


}


/*
 *  This version is almost the same. But here we require that we trace along the FL
 *  a minimum distance before trying to do the convergence on radius.
 */
int Lgm_TraceLine2( Lgm_Vector *u, Lgm_Vector *v, double H0, double MinDist, double sgn, double tol, int AddBminPoint, Lgm_MagModelInfo *Info ) {

    Lgm_Vector	u_scale;
    double	    Htry, Hdid, Hnext, Hmin, Hmax, s, ss;
    double	    Sa=0.0, Sc=0.0, d;
    double	    R0, R, Fa, Fb, Fc, F;
    double	    Ra, Rb, Rc;
    Lgm_Vector	Pa, Pc, P, Bvec, Bcdip;
    int		    done, reset, n, SavePnt;

//printf("\n\n\n*************************\n");
//printf("*************************\n");
//printf("*************************\n");

    if (MinDist < 1e-5) MinDist = 0.0;


    Pa.x = Pa.y = Pa.z = 0.0;
    Pc.x = Pc.y = Pc.z = 0.0;
    P.x  = P.y  = P.z  = 0.0;


    /*
     *  H0 is in km above Earth's surface. Convert to geocentric
     *  radius.
     */
    R0 = H0/Re + 1.0;



    Htry = Info->Hmax; // we want to step with constant increments.
    Hmin = 0.0001;      // This may be necessary to find the endpoint.
    Hmin = 1e-7;
    Hmax = Info->Hmax; // Dont use step bigger than this.
    u_scale.x =  100.0;  u_scale.y = 100.0; u_scale.z = 100.0;
    R = Ra = Rb = Rc = 0.0;
    F = Fa = Fb = Fc = 0.0;


    Pa   = *u;
    Ra   = Lgm_Magnitude( &Pa );
    Sa   = 0.0;



    /*
     * Compute field strength at starting point and save the results in the
     * array.
     */
    n = 0; Info->nPnts = n;
    ss = 0.0;
    Info->Bfield( &Pa, &Bvec, Info );
    Info->s[n]    = ss;                         // save arc length
    Info->Px[n]   = Pa.x;                       // save 3D position vector.
    Info->Py[n]   = Pa.y;                       //
    Info->Pz[n]   = Pa.z;                       //
    Info->Bvec[n] = Bvec;                       // save 3D B-field vector.
    Info->Bmag[n] = Lgm_Magnitude( &Bvec );     // save field strength (and increment counter)
    Lgm_B_cdip( &Pa, &Bcdip, Info );
    Info->BminusBcdip[n] = Info->Bmag[n] - Lgm_Magnitude( &Bcdip );     // save field strength (and increment counter)
    ++n;

    if (n > LGM_MAX_INTERP_PNTS){
	    print("Warning: n > LGM_MAX_INTERP_PNTS (%d)\n", LGM_MAX_INTERP_PNTS);
    }


    /*
     *  Keep stepping along the field line until we drop below the target
     *  radius, R0.  (Or bail if its open). This completes the endpoints of the
     *  bracket pair. We get these points first, because there are field lines
     *  that have more than one local minimum in fabs(R-R0). So find Pa and Pc
     *  first and then complete the triple later. CHECK on this. I changed this
     *  from a minima search to a simple root finder. Is it OK still?
     */
    P       = Pa;
    done    = FALSE;
    reset   = TRUE;
    SavePnt = TRUE;
    while ( !done ) {

        if ( Lgm_MagStep( &P, &u_scale, Htry, &Hdid, &Hnext, 1.0e-7, sgn, &s, &reset, Info->Bfield, Info ) < 0 ) { printf("BAILING 4\n");return(-1);}


        R = Lgm_Magnitude( &P );
	    F =  R - R0;
//printf( "P = %g %g %g    R, R0, F, Htry  = %g %g %g %g    ss, MinDist = %g %g\n", P.x, P.y, P.z, R, R0, F, Htry, ss, MinDist);

	    if (   (P.x > Info->OpenLimit_xMax) || (P.x < Info->OpenLimit_xMin) || (P.y > Info->OpenLimit_yMax) || (P.y < Info->OpenLimit_yMin) 
	        || (P.z > Info->OpenLimit_zMax) || (P.z < Info->OpenLimit_zMin) ) {
	        /*
	         *  Open FL!
	         */
	        v->x = v->y = v->z = 0.0;
	        return(0);
	    } else if ( (ss < MinDist) || ( F > 0.0 ) ) { // keep on stepping
            Pa = P;
            Fa = F;
            Ra = R;
            Sa = 0.0;
            ss += Hdid;
            /*
             * Compute field strength and save the results in the array.
             */
            if ( SavePnt && (ss > Info->s[n-1]) ) {
                Info->Bfield( &P, &Bvec, Info );
                Info->s[n]    = ss;                         // save arc length
                Info->Px[n]   = P.x;                        // save 3D position vector.
                Info->Py[n]   = P.y;                        //
                Info->Pz[n]   = P.z;                        //
                Info->Bvec[n] = Bvec;                       // save 3D B-field vector.
                Info->Bmag[n] = Lgm_Magnitude( &Bvec );     // save field strength (and increment counter)
                Lgm_B_cdip( &P, &Bcdip, Info ); 
                Info->BminusBcdip[n] = Info->Bmag[n] - Lgm_Magnitude( &Bcdip );     // save field strength (and increment counter)
                ++n;
            }
        } else {
            done = TRUE;
            Pc = P;
            Rc = R;
            Fc = F;
            Sc = Sa + Hdid;
        } 
	/*
         *  We have been stepping along at a constant Htry. But when we get
         *  close to the end we do not want to drop below the surface of the
         *  Earth -- some models crash and burn if you do that. To make sure we
         *  bracket the minimum, we need to reduce Htry if F is less than our
         *  current altitude above the Earth.  Note that F may well be < 0, in
         *  which case we already have a bracket and we will bail out of this
         *  loop anyway.
         *
         *  If we have to reduce Htry, then we also want to stop saving the
         *  points.  The final point will be the actual end of the FL. (I.e. we
         *  dont want all of these smaller steps cluttering up our array at the
         *  end.)
	     */
        if ( Htry > (R-1.0) ) {
            Htry = 0.5*(R-1.0);
            SavePnt = FALSE;
        } else {
            SavePnt = TRUE;
        }

    }

    /*
     *  We have a bracket. Now go in for the kill. This just finds where mthe
     *  FL intersects the traget altitude. And this point will become the final
     *  value in our saved arrays.
     *
     *  Use golden-section search to converge toward minimum.  (Sa, Sb, Sc) are
     *  the distances of the triple points along the FL. For convenience, we
     *  maintain Sa = 0.0, so that Sb and Sc are always the distances relative
     *  to Pa. (And so Sc is always the width of the bracketed interval. So
     *  when Sc gets small enough we bail out and take Pb as the min).
     *
     */
    done  = FALSE;
    reset = TRUE;
    while (!done) {

	    d = Sc - Sa;
	    if ( fabs(d) < tol ) {
	        done = TRUE;
	    } else {

            P = Pa; Htry = 0.5*d;
            if ( Lgm_MagStep( &P, &u_scale, Htry, &Hdid, &Hnext, 1.0e-7, sgn, &s, &reset, Info->Bfield, Info ) < 0 ) { printf("BAILING 5\n");return(-1);}
            R = Lgm_Magnitude( &P );
            F =  R - R0;
            if ( F >= 0.0 ) {
                Pa = P; Fa = F; Sa += Hdid;
            } else {
                Pc = P; Fc = F; Sc = Sa + Hdid;
            }
	    }    
    }



    /*
     *  Compute a final value based on what we have. One possibility is to
     *  average Pa and Pc.  This is what we used for a long time, but it may
     *  not do well when it comes to interpolating near the endpoints. Another
     *  way is the take the point we know is beyond the mirror point altitude.
     *  That would (probably?) ensure the real mirror point value is contained
     *  in our interp array.
     */
    //v->x = 0.5*(Pa.x + Pc.x); v->y = 0.5*(Pa.y + Pc.y); v->z = 0.5*(Pa.z + Pc.z);
    *v = Pc; ss += Sc;


    /*
     *  Save final point.
     */
    if ( ss > Info->s[n-1] ) {
        Info->Bfield( v, &Bvec, Info );
        Info->s[n]    = ss;                         // save arc length
        Info->Px[n]   = v->x;                       // save 3D position vector.
        Info->Py[n]   = v->y;                       //
        Info->Pz[n]   = v->z;                       //
        Info->Bvec[n] = Bvec;                       // save 3D B-field vector.
        Info->Bmag[n] = Lgm_Magnitude( &Bvec );     // save field strength (and increment counter)
        Lgm_B_cdip( v, &Bcdip, Info ); 
        Info->BminusBcdip[n] = Info->Bmag[n] - Lgm_Magnitude( &Bcdip );     // save field strength (and increment counter)
        ++n;
    }
    Info->nPnts     = n;                       // set total number of points in the array.
    Info->ds        = Info->Hmax;              // spacing in s for the array -- will help to know this
                                               // when trying to interpolate.




    /*
     *  Add the Smin, Bmin point. Only do if AddBminPoint is TRUE.
     *  This will only make sense if these values are legitimate for this FL.
     *  Perhaps it would be better to force user to do this elesewhere.
     */
    if ( AddBminPoint ) {
        AddNewPoint( Info->Smin, Info->Bmin, &Info->Pmin, Info );
//printf("2) ADDING NEW POINT\n");
//printf("n, nPnts = %d %d\n", n, Info->nPnts);
    }


    return( 1 );


}

void ReplaceFirstPoint( double s, double B, Lgm_Vector *P, Lgm_MagModelInfo *Info ) {

    Lgm_Vector  Bcdip;

    Info->s[0]    = s;
    Info->Bmag[0] = B;
    Info->Px[0]   = P->x;
    Info->Py[0]   = P->y;
    Info->Pz[0]   = P->z;
    Lgm_B_cdip( P, &Bcdip, Info ); 
    Info->BminusBcdip[0] = Info->Bmag[0] - Lgm_Magnitude( &Bcdip );

}


void AddNewPoint( double s, double B, Lgm_Vector *P, Lgm_MagModelInfo *Info ) {



    int         i, i1=0, i2=0;
    Lgm_Vector  Bcdip;

    /*
     * Dont Add if we already have it
     */
    for (i=0; i<Info->nPnts; ++i){
        if (fabs(s-Info->s[i]) < 1e-10) {
//            printf("AddNewPoint: Warning, already have a point very close to this s value. (s,B) = (%g,%g) (Info->s[%d], Info->Bmag[%d]) = (%g, %g). Point not added.\n", 
//                s, B, i, i, Info->s[i], Info->Bmag[i]);
            return;
        }
    }


    /*
     *  This determines the indices that bound the (s, B)  point.
     */
    if ( s < Info->s[0]) {
        i1 = -1;
        i2 =  0;
    } else if ( s > Info->s[Info->nPnts-1]) {
        i1 = Info->nPnts-1;
        i2 = Info->nPnts;
    } else {

        /*
         *  Find where new point should be added
         */
        for (i=0; i<Info->nPnts; i++){
            i1 = i-1;
            i2 = i;
            if ( Info->s[i] > s ) {
                break;
            }
        }
        /*
         *  Make room for new point
         */
        for (i=Info->nPnts-1; i>=i2; i--){
            Info->s[i+1]    = Info->s[i];
            Info->Bmag[i+1] = Info->Bmag[i];
            Info->Px[i+1]   = Info->Px[i];
            Info->Py[i+1]   = Info->Py[i];
            Info->Pz[i+1]   = Info->Pz[i];
        }

    }


    /*
     *  Add new point. Note: we are not adding Bvec at this time so dont rely on using Bvec
     */
    Info->s[i2]    = s;
    Info->Bmag[i2] = B;
    Info->Px[i2]   = P->x;
    Info->Py[i2]   = P->y;
    Info->Pz[i2]   = P->z;
    Lgm_B_cdip( P, &Bcdip, Info ); 
    Info->BminusBcdip[i2] = Info->Bmag[i2] - Lgm_Magnitude( &Bcdip );
    ++Info->nPnts;


}



int InitSpline( Lgm_MagModelInfo *Info ) {



    int i;

    if ( ( Info->nPnts < 2 )||( Info->AllocedSplines ) ) return( 0 );

//    gsl_set_error_handler_off(); // Turn off gsl default error handler

    for (i=0; i<Info->nPnts-1; ++i){
        if (Info->s[i] >= Info->s[i+1]) {
            printf("InitSpline: Error? Info->s[%d] = %g    Info->s[%d] = %g\n", i, Info->s[i], i+1, Info->s[i+1]);
        }
    }
    /*
     *  Now that we have the array of points, lets initialize the gsl spline stuff
     *
     *  To evaluate the spline, use the following:
     *
     *      y = gsl_spline_eval( Info->spline, x, Info->acc );
     *
     *
     */
    Info->acc    = gsl_interp_accel_alloc( );
//    Info->accPx  = gsl_interp_accel_alloc( );
//    Info->accPy  = gsl_interp_accel_alloc( );
//    Info->accPz  = gsl_interp_accel_alloc( );
    if ( Info->nPnts > 2 ){
        Info->spline   = gsl_spline_alloc( GSL_INTERP, Info->nPnts );
//        Info->splinePx = gsl_spline_alloc( GSL_INTERP, Info->nPnts );
//        Info->splinePy = gsl_spline_alloc( GSL_INTERP, Info->nPnts );
//        Info->splinePz = gsl_spline_alloc( GSL_INTERP, Info->nPnts );
    } else {
        Info->spline   = gsl_spline_alloc( gsl_interp_linear, Info->nPnts );
//        Info->splinePx = gsl_spline_alloc( gsl_interp_linear, Info->nPnts );
//        Info->splinePy = gsl_spline_alloc( gsl_interp_linear, Info->nPnts );
//        Info->splinePz = gsl_spline_alloc( gsl_interp_linear, Info->nPnts );
    }
    gsl_spline_init( Info->spline, Info->s, Info->Bmag, Info->nPnts );
//    gsl_spline_init( Info->spline,   Info->s, Info->BminusBcdip, Info->nPnts );
//    gsl_spline_init( Info->splinePx, Info->s, Info->Px, Info->nPnts );
//    gsl_spline_init( Info->splinePy, Info->s, Info->Py, Info->nPnts );
//    gsl_spline_init( Info->splinePz, Info->s, Info->Pz, Info->nPnts );

    Info->AllocedSplines = TRUE;

    return( 1 );

}

int FreeSpline( Lgm_MagModelInfo *Info ) {

    if ( !Info->AllocedSplines ) return(0);

//    gsl_set_error_handler_off(); // Turn off gsl default error handler
    if ( Info->spline != NULL ) gsl_spline_free( Info->spline );
//    if ( Info->splinePx != NULL ) gsl_spline_free( Info->splinePx );
//    if ( Info->splinePy != NULL ) gsl_spline_free( Info->splinePy );
//    if ( Info->splinePz != NULL ) gsl_spline_free( Info->splinePz );
    if ( Info->acc != NULL ) gsl_interp_accel_free( Info->acc );
//    if ( Info->accPx != NULL ) gsl_interp_accel_free( Info->accPx );
//    if ( Info->accPy != NULL ) gsl_interp_accel_free( Info->accPy );
//    if ( Info->accPz != NULL ) gsl_interp_accel_free( Info->accPz );

    Info->AllocedSplines = FALSE;

    return(1);
}




/*
 * Routine to interpolate B along a FL that was traced by the routine above.
 * I.e. compute B(s) using the pre-traced arrays.
 */

double  BofS( double s, Lgm_MagModelInfo *Info ) {

    double      m, ds, B, BminusBcdip, Bcdip;
    int         i1, i2;
    Lgm_Vector  P, Bvec;

    /*
     * Use GSL to compute BminusBcdip(s)
     */
//    BminusBcdip = gsl_spline_eval( Info->spline, s, Info->acc );

    /*
     * Interpolate to get P(s) and Bcdip(P(s))
     */
//    P.x = gsl_spline_eval( Info->splinePx, s, Info->accPx );
//    P.y = gsl_spline_eval( Info->splinePy, s, Info->accPy );
//    P.z = gsl_spline_eval( Info->splinePz, s, Info->accPz );
//    Lgm_B_cdip( &P, &Bvec, Info );
//    Bcdip = Lgm_Magnitude( &Bvec );
//    B = BminusBcdip + Bcdip;



    B = gsl_spline_eval( Info->spline, s, Info->acc );
    return( B );


    // NOT REACHABLE. DECIDE WHAT WE WHAT TO DO HERE...


    /*
     * user needs to make sure s is <= the end point value
     */
    if ( s > Info->s[Info->nPnts-1]+1e-7 ) {
        return(-999.9);
    }


    /*
     * Try to guess which indices we need to interpolate between
     */
    i1 = (int)(s/Info->ds);
    i2 = i1+1;
    if ((i1 > Info->nPnts-1)||(i1 > Info->nPnts-1)){
        return(-9999.9);
    }

    if (i2 > Info->nPnts-1){
        i2 = Info->nPnts-1;
        i1 = i2-1;
    }


    /*
     * If Smin (location of Bmin) is in the interval [s(i1):s(i2)]
     * then we need to do more work...
     */
    if ( (i1==Info->imin1) && (i2==Info->imin2)){
        if ( s > Info->Smin ){
            ds = Info->s[i2] - Info->Smin;
            if ( ds > 1e-7) {
                m  = (Info->Bmag[i2] - Info->Bmin)/ds;
                B  = m*(s-Info->Smin) + Info->Bmin;
            } else {
                // handles possibles cases where Smin is extremely close to one
                // of the grid points.
                B = Info->Bmag[i2];
            }
        } else {
            ds = Info->Smin - Info->s[i1];
            if ( ds > 1e-7) {
                m  = (Info->Bmin - Info->Bmag[i1])/ds;
                B  = m*(s-Info->s[i1]) + Info->Bmag[i1];
            } else {
                B = Info->Bmag[i1];
            }
        }
    } else {
        ds = Info->s[i2] - Info->s[i1];             // run (takes care of last interval where ds != Info->ds)
        if (fabs(ds) > 1e-7){
            m  = (Info->Bmag[i2] - Info->Bmag[i1])/ds;  // slope
            B  = m*(s-Info->s[i1]) + Info->Bmag[i1];    // linearly interped value
        } else {    
            B  = Info->Bmag[i1];   
        }
    }


    return( B );


}




/*
 *   Compute the +/-s vals of the mirror points. Start at s(Bmin) point.
 */
int  SofBm( double Bm, double *ss, double *sn, Lgm_MagModelInfo *Info ) {

    int     done;
    double  s, s1, s2;
    double  B, B1, B2, d;



    // test to see that we have enough points
    if (Info->nPnts < 3 ) {
        return(FALSE);
    }

    // test to see that array brackets Bm
    if ( (Bm > Info->Bmag[0]) || (Bm > Info->Bmag[Info->nPnts-1]) ){
        return(FALSE);
    }


    s = Info->Smin;
    B = BofS( s, Info );

    s1 = s2 = s;
    B1 = B2 = B;

    /*
     * Scan up FL to bracket the northern mirror point
     */
    done = FALSE;
    while( !done && (s <= Info->s[Info->nPnts-1]) ) {
        s += .05;
        if (s > Info->s[Info->nPnts-1]) s = Info->s[Info->nPnts-1];

        B = BofS( s, Info );
        if ( B > Bm ) {
            /*
             * we passed the value of Bm (i.e. Bm-Btmp is negative) , use this
             * s as an upper bracket.
             */
            s2 = s;
            B2 = B;
            done = TRUE;
        } else {
            /*
             * havent passed the value of Bm yet, update s1 to be lower edge
             * of bracket
             */
            s1 = s;
            B1 = B;
        }
    }


    /*
     * Find zero using bi-section
     */
    done = FALSE;
    while ( !done ){
        s = 0.5*(s2+s1); // bisect
        B = BofS( s, Info );
        d = Bm-B;
        if ( fabs(d) < 1e-7)  {
            done = TRUE;        // converged
        } else if ( d > 0.0 ) {
            s1 = s;
            B1 = B;
        } else {
            s2 = s;
            B2 = B;
        }
    }

    *sn = s;





    s = Info->Smin;
    B = BofS( s, Info );

    s1 = s2 = s;
    /*
     * Scan down FL to bracket the southern mirror point
     */
    done = FALSE;
    while( !done && (s <= Info->s[Info->nPnts-1]) ) {
        s -= .05;
        if (s < 0.0) s = 0.0;
        if ( (B = BofS( s, Info )) > Bm ) {
            /*
             * we passed the value of Bm (i.e. Bm-Btmp is negative) , use this
             * s as an upper bracket.
             */
            s2 = s;
            B2 = B;
            done = TRUE;
        } else {
            /*
             * havent passed the value of Bm yet, update s1 to be lower edge
             * of bracket
             */
            s1 = s;
            B1 = B;
        }
    }

    /*
     * Find zero using bi-section
     */
    done = FALSE;
    while ( !done ){
        s = 0.5*(s2+s1); // bisect
        B = BofS( s, Info );
        d = Bm-B;
        if ( fabs(d) < 1e-7)  {
            done = TRUE;        // converged
        } else if ( d > 0.0 ) {
            s1 = s;
            B1 = B;
        } else {
            s2 = s;
            B2 = B;
        }
    }

    *ss = s;


    return( TRUE );

}

