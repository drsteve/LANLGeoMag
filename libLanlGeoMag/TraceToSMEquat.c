/* TraceToSMEquat, Copyright (c) 1999 Michael G. Henderson <mghenderson@lanl.gov>
 *
 *    - Traces down to the SM (solar magnetic) equatorial plane.
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



int Lgm_TraceToSMEquat( Lgm_Vector *u, Lgm_Vector *v, double tol, Lgm_MagModelInfo *Info ) {

    Lgm_Vector	u_scale;
    double	    Htry, Hdid, Hnext, Hmin, Hmax, s, sgn, r2;
    double	    Sa=0.0, Sb=0.0, Sc=0.0, d1, d2;
    double	    Za, Zb, Zc, Z;
    Lgm_Vector	Ztmp;
    Lgm_Vector	Pa, Pb, Pc, P;
    int		    done, reset;



    Hmax = 20.0;
    Hmin = 0.001;
    u_scale.x =  100.0;  u_scale.y = 100.0; u_scale.z = 100.0;
    Za = Zb = Zc = 0.0;



    /*
     *  Bracket the minimum. We want to find two points along
     *  the field line such that the location of equatorial plane is
     *  gauranteed to lie between them. To do this, we need to find
     *  three points; Pa, Pb, and Pc such that;
     *
     *		     Pc > Pb > Pa
     *
     *		and  |z_sm( Pa )|   >   |z_sm( Pb )|
     *		and  |z_sm( Pc )|   >   |z_sm( Pb )|
     *
     *
     */


    /* 
     *  Set the start point, Pa and Za (the absolute value of the SM z coord).
     */
    Pa   = *u;
    Lgm_Convert_Coords( &Pa, &Ztmp, GSM_TO_SM, Info->c );
    Za   = fabs( Ztmp.z );
    Sa   = 0.0;




    /*
     *  If we are above the SM plane (i.e. in the northern hemisphere), lets 
     *  try tracing against the field direction.
     */
    sgn = ( Ztmp.z > 0.0 ) ? -1.0 : 1.0;



    /*
     *  Keep stepping along the field line until we have a |Zsm| less than the start point.
     */
    P    = Pa;
    done = FALSE;
    Htry = Za;
    while ( !done ) {

        if ( Lgm_MagStep( &P, &u_scale, Htry, &Hdid, &Hnext, 1.0e-7, sgn, &s, &reset, Info->Bfield, Info ) < 0 ) return(-1);
        Lgm_Convert_Coords( &P, &Ztmp, GSM_TO_SM, Info->c );
        Z = fabs( Ztmp.z );
            if (   (P.x > Info->OpenLimit_xMax) || (P.x < Info->OpenLimit_xMin) || (P.y > Info->OpenLimit_yMax) || (P.y < Info->OpenLimit_yMin)
                || (P.z > Info->OpenLimit_zMax) || (P.z < Info->OpenLimit_zMin) ) {
	        /*
	         *  Open FL!
	         */
	        v->x = v->y = v->z = 0.0;
	        return(0);
	    } else if ( Z < Za ) {
	        done = TRUE;
	        Pb = P;
	        Zb = Z;
	        Sb += Hdid;
        } else {
            Pa = P;
            Za = Z;
            Sa = 0.0;
        }

	    Htry = Z;

    }



    /*
     *  Keep stepping along the field line until we complete the bracket triple.
     */
    done = FALSE;
    while (!done) {

	    P = Pb;
        if ( Lgm_MagStep( &P, &u_scale, Htry, &Hdid, &Hnext, 1.0e-7, sgn, &s, &reset, Info->Bfield, Info ) < 0 ) return(-1);
        Lgm_Convert_Coords( &P, &Ztmp, GSM_TO_SM, Info->c );
        Z = fabs( Ztmp.z );
            if (   (P.x > Info->OpenLimit_xMax) || (P.x < Info->OpenLimit_xMin) || (P.y > Info->OpenLimit_yMax) || (P.y < Info->OpenLimit_yMin)
                || (P.z > Info->OpenLimit_zMax) || (P.z < Info->OpenLimit_zMin) ) {
	        /*
	         *  Open FL!
	         */
	        v->x = v->y = v->z = 0.0;
	        return(0);
	    } else if ( Z < Zb ) {
	        Pa = Pb; Za = Zb; Sa = 0.0;
	        Pb = P;  Zb = Z;  Sb = Hdid; 
            r2  = P.x*P.x + P.y*P.y + P.z*P.z;
            Htry = (Hnext < Hmax) ? ((Hnext > Hmin) ? Hnext : Hmin) : Hmax;
            if (Htry*Htry > r2) Htry = 0.25*sqrt(r2);
	    } else {
	        Pc = P;  Zc = Z; Sc = Sb + Hdid;
	        done = TRUE;
	    }


    }
	






    /*
     *  We have a bracket. Now go in for the kill.
     *  Use golden-section search to converge toward minimum.
     *  (Sa, Sb, Sc) are the distances of the triple points along
     *  the FL. For convenience, we maintain Sa = 0.0, so that
     *  Sb and Sc are always the distances relative to Pa. (And
     *  so Sc is always the width of the bracketed interval. So when
     *  Sc gets small enough we bail out and take Pb as the min).
     */
    done = FALSE;
    while (!done) {

	    d1 = Sb;
	    d2 = Sc - Sb;
	    if ( Sc < tol ) {
	        done = TRUE;
	    } else if ( d1 > d2 ) {

	        P = Pa; Htry = 0.381966011*d1;
            if ( Lgm_MagStep( &P, &u_scale, Htry, &Hdid, &Hnext, 1.0e-7, sgn, &s, &reset, Info->Bfield, Info ) < 0 ) return(-1);
            Lgm_Convert_Coords( &P, &Ztmp, GSM_TO_SM, Info->c );
            Z = fabs( Ztmp.z );
	        if ( Z < Zb ) {
		        Pb = P; Zb = Z; Sb = Hdid;
	        } else {
		        Pa = P; Za = Z; Sb -= Hdid; Sc -= Hdid;
	        }
	    
	    } else {

	        P = Pb; Htry = 0.381966011*d2;
            if ( Lgm_MagStep( &P, &u_scale, Htry, &Hdid, &Hnext, 1.0e-7, sgn, &s, &reset, Info->Bfield, Info ) < 0 ) return(-1);
            Lgm_Convert_Coords( &P, &Ztmp, GSM_TO_SM, Info->c );
            Z = fabs( Ztmp.z );
	        if ( Z < Zb ) {
		        Pb = P; Zb = Z; Sb += Hdid;
	        } else {
		        Pc = P; Zc = Z; Sc = Sb + Hdid;
	        }

	    }

    }


    /*
     *  Take location of the Min-B surface to be Pb.
     */
    *v = Pb;


    return( 1 );

}

