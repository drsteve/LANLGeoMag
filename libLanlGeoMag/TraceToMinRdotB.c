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



int Lgm_TraceToMinRdotB( Lgm_Vector *u, Lgm_Vector *v, double tol, Lgm_MagModelInfo *Info ) {

    Lgm_Vector	u_scale;
    double	    Htry, Hdid, Hnext, Hmin, Hmax, s, sgn, sgn0;
    double	    Sa, Sb, B, f, r2, z2, r3, L, R, hhh;
    Lgm_Vector	Btmp;
    Lgm_Vector	Pa, Pb, P;
    int		    done, reset=TRUE, bracketed=FALSE;

    Pb.x = Pb.y = Pb.z = 0.0;


    Hmax = 20.0;
    Hmin = 0.01;
    u_scale.x = u_scale.y = u_scale.z = 1.0;


    /* 
     *  Set the start point, Pa and |B(Pa)|;
     */
    Pa   = *u;
    Info->Bfield( &Pa, &Btmp, Info );
    B    = Lgm_Magnitude( &Btmp );
    f    = Lgm_DotProduct( &Pa, &Btmp ) / B;
    Sa   = 0.0;
    sgn0 = ( f  < 0.0 ) ? -1.0 : 1.0;




    /*
     *  Choose a good step size to try.
     */
    Lgm_Convert_Coords( &Pa, &P, GSM_TO_SM, Info->c );
    R  = Lgm_Magnitude( &Pa );
    r2 = R*R;
    z2 = P.z*P.z;
    L  = 9e99;
    if ( fabs( r2 - z2 ) < 1e-2 ) {

        /*
         *  High  L
         */
        Htry = 10.0;

    } else {

        r3   = r2*R;
        L    = r3 / (r2 - z2);
        Htry = ( L < 4.5 ) ? 1.165*L : 4.5;

    }

    if ( Htry < 1e-4 ) Htry = 0.001;



    /* 
     *   Now, bracket minimum. And do a bisection search for the minimum.
     */
    done = FALSE;
    sgn  = sgn0;
    P    = Pa;
    Sb   = -9e99;
    while ( !done ) {

        //if ( Lgm_MagStep( &P, &u_scale, Htry, &Hdid, &Hnext, 1.0e-7, sgn, &s, &reset, Info->Bfield, Info ) < 0 ) return(-1);
        if ( Lgm_MagStep( &P, &u_scale, Htry, &Hdid, &Hnext, sgn, &s, &reset, Info->Bfield, Info ) < 0 ) return(-1);
        Info->Bfield( &P, &Btmp, Info );
        B = Lgm_Magnitude( &Btmp );
	    R = Lgm_Magnitude( &P );

        f = Lgm_DotProduct( &P, &Btmp ) / B;

	if ( fabs( Sb - Sa ) < tol ){

	    done = TRUE;

	} else if ( -f*sgn0 < 0 ){

	    Pa    = P;
	    Sa   += Hdid;
	    sgn   = sgn0;
	    if ( bracketed ) {
		Htry = fabs( Sb - Sa )/2.0;
	    } else if ( L < 4.5 ) {
		hhh = 1.165*L+1.0 - R + 1.0;
		if ( hhh > 1e-4 ) Htry = hhh;
Htry = (R-1.0)/2.0;

	    }
	    /*
printf("A: Sa, Sb, P = %g %g    (%g, %g, %g)   Htry = %f\n", Sa, Sb, P.x, P.y, P.z, Htry);
*/

	} else {

	    Pb  = P;
	    if ( bracketed ) Sb += sgn*sgn0*Hdid;
	    else	     Sb  = Sa + Hdid;
	    
	    bracketed = TRUE;
	    sgn  = -sgn0;
	    Htry = fabs( Sb - Sa )/2.0;
	    /*
printf("B: Sa, Sb, P = %g %g    (%g, %g, %g)   Htry = %f\n", Sa, Sb, P.x, P.y, P.z, Htry);
*/

	}

    }




    /*
     *  Take location of the Min-B surface to be Pb.
     */
    *v = Pb;



    return( 1 );

}

