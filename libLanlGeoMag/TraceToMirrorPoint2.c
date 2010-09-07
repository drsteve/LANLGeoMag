/* TraceToMirrorPoint, Copyright (c) 2004 Michael G. Henderson <mghenderson@lanl.gov>
 *
 *    - Starting at Minimum-B point, and given a value of mu_eq (mu = cos(pitch angle)),
 *	trace up field line to the particle's mirror point. The particle will mirror
 *     	when B = Beq/(1-mu_eq^2)
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
#include <omp.h>



int Lgm_TraceToMirrorPoint( Lgm_Vector *u, Lgm_Vector *v, double *Sm, double Bm, double sgn, double tol, Lgm_MagModelInfo *Info ) {

    Lgm_Vector	u_scale;
    double	    Htry, Hdid, Hnext, Hmin, Hmax, s;
    double	    Sa=0.0, Sb=0.0, d;
    double	    Rlc, R, Fa, Fb, F;
    double	    Ra, Rb, Height;
    Lgm_Vector	w, Pa, Pb, P, Bvec;
    int		    done, reset;


    /*
     *  If particle mirrors below Info->Lgm_LossConeHeight, we are in the loss cone.
     */
    Lgm_Convert_Coords( u, &w, GSM_TO_WGS84, Info->c );
    Lgm_WGS84_to_GeodHeight( &w, &Height );
    if ( Height < Info->Lgm_LossConeHeight ) {
        if ( Info->VerbosityLevel > 1 ) printf("Lgm_TraceToMirrorPoint: Initial Height is below %g km -- LOSS CONE \n", Info->Lgm_LossConeHeight );
        return(-1); // below loss cone height -> particle is in loss cone!
    }




    Hmax = Info->Hmax;
    Hmin = 1e-7;
    u_scale.x =  100.0;  u_scale.y = 100.0; u_scale.z = 100.0;
    R = Ra = Rb = 0.0;
    F = Fa = Fb = 0.0;




    /*
     *  Bracket the minimum first. 
     *  We want to stop when F = 0.0
     *
     *  Set the first point of the bracket to be the start point.
     *  A point of caution here...
     *  One might think that if Fa is identically zero (or very small) already 
     *  then we are done. I.e., if F=0, then the caller was asking for the
     *  mirror points of a locally mirroring particle -- so we are done right?
     *  Not quite.  The problem with this is that unless we are exactly
     *  in the min-B surface already, there are two mirror points: one where we   
     *  are already and another some distance away. So dont bail out yet.
     */
    Pa   = Pb = *u;
    Sa   = 0.0;
    Info->Bfield( &Pa, &Bvec, Info );
    Ra   = Lgm_Magnitude( &Pa );
    Fa   = Lgm_Magnitude( &Bvec ) - Bm;


    /*
     *  Set the step size to be XX% of the (linear) distance to surface of earth
     */
    Htry = 0.2*(Ra-1.0);


    /*
     *  If user doesnt want large steps, limit Htry ...
     */
    if (Htry > Hmax) Htry = Hmax;



    /*
     *  To begin with, B - Bm will be negative (or zero already). So all we need to do is
     *  trace along the F.L. until B - Bm changes sign. That will give us the far side of the bracket.
     *  But dont go below 120km altitude.
     */
    P     = Pa;

    done  = FALSE;
    reset = TRUE;
    while ( !done ) {

        if ( Htry > 0.1 ) Htry = 0.1;

        Lgm_MagStep( &P, &u_scale, Htry, &Hdid, &Hnext, 1.0e-7, sgn, &s, &reset, Info->Bfield, Info );

	    /*
 	     *  Get value of quantity we want to minimize
	     */
	    Info->Bfield( &P, &Bvec, Info );
        R = Lgm_Magnitude( &P );
	    F = Lgm_Magnitude( &Bvec ) - Bm;
        Lgm_Convert_Coords( &P, &w, GSM_TO_WGS84, Info->c );
        Lgm_WGS84_to_GeodHeight( &w, &Height );


        if (   (P.x > Info->OpenLimit_xMax) || (P.x < Info->OpenLimit_xMin) || (P.y > Info->OpenLimit_yMax) || (P.y < Info->OpenLimit_yMin)
                || (P.z > Info->OpenLimit_zMax) || (P.z < Info->OpenLimit_zMin) ) {
	        /*
	         *  Open FL!
	         */
	        v->x = v->y = v->z = 0.0;
            if ( Info->VerbosityLevel > 1 ) printf("TraceToMirrorPoint: FL OPEN\n");
	        return(-2);
	    } else if ( F > 0.0 ) { /* not >= because we want to explore at least a step beyond */
	        done = TRUE;
	        Pb = P;
	        Rb = R;
  	        Fb = F;
	        Sb = Sa + Hdid;
	    } else {
	        Pa = P;
	        Ra = R;
	        Fa = F;
	        Sa += Hdid;
	    }
	    /*
	     *  Go no farther than some small distance above 
	     *  the surface of the Earth.
	     */
	    Htry = fabs(0.2*(R-0.999999));
        if (Htry < 1e-12) done = TRUE;

	    if ( Height < Info->Lgm_LossConeHeight ) {
            if ( Info->VerbosityLevel > 1 ) printf("Lgm_TraceToMirrorPoint: Mirror Height is below %g km -- LOSS CONE \n", Info->Lgm_LossConeHeight );
	        return(-1); /* dropped below loss cone height -> particle is in loss cone! */
	    }

    }



    /*
     *  We have a bracket. Now go in for the kill using bisection.
     *
     *  (Note: If all we wanted to know is whether or not the line hits the
     *  Earth, we could stop here: it must hit the Earth or we wouldnt
     *  have a minimum bracketed.)
     *
     */
    done  = FALSE;
    reset = TRUE;
    while (!done) {

	    d = Sb - Sa;

	    if ( (fabs(d) < tol)  ) {
	        done = TRUE;
	    } else {

	        P = Pa; Htry = 0.5*d;
            Lgm_MagStep( &P, &u_scale, Htry, &Hdid, &Hnext, 1.0e-8, sgn, &s, &reset, Info->Bfield, Info );
	        Info->Bfield( &P, &Bvec, Info );
	        F = Lgm_Magnitude( &Bvec ) - Bm;
	        if ( F >= 0.0 ) {
		        Pb = P; Fb = F; Sb = Sa + Hdid;
	        } else {
		        Pa = P; Fa = F; Sa += Hdid; 
	        }
	    
	    }

    }



    /*
     *  Take average of endpoints as final answer
     */
//    v->x = 0.5*(Pa.x + Pb.x); v->y = 0.5*(Pa.y + Pb.y); v->z = 0.5*(Pa.z + Pb.z);
//    *Sm = 0.5*(Sa+Sb);
    *v = Pb;
    *Sm = Sb;


    return( 1 );

}

