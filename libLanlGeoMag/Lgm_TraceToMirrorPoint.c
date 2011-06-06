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

#define LC_TOL  0.99    // Allow height to be as low as .99*Lgm_LossConeHeight, before we call it "in the loss cone"


double mpFunc( Lgm_Vector *P, double Bm, Lgm_MagModelInfo *Info ){

    Lgm_Vector  Bvec;
    double      F;

    Info->Bfield( P, &Bvec, Info );
    F = Lgm_Magnitude( &Bvec ) - Bm;
    return( F );

}

int Lgm_TraceToMirrorPoint( Lgm_Vector *u, Lgm_Vector *v, double *Sm, double Bm, double sgn, double tol, Lgm_MagModelInfo *Info ) {

    Lgm_Vector	u_scale;
    double	    Htry, Hdid, Hnext, Hmin, Hmax, s;
    double	    Sa=0.0, Sb=0.0, Smin, d;
    double	    Rlc, R, Fa, Fb, F, Fmin, B, Fs, Fn;
    double	    Ra, Rb, Height;
    Lgm_Vector	w, Pa, Pb, P, Bvec, Pmin;
    int		    done, FoundBracket, reset, nIts;

    reset = TRUE;
    Fmin = 9e99;

    if ( Info->VerbosityLevel > 4 ) {
        printf( "\n\n**************** Start Detailed Output From TraceToMirrorPoint (VerbosityLevel = %d) ******************\n", Info->VerbosityLevel );
    }


    *Sm = 0.0;

    /*
     *  If particle mirrors below Info->Lgm_LossConeHeight, we are in the loss cone.
     */
    Lgm_Convert_Coords( u, &w, GSM_TO_WGS84, Info->c );
    Lgm_WGS84_to_GeodHeight( &w, &Height );
    if ( Height < LC_TOL*Info->Lgm_LossConeHeight ) {
        if ( Info->VerbosityLevel > 1 )
        printf("    Lgm_TraceToMirrorPoint: Current Height is below specified loss cone height of %g km. In Loss Cone. (Height = %g) \n", Info->Lgm_LossConeHeight, Height );
        return(-1); // below loss cone height -> particle is in loss cone!
    }


    /*
     * Put in some logic to bail if the initial point is below Bm (i.e. if B is
     * significantly higher than Bm). This ammounts to requiring that the input
     * point is at least on a part of the F.L. accessible to bouncing particles
     * (that have Bm as specified.)
     */





    Hmax = Info->Hmax;
    Hmin = 1e-10;
    u_scale.x =  100.0;  u_scale.y = 100.0; u_scale.z = 100.0;
    R = Ra = Rb = 0.0;
    F = Fa = Fb = 0.0;




    /*
     *  Bracket the root first.
     *  We are attempting to find a root in the direction given to us (via ithe
     *  sgn value).  We want to find where F=0 (i.e. where B = Bm).  Thus we
     *  need a bracket [Pa, Pb] such that Fa is negative and Pb is positive.
     *  Ideally, we would like to start at Bmin and trace up or down to get the
     *  mirror points. However, we could easily have a situation where the
     *  input point is already at (or very close) to one of the mirror points.
     *  For example, if we put in the S/C location, locally mirroring particles
     *  (PA=90) will have a mirror point at exactly that location and another
     *  one at some point beyond the Bmin point. To find the other root, we
     *  first need to get a point Pa between them (so that Fa<0). When the
     *  input point is numerically very close to a mirror point, F will be +/-
     *  a very small number, so it cant be used as the first point of the
     *  bracket. Instead, we need to step a tiny bit to get the first bracket.
     */
    if ( Info->VerbosityLevel > 4 )  {
        Info->Bfield( u, &Bvec, Info );
        R = Lgm_Magnitude( u );
        F = Lgm_Magnitude( &Bvec ) - Bm;
        printf("    TraceToMirrorPoint: Starting Point: %15g %15g %15g   R, F (=B-Bm) = %g %g\n\n", u->x, u->y, u->z, R, F );
    }





    /*
     *  To proceed, we need to find a point to start the bracket on (i.e. a
     *  point where B-Bm < 0).  Then we need to trace until the sign changes
     *  (or we bail due to hitting earth or because its open or whatever.)
     *
     *  First, try a small step in the user-supplied direction.
     */
    Htry = 1e-8; // we probably dont ever need to split the mirror points to any finer precision than this(?).
    P    = *u;
    if ( Lgm_MagStep( &P, &u_scale, Htry, &Hdid, &Hnext, 1.0e-8, sgn, &s, &reset, Info->Bfield, Info ) < 0 ) return( LGM_BAD_TRACE );
    Info->Bfield( &P, &Bvec, Info );
    B = Lgm_Magnitude( &Bvec );
    F = B-Bm;
//printf("MIKE: B, Bm = %.15g %.15g\n", B, Bm);

    if ( F < 0.0 ) {

        /*
         *  If F < 0.0, then it means |B| is still less than Bm in this
         *  direction.  Save this point as the start of a potential bracket.
         *  Then continue tracing in this direction.
         */
        Pa  = P;
        Ra  = Lgm_Magnitude( &P );
        Sa  = Hdid;
        Fa  = F;

    } else {

        /*
         * Even with a small step size in this direction, |B| was not found to
         * be less than Bm. This could mean that we are very close to the Pmin
         * point. Since we are trusting the user that this really is the
         * direction that the other root should have been found in, we must
         * conclude that the roots are very closely separated around Pmin. Just
         * return the input point back as the output point.
         */
//printf("HEREEEEEEEEEEEEEEEEEEEEEEEEEE\n");
         *v  = *u;
         *Sm = 0.0;
         return( 1 );

    }




    // Allow larger step sizes. Also, we can go beyond Bm to get the bracket.
    Hmax = 0.5;
    P    = Pa;

    /*
     *  To begin with, B - Bm will be negative (we wouldnt be here otherwise). So all we need to do is
     *  trace along the F.L. until B - Bm changes sign. That will give us the far side of the bracket.
     *  But dont go below surface of the Earth, etc.
     */
    done         = FALSE;
    FoundBracket = FALSE;
    while ( !done ) {

        if ( Info->VerbosityLevel > 4 ) {
            printf("    TraceToMirrorPoint, Finding Bracket: Starting P: %15g %15g %15g\n", P.x, P.y, P.z );
        }

        /*
         *  If user doesnt want large steps, limit Htry ...
         */
        if (Htry > Hmax) Htry = Hmax;

        if ( Lgm_MagStep( &P, &u_scale, Htry, &Hdid, &Hnext, 1.0e-7, sgn, &s, &reset, Info->Bfield, Info ) < 0 ) return(-1);

        /*
         *  Get value of quantity we want to minimize
         */
        Info->Bfield( &P, &Bvec, Info );
        R = Lgm_Magnitude( &P );
        F = Lgm_Magnitude( &Bvec ) - Bm;
        Lgm_Convert_Coords( &P, &w, GSM_TO_WGS84, Info->c );
        Lgm_WGS84_to_GeodHeight( &w, &Height );


        if (   (P.x > Info->OpenLimit_xMax) || (P.x < Info->OpenLimit_xMin) || (P.y > Info->OpenLimit_yMax) || (P.y < Info->OpenLimit_yMin)
                || (P.z > Info->OpenLimit_zMax) || (P.z < Info->OpenLimit_zMin) || ( Sa > 300.0 ) ) {
            /*
             *  Open FL!
             */
            v->x = v->y = v->z = 0.0;
            if ( Info->VerbosityLevel > 1 ) {
                printf("    Lgm_TraceToMirrorPoint(): FL OPEN\n");
            }
            if ( Info->VerbosityLevel > 4 ) {
                printf( "**************** End Detailed Output From TraceToMirrorPoint (VerbosityLevel = %d) ******************\n\n\n", Info->VerbosityLevel );
            }
            return(-2);

        } else if ( F > 0.0 ) { /* not >= because we want to explore at least a step beyond */

            done = TRUE;
            FoundBracket = TRUE;
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



        // Set Htry adaptively. But make sure we wont descend below the Earth's surface.
        Htry = Hnext;
        if (Htry > (R-1.0)) Htry = 0.95*(R-1.0);


        // If Htry gets to be too small, there is probably something wrong.
        if (Htry < 1e-12) done = TRUE;


        if ( Info->VerbosityLevel > 4 ) {
            printf("                                             Got To: %15g %15g %15g     with Htry of: %g\n", P.x, P.y, P.z, Htry );
            printf("        Pa, Ra, Fa, Sa = (%15g, %15g, %15g) %15g %15g %15g\n", Pa.x, Pa.y, Pa.z, Ra, Fa, Sa  );
            printf("        Pb, Rb, Fb, Sb = (%15g, %15g, %15g) %15g %15g %15g\n", Pb.x, Pb.y, Pb.z, Rb, Fb, Sb  );
            printf("        F = %g, |B| Bm = %g %g FoundBracket = %d  done = %d    Current Height = %g\n\n", F, Lgm_Magnitude( &Bvec ), Bm, FoundBracket, done, Height );
        }

        if ( Height < 0.0 ) {
            if ( Info->VerbosityLevel > 1 ) {
                printf("    Lgm_TraceToMirrorPoint: Current Height is below surface of the Earth. (Height = %g) \n", Height );
            }
            if ( Info->VerbosityLevel > 4 ) {
                printf( "**************** End Detailed Output From TraceToMirrorPoint (VerbosityLevel = %d) ******************\n\n\n", Info->VerbosityLevel );
            }
            return(-1); /* dropped below surface of the Earth */
        }

    }



    /*
     * We have found a potential bracket, but lets just be sure.
     */
    if ( (Fa>=0.0)&&(Fb<=0.0) || (Fa<=0.0)&&(Fb>=0.0) ) FoundBracket = TRUE;
    else FoundBracket = FALSE;

    if ( !FoundBracket && Fmin > 1e-4 ) {

        if ( Info->VerbosityLevel > 1 ) {
            printf("    Lgm_TraceToMirrorPoint: Bracket not found.\n");
            printf("Pmin = %.15lf, %.15lf, %.15lf,     Fmin = %.15lf\n", Pmin.x, Pmin.y, Pmin.z, Fmin);
            printf("Pa = %.15lf, %.15lf, %.15lf,     Fa = %.15lf\n", Pa.x, Pa.y, Pa.z, Fa);
            printf("Pb = %.15lf, %.15lf, %.15lf,     Fb = %.15lf\n", Pb.x, Pb.y, Pb.z, Fb);
        }
        if ( Info->VerbosityLevel > 4 ) {
            printf( "**************** End Detailed Output From TraceToMirrorPoint (VerbosityLevel = %d) ******************\n\n\n", Info->VerbosityLevel );
        }
        return( -3 ); /* We have gone as far as we could without finding a bracket. Bail out. */

    } else if ( !FoundBracket && Fmin <= 1e-4 ) {

        Fb = Fmin;
        Sb = Smin;
        Pb = Pmin;

    } else {



        /*
         *  We have a bracket. Now go in for the kill using bisection.
         *
         *  (Note: If all we wanted to know is whether or not the line hits the
         *  Earth, we could stop here: it must hit the Earth or we wouldnt
         *  have a minimum bracketed.)
         *
         */
    if (0==1){
        done  = FALSE;
        //reset = TRUE;
        if ( Info->VerbosityLevel > 4 ) nIts = 0;
        while (!done) {

            if ( Info->VerbosityLevel > 4 ) {
                printf("    TraceToMirrorPoint, Finding Root: Starting P: %15g %15g %15g  ...\n", P.x, P.y, P.z );
            }

            d = Sb - Sa;

            if ( (fabs(d) < tol)  ) {
                done = TRUE;
            } else {

                //P = Pa; Htry = 0.5*d;
                P = Pa; Htry = LGM_1_OVER_GOLD*d;
                if ( Lgm_MagStep( &P, &u_scale, Htry, &Hdid, &Hnext, 1.0e-8, sgn, &s, &reset, Info->Bfield, Info ) < 0 ) return(-1);
                Info->Bfield( &P, &Bvec, Info );
                F = Lgm_Magnitude( &Bvec ) - Bm;
                if ( F >= 0.0 ) {
                    Pb = P; Fb = F; Sb = Sa + Hdid;
                } else {
                    Pa = P; Fa = F; Sa += Hdid;
                }

            }

            if ( Info->VerbosityLevel > 4 ) {
                printf("                                             Got To: %15g %15g %15g     with Htry of: %g\n", P.x, P.y, P.z, Htry );
                printf("        Pa, Ra, Fa, Sa = (%15g, %15g, %15g) %15g %15g %15g\n", Pa.x, Pa.y, Pa.z, Ra, Fa, Sa  );
                printf("        Pb, Rb, Fb, Sb = (%15g, %15g, %15g) %15g %15g %15g\n", Pb.x, Pb.y, Pb.z, Rb, Fb, Sb  );
                printf("        F = %g, |B| Bm = %g %g FoundBracket = %d  done = %d    Current Height = %g\n\n", F, Lgm_Magnitude( &Bvec ), Bm, FoundBracket, done, Height );
                ++nIts;
            }

        }
    }



    //reset = TRUE;



    if (1==1){
        /*
         *  Use Brent's method
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
        f.func    = &mpFunc;
        f.Val     = Bm;
        Lgm_zBrentP( Sa, Sb, Fa, Fb, Pa, Pb, &f, tol, &Sz, &Fz, &Pz );
        Fb = Fz; Sb = Sz; Pb = Pz;
    //printf("Sa, Sb = %g %g  Fa, Fb = %g %g   tol = %g\n", Sa, Sb, Fa, Fb, tol);
    }

    }



    /*
     *  Take Pb as endpoint
     */
    *v = Pb; *Sm = Sb;


    /*
     * Do a final check on the Height to see if we are in loss cone or not.
     */
    Lgm_Convert_Coords( v, &w, GSM_TO_WGS84, Info->c );
    Lgm_WGS84_to_GeodHeight( &w, &Height );

	if ( Height < LC_TOL*Info->Lgm_LossConeHeight ) {
        // dropped below loss cone height -> particle is in loss cone!
        if ( Info->VerbosityLevel > 1 ) {
            printf("    Lgm_TraceToMirrorPoint: Current Height is below specified loss cone height of %g km. In Loss Cone. (Height = %g) \n", Info->Lgm_LossConeHeight, Height );
        }
        if ( Info->VerbosityLevel > 4 ) {
            printf( "**************** End Detailed Output From TraceToMirrorPoint (VerbosityLevel = %d) ******************\n\n\n", Info->VerbosityLevel );
        }
	    return(-1);
	}


    if ( Info->VerbosityLevel > 4 ) {
        printf("    =============================================================================================================\n" );
        printf("    |                                                                                                           |\n" );
        printf("    |    Total iterations to find root: %3d                                                                     |\n", nIts );
        printf("    |    Total Bfield evaluations find root: %4d                                                                 |\n", Info->Lgm_nMagEvals );
        printf("    |    TraceToMirrorPoint, Final Result: %15g %15g %15g     Sm = %g    |\n", v->x, v->y, v->z, *Sm );
        printf("    |                                                                                                           |\n" );
        printf("    =============================================================================================================\n", v->x, v->y, v->z, *Sm );
        printf( "**************** End Detailed Output From TraceToMirrorPoint (VerbosityLevel = %d) ******************\n\n\n", Info->VerbosityLevel );
    }

    if ( Info->VerbosityLevel > 2 ) printf("Lgm_TraceToMirrorPoint(): Number of Bfield evaluations = %d\n", Info->Lgm_nMagEvals );

    return( 1 );

}

