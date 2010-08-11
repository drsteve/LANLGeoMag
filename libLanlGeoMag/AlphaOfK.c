#include <stdio.h>
#include <math.h>
#include "Lgm/Lgm_QuadPack.h"
#include "Lgm/Lgm_MagModelInfo.h"
#include "Lgm/Lgm_LstarInfo.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <string.h>

#define TRACE_TOL   1e-7
#define GOLD        0.38197



double Func( double Kt, double Alpha,Lgm_LstarInfo *LstarInfo );


/*
 *   This routine returns the pitch angle that corresponds to a given value of K
 *   K = sqrt(Bm)I
 *
 */

double  AlphaOfK( double K,Lgm_LstarInfo *LstarInfo ) {

    double  a0, a1, a;
    double  f0, f1, f;
    int     done;


    /*
     *  Set up low side of bracket.
     */
    a0 = 10.0;
    f0 = Func( K, a0, LstarInfo );
    if ( fabs(f0) < 1e-4 ) return( a0 );


    /*
     *  Set up high side of bracket.
     */
    a1 = 90.0;
    f1 = Func( K, a1, LstarInfo );
    if ( fabs(f1) < 1e-4 ) return( a1 );



    /*
     *  If the vals are not opposite signs, we dont have a proper bracket.
     */
    if ( f0*f1 > 0.0 ) {
        printf("AlphaOfK: [a0:a1] = [%g: %g] does not bracket root?\n", a0, a1);
        return(-9e99);
    }


    done = FALSE;
    while( !done ) {

        /*
         * compute a new PA to test
         */
        //a = (a1-a0)*GOLD + a0;
        a = (a1-a0)*0.5 + a0;
        f = Func( K, a, LstarInfo );

        if ( fabs(a1-a0) < 1e-2 ) {
            done = TRUE;
        } else if ( f1*f < 0.0 ) {
            /*
             *  Root must be in [a:a1]
             *  Reset brackets.
             */
            a0 = a;
            f0 = f;
        } else {
            /*
             *  Root must be in [a0:a]
             *  Reset brackets.
             */
            a1 = a;
            f1 = f;
        }

    }

    /*
     * Take the midpoint of the remaining bracket range as the answer.
     */
    a = 0.5*(a0+a1);

    return( a );



}


/*
 *  This function returns the difference between
 *  the target K and the K(Alpha). I.e.;
 *
 *      Func = K - K(Alpha)
 *
 *  On a plot of f(Alpha) vs Alpha, the value starts out
 *  negative for low pitch angles and crosses zero at some
 *  PA, then goes positive. 
 */
double Func( double Kt, double Alpha,Lgm_LstarInfo *LstarInfo ) {

    double  sa, sa2, Sma, Smb, I, K;


    LstarInfo->PitchAngle = Alpha;
    sa = sin( Alpha*RadPerDeg ); sa2 = sa*sa;
    LstarInfo->mInfo->Bm =LstarInfo->mInfo->Blocal/sa2; 

    /*
     * Trace from Bmin point up to northern mirror point and down to
     * southern mirror point. Sma and Smb are the distances along
     * (i.e. up and down) the FL from the Bmin point.
     */
     if ( Lgm_TraceToMirrorPoint( &(LstarInfo->mInfo->Pmin), &(LstarInfo->mInfo->Pm_South), &Sma, 120.0, LstarInfo->mInfo->Bm, -1.0, LstarInfo->mInfo->Lgm_TraceToMirrorPoint_Tol, LstarInfo->mInfo ) ) {
        if ( Lgm_TraceToMirrorPoint( &(LstarInfo->mInfo->Pm_South), &(LstarInfo->mInfo->Pm_North), &Smb, 120.0, LstarInfo->mInfo->Bm,  1.0, LstarInfo->mInfo->Lgm_TraceToMirrorPoint_Tol, LstarInfo->mInfo ) ) {

            
            /*
             *  Set the limits of integration. Also set tolerances for
             *  Quadpack routines.
             */
            //LstarInfo->mInfo->Sm_South = LstarInfo->mInfo->smin - Sma;
            //LstarInfo->mInfo->Sm_North = LstarInfo->mInfo->smin + Smb;
            LstarInfo->mInfo->Sm_South = 0.0;
            LstarInfo->mInfo->Sm_North = Smb;
//            LstarInfo->mInfo->epsrel = 0.0;
//            LstarInfo->mInfo->epsabs = 1e-1;
//            LstarInfo->mInfo->Lgm_I_Integrator  = DQK21;
//            LstarInfo->mInfo->UseInterpRoutines = TRUE;

            /*
             *  Since the Pm_South and Pm_north points are not likely
             *  to be in the interped array, add them (since we have
             *  them now anyway!). This is particularly important to do
             *  for the Sb integrals because the integrand is high near
             *  the mirror points.  If we dont do it, we could get
             *  slight differences due to it.
             */
            FreeSpline( LstarInfo->mInfo );
            //AddNewPoint( LstarInfo->mInfo->smin-Sma, LstarInfo->mInfo->Bm, &LstarInfo->mInfo->Pm_South, LstarInfo->mInfo );
            //AddNewPoint( LstarInfo->mInfo->smin+Smb, LstarInfo->mInfo->Bm, &LstarInfo->mInfo->Pm_North, LstarInfo->mInfo );
            AddNewPoint( 0.0, LstarInfo->mInfo->Bm, &LstarInfo->mInfo->Pm_South, LstarInfo->mInfo );
            AddNewPoint( Smb, LstarInfo->mInfo->Bm, &LstarInfo->mInfo->Pm_North, LstarInfo->mInfo );
            InitSpline( LstarInfo->mInfo );
//            LstarInfo->mInfo->UseInterpRoutines = TRUE;

            if (  LstarInfo->mInfo->UseInterpRoutines ) {

                /*
                 *  Do interped I integral.
                 */
                I = Iinv_interped( LstarInfo->mInfo  );
                if (LstarInfo->VerbosityLevel >= 2) printf("Iinv (Interped Integral) = %g\n",  I );

            } else {

                /*
                 *  Do full blown I integral. 
                 */
                I = Iinv( LstarInfo->mInfo  );
                if (LstarInfo->VerbosityLevel >= 2) printf("Iinv (Full Integral) = %g\n",  I );

            }


            /*
             *  Compute K(Alpha) (units of G^1/2 Re)
             */
            K = 3.16227766e-3*I*sqrt(LstarInfo->mInfo->Bm);


            /*
             *  return the diff
             */
            return( Kt - K );


        }
     }

     return( -9e99 );


}
