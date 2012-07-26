#include <stdio.h>
#include <math.h>
#include "Lgm/Lgm_QuadPack.h"
#include "Lgm/Lgm_MagModelInfo.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <string.h>

#define TRACE_TOL   1e-7
#define GOLD        0.38197



double Lgm_AlphaOfK_Func( double Kt, double Alpha, Lgm_MagModelInfo *m );



/**
 *  Do initial setup for AlphaOfK(). This involves setting time and doing an
 *  initial field trace to get m->Pmin set up properly.
 *
 *      \param[in]      d   Date/Time to use.
 *      \param[in]      u   Position (in GSM) to use.
 *      \param[in,out]  m   A properly initialized and configured Lgm_MagModelInfo structure.
 *
 *      \returns        TraceFlag, the flag retuend by the Lgm_Trace() call.
 *
 *      \author         Mike Henderson
 *      \date           2010-2011
 *
 */
int  Lgm_Setup_AlphaOfK( Lgm_DateTime *d, Lgm_Vector *u, Lgm_MagModelInfo *m ) {

    double      s;
    int         TraceFlag;
    Lgm_Vector  v1, v2, v3, v4, Bvec;


    /*
     * Set the coordinate transformations for the gievn date/time
     */
    Lgm_Set_Coord_Transforms( d->Date, d->Time, m->c );

    /*
     * Set local B-field magnitude
     */
    m->Bfield( u, &Bvec, m );
    m->Blocal = Lgm_Magnitude( &Bvec );

    /*
     * Trace the field line for the given position. This should be a fast
     * adaptive trace. I.e. -- no points are saved.
     */
    TraceFlag = Lgm_Trace( u, &v1, &v2, &v3, m->Lgm_LossConeHeight, TRACE_TOL, TRACE_TOL, m );
    if ( TraceFlag != LGM_CLOSED ) {
        // problem tracing FL?
        return(-5);
    }
    s = m->Trace_s;


    /*
     * If we are going to use interped routines, we need to pre-trace the line.
     * And we need to init the interp stuff.
     */
    if (  m->UseInterpRoutines ) {

        /*
         * Start at Southern Footpoint and trace to Northern Footpoint.
         */
        m->Hmax = s/200.0;
        //Lgm_TraceLine2( &v1, &v4, m->Lgm_LossConeHeight, s/200.0, 1.0, TRACE_TOL, FALSE, m );

        Lgm_TraceLine2( &v1, &v4, 0.0, s/200.0, 1.0, TRACE_TOL, FALSE, m );


        if ( !InitSpline( m ) ) return(-5);

    } 
    
    return( TraceFlag );

}


/**
 *  Clean up the setup we did in Lgm_Setup_AlphaOfK().
 *
 *      \param[in,out]  m   A Lgm_MagModelInfo structure that has had an InitSpline() performed on it.
 *
 *
 *      \author         Mike Henderson
 *      \date           2010-2011
 *
 */
void  Lgm_TearDown_AlphaOfK( Lgm_MagModelInfo *m ) {

    if ( m->AllocedSplines ) FreeSpline( m );

}


/**
 *   This routine returns the pitch angle that corresponds to a given value of K
 *   \f$ K = I \sqrt{B_m}\f$.
 *
 *      \param[in]      K   The value of the second invariant, K        <b> ( Re G^(1/2) )</b>
 *      \param[in,out]  m   A properly initialized and configured Lgm_MagModelInfo structure.
 *
 *      \returns        Pitch angle, \f$\alpha\f$ implied by K          <b> ( Degrees )</b>
 *
 *      \author         Mike Henderson
 *      \date           2010-2011
 *
 *      \note           You must have setup the Lgm_MagModelInfo structure
 *                      properly first. Then you need to call Lgm_Setup_AlphaOfK() before you
 *                      call this routine.
 *
 */
double  Lgm_AlphaOfK( double K, Lgm_MagModelInfo *m ) {

    double  a0, a1, a, B;
    double  f0, f1, f;
    int     done;


    /*
     *  Set up low side of bracket. The footpoints are at 100-ish km (or
     *  whatever), but the integral invariant is from mirror point to mirror
     *  point which is given by a single Bm.  The B vals at the footpoints are
     *  likely not the same. To get a valid I, we need to use the smallest of
     *  the B's -- i.e. find B at 100km N and S. Then use smallest one as the
     *  mirror field value. From that we can set a0 straight off. May need to
     *  adjust the value a small amount to avoid getting too close to the LC
     *  height.
     */
    if ( m->Ellipsoid_Footprint_Bs < m->Ellipsoid_Footprint_Bn ) {
        B = m->Ellipsoid_Footprint_Bs;
    } else {
        B = m->Ellipsoid_Footprint_Bn;
    }
    a0 = DegPerRad*asin( sqrt( m->Bmin/B ) );
    f0 = Lgm_AlphaOfK_Func( K, a0, m );
    if ( fabs(f0) < 1e-4 ) return( a0 );



    /*
     *  Set up high side of bracket. A PA of 90Deg. is as high as you can get.
     *  And this should give I=0, so no need to compute I
     */
    a1 = 90.0;
    f1 = K - 0.0;
    f1 = Lgm_AlphaOfK_Func( K, a1, m );
    if ( fabs(f1) < 1e-4 ) return( a1 );


    /*
     *  If the vals are not opposite signs, we dont have a proper bracket.
     */
    if ( f0*f1 > 0.0 ) {
        if (m->VerbosityLevel >= 2) printf("Lgm_AlphaOfK(): [a0:a1] = [%g: %g] does not bracket root? Lgm_AlphaOfK_Func( %g, %g ) = %g, Lgm_AlphaOfK_Func( %g, %g ) = %g\n", a0, a1, K, a0, f0, K, a1, f1);
        return(-9e99);
    }


    done = FALSE;
    while( !done ) {

        /*
         * compute a new PA to test
         */
        a = (a1-a0)*GOLD + a0;
        //a = (a1-a0)*0.5 + a0;
        f = Lgm_AlphaOfK_Func( K, a, m );

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


/**
 *  This internal function returns K(Alpha).
 *
 *      \param[in]      Alpha           The pitch angle to compute Kt-K(Alpha) for         <b> ( Degrees )</b>
 *      \param[in,out]  m               A properly initialized and configured Lgm_MagModelInfo structure.
 *
 *      \returns        K(Alpha).                                 <b> ( Re G^(1/2) )</b>
 *
 *      \author         Mike Henderson
 *      \date           2010-2011
 *
 */
double Lgm_KofAlpha( double Alpha, Lgm_MagModelInfo *m ) {

    double  rat, sa, sa2, Sma, Smb, I, K;

    m->PitchAngle = Alpha;
    sa = sin( Alpha*RadPerDeg ); sa2 = sa*sa;
    m->Bm = m->Bmin/sa2; 



    /*
     * Trace from Bmin point up to northern mirror point and down to
     * southern mirror point. Sma and Smb are the distances along
     * (i.e. up and down) the FL from the Bmin point.
     */
     if ( Lgm_TraceToMirrorPoint( &(m->Pmin), &(m->Pm_South), &Sma, m->Bm, -1.0, m->Lgm_TraceToMirrorPoint_Tol, m ) >= 0 ) {
        if ( Lgm_TraceToMirrorPoint( &(m->Pm_South), &(m->Pm_North), &Smb, m->Bm,  1.0, m->Lgm_TraceToMirrorPoint_Tol, m ) >= 0 ) {

            
            /*
             *  Set the limits of integration. Also set tolerances for
             *  Quadpack routines. Note that we have pre-traced the FL already.
             *  In the pre-traced arrays, s=0 is the south foot, and smax is
             *  the north. foot. Also, m->smin is distance from southern foot
             *  to Pmin. Therefore, if we want to use these arrays in the
             *  Iinv_interped routines, we better figure out what [a,b] should
             *  be.
             */


            if ( Smb <= 1e-5 ) {
                // if FL length is small, use an approx expression for I
                rat = m->Bmin/m->Bm;
                if ((1.0-rat) < 0.0) {
                    I = 0.0;
                } else {
                    // Eqn 2.66b in Roederer
                    I = Smb*sqrt(1.0 - m->Bmin/m->Bm);
                }

            } else if (  m->UseInterpRoutines ) {

                /*
                 *  Do interped I integral.
                 */
                m->Sm_South = m->Smin - Sma;
                m->Sm_North = m->Sm_South + Smb;
//printf("m->Smin = %g    Sma, Smb = %g %g\n", m->Smin, Sma, Smb );
//printf("m->Sm_South = %g    m->Sm_North = %g\n", m->Sm_South, m->Sm_North );
//exit(0);
                I = Iinv_interped( m );
                // Compute K(Alpha) (units of G^1/2 Re)
                K = 3.16227766e-3*I*sqrt(m->Bm);
                if (m->VerbosityLevel >= 2) printf("Lgm_AlphaOfK(): Iinv (Interped Integral) = %g   K = %g\n",  I, K );

            } else {

                /*
                 *  Do full blown I integral. 
                 */
                m->Sm_South = 0.0;
                m->Sm_North = Smb;
                I = Iinv( m );
                // Compute K(Alpha) (units of G^1/2 Re)
                K = 3.16227766e-3*I*sqrt(m->Bm);
                if (m->VerbosityLevel >= 2) printf("Lgm_AlphaOfK(): Iinv (Full Integral) = %g   K = %g\n",  I, K );

            }


            /*
             *  return K
             */
            return( K );


        }
     }

     return( -9e99 );


}

/**
 *  This internal function returns the difference between the target K and the
 *  K(Alpha). Its the quanity whose zero we are trying to find. I.e., we are
 *  trying to find the Alpha that makes;
 *
 *      Func = K - K(Alpha) = 0
 *
 *  On a plot of f(Alpha) vs Alpha, the value starts out negative for low pitch
 *  angles and crosses zero at some PA, then goes positive. 
 *
 *      \param[in]      Kt              The target value of the second invariant, K        <b> ( Re G^(1/2) )</b>
 *      \param[in]      Alpha           The pitch angle to compute Kt-K(Alpha) for         <b> ( Degrees )</b>
 *      \param[in,out]  m               A properly initialized and configured Lgm_MagModelInfo structure.
 *
 *      \returns        difference between Kt and K(Alpha).                                 <b> ( Re G^(1/2) )</b>
 *
 *      \author         Mike Henderson
 *      \date           2010-2011
 *
 *      \note           This is an internal routine that the user probably should not call. 
 *                      (Perhaps an exception could be if you wanted to know the final difference 
 *                      follwoing a call to Lgm_AlphaOfK().)
 */
double Lgm_AlphaOfK_Func( double Kt, double Alpha, Lgm_MagModelInfo *m ) {

    double K;

    if ( fabs( Alpha - 90.0 ) < 1e-5 ) return( Kt - 0.0 );

    K = Lgm_KofAlpha( Alpha, m );

    if ( K < 0.0 ) {

        return( -9e99 );

    } else {

        /*
         *  return the diff
         */
        return( Kt - K );

    }

}
