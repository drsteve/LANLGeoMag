#include <Lgm_MagModelInfo.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#define TRACE_TOL   1e-7
#define KP_DEFAULT  1




/*
 *      Input Variables:
 *
 *                      Date:
 *                       UTC:
 *                         u:  Input position vector in GSM
 *                    nAlpha:  Number of Pitch Angles to compute
 *                     Alpha:  Pitch Angles to compute
 *                   Quality:  Quality factor (0-8)
 *
 *      Input/OutPut Variables:
 *
 *              MagEphemInfo:  Structure used to input and output parameters/settings/results to/from routine.
 *
 */
void Lgm_McIlwain_L( long int Date, double UTC, Lgm_Vector *u, double Alpha, Lgm_MagModelInfo *mInfo ) {

    Lgm_Vector      v1, v2, v3, vv1, Bvec;
    double          sa, sa2, Blocal;
    int             i, LS_Flag, nn, tk, TraceFlag;



    /*
     * set coord transformations
     */
    Lgm_Set_Coord_Transforms( Date, UTC, mInfo->c );


    /*
     * Save S/C position to Lgm_MagModelInfo structure and compute Blocal.
     */
    mInfo->P_gsm = *u;
    mInfo->Bfield( u, &Bvec, mInfo );
    Blocal = Lgm_Magnitude( &Bvec );


    /*
     * Set Pitch Angle, sin(Alpha), sin^2(Alpha), and Bmirror
     */
    sa = sin( Alpha*RadPerDeg ); sa2 = sa*sa;




    /*
     *  First do a trace to identify the FL type and some of its critical points.
     */
    if ( Lgm_Trace( u, &v1, &v2, &v3, mInfo->Lgm_LossConeHeight, TRACE_TOL, TRACE_TOL, mInfo ) == LGM_CLOSED ) {

        /*
         * Trace from Bmin point up to northern mirror point and down to
         * southern mirror point. dSa and dSb are the distances along the
         * FL from the starting points point. So dSa is from the Bmin value.
         * And dSb is from the Psouth value.
         */
        if ( Lgm_TraceToMirrorPoint( &(mInfo->Pmin), &(mInfo->Pm_South), &dSa, mInfo->Bm, -1.0, mInfo->Lgm_TraceToMirrorPoint_Tol, mInfo ) >= 0 ) {

            if (VerbosityLevel > 0) {
                printf("\n\t\t%sMin-B  Point Location, Pmin (Re):      < %g, %g, %g >%s\n", PreStr, mInfo->Pmin.x, mInfo->Pmin.y, mInfo->Pmin.z, PostStr);
                printf("\n\t\t%sMirror Point Location, Pm_South (Re):      < %g, %g, %g >%s\n", PreStr, mInfo->Pm_South.x, mInfo->Pm_South.y, mInfo->Pm_South.z, PostStr);
                mInfo->Bfield( &mInfo->Pm_South, &Bvec, mInfo );
                B = Lgm_Magnitude( &Bvec );
                printf("\t\t%sMag. Field Strength, Bm at Pm_South (nT):  %g     (mInfo->Bm = %g)%s\n", PreStr, B, mInfo->Bm, PostStr);
            }

            if ( Lgm_TraceToMirrorPoint( &(mInfo->Pm_South), &(mInfo->Pm_North), &dSb, mInfo->Bm,  1.0, mInfo->Lgm_TraceToMirrorPoint_Tol, mInfo ) >= 0 ) {

                if ( mInfo->UseInterpRoutines ) {
                    Lgm_TraceLine2( &(mInfo->Pm_South), &mInfo->Pm_North, (r-1.0)*Re, 0.5*SS-mInfo->Hmax, 1.0, 1e-7, FALSE, mInfo );
                    ReplaceFirstPoint( 0.0, mInfo->Bm, &mInfo->Pm_South, mInfo );
                    AddNewPoint( SS,  mInfo->Bm, &mInfo->Pm_North, mInfo );
                    InitSpline( mInfo );

                    /*
                     *  Do interped I integral.
                     */
                    I = Iinv_interped( mInfo  );
                    if (mInfo->VerbosityLevel > 0) printf("Lgm_McIlwain_L: Integral Invariant, I (interped):      %g\n",  I );
                    FreeSpline( mInfo );

                } else {

                    /*
                     *  Do full blown I integral. (Integrand is evaluated by tracing to required s-values.)
                     */
                    I = Iinv( mInfo  );
                    if (VerbosityLevel > 0) printf("Lgm_McIlwain_L: Integral Invariant, I (full integral): %g\n",  I );

                }

            }

        }

    }


    return;

}
