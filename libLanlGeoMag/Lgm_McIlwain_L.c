#include "Lgm/Lgm_MagModelInfo.h"                                                                                                                                                                                
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#define TRACE_TOL   1e-7

/*
 *   Lgm_McIlwain_L
 *   -------------
 */


//! Compute McIlwain L-shell parameter for a given date, time, location and pitch angle.
/**
 *            \param[in]        Date        Date in format (e.g. 20101231). 
 *            \param[in]        UTC         Universal Time (Coordinated) in decimal hours (e.g. 23.5).
 *            \param[in]        u           Position to compute L-shell.
 *            \param[in]        Alpha       Pitch angle to compute L for. In degrees.
 *            \param[in]        Type        Flag to indicate which alogorithm to use (0=original McIlwain; else use Hilton's formula).
 *            \param[out]       I           The integral invariant, I that was computed along the way.
 *            \param[out]       Bm          The mirror magnetic field value, Bm that was computed along the way.
 *            \param[out]       M           The dipole magnetic moment used to compute L = f(I, Bm, M)
 *            \param[in,out]    mInfo       Properly initialized Lgm_MagModelInfo structure. (A number of otherm usefull things will have been set in mInfo).
 *
 *            \return           L           McIlwain L-shell parameter (a dimensioless number).
 *
 */
double Lgm_McIlwain_L( long int Date, double UTC, Lgm_Vector *u, double Alpha, int Type, double *I, double *Bm, double *M, Lgm_MagModelInfo *mInfo ) {

    Lgm_Vector      v1, v2, v3, Bvec;
    double          B, sa, sa2, Blocal, dSa, dSb, r, SS, L;

    *I  = -9e99;
    *Bm = -9e99;
    *M  = -9e99;
    L   = -9e99;


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
    mInfo->Bm = Blocal/sa2;



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

            if (mInfo->VerbosityLevel > 0) {
                printf("\n\tMin-B  Point Location, Pmin (Re):      < %g, %g, %g >\n", mInfo->Pmin.x, mInfo->Pmin.y, mInfo->Pmin.z );
                printf("\tMirror Point Location, Pm_South (Re):      < %g, %g, %g >\n", mInfo->Pm_South.x, mInfo->Pm_South.y, mInfo->Pm_South.z );
                mInfo->Bfield( &mInfo->Pm_South, &Bvec, mInfo );
                B = Lgm_Magnitude( &Bvec );
                printf("\tMag. Field Strength, Bm at Pm_South (nT):  %g     (mInfo->Bm = %g)\n", B, mInfo->Bm );
            }



            if ( Lgm_TraceToMirrorPoint( &(mInfo->Pm_South), &(mInfo->Pm_North), &dSb, mInfo->Bm,  1.0, mInfo->Lgm_TraceToMirrorPoint_Tol, mInfo ) >= 0 ) {

                if (mInfo->VerbosityLevel > 0) {
                    printf("\n\tMin-B  Point Location, Pmin (Re):      < %g, %g, %g >\n", mInfo->Pmin.x, mInfo->Pmin.y, mInfo->Pmin.z );
                    printf("\tMirror Point Location, Pm_North (Re):      < %g, %g, %g >\n", mInfo->Pm_North.x, mInfo->Pm_North.y, mInfo->Pm_North.z );
                    mInfo->Bfield( &mInfo->Pm_North, &Bvec, mInfo );
                    B = Lgm_Magnitude( &Bvec );
                    printf("\tMag. Field Strength, Bm at Pm_North (nT):  %g     (mInfo->Bm = %g)\n", B, mInfo->Bm );
                }

                /*
                 *  Set the limits of integration. Define s=0 at the sourthern mirror point. Then, sm_North will just be dSb
                 */
                SS = dSb;
                mInfo->Hmax = SS/200.0;
                r  = Lgm_Magnitude( &mInfo->Pm_North );
                mInfo->Sm_South = 0.0;
                mInfo->Sm_North = SS;

                if ( mInfo->UseInterpRoutines ) {
                    Lgm_TraceLine2( &(mInfo->Pm_South), &mInfo->Pm_North, (r-1.0)*Re, 0.5*SS-mInfo->Hmax, 1.0, 1e-7, FALSE, mInfo );
                    ReplaceFirstPoint( 0.0, mInfo->Bm, &mInfo->Pm_South, mInfo );
                    AddNewPoint( SS,  mInfo->Bm, &mInfo->Pm_North, mInfo );
                    InitSpline( mInfo );

                    /*
                     *  Do interped I integral.
                     */
                    *I = Iinv_interped( mInfo  );
                    if (mInfo->VerbosityLevel > 0) printf("Lgm_McIlwain_L: Integral Invariant, I (interped):      %g\n",  *I );
                    FreeSpline( mInfo );

                } else {

                    /*
                     *  Do full blown I integral. (Integrand is evaluated by tracing to required s-values.)
                     */
                    *I = Iinv( mInfo  );
                    if (mInfo->VerbosityLevel > 0) printf("Lgm_McIlwain_L: Integral Invariant, I (full integral): %g\n",  *I );

                }

            }

        }



        /*
         * Current time-dependant value of dipole moement (derived from first 3 vals of IGRF model)
         */
        *M = mInfo->c->M_cd;

        /*
         *  Bmirror value.
         */
        *Bm = mInfo->Bm;


        /*
         *  McIlwain L, via McIlwain's original tables or via Hilton approx.
         */
        if ( Type == 0 ) {
            L = LFromIBmM_McIlwain( *I, *Bm, *M );
        } else {
            L = LFromIBmM_Hilton( *I, *Bm, *M );
        }


    }


    return( L );

}


