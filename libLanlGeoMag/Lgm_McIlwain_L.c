#include "Lgm/Lgm_MagModelInfo.h"                                                                                                                                                                                
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>


/*
 *   Lgm_McIlwain_L
 *   -------------
 */


//! Compute McIlwain L-shell parameter for a given date, time, location and pitch angle.
/**
 *            \param[in]        Date        Date in format (e.g. 20101231). 
 *            \param[in]        UTC         Universal Time (Coordinated) in decimal hours (e.g. 23.5).
 *            \param[in]        u           Position (in GSM) to compute L-shell.
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

    int             reset;
    Lgm_Vector      v1, v2, v3, Bvec, Bvectmp, Ptmp, u_scale;
    double          rat, B, sa, sa2, Blocal, dSa, dSb, r, SS, L, stmp, Hdid, Hnext, Btmp;

    u_scale.x =  100.0;  u_scale.y = 100.0; u_scale.z = 100.0;

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
    mInfo->Blocal = Blocal;

    /*
     * Set Pitch Angle, sin(Alpha), sin^2(Alpha), and Bmirror
     */
    sa = sin( Alpha*RadPerDeg ); sa2 = sa*sa;
    mInfo->Bm = Blocal/sa2;



    /*
     *  First do a trace to identify the FL type and some of its critical points.
     */
    if ( Lgm_Trace( u, &v1, &v2, &v3, mInfo->Lgm_LossConeHeight, mInfo->Lgm_TraceToEarth_Tol, mInfo->Lgm_TraceToBmin_Tol, mInfo ) == LGM_CLOSED ) {


        /*
         *  Test to see if the S/C is already close to the Bmin point or the
         *  B's are almost the same; And the Pitch angle is close to 90.  If
         *  so, use an approximation to I.
         */
        if ( ( ((SS=Lgm_VecDiffMag( u, &v3 )) < 1e-4) || (fabs( mInfo->Blocal - mInfo->Bmin) < 1e-2) ) && (fabs(90.0-Alpha) < 1e-2)  ) {

            // if FL length is small, use an approx expression for I
            rat = mInfo->Bmin/mInfo->Bm;
            if ((1.0-rat) < 0.0) {
                *I = 0.0;
            } else {
                // Eqn 2.66b in Roederer
                *I = SS*sqrt(1.0 - rat);
            }

        } else {


            /*
             * Trace from Bmin point up to northern mirror point and down to
             * southern mirror point. dSa and dSb are the distances along the
             * FL from the starting points point. So dSa is from the Bmin value.
             * And dSb is from the Psouth value.
             */
            if ( Lgm_TraceToMirrorPoint( &(mInfo->Pmin), &(mInfo->Pm_South), &dSa, mInfo->Bm, -1.0, mInfo->Lgm_TraceToMirrorPoint_Tol, mInfo ) >= 0 ) {

                if (mInfo->VerbosityLevel > 0) {
                    printf("\n\tMin-B  Point Location, Pmin (Re):      < %g, %g, %g >\n", mInfo->Pmin.x, mInfo->Pmin.y, mInfo->Pmin.z );
                    printf("\tMirror Point Location, Pm_South (Re):      < %g, %g, %g >  |Pm_South| = %g\n", mInfo->Pm_South.x, mInfo->Pm_South.y, mInfo->Pm_South.z, Lgm_Magnitude(&mInfo->Pm_South) );
                    mInfo->Bfield( &mInfo->Pm_South, &Bvec, mInfo );
                    B = Lgm_Magnitude( &Bvec );
                    printf("\tMag. Field Strength, Bm at Pm_South (nT):  %g     (mInfo->Bm = %g)\n", B, mInfo->Bm );
                }



                if ( Lgm_TraceToMirrorPoint( &(mInfo->Pmin), &(mInfo->Pm_North), &dSb, mInfo->Bm,  1.0, mInfo->Lgm_TraceToMirrorPoint_Tol, mInfo ) >= 0 ) {

                    if (mInfo->VerbosityLevel > 0) {
                        printf("\n\tMin-B  Point Location, Pmin (Re):      < %g, %g, %g >\n", mInfo->Pmin.x, mInfo->Pmin.y, mInfo->Pmin.z );
                        printf("\tMirror Point Location, Pm_North (Re):      < %g, %g, %g >  |Pm_North| = %g\n", mInfo->Pm_North.x, mInfo->Pm_North.y, mInfo->Pm_North.z, Lgm_Magnitude(&mInfo->Pm_North) );
                        mInfo->Bfield( &mInfo->Pm_North, &Bvec, mInfo );
                        B = Lgm_Magnitude( &Bvec );
                        printf("\tMag. Field Strength, Bm at Pm_North (nT):  %g     (mInfo->Bm = %g)\n", B, mInfo->Bm );
                    }

                    /*
                     *  Set the limits of integration. Define s=0 at the sourthern mirror point. Then, sm_North will just be dSb
                     */
                    //SS = dSb;
                    SS = dSa+dSb;
                    mInfo->Hmax = SS/(double)mInfo->nDivs;
                    if ( mInfo->Hmax > mInfo->MaxDiv ) mInfo->Hmax = mInfo->MaxDiv;
                    r  = Lgm_Magnitude( &mInfo->Pm_North );
                    mInfo->Sm_South = 0.0;
                    mInfo->Sm_North = SS;

                    if ( SS <= 1e-5 ) {

                        // if FL length is small, use an approx expression for I
                        rat = mInfo->Bmin/mInfo->Bm;
                        if ((1.0-rat) < 0.0) {
                            *I = 0.0;
                        } else {
                            // Eqn 2.66b in Roederer
                            *I = SS*sqrt(1.0 - rat);
                        }

                    } else if ( mInfo->UseInterpRoutines ) {
                        if ( Lgm_TraceLine2( &(mInfo->Pm_South), &mInfo->Pm_North, (r-1.0)*Re, 0.5*SS-mInfo->Hmax, 1.0, mInfo->Lgm_TraceToEarth_Tol, FALSE, mInfo ) < 0 ) return(-9e99);
//printf("BEFORE mInfo->nPnts = %d   mInfo->s[0] = %g   mInfo->s[1] = %g    mInfo->s[mInfo->nPnts-2] = %g   mInfo->s[mInfo->nPnts-1] = %g    SS = %g\n", mInfo->nPnts, mInfo->s[0], mInfo->s[1], mInfo->s[mInfo->nPnts-2], mInfo->s[mInfo->nPnts-1], SS );
                        ReplaceFirstPoint( 0.0, mInfo->Bm, &mInfo->Pm_South, mInfo );
                        ReplaceLastPoint( SS, mInfo->Bm, &mInfo->Pm_North, mInfo );
//printf("AFTER1 mInfo->nPnts = %d   mInfo->s[0] = %g   mInfo->s[1] = %g    mInfo->s[mInfo->nPnts-2] = %g   mInfo->s[mInfo->nPnts-1] = %g    SS = %g\n", mInfo->nPnts, mInfo->s[0], mInfo->s[1], mInfo->s[mInfo->nPnts-2], mInfo->s[mInfo->nPnts-1], SS );

                        /*
                         * Make sure we have a small margin before and after so
                         * we dont end up trying to extrapolate if s ever getd
                         * slightly out of bounds.
                         */
                        Ptmp = mInfo->Pm_South; stmp = 0.0; reset = FALSE;
                        if ( Lgm_MagStep( &Ptmp, &u_scale, 0.01, &Hdid, &Hnext, -1.0, &stmp, &reset, mInfo->Bfield, mInfo ) < 0 ) { return(-1); }
                        mInfo->Bfield( &Ptmp, &Bvectmp, mInfo ); Btmp = Lgm_Magnitude( &Bvectmp );
                        //printf("-stmp, Btmp = %g %g\n", -stmp, Btmp );
                        AddNewPoint( -stmp, Btmp, &Ptmp, mInfo );
//printf("AFTER2 mInfo->nPnts = %d   mInfo->s[0] = %g   mInfo->s[1] = %g    mInfo->s[mInfo->nPnts-2] = %g   mInfo->s[mInfo->nPnts-1] = %g    SS = %g\n", mInfo->nPnts, mInfo->s[0], mInfo->s[1], mInfo->s[mInfo->nPnts-2], mInfo->s[mInfo->nPnts-1], SS );

                        Ptmp = mInfo->Pm_North; stmp = 0.0; reset = FALSE;
                        if ( Lgm_MagStep( &Ptmp, &u_scale, 0.01, &Hdid, &Hnext, 1.0, &stmp, &reset, mInfo->Bfield, mInfo ) < 0 ) { return(-1); }
                        mInfo->Bfield( &Ptmp, &Bvectmp, mInfo ); Btmp = Lgm_Magnitude( &Bvectmp );
                        //printf("stmp, Btmp = %g %g\n", SS+stmp, Btmp );
                        AddNewPoint( SS+stmp, Btmp, &Ptmp, mInfo );

                        
                        
                        

//printf("AFTER2 mInfo->nPnts = %d   mInfo->s[0] = %g   mInfo->s[1] = %g    mInfo->s[mInfo->nPnts-2] = %g   mInfo->s[mInfo->nPnts-1] = %g    SS = %g\n", mInfo->nPnts, mInfo->s[0], mInfo->s[1], mInfo->s[mInfo->nPnts-2], mInfo->s[mInfo->nPnts-1], SS );
                        //AddNewPoint( SS,  mInfo->Bm, &mInfo->Pm_North, mInfo );
                        if ( InitSpline( mInfo ) ) {

                            /*
                             *  Do interped I integral.
                             */
                            *I = Iinv_interped( mInfo  );
                            if (mInfo->VerbosityLevel > 0) printf("Lgm_McIlwain_L: Integral Invariant, I (interped):      %g\n",  *I );
                            FreeSpline( mInfo );

                        } else {

                            *I = -9e99;

                        }

                    } else {

                        /*
                         *  Do full blown I integral. (Integrand is evaluated by tracing to required s-values.)
                         */
                        *I = Iinv( mInfo  );
                        if (mInfo->VerbosityLevel > 0) printf("Lgm_McIlwain_L: Integral Invariant, I (full integral): %g\n",  *I );

                    }

                } else {
                    if (mInfo->VerbosityLevel > 0) printf("Could not find northern mirror point.\n");
                }

            } else {
                if (mInfo->VerbosityLevel > 0) printf("Could not find southern mirror point.\n");
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
        if ( *I < 0.0 ){
            L = -9e99;
        } else if ( Type == 0 ) {
            L = LFromIBmM_McIlwain( *I, *Bm, *M );
        } else {
            L = LFromIBmM_Hilton( *I, *Bm, *M );
        }


    }


    return( L );

}


