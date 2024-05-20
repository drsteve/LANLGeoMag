/*! \file  Lgm_ComputeLstarVersusPA.c
 *
 *  \brief Parallel routine for computing L* for multiple pitch angles (these are local pitch angles, not equatorial PA) at once. Uses OpenMP for parallelization.
 *
 *  
 */
#include "Lgm/Lgm_MagModelInfo.h"
#include "Lgm/Lgm_LstarInfo.h"
#include "Lgm/Lgm_MagEphemInfo.h"
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#if USE_OPENMP
#include <omp.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#define TRACE_TOL   1e-7
#define TRACE_TOL2   1e-10
#define KP_DEFAULT  1

// For colorizing text ...
//int Colors[8] = { 26, 202, 77, 63, 185, 207, 124, 46 };
int Colors[9] = { 224, 209, 21, 46, 55, 104, 22, 185, 23 };



/**
 *      Input Variables:
 *
 *                      Date:
 *                       UTC:
 *                         u:  Input position vector in GSM
 *                    nAlpha:  Number of Pitch Angles to compute
 *                     Alpha:  Pitch Angles to compute
 *
 *      Input/OutPut Variables:
 *
 *              MagEphemInfo:  Structure used to input and output parameters/settings/results to/from routine.
 *
 */
void Lgm_ComputeLstarVersusPA( long int Date, double UTC, Lgm_Vector *u, int nAlpha, double *Alpha, int Colorize, Lgm_MagEphemInfo *MagEphemInfo ) {

    Lgm_LstarInfo 	*LstarInfo, *LstarInfo2, *LstarInfo3;
    Lgm_Vector      v1, v2, v3, vv1, Bvec;
    double          sa, sa2, Blocal;
    double          Lam, CosLam, LSimple;
    int             i, k, LS_Flag, nn, tk, TraceFlag;
    int             nShabII, nShabI;
    char            *PreStr, *PostStr;

    LstarInfo = MagEphemInfo->LstarInfo;

    // Save Date, UTC to MagEphemInfo structure
    MagEphemInfo->Date   = Date;
    MagEphemInfo->UTC    = UTC;

    // Save nAlpha, and Alpha array to MagEphemInfo structure
    MagEphemInfo->nAlpha = nAlpha;
    for (i=0; i<MagEphemInfo->nAlpha; i++) MagEphemInfo->Alpha[i] = Alpha[i];

    // Set Tolerances
    Lgm_SetLstarTolerances( MagEphemInfo->LstarQuality, MagEphemInfo->nFLsInDriftShell, LstarInfo );


    // set coord transformation
    Lgm_Set_Coord_Transforms( Date, UTC, LstarInfo->mInfo->c );

    /*
     *  Blocal at sat location.
     */
    MagEphemInfo->P = *u;
    LstarInfo->mInfo->Bfield( u, &Bvec, LstarInfo->mInfo );
    Blocal = Lgm_Magnitude( &Bvec );
    MagEphemInfo->B = Blocal;
//printf("Blocal = %g\n", Blocal);

    for ( i=0; i<MagEphemInfo->nAlpha; i++ ){
        sa = sin( MagEphemInfo->Alpha[i]*RadPerDeg ); sa2 = sa*sa;
        MagEphemInfo->Bm[i] = Blocal/sa2;
    }



    /*
     *  Compute Field-related quantities for each Pitch Angle.
     *
     *
     *  First do a trace to identify the FL type and some of its critical points.
     *
     */
    TraceFlag = Lgm_Trace( u, &v1, &v2, &v3, LstarInfo->mInfo->Lgm_LossConeHeight, TRACE_TOL, TRACE_TOL2, LstarInfo->mInfo );
    MagEphemInfo->FieldLineType = TraceFlag;
//printf("P = %g %g %g   Blocal = %lf Bmin = %lf \n",  u->x, u->y, u->z, Blocal, MagEphemInfo->Bmin );
    if ( TraceFlag > 0 ) {
        MagEphemInfo->Ellipsoid_Footprint_Ps  = v1;
        MagEphemInfo->Ellipsoid_Footprint_Pn  = v2;
        MagEphemInfo->Pmin          = v3;
        MagEphemInfo->Snorth        = LstarInfo->mInfo->Snorth;
        MagEphemInfo->Ssouth        = LstarInfo->mInfo->Ssouth;
        MagEphemInfo->Smin          = LstarInfo->mInfo->Smin;
        MagEphemInfo->Bmin          = LstarInfo->mInfo->Bmin;
//printf("P = %g %g %g   Blocal = %lf Bmin = %lf \n",  u->x, u->y, u->z, Blocal, MagEphemInfo->Bmin );
        //MagEphemInfo->Mref          = LstarInfo->mInfo->c->M_cd_McIlwain;
        MagEphemInfo->Mref          = LstarInfo->mInfo->c->M_cd_2010;
        MagEphemInfo->Mcurr         = LstarInfo->mInfo->c->M_cd;
        MagEphemInfo->Mused         = MagEphemInfo->Mref;
    }
    if ( TraceFlag != 1 ) {
        /*
         * FL not closed.
         */
        for ( i=0; i<MagEphemInfo->nAlpha; i++ ){
            MagEphemInfo->Lstar[i] = LGM_FILL_VALUE;
            MagEphemInfo->I[i]     = LGM_FILL_VALUE;
            MagEphemInfo->K[i]     = LGM_FILL_VALUE;
            MagEphemInfo->Sb[i]    = LGM_FILL_VALUE;
        }
        MagEphemInfo->Sb0     = LGM_FILL_VALUE;
        MagEphemInfo->d2B_ds2 = LGM_FILL_VALUE;
        MagEphemInfo->RofC    = LGM_FILL_VALUE;
    }
    if ( TraceFlag == 1 ) {

           MagEphemInfo->Sb0     = LstarInfo->mInfo->Sb0;     // Sb Integral for equatorially mirroring particles.
           MagEphemInfo->d2B_ds2 = LstarInfo->mInfo->d2B_ds2; // second deriv of B wrt s at equator.
           MagEphemInfo->RofC    = LstarInfo->mInfo->RofC;    // radius of curvature at Bmin point.

            /*
             *  Get a simple measure of how big L is
             */
            Lgm_Convert_Coords( &v1, &vv1, GSM_TO_SM, LstarInfo->mInfo->c );
            Lam = asin( vv1.z/Lgm_Magnitude( &vv1 ) );
            CosLam = cos( Lam );
            //LSimple = (1.0+LstarInfo->mInfo->Lgm_LossConeHeight/WGS84_A)/( CosLam*CosLam );
            LSimple = Lgm_Magnitude( &LstarInfo->mInfo->Pmin );

            { // ***** BEGIN PARALLEL EXECUTION *****

            /*
             *  Do all of the PAs in parallel. To control how many threads get run
             *  use the enironment variable OMP_NUM_THREADS. For example,
             *          setenv OMP_NUM_THREADS 8
             *  will use 8 threads to do the loop in parallel. Be very careful what gets
             *  set private here -- the threads must not interfere with each other.
             */
#if USE_OPENMP
            #pragma omp parallel private(LstarInfo2,LstarInfo3,sa,sa2,LS_Flag,nn,tk,PreStr,PostStr,nShabII,nShabI)
            #pragma omp for schedule(dynamic, 1)
#endif
            for ( i=0; i<MagEphemInfo->nAlpha; i++ ){  // LOOP OVER PITCH ANGLES

                /*
                 * make a local copy of LstarInfo structure -- needed for multi-threading
                 */
                LstarInfo3 = Lgm_CopyLstarInfo( LstarInfo );

                /*
                 * colorize the diagnostic messages.
                 */
                if ( Colorize ){
                    sprintf( LstarInfo3->PreStr, "\033[38;5;%dm", Colors[i%9]); sprintf( LstarInfo3->PostStr, "\033[0m");
                } else {
                    LstarInfo3->PreStr[0] = '\0'; LstarInfo3->PostStr[0] = '\0';
                }
                PreStr = LstarInfo3->PreStr; PostStr = LstarInfo3->PostStr;

                /*
                 *  Set Bmirror
                 */
                LstarInfo3->mInfo->Bm = MagEphemInfo->Bm[i];
                NewTimeLstarInfo( Date, UTC, MagEphemInfo->Alpha[i], LstarInfo3->mInfo->Bfield, LstarInfo3 );

                /*
                 *  Compute L*
                 */
                if ( LSimple < LstarInfo3->LSimpleMax ){

                    LstarInfo2 = Lgm_CopyLstarInfo( LstarInfo3 );

                    LstarInfo2->mInfo->Bm = LstarInfo3->mInfo->Bm;
                    if (LstarInfo3->VerbosityLevel >= 2 ) {
                        printf("\n\n\t\t%sComputing L* for: UTC = %g  (Local) PA = %d  (%g)%s\n", PreStr, UTC, i, MagEphemInfo->Alpha[i], PostStr );
                        //printf("    \t\t%s                  I   = %g PA = %d  (%g)%s\n", PreStr, MagEphemInfo->I[i], i, MagEphemInfo->Alpha[i], PostStr );
                    }
//////////////////////////////NOTE
//////////////////////////////NOTE   We are giving Lstar the Min B point ALREADY!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//////////////////////////////NOTE

                    LS_Flag = Lstar( &v3, LstarInfo2);

                    if (LstarInfo3->VerbosityLevel >= 2 ) {
                        printf("\t\t%sUTC, L*          = %g %g%s\n", PreStr, UTC, LstarInfo2->LS, PostStr );
                        printf("\t\t%sUTC, L*_McIlwain = %g %g%s\n", PreStr, UTC, LstarInfo2->LS_McIlwain_M, PostStr );
                        printf("\t\t%sUTC, LSimple     = %g %g%s\n\n\n", PreStr, UTC, LSimple, PostStr );
                    }
                    MagEphemInfo->Lstar[i] = ( LS_Flag >= 0 ) ? LstarInfo2->LS : LGM_FILL_VALUE;
                    MagEphemInfo->I[i]  = LstarInfo2->I[0]; // I[0] is I for the FL that the sat is on.
                    MagEphemInfo->K[i]  = LstarInfo2->I[0]*sqrt(MagEphemInfo->Bm[i]*1e-5); // Second invariant
                    MagEphemInfo->Sb[i] = LstarInfo2->SbIntegral0; // SbIntegral0 is Sb for the FL that the sat is on.


                    /*
                     *  Determine the type of the orbit
                     */
                    nShabII = nShabI = 0;
                    if ( LS_Flag >= 0 )  {

                        LstarInfo2->DriftOrbitType = LGM_DRIFT_ORBIT_CLOSED;
                        for ( nn=0; nn<LstarInfo2->nPnts; ++nn ) {
                            if  ( LstarInfo2->nMinima[nn] > 1 ) {
                                if ( LstarInfo2->nBounceRegions[nn] > 1 ) {
                                    nShabII++;
                                } else {
                                    nShabI++;
                                }
                            }
                        }
                        if (nShabII > 0) {
                            LstarInfo2->DriftOrbitType = LGM_DRIFT_ORBIT_CLOSED_SHABANSKY_II;
                        } else if (nShabI > 0) {
                            LstarInfo2->DriftOrbitType = LGM_DRIFT_ORBIT_CLOSED_SHABANSKY_I;
                        }

                    } else {

                        LstarInfo2->DriftOrbitType = LGM_DRIFT_ORBIT_CLOSED;
                        for ( nn=0; nn<LstarInfo2->nPnts; ++nn ) {
                            if  ( LstarInfo2->nMinima[nn] > 1 ) {
                                if ( LstarInfo2->nBounceRegions[nn] > 1 ) {
                                    nShabII++;
                                } else {
                                    nShabI++;
                                }
                            }
                        }
                        if (nShabII > 0) {
                            LstarInfo2->DriftOrbitType = LGM_DRIFT_ORBIT_OPEN_SHABANSKY_II;
                        } else if (nShabI > 0) {
                            LstarInfo2->DriftOrbitType = LGM_DRIFT_ORBIT_OPEN_SHABANSKY_I;
                        }

                    }
                    MagEphemInfo->DriftOrbitType[i] = LstarInfo2->DriftOrbitType;


                    //printf("\t    %sL* [ %g Deg. ]: Date: %ld   UTC: %g   Lsimple:%g   L*:%.15g%s\n", PreStr, MagEphemInfo->Alpha[i], Date, UTC, LSimple, LstarInfo2->LS, PostStr );
                    if (LstarInfo3->VerbosityLevel > 0 ) {
		                printf("\t    %sL* [ %g\u00b0 ]: Date: %ld   UTC: %g   Lsimple:%g   L*:%.15g%s\n", PreStr, MagEphemInfo->Alpha[i], Date, UTC, LSimple, LstarInfo2->LS, PostStr );
		            }




                    /*
                     * Save detailed results to the MagEphemInfo structure.
                     */
                    MagEphemInfo->nShellPoints[i] = LstarInfo2->nPnts;
                    for (nn=0; nn<LstarInfo2->nPnts; nn++ ){

                        // This Pmin does not seem to be the right one when we have Shabansky
                        MagEphemInfo->Shell_Pmin[i][nn]  = LstarInfo2->Pmin[nn];

                        MagEphemInfo->Shell_Bmin[i][nn]  = LstarInfo2->Bmin[nn];
                        MagEphemInfo->Shell_GradI[i][nn] = LstarInfo2->GradI[nn];
                        MagEphemInfo->Shell_Vgc[i][nn]   = LstarInfo2->Vgc[nn];


// This does really capture the I/2 bit?
MagEphemInfo->ShellI[i][nn] = LstarInfo2->I[nn];
MagEphemInfo->nBounceRegions[i][nn] = LstarInfo2->nBounceRegions[nn];

                        MagEphemInfo->ShellSphericalFootprint_Pn[i][nn] = LstarInfo2->Spherical_Footprint_Pn[nn];
                        MagEphemInfo->ShellSphericalFootprint_Sn[i][nn] = LstarInfo2->Spherical_Footprint_Sn[nn];
                        MagEphemInfo->ShellSphericalFootprint_Bn[i][nn] = LstarInfo2->Spherical_Footprint_Bn[nn];
                        MagEphemInfo->ShellSphericalFootprint_Ps[i][nn] = LstarInfo2->Spherical_Footprint_Ps[nn];
                        MagEphemInfo->ShellSphericalFootprint_Ss[i][nn] = LstarInfo2->Spherical_Footprint_Ss[nn];
                        MagEphemInfo->ShellSphericalFootprint_Bs[i][nn] = LstarInfo2->Spherical_Footprint_Bs[nn];

                        MagEphemInfo->ShellEllipsoidFootprint_Pn[i][nn] = LstarInfo2->Ellipsoid_Footprint_Pn[nn];
                        MagEphemInfo->ShellEllipsoidFootprint_Sn[i][nn] = LstarInfo2->Ellipsoid_Footprint_Sn[nn];
                        MagEphemInfo->ShellEllipsoidFootprint_Bn[i][nn] = LstarInfo2->Ellipsoid_Footprint_Bn[nn];
                        MagEphemInfo->ShellEllipsoidFootprint_Ps[i][nn] = LstarInfo2->Ellipsoid_Footprint_Ps[nn];
                        MagEphemInfo->ShellEllipsoidFootprint_Ss[i][nn] = LstarInfo2->Ellipsoid_Footprint_Ss[nn];
                        MagEphemInfo->ShellEllipsoidFootprint_Bs[i][nn] = LstarInfo2->Ellipsoid_Footprint_Bs[nn];


                        MagEphemInfo->ShellMirror_Pn[i][nn]    = LstarInfo2->Mirror_Pn[nn];
                        //MagEphemInfo->ShellMirror_Sn[i][nn]    = LstarInfo2->mInfo->Sm_North;
                        MagEphemInfo->ShellMirror_Sn[i][nn]    = LstarInfo2->Mirror_Sn[nn];

                        MagEphemInfo->ShellMirror_Ps[i][nn]    = LstarInfo2->Mirror_Ps[nn];
                        //MagEphemInfo->ShellMirror_Ss[i][nn]    = LstarInfo2->mInfo->Sm_South;
                        MagEphemInfo->ShellMirror_Ss[i][nn]    = LstarInfo2->Mirror_Ss[nn];

                        /*
                         *  Save all of the drift shell FLs in MagEphemInfo structure
                         */
                        MagEphemInfo->nFieldPnts[i][nn] = LstarInfo2->nFieldPnts[nn];
                        for (tk=0; tk<LstarInfo2->nFieldPnts[nn]; tk++){    // loop over points in a FL
                            MagEphemInfo->s_gsm[i][nn][tk] = LstarInfo2->s_gsm[nn][tk];
                            MagEphemInfo->Bmag[i][nn][tk]  = LstarInfo2->Bmag[nn][tk];
                            MagEphemInfo->x_gsm[i][nn][tk] = LstarInfo2->x_gsm[nn][tk];
                            MagEphemInfo->y_gsm[i][nn][tk] = LstarInfo2->y_gsm[nn][tk];
                            MagEphemInfo->z_gsm[i][nn][tk] = LstarInfo2->z_gsm[nn][tk];
                        }
//printf("LstarInfo2->nFieldPnts[%d] = %d\n", nn, LstarInfo2->nFieldPnts[nn]);


                        MagEphemInfo->nMinima[i][nn] = LstarInfo2->nMinima[nn];
                        MagEphemInfo->nMaxima[i][nn] = LstarInfo2->nMaxima[nn];

                    }

                    FreeLstarInfo( LstarInfo2 );

                } else {
                    printf(" Lsimple >= %g  ( Not doing L* calculation )\n", LstarInfo3->LSimpleMax );
                    MagEphemInfo->Lstar[i] = LGM_FILL_VALUE;
//printf("Bm = %g\n", MagEphemInfo->Bm[i]);
                    MagEphemInfo->I[i]     = LGM_FILL_VALUE;
                    MagEphemInfo->K[i]     = LGM_FILL_VALUE;
                    MagEphemInfo->nShellPoints[i] = 0;

                }


                FreeLstarInfo( LstarInfo3 );
            }

        }
        // ***** END PARALLEL EXECUTION *****


    }


    return;

}
