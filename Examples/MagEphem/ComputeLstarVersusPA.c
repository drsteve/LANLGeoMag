#include <Lgm_MagModelInfo.h>
#include <Lgm_LstarInfo.h>
#include <Lgm_MagEphemInfo.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#define TRACE_TOL   1e-7
#define KP_DEFAULT  1

// For colorizing text ...
//int Colors[8] = { 26, 202, 77, 63, 185, 207, 124, 46 };
int Colors[9] = { 224, 209, 21, 46, 55, 104, 22, 185, 23 };



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
void ComputeLstarVersusPA( long int Date, double UTC, Lgm_Vector *u, int nAlpha, double *Alpha, int Quality, Lgm_MagEphemInfo *MagEphemInfo ) {

    Lgm_LstarInfo 	*LstarInfo, *LstarInfo2, *LstarInfo3;
    Lgm_Vector      v1, v2, v3, vv1, Bvec;
    double          sa, sa2, Blocal;
    double          Lam, CosLam, LSimple;
    int             i, LS_Flag, nn, tk, TraceFlag;
    char            *PreStr, *PostStr;

    LstarInfo = MagEphemInfo->LstarInfo;

    // Save Date, UTC to MagEphemInfo structure
    MagEphemInfo->Date   = Date;
    MagEphemInfo->UTC    = UTC;

    // Save nAlpha, and Alpha array to MagEphemInfo structure
    MagEphemInfo->nAlpha = nAlpha;
    for (i=0; i<MagEphemInfo->nAlpha; i++) MagEphemInfo->Alpha[i] = Alpha[i];

    // Set Tolerances
    SetLstarTolerances( MagEphemInfo->LstarQuality, LstarInfo );



    // set coord transformation 
    Lgm_Set_Coord_Transforms( Date, UTC, LstarInfo->mInfo->c );


    /*
     *  Blocal at sat location.
     */
    MagEphemInfo->P_gsm = *u;
    LstarInfo->mInfo->Bfield( u, &Bvec, LstarInfo->mInfo );
    Blocal = Lgm_Magnitude( &Bvec );
    MagEphemInfo->B = Blocal;


    /*
     *  Compute Field-related quantities for each Pitch Angle.
     */
    TraceFlag = Lgm_Trace( u, &v1, &v2, &v3, 120.0, TRACE_TOL, TRACE_TOL, LstarInfo->mInfo );
    MagEphemInfo->FieldLineType = TraceFlag;
    if ( TraceFlag > 0 ) {
        MagEphemInfo->Footprint_Ps  = v1;
        MagEphemInfo->Footprint_Pn  = v2;
        MagEphemInfo->Pmin_gsm      = v3;
        MagEphemInfo->Bmin          = LstarInfo->mInfo->Bmin;
        MagEphemInfo->Mref          = LstarInfo->mInfo->c->M_cd_McIllwain;
        MagEphemInfo->Mcurr         = LstarInfo->mInfo->c->M_cd;
        MagEphemInfo->Mused         = MagEphemInfo->Mcurr;
    }
    if ( TraceFlag == 1 ) {



            /*
             *  Get a simple measure of how big L is
             */
            Lgm_Convert_Coords( &v1, &vv1, GSM_TO_SM, LstarInfo->mInfo->c );
            Lam = asin( vv1.z/Lgm_Magnitude( &vv1 ) );
            CosLam = cos( Lam );
            LSimple = (1.0+120.0/WGS84_A)/( CosLam*CosLam );



            { // ***** BEGIN PARALLEL EXECUTION *****

            /*
             *  Do all of the PAs in parallel. To control how many threads get run
             *  use the enironment variable OMP_NUM_THREADS. For example,
             *          setenv OMP_NUM_THREADS 8
             *  will use 8 threads to do the loop in parallel. Be very careful what gets 
             *  set private here -- the threads must not interfere with each other.
             */
            #pragma omp parallel private(LstarInfo2,LstarInfo3,sa,sa2,LS_Flag,nn,tk,PreStr,PostStr)
            #pragma omp for schedule(dynamic, 1)
            for ( i=0; i<MagEphemInfo->nAlpha; i++ ){  // LOOP OVER PITCH ANGLES

                

                // make a local copy of LstarInfo structure -- needed for multi-threading
                LstarInfo3 = Lgm_CopyLstarInfo( LstarInfo );

                // colorize the diagnostic messages.
                sprintf( LstarInfo3->PreStr, "\033[38;5;%dm", Colors[i%9]); sprintf( LstarInfo3->PostStr, "\033[0m");
                PreStr = LstarInfo3->PreStr; PostStr = LstarInfo3->PostStr;

                /*
                 *  Set Pitch Angle, sin, sin^2, and Bmirror
                 */
                sa = sin( MagEphemInfo->Alpha[i]*RadPerDeg ); sa2 = sa*sa;

                LstarInfo3->mInfo->Bm = Blocal/sa2; 
                NewTimeLstarInfo( Date, UTC, MagEphemInfo->Alpha[i], LstarInfo3->mInfo->Bfield, LstarInfo3 );
                MagEphemInfo->Bm[i] = LstarInfo3->mInfo->Bm;

                /*
                 *  Compute L* 
                 */
                LstarInfo2 = Lgm_CopyLstarInfo( LstarInfo3 );
                if ( LSimple < LstarInfo2->LSimpleMax ){

                    LstarInfo2->mInfo->Bm = LstarInfo3->mInfo->Bm;
                    if (LstarInfo3->VerbosityLevel >= 2 ) {
                        printf("\n\n\t%sComputing L* for: UTC = %g PA = %d  (%g)%s\n", PreStr, UTC, i, MagEphemInfo->Alpha[i], PostStr );
                        printf("    \t%s                  I   = %g PA = %d  (%g)%s\n", PreStr, MagEphemInfo->I[i], i, MagEphemInfo->Alpha[i], PostStr );
                    }
                    LS_Flag = Lstar( &v3, LstarInfo2);
                    if (LstarInfo3->VerbosityLevel >= 2 ) {
                        printf("\t%sUTC, L*          = %g %g%s\n", PreStr, UTC, LstarInfo2->LS, PostStr );
                        printf("\t%sUTC, L*_McIlwain = %g %g%s\n", PreStr, UTC, LstarInfo2->LS_McIlwain_M, PostStr );
                        printf("\t%sUTC, LSimple     = %g %g%s\n\n\n", PreStr, UTC, LSimple, PostStr );
                    }
                    MagEphemInfo->Lstar[i] = LstarInfo2->LS;
                    MagEphemInfo->I[i] = LstarInfo2->I[i];

                    printf("%sL* for Pitch Angle: Alpha[%d] = %g Date: %ld   UTC: %g   Lsimple:%g   L*:%g%s\n", PreStr, i, MagEphemInfo->Alpha[i], Date, UTC, LSimple, LstarInfo2->LS, PostStr );

                    /*
                     * Save results to the MagEphemInfo structure.
                     */
                    MagEphemInfo->nShellPoints[i] = LstarInfo2->nPnts;
                    for (nn=0; nn<LstarInfo2->nPnts; nn++ ){
                        MagEphemInfo->ShellI[i][nn] = LstarInfo2->I[nn];
                        MagEphemInfo->ShellFootprint_Pn[i][nn] = LstarInfo2->Footprint_Pn[nn];
                        MagEphemInfo->ShellFootprint_Ps[i][nn] = LstarInfo2->Footprint_Ps[nn];
                        MagEphemInfo->ShellMirror_Pn[i][nn]    = LstarInfo2->Mirror_Pn[nn];
                        MagEphemInfo->ShellMirror_Ps[i][nn]    = LstarInfo2->Mirror_Ps[nn];
                        MagEphemInfo->ShellMirror_Ss[i][nn]    = LstarInfo2->mInfo->Sm_South;
                        MagEphemInfo->ShellMirror_Sn[i][nn]    = LstarInfo2->mInfo->Sm_North;

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
                    }

                } else {
                    printf(" Lsimple >= %g  ( Not doing L* calculation )\n", LstarInfo2->LSimpleMax );
                }

                FreeLstarInfo( LstarInfo2 );

                FreeLstarInfo( LstarInfo3 );
            }

        } 
        // ***** END PARALLEL EXECUTION *****


    }


    return;

}
