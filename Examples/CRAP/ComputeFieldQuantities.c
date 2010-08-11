#include <Lgm_MagModelInfo.h>
#include <Lgm_LstarInfo.h>
#include "MagEphemInfo.h"
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#define TRACE_TOL   1e-7
#define KP_DEFAULT  1


/*
 *  
 */
void ComputeFieldLineQuantities( long int Date, double UTC, Lgm_Vector *u, int nAlpha, double *Alpha, _MagEphemInfo *MagEphemInfo ) {

    _LstarInfo 	*LstarInfo, *LstarInfo2, *LstarInfo3;
    Lgm_Vector  v, v1, v2, v3, vv1, Bvec;
    double      sa, sa2, Blocal;
    double      Lam, CosLam, LSimple, dSa, dSb;
    int         i, LS_Flag, nn, tk;
    double      Ival, Sbval, r, SS;
    int         Kp=KP_DEFAULT;;
    int         DoPitchAngle[1000];
    char        Filename[128];
    FILE    *fpout;


    // Save Date, UTC to MagEphemInfo structure
    MagEphemInfo->Date   = Date;
    MagEphemInfo->UTC    = UTC;

    // Save nAlpha, and Alpha array to MagEphemInfo structure
    MagEphemInfo->nAlpha = nAlpha;
    for (i=0; i<MagEphemInfo->nAlpha; i++) MagEphemInfo->Alpha[i] = Alpha[i];

    // Use interpolated integrands for calculating I integrals (speeds things up).
    MagEphemInfo->UseInterpRoutines = FALSE;



    /*
     *  Initialize the LstarInfo structure
     */
    LstarInfo = InitLstarInfo( 0 );
    LstarInfo->VerbosityLevel   = 0;
    LstarInfo->SaveShellLines   = TRUE;

    LstarInfo->mInfo->UseInterpRoutines = FALSE;
    LstarInfo->mInfo->UseInterpRoutines = TRUE;
    LstarInfo->mInfo->Lgm_I_integrand_JumpMethod = LGM_ABSOLUTE_JUMP_METHOD;
    LstarInfo->mInfo->Lgm_MagFlux_Integrator_epsabs = 0.0;
    LstarInfo->mInfo->Lgm_MagFlux_Integrator_epsrel = 1e-5;
//LstarInfo->mInfo->Lgm_MagFlux_Integrator_epsrel = 1e-8;

    LstarInfo->mInfo->Lgm_LambdaIntegral_Integrator_epsabs = 0.0;
    LstarInfo->mInfo->Lgm_LambdaIntegral_Integrator_epsrel = 1e-5;
//LstarInfo->mInfo->Lgm_LambdaIntegral_Integrator_epsrel = 1e-8;

    LstarInfo->mInfo->Lgm_I_Integrator = DQAGS;
    LstarInfo->mInfo->Lgm_I_Integrator_epsrel = 0.0;
    LstarInfo->mInfo->Lgm_I_Integrator_epsabs = 1e-5;
    LstarInfo->mInfo->Lgm_I_Integrator_epsabs = 1e-7;
//LstarInfo->mInfo->Lgm_I_Integrator_epsabs = 1e-8;

    LstarInfo->mInfo->Lgm_Sb_Integrator = DQAGP; // not changeable (yet...)
    LstarInfo->mInfo->Lgm_Sb_Integrator_epsrel = 0.0;
    LstarInfo->mInfo->Lgm_Sb_Integrator_epsabs = 1e-3;

    LstarInfo->mInfo->Lgm_FindBmRadius_Tol = 1e-10;
    LstarInfo->mInfo->Lgm_TraceToMirrorPoint_Tol = 1e-10;
    LstarInfo->mInfo->Lgm_FindShellLine_I_Tol = 1e-3;
LstarInfo->mInfo->Lgm_FindShellLine_I_Tol = 1e-6;
LstarInfo->mInfo->Lgm_FindShellLine_I_Tol = 1e-8;
LstarInfo->mInfo->Lgm_FindShellLine_I_Tol = 1e-7;
LstarInfo->mInfo->Lgm_FindShellLine_I_Tol = 1e-3;
LstarInfo->mInfo->Lgm_FindShellLine_I_Tol = 1e-4;
LstarInfo->mInfo->Lgm_FindShellLine_I_Tol = 1e-5;
LstarInfo->mInfo->Lgm_FindShellLine_I_Tol = 1e-8;
LstarInfo->mInfo->Lgm_FindShellLine_I_Tol = 1e-6;
    sprintf( Filename, "results_1.5_%.0e.dat", LstarInfo->mInfo->Lgm_FindShellLine_I_Tol );
    fpout = fopen(Filename, "w");


    // set coord transformation 
    Lgm_Set_Coord_Transforms( Date, UTC, LstarInfo->mInfo->c );





    /*
     * Set model for tracing
     * THIS ALSO SHOULD COME IN ALREADY SET.
     */
    Kp = 1;
    LstarInfo->mInfo->Bfield = Lgm_B_T89;
    LstarInfo->mInfo->Bfield = Lgm_B_cdip;
    LstarInfo->mInfo->InternalModel = LGM_CDIP;
    LstarInfo->mInfo->Kp = ( Kp >= 0 ) ? Kp : KP_DEFAULT;
    if ( LstarInfo->mInfo->Kp > 5 ) LstarInfo->mInfo->Kp = 5;




    /*
     *  Init sat location.
     */
    MagEphemInfo->P_gsm = *u;
    LstarInfo->mInfo->Bfield( u, &Bvec, LstarInfo->mInfo );
    Blocal = Lgm_Magnitude( &Bvec );
    MagEphemInfo->B = Blocal;
    LstarInfo->mInfo->VerbosityLevel = 0;


    /*
     *  Compute Field-related quantities for each Pitch Angle.
     */
    if ( Lgm_Trace( u, &v1, &v2, &v3, 120.0, 0.01, TRACE_TOL, LstarInfo->mInfo ) == 1 ) {

        MagEphemInfo->Pmin_gsm = v3;
        MagEphemInfo->Bmin     = LstarInfo->mInfo->Bmin;

        /*
WHY DO WE DO THIS HERE?
         * Trace out FL from southern footprint to northern footprint. This
         * call saves the FL points in a form that can be interpolated.
         */
        LstarInfo->mInfo->Hmax = 0.5;
        //Lgm_TraceLine( &v1, &v, 120.0, 1.0, TRACE_TOL, TRUE, LstarInfo->mInfo );
        Lgm_TraceLine( &v1, &v, 120.0, 1.0, TRACE_TOL, FALSE, LstarInfo->mInfo );
        InitSpline( LstarInfo->mInfo );
        LstarInfo->mInfo->Hmax = 10.0;
        printf("LstarInfo->mInfo->Bmin = %g\n", LstarInfo->mInfo->Bmin );
        printf("Pmin = %g %g %g\n", LstarInfo->mInfo->Pmin.x, LstarInfo->mInfo->Pmin.y, LstarInfo->mInfo->Pmin.z );
FreeSpline( LstarInfo->mInfo );



            /*
             *  Get a simple measure of how big L is
             */
            Lgm_Convert_Coords( &v1, &vv1, GSM_TO_SM, LstarInfo->mInfo->c );
            Lam = asin( vv1.z/Lgm_Magnitude( &vv1 ) );
            CosLam = cos( Lam );
            LSimple = 1.0/( CosLam*CosLam );



            { // ***** BEGIN PARALLEL EXECUTION *****

            /*
             *  Do all of the PAs in parallel. To control how many threads get run
             *  use the enironment variable OMP_NUM_THREADS. For example,
             *          setenv OMP_NUM_THREADS 8
             *  will use 8 threads to do the loop in parallel. Be very careful what gets 
             *  set private here -- the threads must not interfere with each other.
             */
            #pragma omp parallel private(LstarInfo2,LstarInfo3,sa,sa2,dSa,dSb,LS_Flag,Ival,Sbval,nn,tk)
            #pragma omp for schedule(dynamic, 1)
            for ( i=0; i<MagEphemInfo->nAlpha; i++ ){  // LOOP OVER PITCH ANGLES

                LstarInfo3 = Lgm_CopyLstarInfo( LstarInfo );

                /*
                 *  Set Pitch Angle, sin, sin^2, and Bmirror
                 */
                sa = sin( MagEphemInfo->Alpha[i]*RadPerDeg ); sa2 = sa*sa;

                LstarInfo3->mInfo->Bm = Blocal/sa2; 
                NewTimeLstarInfo( Date, UTC, MagEphemInfo->Alpha[i], LstarInfo3->mInfo->Bfield, LstarInfo3 );
                MagEphemInfo->Bm[i] = LstarInfo3->mInfo->Bm;

printf("\033[38;5;%dm i=%d, sa2 = %g, LstarInfo3->mInfo->c->Lgm_IGRF_g[9][9] = %g\033[0m\n", i, i, sa2, LstarInfo3->mInfo->c->Lgm_IGRF_g[9][9] );

                /*
                 * Trace from Bmin point up to northern mirror point and down to
                 * southern mirror point. dS is the distances along
                 * (i.e. up and down) the FL from the starting point.
                 */
                DoPitchAngle[i] = FALSE;
                if ( Lgm_TraceToMirrorPoint( &(LstarInfo3->mInfo->Pmin), &(LstarInfo3->mInfo->Pm_South), &dSa, 120.0, LstarInfo3->mInfo->Bm, -1.0, TRACE_TOL, LstarInfo3->mInfo ) > 0 ) {
                    if ( Lgm_TraceToMirrorPoint( &(LstarInfo3->mInfo->Pm_South), &(LstarInfo3->mInfo->Pm_North), &dSb, 120.0, LstarInfo3->mInfo->Bm,  1.0, TRACE_TOL, LstarInfo3->mInfo ) > 0 ) {

                       DoPitchAngle[i] = TRUE;


                        /*
                         *  Save mirror points to MagEphemInfo
                         */
                        MagEphemInfo->Pms_gsm[i] = LstarInfo3->mInfo->Pm_South;
                        MagEphemInfo->Pmn_gsm[i] = LstarInfo3->mInfo->Pm_North;





                        /*
                         *  Set the limits of integration. Also set tolerances for
                         *  Quadpack routines.
                         */
                        //LstarInfo3->mInfo->Sm_South = LstarInfo3->mInfo->smin - dSa;
                        //LstarInfo3->mInfo->Sm_North = LstarInfo3->mInfo->Sm_South + dSb;
SS = dSb;
LstarInfo3->mInfo->Hmax = SS/200.0;
r  = Lgm_Magnitude( &LstarInfo3->mInfo->Pm_North );
LstarInfo3->mInfo->Sm_South = 0.0;
LstarInfo3->mInfo->Sm_North = SS;

                        



                        /*
                         *  Since the Pm_South and Pm_north points are not likely
                         *  to be in the interped array, add them (since we have
                         *  them now anyway!). This is particularly important to do
                         *  for the Sb integrals becuase the integrand is high near
                         *  the mirror points.  If we dont do it, we could get
                         *  slight differences due to it.
                         */
//                        AddNewPoint( LstarInfo3->mInfo->Sm_South, LstarInfo3->mInfo->Bm, &LstarInfo3->mInfo->Pm_South, LstarInfo3->mInfo );
//                        AddNewPoint( LstarInfo3->mInfo->Sm_North, LstarInfo3->mInfo->Bm, &LstarInfo3->mInfo->Pm_North, LstarInfo3->mInfo );

                        LstarInfo3->VerbosityLevel = 3;
                        LstarInfo3->VerbosityLevel = 0;
                        if ( MagEphemInfo->UseInterpRoutines ) {

                            Lgm_TraceLine2( &(LstarInfo3->mInfo->Pm_South), &LstarInfo3->mInfo->Pm_North, (r-1.0)*Re, 0.5*SS-LstarInfo3->mInfo->Hmax, 1.0, 1e-7, FALSE, LstarInfo3->mInfo );
                            ReplaceFirstPoint( 0.0, LstarInfo3->mInfo->Bm, &LstarInfo3->mInfo->Pm_South, LstarInfo3->mInfo );
                            AddNewPoint( SS,  LstarInfo3->mInfo->Bm, &LstarInfo3->mInfo->Pm_North, LstarInfo3->mInfo );
                            InitSpline( LstarInfo3->mInfo );


                            /*
                             *  Do interped I integral.
                             */
                            Ival = Iinv_interped( LstarInfo3->mInfo  );
                            if (LstarInfo3->VerbosityLevel >= 2 ) printf("Iinv (Interped Integral) = %g\n",  Ival );

                            /*
                             *  Do interped Sb integral.
                             */
                            Sbval = SbIntegral_interped( LstarInfo3->mInfo  );
                            if (LstarInfo3->VerbosityLevel >= 2 ) printf("Sb (Interped Integral) = %g\n", Sbval );

                            FreeSpline( LstarInfo3->mInfo );


                        } else {

                            /*
                             *  Do full blown I integral. 
                             */
                            Ival = Iinv( LstarInfo3->mInfo  );
                            if (LstarInfo3->VerbosityLevel >= 2 ) printf("Iinv (Full Integral) = %g\n",  Ival );

                            /*
                             *  Do full blown Sb integral. 
                             */
                            Sbval = SbIntegral( LstarInfo3->mInfo  );
                            if (LstarInfo3->VerbosityLevel >= 2 ) printf("Sb (Full Integral) = %g\n", Sbval );

                        }

                        MagEphemInfo->Mcurr = LstarInfo3->mInfo->c->M_cd;
                        MagEphemInfo->Mref  = LstarInfo3->mInfo->c->M_cd_McIllwain;
                        MagEphemInfo->Mused = MagEphemInfo->Mref;

                        if (i<MagEphemInfo->nAlpha){
                            MagEphemInfo->I[i]  = Ival;
                            MagEphemInfo->Sb[i] = Sbval;

                            /*
                             *  Compute bounce period (full return time -- not half
                             *  bounce) for a 1Mev Electron.
                             */
                            //MagEphemInfo->Tb[i] = 2.0*MagEphemInfo->Sb[i]/Ek_to_v( 1.0, ELECTRON ); // seconds

                            /*
                             * Compute K invariant in units of G^.5 Re
                             */
                            MagEphemInfo->K[i] = 3.16227766e-3*MagEphemInfo->I[i]*sqrt(MagEphemInfo->Bm[i]);

                            /*
                             *  Add other misc info to MagEphemInfo structure
                             *     Mcurr     -- this is the current value of the dipole moment.
                             *
                             *     Mref      -- this is a reference value of the dipole
                             *                  moment (supposed to be what McIlwain used
                             *                  originally, but who knows now.).

                             *     Mused     -- this is the value of M used in calculations
                             *                  (itll be one of the above).
                             *
                             *      LMcIlwain -- this is the traditional
                             *                   (pitch-angle-dependent) McIlwain 
                             *                   (dimensionless) L-shell parameter.
                             *
                             *      LHilton  -- this is Hilton's approximation to LMcIlwain
                             *
                             */

                            MagEphemInfo->LMcIlwain[i] = LFromIBmM_McIlwain( MagEphemInfo->I[i], MagEphemInfo->Bm[i], MagEphemInfo->Mused );
                            MagEphemInfo->LHilton[i]   = LFromIBmM_Hilton( MagEphemInfo->I[i], MagEphemInfo->Bm[i], MagEphemInfo->Mused );
                        } 

                    }
                }


                /*
                 *  Compute L* 
                 */
                LstarInfo2 = Lgm_CopyLstarInfo( LstarInfo3 );
                if (DoPitchAngle[i] && (LSimple < 10.0) ){

                    LstarInfo2->mInfo->Bm = LstarInfo3->mInfo->Bm;
                    if (LstarInfo3->VerbosityLevel >= 2 ) {
                        printf("\n\n\tComputing L* for: UTC = %g PA = %d  (%g)\n", UTC, i, MagEphemInfo->Alpha[i]);
                        printf("    \t                  I   = %g PA = %d  (%g)\n", MagEphemInfo->I[i], i, MagEphemInfo->Alpha[i]);
                    }
                    LS_Flag = Lstar( &v3, LstarInfo2);
                    if (LstarInfo3->VerbosityLevel >= 2 ) {
                        printf("\tUTC, L*          = %g %g\n", UTC, LstarInfo2->LS );
                        printf("\tUTC, L*_McIlwain = %g %g\n", UTC, LstarInfo2->LS_McIlwain_M);
                        printf("\tUTC, LSimple     = %g %g\n\n\n", UTC, LSimple);
                    }
                    MagEphemInfo->Lstar[i] = LstarInfo2->LS;


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
                    printf(" Lsimple >= 10.0  ( Not doing L* calculation )\n" );
                }
fprintf(fpout, "%.15lf %.15lf\n", MagEphemInfo->Alpha[i], (1.5-LstarInfo2->LS)/LstarInfo->mInfo->Lgm_FindShellLine_I_Tol );
fflush(fpout);

                FreeLstarInfo( LstarInfo2 );

                FreeLstarInfo( LstarInfo3 );
            }

        } 
        // ***** END PARALLEL EXECUTION *****

//        FreeSpline( LstarInfo->mInfo );
    }


    FreeLstarInfo( LstarInfo );

    return;

}
