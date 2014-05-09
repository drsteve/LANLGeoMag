#include <Lgm_MagModelInfo.h>
#include <Lgm_LstarInfo.h>
#include <Lgm_MagEphemInfo.h>
#include <Lgm_QinDenton.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#define MAIN
#define TRACE_TOL   1e-7
#define KP_DEFAULT  1


/*
 *      Input Variables:
 *
 *                      Date: 
 *                       UTC: 
 *                     brac1:  Radial distance for inner edge of search bracket
 *                     brac2:  Radial distance for outer edge of search bracket
 *                     Alpha:  Equatorial pitch angle to compute 
 *                   Quality:  Quality factor (0-8)
 *
 *      Input/OutPut Variables:
 *                      LCDS:  Last closed drift shell Lstar value
 *                         K:  K value for LCDS at given alpha
 *                 LstarInfo:  LstarInfo structure to specify B-field model, etc.
 *
 *              MagEphemInfo:  Structure used to input and output parameters/settings/results to/from routine.
 *  
 */
int LCDS( long int Date, double UTC, double brac1, double brac2, double Alpha, double tol, int Quality, double *K, Lgm_LstarInfo *LstarInfo ) {
    Lgm_LstarInfo   *LstarInfo_brac1, *LstarInfo_brac2, *LstarInfo_test;
    Lgm_Vector      v1, v2, v3, Bvec, Ptest, Pinner, Pouter;
    double          Blocal, Xtest, sa, sa2, LCDS;
    double          nTtoG = 1.0e-5;
    int             LS_Flag, nn, k;
    int             maxIter = 20;

    //Start by creating necessary structures... need an Lstar info for each bracket, plus test point.
    LstarInfo_brac1 = Lgm_CopyLstarInfo( LstarInfo );
    LstarInfo_brac2 = Lgm_CopyLstarInfo( LstarInfo );
    LstarInfo_test = Lgm_CopyLstarInfo( LstarInfo );

    // Set Tolerances
    SetLstarTolerances( Quality, LstarInfo_brac1 );
    SetLstarTolerances( Quality, LstarInfo_brac2 );
    SetLstarTolerances( Quality, LstarInfo_test );

    // set coord transformation 
    Lgm_Set_Coord_Transforms( Date, UTC, LstarInfo_brac1->mInfo->c );
    Lgm_Set_Coord_Transforms( Date, UTC, LstarInfo_brac2->mInfo->c );
    Lgm_Set_Coord_Transforms( Date, UTC, LstarInfo_test->mInfo->c );

    //Check for closed drift shell at brackets
    Pinner.x = brac1;
    Pinner.y = 0.0; Pinner.z = 0.0;
    /*
     *  Test inner bracket location.
     */
    LstarInfo_brac1->mInfo->Bfield( &Pinner, &Bvec, LstarInfo_brac1->mInfo );
    Blocal = Lgm_Magnitude( &Bvec );
    sa = sin( LstarInfo->PitchAngle*RadPerDeg ); sa2 = sa*sa;
    LstarInfo_brac1->mInfo->Bm = Blocal/sa2;

    //Trace to minimum-B
    if ( Lgm_Trace( &Pinner, &v1, &v2, &v3, 120.0, 0.01, TRACE_TOL, LstarInfo_brac1->mInfo ) == LGM_CLOSED ) {
        //Only continue if bracket 1 is closed FL
        //Get L*
        LS_Flag = Lstar( &v3, LstarInfo_brac1);
        LCDS = LstarInfo_brac1->LS;
        *K = (LstarInfo_brac1->I0)*sqrt(Lgm_Magnitude(&LstarInfo_brac1->Bmin[0]));
        if (LstarInfo->VerbosityLevel > 3) printf("Found valid inner bracket. Pinner, Pmin_GSM = (%g, %g, %g), (%g, %g, %g)\n", Pinner.x, Pinner.y, Pinner.z, LstarInfo_brac1->mInfo->Pmin.x, LstarInfo_brac1->mInfo->Pmin.y, LstarInfo_brac1->mInfo->Pmin.z);
    } else {
        FreeLstarInfo( LstarInfo_brac1 );
        FreeLstarInfo( LstarInfo_brac2 );
        FreeLstarInfo( LstarInfo_test );
        return(-1); //Bracket is bad - bail
    }


    /*
     *  Test outer bracket location.
     */
    Pouter.x = brac2;
    Pouter.y = 0.0; Pouter.z = 0.0;

    LstarInfo_brac2->mInfo->Bfield( &Pouter, &Bvec, LstarInfo_brac2->mInfo );
    Blocal = Lgm_Magnitude( &Bvec );
    LstarInfo_brac2->mInfo->Bm = Blocal/sa2;

    //Trace to minimum-B
    if ( Lgm_Trace( &Pouter, &v1, &v2, &v3, 120.0, 0.01, TRACE_TOL, LstarInfo_brac2->mInfo ) == LGM_CLOSED ) {
        //If bracket 2 is closed FL then check for undefined L*. If L* is defined, we have a problem
        //Get L*
        LS_Flag = Lstar( &v3, LstarInfo_brac2);
        if (LstarInfo_brac2->LS != LGM_FILL_VALUE) {
            FreeLstarInfo( LstarInfo_brac1 );
            FreeLstarInfo( LstarInfo_brac2 );
            FreeLstarInfo( LstarInfo_test );
            return(-1); //Bracket is bad - bail
        }
        if (LstarInfo->VerbosityLevel > 3) printf("Found valid outer bracket. Pouter_GSM, Pmin_GSM = (%g, %g, %g), (%g, %g, %g)\n", Pouter.x, Pouter.y, Pouter.z, LstarInfo_brac2->mInfo->Pmin.x, LstarInfo_brac2->mInfo->Pmin.y, LstarInfo_brac2->mInfo->Pmin.z);
    }

    //if brackets are okay, we've moved on without changing anything except setting initial LCDS value as L* at inner bracket
    nn = 0;
    while (Lgm_Magnitude(&Pouter)-Lgm_Magnitude(&Pinner) > tol){
        if (LstarInfo->VerbosityLevel > 1) printf("Current LCDS iteration, bracket width = %d, %g\n", nn, Lgm_Magnitude(&Pouter)-Lgm_Magnitude(&Pinner));
        if (nn > maxIter) {
            printf("********* EXCEEDED MAXITER\n");
            //free structures before returning
            FreeLstarInfo( LstarInfo_brac1 );
            FreeLstarInfo( LstarInfo_brac2 );
            FreeLstarInfo( LstarInfo_test );
            return(-2); //reached max iterations without achieving tolerance - bail
        }
        Xtest = Pinner.x + (Pouter.x - Pinner.x)/2.0;
        Ptest.x = Xtest;

        LstarInfo_test->mInfo->Bfield( &Ptest, &Bvec, LstarInfo_test->mInfo );
        Blocal = Lgm_Magnitude( &Bvec );
        LstarInfo_test->mInfo->Bm = Blocal/sa2;
        //Trace to minimum-B
        if ( Lgm_Trace( &Ptest, &v1, &v2, &v3, 120.0, 0.01, TRACE_TOL, LstarInfo_test->mInfo ) == LGM_CLOSED ) {
            //If test point is closed FL then check for undefined L*.
            //Get L*
            LS_Flag = Lstar( &v3, LstarInfo_test);
            if (LstarInfo_test->LS != LGM_FILL_VALUE) {
                //Drift shell defined
                Pinner = Ptest;
                LCDS = LstarInfo_test->LS;
                *K = (LstarInfo_test->I0)*sqrt(LstarInfo_test->mInfo->Bm*nTtoG);

                //Determine the type of the orbit
                LstarInfo_test->DriftOrbitType = LGM_DRIFT_ORBIT_CLOSED;
                for ( k=0; k<LstarInfo_test->nMinMax; ++k ) {
                    if ( LstarInfo_test->nMinima[k] > 1 ) LstarInfo_test->DriftOrbitType = LGM_DRIFT_ORBIT_CLOSED_SHABANSKY;
                    if ( LstarInfo_test->nMinima[k] < 1 ) {printf("Less than one minimum defined on field line: Exit due to impossibility."); exit(-1); }
                }

                if (LstarInfo->VerbosityLevel > 0) printf("Current LCDS, K is %g, %g\n", LCDS, *K);
                LstarInfo->mInfo->Bm = LstarInfo_test->mInfo->Bm;
                LstarInfo->DriftOrbitType = LstarInfo_test->DriftOrbitType;
            } else {
                Pouter = Ptest;
            }
        } else {
            //FL open -> drift shell open
            Pouter = Ptest;
        }

    nn++;
    }

    LstarInfo->LS = LCDS;
    //free structures
    FreeLstarInfo( LstarInfo_brac1 );
    FreeLstarInfo( LstarInfo_brac2 );
    FreeLstarInfo( LstarInfo_test );
    return(0);
    }




int main( int argc, char *argv[] ){

    double           UTC, Alpha, brac1, brac2, tol, K, JD;
    long int         Date;
    int              nAlpha, i, Quality, ans;
    //char             Filename2[1024];
    Lgm_LstarInfo    *LstarInfo = InitLstarInfo(0);
    FILE             *fp;
    Lgm_QinDentonOne qd;


    // Date and UTC
    Date       = 20120525;
    Date       = 20121008;
    UTC        = 16.0 + 54.0/60.0 + 9.595436/3600.0;
    UTC        = 20.0 + 54.0/60.0 + 9.595436/3600.0;
    JD = Lgm_Date_to_JD( Date, UTC, LstarInfo->mInfo->c);
    Lgm_Set_Coord_Transforms( Date, UTC, LstarInfo->mInfo->c);

    Lgm_get_QinDenton_at_JD( JD, &qd, 1 );
    Lgm_set_QinDenton( &qd, LstarInfo->mInfo );

    // Bracket Position in GSM
    brac1 = -3.0;
    brac2 = -12.0;

Alpha = 85.0;
nAlpha=1;
    Quality = 3;

    LstarInfo->mInfo->Bfield        = Lgm_B_T89c;
//    LstarInfo->mInfo->Bfield        = Lgm_B_TS04;
    LstarInfo->mInfo->InternalModel = LGM_IGRF;
    LstarInfo->VerbosityLevel       = 0;
    LstarInfo->PitchAngle           = Alpha;
    if ( LstarInfo->mInfo->Kp > 5 ) LstarInfo->mInfo->Kp = 5;

    /*
     * Compute L*s, Is, Bms, etc...
     */
    tol = 0.001;
    tol = 0.01;
    ans = LCDS( Date, UTC, brac1, brac2, Alpha, tol, Quality, &K, LstarInfo );
    if (LstarInfo->DriftOrbitType == 1) printf("Drift Orbit Type: Closed; L* = %g \n", LstarInfo->LS);
    if (LstarInfo->DriftOrbitType == 2) printf("Drift Orbit Type: Shebansky; L* = %g\n", LstarInfo->LS);
    if (ans!=0) printf("**==**==**==** Return value: %d\n", ans);
    fp = fopen("LCDS.dat", "w");
    for ( i=0; i<nAlpha; ++i ) {
        fprintf( fp, "%g %g %g\n", Alpha, LstarInfo->LS, K );
    }
    fclose(fp);

    FreeLstarInfo( LstarInfo );

    return(0);

}
