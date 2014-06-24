#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <Lgm_MagModelInfo.h>
#include <Lgm_LstarInfo.h>
#include <Lgm_MagEphemInfo.h>
#include <Lgm_QinDenton.h>
#include <Lgm_Misc.h>

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
 *                         K:  K value for LCDS at given alpha
 *                 LstarInfo:  LstarInfo structure to specify B-field model, store last valid Lstar, etc.
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

    //Trace to minimum-B
    if ( Lgm_Trace( &Pinner, &v1, &v2, &v3, 120.0, 0.01, TRACE_TOL, LstarInfo_brac1->mInfo ) == LGM_CLOSED ) {
        //Only continue if bracket 1 is closed FL
        //Get L*
        LstarInfo_brac1->mInfo->Bm = LstarInfo_brac1->mInfo->Bmin/sa2;
        LS_Flag = Lstar( &v3, LstarInfo_brac1);
        LCDS = LstarInfo_brac1->LS;
        *K = (LstarInfo_brac1->I0)*sqrt(Lgm_Magnitude(&LstarInfo_brac1->Bmin[0]));
        if (LstarInfo->VerbosityLevel > 0) printf("Found valid inner bracket. Pinner, Pmin_GSM = (%g, %g, %g), (%g, %g, %g)\n", Pinner.x, Pinner.y, Pinner.z, LstarInfo_brac1->mInfo->Pmin.x, LstarInfo_brac1->mInfo->Pmin.y, LstarInfo_brac1->mInfo->Pmin.z);
    } else {
        FreeLstarInfo( LstarInfo_brac1 );
        FreeLstarInfo( LstarInfo_brac2 );
        FreeLstarInfo( LstarInfo_test );
        return(-8); //Inner bracket is bad - bail
    }


    /*
     *  Test outer bracket location.
     */
    Pouter.x = brac2;
    Pouter.y = 0.0; Pouter.z = 0.0;

    LstarInfo_brac2->mInfo->Bfield( &Pouter, &Bvec, LstarInfo_brac2->mInfo );
    Blocal = Lgm_Magnitude( &Bvec );

    //Trace to minimum-B
    if ( Lgm_Trace( &Pouter, &v1, &v2, &v3, 120.0, 0.01, TRACE_TOL, LstarInfo_brac2->mInfo ) == LGM_CLOSED ) {
        //If bracket 2 is closed FL then check for undefined L*. If L* is defined, we have a problem
        //Get L*
        LstarInfo_brac2->mInfo->Bm = LstarInfo_brac2->mInfo->Bmin/sa2;
        LS_Flag = Lstar( &v3, LstarInfo_brac2);
        if (LstarInfo_brac2->LS != LGM_FILL_VALUE) {
            //move outer bracket out and try again
            Pouter.x *= 1.7;
            if ( Lgm_Trace( &Pouter, &v1, &v2, &v3, 120.0, 0.01, TRACE_TOL, LstarInfo_brac2->mInfo ) == LGM_CLOSED ) {
                LS_Flag = Lstar( &v3, LstarInfo_brac2);
                if (LstarInfo_brac2->LS != LGM_FILL_VALUE) {
                    FreeLstarInfo( LstarInfo_brac1 );
                    FreeLstarInfo( LstarInfo_brac2 );
                    FreeLstarInfo( LstarInfo_test );
                    return(-9); //Outer bracket is bad - bail
                }
            }
        }
        if (LstarInfo->VerbosityLevel > 0) printf("Found valid outer bracket. Pouter_GSM, Pmin_GSM = (%g, %g, %g), (%g, %g, %g)\n", Pouter.x, Pouter.y, Pouter.z, LstarInfo_brac2->mInfo->Pmin.x, LstarInfo_brac2->mInfo->Pmin.y, LstarInfo_brac2->mInfo->Pmin.z);
    }

    //if brackets are okay, we've moved on without changing anything except setting initial LCDS value as L* at inner bracket
    nn = 0;
    while (Lgm_Magnitude(&Pouter)-Lgm_Magnitude(&Pinner) > tol){
        if (LstarInfo->VerbosityLevel > 0) printf("Current LCDS iteration, bracket width = %d, %g\n", nn, Lgm_Magnitude(&Pouter)-Lgm_Magnitude(&Pinner));
        if (nn > maxIter) {
            printf("********* EXCEEDED MAXITER\n");
            //free structures before returning
            FreeLstarInfo( LstarInfo_brac1 );
            FreeLstarInfo( LstarInfo_brac2 );
            FreeLstarInfo( LstarInfo_test );
            return(-2); //reached max iterations without achieving tolerance - bail
        }
        Xtest = Pinner.x + (Pouter.x - Pinner.x)/2.0;
        Ptest.x = Xtest; Ptest.y = 0.0; Ptest.z = 0.0;

        LstarInfo_test->mInfo->Bfield( &Ptest, &Bvec, LstarInfo_test->mInfo );
        Blocal = Lgm_Magnitude( &Bvec );
        //Trace to minimum-B
        if ( Lgm_Trace( &Ptest, &v1, &v2, &v3, 120.0, 0.01, TRACE_TOL, LstarInfo_test->mInfo ) == LGM_CLOSED ) {
            //If test point is closed FL then check for undefined L*.
            //Get L*
            LstarInfo_test->mInfo->Bm = LstarInfo_test->mInfo->Bmin/sa2;
            LS_Flag = Lstar( &v3, LstarInfo_test);
            if ( (LS_Flag > 0) || (LstarInfo_test->LS != LGM_FILL_VALUE) ){
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

    double           UTC, brac1, brac2, tol, sJD, eJD, JD, t_cadence;
    double           K[500], LS[500], Alpha[500];
    long int         StartDate, EndDate, Date;
    int              nAlpha, i, Quality, ans, aa, Year, Month, Day, Doy;
    char             Str[128], NewStr[2048];
    char             Filename[1024];
    Lgm_LstarInfo    *LstarInfo = InitLstarInfo(0);
    Lgm_LstarInfo    *LstarInfo3;
    FILE             *fp;
    Lgm_DateTime     DT_UTC;
    Lgm_QinDentonOne qd;

    tol = 0.001;

    nAlpha=18;
    for (i=0; i<nAlpha; i++) {
        Alpha[i] = 5.0+i*5.0;
    }

    // Date and UTC
    StartDate       = 20010101;
    EndDate         = 20010101;
    Lgm_Doy( StartDate, &Year, &Month, &Day, &Doy);
    NewStr[0]       = '\0';
    strcpy(Filename, "%YYYY%MM%DD_LCDS_TS04.txt");
    sprintf( Str, "%4d", Year );   Lgm_ReplaceSubString( NewStr, Filename, "%YYYY", Str );  strcpy( Filename, NewStr );
    sprintf( Str, "%02d", Month ); Lgm_ReplaceSubString( NewStr, Filename, "%MM", Str );   strcpy( Filename, NewStr );
    sprintf( Str, "%02d", Day );   Lgm_ReplaceSubString( NewStr, Filename, "%DD", Str );   strcpy( Filename, NewStr );
    fp = fopen(Filename, "w");

    t_cadence       = 0.5/24.0; //cadence in days
    sJD = Lgm_Date_to_JD( StartDate, 0.0, LstarInfo->mInfo->c);
    eJD = Lgm_Date_to_JD( EndDate, 23.9999, LstarInfo->mInfo->c);

    // Bracket Position in GSM
    brac1 = -3.0;
    brac2 = -13.0;

    Quality = 2;


    /*
     * Write Header
     */
    int nCol = 0;

    fprintf( fp, "# {\n");
    if ( nAlpha > 0 ) {
        fprintf( fp, "#  \"Alpha\":            { \"DESCRIPTION\": \"Pitch Angles.\",\n");
        fprintf( fp, "#                               \"NAME\": \"Pitch Angle\",\n");
        fprintf( fp, "#                              \"TITLE\": \"Pitch Angle\",\n");
        fprintf( fp, "#                              \"LABEL\": \"Pitch Angle\",\n");
        fprintf( fp, "#                          \"DIMENSION\": [ %d ],\n", nAlpha );
        fprintf( fp, "#                             \"VALUES\": [ ");
        for (i=0; i<nAlpha-1; i++) fprintf(fp, "%g, ", Alpha[i] );
        fprintf(fp, "%g ],\n", Alpha[i] ); 

        fprintf( fp, "#                      \"ELEMENT_NAMES\": [ ");
        for (i=0; i<nAlpha-1; i++) fprintf(fp, "\"PA%d\", ", i );
        fprintf(fp, "\"PA%d\" ],\n", i ); 

        fprintf( fp, "#                     \"ELEMENT_LABELS\": [ ");
        for (i=0; i<nAlpha-1; i++) fprintf(fp, "\"%g Deg.\", ", Alpha[i] );
        fprintf(fp, "\"%g Deg.\" ],\n", Alpha[i] ); 
        fprintf( fp, "#                              \"UNITS\": \"Degrees\",\n");
        fprintf( fp, "#                          \"VALID_MIN\":  0.0,\n");
        fprintf( fp, "#                          \"VALID_MAX\": 90.0,\n");
        fprintf( fp, "#                         \"FILL_VALUE\": -1e31 },\n");
        fprintf( fp, "#\n");
    }

    fprintf( fp, "#  \"DateTime\":         { \"DESCRIPTION\": \"The date and time in ISO 8601 compliant format.\",\n");
    fprintf( fp, "#                               \"NAME\": \"IsoDateTime\",\n");
    fprintf( fp, "#                              \"TITLE\": \"ISO DateTime\",\n");
    fprintf( fp, "#                              \"LABEL\": \"Time\",\n");
    fprintf( fp, "#                              \"UNITS\": \"UTC\",\n");
    fprintf( fp, "#                       \"START_COLUMN\": %d },\n", nCol++);
    fprintf( fp, "#\n");

    if ( nAlpha > 0 ) {
        fprintf( fp, "#  \"LCDS\":            { \"DESCRIPTION\": \"Last closed generalized Roederer L-shell value (also known as L*).\",\n");
        fprintf( fp, "#                               \"NAME\": \"LCDS\",\n");
        fprintf( fp, "#                              \"TITLE\": \"LCDS\",\n");
        fprintf( fp, "#                              \"LABEL\": \"LCDS, Dimensionless\",\n");
        fprintf( fp, "#                              \"UNITS\": \"Dimensionless\",\n");
        fprintf( fp, "#                          \"DIMENSION\": [ %d ],\n", nAlpha );
        fprintf( fp, "#                       \"START_COLUMN\": %d,\n", nCol); nCol += nAlpha;
        fprintf( fp, "#                      \"ELEMENT_NAMES\": [ ");
        for (i=0; i<nAlpha-1; i++) fprintf(fp, "\"LCDS_%g\", ", Alpha[i] );
        fprintf(fp, "\"LCDS_%g\" ],\n", Alpha[i] ); 
        fprintf( fp, "#                     \"ELEMENT_LABELS\": [ ");
        for (i=0; i<nAlpha-1; i++) fprintf(fp, "\"LCDS %g!Eo!N\", ", Alpha[i] );
        fprintf(fp, "\"LCDS %g!Eo!N\" ],\n", Alpha[i] ); 
        fprintf( fp, "#                           \"DEPEND_1\": \"Alpha\",\n");
        fprintf( fp, "#                          \"VALID_MIN\": 0.0,\n");
        fprintf( fp, "#                          \"VALID_MAX\": 1000.0,\n");
        fprintf( fp, "#                         \"FILL_VALUE\": -1e31 },\n");
        fprintf( fp, "#\n");
    }
    if ( nAlpha > 0 ) {
        fprintf( fp, "#  \"K\":              { \"DESCRIPTION\": \"Modified second adiabatic invariant, K\",\n");
        fprintf( fp, "#                               \"NAME\": \"K\",\n");
        fprintf( fp, "#                              \"TITLE\": \"K\",\n");
        fprintf( fp, "#                              \"LABEL\": \"K, [R!IE!N G!U1/2!N]\",\n");
        fprintf( fp, "#                              \"UNITS\": \"R!IE!N G!U1/2!N\",\n");
        fprintf( fp, "#                          \"DIMENSION\": [ %d ],\n", nAlpha );
        fprintf( fp, "#                       \"START_COLUMN\": %d,\n", nCol); nCol += nAlpha;
        fprintf( fp, "#                      \"ELEMENT_NAMES\": [ ");
        for (i=0; i<nAlpha-1; i++) fprintf(fp, "\"K_%g\", ", Alpha[i] );
        fprintf(fp, "\"K_%g\" ],\n", Alpha[i] ); 
        fprintf( fp, "#                     \"ELEMENT_LABELS\": [ ");
        for (i=0; i<nAlpha-1; i++) fprintf(fp, "\"K %g!Eo!N\", ", Alpha[i] );
        fprintf(fp, "\"K %g!Eo!N\" ],\n", Alpha[i] ); 
        fprintf( fp, "#                           \"DEPEND_1\": \"Alpha\",\n");
        fprintf( fp, "#                          \"VALID_MIN\": 0.0,\n");
        fprintf( fp, "#                          \"VALID_MAX\": 1000.0,\n");
        fprintf( fp, "#                         \"FILL_VALUE\": -1e31 }\n");
        fprintf( fp, "#\n");
    }
    fprintf( fp, "# } end JSON\n");
    fprintf( fp, "#\n");
    // column header
    fprintf( fp, "# %24s", "Time" );
    for (i=0; i<nAlpha; i++) { sprintf( Str, "L*%d", i ); fprintf(fp, " %8s", Str ); }
    fprintf(fp, "    ");
    for (i=0; i<nAlpha; i++) { sprintf( Str, "K%d", i ); fprintf(fp, " %8s", Str ); }
    fprintf(fp, "    ");
    fprintf(fp, "%s", " \n");

    //loop over date/time at given cadence
    for ( JD = sJD; JD <= eJD; JD += t_cadence ) {
        //set date specific stuff
        Date = Lgm_JD_to_Date( JD, &Year, &Month, &Day, &UTC );
        Lgm_Set_Coord_Transforms( Date, UTC, LstarInfo->mInfo->c);
    
        Lgm_get_QinDenton_at_JD( JD, &qd, 1 );
        Lgm_set_QinDenton( &qd, LstarInfo->mInfo );
    
        //LstarInfo->mInfo->Bfield        = Lgm_B_T89c;
    //LstarInfo->mInfo->Bfield        = Lgm_B_Dungey;
    //LstarInfo->mInfo->InternalModel = LGM_CDIP;
        LstarInfo->mInfo->Bfield        = Lgm_B_TS04;
        LstarInfo->mInfo->Bfield        = Lgm_B_T89c;
        LstarInfo->mInfo->InternalModel = LGM_IGRF;
        LstarInfo->VerbosityLevel       = 1;
        //LstarInfo->mInfo->Kp = 0.3;
    
        /*
         * Compute L*s, Is, Bms, etc...
         */
    
        { // ***** BEGIN PARALLEL EXECUTION *****
    
            /*
             *  Do all of the PAs in parallel. To control how many threads get run
             *  use the enironment variable OMP_NUM_THREADS. For example,
             *          setenv OMP_NUM_THREADS 8
             *  will use 8 threads to do the loop in parallel. Be very careful what gets 
             *  set private here -- the threads must not interfere with each other.
             */
            #pragma omp parallel private(LstarInfo3,ans)
            #pragma omp for schedule(dynamic, 1)
    
            for (aa=0; aa<nAlpha; ++aa) {
                // make a local copy of LstarInfo structure -- needed for multi-threading
                LstarInfo3 = Lgm_CopyLstarInfo( LstarInfo );
                LstarInfo3->PitchAngle = Alpha[aa];
                //printf("Date, UTC, aa, Alpha, tol = %ld, %g, %d, %g, %g\n", Date, UTC, aa, LstarInfo3->PitchAngle, tol);
    
                ans = LCDS( Date, UTC, brac1, brac2, Alpha[aa], tol, Quality, &K[aa], LstarInfo3 );
                if (ans==0) LS[aa] = LstarInfo3->LS;
                if (LstarInfo3->DriftOrbitType == 1) printf("Alpha: %g; Drift Orbit Type: Closed; L* = %g \n", Alpha[aa], LstarInfo3->LS);
                if (LstarInfo3->DriftOrbitType == 2) printf("Drift Orbit Type: Shebansky; L* = %g\n", LstarInfo3->LS);
                if (ans!=0) {
                    K[aa] = LGM_FILL_VALUE;
                    LS[aa] = LGM_FILL_VALUE;
                    printf("**==**==**==** (PA = %g) Return value: %d\n", Alpha[aa], ans);
                }
    
                FreeLstarInfo( LstarInfo3 );
            }
    
        } // ***** END PARALLEL EXECUTION *****
        
        // FIX THIS WRITE OUT
        Lgm_Make_UTC( Date, UTC, &DT_UTC, LstarInfo->mInfo->c );
        Lgm_DateTimeToString( Str, &DT_UTC, 0, 3);
        fprintf(fp, "%24s",     Str );
        for ( i=0; i<nAlpha; ++i ) {
            fprintf( fp, "     %8.6e", LS[i]);
        }
        for ( i=0; i<nAlpha; ++i ) {
            fprintf( fp, "     %8.6e", K[i] );
        }
        fprintf(fp, "%s", " \n");
    
    }
    fclose(fp);
    fflush(fp);

    FreeLstarInfo( LstarInfo );

    return(0);

}
