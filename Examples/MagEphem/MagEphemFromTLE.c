#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <Lgm_CTrans.h>
#include <Lgm_Sgp.h>
#include <Lgm_MagEphemInfo.h>
#include <Lgm_DynamicMemory.h>
#include <Lgm_ElapsedTime.h> 
#include <Lgm_Eop.h>                                                                                                                   
#include <ctype.h>


#define KP_DEFAULT 0



/*
 *  Compute position of S/C from TLE (Two Line Elements). 
 *  Then use this to compute Magnetic Ephemerides...
 */

int main( int argc, char *argv[] ){

    long int        StartDate, EndDate;
    double          tsince, JD, StartUT, EndUT;
    double          GpsTime, StartGpsTime, StopGpsTime, GpsInc;
    Lgm_CTrans      *c = Lgm_init_ctrans( 0 );
    Lgm_Vector      Ugsm, Uteme;
    Lgm_DateTime    UTC;
    Lgm_Eop         *e = Lgm_init_eop( 0 );                                                                                            
    Lgm_EopOne      eop;                      
    int             i, nTLEs; 
    char            Line0[100], Line1[100], Line2[100], *ptr;
    char            *InputFile   = "input.txt";
    char            *OutputFile  = "output.txt";
    char            OutputFilename[1024];
    char            IntModel[20], ExtModel[20], opt, ColorizeStr[20];
    int             AppendMode, UseEop, Colorize, Quality;
    FILE            *fp;
    FILE            *fp_MagEphem;


    double           Alpha[1000], Kp, Dst, Delta, FootpointHeight;
    int              nAlpha;
    //Lgm_MagEphemInfo *MagEphemInfo = Lgm_InitMagEphemInfo(0, nPitchAngles);
    Lgm_MagEphemInfo *MagEphemInfo;
    Lgm_ElapsedTimeInfo  tInfo;

    Lgm_ElapsedTimeInit( &tInfo, 255, 95, 0 );


    /* 
     * Open input file and extract:
     *   1) the 3 TLE lines
     *   2) Start Date and Time
     *   3) End Date and Time
     */
    int sY, sM, sD, sh, sm, ss;
    int eY, eM, eD, eh, em, es;
    if ( (fp = fopen( InputFile, "r" )) != NULL ) {
        fgets( Line0, 99, fp );
        fgets( Line1, 99, fp );
        fgets( Line2, 99, fp );
        fscanf( fp, "%*[^:]:%s", OutputFilename );
        fscanf( fp, "%*[^:]:%4d-%2d-%2dT%2d:%2d:%2d", &sY, &sM, &sD, &sh, &sm, &ss );
        fscanf( fp, "%*[^:]:%4d-%2d-%2dT%2d:%2d:%2d", &eY, &eM, &eD, &eh, &em, &es );
        fscanf( fp, "%*[^:]:%lf", &Delta );
        fscanf( fp, "%*[^:]:%s", IntModel );
        fscanf( fp, "%*[^:]:%s", ExtModel );
        fscanf( fp, "%*[^:]:%lf", &Kp );
        fscanf( fp, "%*[^:]:%lf", &Dst );
        fscanf( fp, "%*[^:]:%d", &nAlpha );
        for (i=0; i<nAlpha; i++ ) fscanf( fp, "%lf", &Alpha[i]);
        fscanf( fp, "%*[^:]:%lf", &FootpointHeight );
        fscanf( fp, "%*[^:]:%s", ColorizeStr );
        if ( !strcmp(ColorizeStr, "True") ) {
            Colorize = TRUE;
        } else {
            Colorize = FALSE;
        }
        fscanf( fp, "%*[^:]:%d", &Quality );
        StartDate = sY*10000 + sM*100 + sD;
        StartUT   = sh + sm/60.0 + ss/3600.0;
        EndDate   = eY*10000 + eM*100 + eD;
        EndUT     = eh + em/60.0 + es/3600.0;
    } else {
        printf( "Couldnt open file %s for reading\n", InputFile );
        exit( 1 );
    }
    fclose( fp );

    printf("OutputFilename  = %s\n", OutputFilename);
    printf("StartDate       = %ld\n", StartDate);
    printf("EndDate         = %ld\n", EndDate);
    printf("StartUTC        = %g\n", StartUT);
    printf("EndUTC          = %g\n", EndUT);
    printf("IntModel        = %s\n", IntModel);
    printf("ExtModel        = %s\n", ExtModel);
    printf("Kp              = %g\n", Kp);
    printf("Dst             = %g\n", Dst);
    printf("FootpointHeight = %g\n", FootpointHeight);
    printf("Colorize        = %d\n", Colorize);
    printf("Quality         = %d\n", Quality);
    //exit(0);
    if ( nAlpha > 0 ){
printf("********** nAlpha = %d\n", nAlpha);
        MagEphemInfo = Lgm_InitMagEphemInfo(0, nAlpha);
    } else {
        // doesnt seem to like allocating zero size...
        MagEphemInfo = Lgm_InitMagEphemInfo(0, 1);
    }
i = 0;
int nn = 0;
printf("MagEphemInfo = %p\n", MagEphemInfo);
printf("MagEphemInfo->Shell_Bmin = %p\n", MagEphemInfo->Shell_Bmin);
printf("MagEphemInfo->Shell_Bmin[%d][%d] = %g %g %g\n", i, nn, MagEphemInfo->Shell_Bmin[i][nn].x, MagEphemInfo->Shell_Bmin[i][nn].y, MagEphemInfo->Shell_Bmin[i][nn].z);                                                                                                      

//exit(0);
    /*
     * Remove any extraneous newline and/or linefeeds at the end of the strings.
     * Probably not needed, but TLEs may have non-linux terminating characters...
     */
    if ( (ptr = strstr(Line0, "\n")) != NULL ) *ptr = '\0'; 
    if ( (ptr = strstr(Line0, "\r")) != NULL ) *ptr = '\0';
    if ( (ptr = strstr(Line1, "\n")) != NULL ) *ptr = '\0'; 
    if ( (ptr = strstr(Line1, "\r")) != NULL ) *ptr = '\0';
    if ( (ptr = strstr(Line2, "\n")) != NULL ) *ptr = '\0'; 
    if ( (ptr = strstr(Line2, "\r")) != NULL ) *ptr = '\0';



    // Settings for Lstar calcs
    MagEphemInfo->ComputeVgc = TRUE;
MagEphemInfo->ComputeVgc = FALSE;
    MagEphemInfo->LstarQuality = Quality;
    MagEphemInfo->SaveShellLines = TRUE;
    MagEphemInfo->LstarInfo->LSimpleMax = 10.0;
    MagEphemInfo->LstarInfo->VerbosityLevel = 0;
    MagEphemInfo->LstarInfo->mInfo->VerbosityLevel = 0;
    MagEphemInfo->LstarInfo->mInfo->Lgm_LossConeHeight = FootpointHeight;

//    Kp = 5;
/*
    MagEphemInfo->LstarInfo->mInfo->Bfield        = Lgm_B_edip;
    MagEphemInfo->LstarInfo->mInfo->Bfield        = Lgm_B_igrf;
    MagEphemInfo->LstarInfo->mInfo->Bfield        = Lgm_B_OP77;
    MagEphemInfo->LstarInfo->mInfo->Bfield        = Lgm_B_cdip;
    MagEphemInfo->LstarInfo->mInfo->InternalModel = LGM_CDIP;
*/
    if ( !strcmp( ExtModel, "T87" ) ){
        MagEphemInfo->LstarInfo->mInfo->Bfield = Lgm_B_T87;
    } else if ( !strcmp( ExtModel, "CDIP" ) ){
        MagEphemInfo->LstarInfo->mInfo->Bfield = Lgm_B_cdip;
    } else if ( !strcmp( ExtModel, "EDIP" ) ){
        MagEphemInfo->LstarInfo->mInfo->Bfield = Lgm_B_edip;
    } else if ( !strcmp( ExtModel, "IGRF" ) ){
        MagEphemInfo->LstarInfo->mInfo->Bfield = Lgm_B_igrf;
    } else { //if ( !strcmp( ExtModel, "T89" ) ){
        // default
        MagEphemInfo->LstarInfo->mInfo->Bfield = Lgm_B_T89;
    }

    if ( !strcmp( IntModel, "CDIP" ) ){
        MagEphemInfo->LstarInfo->mInfo->InternalModel = LGM_CDIP;
    } else if ( !strcmp( IntModel, "EDIP" ) ){
        MagEphemInfo->LstarInfo->mInfo->InternalModel = LGM_EDIP;
    } else {
        // default
        MagEphemInfo->LstarInfo->mInfo->InternalModel = LGM_IGRF;
    }

//    MagEphemInfo->LstarInfo->mInfo->Bfield        = Lgm_B_T89;
    MagEphemInfo->LstarInfo->mInfo->Kp = ( Kp >= 0.0 ) ? (int)(Kp+0.5) : KP_DEFAULT;
    if ( MagEphemInfo->LstarInfo->mInfo->Kp > 5 ) MagEphemInfo->LstarInfo->mInfo->Kp = 5;
printf("MagEphemInfo->LstarInfo->mInfo->Kp = %d\n", MagEphemInfo->LstarInfo->mInfo->Kp);

    // Create array of Pitch Angles to compute
    // Read them in from file now...
    //for (nAlpha=0,a=5.0; a<=90.0; a+=5.0,++nAlpha) {
    //    Alpha[nAlpha] = a ;
    //    MagEphemInfo->Alpha[nAlpha] = a;
    //    printf("Alpha[%d] = %g\n", nAlpha, Alpha[nAlpha]);
    //}

    MagEphemInfo->nAlpha = nAlpha;
    for (i=0; i<nAlpha; i++){
        MagEphemInfo->Alpha[i] = Alpha[i];
        printf("Alpha[%d] = %g\n", i, Alpha[i]);
    }




    /*
     * Alloc some memory for the SgpInfo structure and the TLEs array (here we
     * only have a single element in the array)
     */
    _SgpInfo *s   = (_SgpInfo *)calloc( 1, sizeof(_SgpInfo) );
    _SgpTLE *TLEs = (_SgpTLE *)calloc( 1, sizeof(_SgpTLE) );

    
    /*
     * Read in TLEs from the Line0, Line1 and Line2 strings. nTLEs must be
     * initialized to zero by the user.
     */
    nTLEs = 0;
    LgmSgp_ReadTlesFromStrings( Line0, Line1, Line2, &nTLEs, TLEs, 1 );


    /*
     * All the TLEs have their own epoch times in them. And the propagator
     * (sgp4) uses the "time since (in minutes)". So for a given time of
     * interest, we need to compute the tsince needed. Convert Start/End Dates
     * to Julian dates -- they are easier to loop over contiguously.
     */
    Lgm_Make_UTC( StartDate, StartUT, &UTC, c );
    StartGpsTime = Lgm_UTC_to_GpsSeconds( &UTC, c );


    Lgm_Make_UTC( EndDate, EndUT, &UTC, c );
    StopGpsTime = Lgm_UTC_to_GpsSeconds( &UTC, c );

    GpsInc = Delta; // seconds

    if ( (fp = fopen( OutputFile, "w" )) == NULL ) {
        printf( "Couldnt open file %s for writing\n", OutputFile );
        exit( 1 );
    }

    // init SGP4
    LgmSgp_SGP4_Init( s, &TLEs[0] );
    printf("%%%s\n", TLEs[0].Line0 );
    printf("%%%s\n", TLEs[0].Line1 );
    printf("%%%s\n", TLEs[0].Line2 );
    fprintf(fp, "%%%s\n", TLEs[0].Line0 );
    fprintf(fp, "%%%s\n", TLEs[0].Line1 );
    fprintf(fp, "%%%s\n", TLEs[0].Line2 );


    /*
     * Open Mag Ephem file for writing
     */
    AppendMode = FALSE;
    UseEop     = FALSE;
    while ( (opt = getopt( argc, argv, "a::e::" )) != -1 ) {
        switch( opt ) {
           case 'a':
             AppendMode = TRUE;
             break;
           case 'e':
             UseEop = TRUE;
             break;
        }
    }

    if ( AppendMode ){
        fp_MagEphem = fopen( OutputFilename, "ab" );
    } else {
        fp_MagEphem = fopen( OutputFilename, "wb" );
        Lgm_WriteMagEphemHeader( fp_MagEphem, TLEs[0].Line0, TLEs[0].IdNumber, TLEs[0].IntDesig2, IntModel, ExtModel, MagEphemInfo );
    }

    if ( UseEop ) {
        // Read in the EOP vals
        Lgm_read_eop( e );
    }

    // loop over specified time range
    for ( GpsTime = StartGpsTime; GpsTime <= StopGpsTime; GpsTime += GpsInc ) {


        // Convert the current GpsTime back to Date/UT etc..
        // Need JD to compute tsince
        Lgm_GpsSeconds_to_UTC( GpsTime, &UTC, c ) ;
        JD = Lgm_JD( UTC.Year, UTC.Month, UTC.Day, UTC.Time, LGM_TIME_SYS_UTC, c );

                                                                                                                                       
        if ( UseEop ) {
            // Get (interpolate) the EOP vals from the values in the file at the given Julian Date
            Lgm_get_eop_at_JD( JD, &eop, e );

            // Set the EOP vals in the CTrans structure.
            Lgm_set_eop( &eop, c );
        }

        // Set up the trans matrices
        Lgm_Set_Coord_Transforms( UTC.Date, UTC.Time, c );
    
        // "time since" in minutes (thats what SGP4 wants)
        tsince = (JD - TLEs[0].JD)*1440.0; 

        // Call SGP4. Coords are in TEME. 
        LgmSgp_SGP4( tsince, s );
        Uteme.x = s->X/WGS84_A; Uteme.y = s->Y/WGS84_A; Uteme.z = s->Z/WGS84_A;

        // Example of converting TEME->GSM coords.
        Lgm_Convert_Coords( &Uteme, &Ugsm, TEME_TO_GSM, c );


        /*
         * Compute L*s, Is, Bms, Footprints, etc...
         * These quantities are stored in the MagEphemInfo Structure
         */
        printf("\n\n\nDate, ut = %ld %g   Ugsm = %g %g %g \n", UTC.Date, UTC.Time, Ugsm.x, Ugsm.y, Ugsm.z );
        Lgm_ComputeLstarVersusPA( UTC.Date, UTC.Time, &Ugsm, nAlpha, Alpha, MagEphemInfo->LstarQuality, Colorize, MagEphemInfo );

        Lgm_WriteMagEphemData( fp_MagEphem, IntModel, ExtModel, Kp, Dst, MagEphemInfo );

        if ( nAlpha > 0 ){
            WriteMagEphemInfoStruct( "test.dat", nAlpha, MagEphemInfo );
        }

        Lgm_PrintElapsedTime( &tInfo );

    }
    fclose(fp);
    fclose(fp_MagEphem);

    Lgm_free_ctrans( c );
    Lgm_destroy_eop( e );
    free( s );
    free( TLEs );
    Lgm_FreeMagEphemInfo( MagEphemInfo );



    return(0);
}

