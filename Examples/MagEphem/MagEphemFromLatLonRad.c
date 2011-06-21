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
#include <Lgm_Eop.h>                                                                                                                   
#include <Lgm_QinDenton.h>
#include <ctype.h>


#define KP_DEFAULT 0


/*
 *  Compute ephemerii of S/C from input file that contains S/C position in GEO LAt/Lon/Rad
 */

int main( int argc, char *argv[] ){

    Lgm_CTrans       *c = Lgm_init_ctrans( 0 );
    Lgm_Vector       Ugsm, Ugeo;
    Lgm_DateTime     UTC;
    Lgm_Eop          *e = Lgm_init_eop( 0 );                                                                                            
    Lgm_EopOne       eop;                      
    int              i; 
    char             *InputFile   = "input2.txt";
    char             IsoTimeString[1024];
    char             InputFilename[1024];
    char             OutputFilename[1024];
    char             IntModel[20], ExtModel[20], opt, ColorizeStr[20], Line[5000];
    int              AppendMode, UseEop, Colorize;
    FILE             *fp, *fp_in, *fp_MagEphem;
    double           Alpha[1000], Kp, Dst, FootpointHeight, GeoLat, GeoLon, GeoRad;
    int              nAlpha, Quality;
    Lgm_MagEphemInfo *MagEphemInfo;
    Lgm_QinDentonOne p;






    /* 
     * Open input file and extract info.
     */
    if ( (fp = fopen( InputFile, "r" )) != NULL ) {
        fscanf( fp, "%*[^:]:%s", InputFilename );
        fscanf( fp, "%*[^:]:%s", OutputFilename );
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
    } else {
        printf( "Couldnt open file %s for reading\n", InputFile );
        exit( 1 );
    }
    fclose( fp );

    printf("InputFilename   = %s\n", InputFilename);
    printf("OutputFilename  = %s\n", OutputFilename);
    printf("IntModel        = %s\n", IntModel);
    printf("ExtModel        = %s\n", ExtModel);
    printf("Kp              = %g\n", Kp);
    printf("Dst             = %g\n", Dst);
    printf("FootpointHeight = %g\n", FootpointHeight);
    printf("Quality = %d\n", Quality);
    if ( nAlpha > 0 ){
        MagEphemInfo = Lgm_InitMagEphemInfo(0, nAlpha);
    } else {
        // doesnt seem to like allocating zero size...
        MagEphemInfo = Lgm_InitMagEphemInfo(0, 1);
    }


    // Settings for Lstar calcs
    MagEphemInfo->LstarQuality = Quality;
    MagEphemInfo->SaveShellLines = TRUE;
    MagEphemInfo->LstarInfo->LSimpleMax = 10.0;
    MagEphemInfo->LstarInfo->VerbosityLevel = 0;
    MagEphemInfo->LstarInfo->mInfo->VerbosityLevel = 0;
    MagEphemInfo->LstarInfo->mInfo->Lgm_LossConeHeight = FootpointHeight;

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

    
    

    

//    MagEphemInfo->LstarInfo->mInfo->Kp = ( Kp >= 0.0 ) ? (int)(Kp+0.5) : KP_DEFAULT;
//    if ( MagEphemInfo->LstarInfo->mInfo->Kp > 5 ) MagEphemInfo->LstarInfo->mInfo->Kp = 5;
//    printf("MagEphemInfo->LstarInfo->mInfo->Kp = %d\n", MagEphemInfo->LstarInfo->mInfo->Kp);


    /*
     * Set pitch angles
     */
    MagEphemInfo->nAlpha = nAlpha;
    for (i=0; i<nAlpha; i++){
        MagEphemInfo->Alpha[i] = Alpha[i];
        printf("Alpha[%d] = %g\n", i, Alpha[i]);
    }



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
        Lgm_WriteMagEphemHeader( fp_MagEphem, "FIX ME", 99999, "FIX ME", IntModel, ExtModel, MagEphemInfo );
    }

    if ( UseEop ) {
        // Read in the EOP vals
        Lgm_read_eop( e );
    }



    /*
     * Open input file for reading
     */
    if ( (fp_in = fopen( InputFilename, "r" ) ) == NULL ){
        printf("Could not open file %s for reading\n", InputFilename );
        exit(1);
    }



    /*
     *  Loop over the times/positions given. Here we assume the file contains, Time and GEI_J2000 x, y, z
     */
    while ( fgets( Line, 4096, fp_in ) != NULL ) {

        if ( Line[0] != '#' ) {

            sscanf( Line, "%s %lf %lf %lf", IsoTimeString, &GeoLat, &GeoLon, &GeoRad );

            // convert ISO time to DateTime
            IsoTimeStringToDateTime( IsoTimeString, &UTC, c );

            if ( UseEop ) {
                // Get (interpolate) the EOP vals from the values in the file at the given Julian Date
                Lgm_get_eop_at_JD( UTC.JD, &eop, e );

                // Set the EOP vals in the CTrans structure.
                Lgm_set_eop( &eop, c );
            }

            // Set mag model parameters
            Lgm_get_QinDenton_at_JD( UTC.JD, &p, 0 );
            Lgm_set_QinDenton( &p, MagEphemInfo->LstarInfo->mInfo );

            // Set up the trans matrices
            Lgm_Set_Coord_Transforms( UTC.Date, UTC.Time, c );
        

            // Converting GEI->GSM coords.
            Ugeo.x = GeoRad*cos(GeoLon*RadPerDeg)*cos(GeoLat*RadPerDeg);
            Ugeo.y = GeoRad*sin(GeoLon*RadPerDeg)*cos(GeoLat*RadPerDeg);
            Ugeo.z = GeoRad*sin(GeoLat*RadPerDeg);
            Lgm_Convert_Coords( &Ugeo, &Ugsm, GEI2000_TO_GSM, c );


            /*
             * Compute L*s, Is, Bms, Footprints, etc...
             * These quantities are stored in the MagEphemInfo Structure
             */
            printf("\n\n\nDate, ut = %ld %g   Ugsm = %g %g %g \n", UTC.Date, UTC.Time, Ugsm.x, Ugsm.y, Ugsm.z );
            Lgm_ComputeLstarVersusPA( UTC.Date, UTC.Time, &Ugsm, nAlpha, Alpha, MagEphemInfo->LstarQuality, Colorize, MagEphemInfo );

            Lgm_WriteMagEphemData( fp_MagEphem, IntModel, ExtModel, MagEphemInfo->LstarInfo->mInfo->fKp, MagEphemInfo->LstarInfo->mInfo->Dst, MagEphemInfo );

            if ( nAlpha > 0 ){
                WriteMagEphemInfoStruct( "test.dat", nAlpha, MagEphemInfo );
            }

        }


    }
    fclose(fp_in);
    fclose(fp_MagEphem);

    Lgm_free_ctrans( c );
    Lgm_destroy_eop( e );
    Lgm_FreeMagEphemInfo( MagEphemInfo );



    return(0);
}

