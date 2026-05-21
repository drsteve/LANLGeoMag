#include <argp.h>
#include <dirent.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>
#include <Lgm_CTrans.h>
#include <Lgm_DynamicMemory.h>


const  char *ProgramName = "MagEphemFromFile";
const  char *argp_program_version     = "1.0.0";
const  char *argp_program_bug_address = "<mghenderson@lanl.gov>";
static char doc[] = "Runs MagEphemFromFile over range of dates in an openmp loop";




void StringSplit( char *Str, char *StrArray[], int len, int *n );


/*
 *   Description of options accepted. The fields are as follows;
 *
 *   { NAME, KEY, ARG, FLAGS, DOC } where each of these have the following
 *   meaning;
 *      NAME - Name of option's long argument (can be zero).
 *       KEY - Character used as key for the parser and it's short name.
 *       ARG - Name of the option's argument (zero if there isnt one).
 *     FLAGS - OPTION_ARG_OPTIONAL, OPTION_ALIAS, OPTION_HIDDEN, OPTION_DOC,
 *             or OPTION_NO_USAGE
 */
static struct argp_option Options[] = {
    {"StartDateTime",           'S',    "yyyymmdd[Thh:mm:ss]",      0,        "Start date/time in ISO 8601 format. Seconds will be truncated to integers.", 0 },
    {"EndDateTime",             'E',    "yyyymmdd[Thh:mm:ss]",      0,        "End date/time in ISO 8601 format. Seconds will be truncated to integers.", 0 },
    {"Birds",                   'b',    "\"bird1,bird2,etc\"",      0,        "Birds (sats) to use. E.g., \"ns55,ns70\"."   },

    { 0 }
};

//Mandatory arguments
#define     nArgs   0
static char ArgsDoc[] = "MagEphemFromFile -S Date -E Date InFile OutFile";

struct Arguments {
    char        *args[ nArgs ];

    char        StartDateTimeRaw[80];
    char        EndDateTimeRaw[80];

    long int    StartDate;          // Start date (e.g. 20130201)
    long int    EndDate;            // End date (e.g. 20130228)

    long int    StartSeconds;       // Start time in seconds for the first date.
    long int    EndSeconds;         // End time in seconds for the last date.

    char        StartDateTime[80];
    char        EndDateTime[80];

    char        Birds[4096];

};


/* Parse a single option. */
static error_t parse_opt( int key, char *arg, struct argp_state *state ) {

    double       TaiSecs;
    int          ReturnFlag = 0;
    char         TimeString[40];
    Lgm_DateTime d;
    Lgm_CTrans   *c = Lgm_init_ctrans( 0 );


    /* Get the input argument from argp_parse, which we
      know is a pointer to our arguments structure. */
    struct Arguments *arguments = state->input;
    switch( key ) {
        case 'S': // start date
            strncpy( TimeString, arg, 39 );
            strncpy( arguments->StartDateTimeRaw, arg, 39 ); arguments->StartDateTimeRaw[39] = '\0';
            IsoTimeStringToDateTime( TimeString, &d, c );
            arguments->StartDate    = d.Date;
            arguments->StartSeconds = d.Hour*3600 + d.Minute*60 + (int)d.Second;
            Lgm_DateTimeToString( arguments->StartDateTime, &d, 0, 0 );
            break;
        case 'E': // end date
            strncpy( TimeString, arg, 39 );
            strncpy( arguments->EndDateTimeRaw, arg, 39 ); arguments->EndDateTimeRaw[39] = '\0';
            IsoTimeStringToDateTime( TimeString, &d, c );
            arguments->EndDate    = d.Date;

            // if no explicit time field was given assume user wants whole day.
            if ( strstr(TimeString, "T") == NULL ) {
                // be informative about what the exact time range really will be.
                arguments->EndSeconds = d.DaySeconds;
                TaiSecs = d.DaySeconds + Lgm_UTC_to_TaiSeconds( &d, c );
                Lgm_TaiSeconds_to_UTC( TaiSecs, &d, c );
            } else {
                arguments->EndSeconds = d.Hour*3600 + d.Minute*60 + (int)d.Second;
            }
            Lgm_DateTimeToString( arguments->EndDateTime, &d, 0, 0 );
            break; 
        case 'b': // Birds
            strncpy( arguments->Birds, arg, 4095 );
            break;

    }

    return( ReturnFlag );
}
static struct argp argp = { Options, parse_opt, ArgsDoc, doc };




#define EXE_PATH "/home/mghenderson/git/LanlGeoMag_SAFE/Examples/MagEphem/MagEphemFromFile"



/*
 *  Program to run MAGEIS_L3_to_L3Merged as a pool of tasks using openmp.
 *  The typical way to run this code is like this:
 *
 *     ./MagEphemFromFile -F -S 20200105 -E 20200105 -b "ns55" /data2/GPS_DATA/GPS_REFORMAT/%YYYY/%YYYY%MM%DD_%B_Ephem_v1.09.txt /data2/GPS_DATA/GPS_REFORMAT/%YYYY/%YYYY%MM%DD_%B_MagEphem_v1.09.txt
 *
 */
int main( int argc, char *argv[] ){

    struct dirent **namelist;
    int           i, j, n, nn, nFiles, iSat, ntoken, DateFound;
    struct        stat StatBuf;
    char          Filename[256], FilenameFull[512];
    char          Str[512];
    char          *token, **CommandArray;
    int           nCommandArray;
    long int      Date;

    long int      DatesToProcess[5000];
    int           nDatesToProcess;
    struct        Arguments arguments;
    char          Box[20], OutVersion[80];
    char          EphemFile[1024], ReptFile[1024], HopeFile[1024], MagEphem_File[1024], OutFile[1024];

    int           sYear, sMonth, sDay, sDoy;
    int           eYear, eMonth, eDay, eDoy;
    int           Year, Month, Day, Doy;
    double        sJD, eJD, JD, Hour;

    int              nBirds, nBirdsTmp, iBird;
    char             **Birds, **BirdsTmp, Bird[80];

    Lgm_CTrans   *c = Lgm_init_ctrans( 0 );


    /*
     *  Parse CmdLine arguments and options
     */
    argp_parse( &argp, argc, argv, 0, 0, &arguments );

    Lgm_Doy( arguments.StartDate, &sYear, &sMonth, &sDay, &sDoy);
    sJD = Lgm_JD( sYear, sMonth, sDay, 12.0, LGM_TIME_SYS_UTC, c );

    Lgm_Doy( arguments.EndDate, &eYear, &eMonth, &eDay, &eDoy);
    eJD = Lgm_JD( eYear, eMonth, eDay, 12.0, LGM_TIME_SYS_UTC, c );


    nn = 0;
    for ( JD = sJD; JD <= eJD; JD += 1.0 ) {
        DatesToProcess[nn] = Lgm_JD_to_Date( JD, &Year, &Month, &Day, &Hour );
        ++nn;
    }
    nDatesToProcess = nn;
    
    
    LGM_ARRAY_2D( BirdsTmp, 300, 80, char );
    LGM_ARRAY_2D( Birds,    300, 80, char );
    StringSplit( arguments.Birds, BirdsTmp, 80, &nBirdsTmp );


    nBirds = 0;
    for ( iBird = 0; iBird < nBirdsTmp; iBird++ ) {
        strcpy( Bird, BirdsTmp[ iBird ] );
        strcpy( Birds[ nBirds ], BirdsTmp[ iBird ] );
        ++nBirds;
    }





    nCommandArray = nDatesToProcess*2+10;
    strcpy( OutVersion, "v1.0.0" );

    /*
     * Create an array of commands, so we can parcel them out with system()
     */
    LGM_ARRAY_2D( CommandArray, nCommandArray, 2048, char );
    nn = 0;
    for ( i=0; i<nDatesToProcess; i++ ){
        for ( j=0; j<nBirds; ++j ){

            Year = DatesToProcess[i]/10000;

            sprintf( EphemFile,    "/data2/GPS_DATA/GPS_REFORMAT/%4d/%8ld_%s_Ephem_v1.09.txt", Year, DatesToProcess[i], Birds[j] );
            sprintf( OutFile,      "/data2/GPS_DATA/GPS_REFORMAT/%4d/%8ld_%s_MagEphem_v1.09.txt", Year, DatesToProcess[i], Birds[j] );
            sprintf( &CommandArray[nn][0], "%s -e T89D -F -S %8ld -E %8ld -b %s %s %s", EXE_PATH, DatesToProcess[i], DatesToProcess[i], Birds[j], EphemFile, OutFile );
            printf("Command:%s\n\n", CommandArray[nn] );
            ++nn;


        }
    }


    /*********** PROCESS DATES IN PARALLEL **********/
    #pragma omp parallel 
    #pragma omp for schedule(dynamic, 1)
    for ( i=0; i<nn; i++ ){
            system( CommandArray[i] );
            //printf( "%s\n\n\n", CommandArray[i] );
    }
    /*********** PROCESS DATES IN PARALLEL **********/




    LGM_ARRAY_2D_FREE( CommandArray );
    LGM_ARRAY_2D_FREE( BirdsTmp );
    LGM_ARRAY_2D_FREE( Birds );


    Lgm_free_ctrans( c );





    return( 0 );


}
void StringSplit( char *Str, char **StrArray, int len, int *n ) {

    int         nStr;
    const char  delimiters[] = " ,;";
    char        *token, *ss;

    *n   = 0;
    ss   = Str;
    nStr = strlen( Str );
    while ( ( token = strtok( ss, delimiters ) ) != NULL ) {
        strncpy( StrArray[*n], token, len-1 );
        if ( nStr >= len ) StrArray[*n][len-1] = '\0';
        ++(*n);
        ss = NULL;
    }

}
