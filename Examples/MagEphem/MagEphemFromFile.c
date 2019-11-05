#include <stdlib.h>
#include <glob.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <ctype.h>
#include <argp.h>
#include <time.h>
#include <libgen.h>
#include <Lgm/Lgm_CTrans.h>
#include <Lgm/Lgm_Sgp.h>
#include <Lgm/Lgm_MagEphemInfo.h>
#include <Lgm/Lgm_DynamicMemory.h>
#include <Lgm/Lgm_Eop.h>
#include <Lgm/Lgm_QinDenton.h>
#include <Lgm/Lgm_Misc.h>
#include <Lgm/Lgm_HDF5.h>

void StringSplit( char *Str, char *StrArray[], int len, int *n );


#define KP_DEFAULT 0

const  char *ProgramName = "MagEphemFromFile";
const  char *argp_program_version     = "MagEphemFromFile_1.1";
const  char *argp_program_bug_address = "<mghenderson@lanl.gov>";
static char doc[] = "\nComputes magnetic ephemerii of S/C from input file that contains S/C position"
                    " in GEO Lat/Lon/Rad.\n\n InFile and OutFile are paths to files that may contain"
                    " variables that will be substituted if a time range is given.  The"
                    " available time variables are '%YYYY', '%MM', and '%DD' which correspond"
                    " repectively to 4-digit year, 2-digit month (Jan is 01), and 2-digit day of"
                    " month. %B will also get substituted by the list of birds given in the -b option.\n"
                    " Here is an example using time-variables,\n\n \t./MagEphemFromFile -S 20020901 -E 20020930\n"
                    " \t\t/home/jsmith/EphemData/%YYYY/%YYYY%MM%DD_1989-046_Ephem.txt\n"
                    " \t\t/home/jsmith/MagEphemData/%YYYY/%YYYY%MM%DD_1989-046_MagEphem.txt.\n\n Directories"
                    " in the output file will be created if they don't already exist.\n\n";


// Mandatory arguments
#define     nArgs   2
static char ArgsDoc[] = "InFile OutFile";

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
    {"IntModel",        'i',    "internal_model",             0,                                      "Internal Magnetic Field Model to use. Default is IGRF."    },
    {"ExtModel",        'e',    "external_model",             0,                                      "External Magnetic Field Model to use. Default is T89."    },
    {"Birds",           'b',    "\"bird1, bird2, etc\"",      0,                                      "Birds (sats) to use. E.g., \"LANL-02A, 1989-046, POLAR\"."   },
    {"BirdsPath",       'B',    "\"path\"",                   0,                                      "Path for directory containing Birds (sats) subdirectories. Required if you have wildcards in a satellite name."   },
    {"PitchAngles",     'p',    "\"start_pa, end_pa, npa\"",  0,                                      "Pitch angles to compute. Default is \"5.0, 90, 18\"." },
    {"FootPointHeight", 'f',    "height",                     0,                                      "Footpoint height in km. Default is 100km."                  },
    {"Quality",         'q',    "quality",                    0,                                      "Quality to use for L* calculations. Default is 3."      },
    {"StartDate",       'S',    "yyyymmdd",                   0,                                      "StartDate "                              },
    {"EndDate",         'E',    "yyyymmdd",                   0,                                      "EndDate "                                },
    {"Append",          'a',    0,                            0,                                      "Append to OutFile instead of creating a new one"  },
    {"UseEop",          'e',    0,                            0,                                      "Use Earth Orientation Parameters whn comoputing ephemerii" },
    {"Colorize",        'c',    0,                            0,                                      "Colorize output"                         },
    {"Force",           'F',    0,                            0,                                      "Overwrite output file even if it already exists" },
    {"Coords",          'C',    "coord_system",               0,                                      "Coordinate system used in the input file. Can be: LATLONRAD, SM, GSM, GEI2000 or GSE. Default is LATLONRAD." },
    {"verbose",         'v',    "verbosity",                  0,                                      "Produce verbose output"                  },
    {"visualize",       'V',    0,                            0,                                      "Produce input file for ViewDriftShell visualizer"                  },
    {"silent",          's',    0,                            OPTION_ARG_OPTIONAL | OPTION_ALIAS                                                },
    { 0 }
};

struct Arguments {
    char        *args[ nArgs ];       /* START_PA  END_PA  PA_INC */
    int         silent;
    int         verbose;
    int         visualize;

    double      StartPA;
    double      EndPA;
    int         nPA;

    int         Quality;
    int         Colorize;
    int         Force;
    double      FootPointHeight;

    char        IntModel[80];
    char        ExtModel[80];
    char        CoordSystem[80];

    int         Append;
    int         UseEop;

    char        Birds[4096];
    char        BirdsPath[4096];
    long int    StartDate;
    long int    EndDate;
};


/* Parse a single option. */
static error_t parse_opt (int key, char *arg, struct argp_state *state) {

    /* Get the input argument from argp_parse, which we
      know is a pointer to our arguments structure. */
    struct Arguments *arguments = state->input;
    switch( key ) {
        case 'a': // Append
            arguments->Append = 1;
            break;
        case 'b': // Birds
            strncpy( arguments->Birds, arg, 4095 );
            break;
        case 'B': // BirdsPath
            strncpy( arguments->BirdsPath, arg, 4095 );
            break;
        case 'S': // start date
            sscanf( arg, "%ld", &arguments->StartDate );
            break;
        case 'E': // end date
            sscanf( arg, "%ld", &arguments->EndDate );
            break;
        case 'e': // external model
            strcpy( arguments->ExtModel, arg );
            break;
        case 'i': // inbternal model
            strcpy( arguments->IntModel, arg );
            break;
        case 'p':
            sscanf( arg, "%lf, %lf, %d", &arguments->StartPA, &arguments->EndPA, &arguments->nPA );
            break;
        case 'F':
            arguments->Force = 1;
            break;
        case 'f':
            sscanf( arg, "%lf", &arguments->FootPointHeight );
            break;
        case 'q':
            arguments->Quality = atoi( arg );
            break;
        case 'C':
            strcpy( arguments->CoordSystem, arg );
            break;
        case 'c':
            arguments->Colorize = 1;
            break;
        case 's':
            arguments->silent = 1;
            break;
        case 'v':
            arguments->verbose = atoi( arg );
            break;
        case 'V':
            arguments->visualize = 1;
            break;
        case ARGP_KEY_ARG:
            if (state->arg_num >= nArgs) {
                /* Too many arguments. */
                argp_usage (state);
            }
            arguments->args[state->arg_num] = arg;
            break;
        case ARGP_KEY_END:
            if (state->arg_num < nArgs)
            /* Not enough arguments. */
            argp_usage (state);
            break;
        default:
            return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

/* Our argp parser. */
static struct argp argp = { Options, parse_opt, ArgsDoc, doc };





/*
 *  Compute magnetic ephemerii of S/C from input file that contains S/C position in GEO LAt/Lon/Rad
 */
int main( int argc, char *argv[] ){

    struct Arguments arguments;
    glob_t           globbuf;
    int              glob_n;
    Lgm_CTrans       *c = Lgm_init_ctrans( 0 );
    Lgm_Vector       Ugsm, U, Rgeo, W;
    Lgm_DateTime     UTC;
    Lgm_Eop          *e = Lgm_init_eop( 0 );
    Lgm_EopOne       eop;
    double           X, Y, Z;
    int              i;
    char             IsoTimeString[1024];
    char             InputFilename[1024];
    char             OutputFilename[1024];
    char             IntModel[20], ExtModel[20], CoordSystem[80], Line[5000];
    int              AppendMode, UseEop, Colorize, Force;
    FILE             *fp_in, *fp_MagEphem;
    int              nBirds, nBirdsTmp, iBird;
    int              GotBirdsPath;
    char             **Birds, **BirdsTmp, Bird[80];
    double           Inc, Alpha[1000], FootpointHeight, GeoLat, GeoLon, GeoRad;
    int              nAlpha, Quality;
    long int         StartDate, EndDate, Date __attribute__((unused));
    int              sYear, sMonth, sDay, sDoy, eYear, eMonth, eDay, eDoy, Year, Month, Day;
    double           sJD, eJD, JD, Time;
    Lgm_MagEphemInfo *MagEphemInfo;
    Lgm_MagEphemData *med;
    Lgm_QinDentonOne p;

    hid_t           file;
    hid_t           space;
    hid_t           atype;
    hid_t           DataSet, MemSpace;
    herr_t          status __attribute__((unused));
    hsize_t         Dims[4], Offset[4], SlabSize[4];
    int             iT;
    char            **H5_IsoTimes;
    double          UGSM_ARR[3];
    int             H5_nT;
    double          *H5_Alpha;
    double          **H5_Lstar;
    double          **H5_K;

    double pp, rg, cl, T, vel, Beta, Beta2, s, p2c2, E, Ek, Bmin_mag, Bfn_mag, Bfs_mag, Bsc_mag;
    double GeodLat, GeodLong, GeodHeight, MLAT, MLON, MLT, R;
    Lgm_Vector Bsc_gsm, Bvec, Bvec2, badvec;
    double *ProcessedUTC = NULL;
    int iTime;
    int nTimes = 0;
    int AppendFlag, useTS07, staticModel;

    //Set number of elements in MagEphemData to 86401 (N seconds in day +1) as max number of records
    int NNN = 86401;
    med = Lgm_InitMagEphemData( NNN, 80 );

    //Set dummy vars to fill apogee/perigee stuff so it doesn't break in HDF5 versions <1.8.7
    badvec.x = -999;
    badvec.y = -999;
    badvec.z = -999;
    med->H5_nAscend = 1;
    med->H5_nApogee = 1;
    med->H5_nPerigee = 1;
    strcpy(&med->H5_Perigee_IsoTimes[0][0], "1899-01-01T00:00:00Z");
    strcpy(&med->H5_Apogee_IsoTimes[0][0], "1899-01-01T00:00:00Z");
    strcpy(&med->H5_Ascend_IsoTimes[0][0], "1899-01-01T00:00:00Z");

    /*
     * Default option values.
     */
    arguments.StartPA         = 90;     // start at 90.0 Deg.
    arguments.EndPA           = 5.0;    // stop at 2.5 Deg.
    arguments.nPA             = 18;     // 18 pitch angles
    arguments.silent          = 0;
    arguments.verbose         = 0;
    arguments.visualize       = 0;
    arguments.Quality         = 3;
    arguments.Colorize        = 0;
    arguments.Force           = 0;
    arguments.Append          = 0;
    arguments.UseEop          = 0;
    arguments.StartDate       = -1;
    arguments.EndDate         = -1;
    arguments.FootPointHeight = 100.0; // km
    strcpy( arguments.IntModel, "IGRF" );
    strcpy( arguments.ExtModel, "T89c" );
    strcpy( arguments.CoordSystem, "LATLONRAD" );
    arguments.BirdsPath[0]    = '\0';

    /*
     *  Parse CmdLine arguments and options
     */
    argp_parse (&argp, argc, argv, 0, 0, &arguments);
    if (arguments.silent) {
        arguments.verbose = 0; //silent flag overrides verbosity setting
    }

    /*
     * Define pitch angles to use.
     */
    nAlpha = arguments.nPA;
    Inc    = ( nAlpha > 1 ) ? (arguments.EndPA - arguments.StartPA)/(double)(nAlpha-1) : 0.0;
    for (i=0; i<nAlpha; i++ ) Alpha[i] = arguments.StartPA + i*Inc;

    /*
     *  Set input and output filenames
     */
    strcpy( InputFilename,  arguments.args[0] );
    strcpy( OutputFilename, arguments.args[1] );

    /*
     *  Set other options
     */
    FootpointHeight = arguments.FootPointHeight;
    Quality         = arguments.Quality;
    Colorize        = arguments.Colorize;
    Force           = arguments.Force;
    AppendMode      = arguments.Append;
    UseEop          = arguments.UseEop;
    StartDate       = arguments.StartDate;
    EndDate         = arguments.EndDate;
    strcpy( IntModel,  arguments.IntModel );
    strcpy( ExtModel,  arguments.ExtModel );
    strcpy( CoordSystem,  arguments.CoordSystem );

    LGM_ARRAY_2D( BirdsTmp, 300, 80, char );
    LGM_ARRAY_2D( Birds,    300, 80, char );
    StringSplit( arguments.Birds, BirdsTmp, 80, &nBirdsTmp );


    char *Path       = (char *)calloc( 2056, sizeof( char ) );
    char *BirdsPath  = (char *)calloc( 2056, sizeof( char ) );
    strncpy( BirdsPath,  arguments.BirdsPath, 2055 );
    GotBirdsPath     = ( arguments.BirdsPath[0] != '\0') ? TRUE : FALSE;
    

    /*
     *   Do globbing on birds. I.e. if there are wildcards in any of the birds,
     *   then expand them out. To do this we would have to know the directory
     *   in which to look for the S/C subdirs.
     */
    strcpy( BirdsPath, "/home/mgh/git/dream/Spacecraft" );
    nBirds = 0;
    for ( iBird = 0; iBird < nBirdsTmp; iBird++ ) {

        strcpy( Bird, BirdsTmp[ iBird ] );

        if ( strpbrk( Bird, "*?[]!" ) != NULL ){

            if ( !GotBirdsPath ) {
                printf("Wildcards given in Birds, but no path given (use -B option to specify where the satellite dirs are located.)\n");
                exit(0);
            }

            /*
             * User has provided wildcards in this bird name.
             * Construct path with wild cards and expand it.
             */
            sprintf( Path, "%s/%s", BirdsPath, Bird );
            glob( Path, GLOB_NOCHECK, NULL, &globbuf );

            for ( glob_n=0; glob_n<globbuf.gl_pathc; glob_n++){
                strcpy( Birds[ nBirds ], basename( globbuf.gl_pathv[glob_n] ) ); 
                ++nBirds;
            }
            
        } else {

            // no wildcards -- just copy over
            strcpy( Birds[ nBirds ], BirdsTmp[ iBird ] ); 
            ++nBirds;

        }
        
    }

    if ( nAlpha > 0 ){
        MagEphemInfo = Lgm_InitMagEphemInfo(0, nAlpha);
    } else {
        // doesnt seem to like allocating zero size...
        MagEphemInfo = Lgm_InitMagEphemInfo(0, 1);
    }


    // Settings for Lstar calcs
    Lgm_SetLstarTolerances(Quality, 24, MagEphemInfo->LstarInfo);
    MagEphemInfo->SaveShellLines = TRUE;
    MagEphemInfo->LstarInfo->LSimpleMax = 10.0;
    MagEphemInfo->LstarInfo->VerbosityLevel = arguments.verbose;
    MagEphemInfo->LstarInfo->mInfo->VerbosityLevel = 0;
    MagEphemInfo->LstarInfo->mInfo->Lgm_LossConeHeight = FootpointHeight;

    staticModel = FALSE;
    useTS07 = FALSE;
    if ( !strcmp( ExtModel, "T87" ) ){
        MagEphemInfo->LstarInfo->mInfo->Bfield = Lgm_B_T87;
    } else if ( !strcmp( ExtModel, "CDIP" ) ){
        MagEphemInfo->LstarInfo->mInfo->Bfield = Lgm_B_cdip;
        staticModel = TRUE;
    } else if ( !strcmp( ExtModel, "EDIP" ) ){
        MagEphemInfo->LstarInfo->mInfo->Bfield = Lgm_B_edip;
        staticModel = TRUE;
    } else if ( !strcmp( ExtModel, "DUNGEY" ) ){
        MagEphemInfo->LstarInfo->mInfo->Bfield = Lgm_B_Dungey;
        MagEphemInfo->LstarInfo->mInfo->InternalModel = LGM_DUNGEY;
        strcpy(IntModel, "DUNGEY");
        staticModel = TRUE;
    } else if ( !strcmp( ExtModel, "IGRF" ) ){
        MagEphemInfo->LstarInfo->mInfo->Bfield = Lgm_B_igrf;
        staticModel = TRUE;
    } else if ( !strcmp( ExtModel, "TS04" ) ){
        MagEphemInfo->LstarInfo->mInfo->Bfield = Lgm_B_TS04;
    } else if ( !strcmp( ExtModel, "TS04D" ) ){
        MagEphemInfo->LstarInfo->mInfo->Bfield = Lgm_B_TS04;
    } else if ( !strcmp( ExtModel, "T89" ) ){
        MagEphemInfo->LstarInfo->mInfo->Bfield = Lgm_B_T89;
    } else if ( !strcmp( ExtModel, "T89c" ) ){
        MagEphemInfo->LstarInfo->mInfo->Bfield = Lgm_B_T89c;
    } else if ( !strcmp( ExtModel, "T96" ) ){
        MagEphemInfo->LstarInfo->mInfo->Bfield = Lgm_B_T96;
    } else if ( !strcmp( ExtModel, "T01S" ) ){
        MagEphemInfo->LstarInfo->mInfo->Bfield = Lgm_B_T01S;
    } else if ( !strcmp( ExtModel, "T02" ) ){
        MagEphemInfo->LstarInfo->mInfo->Bfield = Lgm_B_T02;
    } else if ( !strcmp( ExtModel, "OP77" ) ){
        MagEphemInfo->LstarInfo->mInfo->Bfield = Lgm_B_OP77;
        staticModel = TRUE;
    } else if ( !strcmp( ExtModel, "TS07" ) ){
        MagEphemInfo->LstarInfo->mInfo->Bfield = Lgm_B_TS07;
        useTS07 = TRUE;
    } else { //if ( !strcmp( ExtModel, "T89" ) ){
        // default
        strcpy(ExtModel, "T89c");
        MagEphemInfo->LstarInfo->mInfo->Bfield = Lgm_B_T89c;
    }

    if ( !strcmp( IntModel, "CDIP" ) ){
        MagEphemInfo->LstarInfo->mInfo->InternalModel = LGM_CDIP;
    } else if ( !strcmp( IntModel, "EDIP" ) ){
        MagEphemInfo->LstarInfo->mInfo->InternalModel = LGM_EDIP;
    } else {
        // default
        MagEphemInfo->LstarInfo->mInfo->InternalModel = LGM_IGRF;
    }

    /*
     *  Print summary of our options
     */
    if ( !arguments.silent ) {
        printf("\n\n");
        printf( "\t        Program/Version: %s\n", argp_program_version );
        printf( "\t         Bug Reports To: %s\n", argp_program_bug_address );
        printf( "\t             Input File: %s\n", InputFilename );
        printf( "\t            Output File: %s\n", OutputFilename );
        printf( "\t     Input Coord System: %s\n", CoordSystem );
        printf( "\t Number of Pitch Angles: %d\n", nAlpha );
        printf( "\t Pitch Angles [Degrees]:" );
        for (i=0; i<nAlpha; i++ ) printf( " %g", Alpha[i] );
        printf("\n");
        printf( "\t         Internal Model: %s\n", IntModel );
        printf( "\t         External Model: %s\n", ExtModel );
        printf( "\t   FootpointHeight [km]: %g\n", FootpointHeight );
        printf( "\t             L* Quality: %d\n", Quality );
        printf( "\t           Force output: %s\n", Force ? "yes" : "no" );
        printf( "\t   Append to OutputFile: %s\n", AppendMode ? "yes" : "no" );
        printf( "\t                Use Eop: %s\n", UseEop ? "yes" : "no" );
        printf( "\t Colorize Thread Output: %d\n", Colorize );
        printf( "\t                Verbose: %d\n", arguments.verbose); // ? "yes" : "no" );
        printf( "\t                 Silent: %s\n", arguments.silent  ? "yes" : "no" );
        if ( (StartDate > 0)&&(EndDate > 0) ){
            printf( "\t              StartDate: %ld\n", StartDate );
            printf( "\t                EndDate: %ld\n", EndDate );
        }

        if ( nBirds > 1  ){
            printf( "\t                  Birds:" );
            for (iBird=0; iBird<nBirds-1; iBird++) printf(" %s,", Birds[iBird]);
            printf(" %s\n", Birds[iBird]);
        } else if ( nBirds == 1  ){
            printf( "\t                   Bird:" );
            printf(" %s\n", Birds[0]);
        } else {
            printf("At least one bird name must be supplied with -b option.\n");
            exit(0);
        }
    }


    int SubstituteVars = TRUE;
    if ( (StartDate > 0)&&(EndDate > 0) ){
        printf("StartDate = %ld\n", StartDate);
        Lgm_Doy( StartDate, &sYear, &sMonth, &sDay, &sDoy);
        sJD = Lgm_JD( sYear, sMonth, sDay, 12.0, LGM_TIME_SYS_UTC, c );

        Lgm_Doy( EndDate, &eYear, &eMonth, &eDay, &eDoy);
        eJD = Lgm_JD( eYear, eMonth, eDay, 12.0, LGM_TIME_SYS_UTC, c );

        SubstituteVars = TRUE;
    } else {
        // only do a single file
        sJD = eJD = 0;
        SubstituteVars = FALSE;
    }


    char *BaseDir, NewStr[2048], Str[24], Command[4096];
    char *OutFile    = (char *)calloc( 2056, sizeof( char ) );
    char *HdfOutFile = (char *)calloc( 2056, sizeof( char ) );
    char *InFileTmp  = (char *)calloc( 2056, sizeof( char ) );
    char *InFile     = (char *)calloc( 2056, sizeof( char ) );




    /*
     * Set pitch angles in MagEphemInfo structure
     */
    MagEphemInfo->nAlpha = nAlpha;
    med->H5_nAlpha = nAlpha;
    for (i=0; i<nAlpha; i++) {
        MagEphemInfo->Alpha[i]    = Alpha[i];
        med->H5_Alpha[i] = Alpha[i];
    }






    

    for ( JD = sJD; JD <= eJD; JD += 1.0 ) {

        Date = Lgm_JD_to_Date( JD, &Year, &Month, &Day, &Time );



        /*
         * loop over all birds
         */
         for ( iBird = 0; iBird < nBirds; iBird++ ) {

            strcpy( InFile, InputFilename );
            strcpy( OutFile, OutputFilename );
            strcpy( Bird, Birds[ iBird ] );
            strcpy( InFile, InputFilename );

            if ( SubstituteVars ) {

                // Substitute times in the files.
                NewStr[0] = '\0';
                sprintf( Str, "%4d", Year );   Lgm_ReplaceSubString( NewStr, InFile, "%YYYY", Str );  strcpy( InFile, NewStr );
                sprintf( Str, "%02d", Month ); Lgm_ReplaceSubString( NewStr, InFile, "%MM", Str );   strcpy( InFile, NewStr );
                sprintf( Str, "%02d", Day );   Lgm_ReplaceSubString( NewStr, InFile, "%DD", Str );   strcpy( InFile, NewStr );
                Lgm_ReplaceSubString( NewStr, InFile, "%B", Bird );    strcpy( InFile, NewStr );

                sprintf( Str, "%4d", Year );   Lgm_ReplaceSubString( NewStr, OutFile, "%YYYY", Str );  strcpy( OutFile, NewStr );
                sprintf( Str, "%02d", Month ); Lgm_ReplaceSubString( NewStr, OutFile, "%MM", Str );   strcpy( OutFile, NewStr );
                sprintf( Str, "%02d", Day );   Lgm_ReplaceSubString( NewStr, OutFile, "%DD", Str );   strcpy( OutFile, NewStr );
                Lgm_ReplaceSubString( NewStr, OutFile, "%B", Bird );    strcpy( OutFile, NewStr );

            }

            sprintf( HdfOutFile, "%s.h5", OutFile );
            strcat( OutFile, ".txt" );

            if ( SubstituteVars ) {
                printf( "\n\n      Processing Date: %4d-%02d-%02d    Bird: %s\n", Year, Month, Day, Bird );
            } else {
                printf( "\n\n      Processing File: %s\n", InFile );
            }
            printf( "      -------------------------------------------------------------------------------------------\n");
            printf( "           Input File: %s\n", InFile);
            printf( "          Output File: %s\n", OutFile);
            printf( "     HDF5 Output File: %s\n\n", HdfOutFile);


            // Create Base directory if it hasnt been created yet.
            char *dirc = strdup( OutFile );
            BaseDir    = dirname(dirc);
            sprintf( Command, "mkdir -p %s", BaseDir); system( Command );
            free( dirc );


            /*
             *   Check to see if HdfOutFile exists or not.
             */
            int     StatError, FileExists;
            struct stat StatBuf;
            StatError = stat( HdfOutFile, &StatBuf );

            FileExists = FALSE;
            if ( StatError != -1 ) {

                FileExists = TRUE;
                if ( !Force && !AppendMode) {
                    printf("\n\n\tHdfOutfile already exists (use -F option to force processing): %s\n", HdfOutFile );
                } else if (!AppendMode) {
                    printf("\n\n\tWarning. Existing HdfOutfile will be overwritten: %s \n", HdfOutFile );
                } else {
                  printf("\n\n\tExisting HdfOutfile will be appended: %s \n", HdfOutFile );
                }

                if ( StatBuf.st_size < 1000LL ){
                    printf("\t\t             File size:  %lld B\n", (long long)StatBuf.st_size );
                } else if ( StatBuf.st_size < 1000000LL ){
                    printf("\t\t             File size:  %g kB\n", (double)StatBuf.st_size/1.0e3 );
                } else if ( StatBuf.st_size < 1000000000LL ){
                    printf("\t\t             File size:  %g MB\n", (double)StatBuf.st_size/1.0e6 );
                } else {
                    printf("\t\t             File size:  %g GB\n", (double)StatBuf.st_size/1.0e9 );
                }
                printf("\t\t    Last status change:  %s", ctime(&StatBuf.st_ctime));
                printf("\t\t      Last file access:  %s", ctime(&StatBuf.st_atime));
                printf("\t\tLast file modification:  %s", ctime(&StatBuf.st_mtime));
                printf("\n");

                /*
                 *   If HdfOutFile exists, check to see if it is readable.
                 */
                if ( !H5Fis_hdf5( HdfOutFile ) ) {
                    printf("\t  Outfile Not Readable: %s. Forcing regeneration of file.\n\n\n", HdfOutFile );
                    FileExists = FALSE;
                }

            }

            hid_t   file_id=0;
            double *xarray, *yarray, *zarray;
            char **UTCarray, isHDF5=0;
            int idims=0;

            if ( !FileExists || Force || AppendMode ) {
                // Check if the file is HDF5 
                char *ext;
                if (((ext = strrchr(InFile, '.')) != NULL) && (strcmp(ext, ".h5") == 0)) {
                    isHDF5 = 1;
                    file_id = H5Fopen(InFile, H5F_ACC_RDONLY, H5P_DEFAULT);
                    if (file_id < 0) //failed to open
                        continue;
                    xarray = Get_DoubleDataset_1D(file_id, "Latitude", Dims);
                    yarray = Get_DoubleDataset_1D(file_id, "Longitude", Dims);
                    zarray = Get_DoubleDataset_1D(file_id, "Radius", Dims);
                    UTCarray = Get_StringDataset_1D(file_id, "UTC", Dims);
                }
                else if ( (fp_in = fopen( InFile, "r" ) ) == NULL ){
                    printf("\tCould not open file %s for reading\n", InFile );
                    continue;
                } //else { 

                printf( "\t    Reading from file: %s\n", InFile );
                 // Open Mag Ephem file for writing/appending
                if ( !Force && AppendMode && FileExists ) {
                  fp_MagEphem = fopen( OutFile, "a" );
                  printf("\t    Appending to file: %s\n", OutFile );
                } else {
                  fp_MagEphem = fopen( OutFile, "w" );
                  Lgm_WriteMagEphemHeader( fp_MagEphem, argp_program_version,ExtModel, 999, "FIX ME", 99999, "FIX ME", NULL, 0, NULL, NULL, 0, NULL, NULL,0, NULL, NULL, MagEphemInfo );
                  printf("\t      Writing to file: %s\n", OutFile );
                }
                if ( UseEop ) {
                    Lgm_read_eop( e );// Read in the EOP vals
                }
                 // Open MagEphem hdf5 file for writing/appending and Create vars
                if ( !Force && AppendMode && FileExists ) {
                  /* Get the time variable from the HDF5 file to know what has already been processed */
                  file = H5Fopen( HdfOutFile, H5F_ACC_RDWR, H5P_DEFAULT );
                  ProcessedUTC = Get_DoubleDataset_1D( file, "UTC", Dims );
                  nTimes = Dims[0];
                } else {
                    file    = H5Fcreate( HdfOutFile, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteMagEphemHeaderHdf( file, argp_program_version, ExtModel,  999,   "FIX ME",    99999, "FIX ME", "FIX ME",       1,       NULL,      NULL,        1,          NULL,        NULL,       1,           NULL,         NULL, MagEphemInfo, med );
                }

                /* Loop over the times/positions given. Here we assume the file contains, Time and 3 columns: X, Y, Z.
                 *  The meaning of the X, Y, Z columns are dependent upon what CoordSystem is set to;
                 *        CoordSystem          X                    Y               Z
                 *        ------------------------------------------------------------------
                 *          LATLONRAD       Latitude (Deg.)   Longitude (Deg.)  Radius (Re)
                 *             SM              Xsm                 Ysm             Zsm
                 */
                med->H5_nT = 0;
                while (1) {
                    if (isHDF5) {
                        if (idims >= Dims[0])
                            break;
                        X = xarray[idims];
                        Y = yarray[idims];
                        Z = zarray[idims];
                        snprintf(IsoTimeString, sizeof(IsoTimeString), "%s", UTCarray[idims]);
                        idims++;
                    }

                    else {
                        if (fgets( Line, 4096, fp_in ) == NULL)
                            break;
                        if (Line[0] == '#')
                            continue;
                        sscanf(Line, "%s %lf %lf %lf", IsoTimeString, &X, &Y, &Z );
                    }

                    // convert ISO time to DateTime
                    IsoTimeStringToDateTime( IsoTimeString, &UTC, c );

                    // Determine whether time has been processed
                    if ( !Force && AppendMode && FileExists ) {// Assume not processed
                      AppendFlag = 1;
                      for(iTime=0; iTime<nTimes; iTime++) {
                        if ( UTC.Time == ProcessedUTC[iTime]) {
                          // Match found - Do not append
                          AppendFlag = 0;
                        }
                      }
                      // Skip to next timestep
                      if ( !AppendFlag ) {
                        // Increment the data row so that the HDF5 appends the data correctly
                        med->H5_nT++;
                        continue;
                      }
                    }

                    if ( UseEop ) {
                        // Get (interpolate) the EOP vals from the values in the file at the given Julian Date
                        Lgm_get_eop_at_JD( UTC.JD, &eop, e );
                        // Set the EOP vals in the CTrans structure.
                        Lgm_set_eop( &eop, c );
                    }

                    // Set mag model parameters
                    if ((!useTS07) && (!staticModel)) {
                        Lgm_get_QinDenton_at_JD( UTC.JD, &p, 1, 1 );
                        Lgm_set_QinDenton( &p, MagEphemInfo->LstarInfo->mInfo );
                    } 
                    
                    else if (useTS07){ //only TS07 should get here...
                        Lgm_Init_TS07( &(MagEphemInfo->LstarInfo->mInfo->TS07_Info) );
                        Lgm_SetCoeffs_TS07( UTC.Date, UTC.Time, &(MagEphemInfo->LstarInfo->mInfo->TS07_Info) );
                    }
                    
                    //static models don't require QD or other coeffs so should fall through this part without reading files...
                    // Set up the trans matrices
                    Lgm_Set_Coord_Transforms( UTC.Date, UTC.Time, c );

                    if ( !strcmp( CoordSystem, "LATLONRAD" ) ){
                        // Convert GEO->GSM coords.
                        GeoLat = X*RadPerDeg; GeoLon = Y*RadPerDeg; GeoRad = Z;
                        U.x = GeoRad*cos(GeoLon)*cos(GeoLat);
                        U.y = GeoRad*sin(GeoLon)*cos(GeoLat);
                        U.z = GeoRad*sin(GeoLat);
                        Lgm_Convert_Coords( &U, &Ugsm, GEO_TO_GSM, c );
                    } else if ( !strcmp( CoordSystem, "SM" ) ){
                        // Convert SM->GSM coords.
                        U.x = X; U.y = Y; U.z = Z;
                        Lgm_Convert_Coords( &U, &Ugsm, SM_TO_GSM, c );
                    } else if ( !strcmp( CoordSystem, "GSM" ) ){
                        // No conversion needed.
                        Ugsm.x = X; Ugsm.y = Y; Ugsm.z = Z;
                    } else if ( !strcmp( CoordSystem, "GEI2000" ) ){
                        // Convert GEI2000->GSM coords.
                        U.x = X; U.y = Y; U.z = Z;
                        Lgm_Convert_Coords( &U, &Ugsm, GEI2000_TO_GSM, c );
                    } else if ( !strcmp( CoordSystem, "GSE" ) ){
                        // Convert GEI->GSM coords.
                        U.x = X; U.y = Y; U.z = Z;
                        Lgm_Convert_Coords( &U, &Ugsm, GSE_TO_GSM, c );
                    }

                    /*
                     * Compute L*s, Is, Bms, Footprints, etc...
                     * These quantities are stored in the MagEphemInfo Structure
                     */
                    printf("\n\n\t[ %s ]: %s  Bird: %s Ugsm: %g %g %g Re\n", ProgramName, IsoTimeString, Bird, Ugsm.x, Ugsm.y, Ugsm.z );
                    printf("\t--------------------------------------------------------------------------------------------------\n");
                    Lgm_ComputeLstarVersusPA( UTC.Date, UTC.Time, &Ugsm, nAlpha, Alpha, Colorize, MagEphemInfo );
                    Lgm_WriteMagEphemData( fp_MagEphem, IntModel, ExtModel, MagEphemInfo->LstarInfo->mInfo->fKp, MagEphemInfo->LstarInfo->mInfo->Dst, MagEphemInfo );

                    if (( nAlpha > 0 ) && (arguments.visualize)) {
                        WriteMagEphemInfoStruct( "test.dat", nAlpha, MagEphemInfo );
                    }

                    // Fill arrays for dumping out as HDF5 files
                    strcpy( med->H5_IsoTimes[ med->H5_nT ], IsoTimeString );
                    strcpy( med->H5_IntModel[ med->H5_nT ], IntModel );
                    strcpy( med->H5_ExtModel[ med->H5_nT ], ExtModel );
                    switch ( MagEphemInfo->FieldLineType ) {
                        case LGM_OPEN_IMF:
                                            sprintf( med->H5_FieldLineType[ med->H5_nT ], "%s",  "LGM_OPEN_IMF" ); // FL Type
                                            break;
                        case LGM_CLOSED:
                                            sprintf( med->H5_FieldLineType[ med->H5_nT ], "%s",  "LGM_CLOSED" ); // FL Type
                                            break;
                        case LGM_OPEN_N_LOBE:
                                            sprintf( med->H5_FieldLineType[ med->H5_nT ], "%s",  "LGM_OPEN_N_LOBE" ); // FL Type
                                            break;
                        case LGM_OPEN_S_LOBE:
                                            sprintf( med->H5_FieldLineType[ med->H5_nT ], "%s",  "LGM_OPEN_S_LOBE" ); // FL Type
                                            break;
                        case LGM_INSIDE_EARTH:
                                            sprintf( med->H5_FieldLineType[ med->H5_nT ], "%s",  "LGM_INSIDE_EARTH" ); // FL Type
                                            break;
                        case LGM_TARGET_HEIGHT_UNREACHABLE:
                                            sprintf( med->H5_FieldLineType[ med->H5_nT ], "%s",  "LGM_TARGET_HEIGHT_UNREACHABLE" ); // FL Type
                                            break;
                        default:
                                            sprintf( med->H5_FieldLineType[ med->H5_nT ], "%s",  "UNKNOWN FIELD TYPE" ); // FL Type
                                            break;
                    }
                    med->H5_Date[ med->H5_nT ]           = UTC.Date;
                    med->H5_Doy[ med->H5_nT ]            = UTC.Doy;
                    med->H5_UTC[ med->H5_nT ]            = UTC.Time;
                    med->H5_JD[ med->H5_nT ]             = UTC.JD;
                    med->H5_InOut[ med->H5_nT ]          = MagEphemInfo->InOut;
                    med->H5_OrbitNumber[ med->H5_nT ]    = MagEphemInfo->OrbitNumber;
                    med->H5_GpsTime[ med->H5_nT ]        = Lgm_UTC_to_GpsSeconds( &UTC, c );
                    med->H5_TiltAngle[ med->H5_nT ]      = c->psi*DegPerRad;

                    med->H5_Rgsm[ med->H5_nT ][0]        = Ugsm.x;
                    med->H5_Rgsm[ med->H5_nT ][1]        = Ugsm.y;
                    med->H5_Rgsm[ med->H5_nT ][2]        = Ugsm.z;

                    Lgm_Set_Coord_Transforms( UTC.Date, UTC.Time, c );
                    Lgm_Convert_Coords( &Ugsm, &Rgeo, GSM_TO_GEO, c );      Lgm_VecToArr( &Rgeo, &med->H5_Rgeo[ med->H5_nT ][0] );
                    Lgm_Convert_Coords( &Ugsm, &W,    GSM_TO_SM, c );       Lgm_VecToArr( &W,    &med->H5_Rsm[ med->H5_nT ][0] );
                    Lgm_Convert_Coords( &Ugsm, &W,    GSM_TO_GEI2000, c );  Lgm_VecToArr( &W,    &med->H5_Rgei[ med->H5_nT ][0] );
                    Lgm_Convert_Coords( &Ugsm, &W,    GSM_TO_GSE, c );      Lgm_VecToArr( &W,    &med->H5_Rgse[ med->H5_nT ][0] );

                    Lgm_WGS84_to_GEOD( &Rgeo, &GeodLat, &GeodLong, &GeodHeight );
                    Lgm_SetArrElements3( &med->H5_Rgeod[ med->H5_nT ][0],        GeodLat, GeodLong, GeodHeight );
                    Lgm_SetArrElements2( &med->H5_Rgeod_LatLon[ med->H5_nT ][0], GeodLat, GeodLong );
                    med->H5_Rgeod_Height[ med->H5_nT ] = GeodHeight;

                    Lgm_Convert_Coords( &Ugsm, &W, GSM_TO_CDMAG, c );
                    Lgm_CDMAG_to_R_MLAT_MLON_MLT( &W, &R, &MLAT, &MLON, &MLT, c );
                    med->H5_CDMAG_MLAT[ med->H5_nT ] = MLAT;
                    med->H5_CDMAG_MLON[ med->H5_nT ] = MLON;
                    med->H5_CDMAG_MLT[ med->H5_nT ]  = MLT;
                    med->H5_CDMAG_R[ med->H5_nT ]    = R;

                    Lgm_Convert_Coords( &Ugsm, &W, GSM_TO_EDMAG, c );
                    Lgm_EDMAG_to_R_MLAT_MLON_MLT( &W, &R, &MLAT, &MLON, &MLT, c );
                    med->H5_EDMAG_MLAT[ med->H5_nT ] = MLAT;
                    med->H5_EDMAG_MLON[ med->H5_nT ] = MLON;
                    med->H5_EDMAG_MLT[ med->H5_nT ]  = MLT;
                    med->H5_EDMAG_R[ med->H5_nT ]    = R;

                    med->H5_Kp[ med->H5_nT ]             = MagEphemInfo->LstarInfo->mInfo->fKp;
                    med->H5_Dst[ med->H5_nT ]            = MagEphemInfo->LstarInfo->mInfo->Dst;

                    med->H5_S_sc_to_pfn[ med->H5_nT ]    = (MagEphemInfo->Snorth > 0.0) ? MagEphemInfo->Snorth : LGM_FILL_VALUE;
                    med->H5_S_sc_to_pfs[ med->H5_nT ]    = (MagEphemInfo->Ssouth > 0.0) ? MagEphemInfo->Ssouth : LGM_FILL_VALUE;
                    med->H5_S_pfs_to_Bmin[ med->H5_nT ]  = (MagEphemInfo->Smin > 0.0) ? MagEphemInfo->Smin : LGM_FILL_VALUE;
                    med->H5_S_Bmin_to_sc[ med->H5_nT ]   = ((MagEphemInfo->Ssouth>0.0)&&(MagEphemInfo->Smin > 0.0)) ? MagEphemInfo->Ssouth-MagEphemInfo->Smin : LGM_FILL_VALUE;
                    med->H5_S_total[ med->H5_nT ]        = ((MagEphemInfo->Snorth > 0.0)&&(MagEphemInfo->Ssouth > 0.0)) ? MagEphemInfo->Snorth + MagEphemInfo->Ssouth : LGM_FILL_VALUE;

                    med->H5_d2B_ds2[ med->H5_nT ]        = MagEphemInfo->d2B_ds2;
                    med->H5_Sb0[ med->H5_nT ]            = MagEphemInfo->Sb0;
                    med->H5_RadiusOfCurv[ med->H5_nT ]   = MagEphemInfo->RofC;


                    MagEphemInfo->LstarInfo->mInfo->Bfield( &Ugsm, &Bsc_gsm, MagEphemInfo->LstarInfo->mInfo );
                    med->H5_Bsc_gsm[ med->H5_nT ][0] = Bsc_gsm.x;
                    med->H5_Bsc_gsm[ med->H5_nT ][1] = Bsc_gsm.y;
                    med->H5_Bsc_gsm[ med->H5_nT ][2] = Bsc_gsm.z;
                    med->H5_Bsc_gsm[ med->H5_nT ][3] = Lgm_Magnitude( &Bsc_gsm );

                    if ( MagEphemInfo->FieldLineType == LGM_CLOSED ) {
                        med->H5_Pmin_gsm[ med->H5_nT ][0] = MagEphemInfo->Pmin.x;
                        med->H5_Pmin_gsm[ med->H5_nT ][1] = MagEphemInfo->Pmin.y;
                        med->H5_Pmin_gsm[ med->H5_nT ][2] = MagEphemInfo->Pmin.z;

                        MagEphemInfo->LstarInfo->mInfo->Bfield( &MagEphemInfo->Pmin, &Bvec, MagEphemInfo->LstarInfo->mInfo );
                        Bmin_mag = Lgm_Magnitude( &Bvec );
                        med->H5_Bmin_gsm[ med->H5_nT ][0] = Bvec.x;
                        med->H5_Bmin_gsm[ med->H5_nT ][1] = Bvec.y;
                        med->H5_Bmin_gsm[ med->H5_nT ][2] = Bvec.z;
                        med->H5_Bmin_gsm[ med->H5_nT ][3] = Bmin_mag;

                    } else {
                        med->H5_Pmin_gsm[ med->H5_nT ][0] = LGM_FILL_VALUE ;
                        med->H5_Pmin_gsm[ med->H5_nT ][1] = LGM_FILL_VALUE ;
                        med->H5_Pmin_gsm[ med->H5_nT ][2] = LGM_FILL_VALUE ;

                        med->H5_Bmin_gsm[ med->H5_nT ][0] = LGM_FILL_VALUE ;
                        med->H5_Bmin_gsm[ med->H5_nT ][1] = LGM_FILL_VALUE ;
                        med->H5_Bmin_gsm[ med->H5_nT ][2] = LGM_FILL_VALUE ;
                        med->H5_Bmin_gsm[ med->H5_nT ][3] = LGM_FILL_VALUE ;
                        Bmin_mag = LGM_FILL_VALUE;
                    }


                    for (i=0; i<nAlpha; i++){
                        med->H5_Lstar[ med->H5_nT ][i]          = MagEphemInfo->Lstar[i];
                        med->H5_DriftShellType[ med->H5_nT ][i] = MagEphemInfo->DriftOrbitType[i];
                        med->H5_Sb[ med->H5_nT ][i]             = MagEphemInfo->Sb[i];
                        med->H5_I[ med->H5_nT ][i]              = MagEphemInfo->I[i];
                        med->H5_Bm[ med->H5_nT ][i]             = MagEphemInfo->Bm[i];

                        Ek    = 1.0; // MeV
                        E     = Ek + LGM_Ee0; // total energy, MeV
                        p2c2  = Ek*(Ek+2.0*LGM_Ee0); // p^2c^2,  MeV^2
                        Beta2 = p2c2/(E*E); // beta^2 = v^2/c^2   (dimensionless)
                        Beta  = sqrt( Beta2 );
                        vel   = Beta*LGM_c;  // velocity in m/s
                        vel  /= (Re*1000.0); // Re/s
                        T     = ( MagEphemInfo->Sb[i] > 0.0 ) ? 2.0*MagEphemInfo->Sb[i]/vel : LGM_FILL_VALUE;

                        pp    = sqrt(p2c2)*1.60217646e-13/LGM_c;  // mks
                        rg    = sin(MagEphemInfo->Alpha[i]*RadPerDeg)*pp/(LGM_e*Bmin_mag*1e-9); // m. Bmin_mag calced above

                        med->H5_Tb[ med->H5_nT ][i]             = T;
                        med->H5_Kappa[ med->H5_nT ][i]          = sqrt( MagEphemInfo->RofC*Re*1e3/rg );


                        if ( (MagEphemInfo->Bm[i]>0.0)&&(MagEphemInfo->I[i]>=0.0) ) {
                            med->H5_K[ med->H5_nT ][i] = 3.16227766e-3*MagEphemInfo->I[i]*sqrt(MagEphemInfo->Bm[i]);
                        } else {
                            med->H5_K[ med->H5_nT ][i] = LGM_FILL_VALUE;
                        }
                        if (MagEphemInfo->I[i]>=0.0) {
                            med->H5_L[ med->H5_nT ][i] = LFromIBmM_McIlwain(MagEphemInfo->I[i], MagEphemInfo->Bm[i], MagEphemInfo->Mused );
                        } else {
                            med->H5_L[ med->H5_nT ][i] = LGM_FILL_VALUE;
                        }
                    }

                    /*
                     * Compute Lsimple
                     */
                    med->H5_Lsimple[ med->H5_nT ] = ( Bmin_mag > 0.0) ? Lgm_Magnitude( &MagEphemInfo->Pmin ) : LGM_FILL_VALUE;

                    /*
                     * Compute InvLat
                     */
                    if (med->H5_Lsimple[ med->H5_nT ] > 0.0) {
                        med->H5_InvLat[ med->H5_nT ] = DegPerRad*acos(sqrt(1.0/med->H5_Lsimple[ med->H5_nT ]));
                    } else {
                        med->H5_InvLat[ med->H5_nT ] = LGM_FILL_VALUE;
                    }

                    /*
                     * Compute Lm_eq
                     */
                    med->H5_Lm_eq[ med->H5_nT ] = (Bmin_mag > 0.0) ? LFromIBmM_McIlwain( 0.0, Bmin_mag, MagEphemInfo->Mcurr ) : LGM_FILL_VALUE;

                    /*
                     * Compute InvLat_eq
                     */
                    if (med->H5_Lm_eq[ med->H5_nT ] > 0.0) {
                        med->H5_InvLat_eq[ med->H5_nT ] = DegPerRad*acos(sqrt(1.0/med->H5_Lm_eq[ med->H5_nT ]));
                    } else {
                        med->H5_InvLat_eq[ med->H5_nT ] = LGM_FILL_VALUE;
                    }

                    /*
                     * Compute BoverBeq
                     */
                    Bsc_mag = Lgm_Magnitude( &Bsc_gsm );
                    med->H5_BoverBeq[ med->H5_nT ] = ( Bmin_mag > 0.0) ? Bsc_mag / Bmin_mag : LGM_FILL_VALUE;

                    /*
                     * Compute MlatFromBoverBeq
                     */
                    if ( med->H5_BoverBeq[ med->H5_nT ] > 0.0 ) {
                        s = sqrt( 1.0/med->H5_BoverBeq[ med->H5_nT ] );
                        cl = Lgm_CdipMirrorLat( s );
                        if ( fabs(cl) <= 1.0 ){
                            med->H5_MlatFromBoverBeq[ med->H5_nT ] = DegPerRad*acos( cl );
                            if (med->H5_S_Bmin_to_sc[ med->H5_nT ]<0.0) med->H5_MlatFromBoverBeq[ med->H5_nT ] *= -1.0;
                        } else {
                            med->H5_MlatFromBoverBeq[ med->H5_nT ] = LGM_FILL_VALUE;
                        }
                    } else {
                        med->H5_MlatFromBoverBeq[ med->H5_nT ] = LGM_FILL_VALUE;
                    }

                    /*
                     * Save M values
                     */
                    med->H5_M_used[ med->H5_nT ] = MagEphemInfo->Mused;
                    med->H5_M_ref[ med->H5_nT ]  = MagEphemInfo->Mref;
                    med->H5_M_igrf[ med->H5_nT ] = MagEphemInfo->Mcurr;


                    if ( (MagEphemInfo->FieldLineType == LGM_CLOSED) || (MagEphemInfo->FieldLineType == LGM_OPEN_N_LOBE) ) {

                        /*
                         * Save northern Footpoint position in different coord systems.
                         */
                        Lgm_VecToArr( &MagEphemInfo->Ellipsoid_Footprint_Pn, med->H5_Pfn_gsm[ med->H5_nT ] );

                        Lgm_Convert_Coords( &MagEphemInfo->Ellipsoid_Footprint_Pn, &W, GSM_TO_GEO, c );
                        Lgm_VecToArr( &W, med->H5_Pfn_geo[ med->H5_nT ] );

                        Lgm_WGS84_to_GEOD( &W, &GeodLat, &GeodLong, &GeodHeight );
                        Lgm_SetArrElements3( med->H5_Pfn_geod[ med->H5_nT ],        GeodLat, GeodLong, GeodHeight );
                        Lgm_SetArrElements2( med->H5_Pfn_geod_LatLon[ med->H5_nT ], GeodLat, GeodLong );
                        med->H5_Pfn_geod_Height[ med->H5_nT ]    = GeodHeight;

                        Lgm_Convert_Coords( &MagEphemInfo->Ellipsoid_Footprint_Pn, &W, GSM_TO_CDMAG, c );
                        Lgm_CDMAG_to_R_MLAT_MLON_MLT( &W, &R, &MLAT, &MLON, &MLT, c );
                        Lgm_SetArrElements3( med->H5_Pfn_cdmag[ med->H5_nT ], MLAT, MLON, MLT );
                        med->H5_Pfn_CD_MLAT[ med->H5_nT ] = MLAT;
                        med->H5_Pfn_CD_MLON[ med->H5_nT ] = MLON;
                        med->H5_Pfn_CD_MLT[ med->H5_nT ]  = MLT;

                        Lgm_Convert_Coords( &MagEphemInfo->Ellipsoid_Footprint_Pn, &W, GSM_TO_EDMAG, c );
                        Lgm_EDMAG_to_R_MLAT_MLON_MLT( &W, &R, &MLAT, &MLON, &MLT, c );
                        Lgm_SetArrElements3( med->H5_Pfn_edmag[ med->H5_nT ], MLAT, MLON, MLT );
                        med->H5_Pfn_ED_MLAT[ med->H5_nT ] = MLAT;
                        med->H5_Pfn_ED_MLON[ med->H5_nT ] = MLON;
                        med->H5_Pfn_ED_MLT[ med->H5_nT ]  = MLT;


                        /*
                         * Save northern Footpoint B-field values in different coord systems.
                         */
                        MagEphemInfo->LstarInfo->mInfo->Bfield( &MagEphemInfo->Ellipsoid_Footprint_Pn, &Bvec, MagEphemInfo->LstarInfo->mInfo );
                        Lgm_Convert_Coords( &Bvec, &Bvec2, GSM_TO_WGS84, c );
                        Lgm_VecToArr( &Bvec,  &med->H5_Bfn_gsm[ med->H5_nT ][0] ); med->H5_Bfn_gsm[ med->H5_nT ][3] = Lgm_Magnitude( &Bvec  );
                        Lgm_VecToArr( &Bvec2, &med->H5_Bfn_geo[ med->H5_nT ][0] ); med->H5_Bfn_geo[ med->H5_nT ][3] = Lgm_Magnitude( &Bvec2 );


                        /*
                         * Save northern loss cone.
                         */
                        Bfn_mag = Lgm_Magnitude( &Bvec );
                        med->H5_LossConeAngleN[ med->H5_nT ] = asin( sqrt( Bsc_mag/Bfn_mag ) )*DegPerRad;


                    } else {

                        Lgm_SetArrVal3( med->H5_Pfn_gsm[ med->H5_nT ],          LGM_FILL_VALUE );
                        Lgm_SetArrVal3( med->H5_Pfn_geo[ med->H5_nT ],          LGM_FILL_VALUE );
                        Lgm_SetArrVal3( med->H5_Pfn_geod[ med->H5_nT ],         LGM_FILL_VALUE );
                        Lgm_SetArrVal2( med->H5_Pfn_geod_LatLon[ med->H5_nT ],   LGM_FILL_VALUE );

                        med->H5_Pfn_geod_Height[ med->H5_nT ]  = LGM_FILL_VALUE;
                        med->H5_Pfn_CD_MLAT[ med->H5_nT ]      = LGM_FILL_VALUE;
                        med->H5_Pfn_CD_MLON[ med->H5_nT ]      = LGM_FILL_VALUE;
                        med->H5_Pfn_CD_MLT[ med->H5_nT ]       = LGM_FILL_VALUE;
                        med->H5_Pfn_ED_MLAT[ med->H5_nT ]      = LGM_FILL_VALUE;
                        med->H5_Pfn_ED_MLON[ med->H5_nT ]      = LGM_FILL_VALUE;
                        med->H5_Pfn_ED_MLT[ med->H5_nT ]       = LGM_FILL_VALUE;

                        Lgm_SetArrVal3( med->H5_Pfn_cdmag[ med->H5_nT ],        LGM_FILL_VALUE );
                        Lgm_SetArrVal3( med->H5_Pfn_edmag[ med->H5_nT ],        LGM_FILL_VALUE );

                        Lgm_SetArrVal4( med->H5_Bfn_geo[ med->H5_nT ],          LGM_FILL_VALUE );
                        Lgm_SetArrVal4( med->H5_Bfn_gsm[ med->H5_nT ],          LGM_FILL_VALUE );

                        med->H5_LossConeAngleN[ med->H5_nT ] = LGM_FILL_VALUE;

                    }

                    if ( (MagEphemInfo->FieldLineType == LGM_CLOSED) || (MagEphemInfo->FieldLineType == LGM_OPEN_S_LOBE) ) {

                        /*
                         * Save southern Footpoint position in different coord systems.
                         */
                        Lgm_VecToArr( &MagEphemInfo->Ellipsoid_Footprint_Ps, med->H5_Pfs_gsm[ med->H5_nT ] );

                        Lgm_Convert_Coords( &MagEphemInfo->Ellipsoid_Footprint_Ps, &W, GSM_TO_GEO, c );
                        Lgm_VecToArr( &W, med->H5_Pfs_geo[ med->H5_nT ] );

                        Lgm_WGS84_to_GEOD( &W, &GeodLat, &GeodLong, &GeodHeight );
                        Lgm_SetArrElements3( med->H5_Pfs_geod[ med->H5_nT ],        GeodLat, GeodLong, GeodHeight );
                        Lgm_SetArrElements2( med->H5_Pfs_geod_LatLon[ med->H5_nT ], GeodLat, GeodLong );
                        med->H5_Pfs_geod_Height[ med->H5_nT ]    = GeodHeight;

                        Lgm_Convert_Coords( &MagEphemInfo->Ellipsoid_Footprint_Ps, &W, GSM_TO_CDMAG, c );
                        Lgm_CDMAG_to_R_MLAT_MLON_MLT( &W, &R, &MLAT, &MLON, &MLT, c );
                        Lgm_SetArrElements3( med->H5_Pfs_cdmag[ med->H5_nT ], MLAT, MLON, MLT );
                        med->H5_Pfs_CD_MLAT[ med->H5_nT ] = MLAT;
                        med->H5_Pfs_CD_MLON[ med->H5_nT ] = MLON;
                        med->H5_Pfs_CD_MLT[ med->H5_nT ]  = MLT;

                        Lgm_Convert_Coords( &MagEphemInfo->Ellipsoid_Footprint_Ps, &W, GSM_TO_EDMAG, c );
                        Lgm_EDMAG_to_R_MLAT_MLON_MLT( &W, &R, &MLAT, &MLON, &MLT, c );
                        Lgm_SetArrElements3( med->H5_Pfs_edmag[ med->H5_nT ], MLAT, MLON, MLT );
                        med->H5_Pfs_ED_MLAT[ med->H5_nT ] = MLAT;
                        med->H5_Pfs_ED_MLON[ med->H5_nT ] = MLON;
                        med->H5_Pfs_ED_MLT[ med->H5_nT ]  = MLT;


                        /*
                         * Save southern Footpoint B-field values in different coord systems.
                         */
                        MagEphemInfo->LstarInfo->mInfo->Bfield( &MagEphemInfo->Ellipsoid_Footprint_Ps, &Bvec, MagEphemInfo->LstarInfo->mInfo );
                        Lgm_Convert_Coords( &Bvec, &Bvec2, GSM_TO_WGS84, c );
                        Lgm_VecToArr( &Bvec,  &med->H5_Bfs_gsm[ med->H5_nT ][0] ); med->H5_Bfs_gsm[ med->H5_nT ][3] = Lgm_Magnitude( &Bvec  );
                        Lgm_VecToArr( &Bvec2, &med->H5_Bfs_geo[ med->H5_nT ][0] ); med->H5_Bfs_geo[ med->H5_nT ][3] = Lgm_Magnitude( &Bvec2 );


                        /*
                         * Save southern loss cone.
                         */
                        Bfs_mag = Lgm_Magnitude( &Bvec );
                        med->H5_LossConeAngleS[ med->H5_nT ] = asin( sqrt( Bsc_mag/Bfs_mag ) )*DegPerRad;


                    } else {

                        Lgm_SetArrVal3( med->H5_Pfs_gsm[ med->H5_nT ],          LGM_FILL_VALUE );
                        Lgm_SetArrVal3( med->H5_Pfs_geo[ med->H5_nT ],          LGM_FILL_VALUE );
                        Lgm_SetArrVal3( med->H5_Pfs_geod[ med->H5_nT ],         LGM_FILL_VALUE );
                        Lgm_SetArrVal2( med->H5_Pfs_geod_LatLon[ med->H5_nT ],   LGM_FILL_VALUE );

                        med->H5_Pfs_geod_Height[ med->H5_nT ]  = LGM_FILL_VALUE;
                        med->H5_Pfs_CD_MLAT[ med->H5_nT ]      = LGM_FILL_VALUE;
                        med->H5_Pfs_CD_MLON[ med->H5_nT ]      = LGM_FILL_VALUE;
                        med->H5_Pfs_CD_MLT[ med->H5_nT ]       = LGM_FILL_VALUE;
                        med->H5_Pfs_ED_MLAT[ med->H5_nT ]      = LGM_FILL_VALUE;
                        med->H5_Pfs_ED_MLON[ med->H5_nT ]      = LGM_FILL_VALUE;
                        med->H5_Pfs_ED_MLT[ med->H5_nT ]       = LGM_FILL_VALUE;

                        Lgm_SetArrVal3( med->H5_Pfs_cdmag[ med->H5_nT ],        LGM_FILL_VALUE );
                        Lgm_SetArrVal3( med->H5_Pfs_edmag[ med->H5_nT ],        LGM_FILL_VALUE );

                        Lgm_SetArrVal4( med->H5_Bfs_geo[ med->H5_nT ],          LGM_FILL_VALUE );
                        Lgm_SetArrVal4( med->H5_Bfs_gsm[ med->H5_nT ],          LGM_FILL_VALUE );

                        med->H5_LossConeAngleS[ med->H5_nT ] = LGM_FILL_VALUE;

                    }

                    /*
                     * Write a row of data into the hdf5 file
                     */
                    Lgm_WriteMagEphemDataHdf( file, med->H5_nT, med->H5_nT, med );
                    ++(med->H5_nT);

            } //end while loop
                fclose(fp_in);
                fclose(fp_MagEphem);

                H5Fclose( file );
            }
            if (isHDF5)
                H5Fclose(file_id);
        } // end birds loop
    } // end JD loop


    free( BirdsPath );
    free( Path );
    free( OutFile );
    free( InFileTmp );
    free( InFile );

    Lgm_free_ctrans( c );
    Lgm_destroy_eop( e );
    Lgm_FreeMagEphemInfo( MagEphemInfo );
    LGM_ARRAY_2D_FREE( BirdsTmp );
    LGM_ARRAY_2D_FREE( Birds );

    return(0);
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
