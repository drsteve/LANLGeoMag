#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <ctype.h>
#include <argp.h>
#include <time.h>
#include <libgen.h>
#include <Lgm_CTrans.h>
#include <Lgm_Sgp.h>
#include <Lgm_MagEphemInfo.h>
#include <Lgm_DynamicMemory.h>
#include <Lgm_Eop.h>
#include <Lgm_QinDenton.h>
#include <Lgm_Misc.h>
#include <Lgm_HDF5.h>

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
    {"verbose",         'v',    0,                            OPTION_ARG_OPTIONAL,                    "Produce verbose output"                  },
    {"silent",          's',    0,                            OPTION_ARG_OPTIONAL | OPTION_ALIAS                                                },
    { 0 }
};

struct Arguments {
    char        *args[ nArgs ];       /* START_PA  END_PA  PA_INC */
    int         silent;
    int         verbose;

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
    long int    StartDate;
    long int    EndDate;
};


/* Parse a single option. */
static error_t parse_opt (int key, char *arg, struct argp_state *state) {

    /* Get the input argument from argp_parse, which we
      know is a pointer to our arguments structure. */
    struct Arguments *arguments = state->input;
    switch( key ) {
        case 'b': // Birds
            strncpy( arguments->Birds, arg, 4095 );
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
            arguments->verbose = 1;
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
    Lgm_CTrans       *c = Lgm_init_ctrans( 0 );
    Lgm_Vector       Ugsm, U;
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
    int              nBirds, iBird;
    char             **Birds, Bird[80];
    double           Inc, Alpha[1000], FootpointHeight, GeoLat, GeoLon, GeoRad;
    int              nAlpha, Quality;
    long int         StartDate, EndDate, Date;
    int              sYear, sMonth, sDay, sDoy, eYear, eMonth, eDay, eDoy, Year, Month, Day;
    double           sJD, eJD, JD, Time;
    Lgm_MagEphemInfo *MagEphemInfo;
    Lgm_QinDentonOne p;

    hid_t           file;
    hid_t           space;
    hid_t           atype;
    hid_t           DataSet, MemSpace;
    herr_t          status;
    hsize_t         Dims[4], Offset[4], SlabSize[4];
    int             iT;
    char            **H5_IsoTimes;
    double          *H5_JD;
    double          UGSM_ARR[3];
    Lgm_Vector      *H5_Ugsm;
    int             H5_nT;
    int             H5_nAlpha;
    double          *H5_Alpha;
    double          **H5_Lstar;
    double          **H5_K;

    // kludge.
    LGM_ARRAY_2D( H5_IsoTimes, 2000, 80, char );
    LGM_ARRAY_1D( H5_JD, 2000, double );
    LGM_ARRAY_1D( H5_Ugsm, 2000, Lgm_Vector );


    /*
     * Default option values.
     */
    arguments.StartPA         = 90;     // start at 90.0 Deg.
    arguments.EndPA           = 5.0;    // stop at 2.5 Deg.
    arguments.nPA             = 18;     // 18 pitch angles
    arguments.silent          = 0;
    arguments.verbose         = 0;
    arguments.Quality         = 3;
    arguments.Colorize        = 0;
    arguments.Force           = 0;
    arguments.Append          = 0;
    arguments.UseEop          = 0;
    arguments.StartDate       = -1;
    arguments.EndDate         = -1;
    arguments.FootPointHeight = 100.0; // km
    strcpy( arguments.IntModel, "IGRF" );
    strcpy( arguments.ExtModel, "T89" );
    strcpy( arguments.CoordSystem, "LATLONRAD" );


    /*
     *  Parse CmdLine arguments and options
     */
    argp_parse (&argp, argc, argv, 0, 0, &arguments);


    /*
     * Define pitch angles to use.
     */
    nAlpha = arguments.nPA;
    Inc    = (arguments.EndPA - arguments.StartPA)/(double)(nAlpha-1);
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

    LGM_ARRAY_2D( Birds, 20, 80, char );
    StringSplit( arguments.Birds, Birds, 80, &nBirds );

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
        //printf( "\t   Append to OutputFile: %s\n", AppendMode ? "yes" : "no" );
        printf( "\t                Use Eop: %s\n", UseEop ? "yes" : "no" );
        printf( "\t Colorize Thread Output: %d\n", Colorize );
        printf( "\t                Verbose: %s\n", arguments.verbose ? "yes" : "no" );
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

exit(0);

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



    /*
     * Set pitch angles in MagEphemInfo structure
     */
    MagEphemInfo->nAlpha = nAlpha;
    LGM_ARRAY_1D( H5_Alpha, nAlpha, double );
    H5_nAlpha = nAlpha;
    for (i=0; i<nAlpha; i++) {
        MagEphemInfo->Alpha[i] = Alpha[i];
        H5_Alpha[i] = Alpha[i];
    }
    LGM_ARRAY_2D( H5_Lstar, 2000, nAlpha, double );
    LGM_ARRAY_2D( H5_K,     2000, nAlpha, double );



    int SubstituteVars = TRUE;
    if ( (StartDate > 0)&&(EndDate > 0) ){
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
    char *InFile     = (char *)calloc( 2056, sizeof( char ) );

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
                if ( !Force ) {
                    printf("\n\n\tHdfOutfile already exists (use -F option to force processing): %s\n", HdfOutFile );
                } else {
                    printf("\n\n\tWarning. Existing HdfOutfile will be overwritten: %s \n", HdfOutFile );
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

            if ( !FileExists || Force ) {


                /*
                 * Open input file for reading
                 */
                if ( (fp_in = fopen( InFile, "r" ) ) == NULL ){

                    printf("\tCould not open file %s for reading\n", InFile );

                } else {

                    printf( "\t    Reading from file: %s\n", InFile );


                    /*
                     * Open Mag Ephem file for writing
                     */
                    fp_MagEphem = fopen( OutFile, "w" );
                    Lgm_WriteMagEphemHeader( fp_MagEphem, "FIX ME", 99999, "FIX ME", NULL, 0, NULL, NULL, 0, NULL, NULL, MagEphemInfo );
                    printf("\t      Writing to file: %s\n", OutFile );

                    if ( UseEop ) {
                        // Read in the EOP vals
                        Lgm_read_eop( e );
                    }


                    /*
                     *  Loop over the times/positions given. Here we assume the file contains, Time and 3 columns: X, Y, Z.
                     *  The meaning of the X, Y, Z columns are dependent upon what CoordSystem is set to;
                     *
                     *        CoordSystem          X                    Y               Z
                     *        ------------------------------------------------------------------
                     *          LATLONRAD       Latitude (Deg.)   Longitude (Deg.)  Radius (Re)
                     *             SM              Xsm                 Ysm             Zsm
                     *            GSM              Xsm                 Ysm             Zsm
                     *            GEI2000          Xgsm                Ygsm            Zgsm
                     *            GSE              Xgse                Ygse            Zgse
                     *
                     *      Add more as needed....
                     */
                    H5_nT = 0;
                    while ( fgets( Line, 4096, fp_in ) != NULL ) {

                        if ( Line[0] != '#' ) {

                            sscanf( Line, "%s %lf %lf %lf", IsoTimeString, &X, &Y, &Z );

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
                            Lgm_ComputeLstarVersusPA( UTC.Date, UTC.Time, &Ugsm, nAlpha, Alpha, MagEphemInfo->LstarQuality, Colorize, MagEphemInfo );

                            Lgm_WriteMagEphemData( fp_MagEphem, IntModel, ExtModel, MagEphemInfo->LstarInfo->mInfo->fKp, MagEphemInfo->LstarInfo->mInfo->Dst, MagEphemInfo );


                            if ( nAlpha > 0 ){
                                WriteMagEphemInfoStruct( "test.dat", nAlpha, MagEphemInfo );
                            }




                            // Fill arrays for dumping out as HDF5 files
                            strcpy( H5_IsoTimes[H5_nT], IsoTimeString );
                            H5_JD[H5_nT]   = UTC.JD;
                            H5_Ugsm[H5_nT] = Ugsm;
                            for (i=0; i<nAlpha; i++){
                                H5_Lstar[H5_nT][i] = MagEphemInfo->Lstar[i];
                                if (MagEphemInfo->Bm[i]>0.0) {
                                    H5_K[H5_nT][i] = 3.16227766e-3*MagEphemInfo->I[i]*sqrt(MagEphemInfo->Bm[i]);
                                } else {
                                    H5_K[H5_nT][i] = LGM_FILL_VALUE;
                                }
                            }


                            ++H5_nT;

                        }


                    }
                    fclose(fp_in);
                    fclose(fp_MagEphem);


                    // Create HDF5 file
                    file    = H5Fcreate( HdfOutFile, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );

                    // Create IsoTime Dataset
                    Dims[0] = H5_nT;
                    space   = H5Screate_simple( 1, Dims, NULL ); // rank 1
                    atype   = H5Tcopy( H5T_C_S1 );
                    status  = H5Tset_size( atype, strlen(H5_IsoTimes[0]) );
                    status  = H5Tset_strpad( atype, H5T_STR_NULLPAD );
                    status  = H5Tset_cset( atype, H5T_CSET_ASCII );
                    DataSet = H5Dcreate( file, "IsoTime", atype, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );
                    status  = H5Tclose( atype );

                    // Create JD Dataset
                    Dims[0] = H5_nT;
                    space   = H5Screate_simple( 1, Dims, NULL ); // rank 1
                    DataSet = H5Dcreate( file, "JulianDate", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );

                    // Create Ugsm Dataset
                    Dims[0] = H5_nT; Dims[1] = 3;
                    space   = H5Screate_simple( 2, Dims, NULL ); // rank 2
                    DataSet = H5Dcreate( file, "Ugsm", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );

                    // Create Alpha Dataset
                    Dims[0] = H5_nAlpha;
                    space   = H5Screate_simple( 1, Dims, NULL ); // rank 1
                    DataSet = H5Dcreate( file, "Alpha", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );

                    // Create Lstar Dataset
                    Dims[0] = H5_nT; Dims[1] = H5_nAlpha;
                    space   = H5Screate_simple( 2, Dims, NULL ); // rank 2
                    DataSet = H5Dcreate( file, "Lstar", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );

                    // Create K Dataset
                    Dims[0] = H5_nT; Dims[1] = H5_nAlpha;
                    space   = H5Screate_simple( 2, Dims, NULL ); // rank 2
                    DataSet = H5Dcreate( file, "K", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );




                    // Write IsoTime Strings
                    SlabSize[0] = 1;
                    DataSet  = H5Dopen( file, "IsoTime", H5P_DEFAULT );
                    space    = H5Dget_space( DataSet );
                    MemSpace = H5Screate_simple( 1, SlabSize, NULL );
                    atype    = H5Tcopy( H5T_C_S1 );
                    status  = H5Tset_size( atype, strlen(H5_IsoTimes[0]) );
                    status  = H5Tset_strpad( atype, H5T_STR_NULLPAD );
                    status  = H5Tset_cset( atype, H5T_CSET_ASCII );
                    for ( iT=0; iT<H5_nT; iT++ ){
                            Offset[0]   = iT;
                            status = H5Sselect_hyperslab( space, H5S_SELECT_SET, Offset, NULL, SlabSize, NULL );
                            status = H5Dwrite( DataSet, atype, MemSpace, space, H5P_DEFAULT, &H5_IsoTimes[iT][0] );
                    }
                    status = H5Sclose( MemSpace );
                    status = H5Dclose( DataSet );
                    status = H5Sclose( space );
                    status = H5Tclose( atype );

                    // Write JD
                    SlabSize[0] = H5_nT;
                    Offset[0]   = 0;
                    DataSet  = H5Dopen( file, "JulianDate", H5P_DEFAULT );
                    space    = H5Dget_space( DataSet );
                    status   = H5Sselect_hyperslab( space, H5S_SELECT_SET, Offset, NULL, SlabSize, NULL );
                    MemSpace = H5Screate_simple( 1, SlabSize, NULL );
                    status   = H5Dwrite( DataSet, H5T_NATIVE_DOUBLE, MemSpace, space, H5P_DEFAULT, &H5_JD[0] );
                    status   = H5Sclose( MemSpace );
                    status   = H5Sclose( space );
                    status   = H5Dclose( DataSet );

                    // Write Ugsm
                    SlabSize[0] = 1; SlabSize[1] = 3;
                    DataSet  = H5Dopen( file, "Ugsm", H5P_DEFAULT );
                    space    = H5Dget_space( DataSet );
                    MemSpace = H5Screate_simple( 2, SlabSize, NULL );
                    for ( iT=0; iT<H5_nT; iT++ ){
                            UGSM_ARR[0] = H5_Ugsm[iT].x; UGSM_ARR[1] = H5_Ugsm[iT].y; UGSM_ARR[2] = H5_Ugsm[iT].z;
                            Offset[0]   = iT; Offset[1] = 0;
                            status = H5Sselect_hyperslab( space, H5S_SELECT_SET, Offset, NULL, SlabSize, NULL );
                            status = H5Dwrite( DataSet, H5T_NATIVE_DOUBLE, MemSpace, space, H5P_DEFAULT, &UGSM_ARR[0] );
                    }
                    status = H5Sclose( MemSpace );
                    status = H5Dclose( DataSet );
                    status = H5Sclose( space );

                    // Write Alpha
                    SlabSize[0] = H5_nAlpha;
                    Offset[0]   = 0;
                    DataSet  = H5Dopen( file, "Alpha", H5P_DEFAULT );
                    space    = H5Dget_space( DataSet );
                    status   = H5Sselect_hyperslab( space, H5S_SELECT_SET, Offset, NULL, SlabSize, NULL );
                    MemSpace = H5Screate_simple( 1, SlabSize, NULL );
                    status   = H5Dwrite( DataSet, H5T_NATIVE_DOUBLE, MemSpace, space, H5P_DEFAULT, &H5_Alpha[0] );
                    status   = H5Sclose( MemSpace );
                    status   = H5Sclose( space );
                    status   = H5Dclose( DataSet );


                    // Write Lstar
                    SlabSize[0] = H5_nT; SlabSize[1] = H5_nAlpha;
                    Offset[0]   = 0; Offset[1] = 0;
                    DataSet  = H5Dopen( file, "Lstar", H5P_DEFAULT );
                    space    = H5Dget_space( DataSet );
                    status   = H5Sselect_hyperslab( space, H5S_SELECT_SET, Offset, NULL, SlabSize, NULL );
                    MemSpace = H5Screate_simple( 2, SlabSize, NULL );
                    status   = H5Dwrite( DataSet, H5T_NATIVE_DOUBLE, MemSpace, space, H5P_DEFAULT, &H5_Lstar[0][0] );
                    status   = H5Sclose( MemSpace );
                    status   = H5Sclose( space );
                    status   = H5Dclose( DataSet );

                    // Write K
                    SlabSize[0] = H5_nT; SlabSize[1] = H5_nAlpha;
                    Offset[0]   = 0; Offset[1] = 0;
                    DataSet  = H5Dopen( file, "K", H5P_DEFAULT );
                    space    = H5Dget_space( DataSet );
                    status   = H5Sselect_hyperslab( space, H5S_SELECT_SET, Offset, NULL, SlabSize, NULL );
                    MemSpace = H5Screate_simple( 2, SlabSize, NULL );
                    status   = H5Dwrite( DataSet, H5T_NATIVE_DOUBLE, MemSpace, space, H5P_DEFAULT, &H5_K[0][0] );
                    status   = H5Sclose( MemSpace );
                    status   = H5Sclose( space );
                    status   = H5Dclose( DataSet );


                    H5Fclose( file );
                }

            }
        } // end birds loop
    } // end JD loop


    free( OutFile );
    free( InFile );

    Lgm_free_ctrans( c );
    Lgm_destroy_eop( e );
    Lgm_FreeMagEphemInfo( MagEphemInfo );
    LGM_ARRAY_2D_FREE( Birds );

    LGM_ARRAY_2D_FREE( H5_IsoTimes );
    LGM_ARRAY_1D_FREE( H5_JD );
    LGM_ARRAY_1D_FREE( H5_Ugsm );
    LGM_ARRAY_1D_FREE( H5_Alpha );
    LGM_ARRAY_2D_FREE( H5_Lstar );
    LGM_ARRAY_2D_FREE( H5_K );



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
