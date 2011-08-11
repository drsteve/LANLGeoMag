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



#define KP_DEFAULT 0

const  char *argp_program_version     = "MagEphemFromLatLonRad 1.0";
const  char *argp_program_bug_address = "<mghenderson@lanl.gov>";
static char doc[] = "Computes magnetic ephemerii of S/C from input file that contains S/C position in GEO Lat/Lon/Rad";

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
    {"IntModel",        'i',    "internal_model",             0,                                      "Internal Magnetic Field Model to use"    },
    {"ExtModel",        'e',    "external_model",             0,                                      "External Magnetic Field Model to use"    },
    {"PitchAngles",     'p',    "\"start_pa, end_pa, npa\"",  0,                                      "Pitch angles to compute"                 },
    {"FootPointHeight", 'f',    "height",                     0,                                      "Footpoint height in km"                  },
    {"Quality",         'q',    "quality",                    0,                                      "Quality to use for L* calculations"      },
    {"StartDate",       'S',    "yyyymmdd",                   0,                                      "StartDate "                              },
    {"EndDate",         'E',    "yyyymmdd",                   0,                                      "EndDate "                                },
    {"Append",          'a',    0,                            0,                                      "Append to OutFile instead of creating a new one"  },
    {"UseEop",          'e',    0,                            0,                                      "Use Earth Orientation Parameters whn comoputing ephemerii" },
    {"Colorize",        'c',    0,                            0,                                      "Colorize output"                         },
    {"Force",           'F',    0,                            0,                                      "Overwrite output file even if it already exists" },
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

    int         Append;
    int         UseEop;

    long int    StartDate;
    long int    EndDate;
};


/* Parse a single option. */
static error_t parse_opt (int key, char *arg, struct argp_state *state) {

    /* Get the input argument from argp_parse, which we
      know is a pointer to our arguments structure. */
    struct Arguments *arguments = state->input;
    switch( key ) {
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
    Lgm_Vector       Ugsm, Ugeo;
    Lgm_DateTime     UTC;
    Lgm_Eop          *e = Lgm_init_eop( 0 );
    Lgm_EopOne       eop;
    int              i;
    char             IsoTimeString[1024];
    char             InputFilename[1024];
    char             OutputFilename[1024];
    char             IntModel[20], ExtModel[20], Line[5000];
    int              AppendMode, UseEop, Colorize, Force;
    FILE             *fp_in, *fp_MagEphem;
    double           Inc, Alpha[1000], FootpointHeight, GeoLat, GeoLon, GeoRad;
    int              nAlpha, Quality;
    long int         StartDate, EndDate, Date;
    int              sYear, sMonth, sDay, sDoy, eYear, eMonth, eDay, eDoy, Year, Month, Day;
    double           sJD, eJD, JD, Time;
    Lgm_MagEphemInfo *MagEphemInfo;
    Lgm_QinDentonOne p;


    /*
     * Default option values.
     */
    arguments.StartPA         = 90;     // start at 90.0 Deg.
    arguments.EndPA           = 2.5;    // stop at 2.5 Deg.
    arguments.nPA             = 36;     // 36 pitch angles
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


    /*
     *  Print summary of our options
     */
    if ( !arguments.silent ) {
        printf("\n\n");
        printf( "\t        Program/Version: %s\n", argp_program_version );
        printf( "\t         Bug Reports To: %s\n", argp_program_bug_address );
        printf( "\t             Input File: %s\n", InputFilename );
        printf( "\t            Output File: %s\n", OutputFilename );
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
    }


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
    for (i=0; i<nAlpha; i++) MagEphemInfo->Alpha[i] = Alpha[i];


    int SubstituteDates = TRUE;
    if ( (StartDate > 0)&&(EndDate > 0) ){ 
        Lgm_Doy( StartDate, &sYear, &sMonth, &sDay, &sDoy);
        sJD = Lgm_JD( sYear, sMonth, sDay, 12.0, LGM_TIME_SYS_UTC, c );

        Lgm_Doy( EndDate, &eYear, &eMonth, &eDay, &eDoy);
        eJD = Lgm_JD( eYear, eMonth, eDay, 12.0, LGM_TIME_SYS_UTC, c );

        SubstituteDates = TRUE;
    } else {
        // only do a single file
        sJD = eJD = 0;
        SubstituteDates = FALSE;
    }

    char *BaseDir, NewStr[2048], Str[24], Command[4096];
    char *OutFile = (char *)calloc( 2056, sizeof( char ) );
    char *InFile  = (char *)calloc( 2056, sizeof( char ) );

    for ( JD = sJD; JD <= eJD; JD += 1.0 ) {

        Date = Lgm_JD_to_Date( JD, &Year, &Month, &Day, &Time );

        strcpy( InFile, InputFilename );
        strcpy( OutFile, OutputFilename );

        if ( SubstituteDates ) {

            // Substitute times in the files.
            NewStr[0] = '\0';
            sprintf( Str, "%4d", Year );   Lgm_ReplaceSubString( NewStr, InFile, "%YYYY", Str );  strcpy( InFile, NewStr );
            sprintf( Str, "%02d", Month ); Lgm_ReplaceSubString( NewStr, InFile, "%MM", Str );   strcpy( InFile, NewStr );
            sprintf( Str, "%02d", Day );   Lgm_ReplaceSubString( NewStr, InFile, "%DD", Str );   strcpy( InFile, NewStr );

            sprintf( Str, "%4d", Year );   Lgm_ReplaceSubString( NewStr, OutFile, "%YYYY", Str );  strcpy( OutFile, NewStr );
            sprintf( Str, "%02d", Month ); Lgm_ReplaceSubString( NewStr, OutFile, "%MM", Str );   strcpy( OutFile, NewStr );
            sprintf( Str, "%02d", Day );   Lgm_ReplaceSubString( NewStr, OutFile, "%DD", Str );   strcpy( OutFile, NewStr );

        } 



        if ( SubstituteDates ) {
            printf( "\n\n      Processing Date: %4d-%02d-%02d\n", Year, Month, Day );
        } else {
            printf( "\n\n      Processing File: %s\n", InFile );
        }
        printf( "      -------------------------------------------------------------------------------------------\n");
        printf( "           Input File: %s\n", InFile);
        printf( "          Output File: %s\n\n", OutFile);


        // Create Base directory if it hasnt been created yet.
        char *dirc = strdup( OutFile );
        BaseDir    = dirname(dirc);
        sprintf( Command, "mkdir -p %s", BaseDir); system( Command );
        free( dirc );


        /*
         *   Check to see if OutFile exists or not.
         */
        int     StatError, FileExists;
        struct stat StatBuf;
        StatError = stat( OutFile, &StatBuf );

        if ( StatError != -1 ) {

            FileExists = TRUE;
            printf("\n\n\tOutfile already exists: %s\n", OutFile );
            if ( StatBuf.st_size < 1000LL ){
                printf("\t             File size:  %lld B\n", (long long)StatBuf.st_size );
            } else if ( StatBuf.st_size < 1000000LL ){
                printf("\t             File size:  %g kB\n", (double)StatBuf.st_size/1.0e3 );
            } else if ( StatBuf.st_size < 1000000000LL ){
                printf("\t             File size:  %g MB\n", (double)StatBuf.st_size/1.0e6 );
            } else {
                printf("\t             File size:  %g GB\n", (double)StatBuf.st_size/1.0e9 );
            }
            printf("\t    Last status change:  %s", ctime(&StatBuf.st_ctime));
            printf("\t      Last file access:  %s", ctime(&StatBuf.st_atime));
            printf("\tLast file modification:  %s", ctime(&StatBuf.st_mtime));

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
                fp_MagEphem = fopen( OutFile, "wb" );
                Lgm_WriteMagEphemHeader( fp_MagEphem, "FIX ME", 99999, "FIX ME", IntModel, ExtModel, MagEphemInfo );
                printf("\t      Writing to file: %s\n", OutFile );

                if ( UseEop ) {
                    // Read in the EOP vals
                    Lgm_read_eop( e );
                }


                /*
                 *  Loop over the times/positions given. Here we assume the file contains, Time and GEO Lat (Deg.), Lon (Deg.), Rad (Re)
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


                        // Converting GEO->GSM coords.
                        Ugeo.x = GeoRad*cos(GeoLon*RadPerDeg)*cos(GeoLat*RadPerDeg);
                        Ugeo.y = GeoRad*sin(GeoLon*RadPerDeg)*cos(GeoLat*RadPerDeg);
                        Ugeo.z = GeoRad*sin(GeoLat*RadPerDeg);
                        Lgm_Convert_Coords( &Ugeo, &Ugsm, GEO_TO_GSM, c );


                        /*
                         * Compute L*s, Is, Bms, Footprints, etc...
                         * These quantities are stored in the MagEphemInfo Structure
                         */
                        printf("\n\n\t\tDate, UTC: %ld %g   Ugsm: %g %g %g \n", UTC.Date, UTC.Time, Ugsm.x, Ugsm.y, Ugsm.z );
                        printf("\t\t------------------------------------------------------------------\n");
                        Lgm_ComputeLstarVersusPA( UTC.Date, UTC.Time, &Ugsm, nAlpha, Alpha, MagEphemInfo->LstarQuality, Colorize, MagEphemInfo );

                        Lgm_WriteMagEphemData( fp_MagEphem, IntModel, ExtModel, MagEphemInfo->LstarInfo->mInfo->fKp, MagEphemInfo->LstarInfo->mInfo->Dst, MagEphemInfo );

                        if ( nAlpha > 0 ){
                            WriteMagEphemInfoStruct( "test.dat", nAlpha, MagEphemInfo );
                        }

                    }


                }
                fclose(fp_in);
                fclose(fp_MagEphem);
            }

        }
    }

    
    free( OutFile );
    free( InFile );

    Lgm_free_ctrans( c );
    Lgm_destroy_eop( e );
    Lgm_FreeMagEphemInfo( MagEphemInfo );



    return(0);
}

