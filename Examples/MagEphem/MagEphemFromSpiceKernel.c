#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
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
#include "SpiceUsr.h"

#define EARTH_ID     399
#define MOON_ID      301
#define RBSPA_ID    -362
#define RBSPB_ID    -363

void StringSplit( char *Str, char *StrArray[], int len, int *n );


#define KP_DEFAULT 0

const  char *ProgramName = "MagEphemFromFile";
const  char *argp_program_version     = "MagEphemFromFile_1.1";
const  char *argp_program_bug_address = "<mghenderson@lanl.gov>";
static char doc[] =
"Computes the magnetic ephemeris of a S/C from trajectories determined from SPICE kernel files.\n\n"
"The input file is a SPICE kernel description file, which is a file describing all of the SPICE kernels that need to be loaded.  A sample kernel description file is written as;\n\n"

"\t\\begindata\n"
"\tKERNELS_TO_LOAD = ( 'RBSPA_2012_139_2014_213_01_deph.bsp', 'naif0009.tls' )\n"
"\t\\begintext\n\n\n"

"Where 'RBSPA_2012_139_2014_213_01_deph.bsp' is a \"binary SPK\" SPICE kernel, and 'naif0009.tls' is a \"text leap second kernel (LSK)\" file.  In SPICE, SPK kernel files contain information that is used to compute the ephemeris (trajectory) of an object.\n\n"

"The output is written into daily files with filename given by the specified template. OutFile may contain time variables that will be substituted.  The available time variables are '%YYYY', '%MM', and '%DD' which correspond repectively to 4-digit year, 2-digit month (Jan is 01), and 2-digit day of month. \n\n"

"The %B variable will also get substituted by the list of birds given in the -b option.  Here is an example using time-variables.\n\n"

"\t/MagEphemFromSpiceKernelFile -S 20020901 -E 20020930 setup.ker \n"
"\t\t/home/jsmith/MagEphemData/%YYYY/%YYYY%MM%DD_1989-046_MagEphem.txt.  \n"
"\n"

"Directories in the output file will be created if they don't already exist.\n";



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
 * Returns Radial distance of S/C.
 */
typedef struct afInfo {
    long int        Date;
    int             Sgn;
    Lgm_CTrans      *c;
} afInfo;
double ApogeeFunc( double T, double val, void *Info ){

    afInfo          *a;
    Lgm_CTrans      *c;
    Lgm_DateTime    UTC;
    Lgm_Vector      U;
    double          R, et, pos[3], lt;

    a = (afInfo *)Info;
    c = a->c;

    Lgm_Make_UTC( a->Date, T/3600.0, &UTC, c );
    et = Lgm_TDBSecSinceJ2000( &UTC, c );
    spkezp_c( RBSPA_ID,    et,   "J2000",  "NONE", EARTH_ID,  pos,  &lt );
    U.x = pos[0]/WGS84_A; U.y = pos[1]/WGS84_A; U.z = pos[2]/WGS84_A;
    R = Lgm_Magnitude( &U );

    return( R*a->Sgn );

}


/*
 *  Compute magnetic ephemerii of S/C from input file that contains S/C position in GEO LAt/Lon/Rad
 */
int main( int argc, char *argv[] ){

    struct Arguments arguments;
    Lgm_CTrans       *c = Lgm_init_ctrans( 0 );
    Lgm_Vector       Rgsm, W, U;
    Lgm_DateTime     UTC;
    Lgm_Eop          *e = Lgm_init_eop( 0 );
    Lgm_EopOne       eop;
    int              i;
    char             IsoTimeString[1024];
    char             InputFilename[1024];
    char             OutputFilename[1024];
    char             IntModel[20], ExtModel[20], CoordSystem[80];
    int              UseEop, Colorize, Force;
    FILE             *fp_in, *fp_MagEphem;
    int              nBirds, iBird;
    char             **Birds, Bird[80];
    double           Inc, Alpha[1000], FootpointHeight;
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
    double          GeodLat, GeodLong, GeodHeight, MLAT, MLON, MLT;
    char            **H5_IsoTimes;
    long int        *H5_Date;
    int             *H5_Doy;
    double          *H5_UTC;
    double          *H5_JD;
    double          *H5_GpsTime;
    double          *H5_TiltAngle;
    double          U_ARR[4];
    Lgm_Vector      *H5_Rgsm;
    Lgm_Vector      *H5_Bsc_gsm;
    Lgm_Vector      *H5_Bfn_geo;
    Lgm_Vector      *H5_Bfn_gsm;
    Lgm_Vector      *H5_Pfn_geo;
    Lgm_Vector      *H5_Pfn_gsm;
    Lgm_Vector      *H5_Pfn_geod;
    Lgm_Vector      *H5_Pfn_cdmag;
    Lgm_Vector      *H5_Pfn_edmag;
    Lgm_Vector      *H5_Bfs_geo;
    Lgm_Vector      *H5_Bfs_gsm;
    Lgm_Vector      *H5_Pfs_geo;
    Lgm_Vector      *H5_Pfs_gsm;
    Lgm_Vector      *H5_Pfs_geod;
    Lgm_Vector      *H5_Pfs_cdmag;
    Lgm_Vector      *H5_Pfs_edmag;
    double          *H5_Kp;
    double          *H5_Dst;
    int             H5_nT;
    int             H5_nAlpha;
    double          *H5_Alpha;
    double          **H5_Lstar;
    double          **H5_K;
    double          **H5_I;
    double          **H5_L;
    double          **H5_Bm;

    int             done;
    long int        Seconds, Ta, Tb, Tc;
    double          R, Ra, Rb, Rc, Rmin, Tmin;
    BrentFuncInfo   bInfo;
    afInfo          *afi;
    int             nPerigee, nApogee;
    Lgm_DateTime    *Perigee_UTC;
    Lgm_DateTime    *Apogee_UTC;
    Lgm_Vector      *Perigee_U;
    Lgm_Vector      *Apogee_U;
    Lgm_Vector      Bvec, Bvec2;

    // kludge.
    LGM_ARRAY_2D( H5_IsoTimes,  2000, 80,    char );
    LGM_ARRAY_1D( H5_Date,      2000,        long int );
    LGM_ARRAY_1D( H5_Doy,       2000,        int );
    LGM_ARRAY_1D( H5_UTC,       2000,        double );
    LGM_ARRAY_1D( H5_JD,        2000,        double );
    LGM_ARRAY_1D( H5_GpsTime,   2000,        double );
    LGM_ARRAY_1D( H5_TiltAngle, 2000,        double );
    LGM_ARRAY_1D( H5_Rgsm,      2000,        Lgm_Vector );
    LGM_ARRAY_1D( H5_Bsc_gsm,   2000,        Lgm_Vector );
    LGM_ARRAY_1D( H5_Bfn_geo,   2000,        Lgm_Vector );
    LGM_ARRAY_1D( H5_Bfn_gsm,   2000,        Lgm_Vector );
    LGM_ARRAY_1D( H5_Pfn_geo,   2000,        Lgm_Vector );
    LGM_ARRAY_1D( H5_Pfn_gsm,   2000,        Lgm_Vector );
    LGM_ARRAY_1D( H5_Pfn_geod,  2000,        Lgm_Vector );
    LGM_ARRAY_1D( H5_Pfn_cdmag, 2000,        Lgm_Vector );
    LGM_ARRAY_1D( H5_Pfn_edmag, 2000,        Lgm_Vector );
    LGM_ARRAY_1D( H5_Bfs_geo,   2000,        Lgm_Vector );
    LGM_ARRAY_1D( H5_Bfs_gsm,   2000,        Lgm_Vector );
    LGM_ARRAY_1D( H5_Pfs_geo,   2000,        Lgm_Vector );
    LGM_ARRAY_1D( H5_Pfs_gsm,   2000,        Lgm_Vector );
    LGM_ARRAY_1D( H5_Pfs_geod,  2000,        Lgm_Vector );
    LGM_ARRAY_1D( H5_Pfs_cdmag, 2000,        Lgm_Vector );
    LGM_ARRAY_1D( H5_Pfs_edmag, 2000,        Lgm_Vector );
    LGM_ARRAY_1D( H5_Kp,        2000,        double );
    LGM_ARRAY_1D( H5_Dst,       2000,        double );

    LGM_ARRAY_1D( Perigee_UTC, 30, Lgm_DateTime );
    LGM_ARRAY_1D( Apogee_UTC,  30, Lgm_DateTime );
    LGM_ARRAY_1D( Perigee_U, 30, Lgm_Vector );
    LGM_ARRAY_1D( Apogee_U,  30, Lgm_Vector );

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
        printf( "\t               Program/Version: %s\n", argp_program_version );
        printf( "\t                Bug Reports To: %s\n", argp_program_bug_address );
        printf( "\tInput SPICE Kernel Descr. File: %s\n", InputFilename );
        printf( "\t                   Output File: %s\n", OutputFilename );
        printf( "\t            Input Coord System: %s\n", CoordSystem );
        printf( "\t        Number of Pitch Angles: %d\n", nAlpha );
        printf( "\t        Pitch Angles [Degrees]:" );
        for (i=0; i<nAlpha; i++ ) printf( " %g", Alpha[i] );
        printf("\n");
        printf( "\t                Internal Model: %s\n", IntModel );
        printf( "\t                External Model: %s\n", ExtModel );
        printf( "\t          FootpointHeight [km]: %g\n", FootpointHeight );
        printf( "\t                    L* Quality: %d\n", Quality );
        printf( "\t                  Force output: %s\n", Force ? "yes" : "no" );
        printf( "\t                       Use Eop: %s\n", UseEop ? "yes" : "no" );
        printf( "\t        Colorize Thread Output: %d\n", Colorize );
        printf( "\t                       Verbose: %s\n", arguments.verbose ? "yes" : "no" );
        printf( "\t                        Silent: %s\n", arguments.silent  ? "yes" : "no" );
        if ( (StartDate > 0)&&(EndDate > 0) ){
            printf( "\t                     StartDate: %ld\n", StartDate );
            printf( "\t                       EndDate: %ld\n", EndDate );
        }

        if ( nBirds > 1  ){
            printf( "\t                         Birds:" );
            for (iBird=0; iBird<nBirds-1; iBird++) printf(" %s,", Birds[iBird]);
            printf(" %s\n", Birds[iBird]);
        } else if ( nBirds == 1  ){
            printf( "\t                          Bird:" );
            printf(" %s\n", Birds[0]);
        } else {
            printf("At least one bird name must be supplied with -b option.\n");
            exit(0);
        }
    }


SpiceDouble et;
SpiceDouble pos[3], lt;

/*
furnsh_c( InputFilename );

Lgm_DateTime *myUTC;

sprintf(time, "2013-01-01T00:00:00", StartDate );
printf("time = %s\n", time);
myUTC = Lgm_DateTime_Create( 2013, 1, 1, 0.0, LGM_TIME_SYS_UTC, c );
et = Lgm_TDBSecSinceJ2000( myUTC, c );
printf("et = %.15lf\n", et);



str2et_c ( time, &et );
printf("et = %.15lf\n", et);



spkezp_c( RBSPA_ID,    et,   "J2000",  "NONE", EARTH_ID,  pos,  &lt );
printf("pos = %.8lf %.8lf %.8lf\n", pos[0]/WGS84_A, pos[1]/WGS84_A, pos[2]/WGS84_A);
*/







//exit(0);

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
    LGM_ARRAY_2D( H5_I,     2000, nAlpha, double );
    LGM_ARRAY_2D( H5_L,     2000, nAlpha, double );
    LGM_ARRAY_2D( H5_Bm,    2000, nAlpha, double );



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

                    printf("\tCould not open SPICE Kernel Desrciption File ( %s ) for reading\n", InFile );

                } else {

                    fclose( fp_in );

                    /*
                     * Make kernels known to SPICE
                     */
                    furnsh_c( InputFilename );
                    printf( "\t    Using SPICE Kernel Desrciption File: %s\n", InFile );



                    /*
                     * Find Apogees.
                     * Pack aInfo structure with required info.
                     * Find brackets.
                     * Call Lgm_Brent() to get minimum.
                     *
                     */
                    afi = (afInfo *)calloc( 1, sizeof(afInfo) );
                    afi->Date = Date;
                    afi->Sgn  = -1; // +1 => find min (perigee)   -1 => find max (apogee)
                    afi->c    = c;
                    bInfo.Val  = 0.0;
                    bInfo.Info = (void *)afi;
                    bInfo.func = ApogeeFunc;
                    Ta = 0;    Ra = bInfo.func( (double)Ta, 0.0, (void *)afi );
                    Tb = 900;  Rb = bInfo.func( (double)Tb, 0.0, (void *)afi );
                    Tc = 1800; Rc = bInfo.func( (double)Tc, 0.0, (void *)afi );
                    done = FALSE;
                    nApogee = 0;
                    while ( !done ) {
                        if ( (Rb < Ra) && (Rb < Rc) ) {
                            // we have a bracket. Find Min.
                            Lgm_Brent( (double)Ta, (double)Tb, (double)Tc, &bInfo, 1e-3, &Tmin, &Rmin );
                            Lgm_Make_UTC( Date, Tmin/3600.0, &Apogee_UTC[nApogee], c );
                            et = Lgm_TDBSecSinceJ2000( &Apogee_UTC[nApogee], c );
                            spkezp_c( RBSPA_ID,    et,   "J2000",  "NONE", EARTH_ID,  pos,  &lt );
                            Apogee_U[nApogee].x = pos[0]/WGS84_A; Apogee_U[nApogee].y = pos[1]/WGS84_A; Apogee_U[nApogee].z = pos[2]/WGS84_A;
                            ++nApogee;

                            //Lgm_DateTimeToString( IsoTimeString, &UTC, 0, 3 );
                            //printf("Bracket: T = %ld %ld %ld    R = %g %g %g    Tapogee, Rapogee = %s %g\n", Ta, Tb, Tc, Ra, Rb, Rc, IsoTimeString, fabs(Rmin) );
                        }
                        // advance another 900 seconds.
                        Ta = Tb; Ra = Rb;
                        Tb = Tc; Rb = Rc;
                        Tc += 900.0; Rc = bInfo.func( (double)Tc, 0.0, (void *)afi );
                        if (Tc >= 86400) done = TRUE;
                    }


                    /*
                     * Find Perigees. Same as above except set afi->Sgn = 1.
                     */
                    afi->Sgn  = 1; // +1 => find min (perigee)   -1 => find max (apogee)
                    Ta = 0;    Ra = bInfo.func( (double)Ta, 0.0, (void *)afi );
                    Tb = 900;  Rb = bInfo.func( (double)Tb, 0.0, (void *)afi );
                    Tc = 1800; Rc = bInfo.func( (double)Tc, 0.0, (void *)afi );
                    done = FALSE;
                    nPerigee = 0;
                    while ( !done ) {
                        if ( (Rb < Ra) && (Rb < Rc) ) {
                            // we have a bracket. Find Min.
                            Lgm_Brent( (double)Ta, (double)Tb, (double)Tc, &bInfo, 1e-10, &Tmin, &Rmin );
                            Lgm_Make_UTC( Date, Tmin/3600.0, &Perigee_UTC[nPerigee], c );
                            et = Lgm_TDBSecSinceJ2000( &Perigee_UTC[nPerigee], c );
                            spkezp_c( RBSPA_ID,    et,   "J2000",  "NONE", EARTH_ID,  pos,  &lt );
                            Perigee_U[nPerigee].x = pos[0]/WGS84_A; Perigee_U[nPerigee].y = pos[1]/WGS84_A; Perigee_U[nPerigee].z = pos[2]/WGS84_A;
                            ++nPerigee;
                            //Lgm_DateTimeToString( IsoTimeString, &UTC, 0, 3 );
                            //printf("Bracket: T = %ld %ld %ld    R = %g %g %g    Tperigee, Rperigee = %s %g\n", Ta, Tb, Tc, Ra, Rb, Rc, IsoTimeString, fabs(Rmin) );
                        }
                        // advance another 900 seconds.
                        Ta = Tb; Ra = Rb;
                        Tb = Tc; Rb = Rc;
                        Tc += 900.0; Rc = bInfo.func( (double)Tc, 0.0, (void *)afi );
                        if (Tc >= 86400) done = TRUE;
                    }

                    free( afi );





                    /*
                     * Open Mag Ephem file for writing and
                     * Write header to Mag Ephem file
                     */
                    /*
                     */
                    fp_MagEphem = fopen( OutFile, "w" );
printf("nPerigee = %d\n", nPerigee);
                    Lgm_WriteMagEphemHeader( fp_MagEphem, "FIX ME", 99999, "FIX ME", nPerigee, Perigee_UTC, Perigee_U, nApogee, Apogee_UTC, Apogee_U, MagEphemInfo );
                    printf("\t      Writing to file: %s\n", OutFile );
//exit(0);

                    if ( UseEop ) {
                        // Read in the EOP vals
                        Lgm_read_eop( e );
                    }




                    //for ( Seconds=0; Seconds<=86400; Seconds += 900 ) {
                    H5_nT = 0;
                    for ( Seconds=0; Seconds<=86400; Seconds += 60 ) {

                        Lgm_Make_UTC( Date, Seconds/3600.0, &UTC, c );
                        Lgm_DateTimeToString( IsoTimeString, &UTC, 0, 0 );
et = Lgm_TDBSecSinceJ2000( &UTC, c );
spkezp_c( RBSPA_ID,    et,   "J2000",  "NONE", EARTH_ID,  pos,  &lt );
printf("pos = %.8lf %.8lf %.8lf\n", pos[0], pos[1], pos[2]);
U.x = pos[0]/WGS84_A; U.y = pos[1]/WGS84_A; U.z = pos[2]/WGS84_A;


//// convert ISO time to DateTime
//IsoTimeStringToDateTime( IsoTimeString, &UTC, c );

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

Lgm_Convert_Coords( &U, &Rgsm, GEI2000_TO_GSM, c );




                            /*
                             * Compute L*s, Is, Bms, Footprints, etc...
                             * These quantities are stored in the MagEphemInfo Structure
                             */
                            printf("\n\n\t[ %s ]: %s  Bird: %s Rgsm: %g %g %g Re\n", ProgramName, IsoTimeString, Bird, Rgsm.x, Rgsm.y, Rgsm.z );
                            printf("\t--------------------------------------------------------------------------------------------------\n");
                            Lgm_ComputeLstarVersusPA( UTC.Date, UTC.Time, &Rgsm, nAlpha, Alpha, MagEphemInfo->LstarQuality, Colorize, MagEphemInfo );

                            Lgm_WriteMagEphemData( fp_MagEphem, IntModel, ExtModel, MagEphemInfo->LstarInfo->mInfo->fKp, MagEphemInfo->LstarInfo->mInfo->Dst, MagEphemInfo );


                            if ( nAlpha > 0 ){
                                WriteMagEphemInfoStruct( "test.dat", nAlpha, MagEphemInfo );
                            }




                            // Fill arrays for dumping out as HDF5 files
                            strcpy( H5_IsoTimes[H5_nT], IsoTimeString );
                            H5_Date[H5_nT]      = UTC.Date;
                            H5_Doy[H5_nT]       = UTC.Doy;
                            H5_UTC[H5_nT]       = UTC.Time;
                            H5_JD[H5_nT]        = UTC.JD;
                            H5_GpsTime[H5_nT]   = Lgm_UTC_to_GpsSeconds( &UTC, c );
                            H5_TiltAngle[H5_nT] = c->psi*DegPerRad;
                            H5_Rgsm[H5_nT]      = Rgsm;
                            H5_Kp[H5_nT]        = MagEphemInfo->LstarInfo->mInfo->fKp;
                            H5_Dst[H5_nT]       = MagEphemInfo->LstarInfo->mInfo->Dst;
                            MagEphemInfo->LstarInfo->mInfo->Bfield( &Rgsm, &H5_Bsc_gsm[H5_nT], MagEphemInfo->LstarInfo->mInfo );


                            for (i=0; i<nAlpha; i++){
                                H5_Lstar[H5_nT][i] = MagEphemInfo->Lstar[i];
                                H5_I[H5_nT][i]     = MagEphemInfo->I[i];
                                H5_Bm[H5_nT][i]    = MagEphemInfo->Bm[i];
                                if ( (MagEphemInfo->Bm[i]>0.0)&&(MagEphemInfo->I[i]>=0.0) ) {
                                    H5_K[H5_nT][i] = 3.16227766e-3*MagEphemInfo->I[i]*sqrt(MagEphemInfo->Bm[i]);
                                } else {
                                    H5_K[H5_nT][i] = LGM_FILL_VALUE;
                                }
                                if (MagEphemInfo->I[i]>=0.0) {
                                    H5_L[H5_nT][i] = LFromIBmM_McIlwain(MagEphemInfo->I[i], MagEphemInfo->Bm[i], MagEphemInfo->Mused );
                                } else {
                                    H5_L[H5_nT][i] = LGM_FILL_VALUE;
                                }
                            }

                            if ( (MagEphemInfo->FieldLineType == LGM_CLOSED) || (MagEphemInfo->FieldLineType == LGM_OPEN_N_LOBE) ) {

                                /*
                                 * Save northern Footpoint position in different coord systems.
                                 */
                                H5_Pfn_gsm[H5_nT] = MagEphemInfo->Ellipsoid_Footprint_Pn;

                                Lgm_Convert_Coords( &H5_Pfn_gsm[H5_nT], &W, GSM_TO_GEO, c );
                                H5_Pfn_geo[H5_nT] = W;

                                Lgm_WGS84_to_GEOD( &W, &GeodLat, &GeodLong, &GeodHeight );
                                H5_Pfn_geod[H5_nT].x = GeodLat; H5_Pfn_geod[H5_nT].y = GeodLong; H5_Pfn_geod[H5_nT].z = GeodHeight;

                                Lgm_Convert_Coords( &H5_Pfn_gsm[H5_nT], &W, GSM_TO_CDMAG, c );
                                Lgm_CDMAG_to_R_MLAT_MLON_MLT( &W, &R, &MLAT, &MLON, &MLT, c );
                                H5_Pfn_cdmag[H5_nT].x = MLAT; H5_Pfn_cdmag[H5_nT].y = MLON; H5_Pfn_cdmag[H5_nT].z = MLT;

                                Lgm_Convert_Coords( &H5_Pfn_gsm[H5_nT], &W, GSM_TO_EDMAG, c );
                                Lgm_EDMAG_to_R_MLAT_MLON_MLT( &W, &R, &MLAT, &MLON, &MLT, c );
                                H5_Pfn_edmag[H5_nT].x = MLAT; H5_Pfn_edmag[H5_nT].y = MLON; H5_Pfn_edmag[H5_nT].z = MLT;


                                /*
                                 * Save northern Footpoint B-field values in different coord systems.
                                 */
                                MagEphemInfo->LstarInfo->mInfo->Bfield( &MagEphemInfo->Ellipsoid_Footprint_Pn, &Bvec, MagEphemInfo->LstarInfo->mInfo );
                                Lgm_Convert_Coords( &Bvec, &Bvec2, GSM_TO_WGS84, c );
                                H5_Bfn_gsm[H5_nT] = Bvec;
                                H5_Bfn_geo[H5_nT] = Bvec2;


                            } else {

                                H5_Pfn_gsm[H5_nT].x = H5_Pfn_gsm[H5_nT].y = H5_Pfn_gsm[H5_nT].z = LGM_FILL_VALUE;
                                H5_Pfn_geo[H5_nT].x = H5_Pfn_geo[H5_nT].y = H5_Pfn_geo[H5_nT].z = LGM_FILL_VALUE;
                                H5_Pfn_geod[H5_nT].x = H5_Pfn_geod[H5_nT].y = H5_Pfn_geod[H5_nT].z = LGM_FILL_VALUE;
                                H5_Pfn_cdmag[H5_nT].x = H5_Pfn_cdmag[H5_nT].y = H5_Pfn_cdmag[H5_nT].z = LGM_FILL_VALUE;
                                H5_Pfn_edmag[H5_nT].x = H5_Pfn_edmag[H5_nT].y = H5_Pfn_edmag[H5_nT].z = LGM_FILL_VALUE;
                                H5_Bfn_gsm[H5_nT].x = H5_Bfn_gsm[H5_nT].y = H5_Bfn_gsm[H5_nT].z = LGM_FILL_VALUE;
                                H5_Bfn_geo[H5_nT].x = H5_Bfn_geo[H5_nT].y = H5_Bfn_geo[H5_nT].z = LGM_FILL_VALUE;

                            }

                            if ( (MagEphemInfo->FieldLineType == LGM_CLOSED) || (MagEphemInfo->FieldLineType == LGM_OPEN_S_LOBE) ) {

                                /*
                                 * Save Southern Footpoint position in different coord systems.
                                 */
                                H5_Pfs_gsm[H5_nT] = MagEphemInfo->Ellipsoid_Footprint_Pn;

                                Lgm_Convert_Coords( &H5_Pfs_gsm[H5_nT], &W, GSM_TO_GEO, c );
                                H5_Pfs_geo[H5_nT] = W;

                                Lgm_WGS84_to_GEOD( &W, &GeodLat, &GeodLong, &GeodHeight );
                                H5_Pfs_geod[H5_nT].x = GeodLat; H5_Pfs_geod[H5_nT].y = GeodLong; H5_Pfs_geod[H5_nT].z = GeodHeight;

                                Lgm_Convert_Coords( &H5_Pfs_gsm[H5_nT], &W, GSM_TO_CDMAG, c );
                                Lgm_CDMAG_to_R_MLAT_MLON_MLT( &W, &R, &MLAT, &MLON, &MLT, c );
                                H5_Pfs_cdmag[H5_nT].x = MLAT; H5_Pfs_cdmag[H5_nT].y = MLON; H5_Pfs_cdmag[H5_nT].z = MLT;

                                Lgm_Convert_Coords( &H5_Pfs_gsm[H5_nT], &W, GSM_TO_EDMAG, c );
                                Lgm_EDMAG_to_R_MLAT_MLON_MLT( &W, &R, &MLAT, &MLON, &MLT, c );
                                H5_Pfs_edmag[H5_nT].x = MLAT; H5_Pfs_edmag[H5_nT].y = MLON; H5_Pfs_edmag[H5_nT].z = MLT;


                                /*
                                 * Save Southern Footpoint B-field values in different coord systems.
                                 */
                                MagEphemInfo->LstarInfo->mInfo->Bfield( &MagEphemInfo->Ellipsoid_Footprint_Pn, &Bvec, MagEphemInfo->LstarInfo->mInfo );
                                Lgm_Convert_Coords( &Bvec, &Bvec2, GSM_TO_WGS84, c );
                                H5_Bfs_gsm[H5_nT] = Bvec;
                                H5_Bfs_geo[H5_nT] = Bvec2;


                            } else {

                                H5_Pfs_gsm[H5_nT].x = H5_Pfs_gsm[H5_nT].y = H5_Pfs_gsm[H5_nT].z = LGM_FILL_VALUE;
                                H5_Pfs_geo[H5_nT].x = H5_Pfs_geo[H5_nT].y = H5_Pfs_geo[H5_nT].z = LGM_FILL_VALUE;
                                H5_Pfs_geod[H5_nT].x = H5_Pfs_geod[H5_nT].y = H5_Pfs_geod[H5_nT].z = LGM_FILL_VALUE;
                                H5_Pfs_cdmag[H5_nT].x = H5_Pfs_cdmag[H5_nT].y = H5_Pfs_cdmag[H5_nT].z = LGM_FILL_VALUE;
                                H5_Pfs_edmag[H5_nT].x = H5_Pfs_edmag[H5_nT].y = H5_Pfs_edmag[H5_nT].z = LGM_FILL_VALUE;
                                H5_Bfs_gsm[H5_nT].x = H5_Bfs_gsm[H5_nT].y = H5_Bfs_gsm[H5_nT].z = LGM_FILL_VALUE;
                                H5_Bfs_geo[H5_nT].x = H5_Bfs_geo[H5_nT].y = H5_Bfs_geo[H5_nT].z = LGM_FILL_VALUE;

                            }





                            ++H5_nT;

                    }
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
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "UTC" );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    //Lgm_WriteStringAttr( DataSet, "VALIDMIN",   "0.0" );
                    //Lgm_WriteStringAttr( DataSet, "VALIDMAX",   "1E30" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );
                    status  = H5Tclose( atype );

                    // Create Date Dataset
                    Dims[0] = H5_nT;
                    space   = H5Screate_simple( 1, Dims, NULL ); // rank 1
                    DataSet = H5Dcreate( file, "Date", H5T_NATIVE_LONG, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "YYYYMMDD" );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    //Lgm_WriteStringAttr( DataSet, "VALIDMIN",   "0.0" );
                    //Lgm_WriteStringAttr( DataSet, "VALIDMAX",   "1E30" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );

                    // Create Doy Dataset
                    Dims[0] = H5_nT;
                    space   = H5Screate_simple( 1, Dims, NULL ); // rank 1
                    DataSet = H5Dcreate( file, "Doy", H5T_NATIVE_INT, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "DDD" );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    //Lgm_WriteStringAttr( DataSet, "VALIDMIN",   "0.0" );
                    //Lgm_WriteStringAttr( DataSet, "VALIDMAX",   "1E30" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );

                    // Create UTC (hours) Dataset
                    Dims[0] = H5_nT;
                    space   = H5Screate_simple( 1, Dims, NULL ); // rank 1
                    DataSet = H5Dcreate( file, "UTC", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "Hours" );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    //Lgm_WriteStringAttr( DataSet, "VALIDMIN",   "0.0" );
                    //Lgm_WriteStringAttr( DataSet, "VALIDMAX",   "1E30" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );

                    // Create JD Dataset
                    Dims[0] = H5_nT;
                    space   = H5Screate_simple( 1, Dims, NULL ); // rank 1
                    DataSet = H5Dcreate( file, "JulianDate", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "Days" );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    //Lgm_WriteStringAttr( DataSet, "VALIDMIN",   "0.0" );
                    //Lgm_WriteStringAttr( DataSet, "VALIDMAX",   "1E30" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );

                    // Create GpsTime Dataset
                    Dims[0] = H5_nT;
                    space   = H5Screate_simple( 1, Dims, NULL ); // rank 1
                    DataSet = H5Dcreate( file, "GpsTime", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "Seconds" );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    //Lgm_WriteStringAttr( DataSet, "VALIDMIN",   "0.0" );
                    //Lgm_WriteStringAttr( DataSet, "VALIDMAX",   "1E30" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );

                    // Create TiltAngle Dataset
                    Dims[0] = H5_nT;
                    space   = H5Screate_simple( 1, Dims, NULL ); // rank 1
                    DataSet = H5Dcreate( file, "TiltAngle", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "Degrees" );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    //Lgm_WriteStringAttr( DataSet, "VALIDMIN",   "-40.0" );
                    //Lgm_WriteStringAttr( DataSet, "VALIDMAX",   "40.0" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );

                    // Create Rgeo Dataset
                    Dims[0] = H5_nT; Dims[1] = 3;
                    space   = H5Screate_simple( 2, Dims, NULL ); // rank 2
                    DataSet = H5Dcreate( file, "Rgeo", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "Re" );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    //Lgm_WriteStringAttr( DataSet, "VALIDMIN",   "0.0" );
                    //Lgm_WriteStringAttr( DataSet, "VALIDMAX",   "1E30" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );

                    // Create Rgeod Dataset
                    Dims[0] = H5_nT; Dims[1] = 3;
                    space   = H5Screate_simple( 2, Dims, NULL ); // rank 2
                    DataSet = H5Dcreate( file, "Rgeod", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "Deg., Deg., km" );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    //Lgm_WriteStringAttr( DataSet, "VALIDMIN",   "0.0" );
                    //Lgm_WriteStringAttr( DataSet, "VALIDMAX",   "1E30" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );

                    // Create Rgsm Dataset
                    Dims[0] = H5_nT; Dims[1] = 3;
                    space   = H5Screate_simple( 2, Dims, NULL ); // rank 2
                    DataSet = H5Dcreate( file, "Rgsm", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "Re" );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    //Lgm_WriteStringAttr( DataSet, "VALIDMIN",   "0.0" );
                    //Lgm_WriteStringAttr( DataSet, "VALIDMAX",   "1E30" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );

                    // Create Rsm Dataset
                    Dims[0] = H5_nT; Dims[1] = 3;
                    space   = H5Screate_simple( 2, Dims, NULL ); // rank 2
                    DataSet = H5Dcreate( file, "Rsm", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "Re" );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    //Lgm_WriteStringAttr( DataSet, "VALIDMIN",   "0.0" );
                    //Lgm_WriteStringAttr( DataSet, "VALIDMAX",   "1E30" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );

                    // Create Rgei Dataset
                    Dims[0] = H5_nT; Dims[1] = 3;
                    space   = H5Screate_simple( 2, Dims, NULL ); // rank 2
                    DataSet = H5Dcreate( file, "Rgei", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "Re" );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    //Lgm_WriteStringAttr( DataSet, "VALIDMIN",   "0.0" );
                    //Lgm_WriteStringAttr( DataSet, "VALIDMAX",   "1E30" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );

                    // Create Rgse Dataset
                    Dims[0] = H5_nT; Dims[1] = 3;
                    space   = H5Screate_simple( 2, Dims, NULL ); // rank 2
                    DataSet = H5Dcreate( file, "Rgse", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "Re" );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    //Lgm_WriteStringAttr( DataSet, "VALIDMIN",   "0.0" );
                    //Lgm_WriteStringAttr( DataSet, "VALIDMAX",   "1E30" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );

                    // Create CDMAG Dataset
                    Dims[0] = H5_nT; Dims[1] = 4;
                    space   = H5Screate_simple( 2, Dims, NULL ); // rank 2
                    DataSet = H5Dcreate( file, "CDMAG", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "Deg., Deg., Hours, Re" );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    //Lgm_WriteStringAttr( DataSet, "VALIDMIN",   "0.0" );
                    //Lgm_WriteStringAttr( DataSet, "VALIDMAX",   "1E30" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );

                    // Create EDMAG Dataset
                    Dims[0] = H5_nT; Dims[1] = 4;
                    space   = H5Screate_simple( 2, Dims, NULL ); // rank 2
                    DataSet = H5Dcreate( file, "EDMAG", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "Deg., Deg., Hours, Re" );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    //Lgm_WriteStringAttr( DataSet, "VALIDMIN",   "0.0" );
                    //Lgm_WriteStringAttr( DataSet, "VALIDMAX",   "1E30" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );

                    // Create Bsc_gsm Dataset
                    Dims[0] = H5_nT; Dims[1] = 4;
                    space   = H5Screate_simple( 2, Dims, NULL ); // rank 2
                    DataSet = H5Dcreate( file, "Bsc_gsm", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "nT" );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    //Lgm_WriteStringAttr( DataSet, "VALIDMIN",   "0.0" );
                    //Lgm_WriteStringAttr( DataSet, "VALIDMAX",   "1E30" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );

                    // Create Pfn_geo Dataset
                    Dims[0] = H5_nT; Dims[1] = 3;
                    space   = H5Screate_simple( 2, Dims, NULL ); // rank 2
                    DataSet = H5Dcreate( file, "Pfn_geo", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "Re" );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    //Lgm_WriteStringAttr( DataSet, "VALIDMIN",   "0.0" );
                    //Lgm_WriteStringAttr( DataSet, "VALIDMAX",   "1E30" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );

                    // Create Pfn_gsm Dataset
                    Dims[0] = H5_nT; Dims[1] = 3;
                    space   = H5Screate_simple( 2, Dims, NULL ); // rank 2
                    DataSet = H5Dcreate( file, "Pfn_gsm", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "Re" );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    //Lgm_WriteStringAttr( DataSet, "VALIDMIN",   "0.0" );
                    //Lgm_WriteStringAttr( DataSet, "VALIDMAX",   "1E30" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );

                    // Create Pfn_geod Dataset
                    Dims[0] = H5_nT; Dims[1] = 3;
                    space   = H5Screate_simple( 2, Dims, NULL ); // rank 2
                    DataSet = H5Dcreate( file, "Pfn_geod", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "Deg., Deg., km" );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    //Lgm_WriteStringAttr( DataSet, "VALIDMIN",   "0.0" );
                    //Lgm_WriteStringAttr( DataSet, "VALIDMAX",   "1E30" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );

                    // Create Pfn_cdmag Dataset
                    Dims[0] = H5_nT; Dims[1] = 3;
                    space   = H5Screate_simple( 2, Dims, NULL ); // rank 2
                    DataSet = H5Dcreate( file, "Pfn_cdmag", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "Deg., Deg., Hours" );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    //Lgm_WriteStringAttr( DataSet, "VALIDMIN",   "0.0" );
                    //Lgm_WriteStringAttr( DataSet, "VALIDMAX",   "1E30" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );

                    // Create Pfn_edmag Dataset
                    Dims[0] = H5_nT; Dims[1] = 3;
                    space   = H5Screate_simple( 2, Dims, NULL ); // rank 2
                    DataSet = H5Dcreate( file, "Pfn_edmag", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "Deg., Deg., Hours" );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    //Lgm_WriteStringAttr( DataSet, "VALIDMIN",   "0.0" );
                    //Lgm_WriteStringAttr( DataSet, "VALIDMAX",   "1E30" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );

                    // Create Pfs_geo Dataset
                    Dims[0] = H5_nT; Dims[1] = 3;
                    space   = H5Screate_simple( 2, Dims, NULL ); // rank 2
                    DataSet = H5Dcreate( file, "Pfs_geo", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "Re" );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    //Lgm_WriteStringAttr( DataSet, "VALIDMIN",   "0.0" );
                    //Lgm_WriteStringAttr( DataSet, "VALIDMAX",   "1E30" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );

                    // Create Pfs_gsm Dataset
                    Dims[0] = H5_nT; Dims[1] = 3;
                    space   = H5Screate_simple( 2, Dims, NULL ); // rank 2
                    DataSet = H5Dcreate( file, "Pfs_gsm", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "Re" );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    //Lgm_WriteStringAttr( DataSet, "VALIDMIN",   "0.0" );
                    //Lgm_WriteStringAttr( DataSet, "VALIDMAX",   "1E30" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );

                    // Create Pfs_geod Dataset
                    Dims[0] = H5_nT; Dims[1] = 3;
                    space   = H5Screate_simple( 2, Dims, NULL ); // rank 2
                    DataSet = H5Dcreate( file, "Pfs_geod", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "Deg., Deg., km" );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    //Lgm_WriteStringAttr( DataSet, "VALIDMIN",   "0.0" );
                    //Lgm_WriteStringAttr( DataSet, "VALIDMAX",   "1E30" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );

                    // Create Pfs_cdmag Dataset
                    Dims[0] = H5_nT; Dims[1] = 3;
                    space   = H5Screate_simple( 2, Dims, NULL ); // rank 2
                    DataSet = H5Dcreate( file, "Pfs_cdmag", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "Deg., Deg., Hours" );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );

                    // Create Pfs_edmag Dataset
                    Dims[0] = H5_nT; Dims[1] = 3;
                    space   = H5Screate_simple( 2, Dims, NULL ); // rank 2
                    DataSet = H5Dcreate( file, "Pfs_edmag", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "Deg., Deg., Hours" );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );




                    // Create Alpha Dataset
                    Dims[0] = H5_nAlpha;
                    space   = H5Screate_simple( 1, Dims, NULL ); // rank 1
                    DataSet = H5Dcreate( file, "Alpha", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "Degrees" );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );

                    // Create Kp Dataset
                    Dims[0] = H5_nT;
                    space   = H5Screate_simple( 1, Dims, NULL ); // rank 1
                    DataSet = H5Dcreate( file, "Kp", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "Dimensionless" );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );

                    // Create Dst Dataset
                    Dims[0] = H5_nT;
                    space   = H5Screate_simple( 1, Dims, NULL ); // rank 1
                    DataSet = H5Dcreate( file, "Dst", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "nT" );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );



                    // Create Lstar Dataset
                    Dims[0] = H5_nT; Dims[1] = H5_nAlpha;
                    space   = H5Screate_simple( 2, Dims, NULL ); // rank 2
                    DataSet = H5Dcreate( file, "Lstar", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "DEPEND_1",   "Alpha" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "Dimensionless" );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    //Lgm_WriteStringAttr( DataSet, "VALIDMIN",   "1.0" );
                    //Lgm_WriteStringAttr( DataSet, "VALIDMAX",   "20.0" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );

                    // Create K Dataset
                    Dims[0] = H5_nT; Dims[1] = H5_nAlpha;
                    space   = H5Screate_simple( 2, Dims, NULL ); // rank 2
                    DataSet = H5Dcreate( file, "K", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "DEPEND_1",   "Alpha" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "Re G^.5" );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "log" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    Lgm_WriteStringAttr( DataSet, "VALIDMIN",   "0.0" );
                    //Lgm_WriteStringAttr( DataSet, "VALIDMAX",   "1E30" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );

                    // Create L Dataset
                    Dims[0] = H5_nT; Dims[1] = H5_nAlpha;
                    space   = H5Screate_simple( 2, Dims, NULL ); // rank 2
                    DataSet = H5Dcreate( file, "L", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "DEPEND_1",   "Alpha" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "Dimensionless" );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    //Lgm_WriteStringAttr( DataSet, "VALIDMIN",   "1.0" );
                    //Lgm_WriteStringAttr( DataSet, "VALIDMAX",   "20.0" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );

                    // Create I Dataset
                    Dims[0] = H5_nT; Dims[1] = H5_nAlpha;
                    space   = H5Screate_simple( 2, Dims, NULL ); // rank 2
                    DataSet = H5Dcreate( file, "I", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "DEPEND_1",   "Alpha" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "Re" );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    //Lgm_WriteStringAttr( DataSet, "VALIDMIN",   "0.0" );
                    //Lgm_WriteStringAttr( DataSet, "VALIDMAX",   "200.0" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );

                    // Create Bm Dataset
                    Dims[0] = H5_nT; Dims[1] = H5_nAlpha;
                    space   = H5Screate_simple( 2, Dims, NULL ); // rank 2
                    DataSet = H5Dcreate( file, "Bm", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "DEPEND_1",   "Alpha" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "nT" );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "log" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    //Lgm_WriteStringAttr( DataSet, "VALIDMIN",   "0.0" );
                    //Lgm_WriteStringAttr( DataSet, "VALIDMAX",   "1E30" );
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

                    // Write Date
                    SlabSize[0] = H5_nT;
                    Offset[0]   = 0;
                    DataSet  = H5Dopen( file, "Date", H5P_DEFAULT );
                    space    = H5Dget_space( DataSet );
                    status   = H5Sselect_hyperslab( space, H5S_SELECT_SET, Offset, NULL, SlabSize, NULL );
                    MemSpace = H5Screate_simple( 1, SlabSize, NULL );
                    status   = H5Dwrite( DataSet, H5T_NATIVE_LONG, MemSpace, space, H5P_DEFAULT, &H5_Date[0] );
                    status   = H5Sclose( MemSpace );
                    status   = H5Sclose( space );
                    status   = H5Dclose( DataSet );

                    // Write Doy
                    SlabSize[0] = H5_nT;
                    Offset[0]   = 0;
                    DataSet  = H5Dopen( file, "Doy", H5P_DEFAULT );
                    space    = H5Dget_space( DataSet );
                    status   = H5Sselect_hyperslab( space, H5S_SELECT_SET, Offset, NULL, SlabSize, NULL );
                    MemSpace = H5Screate_simple( 1, SlabSize, NULL );
                    status   = H5Dwrite( DataSet, H5T_NATIVE_INT, MemSpace, space, H5P_DEFAULT, &H5_Doy[0] );
                    status   = H5Sclose( MemSpace );
                    status   = H5Sclose( space );
                    status   = H5Dclose( DataSet );

                    // Write UTC (hours)
                    SlabSize[0] = H5_nT;
                    Offset[0]   = 0;
                    DataSet  = H5Dopen( file, "UTC", H5P_DEFAULT );
                    space    = H5Dget_space( DataSet );
                    status   = H5Sselect_hyperslab( space, H5S_SELECT_SET, Offset, NULL, SlabSize, NULL );
                    MemSpace = H5Screate_simple( 1, SlabSize, NULL );
                    status   = H5Dwrite( DataSet, H5T_NATIVE_DOUBLE, MemSpace, space, H5P_DEFAULT, &H5_UTC[0] );
                    status   = H5Sclose( MemSpace );
                    status   = H5Sclose( space );
                    status   = H5Dclose( DataSet );

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

                    // Write GpsTime
                    SlabSize[0] = H5_nT;
                    Offset[0]   = 0;
                    DataSet  = H5Dopen( file, "GpsTime", H5P_DEFAULT );
                    space    = H5Dget_space( DataSet );
                    status   = H5Sselect_hyperslab( space, H5S_SELECT_SET, Offset, NULL, SlabSize, NULL );
                    MemSpace = H5Screate_simple( 1, SlabSize, NULL );
                    status   = H5Dwrite( DataSet, H5T_NATIVE_DOUBLE, MemSpace, space, H5P_DEFAULT, &H5_GpsTime[0] );
                    status   = H5Sclose( MemSpace );
                    status   = H5Sclose( space );
                    status   = H5Dclose( DataSet );

                    // Write TiltAngle
                    SlabSize[0] = H5_nT;
                    Offset[0]   = 0;
                    DataSet  = H5Dopen( file, "TiltAngle", H5P_DEFAULT );
                    space    = H5Dget_space( DataSet );
                    status   = H5Sselect_hyperslab( space, H5S_SELECT_SET, Offset, NULL, SlabSize, NULL );
                    MemSpace = H5Screate_simple( 1, SlabSize, NULL );
                    status   = H5Dwrite( DataSet, H5T_NATIVE_DOUBLE, MemSpace, space, H5P_DEFAULT, &H5_TiltAngle[0] );
                    status   = H5Sclose( MemSpace );
                    status   = H5Sclose( space );
                    status   = H5Dclose( DataSet );



                    // Write Rgeo
                    SlabSize[0] = 1; SlabSize[1] = 3;
                    DataSet  = H5Dopen( file, "Rgeo", H5P_DEFAULT );
                    space    = H5Dget_space( DataSet );
                    MemSpace = H5Screate_simple( 2, SlabSize, NULL );
                    for ( iT=0; iT<H5_nT; iT++ ){
                            Lgm_Set_Coord_Transforms( H5_Date[iT], H5_UTC[iT], c );
                            Lgm_Convert_Coords( &H5_Rgsm[iT], &W, GSM_TO_GEO, c ); // Convert to cartesian GEO coords.
                            U_ARR[0] = W.x; U_ARR[1] = W.y; U_ARR[2] = W.z;
                            Offset[0]   = iT; Offset[1] = 0;
                            status = H5Sselect_hyperslab( space, H5S_SELECT_SET, Offset, NULL, SlabSize, NULL );
                            status = H5Dwrite( DataSet, H5T_NATIVE_DOUBLE, MemSpace, space, H5P_DEFAULT, &U_ARR[0] );
                    }
                    status = H5Sclose( MemSpace );
                    status = H5Dclose( DataSet );
                    status = H5Sclose( space );

                    // Write Rgeod
                    SlabSize[0] = 1; SlabSize[1] = 3;
                    DataSet  = H5Dopen( file, "Rgeod", H5P_DEFAULT );
                    space    = H5Dget_space( DataSet );
                    MemSpace = H5Screate_simple( 2, SlabSize, NULL );
                    for ( iT=0; iT<H5_nT; iT++ ){
                            Lgm_Set_Coord_Transforms( H5_Date[iT], H5_UTC[iT], c );
                            Lgm_Convert_Coords( &H5_Rgsm[iT], &W, GSM_TO_GEO, c ); // Convert to cartesian GEO coords.
                            Lgm_WGS84_to_GEOD( &W, &GeodLat, &GeodLong, &GeodHeight );
                            U_ARR[0] = GeodLat; U_ARR[1] = GeodLong; U_ARR[2] = GeodHeight;
                            Offset[0]   = iT; Offset[1] = 0;
                            status = H5Sselect_hyperslab( space, H5S_SELECT_SET, Offset, NULL, SlabSize, NULL );
                            status = H5Dwrite( DataSet, H5T_NATIVE_DOUBLE, MemSpace, space, H5P_DEFAULT, &U_ARR[0] );
                    }
                    status = H5Sclose( MemSpace );
                    status = H5Dclose( DataSet );
                    status = H5Sclose( space );

                    // Write Rgsm
                    SlabSize[0] = 1; SlabSize[1] = 3;
                    DataSet  = H5Dopen( file, "Rgsm", H5P_DEFAULT );
                    space    = H5Dget_space( DataSet );
                    MemSpace = H5Screate_simple( 2, SlabSize, NULL );
                    for ( iT=0; iT<H5_nT; iT++ ){
                            U_ARR[0] = H5_Rgsm[iT].x; U_ARR[1] = H5_Rgsm[iT].y; U_ARR[2] = H5_Rgsm[iT].z;
                            Offset[0]   = iT; Offset[1] = 0;
                            status = H5Sselect_hyperslab( space, H5S_SELECT_SET, Offset, NULL, SlabSize, NULL );
                            status = H5Dwrite( DataSet, H5T_NATIVE_DOUBLE, MemSpace, space, H5P_DEFAULT, &U_ARR[0] );
                    }
                    status = H5Sclose( MemSpace );
                    status = H5Dclose( DataSet );
                    status = H5Sclose( space );

                    // Write Rsm
                    SlabSize[0] = 1; SlabSize[1] = 3;
                    DataSet  = H5Dopen( file, "Rsm", H5P_DEFAULT );
                    space    = H5Dget_space( DataSet );
                    MemSpace = H5Screate_simple( 2, SlabSize, NULL );
                    for ( iT=0; iT<H5_nT; iT++ ){
                            Lgm_Set_Coord_Transforms( H5_Date[iT], H5_UTC[iT], c );
                            Lgm_Convert_Coords( &H5_Rgsm[iT], &W, GSM_TO_SM, c ); // Convert to cartesian SM coords.
                            U_ARR[0] = W.x; U_ARR[1] = W.y; U_ARR[2] = W.z;
                            Offset[0]   = iT; Offset[1] = 0;
                            status = H5Sselect_hyperslab( space, H5S_SELECT_SET, Offset, NULL, SlabSize, NULL );
                            status = H5Dwrite( DataSet, H5T_NATIVE_DOUBLE, MemSpace, space, H5P_DEFAULT, &U_ARR[0] );
                    }
                    status = H5Sclose( MemSpace );
                    status = H5Dclose( DataSet );
                    status = H5Sclose( space );

                    // Write Rgei
                    SlabSize[0] = 1; SlabSize[1] = 3;
                    DataSet  = H5Dopen( file, "Rgei", H5P_DEFAULT );
                    space    = H5Dget_space( DataSet );
                    MemSpace = H5Screate_simple( 2, SlabSize, NULL );
                    for ( iT=0; iT<H5_nT; iT++ ){
                            Lgm_Set_Coord_Transforms( H5_Date[iT], H5_UTC[iT], c );
                            Lgm_Convert_Coords( &H5_Rgsm[iT], &W, GSM_TO_GEI2000, c ); // Convert to cartesian GEI2000 coords.
                            U_ARR[0] = W.x; U_ARR[1] = W.y; U_ARR[2] = W.z;
                            Offset[0]   = iT; Offset[1] = 0;
                            status = H5Sselect_hyperslab( space, H5S_SELECT_SET, Offset, NULL, SlabSize, NULL );
                            status = H5Dwrite( DataSet, H5T_NATIVE_DOUBLE, MemSpace, space, H5P_DEFAULT, &U_ARR[0] );
                    }
                    status = H5Sclose( MemSpace );
                    status = H5Dclose( DataSet );
                    status = H5Sclose( space );

                    // Write Rgse
                    SlabSize[0] = 1; SlabSize[1] = 3;
                    DataSet  = H5Dopen( file, "Rgse", H5P_DEFAULT );
                    space    = H5Dget_space( DataSet );
                    MemSpace = H5Screate_simple( 2, SlabSize, NULL );
                    for ( iT=0; iT<H5_nT; iT++ ){
                            Lgm_Set_Coord_Transforms( H5_Date[iT], H5_UTC[iT], c );
                            Lgm_Convert_Coords( &H5_Rgsm[iT], &W, GSM_TO_GSE, c ); // Convert to cartesian GSE coords.
                            U_ARR[0] = W.x; U_ARR[1] = W.y; U_ARR[2] = W.z;
                            Offset[0]   = iT; Offset[1] = 0;
                            status = H5Sselect_hyperslab( space, H5S_SELECT_SET, Offset, NULL, SlabSize, NULL );
                            status = H5Dwrite( DataSet, H5T_NATIVE_DOUBLE, MemSpace, space, H5P_DEFAULT, &U_ARR[0] );
                    }
                    status = H5Sclose( MemSpace );
                    status = H5Dclose( DataSet );
                    status = H5Sclose( space );

                    // Write CDMAG
                    SlabSize[0] = 1; SlabSize[1] = 4;
                    DataSet  = H5Dopen( file, "CDMAG", H5P_DEFAULT );
                    space    = H5Dget_space( DataSet );
                    MemSpace = H5Screate_simple( 2, SlabSize, NULL );
                    for ( iT=0; iT<H5_nT; iT++ ){
                            Lgm_Set_Coord_Transforms( H5_Date[iT], H5_UTC[iT], c );
                            Lgm_Convert_Coords( &H5_Rgsm[iT], &W, GSM_TO_CDMAG, c ); // Convert to cartesian CDMAG coords.
                            Lgm_CDMAG_to_R_MLAT_MLON_MLT( &W, &R, &MLAT, &MLON, &MLT, c );
                            U_ARR[0] = MLAT; U_ARR[1] = MLON; U_ARR[2] = MLT; U_ARR[3] = R;
                            Offset[0]   = iT; Offset[1] = 0;
                            status = H5Sselect_hyperslab( space, H5S_SELECT_SET, Offset, NULL, SlabSize, NULL );
                            status = H5Dwrite( DataSet, H5T_NATIVE_DOUBLE, MemSpace, space, H5P_DEFAULT, &U_ARR[0] );
                    }
                    status = H5Sclose( MemSpace );
                    status = H5Dclose( DataSet );
                    status = H5Sclose( space );

                    // Write EDMAG
                    SlabSize[0] = 1; SlabSize[1] = 4;
                    DataSet  = H5Dopen( file, "EDMAG", H5P_DEFAULT );
                    space    = H5Dget_space( DataSet );
                    MemSpace = H5Screate_simple( 2, SlabSize, NULL );
                    for ( iT=0; iT<H5_nT; iT++ ){
                            Lgm_Set_Coord_Transforms( H5_Date[iT], H5_UTC[iT], c );
                            Lgm_Convert_Coords( &H5_Rgsm[iT], &W, GSM_TO_EDMAG, c ); // Convert to cartesian CDMAG coords.
                            Lgm_EDMAG_to_R_MLAT_MLON_MLT( &W, &R, &MLAT, &MLON, &MLT, c );
                            U_ARR[0] = MLAT; U_ARR[1] = MLON; U_ARR[2] = MLT; U_ARR[3] = R;
                            Offset[0]   = iT; Offset[1] = 0;
                            status = H5Sselect_hyperslab( space, H5S_SELECT_SET, Offset, NULL, SlabSize, NULL );
                            status = H5Dwrite( DataSet, H5T_NATIVE_DOUBLE, MemSpace, space, H5P_DEFAULT, &U_ARR[0] );
                    }
                    status = H5Sclose( MemSpace );
                    status = H5Dclose( DataSet );
                    status = H5Sclose( space );

                    // Write Bsc_gsm
                    SlabSize[0] = 1; SlabSize[1] = 4;
                    DataSet  = H5Dopen( file, "Bsc_gsm", H5P_DEFAULT );
                    space    = H5Dget_space( DataSet );
                    MemSpace = H5Screate_simple( 2, SlabSize, NULL );
                    for ( iT=0; iT<H5_nT; iT++ ){
                            U_ARR[0] = H5_Bsc_gsm[iT].x; U_ARR[1] = H5_Bsc_gsm[iT].y; U_ARR[2] = H5_Bsc_gsm[iT].z; U_ARR[3] = Lgm_Magnitude( &H5_Bsc_gsm[iT] );
                            Offset[0]   = iT; Offset[1] = 0;
                            status = H5Sselect_hyperslab( space, H5S_SELECT_SET, Offset, NULL, SlabSize, NULL );
                            status = H5Dwrite( DataSet, H5T_NATIVE_DOUBLE, MemSpace, space, H5P_DEFAULT, &U_ARR[0] );
                    }
                    status = H5Sclose( MemSpace );
                    status = H5Dclose( DataSet );
                    status = H5Sclose( space );

                    // Write Pfn_geo
                    SlabSize[0] = 1; SlabSize[1] = 3;
                    DataSet  = H5Dopen( file, "Pfn_geo", H5P_DEFAULT );
                    space    = H5Dget_space( DataSet );
                    MemSpace = H5Screate_simple( 2, SlabSize, NULL );
                    for ( iT=0; iT<H5_nT; iT++ ){
                            U_ARR[0] = H5_Pfn_geo[iT].x; U_ARR[1] = H5_Pfn_geo[iT].y; U_ARR[2] = H5_Pfn_geo[iT].z;
                            Offset[0]   = iT; Offset[1] = 0;
                            status = H5Sselect_hyperslab( space, H5S_SELECT_SET, Offset, NULL, SlabSize, NULL );
                            status = H5Dwrite( DataSet, H5T_NATIVE_DOUBLE, MemSpace, space, H5P_DEFAULT, &U_ARR[0] );
                    }
                    status = H5Sclose( MemSpace );
                    status = H5Dclose( DataSet );
                    status = H5Sclose( space );

                    // Write Pfn_gsm
                    SlabSize[0] = 1; SlabSize[1] = 3;
                    DataSet  = H5Dopen( file, "Pfn_gsm", H5P_DEFAULT );
                    space    = H5Dget_space( DataSet );
                    MemSpace = H5Screate_simple( 2, SlabSize, NULL );
                    for ( iT=0; iT<H5_nT; iT++ ){
                            U_ARR[0] = H5_Pfn_gsm[iT].x; U_ARR[1] = H5_Pfn_gsm[iT].y; U_ARR[2] = H5_Pfn_gsm[iT].z;
                            Offset[0]   = iT; Offset[1] = 0;
                            status = H5Sselect_hyperslab( space, H5S_SELECT_SET, Offset, NULL, SlabSize, NULL );
                            status = H5Dwrite( DataSet, H5T_NATIVE_DOUBLE, MemSpace, space, H5P_DEFAULT, &U_ARR[0] );
                    }
                    status = H5Sclose( MemSpace );
                    status = H5Dclose( DataSet );
                    status = H5Sclose( space );

                    // Write Pfn_geod
                    SlabSize[0] = 1; SlabSize[1] = 3;
                    DataSet  = H5Dopen( file, "Pfn_geod", H5P_DEFAULT );
                    space    = H5Dget_space( DataSet );
                    MemSpace = H5Screate_simple( 2, SlabSize, NULL );
                    for ( iT=0; iT<H5_nT; iT++ ){
                            U_ARR[0] = H5_Pfn_geod[iT].x; U_ARR[1] = H5_Pfn_geod[iT].y; U_ARR[2] = H5_Pfn_geod[iT].z;
                            Offset[0]   = iT; Offset[1] = 0;
                            status = H5Sselect_hyperslab( space, H5S_SELECT_SET, Offset, NULL, SlabSize, NULL );
                            status = H5Dwrite( DataSet, H5T_NATIVE_DOUBLE, MemSpace, space, H5P_DEFAULT, &U_ARR[0] );
                    }
                    status = H5Sclose( MemSpace );
                    status = H5Dclose( DataSet );
                    status = H5Sclose( space );

                    // Write Pfn_cdmag
                    SlabSize[0] = 1; SlabSize[1] = 3;
                    DataSet  = H5Dopen( file, "Pfn_cdmag", H5P_DEFAULT );
                    space    = H5Dget_space( DataSet );
                    MemSpace = H5Screate_simple( 2, SlabSize, NULL );
                    for ( iT=0; iT<H5_nT; iT++ ){
                            U_ARR[0] = H5_Pfn_cdmag[iT].x; U_ARR[1] = H5_Pfn_cdmag[iT].y; U_ARR[2] = H5_Pfn_cdmag[iT].z;
                            Offset[0]   = iT; Offset[1] = 0;
                            status = H5Sselect_hyperslab( space, H5S_SELECT_SET, Offset, NULL, SlabSize, NULL );
                            status = H5Dwrite( DataSet, H5T_NATIVE_DOUBLE, MemSpace, space, H5P_DEFAULT, &U_ARR[0] );
                    }
                    status = H5Sclose( MemSpace );
                    status = H5Dclose( DataSet );
                    status = H5Sclose( space );

                    // Write Pfn_edmag
                    SlabSize[0] = 1; SlabSize[1] = 3;
                    DataSet  = H5Dopen( file, "Pfn_edmag", H5P_DEFAULT );
                    space    = H5Dget_space( DataSet );
                    MemSpace = H5Screate_simple( 2, SlabSize, NULL );
                    for ( iT=0; iT<H5_nT; iT++ ){
                            U_ARR[0] = H5_Pfn_edmag[iT].x; U_ARR[1] = H5_Pfn_edmag[iT].y; U_ARR[2] = H5_Pfn_edmag[iT].z;
                            Offset[0]   = iT; Offset[1] = 0;
                            status = H5Sselect_hyperslab( space, H5S_SELECT_SET, Offset, NULL, SlabSize, NULL );
                            status = H5Dwrite( DataSet, H5T_NATIVE_DOUBLE, MemSpace, space, H5P_DEFAULT, &U_ARR[0] );
                    }
                    status = H5Sclose( MemSpace );
                    status = H5Dclose( DataSet );
                    status = H5Sclose( space );

                    // Write Pfs_geo
                    SlabSize[0] = 1; SlabSize[1] = 3;
                    DataSet  = H5Dopen( file, "Pfs_geo", H5P_DEFAULT );
                    space    = H5Dget_space( DataSet );
                    MemSpace = H5Screate_simple( 2, SlabSize, NULL );
                    for ( iT=0; iT<H5_nT; iT++ ){
                            U_ARR[0] = H5_Pfs_geo[iT].x; U_ARR[1] = H5_Pfs_geo[iT].y; U_ARR[2] = H5_Pfs_geo[iT].z;
                            Offset[0]   = iT; Offset[1] = 0;
                            status = H5Sselect_hyperslab( space, H5S_SELECT_SET, Offset, NULL, SlabSize, NULL );
                            status = H5Dwrite( DataSet, H5T_NATIVE_DOUBLE, MemSpace, space, H5P_DEFAULT, &U_ARR[0] );
                    }
                    status = H5Sclose( MemSpace );
                    status = H5Dclose( DataSet );
                    status = H5Sclose( space );

                    // Write Pfs_gsm
                    SlabSize[0] = 1; SlabSize[1] = 3;
                    DataSet  = H5Dopen( file, "Pfs_gsm", H5P_DEFAULT );
                    space    = H5Dget_space( DataSet );
                    MemSpace = H5Screate_simple( 2, SlabSize, NULL );
                    for ( iT=0; iT<H5_nT; iT++ ){
                            U_ARR[0] = H5_Pfs_gsm[iT].x; U_ARR[1] = H5_Pfs_gsm[iT].y; U_ARR[2] = H5_Pfs_gsm[iT].z;
                            Offset[0]   = iT; Offset[1] = 0;
                            status = H5Sselect_hyperslab( space, H5S_SELECT_SET, Offset, NULL, SlabSize, NULL );
                            status = H5Dwrite( DataSet, H5T_NATIVE_DOUBLE, MemSpace, space, H5P_DEFAULT, &U_ARR[0] );
                    }
                    status = H5Sclose( MemSpace );
                    status = H5Dclose( DataSet );
                    status = H5Sclose( space );

                    // Write Pfs_geod
                    SlabSize[0] = 1; SlabSize[1] = 3;
                    DataSet  = H5Dopen( file, "Pfs_geod", H5P_DEFAULT );
                    space    = H5Dget_space( DataSet );
                    MemSpace = H5Screate_simple( 2, SlabSize, NULL );
                    for ( iT=0; iT<H5_nT; iT++ ){
                            U_ARR[0] = H5_Pfs_geod[iT].x; U_ARR[1] = H5_Pfs_geod[iT].y; U_ARR[2] = H5_Pfs_geod[iT].z;
                            Offset[0]   = iT; Offset[1] = 0;
                            status = H5Sselect_hyperslab( space, H5S_SELECT_SET, Offset, NULL, SlabSize, NULL );
                            status = H5Dwrite( DataSet, H5T_NATIVE_DOUBLE, MemSpace, space, H5P_DEFAULT, &U_ARR[0] );
                    }
                    status = H5Sclose( MemSpace );
                    status = H5Dclose( DataSet );
                    status = H5Sclose( space );

                    // Write Pfs_cdmag
                    SlabSize[0] = 1; SlabSize[1] = 3;
                    DataSet  = H5Dopen( file, "Pfs_cdmag", H5P_DEFAULT );
                    space    = H5Dget_space( DataSet );
                    MemSpace = H5Screate_simple( 2, SlabSize, NULL );
                    for ( iT=0; iT<H5_nT; iT++ ){
                            U_ARR[0] = H5_Pfs_cdmag[iT].x; U_ARR[1] = H5_Pfs_cdmag[iT].y; U_ARR[2] = H5_Pfs_cdmag[iT].z;
                            Offset[0]   = iT; Offset[1] = 0;
                            status = H5Sselect_hyperslab( space, H5S_SELECT_SET, Offset, NULL, SlabSize, NULL );
                            status = H5Dwrite( DataSet, H5T_NATIVE_DOUBLE, MemSpace, space, H5P_DEFAULT, &U_ARR[0] );
                    }
                    status = H5Sclose( MemSpace );
                    status = H5Dclose( DataSet );
                    status = H5Sclose( space );

                    // Write Pfs_edmag
                    SlabSize[0] = 1; SlabSize[1] = 3;
                    DataSet  = H5Dopen( file, "Pfs_edmag", H5P_DEFAULT );
                    space    = H5Dget_space( DataSet );
                    MemSpace = H5Screate_simple( 2, SlabSize, NULL );
                    for ( iT=0; iT<H5_nT; iT++ ){
                            U_ARR[0] = H5_Pfs_edmag[iT].x; U_ARR[1] = H5_Pfs_edmag[iT].y; U_ARR[2] = H5_Pfs_edmag[iT].z;
                            Offset[0]   = iT; Offset[1] = 0;
                            status = H5Sselect_hyperslab( space, H5S_SELECT_SET, Offset, NULL, SlabSize, NULL );
                            status = H5Dwrite( DataSet, H5T_NATIVE_DOUBLE, MemSpace, space, H5P_DEFAULT, &U_ARR[0] );
                    }
                    status = H5Sclose( MemSpace );
                    status = H5Dclose( DataSet );
                    status = H5Sclose( space );






                    // Write Kp
                    SlabSize[0] = H5_nT;
                    Offset[0]   = 0;
                    DataSet  = H5Dopen( file, "Kp", H5P_DEFAULT );
                    space    = H5Dget_space( DataSet );
                    status   = H5Sselect_hyperslab( space, H5S_SELECT_SET, Offset, NULL, SlabSize, NULL );
                    MemSpace = H5Screate_simple( 1, SlabSize, NULL );
                    status   = H5Dwrite( DataSet, H5T_NATIVE_DOUBLE, MemSpace, space, H5P_DEFAULT, &H5_Kp[0] );
                    status   = H5Sclose( MemSpace );
                    status   = H5Sclose( space );
                    status   = H5Dclose( DataSet );

                    // Write Dst
                    SlabSize[0] = H5_nT;
                    Offset[0]   = 0;
                    DataSet  = H5Dopen( file, "Dst", H5P_DEFAULT );
                    space    = H5Dget_space( DataSet );
                    status   = H5Sselect_hyperslab( space, H5S_SELECT_SET, Offset, NULL, SlabSize, NULL );
                    MemSpace = H5Screate_simple( 1, SlabSize, NULL );
                    status   = H5Dwrite( DataSet, H5T_NATIVE_DOUBLE, MemSpace, space, H5P_DEFAULT, &H5_Dst[0] );
                    status   = H5Sclose( MemSpace );
                    status   = H5Sclose( space );
                    status   = H5Dclose( DataSet );



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

                    // Write I
                    SlabSize[0] = H5_nT; SlabSize[1] = H5_nAlpha;
                    Offset[0]   = 0; Offset[1] = 0;
                    DataSet  = H5Dopen( file, "I", H5P_DEFAULT );
                    space    = H5Dget_space( DataSet );
                    status   = H5Sselect_hyperslab( space, H5S_SELECT_SET, Offset, NULL, SlabSize, NULL );
                    MemSpace = H5Screate_simple( 2, SlabSize, NULL );
                    status   = H5Dwrite( DataSet, H5T_NATIVE_DOUBLE, MemSpace, space, H5P_DEFAULT, &H5_I[0][0] );
                    status   = H5Sclose( MemSpace );
                    status   = H5Sclose( space );
                    status   = H5Dclose( DataSet );

                    // Write Bm
                    SlabSize[0] = H5_nT; SlabSize[1] = H5_nAlpha;
                    Offset[0]   = 0; Offset[1] = 0;
                    DataSet  = H5Dopen( file, "Bm", H5P_DEFAULT );
                    space    = H5Dget_space( DataSet );
                    status   = H5Sselect_hyperslab( space, H5S_SELECT_SET, Offset, NULL, SlabSize, NULL );
                    MemSpace = H5Screate_simple( 2, SlabSize, NULL );
                    status   = H5Dwrite( DataSet, H5T_NATIVE_DOUBLE, MemSpace, space, H5P_DEFAULT, &H5_Bm[0][0] );
                    status   = H5Sclose( MemSpace );
                    status   = H5Sclose( space );
                    status   = H5Dclose( DataSet );

                    // Write L
                    SlabSize[0] = H5_nT; SlabSize[1] = H5_nAlpha;
                    Offset[0]   = 0; Offset[1] = 0;
                    DataSet  = H5Dopen( file, "L", H5P_DEFAULT );
                    space    = H5Dget_space( DataSet );
                    status   = H5Sselect_hyperslab( space, H5S_SELECT_SET, Offset, NULL, SlabSize, NULL );
                    MemSpace = H5Screate_simple( 2, SlabSize, NULL );
                    status   = H5Dwrite( DataSet, H5T_NATIVE_DOUBLE, MemSpace, space, H5P_DEFAULT, &H5_L[0][0] );
                    status   = H5Sclose( MemSpace );
                    status   = H5Sclose( space );
                    status   = H5Dclose( DataSet );


                    H5Fclose( file );

                } //end else
            } // end "if ( !FileExists || Force )" control structure
        } // end birds loop
    } // end JD loop


    //Lgm_DateTime_Destroy( myUTC );
    free( OutFile );
    free( HdfOutFile );
    free( InFile );

    Lgm_free_ctrans( c );
    Lgm_destroy_eop( e );
    Lgm_FreeMagEphemInfo( MagEphemInfo );

    LGM_ARRAY_2D_FREE( H5_K );
    LGM_ARRAY_2D_FREE( H5_Bm );
    LGM_ARRAY_2D_FREE( H5_I );
    LGM_ARRAY_2D_FREE( H5_L );
    LGM_ARRAY_2D_FREE( H5_Lstar );
    LGM_ARRAY_1D_FREE( H5_Alpha );
    LGM_ARRAY_2D_FREE( Birds );
    LGM_ARRAY_2D_FREE( H5_IsoTimes );
    LGM_ARRAY_1D_FREE( H5_Date );
    LGM_ARRAY_1D_FREE( H5_Doy );
    LGM_ARRAY_1D_FREE( H5_UTC );
    LGM_ARRAY_1D_FREE( H5_JD );
    LGM_ARRAY_1D_FREE( H5_GpsTime );
    LGM_ARRAY_1D_FREE( H5_TiltAngle );
    LGM_ARRAY_1D_FREE( H5_Kp );
    LGM_ARRAY_1D_FREE( H5_Dst );
    LGM_ARRAY_1D_FREE( H5_Rgsm );
    LGM_ARRAY_1D_FREE( H5_Bsc_gsm );
    LGM_ARRAY_1D_FREE( H5_Bfn_geo );
    LGM_ARRAY_1D_FREE( H5_Bfn_gsm );
    LGM_ARRAY_1D_FREE( H5_Pfn_geo );
    LGM_ARRAY_1D_FREE( H5_Pfn_gsm );
    LGM_ARRAY_1D_FREE( H5_Pfn_geod );
    LGM_ARRAY_1D_FREE( H5_Pfn_cdmag );
    LGM_ARRAY_1D_FREE( H5_Pfn_edmag );
    LGM_ARRAY_1D_FREE( H5_Bfs_geo );
    LGM_ARRAY_1D_FREE( H5_Bfs_gsm );
    LGM_ARRAY_1D_FREE( H5_Pfs_geo );
    LGM_ARRAY_1D_FREE( H5_Pfs_gsm );
    LGM_ARRAY_1D_FREE( H5_Pfs_geod );
    LGM_ARRAY_1D_FREE( H5_Pfs_cdmag );
    LGM_ARRAY_1D_FREE( H5_Pfs_edmag );
    LGM_ARRAY_1D_FREE( Perigee_UTC );
    LGM_ARRAY_1D_FREE( Apogee_UTC );
    LGM_ARRAY_1D_FREE( Perigee_U );
    LGM_ARRAY_1D_FREE( Apogee_U );



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
