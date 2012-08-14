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
#include <Lgm_ElapsedTime.h>
#include "SpiceUsr.h"
#include <Lgm/qsort.h>
#include <Lgm_Octree.h>

#define EARTH_ID     399
#define MOON_ID      301
#define RBSPA_ID    -362
#define RBSPB_ID    -363

/*
 * Stuff for inbound/outbound determination
 */
typedef struct TimeList {
    double  key;    // set to JD or somesuch thing.
    int     val;    // set to 0 for perigee time 1 for apogee time
} TimeList;

void elt_qsort( struct TimeList *arr, unsigned n ) {
    #define elt_lt(a,b) ((a)->key < (b)->key)
    QSORT( struct TimeList, arr, n, elt_lt );
}

/*
 * This routine is for figuring out wheter we are inbound or outbound at any given time...
 * Takers:
 *    List of JD times of apogee/perigee crossings and current JD.
 *    The TimeList structure contans a key (JD) and a value indicating apogee or perigee.
 * Returns;
 *    0 if not enough points
 *    1 if outbound
 *   -1 if inbound
 */
int InOutBound( TimeList *a, int n, double JD ){

    int     i, ihi, FoundUpper;

    if ( n<1 )             return( 0 );
    if ( JD < a[0].key )   return( (a[0].val == 1)   ?  1 : -1 );
    if ( JD > a[n-1].key ) return( (a[n-1].val == 1) ? -1 :  1 );


    FoundUpper = FALSE;
    ihi = 0;
    for (i=0; i<n; i++){
        if ( a[i].key > JD ) {
            FoundUpper = TRUE;
            ihi = i;
            break;
        }
    }

    return( (a[ihi].val == 1)   ?  1 : -1 );

}


void StringSplit( char *Str, char *StrArray[], int len, int *n );

#define KP_DEFAULT 0

const  char *ProgramName = "MagEphemFromSpiceKernel";
const  char *argp_program_version     = "MagEphemFromSpiceKernel_1.1";
const  char *argp_program_bug_address = "<mghenderson@lanl.gov>";
static char doc[] =
"Computes the magnetic ephemeris of a S/C from trajectories determined from SPICE kernel files.\n\n"
"The input file is a SPICE kernel description file, which is a file describing all of the SPICE kernels that need to be loaded.  A sample kernel description file is written as;\n\n"

"\t\\begindata\n"
"\tKERNELS_TO_LOAD = ( 'RBSPA_2012_139_2014_213_01_deph.bsp', 'RBSPB_2012_139_2014_213_01_deph.bsp', 'naif0009.tls' )\n"
"\t\\begintext\n\n\n"

"Where 'RBSPA_2012_139_2014_213_01_deph.bsp' is a \"binary SPK\" SPICE kernel, and 'naif0009.tls' is a \"text leap second kernel (LSK)\" file.  In SPICE, SPK kernel files contain information that is used to compute the ephemeris (trajectory) of an object.\n\n"

"The output is written into daily files with filename given by the specified template. OutFile may contain time variables that will be substituted.  The available time variables are '%YYYY', '%MM', and '%DD' which correspond repectively to 4-digit year, 2-digit month (Jan is 01), and 2-digit day of month. \n\n"

"The %B variable will also get substituted by the list of birds given in the -b option.  Here is an example using time-variables.\n\n"

"\t./MagEphemFromSpiceKernelFile -S 20020901 -E 20020930 setup.ker \n"
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
    { 0, 0, 0, 0,   "Time Options:", 1},
    {"StartDateTime",   'S',    "yyyymmdd[Thh:mm:ss]",          0,        "Start date/time in ISO 8601 format. Seconds will be truncated to integers.", 0 },
    {"EndDateTime",     'E',    "yyyymmdd[Thh:mm:ss]",          0,        "End date/time in ISO 8601 format. Seconds will be truncated to integers.", 0 },
    {"Delta",           'D',    "delta",                        0,        "Time cadence in (integer) seconds.", 0 },

    { 0, 0, 0, 0,   "Internal Model Options:", 2},
    {"IntModel",        'i',    "model",                        0,        "Internal Magnetic Field Model to use. Can be CDIP, EDIP, IGRF. Default is IGRF.\n\n", 0},

    { 0, 0, 0, 0,   "External Model Options:", 3},
    {"ExtModel",        'e',    "model",                        0,        "External Magnetic Field Model to use. Can be OP77Q, T87Q, T89Q, T87D, T89D, TS04D, TS07D.\n", 0},
    {"Kp",              'K',    "Kp",                           0,        "If set, force Kp to be this value. Use values like 0.7, 1.0, 1.3 for 1-, 1, 1+" },

/*
    { 0, 0, 0, 0,
"OP77Q\n"
"\tOlsen-Pfitser 1977 Quiet model. Does not depend\n"
"\ton any input parameters.\n"
"T87Q\n"
"\tTsyganenko 1987 model with a fixed Kp\n"
"\t(default is 2.) The Kp level can be over-\n"
"\t-ridden with the --Kp option.\n"
"T89Q\n"
"\tTsyganenko 1989 (T89c version) model\n"
"\t(default is 2.) The Kp level can be over-\n"
"\t-ridden with the --Kp option.\n"
"T87D\n"
"\tTsyganenko 1987 model with Qin-Denton input\n"
"\tparameters.\n"
"T89D\n"
"\tTsyganenko 1989 model with Qin-Denton input\n"
"\tparameters.\n"
"TS04D\n"
"\tTsyganenko-Sitnov 2004 model with Qin-Denton\n"
"\tinput parameters.\n"
"TS07D\n"
"\tTsyganenko-Sitnov 2007 model.\n\n", 0 },
*/

    { 0, 0, 0, 0,   "Other Options:", 4},
    {"Birds",           'b',    "\"bird1, bird2, etc\"",      0,        "Birds (sats) to use. E.g., \"LANL-02A, 1989-046, POLAR\".", 0   },
    {"PitchAngles",     'p',    "\"start_pa, end_pa, npa\"",  0,        "Pitch angles to compute. Default is \"5.0, 90, 18\"." },
    {"FootPointHeight", 'f',    "height",                     0,        "Footpoint height in km. Default is 100km."                  },
    {"Quality",         'q',    "quality",                    0,        "Quality to use for L* calculations. Default is 3."      },
    {"UseEop",          'z',    0,                            0,        "Use Earth Orientation Parameters whn comoputing ephemerii" },
    {"Coords",          'C',    "coord_system",               0,        "Coordinate system used in the input file. Can be: LATLONRAD, SM, GSM, GEI2000 or GSE. Default is LATLONRAD." },

    { 0, 0, 0, 0,   "Output Options:", 5},
    {"Force",           'F',    0,                            0,        "Overwrite output file even if it already exists" },
    {"Verbosity",       'v',    "verbosity",                  0,        "Verbosity level to use. (0-4)"      },
    {"DumpShellFiles",  'd',    0,                            0,        "Dump full binary shell files (for use in viaualizing drift shells)." },
    {"Colorize",        'c',    0,                            0,        "Colorize output"                         },
    {"silent",          's',    0,                            OPTION_ARG_OPTIONAL | OPTION_ALIAS                                                },

    { 0 }
};

struct Arguments {
    char        *args[ nArgs ];       /* START_PA  END_PA  PA_INC */
    int         silent;
    int         Verbosity;

    double      StartPA;
    double      EndPA;
    int         nPA;

    int         Quality;
    double      Kp;
    int         Colorize;
    int         Force;
    double      FootPointHeight;

    char        IntModel[80];
    char        ExtModel[80];
    char        CoordSystem[80];

    int         UseEop;
    int         DumpShellFiles;

    char        Birds[4096];

    long int    Delta;              // Time step cadence

    char        StartDateTime[80];
    char        EndDateTime[80];

    long int    StartDate;          // Start date (e.g. 20130201)
    long int    EndDate;            // End date (e.g. 20130228)

    long int    StartSeconds;       // Start time in seconds for the first date.
    long int    EndSeconds;         // End time in seconds for the last date.
};


/* Parse a single option. */
static error_t parse_opt (int key, char *arg, struct argp_state *state) {

    double       TaiSecs;
    int          ReturnFlag = 0, i;
    char         TimeString[40];
    Lgm_DateTime d;
    Lgm_CTrans   *c = Lgm_init_ctrans( 0 );

    //Lgm_LoadLeapSeconds( c );

    /* Get the input argument from argp_parse, which we
      know is a pointer to our arguments structure. */
    struct Arguments *arguments = state->input;
    switch( key ) {
        case 'b': // Birds
            strncpy( arguments->Birds, arg, 4095 );
            arguments->Birds[4094] = '\0'; // just in case
            for ( i=0; i<strlen(arguments->Birds); i++ ) arguments->Birds[i] = tolower(arguments->Birds[i]);
            break;
        case 'D':
            arguments->Delta = atol( arg );
            break;
        case 'S': // start date
            strncpy( TimeString, arg, 39 );
            IsoTimeStringToDateTime( TimeString, &d, c );
            arguments->StartDate    = d.Date;
            arguments->StartSeconds = d.Hour*3600 + d.Minute*60 + (int)d.Second;
            Lgm_DateTimeToString( arguments->StartDateTime, &d, 0, 0 );
            break;
        case 'E': // end date
            strncpy( TimeString, arg, 39 );
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
        case 'z':
            arguments->UseEop = 1;
            break;
        case 'd':
            arguments->DumpShellFiles = 1;
            break;
        case 'f':
            sscanf( arg, "%lf", &arguments->FootPointHeight );
            break;
        case 'q':
            arguments->Quality = atoi( arg );
            break;
        case 'K':
            sscanf( arg, "%lf", &arguments->Kp );
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
            arguments->Verbosity = atoi( arg );
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
            ReturnFlag = ARGP_ERR_UNKNOWN;
            break;
    }

    Lgm_free_ctrans( c );

    return( ReturnFlag );
}

/* Our argp parser. */
static struct argp argp = { Options, parse_opt, ArgsDoc, doc };



/*
 * Returns Radial distance of S/C.
 */
typedef struct afInfo {
    long int        Date;
    int             Sgn;
    int             BODY;
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
    spkezp_c( a->BODY,    et,   "J2000",  "NONE", EARTH_ID,  pos,  &lt );
    U.x = pos[0]/WGS84_A; U.y = pos[1]/WGS84_A; U.z = pos[2]/WGS84_A;
    R = Lgm_Magnitude( &U );

    return( R*a->Sgn );

}


/*
 *  Compute magnetic ephemerii of S/C from input file that contains S/C position in GEO LAt/Lon/Rad
 */
int main( int argc, char *argv[] ){

    Lgm_ElapsedTimeInfo t;

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
    int              DumpShellFiles, UseEop, Colorize, Force;
    FILE             *fp_in, *fp_MagEphem;
    int              nBirds, iBird;
    char             **Birds, Bird[80];
    double           Inc, Alpha[1000], FootpointHeight;
    int              nAlpha, Quality, Verbosity;
    long int         StartDate, EndDate, Date, Delta;
    int              sYear, sMonth, sDay, sDoy, eYear, eMonth, eDay, eDoy, Year, Month, Day;
    double           sJD, eJD, JD, Time, StartSeconds, EndSeconds;
    Lgm_MagEphemInfo *MagEphemInfo;
    Lgm_QinDentonOne p;

    hid_t           file;
    hid_t           space;
    hid_t           atype;
    hid_t           DataSet, MemSpace;
    herr_t          status;
    hsize_t         Dims[4], Offset[4], SlabSize[4];
    int             iT;
    double          Kp, T89Q_Kp;
    double          GeodLat, GeodLong, GeodHeight, MLAT, MLON, MLT;
    char            **H5_IsoTimes;
    long int        *H5_Date;
    int             *H5_Doy;
    int             *H5_InOut;
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
    double          *H5_S_sc_to_pfn;
    double          *H5_S_sc_to_pfs;
    double          *H5_S_pfs_to_Bmin;
    double          *H5_S_Bmin_to_sc;
    double          *H5_S_total;
    int             H5_nT;
    int             H5_nAlpha;
    double          *H5_Alpha;
    double          **H5_Lstar;
    double          **H5_K;
    double          **H5_I;
    double          **H5_L;
    double          **H5_Bm;

    int             done, BODY;
    long int        ss, es, Seconds, Ta, Tb, Tc;
    double          R, Ra, Rb, Rc, Rmin, Tmin;
    BrentFuncInfo   bInfo;
    afInfo          *afi;
    int             nPerigee, nApogee;
    char            **Perigee_IsoTimes;
    char            **Apogee_IsoTimes;
    Lgm_DateTime    *Perigee_UTC;
    Lgm_DateTime    *Apogee_UTC;
    Lgm_Vector      *Perigee_U;
    Lgm_Vector      *Apogee_U;
    double          **Perigee_Geod;
    double          **Apogee_Geod;
    TimeList        *ApoPeriTimeList;
    int             nApoPeriTimeList;
    Lgm_Vector      Bvec, Bvec2, w;
    int             n, OverRideKp;
    char            *CmdLine;


    t.ColorizeText = TRUE;
    Lgm_ElapsedTimeInit( &t, 255, 150, 0 );


    // kludge.
    LGM_ARRAY_2D( H5_IsoTimes,  2000, 80,    char );
    LGM_ARRAY_1D( H5_Date,      2000,        long int );
    LGM_ARRAY_1D( H5_Doy,       2000,        int );
    LGM_ARRAY_1D( H5_InOut,     2000,        int );
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

    LGM_ARRAY_1D( H5_S_sc_to_pfn,   2000,     double );
    LGM_ARRAY_1D( H5_S_sc_to_pfs,   2000,     double );
    LGM_ARRAY_1D( H5_S_pfs_to_Bmin, 2000,     double );
    LGM_ARRAY_1D( H5_S_Bmin_to_sc,  2000,     double );
    LGM_ARRAY_1D( H5_S_total,       2000,     double );

    LGM_ARRAY_1D( ApoPeriTimeList, 20, TimeList );
    LGM_ARRAY_1D( Perigee_UTC, 10, Lgm_DateTime );
    LGM_ARRAY_1D( Apogee_UTC,  10, Lgm_DateTime );
    LGM_ARRAY_2D( Perigee_Geod, 10, 3, double );
    LGM_ARRAY_2D( Apogee_Geod,  10, 3, double );
    LGM_ARRAY_1D( Perigee_U, 10, Lgm_Vector );
    LGM_ARRAY_1D( Apogee_U,  10, Lgm_Vector );
    LGM_ARRAY_2D( Perigee_IsoTimes, 10, 80, char );
    LGM_ARRAY_2D( Apogee_IsoTimes,  10, 80, char );

    /*
     * Default option values.
     */
    arguments.StartPA         = 90;     // start at 90.0 Deg.
    arguments.EndPA           = 5.0;    // stop at 2.5 Deg.
    arguments.nPA             = 18;     // 18 pitch angles
    arguments.silent          = 0;
    arguments.Verbosity       = 0;
    arguments.Quality         = 3;
    arguments.Kp              = -999.9;
    arguments.Delta           = 60;     // 60s default cadence
    arguments.Colorize        = 0;
    arguments.Force           = 0;
    arguments.UseEop          = 0;
    arguments.DumpShellFiles  = 0;
    arguments.StartDate       = -1;
    arguments.EndDate         = -1;
    arguments.FootPointHeight = 100.0; // km
    strcpy( arguments.IntModel, "IGRF" );
    strcpy( arguments.ExtModel, "T89D" );
    strcpy( arguments.CoordSystem, "LATLONRAD" );


    /*
     * Create a string that shows how we we called.
     */
    for (n=0, i=0; i<argc; i++) n += strlen(argv[i]);
    n += argc;
    LGM_ARRAY_1D( CmdLine, n+1, char );
    for (i=0; i<argc; i++) {
        strcat( CmdLine, argv[i]);  // add all argv items to CmdLine string
        strcat( CmdLine, " ");      // pad with spaces
    }



    /*
     *  Parse CmdLine arguments and options
     */
    argp_parse( &argp, argc, argv, 0, 0, &arguments );



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
    Kp              = arguments.Kp;
    Verbosity       = arguments.Verbosity;
    Colorize        = arguments.Colorize;
    Force           = arguments.Force;
    UseEop          = arguments.UseEop;
    DumpShellFiles  = arguments.DumpShellFiles;
    Delta           = arguments.Delta;
    StartDate       = arguments.StartDate;
    EndDate         = arguments.EndDate;
    StartSeconds    = arguments.StartSeconds;
    EndSeconds      = arguments.EndSeconds;
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
        if ( Kp >= 0.0 ) {
            printf( "\t              Forcing Kp to be: %g\n", Kp );
        }
        printf( "\t          FootpointHeight [km]: %g\n", FootpointHeight );
        printf( "\t                    L* Quality: %d\n", Quality );
        printf( "\t                  Force output: %s\n", Force ? "yes" : "no" );
        printf( "\t                       Use Eop: %s\n", UseEop ? "yes" : "no" );
        printf( "\t         Dump Full Shell Files: %s\n", DumpShellFiles ? "yes" : "no" );
        printf( "\t        Colorize Thread Output: %d\n", Colorize );
        printf( "\t               Verbosity Level: %d\n", Verbosity );
        printf( "\t                        Silent: %s\n", arguments.silent  ? "yes" : "no" );
        if ( (StartDate > 0)&&(EndDate > 0) ){
            printf( "\t                 StartDateTime: %s\n", arguments.StartDateTime );
            printf( "\t                   EndDateTime: %s\n", arguments.EndDateTime );
            printf( "\t      Time Increment [Seconds]: %ld\n", Delta);
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
    Lgm_PrintElapsedTime( &t );


    SpiceDouble et;
    SpiceDouble pos[3], lt;


    if ( nAlpha > 0 ){
        MagEphemInfo = Lgm_InitMagEphemInfo(0, nAlpha);
    } else {
        // doesnt seem to like allocating zero size...
        MagEphemInfo = Lgm_InitMagEphemInfo(0, 1);
    }


    // Settings for Lstar calcs
    MagEphemInfo->LstarQuality = Quality;
    MagEphemInfo->SaveShellLines = TRUE;
    MagEphemInfo->LstarInfo->LSimpleMax = 12.0;
    MagEphemInfo->LstarInfo->VerbosityLevel = Verbosity;
    MagEphemInfo->LstarInfo->mInfo->VerbosityLevel = Verbosity;
    MagEphemInfo->LstarInfo->mInfo->Lgm_LossConeHeight = FootpointHeight;

    OverRideKp = FALSE;
    if ( !strcmp( ExtModel, "T87" ) ){
        MagEphemInfo->LstarInfo->mInfo->Bfield = Lgm_B_T87;
    } else if ( !strcmp( ExtModel, "CDIP" ) ){
        MagEphemInfo->LstarInfo->mInfo->Bfield = Lgm_B_cdip;
    } else if ( !strcmp( ExtModel, "EDIP" ) ){
        MagEphemInfo->LstarInfo->mInfo->Bfield = Lgm_B_edip;
    } else if ( !strcmp( ExtModel, "IGRF" ) ){
        MagEphemInfo->LstarInfo->mInfo->Bfield = Lgm_B_igrf;
    } else if ( !strcmp( ExtModel, "OP77Q" ) ){
        MagEphemInfo->LstarInfo->mInfo->Bfield = Lgm_B_T89;
        Kp = 2;
        OverRideKp = TRUE;
    } else if ( !strcmp( ExtModel, "T89Q" ) ){
        MagEphemInfo->LstarInfo->mInfo->Bfield = Lgm_B_T89;
        T89Q_Kp = 2;
        OverRideKp = TRUE;
    } else if ( !strcmp( ExtModel, "T89D" ) ){
        MagEphemInfo->LstarInfo->mInfo->Bfield = Lgm_B_T89;
    } else if ( !strcmp( ExtModel, "TS04D" ) ){
        MagEphemInfo->LstarInfo->mInfo->Bfield = Lgm_B_TS04;
    } else if ( !strcmp( ExtModel, "TS07D" ) ){
        //MagEphemInfo->LstarInfo->mInfo->Bfield = Lgm_B_TS04;
        printf("TS07D not yet implemented\n");
        exit(0);
    } else { //if ( !strcmp( ExtModel, "T89c" ) ){
        // default
        printf("Unknown model. ExtModel: %s\n", ExtModel );
        exit(0);
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
    char *ShellFile  = (char *)calloc( 2056, sizeof( char ) );

///*
// * Read in a BATS-R-US mesh
// */
//FILE *fp;
//double              x, y, z, Bx, By, Bz, B1x, B1y, B1z;
//long int     rmin_n;
//double r,rmin;
//Lgm_Vector          *u, *B, *B1;
//Lgm_Octree          *Octree;
//LGM_ARRAY_1D( u, 2000000 , Lgm_Vector );
//    LGM_ARRAY_1D( B, 2000000 , Lgm_Vector );
//    LGM_ARRAY_1D( B1, 2000000 , Lgm_Vector );
//n = 0;
//Lgm_ElapsedTimeInit( &t, 255, 150, 0 );
//double theta, phi;
//rmin = 100.0;
//rmin_n = -1;
//if (0==1){
//for (r=1.0; r<1.5; r += 0.1){
//    for (i=0; i<10000; i++){
//        theta = -M_PI/2.0 + M_PI*rand()/(double)RAND_MAX;
//        phi =  2.0*M_PI*rand()/(double)RAND_MAX;
//        u[n].x = r*sin(theta)*cos(phi);
//        u[n].y = r*sin(theta)*sin(phi);
//        u[n].z = r*cos(theta);
//        B1[n].x = 0.0; B1[n].y = 0.0; B1[n].z = 0.0;
//
//        ++n;
//    }
//}
//}
//fp = fopen( "bats_r_us.txt", "r" );
//long int rmin_line_num, linenum = 0;
//while ( fscanf( fp, "%lf %lf %lf %lf %lf %lf %lf %lf %lf", &x, &y, &z, &Bx, &By, &Bz, &B1x, &B1y, &B1z ) != EOF ) {
//    if ( ( fabs(B1x) > 1e-9 ) && ( fabs(B1y) >  1e-9 ) && ( fabs(B1y) > 1e-9 ) ) {
//            u[n].x  = x; u[n].y  = y; u[n].z  = z;
//            B[n].x  = Bx; B[n].y  = By; B[n].z  = Bz;
//            B1[n].x = B1x; B1[n].y = B1y; B1[n].z = B1z;
//r = Lgm_Magnitude(&u[n]);
//if (r<rmin){
//rmin = r;
//rmin_n = n;
//rmin_line_num = linenum;
//}
//            ++n;
//    }
//    ++linenum;
//}
//fclose(fp);
//printf("Number of points in Mesh: %ld\n", n);
//Lgm_PrintElapsedTime( &t );
//printf("rmin = %g rmin_n = %ld rmin_line_num = %ld\n", rmin, rmin_line_num, rmin_n);
////exit(0);
///*
//* Create Octree
//*/
//printf("Creating Octree\n");
//Lgm_ElapsedTimeInit( &t, 255, 150, 0 );
//Octree = Lgm_InitOctree( u, B1, n );
//printf("Min, Max, Diff = %g %g %g\n", Octree->Min, Octree->Max, Octree->Diff);
//Lgm_PrintElapsedTime( &t );
//
//Lgm_MagModelInfo_Set_Octree( Octree, 8, MagEphemInfo->LstarInfo->mInfo );
//
//
//Lgm_MagModelInfo_Set_MagModel( LGM_CDIP, LGM_EXTMODEL_SCATTERED_DATA2, MagEphemInfo->LstarInfo->mInfo );
////MagEphemInfo->LstarInfo->mInfo->Lgm_MagStep_Integrator = LGM_MAGSTEP_ODE_BS;
//MagEphemInfo->LstarInfo->mInfo->Lgm_MagStep_RK5_Eps = 1e-1;
//MagEphemInfo->LstarInfo->mInfo->Lgm_MagStep_RK5_Eps = 1e-4;
//MagEphemInfo->LstarInfo->mInfo->Lgm_MagStep_BS_Eps  = 1e-4;
//MagEphemInfo->LstarInfo->mInfo->Lgm_TraceLine_Tol   = 1e-6;
//Lgm_B_FromScatteredData_SetUp( MagEphemInfo->LstarInfo->mInfo );





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

            if        ( !strcmp( Bird, "rbspa" ) ){
                BODY = RBSPA_ID;
            } else if ( !strcmp( Bird, "rbspb" ) ){
                BODY = RBSPB_ID;
            } else {
                // maybe user just has a SPICE BODY number?
                BODY = atoi( Bird );
            }

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
            printf( "      TXT Output File: %s\n", OutFile);
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
                    printf("\n\n\tWarning. Existing HdfOutfile will be overwritten: %s \n\n", HdfOutFile );
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

                    afi->BODY = BODY;
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
                    nApoPeriTimeList = 0;
                    while ( !done ) {
                        if ( (Rb < Ra) && (Rb < Rc) ) {
                            // we have a bracket. Find Min.
                            Lgm_Brent( (double)Ta, (double)Tb, (double)Tc, &bInfo, 1e-3, &Tmin, &Rmin );
                            Lgm_Make_UTC( Date, Tmin/3600.0, &Apogee_UTC[nApogee], c );
                            et = Lgm_TDBSecSinceJ2000( &Apogee_UTC[nApogee], c );
                            spkezp_c( BODY,    et,   "J2000",  "NONE", EARTH_ID,  pos,  &lt );
                            Apogee_U[nApogee].x = pos[0]/WGS84_A; Apogee_U[nApogee].y = pos[1]/WGS84_A; Apogee_U[nApogee].z = pos[2]/WGS84_A;

                            Lgm_Set_Coord_Transforms( Apogee_UTC[nApogee].Date, Apogee_UTC[nApogee].Time, c );
                            Lgm_Convert_Coords( &Apogee_U[nApogee], &w, GEI2000_TO_WGS84, c );
                            Lgm_WGS84_to_GEOD( &w, &Apogee_Geod[nApogee][0], &Apogee_Geod[nApogee][1], &Apogee_Geod[nApogee][2] );

                            Lgm_DateTimeToString( Apogee_IsoTimes[nApogee], &Apogee_UTC[nApogee], 0, 3 );
                            //printf("nApogee: %d Bracket: T = %ld %ld %ld    R = %g %g %g    Tapogee, Rapogee = %s %g\n", nApogee, Ta, Tb, Tc, Ra, Rb, Rc, Apogee_IsoTimes[nApogee], fabs(Rmin) );
                            ApoPeriTimeList[nApoPeriTimeList].key = Apogee_UTC[nApogee].JD;
                            ApoPeriTimeList[nApoPeriTimeList].val = 1; // because its apogee
                            ++nApoPeriTimeList;

                            ++nApogee;

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
                            spkezp_c( BODY,    et,   "J2000",  "NONE", EARTH_ID,  pos,  &lt );
                            Perigee_U[nPerigee].x = pos[0]/WGS84_A; Perigee_U[nPerigee].y = pos[1]/WGS84_A; Perigee_U[nPerigee].z = pos[2]/WGS84_A;

                            Lgm_Set_Coord_Transforms( Perigee_UTC[nPerigee].Date, Perigee_UTC[nPerigee].Time, c );
                            Lgm_Convert_Coords( &Perigee_U[nPerigee], &w, GEI2000_TO_WGS84, c );
                            Lgm_WGS84_to_GEOD( &w, &Perigee_Geod[nPerigee][0], &Perigee_Geod[nPerigee][1], &Perigee_Geod[nPerigee][2] );

                            Lgm_DateTimeToString( Perigee_IsoTimes[nPerigee], &Perigee_UTC[nPerigee], 0, 3 );
                            //printf("nPerigee: %d Bracket: T = %ld %ld %ld    R = %g %g %g    Tperigee, Rperigee = %s %g\n", nPerigee, Ta, Tb, Tc, Ra, Rb, Rc, Perigee_IsoTimes[nPerigee], fabs(Rmin) );

                            ApoPeriTimeList[nApoPeriTimeList].key = Perigee_UTC[nPerigee].JD;
                            ApoPeriTimeList[nApoPeriTimeList].val = 0; // because its perigee
                            ++nApoPeriTimeList;

                            ++nPerigee;
                        }
                        // advance another 900 seconds.
                        Ta = Tb; Ra = Rb;
                        Tb = Tc; Rb = Rc;
                        Tc += 900.0; Rc = bInfo.func( (double)Tc, 0.0, (void *)afi );
                        if (Tc >= 86400) done = TRUE;
                    }

                    free( afi );


                    /*
                     * sort the merged list of apogee/perigee times.
                     */
                    elt_qsort( ApoPeriTimeList, nApoPeriTimeList );



                    /*
                     * Open Mag Ephem file for writing and
                     * Write header to Mag Ephem file
                     */
                    /*
                     */
                    fp_MagEphem = fopen( OutFile, "w" );
//printf("nPerigee = %d\n", nPerigee);
                    Lgm_WriteMagEphemHeader( fp_MagEphem, Bird, 9, "FIX ME", CmdLine, nPerigee, Perigee_UTC, Perigee_U, nApogee, &Apogee_UTC[0], &Apogee_U[0], MagEphemInfo );
                    printf("\t      Writing to file: %s\n", OutFile );

                    if ( UseEop ) {
                        // Read in the EOP vals
                        Lgm_read_eop( e );
                    }




                    /*
                     *  Loop over seconds of the day
                     */
                    ss = (Date == StartDate) ? StartSeconds : 0;
                    es = (Date == EndDate) ? EndSeconds : 86400;
                    H5_nT = 0;
                    Lgm_ElapsedTimeInit( &t, 255, 150, 0 );
                    for ( Seconds=ss; Seconds<=es; Seconds += Delta ) {

                        Lgm_Make_UTC( Date, Seconds/3600.0, &UTC, c );
                        Lgm_DateTimeToString( IsoTimeString, &UTC, 0, 0 );

                        et = Lgm_TDBSecSinceJ2000( &UTC, c );
                        spkezp_c( BODY,    et,   "J2000",  "NONE", EARTH_ID,  pos,  &lt );
                        //printf("pos = %.8lf %.8lf %.8lf\n", pos[0], pos[1], pos[2]);
                        U.x = pos[0]/WGS84_A; U.y = pos[1]/WGS84_A; U.z = pos[2]/WGS84_A;


                        if ( UseEop ) {
                            // Get (interpolate) the EOP vals from the values in the file at the given Julian Date
                            Lgm_get_eop_at_JD( UTC.JD, &eop, e );

                            // Set the EOP vals in the CTrans structure.
                            Lgm_set_eop( &eop, c );
                        }

                        // Set mag model parameters
                        Lgm_get_QinDenton_at_JD( UTC.JD-3153.0, &p, 1 ); // for date 20120616 this puts us back to halloween storm (Oct 29, 2003)
                        //Lgm_get_QinDenton_at_JD( Lgm_JD( 2003, 10, 30, 12.0, LGM_TIME_SYS_UTC, c ), &p, 1 ); // for date 20120616 this puts us back to halloween storm (Oct 29, 2003)
                        Lgm_set_QinDenton( &p, MagEphemInfo->LstarInfo->mInfo );

//MIKE
                        if ( Kp >= 0.0 ) {
                            MagEphemInfo->LstarInfo->mInfo->fKp = Kp;
                            MagEphemInfo->LstarInfo->mInfo->Kp  = (int)(Kp+0.5);
                            if (MagEphemInfo->LstarInfo->mInfo->Kp > 5) MagEphemInfo->LstarInfo->mInfo->Kp = 5;
                            if (MagEphemInfo->LstarInfo->mInfo->Kp < 0 ) MagEphemInfo->LstarInfo->mInfo->Kp = 0;
                        }


                        // Set up the trans matrices
                        Lgm_Set_Coord_Transforms( UTC.Date, UTC.Time, c );

                        Lgm_Convert_Coords( &U, &Rgsm, GEI2000_TO_GSM, c );



                            /*
                             * Compute L*s, Is, Bms, Footprints, etc...
                             * These quantities are stored in the MagEphemInfo Structure
                             */
                            printf("\t\t"); Lgm_PrintElapsedTime( &t ); printf("\n");
                            printf("\n\n\t[ %s ]: %s  Bird: %s Rgsm: %g %g %g Re\n", ProgramName, IsoTimeString, Bird, Rgsm.x, Rgsm.y, Rgsm.z );
                            printf("\t--------------------------------------------------------------------------------------------------\n");
//Lgm_Vector TMPTMP;
//Lgm_Set_Coord_Transforms( UTC.Date, UTC.Time, MagEphemInfo->LstarInfo->mInfo->c );
//Lgm_TraceToMinBSurf( &Rgsm, &TMPTMP, 0.1, 1e-7, MagEphemInfo->LstarInfo->mInfo );
//Lgm_Setup_AlphaOfK( &UTC, &TMPTMP, MagEphemInfo->LstarInfo->mInfo );
//printf("Lgm_AlphaOfK( 0.1 ) = %g\n", Lgm_AlphaOfK( 0.1, MagEphemInfo->LstarInfo->mInfo ) );
//Lgm_TearDown_AlphaOfK( MagEphemInfo->LstarInfo->mInfo );
//exit(0);
                            Lgm_ComputeLstarVersusPA( UTC.Date, UTC.Time, &Rgsm, nAlpha, Alpha, MagEphemInfo->LstarQuality, Colorize, MagEphemInfo );
//                            Lgm_ComputeLstarVersusPA( UTC.Date, UTC.Time, &TMPTMP, nAlpha, Alpha, MagEphemInfo->LstarQuality, Colorize, MagEphemInfo );

                            MagEphemInfo->InOut = InOutBound( ApoPeriTimeList, nApoPeriTimeList, UTC.JD );
                            Lgm_WriteMagEphemData( fp_MagEphem, IntModel, ExtModel, MagEphemInfo->LstarInfo->mInfo->fKp, MagEphemInfo->LstarInfo->mInfo->Dst, MagEphemInfo );


                            if ( DumpShellFiles && (nAlpha > 0) ){

                                sprintf( ShellFile, "%s_%ld.dat", OutFile, Seconds );
                                printf( "Writing Full Shell File: %s\n", ShellFile );
                                WriteMagEphemInfoStruct( ShellFile, nAlpha, MagEphemInfo );
                            }




                            // Fill arrays for dumping out as HDF5 files
                            strcpy( H5_IsoTimes[H5_nT], IsoTimeString );
                            H5_Date[H5_nT]          = UTC.Date;
                            H5_Doy[H5_nT]           = UTC.Doy;
                            H5_UTC[H5_nT]           = UTC.Time;
                            H5_JD[H5_nT]            = UTC.JD;
                            H5_InOut[H5_nT]         = MagEphemInfo->InOut;
                            H5_GpsTime[H5_nT]       = Lgm_UTC_to_GpsSeconds( &UTC, c );
                            H5_TiltAngle[H5_nT]     = c->psi*DegPerRad;
                            H5_Rgsm[H5_nT]          = Rgsm;
                            H5_Kp[H5_nT]            = MagEphemInfo->LstarInfo->mInfo->fKp;
                            H5_Dst[H5_nT]           = MagEphemInfo->LstarInfo->mInfo->Dst;
                            H5_S_sc_to_pfn[H5_nT]   = (MagEphemInfo->Snorth > 0.0) ? MagEphemInfo->Snorth : LGM_FILL_VALUE;
                            H5_S_sc_to_pfs[H5_nT]   = (MagEphemInfo->Ssouth > 0.0) ? MagEphemInfo->Ssouth : LGM_FILL_VALUE;
                            H5_S_pfs_to_Bmin[H5_nT] = (MagEphemInfo->Smin > 0.0) ? MagEphemInfo->Smin : LGM_FILL_VALUE;
                            H5_S_Bmin_to_sc[H5_nT]  = ((MagEphemInfo->Ssouth>0.0)&&(MagEphemInfo->Smin > 0.0)) ? MagEphemInfo->Ssouth-MagEphemInfo->Smin : LGM_FILL_VALUE;
                            H5_S_total[H5_nT]       = ((MagEphemInfo->Snorth > 0.0)&&(MagEphemInfo->Ssouth > 0.0)) ? MagEphemInfo->Snorth + MagEphemInfo->Ssouth : LGM_FILL_VALUE;
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

                    Lgm_PrintElapsedTime( &t );
                    Lgm_SetElapsedTimeStr( &t );
                    sprintf( Command, "sed 's/ELAPSED_TIME/%s/' <%s >%s.new", t.ElapsedTimeStr, OutFile, OutFile); system( Command );
                    sprintf( Command, "mv %s.new %s", OutFile, OutFile); system( Command );




                    // Create HDF5 file
                    file    = H5Fcreate( HdfOutFile, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );

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

                    // Create PerigeeTimes Dataset
                    Dims[0] = nPerigee;
                    space   = H5Screate_simple( 1, Dims, NULL ); // rank 1
                    atype   = H5Tcopy( H5T_C_S1 );
                    status  = H5Tset_size( atype, strlen(Perigee_IsoTimes[0]) );
                    status  = H5Tset_strpad( atype, H5T_STR_NULLPAD );
                    status  = H5Tset_cset( atype, H5T_CSET_ASCII );
                    DataSet = H5Dcreate( file, "PerigeeTimes", atype, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "UTC" );
                    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "The times of Perigee in ISO 8601 compliant format. Defined as smallest geocentric distance." );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );
                    status  = H5Tclose( atype );

                    // Create PerigeePosGeod Dataset
                    Dims[0] = nPerigee; Dims[1] = 3;
                    space   = H5Screate_simple( 2, Dims, NULL ); // rank 2
                    DataSet = H5Dcreate( file, "PerigeePosGeod", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "PerigeeTimes" );
                    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Geodetic position of perigee (lat/lon/rad)." );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "Deg./Deg./km" );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );


                    // Create ApogeeTimes Dataset
                    Dims[0] = nApogee;
                    space   = H5Screate_simple( 1, Dims, NULL ); // rank 1
                    atype   = H5Tcopy( H5T_C_S1 );
                    status  = H5Tset_size( atype, strlen(Apogee_IsoTimes[0]) );
                    status  = H5Tset_strpad( atype, H5T_STR_NULLPAD );
                    status  = H5Tset_cset( atype, H5T_CSET_ASCII );
                    DataSet = H5Dcreate( file, "ApogeeTimes", atype, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "The times of Perigee in ISO 8601 compliant format. Defined as largest geocentric distance." );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "UTC" );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );
                    status  = H5Tclose( atype );

                    // Create ApogeePosGeod Dataset
                    Dims[0] = nApogee; Dims[1] = 3;
                    space   = H5Screate_simple( 2, Dims, NULL ); // rank 2
                    DataSet = H5Dcreate( file, "ApogeePosGeod", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "ApogeeTimes" );
                    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Geodetic position of apogee (lat/lon/rad)." );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "Deg./Deg./km" );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );


                    // Create IsoTime Dataset
                    Dims[0] = H5_nT;
                    space   = H5Screate_simple( 1, Dims, NULL ); // rank 1
                    atype   = H5Tcopy( H5T_C_S1 );
                    status  = H5Tset_size( atype, strlen(H5_IsoTimes[0]) );
                    status  = H5Tset_strpad( atype, H5T_STR_NULLPAD );
                    status  = H5Tset_cset( atype, H5T_CSET_ASCII );
                    DataSet = H5Dcreate( file, "IsoTime", atype, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "The date and time in ISO 8601 compliant format." );
                    Lgm_WriteStringAttr( DataSet, "UNITS",       "UTC" );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",    "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",     "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",    "data" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );
                    status  = H5Tclose( atype );

                    // Create Date Dataset
                    Dims[0] = H5_nT;
                    space   = H5Screate_simple( 1, Dims, NULL ); // rank 1
                    DataSet = H5Dcreate( file, "Date", H5T_NATIVE_LONG, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "The date. In YYYMMDD format." );
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "YYYYMMDD" );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );

                    // Create Doy Dataset
                    Dims[0] = H5_nT;
                    space   = H5Screate_simple( 1, Dims, NULL ); // rank 1
                    DataSet = H5Dcreate( file, "Doy", H5T_NATIVE_INT, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Ordinal Day of Year." );
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "DDD" );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );

                    // Create UTC (hours) Dataset
                    Dims[0] = H5_nT;
                    space   = H5Screate_simple( 1, Dims, NULL ); // rank 1
                    DataSet = H5Dcreate( file, "UTC", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Universal Time (Coordinated). In decimal hours." );
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "Hours" );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );

                    // Create JD Dataset
                    Dims[0] = H5_nT;
                    space   = H5Screate_simple( 1, Dims, NULL ); // rank 1
                    DataSet = H5Dcreate( file, "JulianDate", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Julian Date. In decimal days." );
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "Days" );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );

                    // Create GpsTime Dataset
                    Dims[0] = H5_nT;
                    space   = H5Screate_simple( 1, Dims, NULL ); // rank 1
                    DataSet = H5Dcreate( file, "GpsTime", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Number of SI seconds since 0h Jan 6, 1980 UTC." );
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "Seconds" );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );

                    // Create DipoleTiltAngle Dataset
                    Dims[0] = H5_nT;
                    space   = H5Screate_simple( 1, Dims, NULL ); // rank 1
                    DataSet = H5Dcreate( file, "DipoleTiltAngle", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Angle between Zgsm and Zsm (i.e. between Zgsm and dipole axis direction)." );
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "Degrees" );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );

                    // Create Rgeo Dataset
                    Dims[0] = H5_nT; Dims[1] = 3;
                    space   = H5Screate_simple( 2, Dims, NULL ); // rank 2
                    DataSet = H5Dcreate( file, "Rgeo", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Geocentric Geographic position vector of S/C." );
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "Re" );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );

                    // Create Rgeod_LatLon Dataset
                    Dims[0] = H5_nT; Dims[1] = 2;
                    space   = H5Screate_simple( 2, Dims, NULL ); // rank 2
                    DataSet = H5Dcreate( file, "Rgeod_LatLon", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Geodetic Geographic Latitude and Longitude of S/C." );
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "Deg." );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );

                    // Create Rgeod_Height Dataset
                    Dims[0] = H5_nT;
                    space   = H5Screate_simple( 1, Dims, NULL ); // rank 1
                    DataSet = H5Dcreate( file, "Rgeod_Height", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Geodetic Geographic Height (Above WGS84 Ellipsoid) of S/C." );
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "km" );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );


                    // Create Rgsm Dataset
                    Dims[0] = H5_nT; Dims[1] = 3;
                    space   = H5Screate_simple( 2, Dims, NULL ); // rank 2
                    DataSet = H5Dcreate( file, "Rgsm", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Geocentric Solar Magnetospheric position vector of S/C." );
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "Re" );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );

                    // Create Rsm Dataset
                    Dims[0] = H5_nT; Dims[1] = 3;
                    space   = H5Screate_simple( 2, Dims, NULL ); // rank 2
                    DataSet = H5Dcreate( file, "Rsm", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Geocentric Solar Magnetic position vector of S/C." );
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "Re" );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );

                    // Create Rgei Dataset
                    Dims[0] = H5_nT; Dims[1] = 3;
                    space   = H5Screate_simple( 2, Dims, NULL ); // rank 2
                    DataSet = H5Dcreate( file, "Rgei", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Geocentric Equatorial Inertial position vector of S/C." );
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "Re" );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );

                    // Create Rgse Dataset
                    Dims[0] = H5_nT; Dims[1] = 3;
                    space   = H5Screate_simple( 2, Dims, NULL ); // rank 2
                    DataSet = H5Dcreate( file, "Rgse", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Geocentric Solar Ecliptic position vector of S/C." );
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "Re" );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );



                    // Create CDMAG_MLAT Dataset
                    Dims[0] = H5_nT;
                    space   = H5Screate_simple( 1, Dims, NULL ); // rank 1
                    DataSet = H5Dcreate( file, "CDMAG_MLAT", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Magnetic Latitude of S/C in Centerted Dipole Coordinates." );
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "Deg." );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );

                    // Create CDMAG_MLON Dataset
                    Dims[0] = H5_nT;
                    space   = H5Screate_simple( 1, Dims, NULL ); // rank 1
                    DataSet = H5Dcreate( file, "CDMAG_MLON", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Magnetic Longitude of S/C in Centerted Dipole Coordinates." );
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "Deg." );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );

                    // Create CDMAG_MLT Dataset
                    Dims[0] = H5_nT;
                    space   = H5Screate_simple( 1, Dims, NULL ); // rank 1
                    DataSet = H5Dcreate( file, "CDMAG_MLT", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Magnetic Local Time of S/C in Centerted Dipole Coordinates." );
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "Deg." );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );

                    // Create CDMAG_R Dataset
                    Dims[0] = H5_nT;
                    space   = H5Screate_simple( 1, Dims, NULL ); // rank 1
                    DataSet = H5Dcreate( file, "CDMAG_R", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Radial distance of S/C from center of CDMAG coordinate system." );
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "Deg." );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );


                    // Create EDMAG_MLAT Dataset
                    Dims[0] = H5_nT;
                    space   = H5Screate_simple( 1, Dims, NULL ); // rank 1
                    DataSet = H5Dcreate( file, "EDMAG_MLAT", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Magnetic Latitude of S/C in Eccentric Dipole Coordinates." );
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "Deg." );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );

                    // Create EDMAG_MLON Dataset
                    Dims[0] = H5_nT;
                    space   = H5Screate_simple( 1, Dims, NULL ); // rank 1
                    DataSet = H5Dcreate( file, "EDMAG_MLON", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Magnetic Longitude of S/C in Eccentric Dipole Coordinates." );
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "Deg." );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );

                    // Create EDMAG_MLT Dataset
                    Dims[0] = H5_nT;
                    space   = H5Screate_simple( 1, Dims, NULL ); // rank 1
                    DataSet = H5Dcreate( file, "EDMAG_MLT", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Magnetic Local Time of S/C in Eccentric Dipole Coordinates." );
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "Deg." );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );

                    // Create EDMAG_R Dataset
                    Dims[0] = H5_nT;
                    space   = H5Screate_simple( 1, Dims, NULL ); // rank 1
                    DataSet = H5Dcreate( file, "EDMAG_R", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Radial distance of S/C from center of EDMAG coordinate system." );
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "Deg." );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );




                    // Create IntModel Dataset
                    Dims[0] = H5_nT;
                    space   = H5Screate_simple( 1, Dims, NULL ); // rank 1
                    atype   = H5Tcopy( H5T_C_S1 );
                    status  = H5Tset_size( atype, 30 );
                    status  = H5Tset_strpad( atype, H5T_STR_NULLPAD );
                    status  = H5Tset_cset( atype, H5T_CSET_ASCII );
                    DataSet = H5Dcreate( file, "IntModel", atype, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Internal magnetic field model." );
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );
                    status  = H5Tclose( atype );


                    // Create ExtModel Dataset
                    Dims[0] = H5_nT;
                    space   = H5Screate_simple( 1, Dims, NULL ); // rank 1
                    atype   = H5Tcopy( H5T_C_S1 );
                    status  = H5Tset_size( atype, 30 );
                    status  = H5Tset_strpad( atype, H5T_STR_NULLPAD );
                    status  = H5Tset_cset( atype, H5T_CSET_ASCII );
                    DataSet = H5Dcreate( file, "ExtModel", atype, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "External magnetic field model." );
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );
                    status  = H5Tclose( atype );


                    // Create Kp Dataset
                    Dims[0] = H5_nT;
                    space   = H5Screate_simple( 1, Dims, NULL ); // rank 1
                    DataSet = H5Dcreate( file, "Kp", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Kp index value." );
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
                    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Dst index value." );
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "nT" );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );


                    // Create Bsc_gsm Dataset
                    Dims[0] = H5_nT; Dims[1] = 4;
                    space   = H5Screate_simple( 2, Dims, NULL ); // rank 2
                    DataSet = H5Dcreate( file, "Bsc_gsm", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Magnetic field vector at S/C (in GSM coords)." );
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "nT" );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );



                    // Create FieldLineType Dataset
                    Dims[0] = H5_nT;
                    space   = H5Screate_simple( 1, Dims, NULL ); // rank 1
                    atype   = H5Tcopy( H5T_C_S1 );
                    status  = H5Tset_size( atype, 30 );
                    status  = H5Tset_strpad( atype, H5T_STR_NULLPAD );
                    status  = H5Tset_cset( atype, H5T_CSET_ASCII );
                    DataSet = H5Dcreate( file, "FieldLineType", atype, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Description of the type of field line the S/C is on., Can be one of 4 types: LGM_CLOSED      - FL hits Earth at both ends. LGM_OPEN_N_LOBE - FL is an OPEN field line rooted in the Northern polar cap. LGM_OPEN_S_LOBE - FL is an OPEN field line rooted in the Southern polar cap. LGM_OPEN_IMF    - FL does not hit Earth at eitrher end." );
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );
                    status  = H5Tclose( atype );


                    // Create S_sc_to_pfn Dataset
                    Dims[0] = H5_nT;
                    space   = H5Screate_simple( 1, Dims, NULL ); // rank 2
                    DataSet = H5Dcreate( file, "S_sc_to_pfn", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Distance between S/C and Northern Footpoint along field line." );
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "Re" );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );

                    // Create S_sc_to_pfs Dataset
                    Dims[0] = H5_nT;
                    space   = H5Screate_simple( 1, Dims, NULL ); // rank 2
                    DataSet = H5Dcreate( file, "S_sc_to_pfs", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Distance between S/C and Southern Footpoint along field line." );
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "Re" );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );

                    // Create S_pfs_to_Bmin Dataset
                    Dims[0] = H5_nT;
                    space   = H5Screate_simple( 1, Dims, NULL ); // rank 2
                    DataSet = H5Dcreate( file, "S_pfs_to_Bmin", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Distance between Southern Footpoint and Bmin point along field line.");
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "Re" );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );

                    // Create S_Bmin_to_sc Dataset
                    Dims[0] = H5_nT;
                    space   = H5Screate_simple( 1, Dims, NULL ); // rank 2
                    DataSet = H5Dcreate( file, "S_Bmin_to_sc", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Distance between Bmin point and S/C along field line (positive if north of Bmin).");
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "Re" );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );

                    // Create S_total Dataset
                    Dims[0] = H5_nT;
                    space   = H5Screate_simple( 1, Dims, NULL ); // rank 2
                    DataSet = H5Dcreate( file, "S_total", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Total Field Line length (along field line).");
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "Re" );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );


                    // Create Pfn_geo Dataset
                    Dims[0] = H5_nT; Dims[1] = 3;
                    space   = H5Screate_simple( 2, Dims, NULL ); // rank 2
                    DataSet = H5Dcreate( file, "Pfn_geo", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Location of Northern Footpoint (in GEO coords).");
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "Re" );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );

                    // Create Pfn_gsm Dataset
                    Dims[0] = H5_nT; Dims[1] = 3;
                    space   = H5Screate_simple( 2, Dims, NULL ); // rank 2
                    DataSet = H5Dcreate( file, "Pfn_gsm", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Location of Northern Footpoint (in GSM coords).");
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "Re" );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );



                    // Create Pfn_geod_LatLon Dataset
                    Dims[0] = H5_nT; Dims[1] = 2;
                    space   = H5Screate_simple( 2, Dims, NULL ); // rank 2
                    DataSet = H5Dcreate( file, "Pfn_geod_LatLon", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Geodetic Latitude and Longitude of Northern Footpoint.");
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "Deg." );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );

                    // Create Pfn_geod_Height Dataset
                    Dims[0] = H5_nT;
                    space   = H5Screate_simple( 1, Dims, NULL ); // rank 2
                    DataSet = H5Dcreate( file, "Pfn_geod_Height", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Geodetic Height of Northern Footpoint.");
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "km" );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );


                    // Create Pfn_CD_MLAT Dataset
                    Dims[0] = H5_nT;
                    space   = H5Screate_simple( 1, Dims, NULL ); // rank 1
                    DataSet = H5Dcreate( file, "Pfn_CD_MLAT", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Magnetic Latitude of Northern Footpoint in Centerted Dipole Coordinates.");
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "Deg." );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );

                    // Create Pfn_CD_MLON Dataset
                    Dims[0] = H5_nT;
                    space   = H5Screate_simple( 1, Dims, NULL ); // rank 1
                    DataSet = H5Dcreate( file, "Pfn_CD_MLON", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Magnetic Longitude of Northern Footpoint in Centerted Dipole Coordinates.");
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "Deg." );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );

                    // Create Pfn_CD_MLT Dataset
                    Dims[0] = H5_nT;
                    space   = H5Screate_simple( 1, Dims, NULL ); // rank 1
                    DataSet = H5Dcreate( file, "Pfn_CD_MLT", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Magnetic Local Time of Northern Footpoint in Centerted Dipole Coordinates.");
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "Hours" );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );

                    // Create Pfn_ED_MLAT Dataset
                    Dims[0] = H5_nT;
                    space   = H5Screate_simple( 1, Dims, NULL ); // rank 1
                    DataSet = H5Dcreate( file, "Pfn_ED_MLAT", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Magnetic Latitude of Northern Footpoint in Eccentric Dipole Coordinates.");
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "Deg." );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );

                    // Create Pfn_ED_MLON Dataset
                    Dims[0] = H5_nT;
                    space   = H5Screate_simple( 1, Dims, NULL ); // rank 1
                    DataSet = H5Dcreate( file, "Pfn_ED_MLON", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Magnetic Longitude of Northern Footpoint in Eccentric Dipole Coordinates.");
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "Deg." );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );

                    // Create Pfn_ED_MLT Dataset
                    Dims[0] = H5_nT;
                    space   = H5Screate_simple( 1, Dims, NULL ); // rank 1
                    DataSet = H5Dcreate( file, "Pfn_ED_MLT", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Magnetic Local Time of Northern Footpoint in Eccentric Dipole Coordinates.");
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "Hours" );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );

                    // Create Bfn_geo Dataset
                    Dims[0] = H5_nT; Dims[1] = 4;
                    space   = H5Screate_simple( 2, Dims, NULL ); // rank 2
                    DataSet = H5Dcreate( file, "Bfn_geo", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Magnetic field vector at Northern Footpoint (in GEO coords).");
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "nT" );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );

                    // Create Bfn_gsm Dataset
                    Dims[0] = H5_nT; Dims[1] = 4;
                    space   = H5Screate_simple( 2, Dims, NULL ); // rank 2
                    DataSet = H5Dcreate( file, "Bfn_gsm", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Magnetic field vector at Northern Footpoint (in GSM coords).");
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "nT" );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );

                    // Create Loss_Cone_Alpha_n Dataset
                    Dims[0] = H5_nT;
                    space   = H5Screate_simple( 1, Dims, NULL ); // rank 1
                    DataSet = H5Dcreate( file, "Loss_Cone_Alpha_n", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Value of Northern Loss Cone angle. asin( sqrt(Bsc/Bfn) ).");
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "Deg." );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );







                    // Create Pfs_geo Dataset
                    Dims[0] = H5_nT; Dims[1] = 3;
                    space   = H5Screate_simple( 2, Dims, NULL ); // rank 2
                    DataSet = H5Dcreate( file, "Pfs_geo", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Location of Southern Footpoint (in GEO coords).");
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "Re" );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );

                    // Create Pfs_gsm Dataset
                    Dims[0] = H5_nT; Dims[1] = 3;
                    space   = H5Screate_simple( 2, Dims, NULL ); // rank 2
                    DataSet = H5Dcreate( file, "Pfs_gsm", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Location of Southern Footpoint (in GSM coords).");
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "Re" );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );



                    // Create Pfs_geod_LatLon Dataset
                    Dims[0] = H5_nT; Dims[1] = 2;
                    space   = H5Screate_simple( 2, Dims, NULL ); // rank 2
                    DataSet = H5Dcreate( file, "Pfs_geod_LatLon", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Geodetic Latitude and Longitude of Southern Footpoint.");
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "Deg." );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );

                    // Create Pfs_geod_Height Dataset
                    Dims[0] = H5_nT;
                    space   = H5Screate_simple( 1, Dims, NULL ); // rank 2
                    DataSet = H5Dcreate( file, "Pfs_geod_Height", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Geodetic Height of Southern Footpoint.");
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "km" );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );


                    // Create Pfs_CD_MLAT Dataset
                    Dims[0] = H5_nT;
                    space   = H5Screate_simple( 1, Dims, NULL ); // rank 1
                    DataSet = H5Dcreate( file, "Pfs_CD_MLAT", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Magnetic Latitude of Southern Footpoint in Centerted Dipole Coordinates.");
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "Deg." );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );

                    // Create Pfs_CD_MLON Dataset
                    Dims[0] = H5_nT;
                    space   = H5Screate_simple( 1, Dims, NULL ); // rank 1
                    DataSet = H5Dcreate( file, "Pfs_CD_MLON", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Magnetic Longitude of Southern Footpoint in Centerted Dipole Coordinates.");
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "Deg." );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );

                    // Create Pfs_CD_MLT Dataset
                    Dims[0] = H5_nT;
                    space   = H5Screate_simple( 1, Dims, NULL ); // rank 1
                    DataSet = H5Dcreate( file, "Pfs_CD_MLT", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Magnetic Local Time of Southern Footpoint in Centerted Dipole Coordinates.");
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "Hours" );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );

                    // Create Pfs_ED_MLAT Dataset
                    Dims[0] = H5_nT;
                    space   = H5Screate_simple( 1, Dims, NULL ); // rank 1
                    DataSet = H5Dcreate( file, "Pfs_ED_MLAT", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Magnetic Latitude of Southern Footpoint in Eccentric Dipole Coordinates.");
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "Deg." );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );

                    // Create Pfs_ED_MLON Dataset
                    Dims[0] = H5_nT;
                    space   = H5Screate_simple( 1, Dims, NULL ); // rank 1
                    DataSet = H5Dcreate( file, "Pfs_ED_MLON", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Magnetic Longitude of Southern Footpoint in Eccentric Dipole Coordinates.");
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "Deg." );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );

                    // Create Pfs_ED_MLT Dataset
                    Dims[0] = H5_nT;
                    space   = H5Screate_simple( 1, Dims, NULL ); // rank 1
                    DataSet = H5Dcreate( file, "Pfs_ED_MLT", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Magnetic Local Time of Southern Footpoint in Eccentric Dipole Coordinates.");
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "Hours" );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );

                    // Create Bfs_geo Dataset
                    Dims[0] = H5_nT; Dims[1] = 4;
                    space   = H5Screate_simple( 2, Dims, NULL ); // rank 2
                    DataSet = H5Dcreate( file, "Bfs_geo", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Magnetic field vector at Southern Footpoint (in GEO coords).");
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "nT" );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );

                    // Create Bfs_gsm Dataset
                    Dims[0] = H5_nT; Dims[1] = 4;
                    space   = H5Screate_simple( 2, Dims, NULL ); // rank 2
                    DataSet = H5Dcreate( file, "Bfs_gsm", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Magnetic field vector at Southern Footpoint (in GSM coords).");
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "nT" );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );

                    // Create Loss_Cone_Alpha_s Dataset
                    Dims[0] = H5_nT;
                    space   = H5Screate_simple( 1, Dims, NULL ); // rank 1
                    DataSet = H5Dcreate( file, "Loss_Cone_Alpha_s", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Value of Southern Loss Cone angle. asin( sqrt(Bsc/Bfs) ).");
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "Deg." );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );


                    // Create Pmin_gsm Dataset
                    Dims[0] = H5_nT; Dims[1] = 3;
                    space   = H5Screate_simple( 2, Dims, NULL ); // rank 2
                    DataSet = H5Dcreate( file, "Pmin_gsm", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Location of minimum-|B| point (in GSM coords).");
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "Re" );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );

                    // Create Bmin_gsm Dataset
                    Dims[0] = H5_nT; Dims[1] = 4;
                    space   = H5Screate_simple( 2, Dims, NULL ); // rank 2
                    DataSet = H5Dcreate( file, "Bmin_gsm", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "B-field at minimum-|B| point (in GSM coords).");
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "Re" );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );



                    // Create Lsimple Dataset
                    Dims[0] = H5_nT;
                    space   = H5Screate_simple( 1, Dims, NULL ); // rank 1
                    DataSet = H5Dcreate( file, "Lsimple", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Geocentric distance to Bmin point for FL threading vehicle (i.e. |Pmin|).");
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "Dimless" );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );

                    // Create InvLat Dataset
                    Dims[0] = H5_nT;
                    space   = H5Screate_simple( 1, Dims, NULL ); // rank 1
                    DataSet = H5Dcreate( file, "InvLat", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Invariant latitude of vehicle computed from Lambda=acos(sqrt(1/Lsimple)).");
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "Deg." );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );


                    // Create Lm_eq Dataset
                    Dims[0] = H5_nT;
                    space   = H5Screate_simple( 1, Dims, NULL ); // rank 1
                    DataSet = H5Dcreate( file, "Lm_eq", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "McIlwain L of an eq. mirroring particle on same FL as vehicle (computed from L=Lm_eq, I=0, and Bm=|Bmin_gsm|, M=M_igrf).");
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "Dimless" );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );

                    // Create InvLat_eq Dataset
                    Dims[0] = H5_nT;
                    space   = H5Screate_simple( 1, Dims, NULL ); // rank 1
                    DataSet = H5Dcreate( file, "InvLat_eq", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Invariant latitude of vehicle computed from Lambda=acos(sqrt(1.0/Lm_eq)).");
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "Deg." );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );


                    // Create BoverBeq Dataset
                    Dims[0] = H5_nT;
                    space   = H5Screate_simple( 1, Dims, NULL ); // rank 1
                    DataSet = H5Dcreate( file, "BoverBeq", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Magntiude of Bsc over magnitude of Bmin.");
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "Dimless" );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );

                    // Create MlatFromBoverBeq Dataset
                    Dims[0] = H5_nT;
                    space   = H5Screate_simple( 1, Dims, NULL ); // rank 1
                    DataSet = H5Dcreate( file, "MlatFromBoverBeq", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Dipole latitude where (B/Beq)_dipole == BoverBeq.");
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "Deg." );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );


                    // Create M_used Dataset
                    Dims[0] = H5_nT;
                    space   = H5Screate_simple( 1, Dims, NULL ); // rank 1
                    DataSet = H5Dcreate( file, "M_used", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "The magnetic dipole moment that was used to convert magnetic flux to L*. In units of nT.");
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "nT" );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );

                    // Create M_ref Dataset
                    Dims[0] = H5_nT;
                    space   = H5Screate_simple( 1, Dims, NULL ); // rank 1
                    DataSet = H5Dcreate( file, "M_ref", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "The fixed reference magnetic dipole moment for converting magnetic flux to L*. In units of nT.");
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "nT" );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );

                    // Create M_igrf Dataset
                    Dims[0] = H5_nT;
                    space   = H5Screate_simple( 1, Dims, NULL ); // rank 1
                    DataSet = H5Dcreate( file, "M_igrf", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Time-dependant magnetic dipole moment (probably shouldn't be used for converting magnetic flux to L*, but it sometimes is). In units of nT.");
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
                    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Generalized Roederer L-shell value (also known as L*).");
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "DEPEND_1",   "Alpha" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "Dimensionless" );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );

                    // Create L Dataset
                    Dims[0] = H5_nT; Dims[1] = H5_nAlpha;
                    space   = H5Screate_simple( 2, Dims, NULL ); // rank 2
                    DataSet = H5Dcreate( file, "L", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "McIlwain L-shell value.");
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "DEPEND_1",   "Alpha" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "Dimensionless" );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );

                    // Create Bm Dataset
                    Dims[0] = H5_nT; Dims[1] = H5_nAlpha;
                    space   = H5Screate_simple( 2, Dims, NULL ); // rank 2
                    DataSet = H5Dcreate( file, "Bm", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Magnetic field strength at mirror points for each pitch angle.");
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "DEPEND_1",   "Alpha" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "nT" );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "log" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );

                    // Create I Dataset
                    Dims[0] = H5_nT; Dims[1] = H5_nAlpha;
                    space   = H5Screate_simple( 2, Dims, NULL ); // rank 2
                    DataSet = H5Dcreate( file, "I", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Integral invariant for each pitch angle.");
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "DEPEND_1",   "Alpha" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "Re" );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );

                    // Create K Dataset
                    Dims[0] = H5_nT; Dims[1] = H5_nAlpha;
                    space   = H5Screate_simple( 2, Dims, NULL ); // rank 2
                    DataSet = H5Dcreate( file, "K", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Second Invariant ( I*sqrt(Bm) ) for each pitch angle.");
                    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
                    Lgm_WriteStringAttr( DataSet, "DEPEND_1",   "Alpha" );
                    Lgm_WriteStringAttr( DataSet, "UNITS",      "Re G^.5" );
                    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "log" );
                    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
                    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
                    Lgm_WriteStringAttr( DataSet, "VALIDMIN",   "0.0" );
                    status  = H5Sclose( space );
                    status  = H5Dclose( DataSet );




                    /*
                     * Write Datasets
                     */

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


                    // Write Apogee Time Strings
                    SlabSize[0] = 1;
                    DataSet  = H5Dopen( file, "ApogeeTimes", H5P_DEFAULT );
                    space    = H5Dget_space( DataSet );
                    MemSpace = H5Screate_simple( 1, SlabSize, NULL );
                    atype    = H5Tcopy( H5T_C_S1 );
                    status  = H5Tset_size( atype, strlen(Apogee_IsoTimes[0]) );
                    status  = H5Tset_strpad( atype, H5T_STR_NULLPAD );
                    status  = H5Tset_cset( atype, H5T_CSET_ASCII );
                    for ( iT=0; iT<nApogee; iT++ ){
                            Offset[0]   = iT;
                            status = H5Sselect_hyperslab( space, H5S_SELECT_SET, Offset, NULL, SlabSize, NULL );
                            status = H5Dwrite( DataSet, atype, MemSpace, space, H5P_DEFAULT, &Apogee_IsoTimes[iT][0] );
                    }
                    status = H5Sclose( MemSpace );
                    status = H5Dclose( DataSet );
                    status = H5Sclose( space );
                    status = H5Tclose( atype );

                    // Write ApogeePosGeod
                    SlabSize[0] = 1; SlabSize[1] = 3;
                    DataSet  = H5Dopen( file, "ApogeePosGeod", H5P_DEFAULT );
                    space    = H5Dget_space( DataSet );
                    MemSpace = H5Screate_simple( 2, SlabSize, NULL );
                    for ( iT=0; iT<nApogee; iT++ ){
                            U_ARR[0] = Apogee_Geod[iT][0]; U_ARR[1] = Apogee_Geod[iT][1]; U_ARR[2] = Apogee_Geod[iT][2];
printf("U_ARR = %g %g %g\n", U_ARR[0], U_ARR[1], U_ARR[2]);
                            Offset[0]   = iT; Offset[1] = 0;
                            status = H5Sselect_hyperslab( space, H5S_SELECT_SET, Offset, NULL, SlabSize, NULL );
                            status = H5Dwrite( DataSet, H5T_NATIVE_DOUBLE, MemSpace, space, H5P_DEFAULT, &U_ARR[0] );
                    }
                    status = H5Sclose( MemSpace );
                    status = H5Dclose( DataSet );
                    status = H5Sclose( space );




                    // Write Perigee Time Strings
                    SlabSize[0] = 1;
                    DataSet  = H5Dopen( file, "PerigeeTimes", H5P_DEFAULT );
                    space    = H5Dget_space( DataSet );
                    MemSpace = H5Screate_simple( 1, SlabSize, NULL );
                    atype    = H5Tcopy( H5T_C_S1 );
                    status  = H5Tset_size( atype, strlen(Perigee_IsoTimes[0]) );
                    status  = H5Tset_strpad( atype, H5T_STR_NULLPAD );
                    status  = H5Tset_cset( atype, H5T_CSET_ASCII );
                    for ( iT=0; iT<nPerigee; iT++ ){
                            Offset[0]   = iT;
                            status = H5Sselect_hyperslab( space, H5S_SELECT_SET, Offset, NULL, SlabSize, NULL );
                            status = H5Dwrite( DataSet, atype, MemSpace, space, H5P_DEFAULT, &Perigee_IsoTimes[iT][0] );
                    }
                    status = H5Sclose( MemSpace );
                    status = H5Dclose( DataSet );
                    status = H5Sclose( space );
                    status = H5Tclose( atype );

                    // Write PerigeePosGeod
                    SlabSize[0] = 1; SlabSize[1] = 3;
                    DataSet  = H5Dopen( file, "PerigeePosGeod", H5P_DEFAULT );
                    space    = H5Dget_space( DataSet );
                    MemSpace = H5Screate_simple( 2, SlabSize, NULL );
                    for ( iT=0; iT<nPerigee; iT++ ){
                            U_ARR[0] = Perigee_Geod[iT][0]; U_ARR[1] = Perigee_Geod[iT][1]; U_ARR[2] = Perigee_Geod[iT][2];
                            Offset[0]   = iT; Offset[1] = 0;
                            status = H5Sselect_hyperslab( space, H5S_SELECT_SET, Offset, NULL, SlabSize, NULL );
                            status = H5Dwrite( DataSet, H5T_NATIVE_DOUBLE, MemSpace, space, H5P_DEFAULT, &U_ARR[0] );
                    }
                    status = H5Sclose( MemSpace );
                    status = H5Dclose( DataSet );
                    status = H5Sclose( space );



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
                    DataSet  = H5Dopen( file, "DipoleTiltAngle", H5P_DEFAULT );
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



                    // Write Rgeod_LatLon
                    SlabSize[0] = 1; SlabSize[1] = 2;
                    DataSet  = H5Dopen( file, "Rgeod_LatLon", H5P_DEFAULT );
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

                    // Write Rgeod_Height
                    SlabSize[0] = 1;
                    DataSet  = H5Dopen( file, "Rgeod_Height", H5P_DEFAULT );
                    space    = H5Dget_space( DataSet );
                    MemSpace = H5Screate_simple( 1, SlabSize, NULL );
                    for ( iT=0; iT<H5_nT; iT++ ){
                            Lgm_Set_Coord_Transforms( H5_Date[iT], H5_UTC[iT], c );
                            Lgm_Convert_Coords( &H5_Rgsm[iT], &W, GSM_TO_GEO, c ); // Convert to cartesian GEO coords.
                            Lgm_WGS84_to_GEOD( &W, &GeodLat, &GeodLong, &GeodHeight );
                            U_ARR[0] = GeodLat; U_ARR[1] = GeodLong; U_ARR[2] = GeodHeight;
                            Offset[0]   = iT;
                            status = H5Sselect_hyperslab( space, H5S_SELECT_SET, Offset, NULL, SlabSize, NULL );
                            status = H5Dwrite( DataSet, H5T_NATIVE_DOUBLE, MemSpace, space, H5P_DEFAULT, &U_ARR[2] );
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







                    // Write CDMAG_MLAT
                    SlabSize[0] = 1;
                    DataSet  = H5Dopen( file, "CDMAG_MLAT", H5P_DEFAULT );
                    space    = H5Dget_space( DataSet );
                    MemSpace = H5Screate_simple( 1, SlabSize, NULL );
                    for ( iT=0; iT<H5_nT; iT++ ){
                            Lgm_Set_Coord_Transforms( H5_Date[iT], H5_UTC[iT], c );
                            Lgm_Convert_Coords( &H5_Rgsm[iT], &W, GSM_TO_CDMAG, c ); // Convert to cartesian CDMAG coords.
                            Lgm_CDMAG_to_R_MLAT_MLON_MLT( &W, &R, &MLAT, &MLON, &MLT, c );
                            U_ARR[0] = MLAT; U_ARR[1] = MLON; U_ARR[2] = MLT; U_ARR[3] = R;
                            Offset[0]   = iT;
                            status = H5Sselect_hyperslab( space, H5S_SELECT_SET, Offset, NULL, SlabSize, NULL );
                            status = H5Dwrite( DataSet, H5T_NATIVE_DOUBLE, MemSpace, space, H5P_DEFAULT, &U_ARR[0] );
                    }
                    status = H5Sclose( MemSpace );
                    status = H5Dclose( DataSet );
                    status = H5Sclose( space );

                    // Write CDMAG_MLON
                    SlabSize[0] = 1;
                    DataSet  = H5Dopen( file, "CDMAG_MLON", H5P_DEFAULT );
                    space    = H5Dget_space( DataSet );
                    MemSpace = H5Screate_simple( 1, SlabSize, NULL );
                    for ( iT=0; iT<H5_nT; iT++ ){
                            Lgm_Set_Coord_Transforms( H5_Date[iT], H5_UTC[iT], c );
                            Lgm_Convert_Coords( &H5_Rgsm[iT], &W, GSM_TO_CDMAG, c ); // Convert to cartesian CDMAG coords.
                            Lgm_CDMAG_to_R_MLAT_MLON_MLT( &W, &R, &MLAT, &MLON, &MLT, c );
                            U_ARR[0] = MLAT; U_ARR[1] = MLON; U_ARR[2] = MLT; U_ARR[3] = R;
                            Offset[0]   = iT;
                            status = H5Sselect_hyperslab( space, H5S_SELECT_SET, Offset, NULL, SlabSize, NULL );
                            status = H5Dwrite( DataSet, H5T_NATIVE_DOUBLE, MemSpace, space, H5P_DEFAULT, &U_ARR[1] );
                    }
                    status = H5Sclose( MemSpace );
                    status = H5Dclose( DataSet );
                    status = H5Sclose( space );

                    // Write CDMAG_MLT
                    SlabSize[0] = 1;
                    DataSet  = H5Dopen( file, "CDMAG_MLT", H5P_DEFAULT );
                    space    = H5Dget_space( DataSet );
                    MemSpace = H5Screate_simple( 1, SlabSize, NULL );
                    for ( iT=0; iT<H5_nT; iT++ ){
                            Lgm_Set_Coord_Transforms( H5_Date[iT], H5_UTC[iT], c );
                            Lgm_Convert_Coords( &H5_Rgsm[iT], &W, GSM_TO_CDMAG, c ); // Convert to cartesian CDMAG coords.
                            Lgm_CDMAG_to_R_MLAT_MLON_MLT( &W, &R, &MLAT, &MLON, &MLT, c );
                            U_ARR[0] = MLAT; U_ARR[1] = MLON; U_ARR[2] = MLT; U_ARR[3] = R;
                            Offset[0]   = iT;
                            status = H5Sselect_hyperslab( space, H5S_SELECT_SET, Offset, NULL, SlabSize, NULL );
                            status = H5Dwrite( DataSet, H5T_NATIVE_DOUBLE, MemSpace, space, H5P_DEFAULT, &U_ARR[2] );
                    }
                    status = H5Sclose( MemSpace );
                    status = H5Dclose( DataSet );
                    status = H5Sclose( space );

                    // Write CDMAG_R
                    SlabSize[0] = 1;
                    DataSet  = H5Dopen( file, "CDMAG_R", H5P_DEFAULT );
                    space    = H5Dget_space( DataSet );
                    MemSpace = H5Screate_simple( 1, SlabSize, NULL );
                    for ( iT=0; iT<H5_nT; iT++ ){
                            Lgm_Set_Coord_Transforms( H5_Date[iT], H5_UTC[iT], c );
                            Lgm_Convert_Coords( &H5_Rgsm[iT], &W, GSM_TO_CDMAG, c ); // Convert to cartesian CDMAG coords.
                            Lgm_CDMAG_to_R_MLAT_MLON_MLT( &W, &R, &MLAT, &MLON, &MLT, c );
                            U_ARR[0] = MLAT; U_ARR[1] = MLON; U_ARR[2] = MLT; U_ARR[3] = R;
                            Offset[0]   = iT;
                            status = H5Sselect_hyperslab( space, H5S_SELECT_SET, Offset, NULL, SlabSize, NULL );
                            status = H5Dwrite( DataSet, H5T_NATIVE_DOUBLE, MemSpace, space, H5P_DEFAULT, &U_ARR[3] );
                    }
                    status = H5Sclose( MemSpace );
                    status = H5Dclose( DataSet );
                    status = H5Sclose( space );



                    // Write EDMAG_MLAT
                    SlabSize[0] = 1;
                    DataSet  = H5Dopen( file, "EDMAG_MLAT", H5P_DEFAULT );
                    space    = H5Dget_space( DataSet );
                    MemSpace = H5Screate_simple( 1, SlabSize, NULL );
                    for ( iT=0; iT<H5_nT; iT++ ){
                            Lgm_Set_Coord_Transforms( H5_Date[iT], H5_UTC[iT], c );
                            Lgm_Convert_Coords( &H5_Rgsm[iT], &W, GSM_TO_EDMAG, c ); // Convert to cartesian EDMAG coords.
                            Lgm_EDMAG_to_R_MLAT_MLON_MLT( &W, &R, &MLAT, &MLON, &MLT, c );
                            U_ARR[0] = MLAT; U_ARR[1] = MLON; U_ARR[2] = MLT; U_ARR[3] = R;
                            Offset[0]   = iT;
                            status = H5Sselect_hyperslab( space, H5S_SELECT_SET, Offset, NULL, SlabSize, NULL );
                            status = H5Dwrite( DataSet, H5T_NATIVE_DOUBLE, MemSpace, space, H5P_DEFAULT, &U_ARR[0] );
                    }
                    status = H5Sclose( MemSpace );
                    status = H5Dclose( DataSet );
                    status = H5Sclose( space );

                    // Write EDMAG_MLON
                    SlabSize[0] = 1;
                    DataSet  = H5Dopen( file, "EDMAG_MLON", H5P_DEFAULT );
                    space    = H5Dget_space( DataSet );
                    MemSpace = H5Screate_simple( 1, SlabSize, NULL );
                    for ( iT=0; iT<H5_nT; iT++ ){
                            Lgm_Set_Coord_Transforms( H5_Date[iT], H5_UTC[iT], c );
                            Lgm_Convert_Coords( &H5_Rgsm[iT], &W, GSM_TO_EDMAG, c ); // Convert to cartesian EDMAG coords.
                            Lgm_EDMAG_to_R_MLAT_MLON_MLT( &W, &R, &MLAT, &MLON, &MLT, c );
                            U_ARR[0] = MLAT; U_ARR[1] = MLON; U_ARR[2] = MLT; U_ARR[3] = R;
                            Offset[0]   = iT;
                            status = H5Sselect_hyperslab( space, H5S_SELECT_SET, Offset, NULL, SlabSize, NULL );
                            status = H5Dwrite( DataSet, H5T_NATIVE_DOUBLE, MemSpace, space, H5P_DEFAULT, &U_ARR[1] );
                    }
                    status = H5Sclose( MemSpace );
                    status = H5Dclose( DataSet );
                    status = H5Sclose( space );

                    // Write EDMAG_MLT
                    SlabSize[0] = 1;
                    DataSet  = H5Dopen( file, "EDMAG_MLT", H5P_DEFAULT );
                    space    = H5Dget_space( DataSet );
                    MemSpace = H5Screate_simple( 1, SlabSize, NULL );
                    for ( iT=0; iT<H5_nT; iT++ ){
                            Lgm_Set_Coord_Transforms( H5_Date[iT], H5_UTC[iT], c );
                            Lgm_Convert_Coords( &H5_Rgsm[iT], &W, GSM_TO_EDMAG, c ); // Convert to cartesian EDMAG coords.
                            Lgm_EDMAG_to_R_MLAT_MLON_MLT( &W, &R, &MLAT, &MLON, &MLT, c );
                            U_ARR[0] = MLAT; U_ARR[1] = MLON; U_ARR[2] = MLT; U_ARR[3] = R;
                            Offset[0]   = iT;
                            status = H5Sselect_hyperslab( space, H5S_SELECT_SET, Offset, NULL, SlabSize, NULL );
                            status = H5Dwrite( DataSet, H5T_NATIVE_DOUBLE, MemSpace, space, H5P_DEFAULT, &U_ARR[2] );
                    }
                    status = H5Sclose( MemSpace );
                    status = H5Dclose( DataSet );
                    status = H5Sclose( space );

                    // Write EDMAG_R
                    SlabSize[0] = 1;
                    DataSet  = H5Dopen( file, "EDMAG_R", H5P_DEFAULT );
                    space    = H5Dget_space( DataSet );
                    MemSpace = H5Screate_simple( 1, SlabSize, NULL );
                    for ( iT=0; iT<H5_nT; iT++ ){
                            Lgm_Set_Coord_Transforms( H5_Date[iT], H5_UTC[iT], c );
                            Lgm_Convert_Coords( &H5_Rgsm[iT], &W, GSM_TO_EDMAG, c ); // Convert to cartesian EDMAG coords.
                            Lgm_EDMAG_to_R_MLAT_MLON_MLT( &W, &R, &MLAT, &MLON, &MLT, c );
                            U_ARR[0] = MLAT; U_ARR[1] = MLON; U_ARR[2] = MLT; U_ARR[3] = R;
                            Offset[0]   = iT;
                            status = H5Sselect_hyperslab( space, H5S_SELECT_SET, Offset, NULL, SlabSize, NULL );
                            status = H5Dwrite( DataSet, H5T_NATIVE_DOUBLE, MemSpace, space, H5P_DEFAULT, &U_ARR[3] );
                    }
                    status = H5Sclose( MemSpace );
                    status = H5Dclose( DataSet );
                    status = H5Sclose( space );


                    // Write IntModel





                    // Write ExtModel


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



                    // Write FieldLineType




                    // Write S_sc_to_pfn
                    SlabSize[0] = H5_nT;
                    Offset[0]   = 0;
                    DataSet  = H5Dopen( file, "S_sc_to_pfn", H5P_DEFAULT );
                    space    = H5Dget_space( DataSet );
                    status   = H5Sselect_hyperslab( space, H5S_SELECT_SET, Offset, NULL, SlabSize, NULL );
                    MemSpace = H5Screate_simple( 1, SlabSize, NULL );
                    status   = H5Dwrite( DataSet, H5T_NATIVE_DOUBLE, MemSpace, space, H5P_DEFAULT, &H5_S_sc_to_pfn[0] );
                    status   = H5Sclose( MemSpace );
                    status   = H5Sclose( space );
                    status   = H5Dclose( DataSet );

                    // Write S_sc_to_pfs
                    SlabSize[0] = H5_nT;
                    Offset[0]   = 0;
                    DataSet  = H5Dopen( file, "S_sc_to_pfs", H5P_DEFAULT );
                    space    = H5Dget_space( DataSet );
                    status   = H5Sselect_hyperslab( space, H5S_SELECT_SET, Offset, NULL, SlabSize, NULL );
                    MemSpace = H5Screate_simple( 1, SlabSize, NULL );
                    status   = H5Dwrite( DataSet, H5T_NATIVE_DOUBLE, MemSpace, space, H5P_DEFAULT, &H5_S_sc_to_pfs[0] );
                    status   = H5Sclose( MemSpace );
                    status   = H5Sclose( space );
                    status   = H5Dclose( DataSet );

                    // Write S_pfs_to_Bmin
                    SlabSize[0] = H5_nT;
                    Offset[0]   = 0;
                    DataSet  = H5Dopen( file, "S_pfs_to_Bmin", H5P_DEFAULT );
                    space    = H5Dget_space( DataSet );
                    status   = H5Sselect_hyperslab( space, H5S_SELECT_SET, Offset, NULL, SlabSize, NULL );
                    MemSpace = H5Screate_simple( 1, SlabSize, NULL );
                    status   = H5Dwrite( DataSet, H5T_NATIVE_DOUBLE, MemSpace, space, H5P_DEFAULT, &H5_S_pfs_to_Bmin[0] );
                    status   = H5Sclose( MemSpace );
                    status   = H5Sclose( space );
                    status   = H5Dclose( DataSet );

                    // Write S_Bmin_to_sc
                    SlabSize[0] = H5_nT;
                    Offset[0]   = 0;
                    DataSet  = H5Dopen( file, "S_Bmin_to_sc", H5P_DEFAULT );
                    space    = H5Dget_space( DataSet );
                    status   = H5Sselect_hyperslab( space, H5S_SELECT_SET, Offset, NULL, SlabSize, NULL );
                    MemSpace = H5Screate_simple( 1, SlabSize, NULL );
                    status   = H5Dwrite( DataSet, H5T_NATIVE_DOUBLE, MemSpace, space, H5P_DEFAULT, &H5_S_Bmin_to_sc[0] );
                    status   = H5Sclose( MemSpace );
                    status   = H5Sclose( space );
                    status   = H5Dclose( DataSet );

                    // Write S_total
                    SlabSize[0] = H5_nT;
                    Offset[0]   = 0;
                    DataSet  = H5Dopen( file, "S_total", H5P_DEFAULT );
                    space    = H5Dget_space( DataSet );
                    status   = H5Sselect_hyperslab( space, H5S_SELECT_SET, Offset, NULL, SlabSize, NULL );
                    MemSpace = H5Screate_simple( 1, SlabSize, NULL );
                    status   = H5Dwrite( DataSet, H5T_NATIVE_DOUBLE, MemSpace, space, H5P_DEFAULT, &H5_S_total[0] );
                    status   = H5Sclose( MemSpace );
                    status   = H5Sclose( space );
                    status   = H5Dclose( DataSet );


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

                    // Write Pfn_geod_LatLon
                    SlabSize[0] = 1; SlabSize[1] = 2;
                    DataSet  = H5Dopen( file, "Pfn_geod_LatLon", H5P_DEFAULT );
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

                    // Write Pfn_geod_Height
                    SlabSize[0] = 1;
                    DataSet  = H5Dopen( file, "Pfn_geod_Height", H5P_DEFAULT );
                    space    = H5Dget_space( DataSet );
                    MemSpace = H5Screate_simple( 1, SlabSize, NULL );
                    for ( iT=0; iT<H5_nT; iT++ ){
                            U_ARR[0] = H5_Pfn_geod[iT].x; U_ARR[1] = H5_Pfn_geod[iT].y; U_ARR[2] = H5_Pfn_geod[iT].z;
                            Offset[0]   = iT;
                            status = H5Sselect_hyperslab( space, H5S_SELECT_SET, Offset, NULL, SlabSize, NULL );
                            status = H5Dwrite( DataSet, H5T_NATIVE_DOUBLE, MemSpace, space, H5P_DEFAULT, &U_ARR[2] );
                    }
                    status = H5Sclose( MemSpace );
                    status = H5Dclose( DataSet );
                    status = H5Sclose( space );

                    // Write Pfn_CD_MLAT
                    SlabSize[0] = 1;
                    DataSet  = H5Dopen( file, "Pfn_CD_MLAT", H5P_DEFAULT );
                    space    = H5Dget_space( DataSet );
                    MemSpace = H5Screate_simple( 1, SlabSize, NULL );
                    for ( iT=0; iT<H5_nT; iT++ ){
                            U_ARR[0] = H5_Pfn_cdmag[iT].x; U_ARR[1] = H5_Pfn_cdmag[iT].y; U_ARR[2] = H5_Pfn_cdmag[iT].z;
                            Offset[0]   = iT;
                            status = H5Sselect_hyperslab( space, H5S_SELECT_SET, Offset, NULL, SlabSize, NULL );
                            status = H5Dwrite( DataSet, H5T_NATIVE_DOUBLE, MemSpace, space, H5P_DEFAULT, &U_ARR[0] );
                    }
                    status = H5Sclose( MemSpace );
                    status = H5Dclose( DataSet );
                    status = H5Sclose( space );

                    // Write Pfn_CD_MLON
                    SlabSize[0] = 1;
                    DataSet  = H5Dopen( file, "Pfn_CD_MLON", H5P_DEFAULT );
                    space    = H5Dget_space( DataSet );
                    MemSpace = H5Screate_simple( 1, SlabSize, NULL );
                    for ( iT=0; iT<H5_nT; iT++ ){
                            U_ARR[0] = H5_Pfn_cdmag[iT].x; U_ARR[1] = H5_Pfn_cdmag[iT].y; U_ARR[2] = H5_Pfn_cdmag[iT].z;
                            Offset[0]   = iT;
                            status = H5Sselect_hyperslab( space, H5S_SELECT_SET, Offset, NULL, SlabSize, NULL );
                            status = H5Dwrite( DataSet, H5T_NATIVE_DOUBLE, MemSpace, space, H5P_DEFAULT, &U_ARR[0] );
                    }
                    status = H5Sclose( MemSpace );
                    status = H5Dclose( DataSet );
                    status = H5Sclose( space );

                    // Write Pfn_CD_MLT
                    SlabSize[0] = 1;
                    DataSet  = H5Dopen( file, "Pfn_CD_MLT", H5P_DEFAULT );
                    space    = H5Dget_space( DataSet );
                    MemSpace = H5Screate_simple( 1, SlabSize, NULL );
                    for ( iT=0; iT<H5_nT; iT++ ){
                            U_ARR[0] = H5_Pfn_cdmag[iT].x; U_ARR[1] = H5_Pfn_cdmag[iT].y; U_ARR[2] = H5_Pfn_cdmag[iT].z;
                            Offset[0]   = iT;
                            status = H5Sselect_hyperslab( space, H5S_SELECT_SET, Offset, NULL, SlabSize, NULL );
                            status = H5Dwrite( DataSet, H5T_NATIVE_DOUBLE, MemSpace, space, H5P_DEFAULT, &U_ARR[0] );
                    }
                    status = H5Sclose( MemSpace );
                    status = H5Dclose( DataSet );
                    status = H5Sclose( space );


                    // Write Pfn_ED_MLAT
                    SlabSize[0] = 1;
                    DataSet  = H5Dopen( file, "Pfn_ED_MLAT", H5P_DEFAULT );
                    space    = H5Dget_space( DataSet );
                    MemSpace = H5Screate_simple( 1, SlabSize, NULL );
                    for ( iT=0; iT<H5_nT; iT++ ){
                            U_ARR[0] = H5_Pfn_cdmag[iT].x; U_ARR[1] = H5_Pfn_cdmag[iT].y; U_ARR[2] = H5_Pfn_cdmag[iT].z;
                            Offset[0]   = iT;
                            status = H5Sselect_hyperslab( space, H5S_SELECT_SET, Offset, NULL, SlabSize, NULL );
                            status = H5Dwrite( DataSet, H5T_NATIVE_DOUBLE, MemSpace, space, H5P_DEFAULT, &U_ARR[0] );
                    }
                    status = H5Sclose( MemSpace );
                    status = H5Dclose( DataSet );
                    status = H5Sclose( space );

                    // Write Pfn_ED_MLON
                    SlabSize[0] = 1;
                    DataSet  = H5Dopen( file, "Pfn_ED_MLON", H5P_DEFAULT );
                    space    = H5Dget_space( DataSet );
                    MemSpace = H5Screate_simple( 1, SlabSize, NULL );
                    for ( iT=0; iT<H5_nT; iT++ ){
                            U_ARR[0] = H5_Pfn_cdmag[iT].x; U_ARR[1] = H5_Pfn_cdmag[iT].y; U_ARR[2] = H5_Pfn_cdmag[iT].z;
                            Offset[0]   = iT;
                            status = H5Sselect_hyperslab( space, H5S_SELECT_SET, Offset, NULL, SlabSize, NULL );
                            status = H5Dwrite( DataSet, H5T_NATIVE_DOUBLE, MemSpace, space, H5P_DEFAULT, &U_ARR[0] );
                    }
                    status = H5Sclose( MemSpace );
                    status = H5Dclose( DataSet );
                    status = H5Sclose( space );

                    // Write Pfn_ED_MLT
                    SlabSize[0] = 1;
                    DataSet  = H5Dopen( file, "Pfn_ED_MLT", H5P_DEFAULT );
                    space    = H5Dget_space( DataSet );
                    MemSpace = H5Screate_simple( 1, SlabSize, NULL );
                    for ( iT=0; iT<H5_nT; iT++ ){
                            U_ARR[0] = H5_Pfn_cdmag[iT].x; U_ARR[1] = H5_Pfn_cdmag[iT].y; U_ARR[2] = H5_Pfn_cdmag[iT].z;
                            Offset[0]   = iT;
                            status = H5Sselect_hyperslab( space, H5S_SELECT_SET, Offset, NULL, SlabSize, NULL );
                            status = H5Dwrite( DataSet, H5T_NATIVE_DOUBLE, MemSpace, space, H5P_DEFAULT, &U_ARR[0] );
                    }
                    status = H5Sclose( MemSpace );
                    status = H5Dclose( DataSet );
                    status = H5Sclose( space );


                    // Write Bfn_geo
                    SlabSize[0] = 1; SlabSize[1] = 4;
                    DataSet  = H5Dopen( file, "Bfn_geo", H5P_DEFAULT );
                    space    = H5Dget_space( DataSet );
                    MemSpace = H5Screate_simple( 2, SlabSize, NULL );
                    for ( iT=0; iT<H5_nT; iT++ ){
                            U_ARR[0] = H5_Bfn_geo[iT].x; U_ARR[1] = H5_Bfn_geo[iT].y; U_ARR[2] = H5_Bfn_geo[iT].z; U_ARR[3] = Lgm_Magnitude( &H5_Bfn_geo[iT] );
                            Offset[0]   = iT; Offset[1] = 0;
                            status = H5Sselect_hyperslab( space, H5S_SELECT_SET, Offset, NULL, SlabSize, NULL );
                            status = H5Dwrite( DataSet, H5T_NATIVE_DOUBLE, MemSpace, space, H5P_DEFAULT, &U_ARR[0] );
                    }
                    status = H5Sclose( MemSpace );
                    status = H5Dclose( DataSet );
                    status = H5Sclose( space );

                    // Write Bfn_gsm
                    SlabSize[0] = 1; SlabSize[1] = 4;
                    DataSet  = H5Dopen( file, "Bfn_gsm", H5P_DEFAULT );
                    space    = H5Dget_space( DataSet );
                    MemSpace = H5Screate_simple( 2, SlabSize, NULL );
                    for ( iT=0; iT<H5_nT; iT++ ){
                            U_ARR[0] = H5_Bfn_gsm[iT].x; U_ARR[1] = H5_Bfn_gsm[iT].y; U_ARR[2] = H5_Bfn_gsm[iT].z; U_ARR[3] = Lgm_Magnitude( &H5_Bfn_gsm[iT] );
                            Offset[0]   = iT; Offset[1] = 0;
                            status = H5Sselect_hyperslab( space, H5S_SELECT_SET, Offset, NULL, SlabSize, NULL );
                            status = H5Dwrite( DataSet, H5T_NATIVE_DOUBLE, MemSpace, space, H5P_DEFAULT, &U_ARR[0] );
                    }
                    status = H5Sclose( MemSpace );
                    status = H5Dclose( DataSet );
                    status = H5Sclose( space );



                    // Write Loss_Cone_Alpha_n












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

                    // Write Pfs_geod_LatLon
                    SlabSize[0] = 1; SlabSize[1] = 2;
                    DataSet  = H5Dopen( file, "Pfs_geod_LatLon", H5P_DEFAULT );
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

                    // Write Pfs_geod_Height
                    SlabSize[0] = 1;
                    DataSet  = H5Dopen( file, "Pfs_geod_Height", H5P_DEFAULT );
                    space    = H5Dget_space( DataSet );
                    MemSpace = H5Screate_simple( 1, SlabSize, NULL );
                    for ( iT=0; iT<H5_nT; iT++ ){
                            U_ARR[0] = H5_Pfs_geod[iT].x; U_ARR[1] = H5_Pfs_geod[iT].y; U_ARR[2] = H5_Pfs_geod[iT].z;
                            Offset[0]   = iT;
                            status = H5Sselect_hyperslab( space, H5S_SELECT_SET, Offset, NULL, SlabSize, NULL );
                            status = H5Dwrite( DataSet, H5T_NATIVE_DOUBLE, MemSpace, space, H5P_DEFAULT, &U_ARR[2] );
                    }
                    status = H5Sclose( MemSpace );
                    status = H5Dclose( DataSet );
                    status = H5Sclose( space );

                    // Write Pfs_CD_MLAT
                    SlabSize[0] = 1;
                    DataSet  = H5Dopen( file, "Pfs_CD_MLAT", H5P_DEFAULT );
                    space    = H5Dget_space( DataSet );
                    MemSpace = H5Screate_simple( 1, SlabSize, NULL );
                    for ( iT=0; iT<H5_nT; iT++ ){
                            U_ARR[0] = H5_Pfs_cdmag[iT].x; U_ARR[1] = H5_Pfs_cdmag[iT].y; U_ARR[2] = H5_Pfs_cdmag[iT].z;
                            Offset[0]   = iT;
                            status = H5Sselect_hyperslab( space, H5S_SELECT_SET, Offset, NULL, SlabSize, NULL );
                            status = H5Dwrite( DataSet, H5T_NATIVE_DOUBLE, MemSpace, space, H5P_DEFAULT, &U_ARR[0] );
                    }
                    status = H5Sclose( MemSpace );
                    status = H5Dclose( DataSet );
                    status = H5Sclose( space );

                    // Write Pfs_CD_MLON
                    SlabSize[0] = 1;
                    DataSet  = H5Dopen( file, "Pfs_CD_MLON", H5P_DEFAULT );
                    space    = H5Dget_space( DataSet );
                    MemSpace = H5Screate_simple( 1, SlabSize, NULL );
                    for ( iT=0; iT<H5_nT; iT++ ){
                            U_ARR[0] = H5_Pfs_cdmag[iT].x; U_ARR[1] = H5_Pfs_cdmag[iT].y; U_ARR[2] = H5_Pfs_cdmag[iT].z;
                            Offset[0]   = iT;
                            status = H5Sselect_hyperslab( space, H5S_SELECT_SET, Offset, NULL, SlabSize, NULL );
                            status = H5Dwrite( DataSet, H5T_NATIVE_DOUBLE, MemSpace, space, H5P_DEFAULT, &U_ARR[0] );
                    }
                    status = H5Sclose( MemSpace );
                    status = H5Dclose( DataSet );
                    status = H5Sclose( space );

                    // Write Pfs_CD_MLT
                    SlabSize[0] = 1;
                    DataSet  = H5Dopen( file, "Pfs_CD_MLT", H5P_DEFAULT );
                    space    = H5Dget_space( DataSet );
                    MemSpace = H5Screate_simple( 1, SlabSize, NULL );
                    for ( iT=0; iT<H5_nT; iT++ ){
                            U_ARR[0] = H5_Pfs_cdmag[iT].x; U_ARR[1] = H5_Pfs_cdmag[iT].y; U_ARR[2] = H5_Pfs_cdmag[iT].z;
                            Offset[0]   = iT;
                            status = H5Sselect_hyperslab( space, H5S_SELECT_SET, Offset, NULL, SlabSize, NULL );
                            status = H5Dwrite( DataSet, H5T_NATIVE_DOUBLE, MemSpace, space, H5P_DEFAULT, &U_ARR[0] );
                    }
                    status = H5Sclose( MemSpace );
                    status = H5Dclose( DataSet );
                    status = H5Sclose( space );


                    // Write Pfs_ED_MLAT
                    SlabSize[0] = 1;
                    DataSet  = H5Dopen( file, "Pfs_ED_MLAT", H5P_DEFAULT );
                    space    = H5Dget_space( DataSet );
                    MemSpace = H5Screate_simple( 1, SlabSize, NULL );
                    for ( iT=0; iT<H5_nT; iT++ ){
                            U_ARR[0] = H5_Pfs_cdmag[iT].x; U_ARR[1] = H5_Pfs_cdmag[iT].y; U_ARR[2] = H5_Pfs_cdmag[iT].z;
                            Offset[0]   = iT;
                            status = H5Sselect_hyperslab( space, H5S_SELECT_SET, Offset, NULL, SlabSize, NULL );
                            status = H5Dwrite( DataSet, H5T_NATIVE_DOUBLE, MemSpace, space, H5P_DEFAULT, &U_ARR[0] );
                    }
                    status = H5Sclose( MemSpace );
                    status = H5Dclose( DataSet );
                    status = H5Sclose( space );

                    // Write Pfs_ED_MLON
                    SlabSize[0] = 1;
                    DataSet  = H5Dopen( file, "Pfs_ED_MLON", H5P_DEFAULT );
                    space    = H5Dget_space( DataSet );
                    MemSpace = H5Screate_simple( 1, SlabSize, NULL );
                    for ( iT=0; iT<H5_nT; iT++ ){
                            U_ARR[0] = H5_Pfs_cdmag[iT].x; U_ARR[1] = H5_Pfs_cdmag[iT].y; U_ARR[2] = H5_Pfs_cdmag[iT].z;
                            Offset[0]   = iT;
                            status = H5Sselect_hyperslab( space, H5S_SELECT_SET, Offset, NULL, SlabSize, NULL );
                            status = H5Dwrite( DataSet, H5T_NATIVE_DOUBLE, MemSpace, space, H5P_DEFAULT, &U_ARR[0] );
                    }
                    status = H5Sclose( MemSpace );
                    status = H5Dclose( DataSet );
                    status = H5Sclose( space );

                    // Write Pfs_ED_MLT
                    SlabSize[0] = 1;
                    DataSet  = H5Dopen( file, "Pfs_ED_MLT", H5P_DEFAULT );
                    space    = H5Dget_space( DataSet );
                    MemSpace = H5Screate_simple( 1, SlabSize, NULL );
                    for ( iT=0; iT<H5_nT; iT++ ){
                            U_ARR[0] = H5_Pfs_cdmag[iT].x; U_ARR[1] = H5_Pfs_cdmag[iT].y; U_ARR[2] = H5_Pfs_cdmag[iT].z;
                            Offset[0]   = iT;
                            status = H5Sselect_hyperslab( space, H5S_SELECT_SET, Offset, NULL, SlabSize, NULL );
                            status = H5Dwrite( DataSet, H5T_NATIVE_DOUBLE, MemSpace, space, H5P_DEFAULT, &U_ARR[0] );
                    }
                    status = H5Sclose( MemSpace );
                    status = H5Dclose( DataSet );
                    status = H5Sclose( space );


                    // Write Bfs_geo
                    SlabSize[0] = 1; SlabSize[1] = 4;
                    DataSet  = H5Dopen( file, "Bfs_geo", H5P_DEFAULT );
                    space    = H5Dget_space( DataSet );
                    MemSpace = H5Screate_simple( 2, SlabSize, NULL );
                    for ( iT=0; iT<H5_nT; iT++ ){
                            U_ARR[0] = H5_Bfs_geo[iT].x; U_ARR[1] = H5_Bfs_geo[iT].y; U_ARR[2] = H5_Bfs_geo[iT].z; U_ARR[3] = Lgm_Magnitude( &H5_Bfs_geo[iT] );
                            Offset[0]   = iT; Offset[1] = 0;
                            status = H5Sselect_hyperslab( space, H5S_SELECT_SET, Offset, NULL, SlabSize, NULL );
                            status = H5Dwrite( DataSet, H5T_NATIVE_DOUBLE, MemSpace, space, H5P_DEFAULT, &U_ARR[0] );
                    }
                    status = H5Sclose( MemSpace );
                    status = H5Dclose( DataSet );
                    status = H5Sclose( space );

                    // Write Bfs_gsm
                    SlabSize[0] = 1; SlabSize[1] = 4;
                    DataSet  = H5Dopen( file, "Bfs_gsm", H5P_DEFAULT );
                    space    = H5Dget_space( DataSet );
                    MemSpace = H5Screate_simple( 2, SlabSize, NULL );
                    for ( iT=0; iT<H5_nT; iT++ ){
                            U_ARR[0] = H5_Bfs_gsm[iT].x; U_ARR[1] = H5_Bfs_gsm[iT].y; U_ARR[2] = H5_Bfs_gsm[iT].z; U_ARR[3] = Lgm_Magnitude( &H5_Bfs_gsm[iT] );
                            Offset[0]   = iT; Offset[1] = 0;
                            status = H5Sselect_hyperslab( space, H5S_SELECT_SET, Offset, NULL, SlabSize, NULL );
                            status = H5Dwrite( DataSet, H5T_NATIVE_DOUBLE, MemSpace, space, H5P_DEFAULT, &U_ARR[0] );
                    }
                    status = H5Sclose( MemSpace );
                    status = H5Dclose( DataSet );
                    status = H5Sclose( space );



                    // Write Loss_Cone_Alpha_s












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

                } //end else
            } // end "if ( !FileExists || Force )" control structure


        } // end birds loop
    } // end JD loop


    //Lgm_DateTime_Destroy( myUTC );
    free( ShellFile );
    free( OutFile );
    free( HdfOutFile );
    free( InFile );

    Lgm_free_ctrans( c );
    Lgm_destroy_eop( e );
    Lgm_FreeMagEphemInfo( MagEphemInfo );

    LGM_ARRAY_1D_FREE( CmdLine );
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
    LGM_ARRAY_1D_FREE( H5_InOut );
    LGM_ARRAY_1D_FREE( H5_UTC );
    LGM_ARRAY_1D_FREE( H5_JD );
    LGM_ARRAY_1D_FREE( H5_GpsTime );
    LGM_ARRAY_1D_FREE( H5_TiltAngle );
    LGM_ARRAY_1D_FREE( H5_Kp );
    LGM_ARRAY_1D_FREE( H5_Dst );
    LGM_ARRAY_1D_FREE( H5_S_sc_to_pfn );
    LGM_ARRAY_1D_FREE( H5_S_sc_to_pfs );
    LGM_ARRAY_1D_FREE( H5_S_pfs_to_Bmin );
    LGM_ARRAY_1D_FREE( H5_S_Bmin_to_sc );
    LGM_ARRAY_1D_FREE( H5_S_total );
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
    LGM_ARRAY_2D_FREE( Perigee_IsoTimes );
    LGM_ARRAY_2D_FREE( Apogee_IsoTimes );
    LGM_ARRAY_2D_FREE( Perigee_Geod );
    LGM_ARRAY_2D_FREE( Apogee_Geod );
    LGM_ARRAY_1D_FREE( ApoPeriTimeList );



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
