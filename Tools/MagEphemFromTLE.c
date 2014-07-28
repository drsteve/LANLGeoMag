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
#include <Lgm/qsort.h>
#include <Lgm_Octree.h>

#define EARTH_ID     399
#define MOON_ID      301
#define RBSPA_ID    -362
#define RBSPB_ID    -363
#define  MAXIV      1000
#define  WINSIZ     ( 2 * MAXIV )

#define T89Q_KP  2.0
#define OP77Q_KP 2.0

/*
 * Stuff for inbound/outbound determination
 */
typedef struct TimeList {
    double  key;    // set to JD or somesuch thing.
    int     val;    // set to 0 for perigee time 1 for apogee time
} TimeList;

static void elt_qsort( struct TimeList *arr, unsigned n ) {
    #define elt_lt(a,b) ((a)->key < (b)->key)
    QSORT( struct TimeList, arr, n, elt_lt );
}

/*
 * This routine is for figuring out whether we are inbound or outbound at any given time...
 * Takers:
 *    List of JD times of apogee/perigee crossings and current JD.
 *    The TimeList structure contains a key (JD) and a value indicating apogee or perigee.
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

const  char *ProgramName = "MagEphemFromTLE";
const  char *argp_program_version     = "2.4.0";
const  char *argp_program_bug_address = "<mghenderson@lanl.gov>";
static char doc[] =
"Computes the magnetic ephemeris of a S/C from trajectories determined from Two-Line Element sets.\n\n"
"The input file is an ASCII TLE-set with epochs spanning the time interval of interest.\n\n"

"The output is written into daily files with filename given by the specified template. OutFile may contain time variables that will be substituted.  The available time variables are '%YYYY', '%MM', and '%DD' which correspond respectively to 4-digit year, 2-digit month (Jan is 01), and 2-digit day of month. \n\n"

"The %B variable will also get substituted by the list of birds given in the -b option.  Here is an example using time-variables.\n\n"

"\t./MagEphemFromTLE -S 20020901 -E 20020930 setup.ker \n"
"\t\t/home/jsmith/MagEphemStuff/%YYYY/%YYYY%MM%DD_1989-046_MagEphem.txt.  \n"
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
 *       ARG - Name of the option's argument (zero if there isn't one).
 *     FLAGS - OPTION_ARG_OPTIONAL, OPTION_ALIAS, OPTION_HIDDEN, OPTION_DOC,
 *             or OPTION_NO_USAGE
 */
static struct argp_option Options[] = {
    { 0, 0, 0, 0,   "Time Options:", 1},
    {"StartDateTime",   'S',    "yyyymmdd[Thh:mm:ss]",          0,        "Start date/time in ISO 8601 format. Seconds will be truncated to integers.", 0 },
    {"EndDateTime",     'E',    "yyyymmdd[Thh:mm:ss]",          0,        "End date/time in ISO 8601 format. Seconds will be truncated to integers.", 0 },
    {"ModelDateTime",   'T',    "yyyymmdd[Thh:mm:ss]",          0,        "Force the model to use parameters for the given date/time. This can be used to force the model to be static instead of dynamic. Date/time in ISO 8601 format .", 0 },
    {"Delta",           'D',    "delta",                        0,        "Time cadence in (integer) seconds.", 0 },

    { 0, 0, 0, 0,   "Internal Model Options:", 2},
    {"IntModel",        'i',    "model",                        0,        "Internal Magnetic Field Model to use. Can be CDIP, EDIP, IGRF. Default is IGRF.\n\n", 0},

    { 0, 0, 0, 0,   "External Model Options:", 3},
    {"ExtModel",        'e',    "model",                        0,        "External Magnetic Field Model to use. Can be OP77Q, T87Q, T89Q, T87D, T89D, T96, T02, T01S, TS04D, TS07D. Here, Q stands for Quiet, D for Dynamic.\n", 0},
    {"Kp",              'K',    "Kp",                           0,        "If set, force Kp to be this value. Use values like 0.7, 1.0, 1.3 for 1-, 1, 1+" },
    { 0, 0, 0, 0,   "Other Options:", 4},
    {"Birds",           'b',    "\"bird1, bird2, etc\"",      0,        "Birds (sats) to use. E.g., \"LANL-02A, 1989-046, POLAR\".", 0   },
    {"PitchAngles",     'p',    "\"start_pa, end_pa, npa\"",  0,        "Pitch angles to compute. Default is \"5.0, 90, 18\"." },
    {"FootPointHeight", 'f',    "height",                     0,        "Footpoint height in km. Default is 100km."                  },
    {"Quality",         'q',    "quality",                    0,        "Quality to use for L* calculations. Default is 3."      },
    {"nFLsInDriftShell",'n',    "nFLsInDriftShell",           0,        "Number of Field Lines to use in construction of drift shell. Use values in the range [6,240]. Default is 24." },
    {"UseEop",          'z',    0,                            0,        "Use Earth Orientation Parameters when computing ephemerii" },
    {"Coords",          'C',    "coord_system",               0,        "Coordinate system used in the input file. Can be: LATLONRAD, SM, GSM, GEI2000 or GSE. Default is LATLONRAD." },

    { 0, 0, 0, 0,   "Output Options:", 5},
    {"Force",           'F',    0,                            0,        "Overwrite output file even if it already exists" },
    {"Verbosity",       'v',    "verbosity",                  0,        "Verbosity level to use. (0-4)"      },
    {"DumpShellFiles",  'd',    0,                            0,        "Dump full binary shell files (for use in visualizing drift shells)." },
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
    int         nFLsInDriftShell;
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

    int         FixModelDateTime;
    char        ModelDateTimeStr[80];

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
        case 'T': // end date
            strncpy( arguments->ModelDateTimeStr, arg, 80 );
            arguments->FixModelDateTime = 1;
            break;
        case 'e': // external model
            strcpy( arguments->ExtModel, arg );
            break;
        case 'i': // internal model
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
        case 'n':
            arguments->nFLsInDriftShell = atoi( arg );
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
    _SgpInfo        *sgp;
    _SgpTLE         tle;
    Lgm_CTrans      *c;
} afInfo;
double ApogeeFunc( double T, double val, void *Info ){

    afInfo          *a;
    Lgm_CTrans      *c;
    Lgm_DateTime    UTC;
    Lgm_Vector      U, Uteme;
    double          R, tsince;

    a   = (afInfo *)Info;
    c   = a->c;

    Lgm_Make_UTC( a->Date, T/3600.0, &UTC, c );


    // do the propagation
    tsince = (UTC.JD - a->tle.JD)*1440.0;
    LgmSgp_SGP4( tsince, a->sgp );

    // Convert from TEME -> GEI2000
    Lgm_Set_Coord_Transforms( UTC.Date, UTC.Time, c );
    Uteme.x = a->sgp->X/WGS84_A; Uteme.y = a->sgp->Y/WGS84_A; Uteme.z = a->sgp->Z/WGS84_A;
    Lgm_Convert_Coords( &Uteme, &U, TEME_TO_GEI2000, c );
    R = Lgm_Magnitude( &U );

    return( R*a->Sgn );

}



/*
 * Returns latitude of S/C in GEO.
 */
typedef struct lfInfo {
    long int        Date;
    int             Sgn;
    int             BODY;
    _SgpInfo        *sgp;
    _SgpTLE         tle;
    Lgm_CTrans      *c;
} lfInfo;
double LatitudeFunc( double T, double val, void *Info ){

    lfInfo          *a;
    Lgm_CTrans      *c;
    Lgm_DateTime    UTC;
    Lgm_Vector      U, V, Uteme;
    double          R, tsince;

    a   = (lfInfo *)Info;
    c   = a->c;

    Lgm_Set_Coord_Transforms( a->Date, T/3600.0, c );
    Lgm_Make_UTC( a->Date, T/3600.0, &UTC, c );

    // do the propagation
    tsince = (UTC.JD - a->tle.JD)*1440.0;
    LgmSgp_SGP4( tsince, a->sgp );

    // Convert from TEME -> GEI2000
    Lgm_Set_Coord_Transforms( UTC.Date, UTC.Time, c );
    Uteme.x = a->sgp->X/WGS84_A; Uteme.y = a->sgp->Y/WGS84_A; Uteme.z = a->sgp->Z/WGS84_A;
    Lgm_Convert_Coords( &Uteme, &V, TEME_TO_GEI2000, c );
    R = Lgm_Magnitude( &V );
    return( asin( V.z/R )*DegPerRad );

}


int GetOrbitNumber( Lgm_DateTime *UTC, int nPerigee, Lgm_DateTime *Perigee_UTC, int *PerigeeOrbitNumber ){

    // this will not work in general....
    // just return -999
    return( -999 );

}

void Lgm_WriteMagEphemDataKML( FILE *fp, int i, Lgm_MagEphemData *med ) {

    int j;

    fprintf( fp, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n" );
    fprintf( fp, "<kml xmlns=\"http://www.opengis.net/kml/2.2\" xmlns:gx=\"http://www.google.com/kml/ext/2.2\">\n" );
    fprintf( fp, "  <Document>\n" );
    fprintf( fp, "    <name>RBSP A</name>\n" );
    fprintf( fp, "    <Snippet>Created Wed Jun 2 15:33:39 2010</Snippet>\n" );
    fprintf( fp, "\n" );
    fprintf( fp, "    <!-- Normal track style -->\n" );
    fprintf( fp, "    <LookAt>\n" );
    fprintf( fp, "      <gx:TimeSpan>\n" );
    fprintf( fp, "        <begin>%s</begin>\n", med->H5_IsoTimes[0] );
    fprintf( fp, "        <end>%s</end>\n", med->H5_IsoTimes[med->H5_nT-1] );
    fprintf( fp, "      </gx:TimeSpan>\n" );
    fprintf( fp, "      <longitude>%.10g</longitude>\n", med->H5_Pfs_geod_LatLon[0][1] );
    fprintf( fp, "      <latitude>%.10g</latitude>\n", med->H5_Pfs_geod_LatLon[0][0] );
    //fprintf( fp, "      <range>%.10g</range>\n", med->H5_Pfs_geod_Height[0] );
    fprintf( fp, "      <range>%.10g</range>\n", 20000*1e3 );
    fprintf( fp, "    </LookAt>\n" );
    fprintf( fp, "    <Style id=\"track_n\">\n" );
    fprintf( fp, "      <IconStyle>\n" );
    fprintf( fp, "        <scale>.5</scale>\n" );
    fprintf( fp, "        <Icon>\n" );
    fprintf( fp, "          <href>http://earth.google.com/images/kml-icons/track-directional/track-none.png</href>\n" );
    fprintf( fp, "        </Icon>\n" );
    fprintf( fp, "      </IconStyle>\n" );
    fprintf( fp, "      <LabelStyle>\n" );
    fprintf( fp, "        <scale>0</scale>\n" );
    fprintf( fp, "      </LabelStyle>\n" );
    fprintf( fp, "\n" );
    fprintf( fp, "    </Style>\n" );
    fprintf( fp, "    <!-- Highlighted track style -->\n" );
    fprintf( fp, "    <Style id=\"track_h\">\n" );
    fprintf( fp, "      <IconStyle>\n" );
    fprintf( fp, "        <scale>1.2</scale>\n" );
    fprintf( fp, "        <Icon>\n" );
    fprintf( fp, "          <href>http://earth.google.com/images/kml-icons/track-directional/track-none.png</href>\n" );
    fprintf( fp, "        </Icon>\n" );
    fprintf( fp, "      </IconStyle>\n" );
    fprintf( fp, "    </Style>\n" );
    fprintf( fp, "    <StyleMap id=\"track\">\n" );
    fprintf( fp, "      <Pair>\n" );
    fprintf( fp, "        <key>normal</key>\n" );
    fprintf( fp, "        <styleUrl>#track_n</styleUrl>\n" );
    fprintf( fp, "      </Pair>\n" );
    fprintf( fp, "      <Pair>\n" );
    fprintf( fp, "        <key>highlight</key>\n" );
    fprintf( fp, "        <styleUrl>#track_h</styleUrl>\n" );
    fprintf( fp, "      </Pair>\n" );
    fprintf( fp, "    </StyleMap>\n" );
    fprintf( fp, "    <!-- Normal multiTrack style -->\n" );
    fprintf( fp, "    <Style id=\"multiTrack_n\">\n" );
    fprintf( fp, "      <IconStyle>\n" );
    fprintf( fp, "        <Icon>\n" );
    fprintf( fp, "          <href>http://earth.google.com/images/kml-icons/track-directional/track-0.png</href>\n" );
    fprintf( fp, "        </Icon>\n" );
    fprintf( fp, "      </IconStyle>\n" );
    fprintf( fp, "      <LineStyle>\n" );
    fprintf( fp, "        <color>99ffac59</color>\n" );
    fprintf( fp, "        <width>6</width>\n" );
    fprintf( fp, "      </LineStyle>\n" );
    fprintf( fp, "\n" );
    fprintf( fp, "    </Style>\n" );
    fprintf( fp, "    <!-- Highlighted multiTrack style -->\n" );
    fprintf( fp, "    <Style id=\"multiTrack_h\">\n" );
    fprintf( fp, "      <IconStyle>\n" );
    fprintf( fp, "        <scale>1.2</scale>\n" );
    fprintf( fp, "        <Icon>\n" );
    fprintf( fp, "          <href>http://earth.google.com/images/kml-icons/track-directional/track-0.png</href>\n" );
    fprintf( fp, "        </Icon>\n" );
    fprintf( fp, "      </IconStyle>\n" );
    fprintf( fp, "      <LineStyle>\n" );
    fprintf( fp, "        <color>99ffac59</color>\n" );
    fprintf( fp, "        <width>8</width>\n" );
    fprintf( fp, "      </LineStyle>\n" );
    fprintf( fp, "    </Style>\n" );
    fprintf( fp, "    <StyleMap id=\"multiTrack\">\n" );
    fprintf( fp, "      <Pair>\n" );
    fprintf( fp, "        <key>normal</key>\n" );
    fprintf( fp, "        <styleUrl>#multiTrack_n</styleUrl>\n" );
    fprintf( fp, "      </Pair>\n" );
    fprintf( fp, "      <Pair>\n" );
    fprintf( fp, "        <key>highlight</key>\n" );
    fprintf( fp, "        <styleUrl>#multiTrack_h</styleUrl>\n" );
    fprintf( fp, "      </Pair>\n" );
    fprintf( fp, "    </StyleMap>\n" );
    fprintf( fp, "    <!-- Normal waypoint style -->\n" );
    fprintf( fp, "    <Style id=\"waypoint_n\">\n" );
    fprintf( fp, "      <IconStyle>\n" );
    fprintf( fp, "        <Icon>\n" );
    fprintf( fp, "          <href>http://maps.google.com/mapfiles/kml/pal4/icon61.png</href>\n" );
    fprintf( fp, "        </Icon>\n" );
    fprintf( fp, "      </IconStyle>\n" );
    fprintf( fp, "    </Style>\n" );
    fprintf( fp, "    <!-- Highlighted waypoint style -->\n" );
    fprintf( fp, "    <Style id=\"waypoint_h\">\n" );
    fprintf( fp, "      <IconStyle>\n" );
    fprintf( fp, "        <scale>1.2</scale>\n" );
    fprintf( fp, "        <Icon>\n" );
    fprintf( fp, "          <href>http://maps.google.com/mapfiles/kml/pal4/icon61.png</href>\n" );
    fprintf( fp, "        </Icon>\n" );
    fprintf( fp, "      </IconStyle>\n" );
    fprintf( fp, "    </Style>\n" );
    fprintf( fp, "    <StyleMap id=\"waypoint\">\n" );
    fprintf( fp, "      <Pair>\n" );
    fprintf( fp, "        <key>normal</key>\n" );
    fprintf( fp, "        <styleUrl>#waypoint_n</styleUrl>\n" );
    fprintf( fp, "      </Pair>\n" );
    fprintf( fp, "      <Pair>\n" );
    fprintf( fp, "        <key>highlight</key>\n" );
    fprintf( fp, "        <styleUrl>#waypoint_h</styleUrl>\n" );
    fprintf( fp, "      </Pair>\n" );
    fprintf( fp, "    </StyleMap>\n" );
    fprintf( fp, "    <Style id=\"lineStyle\">\n" );
    fprintf( fp, "      <LineStyle>\n" );
    fprintf( fp, "        <color>99ffac59</color>\n" );
    fprintf( fp, "        <width>6</width>\n" );
    fprintf( fp, "      </LineStyle>\n" );
    fprintf( fp, "    </Style>\n" );
    fprintf( fp, "    <Schema id=\"schema\">\n" );
    fprintf( fp, "      <gx:SimpleArrayField name=\"heartrate\" type=\"int\">\n" );
    fprintf( fp, "        <displayName>Heart Rate</displayName>\n" );
    fprintf( fp, "      </gx:SimpleArrayField>\n" );
    fprintf( fp, "      <gx:SimpleArrayField name=\"cadence\" type=\"int\">\n" );
    fprintf( fp, "        <displayName>Cadence</displayName>\n" );
    fprintf( fp, "      </gx:SimpleArrayField>\n" );
    fprintf( fp, "      <gx:SimpleArrayField name=\"power\" type=\"float\">\n" );
    fprintf( fp, "        <displayName>Power</displayName>\n" );
    fprintf( fp, "      </gx:SimpleArrayField>\n" );
    fprintf( fp, "    </Schema>\n" );
    fprintf( fp, "    <Folder>\n" );
    fprintf( fp, "      <name>Tracks</name>\n" );

    fprintf( fp, "      <Placemark>\n" );
    fprintf( fp, "        <name>%s</name>\n", "Southern Footpoint" );
    fprintf( fp, "        <styleUrl>#multiTrack</styleUrl>\n" );
    fprintf( fp, "        <gx:Track>\n" );
    fprintf( fp, "          <altitudeMode>absolute</altitudeMode>\n");
    fprintf( fp, "          <gx:interpolate>1</gx:interpolate>\n");
    for (j=0; j<med->H5_nT; j++) fprintf( fp, "          <when>%s</when>\n", med->H5_IsoTimes[j] );
    for (j=0; j<med->H5_nT; j++) fprintf( fp, "          <gx:coord>%.10g %.10g %.10g</gx:coord>\n", med->H5_Pfs_geod_LatLon[j][1], med->H5_Pfs_geod_LatLon[j][0], med->H5_Pfs_geod_Height[j]*1e3 );
/*
    fprintf( fp, "          <ExtendedData>\n" );
    fprintf( fp, "            <SchemaData schemaUrl=\"#schema\">\n" );
    fprintf( fp, "              <gx:SimpleArrayData name=\"cadence\">\n" );
    fprintf( fp, "                <gx:value>86</gx:value>\n" );
    fprintf( fp, "                <gx:value>103</gx:value>\n" );
    fprintf( fp, "                <gx:value>108</gx:value>\n" );
    fprintf( fp, "                <gx:value>113</gx:value>\n" );
    fprintf( fp, "                <gx:value>113</gx:value>\n" );
    fprintf( fp, "                <gx:value>113</gx:value>\n" );
    fprintf( fp, "                <gx:value>113</gx:value>\n" );
    fprintf( fp, "              </gx:SimpleArrayData>\n" );
    fprintf( fp, "              <gx:SimpleArrayData name=\"heartrate\">\n" );
    fprintf( fp, "                <gx:value>181</gx:value>\n" );
    fprintf( fp, "                <gx:value>177</gx:value>\n" );
    fprintf( fp, "                <gx:value>175</gx:value>\n" );
    fprintf( fp, "                <gx:value>173</gx:value>\n" );
    fprintf( fp, "                <gx:value>173</gx:value>\n" );
    fprintf( fp, "                <gx:value>173</gx:value>\n" );
    fprintf( fp, "                <gx:value>173</gx:value>\n" );
    fprintf( fp, "              </gx:SimpleArrayData>\n" );
    fprintf( fp, "              <gx:SimpleArrayData name=\"power\">\n" );
    fprintf( fp, "                <gx:value>327.0</gx:value>\n" );
    fprintf( fp, "                <gx:value>177.0</gx:value>\n" );
    fprintf( fp, "                <gx:value>179.0</gx:value>\n" );
    fprintf( fp, "                <gx:value>162.0</gx:value>\n" );
    fprintf( fp, "                <gx:value>166.0</gx:value>\n" );
    fprintf( fp, "                <gx:value>177.0</gx:value>\n" );
    fprintf( fp, "                <gx:value>183.0</gx:value>\n" );
    fprintf( fp, "              </gx:SimpleArrayData>\n" );
    fprintf( fp, "            </SchemaData>\n" );
    fprintf( fp, "          </ExtendedData>\n" );
*/
    fprintf( fp, "        </gx:Track>\n" );
    fprintf( fp, "      </Placemark>\n" );

    fprintf( fp, "      <Placemark>\n" );
    fprintf( fp, "        <name>%s</name>\n", "Northern Footpoint" );
    fprintf( fp, "        <styleUrl>#multiTrack</styleUrl>\n" );
    fprintf( fp, "        <gx:Track>\n" );
    fprintf( fp, "          <altitudeMode>absolute</altitudeMode>\n");
    fprintf( fp, "          <gx:interpolate>1</gx:interpolate>\n");
    for (j=0; j<med->H5_nT; j++) fprintf( fp, "          <when>%s</when>\n", med->H5_IsoTimes[j] );
    for (j=0; j<med->H5_nT; j++) fprintf( fp, "          <gx:coord>%.10g %.10g %.10g</gx:coord>\n", med->H5_Pfn_geod_LatLon[j][1], med->H5_Pfn_geod_LatLon[j][0], med->H5_Pfn_geod_Height[j] );
    fprintf( fp, "        </gx:Track>\n" );
    fprintf( fp, "      </Placemark>\n" );


    fprintf( fp, "    </Folder>\n" );
    fprintf( fp, "  </Document>\n" );
    fprintf( fp, "</kml>\n" );


}



/*
 *  Compute magnetic ephemerii of S/C.
 */
int main( int argc, char *argv[] ){

    Lgm_ElapsedTimeInfo t;
    long int         IdNumber;
    char             IntDesig[10], CommonName[80];
    double           et, pos[3];

    struct Arguments arguments;
    Lgm_CTrans       *c = Lgm_init_ctrans( 0 );
    Lgm_Vector       Rgsm, Rgeo, W, U;
    Lgm_DateTime     UTC, ModelDateTime;
    Lgm_Eop          *e = Lgm_init_eop( 0 );
    Lgm_EopOne       eop;
    int              FixModelDateTime;
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
    int              nAlpha, Quality, nFLsInDriftShell, Verbosity;
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
    double          ForceKp;
    double          GeodLat, GeodLong, GeodHeight, MLAT, MLON, MLT;
    Lgm_MagEphemData *med;


    lfInfo          *lfi;
    double          La, Lb, Lmin;
    int             done, BODY;
    long int        ss, es, Seconds, Ta, Tb, Tc;
    double          R, Ra, Rb, Rc, Rmin, Tmin;
    BrentFuncInfo   bInfo;
    afInfo          *afi;
    int             nPerigee, nApogee, nAscend;
    int             ApogeeOrbitNumber[4];
    int             PerigeeOrbitNumber[4];
    Lgm_DateTime    *Perigee_UTC;
    Lgm_DateTime    *Apogee_UTC;
    Lgm_DateTime    *Ascend_UTC;
    Lgm_Vector      *Perigee_U;
    Lgm_Vector      *Apogee_U;
    Lgm_Vector      *Ascend_U;
    TimeList        *ApoPeriTimeList;
    int             nApoPeriTimeList;
    Lgm_Vector      Bvec, Bvec2, w, Bsc_gsm;
    double          Bsc_mag, Bfn_mag, Bfs_mag;
    double          Bmin_mag, s, cl, Ek, E, pp, p2c2, Beta2, Beta, vel, T, rg;
    int             n, OverRideKp;
    char            *CmdLine, TmpStr[2048];

    double          tsince;
    Lgm_Vector      Uteme, Ugei;
    int             nTle=0; // this is not the number of TLE but the index in the tle array
    _SgpTLE         *tle;  // pointer to a struct
    _SgpInfo        *sgp;




    t.ColorizeText = TRUE;
    Lgm_ElapsedTimeInit( &t, 255, 150, 0 );




    LGM_ARRAY_1D( ApoPeriTimeList, 20, TimeList );
    LGM_ARRAY_1D( Perigee_UTC, 10, Lgm_DateTime );
    LGM_ARRAY_1D( Apogee_UTC,  10, Lgm_DateTime );
    LGM_ARRAY_1D( Ascend_UTC,  10, Lgm_DateTime );
    LGM_ARRAY_1D( Perigee_U, 10, Lgm_Vector );
    LGM_ARRAY_1D( Apogee_U,  10, Lgm_Vector );
    LGM_ARRAY_1D( Ascend_U,  10, Lgm_Vector );

    /*
     * Default option values.
     */
    arguments.StartPA          = 90;     // start at 90.0 Deg.
    arguments.EndPA            = 5.0;    // stop at 2.5 Deg.
    arguments.nPA              = 18;     // 18 pitch angles
    arguments.silent           = 0;
    arguments.Verbosity        = 0;
    arguments.Quality          = 3;
    arguments.nFLsInDriftShell = 24;
    arguments.Kp               = -999.9;
    arguments.Delta            = 60;     // 60s default cadence
    arguments.Colorize         = 0;
    arguments.Force            = 0;
    arguments.UseEop           = 0;
    arguments.DumpShellFiles   = 0;
    arguments.FixModelDateTime =  0;
    arguments.StartDate        = -1;
    arguments.EndDate          = -1;
    arguments.FootPointHeight  = 100.0; // km
    strcpy( arguments.IntModel, "IGRF" );
    strcpy( arguments.ExtModel, "T89D" );
    strcpy( arguments.CoordSystem, "LATLONRAD" );
    strcpy( arguments.ModelDateTimeStr,  "" );


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
    FootpointHeight  = arguments.FootPointHeight;
    Quality          = arguments.Quality;
    nFLsInDriftShell = arguments.nFLsInDriftShell;
    ForceKp          = arguments.Kp;
    Verbosity        = arguments.Verbosity;
    Colorize         = arguments.Colorize;
    Force            = arguments.Force;
    UseEop           = arguments.UseEop;
    DumpShellFiles   = arguments.DumpShellFiles;
    Delta            = arguments.Delta;
    StartDate        = arguments.StartDate;
    EndDate          = arguments.EndDate;
    StartSeconds     = arguments.StartSeconds;
    EndSeconds       = arguments.EndSeconds;
    FixModelDateTime = arguments.FixModelDateTime;
    strcpy( IntModel,  arguments.IntModel );
    strcpy( ExtModel,  arguments.ExtModel );
    strcpy( CoordSystem,  arguments.CoordSystem );
    if ( FixModelDateTime ) IsoTimeStringToDateTime( arguments.ModelDateTimeStr, &ModelDateTime, c );

    LGM_ARRAY_2D( Birds, 20, 80, char );
    StringSplit( arguments.Birds, Birds, 80, &nBirds );







    /*
     *  Print summary of our options
     */
    if ( !arguments.silent ) {
        printf("\n\n");
        printf( "\t               Program/Version: %s\n", argp_program_version );
        printf( "\t                Bug Reports To: %s\n", argp_program_bug_address );
        printf( "\t                Input TLE File: %s\n", InputFilename );
        printf( "\t                   Output File: %s\n", OutputFilename );
        printf( "\t            Input Coord System: %s\n", CoordSystem );
        printf( "\t        Number of Pitch Angles: %d\n", nAlpha );
        printf( "\t        Pitch Angles [Degrees]:" );
        for (i=0; i<nAlpha; i++ ) printf( " %g", Alpha[i] );
        printf("\n");
        printf( "\t                Internal Model: %s\n", IntModel );
        printf( "\t                External Model: %s\n", ExtModel );
        if ( FixModelDateTime ) {
            printf( "\t     Model Parameters Fixed to: %s  ( ", arguments.ModelDateTimeStr ); Lgm_Print_DateTime( &ModelDateTime, 4, 8 );  printf( " )\n" );
        }
        if ( ForceKp >= 0.0 ) {
            printf( "\t              Forcing Kp to be: %g\n", ForceKp );
        }
        printf( "\t          FootpointHeight [km]: %g\n", FootpointHeight );
        printf( "\t                    L* Quality: %d\n", Quality );
        printf( "\t       Num. FLs in drift shell: %d\n", nFLsInDriftShell );
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

// PROBLEM AREA?
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

    Lgm_Doy( StartDate, &sYear, &sMonth, &sDay, &sDoy);
    sJD = Lgm_JD( sYear, sMonth, sDay, 12.0, LGM_TIME_SYS_UTC, c );
    Lgm_Doy( EndDate, &eYear, &eMonth, &eDay, &eDoy);
    eJD = Lgm_JD( eYear, eMonth, eDay, 12.0, LGM_TIME_SYS_UTC, c );

    
    int NNN = (int)( (eJD-sJD+1)*86400.0/(double)Delta + 1.0);
    if (NNN > 86401 ) NNN = 86401;
    med = Lgm_InitMagEphemData( NNN, 80 );


    if ( nAlpha > 0 ){
        MagEphemInfo = Lgm_InitMagEphemInfo(0, nAlpha);
    } else {
        // doesn't seem to like allocating zero size...
        MagEphemInfo = Lgm_InitMagEphemInfo(0, 1);
    }
    MagEphemInfo->PropagatorType = SGP4;


    // Settings for Lstar calcs
    Lgm_SetMagEphemLstarQuality( Quality, nFLsInDriftShell, MagEphemInfo );
    MagEphemInfo->SaveShellLines = TRUE;
    MagEphemInfo->LstarInfo->LSimpleMax = 12.0;
    MagEphemInfo->LstarInfo->VerbosityLevel = Verbosity;
    MagEphemInfo->LstarInfo->mInfo->VerbosityLevel = Verbosity;
    MagEphemInfo->LstarInfo->mInfo->Lgm_LossConeHeight = FootpointHeight;

    // The default model is Lgm_B_T89c
    MagEphemInfo->LstarInfo->mInfo->Bfield = Lgm_B_T89c;

    // Now override the default if user tells us to.
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
        MagEphemInfo->LstarInfo->mInfo->Bfield = Lgm_B_OP77;
        ForceKp    = OP77Q_KP;
        OverRideKp = TRUE;
    } else if ( !strcmp( ExtModel, "T89Q" ) ){
        MagEphemInfo->LstarInfo->mInfo->Bfield = Lgm_B_T89c;
        ForceKp    = T89Q_KP;
        OverRideKp = TRUE;
    } else if ( !strcmp( ExtModel, "T89D" ) ){
        MagEphemInfo->LstarInfo->mInfo->Bfield = Lgm_B_T89c;
    } else if ( !strcmp( ExtModel, "T96" ) ){
        MagEphemInfo->LstarInfo->mInfo->Bfield = Lgm_B_T96;
    } else if ( !strcmp( ExtModel, "T01S" ) ){
        MagEphemInfo->LstarInfo->mInfo->Bfield = Lgm_B_T01S;
    } else if ( !strcmp( ExtModel, "T02" ) ){
        MagEphemInfo->LstarInfo->mInfo->Bfield = Lgm_B_T02;
    } else if ( !strcmp( ExtModel, "TS04D" ) ){
        MagEphemInfo->LstarInfo->mInfo->Bfield = Lgm_B_TS04;
    } else if ( !strcmp( ExtModel, "TS07D" ) ){
//        MagEphemInfo->LstarInfo->mInfo->Bfield = Lgm_B_TS07;
//Lgm_SetCoeffs_TS07( 0, 0, &(MagEphemInfo->LstarInfo->mInfo->TS07_Info) );
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

    med->H5_nAlpha = nAlpha;
    for (i=0; i<nAlpha; i++) {
        MagEphemInfo->Alpha[i]    = Alpha[i];
        med->H5_Alpha[i] = Alpha[i];
    }



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


    // create a place to hold an array of TLEs (1 element)
    tle = (_SgpTLE *)calloc( 1000, sizeof(_SgpTLE) );

    sgp = (_SgpInfo *)calloc( 1, sizeof(_SgpInfo) );



    for ( JD = sJD; JD <= eJD; JD += 1.0 ) {

        Date = Lgm_JD_to_Date( JD, &Year, &Month, &Day, &Time );



        /*
         * loop over all birds
         */
         for ( iBird = 0; iBird < nBirds; iBird++ ) {


            strcpy( InFile, InputFilename );
            strcpy( OutFile, OutputFilename );
            strcpy( Bird, Birds[ iBird ] );

// PROBLEM AREA?
if        ( !strcmp( Bird, "rbspa" ) ){
    BODY     = RBSPA_ID;
    IdNumber = 38752;
    strcpy( IntDesig,   "2012-046A" );
    strcpy( CommonName, "RBSP A" );
} else if ( !strcmp( Bird, "rbspb" ) ){
    BODY     = RBSPB_ID;
    IdNumber = 38753;
    strcpy( IntDesig, "2012-046B" );
    strcpy( CommonName, "RBSP A" );
} else {
    // maybe user just has a SPICE BODY number?
    IdNumber = -999999;
    sprintf( IntDesig, "SPICE Object: %s", Bird );
    sprintf( CommonName, "%s", Bird );
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

            /*
             * Assume that the Outfile name specified on the command line is the basename (i.e. no extensions on the end..)
             * Then we need to add ".txt" and ".h5" to each
             */
            sprintf( HdfOutFile, "%s.h5", OutFile );
            strcat( OutFile, ".txt" );



            if ( SubstituteVars ) {
                printf( "\n\n      Processing Date: %4d-%02d-%02d    Bird: %s\n", Year, Month, Day, Bird );
            } else {
                printf( "\n\n      Processing File: %s\n", InFile );
            }
            printf( "      -------------------------------------------------------------------------------------------\n");
            printf( "           Input File: %s\n", InFile);
            printf( "      TXT Output File: %s\n", OutFile);
            printf( "     HDF5 Output File: %s\n\n", HdfOutFile);


            // Create Base directory if it hasn't been created yet.
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

                    printf("\tCould not open TLE File ( %s ) for reading\n", InFile );

                } else {

                    fclose( fp_in );





// WE NEED TO free this list each time
                    /*
                     * Read TLEs from file
                     */
                    if (!LgmSgp_ReadTlesFromFile( InFile, &nTle, tle, 4)){
                        printf("TLE not parsed!\n");
                        free(tle);
                        return(-1);
                    } else {
                        printf("TLE read and parsed.\n");
                    }
                    LgmSgp_SortListOfTLEs( nTle, tle );


                    printf("Line0: %s\n", tle[10].Line0);
                    printf("Line1: %s\n", tle[10].Line1);
                    printf("Line2: %s\n", tle[10].Line2);
                    printf("\n");



                    LgmSgp_SGP4_Init( sgp, &tle[10] );







                    /*
                     * Find Apogees.
                     * Pack aInfo structure with required info.
                     * Find brackets.
                     * Call Lgm_Brent() to get minimum.
                     *
                     */
                    afi = (afInfo *)calloc( 1, sizeof(afInfo) );
                    afi->Date = Date;



// PROBLEM AREA?
afi->BODY = BODY;
                    // initialize the propagator
                    afi->sgp = (_SgpInfo *)calloc( 1, sizeof(_SgpInfo) );
                    afi->tle = tle[10];
                    LgmSgp_SGP4_Init( afi->sgp, &(afi->tle) );
                    

                    afi->Sgn  = -1; // +1 => find min (perigee)   -1 => find max (apogee)
                    afi->c    = c;
                    bInfo.Val  = 0.0;
                    bInfo.Info = (void *)afi;
                    bInfo.func = ApogeeFunc;
                    Ta = -900; Ra = bInfo.func( (double)Ta, 0.0, (void *)afi );
                    Tb =    0; Rb = bInfo.func( (double)Tb, 0.0, (void *)afi );
                    Tc =  900; Rc = bInfo.func( (double)Tc, 0.0, (void *)afi );
                    done = FALSE;
                    nApogee = 0;
                    nApoPeriTimeList = 0;
                    while ( !done ) {
                        if ( (Rb < Ra) && (Rb < Rc) ) {
                            // we have a bracket. Find Min.
                            Lgm_Brent( (double)Ta, (double)Tb, (double)Tc, &bInfo, 1e-10, &Tmin, &Rmin );
                            if ( (Tmin >=0) && (Tmin < 86400) ) {
                                Lgm_Make_UTC( Date, Tmin/3600.0, &Apogee_UTC[nApogee], c );


                                // do the propagation
                                tsince = (Apogee_UTC[nApogee].JD - tle[10].JD)*1440.0;
                                LgmSgp_SGP4( tsince, afi->sgp );

                                // Convert from TEME -> GEI2000
                                Uteme.x = afi->sgp->X/WGS84_A; Uteme.y = afi->sgp->Y/WGS84_A; Uteme.z = afi->sgp->Z/WGS84_A; 
                                Lgm_Set_Coord_Transforms( Apogee_UTC[nApogee].Date, Apogee_UTC[nApogee].Time, c );
                                Lgm_Convert_Coords( &Uteme, &Ugei, TEME_TO_GEI2000, c );
                                Apogee_U[nApogee].x = Ugei.x; Apogee_U[nApogee].y = Ugei.y; Apogee_U[nApogee].z = Ugei.z;
                                Lgm_Convert_Coords( &Apogee_U[nApogee], &w, GEI2000_TO_WGS84, c );
                                Lgm_WGS84_to_GEOD( &w, &med->H5_Apogee_Geod[nApogee][0], &med->H5_Apogee_Geod[nApogee][1], &med->H5_Apogee_Geod[nApogee][2] );

                                Lgm_DateTimeToString( med->H5_Apogee_IsoTimes[nApogee], &Apogee_UTC[nApogee], 0, 3 );
                                printf("nApogee: %d Bracket: T = %8ld %8ld %8ld    R = %g %g %g    Tapogee, Rapogee = %s %g\n", nApogee, Ta, Tb, Tc, Ra, Rb, Rc, med->H5_Apogee_IsoTimes[nApogee], fabs(Rmin) );
                                ApoPeriTimeList[nApoPeriTimeList].key = Apogee_UTC[nApogee].JD;
                                ApoPeriTimeList[nApoPeriTimeList].val = 1; // because its apogee
                                ++nApoPeriTimeList;

                                ++nApogee;
                            }

                        }
                        // advance another 900 seconds.
                        Ta = Tb; Ra = Rb;
                        Tb = Tc; Rb = Rc;
                        Tc += 900; Rc = bInfo.func( (double)Tc, 0.0, (void *)afi );
                        if (Tc > 86400+900) done = TRUE;
                    }
                    med->H5_nApogee  = nApogee;


                    


                    /*
                     * Find Perigees. Same as above except set afi->Sgn = 1.
                     */
                    printf("\n");
                    afi->Sgn  = 1; // +1 => find min (perigee)   -1 => find max (apogee)
                    Ta = -900; Ra = bInfo.func( (double)Ta, 0.0, (void *)afi );
                    Tb =    0; Rb = bInfo.func( (double)Tb, 0.0, (void *)afi );
                    Tc =  900; Rc = bInfo.func( (double)Tc, 0.0, (void *)afi );
                    done = FALSE;
                    nPerigee = 0;
                    while ( !done ) {
                        if ( (Rb < Ra) && (Rb < Rc) ) {
                            // we have a bracket. Find Min.
                            Lgm_Brent( (double)Ta, (double)Tb, (double)Tc, &bInfo, 1e-10, &Tmin, &Rmin );

                            if ( (Tmin >=0) && (Tmin < 86400) ) {
                                Lgm_Make_UTC( Date, Tmin/3600.0, &Perigee_UTC[nPerigee], c );

                                // do the propagation
                                tsince = (Perigee_UTC[nApogee].JD - tle[10].JD)*1440.0;
                                LgmSgp_SGP4( tsince, afi->sgp );

                                // Convert from TEME -> GEI2000
                                Uteme.x = afi->sgp->X/WGS84_A; Uteme.y = afi->sgp->Y/WGS84_A; Uteme.z = afi->sgp->Z/WGS84_A; 
                                Lgm_Set_Coord_Transforms( Perigee_UTC[nPerigee].Date, Perigee_UTC[nPerigee].Time, c );
                                Lgm_Convert_Coords( &Uteme, &Ugei, TEME_TO_GEI2000, c );
                                Perigee_U[nPerigee].x = Ugei.x; Perigee_U[nPerigee].y = Ugei.y; Perigee_U[nPerigee].z = Ugei.z;

                                Lgm_DateTimeToString( med->H5_Perigee_IsoTimes[nPerigee], &Perigee_UTC[nPerigee], 0, 3 );
                                printf("nPerigee: %d Bracket: T = %8ld %8ld %8ld    R = %g %g %g    Tperigee, Rperigee = %s %g\n", nPerigee, Ta, Tb, Tc, Ra, Rb, Rc, med->H5_Perigee_IsoTimes[nPerigee], fabs(Rmin) );

                                ApoPeriTimeList[nApoPeriTimeList].key = Perigee_UTC[nPerigee].JD;
                                ApoPeriTimeList[nApoPeriTimeList].val = 0; // because its perigee
                                ++nApoPeriTimeList;

                                ++nPerigee;
                            }
                        }
                        // advance another 900 seconds.
                        Ta = Tb; Ra = Rb;
                        Tb = Tc; Rb = Rc;
                        Tc += 900; Rc = bInfo.func( (double)Tc, 0.0, (void *)afi );
                        if (Tc > 86400+900) done = TRUE;
                    }
                    free( afi->sgp );
                    free( afi );
                    med->H5_nPerigee = nPerigee;
                    printf("\n");

                    /*
                     * sort the merged list of apogee/perigee times.
                     */
                    elt_qsort( ApoPeriTimeList, nApoPeriTimeList );




                    /*
                     * Find Times when S/C crosses nodes.
                     */

                    lfi = (lfInfo *)calloc( 1, sizeof(lfInfo) );

                    lfi->sgp = (_SgpInfo *)calloc( 1, sizeof(_SgpInfo) );
                    lfi->tle = tle[10];
                    LgmSgp_SGP4_Init( lfi->sgp, &(lfi->tle) );

                    lfi->Date = Date;
                    lfi->BODY = BODY;
                    lfi->c    = c;
                    bInfo.Val  = 0.0;
                    bInfo.Info = (void *)lfi;
                    bInfo.func = LatitudeFunc;
                    Ta = -900; La = bInfo.func( (double)Ta, 0.0, (void *)lfi );
                    Tb =    0; Lb = bInfo.func( (double)Tb, 0.0, (void *)lfi );
                    done    = FALSE;
                    nAscend = 0;
                    while ( !done ) {
                        if ( (La < 0.0) && (Lb > 0.0) ) {
                            // we have an ascending node bracketed. Find Zero.
                            Lgm_zBrent( (double)Ta, (double)Tb, La, Lb, &bInfo, 1e-10, &Tmin, &Lmin );
                            if ( (Tmin >=0) && (Tmin < 86400) ) {
                                Lgm_Make_UTC( Date, Tmin/3600.0, &Ascend_UTC[nAscend], c );

                                // do the propagation
                                tsince = (Ascend_UTC[nAscend].JD - tle[10].JD)*1440.0;
                                LgmSgp_SGP4( tsince, lfi->sgp );

                                // Convert from TEME -> GEI2000
                                Uteme.x = lfi->sgp->X/WGS84_A; Uteme.y = lfi->sgp->Y/WGS84_A; Uteme.z = lfi->sgp->Z/WGS84_A; 
                                Lgm_Set_Coord_Transforms( Ascend_UTC[nAscend].Date, Ascend_UTC[nAscend].Time, c );
                                Lgm_Convert_Coords( &Uteme, &Ugei, TEME_TO_GEI2000, c );
                                Ascend_U[nAscend].x = Ugei.x; Ascend_U[nAscend].y = Ugei.y; Ascend_U[nAscend].z = Ugei.z;
                                Lgm_Convert_Coords( &Ascend_U[nAscend], &w, GEI2000_TO_WGS84, c );
                                Lgm_WGS84_to_GEOD( &w, &med->H5_Ascend_Geod[nAscend][0], &med->H5_Ascend_Geod[nAscend][1], &med->H5_Ascend_Geod[nAscend][2] );

                                Lgm_DateTimeToString( med->H5_Ascend_IsoTimes[nAscend], &Ascend_UTC[nAscend], 0, 3 );
                                printf("nNode: %d Bracket: T = %ld %ld    L = %g %g    Tnode, Lnode = %s %g\n", nAscend, Ta, Tb, La, Lb, med->H5_Ascend_IsoTimes[nAscend], Lmin );

                                ++nAscend;
                            }
                        }
                        // advance another 900 seconds.
                        Ta = Tb;     La = Lb;
                        Tb += 900; Lb = bInfo.func( (double)Tb, 0.0, (void *)lfi );
                        if (Tb > 86400+900) done = TRUE;
                    }
                    free( lfi->sgp );
                    free( lfi );
                    med->H5_nAscend = nAscend;



//PROBLEM AREA...
//BODY?

                    /*
                     * Open MagEphem txt file for writing and write header.
                     */
                    fp_MagEphem = fopen( OutFile, "w" );
                    Lgm_WriteMagEphemHeader( fp_MagEphem, argp_program_version, ExtModel, BODY, CommonName, IdNumber, IntDesig, CmdLine, nAscend, Ascend_UTC, Ascend_U, nPerigee, Perigee_UTC, Perigee_U, nApogee, &Apogee_UTC[0], &Apogee_U[0], MagEphemInfo );
                    printf("\t      Writing to file: %s\n", OutFile );

                    if ( UseEop ) {
                        // Read in the EOP vals
                        Lgm_read_eop( e );
                    }



                    /*
                     * Open MagEphem hdf5 file for writing and
                     * Create variables.
                     */
                    file    = H5Fcreate( HdfOutFile, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
//PROBLEM AREA...
//BODY?
                    Lgm_WriteMagEphemHeaderHdf( file, argp_program_version, ExtModel, BODY, CommonName, IdNumber, IntDesig, CmdLine, nAscend, Ascend_UTC, Ascend_U, nPerigee, Perigee_UTC, Perigee_U, nApogee, &Apogee_UTC[0], &Apogee_U[0], MagEphemInfo, med );




                    /*
                     *  Loop over seconds of the day
                     */
                    ss = (Date == StartDate) ? StartSeconds : 0;
                    es = (Date == EndDate) ? EndSeconds : 86400;
                    med->H5_nT = 0;
                    Lgm_ElapsedTimeInit( &t, 255, 150, 0 );
                    for ( Seconds=ss; Seconds<=es; Seconds += Delta ) {

                        Lgm_Make_UTC( Date, Seconds/3600.0, &UTC, c );
                        Lgm_DateTimeToString( IsoTimeString, &UTC, 0, 0 );

                        et = Lgm_TDBSecSinceJ2000( &UTC, c );

                        // do the propagation
                        tsince = (UTC.JD - tle[10].JD)*1440.0;
                        LgmSgp_SGP4( tsince, sgp );


// Lets just see what TLE we actually would get by searching...
int tiii = LgmSgp_FindTLEforGivenTime( nTle, tle, 1, UTC.JD, 1 );
if (tiii < 0 ) tiii = 0;
printf("tle index to use: %d  tsince = %g    (prev, next = %g %g)\n", tiii, (UTC.JD - tle[tiii].JD)*1440.0, (UTC.JD - tle[tiii-1].JD)*1440.0, (UTC.JD - tle[tiii+1].JD)*1440.0);





                        // Convert from TEME -> GEI2000
                        Uteme.x = sgp->X/WGS84_A; Uteme.y = sgp->Y/WGS84_A; Uteme.z = sgp->Z/WGS84_A; 
                        Lgm_Set_Coord_Transforms( UTC.Date, UTC.Time, c );
                        Lgm_Convert_Coords( &Uteme, &U, TEME_TO_GEI2000, c );


                        if ( UseEop ) {
                            // Get (interpolate) the EOP vals from the values in the file at the given Julian Date
                            Lgm_get_eop_at_JD( UTC.JD, &eop, e );

                            // Set the EOP vals in the CTrans structure.
                            Lgm_set_eop( &eop, c );
                        }

                        // Set mag model parameters
                        if ( FixModelDateTime ) {
                            Lgm_get_QinDenton_at_JD( ModelDateTime.JD, &p, (Verbosity > 0)? 1 : 0 );
                            Lgm_set_QinDenton( &p, MagEphemInfo->LstarInfo->mInfo );
                        } else {
                            Lgm_get_QinDenton_at_JD( UTC.JD, &p, (Verbosity > 0)? 1 : 0 );
                            Lgm_set_QinDenton( &p, MagEphemInfo->LstarInfo->mInfo );
                        }


                        if ( ForceKp >= 0.0 ) {
                            MagEphemInfo->LstarInfo->mInfo->fKp = ForceKp;
                            MagEphemInfo->LstarInfo->mInfo->Kp  = (int)(ForceKp+0.5);
                            if (MagEphemInfo->LstarInfo->mInfo->Kp > 6) MagEphemInfo->LstarInfo->mInfo->Kp = 6;
                            if (MagEphemInfo->LstarInfo->mInfo->Kp < 0 ) MagEphemInfo->LstarInfo->mInfo->Kp = 0;
                        }


                        // Set up the trans matrices
                        Lgm_Set_Coord_Transforms( UTC.Date, UTC.Time, c );

                        MagEphemInfo->OrbitNumber = GetOrbitNumber( &UTC, nPerigee, Perigee_UTC, PerigeeOrbitNumber );

                        Lgm_Convert_Coords( &U, &Rgsm, GEI2000_TO_GSM, c );

//PROBLEM AREA...
//sce2s_c( BODY,    et, 30, sclkch );

                        /*
                         * Compute L*s, Is, Bms, Footprints, etc...
                         * These quantities are stored in the MagEphemInfo Structure
                         */
                        if ( Verbosity > 0 ) {
                            printf("\t\t"); Lgm_PrintElapsedTime( &t ); printf("\n");
                            printf("\n\n\t[ %s ]: %s  Bird: %s Kp: %g    Rgsm: %g %g %g Re\n", ProgramName, IsoTimeString, Bird, MagEphemInfo->LstarInfo->mInfo->fKp, Rgsm.x, Rgsm.y, Rgsm.z );
                            printf("\t-------------------------------------------------------------------------------------------------------------------\n");
                        }
                        Lgm_ComputeLstarVersusPA( UTC.Date, UTC.Time, &Rgsm, nAlpha, Alpha, Colorize, MagEphemInfo );

                        MagEphemInfo->InOut = InOutBound( ApoPeriTimeList, nApoPeriTimeList, UTC.JD );

                        /*
                         * Write a row of data into the txt file
                         */
                        Lgm_WriteMagEphemData( fp_MagEphem, IntModel, ExtModel, MagEphemInfo->LstarInfo->mInfo->fKp, MagEphemInfo->LstarInfo->mInfo->Dst, MagEphemInfo );




                        if ( DumpShellFiles && (nAlpha > 0) ){

                            sprintf( ShellFile, "%s_%ld.dat", OutFile, Seconds );
                            printf( "Writing Full Shell File: %s\n", ShellFile );
                            WriteMagEphemInfoStruct( ShellFile, nAlpha, MagEphemInfo );
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

                        med->H5_Rgsm[ med->H5_nT ][0]        = Rgsm.x;
                        med->H5_Rgsm[ med->H5_nT ][1]        = Rgsm.y;
                        med->H5_Rgsm[ med->H5_nT ][2]        = Rgsm.z;

                        Lgm_Set_Coord_Transforms( UTC.Date, UTC.Time, c );
                        Lgm_Convert_Coords( &Rgsm, &Rgeo, GSM_TO_GEO, c );      Lgm_VecToArr( &Rgeo, &med->H5_Rgeo[ med->H5_nT ][0] );
                        Lgm_Convert_Coords( &Rgsm, &W,    GSM_TO_SM, c );       Lgm_VecToArr( &W,    &med->H5_Rsm[ med->H5_nT ][0] );
                        Lgm_Convert_Coords( &Rgsm, &W,    GSM_TO_GEI2000, c );  Lgm_VecToArr( &W,    &med->H5_Rgei[ med->H5_nT ][0] );
                        Lgm_Convert_Coords( &Rgsm, &W,    GSM_TO_GSE, c );      Lgm_VecToArr( &W,    &med->H5_Rgse[ med->H5_nT ][0] );

                        Lgm_WGS84_to_GEOD( &Rgeo, &GeodLat, &GeodLong, &GeodHeight );
                        Lgm_SetArrElements3( &med->H5_Rgeod[ med->H5_nT ][0],        GeodLat, GeodLong, GeodHeight );
                        Lgm_SetArrElements2( &med->H5_Rgeod_LatLon[ med->H5_nT ][0], GeodLat, GeodLong );
                        med->H5_Rgeod_Height[ med->H5_nT ] = GeodHeight;

                        Lgm_Convert_Coords( &Rgsm, &W, GSM_TO_CDMAG, c );
                        Lgm_CDMAG_to_R_MLAT_MLON_MLT( &W, &R, &MLAT, &MLON, &MLT, c );
                        med->H5_CDMAG_MLAT[ med->H5_nT ] = MLAT;
                        med->H5_CDMAG_MLON[ med->H5_nT ] = MLON;
                        med->H5_CDMAG_MLT[ med->H5_nT ]  = MLT;
                        med->H5_CDMAG_R[ med->H5_nT ]    = R;

                        Lgm_Convert_Coords( &Rgsm, &W, GSM_TO_EDMAG, c );
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


                        MagEphemInfo->LstarInfo->mInfo->Bfield( &Rgsm, &Bsc_gsm, MagEphemInfo->LstarInfo->mInfo );
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
                        Lgm_WriteMagEphemDataHdf( file, med->H5_nT, med );
                        ++(med->H5_nT);


                    }

                    fclose(fp_MagEphem);
                    H5Fclose( file );

                    printf("DONE.\n");
                    Lgm_PrintElapsedTime( &t );
                    Lgm_SetElapsedTimeStr( &t );
                    sprintf( Command, "sed -i '/ELAPSED_TIME/s++%s+g' %s", t.ElapsedTimeStr, OutFile); system( Command );

// PROBLEM AREA?
//sprintf( Command, "sed -i '/SPICE_KERNEL_FILES_LOADED/s++%s+' %s", SpiceKernelFilesLoaded, OutFile); system( Command );


                    /*
                     * Write KML file -- testing....
                     */
                    FILE *KmlFile = fopen( "Puke.kml", "w" );
                    Lgm_WriteMagEphemDataKML( KmlFile, med->H5_nT, med );
                    fclose( KmlFile );





// PROBLEM AREA?
/*
* Unload spice kernels
*/
//unload_c( InFile );






                } //end else
            } // end "if ( !FileExists || Force )" control structure


        } // end birds loop
    } // end JD loop


    free( tle );
    free( sgp );

    free( ShellFile );
    free( OutFile );
    free( HdfOutFile );
    free( InFile );

    Lgm_free_ctrans( c );
    Lgm_destroy_eop( e );
    Lgm_FreeMagEphemInfo( MagEphemInfo );
    Lgm_FreeMagEphemData( med );

    LGM_ARRAY_1D_FREE( CmdLine );
    LGM_ARRAY_2D_FREE( Birds );
    LGM_ARRAY_1D_FREE( Perigee_UTC );
    LGM_ARRAY_1D_FREE( Apogee_UTC );
    LGM_ARRAY_1D_FREE( Ascend_UTC );
    LGM_ARRAY_1D_FREE( Perigee_U );
    LGM_ARRAY_1D_FREE( Apogee_U );
    LGM_ARRAY_1D_FREE( Ascend_U );
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
