#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <argp.h>
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

void StringSplit( char *Str, char *StrArray[], int len, int *n );


const  char *ProgramName = "LastClosedDriftShell";
const  char *argp_program_version     = "LastClosedDriftShell_0.1";
const  char *argp_program_bug_address = "<smorley@lanl.gov>";
static char doc[] = "\nComputes last closed drift shell for given times, pitch angles and magnetic field model."
                    " \n\n OutFile is a path to the output file(s) that may contain"
                    " variables that will be substituted if a time range is given.  The"
                    " available time variables are '%YYYY', '%MM', and '%DD' which correspond"
                    " repectively to 4-digit year, 2-digit month (Jan is 01), and 2-digit day of"
                    " month.\n"
                    " Here is an example using time variables,\n\n \t./LastClosedDriftShell -S 20020901 -E 20020930\n"
                    " \t\t/home/jsmith/MagEphemData/%YYYY/%YYYY%MM%DD_%EE_MagEphem.txt.\n\n Directories"
                    " in the output file will be created if they don't already exist.\n\n";


// Mandatory arguments
#define     nArgs   1
static char ArgsDoc[] = "OutFile";

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
    {"PitchAngles",     'p',    "\"start_pa, end_pa, npa\"",  0,                                      "Pitch angles to compute. Default is \"5.0, 90, 18\"." },
    {"FootPointHeight", 'f',    "height",                     0,                                      "Footpoint height in km. Default is 100km."                  },
    {"Quality",         'q',    "quality",                    0,                                      "Quality to use for L* calculations. Default is 3."      },
    {"Delta",           'd',    "delta",                      0,                                      "Cadence of LCDS calculation [minutes]. Default is 30."      },
    {"StartDate",       'S',    "yyyymmdd",                   0,                                      "StartDate "                              },
    {"EndDate",         'E',    "yyyymmdd",                   0,                                      "EndDate "                                },
    {"UseEop",          'e',    0,                            0,                                      "Use Earth Orientation Parameters whn comoputing ephemerii" },
    {"Force",           'F',    0,                            0,                                      "Overwrite output file even if it already exists" },
    {"verbose",         'v',    "verbosity",                  0,                                      "Produce verbose output"                  },
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
    int         Force;
    double      FootPointHeight;
    double      Delta;

    char        IntModel[80];
    char        ExtModel[80];

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
        case 'd':
            sscanf( arg, "%lf", &arguments->Delta );
            break;
        case 'q':
            arguments->Quality = atoi( arg );
            break;
        case 's':
            arguments->silent = 1;
            break;
        case 'v':
            arguments->verbose = atoi( arg );
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


int main( int argc, char *argv[] ){

    struct Arguments arguments;
    double           UTC, brac1, brac2, tol, sJD, eJD, JD, jDate, t_cadence;
    double           K[500], LS[500], Alpha[500];
    double           Inc, FootpointHeight;
    int              Force, UseEop;
    long int         StartDate, EndDate, Date, currDate;
    int              nAlpha, i, Quality, ans, aa, Year, Month, Day;
    char             Str[128], NewStr[2048];
    char             IntModel[20], ExtModel[20];
    char             Filename[1024];
    Lgm_LstarInfo    *LstarInfo = InitLstarInfo(0);
    Lgm_LstarInfo    *LstarInfo3;
    FILE             *fp;
    Lgm_DateTime     DT_UTC;
    Lgm_QinDentonOne qd;

    tol = 0.001;


   /*
     * Default option values.
     */
    arguments.StartPA         = 90;     // start at 90.0 Deg.
    arguments.EndPA           = 5.0;    // stop at 2.5 Deg.
    arguments.nPA             = 18;     // 18 pitch angles
    arguments.silent          = 0;
    arguments.verbose         = 0;
    arguments.Quality         = 2;
    arguments.Delta           = 30;
    arguments.Force           = 0;
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
     *  Set other options
     */
    FootpointHeight = arguments.FootPointHeight;
    Quality         = arguments.Quality;
    t_cadence       = arguments.Delta/1440.0; //needs to be in days
    Force           = arguments.Force;
    UseEop          = arguments.UseEop;
    StartDate       = arguments.StartDate;
    EndDate         = arguments.EndDate;
    strcpy( IntModel,  arguments.IntModel );
    strcpy( ExtModel,  arguments.ExtModel );

    // Settings for Lstar calcs
    LstarInfo->VerbosityLevel = arguments.verbose;
    LstarInfo->mInfo->VerbosityLevel = 0;
    LstarInfo->VerbosityLevel       = 0;
    LstarInfo->mInfo->Lgm_LossConeHeight = FootpointHeight;

    if ( !strcmp( ExtModel, "T87" ) ){
        LstarInfo->mInfo->Bfield = Lgm_B_T87;
    } else if ( !strcmp( ExtModel, "CDIP" ) ){
        LstarInfo->mInfo->Bfield = Lgm_B_cdip;
    } else if ( !strcmp( ExtModel, "EDIP" ) ){
        LstarInfo->mInfo->Bfield = Lgm_B_edip;
    } else if ( !strcmp( ExtModel, "DUNGEY" ) ){
        LstarInfo->mInfo->Bfield = Lgm_B_Dungey;
        LstarInfo->mInfo->InternalModel = LGM_DUNGEY;
        strcpy(IntModel, "DUNGEY");
    } else if ( !strcmp( ExtModel, "IGRF" ) ){
        LstarInfo->mInfo->Bfield = Lgm_B_igrf;
    } else if ( !strcmp( ExtModel, "TS04" ) ){
        LstarInfo->mInfo->Bfield = Lgm_B_TS04;
    } else if ( !strcmp( ExtModel, "T89c" ) ){
        LstarInfo->mInfo->Bfield = Lgm_B_T89c;
    } else if ( !strcmp( ExtModel, "T96" ) ){
        LstarInfo->mInfo->Bfield = Lgm_B_T96;
    } else if ( !strcmp( ExtModel, "OP77" ) ){
        LstarInfo->mInfo->Bfield = Lgm_B_OP77;
    } else { 
        // default
        LstarInfo->mInfo->Bfield = Lgm_B_T89c;
    }

    if ( !strcmp( IntModel, "CDIP" ) ){
        LstarInfo->mInfo->InternalModel = LGM_CDIP;
    } else if ( !strcmp( IntModel, "EDIP" ) ){
        LstarInfo->mInfo->InternalModel = LGM_EDIP;
    } else {
        // default
        LstarInfo->mInfo->InternalModel = LGM_IGRF;
    }

    sJD = Lgm_Date_to_JD( StartDate, 0.0, LstarInfo->mInfo->c);
    eJD = Lgm_Date_to_JD( EndDate, 23.9999, LstarInfo->mInfo->c);


    for (jDate = sJD; jDate <= eJD; jDate+= 1.0 ) {
    
        /*
         *  Set output filename, etc.
         */
        currDate = Lgm_JD_to_Date( jDate, &Year, &Month, &Day, &UTC );
        NewStr[0]       = '\0';
        strcpy( Filename, arguments.args[0] );
    
        // Date and UTC
        //Lgm_Doy( StartDate, &Year, &Month, &Day, &Doy);
        NewStr[0]       = '\0';
        sprintf( Str, "%4d", Year );   Lgm_ReplaceSubString( NewStr, Filename, "%YYYY", Str );  strcpy( Filename, NewStr );
        sprintf( Str, "%02d", Month ); Lgm_ReplaceSubString( NewStr, Filename, "%MM", Str );   strcpy( Filename, NewStr );
        sprintf( Str, "%02d", Day );   Lgm_ReplaceSubString( NewStr, Filename, "%DD", Str );   strcpy( Filename, NewStr );
        sprintf( Str, "%s", ExtModel );   Lgm_ReplaceSubString( NewStr, Filename, "%EE", Str );   strcpy( Filename, NewStr );
        fp = fopen(Filename, "w");
    
        // Bracket Position in GSM
        brac1 = -3.0;
        brac2 = -13.0;
    
    
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
        for ( JD = jDate; JD < jDate+1.0; JD += t_cadence ) {
            //set date specific stuff
            Date = Lgm_JD_to_Date( JD, &Year, &Month, &Day, &UTC );
            Lgm_Set_Coord_Transforms( Date, UTC, LstarInfo->mInfo->c);
        
            Lgm_get_QinDenton_at_JD( JD, &qd, 1 );
            Lgm_set_QinDenton( &qd, LstarInfo->mInfo );
        
            if ( !strcmp( ExtModel, "T89c" )) { //truncate Kp at 5 for T89
                if ( LstarInfo->mInfo->Kp > 5 ) LstarInfo->mInfo->Kp = 5;
                }
    
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
        
                    ans = Lgm_LCDS( Date, UTC, brac1, brac2, Alpha[aa], tol, Quality, &K[aa], LstarInfo3 );
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
            
            Lgm_Make_UTC( Date, UTC, &DT_UTC, LstarInfo->mInfo->c );
            Lgm_DateTimeToString( Str, &DT_UTC, 0, 3);
            fprintf(fp, "%24s",     Str );
            for ( i=0; i<nAlpha; ++i ) {
                fprintf( fp, "     %13g", LS[i]);
            }
            for ( i=0; i<nAlpha; ++i ) {
                fprintf( fp, "     %13g", K[i] );
            }
            fprintf(fp, "%s", " \n");
        
        }
        fclose(fp);
    }
    FreeLstarInfo( LstarInfo );

    return(0);

}
