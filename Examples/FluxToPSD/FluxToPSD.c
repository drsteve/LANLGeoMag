/*
CHECK ME FOR MEMORY LEAKS.....
THIS ALSO eeds some argument changes...
e.g. The Mu's K's instead of pitch angles...
*/
#include <argp.h>
#include <ctype.h>
#include <fcntl.h>
#include <libgen.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>
#include <Lgm_CTrans.h>
#include <Lgm_DynamicMemory.h>
#include <Lgm_Eop.h>
#include <Lgm_FluxToPsd.h>
#include <Lgm_HDF5.h>
#include <Lgm_MagModelInfo.h>
#include <Lgm_Misc.h>
#include <Lgm_QinDenton.h>
#include <Lgm_Sgp.h>
#include "CDFUtils.h"



#define KP_DEFAULT 0
void InterpVec( double *JD, Lgm_Vector *Ugsm, int n1, double cJD, Lgm_Vector *v );
void InterpArr( double *JD, double *y, int n, double cJD, double *r ) ;
void WriteHDF5PSDHeader(hid_t file, int nMu, int nK, double *Mu, double *K);
void WriteHDF5PSDData(hid_t file, int i, int nMu, int nK, char **IsoTimes,
                      double *PsdLstar, double *AlphaOfK, double **PSD);

const  char *argp_program_version     = "FluxToPSD 1.0";
const  char *argp_program_bug_address = "<mghenderson@lanl.gov>";
static char doc[] = "Converts Flux to PSD. Blah Blah";

// Mandatory arguments
#define     nArgs   2
static char ArgsDoc[] = "MergedFluxFile StdPsdFile";


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
    {"ModelTime",       'T',    "ISO8601_time_string",        0,                                      "Force the model to use parameters for the given date/time. This can be used to force the model to be static instead of dynamic."    },
    {"Mus",             'M',    "\"start_Mu, end_Mu, nMu\"",  0,                                      "Mu values to compute. Default is \"5, 5000.0, 18\"." },
    {"Ks",              'K',    "\"start_K, end_K, nK\"",     0,                                      "K values to compute. Default is \"0.01, 10.0, 18\"." },
    {"FootPointHeight", 'f',    "height",                     0,                                      "Footpoint height in km. Default is 100km."                  },
    {"Quality",         'q',    "quality",                    0,                                      "Quality to use for L* calculations. Default is 3."      },
    {"FitType",         't',    "fit_type",                   0,                                      "Fit type for PSD energy spectrum: default (0) is 2 rel. Maxwellians, set to 1 for spline fit."},
    {"StartDate",       'S',    "yyyymmdd",                   0,                                      "StartDate "                              },
    {"EndDate",         'E',    "yyyymmdd",                   0,                                      "EndDate "                                },
    {"Append",          'a',    0,                            0,                                      "Append to OutFile instead of creating a new one"  },
//    {"UseEop",          'e',    0,                            0,                                      "Use Earth Orientation Parameters whn comoputing ephemerii" },
    {"Colorize",        'c',    0,                            0,                                      "Colorize output"                         },
    {"Force",           'F',    0,                            0,                                      "Overwrite output file even if it already exists" },
    {"verbose",         'v',    "verbosity",                  0,                                      "Produce verbose output"                  },
    {"silent",          's',    0,                            OPTION_ARG_OPTIONAL | OPTION_ALIAS                                                },
    {"nchannels",       'n',    "num_channels",               0,                                      "Use first n energy channels"                  },
    { 0 }
};

struct Arguments {
    char        *args[ nArgs ];
    int         silent;
    int         verbose;
    int         nchannels;

    double      StartMu;
    double      EndMu;
    int         nMu;

    double      StartK;
    double      EndK;
    int         nK;

    int         Quality;
    int         Colorize;
    int         Force;
    double      FootPointHeight;

    char        IntModel[80];
    char        ExtModel[80];

    int         Append;
    int         UseEop;
    int         FitType;


    long int    StartDate;
    long int    EndDate;

    int         FixModelDateTime;
    char        ModelDateTime[80];
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
        case 'T': // Date/time to use for forcing model parameters to be static.
            strcpy( arguments->ModelDateTime, arg );
            arguments->FixModelDateTime = 1;
            break;
        case 'K':
            sscanf( arg, "%lf, %lf, %d", &arguments->StartK, &arguments->EndK, &arguments->nK );
            break;
        case 'M':
            sscanf( arg, "%lf, %lf, %d", &arguments->StartMu, &arguments->EndMu, &arguments->nMu );
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
            arguments->verbose = atoi( arg );
            break;
        case 'n':
            arguments->nchannels = atoi( arg );
            break;
        case 't':
            arguments->FitType = atoi( arg );
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
 *
 */
int main( int argc, char *argv[] ){

    struct Arguments arguments;

    Lgm_CTrans       *c = Lgm_init_ctrans( 0 );
    Lgm_DateTime     ModelDateTime;
    Lgm_Eop          *e = Lgm_init_eop( 0 );
    Lgm_EopOne       eop;
    int              FixModelDateTime;
    int              i, iii, jjj, kk;
    char             InputFilename[1024];
    char             OutputFilename[1024];
    char             IntModel[20], ExtModel[20];
    int              AppendMode, UseEop, Colorize, Force, nv, Quality;
    double           FootpointHeight;
    long int         StartDate, EndDate, Date __attribute__((unused));
    int              sYear, sMonth, sDay, sDoy, eYear, eMonth, eDay, eDoy, Year, Month, Day;
    int              cYear, cMonth, cDay, cDoy;
    double           sJD, eJD, JD, Time, cJD;
    hsize_t           dims[3];
    double           *Flux_Energy;
    double           *Flux_Alpha;
    int               Flux_nEnergy;
    int               Flux_nAlpha;
    double           **GSM;
    double           **MagEphem_Lstar;
    double           *MagEphem_Alpha;
    double           ***J;
    int              fromHDF, retval __attribute__((unused));


    /*
     * App-specific declarations
     */
//    Lgm_CTrans      *c = Lgm_init_ctrans( 0 );
    Lgm_DateTime     d;
    Lgm_Vector       u, v, v_geo, Ugsm[300];
    Lgm_FluxToPsd   *f2p;
    Lgm_DateTime     UTC_DateTime;
    Lgm_MagModelInfo *mInfo;
    Lgm_QinDentonOne QinDen;
    double           *E, *A, **J_DIFF;
    int              nE, nA, j, n, k, nTimes, MagEphem_nAlpha;
    double           *Mu, *K, *LstarOfK, *AlphaOfK, *PsdLstar;
    int              nMu, nK;
    int              OverRideKp;
    double           T89Q_Kp;
    char             Line[5000], **IsoTimes, buf[512], FillStr[20], fitDesc[25];

    hid_t            file;
    hid_t            atype;
    hid_t            space;

    herr_t           status __attribute__((unused));

    hsize_t          Dims[4], DataSet, MemSpace, SlabSize[4];
    CDFid       id;
    CDFstatus   cdf_status;
    int         Ni, Nj, Nk, numMaxw, fitMaxw;
    long        numDims, dimSizes[CDF_MAX_DIMS], numRecs;
    char        varName[CDF_VAR_NAME_LEN256+1];
    long        dataType, numElems, recVary, dimVarys[CDF_MAX_DIMS];
    double      **Position, ***FEDU;
    double *Epoch = NULL;
    double  *Flux_Energy_Tmp;
    double  *Flux_Alpha_Tmp;
    double **FEDU_Lstar_Tmp;
    char **ProcessedTime = NULL;
    Lgm_DateTime CurrentUTC;
    int nProcTimes;
    int H5_nT;
    int iTime;
    int AppendFlag;


    /*
     * Default option values.
     */
    arguments.StartMu          =  5.0;    // start Mu [MeV/G]
    arguments.EndMu            =  5000.0; // stop Mu [MeV/G]
    arguments.nMu              =  18;     // 18 Mus
    arguments.StartK           =  0.01;   // start K [Re G^.5]
    arguments.EndK             =  10.0;   // stop K [Re G^.5]
    arguments.nK               =  18;     // 18 K
    arguments.silent           =  0;
    arguments.verbose          =  0;
    arguments.nchannels        =  0;
    arguments.Quality          =  3;
    arguments.Colorize         =  0;
    arguments.Force            =  0;
    arguments.Append           =  0;
    arguments.UseEop           =  0;
    arguments.FixModelDateTime =  0;
    arguments.StartDate        = -1;
    arguments.EndDate          = -1;
    arguments.FootPointHeight  =  100.0; // km
    numMaxw = 2;
    arguments.FitType          =  0; //0 for 2 Maxwellians, 1 for spline
    strcpy( arguments.IntModel,  "IGRF" );
    strcpy( arguments.ExtModel,  "T89" );
    strcpy( arguments.ModelDateTime,  "" );
    fitMaxw          = arguments.FitType; 
    if ( fitMaxw==0 ) {
        strcpy( fitDesc, "Double Maxwell-Juttner" );
    } else {
        strcpy( fitDesc, "B-spline" );
    }


    /*
     *  Parse CmdLine arguments and options
     */
    argp_parse (&argp, argc, argv, 0, 0, &arguments);


    /*
     *  Set input and output filenames
     */
    strcpy( InputFilename, arguments.args[0] );
    strcpy( OutputFilename, arguments.args[1] );

    /*
     *  Set other options
     */
    FootpointHeight  = arguments.FootPointHeight;
    Quality          = arguments.Quality;
    Colorize         = arguments.Colorize;
    Force            = arguments.Force;
    AppendMode       = arguments.Append;
    UseEop           = arguments.UseEop;
    StartDate        = arguments.StartDate;
    EndDate          = arguments.EndDate;
    FixModelDateTime = arguments.FixModelDateTime;
    strcpy( IntModel,  arguments.IntModel );
    strcpy( ExtModel,  arguments.ExtModel );

    if ( FixModelDateTime ) IsoTimeStringToDateTime( arguments.ModelDateTime, &ModelDateTime, c );





    /*
     *  Choose a set of Mu's and K's
     *  and compute Phase Space Density at these constant Mu's and K's
     */
    nMu = arguments.nMu;
    nK  = arguments.nK;
    LGM_ARRAY_1D( Mu, nMu, double );
    LGM_ARRAY_1D( K, nK, double );
    LGM_ARRAY_1D( AlphaOfK, nK, double );
    LGM_ARRAY_1D( LstarOfK, nK, double );
    LGM_ARRAY_1D( PsdLstar, nK, double );

    // Geometric progressions for K's and Mu's
    Lgm_GeometricSeq( arguments.StartK,  arguments.EndK,   nK,   K );
    Lgm_GeometricSeq( arguments.StartMu, arguments.EndMu,  nMu,  Mu );




    /*
     *  Print summary of our options
     */
    if ( !arguments.silent ) {
        printf("\n\n");
        printf( "\t         Program/Version: %s\n", argp_program_version );
        printf( "\t          Bug Reports To: %s\n", argp_program_bug_address );
        printf( "\tFlux/MagEphem Input File: %s\n", InputFilename );
        printf( "\t         Output PSD File: %s\n", OutputFilename );
        printf( "\t     Number of Mu Values: %d\n", nMu );
        printf( "\t       Mu Values [MeV/G]:" );
        for (i=0, nv=0; i<nMu; i++ ) {
            printf( " %12g", Mu[i] );
            if ( i==nMu-1 ) printf("\n");
            else if (nv>=8) { nv=0; printf("\n\t\t\t\t"); }
            else { ++nv; };
        }
        printf( "\t     Number of K Values: %d\n", nK );
        printf( "\t     K Values [RE G^.5]:" );
        for (i=0, nv=0; i<nK; i++ ) {
            printf( " %12g", K[i] );
            if ( i==nK-1 ) printf("\n");
            else if (nv>=8) { nv=0; printf("\n\t\t\t\t"); }
            else { ++nv; };
        }
        printf( "\t         Internal Model: %s\n", IntModel );
        printf( "\t         External Model: %s\n", ExtModel );
        printf( "\t           Spectral Fit: %s\n", fitDesc );
        if ( FixModelDateTime ) {
            printf( "\t  Model Params Fixed to: %s\n", arguments.ModelDateTime );
        }
        printf( "\t   FootpointHeight [km]: %g\n", FootpointHeight );
        printf( "\t             L* Quality: %d\n", Quality );
        printf( "\t           Force output: %s\n", Force ? "yes" : "no" );
        printf( "\t   Append to OutputFile: %s\n", AppendMode ? "yes" : "no" );
        printf( "\t                Use Eop: %s\n", UseEop ? "yes" : "no" );
        printf( "\t Colorize Thread Output: %d\n", Colorize );
        printf( "\t                Verbose: %d\n", arguments.verbose);
        printf( "\t                 Silent: %s\n", arguments.silent  ? "yes" : "no" );
        if ( (StartDate > 0)&&(EndDate > 0) ){
            printf( "\t              StartDate: %ld\n", StartDate );
            printf( "\t                EndDate: %ld\n", EndDate );
        }
    }


    mInfo = Lgm_InitMagInfo( );

    OverRideKp = FALSE;
    if ( !strcmp( ExtModel, "T87" ) ){
        mInfo->Bfield = Lgm_B_T87;
    } else if ( !strcmp( ExtModel, "CDIP" ) ){
        mInfo->Bfield = Lgm_B_cdip;
    } else if ( !strcmp( ExtModel, "EDIP" ) ){
        mInfo->Bfield = Lgm_B_edip;
    } else if ( !strcmp( ExtModel, "IGRF" ) ){
        mInfo->Bfield = Lgm_B_igrf;
    } else if ( !strcmp( ExtModel, "TS04" ) ){
        mInfo->Bfield = Lgm_B_TS04;
    } else if ( !strcmp( ExtModel, "OP77Q" ) ){
        mInfo->Bfield = Lgm_B_OP77;
    } else if ( !strcmp( ExtModel, "OP77" ) ){
        mInfo->Bfield = Lgm_B_OP77;
    } else if ( !strcmp( ExtModel, "T89Q" ) ){
        mInfo->Bfield = Lgm_B_T89c;
        T89Q_Kp = 2;
        OverRideKp = TRUE;
    } else if ( !strcmp( ExtModel, "T89D" ) ){
        mInfo->Bfield = Lgm_B_T89c;
    } else if ( !strcmp( ExtModel, "T96" ) ){
        mInfo->Bfield = Lgm_B_T96;
    } else if ( !strcmp( ExtModel, "TS04D" ) ){
        mInfo->Bfield = Lgm_B_TS04;
    } else if ( !strcmp( ExtModel, "T01S" ) ){
        mInfo->Bfield = Lgm_B_T01S;
    } else if ( !strcmp( ExtModel, "T01" ) ){
        mInfo->Bfield = Lgm_B_T01S;
    } else if ( !strcmp( ExtModel, "T02" ) ){
        mInfo->Bfield = Lgm_B_T02;
    } else { //if ( !strcmp( ExtModel, "T89" ) ){
        // default
        mInfo->Bfield = Lgm_B_T89c;
    }

    if ( !strcmp( IntModel, "CDIP" ) ){
        mInfo->InternalModel = LGM_CDIP;
    } else if ( !strcmp( IntModel, "EDIP" ) ){
        mInfo->InternalModel = LGM_EDIP;
    } else {
        // default
        mInfo->InternalModel = LGM_IGRF;
    }




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





    /*
     * Loop over the days that we nbeed to process. This translates into processing a range of files.
     */
    for ( JD = sJD; JD <= eJD; JD += 1.0 ) {


        strcpy( InFile, InputFilename );
        strcpy( OutFile, OutputFilename );

        if ( SubstituteDates ) {

            Date = Lgm_JD_to_Date( JD, &Year, &Month, &Day, &Time );

            // Substitute times in the files.
            NewStr[0] = '\0';
            sprintf( Str, "%4d", Year );   Lgm_ReplaceSubString( NewStr, InFile, "%YYYY", Str );  strcpy( InFile, NewStr );
            sprintf( Str, "%02d", Month ); Lgm_ReplaceSubString( NewStr, InFile, "%MM", Str );    strcpy( InFile, NewStr );
            sprintf( Str, "%02d", Day );   Lgm_ReplaceSubString( NewStr, InFile, "%DD", Str );    strcpy( InFile, NewStr );

            sprintf( Str, "%4d", Year );   Lgm_ReplaceSubString( NewStr, OutFile, "%YYYY", Str );  strcpy( OutFile, NewStr );
            sprintf( Str, "%02d", Month ); Lgm_ReplaceSubString( NewStr, OutFile, "%MM", Str );    strcpy( OutFile, NewStr );
            sprintf( Str, "%02d", Day );   Lgm_ReplaceSubString( NewStr, OutFile, "%DD", Str );    strcpy( OutFile, NewStr );

        }


        if ( SubstituteDates ) {
            printf( "\n\n      Processing Date: %4d-%02d-%02d\n", Year, Month, Day );
        } else {
            printf( "\n\n      Processing Files: %s\n", InFile);
        }
        printf( "      -------------------------------------------------------------------------------------------\n");
        printf( "      Flux/MagEphem Input File: %s\n", InFile);
        printf( "               PSD output File: %s\n\n", OutFile);


        // Create Base directory if it hasnt been created yet.
        char *dirc = strdup( OutFile );
        BaseDir    = dirname(dirc);
        sprintf( Command, "mkdir -p %s", BaseDir); retval = system( Command );
        free( dirc );


        /*
         *   Check to see if OutFile exists or not.
         */
        int     StatError, FileExists, InFileExists;
        struct stat StatBuf;
        StatError = stat( OutFile, &StatBuf );

        FileExists = FALSE;
        if ( StatError != -1 ) {

            FileExists = TRUE;
            if ( !Force && !AppendMode) {
              printf("\n\n\tOutfile already exists (use -F option to force processing): %s\n", OutFile );
            } else if (!AppendMode) {
              printf("\n\n\tWarning. Existing Outfile will be overwritten: %s \n", OutFile );
            } else {
              printf("\n\n\tExisting Outfile will be appended: %s \n", OutFile );
            }

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

            /*
             *   If OutFile exists, check to see if it is readable.
             */
            if ( !H5Fis_hdf5( OutFile ) ) {
                printf("\t  Outfile Not Readable: %s. Forcing regeneration of file.\n\n\n", OutFile );
                FileExists = FALSE;
            }

        }


        InFileExists = ( (StatError = stat( InFile, &StatBuf )) != -1 ) ? TRUE : FALSE;

       /*
        * Read in a combined Flux/MagEphem HDF5 File
        */
        if ( InFileExists ) {
            fromHDF = 0;
            if (H5Fis_hdf5( InFile )){
                fromHDF = 1;
                // Read the InFile
                file = H5Fopen( InFile, H5F_ACC_RDONLY, H5P_DEFAULT );

                //Read ISO time
                IsoTimes = Get_StringDataset_1D( file, "/Epoch", dims);

                // Read GSM and Lstar position
                Position  = Get_DoubleDataset_2D( file, "/Ephemeris/GSM", dims);
                MagEphem_Lstar  = Get_DoubleDataset_2D( file, "/Lstar", dims);

                //Read pitch angles for input flux
                MagEphem_Alpha = Get_DoubleDataset_1D( file, "/Lstar_Alpha", dims);
                MagEphem_nAlpha = dims[0];

                //Read flux array, pitch angle and energies
                Flux_Alpha  = Get_DoubleDataset_1D( file, "/PA", dims);
                Flux_nAlpha = dims[0];
                Flux_Energy = Get_DoubleDataset_1D( file, "/Energy", dims);
                Flux_nEnergy = dims[0];
                J = Get_DoubleDataset_3D( file, "/Flux", dims);
                nTimes = dims[0];
                if (Flux_nEnergy != dims[2]) {
                    printf("Dimensions of input flux array are incorrectly ordered! Exiting.\n");
                    exit(-1);
                }
                status = H5Fclose( file );
                }
            else {
    
                /*
                 *  Read CDF file
                 */
                cdf_status = CDFopenCDF( InFile, &id );
                if ( cdf_status != CDF_OK ) StatusHandler( cdf_status );
    
                /*
                 *  Get FEDU variable.
                 */
                i = CDFgetVarNum( id, "FEDU");
                CDFinquirezVar( id, i, varName, &dataType, &numElems, &numDims, dimSizes, &recVary, dimVarys );
                if ( cdf_status != CDF_OK ) StatusHandler( cdf_status );
                cdf_status = CDFgetzVarNumRecsWritten( id, i, &numRecs );
                if (cdf_status != CDF_OK) StatusHandler( cdf_status );
    
                Ni = numRecs;
                Nj = dimSizes[0];
                Nk = dimSizes[1];
                printf("FEDU dimensions [time][alpha][E]: %d %d %d\n", Ni, Nj, Nk );
    
                LGM_ARRAY_3D( FEDU, Ni, Nj, Nk, double );
                cdf_status = CDFgetzVarAllRecordsByVarID( id, i, &FEDU[0][0][0] );
                if (cdf_status != CDF_OK) StatusHandler( cdf_status );
    
                nTimes       = Ni;
                Flux_nAlpha  = (int)(Nj/2)+1; // assume an odd number like 17. Then 17/2+1 is 9
                //Flux_nEnergy = Nk-1; // drop last integral channel #Integral channel shouldn't be in FEDU variable#
                Flux_nEnergy = (arguments.nchannels != 0) ? arguments.nchannels : Nk-1; //Nk-1 used for compatibility with prev revs
                printf("J dimensions [time][alpha][E]: %d %d %d\n", Ni, Flux_nAlpha, Flux_nEnergy );
                LGM_ARRAY_3D( J, Ni, Flux_nAlpha, Flux_nEnergy, double );
                /*
                 *  Copy over array. But these arrays go 0-180 in PA. So only use half of the PA range.
                 */
                for (i=0; i<Ni; i++){
                    for (j=0; j<Flux_nAlpha; j++){
                        for (k=0; k<Flux_nEnergy; k++){
                            J[i][j][k] = FEDU[i][j][k];
                        }
                    }
                }
                LGM_ARRAY_3D_FREE( FEDU );
    
    
                /*
                 *  Get FEDU_Energy variable.
                 */
                LGM_ARRAY_1D( Flux_Energy_Tmp, Nk, double );
                LGM_ARRAY_1D( Flux_Energy, Flux_nEnergy, double );
                cdf_status = CDFgetzVarAllRecordsByVarID( id, CDFgetVarNum( id, "FEDU_Energy"), &Flux_Energy_Tmp[0] );
                if (cdf_status != CDF_OK) StatusHandler( cdf_status );
                for (k=0; k<Flux_nEnergy; k++ ) {
                    Flux_Energy[k] = Flux_Energy_Tmp[k];
                    printf("Flux_Energy[%d] = %lf\n", k, Flux_Energy[k]);
                }
                LGM_ARRAY_1D_FREE( Flux_Energy_Tmp );
    
    
    
                /*
                 *  Get FEDU_Alpha variable.
                 */
                LGM_ARRAY_1D( Flux_Alpha_Tmp, Nj, double );
                LGM_ARRAY_1D( Flux_Alpha, Flux_nAlpha, double );
                cdf_status = CDFgetzVarAllRecordsByVarID( id, CDFgetVarNum( id, "FEDU_Alpha"), &Flux_Alpha_Tmp[0] );
                if (cdf_status != CDF_OK) StatusHandler( cdf_status );
                for (j=0; j<Flux_nAlpha; j++ ) {
                    Flux_Alpha[j] = Flux_Alpha_Tmp[j];
                    printf("Flux_Alpha[%d] = %lf\n", j, Flux_Alpha[j]);
                }
                LGM_ARRAY_1D_FREE( Flux_Alpha_Tmp );
    
    
                /*
                 *  Get Epoch variable.
                 */
                i = CDFgetVarNum( id, "Epoch");
                LGM_ARRAY_1D( Epoch, nTimes, double );
                cdf_status = CDFgetzVarAllRecordsByVarID( id, i, &Epoch[0] );
                if (cdf_status != CDF_OK) StatusHandler( cdf_status );
                LGM_ARRAY_2D( IsoTimes, nTimes, 30, char );
                for (i=0; i<nTimes; i++){
                    CDF_TT2000_to_UTC_string( CDF_TT2000_from_UTC_EPOCH(Epoch[i]), &IsoTimes[i][0], 3 );
                    //printf("Time = %lf Str = %s\n", Epoch[i], IsoTimes[i] );
                }
    
    
                /*
                 *  Get Position variable. Assumed in km from CDF file.
                 */
                i = CDFgetVarNum( id, "Position");
                LGM_ARRAY_2D( Position, nTimes, 3, double );
                cdf_status = CDFgetzVarAllRecordsByVarID( id, i, &Position[0][0] );
                if (cdf_status != CDF_OK) StatusHandler( cdf_status );
    
    
                /*
                 *  Get FEDU_Lstar
                 */
                MagEphem_nAlpha = Flux_nAlpha;
                LGM_ARRAY_2D( FEDU_Lstar_Tmp, Ni, Nj, double );
                LGM_ARRAY_2D( MagEphem_Lstar, Ni, MagEphem_nAlpha, double );
    printf("Ni, Nj = %d %d MagEphem_nAlpha = %d\n", Ni, Nj, MagEphem_nAlpha);
                LGM_ARRAY_1D( MagEphem_Alpha, MagEphem_nAlpha, double );
                i = CDFgetVarNum( id, "FEDU_Lstar");
                cdf_status = CDFgetzVarAllRecordsByVarID( id, i, &FEDU_Lstar_Tmp[0][0] );
                if (cdf_status != CDF_OK) StatusHandler( cdf_status );
                for (j=0; j<MagEphem_nAlpha; j++){
                    MagEphem_Alpha[j] = Flux_Alpha[j];
                }
    
                for (i=0; i<nTimes; i++){
                    for (j=0; j<MagEphem_nAlpha; j++){
                        if ( FEDU_Lstar_Tmp[i][j] < 0.9 ) {
                            MagEphem_Lstar[i][j] = LGM_FILL_VALUE;
                            }
                        else {
                            MagEphem_Lstar[i][j] = FEDU_Lstar_Tmp[i][j];
                            }
                    }
                }
                LGM_ARRAY_2D_FREE( FEDU_Lstar_Tmp );
    
    
                cdf_status = CDFcloseCDF( id );
                if ( cdf_status != CDF_OK ) StatusHandler( cdf_status );
                }
//exit(0);

                if ( !FileExists || Force || AppendMode ) {


                  /*
                   * Open Output HDF5 file for writing/appending and
                   * Create variables.
                   */
                  if ( !Force && AppendMode && FileExists ) {

                    /* Get the time variable from the HDF5 file to know what has already been processed */
                    file = H5Fopen( OutFile, H5F_ACC_RDWR, H5P_DEFAULT );
                    ProcessedTime = Get_StringDataset_1D( file, "IsoTime", Dims );
                    nProcTimes = Dims[0];

                  } else {

                    file = H5Fcreate( OutFile, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
                    WriteHDF5PSDHeader(file, nMu, nK, Mu, K);

                  }


                  /*
                   *   Allocate  space for J_DIFF
                   Lgm_F2P_SetFlux( J[k], Flux_Energy, Flux_nEnergy, Flux_Alpha, Flux_nAlpha, f2p );
                  */
                  LGM_ARRAY_2D( J_DIFF, Flux_nEnergy, Flux_nAlpha, double );
                  H5_nT = 0;
                  for (k=0; k<nTimes; k++){

                    // Flux -> PSD
                    f2p = Lgm_F2P_CreateFluxToPsd(1);
                    if (fitMaxw==0) {
                        f2p->FitType = LGM_F2P_MAXWELLIAN;
                        f2p->nMaxwellians = numMaxw;
                    } else {
                        f2p->FitType = LGM_F2P_SPLINE;
                    }
                    f2p->Extrapolate  = 1;
                    f2p->DumpDiagnostics = 0;

                    IsoTimeStringToDateTime( IsoTimes[k], &d, c );
                    printf("\tDate/Time: %s\n", &IsoTimes[k][0] );
                    Lgm_Set_Coord_Transforms( d.Date, d.Time, c );
                    Lgm_Doy( d.Date, &cYear, &cMonth, &cDay, &cDoy);
                    cJD = Lgm_JD( cYear, cMonth, cDay, d.Time, LGM_TIME_SYS_UTC, c );


                    // Determine whether time has been processed
                    if ( !Force && AppendMode && FileExists ) {

                      // Assume not processed
                      AppendFlag = 1;
                      // Loop through all times
                      for(iTime=0; iTime<nProcTimes; iTime++) {
                        // Check current time versus iterated processed time
                        if ( strcmp(ProcessedTime[iTime], IsoTimes[k]) == 0) {
                          // Match found - Do not append
                          AppendFlag = 0;
                        }
                      }

                      // Skip to next timestep
                      if ( !AppendFlag ) {

                        // Increment the data row so that the HDF5 appends the data correctly
                        H5_nT++;
                        continue;
                      }

                    }

/*
                    // interp the GSM position to this time. IS THIS CORRECT?
                    v.x = GSM[k][0]; v.y = GSM[k][1]; v.z = GSM[k][2];
                    printf("\tS/C GSM Position: [ %g, %g, %g ] Re\n", v.x, v.y, v.z );
                    Lgm_F2P_SetDateTimeAndPos( &d, &v, f2p );
*/

                    if ( Position[k][0] > LGM_FILL_VALUE ) {
                        if (fromHDF == 1) {
                            v.x = Position[k][0]; v.y = Position[k][1]; v.z = Position[k][2];
                        } else { //FIXME: Kluge to make CDF and HDF5 do the same...
                            v_geo.x = Position[k][0]/Re; v_geo.y = Position[k][1]/Re; v_geo.z = Position[k][2]/Re;
                            Lgm_Convert_Coords( &v_geo, &v, GEO_TO_GSM, c );
                        }
                        printf("\tS/C GSM Position: [ %g, %g, %g ] Re\n", v.x, v.y, v.z );


                        // We need to compute the AlphaOfK array.
                        Lgm_Make_UTC( d.Date, d.Time, &UTC_DateTime, c );

                        if ( FixModelDateTime ) {
                            Lgm_F2P_SetDateTimeAndPos( &ModelDateTime, &v, f2p );
                            Lgm_get_QinDenton_at_JD( ModelDateTime.JD, &QinDen, 1, 1 );
                            Lgm_set_QinDenton( &QinDen, mInfo );
                        } else {
                            Lgm_F2P_SetDateTimeAndPos( &d, &v, f2p );
                            Lgm_get_QinDenton_at_JD( d.JD, &QinDen, 1, 1 );
                            Lgm_set_QinDenton( &QinDen, mInfo );
                        }
                
                        Lgm_Set_Coord_Transforms( d.Date, d.Time, mInfo->c );

                        if ( OverRideKp ) {
                            mInfo->Kp = T89Q_Kp;
                        }


                        if (arguments.verbose > 0) {
                            printf("\tKp: %d\n", mInfo->Kp );
                        }
                        // Add diff. flux data/info to f2p structure.
                        for (iii=0; iii<Flux_nEnergy; iii++ ){
                            for (jjj=0; jjj<Flux_nAlpha; jjj++ ){
                                J_DIFF[iii][jjj] = J[k][jjj][iii]; ////////
                            }
                        }


                        /*
                         * This call will populate the f->AofK[] array.
                         */
                        Lgm_F2P_SetFlux( J_DIFF, Flux_Energy, Flux_nEnergy, Flux_Alpha, Flux_nAlpha, f2p );



                        //  Compute Phase Space Density at the constant Mu's and K's
                        mInfo->UseInterpRoutines = TRUE;
                        Lgm_F2P_GetPsdAtConstMusAndKs( Mu, nMu, K, nK, mInfo, f2p );


                        if (arguments.verbose >= 2) {
                            printf("\tPSD Array (K->Columns, Mu->rows):\n");
                            for ( iii=0; iii<nMu; iii++ ){
                                printf("\t  [");
                                for (jjj=0; jjj<nK; jjj++ ) {
                                    printf("%12g ", f2p->PSD_MK[iii][jjj]);
                                }
                                printf("]\n");
                            }

                            printf("\n\n");
                        }

                        /*
                         * Compute the PsdLstar arrays. We already have the Lstar
                         * values that were computed by MagEphem. But these are not
                         * at the PAs (i.e. Ks) that we need them to be at. So lets
                         * interpolate the MagEphem arrays to the Ks that we need.
                         */
                        //int kkk;
                        for (kk=0; kk<nK; kk++){
                            Lgm_InterpArr( MagEphem_Alpha, MagEphem_Lstar[k], MagEphem_nAlpha,   f2p->AofK[kk], &PsdLstar[kk] );
                            AlphaOfK[kk] = ( f2p->AofK[kk] > 0.0 ) ? f2p->AofK[kk] : LGM_FILL_VALUE;
//printf("f2p->AofK[%d] = %g\n", kk, f2p->AofK[kk] );
                            //for(kkk=0; kkk< MagEphem_nAlpha; ++kkk){
                                //printf("MagEphem_Alpha[%d], MagEphem_Lstar[%d][%d] = %g %g\n", kkk, k, kkk, MagEphem_Alpha[kkk], MagEphem_Lstar[k][kkk] );
                            //}
                            //printf("******** f2p->AofK[%d], PsdLstar[%d] = %g %g\n", kk, f2p->AofK[kk], kk, PsdLstar[kk] );
                        }

                        WriteHDF5PSDData(file, H5_nT, nMu, nK, IsoTimes,
                                         PsdLstar, AlphaOfK, f2p->PSD_MK);


                    } else {

                        for (kk=0; kk<nK; kk++){
                            PsdLstar[kk] = LGM_FILL_VALUE;
                        }

                        // Write PSD FILL VALS.
                        double  **PSD_MK_FILL;
                        LGM_ARRAY_2D( PSD_MK_FILL, nMu, nK, double );
                        for ( iii=0; iii<nMu; iii++ ){
                            for (jjj=0; jjj<nK; jjj++ ) {
                                PSD_MK_FILL[iii][jjj] = LGM_FILL_VALUE;
                            }
                        }

                        WriteHDF5PSDData(file, H5_nT, nMu, nK, IsoTimes,
                                         PsdLstar, AlphaOfK, PSD_MK_FILL);

                        LGM_ARRAY_2D_FREE( PSD_MK_FILL );

                    }

                    ++H5_nT;

                    Lgm_F2P_FreeFluxToPsd( f2p );

                }

                status = H5Fclose( file );


                LGM_ARRAY_2D_FREE( J_DIFF );

            } else { // OutFile exists already
                printf("\t\tOutput file already exists: %s \n", OutFile );
            }

            if ( !Force && AppendMode && FileExists ) {
              LGM_ARRAY_2D_FREE( ProcessedTime );
            }
            LGM_ARRAY_2D_FREE( IsoTimes );
            LGM_ARRAY_1D_FREE( Epoch );
            LGM_ARRAY_2D_FREE( Position );
//            LGM_ARRAY_2D_FREE( GSM );
            LGM_ARRAY_2D_FREE( MagEphem_Lstar );
            LGM_ARRAY_1D_FREE( MagEphem_Alpha );
            LGM_ARRAY_1D_FREE( Flux_Alpha );
            LGM_ARRAY_1D_FREE( Flux_Energy );
            LGM_ARRAY_3D_FREE( J );



        } else { // error on opening InFile
            printf("\t\tCould not open input file for reading: %s \n", InFile );
        }


    } // End date loop

    LGM_ARRAY_1D_FREE( Mu );
    LGM_ARRAY_1D_FREE( K );
    LGM_ARRAY_1D_FREE( AlphaOfK );
    LGM_ARRAY_1D_FREE( LstarOfK );
    LGM_ARRAY_1D_FREE( PsdLstar );



    free( OutFile );
    free( InFile );

    Lgm_free_ctrans( c );
    Lgm_FreeMagInfo( mInfo );
    Lgm_destroy_eop( e );


    return(0);
}

void InterpVec( double *JD, Lgm_Vector *Ugsm, int n1, double cJD, Lgm_Vector *v ) {

    int                 i;
    gsl_interp_accel    *acc;
    gsl_spline          *spline;

    double              *y = (double *)calloc( n1+1, sizeof(double) );


    // X-Comp
    acc    = gsl_interp_accel_alloc( );
    spline = gsl_spline_alloc( gsl_interp_akima, n1 );
    for (i=0; i<n1; i++) y[i] = Ugsm[i].x;
    gsl_spline_init( spline, JD, y, n1 );
    v->x = gsl_spline_eval( spline, cJD, acc );

    // Y-Comp
    acc    = gsl_interp_accel_alloc( );
    spline = gsl_spline_alloc( gsl_interp_akima, n1 );
    for (i=0; i<n1; i++) y[i] = Ugsm[i].y;
    gsl_spline_init( spline, JD, y, n1 );
    v->y = gsl_spline_eval( spline, cJD, acc );

    // Z-Comp
    acc    = gsl_interp_accel_alloc( );
    spline = gsl_spline_alloc( gsl_interp_akima, n1 );
    for (i=0; i<n1; i++) y[i] = Ugsm[i].z;
    gsl_spline_init( spline, JD, y, n1 );
    v->z = gsl_spline_eval( spline, cJD, acc );


    free( y );
    gsl_spline_free( spline );
    gsl_interp_accel_free( acc );

    return;

}


void InterpArr( double *xa, double *ya, int n, double x, double *y ) {

    gsl_interp_accel    *acc;
    gsl_spline          *spline;
    double              *xa2, *ya2;
    int                 i, Flag;

    Flag = 0;
    if ( xa[1] < xa[0] ) {

        xa2 = (double *)calloc( n, sizeof(double) );
        ya2 = (double *)calloc( n, sizeof(double) );

        for (i=0; i<n; i++){
            xa2[i] = xa[n-1-i];
            ya2[i] = ya[n-1-i];
        }

        Flag = 1;

    } else {

        xa2 = xa;
        ya2 = ya;

    }


    acc    = gsl_interp_accel_alloc( );
    spline = gsl_spline_alloc( gsl_interp_akima, n );
    gsl_spline_init( spline, xa2, ya2, n );
    *y = gsl_spline_eval( spline, x, acc );
    gsl_spline_free( spline );
    gsl_interp_accel_free( acc );

    if ( Flag ){
        free( xa2 );
        free( ya2 );
    }

    return;

}


void WriteHDF5PSDHeader(hid_t file, int nMu, int nK, double *Mu, double *K) {

  hid_t atype;
  hid_t DataSet;
  hid_t status;
  hid_t space;
  char buf[80];

  // Create Phase Space Density dataset
  DataSet = CreateExtendableRank3DataSet( file, "PhaseSpaceDensity", nMu, nK,
                                          H5T_NATIVE_DOUBLE, &space );
  Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Phase Space Density.");
  Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
  Lgm_WriteStringAttr( DataSet, "DEPEND_1",   "Mu" );
  Lgm_WriteStringAttr( DataSet, "DEPEND_2",   "K" );
  Lgm_WriteStringAttr( DataSet, "UNITS",      "(c/cm/MeV)^3" );
  Lgm_WriteStringAttr( DataSet, "SCALETYP",   "log" );
  Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
  Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
  Lgm_WriteStringAttr( DataSet, "SCALEMIN",   "1e-10" );
  Lgm_WriteStringAttr( DataSet, "SCALEMAX",   "1.0" );
  Lgm_WriteStringAttr( DataSet, "VALIDMIN",   "1e-30" );
  Lgm_WriteStringAttr( DataSet, "VALIDMAX",   "1.0" );
  status  = H5Sclose( space );
  status  = H5Dclose( DataSet );

  // Create Lstar dataset
  DataSet = CreateExtendableRank2DataSet( file, "Lstar", nK,
                                          H5T_NATIVE_DOUBLE, &space );
  Lgm_WriteStringAttr( DataSet, "DESCRIPTION",
                       "Generalized Roederer L-shell value (also known as L*).");
  Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
  Lgm_WriteStringAttr( DataSet, "DEPEND_1",   "K" );
  Lgm_WriteStringAttr( DataSet, "UNITS",      "Dimensionless" );
  Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
  Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
  Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
  Lgm_WriteStringAttr( DataSet, "SCALEMIN",   "1.0" );
  Lgm_WriteStringAttr( DataSet, "SCALEMAX",   "10.0" );
  Lgm_WriteStringAttr( DataSet, "VALIDMIN",   "1.0" );
  Lgm_WriteStringAttr( DataSet, "VALIDMAX",   "100.0" );
  status  = H5Sclose( space );
  status  = H5Dclose( DataSet );

  // Create AlphaEq dataset
  DataSet = CreateExtendableRank2DataSet( file, "AlphaEq", nK,
                                          H5T_NATIVE_DOUBLE, &space );
  Lgm_WriteStringAttr( DataSet, "DESCRIPTION",
                       "Equatorial pitch angle corresponding to each of the K values.");
  Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
  Lgm_WriteStringAttr( DataSet, "DEPEND_1",   "K" );
  Lgm_WriteStringAttr( DataSet, "UNITS",      "Degrees" );
  Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
  Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
  Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
  Lgm_WriteStringAttr( DataSet, "SCALEMIN",   "0.0" );
  Lgm_WriteStringAttr( DataSet, "SCALEMAX",   "90.0" );
  Lgm_WriteStringAttr( DataSet, "VALIDMIN",   "0.0" );
  Lgm_WriteStringAttr( DataSet, "VALIDMAX",   "90.0" );
  status  = H5Sclose( space );
  status  = H5Dclose( DataSet );

  // Create Times dataset
  atype = CreateStrType(32);
  DataSet = CreateExtendableRank1DataSet( file, "IsoTime", atype, &space);
  status = H5Tclose( atype );
  status = H5Sclose( space );
  status = H5Dclose( DataSet );

  // Create Mu dataset
  DataSet = CreateSimpleRank1DataSet( file, "Mu", nMu,
                                      H5T_NATIVE_DOUBLE, &space );
  Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Mu, 1st adiabatic invariant.");
  Lgm_WriteStringAttr( DataSet, "UNITS",     "MeV/G" );
  Lgm_WriteStringAttr( DataSet, "SCALETYP",  "log" );
  Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
  Lgm_WriteStringAttr( DataSet, "VAR_TYPE",  "data" );
  sprintf(buf, "%g", Mu[0] );     Lgm_WriteStringAttr( DataSet, "VALIDMIN",  buf );
  sprintf(buf, "%g", Mu[nMu-1] ); Lgm_WriteStringAttr( DataSet, "VALIDMAX",  buf );
  sprintf(buf, "%g", Mu[0] );     Lgm_WriteStringAttr( DataSet, "SCALEMIN",  buf );
  sprintf(buf, "%g", Mu[nMu-1] ); Lgm_WriteStringAttr( DataSet, "SCALEMAX",  buf );
  status  = H5Sclose( space );
  status  = H5Dclose( DataSet );
  LGM_HDF5_WRITE_SIMPLE_RANK1_DATASET( file, "Mu", nMu, H5T_NATIVE_DOUBLE, &Mu[0] );

  // Create K dataset
  DataSet = CreateSimpleRank1DataSet( file, "K", nK,
                                      H5T_NATIVE_DOUBLE, &space );
  Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "K, 2nd adiabatic invariant.");
  Lgm_WriteStringAttr( DataSet, "UNITS",     "Re G^.5" );
  Lgm_WriteStringAttr( DataSet, "SCALETYP",  "log" );
  Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
  Lgm_WriteStringAttr( DataSet, "VAR_TYPE",  "data" );
  sprintf(buf, "%g", K[0] );    Lgm_WriteStringAttr( DataSet, "VALIDMIN",  buf );
  sprintf(buf, "%g", K[nK-1] ); Lgm_WriteStringAttr( DataSet, "VALIDMAX",  buf );
  sprintf(buf, "%g", K[0] );    Lgm_WriteStringAttr( DataSet, "SCALEMIN",  buf );
  sprintf(buf, "%g", K[nK-1] ); Lgm_WriteStringAttr( DataSet, "SCALEMAX",  buf );
  status = H5Sclose( space );
  status = H5Dclose( DataSet );
  LGM_HDF5_WRITE_SIMPLE_RANK1_DATASET( file, "K", nK, H5T_NATIVE_DOUBLE, &K[0] );

}


void WriteHDF5PSDData(hid_t file, int i, int nMu, int nK, char **IsoTimes,
                      double *PsdLstar, double *AlphaOfK, double **PSD) {

  // Write IsoTime
  hid_t atype = CreateStrType(32);
  LGM_HDF5_EXTEND_RANK1_DATASET(file, "IsoTime", i,
                                atype, &IsoTimes[i][0]);

  // Write Lstar
  atype = H5T_NATIVE_DOUBLE;
  LGM_HDF5_EXTEND_RANK2_DATASET(file, "Lstar", i, nK, atype, &PsdLstar[0]);

  // Write AlphaEq
  LGM_HDF5_EXTEND_RANK2_DATASET(file, "AlphaEq", i, nK, atype, &AlphaOfK[0]);

  // Write Phase Space Density
  LGM_HDF5_EXTEND_RANK3_DATASET(file, "PhaseSpaceDensity", i, nMu, nK,
                                atype, &PSD[0][0]);

}
