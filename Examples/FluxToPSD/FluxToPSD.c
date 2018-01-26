/*
CHECK ME FOR MEMORY LEAKS.....
THIS ALSO needs some argument changes...
e.g. The Mu's K's instead of pitch angles...
*/
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
#include <Lgm_DynamicMemory.h>
#include <Lgm_MagModelInfo.h>
#include <Lgm_FluxToPsd.h>
#include <Lgm_Eop.h>
#include <Lgm_QinDenton.h>
#include <Lgm_Misc.h>
#include <Lgm_HDF5.h>



#define KP_DEFAULT 0
void InterpVec( double *JD, Lgm_Vector *Ugsm, int n1, double cJD, Lgm_Vector *v );
void InterpArr( double *JD, double *y, int n, double cJD, double *r ) ;

const  char *argp_program_version     = "FluxToPSD 1.0";
const  char *argp_program_bug_address = "<mghenderson@lanl.gov>";
static char doc[] = "Converts Flux to PSD. Blah Blah";

// Mandatory arguments
#define     nArgs   2
static char ArgsDoc[] = "StdFluxFile StdPsdFile";


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
    {"Mus",             'M',    "\"start_Mu, end_Mu, nMu\"",  0,                                      "Mu values to compute. Default is \"1, 2000.0, 18\"." },
    {"Ks",              'K',    "\"start_K, end_K, nK\"",     0,                                      "K values to compute. Default is \"0.01, 10.0, 18\"." },
    {"FootPointHeight", 'f',    "height",                     0,                                      "Footpoint height in km. Default is 100km."                  },
    {"Quality",         'q',    "quality",                    0,                                      "Quality to use for L* calculations. Default is 3."      },
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
    char        *args[ nArgs ];
    int         silent;
    int         verbose;

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
    Lgm_Vector       U;
    Lgm_DateTime     UTC;
    Lgm_Eop          *e = Lgm_init_eop( 0 );
    Lgm_EopOne       eop;
    int              i, iii, jjj;
    char             InputFilename[1024];
    char             OutputFilename[1024];
    char             IntModel[20], ExtModel[20];
    int              AppendMode, UseEop, Colorize, Force, nv;
    double           Inc, Alpha[1000], FootpointHeight;
    int              nAlpha, Quality;
    long int         StartDate, EndDate, Date;
    int              sYear, sMonth, sDay, sDoy, eYear, eMonth, eDay, eDoy, Year, Month, Day;
    int              cYear, cMonth, cDay, cDoy;
    double           sJD, eJD, JD, Time, cJD;
    Lgm_QinDentonOne p;
    hsize_t           dims[3];
    double           *Flux_Energy;
    double           *Flux_Alpha;
    int               Flux_nEnergy;
    int               Flux_nAlpha;
    double           **GSM;
    double           **MagEphem_Lstar;
    double           *MagEphem_Alpha;
    double           ***J;



    /*
     * App-specific declartations
     */
//    Lgm_CTrans      *c = Lgm_init_ctrans( 0 );
    Lgm_DateTime     d;
    Lgm_Vector       u, v, Ugsm[300];
    Lgm_FluxToPsd   *f2p;
    Lgm_PsdToFlux   *p2f;
    Lgm_DateTime     UTC_DateTime;
    Lgm_MagModelInfo *mInfo, *mInfo2;
    Lgm_QinDentonOne QinDen;
    double           *E, *A, **J_DIFF;
    int              nE, nA, j, n, k, n1, nTimes, MagEphem_nAlpha;
    double           *Mu, *K, *LstarOfK, *AlphaOfK, *PsdLstar;
    int              nMu, nK;
    double           f, p2c2, df, sa;
    double           S[100], a, b, r;
    char             Line[5000], **IsoTimes;
    double           MagEphem_JD1[300];
    long int         pos;
    FILE             *fp;

    hid_t            file;
    hid_t            space;

    hid_t            TimeDataSet, TimeMemSpace;
    hsize_t          TimeOffset[1]   = { 0 };
    hsize_t          TimeDims[1]     = { 1 };
    hsize_t          TimeSlabSize[1] = { 1 };

    hid_t            LstarDataSet, LstarMemSpace;
    hsize_t          LstarOffset[2]   = { 0,  0, };
    hsize_t          LstarDims[2]     = { 1, 18, };
    hsize_t          LstarSlabSize[2] = { 1, 18, };

    hid_t            PsdDataSet, PsdMemSpace;
    hsize_t          PsdOffset[3]   = { 0,  0,  0 };
    hsize_t          PsdDims[3]     = { 1, 18, 18 };
    hsize_t          PsdSlabSize[3] = { 1, 18, 18};

    herr_t           status;



    /*
     * Default option values.
     */
    arguments.StartMu         = 1.0;    // start Mu [MeV/G]
    arguments.EndMu           = 2000.0; // stop Mu [MeV/G]
    arguments.nMu             = 18;     // 18 Mus
    arguments.StartK          = 0.01;   // start K [Re G^.5]
    arguments.EndK            = 10.0;   // stop K [Re G^.5]
    arguments.nK              = 18;     // 18 K
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
     *  Set input and output filenames
     */
    strcpy( InputFilename, arguments.args[0] );
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


    mInfo = Lgm_InitMagInfo( );

    if ( !strcmp( ExtModel, "T87" ) ){
        mInfo->Bfield = Lgm_B_T87;
    } else if ( !strcmp( ExtModel, "CDIP" ) ){
        mInfo->Bfield = Lgm_B_cdip;
    } else if ( !strcmp( ExtModel, "EDIP" ) ){
        mInfo->Bfield = Lgm_B_edip;
    } else if ( !strcmp( ExtModel, "IGRF" ) ){
        mInfo->Bfield = Lgm_B_igrf;
    } else { //if ( !strcmp( ExtModel, "T89" ) ){
        // default
        mInfo->Bfield = Lgm_B_T89;
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
        sprintf( Command, "mkdir -p %s", BaseDir); system( Command );
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

        InFileExists = ( (StatError = stat( InFile, &StatBuf )) != -1 ) ? TRUE : FALSE;

       /*
        * Read in a combined Flux/MagEphem HDF5 File
        */
        if ( InFileExists ) {

            // Read the InFile
            file = H5Fopen( InFile, H5F_ACC_RDONLY, H5P_DEFAULT );

            //Read ISO time
            IsoTimes = Get_StringDataset_1D( file, "/Epoch", dims);

            // Read GSM and Lstar position
            GSM  = Get_DoubleDataset_2D( file, "/Ephemeris/GSM", dims);
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


            if ( !FileExists || Force ) {

                printf("\t\tDo work here to create file.\n");

                /*
                 *   Create HDF5 file.
                 */
                PsdDims[0] = nTimes; PsdDims[1] = nMu; PsdDims[2] = nK;
                PsdSlabSize[0] = 1; PsdSlabSize[1] = nMu; PsdSlabSize[2] = nK;
                file       = H5Fcreate( OutFile, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
                space      = H5Screate_simple( 3, PsdDims, NULL ); // rank 3
                PsdDataSet = H5Dcreate( file, "PhaseSpaceDensity", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                status     = H5Sclose( space );

                LstarDims[0] = nTimes; LstarDims[1] = nK;
                LstarSlabSize[0] = 1; LstarSlabSize[1] = nK;
                space        = H5Screate_simple( 2, LstarDims, NULL ); // rank 2
                LstarDataSet = H5Dcreate( file, "Lstar", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                status       = H5Sclose( space );

                TimeDims[0] = nTimes;
                space        = H5Screate_simple( 1, TimeDims, NULL ); // rank 1
                hid_t atype = H5Tcopy (H5T_C_S1);
                status = H5Tset_size( atype, strlen(IsoTimes[0]) );
                status  = H5Tset_strpad( atype, H5T_STR_NULLPAD );
                status  = H5Tset_cset( atype, H5T_CSET_ASCII );
                TimeDataSet  = H5Dcreate( file, "IsoTime", atype, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                status       = H5Sclose( space );



                PsdMemSpace   = H5Screate_simple( 3, PsdSlabSize, NULL );
                LstarMemSpace = H5Screate_simple( 2, LstarSlabSize, NULL );
                TimeMemSpace  = H5Screate_simple( 1, TimeSlabSize, NULL );


                /*
                 *   Allocate  space for J_DIFF
                    Lgm_F2P_SetFlux( J[k], Flux_Energy, Flux_nEnergy, Flux_Alpha, Flux_nAlpha, f2p );
                 */
                LGM_ARRAY_2D( J_DIFF, Flux_nEnergy, Flux_nAlpha, double );
                for (k=0; k<nTimes; k++){

                    // Flux -> PSD
                    f2p = Lgm_F2P_CreateFluxToPsd(1);
                    f2p->Extrapolate = 3;

                    IsoTimeStringToDateTime( IsoTimes[k], &d, c );
                    printf("\t\tDate/Time: %s\n", IsoTimes[k] );
                    Lgm_Set_Coord_Transforms( d.Date, d.Time, c );
                    Lgm_Doy( d.Date, &cYear, &cMonth, &cDay, &cDoy);
                    cJD = Lgm_JD( cYear, cMonth, cDay, d.Time, LGM_TIME_SYS_UTC, c );

                    // interp the GSM position to this time. IS THIS CORRECT?
                    v.x = GSM[k][0]; v.y = GSM[k][1]; v.z = GSM[k][2];
                    printf("\t\tv_gsm = %g %g %g\n", v.x, v.y, v.z );
                    Lgm_F2P_SetDateTimeAndPos( &d, &v, f2p );


                    // We need to compute the AlphaOfK array.
//                    Lgm_Make_UTC( d->Date, d->Time, &UTC_DateTime, c );

                    Lgm_get_QinDenton_at_JD( d.JD, &QinDen, 0 , 0);
                    Lgm_set_QinDenton( &QinDen, mInfo );
                    Lgm_Set_Coord_Transforms( d.Date, d.Time, mInfo->c );
                    Lgm_Setup_AlphaOfK( &d, &v, mInfo );
                    #pragma omp parallel private(mInfo2)
                    #pragma omp for schedule(dynamic, 1)
                    for ( jjj=0; jjj<nK; jjj++ ) {
                        mInfo2 = Lgm_CopyMagInfo( mInfo );  // make a private (per-thread) copy of mInfo
                        AlphaOfK[jjj] = Lgm_AlphaOfK( K[jjj], mInfo2 );
                        printf("AlphaOfK[%d] = %g\n", jjj, AlphaOfK[jjj] );
                        Lgm_FreeMagInfo( mInfo2 ); // free mInfo2
                    }
                    Lgm_TearDown_AlphaOfK( mInfo );

                    // Now, from the AlphaOfK values, we need to find L*
                    // from the L* versus alpha array that was computed by
                    // MagEphem.
                    for ( jjj=0; jjj<nK; jjj++ ) {
                        InterpArr( MagEphem_Alpha, MagEphem_Lstar[k], MagEphem_nAlpha,    AlphaOfK[jjj], &PsdLstar[jjj] );
                        printf("PsdLstar[%d] = %g\n", jjj, PsdLstar[jjj]);
                    }



                    // Add diff. flux data/info to f2p structure.
                    for (iii=0; iii<Flux_nEnergy; iii++ ){
                        for (jjj=0; jjj<Flux_nAlpha; jjj++ ){
                            J_DIFF[iii][jjj] = J[k][jjj][iii];
                        }
                    }
                    //Lgm_F2P_SetFlux( J[k], Flux_Energy, Flux_nEnergy, Flux_Alpha, Flux_nAlpha, f2p );
                    Lgm_F2P_SetFlux( J_DIFF, Flux_Energy, Flux_nEnergy, Flux_Alpha, Flux_nAlpha, f2p );

                    //  Compute Phase Space Density at the constant Mu's and K's
                    Lgm_F2P_GetPsdAtConstMusAndKs( Mu, nMu, K, nK, mInfo, f2p );


                    for ( iii=0; iii<nMu; iii++ ){
                        printf("[");
                        for (jjj=0; jjj<nK; jjj++ ) {
                            printf("%12g ", f2p->PSD_MK[iii][jjj]);
                        }
                        printf("]\n");
                    }
                        /*
                        */

//TODO;
// need Mu and K arrays.
// Others?
                    // Write ISoTime Strings
                    TimeOffset[0] = k; //TimeOffset[1] = 0;
                    space  = H5Dget_space( TimeDataSet );
                    status = H5Sselect_hyperslab( space, H5S_SELECT_SET, TimeOffset, NULL, TimeSlabSize, NULL );
                    status = H5Dwrite( TimeDataSet, atype, TimeMemSpace, space, H5P_DEFAULT, IsoTimes[k] );
                    status = H5Sclose(space);


                    // Write Lstar values
                    LstarOffset[0] = k; LstarOffset[1] = 0;
                    space  = H5Dget_space( LstarDataSet );
                    status = H5Sselect_hyperslab( space, H5S_SELECT_SET, LstarOffset, NULL, LstarSlabSize, NULL );
                    status = H5Dwrite( LstarDataSet, H5T_NATIVE_DOUBLE, LstarMemSpace, space, H5P_DEFAULT, PsdLstar );
                    status = H5Sclose(space);

                    // Write PSD values.
                    PsdOffset[0] = k; PsdOffset[1] = 0; PsdOffset[2] = 0;
                    space  = H5Dget_space( PsdDataSet );
                    status = H5Sselect_hyperslab( space, H5S_SELECT_SET, PsdOffset, NULL, PsdSlabSize, NULL );
                    status = H5Dwrite( PsdDataSet, H5T_NATIVE_DOUBLE, PsdMemSpace, space, H5P_DEFAULT, &f2p->PSD_MK[0][0] );
                    status = H5Sclose(space);


                    Lgm_F2P_FreeFluxToPsd( f2p );

                }


                status = H5Sclose( PsdMemSpace );
                status = H5Dclose( PsdDataSet );
                status = H5Sclose( LstarMemSpace );
                status = H5Dclose( LstarDataSet );
                status = H5Sclose( TimeMemSpace );
                status = H5Dclose( TimeDataSet );
                status = H5Fclose( file );


                LGM_ARRAY_2D_FREE( J_DIFF );

            } else { // OutFile exists already
                printf("\t\tOutput file already exists: %s \n", OutFile );
            }


            LGM_ARRAY_2D_FREE( IsoTimes );
            LGM_ARRAY_2D_FREE( GSM );
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








