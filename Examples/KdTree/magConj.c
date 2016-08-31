#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <libgen.h>
#include <glob.h>
#include <argp.h>
#include <regex.h>
#include <Lgm_DynamicMemory.h>
#include <Lgm_CTrans.h>
#include <Lgm_KdTree.h>
#include <Lgm_Misc.h>
#include <Lgm_HDF5.h>

void StringSplit( char *Str, char *StrArray[], int len, int *n );


const  char *ProgramName = "magConj";
const  char *argp_program_version     = "magConj_0.1";
const  char *argp_program_bug_address = "<smorley@lanl.gov>";
static char doc[] = "\nUsing magnetic ephemerii of S/C from input files, computes magnetic conjunctions"
                    " using a kd-tree storage structure and nearest neighbour algorithm.\n\n"
                    " Satellite names can be left as shell-type wildcards or can be replaced by '%B', both"
                    " together is not currently supported.\n"
                    " Sample usage with wildcards:\n\n\t ./magConj -v 3 -p "
                    "\"/Spacecraft/*/MagEphem/2014/20140126_*_MagEphem_T89c.h5\" "
                    "\n\t\t -T RBSPA /Spacecraft/RBSPA/MagEphem/2014/rbspa_MagEphem_T89D_20140126.h5\n"
                    "\t\t test_magConj.h5"
                    " \n\n"
                    " Sample usage with '%B' token:\n\n\t ./magConj -b \"ns61, ns62, ns63\" -p \n"
                    "\t\t \"/Spacecraft/%B/MagEphem/2014/20140126_%B_MagEphem_T89c.h5\" "
                    "\n\t\t -T RBSPA /Spacecraft/RBSPA/MagEphem/2014/rbspa_MagEphem_T89D_20140126.h5\n"
                    "\t\t test_magConj.h5"
                    " \n\n"
                    " Sample usage with fixed query point (requires dummy input filename):\n\n"
                    "\t ./magConj -Q \"1.0,6.0\" -s 2014-01-21T18:00:00 -e 2014-01-22T22:00:00 -c 20 -p "
                    "\t\t \"/mnt/projects/dream/Spacecraft/*/MagEphem/2014/*2014012[12]_*T89c*.h5\"\n"
                    "\t\t dum fixed_1_6_magConj.h5"
                    " \n\n";

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
    {"Pattern",         'p',    "pattern",           0,         "Glob pattern for finding input files"},
    {"Object",          'o',    "object",            0,         "Name of search object"},
    {"Bird",            'b',    "bird",              0,         "List of satellite names to substitute for %B token"},
    {"Target",          'T',    "target",            0,         "Name of target (query object)"},
    {"Distance",        'd',    "distance",          0,         "Maximum 1-D distance between leaf nodes and points for retrieval."},
    {"QueryPoint",      'Q',    "querypoint",        0,         "Query Point (K, L*), e.g. \"0.11,4.1\" "},
    {"StartDate",       's',    "yyyymmdd",          0,         "ISO StartDate "},
    {"EndDate",         'e',    "yyyymmdd",          0,         "ISO EndDate "},
    {"Cadence",         'c',    "cadence",           0,         "Cadence for query points, when using fixed target+start date+end date"},
    {"Force",           'F',    0,                   0,         "Overwrite output file even if it already exists"},
    {"verbose",         'v',    "verbosity",         0,         "Produce verbose output"},
    { 0 }
};


struct Arguments {
    char        *args[ nArgs ];
    char        Pattern[512];
    char        Object[512];
    char        Target[128];
    char        Birds[512];
    int         verbose;
    int         Force;
    double      Distance;

    char        QueryPoint[128];
    char        StartDate[128];
    char        EndDate[128];
    double      Cadence;
};


/* Parse a single option. */
static error_t parse_opt (int key, char *arg, struct argp_state *state) {

    /* Get the input argument from argp_parse, which we
      know is a pointer to our arguments structure. */
    struct Arguments *arguments = state->input;
    switch( key ) {
        case 's': // start date
            strncpy( arguments->StartDate, arg, 127 );
            break;
        case 'e': // end date
            strncpy( arguments->EndDate, arg, 127 );
            break;
        case 'c':
            sscanf( arg, "%lf", &arguments->Cadence );
            break;
        case 'Q':
            strncpy( arguments->QueryPoint, arg, 127 );
            break;
        case 'b':
            strncpy( arguments->Birds, arg, 511 );
            break;
        case 'p': // glob pattern
            strncpy( arguments->Pattern, arg, 511 );
            break;
        case 'o': // glob pattern
            strncpy( arguments->Object, arg, 511 );
            break;
        case 'd': // distance for search
            sscanf( arg, "%lf", &arguments->Distance );
            break;
        case 'T':
            strncpy( arguments->Target, arg, 127 );
            break;
        case 'F':
            arguments->Force = 1;
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


long int NearestPerObject(Lgm_KdTreeData *kNN, long int N, Lgm_KdTreeData *new_kNN, int verbose) {
    long int i, j, Nout;
    long int Id __attribute__((unused));
    int keep[N], DONE;
    char objectList[256][256] = {{ '\0' }};

    //loop over kNN in dist order (min->max)
    Nout = 0;
    for (i=N-1; i>=0; i--){
    //for (i=0; i<N; i++){
        Id = kNN[i].Id;
        //printf("Testing %dth entry (object = %s, position = %g %g %g)\n", i, kNN[i].Object, kNN[i].Position[0], kNN[i].Position[1]/(30.0*60.0*20.0), kNN[i].Position[2]/(30.0*60.0*4.0));
        j=0;
        DONE=0;
        //if object has member on stack, move on
        while (objectList[j][0] && DONE==0) {
            //printf("Testing %s for presence in list\n", objectList[j]);
            if (!strcmp((char *)kNN[i].Object, &objectList[j][0])) { //string is current item on list
                DONE = 1;
            }
            j++;
        } 
        //if object not on stack, add name to list of members and store point
        if (!DONE) {
             keep[i] = 1;
             strcpy(&objectList[Nout][0], (char *)kNN[i].Object);
             Nout++;
        }
        else {
            keep[i] = 0;
        }
    }
    if (!Nout) {
        if (verbose) printf("****** NearestPerObject: Nothing to keep ******\n");
        exit(-1);
        }
    if (verbose > 1) printf("\nNumber of kept entries = %ld\n", Nout);

    j = 0;
    for (i=0; i<N; i++) {
        if (keep[i]) {
            //printf("Keeping %dth entry (object = %s, position = %g %g %g)\n", i, kNN[i].Object, kNN[i].Position[0], kNN[i].Position[1]/(30.0*60.0*20.0), kNN[i].Position[2]/(30.0*60.0*4.0));
            new_kNN[j] = kNN[i];
            j++;
        }

    }
    return (Nout);
}


long int purgeObject(Lgm_KdTreeData *kNN, long int N, void *Object, Lgm_KdTreeData *new_kNN, int verbose) {
    long int i, j, Nout;
    long int Id __attribute__((unused));
    int keep[N];

    Nout = 0;
    for (i=0; i<N; i++){
        Id = kNN[i].Id;
        //printf("Testing %dth entry (object = %s, position = %g %g %g)\n", i, kNN[i].Object, kNN[i].Position[0], kNN[i].Position[1]/(30.0*60.0*10.0), kNN[i].Position[2]/(30.0*60.0*4.0));
        if (strcmp(kNN[i].Object, Object)) {
             keep[i] = 1;
             Nout++;
        }
        else {
            keep[i] = 0;
            //printf("Match: value (%s) != match (%s)\n", kNN[i].Object, Object);
        }
    }
    if (!Nout) {
        if (verbose) printf("****** purgeObject: Nothing to keep ******\n");
        kNN[0].Dist2=LGM_FILL_VALUE; //set empty entry...
        }
    if (verbose > 1) printf("\nNumber of non-matching entries = %ld\n", Nout);

    j = 0;
    for (i=0; i<N; i++) {
        if (keep[i]) {
            //printf("Keeping %dth entry (object = %s, position = %g %g %g)\n", i, kNN[i].Object, kNN[i].Position[0], kNN[i].Position[1]/(30.0*60.0*10.0), kNN[i].Position[2]/(30.0*60.0*4.0));
            new_kNN[j] = kNN[i];
            j++;
        }

    }
    return (Nout);
}


void SetupHDF5(hid_t *file, int maxNsats, char *target, char *CmdLine) {
    hid_t       space, atype, DataSet;
    herr_t      status __attribute__((unused));

    //Write global attributes
    Lgm_WriteStringAttr( *file, "Target", target );
    Lgm_WriteStringAttr( *file, "CommandLine", CmdLine );

    //Create spaces, etc. for Time, matches (4D: ID number, time, K, L*)
    atype   = CreateStrType(24);
    DataSet = CreateExtendableRank1DataSet(*file, "IsoTime", atype, &space);
    status  = H5Sclose(space);
    status  = H5Dclose(DataSet);

    //Create spaces, etc. for Time, matches (4D: ID number, time, K, L*)
    DataSet = CreateExtendableRank2DataSet(*file, "ConjTime", maxNsats, atype, &space);
    status  = H5Sclose(space);
    status  = H5Dclose(DataSet);

    //create dataset for object identifier for matches
    DataSet = CreateExtendableRank2DataSet(*file, "MatchObject", maxNsats, atype, &space);
    status  = H5Sclose(space);
    status  = H5Dclose(DataSet);
    status  = H5Tclose(atype);

    //Create dataset for query point
    DataSet = CreateExtendableRank1DataSet(*file, "QueryK", H5T_NATIVE_DOUBLE, &space);
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose(space);
    status  = H5Dclose(DataSet);

    //Create dataset for query point
    DataSet = CreateExtendableRank1DataSet(*file, "QueryLstar", H5T_NATIVE_DOUBLE, &space);
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose(space);
    status  = H5Dclose(DataSet);

    //create dataset with distance score
    DataSet = CreateExtendableRank2DataSet(*file, "Distance", maxNsats, H5T_NATIVE_DOUBLE, &space);
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose(space);
    status  = H5Dclose(DataSet);

    //create dataset with K for matches
    DataSet = CreateExtendableRank2DataSet(*file, "MatchK", maxNsats, H5T_NATIVE_DOUBLE, &space);
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose(space);
    status  = H5Dclose(DataSet);

    //create dataset with L* for matches
    DataSet = CreateExtendableRank2DataSet(*file, "MatchLstar", maxNsats, H5T_NATIVE_DOUBLE, &space);
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose(space);
    status  = H5Dclose(DataSet);
}


void WriteToHDF5(hid_t *file, int nTime, char **queryTime, double *queryK, double *queryL, int Nmatchmax, int Nmatch, char ***matchTime, Lgm_KdTreeData *kNN, double *fac) {
    hid_t       atype;
    herr_t      status __attribute__((unused));
    int         i, j;
    char        **matchObject;
    char        **outTime;
    double      *dist, *matchK, *matchL;

    LGM_ARRAY_2D( matchObject, Nmatchmax, 24, char );
    LGM_ARRAY_2D( outTime, Nmatchmax, 24, char );
    LGM_ARRAY_1D_WITH_VAL( dist, Nmatchmax, double, LGM_FILL_VALUE);
    LGM_ARRAY_1D_WITH_VAL( matchK, Nmatchmax, double, LGM_FILL_VALUE);
    LGM_ARRAY_1D_WITH_VAL( matchL, Nmatchmax, double, LGM_FILL_VALUE);

    // Write query IsoTime Strings
    atype = CreateStrType(24);
    LGM_HDF5_EXTEND_RANK1_DATASET( *file, "IsoTime", nTime, atype, queryTime[nTime] );
    LGM_HDF5_EXTEND_RANK1_DATASET( *file, "QueryK", nTime, H5T_NATIVE_DOUBLE, queryK );
    LGM_HDF5_EXTEND_RANK1_DATASET( *file, "QueryLstar", nTime, H5T_NATIVE_DOUBLE, queryL );

    // Write data for current record
    if (Nmatch>0) { //if no matchest then just write fill to the HDF5 file
        for (i=0, j=Nmatch-1; i<Nmatch; i++, j--) { //KdTreeData stores things in reverse order.
            if (kNN[j].Object) {
                strcpy(matchObject[i], (char *)kNN[j].Object);
                strcpy(outTime[i], matchTime[nTime][j]);
                dist[i] = sqrt(kNN[j].Dist2);
                matchK[i] = kNN[j].Position[1]/fac[1];
                matchL[i] = kNN[j].Position[2]/fac[2];
                }
        }
    }
    LGM_HDF5_EXTEND_RANK2_DATASET( *file, "ConjTime", nTime, Nmatchmax, atype, outTime[0] );
    LGM_HDF5_EXTEND_RANK2_DATASET( *file, "MatchObject", nTime, Nmatchmax, atype, matchObject[0] );
    status = H5Tclose(atype);
    LGM_HDF5_EXTEND_RANK2_DATASET( *file, "Distance", nTime, Nmatchmax, H5T_NATIVE_DOUBLE, dist );
    LGM_HDF5_EXTEND_RANK2_DATASET( *file, "MatchK", nTime, Nmatchmax, H5T_NATIVE_DOUBLE, matchK );
    LGM_HDF5_EXTEND_RANK2_DATASET( *file, "MatchLstar", nTime, Nmatchmax, H5T_NATIVE_DOUBLE, matchL );
}


void getObjectName(char *strin, char *strout) {
    char *pattern = "([LAN[:digit:]]{4}-.{3}|ns[[:digit:]]{2}|rbsp[-]?[ab]|th[-]?[abcde])";
    //search pattern captures all LANL geo names, nsXX, rbsp and themis
    int err, nsub, len;

    regex_t compiled;
    regcomp(&compiled, pattern, REG_EXTENDED|REG_ICASE);
    nsub = compiled.re_nsub;
    regmatch_t matchptr[nsub];
    err = regexec(&compiled, strin, nsub, matchptr, 0);
    if (err){
        printf("%s: Did not match any substrings. Satellite name unrecognized or not in filename (%s).\n", __FILE__, strin);
        exit(-1);
    }
    len = matchptr[0].rm_eo - matchptr[0].rm_so;
    strncpy(strout, strin+matchptr[0].rm_so, len);
    strout[len]='\0';
    regfree(&compiled);
}


int main( int argc, char *argv[] ) {
    struct Arguments arguments;
    double              *q, *fac, **u, **MagEphem_K, **MagEphem_Lstar, mjd, sep, **Kquery, **Lquery, jDate;
    char                **Info, **IsoTimes, *pch, *dum, InFile[512], oname[256];
    char                **queryTime, ***matchTime, *temp, **Birds, **BirdsTmp, *CmdLine, **QueryTmp;
    char                fname[128], tStr[128], Bird[128], NewStr[512];
    int                 nBirds, nBirdsTmp, iBird, setBirds, nQueryTmp, fixedTarget;
    long int            Id;
    unsigned long int   n;
    Lgm_CTrans          *c = Lgm_init_ctrans( 0 );
    Lgm_KdTree          *KdTree;
    Lgm_KdTreeData      *kNN;
    Lgm_KdTreeData      *new_kNN;
    Lgm_KdTreeData      *new_kNN2;
    Lgm_DateTime        d;
    int                 K, Kgot, i, ii, j, D=3, t, k, nQueries;
    int                 nPts = 0, ndims, status_n __attribute__((unused));
    size_t              nSats=150; //maximum number of files expected on input
    hsize_t             Dims[3];
    hid_t               file, dataset, dataspace;
    herr_t              status __attribute__((unused));
    int                 cumTimes[nSats+1], cumPts[nSats+1], nTimes[nSats+1], glob_status, verbose;
    glob_t              glob_buffer;
    long int            Nout;

    /* Set pattern for FN matching and match string for filtering given object */
    //const char          *pattern = "/mnt/projects/dream/Spacecraft/ns*/MagEphem/*/201210*h5";
    //const char          *pattern = "/mnt/projects/dream/Spacecraft/\?\?\?\?-\?\?\?/MagEphem/*/2012101[01]*h5";
    char                match[] = "default_match_string";

    /*
     *  Parse CmdLine arguments and options
     */
    strcpy( arguments.Object, match ); //set default
    strcpy( arguments.Target, "Undefined");
    strcpy( arguments.Birds, "\0");
    arguments.verbose = 0;
    arguments.Cadence = 5.0;
    /* Reconstruct the full command line as a string showing how we called this */
    for (n=0, i=0; i<argc; i++) n += strlen(argv[i]);
    n += argc;
    LGM_ARRAY_1D( CmdLine, n+1, char );
    for (i=0; i<argc; i++) {
        strcat( CmdLine, argv[i]);  // add all argv items to CmdLine string
        strcat( CmdLine, " ");      // pad with spaces
    }
    /* Now call argp_parse and make the arguments structure */
    argp_parse (&argp, argc, argv, 0, 0, &arguments);
    verbose = (arguments.verbose) ? arguments.verbose : 0;
    sep = (arguments.Distance) ? arguments.Distance : 1500.0;
    strcpy(&tStr[0], arguments.Target);

    LGM_ARRAY_2D( BirdsTmp, 300, 128, char );
    StringSplit( arguments.Birds, BirdsTmp, 128, &nBirdsTmp );
    setBirds = ( arguments.Birds[0] != '\0') ? TRUE : FALSE;

    fixedTarget = (arguments.QueryPoint) ? TRUE : FALSE;
    if (fixedTarget) {
        strcpy(&tStr[0], "FixedPoint");
        LGM_ARRAY_2D( QueryTmp, 2, 128, char );
        StringSplit( arguments.QueryPoint, QueryTmp, 128, &nQueryTmp );
    }

    //set up scaling for search/storage
    Lgm_DateTime targ, UTC, UTC_start, UTC_end;
    //LGM_ARRAY_2D( queryTime, 2, 24, char );
    LGM_ARRAY_1D( q, D, double );
    LGM_ARRAY_1D( fac, D, double );
    fac[0] = 86400.0;
    fac[1] = 30.0*60.0*8.0;//*10.0;
    fac[2] = 30.0*60.0*1.5;

    /*
     * Read in magephem data from files (just need time, K, L*)
     * UNLESS
     * target specified is a query point, then generate series
     */
    if (fixedTarget) {
        //need to set arrays: queryTime, Kquery, Lquery, nQueries
        //so need start time, end time, cadence and point (K,L)
        if (!arguments.StartDate) exit(-9);
        if (!arguments.EndDate) exit(-9);
        IsoTimeStringToDateTime( arguments.StartDate, &UTC_start, c); //JD is in UTC_start.JD
        IsoTimeStringToDateTime( arguments.EndDate, &UTC_end, c);
        //calculate number of elements between start and end dates at given cadence
        double sJD, eJD;
        sJD = Lgm_JD( UTC_start.Year, UTC_start.Month, UTC_start.Day, UTC_start.Time, LGM_TIME_SYS_UTC, c );
        eJD = Lgm_JD( UTC_end.Year, UTC_end.Month, UTC_end.Day, UTC_end.Time, LGM_TIME_SYS_UTC, c );
        nQueries = (int)((eJD - sJD)/(arguments.Cadence/1440.0));
        //nQueries++;
        //set array sizes
        LGM_ARRAY_2D( Kquery, nQueries, 2, double );
        LGM_ARRAY_2D( Lquery, nQueries, 2, double );
        LGM_ARRAY_2D( queryTime, nQueries, 32, char );
        int ind=0;
        for (jDate = sJD; jDate < eJD; jDate+=arguments.Cadence/1440.0) {
            //stuff t, K and L* into queryTime, Kquery and Lquery
            Lgm_JD_to_DateTime( jDate, &UTC, c );
            Lgm_DateTimeToString( queryTime[ind], &UTC, 0, 3 );
            Kquery[ind][0] = atof(QueryTmp[0]);
            Lquery[ind][0] = atof(QueryTmp[1]);
            ind++;
        }
    }
    else {
        //start with file for input times, positions
        strcpy( &InFile[0], arguments.args[0]);
        if (H5Fis_hdf5( InFile )) {
            // Read the InFile
            printf("Reading %s\n", InFile);
            file      = H5Fopen( InFile, H5F_ACC_RDONLY, H5P_DEFAULT );
            //Read ISO time
            queryTime = Get_StringDataset_1D( file, "IsoTime", Dims);
            if (verbose) printf("Get IsoTime done. Dims are %d\n", (int)Dims[0]);
            // Read GSM and Lstar position
            Kquery    = Get_DoubleDataset_2D( file, "K", Dims);
            if (verbose) printf("Get K done. Dims are %d, %d\n", (int)Dims[0], (int)Dims[1]);
            Lquery    = Get_DoubleDataset_2D( file, "Lstar", Dims);
            nQueries  = Dims[0];
            status    = H5Fclose( file );
            }
        else {
            printf("%s: Failed to read input HDF5 file (%s)\n", __FILE__, InFile);
            exit(-1);
            }
    }

    //if satellites names specified, search/replace on %B
    //TODO: allow %B tokens and wildcards...
    if (setBirds) {
        LGM_ARRAY_2D( Birds, nBirdsTmp, 256, char );
        if (verbose) printf("Substituting tokens in filenames\n");
        for ( iBird = 0; iBird < nBirdsTmp; iBird++ ) {
        NewStr[0] = '\0';
            Lgm_ReplaceSubString( NewStr, arguments.Pattern, "%B", BirdsTmp[iBird] );    strcpy( Birds[iBird], NewStr );
            if (verbose) printf("Adding file %s to list; %d of %d\n", Birds[iBird], iBird+1, nBirdsTmp);
        }
    }
    else {
        //now do wildcard matches for target dataset
        glob_status = glob( arguments.Pattern , 0 , NULL , &glob_buffer );
        if (glob_status != GLOB_NOMATCH) {
            nBirdsTmp = glob_buffer.gl_pathc;
            LGM_ARRAY_2D( Birds, nBirdsTmp, 256, char );
            for (iBird=0; iBird<nBirdsTmp; iBird++) {
                strcpy( Birds[iBird], glob_buffer.gl_pathv[iBird] );
                if (verbose) printf("Found file %s; %d of %d\n", Birds[iBird], iBird+1, nBirdsTmp);
            }
        }
        else {
            printf("%s: Glob couldn't match any files with the specified pattern. Status = %d.\n\n", __FILE__, glob_status);
            exit(-1);
        }
    }
    nSats = nBirdsTmp;

    // Get number of data points across all input files
    LGM_ARRAY_3D( matchTime, nQueries, (int)nSats, 24, char );
    cumTimes[0] = 0;
    cumPts[0] = 0;
    for (n=0; n<nSats; n++) {
        strcpy( &InFile[0], &Birds[n][0]);
        if (H5Fis_hdf5( InFile )) {
            // Read the InFile
            if (verbose > 2) printf("Inspecting %s\n", InFile);
            file          = H5Fopen( InFile, H5F_ACC_RDONLY, H5P_DEFAULT );
            dataset       = H5Dopen( file, "/Lstar", H5P_DEFAULT );
            dataspace     = H5Dget_space( dataset );
            status_n      = H5Sget_simple_extent_dims( dataspace, Dims, NULL );
            nTimes[n]     = Dims[0];
            cumTimes[n+1] = cumTimes[n] + Dims[0];
            cumPts[n+1]   = cumPts[n] + Dims[0]*Dims[1];
            status = H5Fclose( file );
            }
        else {
            continue;
            }
        nPts = cumPts[n+1];
        }
    if (verbose > 1) printf("nPts = %d, nSats=%lu\n", nPts, nSats);

    // alloc mem for 4-D array of state vectors
    LGM_ARRAY_2D( u, D, nPts, double );
    LGM_ARRAY_2D( Info, nPts, 180, char );

    //do this for each given satellite
    j = 0;
    for (n=0; n<nSats; n++) {
        strcpy( &InFile[0], &Birds[n][0]);
        if (H5Fis_hdf5( InFile )) {
            // Read the InFile
            printf("Reading %s\n", InFile);
            file = H5Fopen( InFile, H5F_ACC_RDONLY, H5P_DEFAULT );
            //Read ISO time
            IsoTimes = Get_StringDataset_1D( file, "/IsoTime", Dims);
            // Read GSM and Lstar position
            MagEphem_K     = Get_DoubleDataset_2D( file, "/K", Dims);
            MagEphem_Lstar = Get_DoubleDataset_2D( file, "/Lstar", Dims);
            ndims = Dims[1];
            status = H5Fclose( file );
            }
        else {
            continue;
            }
        printf("Processing file %lu (of %zu), %s\n", n+1, nSats, InFile);
        //loop over times and scale values so 30min ~ 0.1 K ~0.25 L*
        //printf("Looping over times: range %d - %d, values %d-%d\n", cumTimes[n], cumTimes[n+1], 0, nTimes[n]);
        for (t=0; t<nTimes[n]; t++) {
            IsoTimeStringToDateTime( IsoTimes[t], &d, c ); //time system is now LGM_TIME_SYS_UTC
            for (k=0; k<ndims; k++) {
                //first check whether K or Lstar are negative
                if ((MagEphem_K[t][k] < 0) || (MagEphem_Lstar[t][k] < 0)) {
                    continue; //if they are, don't add them to the Kd-tree
                    }
                u[0][j] = fac[0]*Lgm_MJD( d.Year, d.Month, d.Day, d.Time, LGM_TIME_SYS_UTC, c ); //seconds -- 30 minutes is 1800 s
                u[1][j] = fac[1]*MagEphem_K[t][k];
                u[2][j] = fac[2]*MagEphem_Lstar[t][k];
                pch = basename(InFile);
                //now apply satellite name to point
                getObjectName(pch, oname);
                strcpy( Info[j], oname);
                j++;
                }
            //printf("t = %d, j = %d, k = %d, n = %d\n", t, j, k, n);
            }
        }
    
    if (verbose) printf("Creating %d dimensional KdTree with %d points\n", D, j);
    KdTree = Lgm_KdTree_Init( u, (void **)Info, j, D );
    K = 500; //Number of Nearest Neighours to find

    /*
     * Now do nearest neighbour queries.
     *
     * To find efficiently find all NN we should loop over all target times for all(?) satellites
     * and keep only those that fall within a certain range. Store object IDs so that duplicate
     * matches can be weeded out. For each search filter out the current object.
     *
     */
    strcpy(fname, arguments.args[1]);
    file = H5Fcreate( fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
    SetupHDF5( &file, nSats, tStr, CmdLine );

    int nskip=0;
    printf("Making %d queries\n", nQueries);
    for (ii=0; ii<nQueries; ii++) {
        LGM_ARRAY_1D( kNN, K, Lgm_KdTreeData );
        LGM_ARRAY_1D( new_kNN, K, Lgm_KdTreeData );
        LGM_ARRAY_1D( new_kNN2, K, Lgm_KdTreeData );
        /*
         *  Set target ("query") state 
         */
        if (verbose>1) printf("Query %d: %s\n", ii, queryTime[ii]);
        IsoTimeStringToDateTime( queryTime[ii], &targ, c );
        //Kquery = 0.25;
        //Lquery = 4.8;
    
        mjd = Lgm_MJD( targ.Year, targ.Month, targ.Day, targ.Time, LGM_TIME_SYS_UTC, c );
        q[0] = fac[0]*mjd;
        q[1] = fac[1]*Kquery[ii][0]; // <--- Index 0 hardcodes this to lowest K in MagEphem?
        q[2] = fac[2]*Lquery[ii][0];
    
        if (verbose) printf("Finding %d nearest neighbours.\n", K );
        Lgm_KdTree_kNN( q, D, KdTree, K, &Kgot, sep*sep, kNN );
        //drop all kNN where kNN[i].Object matches value
        Nout = purgeObject( kNN, Kgot, arguments.Object, new_kNN, verbose);
        if (verbose) printf("Number of NN (excluding %s) = %ld\n", arguments.Object, Nout);
        if (Nout==0) {
            if (verbose>2) printf("%d\n", nskip++);
                WriteToHDF5( &file, ii, queryTime, &Kquery[ii][0], &Lquery[ii][0], nSats, Nout, matchTime, new_kNN, fac );
                LGM_ARRAY_1D_FREE( kNN );
                LGM_ARRAY_1D_FREE( new_kNN );
                LGM_ARRAY_1D_FREE( new_kNN2 );
                continue; //if there aren't any matches, then skip the rest of the loop
            }
        //drop all kNN for each object beyond the closest match
        int Ncopy = Nout;
        Ncopy = NearestPerObject(new_kNN, Ncopy, new_kNN2, verbose);
    
        //printf("\nQuery Point: q = (" );
        Lgm_DateTimeToString( &queryTime[ii][0], &targ, 0, 2);
        //printf("  %s", &queryTime[ii][0]);
        //for (j=1; j<D; j++) printf("  %lf", q[j]/fac[j] );
        //printf("\n");
        //printf(" )\nNearest neighbour per object\n");
    
        for (i=0; i<Ncopy; i++){ 
            Id = new_kNN2[i].Id;
            //printf("%02d: dist = %g    Id = %ld ( Object: %s ) p = (", Ncopy-i, sqrt(new_kNN2[i].Dist2), Id, (char *)new_kNN2[i].Object);
            Lgm_MJD_to_DateTime( new_kNN2[i].Position[0]/fac[0], &UTC, c  );
            Lgm_DateTimeToString( &matchTime[ii][i][0], &UTC, 0, 2);
            if (verbose > 2) {
                printf("  %s", &matchTime[ii][i][0] ); 
                for (j=1; j<D; j++) printf("  %g", new_kNN2[i].Position[j]/fac[j] );
                printf(" )\n");
            }
        }
        /*
         * Write to HDF5 file
         */
        WriteToHDF5( &file, ii, queryTime, &Kquery[ii][0], &Lquery[ii][0], nSats, Ncopy, matchTime, new_kNN2, fac );

        LGM_ARRAY_1D_FREE( kNN );
        LGM_ARRAY_1D_FREE( new_kNN );
        LGM_ARRAY_1D_FREE( new_kNN2 );
    }
    H5Fclose( file );

    LGM_ARRAY_1D_FREE( q );
    LGM_ARRAY_2D_FREE( u );
    LGM_ARRAY_2D_FREE( queryTime );
    LGM_ARRAY_2D_FREE( matchTime );
    Lgm_free_ctrans( c ); // free the structure
    if (!setBirds) globfree( &glob_buffer );
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
