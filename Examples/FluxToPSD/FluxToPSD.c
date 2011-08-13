#include <stdio.h>
#include <math.h>
#include <Lgm_MagModelInfo.h>
#include <Lgm_FluxToPsd.h>
#include <hdf5.h>

int main( ) {

    Lgm_CTrans     *c = Lgm_init_ctrans( 0 );
    Lgm_DateTime    d;
    Lgm_Vector      u, v;
    Lgm_FluxToPsd  *f2p;
    Lgm_PsdToFlux  *p2f;

    double          ***J, *E, *A, **J_DIFF;
    int             nE, nA, i, j, n, k;

    double          *Mu, *K;
    int             nMu, nK;

    double          f, p2c2, df, sa;

    double          S[100], a, b, r;

    char            Line[5000], ISO_TIME[200][80], Command[512];

    long int        pos;
    FILE            *fp;


    /*
     *  Allocate some memory.
     */
    nE = 10;
    nA = 18;
    LGM_ARRAY_1D( E, nE, double );
    LGM_ARRAY_1D( A, nA, double );
    LGM_ARRAY_3D( J, 200, nE, nA, double );


    /*
     *  Choose a set of Mu's and K's
     *  and compute Phase Space Density at these constant Mu's and K's
     */
    nMu = 36;
    nK  = 36;
    LGM_ARRAY_1D( Mu, nMu, double );
    LGM_ARRAY_1D( K, nK, double );
    for (i=0; i<nMu; i++) Mu[i] = 1.0 + 60.0*i;


    // Geometric progression for K's
    a = 0.01; b = 10.0;
    r = pow( b/a, 1.0/((double)(nK-1)));
    S[0] = a;
    for (j=1; j<nK; j++) S[j] = S[j-1]*r;
    //for (j=0; j<nK; j++) K[j] = S[nK-1-j];
    for (j=0; j<nK; j++) K[j] = S[j];


    // Geometric progression for Mu's
    a = 1.0; b = 2000.0;
    r = pow( b/a, 1.0/((double)(nMu-1)));
    S[0] = a;
    for (j=1; j<nMu; j++) S[j] = S[j-1]*r;
    for (j=0; j<nMu; j++) Mu[j] = S[j];


    /*
     * Set date/time and position
     *  NEED TO READ THIS IN.
     */
    //d = Lgm_DateTime_Create( 2001, 1, 1, 0.0, LGM_TIME_SYS_UTC, c );
    u.x = 6.6*cos(133.0*RadPerDeg);
    u.y = 6.6*sin(133.0*RadPerDeg);
    u.z = 0.0;



    /*
     * Read in a sample GEO file
     */
    n = 0;
    fp = fopen("/data2/DREAM/Satellites/LANL-97A/RawFlux/2002/20021211_LANL-97A_l3_data.txt", "r");
    pos = ftell( fp );
    n = 0;
    while ( fgets( Line, 4096, fp ) != NULL ) {
        if ( Line[0] == '\n' ) {
            // do nothing
        } else if ( Line[0] != '#' ) {
            // backup so the line is available to read again.
            fseek( fp, pos, SEEK_SET );

            // read in data
            for (j=0; j<18; j++) {
                fscanf( fp, "%s", ISO_TIME[n] );
                if (j==0) printf("Read data for Time: %s\n", ISO_TIME[n] );
                fscanf( fp, "%lf", &A[j] );
                for (i=0; i<10; i++) {
                    fscanf( fp, "%lf", &J[n][i][j] );
                }
            }
            ++n;
        } else if ( strstr( Line, "Energy Channels:" ) != NULL ) {
            // extract energy channel info from header.
            sscanf( Line, "%*[^:]:%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &E[0], &E[1], &E[2], &E[3], &E[4], &E[5], &E[6], &E[7], &E[8], &E[9] );
        }
        pos = ftell( fp );
    }
    fclose( fp );



    /*
     *   Create HDF5 file.
     */
    hid_t       file;
    hid_t       dataspace, dataset, dataspace_Time, dataset_Time;
    hid_t       filespace;
    hid_t       cparms;
    hsize_t     dims[3]    = { 1, 36, 36 };
    hsize_t     maxdims[3] = {H5S_UNLIMITED, 36, 36};
    hsize_t     offset[3];
    hsize_t     chunk_dims[3] = {1, 36, 36};
    herr_t      status;


    // Create a new file. If file exists its contents will be overwritten.
    file = H5Fcreate( "myfile.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );

    // Create the data space with unlimited dimensions.
    dataspace = H5Screate_simple( 3, dims, maxdims ); // rank 3

    // Modify dataset creation properties, i.e. enable chunking
    cparms = H5Pcreate( H5P_DATASET_CREATE );
    status = H5Pset_chunk(  cparms, 3, chunk_dims );

    // Create a new dataset within the file using cparms creation properties.
    dataset = H5Dcreate( file, "PhaseSpaceDensity", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, cparms, H5P_DEFAULT );






    /*
     *   Allocate  space for J_DIFF
     */
    LGM_ARRAY_2D( J_DIFF, nE, nA, double );
    //for (k=0; k<n; k++){
    for (k=0; k<1; k++){

            // Flux -> PSD
            f2p = Lgm_F2P_CreateFluxToPsd(1);
            f2p->Extrapolate = 3;

            IsoTimeStringToDateTime( ISO_TIME[k], &d, c );
            printf("Date = %ld   UTC = %g\n", d.Date, d.Time );
            Lgm_Set_Coord_Transforms( d.Date, d.Time, c );

            Lgm_Convert_Coords( &u, &v, WGS84_TO_GSM, c );
            Lgm_F2P_SetDateTimeAndPos( &d, &v, f2p );

            // Add diff. flux data/info to f2p structure.
            Lgm_F2P_SetFlux( J[k], E, nE, A, nA, f2p );

            //  Compute Phase Space Density at the constant Mu's and K's
            Lgm_F2P_GetPsdAtConstMusAndKs( Mu, nMu, K, nK, f2p );

            int iii, jjj;
            for ( iii=0; iii<nMu; iii++ ){
                printf("[");
                for (jjj=0; jjj<nK; jjj++ ) {
                    printf("%8.1g ", f2p->PSD_MK[iii][jjj]);
                }
                printf("]\n");
            }


            // Select a hyperslab
            filespace = H5Dget_space( dataset );
            offset[0] = k; offset[1] = 0; offset[2] = 0;
            status = H5Sselect_hyperslab( filespace, H5S_SELECT_SET, offset, NULL, chunk_dims, NULL );

            // Define memory space.
            dataspace = H5Screate_simple( 3, chunk_dims, NULL );


            // Write the data to the hyperslab
            status = H5Dwrite( dataset, H5T_NATIVE_DOUBLE, dataspace, filespace, H5P_DEFAULT, &f2p->PSD_MK[0][0] );

            if ( k != n-1 ) {
                // extend the dataset (add a new time slot).
                ++dims[0];
                status = H5Dextend( dataset, dims );
            }

            Lgm_F2P_FreeFluxToPsd( f2p );

    }


    status = H5Sclose( filespace );
    status = H5Sclose( dataspace );
    status = H5Dclose( dataset );
    status = H5Fclose( file );


    LGM_ARRAY_1D_FREE( E );
    LGM_ARRAY_1D_FREE( A );
    LGM_ARRAY_3D_FREE( J );
    LGM_ARRAY_2D_FREE( J_DIFF );
    LGM_ARRAY_1D_FREE( Mu );
    LGM_ARRAY_1D_FREE( K );


    Lgm_free_ctrans( c );


    return(0);
}

