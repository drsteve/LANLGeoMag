#include <math.h>
#include <quicksort.h>
#include <stdio.h>
#include <stdlib.h>
#include <Lgm_CTrans.h>
#include <Lgm_DynamicMemory.h>
#include <Lgm_ElapsedTime.h>
#include <Lgm_KdTree.h>
#include <Lgm_Sgp.h>
int main( ) {
    Lgm_ElapsedTimeInfo t;
    double              *q, **u, **B, DstIndex, dum1, dum2, dum3;
    char                **Info;
    double              Time, JDs, JDe, UTC, x, y, z, r, dist, *Dist, delta;
    long int            Id;
    unsigned long int   n, *Idx, nSearches, Id_Old;
    char                Filename[256];
    Lgm_KdTree          *KdTree;
    Lgm_KdTreeData      *kNN;
    int                 K, Kgot, i, j, jj, d, D, k, ny, nm, nd, Keep;
    Lgm_CTrans          *c = Lgm_init_ctrans( 0 );
    FILE                *fp, *fp_out;

    int         Year, Month, Day;
    int         StartYear, StartMonth, StartDay, StartDoy;
    int         EndYear, EndMonth, EndDay, EndDoy;
    long int    Date, StartDate, EndDate;
    double      UT, tsince, JD, StartUT, EndUT, StartJD, EndJD, JDinc;
    Lgm_Vector  Ugse, Uteme;
    
    


    t.ColorizeText = TRUE;
    Lgm_ElapsedTimeInit( &t, 255, 150, 0 );


    /*
     * Read in TLEs from a file
     */
    int      nTLEs = 0;
    char     TleFile[512];
    char     *OutputFile = "output.txt";
    sprintf( TleFile, "catalog.txt" );
    _SgpInfo *s    = (_SgpInfo *)calloc( 1, sizeof(_SgpInfo) );
    _SgpTLE  *TLEs = (_SgpTLE *)calloc( 20000, sizeof(_SgpTLE) );

    LgmSgp_ReadTlesFromFile( TleFile, &nTLEs, TLEs, 1 );
    printf("Found %d Valid TLEs in file %s\n", nTLEs, TleFile );

    StartDate = 20130802; StartUT = 4.0+39.0/60.0;
    EndDate   = 20130802; EndUT   = 4.0+41.0/60.0;


    Lgm_Doy( StartDate, &StartYear, &StartMonth, &StartDay, &StartDoy );
    StartJD = Lgm_JD( StartYear, StartMonth, StartDay, StartUT, LGM_TIME_SYS_UTC, c );

    Lgm_Doy( EndDate, &EndYear, &EndMonth, &EndDay, &EndDoy );
    EndJD = Lgm_JD( EndYear, EndMonth, EndDay, EndUT, LGM_TIME_SYS_UTC, c );

    //JDinc = 1.0/1440.0; // 1-min increment
    JDinc = 1.0/86400.0; // 1-sec increment

    if ( (fp = fopen( OutputFile, "w" )) == NULL ) {
        printf( "Couldnt open file %s for writing\n", OutputFile );
        exit( 1 );
    }



    // alloc mem for 4-D array of state vectors
    D = 4;
    LGM_ARRAY_2D( u, D, 14300*61*5, double );
    LGM_ARRAY_2D( Info, 14300*61*5, 80, char );

    // loop through each TLE found
    for (n=0, i=0; i<nTLEs; i++){

        // init SGP4
        LgmSgp_SGP4_Init( s, &TLEs[i] );
        fprintf(fp, "%%%s\n", TLEs[i].Line0 );
        fprintf(fp, "%%%s\n", TLEs[i].Line1 );
        fprintf(fp, "%%%s\n", TLEs[i].Line2 );
        // loop over specified time range
        for ( JD = StartJD; JD <= EndJD; JD += JDinc ) {

            // Convert the current JD back to Date/UT etc..
            Lgm_jd_to_ymdh ( JD, &Date, &Year, &Month, &Day, &UT );

            // Set up the trans matrices
            Lgm_Set_Coord_Transforms( Date, UT, c );

            // "time since" in minutes (thats what SGP4 wants)
            tsince = (JD - TLEs[i].JD)*1440.0;

            // Call SGP4. Coords are in TEME. 
            LgmSgp_SGP4( tsince, s );
            Uteme.x = s->X/WGS84_A; Uteme.y = s->Y/WGS84_A; Uteme.z = s->Z/WGS84_A;

            // Example of converting TEME->GSM coords.
            Lgm_Convert_Coords( &Uteme, &Ugse, TEME_TO_GSM, c );

            // Dump result to file
            fprintf(fp, "%8ld %.10lf %10.6lf %10.6lf %10.6lf\n", Date, UT, Ugse.x, Ugse.y, Ugse.z  );

            // Save as an array of 4D points
            u[0][n] = UT; u[1][n] = Ugse.x; u[2][n] = Ugse.y; u[3][n] = Ugse.z; 
if (i==14254 ) printf("n=%ld\n", n);
            strcpy( Info[n], TLEs[i].Line0 );
            ++n;
//printf("n = %ld\n", n);
        }
    }

    fclose(fp);



    printf("Creating %d dimensional KdTree with %ld points\n", D, n);
    Lgm_ElapsedTimeInit( &t, 255, 150, 0 );
    KdTree = Lgm_KdTree_Init( u, (void **)Info, n, D );
    Lgm_PrintElapsedTime( &t );


    /*
     *  Set target ("query") state 
     */
    LGM_ARRAY_1D( q, D, double );
long int ii;
for (ii=1724734; ii<=1724854; ii++){

    for ( d=0; d<D; d++ ) {
        q[d] = u[d][ii];
    }
    


    K = 100; 
    printf("Finding %d nearest neigbors.\n", K );
    LGM_ARRAY_1D( kNN, K, Lgm_KdTreeData );
    Lgm_KdTree_kNN( q, D, KdTree, K, &Kgot, .001*.001, kNN );


    printf("Query Point: q = (" );
    for (d=0; d<D; d++) printf("  %g", q[d] );
    printf(" )\n");

    printf("K, Kgot = %d %d\n", K, Kgot );
    for (i=0; i<Kgot; i++){
        Id = kNN[i].Id;
        printf("%02d: dist = %g    Id = %ld ( Object: %s ) p = (", i, sqrt(kNN[i].Dist2), Id, Info[Id] );
        //printf("%02d: dist = %g    Id = %ld ( Object: %s ) p = (", i, sqrt(kNN[i].Dist2), Id, "dummy" );
        printf("  %g", kNN[i].Position[0] );
        for (d=1; d<D; d++) printf("  %g", kNN[i].Position[d] );
        printf(" )\n");
    }
    printf("\n");

}










    LGM_ARRAY_1D_FREE( q );
    LGM_ARRAY_1D_FREE( kNN );
    LGM_ARRAY_2D_FREE( u );


    free( TLEs );
    free( s );
    Lgm_free_ctrans( c ); // free the structure


    return(0);

}

