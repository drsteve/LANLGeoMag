#include <stdlib.h>
#include <math.h>
#include <Lgm_DynamicMemory.h>
#include <Lgm_ElapsedTime.h>
#include <Lgm_CTrans.h>
#include <Lgm_KdTree.h>
#include <quicksort.h>
#include <stdio.h>
int main( ) {
    Lgm_ElapsedTimeInfo t;
    double              *q, **u, **B, *Dst, DstIndex, dum1, dum2, dum3;
    double              Time, JD, JDs, JDe, UTC, x, y, z, r, dist, *Dist, delta;
    long int            Date, KeepList[10000], nKeepList, Id;
    unsigned long int   n, *Idx, nSearches, Id_Old;
    char                Filename[256];
    Lgm_KdTree          *KdTree;
    Lgm_KdTreeData      *kNN;
    int                 K, Kgot, i, j, jj, d, D, k, ny, nm, nd, Keep;
    Lgm_CTrans          *c = Lgm_init_ctrans( 0 );
    FILE                *fp, *fp_out;
    


    t.ColorizeText = TRUE;
    Lgm_ElapsedTimeInit( &t, 255, 150, 0 );

    Dst = (double *)calloc( 5000000, sizeof(double));


    /*
     * Read in Dst data from 1976 -> 2011
     */
    Date = 19580101; UTC = 12.0; JDs = Lgm_Date_to_JD( Date, UTC, c );
    Date = 19760101; UTC = 12.0; JDs = Lgm_Date_to_JD( Date, UTC, c );
    Date = 20110101; UTC = 12.0; JDe = Lgm_Date_to_JD( Date, UTC, c );
    Date = 19760102; UTC = 12.0; JDe = Lgm_Date_to_JD( Date, UTC, c );
    n = 0;
    for ( JD=JDs; JD<=JDe; JD += 1.0 ) {
        Date = Lgm_JD_to_Date( JD, &ny, &nm, &nd, &UTC );
        sprintf( Filename, "/home/mghenderson/Data/Dst/%4d/Dst_%ld.dat", ny, Date );
        printf("reading: %s\n", Filename );
        //if ( (fp = fopen( Filename, "r" ) ) != NULL ) {
        //    while ( fscanf( fp, "%lf %lf %lf %lf %lf\n", &Time, &dum1, &dum2, &dum3, &DstIndex ) != EOF ){
        //        Dst[n] = DstIndex;
        //        if ((Time > 0.0)&&(Time < 1.0)&&(Date == 20011124)) printf("n = %ld\n", n);
        //        if ( n == 120228 ) printf("Date = %ld\n", Date);
        //        if ( n == 250434 ) printf("Date = %ld\n", Date);
        //        ++n;
        //    }
        //    fclose( fp );

        //}
    }
    printf("n = %ld\n", n);


    /*
     *  Now construct state vectors.
     */
    D = 96*2;
    LGM_ARRAY_2D( u, D, n, double );
    for ( j=0; j<n-2*D; j++ ) {
        for ( d=0; d<D; d++ ) {
            u[d][j] = Dst[j+d];
        }
    }
    printf("Creating %d dimensional KdTree with %ld points\n", D, n);
    Lgm_ElapsedTimeInit( &t, 255, 150, 0 );
    KdTree = Lgm_KdTree_Init( u, (void **)B, n, D );
    Lgm_PrintElapsedTime( &t );


    /*
     *  Set target ("query") state 
     */
    LGM_ARRAY_1D( q, D, double );
    for ( d=0; d<D; d++ ) {
        q[d] = u[d][227016];
    }
    


    K = 500; 
    printf("Finding %d nearest neigbors.\n", K );
    LGM_ARRAY_1D( kNN, K, Lgm_KdTreeData );
    Lgm_KdTree_kNN( q, D, KdTree, K, &Kgot, 1000.0*1000.0, kNN );


    printf("Query Point: q = (" );
    for (d=0; d<D; d++) printf("  %g", q[d] );
    printf(" )\n");

    printf("K, Kgot = %d %d\n", K, Kgot );
    for (i=0; i<Kgot; i++){
        printf("%02d: dist = %g    Id = %ld p = (", i, sqrt(kNN[i].Dist2), kNN[i].Id );
        for (d=0; d<D; d++) printf("  %g", kNN[i].Position[d] );
        printf(" )\n");
    }
    printf("\n");











    KeepList[0] = Kgot-1; // keep the nearest one
    nKeepList   = 1;
    for ( i=Kgot-2; i>=0; i-- ){
        Id = kNN[i].Id;

        // loop over the events we have in the keeplist to see if this is a differentm event
        Keep = 1;
        for ( j=0; j<nKeepList; j++ ) {
            //printf("A, B, Diff = %ld %ld %g\n", Id, kNN[ KeepList[j] ].Id, fabs( (double)Id - (double)kNN[ KeepList[j] ].Id ));
            if ( fabs( (double)Id - (double)kNN[ KeepList[j] ].Id ) < 48 ) {
                //printf("i=%d too short: Id, kNN[ KeepList[j] ].Id = %ld %ld\n", i, Id, kNN[ KeepList[j] ].Id);
                Keep = 0;
                break;
            }
        }

        if ( Keep ) {
            KeepList[ nKeepList++ ] = i;
        }

    }
            

    printf("Query Point: q = (" );
    for (d=0; d<D; d++) printf("  %g", q[d] );
    printf(" )\n");

    fp_out = fopen("DstNNs.dat", "w");
if (nKeepList > 4) nKeepList = 4;
    for ( jj=0; jj<nKeepList; jj++ ) {
        j = KeepList[jj];
        fprintf( fp_out, "%d %g 3\n", 0, kNN[j].Position[0] );
        for (d=1; d<D; d++)  fprintf( fp_out, "%d %g 2\n", d, kNN[j].Position[d] );

        printf("%02d: dist = %g    Id = %ld p = (", j, sqrt(kNN[j].Dist2), kNN[j].Id );
        for (d=0; d<D; d++) printf("  %g", kNN[j].Position[d] );
        printf(" )\n");

    }
    fclose(fp_out);

        








    LGM_ARRAY_1D_FREE( q );
    LGM_ARRAY_1D_FREE( kNN );
    LGM_ARRAY_2D_FREE( u );


    Lgm_free_ctrans( c ); // free the structure

    free( Dst );

    return(0);

}

