#include <math.h>
#include <quicksort.h>
#include <stdlib.h>
#include <Lgm_DynamicMemory.h>
#include <Lgm_ElapsedTime.h>
#include <Lgm_KdTree.h>
int main( ) {
    Lgm_ElapsedTimeInfo t;
    double              **q, **u, **B;
    double              Time, JD, x, y, z, r, dist, *Dist, delta;
    long int            Date;
    unsigned long int   n, *Idx, nSearches;

    Lgm_KdTree          *KdTree;
    Lgm_KdTreeData      *kNN;
    int                 K, Kgot, i, j, d, D, k;
    

    t.ColorizeText = TRUE;
    Lgm_ElapsedTimeInit( &t, 255, 150, 0 );


    /*
     * Test kNN algorithm.
     */
    n = 24*365*50;
    D = 24;
    LGM_ARRAY_2D( u, D, n, double );
    LGM_ARRAY_2D( B, D, n, double );
    // generate random data points
    for ( j=0; j<n; j++ ) {
        for ( d=0; d<D; d++ ) u[d][j] = 1.0*rand()/(double)RAND_MAX - 0.0;
    }
    printf("Creating %d dimensional KdTree with %ld points\n", D, n);
    Lgm_ElapsedTimeInit( &t, 255, 150, 0 );
    KdTree = Lgm_KdTree_Init( u, (void **)B, n, D );
    Lgm_PrintElapsedTime( &t );



    printf("Testing kNN (k-Nearest-Neighbor) search.\n");
    nSearches = 1;
    nSearches = 100;
    Lgm_ElapsedTimeInit( &t, 255, 150, 0 );
    // generate random query point
    LGM_ARRAY_2D( q, nSearches, D, double );
    K = 50; 
    LGM_ARRAY_1D( kNN, K, Lgm_KdTreeData );
clock_t StartTime = clock();
    for ( k=0; k<nSearches; k++ ) {
        for ( d=0; d<D; d++ ) q[k][d] = 1.0*rand()/(double)RAND_MAX - 0.0;

/*
        printf("Query Point: q = (" );
        for (d=0; d<D; d++) printf("  %g", q[k][d] );
        printf(" )\n");
*/

        
        Lgm_KdTree_kNN( q[k], D, KdTree, K, &Kgot, 2.0*2.0, kNN );

/*
        printf("K, Kgot = %d %d\n", K, Kgot );
        for (i=0; i<Kgot; i++){
            printf("%02d: dist = %g    p = (", i, sqrt(kNN[i].Dist2) );
            for (d=0; d<D; d++) printf("  %g", kNN[i].Position[d] );
            printf(" )\n");
        }
        printf("\n");
*/

    }
clock_t EndTime = clock();
printf("Sec = %g\n", (EndTime-StartTime)/(double)CLOCKS_PER_SEC);
printf("Searches/Sec = %g\n", nSearches/((EndTime-StartTime)/(double)CLOCKS_PER_SEC));
    Lgm_PrintElapsedTime( &t );


    /*
     * Brute Force
     */
if (0==1){
    Lgm_ElapsedTimeInit( &t, 255, 150, 0 );
    LGM_ARRAY_1D( Dist, n, double );
    LGM_ARRAY_1D( Idx,  n, unsigned long int );
StartTime = clock();
    for ( k=0; k<nSearches; k++ ) {
        for ( j=0; j<n; j++ ) {
            Dist[j] = 0.0;
            for ( d=0; d<D; d++ ) {
                delta = u[d][j] - q[k][d];
                Dist[j] += delta*delta;
            }
            Idx[j] = j;
        }
        quicksort2uli( n, Dist-1, Idx-1 );
        printf("Query Point: q = (" );
        for (d=0; d<D; d++) printf("  %g", q[k][d] );
        printf(" )\n");
        for (i=0; i<K; i++){
            printf("%02d: dist = %g    p = (", i, sqrt(Dist[i]) );
            for (d=0; d<D; d++) printf("  %g", u[d][ Idx[i] ] );
            printf(" )\n");
        }
        printf("\n");

    }
EndTime = clock();
printf("Time = %g\n", (EndTime-StartTime)/(double)CLOCKS_PER_SEC);
}
    LGM_ARRAY_1D_FREE( Dist );
    LGM_ARRAY_1D_FREE( Idx );
    Lgm_PrintElapsedTime( &t );






/*
    printf("Freeing KdTree\n");
    Lgm_ElapsedTimeInit( &t, 255, 150, 0 );
    Lgm_FreeKdTree( KdTree );
    printf("KdTree Freed\n");
    Lgm_PrintElapsedTime( &t );

*/



    LGM_ARRAY_1D_FREE( q );
    LGM_ARRAY_1D_FREE( kNN );
    LGM_ARRAY_2D_FREE( u );
    LGM_ARRAY_2D_FREE( B );

    return(0);

}

