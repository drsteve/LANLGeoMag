#include <stdlib.h>
#include <Lgm_CTrans.h>
#include <Lgm_QinDenton.h>
#include <Lgm_MagModelInfo.h>
#include <Lgm_DynamicMemory.h>
#include <Lgm_ElapsedTime.h>
#include <Lgm_Octree.h>
int main( ) {
    Lgm_ElapsedTimeInfo t;
    Lgm_CTrans          *c = Lgm_init_ctrans( 1 ); // more compact declaration
    Lgm_QinDentonOne    p;
    Lgm_MagModelInfo    *mInfo;
    Lgm_Vector          *u, *B, q, v;
    double              Time, JD, x, y, z, r, d, dist;
    long int            Date, n;

    Lgm_Octree          *Octree;
    Lgm_OctreeData      *kNN;
    int                 K, Kgot, i, j;
    

    t.ColorizeText = TRUE;
    Lgm_ElapsedTimeInit( &t, 255, 150, 0 );



    LGM_ARRAY_1D( u, 5000000, Lgm_Vector );
    LGM_ARRAY_1D( B, 5000000, Lgm_Vector );


    mInfo = Lgm_InitMagInfo( );

    Date = 20020713;                        // August 12, 2004
    Time = 18.0 + 0.0/60.0 + 30.0/3600.0;   // Universal Time Coordinated (in decimal hours)
    JD = Lgm_Date_to_JD( Date, Time, c );    // Compute JD

    // Get (interpolate) the QinDenton vals from the values in the file at the given Julian Date
    Lgm_get_QinDenton_at_JD( JD, &p, 0 );


    Lgm_Set_Coord_Transforms( Date, Time, mInfo->c );
    Lgm_set_QinDenton( &p, mInfo );
    mInfo->Bfield = Lgm_B_T89;


    d = 0.3;
    n = 0;
    Lgm_PrintCurrentTime( &t );
    for ( x = -15.0; x <= 15.0; x += d ){
        for ( y = -15.0; y <= 15.0; y += d ){
            for ( z = -15.0; z <= 15.0; z += d ){
                u[n].x = x; u[n].y = y; u[n].z = z;
                r = Lgm_Magnitude( &u[n] );
                //if (r > 1.1){
                //    mInfo->Bfield( &u[n], &B[n], mInfo );
                    ++n;
                //}
            }
        }
    }
    Lgm_PrintElapsedTime( &t );
    printf("n = %ld\n", n);




    /*
     * Test kNN algorithm.
     */
    printf("Creating Octree\n");
    Lgm_ElapsedTimeInit( &t, 255, 150, 0 );
    Octree = Lgm_InitOctree( u, B, n );
    printf("Min, Max, Diff = %g %g %g\n", Octree->Min, Octree->Max, Octree->Diff);
    Lgm_PrintElapsedTime( &t );

    q.x = 3.2233;
    q.y = 2.4698;
    q.z = 1.35193;

    K = 4; 
    LGM_ARRAY_1D( kNN, K, Lgm_OctreeData );



    printf("\n\nTesting kNN (1000000 times)\n");
    Lgm_ElapsedTimeInit( &t, 255, 150, 0 );
    for (j=0; j<=1000000; j++){
        q.x = 30.0*rand()/(double)RAND_MAX - 15.0;
        q.y = 30.0*rand()/(double)RAND_MAX - 15.0;
        q.z = 30.0*rand()/(double)RAND_MAX - 15.0;
        Lgm_Octree_kNN( &q, Octree, K, &Kgot, 1.0*1.0, kNN );
    }
    Lgm_PrintElapsedTime( &t );



    printf("\n\nKgot = %d   q = %g %g %g\n", Kgot, q.x, q.y, q.z);
    for (i=0; i<Kgot; i++){
        Lgm_OctreeUnScalePosition( &(kNN[i].Position), &v, Octree );
        Lgm_OctreeUnScaleDistance( sqrt(kNN[i].Dist2), &dist, Octree );
        printf("%02d: dist = %g   v = %g %g %g\n", i, dist, v.x, v.y, v.z );
    }





    LGM_ARRAY_1D_FREE( kNN );

    printf("Freeing Octree\n");
    Lgm_ElapsedTimeInit( &t, 255, 150, 0 );
    Lgm_FreeOctree( Octree );
    printf("Octree Freed\n");
    Lgm_PrintElapsedTime( &t );




    Lgm_free_ctrans( c );
    Lgm_FreeMagInfo( mInfo );
    LGM_ARRAY_1D_FREE( u );
    LGM_ARRAY_1D_FREE( B );

    return(0);

}

