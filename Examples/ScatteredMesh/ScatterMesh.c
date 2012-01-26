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
    Lgm_Vector          *u, *B, q, Binterp, Bcd, Bmodel, Diff;
    double              Time, JD, x, y, z, r, d;
    long int            Date, n;

    Lgm_Octree          *Octree;
    

    t.ColorizeText = TRUE;
    Lgm_ElapsedTimeInit( &t, 255, 150, 0 );



    LGM_ARRAY_1D( u, 50000000, Lgm_Vector );
    LGM_ARRAY_1D( B, 50000000, Lgm_Vector );


    mInfo = Lgm_InitMagInfo( );

    Date = 20020713;                        // August 12, 2004
    Time = 18.0 + 0.0/60.0 + 30.0/3600.0;   // Universal Time Coordinated (in decimal hours)
    JD = Lgm_Date_to_JD( Date, Time, c );    // Compute JD

    // Get (interpolate) the QinDenton vals from the values in the file at the given Julian Date
    Lgm_get_QinDenton_at_JD( JD, &p, 0 );


    Lgm_Set_Coord_Transforms( Date, Time, mInfo->c );
    Lgm_set_QinDenton( &p, mInfo );


    d = 0.25;
    n = 0;
    Lgm_PrintCurrentTime( &t );
    for ( x = -15.0; x <= 15.0; x += d ){
        for ( y = -15.0; y <= 15.0; y += d ){
            for ( z = -15.0; z <= 15.0; z += d ){
                q.x = x; q.y = y; q.z = z;
                r = Lgm_Magnitude( &q );
                if (r > 1.1){
                    u[n].x = q.x; u[n].y = q.y; u[n].z = q.z;

                    mInfo->Bfield = Lgm_B_T89;
                    mInfo->Bfield( &u[n], &Bmodel, mInfo );

                    mInfo->Bfield = Lgm_B_edip;
                    mInfo->Bfield( &u[n], &Bcd, mInfo );

                    B[n].x = Bmodel.x - Bcd.x;
                    B[n].y = Bmodel.y - Bcd.y;
                    B[n].z = Bmodel.z - Bcd.z;

                    ++n;
                }
            }
        }
    }
    Lgm_PrintElapsedTime( &t );
    printf("n = %ld\n", n);




    /*
     * Create Octree
     */
    printf("Creating Octree\n");
    Lgm_ElapsedTimeInit( &t, 255, 150, 0 );
    Octree = Lgm_InitOctree( u, B, n );
    printf("Min, Max, Diff = %g %g %g\n", Octree->Min, Octree->Max, Octree->Diff);
    Lgm_PrintElapsedTime( &t );




    /*
     * Test interpolation...
     */
    mInfo->Octree = Octree;
    Lgm_Set_Octree_kNN_k( mInfo, 12 );
    d = 0.01; y = 0.0; z = 0.0;
    y = 0.0; z = 0.0;
    for ( x = -9.0; x <= -2.0; x += d ){
        q.x = x; q.y = y; q.z = z;

        mInfo->Bfield = Lgm_B_FromScatteredData;
        mInfo->Bfield( &q, &Binterp, mInfo );

        mInfo->Bfield = Lgm_B_edip;
        mInfo->Bfield( &q, &Bcd, mInfo );

        Binterp.x += Bcd.x;
        Binterp.y += Bcd.y;
        Binterp.z += Bcd.z;

        mInfo->Bfield = Lgm_B_T89;
        mInfo->Bfield( &q, &Bmodel, mInfo );

        Lgm_VecSub( &Diff, &Bmodel, &Binterp );

        printf("q = %g %g %g   Binterp = %g %g %g  Bmodel = %g %g %g   Diff = %g %g %g   |Diff| = %g  Pcnt Error = %g\n", q.x, q.y, q.z, Binterp.x, Binterp.y, Binterp.z, Bmodel.x, Bmodel.y, Bmodel.z, Diff.x, Diff.y, Diff.z, Lgm_Magnitude(&Diff), Lgm_Magnitude(&Diff)/Lgm_Magnitude(&Bmodel)*100.0);
    }








    printf("Feeing Octree\n");
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

