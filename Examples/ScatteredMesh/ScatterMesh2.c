#include <stdio.h>
#include <stdlib.h>
#include <Lgm_CTrans.h>
#include <Lgm_QinDenton.h>
#include <Lgm_MagModelInfo.h>
#include <Lgm_DynamicMemory.h>
#include <Lgm_ElapsedTime.h>
#include <Lgm_Octree.h>
void quicksort_uli( unsigned long n, unsigned long *arr );
int main( ) {
    Lgm_ElapsedTimeInfo t;
    Lgm_CTrans          *c = Lgm_init_ctrans( 1 ); // more compact declaration
    Lgm_QinDentonOne    p;
    Lgm_MagModelInfo    *mInfo;
    Lgm_Vector          *u, *B, *B1, q, Binterp, Bcd, Bmodel, Diff, v, w, v1 , v2, v3;
    double              Bx, By, Bz, B1x, B1y, B1z;
    double              Time, JD, x, y, z, r, d, Lat;
    long int            Date, n;
    int                 i, Gap, FLAG;
    FILE                *fp, *fp2;

    Lgm_Octree          *Octree;
    

    t.ColorizeText = TRUE;
    Lgm_ElapsedTimeInit( &t, 255, 150, 0 );



    LGM_ARRAY_1D( u, 50000000, Lgm_Vector );
    LGM_ARRAY_1D( B, 50000000, Lgm_Vector );
    LGM_ARRAY_1D( B1, 50000000, Lgm_Vector );


    mInfo = Lgm_InitMagInfo( );


    Date = 20020313;                        // March 3, 2002.
    Time = 12.001;                          // Gives approx zero tilt angle.
    Lgm_Set_Coord_Transforms( Date, Time, mInfo->c );

    JD = Lgm_Date_to_JD( Date, Time, c );    // Compute JD

    // Get (interpolate) the QinDenton vals from the values in the file at the given Julian Date
    Lgm_get_QinDenton_at_JD( JD, &p, 0 );


    Lgm_Set_Coord_Transforms( Date, Time, mInfo->c );
    Lgm_set_QinDenton( &p, mInfo );



    n = 0;
    Lgm_ElapsedTimeInit( &t, 255, 150, 0 );
    fp = fopen( "bats_r_us.txt", "r" );
    while ( fscanf( fp, "%lf %lf %lf %lf %lf %lf %lf %lf %lf", &x, &y, &z, &Bx, &By, &Bz, &B1x, &B1y, &B1z ) != EOF ) {
        if ( ( B1x != 0.0 ) && ( B1y != 0.0 ) && ( B1y != 0.0 ) ) {
//            if ( x*x+y*y+z*z < 20.0*20.0 ) {
                u[n].x  = x; u[n].y  = y; u[n].z  = z;
                B[n].x  = Bx; B[n].y  = By; B[n].z  = Bz;
                B1[n].x = B1x; B1[n].y = B1y; B1[n].z = B1z;
                ++n;
//            }
        }
    }
    fclose(fp);
    printf("n = %ld\n", n);
    Lgm_PrintElapsedTime( &t );






    /*
     * Create Octree
     */
    printf("Creating Octree\n");
    Lgm_ElapsedTimeInit( &t, 255, 150, 0 );
    Octree = Lgm_InitOctree( u, B1, n );
    //Octree = Lgm_InitOctree( u, B, n );
    printf("Min, Max, Diff = %g %g %g\n", Octree->Min, Octree->Max, Octree->Diff);
    Lgm_PrintElapsedTime( &t );

mInfo->Hmax = 0.1;
printf("mInfo->Hmax = %g\n", mInfo->Hmax);


    Lgm_ElapsedTimeInit( &t, 255, 150, 0 );

    /*
     * Test a FL trace
     */
    mInfo->Octree = Octree;
    Lgm_Set_Octree_kNN_k( mInfo, 4 );

    //Lgm_MagModelInfo_Set_MagModel( LGM_CDIP, LGM_EXTMODEL_TS04, mInfo );
    //Lgm_MagModelInfo_Set_MagModel( LGM_CDIP, LGM_EXTMODEL_NULL, mInfo );
    Lgm_MagModelInfo_Set_MagModel( LGM_CDIP, LGM_EXTMODEL_SCATTERED_DATA2, mInfo );
//mInfo->Lgm_MagStep_Integrator = LGM_MAGSTEP_ODE_BS;
mInfo->Lgm_MagStep_RK5_Eps    = 1e-1;
mInfo->Lgm_MagStep_BS_Eps     = 1e-4;
    mInfo->Lgm_TraceLine_Tol  = 1e-6;


    //mInfo->Bfield = Lgm_B_T89;
    fp = fopen("line_xz.txt", "w");
    fp2 = fopen("line_xz2.txt", "w");
    for (Lat = 30.0; Lat<=90.0; Lat += 1.0){
        v.x = -2.0*cos( Lat*RadPerDeg ); v.y = 0.0; v.z = 2.0*sin( Lat*RadPerDeg );
        printf("v = %g %g %g\n", v.x, v.y, v.z);
        Lgm_TraceLine( &v, &w, 120.0, -1.0, 1e-7, FALSE, mInfo );
        //FLAG = Lgm_Trace( &v, &v1, &v2, &v3, 120.0, 1e-2, 1e-4, mInfo );
        //printf("%d   v1, v2, v3  = %g %g %g   %g %g %g    %g %g %g\n", FLAG, v1.x, v1.y, v1.z, v2.x, v2.y, v2.z, v3.x, v3.y, v3.z );
        for (Gap = 3, i=0; i<mInfo->nPnts; i++){
            fprintf( fp, "%g %g %d\n", mInfo->Px[i], mInfo->Pz[i], Gap );
            fprintf( fp2, "%g %g\n", mInfo->Px[i], mInfo->Pz[i] );
            Gap = 2;
        }
    }

    for (Lat = 30.0; Lat<=90.0;  Lat += 1.0){
        v.x = 1.0*cos( Lat*RadPerDeg ); v.y = 0.0; v.z = 1.0*sin( Lat*RadPerDeg );
        printf("v = %g %g %g\n", v.x, v.y, v.z);
        Lgm_TraceLine( &v, &w, 120.0, -1.0, 1e-7, FALSE, mInfo );
        for (Gap = 3, i=0; i<mInfo->nPnts; i++){
            fprintf( fp, "%g %g %d\n", mInfo->Px[i], mInfo->Pz[i], Gap );
            fprintf( fp2, "%g %g\n", mInfo->Px[i], mInfo->Pz[i] );
            Gap = 2;
        }
    }

    for (Lat = -30.0; Lat>=-90.0;  Lat -= 1.0){
        v.x = 1.0*cos( Lat*RadPerDeg ); v.y = 0.0; v.z = 1.0*sin( Lat*RadPerDeg );
        printf("v = %g %g %g\n", v.x, v.y, v.z);
        Lgm_TraceLine( &v, &w, 120.0, 1.0, 1e-7, FALSE, mInfo );
        for (Gap = 3, i=0; i<mInfo->nPnts; i++){
            fprintf( fp, "%g %g %d\n", mInfo->Px[i], mInfo->Pz[i], Gap );
            fprintf( fp2, "%g %g\n", mInfo->Px[i], mInfo->Pz[i] );
            Gap = 2;
        }
    }

    for (Lat = -30.0; Lat>=-90.0;  Lat -= 1.0){
        v.x = -1.0*cos( Lat*RadPerDeg ); v.y = 0.0; v.z = 1.0*sin( Lat*RadPerDeg );
        printf("v = %g %g %g\n", v.x, v.y, v.z);
        Lgm_TraceLine( &v, &w, 120.0, 1.0, 1e-7, FALSE, mInfo );
        for (Gap = 3, i=0; i<mInfo->nPnts; i++){
            fprintf( fp, "%g %g %d\n", mInfo->Px[i], mInfo->Pz[i], Gap );
            fprintf( fp2, "%g %g\n", mInfo->Px[i], mInfo->Pz[i] );
            Gap = 2;
        }
    }

    fclose(fp);
    fclose(fp2);

    Lgm_PrintElapsedTime( &t );

    printf("Number of kNN lookups: %ld\n",  Octree->kNN_Lookups );
    printf("Number of HASH_FIND()s: %ld\n", mInfo->RBF_nHashFinds );
    printf("Number of HASH_ADD_KEYPTR()s: %ld\n", mInfo->RBF_nHashAdds );

//sleep(3000);





    printf("Feeing Octree\n");
    Lgm_FreeOctree( Octree );
    printf("Octree Freed\n");




    Lgm_free_ctrans( c );
    Lgm_FreeMagInfo( mInfo );
    LGM_ARRAY_1D_FREE( u );
    LGM_ARRAY_1D_FREE( B );
    LGM_ARRAY_1D_FREE( B1 );

    return(0);

}

