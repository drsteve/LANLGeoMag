#include <Lgm_MagModelInfo.h>
#include <Lgm_DynamicMemory.h>

/* BAL 04-Jan-2011 modified for more cases */

int main(){


    long int            Date;
    double              UTC;
    double              GradBvec[4][4];
    Lgm_Vector          u, B, CurlB, GradB;
    Lgm_MagModelInfo    **mInfo;
    int                 i, j, n;

    j = 167;
    n = 8640;
    LGM_ARRAY_1D( mInfo, n, Lgm_MagModelInfo * );

    /*
     * Initialize an array of mInfo[i] structures...
     */
    Date = 20050831;
    for ( i=0; i<8640; ++i ) {
        UTC  = (double)i*10.0/3600.0;
        mInfo[i] = Lgm_InitMagInfo( );
        mInfo[i]->Kp = i%5;
        Lgm_Set_Coord_Transforms( Date, UTC, mInfo[i]->c );
        //printf("%8ldT%02d%02d%02d %g\n", mInfo[i]->c->UTC.Date, mInfo[i]->c->UTC.Hour, mInfo[i]->c->UTC.Minute, (int)mInfo[i]->c->UTC.Second, mInfo[i]->c->psi*DegPerRad );
        Lgm_MagModelInfo_Set_MagModel( LGM_IGRF, LGM_EXTMODEL_T89, mInfo[i] );
        mInfo[i]->P      = 4.1011111111111118;
        mInfo[i]->Dst    = 7.7777777777777777;
        mInfo[i]->By     = 3.7244444444444444;
        mInfo[i]->Bz     = -0.12666666666666665;
        mInfo[i]->W[0]   = 0.12244444444444445;
        mInfo[i]->W[1]   = 0.2514;
        mInfo[i]->W[2]   = 0.089266666666666661;
        mInfo[i]->W[3]   = 0.047866666666666668;
        mInfo[i]->W[4]   = 0.22586666666666666;
        mInfo[i]->W[5]   = 1.0461333333333334;
    }


    printf("%5s", "Kp");
    printf(" %10s", "Ux (Re)");
    printf(" %10s", "Uy (Re)");
    printf(" %10s", "Uz (Re)");
    printf(" %10s", "Bx (nT)");
    printf(" %10s", "By (nT)");
    printf(" %10s", "Bz (nT)");
    printf(" %10s", "Bmag (nT)");
    printf("   %15s", "(\u2207B)_x");
    printf("   %15s", "(\u2207B)_y");
    printf("   %15s", "(\u2207B)_z\n");


    for (i=0; i<=5; i++) {
        u.x = -6.6; u.y =  0.3;  u.z =  -2.0;
        mInfo[j]->Bfield( &u, &B, mInfo[j] );
        Lgm_CurlB( &u, &CurlB, LGM_DERIV_SIX_POINT, 1e-3, mInfo[j] );
        Lgm_GradB( &u, &GradB, LGM_DERIV_SIX_POINT, 1e-3, mInfo[j] );
        Lgm_GradBvec2( j, &u, GradBvec, LGM_DERIV_SIX_POINT, 1e-3, n, 10.0, &mInfo[0] );
        
        printf("\n\n");
        printf( "\t Kp        = %5d\n", mInfo[i]->Kp);
        printf( "\t u         = %10g %10g %10g\n", u.x, u.y, u.z );
        printf( "\t B         = %10g %10g %10g\n", B.x, B.y, B.z );
        printf( "\t |B|       = %10g\n", Lgm_Magnitude( &B ) );
        printf( "\t CurlB     = %15g %15g %15g\n", CurlB.x, CurlB.y, CurlB.z );
        printf( "\t GradB     = %15g %15g %15g\n", GradB.x, GradB.y, GradB.z );
        printf( "\t             (%8g %8g %8g %8g)\n", GradBvec[0][0], GradBvec[0][1], GradBvec[0][2], GradBvec[0][3] );
        printf( "\t Grad Bvec = (%8g %8g %8g %8g)\n", GradBvec[1][0], GradBvec[1][1], GradBvec[1][2], GradBvec[1][3] );
        printf( "\t             (%8g %8g %8g %8g)\n", GradBvec[2][0], GradBvec[2][1], GradBvec[2][2], GradBvec[2][3] );
        printf( "\t             (%8g %8g %8g %8g)\n", GradBvec[3][0], GradBvec[3][1], GradBvec[3][2], GradBvec[3][3] );
        double Bmag = Lgm_Magnitude( &B );
        printf("\n\t Show that GradB can be derived from GradBvec...\n");
        printf("\t GradB_T: %g\n", (B.x*GradBvec[1][0] + B.y*GradBvec[2][0] + B.z*GradBvec[3][0])/Bmag );
        printf("\t GradB_X: %g\n", (B.x*GradBvec[1][1] + B.y*GradBvec[2][1] + B.z*GradBvec[3][1])/Bmag );
        printf("\t GradB_Y: %g\n", (B.x*GradBvec[1][2] + B.y*GradBvec[2][2] + B.z*GradBvec[3][2])/Bmag );
        printf("\t GradB_Z: %g\n", (B.x*GradBvec[1][3] + B.y*GradBvec[2][3] + B.z*GradBvec[3][3])/Bmag );
        printf("\n\t Show that CurlB can be derived from GradBvec...\n");
        printf("\t CurlB_X: %g\n", GradBvec[3][2] - GradBvec[2][3] );
        printf("\t CurlB_Y: %g\n", GradBvec[1][3] - GradBvec[3][1] );
        printf("\t CurlB_Z: %g\n", GradBvec[2][1] - GradBvec[1][2] );
    }

/*
    for (j=0; j<100; j++){
        mInfo[i]->Kp = 3;
        for (i=0; i<13; i++) {
            u.x = -1. - (double)i * 0.5;
            u.y =  0.0;  u.z =  0.0;
            mInfo[i]->Bfield( &u, &B, mInfo[i] );
            Lgm_GradB( &u, &GradB, LGM_DERIV_SIX_POINT, 1e-3, mInfo[i] );
            printf( "%5d", mInfo[i]->Kp);
            printf( " %10g %10g %10g", u.x, u.y, u.z );
            printf( " %10g %10g %10g", B.x, B.y, B.z );
            printf( " %10g", Lgm_Magnitude( &B ) );
            printf( " %15g %15g %15g \n", GradB.x, GradB.y, GradB.z );
        }
    }
*/





    for ( i=0; i<8640; ++i ) {
        Lgm_FreeMagInfo( mInfo[i] );
    }
    LGM_ARRAY_1D_FREE( mInfo );


    exit(0);
}
