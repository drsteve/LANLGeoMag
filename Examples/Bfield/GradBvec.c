#include <Lgm_MagModelInfo.h>

/* BAL 04-Jan-2011 modified for more cases */

int main(){


    long int            Date;
    double              UTC;
    double              GradBvec[3][3];
    Lgm_Vector          u, B, GradB;
    Lgm_MagModelInfo    *mInfo;
    int                 i, j;


    Date = 20050831;
    UTC  = 9.0;

    mInfo = Lgm_InitMagInfo( );
    Lgm_Set_Coord_Transforms( Date, UTC, mInfo->c );

    mInfo->Bfield = Lgm_B_T89;
   mInfo->Bfield = Lgm_B_TS04;
//mInfo->Bfield = Lgm_B_OP77;
//mInfo->Bfield = Lgm_B_T89;
    mInfo->Bfield = Lgm_B_TS04;
    mInfo->P      = 4.1011111111111118;
    mInfo->Dst    = 7.7777777777777777;
    mInfo->By     = 3.7244444444444444;
    mInfo->Bz     = -0.12666666666666665;
    mInfo->W[0]   = 0.12244444444444445;
    mInfo->W[1]   = 0.2514;
    mInfo->W[2]   = 0.089266666666666661;
    mInfo->W[3]   = 0.047866666666666668;
    mInfo->W[4]   = 0.22586666666666666;
    mInfo->W[5]   = 1.0461333333333334;


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
        mInfo->Kp = i;
        u.x = -6.6; u.y =  0.0;  u.z =  0.0;
        mInfo->Bfield( &u, &B, mInfo );
        Lgm_GradB( &u, &GradB, LGM_DERIV_SIX_POINT, 1e-3, mInfo );
        Lgm_GradB( &u, GradBvec, LGM_DERIV_SIX_POINT, 1e-3, mInfo );
        printf( "%5d", mInfo->Kp);
        printf( " %10g %10g %10g", u.x, u.y, u.z );
        printf( " %10g %10g %10g", B.x, B.y, B.z );
        printf( " %10g", Lgm_Magnitude( &B ) );
        printf( " %15g %15g %15g \n", GradB.x, GradB.y, GradB.z );
        
        printf("\n\n");
        printf("             = (%8g %8g %8g)\n", GradBvec[0][0], GradBvec[0][1], GradBvec[0][2] );
        printf("   Grad Bvec = (%8g %8g %8g)\n", GradBvec[1][0], GradBvec[1][1], GradBvec[1][2] );
        printf("             = (%8g %8g %8g)\n", GradBvec[2][0], GradBvec[2][1], GradBvec[2][2] );
    }

/*
    for (j=0; j<100; j++){
        mInfo->Kp = 3;
        for (i=0; i<13; i++) {
            u.x = -1. - (double)i * 0.5;
            u.y =  0.0;  u.z =  0.0;
            mInfo->Bfield( &u, &B, mInfo );
            Lgm_GradB( &u, &GradB, LGM_DERIV_SIX_POINT, 1e-3, mInfo );
            printf( "%5d", mInfo->Kp);
            printf( " %10g %10g %10g", u.x, u.y, u.z );
            printf( " %10g %10g %10g", B.x, B.y, B.z );
            printf( " %10g", Lgm_Magnitude( &B ) );
            printf( " %15g %15g %15g \n", GradB.x, GradB.y, GradB.z );
        }
    }
*/





    Lgm_FreeMagInfo( mInfo );


    exit(0);
}
