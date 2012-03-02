#include <Lgm_MagModelInfo.h>

/* BAL 04-Jan-2011 modified for more cases */

int main(){


    long int            Date;
    double              UTC;
    Lgm_Vector          u, B, CurlB;
    Lgm_MagModelInfo    *mInfo;
    int                 i, j;


    Date = 20050831;
    UTC  = 9.0;

    mInfo = Lgm_InitMagInfo( );
    Lgm_Set_Coord_Transforms( Date, UTC, mInfo->c );

   mInfo->Bfield = Lgm_B_TS04;
//mInfo->Bfield = Lgm_B_OP77;
//mInfo->Bfield = Lgm_B_T89;
    mInfo->Bfield = Lgm_B_TS04_opt;
    mInfo->Bfield = Lgm_B_igrf;
    mInfo->Bfield = Lgm_B_T89;
    mInfo->Bfield = Lgm_B_cdip;
    mInfo->Bfield = Lgm_B_edip;
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


    printf("%13s", "Kp");
    printf(" %13s", "Ux (Re)");
    printf(" %13s", "Uy (Re)");
    printf(" %13s", "Uz (Re)");
    printf(" %13s", "Bx (nT)");
    printf(" %13s", "By (nT)");
    printf(" %13s", "Bz (nT)");
    printf(" %13s", "Bmag (nT)");
    printf(" %16s", "CurlB_x (nT/Re)");
    printf(" %16s", "CurlB_y (nT/Re)");
    printf(" %16s", "CurlB_z (nT/Re)\n");


    for (i=0; i<=5; i++) {
        mInfo->Kp = i;
        u.x = -6.6; u.y =  0.0;  u.z =  0.0;
        mInfo->Bfield( &u, &B, mInfo );
        Lgm_CurlB( &u, &CurlB, LGM_DERIV_SIX_POINT, 1e-3, mInfo );
        printf( "%13i", mInfo->Kp);
        printf( " %13g %13g %13g", u.x, u.y, u.z );
        printf( " %13g %13g %13g", B.x, B.y, B.z );
        printf( " %13g", Lgm_Magnitude( &B ) );
        printf( " %16g %16g %16g \n", CurlB.x, CurlB.y, CurlB.z );
    }

    for (j=0; j<100; j++){
        mInfo->Kp = 3;
        for (i=0; i<13; i++) {
            u.x = -1. - (double)i * 0.5;
            u.y =  1.0;  u.z =  -1.0;
            mInfo->Bfield( &u, &B, mInfo );
            Lgm_CurlB( &u, &CurlB, LGM_DERIV_SIX_POINT, 1e-3, mInfo );
            printf( "%13i", mInfo->Kp);
            printf( " %13g %13g %13g", u.x, u.y, u.z );
            printf( " %13g %13g %13g", B.x, B.y, B.z );
            printf( " %13g", Lgm_Magnitude( &B ) );
            printf( " %16g %16g %16g \n", CurlB.x, CurlB.y, CurlB.z );
        }
    }





    Lgm_FreeMagInfo( mInfo );


    exit(0);
}
