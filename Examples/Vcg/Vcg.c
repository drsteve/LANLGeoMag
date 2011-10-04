#include <Lgm_MagModelInfo.h>

/* BAL 04-Jan-2011 modified for more cases */

int main(){


    long int            Date;
    double              UTC, B, T, q, Alpha, SinAlpha, Bm, E0;
    Lgm_Vector          u, w, Bvec, Vcg;
    Lgm_MagModelInfo    *mInfo;
    int                 i;


    Date = 20050831;
    UTC  = 9.0;

    mInfo = Lgm_InitMagInfo( );
    Lgm_Set_Coord_Transforms( Date, UTC, mInfo->c );

    mInfo->Bfield = Lgm_B_TS04;
    mInfo->Bfield = Lgm_B_OP77;
    mInfo->Bfield = Lgm_B_TS04_opt;
    mInfo->Bfield = Lgm_B_cdip;
    mInfo->Bfield = Lgm_B_T89;
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


    printf(" %10s", "usm_x (Re)");
    printf(" %10s", "usm_y (Re)");
    printf(" %10s", "usm_z (Re)");
    printf(" %10s", "Bx (nT)");
    printf(" %10s", "By (nT)");
    printf(" %10s", "Bz (nT)");
    printf(" %10s", "Bmag (nT)");
    printf(" %12s", "      (Vcg)x");
    printf(" %12s", "      (Vcg)y");
    printf(" %12s", "      (Vcg)z\n");


    for (i=0; i<=5; i++) {

        // -6.6 in SM Eq Plane.
        u.x = -2.6-(double)i; u.y = 0.0;  u.z =  0.0;
        Lgm_Convert_Coords( &u, &w, SM_TO_GSM, mInfo->c );
        mInfo->Bfield( &w, &Bvec, mInfo );
        B = Lgm_Magnitude( &Bvec );

        Alpha    = 90.0;    // Degrees
        q        = -LGM_e;   // Electron charge
        E0       = LGM_Ee0; // Electron rest energy
        T        = 0.1;    // 100 keV Kinetic Energy
        SinAlpha = sin(Alpha*RadPerDeg);
        Bm       = B/(SinAlpha*SinAlpha);

        Lgm_GradAndCurvDriftVel( &w, &Vcg,  q, T, E0, Bm, LGM_DERIV_SIX_POINT, 1e-3, mInfo );

        printf( " %10g %10g %10g", u.x, u.y, u.z );
        printf( " %10g %10g %10g", Bvec.x, Bvec.y, Bvec.z );
        printf( " %10g", B );
        printf( " %12g %12g %12g\n",   Vcg.x, Vcg.y, Vcg.z );
    }





    Lgm_FreeMagInfo( mInfo );


    exit(0);
}
