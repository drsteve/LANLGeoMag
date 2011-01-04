#include <Lgm_MagModelInfo.h>

/* BAL 4-Jan2013 modified for more cases */

int main(){


    long int            Date;
    double              UTC;
    Lgm_Vector          u, B;
    Lgm_MagModelInfo    *MagInfo;
    int i, j;


    Date = 20050831;
    UTC  = 9.0;

    MagInfo = Lgm_InitMagInfo( );
    Lgm_Set_Coord_Transforms( Date, UTC, MagInfo->c );

    MagInfo->Bfield = Lgm_B_T89;


    printf("%6s", "Kp");
    printf("%13s", "Ux (Re)");
    printf("%13s", "Uy (Re)");
    printf("%13s", "Uz (Re)");
    printf("%13s", "Bx (nT)");
    printf("%13s", "By (nT)");
    printf("%13s", "Bz (nT)");
    printf("%13s", "Bmag (nT)\n");


    for (i=0; i<13; i++) {
        MagInfo->Kp = i;
        u.x = -6.6; u.y =  0.0;  u.z =  0.0;
        MagInfo->Bfield( &u, &B, MagInfo );
        printf( "%6i", MagInfo->Kp);
        printf( "%13g%13g%13g", u.x, u.y, u.z );
        printf( "%13g%13g%13g", B.x, B.y, B.z );
        printf( "%13g\n", Lgm_Magnitude( &B ) );
    }

    MagInfo->Kp = 3;
    for (i=0; i<13; i++) {
        u.x = -1. - (float)i * 0.5;
        u.y =  0.0;  u.z =  0.0;
        MagInfo->Bfield( &u, &B, MagInfo );
        printf( "%6i", MagInfo->Kp);
        printf( "%13g%13g%13g", u.x, u.y, u.z );
        printf( "%13g%13g%13g", B.x, B.y, B.z );
        printf( "%13g\n", Lgm_Magnitude( &B ) );
    }







    Lgm_FreeMagInfo( MagInfo );


    exit(0);
}
