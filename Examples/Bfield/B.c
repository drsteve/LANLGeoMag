#include <Lgm_MagModelInfo.h>

int main(){


    long int            Date;
    double              UTC;
    Lgm_Vector          u, B;
    Lgm_MagModelInfo    *MagInfo;
    


    Date = 20050831;
    UTC  = 9.0;

    MagInfo = Lgm_InitMagInfo( );
    Lgm_Set_Coord_Transforms( Date, UTC, MagInfo->c );

    MagInfo->Bfield = Lgm_B_T89;
    MagInfo->Kp     = 5;

    u.x = -6.6; u.y =  0.0;  u.z =  0.0;
    MagInfo->Bfield( &u, &B, MagInfo );
    printf( "u = %g %g %g Re\n", u.x, u.y, u.z );
    printf( "B = %g %g %g nT\n", B.x, B.y, B.z );
    printf( "Manitude(B) = %g\n", Lgm_Magnitude( &B ) );


    Lgm_FreeMagInfo( MagInfo );


    exit(0);
}
