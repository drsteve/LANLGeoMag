#include <Lgm_MagModelInfo.h>

int main(){


    long int            Date;
    double              UTC;
    Lgm_Vector          u, v;
    Lgm_MagModelInfo    *MagInfo;
    

    Date = 20100203;
    UTC  = 12.34567;

    MagInfo = Lgm_InitMagInfo( );
//    MagInfo->c->Verbose = 1;
    Lgm_Set_Coord_Transforms( Date, UTC, MagInfo->c );

    MagInfo->Bfield = Lgm_B_T01S;
    MagInfo->Bfield = Lgm_B_TS04;
    MagInfo->Bfield = Lgm_B_T89;
    MagInfo->Bfield = Lgm_B_OP77;
    MagInfo->Bfield = Lgm_B_cdip;
    MagInfo->Bfield = Lgm_B_igrf;
    MagInfo->Kp     = 5;

    MagInfo->VerbosityLevel = 5;



    u.x = -1.6; u.y =  0.0;  u.z =  0.0; // Re
    u.x = -1.601; u.y =  0.0;  u.z =  0.0; // Re
    u.x = -6.6; u.y =  0.0;  u.z =  0.0; // Re
    Lgm_TraceToEarth( &u, &v, 120.0, -1.0, 1e-7, MagInfo );
    printf( "u = %g %g %g Re\n", u.x, u.y, u.z );
    printf( "v = %g %g %g Re\n", v.x, v.y, v.z );
    printf( "Manitude(v) = %g\n", Lgm_Magnitude( &v ) );


    u.x = -2.6; u.y =  0.0;  u.z =  0.0; // Re
    Lgm_TraceToEarth( &u, &v, 120.0, -1.0, 1e-7, MagInfo );
    printf( "u = %g %g %g Re\n", u.x, u.y, u.z );
    printf( "v = %g %g %g Re\n", v.x, v.y, v.z );
    printf( "Manitude(v) = %g\n", Lgm_Magnitude( &v ) );


    Lgm_FreeMagInfo( MagInfo );


    exit(0);
}
