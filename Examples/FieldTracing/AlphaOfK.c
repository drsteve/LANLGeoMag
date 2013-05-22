#include <Lgm_MagModelInfo.h>

int main(){


    long int            Date;
    double              UTC, val;
    Lgm_Vector          u, v;
    Lgm_MagModelInfo    *MagInfo;
    Lgm_DateTime        *d;
    

    Date = 20010101;
    UTC  = 0.0; //12.34567;

    MagInfo = Lgm_InitMagInfo( );
    d = Lgm_DateTime_Create( 2001, 1, 1, 0.0, LGM_TIME_SYS_UTC, MagInfo->c );
//    MagInfo->c->Verbose = 1;
    Lgm_Set_Coord_Transforms( Date, UTC, MagInfo->c );

    MagInfo->Bfield = Lgm_B_T89;
    MagInfo->Kp     = 3;

    MagInfo->VerbosityLevel = 1;

    u.x = -6.6; u.y =  0.0;  u.z =  0.0; // Re
    Lgm_Convert_Coords(&u, &v, SM_TO_GSM, MagInfo->c);
    Lgm_Setup_AlphaOfK( d, &v, MagInfo );
    val = Lgm_AlphaOfK( 0.31, MagInfo );
    printf( "Alpha = %g \n", val );
    Lgm_TearDown_AlphaOfK( MagInfo );

//    u.x = -1.6; u.y =  0.0;  u.z =  0.0; // Re
//    u.x = -1.601; u.y =  0.0;  u.z =  0.0; // Re
//    u.x = -6.6; u.y =  0.0;  u.z =  0.0; // Re
//    Lgm_TraceToEarth( &u, &v, 120.0, -1.0, 1e-7, MagInfo );
//    printf( "u = %g %g %g Re\n", u.x, u.y, u.z );
//    printf( "v = %g %g %g Re\n", v.x, v.y, v.z );
//    printf( "Manitude(v) = %g\n", Lgm_Magnitude( &v ) );
//
//
//    u.x = -2.6; u.y =  0.0;  u.z =  0.0; // Re
//    Lgm_TraceToEarth( &u, &v, 120.0, -1.0, 1e-7, MagInfo );
//    printf( "u = %g %g %g Re\n", u.x, u.y, u.z );
//    printf( "v = %g %g %g Re\n", v.x, v.y, v.z );
//    printf( "Manitude(v) = %g\n", Lgm_Magnitude( &v ) );


    Lgm_FreeMagInfo( MagInfo );


    exit(0);
}
