#include <Lgm_MagModelInfo.h>
#include <Lgm_QinDenton.h>

/* BAL 04-Jan-2011 modified for more cases */

int main(){


    long int            Date;
    double              UTC, JD;
    Lgm_Vector          u, w, B;
    Lgm_MagModelInfo    *mInfo;
    Lgm_QinDentonOne    p;
    int                 i, j;

    mInfo = Lgm_InitMagInfo( );

    Date  = 20040812;                        // August 12, 2004
    UTC   = 16.0 + 27.0/60.0 + 29.99/3600.0;   // Universal Time Coordinated (in decimal hours)
    UTC   = 0.0;
    JD    = Lgm_Date_to_JD( Date, UTC, mInfo->c );    // Compute JD

    Lgm_Set_Coord_Transforms( Date, UTC, mInfo->c );


    // Get (interpolate) the QinDenton vals from the values in the file at the given Julian Date
    Lgm_get_QinDenton_at_JD( JD, &p, 1, 0 );

    // Set params in mInfo structure.
    Lgm_set_QinDenton( &p, mInfo );

    
    u.x = -6.6; u.y = 0.0; u.z = 0.0;
//    u.x = 1.490696277499862; u.y =  0.000943013537158; u.z =  0.051705206498810;
//mInfo->Kp = 3;

    Lgm_MagModelInfo_Set_MagModel( LGM_CDIP, LGM_EXTMODEL_T89, mInfo );
    mInfo->Bfield( &u, &B, mInfo );
    printf( "%20.14lf%20.14lf%20.14lf", u.x, u.y, u.z );
    printf( "%20.14lf%20.14lf%20.14lf", B.x, B.y, B.z );
    printf( "%20.14lf\n", Lgm_Magnitude( &B ) );

    Lgm_MagModelInfo_Set_MagModel( LGM_EDIP, LGM_EXTMODEL_T89, mInfo );
    mInfo->Bfield( &u, &B, mInfo );
    printf( "%20.14lf%20.14lf%20.14lf", u.x, u.y, u.z );
    printf( "%20.14lf%20.14lf%20.14lf", B.x, B.y, B.z );
    printf( "%20.14lf\n", Lgm_Magnitude( &B ) );

    Lgm_MagModelInfo_Set_MagModel( LGM_IGRF, LGM_EXTMODEL_T89, mInfo );
    mInfo->Bfield( &u, &B, mInfo );
    printf( "%20.14lf%20.14lf%20.14lf", u.x, u.y, u.z );
    printf( "%20.14lf%20.14lf%20.14lf", B.x, B.y, B.z );
    printf( "%20.14lf\n", Lgm_Magnitude( &B ) );

    Lgm_MagModelInfo_Set_MagModel( LGM_IGRF, LGM_EXTMODEL_TS04, mInfo );
    mInfo->Bfield( &u, &B, mInfo );
    printf( "%20.14lf%20.14lf%20.14lf", u.x, u.y, u.z );
    printf( "%20.14lf%20.14lf%20.14lf", B.x, B.y, B.z );
    printf( "%20.14lf\n", Lgm_Magnitude( &B ) );

    Lgm_MagModelInfo_Set_MagModel( LGM_IGRF, LGM_EXTMODEL_T96, mInfo );
    mInfo->Bfield( &u, &B, mInfo );
    printf( "%20.14lf%20.14lf%20.14lf", u.x, u.y, u.z );
    printf( "%20.14lf%20.14lf%20.14lf", B.x, B.y, B.z );
    printf( "%20.14lf\n", Lgm_Magnitude( &B ) );

    Lgm_MagModelInfo_Set_MagModel( LGM_IGRF, LGM_EXTMODEL_T01S, mInfo );
    mInfo->Bfield( &u, &B, mInfo );
    printf( "%20.14lf%20.14lf%20.14lf", u.x, u.y, u.z );
    printf( "%20.14lf%20.14lf%20.14lf", B.x, B.y, B.z );
    printf( "%20.14lf\n", Lgm_Magnitude( &B ) );



    printf("\n\n\n");

    /* 
     * At various geodetic points
     */
    Lgm_MagModelInfo_Set_MagModel( LGM_IGRF, LGM_EXTMODEL_NULL, mInfo );
    double  GeodLat = 35.0;
    double  GeodLon = 123.0;
    double  GeodAlt = 0.0;
    Lgm_GEOD_to_WGS84( GeodLat, GeodLon, GeodAlt, &w );

//w.x = 7000.0/Re*cos(123.0*RadPerDeg)*cos(35.0*RadPerDeg);
//////////////////w.y = 7000.0/Re*sin(123.0*RadPerDeg)*cos(35.0*RadPerDeg);
//////////////////w.z = 7000.0/Re*sin(35.0*RadPerDeg);
    Lgm_Convert_Coords( &w, &u, WGS84_TO_GSM, mInfo->c );
    mInfo->Bfield( &u, &B, mInfo );
    printf( "Geodetic (lat,lon,alt): %20.14lf  %20.14lf  %20.14lf\n", GeodLat, GeodLon, GeodAlt );
    printf( "       Bx,By,Bz in GSM: %20.14lf  %20.14lf  %20.14lf\n", B.x, B.y, B.z );
    printf( "                   |B|: %20.14lf\n", Lgm_Magnitude( &B ) );

    printf( "   http://www.ngdc.noaa.gov/geomag-web/#igrfwmm value for |B|: %20.14lf\n",  50371.8 );





    Lgm_FreeMagInfo( mInfo );


    exit(0);
}
