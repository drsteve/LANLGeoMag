#include <Lgm_CTrans.h>

int main( ) {
    long int        sec;
    int             s, m, smax, i;
    Lgm_CTrans      *c = Lgm_init_ctrans( 0 );
    Lgm_DateTime    UTC, TAI, TDB, TT, GPS, GPS2, TAI2;
    double          u0, u1, t, GpsSeconds, TaiSeconds;




    printf("JD( 1985, 6, 30, 12.3456, c  ) = %lf\n", Lgm_JD( 1985, 6, 30, 12.3456, LGM_TIME_SYS_UTC, c ) );




    printf( "Test time conversions forward and back again\n\n" );

    /*
     * Lets go UTC -> TDB and back
     */
    printf("Test UTC -> TDB\n");
    Lgm_Make_UTC( 20040514, 16.0 + 43.0/60.0, &UTC, c );
    printf("\tUTC: "); Lgm_Print_DateTime( UTC, 4, 8 ); printf("\n");

    Lgm_UTC_to_TDB( &UTC, &TDB, c );
    printf("\tTDB: "); Lgm_Print_DateTime( TDB, 4, 8 ); printf("\n");
    
    Lgm_TDB_to_UTC( &TDB, &UTC, c );
    printf("\tUTC: "); Lgm_Print_DateTime( UTC, 4, 8 ); printf("\n");
    
    printf("\n\n");
    
    /*
     * Lets go UTC -> TAI and back
     */
    printf("Test UTC -> TAI\n");
    Lgm_Make_UTC( 20040514, 16.0 + 43.0/60.0, &UTC, c );
    printf("\tUTC: "); Lgm_Print_DateTime( UTC, 4, 8 ); printf("\n");

    Lgm_UTC_to_TAI( &UTC, &TAI, c );
    printf("\tTAI: "); Lgm_Print_DateTime( TAI, 4, 8 ); printf("\n");
    
    Lgm_TAI_to_UTC( &TAI, &UTC, c );
    printf("\tUTC: "); Lgm_Print_DateTime( UTC, 4, 8 ); printf("\n");
    

    printf("\n\n");
    
    /*
     * Lets go UTC -> GPS and back
     */
    printf("Test UTC -> GPS\n");
    Lgm_Make_UTC( 20040514, 16.0 + 43.0/60.0, &UTC, c );
    printf("\tUTC: "); Lgm_Print_DateTime( UTC, 4, 8 ); printf("\n");

    Lgm_UTC_to_GPS( &UTC, &GPS, c );
    printf("\tGPS: "); Lgm_Print_DateTime( GPS, 4, 8 ); printf("\n");
    
    Lgm_GPS_to_UTC( &GPS, &UTC, c );
    printf("\tUTC: "); Lgm_Print_DateTime( UTC, 4, 8 ); printf("\n");
    

    printf("\n\n");
    
    /*
     * Lets go UTC -> TT and back
     */
    printf("Test UTC -> TT\n");
    Lgm_Make_UTC( 20040514, 16.0 + 43.0/60.0, &UTC, c );
    printf("\tUTC: "); Lgm_Print_DateTime( UTC, 4, 8 ); printf("\n");

    Lgm_UTC_to_TT( &UTC, &TT, c );
    printf("\tTT:  "); Lgm_Print_DateTime( TT, 4, 8 ); printf("\n");
    
    Lgm_TT_to_UTC( &TT, &UTC, c );
    printf("\tUTC: "); Lgm_Print_DateTime( UTC, 4, 8 ); printf("\n");


    printf("\n\n");


    /*
     * Lets go UTC -> GPS_Time
     */
    printf("Tests of UTC <-> GpsSeconds and UTC <-> TaiSeconds\n");
    Lgm_Make_UTC( 20100526, 20.0 + 40.0/60.0 + 23.0/3600.0, &UTC, c );
    printf("\tUTC: "); Lgm_Print_DateTime( UTC, 4, 8 ); printf("\n");

    Lgm_UTC_to_TAI( &UTC, &TAI, c );
    printf("\tTAI: "); Lgm_Print_DateTime( TAI, 4, 8 ); printf("\n");

    Lgm_UTC_to_GPS( &UTC, &GPS, c );
    printf("\tGPS: "); Lgm_Print_DateTime( GPS, 4, 8 ); printf("\n");

    printf("\n\n");

    printf("Test UTC <-> GpsSeconds\n");
    GpsSeconds = Lgm_UTC_to_GpsSeconds( &UTC, c );
    printf("\tGpsSeconds: %.8lf Seconds\n", GpsSeconds );

    Lgm_GpsSeconds_to_GPS( GpsSeconds, &GPS2 );
    printf("\tGPS: "); Lgm_Print_DateTime( GPS2, 4, 8 ); printf("\n");

    Lgm_GpsSeconds_to_UTC( GpsSeconds, &UTC, c );
    printf("\tUTC: "); Lgm_Print_DateTime( UTC, 4, 8 ); printf("\n");

    printf("\n\n");


    printf("Test UTC <-> TaiSeconds\n");
    TaiSeconds = Lgm_UTC_to_TaiSeconds( &UTC, c );
    printf("\tTaiSeconds: %.8lf Seconds\n", TaiSeconds );

    Lgm_TaiSeconds_to_TAI( TaiSeconds, &TAI2 );
    printf("\tTAI: "); Lgm_Print_DateTime( TAI2, 4, 8 ); printf("\n");

    Lgm_TaiSeconds_to_UTC( TaiSeconds, &UTC, c );
    printf("\tUTC: "); Lgm_Print_DateTime( UTC, 4, 8 ); printf("\n");

    
    
    // free structures
    Lgm_free_ctrans( c ); 

    return(0);

}

