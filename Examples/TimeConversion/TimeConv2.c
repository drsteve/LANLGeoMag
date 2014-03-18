#include <Lgm_CTrans.h>

int main( ) {
    long int        sec;
    int             s, m, smax, i;
    Lgm_CTrans      *c = Lgm_init_ctrans( 0 );
    Lgm_DateTime    UTC, TAI, TDB, TT, GPS, GPS2, TAI2;
    double          u0, u1, t, GpsTime, TaiSeconds;




    /*
     * Lets go UTC -> GPS_Time
     */
    printf("Tests TaiSeconds <-> UTC over a Leap Seconds Day boundary\n");
    printf("------------------------------------------------------\n\n");
    printf("    Starting with 20081231 23:59:40 UTC\n");
    Lgm_Make_UTC( 20081231, 86380.0/3600.0, &UTC, c );
    printf("\tUTC: "); Lgm_Print_DateTime( &UTC, 0, 8 ); printf("\n\n");

    Lgm_UTC_to_TAI( &UTC, &TAI, c );
    printf("    Results of Lgm_UTC_to_TAI\n");
    printf("\tTAI: "); Lgm_Print_DateTime( &TAI, 0, 8 ); printf("\n\n");

    TaiSeconds = Lgm_UTC_to_TaiSeconds( &UTC, c );
    printf("    Results of Lgm_UTC_to_TaiSeconds\n");
    printf("\tTaiSeconds: %.8lf Seconds\n\n", TaiSeconds );

    Lgm_TaiSeconds_to_TAI( TaiSeconds, &TAI2 );
    printf("    Results of Lgm_TaiSeconds_to_TAI\n");
    printf("\tTAI: "); Lgm_Print_DateTime( &TAI2, 0, 8 ); printf("\n\n");

    Lgm_TaiSeconds_to_UTC( TaiSeconds, &UTC, c );
    printf("    Results of Lgm_TaiSeconds_to_UTC\n");
    printf("\tUTC: "); Lgm_Print_DateTime( &UTC, 0, 8 ); printf("\n\n");




    
    
    // free structures
    Lgm_free_ctrans( c ); 

    return(0);

}

