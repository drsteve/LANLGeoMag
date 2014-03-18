#include <Lgm_CTrans.h>

int main( ) {
    long int        sec, Date, JDN;
    int             s, m, smax, i, Year, Month, Day;
    Lgm_CTrans      *c = Lgm_init_ctrans( 0 );
    Lgm_DateTime    UTC, TAI, TDB, TT, GPS, GPS2, TAI2;
    double          u0, u1, t, GpsSeconds, TaiSeconds, Time, JD;



    // For comparison to SpacePy answer in SpacePy docs...
    printf("Lgm_JD( 2002, 2, 25, 12.0+20.0/60.0+30.0/3600.0, LGM_TIME_SYS_UTC, c  ) = %.15lf\n", (JD = Lgm_JD( 2002, 2, 25, 12.0+(20.0+30.0/60.0)/60.0, LGM_TIME_SYS_UTC, c )) );

    // Try to go backwards
/*
JD = 2452331.51423611111111111;
JD -= 2452331.0;
JDN = (long int) JD;
JD = JD - (double)JDN ;

JD = 52331.51423611111111111;
JD -= 52331.0;

JD *= 24.0;
printf( "%.10lf\n", JD );
JD -= 12.0;
JD *= 60.0;
printf( "%.10lf\n", JD );
JD -= 20.0;
JD *= 60.0;
printf( "%.10lf\n", JD );
exit(0);
*/
    Date = Lgm_JD_to_Date( JD, &Year, &Month, &Day, &Time );
    
    // Dump to UTC structure 
    Lgm_Make_UTC( Date, Time, &UTC, c );

    // Print the time...
    printf("\tUTC: "); Lgm_Print_DateTime( &UTC, 4, 8 ); printf("\n");


    // free structures
    Lgm_free_ctrans( c ); 

    return(0);

}

