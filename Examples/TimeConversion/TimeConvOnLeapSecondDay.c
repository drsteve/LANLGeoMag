#include <Lgm_CTrans.h>

int main( ) {
    long int        sec;
    int             s, m, smax, i;
    Lgm_CTrans      *c = Lgm_init_ctrans( 0 );
    Lgm_DateTime    UTC, TAI, TDB;
    double          u0, u1, t;




    printf("JD( 1985, 6, 30, 12.3456, c  ) = %lf\n", Lgm_JD( 1985, 6, 30, 12.3456, LGM_TIME_SYS_UTC, c ) );





/*
    for ( m=58; m<=59; m++){
        smax = (m==59) ? 59 + Lgm_IsLeapSecondDay( 19850701, c ) : 59;
        for ( s=0; s<=smax; s++){

            Lgm_Make_UTC( 19850701, 23.0 + ((double)m + (double)s/60.0)/60.0, &UTC, c );

            //printf("UTC.Date = %8ld   UTC.Time = %lf\n", UTC.Date, UTC.Time );
            printf("\r%02d %02d   ", m, s ); Lgm_Print_DateTime( UTC, 4, 8 ); fflush(NULL); printf("");
            //printf("TAI.Date = %8ld   TAI.Time = %lf\n", TAI.Date, TAI.Time );
            usleep(200000);
        }
    }
*/

    //Lgm_Make_UTC( 19850701, 23.0 + (23.0 + 34.0/60.0)/60.0, &UTC, c );


    printf( "A leap second was added at the end of the day of June 30, 1985.  This addition\n"
            "increased the leap seconds count from 22 to 23. I.e.  on June 30, 1985 the\n"
            "count is still 22, but on July 1, 1985 it increases to 23. This test confirms\n"
            "that times are properly computed around this date.\n\n");

    printf("UTC Date/Time 0 seconds into June 29, 1985\n");
    Lgm_Make_UTC( 19850629, 0.0/3600.0, &UTC, c );
    printf("\tUTC: "); Lgm_Print_DateTime( &UTC, 4, 8 ); printf("\n");
    Lgm_UTC_to_TAI( &UTC, &TAI, c );
    printf("\tTAI: "); Lgm_Print_DateTime( &TAI, 4, 8 ); printf("\n");
    Lgm_TAI_to_UTC( &TAI, &UTC, c );
    printf("\tUTC: "); Lgm_Print_DateTime( &UTC, 4, 8 ); printf("\n\n");

    printf("UTC Date/Time 43200 seconds into June 29, 1985\n");
    Lgm_Make_UTC( 19850629, 43200.0/3600.0, &UTC, c );
    printf("\tUTC: "); Lgm_Print_DateTime( &UTC, 4, 8 ); printf("\n");
    Lgm_UTC_to_TAI( &UTC, &TAI, c );
    printf("\tTAI: "); Lgm_Print_DateTime( &TAI, 4, 8 ); printf("\n");
    Lgm_TAI_to_UTC( &TAI, &UTC, c );
    printf("\tUTC: "); Lgm_Print_DateTime( &UTC, 4, 8 ); printf("\n\n");

    printf("UTC Date/Time 86399.5 seconds into June 29, 1985\n");
    Lgm_Make_UTC( 19850629, 86399.5/3600.0, &UTC, c );
    printf("\tUTC: "); Lgm_Print_DateTime( &UTC, 4, 8 ); printf("\n");
    Lgm_UTC_to_TAI( &UTC, &TAI, c );
    printf("\tTAI: "); Lgm_Print_DateTime( &TAI, 4, 8 ); printf("\n");
    Lgm_TAI_to_UTC( &TAI, &UTC, c );
    printf("\tUTC: "); Lgm_Print_DateTime( &UTC, 4, 8 ); printf("\n\n");

    printf("UTC Date/Time 86400.0 seconds into June 29, 1985\n");
    Lgm_Make_UTC( 19850629, 86400.0/3600.0, &UTC, c );
    printf("\tUTC: "); Lgm_Print_DateTime( &UTC, 4, 8 ); printf("\n");
    Lgm_UTC_to_TAI( &UTC, &TAI, c );
    printf("\tTAI: "); Lgm_Print_DateTime( &TAI, 4, 8 ); printf("\n");
    Lgm_TAI_to_UTC( &TAI, &UTC, c );
    printf("\tUTC: "); Lgm_Print_DateTime( &UTC, 4, 8 ); printf("\n\n");

    printf("UTC Date/Time 86400.5 seconds into June 29, 1985\n");
    Lgm_Make_UTC( 19850629, 86400.5/3600.0, &UTC, c );
    printf("\tUTC: "); Lgm_Print_DateTime( &UTC, 4, 8 ); printf("\n");
    Lgm_UTC_to_TAI( &UTC, &TAI, c );
    printf("\tTAI: "); Lgm_Print_DateTime( &TAI, 4, 8 ); printf("\n");
    Lgm_TAI_to_UTC( &TAI, &UTC, c );
    printf("\tUTC: "); Lgm_Print_DateTime( &UTC, 4, 8 ); printf("\n\n");

    printf("\n\n");

    printf("UTC Date/Time 0 seconds into June 30, 1985\n");
    Lgm_Make_UTC( 19850630, 0.0/3600.0, &UTC, c );
    printf("\tUTC: "); Lgm_Print_DateTime( &UTC, 4, 8 ); printf("\n");
    Lgm_UTC_to_TAI( &UTC, &TAI, c );
    printf("\tTAI: "); Lgm_Print_DateTime( &TAI, 4, 8 ); printf("\n");
    Lgm_TAI_to_UTC( &TAI, &UTC, c );
    printf("\tUTC: "); Lgm_Print_DateTime( &UTC, 4, 8 ); printf("\n\n");

    printf("UTC Date/Time 43200 seconds into June 30, 1985\n");
    Lgm_Make_UTC( 19850630, 43200.0/3600.0, &UTC, c );
    printf("\tUTC: "); Lgm_Print_DateTime( &UTC, 4, 8 ); printf("\n");
    Lgm_UTC_to_TAI( &UTC, &TAI, c );
    printf("\tTAI: "); Lgm_Print_DateTime( &TAI, 4, 8 ); printf("\n");
    Lgm_TAI_to_UTC( &TAI, &UTC, c );
    printf("\tUTC: "); Lgm_Print_DateTime( &UTC, 4, 8 ); printf("\n\n");

    printf("UTC Date/Time 86399.5 seconds into June 30, 1985\n");
    Lgm_Make_UTC( 19850630, 86399.5/3600.0, &UTC, c );
    printf("\tUTC: "); Lgm_Print_DateTime( &UTC, 4, 8 ); printf("\n");
    Lgm_UTC_to_TAI( &UTC, &TAI, c );
    printf("\tTAI: "); Lgm_Print_DateTime( &TAI, 4, 8 ); printf("\n");
    Lgm_TAI_to_UTC( &TAI, &UTC, c );
    printf("\tUTC: "); Lgm_Print_DateTime( &UTC, 4, 8 ); printf("\n\n");

    printf("UTC Date/Time 86400.0 seconds into June 30, 1985\n");
    Lgm_Make_UTC( 19850630, 86400.0/3600.0, &UTC, c );
    printf("\tUTC: "); Lgm_Print_DateTime( &UTC, 4, 8 ); printf("\n");
    Lgm_UTC_to_TAI( &UTC, &TAI, c );
    printf("\tTAI: "); Lgm_Print_DateTime( &TAI, 4, 8 ); printf("\n");
    Lgm_TAI_to_UTC( &TAI, &UTC, c );
    printf("\tUTC: "); Lgm_Print_DateTime( &UTC, 4, 8 ); printf("\n\n");

    printf("UTC Date/Time 86400.5 seconds into June 30, 1985\n");
    Lgm_Make_UTC( 19850630, 86400.5/3600.0, &UTC, c );
    printf("\tUTC: "); Lgm_Print_DateTime( &UTC, 4, 8 ); printf("\n");
    Lgm_UTC_to_TAI( &UTC, &TAI, c );
    printf("\tTAI: "); Lgm_Print_DateTime( &TAI, 4, 8 ); printf("\n");
Lgm_UTC_to_TDB( &UTC, &TDB,c );
printf("\tTDB: "); Lgm_Print_DateTime( &TDB, 4, 8 ); printf("\n\n");
    Lgm_TAI_to_UTC( &TAI, &UTC, c );
    printf("\tUTC: "); Lgm_Print_DateTime( &UTC, 4, 8 ); printf("\n\n");

    printf("\n\n");

    printf("UTC Date/Time 0 seconds into July 1, 1985\n");
    Lgm_Make_UTC( 19850701, 0.0/3600.0, &UTC, c );
    printf("\tUTC: "); Lgm_Print_DateTime( &UTC, 4, 8 ); printf("\n");
    Lgm_UTC_to_TAI( &UTC, &TAI, c );
    printf("\tTAI: "); Lgm_Print_DateTime( &TAI, 4, 8 ); printf("\n");
    Lgm_TAI_to_UTC( &TAI, &UTC, c );
    printf("\tUTC: "); Lgm_Print_DateTime( &UTC, 4, 8 ); printf("\n\n");

    printf("UTC Date/Time 43200 seconds into July 1, 1985\n");
    Lgm_Make_UTC( 19850701, 43200.0/3600.0, &UTC, c );
    printf("\tUTC: "); Lgm_Print_DateTime( &UTC, 4, 8 ); printf("\n");
    Lgm_UTC_to_TAI( &UTC, &TAI, c );
    printf("\tTAI: "); Lgm_Print_DateTime( &TAI, 4, 8 ); printf("\n");
    Lgm_TAI_to_UTC( &TAI, &UTC, c );
    printf("\tUTC: "); Lgm_Print_DateTime( &UTC, 4, 8 ); printf("\n\n");

    printf("UTC Date/Time 86399.5 seconds into July 1, 1985\n");
    Lgm_Make_UTC( 19850701, 86399.5/3600.0, &UTC, c );
    printf("\tUTC: "); Lgm_Print_DateTime( &UTC, 4, 8 ); printf("\n");
    Lgm_UTC_to_TAI( &UTC, &TAI, c );
    printf("\tTAI: "); Lgm_Print_DateTime( &TAI, 4, 8 ); printf("\n");
    Lgm_TAI_to_UTC( &TAI, &UTC, c );
    printf("\tUTC: "); Lgm_Print_DateTime( &UTC, 4, 8 ); printf("\n\n");

    printf("UTC Date/Time 86400.0 seconds into July 1, 1985\n");
    Lgm_Make_UTC( 19850701, 86400.0/3600.0, &UTC, c );
    printf("\tUTC: "); Lgm_Print_DateTime( &UTC, 4, 8 ); printf("\n");
    Lgm_UTC_to_TAI( &UTC, &TAI, c );
    printf("\tTAI: "); Lgm_Print_DateTime( &TAI, 4, 8 ); printf("\n");
    Lgm_TAI_to_UTC( &TAI, &UTC, c );
    printf("\tUTC: "); Lgm_Print_DateTime( &UTC, 4, 8 ); printf("\n\n");

    printf("UTC Date/Time 86400.5 seconds into July 1, 1985\n");
    Lgm_Make_UTC( 19850701, 86400.5/3600.0, &UTC, c );
    printf("\tUTC: "); Lgm_Print_DateTime( &UTC, 4, 8 ); printf("\n");
    Lgm_UTC_to_TAI( &UTC, &TAI, c );
    printf("\tTAI: "); Lgm_Print_DateTime( &TAI, 4, 8 ); printf("\n");
    Lgm_TAI_to_UTC( &TAI, &UTC, c );
    printf("\tUTC: "); Lgm_Print_DateTime( &UTC, 4, 8 ); printf("\n\n");



    /*
     * Lets go TAI to UTCs...
     */
    Lgm_DateTime *UTC2;
    UTC2 = Lgm_DateTime_Create( 2004, 5, 14, 23.0 + 59.0/60.0 + 43.0/3600.0, LGM_TIME_SYS_UTC, c );
    Lgm_Make_UTC( 20040514, 16.0 + 43.0/60.0, &UTC, c );
    UTC.JD = Lgm_JD( TAI.Year, TAI.Month, TAI.Day, TAI.Time, LGM_TIME_SYS_TAI, c );
    Lgm_TAI_to_UTC( &TAI, &UTC, c );
    printf("\tUTC: "); Lgm_Print_DateTime( &UTC2, 4, 8 ); printf("\n\n");
    


    Lgm_UTC_to_TDB( UTC2, &TDB,c );
    printf("\tTDB: "); Lgm_Print_DateTime( &TDB, 4, 8 ); printf("\n\n");
    
    



    

    
    
    // free structures
    Lgm_free_ctrans( c ); 
    free( UTC2 );

    return(0);

}

