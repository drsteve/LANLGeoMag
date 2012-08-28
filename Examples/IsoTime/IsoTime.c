#include <Lgm_CTrans.h>

/*
 * Shows how to get time into ISO format.
 */
int main( ) {

    Lgm_CTrans      *c = Lgm_init_ctrans( 0 );
    Lgm_DateTime    UTC;
    long int        Date;
    double          Time;
    char            Str[80];


    /*
     *  Fill DateTime object (i.e. a structure) for give date and time
     */
    Date = 20011031;
    Time = 6.5;
    Lgm_Make_UTC( Date, Time, &UTC, c );


    /*
     *  Print IsoTime like this. Note: Lgm_Print_DateTime() doesnt print a newline.
     */
    Lgm_Print_DateTime( &UTC, 0, 3 );
    printf("\n");


    /*
     *  Or like this.
     */
    Lgm_DateTimeToString( Str, &UTC, 0, 3 );
    printf("Date: %ld Time: %lf  in ISO 8601 format is: %s\n", Date, Time, Str );


    /*
     *  Handles leap seconds...
     *  This one should be leap second.
     */
    Date = 20120630;
    Time = 24.0 + 0.5/3600.0;
    Lgm_Make_UTC( Date, Time, &UTC, c );
    Lgm_DateTimeToString( Str, &UTC, 0, 3 );
    printf("Date: %ld Time: %lf  in ISO 8601 format is: %s\n", Date, Time, Str );


    /*
     *  This one should not be leap second.
     */
    Date = 20110630;
    Time = 24.0 + 0.5/3600.0;
    Lgm_Make_UTC( Date, Time, &UTC, c );
    Lgm_DateTimeToString( Str, &UTC, 0, 3 );
    printf("Date: %ld Time: %lf  in ISO 8601 format is: %s\n", Date, Time, Str );


    
    /*
     *  free structures
     */
    Lgm_free_ctrans( c ); 


    return(0);

}

