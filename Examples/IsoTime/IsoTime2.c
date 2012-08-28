#include <Lgm_CTrans.h>

/*
 * Shows how to get time from ISO format using IsoTimeStringToDateTime()
 */
int main( ) {

    Lgm_CTrans      *c = Lgm_init_ctrans( 0 );
    Lgm_DateTime    UTC;
    char            Str[80];



    /*
     * Shows how to get ISO string to DateTime Structure.
     */
    sprintf( Str, "20110630T23:45:51.45678Z" );
    IsoTimeStringToDateTime( Str, &UTC, c );
    printf("ISO 8601 string: %s\n", Str );
    printf("Date: %ld\n", UTC.Date );
    printf("Time: %lf\n\n\n", UTC.Time );

    /*
     * Handles leap seconds
     */
    sprintf( Str, "20120630T23:59:60.5Z" );
    IsoTimeStringToDateTime( Str, &UTC, c );
    printf("ISO 8601 string: %s\n", Str );
    printf("Date: %ld\n", UTC.Date );
    printf("Time: %lf\n\n\n", UTC.Time );


    sprintf( Str, "20110630T23:59:60.5Z" );
    IsoTimeStringToDateTime( Str, &UTC, c );
    printf("ISO 8601 string: %s\n", Str );
    printf("Date: %ld\n", UTC.Date );
    printf("Time: %lf\n\n\n", UTC.Time );

    
    /*
     *  free structures
     */
    Lgm_free_ctrans( c ); 


    return(0);

}

