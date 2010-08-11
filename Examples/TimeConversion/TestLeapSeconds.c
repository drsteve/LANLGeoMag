#include <Lgm_CTrans.h>
int main( ) {
    long int        sec, JDN, Date;
    int             s, m, smax, i, Pass, Year, Month, Day;
    Lgm_CTrans      *c = Lgm_init_ctrans( 0 );
    Lgm_DateTime    UTC, TAI;
    double          u0, u1, t, Time, Seconds;

    /*
     * Confirm that the all days betwen 1972 to present are categorized correctly as leap second
     * or non-leapsecond days.
     */
    printf("Test to see if all JDN between 2441316 and 2454834 are\n"
           "categorized correctly as leap second or non-leapsecond days: ");
    Pass = TRUE;
    for ( JDN = 2441316; JDN <= 2454834; JDN += 1 ){
        Lgm_JD_to_Date( (double)JDN, &Year, &Month, &Day, &Time );
        Date = Year*10000 + Month*100 + Day;

        if      ( (Date == 19711231) && ( Lgm_IsLeapSecondDay( Date, &Seconds, c ) != 1 ) ) { Pass = FALSE; }
        else if ( (Date == 19720630) && ( Lgm_IsLeapSecondDay( Date, &Seconds, c ) != 1 ) ) { Pass = FALSE; }
        else if ( (Date == 19721231) && ( Lgm_IsLeapSecondDay( Date, &Seconds, c ) != 1 ) ) { Pass = FALSE; }
        else if ( (Date == 19731231) && ( Lgm_IsLeapSecondDay( Date, &Seconds, c ) != 1 ) ) { Pass = FALSE; }
        else if ( (Date == 19741231) && ( Lgm_IsLeapSecondDay( Date, &Seconds, c ) != 1 ) ) { Pass = FALSE; }
        else if ( (Date == 19751231) && ( Lgm_IsLeapSecondDay( Date, &Seconds, c ) != 1 ) ) { Pass = FALSE; }
        else if ( (Date == 19761231) && ( Lgm_IsLeapSecondDay( Date, &Seconds, c ) != 1 ) ) { Pass = FALSE; }
        else if ( (Date == 19771231) && ( Lgm_IsLeapSecondDay( Date, &Seconds, c ) != 1 ) ) { Pass = FALSE; }
        else if ( (Date == 19781231) && ( Lgm_IsLeapSecondDay( Date, &Seconds, c ) != 1 ) ) { Pass = FALSE; }
        else if ( (Date == 19791231) && ( Lgm_IsLeapSecondDay( Date, &Seconds, c ) != 1 ) ) { Pass = FALSE; }
        else if ( (Date == 19810630) && ( Lgm_IsLeapSecondDay( Date, &Seconds, c ) != 1 ) ) { Pass = FALSE; }
        else if ( (Date == 19820630) && ( Lgm_IsLeapSecondDay( Date, &Seconds, c ) != 1 ) ) { Pass = FALSE; }
        else if ( (Date == 19830630) && ( Lgm_IsLeapSecondDay( Date, &Seconds, c ) != 1 ) ) { Pass = FALSE; }
        else if ( (Date == 19850630) && ( Lgm_IsLeapSecondDay( Date, &Seconds, c ) != 1 ) ) { Pass = FALSE; }
        else if ( (Date == 19871231) && ( Lgm_IsLeapSecondDay( Date, &Seconds, c ) != 1 ) ) { Pass = FALSE; }
        else if ( (Date == 19891231) && ( Lgm_IsLeapSecondDay( Date, &Seconds, c ) != 1 ) ) { Pass = FALSE; }
        else if ( (Date == 19901231) && ( Lgm_IsLeapSecondDay( Date, &Seconds, c ) != 1 ) ) { Pass = FALSE; }
        else if ( (Date == 19920630) && ( Lgm_IsLeapSecondDay( Date, &Seconds, c ) != 1 ) ) { Pass = FALSE; }
        else if ( (Date == 19930630) && ( Lgm_IsLeapSecondDay( Date, &Seconds, c ) != 1 ) ) { Pass = FALSE; }
        else if ( (Date == 19940630) && ( Lgm_IsLeapSecondDay( Date, &Seconds, c ) != 1 ) ) { Pass = FALSE; }
        else if ( (Date == 19951231) && ( Lgm_IsLeapSecondDay( Date, &Seconds, c ) != 1 ) ) { Pass = FALSE; }
        else if ( (Date == 19970630) && ( Lgm_IsLeapSecondDay( Date, &Seconds, c ) != 1 ) ) { Pass = FALSE; }
        else if ( (Date == 19981231) && ( Lgm_IsLeapSecondDay( Date, &Seconds, c ) != 1 ) ) { Pass = FALSE; }
        else if ( (Date == 20051231) && ( Lgm_IsLeapSecondDay( Date, &Seconds, c ) != 1 ) ) { Pass = FALSE; }
        else if ( (Date == 20081231) && ( Lgm_IsLeapSecondDay( Date, &Seconds, c ) != 1 ) ) { Pass = FALSE; }
        else if ( (Date != 19711231) && (Date != 19720630) && (Date != 19721231) && (Date != 19731231)
                && (Date != 19741231) && (Date != 19751231) && (Date != 19761231) && (Date != 19771231)
                && (Date != 19781231) && (Date != 19791231) && (Date != 19810630) && (Date != 19820630)
                && (Date != 19830630) && (Date != 19850630) && (Date != 19871231) && (Date != 19891231)
                && (Date != 19901231) && (Date != 19920630) && (Date != 19930630) && (Date != 19940630)
                && (Date != 19951231) && (Date != 19970630) && (Date != 19981231) && (Date != 20051231)
                && (Date != 20081231) && ( Lgm_IsLeapSecondDay( Date, &Seconds, c ) != 0 ) ) { Pass = FALSE; }

    }

    if ( Pass ) {
        printf("Passed\n");
    } else {
        printf("Failed\n");
    }
    
    // free structures
    Lgm_free_ctrans( c ); 

    return(0);

}

