#include <Lgm_CTrans.h>
#include <Lgm_QinDenton.h>
int main( ) {
    Lgm_CTrans       *c = Lgm_init_ctrans( 1 ); // more compact declaration
    Lgm_QinDentonOne p;
    long int         Date;
    double           UTC, JD;
    
    Date = 20020713;                        // August 12, 2004
    UTC  = 18.0 + 0.0/60.0 + 30.0/3600.0;   // Universal Time Coordinated (in decimal hours)
    JD = Lgm_Date_to_JD( Date, UTC, c );    // Compute JD

    // Get (interpolate) the QinDenton vals from the values in the file at the given Julian Date
    Lgm_get_QinDenton_at_JD( JD, &p, 1 );

    UTC  = 18.0 + 1.0/60.0 + 30.0/3600.0;   // Universal Time Coordinated (in decimal hours)
    JD = Lgm_Date_to_JD( Date, UTC, c );    // Compute JD
    Lgm_get_QinDenton_at_JD( JD, &p, 1 );

    Lgm_free_ctrans( c ); // free the structure
    return(0);

}

