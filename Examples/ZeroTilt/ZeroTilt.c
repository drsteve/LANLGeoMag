/*
 * Iteratively find time of zero tilt on April 1, 2017. Useful
 * for setting times when you need zero tilt.
 */
#include <Lgm_CTrans.h>
#include <Lgm_Vec.h>

int main( ) {

    Lgm_CTrans  *c = Lgm_init_ctrans( 0 ); 
    int         done;
    long int    Date;
    double      UTa, UTb, UTc;
    double      Psia, Psib, Psic;
    
    Date = 20170401;    // Apr 1, 2017

    // set initial value
    UTa =  0.0; Lgm_Set_Coord_Transforms( Date, UTa, c ); Psia = c->psi;

    // find time when Psi has opposite sign
    UTc = UTa; done = 0;
    while ( !done ) {
        UTc += 1.0; Lgm_Set_Coord_Transforms( Date, UTc, c ); Psic = c->psi;
        if ( Psic*Psia <= 0.0 ) done = 1;
    }

    printf("Initial bracket:\n");
    printf("UTa, PSia = %g %g\n", UTa, Psia );
    printf("UTc, PSic = %g %g\n\n", UTc, Psic );
    
    
    

    done = 0;
    while ( !done ) {

        /*
         * compute Psi at midpoint time.
         */
        UTb = 0.5*(UTa+UTc);
        Lgm_Set_Coord_Transforms( Date, UTb, c ); Psib = c->psi;

        if ( Psib*Psia <= 0.0 ) {
            Psic = Psib;
            UTc  = UTb;
        } else {
            Psia = Psib;
            UTa  = UTb;
        }

        if ( fabs(UTc-UTa) < 1e-8 ) done = 1;

    }
    printf("Final result for Date=%8ld :\n", Date );
    printf("UTa, PSia = %.10lf %.10lf\n", UTa, Psia );
    printf("UTb, PSib = %.10lf %.10lf\n", UTb, Psib );
    printf("UTc, PSic = %.10lf %.10lf\n", UTc, Psic );
    


    Lgm_free_ctrans( c ); // free the structure

    return(0);
}
