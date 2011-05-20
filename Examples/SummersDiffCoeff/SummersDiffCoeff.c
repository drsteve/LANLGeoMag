#include <stdio.h>
#include <Lgm_SummersDiffCoeff.h>
int main( ) {

    double  Daa_ba, Dap_ba, Dpp_ba;

    getDs_ba( 42.0, 1.0, 4.0, 0.0, &Daa_ba, &Dap_ba, &Dpp_ba );

    printf("Daa_ba, Dap_ba, Dpp_ba = %g %g %g\n", Daa_ba, Dap_ba, Dpp_ba);


    return(0);
}

