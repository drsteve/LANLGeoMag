#include <stdio.h>
#include <math.h>
#include <Lgm_MagModelInfo.h>
#include <Lgm_FluxToPsd.h>
int main( ) {

    double  Ek, a, B, Mu;

    Ek = 1.0;   // MeV
    a  = 45.0;  // Degrees
    B  = 100.0;  // nT

    printf("Testing Fwd/Bwd calculations...\n");
    Mu = Lgm_Ek_to_Mu( Ek, a, B, LGM_Ee0 );
    printf("Mu( Ek=%g MeV, a=%g Deg., B=%g nT ) = %g MeV/G\n", Ek, a, B, Mu );

    Ek = Lgm_Mu_to_Ek( Mu, a, B, LGM_Ee0 );
    printf("Ek( Mu=%g MeV/G, a=%g Deg., B=%g nT ) = %g MeV\n", Mu, a, B, Ek );


    printf("\n\n");
    Mu = 700.0; // MeV/G
    B  = 315.0; // nT (B at about L of 4.5 in a dipole)
    a  = 52.0;  // Deg.
    Ek = Lgm_Mu_to_Ek( Mu, a, B, LGM_Ee0 );
    printf("Ek( Mu=%g MeV/G, a=%g Deg., B=%g nT ) = %g MeV\n", Mu, a, B, Ek );


    printf("\n\n");
    Ek = 1.5; // MeV
    B  = 100.0; // nT (B at about L of 4.5 in a dipole)
    a  = 90.0;  // Deg.
    Mu = Lgm_Ek_to_Mu( Ek, a, B, LGM_Ee0 );
    printf("Mu( Ek=%g MeV, a=%g Deg., B=%g nT ) = %g MeV/G\n", Ek, a, B, Mu );

    return(0);
}

