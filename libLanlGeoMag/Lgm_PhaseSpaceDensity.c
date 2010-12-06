#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Lgm/Lgm_PhaseSpaceDensity.h"
#include "Lgm/Lgm_CTrans.h"


/*
 *  Routines for converting Flux <-> PSD.
 */




/*
 * Returns First adiabatic invariant, Mu, given: Particle's kinetic energy,
 * Pitch Angle and the local B-field strength.
 *  
 *  Inputs:
 *          E -- Kinetic Energy of Particle.    ( MeV )
 *          a -- Pitch Angle of Particle.       ( Degrees )
 *          B -- Local magnetic field strength. ( nT )
 *
 *  Returns:
 *          First adiabatic invariant, Mu. ( MeV/nT )
 */
double  Lgm_Energy_to_Mu( double E, double a, double B ) {

    double  sa, sa2;

    if ( B <= 0.0 ) return( -9e99 );

    sa = sin( a*RadPerDeg );    // sin(Alpha)
    sa2 = sa*sa;                // sin^2(Alpha)
    
    return( E*sa2/B );          // Mu = E*sin^2(Alpha)/B
    
}



/*
 * Returns Particle's Kinetic Energy, E, given: Particle's first invariant, Mu,
 * Pitch Angle and the local B-field strength.
 *  
 *  Inputs:
 *          Mu -- Kinetic Energy of Particle.    ( MeV/nT )
 *          a  -- Pitch Angle of Particle.       ( Degrees )
 *          B  -- Local magnetic field strength. ( nT )
 *
 *  Returns:
 *          First adiabatic invariant, Mu. ( MeV )
 */
double  Lgm_Mu_to_Energy( double Mu, double a, double B ) {

    double  sa, sa2;

    
    sa = sin( a*RadPerDeg );    // sin(Alpha)
    sa2 = sa*sa;                // sin^2(Alpha)

    if ( sa2 < 1e-10 ) return( -9e99 );
    else return( Mu*B/sa2 );    // Mu = E*sin^2(Alpha)/B


}






