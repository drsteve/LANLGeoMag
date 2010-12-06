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


/*
 * Compute p^2c^2 given a particle's kinetic energy and rest energy.
 *
 *  Some relativistic equations:
 *  
 *  With,
 *      m     = Particle mass.
 *      m0    = Particle rest mass.
 *      v     = Particle speed.
 *      c     = Speed of Light.
 *      gamma = (1-v^2/c^2)^(-1/2).
 *
 *      m = m0 gamma = m0 (1-v^2/c^2)^(-1/2)
 *
 *  Rearrange this to get,
 *      (mc^2)^2 = (m0c^2)^ + p^2c^2
 *
 *  With,
 *      E  = mc^2  = Total Energy of particle
 *      E0 = m0c^2 = Rest Energy of particle
 *      p  = relativistic momentum of particle
 *
 *  this is,
 *      E^2 = E0^2 + (pc)^2
 *  
 *  so, 
 *      (pc)^2  = E^2 - E0^2
 *  
 *  Let E = Ek+E0 (kinetic energy + rest energy). Then,
 *      p^2c^2  = (Ek+E0)^2 - E0^2
 *              = Ek (Ek+2E0)
 *  
 *
 *  Inputs:
 *          E  -- Kinetic Energy of Particle.    ( MeV )
 *          E0 -- Rest energy of Particle.       ( MeV )
 *
 *  Returns:
 *          p^2c^2 = Ek(Ek+2E0) ( MeV^2 )
 */
double  Lgm_p2c2( double Ek, double E0 ) {
    return( Ek*(Ek+2.0*E0) );    // p^2c^2 in units of MeV^2
}






/*
 * Convert differential flux to phase space density.
 *
 *  The basic relationship is;
 *      f = j/p^2
 *
 *  Multiply top and bottom by c^2 gives,
 *      f = (j c^2)/(p^2c^2)
 *      f = j/c * c^3/(p^2c^2)
 *  
 *  Reason for making it c^3 is that the final units become,
 *  c^3 cm^-3 MeV^-3 or (c/cm/MeV)^3
 *  
 *  Inputs:
 *          j    -- Differential Flux in units of #/cm^2/s/sr/MeV
 *          p2c2 -- (pc)^2 in units of Mev^2
 *  
 *  Output:
 *          Phase space density in units of (c/cm/MeV)^3
 *  
 */
double Lgm_DiffFluxToPsd( double j, double p2c2 ){
    return( j/(p2c2*2.9979e10) ); // f in units of (c/cm/MeV)^3
}


/*
 * Convert phase space density to differential flux.
 *
 *  The basic relationship is;
 *      f = j/p^2
 *
 *  Multiply top and bottom by c^2 gives,
 *      f = (j c^2)/(p^2c^2)
 *      f = j/c * c^3/(p^2c^2)
 *
 *  Reason for making it c^3 is that the final units become,
 *  c^3 cm^-3 MeV^-3 or (c/cm/MeV)^3
 *
 *  Inputs:
 *          f    -- Phase space density in units of (c/cm/MeV)^3
 *          p2c2 -- (pc)^2 in units of Mev^2
 *
 *  Output:
 *          Differential Flux in units of #/cm^2/s/sr/MeV
 *
 */
double Lgm_PsdToDiffFlux( double f, double p2c2 ){
    return( f*2.9979e10/p2c2 ); // j in units of #/cm^2/s/sr/MeV
}
































void Lgm_ComputePsdVersusMuAndAlpha(     DATA STRUCT   ) {


    int nE = 400; 
    int na = 400;



    /*
     * Create: 
     *          1. PSD( E, Alpha ) at native resolution 
     *          2. PSD( E, Alpha ) at high resolution
     *
     * The input data structure has an array containing j(E, Alpha). I.e.
     * differential flux versus energy and pitch angle. This will normally be a
     * fairly small array (e.g. 7 energy channels by 9 pitch angles). Here we
     * want to interpolate this array up to a much larger smooth array.
     */
    ARRAY_2D( PsdArray_Native, d->nE, d->na, 9, double );   // (Energy, Pitch Angle) -> (Row/Col)
    ARRAY_2D( PsdArray_HiRes, nE, na, double );             // (Energy, Pitch Angle) -> (Row/Col)
    for ( j=0; j<d->na; j++ ){
        for ( i=0; i<d->nE; i++ ){
            flux   = d->Flux[i][j];
            fp     = j_to_fp_1( flux, Data->Energy[i] );
            fp = (fp > 0.0) ? log10(fp) : -9e99;
            Image[6-i][j] = fp;
            if ((fp > Max)&&(flux > -1e99))  Max = fp;
            if ((fp < Min)&&(flux > -1e99))  Min = fp;
        }
    }
    DumpGif( "PSD_Versus_E_and_Alpha.gif", 9, 7, Image );
    UpSizeImage( Image, 7, 9, NewImage, 400, 400 );
    ARRAY_2D_FREE( Image );
    DumpGif( "PSD_Versus_E_and_Alpha_NEW.gif", 400, 400, NewImage );




















}














