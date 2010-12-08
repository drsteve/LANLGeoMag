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
 * returns (v/c)^2
 *
 *  (v/c)^2 = m^2v^2/(m^2c^2) 
 *          = m^2v^2c^2/(m^2c^4) 
 *          = p^2c^2/(m^2c^4) 
 *          = p^2c^2/E^2
 *          = Ek(Ek+2E0)/(Ek+E0)^2
 */
double  Lgm_v2overc2( double Ek, double E0 ) {
    double  E = Ek + E0;
    return( Ek*(Ek+2.0*E0)/(E*E) );    // dimensionless
}


/*
 * returns relativistic factor gamma = [ 1 - (v/c)^2 ]^(-1/2)
 */
double  Lgm_gamma( double Ek, double E0 ) {
    double  E = Ek + E0;
    return( 1.0/sqrt( 1.0 - Ek*(Ek+2.0*E0)/(E*E) ) );    // dimensionless
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





Lgm_PhaseSpaceDensity *Lgm_InitPhaseSpaceDensity( double Flux, double E, double A, int nE, int nA, int DumpDiagnostics ) {

    
    int     i, j, nE, na;
    double  flux, p2c2, fp, Min, Max;
    double  **PsdArray_LoRes, **PsdArray_HiRes;

    Lgm_PhaseSpaceDensity *p;
    

    p = (Lgm_PhaseSpaceDensity *) calloc( 1, sizeof(*p) );
    p->DumpDiagnostics = DumpDiagnostics;


    /*
     * Add Flux array info to p
     */
    p->nE1 = nE; 
    p->nA1 = nA; 
    ARRAY_1D( p->E1, p->nE, double );
    ARRAY_1D( p->A1, p->nA, double );
    ARRAY_2D( p->FLUX_EA1, p->nE, p->nA, double );
    for (i=0; i<p->nE; i++) p->E1[i] = E[i];
    for (i=0; i<p->nA; i++) p->A1[i] = A[i];
    for (j=0; j<p->nA; j++) {
        for (i=0; i<p->nE; i++) {
            p->FLUX_EA1[i][j] = Flux[i][j];
        }
    }


    /*
     * Convert Flux array info to PSD
     */
    ARRAY_2D( p->PSD_EA1, p->nE, p->nA, double );
    for (j=0; j<p->nA; j++) {
        for (i=0; i<p->nE; i++) {

            flux   = p->Flux[i][j];
            p2c2   = Lgm_p2c2( p->E1[i], LGM_Ee0 );
            fp     = Lgm_DiffFluxToPsd( flux, p2c2 );
            fp     = (fp > 0.0) ? log10(fp) : LGM_FILL_VALUE;
            p->PSD_EA1[p->nE1-i][j] = fp;
            
        }
    }
    if ( p->DumpDiagnostics ) {
        DumpGif( "PSD_Versus_E_and_A_LoRes.gif", p->nE1, p->nA1, p->PSD_EA1 );
    }




    /*
     * Set desired size of HiRes array.
     */
    p->nE2 = 400; // Energy 
    p->nA2 = 400; // Pitch Angles
    ARRAY_2D( p->PSD_EA2, p->nE2, p->nA2, double );         // (Energy, Pitch Angle) -> (Row/Col)
    UpSizeImage( p->PSD_EA1, p->E1, p->A1, p->nE1, p->nA1, 
                    p->PSD_EA2, p->E2, p->A2, p->nE2, p->nA2 ); // returns p->PSD_EA2, p->E2, p->A2
    if ( p->DumpDiagnostics ) {
        DumpGif( "PSD_Versus_E_and_A_HiRes.gif", p->nE2, p->nA2, p->PSD_EA2 );
    }



    return p;

}




/*
 *  The routine Lgm_ComputePsdVersusEAndAlpha() gives us a Hi-Res array of f(E,
 *  alpha).  BUT, what we really need in the end is f( mu, K ). Although mu is
 *  easy to compute, it is dependant on both E and alpha. K is only dependant
 *  upon alpha, but on the other hand K is not so easy to compute and we dont
 *  want to have to compute lots of K's if we dont have to. So we will use the
 *  following strategy instead,
 * 
 *  Note that f( E, a ) is the same as f( E(mu, a(K)), a(K) ). Thus, for a
 *  given mu and K, we can figure out what E and a they correspond to and then
 *  we can just look up the f value in our HiRes array. The steps are;
 *
 *      1. For each K, compute a(K). We already have this routine ( AlphaOfK() ).
 *
 *      2. Then we compute E from a and the given mu value.
 *
 *      3. Then just look up f(E,a) from the array (interp or fit or whatever).
 *
 *  Inputs:
 *          nMu -- Number of Mu values
 *           Mu -- 1-D array of Mu values
 *           nK -- Number of K values
 *            K -- 1-D array of K values
 *  Outputs:
 *          PSD -- 2-D array of f(Mu, K). The user must alloc memory for this array.
 * 
 *  Usage:
 *          If m->UseInterpRoutines is TRUE, then the user must have pre-traced
 *          a FL with Lgm_TraceLine(). This initializes field line dependant
 *          information that is needed for Lgm_AlphaOfK() to work in the
 *          interpolated mode. Basically Lgm_AlphaOfK() uses bisection to solve
 *          for a(K) and it much faster to use the pre-tracing strategy. E.g.;
 *
 *
 *              m->UseInterpRoutines = TRUE;
 *              Lgm_TraceLine( &u_in, &u_out, m->Lgm_LossConeHeight, -1.0, 1e-8, FALSE, m );
 *              LGM_ARRAY_2D( PSD, nMU, nK, double );
 *              Lgm_PsdAtConstMuAndK( PSD, nMu, MU, nK, K, m );
 *                  ... do stuff with PSD ...
 *              LGM_ARRAY_2D_FREE( PSD );
 * 
 * 
 */
void Lgm_PsdAtConstMuAndK( double **PSD, int nMu, double *Mu, int nK, double *K, Lgm_MagModelInfo *m ) {

    int     i, k;
    double  *a, E;

    
    /*
     * Compute the alpha's -- 1d array
     * parallelize this loop?
     */
    // #pragma etc..
    LGM_ARRAY_1D( a, nK, double );
    for ( k=0; k<nK; k++ ){
        a[k] = Lgm_AlphaOfK( K[k], m );
    }


    /*
     * Compute the PSD's -- 2d array
     */
    for ( i=0; i<nMu; i++ ){
        for ( k=0; k<nK; k++ ){
            E = Lgm_Mu_to_Energy( Mu[i], a[k], B );
            PSD[i][k] = Lgm_Psd( E, a, INFO THAT CONTAINS THE HIRES array );
        }
    }

    LGM_ARRAY_1D_FREE( a );

}






/*
 *  Reverse direction....
 *
 *  Note that f( mu, K ) is the same as f( mu( E, a), K(a) ). Thus, for a given
 *  E and a, we can figure out what mu and K they correspond to and then we
 *  can just look up the f value in our HiRes array. The steps would be;
 *
 *      1. For each a, compute K(a). We already have this routine somewhere.
 *
 *      2. Also compute mu from the given E and a.
 *
 *      3. Then just look f( mu, K ) up from the array (interp or fit or whatever).
 * 
 */












