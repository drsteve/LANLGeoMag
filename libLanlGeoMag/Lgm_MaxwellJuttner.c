#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Lgm/Lgm_FluxToPsd.h"
#include <gsl/gsl_sf_bessel.h>

#define  kb = 8.61734315e-11  //  MeV/K


/*
 *  Maxwell-Boltzman distribution. (AKA Maxwellian).
 *  Maxwell-Juttner distribution. (AKA relativistic Maxwellian).
 */




/*
 * Maxwell-Juttner is (in a dumbed-down form);
 *
 *
 *                 n
 * f = -------------------------  exp( -E/kT ) (1)
 *      4pi kT m0^2 c K2(m0 c^2/kT)
 *
 *
 *  Define a = kT/(m0 c^2). 
 *
 *  Also, we have,
 *      E = sqrt( p^2c^2 + E0^2 )
 *
 *  Thus we can rewrite (1) as,
 *
 *
 *                c^3  n
 *     f = -------------------------  exp( -sqrt(p^2c^2 + E0^2)/(aE0) ) (2)
 *          4pi E0^3 a K2(1/a)
 *
 *
 *   or
 *
 *                c^3  n
 *     f = -------------------------  exp( -(Ek + E0)/(a E0) )    (3)
 *          4pi E0^3 a K2(1/a)
 *
 *
 *
 *  Note that a is dimensionless. If m0 c^2 = E0 is in MeV, then kT is in Mev.
 *  The boltsmann constant k is in MeV/K, so T is in Kelvin. But we want the
 *  temperature to be in MeV to start with, so we will replace kT by T. I.e.,
 *      
 *      a = T/E0
 *
 *  We will also want to end up with units of c^3 cm^-3 MeV^-3, so we will not explicity
 *  multiply the c^3.
 *
 *  Alternate form: Let b = 1/a (=E0/T). And since E = Ek + E0,
 *
 *                c^3  b n
 *     f = -------------------------  exp( -b E / E0 )             (4)
 *          4pi E0^3 K2(b)
 *
 *  For low temperature plasmas, K2(b) ~ sqrt( pi/(2 b) ) exp( -b )
 *  So for b >> 1, (4) becomes:
 *
 *
 *                c^3  b n exp( b )
 *     f = -------------------------  exp( -b E / E0 )             (5)
 *          4pi E0^3 sqrt( pi/(2 b) )
 *
 *  Putting b=E0/T and E0=mc^2 then gives;
 *
 *                n
 *     f = -------------------------  exp( -Ek / T )             (6)
 *          ( 2pi m T )^(3/2) 
 * 
 *  which is just the non-relativistic maxwellian. Or in slightly different form;
 *
 *                n c^3
 *     f = -------------------------  exp( -Ek / T )             (7)
 *          ( 2pi E0 T )^(3/2) 
 * 
 * 
 *
 *  A more rigorous form of Max-Jut is given by (Chacon-Acosta, G., L. Dagdug,
 *  H. A. Morales-Tecotl, "Manifestly covariant Juttner distribution and
 *  equipartition theorem", Phys Rev E, 81, 021126 [2010]);
 *
 *
 *                                                          (d-1)/2
 *                    n                         [ mc Theta ]
 *     f(p) = --------------------------------  [----------]      exp( -Theta_mu p^mu )
 *             2c(mc)^d K_{(d+1)/2}(mc Theta)   [  2pi     ]
 *
 *
 *
 *  where, d is number of spatial dimensions. Theta = c/kT. Theta_mu and p_mu
 *  are four-vectors.  For d=3, in the comoving frame, this is the same as (3).
 *  (In a co-moving frame, Theta_mu = (Theta, 0, 0, 0) and p_mu = (E/c, 0, 0,
 *  0).) 
 *
 *  See also: Lazar, M., A. Stockem, and R. Schlickeiser, "Towards a
 *  relativistically correct characterization of counterstreaming plasma. I.
 *  Distribution functions", The Open Plasma Physics Journal, 3, 138--147,
 *  2010. This paper described the defs and approximations for for "ultra-rel",
 *  "rel", "weak rel", "low temp", "cold".
 *
 *  Inputs:
 *          n  -- Density.        (#/cm^3)
 *          T  -- Temperature.    (keV)
 *          Ek -- Kinetic Energy. (MeV)
 *          E0 -- particle rest energy. (MeV)
 *
 *  Returns:
 *          PSD
 * 
 */
double  Lgm_MaxJut( double n, double T, double Ek, double E0 ) {

    double  p2c2, Theta, K2, E02, E03, f0, f;

    //p2c2 = Lgm_p2c2( Ek, E0 );

    E02 = E0*E0;    // MeV^2
    E03 = E02*E0;   // MeV^3

    T /= 1000.0; // convert T from keV -> MeV
    Theta = E0/T;  // dimensionless (MeV/MeV)

    if ( Theta > 10.0 ) {
        // Low temperature limit. K2 ~ sqrt( pi/(2*Theta) ) exp(-Theta)
        K2 = sqrt( 0.5*M_PI/Theta );
        f0 = n*Theta/( 4.0*M_PI * E03 * K2);    // units of c^3 cm^-3 MeV^-3
        //f = f0*exp( Theta*(1.0-(Ek+E0)/E0) );    // units of c^3 cm^-3 MeV^-3
        f = f0*exp( -Ek/T );    // units of c^3 cm^-3 MeV^-3
    } else {
        // compute Bessel function
        K2 = gsl_sf_bessel_Kn( 2, Theta ); // dimensionless
        f0 = n*Theta/( 4.0*M_PI * E03 * K2);    // units of c^3 cm^-3 MeV^-3
        f = f0*exp( -Theta*(Ek+E0)/E0 );    // units of c^3 cm^-3 MeV^-3
    }
        

    return( f ); 
    
}

/*
 * Derivatives wrt n and T
 */
void Lgm_MaxJut_Derivs( double n, double T, double Ek, double E0, double *dfdn, double *dfdT ) {
    double  K2, E02, E03, f0, f;
    double  sq, ex, b, g, h, q, gnq, dK2db, dfdb, dbdT;
    

    E02 = E0*E0;    // MeV^2
    E03 = E02*E0;   // MeV^3
    b = 1000.0*E0/T;  // dimensionless (MeV/MeV)

    if ( b < 10.0 ) {
        K2 = gsl_sf_bessel_Kn( 2, b ); // dimensionless
        dK2db = -0.5*( gsl_sf_bessel_Kn( 1, b ) + gsl_sf_bessel_Kn( 3, b ) );
        //printf("b < 100: K2, dK2db = %g %g\n", K2, dK2db);
    } else {
        sq = sqrt( 0.5*M_PI/b );
        ex = exp( -b );
        K2 = sq * ex;
        dK2db = -ex * sq * ( 1.0 + 0.5/b );
        //printf("b >= 100: b = %g K2, dK2db = %g %g\n", b, K2, dK2db);
    }

    g  = 1.0/( 4.0*M_PI * E03 );    // units of c^3 cm^-3 MeV^-3
    h  = -(Ek+E0)/E0;
    q  = exp( b*h );
    f  = n*g*b/K2 * q;    // units of c^3 cm^-3 MeV^-3
    //printf("g, h, q, f = %g %g %g %g\n", g, h, q, f);
    
    // df/dn is trivial
    *dfdn = f/n;

    // df/dT is more involved
    gnq = g*n*q;
    dfdb = gnq/K2 + gnq*h*b/K2 - gnq*b/(K2*K2) * dK2db;
    dbdT = -1000.0*E0/(T*T);
    //printf("gnq, dfdb, dbdT = %g %g %g\n", gnq, dfdb, dbdT);
    *dfdT = dfdb*dbdT;
    
}


/*
 * Non-relativistic maxwellian;
 * 
 * 
 *                 n
 *     f = -------------------------  exp( -Ek/(kT) )
 *          (2pi m k T )^(3/2)
 *
 * Note that m = E0/c^2 (E0 is rest energy of particle)
 *
 *                 n c^3
 *     f = -------------------------  exp( -Ek/(kT) )
 *          (2pi E0 k T )^(3/2)
 *
 * And kT -> T if T is in units of MeV already.
 *
 *  Inputs:
 *          n  -- Density.        (#/cm^3)
 *          T  -- Temperature.    (keV)
 *          Ek -- Kinetic Energy. (MeV)
 *          E0 -- particle rest energy. (MeV)
 *
 *  Returns:
 *          PSD
 */
double  Lgm_Maxwellian( double n, double T, double Ek, double E0 ) {

    double f;

    T /= 1000.0; // keV -> MeV
    
    f = n*pow( 1.0/(2.0*M_PI*E0*T), 1.5) * exp( -Ek/T ); // units of c^3 / (cm^3 MeV^3)

    return(f);

}


double  Lgm_Maxwellian_dfdn( double n, double T, double Ek, double E0 ) {

    double dfdn;

    T /= 1000.0; // keV -> MeV
    
    dfdn = pow( 1.0/(2.0*M_PI*E0*T), 1.5) * exp( -Ek/T ); // units of c^3 / (cm^3 MeV^3)

    return(dfdn);
    
}


double  Lgm_Maxwellian_dfdT( double n, double T, double Ek, double E0 ) {

    double dfdT;

    T /= 1000.0; // keV -> MeV
    
    dfdT = n*pow( 1.0/(2.0*M_PI*E0*T), 1.5) * exp( -Ek/T ) * (Ek/(T*T) - 1.5/T); // units of c^3 / (cm^3 MeV^3)

    return(dfdT);
    
}







