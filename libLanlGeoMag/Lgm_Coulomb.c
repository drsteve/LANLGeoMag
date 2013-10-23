#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#if USE_OPENMP
#include <omp.h>
#endif
#include "Lgm/Lgm_Coulomb.h"
#include "Lgm/Lgm_Constants.h"




/**
 *  \brief
 *     Computes the differential cross section for relativistic electron
 *     scattering from an unshielded nucleus or ion (Mott and Massey, 1965).
 *  \details
 *     The differential cross section for relativistic electron scattering from
 *     an unshielded nucleus or ion with atomic number \f$Z_i\f$ is given by
 *     (see Equation (47) in Chapter IX of Mott and Massey, "The Theory of
 *     Atomic Collisions", Third Edition, 1965.);
 *
 *      \f[{d\sigma\over d\Omega} = { Z_i^2 r_e^2\over 4\gamma^2\beta^4\sin^4(\theta / 2) } \left( 1 - \beta^2\sin^2(\theta / 2) \right)\f]
 *
 *      where \f$\theta\f$ is the scattering angle, \f$r_e\f$ is the classical
 *      electron radius (\f$r_e=e^2/(4\pi\epsilon_\circ m_ec^2)\f$),
 *      \f$\beta=v/c\f$, \f$\gamma=(1-\beta^2)^{-1/2}\f$. This form of the Mott
 *      differential scattering cross section is an approximation valid for
 *      light elements such that \f$\alpha=Z/137\f$ is small compared with
 *      unity.
 *
 *
 *
 *      \param[in]      Z       Atomic mass of the scattering nucleus. Z/137 must be << 1.
 *      \param[in]      beta    v/c.       
 *      \param[in]      Theta   scattering angle.  <b>( radians )</b>
 *
 *      \return         Differential scattering cross section.      <b>( m )</b>
 *
 *      \author         Mike Henderson
 *      \date           2013
 */
double  Lgm_MottScattering( double Z, double beta, double Theta ) {

    double  beta2, beta4, gamma2;
    double  SinThetaOver2, SinThetaOver2_2, SinThetaOver2_4;
    double  d_sigma_d_Omega;

    beta2  = beta*beta;         // beta^2
    beta4  = beta2*beta2;       // beta^4
    gamma2 = 1.0/(1.0-beta2);   // gamma^2
    
    SinThetaOver2   = sin( Theta/2.0 );
    SinThetaOver2_2 = SinThetaOver2*SinThetaOver2;
    SinThetaOver2_4 = SinThetaOver2_2*SinThetaOver2_2;
    d_sigma_d_Omega = Z*Z*LGM_ELECTRON_RADIUS*LGM_ELECTRON_RADIUS/(4.0*gamma2*beta4*SinThetaOver2_4) * ( 1.0 - beta2*SinThetaOver2_2 );
    

    return( d_sigma_d_Omega );

}




//void Lgm_MottInit( Lgm_MottInfo *m ) {
//
//    // select random number generator
//    m->RngType = gsl_rng_default;
//    
//    // create rng
//    m->Rng = gsl_rng_alloc( m->RngType );
//
//
//    return;
//
//
//}
//
//
//void Lgm_MottFree( Lgm_MottInfo *m ) {
//
//    // free random number generator
//    gsl_rng_free( m->rng );
//
//
//    return;
//
//}
//
//
//
//
//
//void Lgm_MottGetRandomlySampledAngles( double *Theta, double *Phi ) {
//
//    double u1, u2;
//
//
//    /*
//     *  Get two uniform random variates in [0.0, 1.0)
//     */
//    u1 = gsl_rng_uniform( m->Rng );
//    u2 = gsl_rng_uniform( m->Rng );
//
//
//    *Theta = FunctionOf( u1 );
//    *Phi   = DifferentFunctionOf( u2 );
//    
//    return;
//
//}
//





