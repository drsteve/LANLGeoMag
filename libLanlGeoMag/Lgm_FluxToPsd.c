#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_multifit.h>
#if USE_OPENMP
#include <omp.h>
#endif
#include "Lgm/Lgm_FluxToPsd.h"
#include "Lgm/Lgm_CTrans.h"
#include "Lgm/Lgm_MagModelInfo.h"
#include "Lgm/Lgm_DynamicMemory.h"
#include "Lgm/GPR.h"
//int FLAGFLAG=0;


void praxis( int n, double *x, int *data, double (*funct)(double *, void *data), double *in, double *out);
int WriteGIF( FILE *fp, unsigned char *pic, int ptype, int w, int h, unsigned char *rmap, unsigned char *gmap, unsigned char *bmap, int numcols, int colorstyle, char *comment);
void DumpGif2( char *FilenameBase, double Min, double Max, int W, int H, double **Image );





typedef struct _FitData {

    int     nMaxwellians;
    int     n;
    double  *E;
    double  *g;
    double  *dg;


} _FitData;



// Routines for converting between Flux and Phase Space Density





/**
 *  \brief
 *     Computes relativisitic first adiabatic invariant.
 *  \details
 *     Computes relativisitic first adiabatic invariant, \f$\mu\f$, given: Particle's
 *     kinetic energy \f$E_k\f$, Pitch Angle \f$\alpha\f$, the local B-field strength \f$B\f$, and the
 *     particle's rest energy \f$E_\circ\f$. The relationship is:
 *
 *      \f[\mu = {p_\perp^2\over 2 m_\circ B}\f].
 *
 *      Since,
 *
 *            \f[ p^2c^2 = E_k^2 + 2 E_k E_\circ \f]
 *
 *      we have,
 *
 *          \f{eqnarray*}
 *              \mu = {p_\perp^2 \over 2 m_\circ B } &=& {p^2c^2\over (2 m_\circ c^2 B)} \sin^2(\alpha) \\
 *                                                   &=& {p^2c^2\over (2 E_\circ B)} \sin^2(\alpha) \\
 *                                                   &=& {E_k\over B}\left[1+{E_k \over 2 E_\circ}\right] \sin^2(\alpha)
 *          \f}
 *
 *      \param[in]      Ek  Kinetic Energy of Particle.     <b>( MeV )</b>
 *      \param[in]      a   Pitch Angle of Particle.        <b>( Degrees )</b>
 *      \param[in]      B   Local magnetic field strength.  <b>( nT )</b>
 *      \param[in]      E0  Rest mass of Particle.          <b>( MeV )</b>
 *
 *      \return         First adiabatic invariant, Mu.      <b>( MeV/G )</b>
 *
 *      \author         Mike Henderson
 *      \date           2010-2011
 */
double  Lgm_Ek_to_Mu( double Ek, double a, double B, double E0 ) {

    double  sa, sa2, mu;

    if ( B <= 1e-12 ) return( 9e99 );

    sa = sin( a*RadPerDeg );    // sin(Alpha)
    sa2 = sa*sa;                // sin^2(Alpha)

    mu = Ek*(1.0+0.5*Ek/E0)*sa2*nT_Per_Gauss/B;

    return( mu );

}



/**
 *  \brief
 *     Computes kinetic energy from Mu, pitch angle, B-field strength, and rest energy.
 *  \details
 *     This Returns Particle's Kinetic Energy, \f$E_k\f$, given: Particle's relativistic
 *     first invariant, \f$\alpha\f$, Pitch Angle \f$\alpha\f$ and the local
 *     B-field strength and rest energy. This is the inverse of Lgm_Ek_to_Mu(). See
 *     description of Lgm_Ek_to_Mu() for more details on the equations used.
 *
 *      \param[in]      Mu  First adiabatic invariant.      <b>( MeV/G )</b>
 *      \param[in]      a   Pitch Angle of Particle.        <b>( Degrees )</b>
 *      \param[in]      B   Local magnetic field strength.  <b>( nT )</b>
 *      \param[in]      E0  Rest mass of Particle.          <b>( MeV )</b>
 *
 *      \return         First adiabatic invariant, Mu.      <b>( MeV )</b>
 *
 *      \author         Mike Henderson
 *      \date           2010-2011
 */
double  Lgm_Mu_to_Ek( double Mu, double a, double B, double E0 ) {

    double  sa, sa2, Ek;

    if ( a < 0.0 ) return( LGM_FILL_VALUE );

    sa = sin( a*RadPerDeg );    // sin(Alpha)
    sa2 = sa*sa;                // sin^2(Alpha)

    if ( sa2 < 1e-12 ) {
        return( 9e99 );
    } else {
        Ek = sqrt( 2.0*E0*B*Mu/(sa2*nT_Per_Gauss) + E0*E0) - E0;
        return( Ek );
    }


}


/**
 *  \brief
 *      Compute \f$ p^2c^2 \f$  given a particle's kinetic energy and rest energy.
 *  \details
 *
 *      Some relativistic equations:
 *
 *      With,
 *              \f{eqnarray*}
 *                      m       &=& \mbox{Particle mass.}\\
 *                      m_\circ &=& \mbox{Particle rest mass.}\\
 *                      v       &=& \mbox{Particle speed.}\\
 *                      c       &=& \mbox{Speed of Light.}\\
 *                      \gamma  &=& (1-v^2/c^2)^{-1/2}.\\
 *
 *              \f}
 *
 *          \f[
 *              m = m_\circ \gamma = m_\circ (1-v^2/c^2)^{-1/2}
 *          \f]
 *
 *      Rearrange this to get,
 *          \f[
 *              (mc^2)^2 = (m_\circ c^2)^ + p^2c^2
 *          \f]
 *
 *      With,
 *              \f{eqnarray*}
 *                  E  = mc^2             &=& \mbox{Total Energy of particle} \\
 *                  E_\circ = m_\circ c^2 &=& \mbox{Rest Energy of particle} \\
 *                  p                     &=& \mbox{relativistic momentum of particle}
 *              \f}
 *
 *      this is,
 *          \f[
 *              E^2 = E_\circ^2 + (pc)^2
 *          \f]
 *
 *      so,
 *          \f[
 *              p^2c^2  = E^2 - E_\circ^2
 *          \f]
 *
 *      Let \f$ E = E_k+E_\circ \f$ (kinetic energy + rest energy). Then,
 *              \f{eqnarray*}
 *                  p^2c^2  &=& (E_k+E_\circ)^2 - E_\circ^2 \\
 *                          &=& E_k (E_k+2E_\circ)
 *              \f}
 *
 *
 *      \param[in]      Ek (\f$ = E_k) \f$  Kinetic Energy of Particle.   <b>( MeV )</b>
 *      \param[in]      E0 (\f$ = E_\circ) \f$ Rest energy of Particle.   <b>( MeV )</b>
 *
 *      \return         p2c2 (\f$ = p^2c^2 = E_k(E_k+2E_\circ)\f$)        <b>( MeV^2 )</b>
 *
 *      \author         Mike Henderson
 *      \date           2010-2011
 */
double  Lgm_p2c2( double Ek, double E0 ) {
    return( Ek*(Ek+2.0*E0) );    // p^2c^2 in units of MeV^2
}


/**
 *  \brief
 *      Returns \f$ \beta^2 = (v/c)^2 \f$ as a function of \f$E_k\f$ and \f$E_\circ\f$. 
 *  \details
 *      Uses the following relation,
 *      
 *
 *          \f{eqnarray*}
 *              \beta^2 &=& 1 - {1\over \gamma^2} \\
 *                      &=& 1 - {1\over (1+E_k/E_\circ)^2}
 *          \f}
 *
 *
 *      \param[in]    Ek   Kinetic Energy of Particle.                  <b>( arbitrary units )</b>
 *      \param[in]    E0 (\f$ = E_\circ) \f$ Rest energy of Particle.   <b>( units of Ek )</b>
 *
 *      \return       (\f$ \beta^2 = v^2/c^2)\f$)                      <b>( dimensionless )</b>
 *
 *      \author         Mike Henderson
 *      \date           2010-2011
 */
double  Lgm_v2overc2( double Ek, double E0 ) {
    double  gamma = 1.0 + Ek/E0;
    return( 1.0 - 1.0/(gamma*gamma) );    // dimensionless
}


/**
 *  \brief
 *      Computes relativistic factor \f$ \gamma = [ 1 - (v/c)^2 ]^{-1/2} \f$.
 *  \details
 *       Returns relativistic gamma factor,
 *            \f$\gamma=[1-(v/c)^2]^{-1/2}\f$ .
 *
 *       Since \f$E=\gamma m_\circ c^2 = \gamma E_\circ\f$ and 
 *       \f$E=E_k+E_\circ\f$, we have;
 *
 *              \f[
 *                  \gamma = { E_k\over E_\circ } + 1
 *              \f]
 *
 *
 *      \param[in]    Ek   Kinetic Energy of Particle.                  <b>( arbitrary units )</b>
 *      \param[in]    E0 (\f$ = E_\circ) \f$ Rest energy of Particle.   <b>( units of Ek )</b>
 *
 *      \return       \f$\gamma = [ 1 - (v/c)^2 ]^{-1/2} \f$)          <b>( dimensionless )</b>
 *
 *      \author         Mike Henderson
 *      \date           2010-2011
 */
double  Lgm_gamma( double Ek, double E0 ) {
    return( 1.0 + Ek/E0 );  // dimensionless
}


/**
 *  \brief
 *      Convert differential flux to phase space density.
 *  \details
 *      The basic relationship is;
 *          \f[
 *          f = {j \over p^2}
 *          \f]
 *
 *      Multiplying the top and bottom by \f$ c^2 \f$ gives,
 *          \f{eqnarray*}
 *              f &=& {j c^2 \over p^2c^2 } \\
 *              f &=& { j\over c } {c^3 \over (p^2c^2)}
 *          \f}
 *
 *      The reason for making it \f$ c^3 \f$ is that the final units simplify to,
 *      \f$ c^3 \mbox{cm}^{-3} \mbox{MeV}^{-3}\f$ or \f$ \left[c\over \mbox{cm}
 *      \mbox{MeV}\right]^3 \f$. These are the standard GEM phase space density
 *      units.
 *
 *      \param[in]      j               Differential Flux in units of   <b>#/cm^2/s/sr/MeV</b>
 *      \param[in]      p2c2 (\f$ = p^2 c^2\f$) in units of             <b>Mev^2</b>
 *
 *      \return         f, Phase space density in units of              <b>(c/cm/MeV)^3</b>
 *
 *      \author         Mike Henderson
 *      \date           2010-2011
 *
 */
double Lgm_DiffFluxToPsd( double j, double p2c2 ){
    return( j/(p2c2*2.9979e10) ); // f in units of (c/cm/MeV)^3
}


/**
 *  \brief
 *      Convert differential flux to (non relativistic velocity space) phase space density.
 *  \details
 *      The basic relationship is;
 *          \f[
 *          f = {(m^2\over 2E)} j 
 *          \f]
 *
 *      \param[in]      j               Differential Flux in units of   <b>#/cm^2/s/sr/eV</b>
 *      \param[in]      m               non-relativistic mass           <b>kg</b>
 *      \param[in]      E               non-relativistic Energy         <b>eV</b>
 *
 *      \return         f, Phase space density in units of              <b>s^3/cm^6</b>
 *
 *      \author         Mike Henderson
 *      \date           2021
 *
 */
double Lgm_DiffFluxToPsd2( double j, double m, double E ){

    double eVPerJoule, g, K, f;

    // f = (m^2/(2E)) * j
    // convert mass fro kg to eV^2 s^4/cm^4
    // 1J = 1kg m^2/s^2 = 6.242e18 eV
    eVPerJoule = 6.241509e18; // eV/J
    g = m * eVPerJoule;
    K = 0.5*g*g; // ev^2 s^4/m^4
    K /= 1.0e8; // ev^2 s^4/cm^4

    f = K * j/E;

    return( f ); // f in units of s^3/cm^6
}



/**
 *  \brief
 *      Convert phase space density to differential flux.
 *  \details
 *      The basic relationship is;
 *          \f[
 *          f = j/p^2
 *          \f]
 *
 *      Multiply top and bottom by \f$ c^2 \f$ gives,
 *          \f{eqnarray*}
 *              f &=& {j c^2 \over p^2c^2 } \\
 *              f &=& { j\over c } {c^3 \over (p^2c^2)}
 *          \f}
 *
 *      The reason for making it \f$ c^3 \f$ is that the final units simplify to,
 *      \f$ c^3 \mbox{cm}^{-3} \mbox{MeV}^{-3}\f$ or \f$ \left[c\over \mbox{cm} \mbox{MeV}\right]^3 \f$
 *
 *      \param[in]      f    Phase space density in units of    <b>(c/cm/MeV)^3</b>
 *      \param[in]      p2c2 (\f$ = p^2 c^2\f$) in units of     <b>Mev^2</b>
 *
 *      \return         j, Differential Flux in units of        <b>#/cm^2/s/sr/MeV</b>
 *
 *      \author         Mike Henderson
 *      \date           2010-2011
 *
 */
double Lgm_PsdToDiffFlux( double f, double p2c2 ){
    return( f*2.9979e10*p2c2 ); // j in units of #/cm^2/s/sr/MeV
}



/**
 *  \brief
 *      Returns a pointer to a dynamically allocated Lgm_FluxToPsd structure.
 *  \details
 *      User must destroy this with Lgm_F2P_FreeFluxToPsd() when done.
 *
 *      \param[in]      DumpDiagnostics  Boolean flag to turn on/off dumping of diagnostics.
 *      \return         A pointer to an allocated and initialized Lgm_FluxToPsd stucture.
 *
 *      \author         Mike Henderson
 *      \date           2010-2011
 *
 */
Lgm_FluxToPsd *Lgm_F2P_CreateFluxToPsd( int DumpDiagnostics ) {

    Lgm_FluxToPsd *f;

    /*
     * Allocate memory for a Lgm_PhaseSpaceDensity structure.
     */
    f = (Lgm_FluxToPsd *) calloc( 1, sizeof(*f) );

    /*
     * Set DumpDiagnostics flag to what we got here. This can be changed later as well.
     */
    f->DumpDiagnostics = DumpDiagnostics;

    f->Extrapolate = TRUE;
    f->nMaxwellians = 2;
    f->FitType = LGM_F2P_SPLINE;
    f->FitType = LGM_F2P_MAXWELLIAN;

    f->Alloced1 = FALSE;
    f->Alloced2 = FALSE;

    f->UseModelB = TRUE; // default.

    return f;

}

/**
 *  \brief
 *      Destroy a dynamically allocated Lgm_FluxToPsd structure. 
 *  \details
 *      Destroys a dynamically allocated Lgm_FluxToPsd structure that was
 *      created by Lgm_F2P_CreateFluxToPsd().
 *
 *      \param          f  Pointer to the allocated Lgm_FluxToPsd structure that you want to destroy.
 *
 *      \author         Mike Henderson
 *      \date           2010-2011
 *
 */
void Lgm_F2P_FreeFluxToPsd( Lgm_FluxToPsd *f ) {

    if ( f->Alloced1 ) {
        LGM_ARRAY_1D_FREE( f->E );
        LGM_ARRAY_1D_FREE( f->A );
        LGM_ARRAY_2D_FREE( f->FLUX_EA );
        LGM_ARRAY_2D_FREE( f->PSD_EA );
    }

    if ( f->Alloced2 ) {
        LGM_ARRAY_1D_FREE( f->Mu );
        LGM_ARRAY_1D_FREE( f->K );
        LGM_ARRAY_1D_FREE( f->AofK );
        LGM_ARRAY_2D_FREE( f->EofMu );
        LGM_ARRAY_2D_FREE( f->PSD_MK );
    }

    free( f );

    return;
}


/**
 *  \brief
 *      Set Date/Time and position in the Lgm_FluxToPsd structure.
 *  \details
 *
 *
 *      \param[in]      d   Date/Time of measurement.
 *      \param[in]      u   Position of measurment (in GSM).
 *      \param[in,out]  f   Lgm_FluxToPsd sturcture.
 *
 *      \author         Mike Henderson
 *      \date           2010-2011
 *
 */
void Lgm_F2P_SetDateTimeAndPos( Lgm_DateTime *d, Lgm_Vector *u, Lgm_FluxToPsd *f ) {

    f->DateTime = *d;
    f->Position = *u;

}

/**
 *  \brief
 *      Forces Mu/E conversions to use observed magntitude of B.
 *  \details
 *
 *
 *      \param[in]      |B|   Magnitude of observed B-field (nT).
 *      \param[in,out]   f   Lgm_FluxToPsd sturcture.
 *
 *      \author         Mike Henderson
 *      \date           2014
 *
 */
void Lgm_F2P_SetObservedB( double B_obs, Lgm_FluxToPsd *f ) {

    f->B_obs     = B_obs;
    f->UseModelB = FALSE;

}




/**
 *  \brief
 *      Adds (to a Lgm_FluxToPsd structure) the user-supplied arrays containing
 *      J[Energy][Alpha],  Energy[], Alpha[]. Also computes the PSD (i.e. f as a
 *      function of energy and Alpha). Note that f=j/p^2, so this is not very
 *      much work. The hard work following these steps is to compute f at the
 *      given fixed mu's and K's and that requires interpolation/extrapolation,
 *      etc.
 *  \details
 *
 *      \param[in]      J                 2D array containing the differential flux values as a function of energy and pitch angle.
 *      \param[in]      E                 1D array containing the energy values implied by the first index of Flux[][] array.
 *      \param[in]      nE                number of energies.
 *      \param[in]      A                 1D array containing the pitch angles values implied by the second index of Flux[][] array.
 *      \param[in]      nA                number of pitch angles.
 *      \param[in,out]  f                 Lgm_FluxToPsd sturcture.
 *
 *      \author         Mike Henderson
 *      \date           2010-2011
 *
 */
void Lgm_F2P_SetFlux( double **J, double **dJ, double *E, int nE, double *A, int nA, Lgm_FluxToPsd *f ) {


    int     i, j;
    double  flux, dflux, p2c2, fp, dfp, Min, Max;


    /*
     * If arrays are already alloc'd, free them first
     */
    if ( f->Alloced1 ) {
        LGM_ARRAY_1D_FREE( f->E );
        LGM_ARRAY_1D_FREE( f->A );
        LGM_ARRAY_1D_FREE( f->Aeq );
        LGM_ARRAY_2D_FREE( f->FLUX_EA );
        LGM_ARRAY_2D_FREE( f->dFLUX_EA );
        LGM_ARRAY_2D_FREE( f->PSD_EA );
        LGM_ARRAY_2D_FREE( f->dPSD_EA );
    }


    /*
     * Add Flux array to f structure. Alloc arrays appropriately.
     */
    f->nE = nE; // Number of Energies
    f->nA = nA; // Number of Pitch Angles

    LGM_ARRAY_1D( f->E, f->nE, double );               // Values for the Energies
    LGM_ARRAY_1D( f->A, f->nA, double );               // Values for the Pitch Angles
    LGM_ARRAY_1D( f->Aeq, f->nA, double );             // Values for the Eq. Pitch Angles
    LGM_ARRAY_2D( f->FLUX_EA,  f->nE, f->nA, double ); // Differential Flux as a function of Energy and Pitch Angle.
    LGM_ARRAY_2D( f->dFLUX_EA, f->nE, f->nA, double ); // Unc. in Differential Flux as a function of Energy and Pitch Angle.
    LGM_ARRAY_2D( f->PSD_EA,  f->nE, f->nA, double );  // PSD as a function of Energy and Pitch Angle.
    LGM_ARRAY_2D( f->dPSD_EA, f->nE, f->nA, double );  // Unc. in PSD as a function of Energy and Pitch Angle.
    f->Alloced1 = TRUE;                                // Flag that we have alloced this mem.

    for (i=0; i<f->nE; i++) f->E[i] = E[i];
    for (i=0; i<f->nA; i++) f->A[i] = A[i];

    for (i=0; i<f->nE; i++) {
        for (j=0; j<f->nA; j++) {
            f->FLUX_EA[i][j]  =  J[i][j]; // FLUX_EA is "Differential Flux versus Energy and Pitch Angle".
            f->dFLUX_EA[i][j] = dJ[i][j]; // uncertainty in FLUX_EA 
        }
    }


    if ( f->DumpDiagnostics ) { DumpGif( "Lgm_FluxToPsd_FLUX_EA", f->nA, f->nE, f->FLUX_EA ); }


    /*
     * Convert Flux array into to PSD array.
     * Note, the result here is not PSD at constant Mu and K, it is PSD at the
     * same Es and Alphas we started with.
     */
    for (j=0; j<f->nA; j++) {
        for (i=0; i<f->nE; i++) {

            flux   = f->FLUX_EA[i][j];
            dflux  = f->dFLUX_EA[i][j];

            p2c2   = (f->E[i] >= 0.0 ) ? Lgm_p2c2( f->E[i], LGM_Ee0 ) : LGM_FILL_VALUE;

            fp = Lgm_DiffFluxToPsd( flux, p2c2 );
            f->PSD_EA[i][j] = fp; // PSD_EA is "PSD versus Energy and Pitch Angle".

            dfp = Lgm_DiffFluxToPsd( dflux, p2c2 );
            f->dPSD_EA[i][j] = dfp; // dPSD_EA is Uncertainty in PSD_EA

        }
    }
    if ( f->DumpDiagnostics ) {
        DumpGif( "Lgm_FluxToPsd_PSD_EA", f->nA, f->nE, f->PSD_EA );
    }


    return;

}




/**
 *  \brief
 *      Computes Phase Space Density at user-supplied constant values of \f$\mu\f$
 *      and K.
 *
 *  \details
 *      This routine ( Lgm_F2P_GetPsdAtConstMusAndKs() ) must operate on a
 *      pre-initialized Lgm_FluxToPsd structure.  The routine Lgm_FluxToPsd_SetFlux()
 *      is used to add differential flux data/info to an Lgm_FluxToPsd structure
 *      and it also converts the differential flux to Phase Space Density at
 *      constant E and \f$\alpha\f$ (i.e. it computes \f$ f(E, \alpha) \f$ ).
 *    
 *      However, what we really want is Phase Space Density at constant \f$\mu and
 *      K\f$ (i.e. we want \f$f( \mu, K )\f$). Although \f$mu\f$ is easy to
 *      compute, it is dependant on both E and \f$\alpha\f$. K is only dependant
 *      upon \f$\alpha\f$, but on the other hand K is not so easy to compute from
 *      \f$\alpha\f$.
 *    
 *      To perform the calculation we note that \f$f( E, \alpha )\f$ is the same as
 *      \f$f( E(\mu, \alpha(K)), \alpha(K) )\f$. Thus, for a given \f$\mu\f$ and K,
 *      we can figure out what E and \f$\alpha\f$ they correspond to and then we
 *      can just use the \f$f(E, \alpha)\f$ array to compute the desired f values.
 *      The steps are;
 *    
 *          - For each K, compute \f$\alpha(K)\f$. This is done with the routine
 *            Lgm_AlphaOfK().
 *    
 *          - Then we compute E from \f$\alpha\f$ and the given mu and Alpha
 *            values.
 *    
 *          - Then we look up \f$f(E, \alpha)\f$ from the array (interp or fit or
 *            whatever).
 *    
 *      \param[in]      nMu         Number of Mu values
 *      \param[in]      Mu          1-D array of Mu values
 *      \param[in]      nK          Number of K values
 *      \param[in]      K           1-D array of K values
 *      \param[in]      Extrapolate Turns on/off extrapolation capability
 *      \param[in,out]  f           A properly pre-initialized Lgm_FluxToPsd structure.
 *
 *      \author     Mike Henderson
 *      \date       2011
 *      \warning    Still working on this code. It is not finished.
 *
 */
void Lgm_F2P_GetPsdAtConstMusAndKs( double *Mu, int nMu, double *K, int nK, Lgm_MagModelInfo *mInfo, Lgm_FluxToPsd *f ) {

    int     i, j, k, m, DoIt, nGood, nGood2;
    double  AlphaEq, SqrtBrat, SinA, SinAeq, dPsd;

    int     n, n_star, ngy_sum;
    double  xin[3000], yin[3000], dyin[3000];
    double  x_star[300];
    int     x_star_index[300];
    double  gsig, gymax, gymin, gy_hat, gdy, gy_sum, gy_avg;
    GprInfo *Info;

    Lgm_MagModelInfo    *mInfo2; // This is used to hold thread-safe copies of mInfo.


    /*
     * If arrays are already alloc'd, free them first
     */
    if ( f->Alloced2 ) {
        LGM_ARRAY_1D_FREE( f->Mu );
        LGM_ARRAY_1D_FREE( f->K );
        LGM_ARRAY_1D_FREE( f->AofK );
        LGM_ARRAY_1D_FREE( f->AEqofK );
        LGM_ARRAY_2D_FREE( f->EofMu );
        LGM_ARRAY_2D_FREE( f->PSD_EAeq );
        LGM_ARRAY_2D_FREE( f->dPSD_EAeq );
        LGM_ARRAY_2D_FREE( f->PSD_MK );
        LGM_ARRAY_2D_FREE( f->dPSD_MK );
    }

    /*
     * Alloc arrays, and flag that we have done this.
     */
    f->nMu = nMu;   // Number of Mu values -- specified by user
    f->nK  = nK;    // Number of K values -- specified by user
    LGM_ARRAY_1D( f->Mu,        f->nMu, double ); // Values for the Mu's -- user specified
    LGM_ARRAY_1D( f->K,         f->nK,  double ); // Values for the K's -- user specified
    LGM_ARRAY_1D( f->AofK,      f->nK,  double ); // Array to hold the Local Pitch Angle corresponding to the K's
    LGM_ARRAY_1D( f->AEqofK,    f->nK,  double ); // Array to hold the Equatorial Pitch Angle corresponding to the K's
    LGM_ARRAY_2D( f->EofMu,     f->nMu, f->nK,  double ); // 2D array holding the Energy that corresponds to a Mu,K pair
    LGM_ARRAY_2D( f->PSD_EAeq,  f->nE,  f->nK,  double ); // Phase Space Density as a function of E and Aeq
    LGM_ARRAY_2D( f->dPSD_EAeq, f->nE,  f->nK,  double ); // Unc. in Phase Space Density as a function of E and Aeq
    LGM_ARRAY_2D( f->PSD_MK,    f->nMu, f->nK,  double ); // Phase Space Density as a function of Mu and K
    LGM_ARRAY_2D( f->dPSD_MK,   f->nMu, f->nK,  double ); // Unc. in Phase Space Density as a function of Mu and K
    f->Alloced2 = TRUE; // Flag that e have alloced this mem.



    /*
     * Copy K's (given in the arguments) into f structure.
     * Transform the K's into Alpha's using Lgm_AlphaOfK().
     * Save the results in the f structure.
     */
    if ( Lgm_Setup_AlphaOfK( &(f->DateTime), &(f->Position), mInfo ) > 0 ) {


        f->B   = mInfo->Blocal;
        f->Beq = mInfo->Bmin;

        { // start openmp parallel execution 

#if USE_OPENMP
            #pragma omp parallel private(mInfo2,AlphaEq,SinA)
            #pragma omp for schedule(dynamic, 1)
#endif
            for ( k=0; k<nK; k++ ){

                mInfo2 = Lgm_CopyMagInfo( mInfo );  // make a private (per-thread) copy of mInfo

                f->K[k]    = K[k];
                //printf("\n\nK[%d] = %g   f->DateTime.UTC = %g f->Position = %g %g %g\n", k, K[k], f->DateTime.Time, f->Position.x, f->Position.y, f->Position.z);
                AlphaEq    = Lgm_AlphaOfK( f->K[k], mInfo2 ); // Lgm_AlphaOfK() returns equatorial pitch angle.

                // Save the AlphaEq values
                f->AEqofK[k] = AlphaEq; //These are the a_eq vals that correspond to our desired K values.


                SinA       = sqrt( mInfo2->Blocal/mInfo2->Bmin ) * sin( RadPerDeg*AlphaEq );
                if ( AlphaEq > 0.0 ) {
                    if ( SinA <= 1.0 ) {
                        f->AofK[k] = DegPerRad*asin( SinA );
                    } else {
                        f->AofK[k] = LGM_FILL_VALUE;
                        printf("Particles with Eq. PA of %g mirror below us. (I.e. S/C does not see Ks this low).\n", AlphaEq);
                    }
                } else {
                    f->AofK[k] = LGM_FILL_VALUE;
                    printf("Particles with K of %g mirror below LC height. (I.e. S/C does not see Ks this high).\n", f->K[k]);
                }
                //printf("f->K[k] = %g   AlphaEq = %g SinA = %g f->AofK[k] = %g\n", f->K[k], AlphaEq, SinA, f->AofK[k]);

                Lgm_FreeMagInfo( mInfo2 ); // free mInfo2


            }

        } // end parallel




        Lgm_TearDown_AlphaOfK( mInfo );

    } else {

        // Blocal will have been set in Lgm_Setup_AlphaOfK() even if it returned a value <= 0.
        f->B = mInfo->Blocal;
        for ( k=0; k<nK; k++ ) {
            f->AEqofK[k] = LGM_FILL_VALUE;
            f->AofK[k]   = LGM_FILL_VALUE;
        }

    }

    // how many good values did we get?
    for ( nGood=0, k=0; k<nK; k++ ){
        if ( f->AEqofK[k] > 0.0 ) ++nGood;
    }

    if (  !(f->UseModelB)  ) {
        printf("Bmodel, Bobs = %g %g\n", f->B, f->B_obs );
    }




    /*
     * Create the Flux[E][a_eq] array using GPR interpolation.
     *
     * This was added Oct, 2020.
     * The original algorithm was:
     *      1) for each Mu and K, we figure out what E and Alpha_local we need.
     *      2) Then the problem is just one of interpolating the PSD_EA array
     *         to get the PSD at the impplied E.Alpha_local.
     *
     * This works, but obviously, if the S/C is not at the Bmin point, it
     * cannot see the lower K values (i.e. near 90-deg eq. pitch angles). In
     * other words, if you just take the local PAD and re-label the pitch
     * angles to their corresponding eq. pitch angles, a gap will open up
     * around 90-deg. This is a PITA because even if the gap is small, this
     * method is not ameniable to using interpolation, since the Alpha
     * implied by a given K will be undefined (because Blocal>Bmin).
     *
     * This new algorithm changes this:
     *      1) For each energy, the f( Alpha_local )  PAD is relabeled as f(
     *         Alpha_eq ) which typically opens a gap around 90.
     *      2) Then feed this data into our Gaussian Process Regression Routine
     *         to interp to points at the AlphaEq implied by each K (i.e. instead
     *         of the Alpha_local's which may be undefined.)
     */


    // For each local pitch angle, compute the equatorial pitch angle.
    SqrtBrat = sqrt( mInfo->Bmin/mInfo->Blocal );
    for ( nGood2=0, j=0; j<f->nA; j++ ) { 
        SinAeq = SqrtBrat * sin( RadPerDeg*f->A[j] );
        //if ( SinA <= 1.0 ) {
        if ( SinAeq <= 1.0 ) {
            f->Aeq[j] = DegPerRad*asin( SinAeq );
            ++nGood2;
        } else {
            f->Aeq[j] = LGM_FILL_VALUE;
        }
    }


//FILE *fp_pad, *fp_gpr;
//fp_pad=fopen("PAD.dat","w");
//fp_gpr=fopen("GPR.dat","w");
    for ( i=0; i<f->nE; i++ ) { // for each energy

         // for each PAD
         n = 0; gymax = -9e99; gymin = 9e99;
         ngy_sum = 0; gy_sum = 0.0;

        /*
         * Try adding mirror of first half of PAD. We are adding the main PAD
         * in -90->90. But also adding:
         *      1) the mirror of the first half is placed -180->-90
         *      2) the mirror of the last half is placed 90->180
         *      3) f=0 is forced at -90 and 90. With unc. of the 90-deg. val.
         *         CHECK ON THIS -- we may want to scan to find largest unc.
         *         available (90deg may be undefined.?) Note that GPR wont force the curve through zero!
         */
        for ( j=f->nA/2; j>=0; j-- ) { 
            // add in mirror of first half of PAD
            if ( !isnan(f->Aeq[j]) && (f->PSD_EA[i][j] > 0.0 ) && (f->PSD_EA[i][j] < 1e20 ) && (f->dPSD_EA[i][j] > 0.0) ){
                xin[n]  = -f->Aeq[j] - 90.0;
                yin[n]  = f->PSD_EA[i][j];
                dyin[n] = f->dPSD_EA[i][j];
                ++n;
            }
        }

//        xin[n]  = -90.0;
//        yin[n]  = 0.0;
//        dyin[n] = f->dPSD_EA[i][f->nA/2]; // take unc. from 90deg PA.
//        ++n;

        for ( j=0; j<f->nA; j++ ) { // pitch angle bin
            if ( !isnan(f->Aeq[j]) && (f->PSD_EA[i][j] > 0.0 ) && (f->PSD_EA[i][j] < 1e20 ) && (f->dPSD_EA[i][j] > 0.0) ){
                xin[n]  = f->Aeq[j]-90.0; // Give it values -90 to 90
                yin[n]  = f->PSD_EA[i][j];
                dyin[n] = f->dPSD_EA[i][j];
                if ( yin[n] > gymax )      gymax = yin[n];
                else if ( yin[n] < gymin ) gymin = yin[n];
                ++n;
            }
        }
        for ( j=f->nA-1; j>=0; j-- ) { // pitch angle bin
            if ( !isnan(f->Aeq[j]) && (f->PSD_EA[i][j] > 0.0 ) && (f->PSD_EA[i][j] < 1e20 ) && (f->dPSD_EA[i][j] > 0.0) ){
                xin[n]  = 90.0-f->Aeq[j]; // Give it values 90 to 180
                yin[n]  = f->PSD_EA[i][j];
                dyin[n] = f->dPSD_EA[i][j];
                if ( yin[n] > gymax )      gymax = yin[n];
                else if ( yin[n] < gymin ) gymin = yin[n];
                ++n;
            }
        }

//        xin[n]  = 90.0;
//        yin[n]  = 0.0;
//        dyin[n] = f->dPSD_EA[i][f->nA/2]; // take unc. from 90deg PA.
//        ++n;

        // Try adding mirror of last half of PAD
        for ( j=0; j<=f->nA/2; j++ ) { // pitch angle bin
            if ( !isnan(f->Aeq[j]) && (f->PSD_EA[i][j] > 0.0 ) && (f->PSD_EA[i][j] < 1e20 ) && (f->dPSD_EA[i][j] > 0.0) ){
                xin[n]  = 90.0+f->Aeq[j]; // Give it values 180 to 270
                yin[n]  = f->PSD_EA[i][j];
                dyin[n] = f->dPSD_EA[i][j];
                ++n;
            }
        }
/*
for (j=0; j<n; j++){
//fprintf(fp_pad, "%g %g %g\n", xin[j], yin[j], dyin[j]);
if ( dyin[j]/yin[j]*100.0 < 1.0){
printf("i,j=%d,%d: Energy: %g  xin, yin, dyin, rel error = %g %g %g %g %\n", i, j, f->E[i],  xin[j], yin[j], dyin[j], dyin[j]/yin[j]*100.0);
}
}
*/
//fflush(fp_pad);

//printf("E. Lgm_F2P_GetPsdAtConstMusAndKs %d %d\n", n, nGood);
//printf("n, nGood=%d %d\n", n, nGood);
//exit(0);

        // its possible that some of the AEqofK may not have computed correctly...
        // so only use good values. But also remember what indices they are coming from 
        n_star = 0;
        for (j=0; j<f->nK; j++){
            if ( f->AEqofK[j] > 0.0 ) {
                x_star[n_star] = f->AEqofK[j]-90.0;
                x_star_index[n_star] = j; // remember what index this is from...
                ++n_star;
            }
        }

        if ( (n>4)&&(nGood>2) ) {
            /*
             * OK, now we need to define the array of points at which we want to
             * evaluate the function at.  For visualization purposes, we could do
             * this over a lot of regularly spaced points. But we really only need
             * them at the f->AEqofK[k] points we calculated above.
             */


            // Initialize the Info structure.
            Info = InitGprInfo( n, n_star, 1 ); // The "1" means force symmetryo
            for ( j=0; j<n_star; j++ ) { 
                gsl_vector_set( Info->x_star, j, x_star[j] );
            }

            for (j=0; j<n; j++ ) {
                gsl_vector_set( Info->x, j, xin[j] );
                gsl_vector_set( Info->y, j, yin[j]/gymax - 0.5 );
                gsig = dyin[j]/gymax;
                //printf("xin, yin = %g %g gsig = %g   gymax = %g\n", xin[j], yin[j], gsig, gymax);
                gsl_vector_set( Info->sigma_n_vec, j, gsig );
                gsl_vector_set( Info->sigma_n_2_vec, j, gsig*gsig );
            }

            // Do regression
            //printf("i=%d n=%d\n", i, n);
            GPR( Info );
            

            /*
             * Save results. The FLUX_EAeq array should be interpolated to all of
             * the required AlphaEq (corresponding to the Ks we want).
             */
            for ( j=0; j<f->nK; j++ ) f->PSD_EAeq[i][j] = LGM_FILL_VALUE; // Initialize to FILLs
            for ( j=0; j<n_star; j++ ) { 
                gy_hat = (gsl_vector_get( Info->y_hat, j ) +0.5)*gymax;
                gdy    = (gsl_vector_get( Info->y_cred_hi, j )+0.5)*gymax - gy_hat;
//printf("here   gy_hat = %g   gdy = %g\n", gy_hat, gdy);
                f->PSD_EAeq[i][ x_star_index[j] ]  = gy_hat;
                f->dPSD_EAeq[i][ x_star_index[j] ] = gdy;
            }
        

            // Clean up memory
            FreeGprInfo( Info );

        } else {
            for ( j=0; j<n_star; j++ ) { 
                f->PSD_EAeq[i][j]  = LGM_FILL_VALUE;
                f->dPSD_EAeq[i][j] = LGM_FILL_VALUE;
            }
        }
//printf("F. Lgm_F2P_GetPsdAtConstMusAndKs %d %d\n", n, nGood);
//exit(0);


    } // end energy loop
//fclose(fp_pad);
//fclose(fp_gpr);
//exit(0);



    /*
     * Copy Mu's (given in the arguments) into f structure.
     * Transform the Mu's into (Kinetic) Energies.
     * Save the results in the f structure.
     * Note that since this conversion involves Mu and Alpha, the result is 2D.
assumes electrons -- generalize this...
     */

    for ( m=0; m<nMu; m++ ){
        f->Mu[m] = Mu[m];
        for ( k=0; k<nK; k++ ){
            if ( f->UseModelB ) {
                //f->EofMu[m][k] = Lgm_Mu_to_Ek( f->Mu[m], f->AofK[k], f->B, LGM_Ee0 );
                f->EofMu[m][k] = Lgm_Mu_to_Ek( f->Mu[m], f->AEqofK[k], f->Beq, LGM_Ee0 );

            } else {
                f->EofMu[m][k] = Lgm_Mu_to_Ek( f->Mu[m], f->AofK[k], f->B_obs, LGM_Ee0 );
            }
            //printf("f->Mu[%d], f->K[%d], f->AofK[%d], f->B, f->EofMu[%d][%d] = %g %g %g %g %g\n", 
            //       m, k, k, m, k, f->Mu[m], f->K[k], f->AofK[k], f->B, f->EofMu[m][k]);
        }
    }



    /*
     * OLD ALGORITHM: 
     *      Now, from the PSD[E][a] array, get PSD at the E's and Alpha's we just computed.
     *      The result will be the same as PSD at the given Mu's and K's
     * 
     *  NEW ALGORITHM:
     *      Grab the PSD from the PSD_EAeq array instead (that we just computed above).
     */
    for ( m=0; m<nMu; m++ ){
        for ( k=0; k<nK; k++ ){
            DoIt = FALSE;
            //printf("f->EofMu[m][k] %g, f->E[0] %g \n", f->EofMu[m][k], f->E[0]);
f->Extrapolate=1;
f->Extrapolate=0;
f->Extrapolate=2;
            if ( f->Extrapolate > 2 ){ // extrapolate above and below

                DoIt = TRUE;

            } else if ( f->Extrapolate == 2) {

                if (f->EofMu[m][k] >= f->E[0] ) DoIt = TRUE; // extrapolate above

            } else if ( f->Extrapolate == 1) {

                if (f->EofMu[m][k] <= f->E[f->nE-1] ) DoIt = TRUE; // extrapolate below

            } else if ( f->Extrapolate == 0) {

                if ((f->EofMu[m][k] >= f->E[0])&&(f->EofMu[m][k] <= f->E[f->nE-1])) DoIt = TRUE; // interp only


            }

            if (DoIt) {

                // This call is where the OLD interpolation would get done.
                // Note that the one on pitch angle will no longer be needed....
                // The routine Lgm_F2P_GetPsdAtEandAlpha2() was added for this NEW way.
                //f->PSD_MK[m][k] =  Lgm_F2P_GetPsdAtEandAlpha( m, k, f->EofMu[m][k], f->AofK[k], &dPsd, f );
                f->PSD_MK[m][k] =  Lgm_F2P_GetPsdAtEandAlpha2( m, f->EofMu[m][k], k, &dPsd, f );
                f->dPSD_MK[m][k] = dPsd;

            } else {

                f->PSD_MK[m][k]  = LGM_FILL_VALUE;
                f->dPSD_MK[m][k] = LGM_FILL_VALUE;

            }

        }
    }

    if ( f->DumpDiagnostics ) { DumpGif( "Lgm_FluxToPsd_PSD_MK", f->nK, f->nMu, f->PSD_MK ); }

    return;

}




double  Model( double *x, int n, double E ) {

    int i;
    double  nn, TT, val;

    // Sum of maxwellians
    for ( val=0.0, i=0; i<n; i++){
        nn = pow( 10.0,  x[2*i+1] );
        TT = fabs( x[2*i+2] );
        val += Lgm_MaxJut( nn, TT, E, LGM_Ee0 );
    }

    return( val );

}

double Cost( double *x, void *data ){

    _FitData    *FitData;
    int         i;
    double      g_model, d, sum;

    FitData = (_FitData *)data;

    // These constraints wont work for nMaxwellians>2
/*
    if ( (x[1] > 2.0) || ( x[1] < -30.0) ) return( 9e99 );
    if ( (fabs( x[2] ) > 1000.0) || (fabs( x[2] ) < 1.0) ) return( 9e99 );
    if ( FitData->nMaxwellians > 1 ) {
        if ( (x[3] > 2.0) || ( x[3] < -30.0) ) return( 9e99 );
        if ( (fabs( x[4] ) > 1000.0) || (fabs( x[4] ) < 1.0) ) return( 9e99 );
    }
*/





    for ( sum = 0.0, i=0; i<FitData->n; ++i ){

        g_model = Model( x, FitData->nMaxwellians, FitData->E[i] ) ;
        d = log10( FitData->g[i]) - log10( g_model);

        sum += d*d;
        //sum += fabs(d);
        ///d = (d >= 0.0) ? d : -d;
        ///sum += d;
if (isinf(sum)) {
//    printf("Cost, INF: g_model, g = %g %g %g     x[1], x[2] = %g %g\n", g_model, FitData->g[i], log10( FitData->g[i] ), x[1], x[2]);
    return( 9e99 );
}
if (isnan(sum)) {
//    printf("Cost, NaN: g_model, g = %g %g %g     x[1], x[2] = %g %g\n", g_model, FitData->g[i], log10( FitData->g[i] ), x[1], x[2]);
//    printf("FitData->n = %d\n", FitData->n);
    return( 9e99 );
}

    }


    return( sum );

}



/**
 * The f structure should have an initialized PSD[E][a] array in it.
 * This routine computes psd given a value of E and a. dPsd is the uncertainty in Psd (return value).
 */
double  Lgm_F2P_GetPsdAtEandAlpha( int iMu, int iK, double E, double a, double *dpsd, Lgm_FluxToPsd *f ) {

//NEED to get dPsd computed!!!

    int         j, i, i0, i1, nn;
    double      a0, a1, y0, y1, slp, psd, g;
    _FitData    *FitData;

    // if a < 0, we should return fill value.
    if ( a < 0.0 ) return(LGM_FILL_VALUE);

    FitData = (_FitData *) calloc( 1, sizeof( _FitData ) );
    FitData->nMaxwellians = f->nMaxwellians;

//FROM HERE 
    /*
     * Since pitch angle, a is bounded (here its constrained to be between 0
     * and 90), we will interpolate on that first to produce a 1D array of
     * f(E).
     */
    if ( a < f->A[0] ) {
        return(LGM_FILL_VALUE);
        //i0 = 0; i1 = 1;
    } else if ( a > f->A[f->nA - 1] ) {
        return(LGM_FILL_VALUE);
        //i0 = f->nA - 2; i1 = f->nA - 1;
    } else {
        for (i=1; i<f->nA; i++) {
            if ( a < f->A[i] ) {
                i0 = i-1; i1 = i;
                break;
            }
        }
    }
    //printf("i0, i1 = %d %d\n", i0, i1);


    /************************
     *  interpolate PA 
     *  For each energy, we interpolate the PA profile to get f at the user-specified alpha.
     ************************/
    LGM_ARRAY_1D( FitData->E, f->nE, double );
    LGM_ARRAY_1D( FitData->g, f->nE, double );
    FitData->n = 0;
    for (j=0; j<f->nE; ++j){

        a0   = f->A[i0];
        a1   = f->A[i1];
        y0   = f->PSD_EA[j][i0];
        y1   = f->PSD_EA[j][i1];

        if ( (y0>0.0)&&(y1>0.0) ){ // if one of the 'y' vals is zero, dont use this point.
            slp  = (y1-y0)/(a1-a0);
            g = slp*(a-a1) + y1;
            if ( g > 1e-40 ) { // dont use if it looks bogus
                FitData->g[ FitData->n ] = g;
                FitData->E[ FitData->n ] = f->E[j];
//printf("j=%d a = %g E=%g i0, i1 = %d %d   a0, a1 = %g %g   y0, y1, slp = %g %g %g    a = %g, FitData->E[%d] = %g FitData->g[%d] = %g\n", j, a, E, i0, i1, a0, a1, y0, y1, slp, a, FitData->n, FitData->E[ FitData->n ], FitData->n, FitData->g[FitData->n]);
                ++(FitData->n);
            }
        }


    }
// TO HERE 
//The above could just be replaced by picking out the right index to the f->PSD_EAeq[][] array now...
// Basically the g value is just read out of the array -- the indices are already at the right alpha val.
// E.g.
// g = f->PSD_EAeq[j][i] where i is given to us. I.e. we dont need all the crap to do the interp above...





    /*
     * Now we have f versus E (at the specified alpha).
     */

    if ( f->FitType == LGM_F2P_SPLINE ) {

        /*
         * Use smoothing spline.
         */
        int n;
        //int ncoeffs = 12;
        int ncoeffs = 8;
        int nbreak  = ncoeffs - 2;

        // determine number of points
        for (n=0, j=0; j<FitData->n; ++j){ 
            //if ( (f->E[j] > 1.0/1000.0) && ( FitData->g[j] > 0.0 ) ) {
            if ( (f->E[j] > 0.0/1000.0) && ( FitData->g[j] > 0.0 ) ) {
                //printf("E[%d] = %g\n", j, f->E[j] );
                ++n;
            }
        }

        if ( n > 20 ) {
            ncoeffs = 12;
        } else if ( n > 5 ) {
            ncoeffs = 4;
        } else {
            ncoeffs = n/2;
        }
        nbreak  = ncoeffs-2;
        
//if ((iK==4)&&(iMu==7)) printf("n, ncoeffs, nbreak = %d %d %d\n", n, ncoeffs, nbreak );

        if ( (n >= 4) && (nbreak>=2) ) {

            // allocate a cubic bspline workspace (k = 4)
            gsl_bspline_workspace *bs_bw;
            gsl_vector            *bs_B;
            bs_bw = gsl_bspline_alloc( 4, nbreak );
            bs_B  = gsl_vector_alloc( ncoeffs );



            gsl_vector *bs_x, *bs_y, *bs_c, *bs_w;
            gsl_matrix *bs_X, *bs_cov;
            gsl_multifit_linear_workspace *bs_mw;
            double  bs_chisq, xi, yi, yierr, sigma;
            
            
            bs_x = gsl_vector_alloc( n );
            bs_y = gsl_vector_alloc( n );
            bs_X = gsl_matrix_alloc( n, ncoeffs );
            bs_c = gsl_vector_alloc( ncoeffs );
            bs_w = gsl_vector_alloc( n );
            bs_cov = gsl_matrix_alloc( ncoeffs, ncoeffs );
            bs_mw = gsl_multifit_linear_alloc( n, ncoeffs );

            // set up data arrays.
            for ( nn=0, j=0; j<FitData->n; j++ ){

                //if ( (f->E[j] > 1.0/1000.0) && ( FitData->g[j] > 0.0 ) ) {
                if ( (FitData->E[j] > 0.0/1000.0) && ( FitData->g[j] > 0.0 ) && ( FitData->dg[j] > 0.0 )) {
                    xi    = log10( FitData->E[j] );
                    yi    = log10( FitData->g[j] );
                    sigma = .434*( FitData->dg[j]/FitData->g[j] );

                    gsl_vector_set( bs_x, nn, xi );
                    gsl_vector_set( bs_y, nn, yi );
                    gsl_vector_set( bs_w, nn, 1.0/(sigma*sigma) );

                    ++nn;
                }

            }

            // use uniform breakpoints on defined data interval
            gsl_bspline_knots_uniform( gsl_vector_get( bs_x, 0 ), gsl_vector_get( bs_x, n-1 ), bs_bw );

            // construct the fit matrix X
            for ( i=0; i<n; i++ ) {
                xi = gsl_vector_get( bs_x, i );

                // compute B_j(xi) for all j
                gsl_bspline_eval( xi, bs_B, bs_bw );

                // fill in row i of X
                for ( j=0; j<ncoeffs; j++ ) {
                    double Bj = gsl_vector_get( bs_B, j );
                    gsl_matrix_set( bs_X, i, j, Bj );
                }
            }

            /* do the fit */
            gsl_multifit_wlinear( bs_X, bs_w, bs_y, bs_c, bs_cov, &bs_chisq, bs_mw );


            /*
            FILE    *fp;
            printf("E = %g\n", E);
            fp = fopen("data.txt", "w");
            for (j=0; j<FitData->n; ++j){
                fprintf(fp, "%g %g\n", log10(f->E[j]), log10(FitData->g[j]));
            }
            fclose(fp);
            fp = fopen("fit.txt", "w");
            for (xi = gsl_vector_get( bs_x, 0 ); xi < gsl_vector_get( bs_x, n-1 ); xi += 0.01) {
                gsl_bspline_eval( xi, bs_B, bs_bw );
                gsl_multifit_linear_est( bs_B, bs_c, bs_cov, &yi, &yierr );
                fprintf(fp, "%g %g\n", xi, yi );
            }
            fclose(fp);
            exit(0);
            */

            xi = log10( E );

            if ( (xi >= gsl_vector_get( bs_x, 0 )) && (xi <= gsl_vector_get( bs_x, n-1 )) ) {
                gsl_bspline_eval( xi, bs_B, bs_bw );
                gsl_multifit_linear_est( bs_B, bs_c, bs_cov, &yi, &yierr );
                psd  = pow( 10.0, yi );
                *dpsd = 2.303*psd*yierr;
            } else {
                psd  = LGM_FILL_VALUE;
                *dpsd = LGM_FILL_VALUE;
            }

            gsl_bspline_free( bs_bw );
            gsl_vector_free( bs_B );
            gsl_vector_free( bs_x );
            gsl_vector_free( bs_y );
            gsl_matrix_free( bs_X );
            gsl_vector_free( bs_c );
            gsl_vector_free( bs_w );
            gsl_matrix_free( bs_cov );
            gsl_multifit_linear_free( bs_mw );



        } else {

            psd = LGM_FILL_VALUE;
            *dpsd = LGM_FILL_VALUE;

        }

        // spline fit to energy spectra is done.


    } else if ( f->FitType == LGM_F2P_MAXWELLIAN ) {
//Uncertainties are not propery set for thej maxwellian fits.
// Currently we use praxis -- a line minimization method, but we really should go to a differentm one
// in order to be able to get unc as well.
        if ( FitData->n > 2 ) {


            // interpolate/fit E
            // for now just do a linear interp.
            // no lets try a fit...
            double  in[10], out[7], x[10];
            in[0] = 1e-8;
            in[1] = in[2] = 1e-9; //Info->Praxis_Tolerance;
            in[5] = 30000.0; //(double)Info->Praxis_Max_Function_Evals;
            in[6] = 10.0; //Info->Praxis_Maximum_Step_Size;
            in[7] = 10.0; //Info->Praxis_Bad_Scale_Paramater;
            in[8] = 4.0; //(double)Info->Praxis_Max_Its_Without_Improvement;
            in[9] = 1.0; //(double)Info->Praxis_Ill_Conditioned_Problem;
            x[0] = 0.0;
            x[1] = -1.0;
            x[2] = 25.0;
            x[3] = -2.0;
            x[4] = 200.0;
            praxis( 2*FitData->nMaxwellians, x, (void *)FitData, Cost, in, out);
            /*
            printf("out[0] = %g\n", out[0]);
            printf("out[1] = %g\n", out[1]);
            printf("out[2] = %g\n", out[2]);
            printf("out[3] = %g\n", out[3]);
            printf("out[4] = %g\n", out[4]);
            printf("out[5] = %g\n", out[5]);
            printf("out[6] = %g\n", out[6]);
            */
            //printf("x[1] = %g   x[2] = %g   Cost = %g\n", x[1], x[2], out[6]);
 //           FILE *fp;
 //           printf("E = %g\n", E);
 //           fp = fopen("data.txt", "w");
 //           for (j=0; j<FitData->n; ++j){
 //               fprintf(fp, "%g %g\n", log10(f->E[j]), log10(FitData->g[j]));
 //           }
 //           fclose(fp);

 //           fp = fopen("fit.txt", "w");
 //           for (j=0; j<FitData->n; ++j){
 //               psd = Model( x,  FitData->nMaxwellians, f->E[j] );
 //               fprintf(fp, "%g %g\n", log10(f->E[j]), log10( psd ) );
 //           }
 //           fclose(fp);

//            exit(0);
            /*
            */


            //x[2] = 200.0;
            //printf("x - %g %g %g %g\n", x[1], x[2], x[3], x[4]);
            psd = Model( x,  FitData->nMaxwellians, E );
*dpsd = 0.0;
//dpsd = ?;
            //psd = (double)a;

            //printf("E, a = %g %g  x = %g %g psd = %g\n", E, a, x[1], x[2], psd);
        } else {

            psd  = LGM_FILL_VALUE;
            *dpsd = LGM_FILL_VALUE;

        }

    } else {
        //FitType not valid...
        psd  = LGM_FILL_VALUE;
        *dpsd = LGM_FILL_VALUE;
    }

    LGM_ARRAY_1D_FREE( FitData->E );
    LGM_ARRAY_1D_FREE( FitData->g );
    free( FitData );


    return( psd );

}


/**
 * This version uses the GP-interped array
 * Here we just need to know what index to grab from the AlphaEq slot...
 * The f structure should have an initialized PSD[E][a] array in it.
 * This routine computes psd given a value of E and a. dPsd is the uncertainty in Psd (return value).
 *
 */
double  Lgm_F2P_GetPsdAtEandAlpha2( int iMu, double E, int iAEq, double *dpsd, Lgm_FluxToPsd *f ) {


    int         j, i, nn;
    double      psd, g, dg;
    _FitData    *FitData;
    double MinEnergy, MaxEnergy;

    FitData = (_FitData *) calloc( 1, sizeof( _FitData ) );
    FitData->nMaxwellians = f->nMaxwellians;

    MinEnergy = 0.0/1000.0/1000.0; // 5eV in units of MeV
MinEnergy = 1.0;
    MaxEnergy = 20.0; // MeV

    /************************
     *  Grab the energy spectra at the given PA index
     ************************/
    LGM_ARRAY_1D( FitData->E, f->nE, double );
    LGM_ARRAY_1D( FitData->g, f->nE, double );
    LGM_ARRAY_1D( FitData->dg, f->nE, double );
    FitData->n = 0;
    for (j=0; j<f->nE; ++j){
        if ( (f->E[j] > MinEnergy) && (f->E[j] <= MaxEnergy) && ( f->PSD_EAeq[j][iAEq] > 0.0 ) && (f->dPSD_EAeq[j][iAEq] > 0.0) ) {
            g  = f->PSD_EAeq[j][iAEq];
            dg = f->dPSD_EAeq[j][iAEq];
//            if ( g > 1e-40 ) { // dont use if it looks bogus
                FitData->g[ FitData->n ]  = g;
                FitData->dg[ FitData->n ] = dg;
                FitData->E[ FitData->n ] = f->E[j];
                ++(FitData->n);
//printf("A. %g %g %g\n", g, dg, f->E[j]);
 //           }
        }
    }



    /*
     * Now we have f versus E (at the specified alpha).
     */

    if ( f->FitType == LGM_F2P_SPLINE ) {

        /*
         * Use smoothing spline.
         */
        int n;
        //int ncoeffs = 12;
        int ncoeffs = 8;
        int nbreak  = ncoeffs - 2;

        // determine number of points
        for (n=0, j=0; j<FitData->n; ++j){ 
            if ( (FitData->E[j] > MinEnergy) && (f->E[j] <= MaxEnergy) && ( FitData->g[j] > 0.0 ) && ( FitData->dg[j] > 0.0 )) {
                //printf("E[%d] = %g  g[%d] = %g\n", j, f->E[j], j, FitData->g[j] );
                ++n;
            }
        }

        if ( n > 20 ) {
            ncoeffs = 12;
        } else if ( n > 5 ) {
            ncoeffs = 4;
        } else {
            ncoeffs = n/2;
        }
        nbreak  = ncoeffs-2;
        
//if ((iK==4)&&(iMu==7)) printf("n, ncoeffs, nbreak = %d %d %d\n", n, ncoeffs, nbreak );

        if ( (n >= 4) && (nbreak>=2) ) {

            // allocate a cubic bspline workspace (k = 4)
            gsl_bspline_workspace *bs_bw;
            gsl_vector            *bs_B;
            bs_bw = gsl_bspline_alloc( 4, nbreak );
            bs_B  = gsl_vector_alloc( ncoeffs );



            gsl_vector *bs_x, *bs_y, *bs_c, *bs_w;
            gsl_matrix *bs_X, *bs_cov;
            gsl_multifit_linear_workspace *bs_mw;
            double  bs_chisq, xi, yi, yierr, sigma;
            
            
            bs_x = gsl_vector_alloc( n );
            bs_y = gsl_vector_alloc( n );
            bs_X = gsl_matrix_alloc( n, ncoeffs );
            bs_c = gsl_vector_alloc( ncoeffs );
            bs_w = gsl_vector_alloc( n );
            bs_cov = gsl_matrix_alloc( ncoeffs, ncoeffs );
            bs_mw = gsl_multifit_linear_alloc( n, ncoeffs );

            // set up data arrays.
            for ( nn=0, j=0; j<FitData->n; j++ ){

                if ( (FitData->E[j] > MinEnergy) && (f->E[j] <= MaxEnergy) && ( FitData->g[j] > 0.0 ) && ( FitData->dg[j] > 0.0 )) {
                    xi    = log10( FitData->E[j] );
                    yi    = log10( FitData->g[j] );
if (FitData->dg[j]/FitData->g[j]*100.0 < 10.0){
                    sigma = .0434;
} else {
                    sigma = .434*( FitData->dg[j]/FitData->g[j] );
}
                    //sigma = 0.2*yi;
//printf("xi, yi, 0.2*yi, .434*( FitData->dg[j]/FitData->g[j] ) = %g %g %g %g   g, dg = %g %g\n", xi, yi, 0.2*yi, .434*( FitData->dg[j]/FitData->g[j] ), FitData->g[j], FitData->dg[j]);
//printf("%g %g %g       %g %g %g\n", xi, yi, sigma, FitData->E[j], FitData->g[j], FitData->dg[j]/FitData->g[j]);

                    gsl_vector_set( bs_x, nn, xi );
                    gsl_vector_set( bs_y, nn, yi );
                    gsl_vector_set( bs_w, nn, 1.0/(sigma*sigma) );

                    ++nn;
                }

            }

            // use uniform breakpoints on defined data interval
            //printf("gsl_vector_get( bs_x, 0 ), gsl_vector_get( bs_x, n-1 ) = %g %g\n", gsl_vector_get( bs_x, 0 ), gsl_vector_get( bs_x, n-1 ));
            gsl_bspline_knots_uniform( gsl_vector_get( bs_x, 0 ), gsl_vector_get( bs_x, n-1 ), bs_bw );
            //printf("B. here\n");

            // construct the fit matrix X
//printf("n=%d nn=%d   ", n, nn);
//printf("n, ncoeffs, nbreak = %d %d %d\n", n, ncoeffs, nbreak );
//if (nn==92) FLAGFLAG = 1;
//if ((nn==93)&&(FLAGFLAG==1))exit(0);
            for ( i=0; i<n; i++ ) {
                xi = gsl_vector_get( bs_x, i );

                // compute B_j(xi) for all j
                gsl_bspline_eval( xi, bs_B, bs_bw );

                // fill in row i of X
                for ( j=0; j<ncoeffs; j++ ) {
                    double Bj = gsl_vector_get( bs_B, j );
//printf("Bj = %g\n", Bj);
                    gsl_matrix_set( bs_X, i, j, Bj );
                }
            }

            /* do the fit */
            gsl_multifit_wlinear( bs_X, bs_w, bs_y, bs_c, bs_cov, &bs_chisq, bs_mw );
//printf("bs_chisq = %g\n", bs_chisq);



            /*
            FILE    *fp;
            printf("E = %g\n", E);
            fp = fopen("data.txt", "w");
            for (j=0; j<FitData->n; ++j){
                fprintf(fp, "%g %g\n", log10(f->E[j]), log10(FitData->g[j]));
            }
            fclose(fp);
            fp = fopen("fit.txt", "w");
            for (xi = gsl_vector_get( bs_x, 0 ); xi < gsl_vector_get( bs_x, n-1 ); xi += 0.01) {
                gsl_bspline_eval( xi, bs_B, bs_bw );
                gsl_multifit_linear_est( bs_B, bs_c, bs_cov, &yi, &yierr );
                fprintf(fp, "%g %g\n", xi, yi );
            }
            fclose(fp);
            exit(0);
            */

            xi = log10( E );

            if ( (xi >= gsl_vector_get( bs_x, 0 )) && (xi <= gsl_vector_get( bs_x, n-1 )) ) {
                gsl_bspline_eval( xi, bs_B, bs_bw );
                gsl_multifit_linear_est( bs_B, bs_c, bs_cov, &yi, &yierr );
                psd  = pow( 10.0, yi );
                *dpsd = 2.303*psd*yierr;

printf("E = %g  xi = %g   yi, yierr = %g %g   psd, dpsd = %g %g\n", E, xi, yi, yierr, psd, *dpsd);
            } else {
                psd  = LGM_FILL_VALUE;
                *dpsd = LGM_FILL_VALUE;
            }

            gsl_bspline_free( bs_bw );
            gsl_vector_free( bs_B );
            gsl_vector_free( bs_x );
            gsl_vector_free( bs_y );
            gsl_matrix_free( bs_X );
            gsl_vector_free( bs_c );
            gsl_vector_free( bs_w );
            gsl_matrix_free( bs_cov );
            gsl_multifit_linear_free( bs_mw );



        } else {

            psd = LGM_FILL_VALUE;
            *dpsd = LGM_FILL_VALUE;

        }

        // spline fit to energy spectra is done.


    } else if ( f->FitType == LGM_F2P_MAXWELLIAN ) {
//Uncertainties are not propery set for thej maxwellian fits.
// Currently we use praxis -- a line minimization method, but we really should go to a differentm one
// in order to be able to get unc as well.
        if ( FitData->n > 2 ) {


            // interpolate/fit E
            // for now just do a linear interp.
            // no lets try a fit...
            double  in[10], out[7], x[10];
            in[0] = 1e-8;
            in[1] = in[2] = 1e-9; //Info->Praxis_Tolerance;
            in[5] = 30000.0; //(double)Info->Praxis_Max_Function_Evals;
            in[6] = 10.0; //Info->Praxis_Maximum_Step_Size;
            in[7] = 10.0; //Info->Praxis_Bad_Scale_Paramater;
            in[8] = 4.0; //(double)Info->Praxis_Max_Its_Without_Improvement;
            in[9] = 1.0; //(double)Info->Praxis_Ill_Conditioned_Problem;
            x[0] = 0.0;
            x[1] = -1.0;
            x[2] = 25.0;
            x[3] = -2.0;
            x[4] = 200.0;
            praxis( 2*FitData->nMaxwellians, x, (void *)FitData, Cost, in, out);
            /*
            printf("out[0] = %g\n", out[0]);
            printf("out[1] = %g\n", out[1]);
            printf("out[2] = %g\n", out[2]);
            printf("out[3] = %g\n", out[3]);
            printf("out[4] = %g\n", out[4]);
            printf("out[5] = %g\n", out[5]);
            printf("out[6] = %g\n", out[6]);
            */
            //printf("x[1] = %g   x[2] = %g   Cost = %g\n", x[1], x[2], out[6]);
 //           FILE *fp;
 //           printf("E = %g\n", E);
 //           fp = fopen("data.txt", "w");
 //           for (j=0; j<FitData->n; ++j){
 //               fprintf(fp, "%g %g\n", log10(f->E[j]), log10(FitData->g[j]));
 //           }
 //           fclose(fp);

 //           fp = fopen("fit.txt", "w");
 //           for (j=0; j<FitData->n; ++j){
 //               psd = Model( x,  FitData->nMaxwellians, f->E[j] );
 //               fprintf(fp, "%g %g\n", log10(f->E[j]), log10( psd ) );
 //           }
 //           fclose(fp);

//            exit(0);
            /*
            */


            //x[2] = 200.0;
            //printf("x - %g %g %g %g\n", x[1], x[2], x[3], x[4]);
            psd = Model( x,  FitData->nMaxwellians, E );
*dpsd = 0.2*psd;
//dpsd = ?;
            //psd = (double)a;

            //printf("E, a = %g %g  x = %g %g psd = %g\n", E, a, x[1], x[2], psd);
        } else {

            psd  = LGM_FILL_VALUE;
            *dpsd = LGM_FILL_VALUE;

        }

    } else {
        //FitType not valid...
        psd  = LGM_FILL_VALUE;
        *dpsd = LGM_FILL_VALUE;
    }

    LGM_ARRAY_1D_FREE( FitData->E );
    LGM_ARRAY_1D_FREE( FitData->g );
    LGM_ARRAY_1D_FREE( FitData->dg );
    free( FitData );


    return( psd );

}



/**
 *  Returns a pointer to a dynamically allocated Lgm_PsdToFlux structure.
 *  User must destroy this with Lgm_P2F_FreePsdToFlux() when done.
 *
 *      \param[in]      DumpDiagnostics  Boolean flag to turn on/off dumping of diagnostics.
 *      \return         A pointer to an allocated and initialized Lgm_PsdToFlux stucture.
 *
 *      \author         Mike Henderson
 *      \date           2010-2011
 *
 */
Lgm_PsdToFlux *Lgm_P2F_CreatePsdToFlux( int DumpDiagnostics ) {

    Lgm_PsdToFlux *p;

    /*
     * Allocate memory for a Lgm_PsdToFlux structure.
     */
    p = (Lgm_PsdToFlux *) calloc( 1, sizeof(*p) );

    /*
     * Set DumpDiagnostics flag to what we got here. This can be changed later as well.
     */
    p->DumpDiagnostics = DumpDiagnostics;

    p->Extrapolate  = TRUE;
    p->nMaxwellians = 2;
    p->FitType = LGM_P2F_SPLINE;

    p->Alloced1 = FALSE;
    p->Alloced2 = FALSE;


    return p;

}

/**
 * Destroy a dynamically allocated Lgm_PsdToFlux structure. (E.g. one that was
 * created by Lgm_P2F_CreatePsdToFlux().)
 *
 *      \param          p  Pointer to the allocated Lgm_PsdToFlux structure that you want to destroy.
 *
 *      \author         Mike Henderson
 *      \date           2010-2011
 *
 */
void Lgm_P2F_FreePsdToFlux( Lgm_PsdToFlux *p ) {

    if ( p->Alloced1 ) {
        LGM_ARRAY_1D_FREE( p->Mu );
        LGM_ARRAY_1D_FREE( p->K );
        LGM_ARRAY_2D_FREE( p->PSD_MK );
        LGM_ARRAY_3D_FREE( p->PSD_LMK );
    }

    if ( p->Alloced2 ) {
        LGM_ARRAY_1D_FREE( p->E );
        LGM_ARRAY_1D_FREE( p->A );
        LGM_ARRAY_1D_FREE( p->KofA );
        LGM_ARRAY_1D_FREE( p->LstarOfA );
        LGM_ARRAY_2D_FREE( p->MuofE );
        LGM_ARRAY_2D_FREE( p->PSD_EA );
        LGM_ARRAY_2D_FREE( p->FLUX_EA );
    }

    free( p );

    return;
}


/**
 *  \brief
 *      Set Date/Time and position in the Lgm_PsdToFlux structure.
 *  \details
 *
 *
 *      \param[in]      d   Date/Time of measurement.
 *      \param[in]      u   Position of measurment (in GSM).
 *      \param[in,out]  p   Lgm_PsdToFlux sturcture.
 *
 *      \author         Mike Henderson
 *      \date           2010-2011
 *
 */
void Lgm_P2F_SetDateTimeAndPos( Lgm_DateTime *d, Lgm_Vector *u, Lgm_PsdToFlux *p ) {

    p->DateTime = *d;
    p->Position = *u;

}


/**
 *  \brief
 *      Adds (to a Lgm_PsdToFlux structure) the user-supplied arrays containing PSD[Mu][K],  Mu[], K[]
 *  \detail
 *
 *      \param[in]      P                 3D array containing the Phase Space Density as a function of Mu and K.
 *      \param[in]      L                 1D array containing the Lstar values implied by the first index of PSD[][][] array.
 *      \param[in]      nL                number of Lstar values.
 *      \param[in]      Mu                1D array containing the energy values implied by the second index of PSD[][][] array.
 *      \param[in]      nMu               number of energies.
 *      \param[in]      K                 1D array containing the pitch angles values implied by the third index of PSD[][][] array.
 *      \param[in]      nK                number of pitch angles.
 *      \param[in,out]  p                 Lgm_FluxToPsd structure.
 *
 *      \author         Mike Henderson
 *      \date           2010-2011
 *
 */
void Lgm_P2F_SetPsd( double ***P, double *L, int nL, double *Mu, int nMu, double *K, int nK, Lgm_PsdToFlux *p ) {


    int     i, j, k;


    /*
     * If arrays are already alloc'd, free them first
     */
    if ( p->Alloced1 ) {
        LGM_ARRAY_1D_FREE( p->Mu );
        LGM_ARRAY_1D_FREE( p->K );
        LGM_ARRAY_1D_FREE( p->L );
        LGM_ARRAY_3D_FREE( p->PSD_LMK );
        LGM_ARRAY_2D_FREE( p->PSD_MK );
    }


    /*
     * Add Psd array to p structure. Alloc arrays appropriately.
     */
    p->nL  = nL;
    p->nMu = nMu;
    p->nK  = nK;
    LGM_ARRAY_1D( p->L, p->nL, double );
    LGM_ARRAY_1D( p->Mu, p->nMu, double );
    LGM_ARRAY_1D( p->K, p->nK, double );
    LGM_ARRAY_3D( p->PSD_LMK, p->nL, p->nMu, p->nK, double );
    LGM_ARRAY_2D( p->PSD_MK, p->nMu, p->nK, double );

    for (i=0; i<p->nL; i++)  p->L[i]  = L[i];
    for (i=0; i<p->nMu; i++) p->Mu[i] = Mu[i];
    for (i=0; i<p->nK; i++)  p->K[i]  = K[i];

    for (i=0; i<p->nL; i++) {
        for (j=0; j<p->nMu; j++) {
            for (k=0; k<p->nK; k++) {
                p->PSD_LMK[i][j][k] = P[i][j][k]; // PSD_LMK is "PSD versus L*, Mu and K".
            }
        }
    }

    p->Alloced1 = TRUE;

    return;

}


/**
 *  \brief
 *      Computes Flux at user-supplied constant values of E
 *      and \f$\alpha\f$.
 *  \details
 *      This routine ( Lgm_P2F_GetFluxAtConstEsAndAs() ) must operate on a
 *      pre-initialized Lgm_PsdToFlux structure.  The routine Lgm_P2F_SetPsd()
 *      is used to add PSD data/info to an Lgm_PsdToFlux structure.
 *    
 *      We want Flux at constant E and \f$\alpha \f$ (i.e. we want J\f$( E, \alpha
 *      )\f$).
 *    
 *      To perform the calculation we note that \f$f( \mu, K )\f$ is the same as
 *      \f$f( \mu(E, \alpha), K(\alpha) )\f$. Thus, for a given E and \f$\alpha\f$,
 *      we can figure out what \f$\mu\f$ and K  they correspond to and then we can
 *      just use the \f$f(\mu, K)\f$ array to compute the desired f values. Then
 *      finally convert f to J.
 *      The steps are;
 *    
 *          - For each \f$\alpha\f$ compute \f$K(\alpha)\f$. This is done with the routine
 *            Lgm_KofAlpha().
 *    
 *          - Then we compute \f$\mu\f$ from E and \f$\alpha\f$.
 *    
 *          - Then we look up \f$f(\mu, K)\f$ from the array (interp or fit or
 *            whatever).
 *    
 *      \param[in]      nE          Number of E values
 *      \param[in]      E           1-D array of E values
 *      \param[in]      nA          Number of Alpha values
 *      \param[in]      A           1-D array of Alpha values
 *      \param[in]      Larr        Array of precomputed L* versus Alpha values.
 *      \param[in]      Karr        Array of precomputed K versus Alpha values.
 *      \param[in]      Aarr        Array of Alpha values for Larr and Karr.
 *      \param[in]      narr        Number of values in Larr, Karr, Aarr arrays.
 *      \param[in,out]  p           A properly pre-initialized Lgm_PsdToFlux structure.
 *
 *      \author     Mike Henderson
 *      \date       2011
 *      \warning    Still working on this code. It is not finished.
 *
 */
void Lgm_P2F_GetFluxAtConstEsAndAs( double *E, int nE, double *A, int nA, double *Larr, double *Karr, double *Aarr, int narr, Lgm_MagModelInfo *mInfo, Lgm_PsdToFlux *p ) {

    int                 k, m, DoIt, i, iL, iMu, iK, done;
    double              SinAlphaEq, AlphaEq, p2c2, Lstar;
    Lgm_MagModelInfo    *mInfo2;
    Lgm_Vector          Bvec;


    /*
     * If arrays are already alloc'd, free them first
     */
    if ( p->Alloced2 ) {
        LGM_ARRAY_1D_FREE( p->E );
        LGM_ARRAY_1D_FREE( p->A );
        LGM_ARRAY_1D_FREE( p->KofA );
        LGM_ARRAY_1D_FREE( p->LstarOfA );
        LGM_ARRAY_2D_FREE( p->MuofE );
        LGM_ARRAY_2D_FREE( p->PSD_EA );
        LGM_ARRAY_2D_FREE( p->FLUX_EA );
    }

    /*
     * Alloc arrays
     */
    p->nE = nE;
    p->nA = nA;
    LGM_ARRAY_1D( p->E,     p->nE, double );
    LGM_ARRAY_1D( p->A,     p->nA,  double );
    LGM_ARRAY_1D( p->KofA,  p->nA,  double );
    LGM_ARRAY_1D( p->LstarOfA,  p->nA,  double );
    LGM_ARRAY_2D( p->MuofE, p->nE,  p->nA,  double );



    /*
     * Copy A's (given in the arguments) into p structure.
     * Transform the A's into K's using Lgm_KofAlpha().
     * Save the results in the p structure.
     * Also copy the Lstars into p structure.
     */
// If we arent using Lgm_KofAlpha(), why do we need the setup/teardown?
// This has got to be wasteful....
//    Lgm_Setup_AlphaOfK( &(p->DateTime), &(p->Position), mInfo );
    Lgm_Set_Coord_Transforms( p->DateTime.Date, p->DateTime.Time, mInfo->c );
    mInfo->Bfield( &(p->Position), &Bvec, mInfo );
    mInfo->Blocal = Lgm_Magnitude( &Bvec );
    p->B = mInfo->Blocal;
    {   // start parallel
#if USE_OPENMP
        #pragma omp parallel private(mInfo2,SinAlphaEq,AlphaEq)
        #pragma omp for schedule(dynamic, 1)
#endif
        for ( k=0; k<nA; k++ ){

            mInfo2 = Lgm_CopyMagInfo( mInfo );  // make a private (per-thread) copy of mInfo

            p->A[k]    = A[k]; // A is local Pitch Angle
            SinAlphaEq = sqrt( mInfo2->Bmin/mInfo2->Blocal ) * sin( RadPerDeg*p->A[k] );
            AlphaEq    = DegPerRad*asin( SinAlphaEq );

// REALLY SHOULD ASSUME WE HAVE THESE ALREADY. I.E. from MahEphemInfo pre-processing.
            //p->KofA[k] = Lgm_KofAlpha( AlphaEq, mInfo2 );
//            Lgm_InterpArr( Aarr, Karr, narr,   AlphaEq, &p->KofA[k] );
//            Lgm_InterpArr( Aarr, Larr, narr,   AlphaEq, &p->LstarOfA[k] );
//GUARD AGAINST BAD VALS HERE!!!? i.e. negative vals
            //printf("KofA[%d] = %g    LstarOfA[%d] = %g\n", k, p->KofA[k], k, p->LstarOfA[k]);


            Lgm_InterpArr( Aarr, Karr, narr,   A[k], &p->KofA[k] );
            Lgm_InterpArr( Aarr, Larr, narr,   A[k], &p->LstarOfA[k] );


            Lgm_FreeMagInfo( mInfo2 ); // free mInfo2

        }
    }   // end parallel
//    Lgm_TearDown_AlphaOfK( mInfo );




    /*
     * Copy E's (given in the arguments) into p structure.
     * Transform the E's into Mu's.
     * Save the results in the p structure.
     * Note that since this conversion involves E and Alpha, the result is 2D.
assumes electrons -- generalize this...
     */
    for ( m=0; m<nE; m++ ){ // energy loop
        p->E[m] = E[m];
        for ( k=0; k<nA; k++ ){ // pitch angle loop
            p->MuofE[m][k] = Lgm_Ek_to_Mu( p->E[m], p->A[k], p->B, LGM_Ee0 );
            //printf("p->E[%d], p->A[%d], p->KofA[%d], p->B, f->MuofE[%d][%d] = %g %g %g %g %g\n", m, k, k, m, k, p->E[m], p->A[k], p->KofA[k], p->B, p->MuofE[m][k]);
        }
    }





    /*
     * Now, from the PSD[Mu][K] array, get PSD at the Mu's and K's we just computed.
     * The result will be the same as PSD at the given E's and A's
     */
    LGM_ARRAY_2D( p->PSD_EA,  p->nE,  p->nA,  double );
    LGM_ARRAY_2D( p->FLUX_EA, p->nE,  p->nA,  double );
    for ( k=0; k<nA; k++ ){ // loop over pitch angles


        /*
         * Get L index and MK array
         */
        Lstar = p->LstarOfA[k];
        i=0; done = FALSE;
        while ( !done ){
            if (i >= p->nL ) {
                done = TRUE;
            } else if ( p->L[i] < Lstar ) {
                ++i;
            } else {
                done = TRUE;
            }
        }
        iL = i; if (iL<0)iL=0; if (iL>=p->nL) iL=p->nL-1;



        /*
         * Grab the MK subarray at iL (this is PSD as a func of mu and K at the Lstar implied for this PA.)
         */
        for ( iMu=0; iMu<p->nMu; iMu++ ){
            for ( iK=0; iK<p->nK; iK++ ){
                p->PSD_MK[iMu][iK] = p->PSD_LMK[iL][iMu][iK];
            }
        }



        for ( m=0; m<nE; m++ ){ // loop over energy
            p2c2 = Lgm_p2c2( p->E[m], LGM_Ee0 );

            DoIt = FALSE;

            if ( p->Extrapolate > 2 ){ // extrapolate above and below

                DoIt = TRUE;

            } else if ( p->Extrapolate == 2) {

                if (p->MuofE[m][k] >= p->Mu[0] ) DoIt = TRUE; // extrapolate above

            } else if ( p->Extrapolate == 1) {

                if (p->MuofE[m][k] <= p->Mu[p->nMu-1] ) DoIt = TRUE; // extrapolate below

            } else if ( p->Extrapolate == 0) {

                if ((p->MuofE[m][k] >= p->Mu[0])&&(p->MuofE[m][k] <= p->Mu[p->nMu-1])) DoIt = TRUE; // interp only

            }

            if (DoIt) {
                p->PSD_EA[m][k]  = Lgm_P2F_GetPsdAtMuAndK( p->MuofE[m][k], p->KofA[k], p->A[k], p );
                // Now do conversion from units of PSD to Flux
                if ( p->PSD_EA[m][k] < 0.0 ) {
//printf("p->PSD_EA[%d][%d] = %g\n", m, k, p->PSD_EA[m][k]);
                    p->FLUX_EA[m][k] = LGM_FILL_VALUE;
                } else {
                    p->FLUX_EA[m][k] = Lgm_PsdToDiffFlux( p->PSD_EA[m][k], p2c2 );
                }
//CRAP
//if (m==7)
//printf("m=%d,k=%d, p->MuofE[m][k] = %g p->KofA[k], p->A[k] = %g %g       p->PSD_EA[m][k] = %g  p->FLUX_EA[m][k] = %g\n", m,k,p->MuofE[m][k], p->KofA[k], p->A[k], p->PSD_EA[m][k], p->FLUX_EA[m][k]);
            } else {
                p->PSD_EA[m][k] = 0.0;
            }


        } // end energy loop (m index)
    } // end pitch angle loop (k index)

    if ( p->DumpDiagnostics ) {
        DumpGif( "Lgm_PsdToFlux_PSD_EA", p->nA, p->nE, p->PSD_EA );
        DumpGif( "Lgm_PsdToFlux_FLUX_EA", p->nA, p->nE, p->FLUX_EA );
    }


    p->Alloced2 = TRUE;

    return;

}



/**
 *  \brief
 *      Computes PSD at a given mu and K
 *  \details
 *     The p structure should have an initialized PSD[L][Mu][K] array in it (i.e.
 *     as added by Lgm_P2F_SetPsd()).  When transforming from PSD at constant Mu
 *     and K, to Flux at constant E and Alpha, we note that a given Alpha implies
 *     both a value for K and a value for L*. First, since the L* dimension is
 *     usually large (because these types of arrays typically come from diffusion
 *     codes), we can just use the nearest bin for iL. However, that wont
 *     necessarily work so well for the K dimension because there are usually far
 *     fewer of those bins. For K, we will do a full interpolation. After that we
 *     will have an array of PSD(Mu) at the pitch angle we are interested in. This
 *     is easily transformed into PSD(E) which we then fit with a relativistic
 *     Maxwellian and evaluate at the E we desire.
 *
 *     In this routine, we assume we have already done the iL step to reduce the 3D
 *     array down to a PSD[Mu][K] array.
 */
double  Lgm_P2F_GetPsdAtMuAndK( double Mu, double K, double A, Lgm_PsdToFlux *p ) {

    int         j, i, i0, i1, nn;
    double      K0, K1, y0, y1, slp, psd, E, g;
    _FitData    *FitData;

    // if K < 0, we should return fill value.
    if ( K < 0.0 ) return(LGM_FILL_VALUE);

    FitData = (_FitData *) calloc( 1, sizeof( _FitData ) );
    FitData->nMaxwellians = p->nMaxwellians;


    /*
     * Interpolate on K first to get a 1D array of f(mu).
     */
    if ( K > p->K[p->nK - 1] ) {
        //printf("Lgm_P2F_GetPsdAtMuAndK: (A) K >  p->K[%d] = %d %g\n", K, p->nK - 1, p->K[p->nK - 1] );
        return(LGM_FILL_VALUE);
        i0 = p->nK - 2; i1 = p->nK - 1;
    } else if ( K < p->K[0] ) {
        //printf("Lgm_P2F_GetPsdAtMuAndK: (B) K <  p->K[0] = %g %g\n", K, p->K[0] );
        return(LGM_FILL_VALUE);
        i0 = 0; i1 = 1;
    } else {
        for (i=1; i<p->nK; i++) {
            if ( K < p->K[i] ) {
                i0 = i-1; i1 = i;
                break;
            }
        }
    }

    // interpolate K
    LGM_ARRAY_1D( FitData->E,  p->nMu, double );
    LGM_ARRAY_1D( FitData->g,  p->nMu, double );
    LGM_ARRAY_1D( FitData->dg, p->nMu, double );
    FitData->n = 0;
    for (j=0; j<p->nMu; ++j){
        K0   = p->K[i0];
        K1   = p->K[i1];
        y0   = p->PSD_MK[j][i0];
        y1   = p->PSD_MK[j][i1];
        if ( (y0>0.0)&&(y1>0.0) ){
            slp  = (y1-y0)/(K1-K0);
            g = slp*(K-K1) + y1;
            if ( g > 1e-40 ) {
                FitData->g[ FitData->n ] = g;
                FitData->E[ FitData->n ] = Lgm_Mu_to_Ek( p->Mu[j], A, p->B, LGM_Ee0 );
                ++(FitData->n);
            }
        }
    }

    if ( p->FitType == LGM_P2F_SPLINE ) {

        /*
         * Use smoothing spline.
         */
        int n;
        int ncoeffs = 8;
        int nbreak  = ncoeffs - 2;

        // determine number of points
        for (n=0, j=0; j<FitData->n; ++j){
            //if ( (FitData->E[j] > 1.0/1000.0) && ( FitData->g[j] > 0.0 ) ) {
            if ( (FitData->E[j] > 0.0/1000.0) && ( FitData->g[j] > 0.0 ) ) {
                ++n;
            }
        }

        if ( n > 20 ) {
            ncoeffs = 12;
        } else if ( n > 5 ) {
            ncoeffs = 4;
        } else {
            ncoeffs = n/2;
        }
        nbreak  = ncoeffs-2;

        if ( (n >= 4) && (nbreak>=2) ) {

            // allocate a cubic bspline workspace (k = 4)
            gsl_bspline_workspace *bs_bw;
            gsl_vector            *bs_B;
            bs_bw = gsl_bspline_alloc( 4, nbreak );
            bs_B  = gsl_vector_alloc( ncoeffs );

            gsl_vector *bs_x, *bs_y, *bs_c, *bs_w;
            gsl_matrix *bs_X, *bs_cov;
            gsl_multifit_linear_workspace *bs_mw;
            double  bs_chisq, xi, yi, yierr, sigma;

            bs_x = gsl_vector_alloc( n );
            bs_y = gsl_vector_alloc( n );
            bs_X = gsl_matrix_alloc( n, ncoeffs );
            bs_c = gsl_vector_alloc( ncoeffs );
            bs_w = gsl_vector_alloc( n );
            bs_cov = gsl_matrix_alloc( ncoeffs, ncoeffs );
            bs_mw = gsl_multifit_linear_alloc( n, ncoeffs );

            // set up data arrays.
            for ( nn=0, j=0; j<FitData->n; j++ ){

                //if ( (FitData->E[j] > 1.0/1000.0) && ( FitData->g[j] > 0.0 ) ) {
                if ( (FitData->E[j] > 0.0/1000.0) && ( FitData->g[j] > 0.0 ) && ( FitData->dg[j] > 0.0 )) {
                    xi    = log10( FitData->E[j] );
                    yi    = log10( FitData->g[j] );
                    sigma = .434*( FitData->dg[j]/FitData->g[j] );

                    gsl_vector_set( bs_x, nn, xi );
                    gsl_vector_set( bs_y, nn, yi );
                    gsl_vector_set( bs_w, nn, 1.0/(sigma*sigma) );

                    ++nn;
                }

            }

            // use uniform breakpoints on defined data interval
            gsl_bspline_knots_uniform( gsl_vector_get( bs_x, 0 ), gsl_vector_get( bs_x, n-1 ), bs_bw );

            // construct the fit matrix X
            for ( i=0; i<n; i++ ) {
                xi = gsl_vector_get( bs_x, i );

                // compute B_j(xi) for all j
                gsl_bspline_eval( xi, bs_B, bs_bw );

                // fill in row i of X
                for ( j=0; j<ncoeffs; j++ ) {
                    double Bj = gsl_vector_get( bs_B, j );
                    gsl_matrix_set( bs_X, i, j, Bj );
                }
            }

            /* do the fit */
            gsl_multifit_wlinear( bs_X, bs_w, bs_y, bs_c, bs_cov, &bs_chisq, bs_mw );

            E = Lgm_Mu_to_Ek( Mu, A, p->B, LGM_Ee0 );
            xi = log10( E );

            if ( (xi >= gsl_vector_get( bs_x, 0 )) && (xi <= gsl_vector_get( bs_x, n-1 )) ) {
                gsl_bspline_eval( xi, bs_B, bs_bw );
                gsl_multifit_linear_est( bs_B, bs_c, bs_cov, &yi, &yierr );
                psd = pow( 10.0, yi );
            } else {
                psd = LGM_FILL_VALUE;
            }

            gsl_bspline_free( bs_bw );
            gsl_vector_free( bs_B );
            gsl_vector_free( bs_x );
            gsl_vector_free( bs_y );
            gsl_matrix_free( bs_X );
            gsl_vector_free( bs_c );
            gsl_vector_free( bs_w );
            gsl_matrix_free( bs_cov );
            gsl_multifit_linear_free( bs_mw );



        } else {

            psd = LGM_FILL_VALUE;

        }

    } else if ( p->FitType == LGM_P2F_MAXWELLIAN ) {

      if ( FitData->n >=  2 ) {


        // interpolate/fit E
        // for now just do a linear interp.
        // no lets try a fit...
        double  in[10], out[7], x[5];
        in[0] = 1e-8;
        in[1] = in[2] = 1e-9; //Info->Praxis_Tolerance;
        in[5] = 30000.0; //(double)Info->Praxis_Max_Function_Evals;
        in[6] = 10.0; //Info->Praxis_Maximum_Step_Size;
        in[7] = 10.0; //Info->Praxis_Bad_Scale_Paramater;
        in[8] = 4.0; //(double)Info->Praxis_Max_Its_Without_Improvement;
        in[9] = 1.0; //(double)Info->Praxis_Ill_Conditioned_Problem;
        x[0] = 0.0;
        x[1] = -1.0;
        x[2] = 25.0;
        x[3] = -2.0;
        x[4] = 200.0;
        praxis( 2*FitData->nMaxwellians, x, (void *)FitData, Cost, in, out);

        //printf("PSD_TO_FLUX: x - %g %g %g %g\n", x[1], x[2], x[3], x[4]);
        E = Lgm_Mu_to_Ek( Mu, A, p->B, LGM_Ee0 );
        psd = Model( x,  FitData->nMaxwellians, E );
        //psd = (double)a;

        //printf("E, A = %g %g  x = %g %g psd = %g\n", E, A, x[1], x[2], psd);
      } else {

        psd = LGM_FILL_VALUE;

      }

    } else {
      //FitType not valid...
      psd = LGM_FILL_VALUE;
    }



    LGM_ARRAY_1D_FREE( FitData->E );
    LGM_ARRAY_1D_FREE( FitData->g );
    LGM_ARRAY_1D_FREE( FitData->dg );
    free( FitData );


    return( psd );

}













/**
 * Routine to increase the size of an image smoothly. Uses an area-weighting
 * interpolation scheme.
 */
void UpSizeImage( double **Orig, double *X, double *Y, int M, int N, double **New, double *U, double *V, int K, int J ) {

    double  DX, DY, DXO2, DYO2, A;
    double  x, y, xe, ye, xx, yy, xp, yp;
    double  f[4], w0, w1, w2, w3;
    double  dxm, dxp, dym, dyp;
    int     j, k, np, mp, ne, me;
    int     n0, n1, n2, n3;
    int     m0, m1, m2, m3;



    /*
     * compute size and area of a pixel in orig image.
     */
    DX = 1.0/(double)N;
    DY = 1.0/(double)M;
    DXO2 = 0.5*DX;
    DYO2 = 0.5*DY;
    A    = DX*DY;
    printf("A = %g\n", A);


    for (j=0; j<J; j++ ){
        for (k=0; k<K; k++ ){


            /*
             *  Compute the normalized coords (x,y) for the center of this
             *  pixel.
             */
            x = ((double)j+0.5)/(double)J;
            y = ((double)k+0.5)/(double)K;



            /*
             *  Determine which pixel of the original image we are in.
             */
            np = (int)(x*N);
            mp = (int)(y*M);
            xp = ((double)np+0.5)/(double)N;
            yp = ((double)mp+0.5)/(double)M;



            /*
             *  Find the edges in the original image that this pixel is closest to.
             */
            ne = (int)(x*N + 0.5);
            me = (int)(y*M + 0.5);
            xe = (double)ne/(double)N;
            ye = (double)me/(double)M;




            /*
             * Compute coords of the four closest pixels
             */
            //if        ( ( np < ne ) && ( mp < me ) ){  // we are in upper left  Pix #0
            if        ( ( x > xp ) && ( y > yp ) ){  // we are in upper left  Pix #0
                n0 = np;   m0 = mp;
                n1 = np+1; m1 = mp;
                n2 = np+1; m2 = mp+1;
                n3 = np;   m3 = mp+1;
            //} else if ( ( np > ne ) && ( mp < me ) ){  // we are in upper right Pix #1
            } else if ( ( x < xp ) && ( y > yp ) ){  // we are in upper right Pix #1
                n0 = np-1; m0 = mp;
                n1 = np;   m1 = mp;
                n2 = np;   m2 = mp+1;
                n3 = np-1; m3 = mp+1;
            //} else if ( ( np > ne ) && ( mp > me ) ){  // we are in lower right Pix #2
            } else if ( ( x < xp ) && ( y < yp ) ){  // we are in lower right Pix #2
                n0 = np-1; m0 = mp-1;
                n1 = np;   m1 = mp-1;
                n2 = np;   m2 = mp;
                n3 = np-1; m3 = mp;
            } else {                                   // we are in lower left  Pix #3
                n0 = np;   m0 = mp-1;
                n1 = np+1; m1 = mp-1;
                n2 = np+1; m2 = mp;
                n3 = np;   m3 = mp;
            }

  xe = 0.5*(n0+n1+1.0)/(double)N;
  ye = 0.5*(m1+m2+1.0)/(double)M;



            /*
             *  Check for edge problems
             */
            if (n0 < 0)   n0 = 0;
            if (n3 < 0)   n3 = 0;
            if (n1 > N-1) n1 = N-1;
            if (n2 > N-1) n2 = N-1;

            if (m0 < 0)     m0 = 0;
            if (m1 < 0)     m1 = 0;
            if (m2 > M-1)   m2 = M-1;
            if (m3 > M-1)   m3 = M-1;



            /*
             * Extract the pixel vals.
             */
            f[0] = Orig[ m0 ][ n0 ];
            f[1] = Orig[ m1 ][ n1 ];
            f[2] = Orig[ m2 ][ n2 ];
            f[3] = Orig[ m3 ][ n3 ];



            /*
             *  Compute (x',y') (i.e. the coords relative to center of 4 pixels.)
             */
            xx = x - xe;
            yy = y - ye;



            /*
             *  Compute weights (i.e. overlapping areas)
             */
            dxm = DXO2 - xx;
            dxp = DXO2 + xx;
            dym = DYO2 - yy;
            dyp = DYO2 + yy;
//if (dxm < 0.0) { printf("dxm = %g\n", dxm); exit(0); }
//if (dxp < 0.0) { printf("dxp = %g\n", dxp); exit(0); }
//if (dym < 0.0) { printf("dym = %g    DYO2 = %g  x, y = %g %g   xx, yy = %g %g\n", dym, DYO2, x, y, xx, yy); exit(0); }
//if (dyp < 0.0) { printf("dyp = %g\n", dyp); exit(0); }

            w0 = dxm*dym;
            w1 = dxp*dym;
            w2 = dxp*dyp;
            w3 = dxm*dyp;
            /*
            w0=0.1;
            w1=0.2;
            w2=0.3;
            w3=0.4;
            */

if ((j==189)&&(k==118)) {
                printf("np,mp = %d %d    pixels:  %d %d   %d %d   %d %d   %d %d\n", np, mp, n0, m0, n1, m1, n2, m2, n3, m3);
                printf("j,k = %d %d     f = %g %g %g %g     w = %g %g %g %g     xp, yp = %g %g   x,y = %g %g    xx, yy = %g %g\n", j, k, f[0], f[1], f[2], f[3], w0, w1, w2, w3, xp, yp, x, y, xx, yy);
}

            /*
             *  Finally compute weighted average and assign to current pixel of new image.
             *  A is the total area of an original pixel. w's are the partial overlap areas.
             */
            New[k][j] = (f[0]*w0 + f[1]*w1 + f[2]*w2 + f[3]*w3)/A;
            //New[k][j] = (f[0]*w0 + f[1]*w1 )/A;
            //New[k][j] = (f[0]*w0 )/A;
            //New[k][j] = (double)k;


        }
    }






}

/**
 *   Routine to write out a GIF image
 */
void DumpGif( char *FilenameBase, int W, int H, double **Image ){

    double           Val, Val2, Min, Max, dVal;
    int              w, h;
    unsigned char   *uImage, uVal, Filename[1024];
    FILE            *fp_gif, *fp_info;

    int             LogScale;

    LogScale = FALSE;
    LogScale = TRUE;


    // Determine Min/Max values...
    Min =  9e99;
    Max = LGM_FILL_VALUE;
    for ( w=0; w<W; w++ ){
        for ( h=0; h<H; h++ ) {

            Val = Image[h][w];
            if ( LogScale ) {
                Val2 = (Val > 0.0) ? log10( Val ) : LGM_FILL_VALUE;
            } else {
                Val2 = Image[h][w];
            }
            if (Val2 > Max) Max = Val2;
            if ((Val2 < Min)&&(Val > 0.0)) Min = Val2;

        }
    }

    //printf("Min, Max = %g %g\n", Min, Max);

    sprintf( Filename, "%s.info", FilenameBase);
    fp_info = fopen( Filename, "w" );
    fprintf( fp_info, "Min: %g\n", Min );
    fprintf( fp_info, "Max: %g\n", Max );
    fclose( fp_info );



    uImage = (unsigned char *)calloc( W*H, sizeof(unsigned char) );

    for ( w=0; w<W; w++ ){
        for ( h=0; h<H; h++ ) {

            if ( LogScale ) {
                Val = Image[h][w] > 0.0 ? log10( Image[h][w] ) : LGM_FILL_VALUE;
            } else {
                Val = Image[h][w];
            }
            if ( Val <= LGM_FILL_VALUE ) {
                uVal = 0;
            } else {
                dVal = (Val - Min)/(Max-Min)*255.0;
                uVal = (dVal > 255.0) ? 255 : (unsigned char)dVal;
            }
            *(uImage + W*(H-1-h) + w) = uVal;

        }
    }

    sprintf( Filename, "%s.gif", FilenameBase);
    fp_gif = fopen(Filename, "w");
    WriteGIF(fp_gif, uImage, 0, W, H, Rainbow3_Red, Rainbow3_Grn, Rainbow3_Blu, 256, 0, "");
    fclose(fp_gif);

    free( uImage );



    // dump a colorbar image
    W = 10; H = 256;
    uImage = (unsigned char *)calloc( W*H, sizeof(unsigned char) );
    for ( w=0; w<W; w++ ){
        for ( h=0; h<H; h++ ) {
            *(uImage + W*(H-1-h) + w) = h;
        }
    }
    sprintf( Filename, "%s_Bar.gif", FilenameBase);
    fp_gif = fopen(Filename, "w");
    WriteGIF(fp_gif, uImage, 0, W, H, Rainbow3_Red, Rainbow3_Grn, Rainbow3_Blu, 256, 0, "");
    fclose(fp_gif);
    free( uImage );

}

/**
 *   Routine to write out a GIF image
 */
void DumpGif2( char *FilenameBase, double Min, double Max, int W, int H, double **Image ){

    double           Val, dVal;
    int              w, h;
    unsigned char   *uImage, uVal, Filename[1024];
    FILE            *fp_gif, *fp_info;

    int             LogScale;

    LogScale = TRUE;
    LogScale = FALSE;


    //printf("Min, Max = %g %g\n", Min, Max);

    sprintf( Filename, "%s.info", FilenameBase);
    fp_info = fopen( Filename, "w" );
    fprintf( fp_info, "Min: %g\n", Min );
    fprintf( fp_info, "Max: %g\n", Max );
    fclose( fp_info );



    uImage = (unsigned char *)calloc( W*H, sizeof(unsigned char) );

    for ( w=0; w<W; w++ ){
        for ( h=0; h<H; h++ ) {

            if ( LogScale ) {
                Val = Image[h][w] > 0.0 ? log10( Image[h][w] ) : LGM_FILL_VALUE;
            } else {
                Val = Image[h][w];
            }
            if ( Val <= LGM_FILL_VALUE  ) {
                uVal = 0;
            } else {
                dVal = (Val - Min)/(Max-Min)*255.0;
                if (dVal > 255.0) {
                    uVal = 255;
                } else if ( dVal < 0.0 ) {
                    uVal = 0;
                } else {
                    uVal = (unsigned char)dVal;
                }
            }
            *(uImage + W*(H-1-h) + w) = uVal;

        }
    }

    sprintf( Filename, "%s.gif", FilenameBase);
    fp_gif = fopen(Filename, "w");
    WriteGIF(fp_gif, uImage, 0, W, H, Rainbow3_Red, Rainbow3_Grn, Rainbow3_Blu, 256, 0, "");
    fclose(fp_gif);

    free( uImage );



    // dump a colorbar image
    W = 10; H = 256;
    uImage = (unsigned char *)calloc( W*H, sizeof(unsigned char) );
    for ( w=0; w<W; w++ ){
        for ( h=0; h<H; h++ ) {
            *(uImage + W*(H-1-h) + w) = h;
        }
    }
    sprintf( Filename, "%s_Bar.gif", FilenameBase);
    fp_gif = fopen(Filename, "w");
    WriteGIF(fp_gif, uImage, 0, W, H, Rainbow3_Red, Rainbow3_Grn, Rainbow3_Blu, 256, 0, "");
    fclose(fp_gif);
    free( uImage );

}



/**
 *  \brief
 *     Fills the given array with a geomteric sequence of values.
 *
 *  \details
 *     Blah blah
 *
 *      \param[in]      a  Starting value.
 *      \param[in]      b  Ending value.
 *      \param[in]      n  Number of values (must be 2 or more.)
 *      \param[out]     G  User supply array (or size >= n ).
 *
 *      \return         1 if successful, 0 otherwise.
 *
 *  \notes
 *      a and be must be greater than zero and n must be greater than or equal to 2.
 *
 *      \author         Mike Henderson
 *      \date           2010-2011
 */
int Lgm_GeometricSeq( double a, double b, int n, double *G ) {

    if ( ( a <= 0.0 ) || ( b <= 0.0 ) || ( n<2 ) ) return(0);

    int     j;
    double  r, *S = (double *)calloc( n, sizeof(double) );

    r = pow( b/a, 1.0/((double)(n-1)));


    S[0] = a;
    for (j=1; j<n; j++) S[j] = S[j-1]*r;
    for (j=0; j<n; j++) G[j] = S[j];

    free( S );

    return(1);
}



int Lgm_InterpArr( double *xa, double *ya, int n, double x, double *y ) {

    gsl_interp_accel    *acc;
    gsl_spline          *spline;
    double              *xa2, *ya2;
    int                 i, n2, RetVal;

    xa2 = (double *)calloc( n, sizeof(double) );
    ya2 = (double *)calloc( n, sizeof(double) );
    
    /*
     *  Put array into ascending order if it isnt already. (gsl needs this).
     */
    if ( xa[1] < xa[0] ) {

        for (n2=0, i=0; i<n; i++){
            if ( ya[n-1-i] > LGM_FILL_VALUE ) {
                xa2[n2] = xa[n-1-i];
                ya2[n2] = ya[n-1-i];
                ++n2;
            }
        }

    } else {
        for (n2=0, i=0; i<n; i++){
            if ( ya[i] > LGM_FILL_VALUE ) {
                xa2[n2] = xa[i];
                ya2[n2] = ya[i];
                ++n2;
            }
        }
    }


    /*
     * Check to see if x would cause an extrapolation instead of an interp.
     */
    if ( (x<xa2[0]) || (x>xa2[n2-1]) ){
//printf("here:  x = %lf     xa2[0] = %lf    n2 = %d  xa2[n2-1] = %lf\n", x, xa2[0], n2, xa2[n2-1]);
        *y = LGM_FILL_VALUE;
        RetVal = -1;
    } else if ( n2 > 4 ) { //think Akima spline needs at least 5 points
        acc    = gsl_interp_accel_alloc( );
        spline = gsl_spline_alloc( gsl_interp_akima, n2 );
        gsl_spline_init( spline, xa2, ya2, n2 );
        *y = gsl_spline_eval( spline, x, acc );
        gsl_spline_free( spline );
        gsl_interp_accel_free( acc );
        RetVal = 1;
    } else if ( n2 > 1) { //fall back to linear
        acc    = gsl_interp_accel_alloc( );
        spline = gsl_spline_alloc( gsl_interp_linear, n2 );
        gsl_spline_init( spline, xa2, ya2, n2 );
        *y = gsl_spline_eval( spline, x, acc );
        gsl_spline_free( spline );
        gsl_interp_accel_free( acc );
        RetVal = 1;
    } else {
       *y = LGM_FILL_VALUE;
        RetVal = -1; 
    }


    /*
     * Free new vars
     */
    free( xa2 );
    free( ya2 );

    return( RetVal );

}


