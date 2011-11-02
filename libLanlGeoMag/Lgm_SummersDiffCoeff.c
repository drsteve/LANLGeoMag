#include "Lgm/Lgm_SummersDiffCoeff.h"

#include <gsl/gsl_poly.h>


/**
 *  \brief
 *      Compute the electron plasma frequency, \f$\omega_{pe}\f$.
 *
 *  \details
 *
 *      \param[in]      Density     Electron number density in units of cm^-3.
 *
 *      \return         Electron plasma frequency in units of s^-1
 *
 *      \author         M. Henderson
 *      \date           2011
 *
 *
 */
double  Lgm_ePlasmaFreq( double Density ) {

    /*
     *  In SI units,
     *      wpe = sqrt( (Density*LGM_e*LGM_e)/(LGM_me*LGM_Eps0) );
     *
     *  Putting Density in units of #/cm^2, we have,
     *      wpe = 1000 * sqrt( (Density*LGM_e*LGM_e)/(LGM_me*LGM_Eps0) );
     *      wpe = 56414.601925 * sqrt( Density );
     */

    return( 56414.601925*sqrt( Density ) );

}

/**
 *  \brief
 *      Compute the gyro-frequency of a particle with change q and mass m, \f$\Omega_\sigma\f$.
 *
 *  \details
 *
 *      \param[in]      Charge      Electric charge of the particle. in Coulombs
 *      \param[in]      B           Magnetic field in units of nT.
 *      \param[in]      mass        Mass of the particle in units of kg.
 *
 *      \return         Gyro-frequency of particle in units of s^-1
 *
 *      \author         M. Henderson
 *      \date           2011
 *
 *
 */
double  Lgm_GyroFreq( double q, double B, double m ) {

    /*
     *  In SI units,
     *      Omega = (q*B)/(m*c)
     *
     *  Putting B in units of nT, we have,
     *      Omega = 1e-9 * q*B/(m*c)
     */


    return( 1e-9*B*q/m );

}


/**
 *  \brief
 *      Computes the bounce-averaged Summer's [2005] diffusion coefficients
 *      in a pure centered dipole field.
 *
 *  \details
 *      The bounce average of a quantity \f$Q\f$ is;
 *
 *          \f[ <Q> = { 1\over S_b} \int_{s_{sm}}^{s_{nm}} {Q\;ds\over [1-{B/B_m}]^{1/2}} \f]
 *
 *      where,
 *          \f[ S_b = \int_{s_{sm}}^{s_{nm}} {ds\over [1-{B/B_m}]^{1/2}} \f]
 *
 *      and \f$s\f$ is the distance along the field line and \f$s_{sm}\f$ and
 *      \f$s_{nm}\f$ are the southern and northern mirror points respectively.
 *
 *      In a pure dipole field, the integrals can be recast as integrals over
 *      latitude using;
 *
 *          \f[ ds = L\cos\lambda[4-3\cos^2\lambda]^{1/2} d\lambda \f]
 *          \f[ B/B_m = { \sin^2\alpha_\circ [4-3\cos^2\lambda]^{1/2}\over cos^6\lambda } \f]
 *
 *      Since \f$L\f$ is constant and is in the numerator and denominator it
 *      cancels (so we can ignore it).
 *
 *      The final result is;
 *
 *          \f[ <Q> = { 1\over S_b^\prime}\displaystyle
 *          \int_{\lambda_{sm}}^{\lambda_{nm}} {Q\;
 *          \cos\lambda[4-3\cos^2\lambda]^{1/2} d\lambda  \over \left[1-
 *          {\displaystyle \sin^2\alpha_\circ
 *          [4-3\cos^2\lambda]^{1/2}\over\displaystyle cos^6\lambda }
 *          \right]^{1/2}} \f]
 *
 *      If we already have B/Bm, then its easier to compute this as;
 *
 *          \f[ <Q> = { 1\over S_b}\displaystyle
 *          \int_{\lambda_{sm}}^{\lambda_{nm}} {Q\;
 *          \cos\lambda[4-3\cos^2\lambda]^{1/2} d\lambda  \over \left[1-
 *          B/B_m\right]^{1/2}} \f]
 *
 *      where,
 *
 *          \f[ S_b^\prime = \displaystyle \int_{\lambda_{sm}}^{\lambda_{nm}}
 *          {\cos\lambda[4-3\cos^2\lambda]^{1/2} d\lambda  \over \left[1-
 *          B/B_m\right]^{1/2}} \f]
 *
 *      Note that \f$ S_b^\prime = S_b/L \f$, is a normalized \f$ S_b \f$ .
 *
 *      This routine computes \f$ <D_{\alpha_\circ\alpha_\circ}>,
 *      <D_{\alpha_\circ p}>, <D_{pp}> \f$.
 *
 *
 *      \param[in]      Version     Version of Summers Model to use. Can be LGM_SUMMERS_2005 or LGM_SUMMERS_2007.
 *      \param[in]      Alpha0      Equatoria PA in Degrees.
 *      \param[in]      Ek          Kinetic energy in MeV.
 *      \param[in]      L           L-shell parameter (dimensionless).
 *      \param[in]      BwFuncData  A (void *) pointer to data that the user may need to use in BwFunc().
 *      \param[in]      BwFunc      A user-defined function that returns Bw as a function of Latitude in units of nT. Prototype is "double BwFunc( double Lat, void *Data )".
 *      \param[in]      n1          Ratio of Hydrogen number density to Electron number density (note that n1+n2+n3=1).
 *      \param[in]      n2          Ratio of Helium number density to Electron number density (note that n1+n2+n3=1).
 *      \param[in]      n3          Ratio of Oxygen number density to Electron number density (note that n1+n2+n3=1).
 *      \param[in]      aStarEq     Equatorial value of the cold plasma parameter \f$\Omega_\sigma/\omega_{pe}\f$.
 *      \param[in]      w1          Lower limit of freq. band, [Hz].
 *      \param[in]      w2          Lower limit of freq. band, [Hz].
 *      \param[in]      wm          Midpoint of freq. band, [Hz].
 *      \param[in]      dw          Width of freq. band, [Hz].
 *      \param[in]      WaveMode    Mode of the wave. Can be LGM_R_MODE_WAVE or LGM_L_MODE_WAVE
 *      \param[in]      Species     Particle species. Can be LGM_ELECTRONS, LGM_PROTONS, or ....?
 *      \param[in]      MaxWaveLat  Latitude (+/-) that waves exist up to. [Degrees].
 *      \param[out]     Daa_ba      Bounce-averaged value of Daa.
 *      \param[out]     Dap_ba      Bounce-averaged value of Dap.
 *      \param[out]     Dpp_ba      Bounce-averaged value of Dpp.
 *
 *      \return         0 if successful. <0 otherwise.
 *
 *      \author         M. Henderson
 *      \date           2011
 *
 */
int Lgm_SummersDxxBounceAvg( int Version, double Alpha0,  double Ek,  double L,  void *BwFuncData, double (*BwFunc)(), double n1, double n2, double n3, double aStarEq,  double w1, double w2, double wm, double dw, int WaveMode, int Species, double MaxWaveLat, double *Daa_ba,  double *Dap_ba,  double *Dpp_ba) {

    double           T, a, b, E0, Omega_eEq, Omega_SigEq, Beq, Rho;
    double           epsabs, epsrel, abserr, work[2002], points[10];
    double           ySing, a_new, b_new, B;
    int              npts2=4, limit=500, lenw=4*limit, iwork[502], last, ier, neval;
    Lgm_SummersInfo *si=(Lgm_SummersInfo *)calloc( 1, sizeof(Lgm_SummersInfo));

    si->Version = Version;
    si->n1 = n1;
    si->n2 = n2;
    si->n3 = n3;
    if ( Version == LGM_SUMMERS_2007 ){
        if ( fabs(1.0-(n1+n2+n3)) > 1e-10 ) {
            printf("Lgm_SummersDxxBounceAvg: n1+n2+n3 is not 1. Got n1+n2+n3 = %g\n", n1+n2+n3 );
            return(-1);
        }
    }




    if ( (fabs(Alpha0-90.0) < 1e-8)||(fabs(Alpha0) < 1e-8) ) {

        *Daa_ba = 0.0;
        *Dap_ba = 0.0;
        *Dpp_ba = 0.0;
        return( 0 );

    }

    /*
     *  Set Rho. See Eqn(3) of Summers2007.
     *      Rho = sqrt(PI)/2.0( erf( (wm-w1)/dw ) + erf( (w2-wm)/dw ) )
     *  Note: If wm-w1 = w2-wm, then this is the same as the Summers2005 paper.
     */
    Rho = ( erf( (wm-w1)/dw ) + erf( (w2-wm)/dw ) )/M_2_SQRTPI;

    /*
     * convert w1, w2, wm, dw to radians/second
     */

    w1 = w1*M_2PI ;
    w2 = w2*M_2PI ;
    wm = wm*M_2PI ;
    dw = dw*M_2PI ;

    /*
     *  Quadpack integration tolerances.
     */
    epsabs = 0.0;
    epsrel = 1e-5;


    /*
     *  Pack integrand parameters into structure
     */
    if ( WaveMode == LGM_R_MODE_WAVE ) {
        si->s = 1;
    } else if ( WaveMode == LGM_L_MODE_WAVE ) {
        si->s = -1;
    } else {
        printf("Lgm_SummersDxxBounceAvg(): Unknown wavemode = %d\n", WaveMode );
        return(0);
    }


    Beq = M_CDIP/(L*L*L);
    if ( Species == LGM_ELECTRONS ) {
        si->Lambda  = -1.0;
        E0          = LGM_Ee0;      // set rest energy to electron rest energy (MeV)
        Omega_SigEq = Lgm_GyroFreq( -LGM_e, Beq, LGM_ELECTRON_MASS );
    } else if ( Species == LGM_PROTONS ) {
        si->Lambda = LGM_EPS;
        E0         = LGM_Ep0;       // set rest energy to proton rest energy (MeV)
        Omega_SigEq = Lgm_GyroFreq( LGM_e, Beq, LGM_PROTON_MASS );
    } else {
        printf("Lgm_SummersDxxBounceAvg(): Unknown species = %d\n", Species );
        return(0);
    }
    Omega_eEq = Lgm_GyroFreq( -LGM_e, Beq, LGM_ELECTRON_MASS );
    si->Omega_eEq   = Omega_eEq;    // Equatorial electron gyro-frequency.
    si->Omega_SigEq = Omega_SigEq;  // Equatorial gyro-frequency for given species.

    si->E           = Ek/E0;        // Dimensionless energy Ek/E0 (kinetic energy over rest energy).
    si->Alpha0      = Alpha0*M_PI/180.0;
    si->SinAlpha0   = sin(si->Alpha0);
    si->SinAlpha02  = si->SinAlpha0*si->SinAlpha0;
    si->CosAlpha0   = cos(si->Alpha0);
    si->CosAlpha02  = si->CosAlpha0*si->CosAlpha0;
    si->TanAlpha0   = si->SinAlpha0/si->CosAlpha0;
    si->TanAlpha02  = si->TanAlpha0*si->TanAlpha0;
    si->L           = L;
    si->aStarEq     = aStarEq;
//    si->dB          = dB;
    si->BwFuncData  = BwFuncData;
    si->BwFunc      = BwFunc;
    si->w1          = w1;           // Lower freq cuttof.
    si->w2          = w2;           // Upper freq cuttof.
    si->wm          = wm;           // Frequency of max wave power.
    si->dw          = dw;           // Frequency bandwidth (at equator).
    si->MaxWaveLat  = MaxWaveLat*RadPerDeg;   // Assume there are no waves at +/-MaxWaveLat.
    si->Rho         = Rho;


    /*
     *  Set integration limits
     */
    a = 0.0;                                            // radians
    B = acos( Lgm_CdipMirrorLat( si->SinAlpha0 ) );     // radians
    b = ( B < si->MaxWaveLat ) ? B : si->MaxWaveLat;
    if ( fabs(a-b) < 1e-9 ) {

        *Daa_ba = 0.0;
        *Dap_ba = 0.0;
        *Dpp_ba = 0.0;
        return( 0 );

    } else {


        /*
         *  Perform integrations. The points[] array contains a list of point
         *  in the integrand where integrable singularities occur. dqagp() uses
         *  these to break up the integral at these points (dqagp() is an
         *  extrapolation algorithm so it handles singularities very well.)
         *  Note: the value of npts2 must be 2+ the number of points added.
         *  This is because it adds the endpoints internally (so dont add those
         *  here).
         *
         *  An approximate form of T( SinAlpha0 ), related to bounce period is given in
         *  Schultz and Lanzeroti (Eqns 1.28a-d). Here we do it exactly.
         *      T0 = 1.38017299815047317375;
         *      T1 = 0.74048048969306104115;
         *      T  = T0 - 0.5*(T0-T1)*(si->SinAlpha0 + sqrt(si->SinAlpha0));
         */
        npts2 = 2;
        //dqagp( CdipIntegrand_Sb, (_qpInfo *)si, a, b, npts2, points, epsabs, epsrel, &T, &abserr, &neval, &ier, limit, lenw, &last, iwork, work );
        dqags( CdipIntegrand_Sb, (_qpInfo *)si, a, B, epsabs, epsrel, &T, &abserr, &neval, &ier, limit, lenw, &last, iwork, work );


        Lgm_SummersFindCutoffs2( SummersIntegrand_Gaa, (_qpInfo *)si, TRUE, a, b, &a_new, &b_new );
        npts2 = 2 + Lgm_SummersFindSingularities( SummersIntegrand_Gaa, (_qpInfo *)si, TRUE, a_new, b_new, &points[1], &ySing );
        if ( b_new > a_new ) {
            if ( npts2 > 2 ) {
                dqagp( SummersIntegrand_Gaa, (_qpInfo *)si, a_new, b_new, npts2, points, epsabs, epsrel, Daa_ba, &abserr, &neval, &ier, limit, lenw, &last, iwork, work );
            } else {
                dqags( SummersIntegrand_Gaa, (_qpInfo *)si, a_new, b_new, epsabs, epsrel, Daa_ba, &abserr, &neval, &ier, limit, lenw, &last, iwork, work );
            }
        } else {
            *Daa_ba = 0.0;
        }


        Lgm_SummersFindCutoffs2( SummersIntegrand_Gap, (_qpInfo *)si, TRUE, a, b, &a_new, &b_new );
        npts2 = 2; npts2 += Lgm_SummersFindSingularities( SummersIntegrand_Gap, (_qpInfo *)si, TRUE, a_new, b_new, &points[1], &ySing );
        if ( b_new > a_new ) {
            if ( npts2 > 2 ) {
                dqagp( SummersIntegrand_Gap, (_qpInfo *)si, a_new, b_new, npts2, points, epsabs, epsrel, Dap_ba, &abserr, &neval, &ier, limit, lenw, &last, iwork, work );
            } else {
                dqags( SummersIntegrand_Gap, (_qpInfo *)si, a_new, b_new, epsabs, epsrel, Dap_ba, &abserr, &neval, &ier, limit, lenw, &last, iwork, work );
            }
        } else {
            *Dap_ba = 0.0;
        }

        Lgm_SummersFindCutoffs2( SummersIntegrand_Gpp, (_qpInfo *)si, TRUE, a, b, &a_new, &b_new );
        npts2 = 2 + Lgm_SummersFindSingularities( SummersIntegrand_Gpp, (_qpInfo *)si, TRUE, a_new, b_new, &points[1], &ySing );
        if ( b_new > a_new ) {
            if ( npts2 > 2 ) {
                dqagp( SummersIntegrand_Gpp, (_qpInfo *)si, a_new, b_new, npts2, points, epsabs, epsrel, Dpp_ba, &abserr, &neval, &ier, limit, lenw, &last, iwork, work );
            } else {
                dqags( SummersIntegrand_Gpp, (_qpInfo *)si, a_new, b_new, epsabs, epsrel, Dpp_ba, &abserr, &neval, &ier, limit, lenw, &last, iwork, work );
            }
        } else {
            *Dpp_ba = 0.0;
        }

        *Daa_ba /= T;
        *Dap_ba /= T;
        *Dpp_ba /= T;

    }


    return(1);

}

/**
 *  \brief
 *      Computes derivatives of the bounce-averaged Summer's [2005, 2007] diffusion coefficients
 *      in a pure centered dipole field.
 *
 *  \details
 *
 *
 *      \param[in]      Version     Version of Summers Model to use. Can be LGM_SUMMERS_2005 or LGM_SUMMERS_2007.
 *      \param[in]      Alpha0      Equatoria PA in Degrees.
 *      \param[in]      Ek          Kinetic energy in MeV.
 *      \param[in]      L           L-shell parameter (dimensionless).
 *      \param[in]      BwFuncData  A (void *) pointer to data that the user may need to use in BwFunc().
 *      \param[in]      BwFunc      A user-defined function that returns Bw as a function of Latitude in units of nT. Prototype is "double BwFunc( double Lat, void *Data )".
 *      \param[in]      n1          Ratio of Hydrogen number density to Electron number density (note that n1+n2+n3=1).
 *      \param[in]      n2          Ratio of Helium number density to Electron number density (note that n1+n2+n3=1).
 *      \param[in]      n3          Ratio of Oxygen number density to Electron number density (note that n1+n2+n3=1).
 *      \param[in]      aStarEq     Equatorial value of the cold plasma parameter \f$\Omega_\sigma/\omega_{pe}\f$.
 *      \param[in]      w1          Lower limit of freq. band, [Hz].
 *      \param[in]      w2          Lower limit of freq. band, [Hz].
 *      \param[in]      wm          Midpoint of freq. band, [Hz].
 *      \param[in]      dw          Width of freq. band, [Hz].
 *      \param[in]      WaveMode    Mode of the wave. Can be LGM_R_MODE_WAVE or LGM_L_MODE_WAVE
 *      \param[in]      Species     Particle species. Can be LGM_ELECTRONS, LGM_PROTONS, or ....?
 *      \param[in]      MaxWaveLat  Latitude (+/-) that waves exist up to. [Degrees].
 *      \param[out]     Daa_ba      Bounce-averaged value of Daa.
 *      \param[out]     Dap_ba      Bounce-averaged value of Dap.
 *      \param[out]     Dpp_ba      Bounce-averaged value of Dpp.
 *
 *      \return         0 if successful. <0 otherwise.
 *
 *      \author         M. Henderson
 *      \date           2011
 *
 */
int Lgm_SummersDxxDerivsBounceAvg( int DerivScheme, double ha, int Version, double Alpha0,  double Ek,  double L,  void *BwFuncData, double (*BwFunc)(), double n1, double n2, double n3, double aStarEq,  double w1, double w2, double wm, double dw, int WaveMode, int Species, double MaxWaveLat, double *dDaa,  double *dDap) {

    double  a, h, H, faa[7], fap[7], Daa_ba, Dap_ba, Dpp_ba;
    int     i, N;

    /*
     * Select the derivative scheme to use.
     *      f_0^(1) = 1/(60h)  ( f_3 - 9f_2 + 45f_1  - 45f_-1 + 9f_-2 - f_-3 )
     * See page 450 of CRC standard Math tables 28th edition.
     */
    switch ( DerivScheme ) {
        case LGM_DERIV_SIX_POINT:
            N = 3;
            break;
        case LGM_DERIV_FOUR_POINT:
            N = 2;
            break;
        case LGM_DERIV_TWO_POINT:
            N = 1;
            break;
    }

    /*
     * Compute points in alpha grid use ha as spacing
     */
    h = ha;
    for (i=-N; i<=N; ++i){
        a = Alpha0; H = (double)i*h; a += H;
        Lgm_SummersDxxBounceAvg( Version, a, Ek, L, BwFuncData, BwFunc, n1, n2, n3, aStarEq, w1, w2, wm, dw, WaveMode, Species, MaxWaveLat, &Daa_ba, &Dap_ba, &Dpp_ba );
        faa[i+N] = Daa_ba;
        fap[i+N] = Dap_ba;
    }

    if (DerivScheme == LGM_DERIV_SIX_POINT){
        *dDaa = (faa[6] - 9.0*faa[5] + 45.0*faa[4] - 45.0*faa[2] + 9.0*faa[1] - faa[0])/(60.0*h);
        *dDap = (fap[6] - 9.0*fap[5] + 45.0*fap[4] - 45.0*fap[2] + 9.0*fap[1] - fap[0])/(60.0*h);
    } else if (DerivScheme == LGM_DERIV_FOUR_POINT){
        *dDaa = (-faa[4] + 8.0*faa[3] - 8.0*faa[1] + faa[0])/(12.0*h);
        *dDap = (-fap[4] + 8.0*fap[3] - 8.0*fap[1] + fap[0])/(12.0*h);
    } else if (DerivScheme == LGM_DERIV_TWO_POINT){
        *dDaa = (faa[2] - faa[0])/(2.0*h);
        *dDap = (fap[2] - fap[0])/(2.0*h);
    }


    return(1);

}




/**
 *  \brief
 *      Integrand for computing \f$ S_b \f$ in a pure centered dipole field.
 *
 *  \details
 *
 *
 *      This routine computes the integrand of the \f$ S_b \f$ integral in a
 *      pure centered dipole field;
 *
 *          \f[ S_b = \displaystyle \int_{\lambda_{sm}}^{\lambda_{nm}}
 *          {\cos\lambda[4-3\cos^2\lambda]^{1/2} d\lambda  \over \left[1-
 *          B/B_m\right]^{1/2}} \f]
 *
 *      Note that \f$ B/B_m \f$ depends on lattitude and on equatorial pitch
 *      angle and is obtained as follows in the code;
 *
 *          \f[ B/B_m = { \sin^2\alpha_\circ [4-3\cos^2\lambda]^{1/2}\over
 *          cos^6\lambda } \f]
 *
 *      \param[in]      Lat         Latitude in radians
 *      \param[in]      qpInfo      Structure contaning additional information needed to compute the integrand.
 *
 *      \return         The value of the integrand
 *
 *      \author         M. Henderson
 *      \date           2010-2011
 */
double  CdipIntegrand_Sb( double Lat, _qpInfo *qpInfo ) {

    double  CosLat, CosLat2, CosLat3, CosLat6, v, arg;
    double  BoverBeq, BoverBm, f;
    Lgm_SummersInfo *si;

    /*
     * Get pointer to parameters for integrand from Lgm_SummersInfo structure.
     *  si->SinAlpha02      pre-computed sin^2( Alpha0 )
     */
    si  = (Lgm_SummersInfo *)qpInfo;

    CosLat  = cos( Lat );     CosLat2 = CosLat*CosLat;
    CosLat3 = CosLat2*CosLat; CosLat6 = CosLat3*CosLat3;
    arg = 4.0 - 3.0*CosLat2;
    v = (arg < 0.0 ) ? 0.0 : sqrt( arg );

    BoverBeq  = v/CosLat6;               // B/Beq
    BoverBm   = BoverBeq*si->SinAlpha02; // B/Bm = B/Beq * sin^2(Alpha0)

    /*
     * We dont need to worry about the possibility that 1-B/Bm is 0 which gives
     * f of infinity. This singularity is not a problem because it is
     * integrable and will never get evaluated when using the proper quadpack
     * integrator.
     */
    f = CosLat*v/sqrt(1.0-BoverBm);

//printf("%g %g\n", Lat*DegPerRad, f);
    return( f );

}

/**
 *  \brief
 *      Integrand for computing \f$ <D_{\alpha\circ\alpha_\circ} > \f$ in a
 *      pure centered dipole field.
 *
 *      This routine computes the integrand;
 *
 *         \f[  g_{\alpha_\circ\alpha_\circ}(\lambda) =
 *         {D_{\alpha_\circ\alpha_\circ}\; \cos\lambda[4-3\cos^2\lambda]^{1/2}
 *         \over \left[1- B/B_m\right]^{1/2}} \f]
 *
 *      In the code, \f$ D_{\alpha_\circ\alpha_\circ} \f$ is determined as
 *      follows;
 *
 *         \f[  D_{\alpha_\circ\alpha_\circ} =
 *         D_{\alpha\alpha}\left({\partial\alpha_\circ\over\partial\alpha}\right)^2
 *         \f]
 *
 *      Because \f$ \sin^2\alpha/B = \sin^2\alpha_{eq}/B_{eq} \f$, we have;
 *
 *         \f[  {\partial\alpha_\circ\over\partial\alpha} =
 *         \tan(\alpha_\circ)/\tan(\alpha) \f]
 *
 */
double  SummersIntegrand_Gaa( double Lat, _qpInfo *qpInfo ) {

    double  CosLat, CosLat2, CosLat3, CosLat6, v, Omega_e, Omega_Sig, x1, x2, xm, dx;
    double  BoverBeq, BoverBm, SinAlpha2, CosAlpha2, TanAlpha2;
    double  Daa, Da0a0, f, aStar, dB, B, dBoverB2;
    Lgm_SummersInfo *si;


    /*
     * Get a pointer toPull parameters for integrand from Lgm_SummersInfo structure.
     *
     *  si->SinAlpha02  pre-computed sin^2( Alpha0 )
     *  si->TanAlpha02  pre-computed tan^2( Alpha0 )
     *  si->E           Dimensionless energy Ek/E0 (kinetic energy over rest energy)
     *  si->L           Dimensionless L-shell parameter (i.e. radius in eq plane).
     *  si->aStarEq     Equatorial value of Summers' "AlphaStar" parameter.
     *  si->Omega_eEq   Equatorial electron gyro-frequency.
     *  si->Omega_SigEq Equatorial gyro-frequency for given species.
     *  si->w1          lower cutoff frequency.
     *  si->w2          upper cutoff frequency.
     *  si->wm          frequency of max power.
     *  si->dw          frequency bandwidth. (Semi bandwidth is sigma*dw).
     *  si->Lambda      -1 for electrons LGM_EPS (=me/mp) for protons.
     *  si->s           +1 for R-mode -1 for L-mode.
     *  si->Rho         Rho = sqrt(PI)/2.0( erf( (wm-w1)/dw ) + erf( (w2-wm)/dw ) )   Eq (3) in Summers2007
     *  si->Sig         Sig -- see Summers [2005], text above eqn (30).
     *  si->n1          nH/Ne
     *  si->n2          nHe/Ne
     *  si->n3          nO/Ne
     *  si->Version     Which version to use. (LGM_SUMMERS_2005 or LGM_SUMMERS_2007).
     */
    si  = (Lgm_SummersInfo *)qpInfo;


    // Return 0.0 if Lat is greater than Latitudinal cutoff for waves.
    if ( fabs(Lat) > si->MaxWaveLat ) return(0.0);


    // Compute dipole-related params
    CosLat  = cos( Lat );     CosLat2 = CosLat*CosLat;
    CosLat3 = CosLat2*CosLat; CosLat6 = CosLat3*CosLat3;
    v = sqrt(4.0 - 3.0*CosLat2);
    BoverBeq  = v/CosLat6;               // B/Beq
    BoverBm   = BoverBeq*si->SinAlpha02; // B/Bm = B/Beq * sin^2(Alpha0)
    SinAlpha2 = BoverBm;                 // sin^2(Alpha) = B/Bm
    if (SinAlpha2 > 1.0) { SinAlpha2 = BoverBm = 1.0; }
    CosAlpha2 = 1.0-SinAlpha2;
    if (CosAlpha2 < 0.0) { CosAlpha2 = 0.0; } // Somewhat moot, because this case will lead to TanAlpha2 = infty
    TanAlpha2 = SinAlpha2/CosAlpha2;

    /*
     * Compute local parameters from the equatorial ones that are specified.
     * Compute local parameters from the equatorial ones that are specified.
     * Assume wave frequency distribution (i.e. defined by the 'w' quantities,
     * e.g. w1, w2, wm, dw) is constant along the field lines.
     */
    B         = BoverBeq*M_CDIP/(si->L*si->L*si->L); // local value of B (warning, hard-coded M value).
    dB        = si->BwFunc( Lat, si->BwFuncData );   // Compute dB from user-supplied function.
    dBoverB2  = dB*dB/(B*B);                         // Local value of R = dB^2/B^2.
    aStar     = si->aStarEq*BoverBeq*BoverBeq;       // Local aStar value.
    Omega_e   = fabs(Lgm_GyroFreq( -LGM_e, B, LGM_ELECTRON_MASS ));  // Local electron gyro-frequency.
    Omega_Sig = si->Omega_SigEq*BoverBeq;    // Local gyro-frequency for given species.
    x1        = si->w1/Omega_e;        // Local value of x1.
    x2        = si->w2/Omega_e;        // Local value of x1.
    xm        = si->wm/Omega_e;        // Local value of xm.
    dx        = si->dw/Omega_e;        // Local value of dx.


    /*
     * Compute the local Daa
     */
    if ( si->Version == LGM_SUMMERS_2005 ) {

        Daa = Lgm_SummersDaaLocal( SinAlpha2, si->E, dBoverB2, BoverBeq, Omega_e, Omega_Sig,
                   si->Rho, si->Sig, x1, x2, xm, dx, si->Lambda, si->s, aStar );

    } else if ( si->Version == LGM_SUMMERS_2007 ) {

        Daa = Lgm_SummersDaaLocal_2007( SinAlpha2, si->E, dBoverB2, BoverBeq, Omega_e, Omega_Sig,
                    si->Rho, si->Sig, x1, x2, xm, dx, si->Lambda, si->s, si->n1, si->n2, si->n3, aStar );

    } else {

        printf("SummersIntegrand_Gaa: Unknown version of Summers Diff. Coeff Model. Version = %d\n", si->Version );

    }


    /*
     * Compute Da0a0. The TanAlpha02/TanAlpha2 term (which is (dAlphaeq/dAlpha)^2 )
     * converts D_Alpha,Alpha to D_Alpha_eq,Alpha_eq.
     */
    Da0a0 = Daa*si->TanAlpha02/TanAlpha2;


    /*
     * Finally the integrand.
     */
    f = Da0a0 * CosLat*v/sqrt(1.0-BoverBm);

//printf("%e %e\n", Lat*DegPerRad, f);
//printf("%g %g\n", Lat*DegPerRad, f);
    return( f );

}

/**
 *  \brief
 *      Integrand for computing \f$ <D_{\alpha\_circ p} > \f$ in a
 *      pure centered dipole field.
 *
 *      This routine computes the integrand;
 *
 *         \f[  g_{\alpha_\circ p }(\lambda) =
 *         {D_{\alpha_\circ p }\; \cos\lambda[4-3\cos^2\lambda]^{1/2}
 *         \over \left[1- B/B_m\right]^{1/2}} \f]
 *
 *      In the code, \f$ D_{\alpha_\circ p } \f$ is determined as
 *      follows;
 *
 *         \f[  D_{\alpha_\circ p } =
 *         D_{\alpha p }{\partial\alpha_\circ\over\partial\alpha}
 *         \f]
 *
 *      Because \f$ \sin^2\alpha/B = \sin^2\alpha_{eq}/B_{eq} \f$, we have;
 *
 *         \f[  {\partial\alpha_\circ\over\partial\alpha} =
 *         \tan(\alpha_\circ)/\tan(\alpha) \f]
 *
 */
double  SummersIntegrand_Gap( double Lat, _qpInfo *qpInfo ) {

    double  CosLat, CosLat2, CosLat3, CosLat6, v, Omega_e, Omega_Sig, x1, x2, xm, dx;
    double  BoverBeq, BoverBm, SinAlpha2, CosAlpha2, TanAlpha;
    double  Dap, Da0p, f, aStar, dB, B, dBoverB2;
    Lgm_SummersInfo *si;

    /*
     * Get a pointer toPull parameters for integrand from Lgm_SummersInfo structure.
     *
     *  si->SinAlpha02  pre-computed sin^2( Alpha0 )
     *  si->TanAlpha0   pre-computed tan( Alpha0 )
     *  si->E           Dimensionless energy Ek/E0 (kinetic energy over rest energy)
     *  si->L           Dimensionless L-shell parameter (i.e. radius in eq plane).
     *  si->aStarEq     Equatorial value of Summers' "AlphaStar" parameter.
     *  si->Omega_eEq   Equatorial electron gyro-frequency.
     *  si->Omega_SigEq Equatorial gyro-frequency for given species.
     *  si->w1          lower cutoff frequency.
     *  si->w2          upper cutoff frequency.
     *  si->wm          frequency of max power.
     *  si->dw          frequency bandwidth. (Semi bandwidth is sigma*dw).
     *  si->Lambda      -1 for electrons LGM_EPS (=me/mp) for protons.
     *  si->s           +1 for R-mode -1 for L-mode.
     *  si->Rho         Rho = sqrt(PI)/2.0( erf( (wm-w1)/dw ) + erf( (w2-wm)/dw ) )   Eq (3) in Summers2007
     *  si->Sig         Sig -- see Summers [2005], text above eqn (30).
     *  si->n1          nH/Ne
     *  si->n2          nHe/Ne
     *  si->n3          nO/Ne
     *  si->Version     Which version to use. (LGM_SUMMERS_2005 or LGM_SUMMERS_2007).
     */
    si  = (Lgm_SummersInfo *)qpInfo;


    // Return 0.0 if Lat is greater than Latitudinal cutoff for waves.
    if ( fabs(Lat) > si->MaxWaveLat ) return(0.0); // Return 0 if beyond cutoff latitude.



    // Compute dipole-related params
    CosLat  = cos( Lat );     CosLat2 = CosLat*CosLat;
    CosLat3 = CosLat2*CosLat; CosLat6 = CosLat3*CosLat3;
    v = sqrt(4.0 - 3.0*CosLat2);
    BoverBeq  = v/CosLat6;               // B/Beq
    BoverBm   = BoverBeq*si->SinAlpha02; // B/Bm = B/Beq * sin^2(Alpha0)
    SinAlpha2 = BoverBm;                 // sin^2(Alpha) = B/Bm
    if (SinAlpha2 > 1.0) { SinAlpha2 = BoverBm = 1.0; }
    CosAlpha2 = 1.0-SinAlpha2;
    if (CosAlpha2 < 0.0) { CosAlpha2 = 0.0; } // Somewhat moot, because this case will lead to TanAlpha2 = infty
    TanAlpha  = sqrt(SinAlpha2/CosAlpha2);

    /*
     * Compute local parameters from the equatorial ones that are specified.
     * Compute local parameters from the equatorial ones that are specified.
     * Assume wave frequency distribution (i.e. defined by the 'w' quantities,
     * e.g. w1, w2, wm, dw) is constant along the field lines.
     */
    B         = BoverBeq*M_CDIP/(si->L*si->L*si->L); // local value of B (warning, hard-coded M value).
    dB        = si->BwFunc( Lat, si->BwFuncData );   // Compute dB from user-supplied function.
    dBoverB2  = dB*dB/(B*B);                         // Local value of R = dB^2/B^2.
    aStar     = si->aStarEq*BoverBeq*BoverBeq;       // Local aStar value.
    Omega_e   = fabs(Lgm_GyroFreq( -LGM_e, B, LGM_ELECTRON_MASS ));  // Local electron gyro-frequency.
    Omega_Sig = si->Omega_SigEq*BoverBeq;    // Local gyro-frequency for given species.
    x1        = si->w1/Omega_e;        // Local value of x1.
    x2        = si->w2/Omega_e;        // Local value of x1.
    xm        = si->wm/Omega_e;        // Local value of xm.
    dx        = si->dw/Omega_e;        // Local value of dx.

    /*
     * Compute the local Dap
     */
    if ( si->Version == LGM_SUMMERS_2005 ) {

        Dap = Lgm_SummersDapLocal( SinAlpha2, si->E, dBoverB2, BoverBeq, Omega_e, Omega_Sig,
                   si->Rho, si->Sig, x1, x2, xm, dx, si->Lambda, si->s, aStar );

    } else if ( si->Version == LGM_SUMMERS_2007 ) {

        Dap = Lgm_SummersDapLocal_2007( SinAlpha2, si->E, dBoverB2, BoverBeq, Omega_e, Omega_Sig,
                    si->Rho, si->Sig, x1, x2, xm, dx, si->Lambda, si->s, si->n1, si->n2, si->n3, aStar );

    } else {

        printf("SummersIntegrand_Gap: Unknown version of Summers Diff. Coeff Model. Version = %d\n", si->Version );
        return(-1);

    }


    /*
     * Compute Da0p. The TanAlpha0/TanAlpha term (which is dAlphaeq/dAlpha)
     * converts D_Alpha,p to D_Alpha_eq,p.
     */
    Da0p = Dap*si->TanAlpha0/TanAlpha;


    /*
     * Finally the integrand.
     */
    f = Da0p * CosLat*v/sqrt(1.0-BoverBm);

    return( f );

}

/**
 *  \brief
 *      Integrand for computing \f$ <D_{pp}> \f$ in a
 *      pure centered dipole field.
 *
 *      This routine computes the integrand;
 *
 *         \f[  g_{pp}(\lambda) =
 *         {D_{pp}\; \cos\lambda[4-3\cos^2\lambda]^{1/2}
 *         \over \left[1- B/B_m\right]^{1/2}} \f]
 *
 *
 */
double  SummersIntegrand_Gpp( double Lat, _qpInfo *qpInfo ) {

    double  CosLat, CosLat2, CosLat3, CosLat6, v, Omega_e, Omega_Sig, x1, x2, xm, dx;
    double  BoverBeq, BoverBm, SinAlpha2;
    double  Dpp, f, aStar, dB, B, dBoverB2;
    Lgm_SummersInfo *si;

    /*
     * Get a pointer toPull parameters for integrand from Lgm_SummersInfo structure.
     *
     *  si->SinAlpha02  pre-computed sin^2( Alpha0 )
     *  si->E           Dimensionless energy Ek/E0 (kinetic energy over rest energy)
     *  si->L           Dimensionless L-shell parameter (i.e. radius in eq plane).
     *  si->aStarEq     Equatorial value of Summers' "AlphaStar" parameter.
     *  si->Omega_eEq   Equatorial electron gyro-frequency.
     *  si->Omega_SigEq Equatorial gyro-frequency for given species.
     *  si->w1          lower cutoff frequency.
     *  si->w2          upper cutoff frequency.
     *  si->wm          frequency of max power.
     *  si->dw          frequency bandwidth. (Semi bandwidth is sigma*dw).
     *  si->Lambda      -1 for electrons LGM_EPS (=me/mp) for protons.
     *  si->s           +1 for R-mode -1 for L-mode.
     *  si->Rho         Rho = sqrt(PI)/2.0( erf( (wm-w1)/dw ) + erf( (w2-wm)/dw ) )   Eq (3) in Summers2007
     *  si->Sig         Sig -- see Summers [2005], text above eqn (30).
     *  si->n1          nH/Ne
     *  si->n2          nHe/Ne
     *  si->n3          nO/Ne
     *  si->Version     Which version to use. (LGM_SUMMERS_2005 or LGM_SUMMERS_2007).
     */
    si  = (Lgm_SummersInfo *)qpInfo;

    // Return 0.0 if Lat is greater than Latitudinal cutoff for waves.
    if ( fabs(Lat) > si->MaxWaveLat ) return(0.0);


    // Compute dipole-related params
    CosLat  = cos( Lat );     CosLat2 = CosLat*CosLat;
    CosLat3 = CosLat2*CosLat; CosLat6 = CosLat3*CosLat3;
    v = sqrt(4.0 - 3.0*CosLat2);
    BoverBeq  = v/CosLat6;               // B/Beq
    BoverBm   = BoverBeq*si->SinAlpha02; // B/Bm = B/Beq * sin^2(Alpha0)
    SinAlpha2 = BoverBm;                 // sin^2(Alpha) = B/Bm
    if (SinAlpha2 > 1.0) { SinAlpha2 = BoverBm = 1.0; }


    /*
     * Compute local parameters from the equatorial ones that are specified.
     * Compute local parameters from the equatorial ones that are specified.
     * Assume wave frequency distribution (i.e. defined by the 'w' quantities,
     * e.g. w1, w2, wm, dw) is constant along the field lines.
     */
    B         = BoverBeq*M_CDIP/(si->L*si->L*si->L); // local value of B (warning, hard-coded M value).
    dB        = si->BwFunc( Lat, si->BwFuncData );   // Compute dB from user-supplied function.
    dBoverB2  = dB*dB/(B*B);                         // Local value of R = dB^2/B^2.
    aStar     = si->aStarEq*BoverBeq*BoverBeq;       // Local aStar value.
    Omega_e   = fabs(Lgm_GyroFreq( -LGM_e, B, LGM_ELECTRON_MASS ));  // Local electron gyro-frequency.
    Omega_Sig = si->Omega_SigEq*BoverBeq;    // Local gyro-frequency for given species.
    x1        = si->w1/Omega_e;        // Local value of x1.
    x2        = si->w2/Omega_e;        // Local value of x1.
    xm        = si->wm/Omega_e;        // Local value of xm.
    dx        = si->dw/Omega_e;        // Local value of dx.

    /*
     * Compute the local Dpp
     */
    if ( si->Version == LGM_SUMMERS_2005 ) {

        Dpp = Lgm_SummersDppLocal( SinAlpha2, si->E, dBoverB2, BoverBeq, Omega_e, Omega_Sig,
                   si->Rho, si->Sig, x1, x2, xm, dx, si->Lambda, si->s, aStar );

    } else if ( si->Version == LGM_SUMMERS_2007 ) {

        Dpp = Lgm_SummersDppLocal_2007( SinAlpha2, si->E, dBoverB2, BoverBeq, Omega_e, Omega_Sig,
                    si->Rho, si->Sig, x1, x2, xm, dx, si->Lambda, si->s, si->n1, si->n2, si->n3, aStar );

    } else {

        printf("SummersIntegrand_Gpp: Unknown version of Summers Diff. Coeff Model. Version = %d\n", si->Version );
        return(-1);

    }


    /*
     * Finally the integrand.
     */
    f = Dpp * CosLat*v/sqrt(1.0-BoverBm);

//printf("%g %g\n", Lat*DegPerRad, f);
    return( f );

}



/**
 *  \brief
 *      Computes the local Summer's [2005] Daa diffusion coefficient.
 *
 *  \details
 *
 *
 *
 *
 *      \param[in]      SinAlpha2   \f$\sin^2(\alpha)\f$, where \f$\alpha\f$ is the local particle pitch angle.
 *      \param[in]      E           Dimensionless energy Ek/E0 (kinetic energy over rest mass).
 *      \param[in]      dBoverB2    Ratio of wave amplitude, dB to local background field, B.
 *      \param[in]      BoverBeq    Ratio of local B to Beq.
 *      \param[in]      Omega_e     Local gyro-frequency of electrons.
 *      \param[in]      Omega_Sig   Local gyro-frequency of particle species we are interested in.
 *      \param[in]      Rho         This is \f$0.5 \sqrt(\pi) ( \mbox{erf}((wm-w1)/dw) + \mbox{erf}(wm-w1)/dw) )\f$.
 *      \param[in]      Sig         Blah.
 *      \param[in]      xm          Blah.
 *      \param[in]      dx          Blah.
 *      \param[in]      Lambda      Particle species. Can be LGM_ELECTRONS or LGM_PROTONS.
 *      \param[in]      s           Mode of the wave. Can be LGM_R_MODE_WAVE ( s = -1 ) or LGM_L_MODE_WAVE ( s = +1 ).
 *      \param[in]      aStar       Local value of the aStar parameter ( \f$ \alpha^* \f$ ) in the Summer's papers.
 *
 *      \return         Local value of Daa
 *
 *      \author         M. Henderson
 *      \date           2010-2011
 *
 */
double Lgm_SummersDaaLocal( double SinAlpha2, double E, double dBoverB2, double BoverBeq, double Omega_e, double Omega_Sig, double Rho, double Sig, double xl, double xh, double xm, double dx, double Lambda, int s, double aStar ) {

    int             nReal, nRoots, n;
    double          Gamma, Gamma2, Beta, Beta2, Mu, Mu2, BetaMu, BetaMu2, OneMinusBetaMu2;
    double          sEpsMinusOne, a1, a2, a3, a4, a, aa, apa, b;
    double complex  z1, z2, z3, z4, z[4];
    double          Daa, R, x0, y0, Ep1, Ep12, xms, xpse, u, arg, c1, c2, c3, c4;
    double          x, x2, x3, x4, y, g, F, DD, fac;

    /*
     * Solve for resonant roots.
     */
    Gamma = E+1.0; Gamma2 = Gamma*Gamma;
    Beta2 = E*(E+2.0) / Gamma2;
    Mu2   = 1.0-SinAlpha2; if (Mu2<0.0) Mu2 = 0.0; // Mu2 is cos^2(Alpha)
    a     = s*Lambda/Gamma; aa = a*a; apa = a+a;
    b     = (1.0 + LGM_EPS)/aStar;

    BetaMu2 = Beta2*Mu2;
    OneMinusBetaMu2 = 1.0 - BetaMu2;

    sEpsMinusOne = s*(LGM_EPS - 1.0);

    a1 = ( apa + sEpsMinusOne*OneMinusBetaMu2 ) / OneMinusBetaMu2;
    a2 = ( aa + apa*sEpsMinusOne - LGM_EPS + BetaMu2*(b+LGM_EPS) ) / OneMinusBetaMu2;
    a3 = ( aa*sEpsMinusOne - apa*LGM_EPS ) / OneMinusBetaMu2;
    a4 = -aa*LGM_EPS / OneMinusBetaMu2;

    nReal = Lgm_QuarticRoots( a1, a2, a3, a4, &z1, &z2, &z3, &z4 );

    R  = dBoverB2;   // The ratio (dB/B)^2


    // This is another hard-wired input parameter. Fix.
    int sum_res = 0;


    /*
     * Gather applicable roots together into the z[] array
     */
    nRoots = 0;
    if ( ( fabs(cimag(z1)) < 1e-10 ) && ( creal(z1) > 0.0 ) ) z[nRoots++] = z1;
    if ( ( fabs(cimag(z2)) < 1e-10 ) && ( creal(z2) > 0.0 ) ) z[nRoots++] = z2;
    if ( ( fabs(cimag(z3)) < 1e-10 ) && ( creal(z3) > 0.0 ) ) z[nRoots++] = z3;
    if ( ( fabs(cimag(z4)) < 1e-10 ) && ( creal(z4) > 0.0 ) ) z[nRoots++] = z4;


    if ( ( nRoots == 0 ) && ( Mu2 > 1e-16 ) ){

        /*
         *  No resonant roots were found
         */
        Daa = 0.0;


    } else if ( Mu2 <= 1e-16 ) {

        /*
         *  If cos(alpha) is zero, then the resonance condition becomnes x = -a = s*Lambda/Gamma.
         *  And x>0 in this case will only happen for (R-MODE with electrons) and (L-MODE with protons)
         */
        if          ( (s ==  1) && ( Lambda < 0.0 ) ) {  // R-MODE with electrons
            x0 = 1.0/Gamma;
        } else if   ( (s == -1) && ( Lambda > 0.0 ) ) { // L-MODE with protons
            x0 = Lambda/Gamma;
        } else {
            x0 = 0.0;
        }

        if ( x0 > 1e-12 ) {
            y0 = x0 * sqrt( 1.0 + b*Gamma2/((Gamma-1.0)*(1.0+LGM_EPS*Gamma)) );
            Ep1 = E+1.0; Ep12 = Ep1*Ep1;
            arg = (x0-xm)/dx;
            DD = M_PI_2/Rho * Omega_Sig*Omega_Sig/fabs(Omega_e) * R/(Ep12*dx) * exp( -arg*arg );
            Daa = ( sum_res ) ? 2.0*DD : DD;
        } else {
            Daa = 0.0;
        }


    } else {

        /*
         *  We have resonant roots and we are not at 90Deg. Use equations
         *  (33)-(35) of Summers2005.
         */
        c1  = 2.0*sEpsMinusOne;
        c2  = 1.0 - 4.0*LGM_EPS + LGM_EPS*LGM_EPS;
        c3  = -sEpsMinusOne * (b + 4.0*LGM_EPS)/2.0;
        c4  = LGM_EPS*(b + LGM_EPS);
        Ep1 = E+1.0; Ep12 = Ep1*Ep1;
        fac = M_PI_2/Rho * Omega_Sig*Omega_Sig/fabs(Omega_e) * R/Ep12;

        /*
         * Execute sum over resonant roots.
         */
        BetaMu = sqrt(BetaMu2);
        Beta   = sqrt(Beta2);
        Mu     = sqrt(Mu2);
        for ( Daa=0.0, n=0; n<nRoots; n++ ){


            x = creal( z[n] ); x2 = x*x; x3 = x2*x; x4 = x2*x2;
            y = ( x+a )/BetaMu;
            if ((x>xl)&&(x<xh)) {
                g = x4 + c1*x3 + c2*x2 + c3*x + c4;
                xms  = x - s;
                xpse = x + s*LGM_EPS;
                F = y*xms*xms*xpse*xpse/(x*g);

                u = 1.0 - x*Mu/(y*Beta);
                arg = (x-xm)/dx;
                Daa += u*u *fabs(F) / ( dx * fabs( BetaMu - F ) ) * exp( -arg*arg );
            }

        }
        Daa *= fac;

    }

    return(Daa);

}

/**
 *  \brief
 *      Computes the local Summer's [2005] Dap diffusion coefficient. (Actually computes Dap/p -- see eqn (34).)
 *
 *  \details
 *
 *
 *
 *
 *      \param[in]      SinAlpha2   \f$\sin^2(\alpha)\f$, where \f$\alpha\f$ is the local particle pitch angle.
 *      \param[in]      E           Dimensionless energy Ek/E0 (kinetic energy over rest mass).
 *      \param[in]      dBoverB2    Ratio of wave amplitude, dB to local background field, B.
 *      \param[in]      BoverBeq    Ratio of local B to Beq.
 *      \param[in]      Omega_e     Local gyro-frequency of electrons.
 *      \param[in]      Omega_Sig   Local gyro-frequency of particle species we are interested in.
 *      \param[in]      Rho         This is \f$0.5 \sqrt(\pi) ( \mbox{erf}((wm-w1)/dw) + \mbox{erf}(wm-w1)/dw) )\f$.
 *      \param[in]      Sig         Blah.
 *      \param[in]      xm          Blah.
 *      \param[in]      dx          Blah.
 *      \param[in]      Lambda      Particle species. Can be LGM_ELECTRONS or LGM_PROTONS.
 *      \param[in]      s           Mode of the wave. Can be LGM_R_MODE_WAVE ( s = -1 ) or LGM_L_MODE_WAVE ( s = +1 ).
 *      \param[in]      aStar       Local value of the aStar parameter ( \f$ \alpha^* \f$ ) in the Summer's papers.
 *
 *      \return         Local value of Daa
 *
 *      \author         M. Henderson
 *      \date           2010-2011
 *
 */
double Lgm_SummersDapLocal( double SinAlpha2, double E, double dBoverB2, double BoverBeq, double Omega_e, double Omega_Sig, double Rho, double Sig, double xl, double xh, double xm, double dx, double Lambda, int s, double aStar ) {

    int             nReal, nRoots, n;
    double          Gamma, Gamma2, Beta, Beta2, Mu, Mu2, BetaMu, BetaMu2, OneMinusBetaMu2;
    double          sEpsMinusOne, a1, a2, a3, a4, a, aa, apa, b, SinAlpha;
    double complex  z1, z2, z3, z4, z[4];
    double          Dap, R, x0, y0, Ep1, Ep12, xms, xpse, u, arg, c1, c2, c3, c4;
    double          x, x2, x3, x4, y, g, F, DD, fac;

    /*
     * Solve for resonant roots.
     */
    Gamma = E+1.0; Gamma2 = Gamma*Gamma;
    Beta2 = E*(E+2.0) / Gamma2;
    Mu2   = 1.0-SinAlpha2; if (Mu2<0.0) Mu2 = 0.0; // Mu2 is cos^2(Alpha)
    a     = s*Lambda/Gamma; aa = a*a; apa = a+a;
    b     = (1.0 + LGM_EPS)/aStar;

    BetaMu2 = Beta2*Mu2;
    OneMinusBetaMu2 = 1.0 - BetaMu2;

    sEpsMinusOne = s*(LGM_EPS - 1.0);

    a1 = ( apa + sEpsMinusOne*OneMinusBetaMu2 ) / OneMinusBetaMu2;
    a2 = ( aa + apa*sEpsMinusOne - LGM_EPS + BetaMu2*(b+LGM_EPS) ) / OneMinusBetaMu2;
    a3 = ( aa*sEpsMinusOne - apa*LGM_EPS ) / OneMinusBetaMu2;
    a4 = -aa*LGM_EPS / OneMinusBetaMu2;

    nReal = Lgm_QuarticRoots( a1, a2, a3, a4, &z1, &z2, &z3, &z4 );

    R  = dBoverB2;   // The ratio (dB/B)^2

    // This is another hard-wired input parameter. Fix.
    int sum_res = 0;


    /*
     * Gather applicable roots together into the z[] array
     */
    nRoots = 0;
    if ( ( fabs(cimag(z1)) < 1e-10 ) && ( fabs(creal(z1)) > 0.0 ) ) z[nRoots++] = z1;
    if ( ( fabs(cimag(z2)) < 1e-10 ) && ( fabs(creal(z2)) > 0.0 ) ) z[nRoots++] = z2;
    if ( ( fabs(cimag(z3)) < 1e-10 ) && ( fabs(creal(z3)) > 0.0 ) ) z[nRoots++] = z3;
    if ( ( fabs(cimag(z4)) < 1e-10 ) && ( fabs(creal(z4)) > 0.0 ) ) z[nRoots++] = z4;


    if ( ( nRoots == 0 ) && ( Mu2 > 1e-16 ) ){

        /*
         *  No resonant roots were found
         */
        Dap = 0.0;


    } else if ( Mu2 <= 1e-16 ) {

        /*
         *  If cos(alpha) is zero, then the resonance condition becomnes x = -a = s*Lambda/Gamma.
         *  And x>0 in this case will only happen for (R-MODE with electrons) and (L-MODE with protons)
         */
        if          ( (s ==  1) && ( Lambda < 0.0 ) ) {  // R-MODE with electrons
            x0 = 1.0/Gamma;
        } else if   ( (s == -1) && ( Lambda > 0.0 ) ) { // L-MODE with protons
            x0 = Lambda/Gamma;
        } else {
            x0 = 0.0;
        }

        if ( x0 > 1e-12 ) {
            y0 = x0 * sqrt( 1.0 + b*Gamma2/((Gamma-1.0)*(1.0+LGM_EPS*Gamma)) );
            Beta  = sqrt(Beta2);
            Ep1   = E+1.0; Ep12 = Ep1*Ep1; arg   = (x0-xm)/dx;
            DD = -M_PI_2/Rho * Omega_Sig*Omega_Sig/fabs(Omega_e) * R*x0/(Ep12*dx*y0*Beta)  * exp( -arg*arg );
            Dap   = ( sum_res ) ? 2.0*DD : DD;
        } else {
            Dap = 0.0;
        }


    } else {

        /*
         *  We have resonant roots and we are not at 90Deg. Use equations
         *  (33)-(35) of Summers2005.
         */
        Beta   = sqrt(Beta2);
        if (SinAlpha2 < 0.0) SinAlpha = 0.0;
        else if ( SinAlpha2 > 1.0) SinAlpha = 1.0;
        else SinAlpha = sqrt( SinAlpha2 );
        c1  = 2.0*sEpsMinusOne;
        c2  = 1.0 - 4.0*LGM_EPS + LGM_EPS*LGM_EPS;
        c3  = -sEpsMinusOne * (b + 4.0*LGM_EPS)/2.0;
        c4  = LGM_EPS*(b + LGM_EPS);
        Ep1 = E+1.0; Ep12 = Ep1*Ep1;
        fac = -M_PI_2/Rho * Omega_Sig*Omega_Sig/fabs(Omega_e) * R*SinAlpha/Ep12 /Beta;

        /*
         * Execute sum over resonant roots.
         */
        BetaMu = sqrt(BetaMu2);
        Mu     = sqrt(Mu2);
        for ( Dap=0.0, n=0; n<nRoots; n++ ){


            x = creal( z[n] ); x2 = x*x; x3 = x2*x; x4 = x2*x2;
            y = ( x+a )/BetaMu;
            if ((x>xl)&&(x<xh)) {
                g = x4 + c1*x3 + c2*x2 + c3*x + c4;
                xms  = x - s;
                xpse = x + s*LGM_EPS;
                F = y*xms*xms*xpse*xpse/(x*g);

                u = 1.0 - x*Mu/(y*Beta);
                arg = (x-xm)/dx;
                Dap += x/y * u *fabs(F) / ( dx * fabs( BetaMu - F ) ) * exp( -arg*arg );
            }

        }
        Dap *= fac;

    }

    return(Dap);

}



/**
 *  \brief
 *      Computes the local Summer's [2005] Dpp diffusion coefficient. (Actually computes Dpp/p^2 -- see eqn (34).)
 *
 *  \details
 *
 *
 *
 *
 *      \param[in]      SinAlpha2   \f$\sin^2(\alpha)\f$, where \f$\alpha\f$ is the local particle pitch angle.
 *      \param[in]      E           Dimensionless energy Ek/E0 (kinetic energy over rest mass).
 *      \param[in]      dBoverB2    Ratio of wave amplitude, dB to local background field, B.
 *      \param[in]      BoverBeq    Ratio of local B to Beq.
 *      \param[in]      Omega_e     Local gyro-frequency of electrons.
 *      \param[in]      Omega_Sig   Local gyro-frequency of particle species we are interested in.
 *      \param[in]      Rho         This is \f$0.5 \sqrt(\pi) ( \mbox{erf}((wm-w1)/dw) + \mbox{erf}(wm-w1)/dw) )\f$.
 *      \param[in]      Sig         Blah.
 *      \param[in]      xm          Blah.
 *      \param[in]      dx          Blah.
 *      \param[in]      Lambda      Particle species. Can be LGM_ELECTRONS or LGM_PROTONS.
 *      \param[in]      s           Mode of the wave. Can be LGM_R_MODE_WAVE ( s = -1 ) or LGM_L_MODE_WAVE ( s = +1 ).
 *      \param[in]      aStar       Local value of the aStar parameter ( \f$ \alpha^* \f$ ) in the Summer's papers.
 *
 *      \return         Local value of Daa
 *
 *      \author         M. Henderson
 *      \date           2010-2011
 *
 */
double Lgm_SummersDppLocal( double SinAlpha2, double E, double dBoverB2, double BoverBeq, double Omega_e, double Omega_Sig, double Rho, double Sig, double xl, double xh, double xm, double dx, double Lambda, int s, double aStar ) {

    int             nReal, nRoots, n;
    double          Gamma, Gamma2, Beta, Beta2, Mu, Mu2, BetaMu, BetaMu2, OneMinusBetaMu2;
    double          sEpsMinusOne, a1, a2, a3, a4, a, aa, apa, b;
    double complex  z1, z2, z3, z4, z[4];
    double          Dpp, R, x0, y0, Ep1, Ep12, xms, xpse, arg, c1, c2, c3, c4;
    double          x, x2, x3, x4, y, g, F, DD, fac;

    /*
     * Solve for resonant roots.
     */
    Gamma = E+1.0; Gamma2 = Gamma*Gamma;
    Beta2 = E*(E+2.0) / Gamma2;
    Mu2   = 1.0-SinAlpha2; if (Mu2<0.0) Mu2 = 0.0; // Mu2 is cos^2(Alpha)
    a     = s*Lambda/Gamma; aa = a*a; apa = a+a;
    b     = (1.0 + LGM_EPS)/aStar;

    BetaMu2 = Beta2*Mu2;
    OneMinusBetaMu2 = 1.0 - BetaMu2;

    sEpsMinusOne = s*(LGM_EPS - 1.0);

    a1 = ( apa + sEpsMinusOne*OneMinusBetaMu2 ) / OneMinusBetaMu2;
    a2 = ( aa + apa*sEpsMinusOne - LGM_EPS + BetaMu2*(b+LGM_EPS) ) / OneMinusBetaMu2;
    a3 = ( aa*sEpsMinusOne - apa*LGM_EPS ) / OneMinusBetaMu2;
    a4 = -aa*LGM_EPS / OneMinusBetaMu2;
    //printf("a1, a2, a3, a4 = %g %g %g %g\n", a1, a2, a3, a4 );

    nReal = Lgm_QuarticRoots( a1, a2, a3, a4, &z1, &z2, &z3, &z4 );

    R  = dBoverB2;   // The ratio (dB/B)^2

    // This is another hard-wired input parameter. Fix.
    int sum_res = 1;

    /*
     * Gather applicable roots together into the z[] array
     */
    nRoots = 0;
    if ( ( fabs(cimag(z1)) < 1e-10 ) && ( fabs(creal(z1)) > 0.0 ) ) z[nRoots++] = z1;
    if ( ( fabs(cimag(z2)) < 1e-10 ) && ( fabs(creal(z2)) > 0.0 ) ) z[nRoots++] = z2;
    if ( ( fabs(cimag(z3)) < 1e-10 ) && ( fabs(creal(z3)) > 0.0 ) ) z[nRoots++] = z3;
    if ( ( fabs(cimag(z4)) < 1e-10 ) && ( fabs(creal(z4)) > 0.0 ) ) z[nRoots++] = z4;


    if ( ( nRoots == 0 ) && ( Mu2 > 1e-16 ) ){

        /*
         *  No resonant roots were found
         */
        Dpp = 0.0;


    } else if ( Mu2 <= 1e-16 ) {

        /*
         *  If cos(alpha) is zero, then the resonance condition becomnes x = -a = s*Lambda/Gamma.
         *  And x>0 in this case will only happen for (R-MODE with electrons) and (L-MODE with protons)
         */
        if          ( (s ==  1) && ( Lambda < 0.0 ) ) {  // R-MODE with electrons
            x0 = 1.0/Gamma;
        } else if   ( (s == -1) && ( Lambda > 0.0 ) ) { // L-MODE with protons
            x0 = Lambda/Gamma;
        } else {
            x0 = 0.0;
        }

        if ( x0 > 1e-12 ) {
            y0 = x0 * sqrt( 1.0 + b*Gamma2/((Gamma-1.0)*(1.0+LGM_EPS*Gamma)) );
            Ep1   = E+1.0; Ep12 = Ep1*Ep1; arg   = (x0-xm)/dx;
            DD = M_PI_2/Rho * Omega_Sig*Omega_Sig/fabs(Omega_e) * R*x0*x0/(Ep12*dx*y0*y0*Beta2)  * exp( -arg*arg );
            Dpp   = ( sum_res ) ? 2.0*DD : DD;
        } else {
            Dpp = 0.0;
        }


    } else {

        /*
         *  We have resonant roots and we are not at 90Deg. Use equations
         *  (33)-(35) of Summers2005.
         */
        c1  = 2.0*sEpsMinusOne;
        c2  = 1.0 - 4.0*LGM_EPS + LGM_EPS*LGM_EPS;
        c3  = -sEpsMinusOne * (b + 4.0*LGM_EPS)/2.0;
        c4  = LGM_EPS*(b + LGM_EPS);
        Ep1 = E+1.0; Ep12 = Ep1*Ep1;
        fac = M_PI_2/Rho * Omega_Sig*Omega_Sig/fabs(Omega_e) * R*SinAlpha2/Ep12 /Beta2;

        /*
         * Execute sum over resonant roots.
         */
        Beta   = sqrt(Beta2);
        BetaMu = sqrt(BetaMu2);
        Mu     = sqrt(Mu2);
        for ( Dpp=0.0, n=0; n<nRoots; n++ ){


            x = creal( z[n] ); x2 = x*x; x3 = x2*x; x4 = x2*x2;
            y = ( x+a )/BetaMu;
            if ((x>xl)&&(x<xh)) {
                g = x4 + c1*x3 + c2*x2 + c3*x + c4;
                xms  = x - s;
                xpse = x + s*LGM_EPS;
                F = y*xms*xms*xpse*xpse/(x*g);

                arg = (x-xm)/dx;
                Dpp += x*x/(y*y) *fabs(F) / ( dx * fabs( BetaMu - F ) ) * exp( -arg*arg );
            }

        }
        Dpp *= fac;

    }

    return(Dpp);

}


/**
 *  \brief
 *      Computes the local Summer's [2007] Daa diffusion coefficient.
 *
 *  \details
 *      The derivative of the dispersion relation is not given by Summers et al. [2007], so here it is;
 *
 *          \f[
 *              2 y(x) y'(x)=\left(\frac{\frac{e  n_1 }{(e
 *                  s+x)^2}+\frac{1}{(s-x)^2}+\frac{4 e  n_2 }{(e s+4 x)^2}+\frac{16 e
 *                   n_3 }{(e s+16 x)^2}}{ \alpha^*  x}-\frac{-\frac{e  n_1 }{e
 *                  s+x}+\frac{1}{s-x}-\frac{e  n_2 }{e s+4 x}-\frac{e  n_3 }{e s+16
 *                  x}}{ \alpha^*  x^2}\right) x^2+2 \left(\frac{-\frac{e  n_1 }{e
 *                  s+x}+\frac{1}{s-x}-\frac{e  n_2 }{e s+4 x}-\frac{e  n_3 }{e s+16
 *                  x}}{ \alpha^*  x}+1\right) x
 *          \f]
 *
 *      Let, \f$ 2 y(x) y'(x) = Q\f$ so that,
 *
 *          \f[
 *              Q = \left(\frac{\frac{e  n_1 }{(e s+x)^2}+\frac{1}{(s-x)^2}+\frac{4 e
 *                   n_2 }{(e s+4 x)^2}+\frac{16 e  n_3 }{(e s+16 x)^2}}{ \alpha^*
 *                  x}-\frac{-\frac{e  n_1 }{e s+x}+\frac{1}{s-x}-\frac{e  n_2 }{e s+4
 *                  x}-\frac{e  n_3 }{e s+16 x}}{ \alpha^*  x^2}\right) x^2+2
 *                  \left(\frac{-\frac{e  n_1 }{e s+x}+\frac{1}{s-x}-\frac{e  n_2 }{e
 *                  s+4 x}-\frac{e  n_3 }{e s+16 x}}{ \alpha^*  x}+1\right) x
 *          \f]
 *
 *      Then, \f$ W=1/Q \f$ gives,
 *
 *          \f[
 *              W = \frac{ \alpha^*  (s-x)^2 (e s+x)^2 (e s+4 x)^2 (e s+16 x)^2}{s^5 \left(2
 *                   \alpha^*  x s^3-\left(4  \alpha^*
 *                  x^2+ n_1 + n_2 + n_3 -1\right) s^2+2 x \left( \alpha^*
 *                  x^2+ n_1 + n_2 + n_3 \right) s-( n_1 + n_2 + n_3 )
 *                  x^2\right) e^6+2 s^4 x \left(42  \alpha^*  x s^3-\left(84  \alpha^*
 *                  x^2+20  n_1 +17  n_2 +5  n_3 -21\right) s^2+2 x \left(21
 *                   \alpha^*  x^2+20  n_1 +17  n_2 +5  n_3 \right) s-(20
 *                   n_1 +17  n_2 +5  n_3 ) x^2\right) e^5+3 s^3 x^2 \left(406
 *                   \alpha^*  x s^3-\left(812  \alpha^*  x^2+176  n_1 +107
 *                   n_2 +11  n_3 -203\right) s^2+2 x \left(203  \alpha^*  x^2+176
 *                   n_1 +107  n_2 +11  n_3 \right) s-(176  n_1 +107
 *                   n_2 +11  n_3 ) x^2\right) e^4+8 s^2 x^3 \left(914  \alpha^*  x
 *                  s^3-\left(1828  \alpha^*  x^2+320  n_1 +68  n_2 +5
 *                   n_3 -457\right) s^2+2 x \left(457  \alpha^*  x^2+320  n_1 +68
 *                   n_2 +5  n_3 \right) s-(320  n_1 +68  n_2 +5  n_3 )
 *                  x^2\right) e^3+16 s x^4 \left(1218  \alpha^*  x s^3-\left(2436
 *                   \alpha^*  x^2+256  n_1 +16  n_2 + n_3 -609\right) s^2+2 x
 *                  \left(609  \alpha^*  x^2+256  n_1 +16  n_2 + n_3 \right)
 *                  s-(256  n_1 +16  n_2 + n_3 ) x^2\right) e^2+10752 s x^5
 *                  \left(2  \alpha^*  x^3-4  \alpha^*  s x^2+2  \alpha^*  s^2
 *                  x+s\right) e+4096 x^6 \left(2  \alpha^*  x^3-4  \alpha^*  s x^2+2
 *                   \alpha^*  s^2 x+s\right)}
 *          \f]
 *
 *      Define, \f$ P =  \alpha^*  (s-x)^2 (e s+x)^2 (e s+4 x)^2 (e s+16 x)^2 \f$, so that \f$ W = P/H \f$, where,
 *
 *          \f[
 *              H = s^5 \left(2  \alpha^*  x s^3-\left(4  \alpha^*
 *                  x^2+ n_1 + n_2 + n_3 -1\right) s^2+2 x \left( \alpha^*
 *                  x^2+ n_1 + n_2 + n_3 \right) s-( n_1 + n_2 + n_3 )
 *                  x^2\right) e^6+2 s^4 x \left(42  \alpha^*  x s^3-\left(84  \alpha^*
 *                  x^2+20  n_1 +17  n_2 +5  n_3 -21\right) s^2+2 x \left(21
 *                   \alpha^*  x^2+20  n_1 +17  n_2 +5  n_3 \right) s-(20
 *                   n_1 +17  n_2 +5  n_3 ) x^2\right) e^5+3 s^3 x^2 \left(406
 *                   \alpha^*  x s^3-\left(812  \alpha^*  x^2+176  n_1 +107
 *                   n_2 +11  n_3 -203\right) s^2+2 x \left(203  \alpha^*  x^2+176
 *                   n_1 +107  n_2 +11  n_3 \right) s-(176  n_1 +107
 *                   n_2 +11  n_3 ) x^2\right) e^4+8 s^2 x^3 \left(914  \alpha^*  x
 *                  s^3-\left(1828  \alpha^*  x^2+320  n_1 +68  n_2 +5
 *                   n_3 -457\right) s^2+2 x \left(457  \alpha^*  x^2+320  n_1 +68
 *                   n_2 +5  n_3 \right) s-(320  n_1 +68  n_2 +5  n_3 )
 *                  x^2\right) e^3+16 s x^4 \left(1218  \alpha^*  x s^3-\left(2436
 *                   \alpha^*  x^2+256  n_1 +16  n_2 + n_3 -609\right) s^2+2 x
 *                  \left(609  \alpha^*  x^2+256  n_1 +16  n_2 + n_3 \right)
 *                  s-(256  n_1 +16  n_2 + n_3 ) x^2\right) e^2+10752 s x^5
 *                  \left(2  \alpha^*  x^3-4  \alpha^*  s x^2+2  \alpha^*  s^2
 *                  x+s\right) e+4096 x^6 \left(2  \alpha^*  x^3-4  \alpha^*  s x^2+2
 *                   \alpha^*  s^2 x+s\right)
 *          \f]
 *      Collect terms in powers of \f$x\f$, and note that \f$s=\pm 1\f$. This gives,
 *
 *          \f{eqnarray*}{
 *              H &=& 8192  \alpha^*  x^9\\
 *                &+& (21504  \alpha^*  e s-16384  \alpha^*  s) x^8\\
 *                &+& \left(19488  \alpha^*  e^2-43008  \alpha^*  e+8192 \alpha^* \right) x^7\\
 *                &+& \left(7312  \alpha^*  s e^3-38976  \alpha^*  s e^2-16 (256  n_1 +16  n_2 + n_3 ) s e^2+21504  \alpha^*  s e+4096 s\right) x^6\\
 *                &+& \left(1218  \alpha^*  e^4-14624  \alpha^*  e^3-8 (320  n_1 +68  n_2 +5  n_3 ) e^3+19488  \alpha^*  e^2+8192 n_1  e^2+512  n_2  e^2+32  n_3  e^2+10752 e\right) x^5\\
 *                &+& \left(84  \alpha^*  s e^5-2436  \alpha^*  s e^4-3 (176 n_1 +107  n_2 +11  n_3 ) s e^4+7312  \alpha^*  s e^3+5120 n_1  s e^3+1088  n_2  s e^3+80  n_3  s e^3-4096  n_1  s e^2-256  n_2  s e^2-16  n_3  s e^2+9744 s e^2\right) x^4\\
 *                &+& \left(2 \alpha^*  e^6-168  \alpha^*  e^5-2 (20  n_1 +17  n_2 +5 n_3 ) e^5+1218  \alpha^*  e^4+1056  n_1  e^4+642  n_2 e^4+66  n_3  e^4-2560  n_1  e^3-544  n_2  e^3-40  n_3 e^3+3656 e^3\right) x^3\\
 *                &+& \left(-4  \alpha^*  s e^6-( n_1 + n_2 + n_3 ) s e^6+84  \alpha^*  s e^5+80 n_1  s e^5+68  n_2  s e^5+20  n_3  s e^5-528  n_1  s e^4-321  n_2  s e^4-33  n_3  s e^4+609 s e^4\right) x^2\\
 *                &+& \left(2 \alpha^*  e^6+2  n_1  e^6+2  n_2  e^6+2  n_3  e^6-40 n_1  e^5-34  n_2  e^5-10  n_3  e^5+42 e^5\right) x\\
 *                &+& e^6 s-e^6 n_1  s-e^6  n_2  s-e^6  n_3  s
 *          \f}
 *
 *      Let \f$ H =  \Sigma_{i=0}^{8} h_i x^i  \f$. This gives,
 *
 *          \f{eqnarray*}{
 *              h_0 &=& -e^6 ( n_1 + n_2 + n_3 -1) s\\
 *              h_1 &=& 2 e^5 ( \alpha^*  e+ n_2  e+ n_3  e+(e-20)  n_1 -17 n_2 -5  n_3 +21)\\
 *              h_2 &=& e^4 \left(-4  \alpha^*  e^2-( n_1 + n_2 + n_3 ) e^2+84 \alpha^*  e+80  n_1  e+68  n_2  e+20  n_3  e-528 n_1 -321  n_2 -33  n_3 +609\right) s\\
 *              h_3 &=& 2 e^3 \left(-17  n_2  e^2-5  n_3  e^2+ \alpha^*  \left(e^2-84 e+609\right) e+321  n_2  e+33  n_3  e-4 \left(5 e^2-132 e+320\right)  n_1 -272  n_2 -20  n_3 +1828\right)\\
 *              h_4 &=& e^2 \left(-321  n_2  e^2-33  n_3  e^2+4  \alpha^*  \left(21 e^2-609 e+1828\right) e+1088  n_2  e+80  n_3  e-16 \left(33 e^2-320 e+256\right)  n_1 -256  n_2 -16  n_3 +9744\right) s\\
 *              h_5 &=& 2 e \left( \alpha^*  e \left(609 e^2-7312 e+9744\right)-4 \left((320 n_1 +68  n_2 +5  n_3 ) e^2-4 (256  n_1 +16 n_2 + n_3 ) e-1344\right)\right)\\
 *              h_6 &=& 16 \left(-(256  n_1 +16  n_2 + n_3 ) e^2+ \alpha^*  \left(457 e^2-2436 e+1344\right) e+256\right) s\\
 *              h_7 &=& 32  \alpha^*  \left(609 e^2-1344 e+256\right)\\
 *              h_8 &=& 1024  \alpha^*  (21 e-16) s\\
 *              h_9 &=& 8192  \alpha^*
 *          \f}
 *
 *
 *  Then, collecting terms in powers of \f$ e\f$, and defining,
 *
 *          \f{eqnarray*}{
 *              q0 &=& 16  n_1 +4  n_2 + n_3 \\
 *              q1 &=&  n_1 + n_2 + n_3 \\
 *              q2 &=& 20  n_1 +17  n_2 +5  n_3 \\
 *              q3 &=& 176  n_1 +107  n_2 +11  n_3 \\
 *              q4 &=& 256  n_1 +16  n_2 + n_3 \\
 *              q5 &=& 320  n_1 +68  n_2 +5  n_3
 *          \f}
 *
 *  We have,
 *
 *
 *          \f{eqnarray*}{
 *              h_0 &=& e^6 (1- q_1 ) s\\
 *              h_1 &=& 2 e^5 (e ( \alpha^* + q_1 )- q_2 +21)\\
 *              h_2 &=& e^4 \left((-4  \alpha^* - q_1 ) e^2+4 (21  \alpha^* + q_2 ) e+3 (203- q_3 )\right) s\\
 *              h_3 &=& 2 e^3 \left( \alpha^*  e^3+(-84  \alpha^* - q_2 ) e^2+3 (203  \alpha^* + q_3 ) e+4 (457- q_5 )\right)\\
 *              h_4 &=& e^2 \left(84  \alpha^*  e^3+3 (-812  \alpha^* - q_3 ) e^2+16 (457  \alpha^* + q_5 ) e+16 (609- q_4 )\right) s\\
 *              h_5 &=& e \left(1218  \alpha^*  e^3-8 (1828  \alpha^* + q_5 ) e^2+32 (609  \alpha^* + q_4 ) e+10752\right)\\
 *              h_6 &=& 16 \left(457  \alpha^*  e^3+(-2436  \alpha^* - q_4 ) e^2+1344  \alpha^*  e+256\right) s\\
 *              h_7 &=& 32  \alpha^*  \left(609 e^2-1344 e+256\right)\\
 *              h_8 &=& 1024  \alpha^*  (21 e-16) s\\
 *              h_9 &=& 8192  \alpha^*
 *          \f}
 *
 *  But, note that by definition, \f$ q_1 = 0\f$. Therefore \f$ h_0 = 0\f$. Define \f$ H = xG = x\Sigma_{i=0}^{8} g_i x^i \f$, so that;
 *
 *          \f{eqnarray*}{
 *              g_0 &=& 2 e^5 ( \alpha^*  e- q_2 +21)\\
 *              g_1 &=& e^4 \left(-4  \alpha^*  e^2+4 (21  \alpha^* + q_2 ) e+3 (203- q_3 )\right) s\\
 *              g_2 &=& 2 e^3 \left( \alpha^*  e^3+(-84  \alpha^* - q_2 ) e^2+3 (203  \alpha^* + q_3 ) e+4 (457- q_5 )\right)\\
 *              g_3 &=& e^2 \left(84  \alpha^*  e^3+3 (-812  \alpha^* - q_3 ) e^2+16 (457  \alpha^* + q_5 ) e+16 (609- q_4 )\right) s\\
 *              g_4 &=& e \left(1218  \alpha^*  e^3-8 (1828  \alpha^* + q_5 ) e^2+32 (609  \alpha^* + q_4 ) e+10752\right)\\
 *              g_5 &=& 16 \left(457  \alpha^*  e^3+(-2436  \alpha^* - q_4 ) e^2+1344  \alpha^*  e+256\right) s\\
 *              g_6 &=& 32  \alpha^*  \left(609 e^2-1344 e+256\right)\\
 *              g_7 &=& 1024  \alpha^*  (21 e-16) s\\
 *              g_8 &=& 8192  \alpha^*
 *          \f}
 *
 *
 *  Putting all of this together we have,
 *
 *      \f{eqnarray*}{
 *          {dx\over dy} &=& {2 y P\over xG } \\
 *          {dx\over dy} &=& {2 y \alpha^*  (s-x)^2 (e s+x)^2 (e s+4 x)^2 (e s+16 x)^2   \over x \Sigma_{i=0}^{8} g_i x^i }
 *      \f}
 *
 *
 *  Also, with these defintions, the dispersion relation can be simplified to (see equation (18) in Summers et al. [2007]),
 *
 *      \f[
 *          \Sigma_{i=0}^{6} A_i x^i = 0
 *      \f]
 *
 *
 *  where,
 *
 *      \f{eqnarray*}{
 *          A_1 &=& \frac{128 a+4  (1-\xi^2)   p_1  s}{ (64(1-\xi^2))} \\
 *          A_2 &=& \frac{64  a^2 +21 e  (1-\xi^2)   p_2 +8 a  p_1  s+\frac{A  \xi^2 }{ \alpha^* }}{ (64(1-\xi^2))}\\
 *          A_3 &=& \frac{42 a e  p_2 +4  a^2   p_1  s+ e^2   (1-\xi^2)   p_3  s+\frac{B  \xi^2 }{ \alpha^* }}{ (64(1-\xi^2))}\\
 *          A_4 &=& \frac{- e^3   (1-\xi^2) +21  a^2  e  p_2 +2 a  e^2   p_3  s+\frac{C  \xi^2 }{ \alpha^* }}{ (64(1-\xi^2))}\\
 *          A_5 &=& \frac{a  e^2  (a  p_3  s-2 e)}{ (64(1-\xi^2))}\\
 *          A_6 &=& -\frac{ a^2   e^3 }{ (64(1-\xi^2))}
 *      \f}
 *
 *
 *
 *
 *
 *  and,
 *
 *      \f{eqnarray*}{
 *          A &=& 4 e q_0 + 64\\
 *          B &=& e (4 (21-q_0) + e  q_2 ) s\\
 *          C &=& e^2 (21- q_2  ) \\
 *          p_1 &=& 21 e - 16 \\
 *          p_2 &=&  e - 4 \\
 *          p_3 &=&  e - 21
 *      \f}
 *
 *
 *
 *      \param[in]      SinAlpha2   \f$\sin^2(\alpha)\f$, where \f$\alpha\f$ is the local particle pitch angle.
 *      \param[in]      E           Dimensionless energy Ek/E0 (kinetic energy over rest mass).
 *      \param[in]      dBoverB2    Ratio of wave amplitude, dB to local background field, B.
 *      \param[in]      BoverBeq    Ratio of local B to Beq.
 *      \param[in]      Omega_e     Local gyro-frequency of electrons.
 *      \param[in]      Omega_Sig   Local gyro-frequency of particle species we are interested in.
 *      \param[in]      Rho         This is \f$0.5 \sqrt(\pi) ( \mbox{erf}((wm-w1)/dw) + \mbox{erf}(wm-w1)/dw) )\f$.
 *      \param[in]      Sig         Blah.
 *      \param[in]      xm          Blah.
 *      \param[in]      dx          Blah.
 *      \param[in]      Lambda      Particle species. Can be LGM_ELECTRONS or LGM_PROTONS.
 *      \param[in]      s           Mode of the wave. Can be LGM_R_MODE_WAVE ( s = -1 ) or LGM_L_MODE_WAVE ( s = +1 ).
 *      \param[in]      n1          Ratio of H+ Number Density to Electron Number Density. (n1+n2+n3 must add to 1).
 *      \param[in]      n2          Ratio of He+ Number Density to Electron Number Density. (n1+n2+n3 must add to 1).
 *      \param[in]      n3          Ratio of O+ Number Density to Electron Number Density. (n1+n2+n3 must add to 1).
 *      \param[in]      aStar       Local value of the aStar parameter ( \f$ \alpha^* \f$ ) in the Summer's papers.
 *
 *      \return         Local value of Daa
 *
 *      \author         M. Henderson
 *      \date           2010-2011
 *
 */
double Lgm_SummersDaaLocal_2007( double SinAlpha2, double E, double dBoverB2, double BoverBeq, double Omega_e, double Omega_Sig, double Rho, double Sig, double xl, double xh, double xm, double dx, double Lambda, int s, double n1, double n2, double n3, double aStar ) {

    int             nRoots, n, i, gsl_err;
    double          q0, q2, q3, q4, q5;
    double          p1, p2, p3, p4;
    double          Gamma, Gamma2, Beta, Beta2, Mu, Mu2, BetaMu;
    double          A, B, C, a, aa, apa, Coeff[7];
    double complex  z[7], zz[7];
    double          gsl_z[20];
    double          Daa, R, x0, y0, Ep1, Ep12, sss, u, arg;
    double          x, x2, x3, x4, x5, x6, x7, x8;
    double          g, g0, g1, g2, g3, g4, g5, g6, g7, g8;
    double          xi2, xi2OveraStar, OneMinusxi2, OneMinusxi2Times64;
    double          y, F, DD, fac;

    /*
     * Solve for resonant roots.
     */
    Gamma = E+1.0; Gamma2 = Gamma*Gamma;
    Beta2 = E*(E+2.0) / Gamma2;
    Mu2   = 1.0-SinAlpha2; if (Mu2<0.0) Mu2 = 0.0; // Mu2 is cos^2(Alpha)
    xi2  = Beta2*Mu2;
    OneMinusxi2 = (1.0 - xi2);
    OneMinusxi2Times64 = 64.0*OneMinusxi2;
    a  = s*Lambda/Gamma; aa = a*a; apa = a+a;

    q0 = 16.0*n1 + 4.0*n2 + n3;
    q2 = 20.0*n1 + 17.0*n2 + 5.0*n3;

    p1 = 21.0*LGM_EPS - 16.0;
    p2 = LGM_EPS - 4.0;
    p3 = LGM_EPS - 21.0;
    p4 = LGM_EPS + 21.0;

    A = 64.0 + 4.0*LGM_EPS*q0;
    B = LGM_EPS*s*( 4.0*(21.0-q0) + LGM_EPS*q2 );
    C = LGM_EPS2*( p4 - q2 );

    xi2OveraStar = xi2/aStar;

    Coeff[0] = -aa*LGM_EPS3 / OneMinusxi2Times64;                                                                        // A6  x^0 Coeff
    Coeff[1] = a*LGM_EPS2*( a*s*p3 - 2.0*LGM_EPS ) / OneMinusxi2Times64;                                                 // A5  x^1 Coeff
    Coeff[2] = ( 21.0*aa*LGM_EPS*p2 + 2.0*a*LGM_EPS2*s*p3 - LGM_EPS3*OneMinusxi2 + C*xi2OveraStar) / OneMinusxi2Times64; // A4  x^2 Coeff
    Coeff[3] = ( 42.0*a*LGM_EPS*p2 + 4.0*aa*s*p1 + LGM_EPS2*s*p3*OneMinusxi2 + B*xi2OveraStar  ) / OneMinusxi2Times64;   // A3  x^3 Coeff
    Coeff[4] = (64.0*aa + 21.0*LGM_EPS*p2*OneMinusxi2 + 8.0*a*s*p1 + A*xi2OveraStar ) / OneMinusxi2Times64;              // A2  x^4 Coeff
    Coeff[5] = ( 128.0*a + 4.0*s*OneMinusxi2*p1 ) / OneMinusxi2Times64;                                                  // A1  x^5 Coeff
    Coeff[6] = 1.0;                                                                                                      // A0  x^6 Coeff

    /*
     * Solve for resonant roots.
     */
    gsl_poly_complex_workspace *w = gsl_poly_complex_workspace_alloc( 7 );
    gsl_err = gsl_poly_complex_solve( Coeff, 7, w, gsl_z );
    //if (gsl_err != GSL_SUCCESS ) printf("CRAP!\n");
    gsl_poly_complex_workspace_free( w );
    for (i=0; i<6; i++) zz[i] = gsl_z[2*i] + gsl_z[2*i+1]*I;


    /*
     * Gather applicable roots together into the z[] array
     */
    nRoots = 0;
    for ( i=0; i<6; i++ ) {
        if ( ( fabs(cimag(zz[i])) < 1e-10 ) && ( creal(zz[i]) > 0.0 ) ) z[nRoots++] = zz[i];
    }


    // This is another hard-wired input parameter. Fix.
    int sum_res = 0;

    R  = dBoverB2;   // The ratio (dB/B)^2

    if ( ( nRoots == 0 ) && ( Mu2 > 1e-16 ) ){

        /*
         *  No resonant roots were found
         */
        Daa = 0.0;


    } else if ( Mu2 <= 1e-16 ) {

        /*
         *  If cos(alpha) is zero, then the resonance condition becomnes x = -a = s*Lambda/Gamma.
         *  And x>0 in this case will only happen for (R-MODE with electrons) and (L-MODE with protons)
         */
        if          ( (s ==  1) && ( Lambda < 0.0 ) ) {  // R-MODE with electrons
            x0 = 1.0/Gamma;
        } else if   ( (s == -1) && ( Lambda > 0.0 ) ) { // L-MODE with protons
            x0 = Lambda/Gamma;
        } else {
            x0 = 0.0;
        }

        if ( x0 > 1e-12 ) {
            y0 = x0 * sqrt( 1.0 + ( 1.0/(s-x0) - LGM_EPS*n1/(x0+s*LGM_EPS) - LGM_EPS*n2/(4.0*x0+s*LGM_EPS) - LGM_EPS*n3/(16.0*x0+s*LGM_EPS) )/(aStar*x0) );
            Ep1 = E+1.0; Ep12 = Ep1*Ep1; arg = (x0-xm)/dx;
            DD = M_PI_2/Rho * Omega_Sig*Omega_Sig/fabs(Omega_e) * R/(Ep12*dx) * exp( -arg*arg );
            Daa = ( sum_res ) ? 2.0*DD : DD;
        } else {
            Daa = 0.0;
        }

    } else {

        /*
         *  We have resonant roots and we are not at 90Deg.
         */
        q3 = 176.0*n1 + 107.0*n2 + 11.0*n3;
        q4 = 256.0*n1 +  16.0*n2 +    n3;
        q5 = 320.0*n1 +  68.0*n2 +  5.0*n3;

        g0 = 2.0*LGM_EPS5*( 21.0 - q2 + LGM_EPS*aStar );
        g1 = LGM_EPS4*s*( 3.0*(203.0 - q3) - 4.0*LGM_EPS2*aStar + 4.0*LGM_EPS*(21.0*aStar + q2) );
        g2 = 2.0*LGM_EPS3*( aStar*LGM_EPS3 + 4.0*(457.0 - q5) + LGM_EPS2*(-84.0*aStar - q2) + 3.0*LGM_EPS*(203.0*aStar + q3) );
        g3 = LGM_EPS2*s*( 84.0*aStar*LGM_EPS3 + 16.0*(609.0 - q4) + 16.0*LGM_EPS*(457.0*aStar + q5) + 3.0*LGM_EPS2*(-812.0*aStar - q3 ) );
        g4 = LGM_EPS*( 10752.0 + 1218.0*aStar*LGM_EPS3 + 32.0*LGM_EPS*(609.0*aStar + q4) - 8.0*LGM_EPS2*(1828.0*aStar + q5) );
        g5 = 16.0*s*( 256.0 + 1344.0*aStar*LGM_EPS + 457.0*aStar*LGM_EPS3  + LGM_EPS2*(-2436.0*aStar - q4 ) );
        g6 = 32.0*aStar*( 256.0 - 1344.0*LGM_EPS + 609.0*LGM_EPS2 );
        g7 = 1024.0*aStar*s*p1;
        g8 = 8192.0*aStar;

        Ep1 = E+1.0; Ep12 = Ep1*Ep1;
        fac = M_PI_2/Rho * Omega_Sig*Omega_Sig/fabs(Omega_e) * R/Ep12;

        /*
         * Execute sum over resonant roots.
         */
        BetaMu = sqrt(xi2);
        Beta   = sqrt(Beta2);
        Mu     = sqrt(Mu2);
        for ( Daa=0.0, n=0; n<nRoots; n++ ){

            x = creal( z[n] ); x2 = x*x; x3 = x2*x; x4 = x2*x2; x5 = x3*x2; x6 = x3*x3; x7 = x4*x3; x8 = x4*x4;
            y = ( x+a )/BetaMu;
            if ((x>xl)&&(x<xh)) {
                g = g0 + g1*x + g2*x2 + g3*x3 + g4*x4 + g5*x5 + g6*x6 + g7*x7 + g8*x8;
                sss = (x - s)*(x + s*LGM_EPS)*(4.0*x + s*LGM_EPS)*(16.0*x + s*LGM_EPS);

                F = 2.0*y*aStar*sss*sss/(x*g);

                u = 1.0 - x*Mu/(y*Beta);
                arg = (x-xm)/dx;
                Daa += u*u *fabs(F) / ( dx * fabs( BetaMu - F ) ) * exp( -arg*arg );
            }

        }
        Daa *= fac;



    }

    return(Daa);

}

/**
 *  \brief
 *      Computes the local Summer's [2007] Dap diffusion coefficient.
 *
 *  \details
 *
 *
 *
 *
 *      \param[in]      SinAlpha2   \f$\sin^2(\alpha)\f$, where \f$\alpha\f$ is the local particle pitch angle.
 *      \param[in]      E           Dimensionless energy Ek/E0 (kinetic energy over rest mass).
 *      \param[in]      dBoverB2    Ratio of wave amplitude, dB to local background field, B.
 *      \param[in]      BoverBeq    Ratio of local B to Beq.
 *      \param[in]      Omega_e     Local gyro-frequency of electrons.
 *      \param[in]      Omega_Sig   Local gyro-frequency of particle species we are interested in.
 *      \param[in]      Rho         This is \f$0.5 \sqrt(\pi) ( \mbox{erf}((wm-w1)/dw) + \mbox{erf}(wm-w1)/dw) )\f$.
 *      \param[in]      Sig         Blah.
 *      \param[in]      xm          Blah.
 *      \param[in]      dx          Blah.
 *      \param[in]      Lambda      Particle species. Can be LGM_ELECTRONS or LGM_PROTONS.
 *      \param[in]      s           Mode of the wave. Can be LGM_R_MODE_WAVE ( s = -1 ) or LGM_L_MODE_WAVE ( s = +1 ).
 *      \param[in]      n1          Ratio of H+ Number Density to Electron Number Density. (n1+n2+n3 must add to 1).
 *      \param[in]      n2          Ratio of He+ Number Density to Electron Number Density. (n1+n2+n3 must add to 1).
 *      \param[in]      n3          Ratio of O+ Number Density to Electron Number Density. (n1+n2+n3 must add to 1).
 *      \param[in]      aStar       Local value of the aStar parameter ( \f$ \alpha^* \f$ ) in the Summer's papers.
 *
 *      \return         Local value of Dap
 *
 *      \author         M. Henderson
 *      \date           2010-2011
 *
 */
double Lgm_SummersDapLocal_2007( double SinAlpha2, double E, double dBoverB2, double BoverBeq, double Omega_e, double Omega_Sig, double Rho, double Sig, double xl, double xh, double xm, double dx, double Lambda, int s, double n1, double n2, double n3, double aStar ) {

    int             nRoots, n, i;
    double          q0, q2, q3, q4, q5;
    double          p1, p2, p3, p4;
    double          Gamma, Gamma2, Beta, Beta2, Mu, Mu2, BetaMu;
    double          A, B, C, a, aa, apa, Coeff[7];
    double complex  z[7], zz[7];
    double          gsl_z[20];
    double          Dap, R, x0, y0, Ep1, Ep12, sss, u, arg, SinAlpha;
    double          x, x2, x3, x4, x5, x6, x7, x8;
    double          g, g0, g1, g2, g3, g4, g5, g6, g7, g8;
    double          xi2, xi2OveraStar, OneMinusxi2, OneMinusxi2Times64;
    double          y, F, DD, fac;

    /*
     * Solve for resonant roots.
     */
    Gamma = E+1.0; Gamma2 = Gamma*Gamma;
    Beta2 = E*(E+2.0) / Gamma2;
    Mu2   = 1.0-SinAlpha2; if (Mu2<0.0) Mu2 = 0.0; // Mu2 is cos^2(Alpha)
    xi2  = Beta2*Mu2;
    OneMinusxi2 = (1.0 - xi2);
    OneMinusxi2Times64 = 64.0*OneMinusxi2;
    a  = s*Lambda/Gamma; aa = a*a; apa = a+a;

    q0 = 16.0*n1 + 4.0*n2 + n3;
    q2 = 20.0*n1 + 17.0*n2 + 5.0*n3;

    p1 = 21.0*LGM_EPS - 16.0;
    p2 = LGM_EPS - 4.0;
    p3 = LGM_EPS - 21.0;
    p4 = LGM_EPS + 21.0;

    A = 64.0 + 4.0*LGM_EPS*q0;
    B = LGM_EPS*s*( 4.0*(21.0-q0) + LGM_EPS*q2 );
    C = LGM_EPS2*( p4 - q2 );

    xi2OveraStar = xi2/aStar;

    Coeff[0] = -aa*LGM_EPS3 / OneMinusxi2Times64;                                                                        // A6  x^0 Coeff
    Coeff[1] = a*LGM_EPS2*( a*s*p3 - 2.0*LGM_EPS ) / OneMinusxi2Times64;                                                 // A5  x^1 Coeff
    Coeff[2] = ( 21.0*aa*LGM_EPS*p2 + 2.0*a*LGM_EPS2*s*p3 - LGM_EPS3*OneMinusxi2 + C*xi2OveraStar) / OneMinusxi2Times64; // A4  x^2 Coeff
    Coeff[3] = ( 42.0*a*LGM_EPS*p2 + 4.0*aa*s*p1 + LGM_EPS2*s*p3*OneMinusxi2 + B*xi2OveraStar  ) / OneMinusxi2Times64;   // A3  x^3 Coeff
    Coeff[4] = (64.0*aa + 21.0*LGM_EPS*p2*OneMinusxi2 + 8.0*a*s*p1 + A*xi2OveraStar ) / OneMinusxi2Times64;              // A2  x^4 Coeff
    Coeff[5] = ( 128.0*a + 4.0*s*OneMinusxi2*p1 ) / OneMinusxi2Times64;                                                  // A1  x^5 Coeff
    Coeff[6] = 1.0;                                                                                                      // A0  x^6 Coeff

    /*
     * Solve for resonant roots.
     */
    gsl_poly_complex_workspace *w = gsl_poly_complex_workspace_alloc( 7 );
    gsl_poly_complex_solve( Coeff, 7, w, gsl_z );
    gsl_poly_complex_workspace_free( w );
    for (i=0; i<6; i++) zz[i] = gsl_z[2*i] + gsl_z[2*i+1]*I;

    /*
     * Gather applicable roots together into the z[] array
     */
    nRoots = 0;
    for ( i=0; i<6; i++ ) {
        if ( ( fabs(cimag(zz[i])) < 1e-10 ) && ( creal(zz[i]) > 0.0 ) ) z[nRoots++] = zz[i];
    }


    R  = dBoverB2;   // The ratio (dB/B)^2


    // This is another hard-wired input parameter. Fix.
    int sum_res = 0;


    if ( ( nRoots == 0 ) && ( Mu2 > 1e-16 ) ){

        /*
         *  No resonant roots were found
         */
        Dap = 0.0;


    } else if ( Mu2 <= 1e-16 ) {

        /*
         *  If cos(alpha) is zero, then the resonance condition becomnes x = -a = s*Lambda/Gamma.
         *  And x>0 in this case will only happen for (R-MODE with electrons) and (L-MODE with protons)
         */
        if          ( (s ==  1) && ( Lambda < 0.0 ) ) {  // R-MODE with electrons
            x0 = 1.0/Gamma;
        } else if   ( (s == -1) && ( Lambda > 0.0 ) ) { // L-MODE with protons
            x0 = Lambda/Gamma;
        } else {
            x0 = 0.0;
        }

        if ( x0 > 1e-12 ) {
            y0 = x0 * sqrt( 1.0 + ( 1.0/(s-x0) - LGM_EPS*n1/(x0+s*LGM_EPS) - LGM_EPS*n2/(4.0*x0+s*LGM_EPS) - LGM_EPS*n3/(16.0*x0+s*LGM_EPS) )/(aStar*x0) );
            Beta  = sqrt(Beta2);
            Ep1   = E+1.0; Ep12 = Ep1*Ep1;
            arg   = (x0-xm)/dx;
            DD = -M_PI_2/Rho * Omega_Sig*Omega_Sig/fabs(Omega_e) * R*x0/(Ep12*dx*y0*Beta)  * exp( -arg*arg );
            Dap   = ( sum_res ) ? 2.0*DD : DD;
        } else {
            Dap = 0.0;
        }

    } else {

        /*
         *  We have resonant roots and we are not at 90Deg.
         */
        q3 = 176.0*n1 + 107.0*n2 + 11.0*n3;
        q4 = 256.0*n1 +  16.0*n2 +    n3;
        q5 = 320.0*n1 +  68.0*n2 +  5.0*n3;

        g0 = 2.0*LGM_EPS5*( 21.0 - q2 + LGM_EPS*aStar );
        g1 = LGM_EPS4*s*( 3.0*(203.0 - q3) - 4.0*LGM_EPS2*aStar + 4.0*LGM_EPS*(21.0*aStar + q2) );
        g2 = 2.0*LGM_EPS3*( aStar*LGM_EPS3 + 4.0*(457.0 - q5) + LGM_EPS2*(-84.0*aStar - q2) + 3.0*LGM_EPS*(203.0*aStar + q3) );
        g3 = LGM_EPS2*s*( 84.0*aStar*LGM_EPS3 + 16.0*(609.0 - q4) + 16.0*LGM_EPS*(457.0*aStar + q5) + 3.0*LGM_EPS2*(-812.0*aStar - q3 ) );
        g4 = LGM_EPS*( 10752.0 + 1218.0*aStar*LGM_EPS3 + 32.0*LGM_EPS*(609.0*aStar + q4) - 8.0*LGM_EPS2*(1828.0*aStar + q5) );
        g5 = 16.0*s*( 256.0 + 1344.0*aStar*LGM_EPS + 457.0*aStar*LGM_EPS3  + LGM_EPS2*(-2436.0*aStar - q4 ) );
        g6 = 32.0*aStar*( 256.0 - 1344.0*LGM_EPS + 609.0*LGM_EPS2 );
        g7 = 1024.0*aStar*s*p1;
        g8 = 8192.0*aStar;


        /*
         * Execute sum over resonant roots.
         */
        Beta   = sqrt(Beta2);
        if (SinAlpha2 < 0.0) SinAlpha = 0.0;
        else if ( SinAlpha2 > 1.0) SinAlpha = 1.0;
        else SinAlpha = sqrt( SinAlpha2 );

        Ep1 = E+1.0; Ep12 = Ep1*Ep1;
        fac = -M_PI_2/Rho * Omega_Sig*Omega_Sig/fabs(Omega_e) * R*SinAlpha/Ep12 /Beta;
        BetaMu = sqrt(xi2);
        Mu     = sqrt(Mu2);

        for ( Dap=0.0, n=0; n<nRoots; n++ ){

            x = creal( z[n] ); x2 = x*x; x3 = x2*x; x4 = x2*x2; x5 = x3*x2; x6 = x3*x3; x7 = x4*x3; x8 = x4*x4;
            y = ( x+a )/BetaMu;
            if ((x>xl)&&(x<xh)) {
                g = g0 + g1*x + g2*x2 + g3*x3 + g4*x4 + g5*x5 + g6*x6 + g7*x7 + g8*x8;
                sss = (x - s)*(x + s*LGM_EPS)*(4.0*x + s*LGM_EPS)*(16.0*x + s*LGM_EPS);

                F = 2.0*y*aStar*sss*sss/(x*g);

                u = 1.0 - x*Mu/(y*Beta);
                arg = (x-xm)/dx;
                // if fabs(BetaMu - F) gets too small, we will get a singularity.
                Dap += x/y * u *fabs(F) / ( dx * fabs( BetaMu - F ) ) * exp( -arg*arg );
            }

        }
        Dap *= fac;

    }

    return(Dap);

}

/**
 *  \brief
 *      Computes the local Summer's [2007] Dpp diffusion coefficient.
 *
 *  \details
 *
 *
 *
 *
 *      \param[in]      SinAlpha2   \f$\sin^2(\alpha)\f$, where \f$\alpha\f$ is the local particle pitch angle.
 *      \param[in]      E           Dimensionless energy Ek/E0 (kinetic energy over rest mass).
 *      \param[in]      dBoverB2    Ratio of wave amplitude, dB to local background field, B.
 *      \param[in]      BoverBeq    Ratio of local B to Beq.
 *      \param[in]      Omega_e     Local gyro-frequency of electrons.
 *      \param[in]      Omega_Sig   Local gyro-frequency of particle species we are interested in.
 *      \param[in]      Rho         This is \f$0.5 \sqrt(\pi) ( \mbox{erf}((wm-w1)/dw) + \mbox{erf}(wm-w1)/dw) )\f$.
 *      \param[in]      Sig         Blah.
 *      \param[in]      xm          Blah.
 *      \param[in]      dx          Blah.
 *      \param[in]      Lambda      Particle species. Can be LGM_ELECTRONS or LGM_PROTONS.
 *      \param[in]      s           Mode of the wave. Can be LGM_R_MODE_WAVE ( s = -1 ) or LGM_L_MODE_WAVE ( s = +1 ).
 *      \param[in]      n1          Ratio of H+ Number Density to Electron Number Density. (n1+n2+n3 must add to 1).
 *      \param[in]      n2          Ratio of He+ Number Density to Electron Number Density. (n1+n2+n3 must add to 1).
 *      \param[in]      n3          Ratio of O+ Number Density to Electron Number Density. (n1+n2+n3 must add to 1).
 *      \param[in]      aStar       Local value of the aStar parameter ( \f$ \alpha^* \f$ ) in the Summer's papers.
 *
 *      \return         Local value of Dpp
 *
 *      \author         M. Henderson
 *      \date           2010-2011
 *
 */
double Lgm_SummersDppLocal_2007( double SinAlpha2, double E, double dBoverB2, double BoverBeq, double Omega_e, double Omega_Sig, double Rho, double Sig, double xl, double xh, double xm, double dx, double Lambda, int s, double n1, double n2, double n3, double aStar ) {

    int             nRoots, n, i;
    double          q0, q2, q3, q4, q5;
    double          p1, p2, p3, p4;
    double          Gamma, Gamma2, Beta, Beta2, Mu, Mu2, BetaMu;
    double          A, B, C, a, aa, apa, Coeff[7];
    double complex  z[6], zz[6];
    double          gsl_z[20];
    double          Dpp, R, x0, y0, Ep1, Ep12, sss, u, arg;
    double          x, x2, x3, x4, x5, x6, x7, x8;
    double          g, g0, g1, g2, g3, g4, g5, g6, g7, g8;
    double          xi2, xi2OveraStar, OneMinusxi2, OneMinusxi2Times64;
    double          y, F, DD, fac;

    /*
     * Solve for resonant roots.
     */
    Gamma = E+1.0; Gamma2 = Gamma*Gamma;
    Beta2 = E*(E+2.0) / Gamma2;
    Mu2   = 1.0-SinAlpha2; if (Mu2<0.0) Mu2 = 0.0; // Mu2 is cos^2(Alpha)
    xi2  = Beta2*Mu2;
    OneMinusxi2 = (1.0 - xi2);
    OneMinusxi2Times64 = 64.0*OneMinusxi2;
    a  = s*Lambda/Gamma; aa = a*a; apa = a+a;

    q0 = 16.0*n1 + 4.0*n2 + n3;
    q2 = 20.0*n1 + 17.0*n2 + 5.0*n3;

    p1 = 21.0*LGM_EPS - 16.0;
    p2 = LGM_EPS - 4.0;
    p3 = LGM_EPS - 21.0;
    p4 = LGM_EPS + 21.0;

    A = 64.0 + 4.0*LGM_EPS*q0;
    B = LGM_EPS*s*( 4.0*(21.0-q0) + LGM_EPS*q2 );
    C = LGM_EPS2*( p4 - q2 );

    xi2OveraStar = xi2/aStar;

    Coeff[0] = -aa*LGM_EPS3 / OneMinusxi2Times64;                                                                        // A6  x^0 Coeff
    Coeff[1] = a*LGM_EPS2*( a*s*p3 - 2.0*LGM_EPS ) / OneMinusxi2Times64;                                                 // A5  x^1 Coeff
    Coeff[2] = ( 21.0*aa*LGM_EPS*p2 + 2.0*a*LGM_EPS2*s*p3 - LGM_EPS3*OneMinusxi2 + C*xi2OveraStar) / OneMinusxi2Times64; // A4  x^2 Coeff
    Coeff[3] = ( 42.0*a*LGM_EPS*p2 + 4.0*aa*s*p1 + LGM_EPS2*s*p3*OneMinusxi2 + B*xi2OveraStar  ) / OneMinusxi2Times64;   // A3  x^3 Coeff
    Coeff[4] = (64.0*aa + 21.0*LGM_EPS*p2*OneMinusxi2 + 8.0*a*s*p1 + A*xi2OveraStar ) / OneMinusxi2Times64;              // A2  x^4 Coeff
    Coeff[5] = ( 128.0*a + 4.0*s*OneMinusxi2*p1 ) / OneMinusxi2Times64;                                                  // A1  x^5 Coeff
    Coeff[6] = 1.0;                                                                                                      // A0  x^6 Coeff

    /*
     * Solve for resonant roots.
     */
    gsl_poly_complex_workspace *w = gsl_poly_complex_workspace_alloc( 7 );
    gsl_poly_complex_solve( Coeff, 7, w, gsl_z );
    gsl_poly_complex_workspace_free( w );
    for (i=0; i<6; i++) zz[i] = gsl_z[2*i] + gsl_z[2*i+1]*I;


    /*
     * Gather applicable roots together into the z[] array
     */
    nRoots = 0;
    for ( i=0; i<6; i++ ) {
        if ( ( fabs(cimag(zz[i])) < 1e-10 ) && ( creal(zz[i]) ) > 0.0  ) z[nRoots++] = zz[i];
    }


    // This is another hard-wired input parameter. Fix.
    int sum_res = 0;

    R  = dBoverB2;   // The ratio (dB/B)^2

    if ( ( nRoots == 0 ) && ( Mu2 > 1e-16 ) ){

        /*
         *  No resonant roots were found
         */
        Dpp = 0.0;


    } else if ( Mu2 <= 1e-16 ) {

        /*
         *  If cos(alpha) is zero, then the resonance condition becomnes x = -a = s*Lambda/Gamma.
         *  And x>0 in this case will only happen for (R-MODE with electrons) and (L-MODE with protons)
         */
        if          ( (s ==  1) && ( Lambda < 0.0 ) ) {  // R-MODE with electrons
            x0 = 1.0/Gamma;
        } else if   ( (s == -1) && ( Lambda > 0.0 ) ) { // L-MODE with protons
            x0 = Lambda/Gamma;
        } else {
            x0 = 0.0;
        }

        if ( x0 > 1e-12 ) {
            y0 = x0 * sqrt( 1.0 + ( 1.0/(s-x0) - LGM_EPS*n1/(x0+s*LGM_EPS) - LGM_EPS*n2/(4.0*x0+s*LGM_EPS) - LGM_EPS*n3/(16.0*x0+s*LGM_EPS) )/(aStar*x0) );
            Ep1   = E+1.0; Ep12 = Ep1*Ep1;
            arg   = (x0-xm)/dx;
            DD = M_PI_2/Rho * Omega_Sig*Omega_Sig/fabs(Omega_e) * R*x0*x0/(Ep12*dx*y0*y0*Beta2)  * exp( -arg*arg );
            Dpp   = ( sum_res ) ? 2.0*DD : DD;
        } else {
            Dpp = 0.0;
        }

    } else {

        /*
         *  We have resonant roots and we are not at 90Deg.
         *  These could be cached.
         */
        q3 = 176.0*n1 + 107.0*n2 + 11.0*n3;
        q4 = 256.0*n1 +  16.0*n2 +      n3;
        q5 = 320.0*n1 +  68.0*n2 +  5.0*n3;

        /*
         * A bunch of these terms could be cached or pre-computed.
         */
        g0 = 2.0*LGM_EPS5*( 21.0 - q2 + LGM_EPS*aStar );
        g1 = LGM_EPS4*s*( 3.0*(203.0 - q3) - 4.0*LGM_EPS2*aStar + 4.0*LGM_EPS*(21.0*aStar + q2) );
        g2 = 2.0*LGM_EPS3*( aStar*LGM_EPS3 + 4.0*(457.0 - q5) + LGM_EPS2*(-84.0*aStar - q2) + 3.0*LGM_EPS*(203.0*aStar + q3) );
        g3 = LGM_EPS2*s*( 84.0*aStar*LGM_EPS3 + 16.0*(609.0 - q4) + 16.0*LGM_EPS*(457.0*aStar + q5) + 3.0*LGM_EPS2*(-812.0*aStar - q3 ) );
        g4 = LGM_EPS*( 10752.0 + 1218.0*aStar*LGM_EPS3 + 32.0*LGM_EPS*(609.0*aStar + q4) - 8.0*LGM_EPS2*(1828.0*aStar + q5) );
        g5 = 16.0*s*( 256.0 + 1344.0*aStar*LGM_EPS + 457.0*aStar*LGM_EPS3  + LGM_EPS2*(-2436.0*aStar - q4 ) );
        g6 = 32.0*aStar*( 256.0 - 1344.0*LGM_EPS + 609.0*LGM_EPS2 );
        g7 = 1024.0*aStar*s*p1;
        g8 = 8192.0*aStar;

        Ep1 = E+1.0; Ep12 = Ep1*Ep1;
        fac = M_PI_2/Rho * Omega_Sig*Omega_Sig/fabs(Omega_e) * R*SinAlpha2/Ep12 /Beta2;

        /*
         * Execute sum over resonant roots.
         */
        Beta   = sqrt(Beta2);
        BetaMu = sqrt(xi2);
        Mu     = sqrt(Mu2);
        for ( Dpp=0.0, n=0; n<nRoots; n++ ){

            x = creal( z[n] ); x2 = x*x; x3 = x2*x; x4 = x2*x2; x5 = x3*x2; x6 = x3*x3; x7 = x4*x3; x8 = x4*x4;
            y = ( x+a )/BetaMu;
            if ((x>xl)&&(x<xh)) {
                g = g0 + g1*x + g2*x2 + g3*x3 + g4*x4 + g5*x5 + g6*x6 + g7*x7 + g8*x8;
                sss = (x - s)*(x + s*LGM_EPS)*(4.0*x + s*LGM_EPS)*(16.0*x + s*LGM_EPS);

                F = 2.0*y*aStar*sss*sss/(x*g);

                u = 1.0 - x*Mu/(y*Beta);
                arg = (x-xm)/dx;
                Dpp += x*x/(y*y) *fabs(F) / ( dx * fabs( BetaMu - F ) ) * exp( -arg*arg );
            }

        }
        Dpp *= fac;

    }

    return(Dpp);

}



/**
 *  \brief
 *      Finds singularities in the integrands
 *
 *  \details
 */
int Lgm_SummersFindSingularities( double  (*f)( double, _qpInfo *), _qpInfo *qpInfo, int Verbose, double Lat0, double Lat1, double *x, double *y ) {

    int     done, founda, foundb, foundc, Found;
    double  a, b, c, d, fa, fb, fc, fd, inc, D1, D2, DELTA;

    inc = 0.01*RadPerDeg;


    /*
     * Find a non-zero starting point
     */
    a = Lat0; done = FALSE; founda = FALSE;
    while ( !done ) {
        a += inc;
        if ( a > Lat1 ){
            done = TRUE;
        } else {
            fa = fabs( f( a, qpInfo ) );
            if ( fa > 0.0 ) {
                founda = TRUE;
                done   = TRUE;
            }
        }
    }


    /*
     * Keep stepping to get a higher point
     */
    b = a; done = FALSE; foundb = FALSE;
    while ( !done ) {
        b += inc;
        if ( b > Lat1 ){
            done = TRUE;
        } else {
            fb = fabs( f( b, qpInfo ) );
            if ( fb > fa ) {
                foundb = TRUE;
                done   = TRUE;
            }
        }
    }

    /*
     * Keep stepping to get a lower point than b
     */
    c = b; done = FALSE; foundc = FALSE;
    while ( !done ) {
        c += inc;
        if ( c > Lat1 ){
            done = TRUE;
        } else {
            fc = fabs( f( c, qpInfo ) );
            if ( fc < fb ) {
                foundc = TRUE;
                done   = TRUE;
            }
        }
    }

    //printf("Found Bracket: %g %g %g %g %g %g\n", a, b, c, fa, fb, fc );


    done = FALSE;
    while (!done) {

        D1 = b-a;
        D2 = c-b;

        if ( fabs(c-a) < 1e-14 ) {
            done = TRUE;
        } else {

            // bisect largest interval
            if ( D1 > D2 ) {
                d = a+0.5*D1; fd = fabs( f( d, qpInfo ) );

                if ( (fd>fb) ) {
                    b = d; fb = fd;
                } else {
                    a = d; fa = fd;
                }
            } else {
                d = b+0.5*D2; fd = fabs( f( d, qpInfo ) );

                if ( (fd>fb) ) {
                    b = d; fb = fd;
                } else {
                    c = d; fc = fd;
                }
            }
            //printf("a, b, c, d = %g %g %g %g\n", a, b, c, d);

        }

    }

    *x = b;
    *y = fb;

    if ( isinf( fb ) ) {
        if ( ( *x <= Lat0 ) || ( *x >= Lat1 )  ) {
            // singularity at endpoint (dqagp already deals with this case).
            Found = FALSE;
        } else {
            Found = TRUE;
        }
    } else {
        DELTA = fb - fabs(f( b-1e-4, qpInfo ));
        Found = ( DELTA > 1.0 ) ?  1 : 0;
    }

    if ( Found ) {
        if (Verbose) printf("Singularity at Lat = %g  val = %g\n", (*x)*DegPerRad, *y);
        return(1);
    } else {
        return(0);
    }

}



/**
 *  \brief
 *      Finds where integrand is non-zero on eaither end.
 *
 *  \details
 */
int Lgm_SummersFindCutoffs( double  (*f)( double, _qpInfo *), _qpInfo *qpInfo, int Verbose, double inc, double Lat0, double Lat1, double *a, double *b ) {

    int     done, founda, foundb;
    double  x, x0, x1, fx, fx0, fx1, D;

    *a = Lat0;
    *b = Lat1;

    founda = FALSE;
    foundb = FALSE;


    // convert inc from deg to rads.
    inc *= RadPerDeg;

    fx = fabs( f( Lat0, qpInfo ) );
    if ( fx != 0.0 ) {

        *a = Lat0;
        founda = TRUE;

    } else {

        /*
         * Find a non-zero starting point
         */
        x0 = x1 = Lat0; done = FALSE;
        while ( !done ) {
            if ( x1 < Lat1 ){
                fx1 = fabs( f( x1, qpInfo ) );
                if ( fx1 > 0.0 ) {
                    done = TRUE;
                    founda = TRUE;
                } else {
                    x0 = x1; fx0 = fx1;
                    x1 += inc;
                }
            } else {
                done = TRUE;
            }
        }

        if ( founda ) {
            /*
             * Refine the point
             */
            done = FALSE;
            while( !done ) {
                D = x1-x0;
                if ( fabs(D) < 1e-8 ) {
                    done = TRUE;
                } else {
                    x = x0 +0.5*D;
                    fx = fabs( f( x, qpInfo ) );
                    if ( fx > 0.0 ) {
                        x1 = x; fx1 = fx;
                    } else {
                        x0 = x; fx0 = fx;
                    }
                }
            }
            *a = x1;
        } else {
            /*
             * We scanned the whole interval and found nothing. We may need to
             * try again with a smaller inc.
             */
            *a = Lat1;
            return(0);
        }
    }


    fx = fabs( f( Lat1, qpInfo ) );
    if ( (fx != 0.0)  ) {

        *b = Lat1;
        foundb = TRUE;

    } else {

        /*
         * Find a non-zero ending point
         */
        x0 = x1 = Lat1; done = FALSE;
        while ( !done ) {
            if ( x0 > Lat0 ){
                fx0 = fabs( f( x0, qpInfo ) );
                if ( fx0 > 0.0 ) {
                    done = TRUE;
                    foundb = TRUE;
                } else {
                    x1 = x0; fx1 = fx0;
                    x0 -= inc;
                }
            } else {
                done = TRUE;
            }
        }

        if ( foundb ) {

            /*
             * Refine the point
             */
            done = FALSE;
            while( !done ) {
                D = x1-x0;
                if ( fabs(D) < 1e-8 ) {
                    done = TRUE;
                } else {
                    x = x0 + 0.5*D;
                    fx = fabs( f( x, qpInfo ) );
                    if ( fx > 0.0 ) {
                        x0 = x; fx0 = fx;
                    } else {
                        x1 = x; fx1 = fx;
                    }
                }
            }
            *b = x0;
        } else {
            /*
             * We scanned the whole interval and found nothing. We may need to
             * try again with a smaller inc. Can we ever get here?
             */
            *b = Lat0;
            return(0);
        }

    }

    return(1);


}


/**
 *  \brief
 *      Finds where integrand is non-zero on eaither end.
 *
 *  \details
 */
int Lgm_SummersFindCutoffs2( double  (*f)( double, _qpInfo *), _qpInfo *qpInfo, int Verbose, double Lat0, double Lat1, double *a, double *b ) {

    /*
     * Make a series of attempts at progressively finer increments.  Another
     * way to accellerate this would be to keep track of the cutoffs we already
     * have and use them as a guide for future searches.
     */
    if ( Lgm_SummersFindCutoffs( f, qpInfo, Verbose, 2.00, Lat0, Lat1, a, b ) == TRUE ) return(1);
    if ( Lgm_SummersFindCutoffs( f, qpInfo, Verbose, 0.10, Lat0, Lat1, a, b ) == TRUE ) return(1);
    if ( Lgm_SummersFindCutoffs( f, qpInfo, Verbose, 0.01, Lat0, Lat1, a, b ) == TRUE ) return(1);
    //if ( Lgm_SummersFindCutoffs( f, qpInfo, Verbose, 0.001, Lat0, Lat1, a, b ) == TRUE ) return(1);

    return(0);

}












