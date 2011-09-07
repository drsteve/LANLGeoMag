#include "Lgm/Lgm_SummersDiffCoeff.h"


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


    return( 1e-9*B*LGM_e/m );

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
 *  \details
 *
 *      \param[in]      Alpha0      Equatoria PA in Degrees.
 *      \param[in]      Ek          Kinetic energy in MeV.
 *      \param[in]      L           L-shell parameter (dimensionless).
 *      \param[in]      BwFuncData  A (void *) pointer to data that the user may need to use in BwFunc().
 *      \param[in]      BwFunc      A user-defined function that returns Bw as a function of Latitude in units of nT. Prototype is "double BwFunc( double Lat, void *Data )".
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
 *      \return         void
 *
 *      \author         M. Henderson
 *      \date           2011
 *
 */
int Lgm_SummersDxxBounceAvg( double Alpha0,  double Ek,  double L,  void *BwFuncData, double (*BwFunc)(), double aStarEq,  double w1, double w2, double wm, double dw, int WaveMode, int Species, double MaxWaveLat, double *Daa_ba,  double *Dap_ba,  double *Dpp_ba) {

    double           T0, T1, T, a, b, E, s, Lambda, E0, Omega_eEq, Omega_SigEq, Beq, Rho;
    double           epsabs, epsrel, result, abserr, resabs, resasc, work[2001], points[3];
    int              npts=4, key=6, limit=500, lenw=4*limit, iwork[501], last, ier, neval;
    Lgm_SummersInfo *si=(Lgm_SummersInfo *)calloc( 1, sizeof(Lgm_SummersInfo));

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
     *  Quadpack integration tolerances.
1e-10 seems exccessive.
     */
    epsabs = 1e-10;
    epsrel = 1e-3;


    /*
     *  Pack integrand parameters into structure
     */
    if ( WaveMode == LGM_R_MODE_WAVE ) {
        si->s = 1.0;
    } else if ( WaveMode == LGM_L_MODE_WAVE ) {
        si->s = -1.0;
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
    b = acos( Lgm_CdipMirrorLat( si->SinAlpha0 ) );     // radians
    if ( fabs(a-b) < 1e-9 ) {

        *Daa_ba = 0.0;
        *Dap_ba = 0.0;
        *Dpp_ba = 0.0;
        return( 0 );

    } else {


        /*
         *  Perform integrations. The points[] array contains a list of point in
         *  the integrand where integrable singularities occur. dqagp() will avoid
         *  computing integrand at exactly these points (dqagp() is an
         *  extrapolation algorithm so it handles singularities very well.)
         */
        points[1] = b; npts=2;
        dqagp( CdipIntegrand_Sb, (_qpInfo *)si, a, b, npts, points, epsabs, epsrel, &T, &abserr, &neval, &ier, limit, lenw, &last, iwork, work );
        dqagp( SummersIntegrand_Gaa, (_qpInfo *)si, a, b, npts, points, epsabs, epsrel, Daa_ba, &abserr, &neval, &ier, limit, lenw, &last, iwork, work );
        dqagp( SummersIntegrand_Gap, (_qpInfo *)si, a, b, npts, points, epsabs, epsrel, Dap_ba, &abserr, &neval, &ier, limit, lenw, &last, iwork, work );
        dqagp( SummersIntegrand_Gpp, (_qpInfo *)si, a, b, npts, points, epsabs, epsrel, Dpp_ba, &abserr, &neval, &ier, limit, lenw, &last, iwork, work );


        /*
         *  Calculate T( SinAlpha0 ), related to bounce period.
         *  Schultz and Lanzeroti Eqns 1.28a-d)
         */
        //T0 = 1.38017299815047317375;
        //T1 = 0.74048048969306104115;
        //printf("T = %g\n", T);
        //T  = T0 - 0.5*(T0-T1)*(si->SinAlpha0 + sqrt(si->SinAlpha0));
        //printf("T = %g\n", T);


        *Daa_ba /= T;
        *Dap_ba /= T;
        *Dpp_ba /= T;

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
    double  BoverBeq, BoverBm, SinAlpha2;
    double  SinAlpha02, Ek_in, L, aStarEq;
    double  f;
    Lgm_SummersInfo *si;

    /*
     * Pull parameters for integrand from Lgm_SummersInfo structure.
     */
    si  = (Lgm_SummersInfo *)qpInfo;
    SinAlpha02 = si->SinAlpha02;    // pre-computed sin^2( Alpha0 )

    CosLat  = cos( Lat );     CosLat2 = CosLat*CosLat;
    CosLat3 = CosLat2*CosLat; CosLat6 = CosLat3*CosLat3;
    arg = 4.0 - 3.0*CosLat2;
    v = (arg < 0.0 ) ? 0.0 : sqrt( arg );

    BoverBeq  = v/CosLat6;              // B/Beq
    BoverBm   = BoverBeq*SinAlpha02;    // B/Bm = B/Beq * sin^2(Alpha0)
    //SinAlpha2 = BoverBm;                // sin^2(Alpha) = B/Bm
    //if (SinAlpha2 > 1.0) { SinAlpha2 = BoverBm = 1.0; }

    /*
     * We dont need to worry about the possibility that 1-B/Bm is 0 which gives
     * f of infinity. This singularity is not a problem because it is
     * integrable and will never get evaluated when using the proper quadpack
     * integrator.
     */
    f = CosLat*v/sqrt(1.0-BoverBm);
//printf("Lat = %.15g    BoverBm = %.15g    f = %.15g\n", Lat*DegPerRad, BoverBm, f );

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
    double  BoverBeq, BoverBm, SinAlpha2, Omega_eEq, Omega_SigEq, w1, w2, wm, dw, MaxWaveLat;
    double  SinAlpha02, CosAlpha2, TanAlpha2, TanAlpha02, Ek_in, L, aStarEq;
    double  Daa, Da0a0, f, aStar, E, dB, B, dBoverB2, s, Lambda, Rho, Sig;
    Lgm_SummersInfo *si;


    /*
     * Pull parameters for integrand from Lgm_SummersInfo structure.
     */
    si  = (Lgm_SummersInfo *)qpInfo;
    MaxWaveLat  = si->MaxWaveLat;    // Latitudinal cutoff for waves.
    if ( fabs(Lat) > MaxWaveLat ) return(0.0);
    SinAlpha02  = si->SinAlpha02;    // pre-computed sin^2( Alpha0 )
    TanAlpha02  = si->TanAlpha02;    // pre-computed tan^2( Alpha0 )
    E           = si->E;             // Dimensionless energy Ek/E0 (kinetic energy over rest energy).
    L           = si->L;
    aStarEq     = si->aStarEq;
//    dB          = si->dB;
    dB          = si->BwFunc( Lat, si->BwFuncData );
    Omega_eEq   = si->Omega_eEq;     // Equatorial electron gyro-frequency.
    Omega_SigEq = si->Omega_SigEq;   // Equatorial gyro-frequency for given species.
    w1          = si->w1;            // lower cutoff frequency.
    w2          = si->w2;            // upper cutoff frequency.
    wm          = si->wm;            // frequency of max power.
    dw          = si->dw;            // frequency bandwidth. (Semi bandwidth is sigma*dw).
    Lambda      = si->Lambda;        //
    s           = si->s;             //
    Rho         = si->Rho;           // Rho = sqrt(PI)/2.0( erf( (wm-w1)/dw ) + erf( (w2-wm)/dw ) )   Eq (3) in Summers2007
    Sig         = si->Sig;           // Sig -- see Summers [2005], text above eqn (30).

    CosLat  = cos( Lat );     CosLat2 = CosLat*CosLat;
    CosLat3 = CosLat2*CosLat; CosLat6 = CosLat3*CosLat3;
    v = sqrt(4.0 - 3.0*CosLat2);

    BoverBeq  = v/CosLat6;              // B/Beq
    BoverBm   = BoverBeq*SinAlpha02;    // B/Bm = B/Beq * sin^2(Alpha0)
    SinAlpha2 = BoverBm;                // sin^2(Alpha) = B/Bm
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
    B         = BoverBeq*M_CDIP/(L*L*L);     // local value of B (warning, hard-coded M value).
    dBoverB2  = dB*dB/(B*B);                 // Local value of R = dB^2/B^2.
    aStar     = aStarEq*BoverBeq*BoverBeq;   // Local aStar value.
    Omega_e   = Lgm_GyroFreq( -LGM_e, B, LGM_ELECTRON_MASS );   // Local electron gyro-frequency.
    Omega_Sig = Omega_SigEq*BoverBeq*BoverBeq;                  // Local gyro-frequency for given species.
    x1        = w1/fabs(Omega_e);                               // Local value of x1.
    x2        = w2/fabs(Omega_e);                               // Local value of x1.
    xm        = wm/fabs(Omega_e);                               // Local value of xm.
    dx        = dw/fabs(Omega_e);                               // Local value of dx.


    /*
     * Compute the local Daa
     */
    Daa = Lgm_SummersDaaLocal( SinAlpha2, E, dBoverB2, BoverBeq, Omega_e, Omega_Sig, Rho, Sig, x1, x2, xm, dx, Lambda, s, aStar );


    /*
     * Compute Da0a0. The TanAlpha02/TanAlpha2 term (which is (dAlphaeq/dAlpha)^2 )
     * converts D_Alpha,Alpha to D_Alpha_eq,Alpha_eq.
     */
    Da0a0 = Daa*TanAlpha02/TanAlpha2;


    /*
     * Finally the integrand.
     */
    f = Da0a0 * CosLat*v/sqrt(1-BoverBm);

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
    double  BoverBeq, BoverBm, SinAlpha2, Omega_eEq, Omega_SigEq, w1, w2, wm, dw, MaxWaveLat;
    double  SinAlpha02, CosAlpha2, TanAlpha, TanAlpha0, Ek_in, L, aStarEq;
    double  Dap, Da0p, f, aStar, E, dB, B, dBoverB2, s, Lambda, Rho, Sig;
    Lgm_SummersInfo *si;


    /*
     * Pull parameters for integrand from Lgm_SummersInfo structure.
     */
    si  = (Lgm_SummersInfo *)qpInfo;
    MaxWaveLat  = si->MaxWaveLat;    // Latitudinal cutoff for waves.
    if ( fabs(Lat) > MaxWaveLat ) return(0.0); // Return 0 if beyond cutoff latitude.
    SinAlpha02  = si->SinAlpha02;    // pre-computed sin^2( Alpha0 )
    TanAlpha0   = si->TanAlpha0;     // pre-computed tan( Alpha0 )
    E           = si->E;             // Dimensionless energy Ek/E0 (kinetic energy over rest energy).
    L           = si->L;
    aStarEq     = si->aStarEq;
//    dB          = si->dB;
    dB          = si->BwFunc( Lat, si->BwFuncData );
    Omega_eEq   = si->Omega_eEq;     // Equatorial electron gyro-frequency.
    Omega_SigEq = si->Omega_SigEq;   // Equatorial gyro-frequency for given species.
    w1          = si->w1;            // lower cutoff frequency.
    w2          = si->w2;            // upper cutoff frequency.
    wm          = si->wm;            // frequency of max power.
    dw          = si->dw;            // frequency bandwidth. (Semi bandwidth is sigma*dw).
    Lambda      = si->Lambda;        //
    s           = si->s;             //
    Rho         = si->Rho;           // Rho = sqrt(PI)/2.0( erf( (wm-w1)/dw ) + erf( (w2-wm)/dw ) )   Eq (3) in Summers2007
    Sig         = si->Sig;           // Sig -- see Summers [2005], text above eqn (30).

    CosLat  = cos( Lat );     CosLat2 = CosLat*CosLat;
    CosLat3 = CosLat2*CosLat; CosLat6 = CosLat3*CosLat3;
    v = sqrt(4.0 - 3.0*CosLat2);

    BoverBeq  = v/CosLat6;              // B/Beq
    BoverBm   = BoverBeq*SinAlpha02;    // B/Bm = B/Beq * sin^2(Alpha0)
    SinAlpha2 = BoverBm;                // sin^2(Alpha) = B/Bm
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
    B         = BoverBeq*M_CDIP/(L*L*L);     // local value of B (warning, hard-coded M value).
    dBoverB2  = dB*dB/(B*B);                 // Local value of R = dB^2/B^2.
    aStar     = aStarEq*BoverBeq*BoverBeq;   // Local aStar value.
    Omega_e   = Lgm_GyroFreq( -LGM_e, B, LGM_ELECTRON_MASS );   // Local electron gyro-frequency.
    Omega_Sig = Omega_SigEq*BoverBeq*BoverBeq;                  // Local gyro-frequency for given species.
    x1        = w1/fabs(Omega_e);                               // Local value of x1.
    x2        = w2/fabs(Omega_e);                               // Local value of x1.
    xm        = wm/fabs(Omega_e);                               // Local value of xm.
    dx        = dw/fabs(Omega_e);                               // Local value of dx.

    /*
     * Compute the local Daa
     */
    Dap = Lgm_SummersDapLocal( SinAlpha2, E, dBoverB2, BoverBeq, Omega_e, Omega_Sig, Rho, Sig, x1, x2, xm, dx, Lambda, s, aStar );


    /*
     * Compute Da0p. The TanAlpha0/TanAlpha term (which is dAlphaeq/dAlpha)
     * converts D_Alpha,p to D_Alpha_eq,p.
     */
    Da0p = Dap*TanAlpha0/TanAlpha;


    /*
     * Finally the integrand.
     */
    f = Da0p * CosLat*v/sqrt(1-BoverBm);

//printf("Lat, f = %.15g %.15g\n", Lat*DegPerRad, log10(fabs(f))  );
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
    double  BoverBeq, BoverBm, SinAlpha2, Omega_eEq, Omega_SigEq, w1, w2, wm, dw, MaxWaveLat;
    double  SinAlpha02, CosAlpha2, TanAlpha2, TanAlpha02, Ek_in, L, aStarEq;
    double  Dpp, f, aStar, E, dB, B, dBoverB2, Lambda, s, Rho, Sig;
    Lgm_SummersInfo *si;


    /*
     * Pull parameters for integrand from Lgm_SummersInfo structure.
     */
    si  = (Lgm_SummersInfo *)qpInfo;
    MaxWaveLat  = si->MaxWaveLat;    // Latitudinal cutoff for waves.
    if ( fabs(Lat) > MaxWaveLat ) return(0.0);
    SinAlpha02  = si->SinAlpha02;    // pre-computed sin^2( Alpha0 )
    E           = si->E;             // Dimensionless energy Ek/E0 (kinetic energy over rest energy).
    L           = si->L;
    aStarEq     = si->aStarEq;
//    dB          = si->dB;
    dB          = si->BwFunc( Lat, si->BwFuncData );
    Omega_eEq   = si->Omega_eEq;     // Equatorial electron gyro-frequency.
    Omega_SigEq = si->Omega_SigEq;   // Equatorial gyro-frequency for given species.
    w1          = si->w1;            // lower cutoff frequency.
    w2          = si->w2;            // upper cutoff frequency.
    wm          = si->wm;            // frequency of max power.
    dw          = si->dw;            // frequency bandwidth. (Semi bandwidth is sigma*dw).
    Lambda      = si->Lambda;        //
    s           = si->s;             //
    Rho         = si->Rho;           // Rho = sqrt(PI)/2.0( erf( (wm-w1)/dw ) + erf( (w2-wm)/dw ) )   Eq (3) in Summers2007
    Sig         = si->Sig;           // Sig -- see Summers [2005], text above eqn (30).

    CosLat  = cos( Lat );     CosLat2 = CosLat*CosLat;
    CosLat3 = CosLat2*CosLat; CosLat6 = CosLat3*CosLat3;
    v = sqrt(4.0 - 3.0*CosLat2);

    BoverBeq  = v/CosLat6;              // B/Beq
    BoverBm   = BoverBeq*SinAlpha02;    // B/Bm = B/Beq * sin^2(Alpha0)
    SinAlpha2 = BoverBm;                // sin^2(Alpha) = B/Bm
    if (SinAlpha2 > 1.0) { SinAlpha2 = BoverBm = 1.0; }

    /*
     * Compute local parameters from the equatorial ones that are specified.
     * Compute local parameters from the equatorial ones that are specified.
     * Assume wave frequency distribution (i.e. defined by the 'w' quantities,
     * e.g. w1, w2, wm, dw) is constant along the field lines.
     */
    B         = BoverBeq*M_CDIP/(L*L*L);     // local value of B (warning, hard-coded M value).
    dBoverB2  = dB*dB/(B*B);                 // Local value of R = dB^2/B^2.
    aStar     = aStarEq*BoverBeq*BoverBeq;   // Local aStar value.
    Omega_e   = Lgm_GyroFreq( -LGM_e, B, LGM_ELECTRON_MASS );   // Local electron gyro-frequency.
    Omega_Sig = Omega_SigEq*BoverBeq*BoverBeq;                  // Local gyro-frequency for given species.
    x1        = w1/fabs(Omega_e);                               // Local value of x1.
    x2        = w2/fabs(Omega_e);                               // Local value of x1.
    xm        = wm/fabs(Omega_e);                               // Local value of xm.
    dx        = dw/fabs(Omega_e);                               // Local value of dx.

    /*
     * Compute the local Daa
     */
    Dpp = Lgm_SummersDppLocal( SinAlpha2, E, dBoverB2, BoverBeq, Omega_e, Omega_Sig, Rho, Sig, x1, x2, xm, dx, Lambda, s, aStar );


    /*
     * Finally the integrand.
     */
    f = Dpp * CosLat*v/sqrt(1-BoverBm);


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
double Lgm_SummersDaaLocal( double SinAlpha2, double E, double dBoverB2, double BoverBeq, double Omega_e, double Omega_Sig, double Rho, double Sig, double xl, double xh, double xm, double dx, double Lambda, double s, double aStar ) {

    int             nReal, nRoots, n;
    double          Gamma, Gamma2, Beta, Beta2, Mu, Mu2, BetaMu, BetaMu2, OneMinusBetaMu2;
    double          sEpsMinusOne, a1, a2, a3, a4, a, aa, apa, b;
    double complex  z1, z2, z3, z4, z[4];
    double          Daa, R, x0, y0, Ep1, Ep12, xms, xpse, u, arg, c1, c2, c3, c4;
    double          x, x2, x3, x4, y, g, F, Dcore, fac;

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


/*
double aaa[4];
int    nn=4;
double complex zz[4];
int num = 0;
aaa[0] = 1.0;
aaa[1] = a1;
aaa[2] = a2;
aaa[3] = a3;
aaa[4] = a4;
num = Lgm_PolyRoots( aaa, nn, zz );
*/




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
/*
int ii;
nRoots = 0;
for (ii=0; ii<num; ii++){
    if ( ( fabs(cimag(zz[ii])) < 1e-10 ) && ( fabs(creal(zz[ii])) > 0.0 ) ) z[nRoots++] = zz[ii];
}
*/


    if ( ( nRoots == 0 ) && ( Mu2 > 1e-16 ) ){

        /*
         *  No resonant roots were found
         */
        Daa = 0.0;


    } else if ( Mu2 <= 1e-16 ) {

        /*
         *  Alpha is essentially 90Deg. Use the limiting forms given by Eqns
         *  (36)-(38) of Summers2005. This is Eqn (36) for Daa.
         */
        x0 = 1.0/Gamma;
        y0 = x0 * sqrt( 1.0 + b*Gamma2/((Gamma-1.0)*(1.0+LGM_EPS*Gamma)) );

        Ep1 = E+1.0; Ep12 = Ep1*Ep1;
        arg = (x0-xm)/dx;
        Dcore = M_PI_2/Rho * Omega_Sig*Omega_Sig/fabs(Omega_e) * R/(Ep12*dx) * exp( -arg*arg );

        Daa = ( sum_res ) ? 2.0*Dcore : Dcore;

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
//if (y<0){
            g = x4 + c1*x3 + c2*x2 + c3*x + c4;
            xms  = x - s;
            xpse = x + s*LGM_EPS;
            F = y*xms*xms*xpse*xpse/(x*g);

            u = 1.0 - x*Mu/(y*Beta);
            arg = (x-xm)/dx;
            Daa += u*u *fabs(F) / ( dx * fabs( BetaMu - F ) ) * exp( -arg*arg );
//}
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
double Lgm_SummersDapLocal( double SinAlpha2, double E, double dBoverB2, double BoverBeq, double Omega_e, double Omega_Sig, double Rho, double Sig, double xl, double xh, double xm, double dx, double Lambda, double s, double aStar ) {

    int             nReal, nRoots, n;
    double          Gamma, Gamma2, Beta, Beta2, Mu, Mu2, BetaMu, BetaMu2, OneMinusBetaMu2;
    double          sEpsMinusOne, a1, a2, a3, a4, a, aa, apa, b, SinAlpha;
    double complex  z1, z2, z3, z4, z[4];
    double          Dap, R, x0, y0, Ep1, Ep12, xms, xpse, u, arg, c1, c2, c3, c4;
    double          x, x2, x3, x4, y, g, F, Dcore, fac;

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
        Dap = 0.0;


    } else if ( Mu2 <= 1e-16 ) {

        /*
         *  Alpha is essentially 90Deg. Use the limiting forms given by Eqns
         *  (36)-(38) of Summers2005. This is Eqn (37) for Dap/p.
CHECK THAT THIS IS USING (A2) and (A3) CORRECTLY.
         */
        x0 = 1.0/Gamma;
        y0 = x0 * sqrt( 1.0 + b*Gamma2/((Gamma-1.0)*(1.0+LGM_EPS*Gamma)) );

        Beta  = sqrt(Beta2);
        Ep1   = E+1.0; Ep12 = Ep1*Ep1;
        arg   = (x0-xm)/dx;
        Dcore = -M_PI_2/Rho * Omega_Sig*Omega_Sig/fabs(Omega_e) * R*x0/(Ep12*dx*y0*Beta)  * exp( -arg*arg );
        Dap   = ( sum_res ) ? 2.0*Dcore : Dcore;

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
//if (y<0){
            g = x4 + c1*x3 + c2*x2 + c3*x + c4;
            xms  = x - s;
            xpse = x + s*LGM_EPS;
            F = y*xms*xms*xpse*xpse/(x*g);

            u = 1.0 - x*Mu/(y*Beta);
            arg = (x-xm)/dx;
            Dap += x/y * u *fabs(F) / ( dx * fabs( BetaMu - F ) ) * exp( -arg*arg );
//}
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
double Lgm_SummersDppLocal( double SinAlpha2, double E, double dBoverB2, double BoverBeq, double Omega_e, double Omega_Sig, double Rho, double Sig, double xl, double xh, double xm, double dx, double Lambda, double s, double aStar ) {

    int             nReal, nRoots, n;
    double          Gamma, Gamma2, Beta, Beta2, Mu, Mu2, BetaMu, BetaMu2, OneMinusBetaMu2;
    double          sEpsMinusOne, a1, a2, a3, a4, a, aa, apa, b, SinAlpha;
    double complex  z1, z2, z3, z4, z[4];
    double          Dpp, R, x0, y0, Ep1, Ep12, xms, xpse, u, arg, c1, c2, c3, c4;
    double          x, x2, x3, x4, y, g, F, Dcore, fac;

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
    if ( ( fabs(cimag(z1)) < 1e-10 ) && ( creal(z1) > 0.0 ) ) z[nRoots++] = z1;
    if ( ( fabs(cimag(z2)) < 1e-10 ) && ( creal(z2) > 0.0 ) ) z[nRoots++] = z2;
    if ( ( fabs(cimag(z3)) < 1e-10 ) && ( creal(z3) > 0.0 ) ) z[nRoots++] = z3;
    if ( ( fabs(cimag(z4)) < 1e-10 ) && ( creal(z4) > 0.0 ) ) z[nRoots++] = z4;


    if ( ( nRoots == 0 ) && ( Mu2 > 1e-16 ) ){

        /*
         *  No resonant roots were found
         */
        Dpp = 0.0;


    } else if ( Mu2 <= 1e-16 ) {

        /*
         *  Alpha is essentially 90Deg. Use the limiting forms given by Eqns
         *  (36)-(38) of Summers2005. This is Eqn (37) for Dpp/p^2.
CHECK THAT THIS IS USING (A2) and (A3) CORRECTLY.
         */
        x0 = 1.0/Gamma;
        y0 = x0 * sqrt( 1.0 + b*Gamma2/((Gamma-1.0)*(1.0+LGM_EPS*Gamma)) );

        Ep1   = E+1.0; Ep12 = Ep1*Ep1;
        arg   = (x0-xm)/dx;
        Dcore = M_PI_2/Rho * Omega_Sig*Omega_Sig/fabs(Omega_e) * R*x0*x0/(Ep12*dx*y0*y0*Beta2)  * exp( -arg*arg );
        Dpp   = ( sum_res ) ? 2.0*Dcore : Dcore;

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
//if (y<0){
            g = x4 + c1*x3 + c2*x2 + c3*x + c4;
            xms  = x - s;
            xpse = x + s*LGM_EPS;
            F = y*xms*xms*xpse*xpse/(x*g);

            arg = (x-xm)/dx;
            Dpp += x*x/(y*y) *fabs(F) / ( dx * fabs( BetaMu - F ) ) * exp( -arg*arg );
//}
}

        }
        Dpp *= fac;



    }

    return(Dpp);

}
