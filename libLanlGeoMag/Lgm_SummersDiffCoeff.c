#include "Lgm/Lgm_SummersDiffCoeff.h"


/**
 *  \brief
 *      Compute the electron plasma frequency, /f$\omega_{pe}/f$.
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
 *      Compute the gyro-frequency of a particle with change q and mass m, /f$\Omega_\sigma/f$.
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
 *      Computes the bounce-averaged Summer's [2005] diffusion coefficients.
 *
 *  \details
 *
 *      \param[in]      Alpha0      equatoria PA in radians.
 *      \param[in]      Ek_in       kinetic energy in MeV.
 *      \param[in]      L           L-shell parameter (dimensionless).
 *      \param[in]      aStarEq     Equatorial value of the cold plasma parameter /f$\Omega_\sigma/\omega_{pe}/f$.
 *      \param[out]     Daa_ba      Bounce-averaged value of Daa.
 *      \param[out]     Dap_ba      Bounce-averaged value of Dap.
 *      \param[out]     Dpp_ba      Bounce-averaged value of Dpp.
 *
 *      \return         void
 *
 *      \author         M. Henderson
 *      \date           2010-2011
 *
 */
int Lgm_SummersDxxBounceAvg( double Alpha0,  double Ek,  double L,  double dB, double aStarEq,  double w1, double w2, double wm, double dw, int WaveMode, int Species, double *Daa_ba,  double *Dap_ba,  double *Dpp_ba) {

    double          T0, T1, T, a, b, E, s, Lambda, E0, Omega_eEq, Omega_SigEq, Beq, Rho;
    double          epsabs, epsrel, result, abserr, resabs, resasc, work[2001], points[3];
    int             npts=4, key=6, limit=500, lenw=4*limit, iwork[501], last, ier, neval;
    Lgm_SummersInfo *si=(Lgm_SummersInfo *)calloc( 1, sizeof(Lgm_SummersInfo));


    /*
     *  Set Rho. See Eqn(3) of Summers2007.
     *      Rho = sqrt(PI)/2.0( erf( (wm-w1)/dw ) + erf( (w2-wm)/dw ) )
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
        E0          = LGM_Ee0;      // set rest energy to proton rest energy (MeV)
        Omega_SigEq = Lgm_GyroFreq( -LGM_e, Beq, LGM_ELECTRON_MASS );
    } else if ( Species == LGM_PROTONS ) {
        si->Lambda = LGM_EPS;
        E0         = LGM_Ep0;       // set rest energy to electron rest energy (MeV)
        Omega_SigEq = Lgm_GyroFreq( LGM_e, Beq, LGM_PROTON_MASS );
    } else {
        printf("Lgm_SummersDxxBounceAvg(): Unknown species = %d\n", Species );
        return(0);
    }
    Omega_eEq = Lgm_GyroFreq( -LGM_e, Beq, LGM_ELECTRON_MASS );
    si->E           = Ek/E0;        // Dimensionless energy Ek/E0 (kinetic energy over rest energy).
    si->Omega_eEq   = Omega_eEq;    // Equatorial electron gyro-frequency.
    si->Omega_SigEq = Omega_SigEq;  // Equatorial gyro-frequency for given species.

    si->Alpha0      = Alpha0*M_PI/180.0;
    si->SinAlpha0   = sin(si->Alpha0);
    si->SinAlpha02  = si->SinAlpha0*si->SinAlpha0;
    si->CosAlpha0   = cos(si->Alpha0);
    si->CosAlpha02  = si->CosAlpha0*si->CosAlpha0;
    si->TanAlpha0   = si->SinAlpha0/si->CosAlpha0;
    si->TanAlpha02  = si->TanAlpha0*si->TanAlpha0;
    si->L           = L;
    si->aStarEq     = aStarEq;
    si->dB          = dB;
    si->wm          = wm;            // center frequency (at equator).
    si->dw          = dw;            // frequency bandwidth (at equator).
    si->Rho         = Rho;


    /*
     *  Set integration limits
     */
    a = 0.0;                                        // radians
    b = acos( Lgm_CdipMirrorLat( si->SinAlpha0 ) );  // radians
    //lambda_m = MIN( lambda_m, 25.*RadPerDeg );    // equatatorial confinement of waves.

    /*
     *  Perform integrations. The points[] array contains a list of point in
     *  the integrand where integrable singularities occur. dqagp() will avoid
     *  computing integrand at exactly these points (dqagp() is an
     *  extrapolation algorithm so it handles singularities very well.)
     */
    points[1] = a; points[2] = b; npts=2;
    dqagp( CdipIntegrand_Sb, (_qpInfo *)si, a, b, npts, points, epsabs, epsrel, &T, &abserr, &neval, &ier, limit, lenw, &last, iwork, work );
    //printf("neval = %d\n", neval);

    dqagp( SummersIntegrand_Gaa, (_qpInfo *)si, a, b, npts, points, epsabs, epsrel, Daa_ba, &abserr, &neval, &ier, limit, lenw, &last, iwork, work );
    //printf("neval = %d\n", neval);

    dqagp( SummersIntegrand_Gap, (_qpInfo *)si, a, b, npts, points, epsabs, epsrel, Dap_ba, &abserr, &neval, &ier, limit, lenw, &last, iwork, work );
    //printf("neval = %d\n", neval);

    dqagp( SummersIntegrand_Gpp, (_qpInfo *)si, a, b, npts, points, epsabs, epsrel, Dpp_ba, &abserr, &neval, &ier, limit, lenw, &last, iwork, work );
    //printf("neval = %d\n", neval);


    /*
     *  Calculate T( SinAlpha0 ), related to bounce period.
     *  Schultz and lanzeroti Eqns 1.28a-d)
     */
    //T0 = 1.38017299815047317375;
    //T1 = 0.74048048969306104115;
    //printf("T = %g\n", T);
    //T  = T0 - 0.5*(T0-T1)*(si->SinAlpha0 + sqrt(si->SinAlpha0));
    //printf("T = %g\n", T);


    *Daa_ba /= T;
    *Dap_ba /= T;
    *Dpp_ba /= T;


    return(1);

}





/**
 *  \brief
 *      Integrand for computing bounce-averaged Summer's [2005] Daa diffusion coefficient.
 *
 *  \details
 *
 *      The bounce average of a quantity \f$Q\f$ is;
 *
 *          \f[ <Q> = { 1\over S_b} \int_{s_{sm}}^{s_{nm}} {Q\;ds\over [1-{B/B_m}]^{1/2}} \f]
 *
 *      where,
 *          \f[ S_b = \int_{s_{sm}}^{s_{nm}} {ds\over [1-{B/B_m}]^{1/2}} \f]
 *
 *      and \f$s\f$ is the distance along the fieldl line and \f$s_{sm}\f$ and
 *      \f$s_{nm}\f$ are the southern and northern mirror points respectively.
 *
 *      The integrals can be recast as integrals over latitude using;
 *
 *          \f[ ds = L\cos\lambda[4-3\cos^2\lambda]^{1/2} d\lambda \f]
 *          \f[ B/B_m = { \sin^2\alpha_\circ [4-3\cos^2\lambda]^{1/2}\over cos^6\lambda } \f]
 *
 *      Since \f$L\f$ is constant and is in the numerator and denominator it cancels (so we can ignore it).
 *      The final result is;
 *
 *          \f[ <Q> = { 1\over S_b}\displaystyle \int_{\lambda_{sm}}^{\lambda_{nm}} {Q\;  \cos\lambda[4-3\cos^2\lambda]^{1/2} d\lambda  \over \left[1- {\displaystyle \sin^2\alpha_\circ [4-3\cos^2\lambda]^{1/2}\over\displaystyle cos^6\lambda }   \right]^{1/2}} \f]
 *
 *
 *
 *
 *      \param[in]      Lat         LAtitude in radians
 *      \param[in]      qpInfo      Structure contaning additional information needed to compute the integrand.
 *
 *      \return         The value of the integrand
 *
 *      \author         M. Henderson
 *      \date           2010-2011
 */
double  CdipIntegrand_Sb( double Lat, _qpInfo *qpInfo ) {

    double  CosLat, CosLat2, CosLat3, CosLat6, v;
    double  BoverBeq, BoverBm, SinAlpha2;
    double  SinAlpha02, Ek_in, L, aStarEq;
    double  f;
    Lgm_SummersInfo *si;

    /*
     * Pull parameters for integrand from Lgm_SummersInfo structure.
     */
    si  = (Lgm_SummersInfo *)qpInfo;
    SinAlpha02 = si->SinAlpha02;    // pre-computed sin( Alpha0 )

    CosLat  = cos( Lat );     CosLat2 = CosLat*CosLat;
    CosLat3 = CosLat2*CosLat; CosLat6 = CosLat3*CosLat3;
    v = sqrt(4.0 - 3.0*CosLat2);

    BoverBeq  = v/CosLat6;          // B/Beq
    BoverBm   = BoverBeq*SinAlpha02; // B/Bm = B/Beq * sin^2(Alpha0)
    SinAlpha2 = BoverBm;            // sin^2(Alpha) = B/Bm
    if (SinAlpha2 > 1.0) { SinAlpha2 = BoverBm = 1.0; }

    f = CosLat*v/sqrt(1.0-BoverBm);

    return( f );

}

double  SummersIntegrand_Gaa( double Lat, _qpInfo *qpInfo ) {

    double  CosLat, CosLat2, CosLat3, CosLat6, v, Omega_e, Omega_Sig, xm, dx;
    double  BoverBeq, BoverBm, SinAlpha2, Omega_eEq, Omega_SigEq, wm, dw;
    double  SinAlpha02, CosAlpha2, TanAlpha2, TanAlpha02, Ek_in, L, aStarEq;
    double  Daa, f, aStar, E, dB, B, dBoverB2, s, Lambda, Rho, Sig;
    Lgm_SummersInfo *si;

    /*
     * Pull parameters for integrand from Lgm_SummersInfo structure.
     */
    si  = (Lgm_SummersInfo *)qpInfo;
    SinAlpha02  = si->SinAlpha02;    // pre-computed sin^2( Alpha0 )
    TanAlpha02  = si->TanAlpha02;    // pre-computed tan^2( Alpha0 )
    E           = si->E;             // Dimensionless energy Ek/E0 (kinetic energy over rest energy).
    L           = si->L;
    aStarEq     = si->aStarEq;
    dB          = si->dB;
    Omega_eEq   = si->Omega_eEq;     // Equatorial electron gyro-frequency.
    Omega_SigEq = si->Omega_SigEq;   // Equatorial gyro-frequency for given species.
    wm          = si->wm;            // center frequency (at equator)
    dw          = si->dw;            // frequency bandwidth (at equator)
    Lambda      = si->Lambda;        // 
    s           = si->s;             // 
    Rho         = si->Rho;           // Rho = sqrt(PI)/2.0( erf( (wm-w1)/dw ) + erf( (w2-wm)/dw ) )   Eq (3) in Summers2007
    Sig         = si->Sig;           // Sig -- seetext above eqn (30).

    CosLat  = cos( Lat );     CosLat2 = CosLat*CosLat;
    CosLat3 = CosLat2*CosLat; CosLat6 = CosLat3*CosLat3;
    v = sqrt(4.0 - 3.0*CosLat2);

    BoverBeq  = v/CosLat6;              // B/Beq
    BoverBm   = BoverBeq*SinAlpha02;    // B/Bm = B/Beq * sin^2(Alpha0)
    SinAlpha2 = BoverBm;                // sin^2(Alpha) = B/Bm
    if (SinAlpha2 > 1.0) { SinAlpha2 = BoverBm = 1.0; }
    CosAlpha2 = 1.0-SinAlpha2;
    if (CosAlpha2 < 0.0) { CosAlpha2 = 0.0; }
    TanAlpha2 = SinAlpha2/CosAlpha2;

    /*
     * Compute local parameters from the equatorial ones that are specified.
     */
    B         = BoverBeq*M_CDIP/(L*L*L);     // local value of B (warning, hard-coded M value).
    dBoverB2  = dB*dB/(B*B);                 // Local value of R = dB^2/B^2.
    aStar     = aStarEq*BoverBeq*BoverBeq;   // Local aStar value.
    Omega_e   = Omega_eEq*BoverBeq;          // Local electron gyro-frequency.
    Omega_Sig = Omega_SigEq*BoverBeq;        // Local gyro-frequency for given species.
    xm        = wm/fabs(Omega_e);            // Local value of xm.
    dx        = dw/fabs(Omega_e);            // Local value of dx.


    /*
     * Compute the local Daa
     */
    Daa = Lgm_SummersDaaLocal( SinAlpha2, E, dBoverB2, BoverBeq, Omega_e, Omega_Sig, Rho, Sig, xm, dx, Lambda, s, aStar );

    /*
     * Finally the integrand. The TanAlpha02/TanAlpha2 term (which is dAlphaeq/dAlpha) converts D_Alpha,Alpha to
     * D_Alpha_eq,Alpha_eq.
     */
    f = (Daa*TanAlpha02/TanAlpha2)   * CosLat*v/sqrt(1-BoverBm);

    return( f );

}

double  SummersIntegrand_Gap( double Lat, _qpInfo *qpInfo ) {

    double  CosLat, CosLat2, CosLat3, CosLat6, v, Omega_e, Omega_Sig, xm, dx;
    double  BoverBeq, BoverBm, SinAlpha2, Omega_eEq, Omega_SigEq, wm, dw;
    double  SinAlpha02, CosAlpha2, TanAlpha2, TanAlpha02, Ek_in, L, aStarEq;
    double  Dap, f, aStar, E, dB, B, dBoverB2, Lambda, s, Rho, Sig;
    Lgm_SummersInfo *si;

    /*
     * Pull parameters for integrand from Lgm_SummersInfo structure.
     */
    si  = (Lgm_SummersInfo *)qpInfo;
    SinAlpha02  = si->SinAlpha02;    // pre-computed sin^2( Alpha0 )
    TanAlpha02  = si->TanAlpha02;    // pre-computed tan^2( Alpha0 )
    E           = si->E;             // Dimensionless energy Ek/E0 (kinetic energy over rest energy).
    L           = si->L;
    aStarEq     = si->aStarEq;
    dB          = si->dB;
    Omega_eEq   = si->Omega_eEq;     // Equatorial electron gyro-frequency.
    Omega_SigEq = si->Omega_SigEq;   // Equatorial gyro-frequency for given species.
    wm          = si->wm;            // center frequency (at equator)
    dw          = si->dw;            // frequency bandwidth (at equator)
    Lambda      = si->Lambda;
    s           = si->s;
    Rho         = si->Rho;           // Rho = sqrt(PI)/2.0( erf( (wm-w1)/dw ) + erf( (w2-wm)/dw ) )   Eq (3) in Summers2007
    Sig         = si->Sig;           // Sig -- seetext above eqn (30).

    CosLat  = cos( Lat );     CosLat2 = CosLat*CosLat;
    CosLat3 = CosLat2*CosLat; CosLat6 = CosLat3*CosLat3;
    v = sqrt(4.0 - 3.0*CosLat2);

    BoverBeq  = v/CosLat6;              // B/Beq
    BoverBm   = BoverBeq*SinAlpha02;    // B/Bm = B/Beq * sin^2(Alpha0)
    SinAlpha2 = BoverBm;                // sin^2(Alpha) = B/Bm
    if (SinAlpha2 > 1.0) { SinAlpha2 = BoverBm = 1.0; }
    CosAlpha2 = 1.0-SinAlpha2;
    if (CosAlpha2 < 0.0) { CosAlpha2 = 0.0; }
    TanAlpha2 = SinAlpha2/CosAlpha2;

    /*
     * Compute local parameters from the equatorial ones that are specified.
     */
    B         = BoverBeq*M_CDIP/(L*L*L);     // local value of B (warning, hard-coded M value).
    dBoverB2  = dB*dB/(B*B);                 // Local value of R = dB^2/B^2.
    aStar     = aStarEq*BoverBeq*BoverBeq;   // Local aStar value.
    Omega_e   = Omega_eEq*BoverBeq;          // Local electron gyro-frequency.
    Omega_Sig = Omega_SigEq*BoverBeq;        // Local gyro-frequency for given species.
    xm        = wm/fabs(Omega_e);            // Local value of xm.
    dx        = dw/fabs(Omega_e);            // Local value of dx.


    /*
     * Compute the local Daa
     */
//FIX
    Dap = Lgm_SummersDaaLocal( SinAlpha2, E, dBoverB2, BoverBeq, Omega_e, Omega_Sig, Rho, Sig, xm, dx, Lambda, s, aStar );

    /*
     * Finally the integrand. The TanAlpha02/TanAlpha2 term (which is dAlphaeq/dAlpha) converts D_Alpha,Alpha to
     * D_Alpha_eq,Alpha_eq.
     */
    f = (Dap*TanAlpha02/TanAlpha2)   * CosLat*v/sqrt(1-BoverBm);

    return( f );

}

double  SummersIntegrand_Gpp( double Lat, _qpInfo *qpInfo ) {

    double  CosLat, CosLat2, CosLat3, CosLat6, v, Omega_e, Omega_Sig, xm, dx;
    double  BoverBeq, BoverBm, SinAlpha2, Omega_eEq, Omega_SigEq, wm, dw;
    double  SinAlpha02, CosAlpha2, TanAlpha2, TanAlpha02, Ek_in, L, aStarEq;
    double  Dpp, f, aStar, E, dB, B, dBoverB2, Lambda, s, Rho, Sig;
    Lgm_SummersInfo *si;

    /*
     * Pull parameters for integrand from Lgm_SummersInfo structure.
     */
    si  = (Lgm_SummersInfo *)qpInfo;
    SinAlpha02  = si->SinAlpha02;    // pre-computed sin^2( Alpha0 )
    TanAlpha02  = si->TanAlpha02;    // pre-computed tan^2( Alpha0 )
    E           = si->E;             // Dimensionless energy Ek/E0 (kinetic energy over rest energy).
    L           = si->L;
    aStarEq     = si->aStarEq;
    dB          = si->dB;
    Omega_eEq   = si->Omega_eEq;     // Equatorial electron gyro-frequency.
    Omega_SigEq = si->Omega_SigEq;   // Equatorial gyro-frequency for given species.
    wm          = si->wm;            // center frequency (at equator)
    dw          = si->dw;            // frequency bandwidth (at equator)
    Lambda      = si->Lambda;
    s           = si->s;
    Rho         = si->Rho;           // Rho = sqrt(PI)/2.0( erf( (wm-w1)/dw ) + erf( (w2-wm)/dw ) )   Eq (3) in Summers2007
    Sig         = si->Sig;           // Sig -- seetext above eqn (30).

    CosLat  = cos( Lat );     CosLat2 = CosLat*CosLat;
    CosLat3 = CosLat2*CosLat; CosLat6 = CosLat3*CosLat3;
    v = sqrt(4.0 - 3.0*CosLat2);

    BoverBeq  = v/CosLat6;              // B/Beq
    BoverBm   = BoverBeq*SinAlpha02;    // B/Bm = B/Beq * sin^2(Alpha0)
    SinAlpha2 = BoverBm;                // sin^2(Alpha) = B/Bm
    if (SinAlpha2 > 1.0) { SinAlpha2 = BoverBm = 1.0; }
    CosAlpha2 = 1.0-SinAlpha2;
    if (CosAlpha2 < 0.0) { CosAlpha2 = 0.0; }
    TanAlpha2 = SinAlpha2/CosAlpha2;

    /*
     * Compute local parameters from the equatorial ones that are specified.
     */
    B         = BoverBeq*M_CDIP/(L*L*L);     // local value of B (warning, hard-coded M value).
    dBoverB2  = dB*dB/(B*B);                 // Local value of R = dB^2/B^2.
    aStar     = aStarEq*BoverBeq*BoverBeq;   // Local aStar value.
    Omega_e   = Omega_eEq*BoverBeq;          // Local electron gyro-frequency.
    Omega_Sig = Omega_SigEq*BoverBeq;        // Local gyro-frequency for given species.
    xm        = wm/fabs(Omega_e);            // Local value of xm.
    dx        = dw/fabs(Omega_e);            // Local value of dx.


    /*
     * Compute the local Daa
     */
//FIX
    Dpp = Lgm_SummersDaaLocal( SinAlpha2, E, dBoverB2, BoverBeq, Omega_e, Omega_Sig, Rho, Sig, xm, dx, Lambda, s, aStar );

    /*
     * Finally the integrand. The TanAlpha02/TanAlpha2 term (which is dAlphaeq/dAlpha) converts D_Alpha,Alpha to
     * D_Alpha_eq,Alpha_eq.
     */
    f = (Dpp*TanAlpha02/TanAlpha2)   * CosLat*v/sqrt(1-BoverBm);

    return( f );

}



// FLESH this out.
/**
 *  \brief
 *      Computes the local Summer's [2005] diffusion coefficients.
 *
 *  \details
 *
 *
 *
 *
 *      \param[in]      Alpha0      equatoria PA in radians.
 *      \param[in]      E           Dimensionless energy Ek/E0 (kinetic energy over rest mass).
 *      \param[in]      dBoverB2    Ratio of wave amplitude, dB to local background field, B.
 *      \param[in]      BoverBeq    Ratio of local B to Beq.
 *      \param[in]      aStar  Local value of the aStar parameter.
 *
 *      \return         Local value of Daa
 *
 *      \author         M. Henderson
 *      \date           2010-2011
 *
 */
double Lgm_SummersDaaLocal( double SinAlpha2, double E, double dBoverB2, double BoverBeq, double Omega_e, double Omega_Sig, double Rho, double Sig, double xm, double dx, double Lambda, double s, double aStar ) {

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
    //printf("a1, a2, a3, a4 = %g %g %g %g\n", a1, a2, a3, a4 );

    nReal = Lgm_QuarticRoots( a1, a2, a3, a4, &z1, &z2, &z3, &z4 );

    //printf("nReal = %d\n", nReal );
    //printf("z1    = %g + %g i\n", creal(z1), cimag(z1) );
    //printf("z2    = %g + %g i\n", creal(z2), cimag(z2) );
    //printf("z3    = %g + %g i\n", creal(z3), cimag(z3) );
    //printf("z4    = %g + %g i\n", creal(z4), cimag(z4) );



    R  = dBoverB2;   // The ratio (dB/B)^2


int sum_res = 0;



    /*
     * Gather applicable roots together into and array
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
            g = x4 + c1*x3 + c2*x2 + c3*x + c4;
            xms  = x - s;
            xpse = x + s*LGM_EPS;
            F = y*xms*xms*xpse*xpse/(x*g);

            u = 1.0 - x*Mu/(y*Beta);
            arg = (x-xm)/dx;
            Daa += u*u *fabs(F) / ( dx * fabs( BetaMu - F ) ) * exp( -arg*arg );

        }
        Daa *= fac;

//printf("%g  %g\n", 180.0/M_PI*acos(Mu), Daa);


    }







    return(Daa);

}
