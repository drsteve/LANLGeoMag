#include "Lgm/Lgm_SummersDiffCoeff.h"


/**
 *  \brief
 *      Computes the bounce-averaged Summer's [2005] diffusion coefficients.
 *
 *  \details
 *
 *      \param[in]      Alpha0      equatoria PA in radians.
 *      \param[in]      Ek_in       kinetic energy in MeV.
 *      \param[in]      L           L-shell parameter (dimensionless).
 *      \param[in]      astar       Cold plasma parameter.
 *      \param[out]     Daa_ba      Bounce-averaged value of Daa.
 *      \param[out]     Dap_ba      Bounce-averaged value of Dap.
 *      \param[out]     Dpp_ba      Bounce-averaged value of Dpp.
 *
 *      \return         void
 *
 *      \author         J. Koller, S. Zaharia, M. Henderson
 *      \date           2010-2011
 *
 */
void getDs_ba( double Alpha0,  double Ek_in,  double L,  double astar,  double *Daa_ba,  double *Dap_ba,  double *Dpp_ba) {

    double          T0, T1, T, a, b;
    double          epsabs, epsrel, result, abserr, resabs, resasc, work[2001], points[3];
    int             npts=4, key=6, limit=500, lenw=4*limit, iwork[501], last, ier, neval;
    Lgm_SummersInfo *s=(Lgm_SummersInfo *)calloc( 1, sizeof(Lgm_SummersInfo));


    /*
     *  Quadpack integration tolerances.
1e-10 seems exccessive.
     */
    epsabs = 1e-10;
    epsrel = 1e-3;






    /*
     *  Pack integrand parameters into structure
     */
    s->Alpha0     = Alpha0*M_PI/180.0;
    s->SinAlpha0  = sin(s->Alpha0);
    s->SinAlpha02 = s->SinAlpha0*s->SinAlpha0;
    s->CosAlpha0  = cos(s->Alpha0);
    s->CosAlpha02 = s->CosAlpha0*s->CosAlpha0;
    s->TanAlpha0  = s->SinAlpha0/s->CosAlpha0;
    s->TanAlpha02 = s->TanAlpha0*s->TanAlpha0;
    s->Ek_in      = Ek_in;
    s->L          = L;
    s->astar      = astar;


    /*
     *  Set integration limits
     */
    a = 0.0;                                        // radians
    b = acos( Lgm_CdipMirrorLat( s->SinAlpha0 ) );  // radians
    //lambda_m = MIN( lambda_m, 25.*RadPerDeg );  // equatatorial confinement of waves.

    /*
     *  Perform integrations.
     */
    points[1] = a; points[2] = b; npts=2;
    dqagp( CdipIntegrand_Sb, (_qpInfo *)s, a, b, npts, points, epsabs, epsrel, &T, &abserr, &neval, &ier, limit, lenw, &last, iwork, work );
printf("neval = %d\n", neval);
    dqagp( SummersIntegrand_Gaa, (_qpInfo *)s, a, b, npts, points, epsabs, epsrel, Daa_ba, &abserr, &neval, &ier, limit, lenw, &last, iwork, work );
printf("neval = %d\n", neval);
    dqagp( SummersIntegrand_Gap, (_qpInfo *)s, a, b, npts, points, epsabs, epsrel, Dap_ba, &abserr, &neval, &ier, limit, lenw, &last, iwork, work );
printf("neval = %d\n", neval);
    dqagp( SummersIntegrand_Gpp, (_qpInfo *)s, a, b, npts, points, epsabs, epsrel, Dpp_ba, &abserr, &neval, &ier, limit, lenw, &last, iwork, work );
printf("neval = %d\n", neval);


    /*
     *  Calculate T( SinAlpha0 ), related to bounce period.
     *  Schultz and lanzeroti Eqns 1.28a-d)
     */
    //T0 = 1.38017299815047317375;
    //T1 = 0.74048048969306104115;
    //T  = T0 - 0.5*(T0-T1)*(s->SinAlpha0 + sqrt(s->SinAlpha0));


    *Daa_ba /= T;
    *Dap_ba /= T;
    *Dpp_ba /= T;


    return;
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
 *      \author         J. Koller, S. Zaharia, M. Henderson
 *      \date           2010-2011
 */
double  CdipIntegrand_Sb( double Lat, _qpInfo *qpInfo ) {

    double  CosLat, CosLat2, CosLat3, CosLat6, v;
    double  BoverB0, BoverBm, SinAlpha2;
    double  SinAlpha02, Ek_in, L, astar;
    double  f;
    Lgm_SummersInfo *s;

    /*
     * Pull parameters for integrand from Lgm_SummersInfo structure.
     */
    s  = (Lgm_SummersInfo *)qpInfo;
    SinAlpha02 = s->SinAlpha02;    // pre-computed sin( Alpha0 )

    CosLat  = cos( Lat );     CosLat2 = CosLat*CosLat;
    CosLat3 = CosLat2*CosLat; CosLat6 = CosLat3*CosLat3;
    v = sqrt(4.0 - 3.0*CosLat2);

    BoverB0   = v/CosLat6;          // B/B0
    BoverBm   = BoverB0*SinAlpha02; // B/Bm = B/B0 * sin^2(Alpha0)
    SinAlpha2 = BoverBm;            // sin^2(Alpha) = B/Bm

    f = CosLat*v/sqrt(1-BoverBm);

    return( f );

}
double  SummersIntegrand_Gaa( double Lat, _qpInfo *qpInfo ) {

    double  CosLat, CosLat2, CosLat3, CosLat6, v;
    double  BoverB0, BoverBm, SinAlpha2;
    double  SinAlpha02, CosAlpha2, TanAlpha2, TanAlpha02, Ek_in, L, astar;
    double  Daa, Dap, Dpp, f;
    Lgm_SummersInfo *s;

    /*
     * Pull parameters for integrand from Lgm_SummersInfo structure.
     */
    s  = (Lgm_SummersInfo *)qpInfo;
    SinAlpha02 = s->SinAlpha02;    // pre-computed sin^2( Alpha0 )
    TanAlpha02 = s->TanAlpha02;    // pre-computed tan^2( Alpha0 )
    Ek_in      = s->Ek_in;
    L          = s->L;
    astar      = s->astar;

    CosLat  = cos( Lat );     CosLat2 = CosLat*CosLat;
    CosLat3 = CosLat2*CosLat; CosLat6 = CosLat3*CosLat3;
    v = sqrt(4.0 - 3.0*CosLat2);

    BoverB0   = v/CosLat6;          // B/B0
    BoverBm   = BoverB0*SinAlpha02; // B/Bm = B/B0 * sin^2(Alpha0)
    SinAlpha2 = BoverBm;            // sin^2(Alpha) = B/Bm
    CosAlpha2 = 1.0-SinAlpha2;
    TanAlpha2 = SinAlpha2/CosAlpha2;

    getDs( SinAlpha2, Ek_in, L, Lat, astar, &Daa, &Dap, &Dpp );

    f = (Daa*TanAlpha02/TanAlpha2)   * CosLat*v/sqrt(1-BoverBm);

    return( f );

}

double  SummersIntegrand_Gap( double Lat, _qpInfo *qpInfo ) {

    double  CosLat, CosLat2, CosLat3, CosLat6, v;
    double  BoverB0, BoverBm, SinAlpha2, CosAlpha2, TanAlpha2;
    double  SinAlpha02, TanAlpha02, Ek_in, L, astar;
    double  Daa, Dap, Dpp, f;
    Lgm_SummersInfo *s;

    /*
     * Pull parameters for integrand from Lgm_SummersInfo structure.
     */
    s  = (Lgm_SummersInfo *)qpInfo;
    SinAlpha02 = s->SinAlpha02;    // pre-computed sin^2( Alpha0 )
    TanAlpha02 = s->SinAlpha02;    // pre-computed tan^2( Alpha0 )
    Ek_in      = s->Ek_in;
    L          = s->L;
    astar      = s->astar;

    CosLat  = cos( Lat );     CosLat2 = CosLat*CosLat;
    CosLat3 = CosLat2*CosLat; CosLat6 = CosLat3*CosLat3;
    v = sqrt(4.0 - 3.0*CosLat2);

    BoverB0   = v/CosLat6;          // B/B0
    BoverBm   = BoverB0*SinAlpha02; // B/Bm = B/B0 * sin^2(Alpha0)
    SinAlpha2 = BoverBm;            // sin^2(Alpha) = B/Bm
    CosAlpha2 = 1.0-SinAlpha2;
    TanAlpha2 = SinAlpha2/CosAlpha2;

    getDs( SinAlpha2, Ek_in, L, Lat, astar, &Daa, &Dap, &Dpp );

    f = (Dap*sqrt(TanAlpha02/TanAlpha2)) * CosLat*v/sqrt(1-BoverBm);

    return( f );

}

double  SummersIntegrand_Gpp( double Lat, _qpInfo *qpInfo ) {

    double  CosLat, CosLat2, CosLat3, CosLat6, v;
    double  BoverB0, BoverBm, SinAlpha2;
    double  SinAlpha02, Ek_in, L, astar;
    double  Daa, Dap, Dpp, f;
    Lgm_SummersInfo *s;

    /*
     * Pull parameters for integrand from Lgm_SummersInfo structure.
     */
    s  = (Lgm_SummersInfo *)qpInfo;
    SinAlpha02 = s->SinAlpha02;    // pre-computed sin( Alpha0 )
    Ek_in      = s->Ek_in;
    L          = s->L;
    astar      = s->astar;

    CosLat  = cos( Lat );     CosLat2 = CosLat*CosLat;
    CosLat3 = CosLat2*CosLat; CosLat6 = CosLat3*CosLat3;
    v = sqrt(4.0 - 3.0*CosLat2);

    BoverB0   = v/CosLat6;          // B/B0
    BoverBm   = BoverB0*SinAlpha02; // B/Bm = B/B0 * sin^2(Alpha0)
    SinAlpha2 = BoverBm;            // sin^2(Alpha) = B/Bm

    getDs( SinAlpha2, Ek_in, L, Lat, astar, &Daa, &Dap, &Dpp );

    f = Dpp * CosLat*v/sqrt(1-BoverBm);

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
 *      \param[in]      Ek_in       kinetic energy in MeV.
 *      \param[in]      L           L-shell parameter (dimensionless).
 *      \param[in]      astar       Cold plasma parameter.
 *      \param[out]     Daa         Local value of Daa.
 *      \param[out]     Dap         Local value of Dap.
 *      \param[out]     Dpp         Local value of Dpp.
 *
 *      \return         void
 *
 *      \author         J. Koller, S. Zaharia, M. Henderson
 *      \date           2010-2011
 *
 */
void getDs( double alpha, double Ek_in, double L, double Lat, double astar, double *Daa, double *Dap, double *Dpp ) {

    // test
    *Daa = 1.0;
    *Dap = 1.0;
    *Dpp = 1.0;

    return;

}
