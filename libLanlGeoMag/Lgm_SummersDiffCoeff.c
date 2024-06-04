#include "Lgm/Lgm_SummersDiffCoeff.h"

#include <gsl/gsl_poly.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_errno.h>

struct inequalityParameters {double D; double E1; double E2; double F1; double F2;};
struct normalizerForWavePowerSpectrumFunctionParams {double x; double aStar; int numberOfWaveNormalAngleDistributions; double *xmArray; double *dxArray; double *weightsOnWaveNormalAngleDistributions; double xmin; double xmax;};
struct localDiffusionCoefficientAtSpecificThetaGlauertAndHorneFunctionParams {int tensorFlag; int s; double KE; double aStar; double wce; double xmin; double xmax; int numberOfWaveNormalAngleDistributions; double *xmArray; double *dxArray; double *weightsOnWaveNormalAngleDistributions; double wm; double dw; double wlc; double wuc; double Bw; double alpha; int nCyclotronLow; int nCyclotronHigh; double *Nw; int nNw; int Dir;};
double normalizerForWavePowerSpectrumFunction(double tanTheta, void *p);
int computeResonantRootsUsingColdElectronPlasma(double KE, double theta, double alpha, int s, int nCyclotron, double aStar, double xlc, double xuc, double *resonantRoots, int *nRoots);
double besselFunctionNormalizer(int s, double p, double alpha, int nCyclotronNumber, double x, double y, double P, double R, double L, double S, double theta);
double electronColdPlasmaGroupVelocity(double aStar, double tanTheta2, double x, double y);
double localDiffusionCoefficientAtSpecificThetaGlauertAndHorneWeightedByTanTheta(double tanTheta, struct localDiffusionCoefficientAtSpecificThetaGlauertAndHorneFunctionParams *p);
double localDiffusionCoefficientAtSpecificThetaGlauertAndHorne(double tanTheta, struct localDiffusionCoefficientAtSpecificThetaGlauertAndHorneFunctionParams *p);
//old calling style double localDiffusionCoefficientAtSpecificThetaGlauertAndHorne(int tensorFlag, int s, double KE, double aStar, double wce, double xm, double xmin, double xmax, double dx, double wm, double dw, double wlc, double wuc, double Bw, double alpha, double tanTheta, int nCyclotronLow, int nCyclotronHigh, double *Nw, int nNw, int Directions);

int Lgm_SimpleRiemannSum( double (*f)( double, _qpInfo * ), _qpInfo *qpInfo, double xleft, double xright, double *result, int VerbosityLevel );
int computeInequalityParameters(double x, double aStar, double vpar, struct inequalityParameters *params);
int computeWaveNormalAngleIntervalsWhereResonantRootsMayExist(  double wn, double wlc, double wuc, struct inequalityParameters *inequalityParamsLC, 
                                                                struct inequalityParameters *inequalityParamsUC, 
                                                                int mode, int *nIntervals, double *thetaStart, double *thetaEnd);

int findSingleIntervalInTheta(struct inequalityParameters *params, double dw, int inequalitySign, int *nIntervals, double *thetaStart, double *thetaEnd);

int mergeOneIntervalWithOneInterval(   double *thisThetaStart, double *thisThetaEnd, double *thetaStart, double *thetaEnd );
int mergeTwoIntervalsWithOneInterval(  double *thisThetaStart, double *thisThetaEnd, double *thetaStart, double *thetaEnd );
int mergeOneIntervalWithTwoIntervals(  double *thisThetaStart, double *thisThetaEnd, double *thetaStart, double *thetaEnd );
int mergeTwoIntervalsWithTwoIntervals( double *thisThetaStart, double *thisThetaEnd, double *thetaStart, double *thetaEnd );

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
int Lgm_SummersDxxBounceAvg( int Version, double Alpha0,  double Ek,  double L,  void *BwFuncData, double (*BwFunc)(), double n1, double n2, double n3, double aStarEq,  int Directions, double w1, double w2, double wm, double dw, int WaveMode, int Species, double MaxWaveLat, double *Daa_ba,  double *Dap_ba,  double *Dpp_ba) {

    double           T, a, b, E0, Omega_eEq, Omega_SigEq, Beq, Rho;
    double           epsabs, epsrel, abserr, work[2002], points[10];
    double           ySing, a_new, b_new, B;
    int              npts2=4, limit=500, lenw=4*limit, iwork[502], last, ier, neval;
    int              VerbosityLevel = 0;
    Lgm_SummersInfo  si;
    si.Version = Version;
    si.n1 = n1;
    si.n2 = n2;
    si.n3 = n3;
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
        si.s = 1;
    } else if ( WaveMode == LGM_L_MODE_WAVE ) {
        si.s = -1;
    } else {
        printf("Lgm_SummersDxxBounceAvg(): Unknown wavemode = %d\n", WaveMode );
        return(0);
    }


    Beq = M_CDIP/(L*L*L);
    if ( Species == LGM_ELECTRONS ) {
        si.Lambda  = -1.0;
        E0          = LGM_Ee0;      // set rest energy to electron rest energy (MeV)
        Omega_SigEq = Lgm_GyroFreq( -LGM_e, Beq, LGM_ELECTRON_MASS );
    } else if ( Species == LGM_PROTONS ) {
        si.Lambda = LGM_EPS;
        E0         = LGM_Ep0;       // set rest energy to proton rest energy (MeV)
        Omega_SigEq = Lgm_GyroFreq( LGM_e, Beq, LGM_PROTON_MASS );
    } else {
        printf("Lgm_SummersDxxBounceAvg(): Unknown species = %d\n", Species );
        return(0);
    }
    Omega_eEq = Lgm_GyroFreq( -LGM_e, Beq, LGM_ELECTRON_MASS );
    si.Omega_eEq   = Omega_eEq;    // Equatorial electron gyro-frequency.
    si.Omega_SigEq = Omega_SigEq;  // Equatorial gyro-frequency for given species.

    si.E           = Ek/E0;        // Dimensionless energy Ek/E0 (kinetic energy over rest energy).
    si.Alpha0      = Alpha0*M_PI/180.0;
    si.SinAlpha0   = sin(si.Alpha0);
    si.SinAlpha02  = si.SinAlpha0*si.SinAlpha0;
    si.CosAlpha0   = cos(si.Alpha0);
    si.CosAlpha02  = si.CosAlpha0*si.CosAlpha0;
    si.TanAlpha0   = si.SinAlpha0/si.CosAlpha0;
    si.TanAlpha02  = si.TanAlpha0*si.TanAlpha0;
    si.L           = L;
    si.aStarEq     = aStarEq;
//    si.dB          = dB;
    si.BwFuncData  = BwFuncData;
    si.BwFunc      = BwFunc;
    si.w1          = w1;           // Lower freq cuttof.
    si.w2          = w2;           // Upper freq cuttof.
    si.wm          = wm;           // Frequency of max wave power.
    si.dw          = dw;           // Frequency bandwidth (at equator).
    si.MaxWaveLat  = MaxWaveLat*RadPerDeg;   // Assume there are no waves at +/-MaxWaveLat.
    si.Rho         = Rho;
    si.Directions  = Directions;


    /*
     *  Set integration limits
     */
    a = 0.0;                                            // radians
    double tempb;
    tempb = acos(Lgm_CdipMirrorLat( si.SinAlpha0 ));     // radians
    if (acos(sqrt(1.0/L)) < tempb) {B=acos(sqrt(1.0/L));} else {B=tempb;} // Modified by Greg Cunningham on 10/14/2014 so that integral doesn't extend to inside earth for low L
    b = ( B < si.MaxWaveLat ) ? B : si.MaxWaveLat;
    if ( fabs(a-b) < 1e-9 ) {

        *Daa_ba = 0.0;
        *Dap_ba = 0.0;
        *Dpp_ba = 0.0;
        return( 0 );

    } else {


        /*
         *  Perform integrations. The points[] array contains a list of points
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
         *      T  = T0 - 0.5*(T0-T1)*(si.SinAlpha0 + sqrt(si.SinAlpha0));
         */
        npts2 = 2;
        //dqagp( CdipIntegrand_Sb, (_qpInfo *)&si, a, b, npts2, points, epsabs, epsrel, &T, &abserr, &neval, &ier, limit, lenw, &last, iwork, work );
        dqags( CdipIntegrand_Sb, (_qpInfo *)&si, a, B, epsabs, epsrel, &T, &abserr, &neval, &ier, limit, lenw, &last, iwork, work, VerbosityLevel );


        Lgm_SummersFindCutoffs2( SummersIntegrand_Gaa, (_qpInfo *)&si, FALSE, a, b, &a_new, &b_new );
        npts2 = 2 + Lgm_SummersFindSingularities( SummersIntegrand_Gaa, (_qpInfo *)&si, FALSE, a_new, b_new, &points[1], &ySing );
        if ( b_new > a_new ) {
            if ( npts2 > 2 ) {
//printf("There is an interior singularity in Daa: Alpha0=%g, Ek=%g\n", Alpha0, Ek);
                /* Approach #1 from Mike: use dqagp to break up the integral into regions that have singularities.  */
                dqagp( SummersIntegrand_Gaa, (_qpInfo *)&si, a_new, b_new, npts2, points, epsabs, epsrel, Daa_ba, &abserr, &neval, &ier, limit, lenw, &last, iwork, work, VerbosityLevel );
		/* Approach #2: break up the integral into regions without interior singularities by hand and use dqags to ingegrate 
                dqags( SummersIntegrand_Gaa, (_qpInfo *)&si, a_new, points[1]-0.001, epsabs, epsrel, Daa_ba, &abserr, &neval, &ier, limit, lenw, &last, iwork, work, VerbosityLevel );
                dqags( SummersIntegrand_Gaa, (_qpInfo *)&si, points[1]+0.001, b_new, epsabs, epsrel, Daa_ba, &abserr, &neval, &ier, limit, lenw, &last, iwork, work, VerbosityLevel );
 		*/
		/* Approach #3: use GSL's implementation of quadpack to do the integral in parts
		gsl_integration_workspace * w =gsl_integration_workspace_alloc( 1000 );
		gsl_function F;
		F.function = SummersIntegrand_Gaa;
		F.params = &si;
		gsl_integration_qags(&F, a_new, points[1]-0.001, epsabs, epsrel, (size_t) 1000, w, Daa_ba, &abserr);
		gsl_integration_qags(&F, points[1]+0.001, b_new, epsabs, epsrel, (size_t) 1000, w, Daa_ba, &abserr);
		gsl_integration_workspace_free(w);
 		*/
		/* Approach #4: use GSL's implementation of quadpack to do the integral all at once
		gsl_integration_workspace * w =gsl_integration_workspace_alloc( 1000 );
		gsl_function F;
		F.function = SummersIntegrand_Gaa;
		F.params = &si;
		points[0] = a_new; points[2] = b_new;
		gsl_integration_qagp(&F, points, (size_t) 3, epsabs, epsrel, (size_t) 1000, w, Daa_ba, &abserr);
		gsl_integration_workspace_free(w);
 		* */
            } else {
                dqags( SummersIntegrand_Gaa, (_qpInfo *)&si, a_new, b_new, epsabs, epsrel, Daa_ba, &abserr, &neval, &ier, limit, lenw, &last, iwork, work, VerbosityLevel );
		/* Approach #2: use GSL's implementation of the simple quadpack routine dqag to do the integral 
		gsl_integration_workspace * w =gsl_integration_workspace_alloc( 1000 );
		gsl_function F;
		F.function = SummersIntegrand_Gaa;
		F.params = &si;
		gsl_integration_qag(&F, a_new, b_new, epsabs, epsrel, (size_t) 1000, GSL_INTEG_GAUSS15, w, Daa_ba, &abserr);
		gsl_integration_workspace_free(w);
 		*/
            }
        } else {
            *Daa_ba = 0.0;
        }

        Lgm_SummersFindCutoffs2( SummersIntegrand_Gap, (_qpInfo *)&si, FALSE, a, b, &a_new, &b_new );
        npts2 = 2 + Lgm_SummersFindSingularities( SummersIntegrand_Gap, (_qpInfo *)&si, FALSE, a_new, b_new, &points[1], &ySing );
        if ( b_new > a_new ) {
            if ( npts2 > 2 ) {
		/* Approach 1 */
                dqagp( SummersIntegrand_Gap, (_qpInfo *)&si, a_new, b_new, npts2, points, epsabs, epsrel, Dap_ba, &abserr, &neval, &ier, limit, lenw, &last, iwork, work, VerbosityLevel );
		/* Approach 2 
                dqags( SummersIntegrand_Gap, (_qpInfo *)&si, a_new, points[1]-0.001, epsabs, epsrel, Dap_ba, &abserr, &neval, &ier, limit, lenw, &last, iwork, work, VerbosityLevel );
                dqags( SummersIntegrand_Gap, (_qpInfo *)&si, points[1]+0.001, b_new, epsabs, epsrel, Dap_ba, &abserr, &neval, &ier, limit, lenw, &last, iwork, work, VerbosityLevel );
 		*/
            } else {
                dqags( SummersIntegrand_Gap, (_qpInfo *)&si, a_new, b_new, epsabs, epsrel, Dap_ba, &abserr, &neval, &ier, limit, lenw, &last, iwork, work, VerbosityLevel );
		/* Approach #2: use GSL's implementation of the simple quadpack routine dqag to do the integral 
		gsl_integration_workspace * w =gsl_integration_workspace_alloc( 1000 );
		gsl_function F;
		F.function = SummersIntegrand_Gap;
		F.params = &si;
		gsl_integration_qag(&F, a_new, b_new, epsabs, epsrel, (size_t) 1000, GSL_INTEG_GAUSS15, w, Dap_ba, &abserr);
		gsl_integration_workspace_free(w);
 		*/
            }
        } else {
            *Dap_ba = 0.0;
        }

        Lgm_SummersFindCutoffs2( SummersIntegrand_Gpp, (_qpInfo *)&si, FALSE, a, b, &a_new, &b_new );
        npts2 = 2 + Lgm_SummersFindSingularities( SummersIntegrand_Gpp, (_qpInfo *)&si, FALSE, a_new, b_new, &points[1], &ySing );
        if ( b_new > a_new ) {
            if ( npts2 > 2 ) {
		/* Approach 1 */
                dqagp( SummersIntegrand_Gpp, (_qpInfo *)&si, a_new, b_new, npts2, points, epsabs, epsrel, Dpp_ba, &abserr, &neval, &ier, limit, lenw, &last, iwork, work, VerbosityLevel );
		/* Approach 2 
                dqags( SummersIntegrand_Gpp, (_qpInfo *)&si, a_new, points[1]-0.001, epsabs, epsrel, Dpp_ba, &abserr, &neval, &ier, limit, lenw, &last, iwork, work, VerbosityLevel );
                dqags( SummersIntegrand_Gpp, (_qpInfo *)&si, points[1]+0.001, b_new, epsabs, epsrel, Dpp_ba, &abserr, &neval, &ier, limit, lenw, &last, iwork, work, VerbosityLevel );
 		*/
            } else {
                dqags( SummersIntegrand_Gpp, (_qpInfo *)&si, a_new, b_new, epsabs, epsrel, Dpp_ba, &abserr, &neval, &ier, limit, lenw, &last, iwork, work, VerbosityLevel );
		/* Approach #2: use GSL's implementation of the simple quadpack routine dqag to do the integral 
		gsl_integration_workspace * w =gsl_integration_workspace_alloc( 1000 );
		gsl_function F;
		F.function = SummersIntegrand_Gpp;
		F.params = &si;
		gsl_integration_qag(&F, a_new, b_new, epsabs, epsrel, (size_t) 1000, GSL_INTEG_GAUSS15, w, Dpp_ba, &abserr);
		gsl_integration_workspace_free(w);
 		*/
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

/*
 *  Following routine is identical to Lgm_SummersDxxBounceAvg except that there are additional parameters that describe the wave normal angle distribution, which in this routine, 
 *  consists of a sum of gaussians in the tangent of the wave normal angle.  
 *
 *   	x1					tangent of the minimum wave normal angle (the angle between the magnetic field and the k-vector of the wave)
 *   	x2					tangent of the maximum wave normal angle 
 *   	numberOfWaveNormalAngleDistributions	the number of gaussian distributions in the tangent of the wave normal angle that are weighted and summed together	
 *   	xmArray					an array with at least numberOfWaveNormalAngleDistributions elements, each element contains the	tangent of the 
 *   							mean wave normal angle for that component
 *   	dxArray					an array with at least numberOfWaveNormalAngleDistributions elements, each element contains the	tangent of the 
 *   	 						spread in wave normal angles for that component
 */
 
int Lgm_GlauertAndHorneDxxBounceAvg( int Version, double Alpha0,  double Ek,  double L,  void *BwFuncData, double (*BwFunc)(), double n1, double n2, double n3, double aStarEq,  int Directions, double w1, double w2, double wm, double dw, double x1, double x2, int numberOfWaveNormalAngleDistributions, double *xmArray, double *dxArray, double *weightsOnWaveNormalAngleDistributions, int WaveMode, int Species, double MaxWaveLat, int nNw, int nPlasmaParameters, double aStarMin, double aStarMax, double *Nw, double *Daa_ba,  double *Dap_ba,  double *Dpp_ba) {


    double           T, a, b, E0, Omega_eEq, Omega_SigEq, Beq, Rho;
    double           epsabs, epsrel, abserr, work[2002], points[10];
    double           ySing, a_new, b_new, B;
    int              npts2=4, limit=500, lenw=4*limit, iwork[502], last, ier, neval;
    int              VerbosityLevel = 1;

    Lgm_SummersInfo  si;
printf("I'm really in Lgm_GlauertAndHorneDxxBounceAvg\n");

/*  Added by Greg Cunningham to enable precomputation of the normalizing function N(w) that integrates over tan(theta) interval [xmin,xmax] and depends only on aStar */
    si.aStarMin = aStarMin;
    si.aStarMax = aStarMax;
    si.nNw = nNw;
    si.nPlasmaParameters= nPlasmaParameters;
    si.Nw = (double *) Nw;

    si.Version = Version;
    si.n1 = n1;
    si.n2 = n2;
    si.n3 = n3;
    if ( Version == LGM_SUMMERS_2007 ){
        if ( fabs(1.0-(n1+n2+n3)) > 1e-10 ) {
            printf("Lgm_GlauertAndHorneDxxBounceAvg: n1+n2+n3 is not 1. Got n1+n2+n3 = %g\n", n1+n2+n3 );
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
        si.s = 1;
    } else if ( WaveMode == LGM_L_MODE_WAVE ) {
        si.s = -1;
    } else {
        printf("Lgm_SummersDxxBounceAvg(): Unknown wavemode = %d\n", WaveMode );
        return(0);
    }


    Beq = M_CDIP/(L*L*L);
    if ( Species == LGM_ELECTRONS ) {
        si.Lambda  = -1.0;
        E0          = LGM_Ee0;      // set rest energy to electron rest energy (MeV)
        Omega_SigEq = Lgm_GyroFreq( -LGM_e, Beq, LGM_ELECTRON_MASS );
    } else if ( Species == LGM_PROTONS ) {
        si.Lambda = LGM_EPS;
        E0         = LGM_Ep0;       // set rest energy to proton rest energy (MeV)
        Omega_SigEq = Lgm_GyroFreq( LGM_e, Beq, LGM_PROTON_MASS );
    } else {
        printf("Lgm_SummersDxxBounceAvg(): Unknown species = %d\n", Species );
        return(0);
    }
    Omega_eEq = Lgm_GyroFreq( -LGM_e, Beq, LGM_ELECTRON_MASS );
    si.Omega_eEq   = Omega_eEq;    // Equatorial electron gyro-frequency.
    si.Omega_SigEq = Omega_SigEq;  // Equatorial gyro-frequency for given species.

    si.E           = Ek/E0;        // Dimensionless energy Ek/E0 (kinetic energy over rest energy).
    si.Alpha0      = Alpha0*M_PI/180.0;
    si.SinAlpha0   = sin(si.Alpha0);
    si.SinAlpha02  = si.SinAlpha0*si.SinAlpha0;
    si.CosAlpha0   = cos(si.Alpha0);
    si.CosAlpha02  = si.CosAlpha0*si.CosAlpha0;
    si.TanAlpha0   = si.SinAlpha0/si.CosAlpha0;
    si.TanAlpha02  = si.TanAlpha0*si.TanAlpha0;
    si.L           = L;
    si.aStarEq     = aStarEq;
//    si.dB          = dB;
    si.BwFuncData  = BwFuncData;
    si.BwFunc      = BwFunc;
    si.w1          = w1;           // Lower freq cuttof.
    si.w2          = w2;           // Upper freq cuttof.
    si.wm          = wm;           // Frequency of max wave power.
    si.dw          = dw;           // Frequency bandwidth (at equator).
    si.x1          = x1;           // Lower tan(wave normal angle) cutoff.
    si.x2          = x2;           // Upper tan(wave normal angle) cutoff.
    si.numberOfWaveNormalAngleDistributions = (int) numberOfWaveNormalAngleDistributions; 
    si.xm          = (double *) xmArray;           // tan(wave normal angle) with max wave power.
    si.dx          = (double *) dxArray;           // spread in tan(wave normal angle).
    si.weightsOnWaveNormalAngleDistributions = (double *) weightsOnWaveNormalAngleDistributions;           
    si.MaxWaveLat  = MaxWaveLat*RadPerDeg;   // Assume there are no waves at +/-MaxWaveLat.
    si.Rho         = Rho;
    si.Directions  = Directions;

    /*
     *  Set integration limits
     */
    a = 0.0;                                            // radians
    double tempb;
    tempb = acos(Lgm_CdipMirrorLat( si.SinAlpha0 ));     // radians
    if (acos(sqrt(1.0/L)) < tempb) {B=acos(sqrt(1.0/L));} else {B=tempb;}  // Modified by Greg Cunningham on 10/14/2014 so that integral doesn't extend to inside earth for low L
    b = ( B < si.MaxWaveLat ) ? B : si.MaxWaveLat;
    if ( fabs(a-b) < 1e-9 ) {

        *Daa_ba = 0.0;
        *Dap_ba = 0.0;
        *Dpp_ba = 0.0;
        return( 0 );

    } else {

        /*
         *  Perform integrations. The points[] array contains a list of points
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
         *      T  = T0 - 0.5*(T0-T1)*(si.SinAlpha0 + sqrt(si.SinAlpha0));
         */
        npts2 = 2;
        //dqagp( CdipIntegrand_Sb, (_qpInfo *)&si, a, b, npts2, points, epsabs, epsrel, &T, &abserr, &neval, &ier, limit, lenw, &last, iwork, work );
        dqags( CdipIntegrand_Sb, (_qpInfo *)&si, a, B, epsabs, epsrel, &T, &abserr, &neval, &ier, limit, lenw, &last, iwork, work, VerbosityLevel );
	Lgm_SimpleRiemannSum( SummersIntegrand_Gaa, (_qpInfo *)&si, a, b, Daa_ba, VerbosityLevel );
/*	Comment all the old stuff out since we are just doing a very simple Riemann sum to compare with old results
        Lgm_SummersFindCutoffs2( SummersIntegrand_Gaa, (_qpInfo *)&si, FALSE, a, b, &a_new, &b_new );
printf("done finding cutoffs\n");
        npts2 = 2 + Lgm_SummersFindSingularities( SummersIntegrand_Gaa, (_qpInfo *)&si, FALSE, a_new, b_new, &points[1], &ySing );
printf("done finding singularity\n");
        if ( b_new > a_new ) {
            if ( npts2 > 2 ) {
printf("There is an interior singularity: Alpha0=%g, Ek=%g\n", Alpha0, Ek);
                dqagp( SummersIntegrand_Gaa, (_qpInfo *)&si, a_new, b_new, npts2, points, epsabs, epsrel, Daa_ba, &abserr, &neval, &ier, limit, lenw, &last, iwork, work, VerbosityLevel );
printf("done finding integral\n");
            } else {
                dqags( SummersIntegrand_Gaa, (_qpInfo *)&si, a_new, b_new, epsabs, epsrel, Daa_ba, &abserr, &neval, &ier, limit, lenw, &last, iwork, work, VerbosityLevel );
printf("done finding integral\n");
            }
        } else {
            *Daa_ba = 0.0;
        }
*/

	Lgm_SimpleRiemannSum( SummersIntegrand_Gap, (_qpInfo *)&si, a, b, Dap_ba, VerbosityLevel );
/*	Comment all the old stuff out since we are just doing a very simple Riemann sum to compare with old results
        Lgm_SummersFindCutoffs2( SummersIntegrand_Gap, (_qpInfo *)&si, FALSE, a, b, &a_new, &b_new );
        npts2 = 2; npts2 += Lgm_SummersFindSingularities( SummersIntegrand_Gap, (_qpInfo *)&si, FALSE, a_new, b_new, &points[1], &ySing );
        if ( b_new > a_new ) {
            if ( npts2 > 2 ) {
                dqagp( SummersIntegrand_Gap, (_qpInfo *)&si, a_new, b_new, npts2, points, epsabs, epsrel, Dap_ba, &abserr, &neval, &ier, limit, lenw, &last, iwork, work, VerbosityLevel );
//printf("There is an interior singularity: Alpha0=%g, Ek=%g\n", Alpha0, Ek);
            } else {
                dqags( SummersIntegrand_Gap, (_qpInfo *)&si, a_new, b_new, epsabs, epsrel, Dap_ba, &abserr, &neval, &ier, limit, lenw, &last, iwork, work, VerbosityLevel );
            }
        } else {
            *Dap_ba = 0.0;
        }
*/

	Lgm_SimpleRiemannSum( SummersIntegrand_Gpp, (_qpInfo *)&si, a, b, Dpp_ba, VerbosityLevel );
/*	Comment all the old stuff out since we are just doing a very simple Riemann sum to compare with old results
        Lgm_SummersFindCutoffs2( SummersIntegrand_Gpp, (_qpInfo *)&si, FALSE, a, b, &a_new, &b_new );
        npts2 = 2 + Lgm_SummersFindSingularities( SummersIntegrand_Gpp, (_qpInfo *)&si, FALSE, a_new, b_new, &points[1], &ySing );
        if ( b_new > a_new ) {
            if ( npts2 > 2 ) {
                dqagp( SummersIntegrand_Gpp, (_qpInfo *)&si, a_new, b_new, npts2, points, epsabs, epsrel, Dpp_ba, &abserr, &neval, &ier, limit, lenw, &last, iwork, work, VerbosityLevel );
//printf("There is an interior singularity: Alpha0=%g, Ek=%g\n", Alpha0, Ek);
            } else {
                dqags( SummersIntegrand_Gpp, (_qpInfo *)&si, a_new, b_new, epsabs, epsrel, Dpp_ba, &abserr, &neval, &ier, limit, lenw, &last, iwork, work, VerbosityLevel );
            }
        } else {
            *Dpp_ba = 0.0;
        }
*/
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
int Lgm_SummersDxxDerivsBounceAvg( int DerivScheme, double ha, int Version, double Alpha0,  double Ek,  double L,  void *BwFuncData, double (*BwFunc)(), double n1, double n2, double n3, double aStarEq,  int Directions, double w1, double w2, double wm, double dw, int WaveMode, int Species, double MaxWaveLat, double *dDaa,  double *dDap) {

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
        Lgm_SummersDxxBounceAvg( Version, a, Ek, L, BwFuncData, BwFunc, n1, n2, n3, aStarEq, Directions, w1, w2, wm, dw, WaveMode, Species, MaxWaveLat, &Daa_ba, &Dap_ba, &Dpp_ba );
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
     *  si->x1          lower cutoff in tan(wave normal angle).
     *  si->x2          upper cutoff in tan(wave normal angle).
     *  si->xm          tan(wave normal angle) with max power.
     *  si->dx          spread in tan(wave normal angle) (Semi bandwidth is dx).
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
                   si->Rho, si->Sig, x1, x2, xm, dx, si->Lambda, si->s, aStar, si->Directions );

    } else if ( si->Version == LGM_SUMMERS_2007 ) {

        Daa = Lgm_SummersDaaLocal_2007( SinAlpha2, si->E, dBoverB2, BoverBeq, Omega_e, Omega_Sig,
                    si->Rho, si->Sig, x1, x2, xm, dx, si->Lambda, si->s, si->n1, si->n2, si->n3, aStar, si->Directions );

    } else if ( si->Version == LGM_GLAUERT_AND_HORNE_HIGH_FREQ ) {

	    Daa = Lgm_GlauertAndHorneHighFrequencyDiffusionCoefficients_Local(SinAlpha2, si->E, dBoverB2, BoverBeq, Omega_e, Omega_Sig, 
		    si->Rho, si->Sig, x1, x2, xm, dx, si->x1, si->x2, si->numberOfWaveNormalAngleDistributions, 
            	    (double *) si->xm, (double *) si->dx, (double *) si->weightsOnWaveNormalAngleDistributions, 
            	    si->Lambda, si->s, aStar, si->Directions, (int) 0, si->nNw, si->nPlasmaParameters, si->aStarMin, si->aStarMax, si->Nw);

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

//printf("%20.11e %20.11e\n", Lat*DegPerRad, f);
//printf("%g %g\n", Lat*DegPerRad, f);
//printf("%g %g %g %g\n", Lat*DegPerRad, Daa, TanAlpha2, si->TanAlpha02);
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
                   si->Rho, si->Sig, x1, x2, xm, dx, si->Lambda, si->s, aStar, si->Directions );

    } else if ( si->Version == LGM_SUMMERS_2007 ) {

        Dap = Lgm_SummersDapLocal_2007( SinAlpha2, si->E, dBoverB2, BoverBeq, Omega_e, Omega_Sig,
                    si->Rho, si->Sig, x1, x2, xm, dx, si->Lambda, si->s, si->n1, si->n2, si->n3, aStar, si->Directions );

    } else if ( si->Version == LGM_GLAUERT_AND_HORNE_HIGH_FREQ ) {

	    Dap = Lgm_GlauertAndHorneHighFrequencyDiffusionCoefficients_Local(SinAlpha2, si->E, dBoverB2, BoverBeq, Omega_e, Omega_Sig, 
		    si->Rho, si->Sig, x1, x2, xm, dx, si->x1, si->x2, si->numberOfWaveNormalAngleDistributions, 
		    (double *) si->xm, (double *) si->dx, (double *) si->weightsOnWaveNormalAngleDistributions, 
                    si->Lambda, si->s, aStar, si->Directions, (int) 1, si->nNw, si->nPlasmaParameters, si->aStarMin, si->aStarMax, si->Nw);

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
                   si->Rho, si->Sig, x1, x2, xm, dx, si->Lambda, si->s, aStar, si->Directions );

    } else if ( si->Version == LGM_SUMMERS_2007 ) {

        Dpp = Lgm_SummersDppLocal_2007( SinAlpha2, si->E, dBoverB2, BoverBeq, Omega_e, Omega_Sig,
                    si->Rho, si->Sig, x1, x2, xm, dx, si->Lambda, si->s, si->n1, si->n2, si->n3, aStar, si->Directions );

    } else if ( si->Version == LGM_GLAUERT_AND_HORNE_HIGH_FREQ ) {

	Dpp = Lgm_GlauertAndHorneHighFrequencyDiffusionCoefficients_Local(SinAlpha2, si->E, dBoverB2, BoverBeq, Omega_e, Omega_Sig, 
		    si->Rho, si->Sig, x1, x2, xm, dx, si->x1, si->x2, si->numberOfWaveNormalAngleDistributions, 
		    (double *) si->xm, (double *) si->dx, (double *) si->weightsOnWaveNormalAngleDistributions, si->Lambda, si->s, aStar, si->Directions, (int) 2,
		    si->nNw, si->nPlasmaParameters, si->aStarMin, si->aStarMax, si->Nw);

    } else {

        printf("SummersIntegrand_Gpp: Unknown version of Summers Diff. Coeff Model. Version = %d\n", si->Version );
        return(-1);

    }


    /*
     * Finally the integrand.
     */
    f = Dpp * CosLat*v/sqrt(1.0-BoverBm);
//double localAlpha;
//localAlpha = asin(sqrt(SinAlpha2));
//printf("Lat=%g Dpp=%g, localAlpha=%g\n", Lat, f, localAlpha);
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
double Lgm_SummersDaaLocal( double SinAlpha2, double E, double dBoverB2, double BoverBeq, double Omega_e, double Omega_Sig, double Rho, double Sig, double xl, double xh, double xm, double dx, double Lambda, int s, double aStar, int Directions ) {

    int             nReal, nRoots, n;
    double          Gamma, Gamma2, Beta, Beta2, Mu, Mu2, BetaMu, BetaMu2, OneMinusBetaMu2;
    double          sEpsMinusOne, a1, a2, a3, a4, a, aa, apa, b;
    double complex  z1, z2, z3, z4, z[4];
    double          Daa, R, x0, y0, Ep1, Ep12, xms, xpse, u, arg, c1, c2, c3, c4;
    double          x, x2, x3, x4, y, g, F, DD, fac;
    double 	    regularizerForSingularity=0.0;

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
            Daa = ( Directions == LGM_FRWD_BKWD ) ? 2.0*DD : DD;
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
                if (  (Directions == LGM_FRWD_BKWD) || ((Directions == LGM_FRWD)&&(y>=0.0)) || ((Directions == LGM_BKWD)&&(y<=0.0)) ) {
                    g = x4 + c1*x3 + c2*x2 + c3*x + c4;
                    xms  = x - s;
                    xpse = x + s*LGM_EPS;
                    F = y*xms*xms*xpse*xpse/(x*g);
                    u = 1.0 - x*Mu/(y*Beta);
                    arg = (x-xm)/dx;
		    /* modified fabs(BetaMu-F) to take away the singularity (Cunningham 1/30/2013) */
		    if (fabs(BetaMu-F) < regularizerForSingularity) { Daa += u*u *fabs(F) / ( dx * regularizerForSingularity ) * exp( -arg*arg ); }
		    else { Daa += u*u *fabs(F) / ( dx * fabs(BetaMu - F) ) * exp( -arg*arg ); }
                }
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
 *
 *      \return         Local value of Daa
 *
 *      \author         M. Henderson
 *      \date           2010-2011
 *
 */
double Lgm_SummersDapLocal( double SinAlpha2, double E, double dBoverB2, double BoverBeq, double Omega_e, double Omega_Sig, double Rho, double Sig, double xl, double xh, double xm, double dx, double Lambda, int s, double aStar, int Directions ) {

    int             nReal, nRoots, n;
    double          Gamma, Gamma2, Beta, Beta2, Mu, Mu2, BetaMu, BetaMu2, OneMinusBetaMu2;
    double          sEpsMinusOne, a1, a2, a3, a4, a, aa, apa, b, SinAlpha;
    double complex  z1, z2, z3, z4, z[4];
    double          Dap, R, x0, y0, Ep1, Ep12, xms, xpse, u, arg, c1, c2, c3, c4;
    double          x, x2, x3, x4, y, g, F, DD, fac;
    double 	    regularizerForSingularity=0.0;

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
            Dap   = ( Directions == LGM_FRWD_BKWD ) ? 2.0*DD : DD;
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
                if (  (Directions == LGM_FRWD_BKWD) || ((Directions == LGM_FRWD)&&(y>=0.0)) || ((Directions == LGM_BKWD)&&(y<=0.0)) ) {
                    g = x4 + c1*x3 + c2*x2 + c3*x + c4;
                    xms  = x - s;
                    xpse = x + s*LGM_EPS;
                    F = y*xms*xms*xpse*xpse/(x*g);

                    u = 1.0 - x*Mu/(y*Beta);
                    arg = (x-xm)/dx;
		    /* modified fabs(BetaMu-F) to take away the singularity (Cunningham 1/30/2013) */
		    if (fabs(BetaMu-F) < regularizerForSingularity) { Dap += x/y * u *fabs(F) / ( dx * regularizerForSingularity) * exp( -arg*arg ); }
		    else { Dap += x/y * u *fabs(F) / ( dx * fabs(BetaMu - F) ) * exp( -arg*arg ); }
                }
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
double Lgm_SummersDppLocal( double SinAlpha2, double E, double dBoverB2, double BoverBeq, double Omega_e, double Omega_Sig, double Rho, double Sig, double xl, double xh, double xm, double dx, double Lambda, int s, double aStar, int Directions ) {

    int             nReal, nRoots, n;
    double          Gamma, Gamma2, Beta, Beta2, Mu, Mu2, BetaMu, BetaMu2, OneMinusBetaMu2;
    double          sEpsMinusOne, a1, a2, a3, a4, a, aa, apa, b;
    double complex  z1, z2, z3, z4, z[4];
    double          Dpp, R, x0, y0, Ep1, Ep12, xms, xpse, arg, c1, c2, c3, c4;
    double          x, x2, x3, x4, y, g, F, DD, fac;
    double 	    regularizerForSingularity=0.0;

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
            Dpp   = ( Directions == LGM_FRWD_BKWD ) ? 2.0*DD : DD;
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
//BetaMu *= -1.0;
        Mu     = sqrt(Mu2);
//Mu *= -1.0;
        for ( Dpp=0.0, n=0; n<nRoots; n++ ){


            x = creal( z[n] ); x2 = x*x; x3 = x2*x; x4 = x2*x2;
            y = ( x+a )/BetaMu;
            if ((x>xl)&&(x<xh)) {
                if (  (Directions == LGM_FRWD_BKWD) || ((Directions == LGM_FRWD)&&(y>=0.0)) || ((Directions == LGM_BKWD)&&(y<=0.0)) ) {
                    g = x4 + c1*x3 + c2*x2 + c3*x + c4;
                    xms  = x - s;
                    xpse = x + s*LGM_EPS;
                    F = y*xms*xms*xpse*xpse/(x*g);

                    arg = (x-xm)/dx;
		    /* modified fabs(BetaMu-F) to take away the singularity (Cunningham 1/30/2013) */
		    if (fabs(BetaMu-F) < regularizerForSingularity) {Dpp += x*x/(y*y) *fabs(F) / ( dx * regularizerForSingularity) * exp( -arg*arg );}
		    else {Dpp += x*x/(y*y) *fabs(F) / ( dx * fabs(BetaMu - F) ) * exp( -arg*arg );}
                }
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
 *          \f{eqnarray*}
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
 *          \f{eqnarray*}
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
 *          \f{eqnarray*}
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
 *          \f{eqnarray*}
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
 *  But, note that by definition, \f$ q_1 = 1\f$. Therefore \f$ h_0 = 0\f$. Define \f$ H = xG = x\Sigma_{i=0}^{8} g_i x^i \f$, so that;
 *
 *          \f{eqnarray*}
 *              g_0 &=& 2 e^5 ( 21 - q_2 + e \alpha^* )\\
 *              g_1 &=& e^4 s \left( 3*(203 - q_3) - 4*e^2*\alpha^* + 4*e*(21*\alpha^* + q_2) \right)\\
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
 *      \f{eqnarray*}
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
 *      \f{eqnarray*}
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
 *      \f{eqnarray*}
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
double Lgm_SummersDaaLocal_2007( double SinAlpha2, double E, double dBoverB2, double BoverBeq, double Omega_e, double Omega_Sig, double Rho, double Sig, double xl, double xh, double xm, double dx, double Lambda, int s, double n1, double n2, double n3, double aStar, int Directions ) {

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
            Daa = ( Directions == LGM_FRWD_BKWD ) ? 2.0*DD : DD;
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
                if (  (Directions == LGM_FRWD_BKWD) || ((Directions == LGM_FRWD)&&(y>=0.0)) || ((Directions == LGM_BKWD)&&(y<=0.0)) ) {
                    g = g0 + g1*x + g2*x2 + g3*x3 + g4*x4 + g5*x5 + g6*x6 + g7*x7 + g8*x8;
                    sss = (x - s)*(x + s*LGM_EPS)*(4.0*x + s*LGM_EPS)*(16.0*x + s*LGM_EPS);

                    F = 2.0*y*aStar*sss*sss/(x*g);

                    u = 1.0 - x*Mu/(y*Beta);
                    arg = (x-xm)/dx;
                    Daa += u*u *fabs(F) / ( dx * fabs( BetaMu - F ) ) * exp( -arg*arg );
                }
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
double Lgm_SummersDapLocal_2007( double SinAlpha2, double E, double dBoverB2, double BoverBeq, double Omega_e, double Omega_Sig, double Rho, double Sig, double xl, double xh, double xm, double dx, double Lambda, int s, double n1, double n2, double n3, double aStar, int Directions ) {

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
            Dap   = ( Directions == LGM_FRWD_BKWD ) ? 2.0*DD : DD;
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
                if (  (Directions == LGM_FRWD_BKWD) || ((Directions == LGM_FRWD)&&(y>=0.0)) || ((Directions == LGM_BKWD)&&(y<=0.0)) ) {
                    g = g0 + g1*x + g2*x2 + g3*x3 + g4*x4 + g5*x5 + g6*x6 + g7*x7 + g8*x8;
                    sss = (x - s)*(x + s*LGM_EPS)*(4.0*x + s*LGM_EPS)*(16.0*x + s*LGM_EPS);

                    F = 2.0*y*aStar*sss*sss/(x*g);

                    u = 1.0 - x*Mu/(y*Beta);
                    arg = (x-xm)/dx;
                    // if fabs(BetaMu - F) gets too small, we will get a singularity.
                    Dap += x/y * u *fabs(F) / ( dx * fabs( BetaMu - F ) ) * exp( -arg*arg );
                }
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
double Lgm_SummersDppLocal_2007( double SinAlpha2, double E, double dBoverB2, double BoverBeq, double Omega_e, double Omega_Sig, double Rho, double Sig, double xl, double xh, double xm, double dx, double Lambda, int s, double n1, double n2, double n3, double aStar, int Directions ) {

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
            Dpp   = ( Directions == LGM_FRWD_BKWD ) ? 2.0*DD : DD;
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
        //Mu     = sqrt(Mu2);
        for ( Dpp=0.0, n=0; n<nRoots; n++ ){

            x = creal( z[n] ); x2 = x*x; x3 = x2*x; x4 = x2*x2; x5 = x3*x2; x6 = x3*x3; x7 = x4*x3; x8 = x4*x4;
            y = ( x+a )/BetaMu;
            if ((x>xl)&&(x<xh)) {
                if (  (Directions == LGM_FRWD_BKWD) || ((Directions == LGM_FRWD)&&(y>=0.0)) || ((Directions == LGM_BKWD)&&(y<=0.0)) ) {
                    g = g0 + g1*x + g2*x2 + g3*x3 + g4*x4 + g5*x5 + g6*x6 + g7*x7 + g8*x8;
                    sss = (x - s)*(x + s*LGM_EPS)*(4.0*x + s*LGM_EPS)*(16.0*x + s*LGM_EPS);

                    F = 2.0*y*aStar*sss*sss/(x*g);

                    //u = 1.0 - x*Mu/(y*Beta);
                    arg = (x-xm)/dx;
                    Dpp += x*x/(y*y) *fabs(F) / ( dx * fabs( BetaMu - F ) ) * exp( -arg*arg );
                }
            }

        }
        Dpp *= fac;

    }

    return(Dpp);

}


/**
 *  \brief
 *      Attempts to find the location of the maximum of the function, f, that is closest to the left-hand edge of the specified range in latitude.
 *
 *  \details
 *      If the smallest increment in latitude (0.01 degrees) is larger than the size of the interval, the function is identically zero, or the 
 *      maximum of the function is at the left- or right-hand edge, then the routine returns '0' and stores Lat1 in *x.  Otherwise, the location 
 *      of the maximum is stored in *x and the value of f at the maximum is stored in *y.  If the maximum is a singularity, then the routine 
 *      returns a '1'; otherwise the routine returns a '0'.
 */
int Lgm_SummersFindFirstMaximumFromLeft( double  (*f)( double, _qpInfo *), _qpInfo *qpInfo, int Verbose, double Lat0, double Lat1, double *x, double *y ) {

    int     done, founda, foundb, foundc, Found;
    double  a, b, c, d, fa, fb, fc, fd, fbb, inc, D1, D2, DELTA;

    inc = 0.01*RadPerDeg;		//alternatively, tried 0.1*RadPerDeg to try to speed things up, but this misses local maxima

//printf("In Lgm_SummersFindSingularities: Lat0=%g, Lat1=%g\n", Lat0, Lat1);
    /*
     * Find a non-zero starting point
     */
    a = Lat0; done = FALSE; founda = FALSE;
    while ( !done ) {
        a += inc;
        if ( a > Lat1 ){
            if (Verbose) printf("Warning: Line %d in %s Lat1-Lat0 too small?\n", __LINE__, __FILE__);
	    (*x) = Lat1; 		//added by Greg Cunningham on 1-28-2013 so that *x contains the right-hand side of the bracket
            return(FALSE);
            //done = TRUE;
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
            if (Verbose) printf("Warning: Line %d in %s Function monotonically non-increasing over entire interval\n", __LINE__, __FILE__);
	        (*x) = Lat1; 		//added by Greg Cunningham on 1-28-2013 so that *x contains the right-hand side of the bracket
            return(FALSE);
            //done = TRUE;
        } else {
            fb = fabs( f( b, qpInfo ) );
            if ( fb > fa ) {
                foundb = TRUE;
                done   = TRUE;
            } else {		//added by Greg Cunningham on 1-28-2013; since f is monotonically non-increasing, update (a,fa) so that it represents the minimum 
		        a = b; fa=fb; 
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
            if (Verbose) printf("Warning: Line %d in %s Function monotonically non-decreasing over entire interval\n", __LINE__, __FILE__);
	        (*x) = Lat1;		//added by Greg Cunningham on 1-28-2013 so that *x contains the right-hand side of the bracket
            return(FALSE);
            //done = TRUE;
        } else {
            fc = fabs( f( c, qpInfo ) );
            if ( fc < fb ) {
                foundc = TRUE;
                done   = TRUE;
            } else {		//added by Greg Cunningham on 1-28-2013; since f is monotonically non-decreasing, update (b,fb) so that it represents the peak 
		        b = c; fb = fc;  
		    }
        }
    }

if(Verbose) {printf("Found Bracket: %20.11e %20.11e %20.11e %20.11e %20.11e %20.11e\n", a, b, c, fa, fb, fc );}


    done = FALSE;
    while (!done) {

        D1 = b-a;
        D2 = c-b;

        if ( fabs(c-a) < 1e-14 ) {
        //if ( fabs(c-a) < 1e-6 ) { an alternative that Greg Cunningham tried on 1-28-2013 to speed things up and be consistent with exclusion of the singularity
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
if(Verbose) {printf("a, b, c, d = %20.11e %20.11e %20.11e %20.11e; fa, fb, fc, fd=%20.11e %20.11e %20.11e %20.11e\n", a, b, c, d, fa, fb, fc, fd);}

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
	fbb = fabs(f( b-1e-4, qpInfo ));
        DELTA = fb - fbb;
        Found = ( DELTA > 1.0 ) ?  1 : 0;  
        //Found = ( DELTA > (0.1*fmax(fb, fbb) )) ?  1 : 0;  an alternative relative test that Greg Cunningham tried on 1-28-2013 to speed things up
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
 *      Finds all of the singularities in a specified latitude range.
 *
 *  \details
 */

int Lgm_SummersFindSingularities( double  (*f)( double, _qpInfo *), _qpInfo *qpInfo, int Verbose, double Lat0, double Lat1, double *x, double *y ) {
	int 	thisn, npts;

	npts = 0;
	while ((Lat0 < Lat1)&&(npts<8))
		{
		thisn = Lgm_SummersFindFirstMaximumFromLeft( f, qpInfo, Verbose, Lat0, Lat1, x+npts, y );
		Lat0 = (double) (*(x+npts));
		npts += thisn;
		}
	return(npts);
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
                //if ( fabs(x1-x0) < 0.001 ) {  an alternative tried by Greg Cunningham on 1-28-2013 to improve speed and be consistent with precision on integral
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
                //if ( fabs(x1-x0) < 0.001 ) {  an alternative tried by Greg Cunningham on 1-28-2013 to improve speed and be consistent with precision on integral
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

	double 	inc;
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

/*
 *  	The following routine is written to have the same interface as Lgm_SummersDxxLocal so that the code can be easily inserted alongside the Summers' 
 *  	approach.   Some of the parameters needed for the Summers' approach mean something different for the Glauert and Horne approach, e.g. s=+1 means 
 *  	R-mode in Summers but s is used to define the sign on the cyclotron freq in Glauert and Horne.  Currently, the input 's' parameter is ignored since 
 *  	only R-mode waves are used below.  The Summers' Lambda parameter is used to define the 's' needed for Glauert and Horne.  
 *  	Other paramters are ignored because my Glauert and Horne code only works for electrons, e.g. Sig and Omega_Sig; need to extend the code
 *  	to handle protons in addition to electrons.
 */

double Lgm_GlauertAndHorneHighFrequencyDiffusionCoefficients_Local( double SinAlpha2, double E, double dBoverB2, double BoverBeq, 
        double Omega_e, double Omega_Sig, double Rho, double Sig, double wxl, double wxh, double wxm, double wdx, double xmin, double xmax, 
        int numberOfWaveNormalAngleDistributions, double *xmArray, double *dxArray, double *weightsOnWaveNormalAngleDistributions, 
        double Lambda, int sIn, double aStar, int Directions, int tensorFlag, int nNw, int nPlasmaParameters, double aStarMin, double aStarMax, double *Nwglobal) {

	double 				alpha, KE, wce, wlc, wuc, wm, dw, Bw, B;
	int	    			gsl_err, iTheta, i, s;
	double				theta, wpe, totalDaa, dTanTheta, tanTheta, localDaa, Daa, tempDaa, tempxmin, tempxmax;
	struct inequalityParameters 	inequalityParamsLC, inequalityParamsUC;
	double				gamma, wn, vpar, thetaStart[4], thetaEnd[4];
	int				indexIntoNormalizer,j,nIntervals;
	double				*Nw;

//printf("In Lgm_GlauertAndHorneHighFrequencyDiffusionCoefficients_Local: SinAlpha2=%g, E=%g, dBoverB2=%g, BoverBeq=%g, Omega_e=%g, Omega_Sig=%g, Rho=%g, Sig=%g, wxl=%g, wxh=%g, wxm=%g, wdx=%g, xmin=%g, xmax=%g, Lambda=%g, s=%d, aStar=%g, Directions=%d, tensorFlag=%d\n", SinAlpha2, E, dBoverB2, BoverBeq, Omega_e, Omega_Sig, Rho, Sig, wxl, wxh, wxm, wdx, xmin, xmax, Lambda, sIn, aStar, Directions, tensorFlag);
//printf("                       :number of wave normal angle distributions to sum together is %d, with the following parameters (xm, dx, weight)\n", numberOfWaveNormalAngleDistributions);
//for (i=0; i<numberOfWaveNormalAngleDistributions; i++) {printf("                      %g %g %g\n", xmArray[i], dxArray[i], weightsOnWaveNormalAngleDistributions[i]);}

/*	Convert from 'Summers' parameters to 'Glauert and Horne' parameters  */
	if ((SinAlpha2>=0.0) && (SinAlpha2<=1.0)) {alpha = asin(sqrt(SinAlpha2));}//FIXME: this means we can only handle particles with parallel momentum in the same direction as the field
	if (Lambda == -1) {KE = E*LGM_Ee0; s=-1;}				  // Lambda is set to -1 for electrons in Lgm_SummersDxxBounceAvg
	if (Lambda == LGM_EPS) {KE = E*LGM_Ep0; s=1;} 				  // Lambda is set to LGM_ESP for protons in Lgm_SummersDxxBounceAvg
	wce=Omega_e/(2.0*M_PI);							  // Mike uses all frequencies in units of rads/second; I used Hz for frequencies in my implementation of Glauert and Horne
	wlc=wxl*wce;								  // wxl is the lower frequency cutoff normalized by local electron gyrofrequency; wlc is the lower frequency cutoff in Hz
	wuc=wxh*wce;								  // wxh is the upper frequency cutoff normalized by local electron gyrofrequency; wuc is the upper frequency cutoff in Hz
	wm=wxm*wce;								  // wxm is the mean frequency normalized by local electron gyrofrequency; wm is the mean frequency in Hz
	dw=wdx*wce;								  // wdx is the frequency range normalized by local electron gyrofrequency; dw is the frequency range in Hz
	Bw = sqrt(dBoverB2*pow(Omega_e*LGM_ELECTRON_MASS*1e12/LGM_e, 2.0));	  // Bw is wave amplitude in pT, so we need to multiply dBoverB2 by B^2, where B is the field in pT

	indexIntoNormalizer = (int) trunc(nPlasmaParameters*log(aStar/aStarMin)/log(aStarMax/aStarMin));
if ((aStar<aStarMin)||(aStar>aStarMax)) {printf("aStar=%g, aStarMin=%g, aStarMax=%g, nPlasmaParameters=%d, indexIntoNormalizer=%d\n", aStar, aStarMin, aStarMax, nPlasmaParameters, indexIntoNormalizer);}
	if (indexIntoNormalizer<0) {indexIntoNormalizer =0;}
	if (indexIntoNormalizer>nPlasmaParameters-1) {indexIntoNormalizer=nPlasmaParameters-1;}
	Nw=Nwglobal+(nNw*indexIntoNormalizer);
 
/*	Need to insert Jay's test here and expand it to include the other ranges of wn (not just wn<wlc, which is what I have now).  May need to change how integration is done. */
/*	Perform the integration over theta: currently done naively, need to port to GSL? */
	struct localDiffusionCoefficientAtSpecificThetaGlauertAndHorneFunctionParams p2;
	p2.tensorFlag=(int) tensorFlag; p2.s=(int) s; p2.KE=(double) KE; p2.aStar=(double) aStar; p2.wce=(double) wce; 
	p2.wm=(double) wm; p2.dw=(double) dw; p2.wlc=(double) wlc; p2.wuc=(double) wuc; p2.xmin=(double) xmin; p2.xmax=(double) xmax; 
	p2.xmArray=(double *) xmArray; p2.dxArray=(double *) dxArray; p2.weightsOnWaveNormalAngleDistributions=(double *) weightsOnWaveNormalAngleDistributions; 
	p2.numberOfWaveNormalAngleDistributions=(int) numberOfWaveNormalAngleDistributions; 
	p2.Bw=(double) Bw; p2.alpha=(double) alpha; p2.Nw=(double *) Nw; p2.nNw=(int) nNw; p2.Dir=(int) Directions;
	totalDaa = (double) 0.0;
	gamma = 1.0+((p2.KE)/LGM_Ee0);
	vpar = cos(alpha)*sqrt(1-(1.0/(gamma*gamma)));
	j = computeInequalityParameters(wxl, p2.aStar, vpar, &inequalityParamsLC);
	j = computeInequalityParameters(wxh, p2.aStar, vpar, &inequalityParamsUC);

	for (i=-5; i<=5; i++) {

 		p2.nCyclotronLow=(int) i; p2.nCyclotronHigh=(int) i; 
//printf("iCyclotron number=%d\n", i);
/*		
 *		wn	doppler shift, s*n*wce/gamma, normalized by the gyrofrequency, wce, to give s*n/gamma
 *		wlc	lower cutoff frequency normalized by the gyrofrequency
 *		wuc	uppper cutoff frequency normalized by the gyrofrequency
 *		aStar	the squared ratio of the electron gyrofrequency to the plasma frequency, wpe
 *		mode	FIXME: should be =1 for R-mode, =-1 for L-mode but p2.x is -1, which is opposite what it should be. need to iron this out.
 *		vpar	the parallel velocity divided by the speed of light
 */
		double wn;

		wn = (p2.s)*i/gamma;
		j = computeWaveNormalAngleIntervalsWhereResonantRootsMayExist(wn, wxl, wxh, &inequalityParamsLC, &inequalityParamsUC, (int) 1, &nIntervals, thetaStart, thetaEnd);
		tempDaa=0.0;

		if (nIntervals !=0) {

			if (nIntervals == 1) {

				tempxmin = tan(thetaStart[0]); if (tempxmin < xmin) { tempxmin = xmin; }
				tempxmax = tan(thetaEnd[0]); if (tempxmax > xmax) { tempxmax = xmax; }

				if (tempxmax>tempxmin) {
					j = Lgm_SimpleRiemannSum( localDiffusionCoefficientAtSpecificThetaGlauertAndHorneWeightedByTanTheta,(_qpInfo *) &p2,(double) tempxmin,(double) tempxmax, &tempDaa, (int) 0);
			    }

	        } else {

				j = Lgm_SimpleRiemannSum( localDiffusionCoefficientAtSpecificThetaGlauertAndHorneWeightedByTanTheta, (_qpInfo *) &p2,(double) xmin,(double) xmax, &tempDaa, (int) 0);

            }

			totalDaa += tempDaa;

        }
    }
	return(totalDaa);
}

// these long names are pretty excessive -- it makes the code structure less readable
double localDiffusionCoefficientAtSpecificThetaGlauertAndHorneWeightedByTanTheta( double tanTheta, struct localDiffusionCoefficientAtSpecificThetaGlauertAndHorneFunctionParams *p2 ) {

	double 	unweightedDaa;
		
	unweightedDaa = (double) localDiffusionCoefficientAtSpecificThetaGlauertAndHorne( tanTheta, p2 );
	return( unweightedDaa*tanTheta );

}

/*
 *	The following routine computes the diffusion coefficient Daa (or Dap or Dpp) at a specific value of the tangent of the wave normal angle using the derivation in Glauert and Horne 2005.
 *	Note that the frequencies in Glauert and Horne use units of radians/second, whereas I use units of cycles/second (Hz).  The difference is mostly irrelevant since the 
 *	solution to the dispersion relation is done using normalized frequencies (frequency divided by gyrofrequency).  The appearance of w^2 in the weighting coefficient, where
 *	w_i is a resonant frequency, is effectively divided out by the k^2 in the N(w) term, so the units here also cancel out.  All of the terms dw/dk similarly cancel out in 
 *	terms of whether or not one uses frequency in Hz or radians/second.  The only place it remains, then, is in the term dw in the denominator of Asquared.  The term dw in Hz
 *	must be multiplied by 2*pi in order to match Glauert and Horne 2005.
 *
 *	tensorFlag 	specifies calculation of Daa (0), Dap (1) or Dpp (2).
 *
 *	Modified by Greg Cunningham on 6-1-2013 to eliminate resonant roots for which the wave normal angle is larger than the resonance cone angle.  This change was made at the same
 *	time as the change above to the normalizer routine, normalizerForWavePowerSpectrumFunction, that computes N(w).
 */

//double localDiffusionCoefficientAtSpecificThetaGlauertAndHorne(int tensorFlag, int s, double KE, double aStar, double wce, double xm, double xmin, double xmax, double dx, double wm, double dw, double wlc, double wuc, double Bw, double alpha, double tanTheta, int nCyclotronLow, int nCyclotronHigh, double *Nw, int nNw, int Dir) 
double localDiffusionCoefficientAtSpecificThetaGlauertAndHorne(double tanTheta, struct localDiffusionCoefficientAtSpecificThetaGlauertAndHorneFunctionParams *p)
{
    	int    	nRoots, i, j, nCyclotron, indexIntoNw, test;
    	double 	Gamma, Gamma2, Beta, Beta2, Mu2;
    	double 	a, a2;
    	double 	cosAlpha, SinAlpha, SinAlpha2,theta,cosTheta,cosTheta2,tanTheta2;
    	double 	z[9];
    	double 	x,y,nSquared,R,L,S,D,P;
    	double 	momentum;
	double 	term1, term2, term3a, term3b, term4, localDaa, Bsquared, Asquared, thisg, g, vpar, dwdk, phink2;
	double 	wpe, wi, E, result;
	int 	tensorFlag, s, nCyclotronLow, nCyclotronHigh, nNw, Dir;
	double 	KE, aStar, wce, xm, xmin, xmax, dx, weight, wm, dw, wlc, wuc, Bw, alpha, *Nw;

	tensorFlag=(int) p->tensorFlag; s=(int) p->s; KE=(double) p->KE; aStar=(double) p->aStar; wce=(double) p->wce; 
	wm=(double) p->wm; dw=(double) p->dw; wlc=(double) p->wlc; wuc=(double) p->wuc; xmin=(double) p->xmin; xmax=(double) p->xmax; 
	Bw=(double) p->Bw; alpha=(double) p->alpha; nCyclotronLow=(int) p->nCyclotronLow; nCyclotronHigh=(int) p->nCyclotronHigh; Nw=(double *) p->Nw; nNw=(int) p->nNw; Dir=(int) p->Dir;
	//printf("tanTheta=%g, tensorFlag=%d, s=%d, KE=%g, aStar=%g, wce=%g, wm=%g, dw=%g, wlc=%g, wuc=%g, xmin=%g, xmax=%g, Bw=%g, alpha=%g, nCyclotronLow=%d, nCyclotronHigh=%d, nNw=%d, Dir=%d\n", tanTheta, tensorFlag, s, KE, aStar, wce, wm, dw, wlc, wuc, xmin, xmax, Bw, alpha, nCyclotronLow, nCyclotronHigh, nNw, Dir); 
	if (fabs(alpha) < 1e-3) {alpha = 1e-3;} //FIXME: alpha=0 causes a divide by zero besselFunctionNormalizer; need to fix this special case at some point
	if (fabs(alpha-(M_PI/2.0)) < 1e-3) {alpha = M_PI/2.0-1e-3;} //FIXME: alpha=M_PI/2 causes a divide by zero in term2 and computeResonantRoots; need to fix this special case at some point
	if (fabs(tanTheta) < 1e-3) {tanTheta=1e-3;} //FIXME: tanTheta=0 is the special case of parallel propagation; should use simpler resonant root finder here; may get degenerate roots?
	if ((tensorFlag!=0)&&(tensorFlag!=1)&&(tensorFlag!=2)) 
		{
		printf("Warning: in localDaaAtSpecificThetaGlauertAndHorne, tensorFlag does not have a value of 0, 1 or 2\n");
		return(0.0); // no valid specification of tensorFlag
		}
	if ((tanTheta<xmin)||(tanTheta>xmax)) return(0.0); 		// no need to perform the calculation since the input wave normal angle is outside the limited range of applicability

/*	Define constants that are needed */
	theta = atan(tanTheta); 					//FIXME: this means that the function is only valid for wave normal angles between -pi/2 and pi/2
    	SinAlpha=sin(alpha);             
    	SinAlpha2=SinAlpha*SinAlpha;             
    	cosAlpha=cos(alpha);             
    	Mu2=1.0-SinAlpha2;
    	E=KE/0.511;               					//input normalized energy (KE/m0c2)
    	Gamma = E+1.0; Gamma2 = Gamma*Gamma;
    	momentum=sqrt((KE+LGM_Ee0)*(KE+LGM_Ee0)-(LGM_Ee0*LGM_Ee0));	//momentum in MeV/c
    	tanTheta2=tanTheta*tanTheta;
    	cosTheta=cos(theta);
    	cosTheta2=cosTheta*cosTheta;					
    	Beta2 = E*(E+2.0)*Mu2*cosTheta2/Gamma2;  // square of (v_par/c)*cos(theta)
	Beta=sqrt(Beta2);
	localDaa=0.0;
    	for (nCyclotron=nCyclotronLow; nCyclotron<=nCyclotronHigh; nCyclotron++)   		//input cyclotron mode
		{
/*		resonance condition depends on nCyclotron */
    		a  = s*nCyclotron/Gamma; a2 = a*a; 			//assume electrons so lambda absorbed in basic formula; s is 1 for R-mode, -1 for L-mode; FIXME: is this right?

/* 		Solve for resonant roots. Only the non-evanescent roots with specified polarization in the specified frequency passband are retained. */
		i=(int) computeResonantRootsUsingColdElectronPlasma(KE, theta, alpha, s, nCyclotron, aStar, wlc/wce, wuc/wce, &z[0], &nRoots);  
//printf("For nCyclotron=%d, number of resonant roots is %d\n", nCyclotron, nRoots);

/*		Sum over the resonant roots, weighting each one */
    		for ( i=0; i<nRoots; i++ ) 
			{
			x=z[i];
			y=(x-a)/Beta;											
			if (((y<0)&&((Dir == LGM_BKWD)||(Dir== LGM_FRWD_BKWD)))||((y>0)&&((Dir == LGM_FRWD)||(Dir == LGM_FRWD_BKWD))))// to make consistent with Summers
			  {
			  nSquared=(y/x)*(y/x);
			  R=1.0-((1.0/(aStar*x*x))*(x/(x-1.0)));
			  L=1.0-((1.0/(aStar*x*x))*(x/(x+1.0)));
			  S=0.5*(R+L);
			  D=0.5*(R-L);
			  P=1.0-(1.0/(aStar*x*x));

/*			  Modification by Greg Cunningham on 6-1-2013: test to see if the resonant root has a resonance cone angle that is less than the given wave normal angle */ 
			  test = ((-1.0*P/S*0.99)<(tanTheta2))&&(-1.0*P/S > 0.0);
			  if (test == 0) 
				{
/*			  	Evaluate the terms needed to weight this particular resonance */
			  	wi=x*wce;										//resonant frequency in Hz
			  	indexIntoNw = (int) (nNw*x);								// modified on 6/6/2014 to assume that the normalized frequency interval [0,1] is spanned
			  	if (indexIntoNw < 0) {indexIntoNw=0;}
			  	if (indexIntoNw > (nNw-1)) {indexIntoNw=nNw-1;}
			  	result=Nw[indexIntoNw]*wce*wce;								//wce*wce added on 6/6/2014 to compensate for fact that Nw[] doesn't include the wce^2 term
			  	phink2=(double) besselFunctionNormalizer(s, momentum, alpha, nCyclotron, x, y, P, R, L, S, theta);  	
			  	dwdk=(double) electronColdPlasmaGroupVelocity(aStar, tanTheta2, x, y);	//the group velocity
			  	vpar = (double) (momentum*1e6*Joules_Per_eV*cosAlpha)/(Gamma*LGM_c*LGM_ELECTRON_MASS);  //vpar in m/s; note that momentum is in MeV/c, so need to convert to kg*m/s
			  	g=0.0;
			  	for (j=0; j<p->numberOfWaveNormalAngleDistributions; j++)
					{
					xm=(double) p->xmArray[j]; dx=(double) p->dxArray[j]; weight = (double) p->weightsOnWaveNormalAngleDistributions[j];
			  		thisg=weight*exp(-1.0*(tanTheta-xm)*(tanTheta-xm)/(dx*dx));				//g(tanTheta)
					g += thisg;
					}	
			  	Asquared = (pow(Bw*1e-12, 2.0)/(2.0*M_PI*dw))*(2/sqrt(M_PI))/(gsl_sf_erf((wm-wlc)/dw)+gsl_sf_erf((wuc-wm)/dw));	//wave spectral intensity in T^2/Hz (actually, these units are somewhat of a misnomer since dw (which is in Hz) has to be multiplied by 2*pi in order to get it in units of radians/second instead of cycles/sec
			  	Bsquared=Asquared*exp(-1.0*(wi-wm)*(wi-wm)/(dw*dw));
			  	term1 = (LGM_e*LGM_e*wi*wi/(4.0*M_PI*(1.0+tanTheta2)*result));
			  	term2 = pow(((s*nCyclotron*wce/(Gamma*wi))-SinAlpha2)/cosAlpha, 2.0);			//FIXME: using wce assumes resonance condition for electrons; check for protons
			  	term3a = Bsquared*g*phink2;
			  	term3b = 1.0/fabs(vpar-(dwdk/cosTheta));						//FIXME: don't I need to check to make sure vpar is ne to dwdk/cosTheta?  
			  	//if (term3b > 3e-5) {term3b = 3e-5;}                                                   //regularize the singularity
			  	term4 = SinAlpha*cosAlpha/((s*nCyclotron*wce/(Gamma*wi))-SinAlpha2); 			//FIXME: using wce assumes resonance condition for electrons; check for protons
//printf("cyclotron number=%d resonant root number=%d root=%g Term1=%g term2=%g Bsquared=%g g=%g phink2=%g term3b=%g term4=%g\n", nCyclotron, i, wi, term1, term2, Bsquared, g, phink2, term3b, term4);
			  	if (tensorFlag == 0) { localDaa += term1*term2*term3a*term3b;}				// compute Daa
			  	if (tensorFlag == 1) { localDaa += term1*term2*term3a*term3b*term4;}			// compute Dap
			  	if (tensorFlag == 2) { localDaa += term1*term2*term3a*term3b*term4*term4;}		// compute Dpp
				}
			  }
			}
		}
	localDaa /= pow(momentum*1e6*Joules_Per_eV/LGM_c, 2.0);  							//divide Glauert/Horne Daa by momentum^2 in (kg*m/s)^2 to get Summers
//printf("localDaa=%g\n", localDaa);
	if (fpclassify(localDaa) != FP_NORMAL) {localDaa = 0.0;}							//eliminate the possible of returning a number that is not a number
	return(localDaa);
}

/*
 *	The following routine computes all 3 components of the diffusion tensor
 *	(Daa, Dap and Dpp) at a specific value of the tangent of the wave normal
 *	angle using the derivation in Glauert and Horne 2005.  Note that the
 *	frequencies in Glauert and Horne use units of radians/second, whereas I use
 *	units of cycles/second (Hz).  The difference is mostly irrelevant since the
 *	solution to the dispersion relation is done using normalized frequencies
 *	(frequency divided by gyrofrequency).  The appearance of w^2 in the
 *	weighting coefficient, where w_i is a resonant frequency, is effectively
 *	divided out by the k^2 in the N(w) term, so the units here also cancel out.
 *	All of the terms dw/dk similarly cancel out in terms of whether or not one
 *	uses frequency in Hz or radians/second.  The only place it remains, then,
 *	is in the term dw in the denominator of Asquared.  The term dw in Hz must
 *	be multiplied by 2*pi in order to match Glauert and Horne 2005.
 *
 *	Modified by Greg Cunningham on 6-1-2013 to eliminate resonant roots for
 *	which the wave normal angle is larger than the resonance cone angle.  This
 *	change was made at the same time as the change above to the normalizer
 *	routine, normalizerForWavePowerSpectrumFunction, that computes N(w).
 */
int localDiffusionCoefficientAtSpecificThetaGlauertAndHorneReturnTensorInPlace(double tanTheta, struct localDiffusionCoefficientAtSpecificThetaGlauertAndHorneFunctionParams *p, double *Daa, double *Dap, double *Dpp)
{
    	int    	nRoots, i, j, nCyclotron, indexIntoNw, test;
    	double 	Gamma, Gamma2, Beta, Beta2, Mu2;
    	double 	a, a2;
    	double 	cosAlpha, SinAlpha, SinAlpha2,theta,cosTheta,cosTheta2,tanTheta2;
    	double 	z[9];
    	double 	x,y,nSquared,R,L,S,D,P;
    	double 	momentum;
	double 	term1, term2, term3a, term3b, term4, localDaa, Bsquared, Asquared, thisg, g, vpar, dwdk, phink2;
	double 	wpe, wi, E, result;
	int 	s, nCyclotronLow, nCyclotronHigh, nNw, Dir;
	double 	KE, aStar, wce, xm, xmin, xmax, dx, weight, wm, dw, wlc, wuc, Bw, alpha, *Nw;
	double 	term0, commonFactor;

	s=(int) p->s; KE=(double) p->KE; aStar=(double) p->aStar; wce=(double) p->wce; 
	wm=(double) p->wm; dw=(double) p->dw; wlc=(double) p->wlc; wuc=(double) p->wuc; xmin=(double) p->xmin; xmax=(double) p->xmax; 
	Bw=(double) p->Bw; alpha=(double) p->alpha; nCyclotronLow=(int) p->nCyclotronLow; nCyclotronHigh=(int) p->nCyclotronHigh; Nw=(double *) p->Nw; nNw=(int) p->nNw; Dir=(int) p->Dir;
	//printf("tanTheta=%g, s=%d, KE=%g, aStar=%g, wce=%g, wm=%g, dw=%g, wlc=%g, wuc=%g, xmin=%g, xmax=%g, Bw=%g, alpha=%g, nCyclotronLow=%d, nCyclotronHigh=%d, nNw=%d, Dir=%d\n", tanTheta, s, KE, aStar, wce, wm, dw, wlc, wuc, xmin, xmax, Bw, alpha, nCyclotronLow, nCyclotronHigh, nNw, Dir); 
	if (fabs(alpha) < 1e-3) {alpha = 1e-3;} //FIXME: alpha=0 causes a divide by zero besselFunctionNormalizer; need to fix this special case at some point
	if (fabs(alpha-(M_PI/2.0)) < 1e-3) {alpha = M_PI/2.0-1e-3;} //FIXME: alpha=M_PI/2 causes a divide by zero in term2 and computeResonantRoots; need to fix this special case at some point
	if (fabs(tanTheta) < 1e-3) {tanTheta=1e-3;} //FIXME: tanTheta=0 is the special case of parallel propagation; should use simpler resonant root finder here; may get degenerate roots?
	if ((tanTheta<xmin)||(tanTheta>xmax)) {*Daa = 0.0; *Dap=0.0; *Dpp=0.0; return(-1);} 		// no need to perform the calculation since the input wave normal angle is outside the limited range of applicability

/*	Define constants that are needed */
	theta = atan(tanTheta); 					//FIXME: this means that the function is only valid for wave normal angles between -pi/2 and pi/2
    SinAlpha=sin(alpha);             
    SinAlpha2=SinAlpha*SinAlpha;             
    cosAlpha=cos(alpha);             
    Mu2=1.0-SinAlpha2;
    E=KE/0.511;               					//input normalized energy (KE/m0c2)
    Gamma = E+1.0; Gamma2 = Gamma*Gamma;
    momentum=sqrt((KE+LGM_Ee0)*(KE+LGM_Ee0)-(LGM_Ee0*LGM_Ee0));	//momentum in MeV/c
    tanTheta2=tanTheta*tanTheta;
    cosTheta=cos(theta);
    cosTheta2=cosTheta*cosTheta;					
    Beta2 = E*(E+2.0)*Mu2*cosTheta2/Gamma2;  // square of (v_par/c)*cos(theta)
	Beta=sqrt(Beta2);
	*Daa=0.0;
	*Dap=0.0;
	*Dpp=0.0;
    for (nCyclotron=nCyclotronLow; nCyclotron<=nCyclotronHigh; nCyclotron++) {  		//input cyclotron mode
/*		resonance condition depends on nCyclotron */
    		a  = s*nCyclotron/Gamma; a2 = a*a; 			//assume electrons so lambda absorbed in basic formula; s is 1 for R-mode, -1 for L-mode; FIXME: is this right?

/* 		Solve for resonant roots. Only the non-evanescent roots with specified polarization in the specified frequency passband are retained. */
		i=(int) computeResonantRootsUsingColdElectronPlasma(KE, theta, alpha, s, nCyclotron, aStar, wlc/wce, wuc/wce, &z[0], &nRoots);  
//printf("For nCyclotron=%d, number of resonant roots is %d\n", nCyclotron, nRoots);

/*		Sum over the resonant roots, weighting each one */
    		for ( i=0; i<nRoots; i++ ) 
			{
			x=z[i];
			y=(x-a)/Beta;											
			if (((y<0)&&((Dir == LGM_BKWD)||(Dir== LGM_FRWD_BKWD)))||((y>0)&&((Dir == LGM_FRWD)||(Dir == LGM_FRWD_BKWD))))// to make consistent with Summers
			  {
			  nSquared=(y/x)*(y/x);
			  R=1.0-((1.0/(aStar*x*x))*(x/(x-1.0)));
			  L=1.0-((1.0/(aStar*x*x))*(x/(x+1.0)));
			  S=0.5*(R+L);
			  D=0.5*(R-L);
			  P=1.0-(1.0/(aStar*x*x));

/*			  Modification by Greg Cunningham on 6-1-2013: test to see if the resonant root has a resonance cone angle that is less than the given wave normal angle */ 
			  test = ((-1.0*P/S*0.99)<(tanTheta2))&&(-1.0*P/S > 0.0);
			  if (test == 0) 
				{
/*			  	Evaluate the terms needed to weight this particular resonance */
			  	wi=x*wce;										//resonant frequency in Hz
			  	indexIntoNw = (int) (nNw*x);								// modified on 6/6/2014 to assume that the normalized frequency interval [0,1] is spanned
			  	if (indexIntoNw < 0) {indexIntoNw=0;}
			  	if (indexIntoNw > (nNw-1)) {indexIntoNw=nNw-1;}
			  	result=Nw[indexIntoNw]*wce*wce;								//wce*wce added on 6/6/2014 to compensate for fact that Nw[] doesn't include the wce^2 term
			  	phink2=(double) besselFunctionNormalizer(s, momentum, alpha, nCyclotron, x, y, P, R, L, S, theta);  	
			  	dwdk=(double) electronColdPlasmaGroupVelocity(aStar, tanTheta2, x, y);	//the group velocity
			  	vpar = (double) (momentum*1e6*Joules_Per_eV*cosAlpha)/(Gamma*LGM_c*LGM_ELECTRON_MASS);  //vpar in m/s; note that momentum is in MeV/c, so need to convert to kg*m/s
			  	g=0.0;
			  	for (j=0; j<p->numberOfWaveNormalAngleDistributions; j++)
					{
					xm=(double) p->xmArray[j]; dx=(double) p->dxArray[j]; weight = (double) p->weightsOnWaveNormalAngleDistributions[j];
			  		thisg=weight*exp(-1.0*(tanTheta-xm)*(tanTheta-xm)/(dx*dx));				//g(tanTheta)
					g += thisg;
					}	
			  	Asquared = (pow(Bw*1e-12, 2.0)/(2.0*M_PI*dw))*(2/sqrt(M_PI))/(gsl_sf_erf((wm-wlc)/dw)+gsl_sf_erf((wuc-wm)/dw));	//wave spectral intensity in T^2/Hz (actually, these units are somewhat of a misnomer since dw (which is in Hz) has to be multiplied by 2*pi in order to get it in units of radians/second instead of cycles/sec
			  	Bsquared=Asquared*exp(-1.0*(wi-wm)*(wi-wm)/(dw*dw));
				term0 = cosAlpha/((s*nCyclotron*wce/(Gamma*wi))-SinAlpha2); 				//FIXME: using wce assumes resonance condition for electrons; check for protons
			  	term1 = (LGM_e*LGM_e*wi*wi/(4.0*M_PI*(1.0+tanTheta2)*result));
			  	term2 = pow(1.0/term0, 2.0);								//FIXME: using wce assumes resonance condition for electrons; check for protons
			  	term3a = Bsquared*g*phink2;
			  	term3b = 1.0/fabs(vpar-(dwdk/cosTheta));						//FIXME: don't I need to check to make sure vpar is ne to dwdk/cosTheta?  
			  	//if (term3b > 3e-5) {term3b = 3e-5;}                                                   //regularize the singularity
			  	//term4 = SinAlpha*term0; 								//FIXME: using wce assumes resonance condition for electrons; check for protons
				commonFactor = term1*term2*term3a*term3b;
if (commonFactor < 0.0) {printf("commonFactor<0.0: x=%g, indexIntoNw=%d, result=%g, term1=%g, term2=%g, term3a=%g, term3b=%g, term4=%g\n", x,indexIntoNw, result, term1, term2, term3a, term3b, term4);}
//printf("cyclotron number=%d resonant root number=%d root=%g Term1=%g term2=%g Bsquared=%g g=%g phink2=%g term3b=%g term4=%g\n", nCyclotron, i, wi, term1, term2, Bsquared, g, phink2, term3b, term4);
			  	*Daa += term1*term2*term3a*term3b;				// compute Daa
			  	*Dap += term1*term3a*term3b*SinAlpha/term0;			// used to be term1*term3a*term3b*term4 compute Dap
			  	*Dpp += term1*term3a*term3b*SinAlpha*SinAlpha;			// used to be term1*term3a*term3b*term4*term4 compute Dpp
//printf("x=%g, indexIntoNw=%d, result=%g, term1=%g, term2=%g, term3a=%g, term3b=%g, term4=%g\n", x,indexIntoNw, result, term1, term2, term3a, term3b, term4);
				}
			  }
			}
		}
	*Daa /= pow(momentum*1e6*Joules_Per_eV/LGM_c, 2.0);  							//divide Glauert/Horne Daa by momentum^2 in (kg*m/s)^2 to get Summers
	*Dap /= pow(momentum*1e6*Joules_Per_eV/LGM_c, 2.0);  							//divide Glauert/Horne Daa by momentum^2 in (kg*m/s)^2 to get Summers
	*Dpp /= pow(momentum*1e6*Joules_Per_eV/LGM_c, 2.0);  							//divide Glauert/Horne Daa by momentum^2 in (kg*m/s)^2 to get Summers
//printf("localDaa=%g, localDap=%g, localDpp=%g\n", *Daa, *Dap, *Dpp);
/*	eliminate the possible of returning a number that is not a number, or a tensor that has negative determinant */
	if ((fpclassify(*Daa) != FP_NORMAL)||(fpclassify(*Dap) != FP_NORMAL)||(fpclassify(*Dpp) != FP_NORMAL)) 
		{
		*Daa = 0.0;
		*Dap = 0.0; 
		*Dpp = 0.0;
		}							
	return(0);
}

/*
 *	The following routine computes the frequencies, w_i, at which the dispersion relation and resonance condition are both simultaneously satisfied.
 *	The dispersion relation that is used is an approximation which assumes that the ion contribution is minimal (frequency is greater than 0.1 times
 *	the electron gyrofrequency); this may be called the cold electron magnetoplasma dispersion relation.  The inputs to this routine are 
 *	
 *		KE		the kinetic energy of the particle (in MeV), needed to compute the parallel velocity for the resonance condition
 *		theta		wave normal angle (in radians) FIXME: don't I need to do something special here if theta=0 or pi/2?
 *		alpha		particle pitch angle (in radians), needed to compute the parallel velocity for the resonance condition
 *		s		particle type (-1=electron, +1=proton)
 *		aStar		the 'cold plasma parameter': square of the ratio of the electron gyrofrequency to the electron plasmafrequency
 *		xlc		the lower normalized frequency (w/wce) limit of the 'passband' in which roots must be
 *		xuc		the higher normalized frequency (w/wce) limit of the 'passband' in which roots must be
 *		resonantRoots	an array pre-allocated to have at most 9 entries, into which the real part of the roots in normalized frequency are stored
 *		nRoots		the number of roots that are not evanescent (complex part=0), have non-zero real component, have the correct polarization, 
 *					and are in the specified normalized frequency passband
 *
 */
int computeResonantRootsUsingColdElectronPlasma(double KE, double theta, double alpha, int s, int nCyclotron, double aStar, double xlc, double xuc, double *resonantRoots, int *nRoots)
{
	double		E, Gamma, Gamma2, Beta, Beta2, Beta4, tanTheta, tanTheta2, cosTheta, cosTheta2, cosAlpha, sinAlpha, sinAlpha2, Mu2, a, a2, a3, a4, Coeff[9];
    	double complex  zz[9];
    	double          gsl_z[30];
	int		i,gsl_err;
	double 		x,y,nSquared,P,R,S,L,D,changeInPolarization,polarization;

/* 	Define constants.  */
    	E=KE/LGM_Ee0;               			// the dimensionless energy equal to the kinetic energy (in MeV) divided by the rest mass energy of the electron
    	Gamma = E+1.0; Gamma2 = Gamma*Gamma;		// the Lorentz factor and its square
    	tanTheta=tan(theta); tanTheta2=tanTheta*tanTheta;
    	cosTheta=cos(theta); cosTheta2=cosTheta*cosTheta;
    	sinAlpha=sin(alpha); sinAlpha2=sinAlpha*sinAlpha;             
	cosAlpha=sqrt(1.0-sinAlpha2);
    	Mu2=1.0-sinAlpha2;
    	Beta2 = E*(E+2.0)*Mu2*cosTheta2/Gamma2;  	// square of (v_par/c)*cos(theta)
    	Beta4 = Beta2*Beta2;
    	a  = s*nCyclotron/Gamma; a2 = a*a; a3=a2*a; a4=a2*a2; //s=-1 for electrons, s=+1 for protons

/*	
 *	Define the coefficients of the 8th degree polynomial whose roots in x=(wave freq)/(electron gyrofrequency) are solutions to the resonance condition 
 *	y=(x-a)/beta and the dispersion equation D(x)=0 for a cold electron magnetoplasma (the wave frequency is assumed to be between 0.1 and 0.9 times the 
 *	electron gyrofrequency and so the ion terms can be ignored).  After finding the roots, we throw away the evanescent terms (these waves do not 
 *	propagate freely because they have a complex component) and the polarizations that we don't want.
 *	
 *	The coefficients defined below were determined using the following command in Mathematica:
 *
 *	Expand[(((x^2*(x^2-1)-(x^2/aStar))*((x-a)/beta)^2-(x^2*(x+1)-(x/aStar))*(x^2*(x-1)-(x/aStar)))*(((x-a)/beta)^2-(x^2-(1/aStar)))*tanTheta^2) 
 *	+ ((x^2-(1/aStar))*(((x-a)/beta)^2*(x -1)-(x^2*(x-1)-(x/aStar)))*(((x-a)/beta)^2*(x+1)-(x^2*(x+1)-(x/aStar))))]
 *
 * 	Modified by Greg Cunningham on 6-12-2014 to get rid of the cos(alpha) terms in the denominator by multiplying through by Beta4.  The divide by zero terms, which
 * 	occur when Alpha=90, are in the original version below:
 *
    	Coeff[0] = a4/(Beta4*aStar); 		// x^0 term
    	Coeff[1] = -4.0*a3/(Beta4*aStar); 		// x^1 term
    	Coeff[2] = -1.0/pow(aStar,3.0)-(a4/Beta4)+(6.0*a2/(aStar*Beta4))-(a4/(aStar*Beta4))-(2.0*a2/(aStar*aStar*Beta2))-(2.0*a2/(aStar*Beta2))-(tanTheta2/(aStar*aStar*aStar))-(a4*tanTheta2/Beta4)-(a4*tanTheta2/(aStar*Beta4))-(2.0*a2*tanTheta2/(aStar*aStar*Beta2))-(a2*tanTheta2/(aStar*Beta2));
    	Coeff[3] = 4.0*a3/Beta4 - (4.0*a/(aStar*Beta4)) + (4.0*a3/(aStar*Beta4)) + (4.0*a/(aStar*aStar*Beta2)) + (4.0*a/(aStar*Beta2))+(4.0*a3*tanTheta2/Beta4)+(4.0*a3*tanTheta2/(aStar*Beta4))+(4.0*a*tanTheta2/(aStar*aStar*Beta2))+(2.0*a*tanTheta2/(aStar*Beta2));
    	Coeff[4] = (3.0/(aStar*aStar))+(1.0/aStar)-(6.0*a2/Beta4)+(a4/Beta4)+(1.0/(aStar*Beta4))-(6.0*a2/(aStar*Beta4))+(2.0*a2/Beta2)- (2.0/(aStar*aStar*Beta2))-(2.0/(aStar*Beta2))+(4.0*a2/(aStar*Beta2))+(3.0*tanTheta2/(aStar*aStar))+(tanTheta2/aStar)-(6.0*a2*tanTheta2/Beta4)+(a4*tanTheta2/Beta4)-(6.0*a2*tanTheta2/(aStar*Beta4))+(2.0*a2*tanTheta2/Beta2)-(2.0*tanTheta2/(aStar*aStar*Beta2))-(tanTheta2/(aStar*Beta2))+(4.0*a2*tanTheta2/(aStar*Beta2));
    	Coeff[5] = (4.0*a/Beta4)-(4.0*a3/Beta4)+(4.0*a/(aStar*Beta4))-(4.0*a/Beta2)-(8.0*a/(aStar*Beta2))+(4.0*a*tanTheta2/Beta4)-(4.0*a3*tanTheta2/Beta4)+(4.0*a*tanTheta2/(aStar*Beta4))-(4.0*a*tanTheta2/Beta2)-(8.0*a*tanTheta2/(aStar*Beta2));
    	Coeff[6] = -1.0-(3.0/aStar)-(1.0/Beta4)+(6.0*a2/Beta4)-(1.0/(aStar*Beta4))+(2.0/Beta2)-(2.0*a2/Beta2)+(4.0/(aStar*Beta2))-(tanTheta2)+(3.0*tanTheta2/aStar)-(tanTheta2/Beta4)+(6.0*a2*tanTheta2/Beta4)-(tanTheta2/(aStar*Beta4))+(2.0*tanTheta2/Beta2)-(2.0*a2*tanTheta2/Beta2)+(4.0*tanTheta2/(aStar*Beta2));
    	Coeff[7] = (-4.0*a/Beta4)+(4.0*a/Beta2)-(4.0*a*tanTheta2/Beta4)+(4.0*a*tanTheta2/Beta2);
    	Coeff[8] = 1.0+(1.0/Beta4)-(2/Beta2)+tanTheta2+(tanTheta2/Beta4)-(2.0*tanTheta2/Beta2);
 */
/*	
 *	For the special case nCyclotron=0, the coefficients of all of the odd terms in x are zero, and so we have something like c2*x^2+c4*x^4+c6*x^6+c8*x^8=0, which, 
 *	if we divide by x^2 and then substitute y=x^2, we have a cubic equation c2+c4*y+c6*y^2+c8*y^3 that has a closed form root-finding solution. Since we are only 
 *	interested in positive x (is that right?), we can ignore the negative roots, but put them in there anyways and then fill out the vector by putting zero as a root
 *	in the list of roots, twice. 
 */
	if (nCyclotron == 0)
		{
    		Coeff[0] = -1.0*Beta4/pow(aStar,3.0)-(tanTheta2*Beta4/(aStar*aStar*aStar));  	   // the coefficient on the x^2 term in the original expansion, which will now be the constant offset
    		Coeff[1] = (3.0*Beta4/(aStar*aStar))+(1.0*Beta4/aStar)+(1.0/(aStar))-(2.0*Beta2/(aStar*aStar))-(2.0*Beta2/(aStar))+(3.0*tanTheta2*Beta4/(aStar*aStar))+(tanTheta2*Beta4/aStar)-(2.0*tanTheta2*Beta2/(aStar*aStar))-(tanTheta2*Beta2/(aStar));										// the coefficient on the x^4 term in the original expansion, which will now be the coefficient on the y term
    		Coeff[2] =-1.0*Beta4-(3.0*Beta4/aStar)-(1.0)-(1.0/(aStar))+(2.0*Beta2)+(4.0*Beta2/(aStar))-(tanTheta2*Beta4)+(3.0*tanTheta2*Beta4/aStar)-(tanTheta2)-(tanTheta2/(aStar))+(2.0*tanTheta2*Beta2)+(4.0*tanTheta2*Beta2/(aStar));												// the coefficient on the x^6 term in the original expansion, which will now be the y^2 term
    		Coeff[3] = 1.0*Beta4+(1.0)-(2.0*Beta2)+(tanTheta2*Beta4)+(tanTheta2)-(2.0*tanTheta2*Beta2); // the coefficient on the x^8 term in the original expansion, which will now be the y^3 term
/* 		Solve for resonant roots and put into a complex array for further processing.  */
    		gsl_poly_complex_workspace *w = gsl_poly_complex_workspace_alloc( 4 );
    		gsl_err = gsl_poly_complex_solve( Coeff, 4, w, gsl_z );
    		gsl_poly_complex_workspace_free( w );
    		(*nRoots) = 0;
    		if (gsl_err != GSL_SUCCESS ) 
			{
			//printf("Problem with finding resonant roots!  Returning with nRoots=0\n");
			return((int) (-1));
			}	
    		for (i=0; i<3; i++) 
			{
			if (gsl_z[2*i] > 0.0) 
				{
				zz[i] = sqrt(gsl_z[2*i]) + gsl_z[2*i+1]*I;  // need to factor x=sqrt(y) into the root; don't care about imaginary part because we only count this root if it is real
				}
			else
				{
				zz[i] = 0.0+0.0*I;  			// the square root of a negative number has a complex component that is significant and so the root will be discarded as evanescent anyway, so just set to 0
				}
			zz[i+3] = 0.0+0.0*I;                             // there is a negative root, but it won't count toward anything so don't need this anymore: -1.0*gsl_z[2*i] - gsl_z[2*i+1]*I;
			}
		zz[7] = 0.0+0.0*I; zz[8] = 0.0+0.0*I;

		/* quick lookup from Ripoll 
		zz[0] = ((1.0/aStar)*(Gamma*Gamma-1.0)*cosAlpha*cosAlpha/(Gamma*Gamma))+(0.0*I); 
    		for (i=1; i<9; i++) {zz[i]=0.0+0.0*I;}
 		*/
		}
	else
		{
    		Coeff[0] = a4/aStar; 		// x^0 term
    		Coeff[1] = -4.0*a3/aStar; 		// x^1 term
    		Coeff[2] = -1.0*Beta4/pow(aStar,3.0)-(a4)+(6.0*a2/(aStar))-(a4/(aStar))-(2.0*a2*Beta2/(aStar*aStar))-(2.0*a2*Beta2/(aStar))-(tanTheta2*Beta4/(aStar*aStar*aStar))-(a4*tanTheta2)-(a4*tanTheta2/(aStar))-(2.0*a2*tanTheta2*Beta2/(aStar*aStar))-(a2*tanTheta2*Beta2/(aStar));
    		Coeff[3] = 4.0*a3- (4.0*a/(aStar)) + (4.0*a3/(aStar)) + (4.0*a*Beta2/(aStar*aStar)) + (4.0*a*Beta2/(aStar))+(4.0*a3*tanTheta2)+(4.0*a3*tanTheta2/(aStar))+(4.0*a*tanTheta2*Beta2/(aStar*aStar))+(2.0*a*tanTheta2*Beta2/(aStar));
    		Coeff[4] = (3.0*Beta4/(aStar*aStar))+(1.0*Beta4/aStar)-(6.0*a2)+(a4)+(1.0/(aStar))-(6.0*a2/(aStar))+(2.0*a2*Beta2)- (2.0*Beta2/(aStar*aStar))-(2.0*Beta2/(aStar))+(4.0*a2*Beta2/(aStar))+(3.0*tanTheta2*Beta4/(aStar*aStar))+(tanTheta2*Beta4/aStar)-(6.0*a2*tanTheta2)+(a4*tanTheta2)-(6.0*a2*tanTheta2/(aStar))+(2.0*a2*tanTheta2*Beta2)-(2.0*tanTheta2*Beta2/(aStar*aStar))-(tanTheta2*Beta2/(aStar))+(4.0*a2*tanTheta2*Beta2/(aStar));
    		Coeff[5] = (4.0*a)-(4.0*a3)+(4.0*a/(aStar))-(4.0*a*Beta2)-(8.0*a*Beta2/(aStar))+(4.0*a*tanTheta2)-(4.0*a3*tanTheta2)+(4.0*a*tanTheta2/(aStar))-(4.0*a*tanTheta2*Beta2)-(8.0*a*tanTheta2*Beta2/(aStar));
    		Coeff[6] = -1.0*Beta4-(3.0*Beta4/aStar)-(1.0)+(6.0*a2)-(1.0/(aStar))+(2.0*Beta2)-(2.0*a2*Beta2)+(4.0*Beta2/(aStar))-(tanTheta2*Beta4)+(3.0*tanTheta2*Beta4/aStar)-(tanTheta2)+(6.0*a2*tanTheta2)-(tanTheta2/(aStar))+(2.0*tanTheta2*Beta2)-(2.0*a2*tanTheta2*Beta2)+(4.0*tanTheta2*Beta2/(aStar));
    		Coeff[7] = (-4.0*a)+(4.0*a*Beta2)-(4.0*a*tanTheta2)+(4.0*a*tanTheta2*Beta2);
    		Coeff[8] = 1.0*Beta4+(1.0)-(2.0*Beta2)+(tanTheta2*Beta4)+(tanTheta2)-(2.0*tanTheta2*Beta2);

/* 		Solve for resonant roots and put into a complex array for further processing.  */
    		gsl_poly_complex_workspace *w = gsl_poly_complex_workspace_alloc( 9 );
    		gsl_err = gsl_poly_complex_solve( Coeff, 9, w, gsl_z );
    		gsl_poly_complex_workspace_free( w );
    		(*nRoots) = 0;
    		if (gsl_err != GSL_SUCCESS ) 
			{
			//printf("Problem with finding resonant roots!  Returning with nRoots=0\n");
			return((int) (-1));
			}	
    		for (i=0; i<8; i++) zz[i] = gsl_z[2*i] + gsl_z[2*i+1]*I;
		}

/* 	Gather applicable roots together into the resonantRoot array. */
    	Beta=sqrt(Beta2);
    	for ( i=0; i<8; i++ ) 
		{
        	if ( ( fabs(cimag(zz[i])) < 1e-10 ) && ( creal(zz[i]) > 0.0 ) )   //not an evanescent wave and index of refraction<infinity
			{
			x=creal(zz[i]);
			y=(x-a)/Beta;
			if ((y<0.0)||(y>=0.0)) 						//FIXME: should I be including waves in both directions or not?
				{
				nSquared=(y/x)*(y/x);
				R=1.0-((1.0/(aStar*x*x))*(x/(x-1.0)));
				L=1.0-((1.0/(aStar*x*x))*(x/(x+1.0)));
				S=0.5*(R+L);
				D=0.5*(R-L);
				P=1.0-(1.0/(aStar*x*x));
				if ((P/S > 0.0)&&(P/S < 1.0)) { changeInPolarization=asin(sqrt(P/S))*180.0/314159; } else { changeInPolarization=90.0;}  
				if ((changeInPolarization>0.01)&&(changeInPolarization<89.99)) {printf("Polarization label can change as a function of wave normal angle\n");}
				polarization = (nSquared-S)/D;
/*		
 *				Only keep the roots that correspond to R-mode polarization and are in the specified 'passband' in normalized frequency.  In Summers' formulation, 
 *				he explicitly only computes the right-hand or left-hand polarized waves because the two roots factor easily for parallel-propagating waves, but in the 
 *				general case we have to compute all the roots and throw out the ones that don't have the specified polarization.  This is complicated by the fact that 
 *				the polarization label may change as a function of wave normal angle (as indicated above).  Not sure what to do in this case; ideally one would track the 
 *				wave branch back to parallel-propagating and use that label...  FIXME: I'm currently assuming only use R-mode with electrons and L-mode with protons.
 */
				if (((s==1)&&(polarization < 0.0))||((s==-1)&&(polarization > 0.0))) 
					{ 
					if ((x>xlc)&&(x<xuc)) 
						{
						resonantRoots[(*nRoots)]= x; (*nRoots)+=1;
						}
					}
				}
			}
		}
return(0);
}

/* 
 *	The value of group velocity dw/dk for a cold electron plasma is obtained by recognizing that the dispersion relation is D(w,k,theta)=0 so that dD/dw=0 
 *	(since D is zero for all w).  Carrying through with partial derivatives wrt to w and k we find that dD/dw+dD/dk*dk/dw=0, and hence dk/dw=-(dD/dw)/(dD/dk).  
 *	Using the identities x=w/wce and y=ck/wce we find that -(dD/dx)/(dD/dy)=-(wce*dD/dw)/((wce/c)*(dD/dy))=-c*(dD/dw)/(dD/dk)=c*dk/dw, or dk/dw=1/c*[-(dD/dx)/(dD/dy)], 
 *	or dw/dk=c*[-(dD/dx)/(dD/dy)]^-1.  The values of dD/dx and dD/dy are found in mathematica by evaluating the following 3 expressions: 
 *
 *		A[x, y] := wce^2*x^2*(y^2/x^2-(1-(wpe^2/(wce^2*x^2))))*((wce^2*x^2*(x^2-1)-(wpe^2*x^2))*wce^2*y^2-((wce^2*x^2*(x+1)-(wpe^2*x))*(wce^2*x^2*(x-1)-(wpe^2*x))))
 *		B[x, y] := -1.0*(wce^2*x^2-wpe^2)*(y^2*wce^2*(x-1)-(wce^2*x^2*(x-1)-(wpe^2*x)))*(y^2*wce^2*(x+1)-(wce^2*x^2*(x+1)-(wpe^2*x)))
 *		Simplify[-1.0*D[A[x, y]*tanTheta2 - B[x, y], x]/D[A[x, y]*tanTheta2 - B[x, y], y]]	
 *
 *	where 
 *		aStar=wce^2/wpe^2, where wce=electron gyrofrequency and wpe=electron plasmafrequency
 *		x=normalized wave frequency (w/wce)
 *		y=normalized wave number (ck/wce)
 *		tanTheta2=tangent squared of the wave normal angle
 *	
 *	returns the group velocity in m/s.
 */
double electronColdPlasmaGroupVelocity(double aStar, double tanTheta2, double x, double y)
{
	double 	x2, x4, x6, y2, y4, dkdwnumerator, dkdwdenominator, dkdwTimesC, groupVelocity, aStar2, aStar3;

	x2 = x*x; x4=x2*x2; x6=x4*x2; y2 = y*y; y4=y2*y2;

/*	This version is not numerically stable
	dkdwnumerator = (2.0+(2.0*tanTheta2))*wpe6*x+(wce2*wpe4*x*((-12.0-(12.0*tanTheta2))*x2+((4.0+(4.0*tanTheta2))*y2)))+ 
   		(wce4*wpe2*x*(((18.0 + (18.0*tanTheta2))*x4)+(x2*(-4.0-(4.0*tanTheta2)-(16.0*y2)-(16.0*tanTheta2*y2)))+(y2*(4.0+(2.0*tanTheta2)+(2.0*y2)+(2.0*tanTheta2*y2)))))+ 
   		(wce6*x*((-8.0-(8.0*tanTheta2))*x6+((2.0+(2.0*tanTheta2))*y4)+(x2*y2*((-8.0-8.0*tanTheta2)-4.0*y2-(4.0*tanTheta2*y2)))+(x4*(6.0+(6.0*tanTheta2)+(12.0*y2)+(12.0*tanTheta2*y2)))));
	dkdwdenominator = (-4.0-(4.0*tanTheta2))*wce2*wpe4*x2*y+ 
   		(wce4*wpe2*y*((8.0+(8.0*tanTheta2))*x4+(4.0*y2)+(x2*(-4.0-(2.0*tanTheta2)-(4.0*y2)-(4.0*tanTheta2*y2)))))+ 
   		(wce6*x2*y*((-4.0-(4.0*tanTheta2))*x4+((-4.0-(4.0*tanTheta2))*y2)+(x2*(4.0+(4.0*tanTheta2)+(4.0*y2)+(4.0*tanTheta2*y2)))));
*/
/* 	In the following version, divide both numerator and denominator by wpe^6 in order to make it more numerically stable */
	aStar2=aStar*aStar; aStar3=aStar2*aStar;
	dkdwnumerator = (2.0+(2.0*tanTheta2))*x+(aStar*x*((-12.0-(12.0*tanTheta2))*x2+((4.0+(4.0*tanTheta2))*y2)))+ 
   		(aStar2*x*(((18.0 + (18.0*tanTheta2))*x4)+(x2*(-4.0-(4.0*tanTheta2)-(16.0*y2)-(16.0*tanTheta2*y2)))+(y2*(4.0+(2.0*tanTheta2)+(2.0*y2)+(2.0*tanTheta2*y2)))))+ 
   		(aStar3*x*((-8.0-(8.0*tanTheta2))*x6+((2.0+(2.0*tanTheta2))*y4)+(x2*y2*((-8.0-8.0*tanTheta2)-4.0*y2-(4.0*tanTheta2*y2)))+(x4*(6.0+(6.0*tanTheta2)+(12.0*y2)+(12.0*tanTheta2*y2)))));
	dkdwdenominator = (-4.0-(4.0*tanTheta2))*aStar*x2*y+ 
   		(aStar2*y*((8.0+(8.0*tanTheta2))*x4+(4.0*y2)+(x2*(-4.0-(2.0*tanTheta2)-(4.0*y2)-(4.0*tanTheta2*y2)))))+ 
   		(aStar3*x2*y*((-4.0-(4.0*tanTheta2))*x4+((-4.0-(4.0*tanTheta2))*y2)+(x2*(4.0+(4.0*tanTheta2)+(4.0*y2)+(4.0*tanTheta2*y2)))));
	dkdwTimesC = dkdwnumerator/dkdwdenominator;
	groupVelocity=(1.0/dkdwTimesC)*LGM_c;
//printf("for wce2=%g wpe2=%g tanTheta2=%g x=%g y=%g group velocity is %g m/s\n", wce2, wpe2, tanTheta2, x, y, groupVelocity);
	return(groupVelocity);
}

/*	
 *	Set up the function that is integrated over tan(theta) to produce the normalizer N(w) and evaluate it for a range of aStar's; the normalizer is 
 *	only valid for the local plasma environment, so it must be re-evaluated at each point along a field line, e.g., althoug the local environment
 *	is specified solely by aStar. 
 *
 *	Outputs		Nw[]	contains the normalizers for nNw values of normalized frequencies (w/wce) in the range [wxl,wxh] for a set of 
 *				plasma parameters of size numberOfPlasmaParameters that span the range [aStarMin,aStarMax] and are sampled geometrically.
 *				Must be pre-allocated to be of size nNw*numberOfPlasmaParameters.
 */
int computeNormalizerForWavePowerSpectrumFunctionForRangeOfPlasmaParameters(double wxl, double wxh, double xmin, double xmax, int numberOfWaveNormalAngleDistributions, double *xmArray, double *dxArray, double *weightsOnWaveNormalAngleDistributions, int nNw, int numberOfPlasmaParameters, double aStarMin, double aStarMax, double *Nw) 
{
	double 	result, error, aStar, gsl_err;
    	struct	normalizerForWavePowerSpectrumFunctionParams p;
	int	i, iPlasma;

	p.xmin=xmin; p.xmax=xmax; p.numberOfWaveNormalAngleDistributions=numberOfWaveNormalAngleDistributions; 
	p.xmArray=(double *) xmArray; p.dxArray=(double *) dxArray; p.weightsOnWaveNormalAngleDistributions=(double *) weightsOnWaveNormalAngleDistributions;
	gsl_integration_workspace * work=gsl_integration_workspace_alloc(nNw);
	gsl_function F;
	gsl_set_error_handler_off(); // so that we don't core dump if the integration routine can't converge on an answer
	F.function = &normalizerForWavePowerSpectrumFunction;
	F.params=&p;
	for (iPlasma=0; iPlasma<numberOfPlasmaParameters; iPlasma++)
		{
		aStar = aStarMin*pow(aStarMax/aStarMin, (iPlasma+0.5)/((double) numberOfPlasmaParameters));
		p.aStar=aStar; 
		for (i=0; i<nNw; i++)
			{
			p.x = wxl+((wxh-wxl)*(i+0.5)/((double) nNw));
			result = 0.0;
			gsl_err = (int) gsl_integration_qags(&F, xmin, xmax, 0, 1e-4, nNw, work, &result, &error);  		
    			if (gsl_err != GSL_SUCCESS ) 
				{
				//printf("gsl_integration_qags failed to produce N(w) for aStar=%g, w/wce=%g, i=%d of %d with errno=%d, result=%g, error=%g; setting to 1.0\n", aStar, p.x, i, nNw, gsl_err, result, error);
				Nw[i+(iPlasma*nNw)] = 1.0;
				}	
			Nw[i+(iPlasma*nNw)] = result;
			}
		}
	gsl_integration_workspace_free(work);
	return(0);
}

/*
 *	Compute N(w), a normalization that ensures that the energy per unit
 *	frequency is given by B^2(w) for any wave normal angle distribution
 *	consistent with the dispersion relation.  I believe that the idea here is
 *	that the wave power is partitioned into wave normal angles as well as
 *	frequency; however, for a specified wave normal angle and frequency, there
 *	are only two wave number magnitudes, |k|, that satisfy the dispersion
 *	relation and presumably only one that is the right 'mode'.  This is the
 *	wave number magnitude we use in the equations from Glauert and Horne.
 *	Thus, we need to perform an integral of the function
 *	h(X)=1/2pi^2[g(X)k^2/(1+X^2)^3/2][dk/dw]X over X.  First, we define the
 *	function h(X) and then we integrate over it, using quadpack since this is
 *	what Mike uses to do the Summers' integrals.
 *
 *	Modified by Greg Cunningham on 6-1-2013 to account for the fact that, for
 *	each frequency w (which is contained in the struct element p.w), there is
 *	potentially a resonance cone angle, theta, in the interval 0<=theta<=pi/2
 *	that limits the range of theta over which resonant roots should be accepted
 *	to [0,theta].  If the input argument, tanTheta, is not within this range
 *	then the routine should return zero since we will not be accepting resonant
 *	roots for values of theta outside this range.
 *
 *	Modified by Greg Cunningham on 6-6-2014 to use only normalized frequencies
 *	and aStar=wce^2/wpe^2 instead of wce, wpe, w, etc. Unfortunately, this
 *	means that the result that is returned is 1/wce^2 times the old result and
 *	so callers of this routine have to multiply the result by wce^2.  This
 *	change was made so that the array Nw[i] that is used to normalize the
 *	spectrum at each frequency can be precomputed for a range of aStar's and
 *	then the column that corresponds to the correct aStar can be selected.
 */
double normalizerForWavePowerSpectrumFunction(double tanTheta, void *p) {
	double	aStar,x2,tanTheta2,sinTheta2,sinTheta4,cosTheta2;
	double 	P,R,L,S,D,P2,D2,A,B,F,F2,n2_1,n2_2,n2,polarization1,polarization2,k,k2,c2,x,y,h,thisg,g,dwdk,dkdw;
	double	xmean, xmin, xmax, dx, weight, test;
	int	j;

	struct 	normalizerForWavePowerSpectrumFunctionParams *params;
	params =(struct normalizerForWavePowerSpectrumFunctionParams *) p;
	aStar = (double) params->aStar;
	x = (double) params->x; x2=x*x;
	xmin = (double) params->xmin;
	xmax = (double) params->xmax;
	tanTheta2=tanTheta*tanTheta;
	sinTheta2 = tanTheta2/(1.0+tanTheta2); sinTheta4=sinTheta2*sinTheta2;
	cosTheta2 = 1.0/(1.0+tanTheta2);
	P = 1.0-(1.0/(aStar*x2)); P2=P*P;
	L = 1.0-((1.0/(aStar*x2))*(x/(x+1.0)));  
	R = 1.0-((1.0/(aStar*x2))*(x/(x-1.0)));  
	S = 0.5*(R+L);
	D = 0.5*(R-L); D2=D*D;
	A = S*sinTheta2+(P*cosTheta2);
	B = R*L*sinTheta2+(P*S*(1.0+cosTheta2));
	F2 = (R*L-(P*S))*(R*L-(P*S))*sinTheta4+(4.0*P2*D2*cosTheta2);
	test = ((-1.0*P/S*0.99)<(tanTheta2))&&(-1.0*P/S > 0.0);
	if (test) 
		{
		//FIXME: may want to print this warning out printf("In normalizerForWavePowerSpectrumFunction, resonance cone angle is less than the input theta: tanTheta=%g x=%g, aStar=%g, P=%g, L=%g, R=%g, S=%g, D=%g, A=%g, B=%g, F2=%g\n", tanTheta, x, aStar, P, L, R, S, D, A, B, F2);
		return(0.0);
		}
	if (F2 < 0.0) 
		{
		//FIXME: may want to print this warning out printf("warning: F2<0 in normalizerForWavePowerSpectrumFunction, function evaluates to 0\n"); 
		return(0.0);
		}
	else
		{
		F = sqrt(F2);
		n2_1 = (B+F)/(2.0*A);
		n2_2 = (B-F)/(2.0*A);
		}
/*	Find the solution for which the polarization is R-mode */
	polarization1 = (n2_1-S)/D;
	polarization2 = (n2_2-S)/D;
	if (polarization1 > 0.0)   //FIXME: this is fixed to use R-mode and electrons only; need to generalize
		{
		if (polarization2 > 0.0)
			{
			//FIXME: may want to print this out printf("warning: both solutions to dispersion relation in normalizerForWavePowerSpectrumFunction are R-mode: pick first one\n");
			}
		n2 = n2_1;
		}
	else 
		{
		if (polarization2 < 0.0)
			{
			//FIXME: may want to print this out printf("warning: both solutions to dispersion relation in normalizerForWavePowerSpectrumFunction are not R-mode: function evaluates to 0\n");
			if (test) { printf("\n");}
			return(0.0);
			}
		n2 = n2_2;
		}
/*
 *	Compute k2 from w2 and n2=(ck)^2/w^2
 */
	c2 = (LGM_c*LGM_c);  // squared speed of light in m^2/s^2
	k2 = n2*x2/c2; k=sqrt(k2);  //k is actually k/wce since x=w/wce
/*	Compute the group velocity */
	y=LGM_c*k;
	dwdk = (double) electronColdPlasmaGroupVelocity(aStar, tanTheta2, x, y);
	dkdw = (1.0/dwdk);
/*	Evaluate the function at this w */
	if ((tanTheta>=xmin)&&(tanTheta<=xmax))
		{
		g = 0.0;
		for (j=0; j<params->numberOfWaveNormalAngleDistributions; j++)
			{
			xmean = (double) params->xmArray[j];
			dx = (double) params->dxArray[j];
			weight = (double) params->weightsOnWaveNormalAngleDistributions[j];
			thisg = weight*exp(-1.0*(tanTheta-xmean)*(tanTheta-xmean)/(dx*dx));
			g += thisg;
			}
 		h=(1.0/(2.0*M_PI*M_PI))*(g*k2/pow(1.0+tanTheta2,1.5))*dkdw*tanTheta;		// h will have to be multiplied outside this subroutine by wce^2 to get correct answer
		}
	else {h=0.0;}  // tanTheta is outside range of valid wave normal angles, so integrand is zero
//printf(" for tanTheta=%g, x=%g, y=%g, dwdk=%g, h=%g\n", tanTheta, x, y, dwdk, h);
	if (test) { printf(" n2=%g, dkdw=%g, h=%g\n", n2,dkdw,h);}
	return(h);
}

/*
 *	The following function evaluates the normalizing constant |/phi_{nk}|^2 in
 *	equation A13 of Glauert and Horne 2005.  The Bessel functions are evaluated
 *	using the argument \frac{k_{/perp}p_{/perp}}{m_{/sigma}/Omega_{/sigma}}
 *	where /sigma denotes the particle type (electron or proton) and the
 *	gyrofreq is signed.  Because p is specified in units of MeV/c, and
 *	y=ck/wce, the product p*y/m_e has units of MeV/kg; need to convert MeV to
 *	kg*m/sec so that the units work out properly.  Since 1 MeV=1e6 eV=1.6e-13
 *	Joules and c=3e8 m/s, we have p(MeV/c)*1.6e-13/3e8 has units of kg*m/s.
 *	Meanwhile, y=ck/wce must be divided by c to get k/wce, which has units of
 *	1/(m/s) so that y*p*sin(theta)/m_e, where m_e is the electron mass in kg,
 *	is unitless.
 *
 *	s			particle type (-1=electron, +1=proton)
 *	p			momentum (in units of MeV/c)
 *	alpha			local pitch angle (in radians)
 *	nCyclotronNumber	the 'n' in the resonance condition
 *	x			normalized wave frequency w/wce, where wce is the electron gyrofrequency
 *	y			normalized wave number ck/wce, where c is the speed of light
 *	P, R, L, S		dispersion relation values as labeled in Stix 1962 page 11
 *	theta			wave normal angle
 *
 */
double besselFunctionNormalizer(int s, double p, double alpha, int nCyclotronNumber, double x, double y, double P, double R, double L, double S, double theta)
{
	double 	n2,arg,sinTheta,sinTheta2,cosTheta,phink2;
	double	J[3];
	int	i,thisN;

	n2 = (y/x)*(y/x);		//square of the index of refraction, n=ck/w
	sinTheta=sin(theta); sinTheta2=sinTheta*sinTheta; 
	cosTheta=cos(theta);
	arg = s*(y/LGM_c)*sinTheta*p*sin(alpha)*(1e6*Joules_Per_eV/LGM_c)/LGM_ELECTRON_MASS; //the argument
	thisN = nCyclotronNumber-1;
	if (thisN < 0) 
		{
		thisN = -1*thisN;
		i=gsl_sf_bessel_Jn_array(thisN, thisN, (double) (-1.0*arg),&J[0]);
		}
	else
		{
		i=gsl_sf_bessel_Jn_array(thisN, thisN, arg, &J[0]);
		}
	thisN = nCyclotronNumber;
	if (thisN < 0) 
		{
		thisN = -1*thisN;
		i=gsl_sf_bessel_Jn_array(thisN, thisN, (double) (-1.0*arg),&J[1]);
		}
	else
		{
		i=gsl_sf_bessel_Jn_array(thisN, thisN, arg, &J[1]);
		}
	thisN = nCyclotronNumber+1;
	if (thisN < 0) 
		{
		thisN = -1*thisN;
		i=gsl_sf_bessel_Jn_array(thisN, thisN, (double) (-1.0*arg),&J[2]);
		}
	else
		{
		i=gsl_sf_bessel_Jn_array(thisN, thisN, arg, &J[2]);
		}
	phink2=pow((((n2-L)/(n2-S))*J[2]+(((n2-R)/(n2-S))*J[0]))*((n2*sinTheta2-P)/(2.0*n2))+(sinTheta*cosTheta*J[1]/tan(alpha)), 2.0)/(pow((R-L)/(2.0*(n2-S)),2.0)*pow((P-(n2*sinTheta2))/n2, 2.0)+pow(P*cosTheta/n2, 2.0));
	return(phink2);
}

/*
 *	The following routine finds intervals in theta for which it is possible
 *	that there is a resonant root at the intersection between the resonance
 *	condition and the dispersion relation.  Inputs are 
 *		wn	doppler shift, s*n*wce/gamma, normalized by the gyrofrequency, wce, to give s*n/gamma
 *		wlc	lower cutoff frequency normalized by the gyrofrequency
 *		wuc	uppper cutoff frequency normalized by the gyrofrequency
 *		aStar	the squared ratio of the electron gyrofrequency to the plasma frequency, wpe
 *		mode	=1 for R-mode, =-1 for L-mode
 *		vpar	the parallel velocity divided by the speed of light
 *
 *	The code determines the intervals in theta for which both of two inequality
 *	conditions is met (the minimum of one monotonic function is less than the
 *	maximum of the other, and vice versa).  A list of intervals in theta that
 *	satisfy each of the two conditions is returned by the code and then the
 *	intersection of these intervals is computed and returned. The process may
 *	result in 0, 1 or as many as 4 intervals being returned.
 *
 */
int computeWaveNormalAngleIntervalsWhereResonantRootsMayExist(double wn, double wlc, double wuc, struct inequalityParameters *inequalityParamsLC, struct inequalityParameters *inequalityParamsUC, int mode, int *nIntervals, double *thetaStart, double *thetaEnd) {
	double	x, dw, thisThetaStart[2], thisThetaEnd[2];
	int	thisN, i;

	*nIntervals = 0;
	thisN = 0;
	/* 	x is either wuc or wlc depending on the inequality */
	if (wn < wlc) 
		{
		/* condition (i) from Jay Albert's 2007 paper: if either condition is satisfied, there can be a root */
		dw = wlc-wn;
		i = findSingleIntervalInTheta(inequalityParamsUC, dw, (int) mode, nIntervals, thetaStart, thetaEnd);

		/* condition (ii) from Jay Albert's 2007 paper */
		dw = wuc-wn;
		i = findSingleIntervalInTheta(inequalityParamsLC, dw, (int) -1*mode, &thisN, thisThetaStart, thisThetaEnd);
		}
	if (wn > wuc)
		{
		/* condition (i) from Jay Albert's 2007 paper: if either condition is satisfied, there can be a root */
		dw = wuc-wn;
		i = findSingleIntervalInTheta(inequalityParamsUC, dw, (int) mode, nIntervals, thetaStart, thetaEnd);

		/* condition (ii) from Jay Albert's 2007 paper */
		dw = wlc-wn;
		i = findSingleIntervalInTheta(inequalityParamsLC, dw, (int) -1*mode, &thisN, thisThetaStart, thisThetaEnd);
		}
	if ((wn <= wuc) && (wn >= wlc))
		{
		/* condition (i) from Jay Albert's 2007 paper automatically satisifed */
		*nIntervals = 1;
		thetaStart[0] = 0.0;
		thetaEnd[0] = M_PI/2.0;
		/* condition (ii) from Jay Albert's 2007 paper */
		if (wn-wlc > (wuc-wn)) {dw = wlc-wn;} else {dw=wuc-wn;}
		i = findSingleIntervalInTheta(inequalityParamsLC, dw, (int) -1*mode, &thisN, thisThetaStart, thisThetaEnd);
		}
/*
printf("In computeWaveNormalAngleIntervalsWhereResonantRootsMayExist: \n");
printf("	from condition (i), have nIntervals=%d\n", *nIntervals);
for (i=0; i<*nIntervals; i++)
	{
	printf("	interval %d is [%g, %g]\n", i, thetaStart[i], thetaEnd[i]);
	}
printf("	from condition (ii), have nIntervals=%d\n", thisN);
for (i=0; i<thisN; i++)
	{
	printf("	interval %d is [%g, %g]\n", i, thisThetaStart[i], thisThetaEnd[i]);
	}
*/
	/* combine intervals */
	if ((*nIntervals) == 0) {return(0);}				/* means that condition (i) is not satisified for any theta */
	if (thisN == 0) {(*nIntervals) = 0; return(0);}			/* means that condition (ii) is not satisified for any theta */
	if ((*nIntervals) == 1)
		{
		if (thisN == 1)
			{
			i = mergeOneIntervalWithOneInterval(thisThetaStart, thisThetaEnd, thetaStart, thetaEnd);
			*nIntervals = i;
			}
		if (thisN == 2)
			{
			i = mergeTwoIntervalsWithOneInterval(thisThetaStart, thisThetaEnd, thetaStart, thetaEnd);
			*nIntervals = i;
			}
		}
	if ((*nIntervals) == 2)
		{
		if (thisN == 1)
			{
			i = mergeOneIntervalWithTwoIntervals(thisThetaStart, thisThetaEnd, thetaStart, thetaEnd);
			*nIntervals = i;
			}
		if (thisN == 2)
			{
			i = mergeTwoIntervalsWithTwoIntervals(thisThetaStart, thisThetaEnd, thetaStart, thetaEnd);
			*nIntervals = i;
			}
		
		}
/*
printf("	number of combined intervals is %d\n", *nIntervals); 
for (i=0; i<*nIntervals; i++) {printf("		combined interval %d is [%g,%g]\n", i, thetaStart[i], thetaEnd[i]); }
*/
	return(0);
}

int mergeOneIntervalWithOneInterval(double *thisThetaStart, double *thisThetaEnd, double *thetaStart, double *thetaEnd) {
	if (thisThetaStart[0] < thetaStart[0])
		{
		if (thisThetaEnd[0] < thetaStart[0]) {return(0);}	// the two intervals do not overlap; return 0 intervals
		else 							// the intervals overlap; only question is what is endpoint of intersection
			{
			if (thisThetaEnd[0] < thetaEnd[0]) { thetaEnd[0] = thisThetaEnd[0]; } 
			return(1);
			}
		}
	else 
		{
		if (thisThetaStart[0] > thetaEnd[0]) {return(0);}	// the two intervals do not overlap; return 0 intervals
		else 							// the intervals overlap; only question is what is endpoint
			{
			thetaStart[0] = thisThetaStart[0];
			if (thisThetaEnd[0] < thetaEnd[0]) { thetaEnd[0] = thisThetaEnd[0]; } 
			return(1);
			}
		}
	return(0);
}

/* 
 *	In the following code it is assumed we have two collections of intervals.  Both collections have two intervals and we are interested in the
 *	intersection of the two intervals.  The first interval of each collection starts at 0 and the second interval of each collection ends at pi/2.  
 *	Return two intervals that are the intersection of the contributing intervals.
 */

int mergeTwoIntervalsWithTwoIntervals(double *thisThetaStart, double *thisThetaEnd, double *thetaStart, double *thetaEnd) {
	if (thisThetaEnd[0] < thetaEnd[0]) 			//left-hand interval of one collection is subsumed by the other
		{
		if (thisThetaStart[1] < thetaEnd[0])		//right-hand interval of this same collection intersects left-hand interval of the other; have 3 intervals
			{
			thetaStart[2] = thetaStart[1];
			thetaEnd[2] = thetaEnd[1];
			thetaStart[1] = thisThetaStart[1];
			thetaEnd[1] = thetaEnd[0];
			thetaEnd[0] = thisThetaEnd[0];
			return(3);
			}
		else
			{
			thetaEnd[0] = thisThetaEnd[0];
			if (thisThetaStart[1] > thetaStart[1])	// right-hand interval from thisTheta* is subsumed by right-hand interval from theta*, so redefine left edge
				{
				thetaStart[1] = thisThetaStart[1];
				return(2);
				}
			else { return(2); }			
			}
		}
	else
		{
		if (thisThetaEnd[0] > thetaStart[1])		// have 3 intervals
			{
			thetaStart[2] = thisThetaStart[1];
			thetaEnd[2] = thetaEnd[1];
			thetaEnd[1] = thisThetaEnd[0];
			return(3);
			}
		else
			{
			if (thisThetaStart[1] > thetaStart[1])
				{
				thetaStart[1] = thisThetaStart[1];
				return(2);
				}
			else {return(2);}			// the two intervals described by thisTheta* subsume the two intervals described by theta*, so no change needed
			}
		}
	return(0);
}

/*
 *	Merge an interval of the sort (a,b) with two intervals of the sort [0,c)+(d,pi/2]
 */

int mergeOneIntervalWithTwoIntervals(double *thisThetaStart, double *thisThetaEnd, double *thetaStart, double *thetaEnd) {
	if (thisThetaStart[0] < thetaEnd[0])
		{
		if (thisThetaEnd[0] < thetaEnd[0]) 		// interval (a,b) subsumed by left-hand interval described by theta*, return (a,b) interval only
			{
			thetaStart[0] = thisThetaStart[0];
			thetaEnd[0] = thisThetaEnd[0];
			return(1);
			}
		else
			{
			if (thisThetaEnd[0] > thetaStart[1]) 	// two intervals 
				{
				thetaStart[0] = thisThetaStart[0];
				thetaEnd[1] = thisThetaEnd[0];
				return(2);
				}
			else					// one interval
				{
				thetaStart[0] = thisThetaStart[0];
				return(1);
				}
			}
		}
	else
		{
		if (thisThetaStart[0] > thetaStart[1])		// interval (a,b) subsumed by right-hand interval described by theta*, return (a,b) only
			{
			thetaStart[0] = thisThetaStart[0];
			thetaEnd[0] = thisThetaEnd[0];
			return(1);
			}
		else
			{
			if (thisThetaEnd[0] < thetaStart[1]) {return(0);}	// no intersection, return 0
			else 
				{
				thetaStart[0] = thetaStart[1];
				thetaEnd[0] = thisThetaEnd[0];
				return(1);
				}	
			}
		}
	return(0);
}

int mergeTwoIntervalsWithOneInterval(double *thisThetaStart, double *thisThetaEnd, double *thetaStart, double *thetaEnd) {
	double	newThetaStart[3], newThetaEnd[3];
	int	i, j;

	newThetaStart[0] = thisThetaStart[0];
	newThetaStart[1] = thisThetaStart[1];
	newThetaEnd[0] = thisThetaEnd[0];
	newThetaEnd[1] = thisThetaEnd[1];
	i= mergeOneIntervalWithTwoIntervals(thetaStart, thetaEnd, newThetaStart, newThetaEnd);
	for (j=0; j<i; j++)
		{
		thetaStart[j] = newThetaStart[j];
		thetaEnd[j] = newThetaEnd[j];
		}
	return(i);
}

/*
 * 	Following routine computes the parameters that are re-usable for all harmonic numbers in order to minimize the cost of doing the inequality test.
 */

int computeInequalityParameters(double x, double aStar, double vpar, struct inequalityParameters *params) {
	double	R, L, S, P;

	R = 1.0-(x/(aStar*x*x*(x-1)));
	L = 1.0-(x/(aStar*x*x*(x+1)));
	S = (R+L)/2.0;
	P = 1.0-(1.0/(aStar*x*x));
	params->D = S;
	params->E1 = -1.0*(R*L+(P*S))*pow(vpar, 2.0)*pow(x, 2.0);
	params->E2 = P-S;
	params->F1 = -1.0*((P*S-(R*L))*pow(vpar, 2.0)*pow(x, 2.0));
	params->F2 = P*R*L*pow(vpar, 4.0)*pow(x, 4.0);

	return(0);
}

int findSingleIntervalInTheta(struct inequalityParameters *params, double dw, int inequalitySign, int *nIntervals, double *thetaStart, double *thetaEnd) {
	double	dw2,dw4,D,E,F,arg,root1,root2;

	dw2 = dw*dw;
	dw4 = dw2*dw2;
	D = (params->D)*dw4;
	E = (params->E1)*dw2+((params->E2)*dw4);
	F = (params->F1)*dw2+(params->F2);
//printf(" D=%g, E=%g, F=%g\n", D, E, F);
	arg = E*E-(4.0*D*F);
	if (arg < 0.0) 	
		{ 					/* then entire interval is either satisfied or not */
		if ((inequalitySign*(D+E+F)) > 0.0)  	/* entire interval does satisify the inequality, so can return whole interval, otherwise go on to condition (ii) */
			{
			*nIntervals = (int) 1;
			thetaStart[0] = (double) 0.0;
			thetaEnd[0] = (double) (M_PI/2.0);
			return(0); 
			}
		else
			{
			*nIntervals = (int) 0;		/* Added on 5-29-2014; if the inequality is not satisfied at theta=0, then there is no interval that satisfies the inequality */
			}
		}
	else 			
		{			/* there are roots of the quadratic D+Ex+Fx^2, where x=cos^2(theta)*/
		if (F > 0.0)
			{
			root1 = (-1.0*E+sqrt(arg))/(2.0*F);  		/* because taking the inverse cos of the sqrt of the root will reverse the order of the roots */
			root2 = (-1.0*E-sqrt(arg))/(2.0*F);
			}
		else
			{
			root2 = (-1.0*E+sqrt(arg))/(2.0*F);  		/* because taking the inverse cos of the sqrt of the root will reverse the order of the roots */
			root1 = (-1.0*E-sqrt(arg))/(2.0*F);
			}
		if (root1 < 0.0) {root1 = 0.0;}
		if (root1 > 1.0) {root1 = 1.0;}
		if (root2 < 0.0) {root2 = 0.0;}
		if (root2 > 1.0) {root2 = 1.0;}
//printf("		root1=%g root2=%g\n", root1, root2);
		if ((root1-root2) > 1e-8) 
			{
			if ((inequalitySign*F) < 0.0) 	/* 	Either the quadratic is concave up and the inequality is < or the quadratic is concave down and the inequality is >. 
								In either case, the interval in cos^2(theta) that can have roots is [root1, root2] */
				{
				*nIntervals = (int) 1;
				thetaStart[0] = (double) acos(sqrt(root1)); 
				thetaEnd[0] = (double) acos(sqrt(root2));
				}
			else 		/* 	Either the quadratic is concave down and the inequality is < or the quadratic is concave up and the inequality is >.
						In either case, the intervals in cos^2(theta) that can have roots are [0,root1] and [root2,1]; get rid of degenerate intervals */
				{
				if (root1 == 1)
					{
					if (root2 ==0)
						{
						*nIntervals = 0;
						}
					else 
						{
						*nIntervals = 1;
						thetaStart[0] = (double) acos(sqrt(root2));
						thetaEnd[0] = (double) (M_PI/2.0);
						}
					}
				else
					{
					thetaStart[0] = (double) 0.0;
					thetaEnd[0] = (double) acos(sqrt(root1));
					if (root2 ==0)
						{
						*nIntervals = (int) 1;
						}
					else 
						{
						*nIntervals = (int) 2;
						thetaStart[1] = (double) acos(sqrt(root2));
						thetaEnd[1] = (double) (M_PI/2.0);
						}
					}
				}
			}
		else			/* this is the case when the quadratic effectively touches zero and so the whole interval can have roots or not */
			{
			if ((inequalitySign*(D+E+F)) > 0.0)  	/* entire interval satisifies the inequality, so return whole interval */
				{
				*nIntervals = (int) 1;
				thetaStart[0] = (double) 0.0;
				thetaEnd[0] = (double) (M_PI/2.0);
				return(0); 
				}
			}
		}
	return(0);
}


/*
 *  Following routine is just a very simple Riemann sum to evaluate the integral of a function, f, over the range [xleft,xright].  
 *  The function is parameterized by args.  The result of the integration is put into result.
 */

int Lgm_SimpleRiemannSum( double (*f)(double, _qpInfo *), _qpInfo *args, double xleft, double xright, double *result, int VerbosityLevel ) {

	int	nIntervals, nSubIntervals, i;
	double	x,sum,dx, epsabs=0.0, epsrel=0.1,abserr;
	size_t 	neval;

/*	Approach #1: use a simple Riemann sum
	nIntervals = 100;

	dx = (xright-xleft)/((double) nIntervals);
	sum = 0.0;
	for (i=0; i<nIntervals; i++)
		{
		x = xleft+((i+0.5)*dx);
		sum += dx*f(x,args);
		}
	*result = sum;
*/

/* 	Approach #2: use GSL's implementation of the simple quadpack routine qng (non-adaptive Gauss-Kronrod quadrature) to do the integral 
	gsl_function F;
	F.function = f;
	F.params = args;
 	gsl_integration_qng (&F, xleft, xright, epsabs, epsrel, result, &abserr, &neval);
	printf("Number of evaluations needed by non-adaptive Gauss-Kronrod quadrature to achieve %g relative accuracy is %d\n", epsrel, neval);
*/

/* 	Approach #3: use GSL's implementation of more complicated quadpack routine qag (adaptive quadrature with no singularities) to do the integral */
	nSubIntervals = 20;  
	gsl_integration_workspace * w =gsl_integration_workspace_alloc( nSubIntervals );
	gsl_function F;
	F.function = f;
	F.params = args;
	*result=0.0;
	gsl_integration_qag(&F, xleft, xright, epsabs, epsrel, (size_t) nSubIntervals, GSL_INTEG_GAUSS21, w, result, &abserr);
//	printf("Number of evaluations needed by adaptive Gauss-Kronrod 21-point quadrature to achieve %g relative accuracy is %d, abserr is %g and integral is %g\n", epsrel, w->size, abserr, *result);
	gsl_integration_workspace_free(w);

	return(w->size);
}
