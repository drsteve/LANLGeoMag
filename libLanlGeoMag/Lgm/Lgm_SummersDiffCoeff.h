#ifndef LGM_SUMMERS_DIFF_COEFF_H
#define LGM_SUMMERS_DIFF_COEFF_H
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_poly.h>
#include "Lgm/Lgm_CTrans.h"
#include "Lgm/Lgm_QuadPack.h"
#include "Lgm/Lgm_PolyRoots.h"
#include "Lgm/Lgm_Constants.h"


#ifndef M_CDIP
#define M_CDIP  30040.680710                // dipole moment (in nT) for center dipole field (computed for Aug 12, 2004  12.34567 UTC) Need to generalize this.
#endif


#define LGM_FRWD         1
#define LGM_FRWD_BKWD    0
#define LGM_BKWD        -1


#define LGM_R_MODE_WAVE     0
#define LGM_L_MODE_WAVE     1

#define LGM_ELECTRONS       0
#define LGM_PROTONS         1

#define LGM_SUMMERS_2005    		2005
#define LGM_SUMMERS_2007    		2007
#define LGM_GLAUERT_AND_HORNE_HIGH_FREQ    1 

// Derivative schemes
#ifndef LGM_DERIV_TWO_POINT
#define LGM_DERIV_TWO_POINT     0
#endif
#ifndef LGM_DERIV_FOUR_POINT
#define LGM_DERIV_FOUR_POINT    1
#endif
#ifndef LGM_DERIV_SIX_POINT
#define LGM_DERIV_SIX_POINT     2
#endif

typedef struct Lgm_SummersInfo {

    int     Version;        //!< Vewrsion of Summers model to use. Can be LGM_SUMMERS_2005 or LGM_SUMMERS_2007.
    double  Alpha0;         //!< Equatorial pitch angle, \f$ \alpha_\circ \f$, [radians].

    double  SinAlpha0;      //!< \f$ \sin\alpha_\circ \f$
    double  CosAlpha0;      //!< \f$ \cos\alpha_\circ \f$

    double  SinAlpha02;     //!< \f$ \sin^2\alpha_\circ \f$
    double  CosAlpha02;     //!< \f$ \cos^2\alpha_\circ \f$

    double  TanAlpha0;      //!< \f$ \tan\alpha_\circ \f$
    double  TanAlpha02;     //!< \f$ \tan^2\alpha_\circ \f$

    double  E;              //!< Normalized dimensionless energy, \f$ E_k/E_\sigma\circ \f$, where \f$ E_\sigma\circ \f$ is rest energy of species.
    double  L;              //!< Dipole L-shell value
    double  n1;             //!< Ratio of Hydrogen number density to Electron number density (note that n1+n2+n3=1).
    double  n2;             //!< Ratio of Helium number density to Electron number density (note that n1+n2+n3=1).
    double  n3;             //!< Ratio of Oxygen number density to Electron number density (note that n1+n2+n3=1).
    double  aStarEq;        //!< Summer's \f$ \alpha^* \f$ value which is \f$ \Omega_e/\omega^2_{pe} \f$.
//    double  dB;           //!< Value of wave amplitude [nT].
    void   *BwFuncData;     //!< Pointer to data that may be needed by BwFunc()
    double (*BwFunc)();     //!< Function to return Bw as a function of latitude.
    double  Omega_eEq;      //!< Equatorial gyrofrequency of electrons [Hz].
    double  Omega_SigEq;    //!< Equatorial gyrofrequency of species [Hz].
    double  w1;             //!< Lower frequency cutoff [Hz].
    double  w2;             //!< Upper frequency cutoff [Hz].
    double  wm;             //!< Frequancy of maximum wave power.
    double  dw;             //!< Measure of width of guassian wavre frequency distribution.
    double  x1;		    //!< Lower cutoff in tangent of the wave normal angle
    double  x2;		    //!< Upper cutoff in tangent of the wave normal angle
    int     numberOfWaveNormalAngleDistributions; //!< Number of gaussians in the tangent of wave normal angle that are weighted and summed together to produce the total distribution
    double  *xm;	    //!< Array of mean tangent of the wave normal angle for each component of the total distribution
    double  *dx;	    //!< Array of widths of the tangent wave normal angle for each component of the total distribution
    double  *weightsOnWaveNormalAngleDistributions; //!< Weight that is applied to each component of the total distribution; must add to 1.0
    double  MaxWaveLat;     //!< Latitudinal cuttoff for waves. I.e. assume no waves at lats greater than +MaxWaveLat or less than -MaxWaveLat.
    double  Sig;            //!< Additional paramter to define 'semi-bandwidth' which is equal to \f$ \sigma d\omega \f$.
    int     s;              //!< Defines wave mode (s=1 for R-mode, s=-1 for L-mode).
    double  Lambda;         //!< Defines sepcies (\f$ \lambda=-1 \f$ for electrons \f$ \lambda = \epsilon \f$ for protons (\f$\epsilon =m_e/m_p \f$).
    double  Rho;            //!< Rho = sqrt(PI)/2.0( erf( (wm-w1)/dw ) + erf( (w2-wm)/dw )
    int     Directions;     //!< Wave directionsm to include. Can be LGM_FRWD, LGM_BKWD, LGM_FRWD_BKWD
    double  aStarMin;	    //!< minimum value of aStar=wce^2/wpe^2 specifying the grid over which the normalizer needed for Glauert and Horne is computed.
    double  aStarMax;	    //!< maximum value of aStar=wce^2/wpe^2 specifying the grid over which the normalizer needed for Glauert and Horne is computed.
    int	    nNw; 	    //!< number of grid points in normalized frequency=w/wce spanning [0,1] on which the normalizer needed for Glauert and Horne is computed.
    int	    nPlasmaParameters; //!< number of logarithmic grid points in aStar on which the normalizer for Glauert and Horne is computed
    double  *Nw;	    //!< pointer to the gridded evaluations of the normalizer needed for Glauert and Horne.

} Lgm_SummersInfo;


int Lgm_SummersDxxBounceAvg( int Version, double Alpha0,  double Ek,  double L,  void *BwFuncData, double (*BwFunc)( double, void * ), double n1, double n2, double n3, double aStarEq,  int Directions, double w1, double w2, double wm, double dw, int WaveMode, int Species, double MaxWaveLat, double *Daa_ba,  double *Dap_ba,  double *Dpp_ba);
//int Lgm_GlauertAndHorneDxxBounceAvg( int Version, double Alpha0,  double Ek,  double L,  void *BwFuncData, double (*BwFunc)( double, void * ), double n1, double n2, double n3, double aStarEq,  int Directions, double w1, double w2, double wm, double dw, double x1, double x2, int numberOfWaveNormalAngleDistributions, double *xm, double *dx, double *weightsOnWaveNormalAngleDistributions,int WaveMode, int Species, double MaxWaveLat, double *Daa_ba,  double *Dap_ba,  double *Dpp_ba);
int Lgm_GlauertAndHorneDxxBounceAvg( int Version, double Alpha0,  double Ek,  double L,  void *BwFuncData, double (*BwFunc)( double, void * ), double n1, double n2, double n3, double aStarEq,  int Directions, double w1, double w2, double wm, double dw, double x1, double x2, int numberOfWaveNormalAngleDistributions, double *xm, double *dx, double *weightsOnWaveNormalAngleDistributions,int WaveMode, int Species, double MaxWaveLat, int nNw, int nPlasmaParameters, double aStarMin, double aStarMax, double *Nw, double *Daa_ba,  double *Dap_ba,  double *Dpp_ba);
int Lgm_SummersDxxDerivsBounceAvg( int DerivScheme, double ha, int Version, double Alpha0,  double Ek,  double L,  void *BwFuncData, double (*BwFunc)(), double n1, double n2, double n3, double aStarEq,  int Directions, double w1, double w2, double wm, double dw, int WaveMode, int Species, double MaxWaveLat, double *dDaa,  double *dDap);
double Lgm_ePlasmaFreq( double Density );
double  Lgm_GyroFreq( double q, double B, double m );
double CdipIntegrand_Sb( double Lat, _qpInfo *qpInfo );
double SummersIntegrand_Gaa( double Lat, _qpInfo *qpInfo );
double SummersIntegrand_Gap( double Lat, _qpInfo *qpInfo );
double SummersIntegrand_Gpp( double Lat, _qpInfo *qpInfo );
//double Lgm_GlauertAndHorneHighFrequencyDiffusionCoefficients_Local(double SinAlpha2, double E, double dBoverB2, double BoverBeq, double Omega_e, double Omega_Sig, double Rho, double Sig, double wxl, double wxh, double wxm, double wdx, double xmin, double xmax, int numberOfWaveNormalAngleDistributions, double *xmArray, double *dxArray, double *weightsOnWaveNormalAngleDistributions, double Lambda, int s, double aStar, int Directions, int tensorFlag);
double Lgm_GlauertAndHorneHighFrequencyDiffusionCoefficients_Local(double SinAlpha2, double E, double dBoverB2, double BoverBeq, double Omega_e, double Omega_Sig, double Rho, double Sig, double wxl, double wxh, double wxm, double wdx, double xmin, double xmax, int numberOfWaveNormalAngleDistributions, double *xmArray, double *dxArray, double *weightsOnWaveNormalAngleDistributions, double Lambda, int s, double aStar, int Directions, int tensorFlag, int nNw, int nPlasmaParameters, double aStarMin, double aStarMax, double *Nw);
double Lgm_SummersDaaLocal( double SinAlpha2, double E, double dBoverB2, double BoverBeq, double Omega_e, double Omega_Sig, double Rho, double Sig, double xl, double xh, double xm, double dx, double Lambda, int s, double aStar, int Directions );
double Lgm_SummersDapLocal( double SinAlpha2, double E, double dBoverB2, double BoverBeq, double Omega_e, double Omega_Sig, double Rho, double Sig, double xl, double xh, double xm, double dx, double Lambda, int s, double aStar, int Directions );
double Lgm_SummersDppLocal( double SinAlpha2, double E, double dBoverB2, double BoverBeq, double Omega_e, double Omega_Sig, double Rho, double Sig, double xl, double xh, double xm, double dx, double Lambda, int s, double aStar, int Directions );
double Lgm_SummersDaaLocal_2007( double SinAlpha2, double E, double dBoverB2, double BoverBeq, double Omega_e, double Omega_Sig, double Rho, double Sig, double xl, double xh, double xm, double dx, double Lambda, int s, double n1, double n2, double n3, double aStar, int Directions );
double Lgm_SummersDapLocal_2007( double SinAlpha2, double E, double dBoverB2, double BoverBeq, double Omega_e, double Omega_Sig, double Rho, double Sig, double xl, double xh, double xm, double dx, double Lambda, int s, double n1, double n2, double n3, double aStar, int Directions );
double Lgm_SummersDppLocal_2007( double SinAlpha2, double E, double dBoverB2, double BoverBeq, double Omega_e, double Omega_Sig, double Rho, double Sig, double xl, double xh, double xm, double dx, double Lambda, int s, double n1, double n2, double n3, double aStar, int Directions );
int Lgm_SummersFindSingularities( double  (*f)( double, _qpInfo *), _qpInfo *qpInfo, int Verbose, double Lat0, double Lat1, double *x, double *y );
int Lgm_SummersFindCutoffs( double  (*f)( double, _qpInfo *), _qpInfo *qpInfo, int Verbose, double inc, double Lat0, double Lat1, double *x1, double *x2 );
int Lgm_SummersFindCutoffs2( double  (*f)( double, _qpInfo *), _qpInfo *qpInfo, int Verbose, double Lat0, double Lat1, double *a, double *b );
#endif
