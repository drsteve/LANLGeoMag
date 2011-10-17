#ifndef LGM_SUMMERS_DIFF_COEFF_H
#define LGM_SUMMERS_DIFF_COEFF_H
#include <math.h>
#include "Lgm/Lgm_CTrans.h"
#include "Lgm/Lgm_QuadPack.h"
#include "Lgm/Lgm_PolyRoots.h"
#include "Lgm/Lgm_Constants.h"


#ifndef M_CDIP
#define M_CDIP  30040.680710                // dipole moment (in nT) for center dipole field (computed for Aug 12, 2004  12.34567 UTC) Need to generalize this.
#endif




#define LGM_R_MODE_WAVE     0
#define LGM_L_MODE_WAVE     1

#define LGM_ELECTRONS       0
#define LGM_PROTONS         1

#define LGM_SUMMERS_2005    2005
#define LGM_SUMMERS_2007    2007

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
    double  MaxWaveLat;     //!< Latitudinal cuttoff for waves. I.e. assume no waves at lats greater than +MaxWaveLat or less than -MaxWaveLat.
    double  Sig;            //!< Additional paramter to define 'semi-bandwidth' which is equal to \f$ \sigma d\omega \f$.
    double  s;              //!< Defines wave mode (s=1 for R-mode, s=-1 for L-mode).
    double  Lambda;         //!< Defines sepcies (\f$ \lambda=-1 \f$ for electrons \f$ \lambda = \epsilon \f$ for protons (\f$\epsilon =m_e/m_p \f$).
    double  Rho;            //!< Rho = sqrt(PI)/2.0( erf( (wm-w1)/dw ) + erf( (w2-wm)/dw )

} Lgm_SummersInfo;

int Lgm_SummersDxxBounceAvg( int Version, double Alpha0,  double Ek,  double L,  void *BwFuncData, double (*BwFunc)( double, void * ), double n1, double n2, double n3, double aStarEq,  double w1, double w2, double wm, double dw, int WaveMode, int Species, double MaxWaveLat, double *Daa_ba,  double *Dap_ba,  double *Dpp_ba);             
double Lgm_ePlasmaFreq( double Density );
double  Lgm_GyroFreq( double q, double B, double m );
double CdipIntegrand_Sb( double Lat, _qpInfo *qpInfo );
double SummersIntegrand_Gaa( double Lat, _qpInfo *qpInfo );
double SummersIntegrand_Gap( double Lat, _qpInfo *qpInfo );
double SummersIntegrand_Gpp( double Lat, _qpInfo *qpInfo );
double Lgm_SummersDaaLocal( double SinAlpha2, double E, double dBoverB2, double BoverBeq, double Omega_e, double Omega_Sig, double Rho, double Sig, double xl, double xh, double xm, double dx, double Lambda, double s, double aStar );
double Lgm_SummersDapLocal( double SinAlpha2, double E, double dBoverB2, double BoverBeq, double Omega_e, double Omega_Sig, double Rho, double Sig, double xl, double xh, double xm, double dx, double Lambda, double s, double aStar );
double Lgm_SummersDppLocal( double SinAlpha2, double E, double dBoverB2, double BoverBeq, double Omega_e, double Omega_Sig, double Rho, double Sig, double xl, double xh, double xm, double dx, double Lambda, double s, double aStar );
double Lgm_SummersDaaLocal_2007( double SinAlpha2, double E, double dBoverB2, double BoverBeq, double Omega_e, double Omega_Sig, double Rho, double Sig, double xl, double xh, double xm, double dx, double Lambda, double s, double n1, double n2, double n3, double aStar );
double Lgm_SummersDapLocal_2007( double SinAlpha2, double E, double dBoverB2, double BoverBeq, double Omega_e, double Omega_Sig, double Rho, double Sig, double xl, double xh, double xm, double dx, double Lambda, double s, double n1, double n2, double n3, double aStar );
double Lgm_SummersDppLocal_2007( double SinAlpha2, double E, double dBoverB2, double BoverBeq, double Omega_e, double Omega_Sig, double Rho, double Sig, double xl, double xh, double xm, double dx, double Lambda, double s, double n1, double n2, double n3, double aStar );
int Lgm_SummersFindSingularities( double  (*f)( double, _qpInfo *), _qpInfo *qpInfo, int Verbose, double Lat0, double Lat1, double *x, double *y );
int Lgm_SummersFindCutoffs( double  (*f)( double, _qpInfo *), _qpInfo *qpInfo, int Verbose, double Lat0, double Lat1, double *x1, double *x2 );
#endif
