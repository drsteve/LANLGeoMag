#ifndef LGM_SUMMERS_DIFF_COEFF_H
#define LGM_SUMMERS_DIFF_COEFF_H
#include <math.h>
#include "Lgm/Lgm_CTrans.h"
#include "Lgm/Lgm_QuadPack.h"
#include "Lgm/Lgm_PolyRoots.h"

#ifndef LGM_c
#define LGM_c           2.99792458e8                // Spped of light  m/s
#endif

#ifndef LGM_e
#define LGM_e           1.60217656535e-19           // Elementary charge, Coulombs (or s·A)
#endif


#ifndef LGM_Eps0
#define LGM_Eps0        8.854187817620e-12          // Permittivity of free space, A^2·s^4·kg^−1·m^−3
#endif

#ifndef LGM_U
#define LGM_U           (1.66053878283e-27)         // Atomic mass unit, kg
#endif

#ifndef LGM_ELECTRON_MASS
#define LGM_ELECTRON_MASS   9.1093821545e-31    // Rest mass of electron, kg
#endif

#ifndef LGM_PROTON_MASS
#define LGM_PROTON_MASS  (1.0072764667710*LGM_U)    // Proton rest mass,  kg
#endif

#ifndef LGM_EPS
#define LGM_EPS  0.000544617021983          // me/mp
#endif

#ifndef LGM_OXYGEN_MASS
#define LGM_OXYGEN_MASS (15.9994*LGM_U)     // Oxygen rest mass,  kg
#endif

#ifndef LGM_Ee0
#define LGM_Ee0     0.510998910             // Electron rest energy in MeV
#endif

#ifndef LGM_Ep0
#define LGM_Ep0   938.27201323              // Proton rest energy in MeV
#endif

#ifndef M_CDIP
#define M_CDIP  30040.680710                // dipole moment (in nT) for center dipole field (computed for Aug 12, 2004  12.34567 UTC) Need to generalize this.
#endif




#define LGM_R_MODE_WAVE     0
#define LGM_L_MODE_WAVE     1

#define LGM_ELECTRONS       0
#define LGM_PROTONS         1

typedef struct Lgm_SummersInfo {

    double  Alpha0;         //!< Equatorial pitch angle, \f$ \alpha_\circ \f$, [radians].

    double  SinAlpha0;      //!< \f$ \sin\alpha_\circ \f$
    double  CosAlpha0;      //!< \f$ \cos\alpha_\circ \f$

    double  SinAlpha02;     //!< \f$ \sin^2\alpha_\circ \f$
    double  CosAlpha02;     //!< \f$ \cos^2\alpha_\circ \f$

    double  TanAlpha0;      //!< \f$ \tan\alpha_\circ \f$
    double  TanAlpha02;     //!< \f$ \tan^2\alpha_\circ \f$

    double  E;              //!< Normalized dimensionless energy, \f$ E_k/E_\sigma\circ \f$, where \f$ E_\sigma\circ \f$ is rest energy of species.
    double  L;              //!< Dipole L-shell value
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

int Lgm_SummersDxxBounceAvg( double Alpha0,  double Ek,  double L,  void *BwFuncData, double (*BwFunc)( double, void * ), double aStarEq,  double w1, double w2, double wm, double dw, int WaveMode, int Species, double MaxWaveLat, double *Daa_ba,  double *Dap_ba,  double *Dpp_ba);             
double Lgm_ePlasmaFreq( double Density );
double  Lgm_GyroFreq( double q, double B, double m );
double CdipIntegrand_Sb( double Lat, _qpInfo *qpInfo );
double SummersIntegrand_Gaa( double Lat, _qpInfo *qpInfo );
double SummersIntegrand_Gap( double Lat, _qpInfo *qpInfo );
double SummersIntegrand_Gpp( double Lat, _qpInfo *qpInfo );
double Lgm_SummersDaaLocal( double SinAlpha2, double E, double dBoverB2, double BoverBeq, double Omega_e, double Omega_Sig, double Rho, double Sig, double xl, double xh, double xm, double dx, double Lambda, double s, double aStar );
double Lgm_SummersDapLocal( double SinAlpha2, double E, double dBoverB2, double BoverBeq, double Omega_e, double Omega_Sig, double Rho, double Sig, double xl, double xh, double xm, double dx, double Lambda, double s, double aStar );
double Lgm_SummersDppLocal( double SinAlpha2, double E, double dBoverB2, double BoverBeq, double Omega_e, double Omega_Sig, double Rho, double Sig, double xl, double xh, double xm, double dx, double Lambda, double s, double aStar );

#endif
