#ifndef LGM_PHASE_SPACE_DENSITY_H
#define LGM_PHASE_SPACE_DENSITY_H
#include "Lgm_DynamicMemory.h"

#define LGM_Ee0     0.510998910  // Electron rest energy in MeV


typedef struct Lgm_PhaseSpaceDensity {


    /*
     * Params for measured Differential Flux versus Energy and Pitch Angle -- native or Lo-Res.
     */
    int         nE1;                //!< Number of energy bins in Flux array.

    double      *E1;                //!< Array of energies corresponding to energy bins in flux array.

    int         nA1;                //!< Number of pitch angle bins in Flux array.

    double      *A1;                //!< Array of pitch angles corresponding to pitch angle bins in flux array.

    double      **PSD_EA1;          //!< Array of PSD at native resolution, PSD[E1][A1].

    double      **FLUX_EA1;         //!< Array of measured differential flux (i.e. at antive resolution), Flux[E1][A1].



    /*
     * Params for Phase Space Density versus Energy and Pitch Angle - Hi-Res.
     */
    int         nE2;                //!< Number of energy bins in Hi-Res PSD array.

    double      *E2;                //!< Array of energies corresponding to energy bins in Hi-Res PSD array.

    int         nA2;                //!< Number of pitch angle bins in Hi-Res PSD array.

    double      *A2;                //!< Array of pitch angles corresponding to pitch angle bins in Hi-Res PSD array.

    double      **PSD_EA2;          //!< Hi-Res Array of PSD,  PSD[E2][A2].

    


    /*
     * Params for Phase Space Density versus Mu and K -- native or Lo-Res.
     */
    int         nM1;                //!< Number of energy bins in Flux array.

    double      *M1;                //!< Array of energies corresponding to energy bins in flux array.

    int         nK1;                //!< Number of pitch angle bins in Flux array.

    double      *K1;                //!< Array of pitch angles corresponding to pitch angle bins in flux array.

    double      **PSD_MK1;          //!< Array of measured differential flux. Flux[E][alpha]
    


    /*
     * Params for Phase Space Density versus Mu and K -- Hi-Res.
     */
    int         nM2;                //!< Number of energy bins in Flux array.

    double      *M2;                //!< Array of energies corresponding to energy bins in flux array.

    int         nK2;                //!< Number of pitch angle bins in Flux array.

    double      *K2;                //!< Array of pitch angles corresponding to pitch angle bins in flux array.

    double      **PSD_MK2;          //!< Array of measured differential flux. Flux[E][alpha]
    


    /*
     * Other things..
     */
    int         DumpDiagnostics;    //!< If true, some diagnostics (images, etc) may get dumped out.


} Lgm_PhaseSpaceDensity;



Lgm_PhaseSpaceDensity *Lgm_InitPhaseSpaceDensity( double **J, double *E, double *A, int nE, int nA, int DumpDiagnostics );
void Lgm_FreePhaseSpaceDensity( Lgm_PhaseSpaceDensity *p );
void UpSizeImage( double **Orig, double *X, double *Y, int M, int N, double **New, double *U, double *V, int K, int J );





double Lgm_Energy_to_Mu( double E, double a, double B );
double Lgm_Mu_to_Energy( double Mu, double a, double B );
double Lgm_p2c2( double Ek, double E0 );
double Lgm_v2overc2( double Ek, double E0 );
double Lgm_gamma( double Ek, double E0 );
double Lgm_PsdToDiffFlux( double f, double p2c2 );
double Lgm_DiffFluxToPsd( double j, double p2c2 );



#endif
