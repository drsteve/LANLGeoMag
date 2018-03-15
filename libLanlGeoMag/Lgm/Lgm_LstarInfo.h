#ifndef LGM_LSTAR_INFO_H
#define LGM_LSTAR_INFO_H

#include <math.h>
#include "Lgm_QuadPack.h"
#include "Lgm_MagModelInfo.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#define ELECTRON_MASS   (9.10938188e-31)    /* Electron Mass , kg */
#define AMU             (1.660538e-27)      /* kg */
#define PROTON_MASS     (1.00794*AMU)       /* kg */
#define OXYGEN_MASS     (15.9994*AMU)       /* kg */
#define RE              (6378.135e3)        /* Earth Radius, m */
#define CC              (2.99792458e8)      /* Speed of Light, m/s */

#define LGM_DRIFT_ORBIT_CLOSED              1
#define LGM_DRIFT_ORBIT_CLOSED_SHABANSKY    2
#define LGM_DRIFT_ORBIT_OPEN                3
#define LGM_DRIFT_ORBIT_OPEN_SHABANSKY      4

#define LGM_SHABANSKY_IGNORE        0
#define LGM_SHABANSKY_HALVE_I       1
#define LGM_SHABANSKY_REJECT        2

#define LGM_LSTARINFO_MAX_FL        300
#define LGM_LSTARINFO_MAX_MINIMA    300


typedef struct Lgm_LstarInfo {

    int         nFLsInDriftShell;   //!< Number of Field Lines to use when constructing Drift Shell.
    int         LstarQuality;       //!< Quality factor to use [0,8] -- higher gives more precise results.
    int         ShabanskyHandling;  //!< Set how Lstar calculations are to handle Shabansky orbits

    double      KineticEnergy;      //!< Particle kinetic energy (only for energy dep. quantities.)
    double      Mass;               //!< Particle mass
    double      PitchAngle;         //!< Particle Pitch Angle
    double      LSimpleMax;         //!< Threshold for doing drift-shell calculation.


    Lgm_MagModelInfo	*mInfo;

    /*
     *  Variables to hold info on field lines defining the Drift Shell
     */
    int         FindShellPmin;      //!< Find the Bmin location on each FL
    int         ComputeVgc;         //!< Compute the gradient of I and Vgc
    int         SaveShellLines;     //!< only save them if this is true
    int         nFieldPnts[ LGM_LSTARINFO_MAX_FL ];    //!< number of points in each FL.
    double      s_gsm[ LGM_LSTARINFO_MAX_FL ][1000];   //!< distance along FL.
    double      Bmag[ LGM_LSTARINFO_MAX_FL ][1000];    //!< Field magnitude
    double      x_gsm[ LGM_LSTARINFO_MAX_FL ][1000];
    double      y_gsm[ LGM_LSTARINFO_MAX_FL ][1000];
    double      z_gsm[ LGM_LSTARINFO_MAX_FL ][1000];


    /*
     *  Variables to hold info on footprint of Drift Shell
     */
    int                 nPnts;
    double              MLT[ LGM_LSTARINFO_MAX_FL ], mlat[ LGM_LSTARINFO_MAX_FL ];

    Lgm_Vector          Spherical_Footprint_Pn[ LGM_LSTARINFO_MAX_FL ];
    double              Spherical_Footprint_Sn[ LGM_LSTARINFO_MAX_FL ];
    double              Spherical_Footprint_Bn[ LGM_LSTARINFO_MAX_FL ];

    Lgm_Vector          Spherical_Footprint_Ps[ LGM_LSTARINFO_MAX_FL ];
    double              Spherical_Footprint_Ss[ LGM_LSTARINFO_MAX_FL ];
    double              Spherical_Footprint_Bs[ LGM_LSTARINFO_MAX_FL ];

    Lgm_Vector          Ellipsoid_Footprint_Pn[ LGM_LSTARINFO_MAX_FL ];
    double              Ellipsoid_Footprint_Sn[ LGM_LSTARINFO_MAX_FL ];
    double              Ellipsoid_Footprint_Bn[ LGM_LSTARINFO_MAX_FL ];

    Lgm_Vector          Ellipsoid_Footprint_Ps[ LGM_LSTARINFO_MAX_FL ];
    double              Ellipsoid_Footprint_Ss[ LGM_LSTARINFO_MAX_FL ];
    double              Ellipsoid_Footprint_Bs[ LGM_LSTARINFO_MAX_FL ];

    Lgm_Vector          Mirror_Pn[ LGM_LSTARINFO_MAX_FL ];
    double              Mirror_Sn[ LGM_LSTARINFO_MAX_FL ];

    Lgm_Vector          Mirror_Ps[ LGM_LSTARINFO_MAX_FL ];
    double              Mirror_Ss[ LGM_LSTARINFO_MAX_FL ];

    double              PhiVal[ LGM_LSTARINFO_MAX_FL ], AngularVelocity[ LGM_LSTARINFO_MAX_FL ];

    double              Sb0;        // Equatorial value of Sb Integral.
    double              d2B_ds2;    // second derivative of B wrt s at equator.
    double              RofC;       // radius of curvature at Bmin point.

    double              I0;
    double              I[ LGM_LSTARINFO_MAX_FL ];

    int                 ComputeSbIntegral;
    double              SbIntegral0; // Sb Integral on initial FL (i.e. one with Sat through it.)
    //double              SbIntegral[ LGM_LSTARINFO_MAX_FL ]; // havent added code to compute these yet (do we need to?)

    Lgm_Vector          Bmin[ LGM_LSTARINFO_MAX_FL ];
    Lgm_Vector          Pmin[ LGM_LSTARINFO_MAX_FL ];
    Lgm_Vector          GradI[ LGM_LSTARINFO_MAX_FL ];
    Lgm_Vector          Vgc[ LGM_LSTARINFO_MAX_FL ];

    int                 DriftOrbitType;         // e.g. Open, Closed, Shabansky
    int                 nMinMax;                // Number of valid FLs represented in nMinima[] and nMaxima[] (we may have bailed early)
    int                 nMinima[ LGM_LSTARINFO_MAX_MINIMA ];           // # of minima on FL
    int                 nMaxima[ LGM_LSTARINFO_MAX_MINIMA ];           // # of maxima on FL (not including endpoints

    int                 nSplnPnts;
    double              xa [ 3*LGM_LSTARINFO_MAX_FL ], ya[ 3*LGM_LSTARINFO_MAX_FL ], y2[ 3*LGM_LSTARINFO_MAX_FL ];

    double	            Phi;

    gsl_interp_accel    *acc;
    gsl_interp          *pspline;



    /*
     *  Other stuff
     */
    int		            VerbosityLevel;
    char                PreStr[64], PostStr[64];
    int                 ISearchMethod;

    double	            LS;
    double	            LS_dip_approx;
    double	            LS_McIlwain_M;

    int                 m;
    double              xma[ 3*LGM_LSTARINFO_MAX_FL ], yma[ 3*LGM_LSTARINFO_MAX_FL ], ym2[ 3*LGM_LSTARINFO_MAX_FL ];


    /*
     * Variables to keep track of the I-I0 versus mlat variations that result
     * in search attempts in FindShellLine(). These vars are used in Lstar() to
     * help figure out ranges of mlat to search over.
     */
    double  MLATarr[ 3*LGM_LSTARINFO_MAX_FL ]; // mlat's
    double  ImI0arr[ 3*LGM_LSTARINFO_MAX_FL ]; // I-I0's
    double  Earr[ 3*LGM_LSTARINFO_MAX_FL ];    // nominally the error on I-I0 (typically set to const).
    int     nImI0;        // number of vals stored.



    /*
     *  variables for keeping track of particles
     */
    int                 nParticles;
    Lgm_Vector          Particles[5000];


} Lgm_LstarInfo;


void        Lgm_SetLstarTolerances( int Quality, int nFLsInDriftShell, Lgm_LstarInfo *LstarInfo );
Lgm_LstarInfo  *InitLstarInfo( int VerbosityLevel );
//void Lgm_InitMagInfoDefaults( Lgm_MagModelInfo  * );
void Lgm_InitLstarInfoDefaults( Lgm_LstarInfo   *LstarInfo );

void FreeLstarInfo( Lgm_LstarInfo *LstarInfo );
Lgm_LstarInfo *Lgm_CopyLstarInfo( Lgm_LstarInfo *s );

int         Grad_I( Lgm_Vector *vin, Lgm_Vector *GradI, Lgm_LstarInfo *LstarInfo );
int         ComputeVcg( Lgm_Vector *vin, Lgm_Vector *Vcg, Lgm_LstarInfo *LstarInfo );
int 	    FindBmRadius( double Bm, double MLT, double mlat, double *r, double tol, Lgm_LstarInfo *LstarInfo );
int 	    FindShellLine( double I0, double *Ifound, double Bm, double MLT, double *mlat, double *rad, double mlat0, double mlat1, double mlat2, int *Iterations, Lgm_LstarInfo *LstarInfo );
double      ComputeI_FromMltMlat(  double Bm, double MLT, double mlat, double *r, double I0, Lgm_LstarInfo *LstarInfo );
double      ComputeI_FromMltMlat1( double Bm, double MLT, double mlat, double *r, double I0, Lgm_LstarInfo *LstarInfo );
double      ComputeI_FromMltMlat2( double Bm, double MLT, double mlat, double *r, double I0, Lgm_LstarInfo *LstarInfo );
void 	    spline( double *x, double *y, int n, double yp1, double ypn, double *y2);
void 	    splint( double *xa, double *ya, double *y2a, int n, double x, double *y);
void 	    quicksort( unsigned long n, double *arr );
void 	    quicksort2(  unsigned long n, double *arr, double *brr );
void 	    quicksort3(  unsigned long n, double *arr, double *brr, double *crr );
Lgm_MagModelInfo     *init_info( ) ;
void        NewTimeLstarInfo( long int Date, double UT, double PitchAngle, int (*Mag)(Lgm_Vector *, Lgm_Vector *, Lgm_MagModelInfo *), Lgm_LstarInfo *LstarInfo );
int         Lstar( Lgm_Vector *vin, Lgm_LstarInfo *LstarInfo );
double      MagFlux( Lgm_LstarInfo *LstarInfo );
double      MagFlux2( Lgm_LstarInfo *LstarInfo );
double      MagFluxIntegrand( double Phi, _qpInfo *qpInfo ) ;
double      MagFluxIntegrand2( double Phi, _qpInfo *qpInfo ) ;
double      LambdaIntegrand( double Lambda, _qpInfo *qpInfo ) ;
double      LambdaIntegral( Lgm_LstarInfo *LstarInfo ) ;
double      AngVelInv( double Phi );
int         Lgm_LCDS( long int Date, double UTC, double brac1, double brac2, double Alpha, double LT, double tol, int Quality, int nFLsInDriftShell, double *K, Lgm_LstarInfo *LstarInfo );
 

#endif
