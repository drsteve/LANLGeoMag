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




typedef struct Lgm_LstarInfo {

    double              KineticEnergy;  // Particle kinetic energy
    double              Mass;           // Particle mass
    double              PitchAngle;     // Particle Pitch Angle
    double              LSimpleMax;     // Threshold for doing
                                        // drift-shell calculation.

    Lgm_MagModelInfo	*mInfo;

    /*
     *  Variables to hold info on field lines defining the Drift Shell
     */
    int         FindShellPmin;      // Find the Bmin location on each FL
    int         ComputeVgc;         // Compute the gradient of I and Vgc
    int         SaveShellLines;     // only save them if this is true
    int         nFieldPnts[100];    // number of points in each FL.
    double      s_gsm[100][1000];   // distance along FL.
    double      Bmag[100][1000];    // Field magnitude
    double      x_gsm[100][1000];
    double      y_gsm[100][1000];
    double      z_gsm[100][1000];


    /*
     *  Variables to hold info on footprint of Drift Shell
     */
    int                 nPnts;
    double              MLT[100], mlat[100];

    Lgm_Vector          Spherical_Footprint_Pn[100];
    double              Spherical_Footprint_Sn[100];
    double              Spherical_Footprint_Bn[100];

    Lgm_Vector          Spherical_Footprint_Ps[100];
    double              Spherical_Footprint_Ss[100];
    double              Spherical_Footprint_Bs[100];

    Lgm_Vector          Ellipsoid_Footprint_Pn[100];
    double              Ellipsoid_Footprint_Sn[100];
    double              Ellipsoid_Footprint_Bn[100];

    Lgm_Vector          Ellipsoid_Footprint_Ps[100];
    double              Ellipsoid_Footprint_Ss[100];
    double              Ellipsoid_Footprint_Bs[100];

    Lgm_Vector          Mirror_Pn[100];
    double              Mirror_Sn[100];

    Lgm_Vector          Mirror_Ps[100];
    double              Mirror_Ss[100];

    double              PhiVal[100], AngularVelocity[100];

    double              Sb0;        // Equatorial value of Sb Integral.
    double              d2B_ds2;    // second derivative of B wrt s at equator.
    double              RofC;       // radius of curvature at Bmin point.

    double              I0;
    double              I[100];

    int                 ComputeSbIntegral;
    double              SbIntegral0; // Sb Integral on initial FL (i.e. one with Sat through it.)
    //double              SbIntegral[100]; // havent added code to compute these yet (do we need to?)

    Lgm_Vector          Bmin[100];
    Lgm_Vector          Pmin[100];
    Lgm_Vector          GradI[100];
    Lgm_Vector          Vgc[100];

    int                 DriftOrbitType;         // e.g. Open, Closed, Shabansky
    int                 nMinMax;                // Number of valid FLs represented in nMinima[] and nMaxima[] (we may have bailed early)
    int                 nMinima[100];           // # of minima on FL
    int                 nMaxima[100];           // # of maxima on FL (not including endpoints

    int                 nSplnPnts;
    double              xa[500], ya[500], y2[500];

    double	            Phi;

    gsl_interp_accel    *acc;
    gsl_interp          *pspline;



    /*
     *  Other stuff
     */
    int		            VerbosityLevel;
    char                PreStr[64], PostStr[64];

    double	            LS;
    double	            LS_dip_approx;
    double	            LS_McIlwain_M;

    int                 m;
    double              xma[500], yma[500], ym2[500];


    /*
     * Variables to keep track of the I-I0 versus mlat variations that result
     * in search attempts in FindShellLine(). These vars are used in Lstar() to
     * help figure out ranges of mlat to search over.
     */
    double  MLATarr[500]; // mlat's
    double  ImI0arr[500]; // I-I0's
    double  Earr[500];    // nominally the error on I-I0 (typically set to const).
    int     nImI0;        // number of vals stored.



    /*
     *  variables for keeping track of particles
     */
    int                 nParticles;
    Lgm_Vector          Particles[5000];


} Lgm_LstarInfo;


void        SetLstarTolerances( int Quality, Lgm_LstarInfo *LstarInfo );
Lgm_LstarInfo  *InitLstarInfo( int VerbosityLevel );
//void Lgm_InitMagInfoDefaults( Lgm_MagModelInfo  * );
void Lgm_InitLstarInfoDefaults( Lgm_LstarInfo   *LstarInfo );

void FreeLstarInfo( Lgm_LstarInfo *LstarInfo );
Lgm_LstarInfo *Lgm_CopyLstarInfo( Lgm_LstarInfo *s );

int         Grad_I( Lgm_Vector *vin, Lgm_Vector *GradI, Lgm_LstarInfo *LstarInfo );
int         ComputeVcg( Lgm_Vector *vin, Lgm_Vector *Vcg, Lgm_LstarInfo *LstarInfo );
int 	    FindBmRadius( double Bm, double MLT, double mlat, double *r, double tol, Lgm_LstarInfo *LstarInfo );
int 	    FindShellLine( double I0, double *Ifound, double Bm, double MLT, double *mlat, double *rad, double mlat0, double mlat1, double mlat2, int *Iterations, Lgm_LstarInfo *LstarInfo );
double      ComputeI_FromMltMlat( double Bm, double MLT, double mlat, double *r, double I0, Lgm_LstarInfo *LstarInfo );
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


#endif
