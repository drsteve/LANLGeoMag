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
#define EE              (1.6022e-19)        /* electron charge, C */
                                                                                                                                                                                              
                                                                                                                                                                                              



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
    double              I[100];
                                                                                                                                                  
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
     *  variables for keeping track of particles
     */
    int                 nParticles;
    Lgm_Vector          Particles[5000];


} Lgm_LstarInfo;


void        SetLstarTolerances( int Quality, Lgm_LstarInfo *LstarInfo );
Lgm_LstarInfo  *InitLstarInfo( int VerbosityLevel );
void FreeLstarInfo( Lgm_LstarInfo *LstarInfo );
Lgm_LstarInfo *Lgm_CopyLstarInfo( Lgm_LstarInfo *s );

double      AlphaOfK( double K, Lgm_LstarInfo *LstarInfo );
int         Grad_I( Lgm_Vector *vin, Lgm_Vector *GradI, Lgm_LstarInfo *LstarInfo );
int         ComputeVcg( Lgm_Vector *vin, Lgm_Vector *Vcg, Lgm_LstarInfo *LstarInfo );
int 	    FindBmRadius( double Bm, double MLT, double mlat, double *r, double tol, Lgm_LstarInfo *LstarInfo );
int 	    FindShellLine( double I0, double *Ifound, double Bm, double MLT, double *mlat, double *rad, double mlat0, double mlat1, Lgm_LstarInfo *LstarInfo );
void 	    spline( double *x, double *y, int n, double yp1, double ypn, double *y2);
void 	    splint( double *xa, double *ya, double *y2a, int n, double x, double *y);
void 	    quicksort( unsigned long n, double *arr );
void 	    quicksort2(  unsigned long n, double *arr, double *brr );
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
