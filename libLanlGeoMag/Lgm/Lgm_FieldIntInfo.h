#include <math.h>
#include "Lgm_QuadPack.h"
#include "gsl/gsl_errno.h"
#include "gsl/gsl_spline.h"

#define ELECTRON_MASS   (9.10938188e-31)    /* Electron Mass , kg */
#define AMU             (1.660538e-27)      /* kg */
#define PROTON_MASS     (1.00794*AMU)       /* kg */
#define OXYGEN_MASS     (15.9994*AMU)       /* kg */
#define RE              (6378.135e3)        /* Earth Radius, m */
#define CC              (2.99792458e8)      /* Speed of Light, m/s */
#define EE              (1.6022e-19)        /* electron charge, C */
                                                                                                                                                                                              
                                                                                                                                                                                              



typedef struct Lgm_FieldIntInfo {

    double              KineticEnergy;  // Particle kinetic energy
    double              Mass;           // Particle mass
    double              PitchAngle;     // Particle Pitch Angle

    Lgm_Vector	        Pm_South;
    Lgm_Vector	        Pm_North;
    double	            Bm, Sm_South, Sm_North;
    int		            FirstCall;
    int		            n_I_integrand_Calls;
    int		            n_Sb_integrand_Calls;
    double	            epsabs, epsrel;

    /*
     * Arrays containing FL vals
     */
    double              s[1000];    // distance along FL
    Lgm_Vector          P[1000];    // 3D points along FL  (in GSM)
    Lgm_Vector          Bvec[1000]; // 3D B-field vector   (in GSM)
    double              Bmag[1000]; // magnitude of B
    int                 nPnts;      // actual number of points defined
    

    /*
     *  Other stuff
     */
    int		            VerbosityLevel;


} Lgm_FieldIntInfo;



double 	    Iinv( Lgm_FieldIntInfo *Lgm_FieldIntInfo );
double 	    I_integrand( double s, _qpInfo *qpInfo );
double 	    SbIntegral( Lgm_FieldIntInfo *Lgm_FieldIntInfo );
double 	    Sb_integrand( double s, _qpInfo *qpInfo );
void        ratint( double *xa, double *ya, int n, double x, double *y, double *dy );
void        polint(double *xa, double *ya, int n, double x, double *y, double *dy);
void        Interp( double xa[],  double ya[], long int n, double x, double *y);
void        Interp2( double xa[],  double ya[], long int n, double x, double *y);
double      LFromIBmM_Hilton( double I, double Bm, double M );
double      IFromLBmM_Hilton( double L, double Bm, double M );
double      LFromIBmM_McIlwain( double I, double Bm, double M );
double      IFromLBmM_McIlwain( double L, double Bm, double M );




/*
 *    $Id$
 */

