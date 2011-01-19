#ifndef LGM_MAG_MODEL_INFO_H
#define LGM_MAG_MODEL_INFO_H

#include <math.h>
#include "Lgm_QuadPack.h"
#include "Lgm_CTrans.h"
#include "Lgm_Octree.h"
#include "gsl/gsl_errno.h"
#include "gsl/gsl_spline.h"

#define LGM_ELECTRON_MASS   (9.10938188e-31)    /* Electron Mass , kg */
#define LGM_AMU             (1.660538e-27)      /* kg */
#define LGM_PROTON_MASS     (1.00794*AMU)       /* kg */
#define LGM_OXYGEN_MASS     (15.9994*AMU)       /* kg */
#define RE                  (6378.135e3)        /* Earth Radius, m */
#define LGM_CC              (2.99792458e8)      /* Speed of Light, m/s */
#define LGM_EE              (1.6022e-19)        /* electron charge, C */


#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif

/*
 * Field line types deduced by Lgm_Trace()
 */
#define		LGM_OPEN_IMF		           0
#define		LGM_CLOSED		               1
#define		LGM_OPEN_N_LOBE		           2
#define		LGM_OPEN_S_LOBE		           3
#define		LGM_INSIDE_EARTH           	  -1
#define     LGM_TARGET_HEIGHT_UNREACHABLE -2


#define 	LGM_MAGSTEP_KMAX	9
#define 	LGM_MAGSTEP_IMAX	(LGM_MAGSTEP_KMAX+1)
#define 	LGM_MAGSTEP_JMAX	(LGM_MAGSTEP_KMAX+2)
#define 	LGM_MAGSTEP_REDMAX 	1.0e-5
#define 	LGM_MAGSTEP_REDMIN 	0.7
#define 	LGM_MAGSTEP_SCLMAX 	0.1
#define 	LGM_MAGSTEP_SAFE1	0.25
#define 	LGM_MAGSTEP_SAFE2	0.70

#define     DQAGS   0   // Quadpack routine (works well for I. Dont use for Sb)
#define     DQAGP   1   // Quadpack routine (works well for Sb. Overkill for I)
#define     DQK21   2   // Quadpack simple routine (not terribly accurate).


/*
 * these are used in B_FromScatteredData()
 * They control which interpolation scheme to use.
 * QUADRATIC_DFI doesnt seem to work very well...
 * The best is probably LINEAR_DFI
 */
#define LINEAR          0
#define LINEAR_DFI      1
#define QUADRATIC       2
#define QUADRATIC_DFI   3
#define NEWTON_INTERP   4

#define LGM_CDIP        0
#define LGM_EDIP        1
#define LGM_IGRF        2

#define LGM_MAX_INTERP_PNTS 10000

#define LGM_RELATIVE_JUMP_METHOD 0
#define LGM_ABSOLUTE_JUMP_METHOD 1


typedef struct Lgm_MagModelInfo {

    Lgm_CTrans  *c;                 /* This contains all time info and a bunch more stuff */
    long int	nFunc;
    int 		(*Bfield)();
    int			SavePoints;
    double		Hmax;
    FILE	    *fp;
    double      W[6];
    double      G1, G2;
    int			Kp;
    double      Dst;
    double      P;
    double      Bx, By, Bz;
    double      T96MOD_V[11];       /* free params for T96_MOD */

    double      Trace_s;

    double      B0, B1, B2;             // params for simplified mead field


    /*
     * Variable to control which internal model to use
     */
    int         InternalModel;          // Can be LGM_CDIP, LGM_EDIP or LGM_IGRF



    /*
     * Stuff that may be neeeded for field integrals, etc.
     */
    double              KineticEnergy;  // Particle kinetic energy
    double              Mass;           // Particle mass
    double              PitchAngle;     // Particle Pitch Angle

    Lgm_Vector          Pm_South;
    Lgm_Vector          Pm_North;
    double              Bm, Sm_South, Sm_North, Blocal;
    int                 FirstCall;
//    double              epsabs, epsrel;

    /*
     * Arrays containing FL vals
     */
    double              s[LGM_MAX_INTERP_PNTS];    // distance along FL
    double              Px[LGM_MAX_INTERP_PNTS];   // Px along FL  (in GSM)
    double              Py[LGM_MAX_INTERP_PNTS];   // Py along FL  (in GSM)
    double              Pz[LGM_MAX_INTERP_PNTS];   // Pz along FL  (in GSM)
    Lgm_Vector          Bvec[LGM_MAX_INTERP_PNTS]; // 3D B-field vector   (in GSM)
    double              Bmag[LGM_MAX_INTERP_PNTS]; // magnitude of B
    double              BminusBcdip[LGM_MAX_INTERP_PNTS]; // magnitude of B minus magnitude of Cent. Dipole
    int                 nPnts;      // actual number of points defined
    double              ds;         // spacing in s (dist. along FL)
                                    // this will help in seacrhing the
                                    // arrays (e.g. for interpolation).

    Lgm_Vector      P_gsm;          //< S/C position in GSM
    double          S;              //< Distance along FL from southern footpoint to S/C location in Re.
    double          B;              //< Local (model) B-field magnitude (i.e. at S/C position)

    Lgm_Vector      Pmin;           //< position of minimum |B| in GSM
    Lgm_Vector      Bvecmin;        //< value of Bvecmin
    double          Bmin;           //< Value of |Bmin|
    double          Smin;           //< Distance from southern footpoint to Pmin along FL.

    Lgm_Vector      Spherical_Footprint_Pn;   //< position of northern footpoint (at 120km)
    double          Spherical_Footprint_Sn;   //< Distance along FL from southern foorpoint in Re
    double          Spherical_Footprint_Bn;   //< Value of |B| at Footprint_Pn

    Lgm_Vector      Spherical_Footprint_Ps;   //< position of southern footpoint (at 120km)
    double          Spherical_Footprint_Ss;   //< Distance along FL from southern foorpoint in Re (i.e. this one is zero by definition)
    double          Spherical_Footprint_Bs;   //< Value of |B| at Footprint_Ps

    Lgm_Vector      Ellipsoid_Footprint_Pn;   //< position of northern footpoint (at 120km above surface of WGS84 ellipsoid)
    double          Ellipsoid_Footprint_Sn;   //< Distance along FL from southern footpoint in Re
    double          Ellipsoid_Footprint_Bn;   //< Value of |B| at Footprint_Pn

    Lgm_Vector      Ellipsoid_Footprint_Ps;   //< position of southern footpoint (at 120km)
    double          Ellipsoid_Footprint_Ss;   //< Distance along FL from southern foorpoint in Re (i.e. this one is zero by definition)
    double          Ellipsoid_Footprint_Bs;   //< Value of |B| at Footprint_Ps

    int             FieldLineType;  //< Field line type. (I.e., LGM_OPEN_IMF, LGM_CLOSED, LGM_OPEN_N_LOBE, LGM_OPEN_S_LOBE, LGM_INSIDE_EARTH, LGM_TARGET_HEIGHT_UNREACHABLE)




    double              d2B_ds2;    // second derivative of B wrt to s at smin.
    double              Sb0;        // value of Sb integral for eq. mirroring particles.
    int                 imin1;      // imin1 and imin2 are the indices in the
    int                 imin2;      //   array between which smin is located.



    /*
     *  GSL defs for spline interpolation
     */
    gsl_interp_accel    *acc;       // accelerator
    gsl_interp_accel    *accPx;     // accelerator
    gsl_interp_accel    *accPy;     // accelerator
    gsl_interp_accel    *accPz;     // accelerator
    gsl_spline          *spline;    // spline object
    gsl_spline          *splinePx;  // spline object
    gsl_spline          *splinePy;  // spline object
    gsl_spline          *splinePz;  // spline object



    /*
     *  Other stuff
     */
    int                 VerbosityLevel;
    int                 UseInterpRoutines; // whether to use fast I and Sb routines.


    /*
     *  These variables are needed to make Lgm_MagStep() reentrant/thread-safe.
     *  They basically used to be static declarations within Lgm_MagStep()
     */
    double  Lgm_MagStep_eps_old;
    int	    Lgm_MagStep_FirstTimeThrough;
    int     Lgm_MagStep_kmax;
    int     Lgm_MagStep_kopt;
    double  Lgm_MagStep_snew;
    double  Lgm_MagStep_A[LGM_MAGSTEP_JMAX+1];
    double  Lgm_MagStep_alpha[LGM_MAGSTEP_IMAX+1][LGM_MAGSTEP_IMAX+1];
    double  Lgm_MagStep_d[LGM_MAGSTEP_JMAX][LGM_MAGSTEP_JMAX];
    double  Lgm_MagStep_x[LGM_MAGSTEP_JMAX];


    /*
     *  These variables are needed to make I_integrand() reentrant/thread-safe.
     *  They basically used to be static declarations.
     */
    int         Lgm_I_integrand_FirstCall;
    int         Lgm_I_integrand_JumpMethod;
    double      Lgm_I_integrand_S;
    Lgm_Vector  Lgm_I_integrand_P;
    Lgm_Vector  Lgm_I_integrand_u_scale;
    int         Lgm_n_I_integrand_Calls;
    int         Lgm_I_Integrator;
    double      Lgm_I_Integrator_epsrel;
    double      Lgm_I_Integrator_epsabs;

    /*
     *  These variables are needed to make Sb_integrand() reentrant/thread-safe.
     *  They basically used to be static declarations.
     */
    int         Lgm_Sb_integrand_FirstCall;
    double      Lgm_Sb_integrand_S;
    Lgm_Vector  Lgm_Sb_integrand_P;
    Lgm_Vector  Lgm_Sb_integrand_u_scale;
    int         Lgm_n_Sb_integrand_Calls;
    int         Lgm_Sb_Integrator;              // This variable is not currently used (for now I force
                                                // Sb to use DQAGP. )
    double      Lgm_Sb_Integrator_epsrel;
    double      Lgm_Sb_Integrator_epsabs;


    /*
     * Variables to control MagFlux integration
     */
    int         Lgm_MagFlux_Integrator;         // This variable is not currently used (for now I force
                                                // MagFlux to use DQAGS. )
    double      Lgm_MagFlux_Integrator_epsrel;
    double      Lgm_MagFlux_Integrator_epsabs;


    /*
     * Variables to control LambdaIntegral integration
     */
    int         Lgm_LambdaIntegral_Integrator;  // This variable is not currently used (for now I force
                                                // LambdaIntegral to use DQAGS. )
    double      Lgm_LambdaIntegral_Integrator_epsrel;
    double      Lgm_LambdaIntegral_Integrator_epsabs;


    /*
     * Some other tolerances
     */
    double      Lgm_FindBmRadius_Tol;
    double      Lgm_FindShellLine_I_Tol;
    double      Lgm_TraceToMirrorPoint_Tol;



    /*
     * Variables for defining Octree stuff
     */
    Lgm_OctreeCell  *OctreeRoot;
    int             Octree_kNN_k;
    int             Octree_kNN_InterpMethod;  // 0 = linear Div Feee; 1 = Quadratic Div Free; 2 = Newton 4-point
    double          Octree_kNN_MaxDist;       // in physical units
    double          OctreeScaleMin;
    double          OctreeScaleMax;
    double          OctreeScaleDiff;


    /*
     * limits for tracing. If you go outside this box, we consider it open
     */
    double OpenLimit_xMin;
    double OpenLimit_xMax;
    double OpenLimit_yMin;
    double OpenLimit_yMax;
    double OpenLimit_zMin;
    double OpenLimit_zMax;

    double  Lgm_LossConeHeight;


    /*
     *  Globals for OP77 Model
     */
    double OP77_TILTL;
    double OP77_A[65], OP77_B[65], OP77_C[45], OP77_D[45], OP77_E[65], OP77_F[65], OP77_TT[5];




} Lgm_MagModelInfo;


Lgm_MagModelInfo *Lgm_InitMagInfo( );
void Lgm_FreeMagInfo( Lgm_MagModelInfo  *Info );
Lgm_MagModelInfo *Lgm_CopyMagInfo( Lgm_MagModelInfo *s );

int Lgm_Trace( Lgm_Vector *u, Lgm_Vector *v1, Lgm_Vector *v2, Lgm_Vector *v3, double Height, double TOL1, double TOL2, Lgm_MagModelInfo *Info );
int Lgm_TraceToMinBSurf( Lgm_Vector *, Lgm_Vector *, double, double, Lgm_MagModelInfo * );
int Lgm_TraceToSMEquat(  Lgm_Vector *, Lgm_Vector *, double, Lgm_MagModelInfo * );
int Lgm_TraceToEarth(  Lgm_Vector *, Lgm_Vector *, double, double, double, Lgm_MagModelInfo * );
int Lgm_TraceToSphericalEarth(  Lgm_Vector *, Lgm_Vector *, double, double, double, Lgm_MagModelInfo * );
int Lgm_TraceLine(  Lgm_Vector *, Lgm_Vector *, double, double, double, int, Lgm_MagModelInfo * );
int Lgm_TraceLine2(  Lgm_Vector *, Lgm_Vector *, double, double, double, double, int, Lgm_MagModelInfo * );
void ReplaceFirstPoint( double s, double B, Lgm_Vector *P, Lgm_MagModelInfo *Info );
void AddNewPoint( double s, double B, Lgm_Vector *P, Lgm_MagModelInfo *Info );
void InitSpline( Lgm_MagModelInfo *Info );
void FreeSpline( Lgm_MagModelInfo *Info );
int Lgm_TraceToMinRdotB( Lgm_Vector *, Lgm_Vector *, double, Lgm_MagModelInfo * );
int Lgm_TraceIDL( int, void *argv[] );
int Lgm_TraceToMirrorPoint( Lgm_Vector *u, Lgm_Vector *v, double *Sm, double Bm, double sgn, double tol, Lgm_MagModelInfo *Info );




void Lgm_ModMid( Lgm_Vector *, Lgm_Vector *, double, int, double,
	     int (*Mag)(Lgm_Vector *, Lgm_Vector *, Lgm_MagModelInfo *), Lgm_MagModelInfo * );
void Lgm_RatFunExt( int, double, Lgm_Vector *, Lgm_Vector *, Lgm_Vector *, Lgm_MagModelInfo * );
int  Lgm_MagStep( Lgm_Vector *, Lgm_Vector *, double, double *, double *, double, double, double *, int *,
              int (*Mag)(Lgm_Vector *, Lgm_Vector *, Lgm_MagModelInfo *), Lgm_MagModelInfo * );



/*
 *  B_internal
 *
 *  Function Prototypes for internal models
 */
int Lgm_B_igrf(Lgm_Vector *, Lgm_Vector *, Lgm_MagModelInfo *);
int Lgm_B_cdip(Lgm_Vector *, Lgm_Vector *, Lgm_MagModelInfo *);
int Lgm_B_edip(Lgm_Vector *, Lgm_Vector *, Lgm_MagModelInfo *);



/*
 *
 *  Function Prototypes for Olsen Pfitzer 1977  Model
 *
 */
int Lgm_B_OP77( Lgm_Vector *v, Lgm_Vector *B, Lgm_MagModelInfo *Info );
void OlsenPfitzerStatic( double XX[], double BF[], double TILT, Lgm_MagModelInfo *m );



/*
 *  T87
 *
 *  Function Prototypes for T87 model
 */
int Lgm_B1_T87( Lgm_Vector *, Lgm_Vector *, Lgm_MagModelInfo * );
int Lgm_B2_T87( Lgm_Vector *, Lgm_Vector *, Lgm_MagModelInfo * );
int Lgm_B3_T87( Lgm_Vector *, Lgm_Vector *, Lgm_MagModelInfo * );
int Lgm_B_T87( Lgm_Vector *, Lgm_Vector *, Lgm_MagModelInfo * );


/*
 *  T89
 *
 *  Function Prototypes for T89 model
 */
int Lgm_BM_T89( Lgm_Vector *, Lgm_Vector *, Lgm_MagModelInfo * );
int Lgm_BT_T89( Lgm_Vector *, Lgm_Vector *, Lgm_MagModelInfo * );
int Lgm_BRC_T89( Lgm_Vector *, Lgm_Vector *, Lgm_MagModelInfo * );
int Lgm_BC_T89( Lgm_Vector *, Lgm_Vector *, Lgm_MagModelInfo * );
int Lgm_B_T89( Lgm_Vector *, Lgm_Vector *, Lgm_MagModelInfo * );


/*
 *  T96MOD
 *
 *  Function Prototypes for T96MOD model
 */
int Lgm_B_T96MOD_MGH( Lgm_Vector *v, Lgm_Vector *B, Lgm_MagModelInfo *Info );
void lgm_field_t96mod_mgh_(double *, double *, int *IYEAR,int *IDAY,int *IH,int *IM,double *SEC,double *X,double *Y,double *Z,double *BX,double *BY,double *BZ);
void lgm_field_t96mod_mgh__(double *PARMOD, double *AMDF, int *IYEAR,int *IDAY,int *IH,int *IM,double *SEC,double *X,double *Y,double *Z,double *BX,double *BY,double *BZ);
void lgm_field_t96mod_( int *, int *, int *, int *, double *, double *, double *, double *, double *, double *, double * );


/*
 *  T01S
 *
 *  Function Prototypes for TS04 model
 */
int  Lgm_B_TS04( Lgm_Vector *v, Lgm_Vector *B, Lgm_MagModelInfo *Info );
void Lgm_ComputeW( double W[], int i, double Nk[], double Vk[], double Bsk[], int nk );
void Tsyg_TS04( int IOPT, double *PARMOD, double PS, double SINPS, double COSPS, double X, double Y, double Z, double *BX, double *BY, double *BZ);
void TS04_EXTERN( int IOPGEN, int IOPT, int IOPB, int IOPR, double *A, int NTOT, double PDYN, double DST, double BXIMF, double BYIMF,
                double BZIMF, double W1, double W2, double W3, double W4, double W5, double W6, double PS,
                double X, double Y, double Z, double *BXCF, double *BYCF, double *BZCF, double *BXT1, double *BYT1,
                double *BZT1, double *BXT2, double *BYT2, double *BZT2, double *BXSRC, double *BYSRC, double *BZSRC,
                double *BXPRC, double *BYPRC, double *BZPRC,  double *BXR11, double *BYR11, double *BZR11, double *BXR12,
                double *BYR12, double *BZR12, double *BXR21, double *BYR21, double *BZR21, double *BXR22, double *BYR22,
                double *BZR22, double *HXIMF, double *HYIMF, double *HZIMF, double *BBX, double *BBY, double *BBZ );



/*
 *  TS04
 *
 *  Function Prototypes for TS04 model
 */
int  Lgm_B_T01S( Lgm_Vector *v, Lgm_Vector *B, Lgm_MagModelInfo *Info );
void Tsyg_T01S( int IOPT, double *PARMOD, double PS, double SINPS, double COSPS, double X, double Y, double Z, double *BX, double *BY, double *BZ);
void T01S_EXTALL( int IOPGEN, int IOPT, int IOPB, int IOPR, double *A, int NTOT, double PDYN, double DST, double BYIMF,
                    double BZIMF, double VBIMF1, double VBIMF2, double PS, double X, double Y, double Z,
                    double *BXCF, double *BYCF, double *BZCF, double *BXT1, double *BYT1, double *BZT1,
                    double *BXT2, double *BYT2, double *BZT2, double *BXSRC, double *BYSRC, double *BZSRC,
                    double *BXPRC, double *BYPRC, double *BZPRC,  double *BXR11, double *BYR11, double *BZR11,
                    double *BXR12, double *BYR12, double *BZR12, double *BXR21, double *BYR21, double *BZR21,
                    double *BXR22, double *BYR22, double *BZR22, double *HXIMF, double *HYIMF, double *HZIMF,
                    double *BX, double *BY, double *BZ);

/*
 *  Computing B from scattered data -- (e.g. an irregular mesh)
 *
 */

int Lgm_B_FromScatteredData( Lgm_Vector *v, Lgm_Vector *B, Lgm_MagModelInfo *Info );


/*
 * Simplified Mead field.
 */
int Lgm_SimplifiedMead(Lgm_Vector *v, Lgm_Vector *B, Lgm_MagModelInfo *Info);


/*
 * routines/functions for field integrals, invariants, etc.
 */
double      Iinv( Lgm_MagModelInfo *fInfo );
double      I_integrand( double s, _qpInfo *qpInfo );
double      Iinv_interped( Lgm_MagModelInfo *fInfo );
double      I_integrand_interped( double s, _qpInfo *qpInfo );
double      SbIntegral( Lgm_MagModelInfo *fInfo );
double      Sb_integrand( double s, _qpInfo *qpInfo );
double      SbIntegral_interped( Lgm_MagModelInfo *fInfo );
double      Sb_integrand_interped( double s, _qpInfo *qpInfo );
void        ratint( double *xa, double *ya, int n, double x, double *y, double *dy );
void        polint(double *xa, double *ya, int n, double x, double *y, double *dy);
void        Interp( double xa[],  double ya[], long int n, double x, double *y);
void        Interp2( double xa[],  double ya[], long int n, double x, double *y);
double      LFromIBmM_Hilton( double I, double Bm, double M );
double      IFromLBmM_Hilton( double L, double Bm, double M );
double      LFromIBmM_McIlwain( double I, double Bm, double M );
double      IFromLBmM_McIlwain( double L, double Bm, double M );

double      BofS( double s, Lgm_MagModelInfo *Info );
int         SofBm( double Bm, double *ss, double *sn, Lgm_MagModelInfo *Info );
double      Lgm_AlphaOfK( double K, Lgm_MagModelInfo *Info );
int         Lgm_Init_AlphaOfK( Lgm_DateTime *d, Lgm_Vector *u, Lgm_MagModelInfo *m );
int         Lgm_Grad_I( Lgm_Vector *vin, Lgm_Vector *GradI, Lgm_MagModelInfo *Info );
//int         ComputeVcg( Lgm_Vector *vin, Lgm_Vector *Vcg, Lgm_LstarInfo *LstarInfo );


/*
 * Getter/Setter functions. May not want to do this...
 *  Just testing for now..
 */
void Lgm_MagModelInfo_Set_Psw( double Psw, Lgm_MagModelInfo *m );
void Lgm_MagModelInfo_Set_Kp( double Kp, Lgm_MagModelInfo *m );
void Lgm_Set_Octree_kNN_InterpMethod( Lgm_MagModelInfo *m, int Method );
void Lgm_Set_Octree_kNN_k( Lgm_MagModelInfo *m, int k );
void Lgm_Set_Octree_kNN_MaxDist( Lgm_MagModelInfo *m, double MaxDist );
void Lgm_Set_Open_Limits( Lgm_MagModelInfo *m, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax );
void Lgm_Set_LossConeHeight( Lgm_MagModelInfo *m, double LossConeHeight );

/*
 * Added for Python wrapping, the function pointer is hard to impossible to
 * deal with setting
 */
void Lgm_Set_Lgm_B_igrf(Lgm_MagModelInfo *MagInfo);
void Lgm_Set_Lgm_B_T01S(Lgm_MagModelInfo *MagInfo);
void Lgm_Set_gm_B_TS04(Lgm_MagModelInfo *MagInfo);
void Lgm_Set_Lgm_B_T89(Lgm_MagModelInfo *MagInfo);
void Lgm_Set_Lgm_B_OP77(Lgm_MagModelInfo *MagInfo);






#endif


/*
 *  $Id$
 */
