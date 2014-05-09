#ifndef LGM_VEL_STEP_INFO_H
#define LGM_VEL_STEP_INFO_H

#include <math.h>
#include "Lgm/Lgm_Vec.h"
#include "Lgm/Lgm_QuadPack.h"
#include "Lgm/Lgm_Constants.h"


#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif

#define 	LGM_VELSTEP_KMAX	16
#define 	LGM_VELSTEP_IMAX	(LGM_VELSTEP_KMAX+1)
#define 	LGM_VELSTEP_JMAX	(LGM_VELSTEP_KMAX+2)
#define 	LGM_VELSTEP_REDMAX 	1.0e-5
#define 	LGM_VELSTEP_REDMIN 	0.7
#define 	LGM_VELSTEP_SCLMAX 	0.1
#define 	LGM_VELSTEP_SAFE1	0.25
#define 	LGM_VELSTEP_SAFE2	0.70



typedef struct Lgm_VelStepInfo {

    void    *data; // to carry arbitrary data

    int     VerbosityLevel;

    /*
     *  Vars for Bulirsch-Stoer ODE solver. Some of these variables are needed
     *  to make Lgm_VelStep() reentrant/thread-safe.  They basically used to be
     *  static declarations within Lgm_VelStep()
     */
    double      Lgm_VelStep_BS_eps_old;
    int         Lgm_VelStep_BS_FirstTimeThrough;
    int         Lgm_VelStep_BS_kmax;
    int         Lgm_VelStep_BS_kopt;
    double      Lgm_VelStep_BS_snew;
    double      Lgm_VelStep_BS_A[LGM_VELSTEP_JMAX+1];
    double      Lgm_VelStep_BS_alpha[LGM_VELSTEP_IMAX+1][LGM_VELSTEP_IMAX+1];
    double      Lgm_VelStep_BS_d[LGM_VELSTEP_JMAX][LGM_VELSTEP_JMAX];
    double      Lgm_VelStep_BS_x[LGM_VELSTEP_JMAX];
    Lgm_Vector  Lgm_VelStep_BS_u[LGM_VELSTEP_JMAX];
    double      Lgm_VelStep_BS_Eps;        // Eps parameter used in BS method. Influences speed greatly.
    

    /*
     *  These variables are needed to make Lgm_MagStep2() reentrant/thread-safe.
     *  They basically used to be static declarations within Lgm_MagStep2()
     */
    long int    Lgm_nVelEvals; // records number of Bfield evals between resets
    int         Lgm_VelStep_Integrator; // ODE solver to use ( LGM_VELSTEP_ODE_BS or LGM_VELSTEP_ODE_RK5)
    double      Lgm_VelStep_eps_old;
    int	        Lgm_VelStep_FirstTimeThrough;
    int         Lgm_VelStep_kmax;
    int         Lgm_VelStep_kopt;
    double      Lgm_VelStep_snew;
    double      Lgm_VelStep_A[LGM_VELSTEP_JMAX+1];
    double      Lgm_VelStep_alpha[LGM_VELSTEP_IMAX+1][LGM_VELSTEP_IMAX+1];
    double      Lgm_VelStep_d[LGM_VELSTEP_JMAX][LGM_VELSTEP_JMAX];
    double      Lgm_VelStep_x[LGM_VELSTEP_JMAX];
    double      Lgm_VelStep_Tol;        // tolerance for Magstep (ODE solver).
    double      Lgm_VelStep_Alpha;
    double      Lgm_VelStep_SinAlpha;
    double      Lgm_VelStep_q;
    double      Lgm_VelStep_T;
    double      Lgm_VelStep_E0;
    double      Lgm_VelStep_Bm;
    double      Lgm_VelStep_h;
    int         Lgm_VelStep_DerivScheme;


} Lgm_VelStepInfo;



int             Lgm_ModMid2( Lgm_Vector *, Lgm_Vector *, Lgm_Vector *, double, int, double, int (*Velocity)(Lgm_Vector *, Lgm_Vector *, Lgm_VelStepInfo *), Lgm_VelStepInfo * );
int             Lgm_VelStep( Lgm_Vector *, Lgm_Vector *, double, double *, double *, double, double, double *, int *, int (*Velocity)(Lgm_Vector *, Lgm_Vector *, Lgm_VelStepInfo *), Lgm_VelStepInfo * );
Lgm_VelStepInfo *Lgm_InitVelInfo( );
void            Lgm_FreeVelInfo( Lgm_VelStepInfo * );
void            Lgm_VelPolFunExt( int k, double x_k, Lgm_Vector *u_k, Lgm_Vector *w, Lgm_Vector *dw, Lgm_VelStepInfo *Info );
void            Lgm_VelRatFunExt( int k, double x_k, Lgm_Vector *u_k, Lgm_Vector *w, Lgm_Vector *dw, Lgm_VelStepInfo *Info );

#endif
