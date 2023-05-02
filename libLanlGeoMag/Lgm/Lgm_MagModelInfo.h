#ifndef LGM_MAG_MODEL_INFO_H
#define LGM_MAG_MODEL_INFO_H

#include <math.h>
#include "Lgm_QuadPack.h"
#include "Lgm_CTrans.h"
#include "Lgm_Octree.h"
#include "Lgm_KdTree.h"
#include "Lgm_Constants.h"
#include "Lgm_RBF.h"
#include "Lgm_Tsyg1996.h"
#include "Lgm_Tsyg2001.h"
#include "Lgm_Tsyg2004.h"
#include "Lgm_Tsyg2007.h"
#include "Lgm_TA2016.h"
#include "gsl/gsl_errno.h"
#include "gsl/gsl_spline.h"
//#include "Lgm/Lgm_FastPowPoly.h"


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
#define     LGM_BAD_TRACE                 -3


//#define 	LGM_MAGSTEP_KMAX	16
#define 	LGM_MAGSTEP_KMAX	8
#define 	LGM_MAGSTEP_IMAX	(LGM_MAGSTEP_KMAX+1)
#define 	LGM_MAGSTEP_JMAX	(LGM_MAGSTEP_KMAX+2)
#define 	LGM_MAGSTEP_REDMAX 	1.0e-5
#define 	LGM_MAGSTEP_REDMIN 	0.7
#define 	LGM_MAGSTEP_SCLMAX 	0.1
#define 	LGM_MAGSTEP_SAFE1	0.25
#define 	LGM_MAGSTEP_SAFE2	0.70

#define 	LGM_VELSTEP_KMAX	16
#define 	LGM_VELSTEP_IMAX	(LGM_VELSTEP_KMAX+1)
#define 	LGM_VELSTEP_JMAX	(LGM_VELSTEP_KMAX+2)
#define 	LGM_VELSTEP_REDMAX 	1.0e-5
#define 	LGM_VELSTEP_REDMIN 	0.7
#define 	LGM_VELSTEP_SCLMAX 	0.1
#define 	LGM_VELSTEP_SAFE1	0.25
#define 	LGM_VELSTEP_SAFE2	0.70

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

#define LGM_CDIP        	0
#define LGM_EDIP        	1
#define LGM_IGRF        	2
#define LGM_DUNGEY      	3
#define LGM_JENSENCAIN1960      4

#define LGM_MAX_INTERP_PNTS 50000

#define LGM_RELATIVE_JUMP_METHOD 0
#define LGM_ABSOLUTE_JUMP_METHOD 1

// External Mag models
#define LGM_EXTMODEL_NULL              -1
#define LGM_EXTMODEL_T87                0
#define LGM_EXTMODEL_T89                1
#define LGM_EXTMODEL_T89c               2
#define LGM_EXTMODEL_T96                3
#define LGM_EXTMODEL_T01S               4
#define LGM_EXTMODEL_T02                5
#define LGM_EXTMODEL_TS04               6
#define LGM_EXTMODEL_TS07               7
#define LGM_EXTMODEL_OP77               8
#define LGM_EXTMODEL_SCATTERED_DATA     9
#define LGM_EXTMODEL_SCATTERED_DATA2    10
#define LGM_EXTMODEL_SCATTERED_DATA3    11
#define LGM_EXTMODEL_SCATTERED_DATA4    12
#define LGM_EXTMODEL_SCATTERED_DATA5    13
#define LGM_EXTMODEL_TU82               14
#define LGM_EXTMODEL_OP88               15
#define LGM_EXTMODEL_TA16               16

#define LGM_EXTMODEL_SCATTERED_DATA6    23


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


/*
 *  ODE solvers for FL tracing
 */
#ifndef LGM_MAGSTEP_ODE_BS
#define LGM_MAGSTEP_ODE_BS      0
#endif
#ifndef LGM_MAGSTEP_ODE_RK5
#define LGM_MAGSTEP_ODE_RK5     1
#endif

typedef struct CircularBuffer {

    int oldest_i;
    int newest_i;
    int n;
    int N;
    int nEntries;
    Lgm_Vec_RBF_Info **Buf1; // array of pointers to rbf's
    Lgm_Vec_RBF_Info **Buf2; // array of pointers to rbf's

} CircularBuffer;

typedef struct Lgm_MagModelInfo {

    Lgm_CTrans  *c;                 /* This contains all time info and a bunch more stuff */
    long int	nFunc;
    int 		(*Bfield)();
    int			SavePoints;
    double		Hmax;
    FILE	    *fp;
    double      W[6];
    double      G1, G2, G3;
    int			Kp, aKp3;
    double      fKp, Dst;
    double      V, Den, P;
    double      Bx, By, Bz;
    double      T96MOD_V[11];       /* free params for T96_MOD */

    double      Trace_s;

    double      B0, B1, B2;             // params for simplified mead field

    double      M_Dungey, dB_Dungey;    // params for simple Dungey model


    /*
     * Variable to control which internal model to use
     */
    int         InternalModel;          // Can be LGM_CDIP, LGM_EDIP or LGM_IGRF
    char        IntMagModelStr1[80];
    char        IntMagModelStr2[1024];
    char        IntMagModelStr3[1024];
    char        IntMagModelStr4[1024];

    int         ExternalModel;          // Can be from list above (e.g. LGM_EXTMODEL_T89)
    char        ExtMagModelStr1[80];
    char        ExtMagModelStr2[1024];
    char        ExtMagModelStr3[1024];
    char        ExtMagModelStr4[1024];

    /*
     * Temporary variable to hold a generic position
     */
    Lgm_Vector  Ptmp;


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
    double              MaxDiv;     // Dont subdivide the FL length with steps bigger than this.
    int                 nDivs;      // Number of divisions of FL length to try to make (actual number of points defined may be different as MAxDiv mux be respected.)
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
    double          Snorth;         //< Distance from initial point to northern footpoint along FL.
    double          Ssouth;         //< Distance from initial point to southern footpoint along FL.
    double          Stotal;         //< Distance from southern footpoint to northern footpoint along FL.
    double          Smin;           //< Distance from southern footpoint to Bmin point along FL.

    Lgm_Vector      Spherical_Footprint_Pn;      //< position of northern footpoint (at 120km)
    double          Spherical_Footprint_Sn;      //< Distance along FL from southern foorpoint in Re
    double          Spherical_Footprint_Bn;      //< Value of |B| at Sph. Footprint_Pn
    Lgm_Vector      Spherical_Footprint_Bvecn;   //< Value B of  at Sph. Footprint_Pn

    Lgm_Vector      Spherical_Footprint_Ps;      //< position of southern footpoint (at 120km)
    double          Spherical_Footprint_Ss;      //< Distance along FL from southern foorpoint in Re (i.e. this one is zero by definition)
    double          Spherical_Footprint_Bs;      //< Value of |B| at Sph. Footprint_Ps
    Lgm_Vector      Spherical_Footprint_Bvecs;   //< Value B of  at Sph. Footprint_Ps

    Lgm_Vector      Ellipsoid_Footprint_Pn;      //< position of northern footpoint (at 120km above surface of WGS84 ellipsoid)
    double          Ellipsoid_Footprint_Sn;      //< Distance along FL from southern footpoint in Re
    double          Ellipsoid_Footprint_Bn;      //< Value of |B| at Ellips. Footprint_Pn
    Lgm_Vector      Ellipsoid_Footprint_Bvecn;   //< Value B of  at Ellips. Footprint_Pn

    Lgm_Vector      Ellipsoid_Footprint_Ps;      //< position of southern footpoint (at 120km)
    double          Ellipsoid_Footprint_Ss;      //< Distance along FL from southern foorpoint in Re (i.e. this one is zero by definition)
    double          Ellipsoid_Footprint_Bs;      //< Value of |B| at Ellips. Footprint_Ps
    Lgm_Vector      Ellipsoid_Footprint_Bvecs;   //< Value B of  at Ellips. Footprint_Ps

    int             FieldLineType;  //< Field line type. (I.e., LGM_OPEN_IMF, LGM_CLOSED, LGM_OPEN_N_LOBE, LGM_OPEN_S_LOBE, LGM_INSIDE_EARTH, LGM_TARGET_HEIGHT_UNREACHABLE)




    int                 ComputeSb0;     // Flag to determine whether or not to compute d2B_ds2 and Sb0 (default is FALSE).
    double              d2B_ds2;        // second derivative of B wrt to s at smin.
    double              Sb0;            // value of Sb integral for eq. mirroring particles.
    double              Kappa, RofC;    // Curvature and Radius of Curvature at Pmin.
    int                 imin1;          // imin1 and imin2 are the indices in the
    int                 imin2;          //   array between which smin is located.


    Lgm_Vector          v1_final;       //< Place to store the final value of v1 when tryingm to trace to N foot.
    Lgm_Vector          v2_final;       //< Place to store the final value of v2 when tryingm to trace to S foot.
    Lgm_Vector          v3_final;       //< Place to store the final value of v3 when tryingm to trace to Min B.



    /*
     *  GSL defs for spline interpolation
     */
    int                 AllocedSplines;     // Flag to indicate that we have alloced splines (via gsl)
    gsl_interp_accel    *acc;               // accelerator
    gsl_interp_accel    *accPx;             // accelerator
    gsl_interp_accel    *accPy;             // accelerator
    gsl_interp_accel    *accPz;             // accelerator
    gsl_spline          *spline;            // spline object
    gsl_spline          *splinePx;          // spline object
    gsl_spline          *splinePy;          // spline object
    gsl_spline          *splinePz;          // spline object


    /*
     *  Other stuff
     */
    int                 VerbosityLevel;
    int                 UseInterpRoutines; // whether to use fast I and Sb routines.


    long int    Lgm_nMagEvals;          // records number of Bfield evals between resets
    int         Lgm_MagStep_Integrator; // ODE solver to use ( LGM_MAGSTEP_ODE_BS or LGM_MAGSTEP_ODE_RK5)

    /*
     *  Vars for Bulirsch-Stoer ODE solver. Some of these variables are needed
     *  to make Lgm_MagStep() reentrant/thread-safe.  They basically used to be
     *  static declarations within Lgm_MagStep()
     */
    double      Lgm_MagStep_BS_eps_old;
    int	        Lgm_MagStep_BS_FirstTimeThrough;
    int         Lgm_MagStep_BS_kmax;
    int         Lgm_MagStep_BS_kopt;
    double      Lgm_MagStep_BS_snew;
    double      Lgm_MagStep_BS_A[LGM_MAGSTEP_JMAX+1];
    double      Lgm_MagStep_BS_alpha[LGM_MAGSTEP_IMAX+1][LGM_MAGSTEP_IMAX+1];
    double      Lgm_MagStep_BS_d[LGM_MAGSTEP_JMAX][LGM_MAGSTEP_JMAX];
    double      Lgm_MagStep_BS_x[LGM_MAGSTEP_JMAX];
    Lgm_Vector  Lgm_MagStep_BS_u[LGM_MAGSTEP_JMAX];
    double      Lgm_MagStep_BS_Eps;        // Eps parameter used in BS method. Influences speed greatly.


    /*
     * For new BS ODE solver...
     */
    //double      Lgm_MagStep_BS_coeff[LGM_MAGSTEP_IMAX+1][LGM_MAGSTEP_IMAX+1];
    double      Lgm_MagStep_BS_coeff[LGM_MAGSTEP_IMAX][LGM_MAGSTEP_IMAX];
    int         Lgm_MagStep_BS_first_step;
    int         Lgm_MagStep_BS_last_step;
    int         Lgm_MagStep_BS_reject;
    int         Lgm_MagStep_BS_prev_reject;
    int         Lgm_MagStep_BS_k_targ;
    double      Lgm_MagStep_BS_atol;
    double      Lgm_MagStep_BS_rtol;
    //double      Lgm_MagStep_BS_cost[LGM_MAGSTEP_KMAX+1];
    double      Lgm_MagStep_BS_cost[LGM_MAGSTEP_IMAX];



    /*
     *  These variables are needed to make Lgm_MagStep2() reentrant/thread-safe.
     *  They basically used to be static declarations within Lgm_MagStep()
     */
    int	        Lgm_MagStep_RK5_FirstTimeThrough;
    double      Lgm_MagStep_RK5_snew;
    double      Lgm_MagStep_RK5_MaxCount;
    double      Lgm_MagStep_RK5_Safety;
    double      Lgm_MagStep_RK5_pGrow;
    double      Lgm_MagStep_RK5_pShrnk;
    double      Lgm_MagStep_RK5_ErrCon;
    double      Lgm_MagStep_RK5_Eps;        // Eps parameter used in RK5 method. Influences speed greatly.

    /*
     *  These variables are needed to make Lgm_MagStep2() reentrant/thread-safe.
     *  They basically used to be static declarations within Lgm_MagStep2()
     */
    long int    Lgm_nVelEvals; // records number of Bfield evals between resets
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

    double      Lgm_I_Integrator_epsrel;        // Quadpack epsrel tolerance for I_integrator
    double      Lgm_I_Integrator_epsabs;        // Quadpack epsabs tolerance for I_integrator


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
    double      Lgm_Sb_Integrator_epsrel;       // Quadpack epsrel tolerance for Sb_integrator
    double      Lgm_Sb_Integrator_epsabs;       // Quadpack epsabs tolerance for Sb_integrator



    /*
     * Variables to control FluxTubeVolume integration
     */
    int         Lgm_n_V_integrand_Calls;
    int         Lgm_V_Integrator;

    double      Lgm_V_Integrator_epsrel;        // Quadpack epsrel tolerance for V_integrator
    double      Lgm_V_Integrator_epsabs;        // Quadpack epsabs tolerance for V_integrator


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
    double		Lgm_TraceToEarth_Tol;    // tolerance for deciding when you have converged to a footpoint. (Even for low Tol, the points obtained should still be on the FL to high precision.)
    double		Lgm_TraceToBmin_Tol;     // tolerance for deciding when you have converged to the Bmin point along a FL. (Even for low Tol, the points obtained should still be on the FL to high precision.)
    double		Lgm_TraceLine_Tol;       // tolerance for magstep in Lgm_TraceLine() routines.



    /*
     * Variables for defining Octree stuff
     */
    Lgm_Octree     *Octree;
    int             Octree_Alloced;     // TRUE if memory is alloced for the octree, FALSE otherwise.
    int             Octree_kNN_InterpMethod;
    int             Octree_kNN_k;
    double          Octree_kNN_MaxDist2;
    Lgm_OctreeData *Octree_kNN;
    int             Octree_kNN_Alloced; // number of elements allocated. (0 if unallocated).

    /*
     * Variables for defining KdTree stuff
     */
    Lgm_KdTree     *KdTree;
    int             KdTree_Alloced;     // TRUE if memory is alloced for the KdTree, FALSE otherwise.
    int             KdTree_kNN_InterpMethod;
    int             KdTree_kNN_k;
    double          KdTree_kNN_MaxDist2;
    Lgm_KdTreeData *KdTree_kNN;
    int             KdTree_kNN_Alloced; // number of elements allocated. (0 if unallocated).
    int             KdTreeCopy;         // If set, then we assume the pointer
                                        // to the tree data is a copy and we
                                        // should not free it when calling
                                        // Lgm_FreeMagInfo().


    /*
     *  hash table, etc.  used in Lgm_B_FromScatteredData*()
     */
    Lgm_DFI_RBF_Info   *rbf_ht;             // hash table (uthash)
    double             dfi_rbf_ht_size;     // hash table size in MB
    double             dfi_rbf_ht_maxsize;  // hash table max size in MB
    CircularBuffer     RBF_DFI_CB;

    Lgm_Vec_RBF_Info   *vec_rbf_ht;         // hash table (uthash)
    double             vec_rbf_ht_size;     // hash table size in MB
    double             vec_rbf_ht_maxsize;  // hash table max size in MB
    CircularBuffer     RBF_CB;

    Lgm_Vec_RBF_Info   *vec_rbf_e_ht;       // hash table (uthash)
    double             vec_rbf_e_ht_size;   // hash table size in MB
    double             vec_rbf_e_ht_maxsize;// hash table max size in MB
    CircularBuffer     RBF_E_CB;

    int                 rbf_ht_alloced;     // Flag to indicate whether or not rbf_ht is allocated with data.
    long int            RBF_nHashFinds;     // Number of HASH_FIND()'s performed.
    long int            RBF_nHashAdds;      // Number of HASH_ADD_KEYPTR()'s performed.
    Lgm_Vector          RBF_dBdx;           // deriv of B-vec wrt x computed using RBF.
    Lgm_Vector          RBF_dBdy;           // deriv of B-vec wrt y computed using RBF.
    Lgm_Vector          RBF_dBdz;           // deriv of B-vec wrt z computed using RBF.
    Lgm_Vector          RBF_E;              // Electric field computed using RBF.
    Lgm_Vector          RBF_dEdx;           // deriv of B-vec wrt x computed using RBF.
    Lgm_Vector          RBF_dEdy;           // deriv of B-vec wrt y computed using RBF.
    Lgm_Vector          RBF_dEdz;           // deriv of B-vec wrt z computed using RBF.
    int                 RBF_CompGradAndCurl;
    Lgm_Vector          RBF_Grad_B;         // Grad_B
    Lgm_Vector          RBF_Curl_B;         // Curl_B
    Lgm_Vector          RBF_Curl_E;         // Curl_B
    int                 RBF_Type;           // Type of RBF to use
    int                 RBF_DoPoly;         // Flag to simultaneously fit a linear polynomail ( i.e. Sum_ijk{ a_ijl x^i y^j z^k }) as well.
    double              RBF_Eps;            // Eps value to use


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


    /*
     *  Info structure for T96
     */
    LgmTsyg1996_Info    T96_Info;

    /*
     *  Info structure for T01S
     */
    LgmTsyg2001_Info    T01_Info;

    /*
     *  Info structure for TS04
     */
    LgmTsyg2004_Info    TS04_Info;


    /*
     *  Info structure for TS07
     */
    LgmTsyg2007_Info    TS07_Info;

    /*
     *  Info structure for TA16
     */
    LgmTA16_Info    TA16_Info;


    /*
     *  hash table used in Lgm_FastPow()
     */
    //Lgm_FastPow   *fastpow_ht;          // hash table (uthash)
    //int            fastpow_ht_alloced;  // Flag to indicate whether or not fastpow_ht is allocated with data.
    //long int       FP_nHashFinds;       // Number of HASH_FIND()'s performed.
    //long int       FP_nHashAdds;        // Number of HASH_ADD_KEYPTR()'s performed.

    /*
     * Transformation matrix from GSM to PQB and (reverse) 
     */
    double Agsm_to_pqb[3][3];
    double Apqb_to_gsm[3][3];
    int    Agsm_to_pqb_set;


    void        *Data;


} Lgm_MagModelInfo;
void Lgm_GSM_TO_PQB( Lgm_Vector *u_gsm, Lgm_Vector *u_pqb, Lgm_MagModelInfo *m );
void Lgm_Set_GSM_TO_PQB( Lgm_Vector *Position, Lgm_MagModelInfo *m );
void Lgm_PQB_TO_GSM( Lgm_Vector *u_pqb, Lgm_Vector *u_gsm, Lgm_MagModelInfo *m );



typedef struct BrentFuncInfoP {

    Lgm_Vector          u_scale;
    double              Htry, sgn;
    int                 reset;
    Lgm_MagModelInfo    *Info;
    double              Val;
    double 		        (*func)( Lgm_Vector *P, double Val, Lgm_MagModelInfo *Info );

} BrentFuncInfoP;
int Lgm_BrentP(double Sa, double Sb, double Sc, double Bb, Lgm_Vector Pa, Lgm_Vector Pb, Lgm_Vector Pc, BrentFuncInfoP *fInfo, double tol, double *Smin, double *Bmin, Lgm_Vector *Pmin );
int Lgm_zBrentP(double S1, double S2, double F1, double F2, Lgm_Vector P1, Lgm_Vector P2, BrentFuncInfoP *fInfo, double tol, double *Sz, double *Fz, Lgm_Vector *Pz );

typedef struct BrentFuncInfo {

    double              Val;
    void                *Info;
    double 		        (*func)( double x, double Val, void *Info );

} BrentFuncInfo;
int Lgm_Brent(double xa, double xb, double xc, BrentFuncInfo *fInfo, double tol, double *xmin, double *fmin );
int Lgm_zBrent(double x1, double x2, double f1, double f2, BrentFuncInfo *fInfo, double tol, double *x, double *f );


Lgm_MagModelInfo *Lgm_InitMagInfo( );
void Lgm_InitMagInfoDefaults( Lgm_MagModelInfo  *MagInfo );

void Lgm_FreeMagInfo_children( Lgm_MagModelInfo  *Info );
void Lgm_FreeMagInfo( Lgm_MagModelInfo  *Info );
Lgm_MagModelInfo *Lgm_CopyMagInfo( Lgm_MagModelInfo *s );

int  Lgm_Trace( Lgm_Vector *u, Lgm_Vector *v1, Lgm_Vector *v2, Lgm_Vector *v3, double Height, double TOL1, double TOL2, Lgm_MagModelInfo *Info );
int  Lgm_TraceToMinBSurf( Lgm_Vector *, Lgm_Vector *, double, double, Lgm_MagModelInfo * );
int  Lgm_TraceToSMEquat(  Lgm_Vector *, Lgm_Vector *, double, Lgm_MagModelInfo * );
int  Lgm_TraceToEarth(  Lgm_Vector *, Lgm_Vector *, double, double, double, Lgm_MagModelInfo * );
int  Lgm_TraceToSphericalEarth(  Lgm_Vector *, Lgm_Vector *, double, double, double, Lgm_MagModelInfo * );
int  Lgm_TraceLine(  Lgm_Vector *, Lgm_Vector *, double, double, double, int, Lgm_MagModelInfo * );
int  Lgm_TraceLine2(  Lgm_Vector *, Lgm_Vector *, double, double, double, double, int, Lgm_MagModelInfo * );
int  Lgm_TraceLine3( Lgm_Vector *u, double S, int N, double sgn, double tol, int AddBminPoint, Lgm_MagModelInfo *Info );
int  Lgm_TraceLine4( Lgm_Vector *Pm_s, Lgm_Vector *Pm_n, double dSa, double dSb, int N, int AddBminPoint, Lgm_MagModelInfo *Info );
int   Lgm_TraceToYZPlane( Lgm_Vector *u, Lgm_Vector *v, double Xtarget, double sgn_in, double tol, Lgm_MagModelInfo *Info );

int  ReplaceFirstPoint( double s, double B, Lgm_Vector *P, Lgm_MagModelInfo *Info );
int  ReplaceLastPoint( double s, double B, Lgm_Vector *P, Lgm_MagModelInfo *Info );
void AddNewPoint( double s, double B, Lgm_Vector *P, Lgm_MagModelInfo *Info );
int  InitSpline( Lgm_MagModelInfo *Info );
int  FreeSpline( Lgm_MagModelInfo *Info );
int  Lgm_TraceToMinRdotB( Lgm_Vector *, Lgm_Vector *, double, Lgm_MagModelInfo * );
int  Lgm_TraceIDL( int, void *argv[] );
int  Lgm_TraceToMirrorPoint( Lgm_Vector *u, Lgm_Vector *v, double *Sm, double Bm, double sgn, double tol, Lgm_MagModelInfo *Info );


int  Lgm_MagStep( Lgm_Vector *, Lgm_Vector *, double, double *, double *, double, double *, int *,
              int (*Mag)(Lgm_Vector *, Lgm_Vector *, Lgm_MagModelInfo *), Lgm_MagModelInfo * );
int  Lgm_MagStep_BS( Lgm_Vector *, Lgm_Vector *, double, double *, double *, double, double, double *, int *,
              int (*Mag)(Lgm_Vector *, Lgm_Vector *, Lgm_MagModelInfo *), Lgm_MagModelInfo * );
int  Lgm_MagStep_RK5( Lgm_Vector *, Lgm_Vector *, double, double *, double *, double, double, double *, int *,
              int (*Mag)(Lgm_Vector *, Lgm_Vector *, Lgm_MagModelInfo *), Lgm_MagModelInfo * );

/*
 * For Bulirsch-Stoer FL tracer
 */
int Lgm_ModMid( Lgm_Vector *, Lgm_Vector *, Lgm_Vector *, double, int, double,
	     int (*Mag)(Lgm_Vector *, Lgm_Vector *, Lgm_MagModelInfo *), Lgm_MagModelInfo * );
void Lgm_RatFunExt( int, double, Lgm_Vector *, Lgm_Vector *, Lgm_Vector *, Lgm_MagModelInfo * );
void Lgm_PolFunExt( int, double, Lgm_Vector *, Lgm_Vector *, Lgm_Vector *, Lgm_MagModelInfo * );




//int Lgm_ModMid2( Lgm_Vector *, Lgm_Vector *, Lgm_Vector *, double, int, double,
//	     int (*Velocity)(Lgm_Vector *, Lgm_Vector *, Lgm_MagModelInfo *), Lgm_MagModelInfo * );
//int  Lgm_VelStep( Lgm_Vector *, Lgm_Vector *, double, double *, double *, double, double, double *, int *,
//              int (*Velocity)(Lgm_Vector *, Lgm_Vector *, Lgm_MagModelInfo *), Lgm_MagModelInfo * );


/*
 * for RK5 FL tracer.
 */
int Lgm_RKCK( Lgm_Vector *u0, Lgm_Vector *b0, Lgm_Vector *v, double h, double sgn, double *yerr,
        int (*Mag)(Lgm_Vector *, Lgm_Vector *, Lgm_MagModelInfo *), Lgm_MagModelInfo *Info );



/*
 *  B_internal
 *
 *  Function Prototypes for internal models
 */
int Lgm_B_igrf(Lgm_Vector *, Lgm_Vector *, Lgm_MagModelInfo *);
int Lgm_B_cdip(Lgm_Vector *, Lgm_Vector *, Lgm_MagModelInfo *);
int Lgm_B_edip(Lgm_Vector *, Lgm_Vector *, Lgm_MagModelInfo *);
int Lgm_B_JensenCain1960(Lgm_Vector *, Lgm_Vector *, Lgm_MagModelInfo *);



/*
 *
 *  Function Prototypes for Olsen Pfitzer 1977  Model
 *
 */
int Lgm_B_OP77( Lgm_Vector *v, Lgm_Vector *B, Lgm_MagModelInfo *Info );
void OlsenPfitzerStatic( double XX[], double BF[], double TILT, Lgm_MagModelInfo *m );


/*
 *
 *  Function Prototypes for Olsen Pfitzer 1988  Model
 *
 */
int     Lgm_B_OP88( Lgm_Vector *v, Lgm_Vector *B, Lgm_MagModelInfo *Info );
void    Lgm_OP88_BDYN( double DEN, double VEL, double DST, double X, double Y, double Z, double *BX, double *BY, double *BZ );
void    Lgm_OP88_BDYNAM( double *XX, double *BB, double SOFFD, double SRING, double STAIL );
void    Lgm_OP88_BFMAGP( double *XX, double *BB );
void    Lgm_OP88_BFTAIL( double *XX,  double *BB );
void    Lgm_OP88_BFRING( double *XX,  double *BB );
double  Lgm_OP88_RINGST( double SOFFD,  double DST );
double  Lgm_OP88_STDOFF( double VEL,  double DEN );



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
 *  TU82
 *
 *  Function Prototypes for TU82 model
 */
int Lgm_Brc_TU82( Lgm_Vector *, Lgm_Vector *, Lgm_MagModelInfo * );
int Lgm_Bt_TU82( Lgm_Vector *, Lgm_Vector *, Lgm_MagModelInfo * );
int Lgm_Bmp_TU82( Lgm_Vector *, Lgm_Vector *, Lgm_MagModelInfo * );
int Lgm_B_TU82( Lgm_Vector *, Lgm_Vector *, Lgm_MagModelInfo * );


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
 *  T89c
 *
 *  Function Prototypes for T89c model
 */
void T89c( int IOPT, double *PARMOD, double PS, double X, double Y, double Z, double *BX, double *BY, double *BZ, Lgm_MagModelInfo *Info );
void T89c_Field( int ID, double *A, double *XI, double *F, Lgm_MagModelInfo *Info );
int Lgm_B_T89c( Lgm_Vector *, Lgm_Vector *, Lgm_MagModelInfo * );


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
 *  T96 
 *
 *  Function Prototypes for T96 model
 */
int  Lgm_B_T96( Lgm_Vector *v, Lgm_Vector *B, Lgm_MagModelInfo *Info );
void Tsyg_T96( int IOPT, double *PARMOD, double PS, double SINPS, double COSPS, double X, double Y, double Z, double *BX, double *BY, double *BZ, LgmTsyg1996_Info *tInfo );

/*
 *  TS04 Optimized
 *
 *  Function Prototypes for TS04 model
 */
int  Lgm_B_TS04( Lgm_Vector *v, Lgm_Vector *B, Lgm_MagModelInfo *Info );
void Tsyg_TS04( int IOPT, double *PARMOD, double PS, double SINPS, double COSPS, double X, double Y, double Z, double *BX, double *BY, double *BZ, LgmTsyg2004_Info *tInfo );

/*
 *  TS07 Optimized
 *
 *  Function Prototypes for TS07 model
 */
void Lgm_DeAllocate_TS07( LgmTsyg2007_Info *tInfo );
int Lgm_B_TS07( Lgm_Vector *v, Lgm_Vector *B, Lgm_MagModelInfo *Info );
void Tsyg_TS07( int IOPT, double *PARMOD, double PS, double SINPS, double COSPS, double X, double Y, double Z, double *BX, double *BY, double *BZ, LgmTsyg2007_Info *tInfo );

/*
 * TA16
 *
 * Function Prototypes for TA16 model
 */
int Lgm_B_TA16( Lgm_Vector *v, Lgm_Vector *B, Lgm_MagModelInfo *Info );

/*
 *  T01S
 *
 *  Function Prototypes for T01S model
 */
double mypow( double, double );
void   Lgm_Init_T01S( LgmTsyg2001_Info *t );
int    Lgm_B_T01S( Lgm_Vector *v, Lgm_Vector *B, Lgm_MagModelInfo *Info );
void   Tsyg_T01S( int IOPT, double *PARMOD, double PS, double SINPS, double COSPS, double X, double Y, double Z, double *BX, double *BY, double *BZ, LgmTsyg2001_Info *tInfo );
int    Lgm_B_T02( Lgm_Vector *v, Lgm_Vector *B, Lgm_MagModelInfo *Info );
void   Tsyg_T02( int IOPT, double *PARMOD, double PS, double SINPS, double COSPS, double X, double Y, double Z, double *BX, double *BY, double *BZ, LgmTsyg2001_Info *tInfo );

/*
 *  Computing B from scattered data -- (e.g. an irregular mesh)
 *
 */

int Lgm_B_FromScatteredData( Lgm_Vector *v, Lgm_Vector *B, Lgm_MagModelInfo *Info );
int Lgm_B_FromScatteredData2( Lgm_Vector *v, Lgm_Vector *B, Lgm_MagModelInfo *Info );
int Lgm_B_FromScatteredData3( Lgm_Vector *v, Lgm_Vector *B, Lgm_MagModelInfo *Info );
int Lgm_B_FromScatteredData4( Lgm_Vector *v, Lgm_Vector *B, Lgm_MagModelInfo *Info );
int Lgm_B_FromScatteredData5( Lgm_Vector *v, Lgm_Vector *B, Lgm_MagModelInfo *Info );
int Lgm_B_FromScatteredData6( Lgm_Vector *v, Lgm_Vector *B, Lgm_MagModelInfo *Info );
void Lgm_B_FromScatteredData_SetUp( Lgm_MagModelInfo *Info );
void Lgm_B_FromScatteredData_TearDown( Lgm_MagModelInfo *Info );
void Lgm_B_FromScatteredData4_TearDown( Lgm_MagModelInfo *Info ); // unify the structs to avoid having this

void Lgm_B_FromScatteredData5_SetUp( Lgm_MagModelInfo *Info );
void Lgm_B_FromScatteredData5_TearDown( Lgm_MagModelInfo *Info ); // I dont like this proliferation of routines here..
void Lgm_B_FromScatteredData6_SetUp( Lgm_MagModelInfo *Info );
void Lgm_B_FromScatteredData6_TearDown( Lgm_MagModelInfo *Info ); // I dont like this proliferation of routines here..


/*
 * Simplified Mead field.
 */
int Lgm_SimplifiedMead(Lgm_Vector *v, Lgm_Vector *B, Lgm_MagModelInfo *Info);

/*
 * Simple Dungey-Type open model. I.e. centered dipole plus constant dB.
 * This model uses Inf0->M_Dungey and Info->dB_Dungey
 */
int Lgm_B_Dungey(Lgm_Vector *v, Lgm_Vector *B, Lgm_MagModelInfo *Info);


/*
 * routines/functions for field integrals, invariants, etc.
 */
double      FluxTubeVolume( Lgm_MagModelInfo *fInfo );
double      V_integrand( double s, _qpInfo *qpInfo );
double      Iinv( Lgm_MagModelInfo *fInfo );
double      I_integrand( double s, _qpInfo *qpInfo );
double      Iinv_interped( Lgm_MagModelInfo *fInfo );
double      I_integrand_interped( double s, _qpInfo *qpInfo );
double      SbIntegral( Lgm_MagModelInfo *fInfo );
double      Sb_integrand( double s, _qpInfo *qpInfo );
double      SbIntegral_interped( Lgm_MagModelInfo *fInfo );
double      SbIntegral_interped2( Lgm_MagModelInfo *fInfo, double a , double b );
double      Sb_integrand_interped( double s, _qpInfo *qpInfo );
void        ratint( double *xa, double *ya, int n, double x, double *y, double *dy );
void        polint(double *xa, double *ya, int n, double x, double *y, double *dy);
void        Interp( double xa[],  double ya[], long int n, double x, double *y);
void        Interp2( double xa[],  double ya[], long int n, double x, double *y);
double      LFromIBmM_Hilton( double I, double Bm, double M );
double      IFromLBmM_Hilton( double L, double Bm, double M );
double      LFromIBmM_McIlwain( double I, double Bm, double M );
double      IFromLBmM_McIlwain( double L, double Bm, double M );
double      Lgm_McIlwain_L( long int Date, double UTC, Lgm_Vector *u, double Alpha, int Type, double *I, double *Bm, double *M, Lgm_MagModelInfo *mInfo );


double      BofS( double s, Lgm_MagModelInfo *Info );
int         SofBm( double Bm, double *ss, double *sn, Lgm_MagModelInfo *Info );
double      Lgm_AlphaOfK( double K, Lgm_MagModelInfo *Info );
double      Lgm_KofAlpha( double Alpha, Lgm_MagModelInfo *Info );
int         Lgm_Setup_AlphaOfK( Lgm_DateTime *d, Lgm_Vector *u, Lgm_MagModelInfo *m );
void        Lgm_TearDown_AlphaOfK( Lgm_MagModelInfo *m );
int         Lgm_Grad_I( Lgm_Vector *vin, Lgm_Vector *GradI, Lgm_MagModelInfo *Info );
//int         ComputeVcg( Lgm_Vector *vin, Lgm_Vector *Vcg, Lgm_LstarInfo *LstarInfo );


/*
 * Getter/Setter functions. May not want to do this...
 *  Just testing for now..
 */
void Lgm_MagModelInfo_Set_Psw( double Psw, Lgm_MagModelInfo *m );
void Lgm_MagModelInfo_Set_Kp( double Kp, Lgm_MagModelInfo *m );
void Lgm_MagModelInfo_Set_Octree( Lgm_Octree *Octree, int k, Lgm_MagModelInfo *m );
void Lgm_Set_Octree_kNN_InterpMethod( Lgm_MagModelInfo *m, int Method );
void Lgm_Set_Octree_kNN_k( Lgm_MagModelInfo *m, int k );
void Lgm_Set_Octree_kNN_MaxDist2( Lgm_MagModelInfo *m, double MaxDist2 );
void Lgm_Set_KdTree( Lgm_KdTree *KdTree, int k, double d2, Lgm_MagModelInfo *m );
void Lgm_Set_KdTree_kNN_k( Lgm_MagModelInfo *m, int k );
void Lgm_Set_KdTree_kNN_MaxDist2( Lgm_MagModelInfo *m, double MaxDist2 );
void Lgm_Set_Open_Limits( Lgm_MagModelInfo *m, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax );
void Lgm_Set_LossConeHeight( Lgm_MagModelInfo *m, double LossConeHeight );
void Lgm_MagModelInfo_Set_MagModel( int InternalModel, int ExternalModel, Lgm_MagModelInfo *m );
void Lgm_Set_MagModel( int InternalModel, int ExternalModel, Lgm_MagModelInfo *m );
void Lgm_Get_MagModel( int *InternalModel, int *ExternalModel, Lgm_MagModelInfo *m );
void Lgm_Get_IntMagModelStrings( char **s1, char **s2, char **s3, char **s4, Lgm_MagModelInfo *m );
void Lgm_Get_ExtMagModelStrings( char **s1, char **s2, char **s3, char **s4, Lgm_MagModelInfo *m );

void Lgm_B_FromScatteredData_SetRbf( Lgm_MagModelInfo *Info, double eps, int RbfType );

/*
 * Added for Python wrapping, the function pointer is hard to impossible to
 * deal with setting.
 *
 *   (Note: The routine Lgm_MagModelInfo_Set_MagModel() sets internal and
 *   external models in a single call for all models. MGH - 20130227)
 */
void Lgm_Set_Lgm_B_cdip(Lgm_MagModelInfo *MagInfo);
void Lgm_Set_Lgm_B_edip(Lgm_MagModelInfo *MagInfo);
void Lgm_Set_Lgm_B_igrf(Lgm_MagModelInfo *MagInfo);
void Lgm_Set_Lgm_B_Dungey(Lgm_MagModelInfo *MagInfo);
void Lgm_Set_Lgm_B_T01S(Lgm_MagModelInfo *MagInfo);
void Lgm_Set_Lgm_B_TS04(Lgm_MagModelInfo *MagInfo);
void Lgm_Set_Lgm_B_TA16(Lgm_MagModelInfo *MagInfo);
void Lgm_Set_Lgm_B_TS04_opt(Lgm_MagModelInfo *MagInfo);
void Lgm_Set_Lgm_B_T89(Lgm_MagModelInfo *MagInfo);
void Lgm_Set_Lgm_B_T89c(Lgm_MagModelInfo *MagInfo);
void Lgm_Set_Lgm_B_T96(Lgm_MagModelInfo *MagInfo);
void Lgm_Set_Lgm_B_OP77(Lgm_MagModelInfo *MagInfo);
void Lgm_Set_Lgm_B_cdip_InternalModel(Lgm_MagModelInfo *MagInfo);
void Lgm_Set_Lgm_B_edip_InternalModel(Lgm_MagModelInfo *MagInfo);
void Lgm_Set_Lgm_B_IGRF_InternalModel(Lgm_MagModelInfo *MagInfo);


/*
 * coordinate conversions to/from corrected Geomag.
 */
int  Lgm_GSM_TO_CBM( Lgm_Vector *u, double *CgmLat, double *CgmLon, double *CgmRad, Lgm_MagModelInfo *m );
int  Lgm_CBM_TO_GSM( double CgmLat, double CgmLon, double CgmRad, Lgm_Vector *u, Lgm_MagModelInfo *m );
int  Lgm_GEOD_TO_CGM( double geoLat, double geoLon, double geoAlt, double *CgmLat, double *CgmLon, double *CgmRad, Lgm_MagModelInfo *m );
int  Lgm_CGM_TO_GEOD( double CgmLat, double CgmLon, double CgmRadi, double *geoLat, double *geoLon, double *geoAlt, Lgm_MagModelInfo *m );



/*
 * Various B-related routines
 */
void    Lgm_GradB( Lgm_Vector *u0, Lgm_Vector *GradB, int DerivScheme, double h, Lgm_MagModelInfo *m );
void    Lgm_Grad_Divb( Lgm_Vector *u0, Lgm_Vector *Grad_Divb, int DerivScheme, double h, Lgm_MagModelInfo *m );
void    Lgm_GradB2( Lgm_Vector *u0, Lgm_Vector *GradB, Lgm_Vector *GradB_para, Lgm_Vector *GradB_perp, int DerivScheme, double h, Lgm_MagModelInfo *m );
void    Lgm_GradBvec( Lgm_Vector *u0, double GradBvec[3][3], int DerivScheme, double h, Lgm_MagModelInfo *m );
void    Lgm_GradBvec2( int j, Lgm_Vector *u0, double GradBvec[4][4], int DerivScheme, double h, int n, double dt, Lgm_MagModelInfo **m );
void    Lgm_CurlB( Lgm_Vector *u0, Lgm_Vector *CurlB, int DerivScheme, double h, Lgm_MagModelInfo *m );
void    Lgm_Curlb( Lgm_Vector *u0, Lgm_Vector *Curlb, int DerivScheme, double h, Lgm_MagModelInfo *m );
void    Lgm_CurlB2( Lgm_Vector *u0, Lgm_Vector *CurlB, Lgm_Vector *CurlB_para, Lgm_Vector *CurlB_perp, int DerivScheme, double h, Lgm_MagModelInfo *m );
void    Lgm_B_Cross_GradB_Over_B( Lgm_Vector *u0, Lgm_Vector *A, int DerivScheme, double h, Lgm_MagModelInfo *m );
void    Lgm_DivB( Lgm_Vector *u0, double *DivB, int DerivScheme, double h, Lgm_MagModelInfo *m );
void    Lgm_Divb( Lgm_Vector *u0, double *Divb, int DerivScheme, double h, Lgm_MagModelInfo *m );
int     Lgm_GradAndCurvDriftVel( Lgm_Vector *u0, Lgm_Vector *Vel, Lgm_MagModelInfo *m );
void    Lgm_dBdcomp( Lgm_Vector *u0, int comp, double *dBdcomp, int DerivScheme, double h, Lgm_MagModelInfo *m );
void    Lgm_dBdx( Lgm_Vector *u0, double *dBx, int DerivScheme, double h, Lgm_MagModelInfo *m );
void    Lgm_dBdy( Lgm_Vector *u0, double *dBy, int DerivScheme, double h, Lgm_MagModelInfo *m );
void    Lgm_dBdz( Lgm_Vector *u0, double *dBz, int DerivScheme, double h, Lgm_MagModelInfo *m );
void    Lgm_LaplacianB( Lgm_Vector *u0, double *LaplacianB, int DerivScheme, double h, Lgm_MagModelInfo *m );
void    Lgm_dDivbdcomp( Lgm_Vector *u0, int comp, double *dDivbdcomp, int DerivScheme, double h, Lgm_MagModelInfo *m );
void    Lgm_Curvature( Lgm_Vector *q, Lgm_Vector *Rc, int DerivScheme, double h, Lgm_MagModelInfo *m );
void    Lgm_Curvature2( Lgm_Vector *q, Lgm_Vector *Rc, int DerivScheme, double h, Lgm_MagModelInfo *m );
int     Lgm_Vdrift_GradB_Curv( Lgm_Vector *u0, Lgm_Vector *v_full, double Mass, double Charge, Lgm_Vector *v_GradB, Lgm_Vector *v_CurvB, Lgm_Vector *Rc_vec, Lgm_MagModelInfo *m );


void quicksort_uli( unsigned long n, unsigned long *arr );


/*
 * Rotuines for computing the first invariant up to second order
 */
double  Lgm_DoubleDot( double A[3][3], double B[3][3] );
void    Lgm_Set_Particle_Frames( Lgm_Vector *q, double Mass, double RestEnergy, double Energy, double Phi, double Delta, 
                Lgm_Vector *a, Lgm_Vector *c, Lgm_Vector *b, Lgm_Vector *v, Lgm_MagModelInfo *mInfo );
void    Lgm_Grad_Cb( Lgm_Vector *u0, double Grad_Cb[3][3], int DerivScheme, double h, Lgm_MagModelInfo *m );
void    Lgm_Grad_b( Lgm_Vector *u0, double Grad_b[3][3], int DerivScheme, double h, Lgm_MagModelInfo *m );
void    Lgm_Grad_GradB( Lgm_Vector *u0, double Grad_GradB[3][3], int DerivScheme, double h, Lgm_MagModelInfo *m );
void    Lgm_dbdotvdcomp( Lgm_Vector *u0, Lgm_Vector *v, int comp, double *dbdotvdcomp, int DerivScheme, double h, Lgm_MagModelInfo *m );
void    Lgm_dbdotvdx( Lgm_Vector *u0, Lgm_Vector *v, double *dbdotvdx, int DerivScheme, double h, Lgm_MagModelInfo *m );
void    Lgm_dbdotvdy( Lgm_Vector *u0, Lgm_Vector *v, double *dbdotvdy, int DerivScheme, double h, Lgm_MagModelInfo *m );
void    Lgm_dbdotvdz( Lgm_Vector *u0, Lgm_Vector *v, double *dbdotvdz, int DerivScheme, double h, Lgm_MagModelInfo *m );
void    Lgm_Laplacian_bdotv( Lgm_Vector *u0, Lgm_Vector *v, double *Laplacian_bdotv, int DerivScheme, double h, Lgm_MagModelInfo *m );
void    Lgm_Grad_bdotv( Lgm_Vector *u0, Lgm_Vector *v, Lgm_Vector *Grad_bdotv, int DerivScheme, double h, Lgm_MagModelInfo *m );
void    Lgm_Grad_Grad_bdotv( Lgm_Vector *u0, Lgm_Vector *v, double Grad_Grad_bdotv[3][3], int DerivScheme, double h, Lgm_MagModelInfo *m );
void    Lgm_Dyad( Lgm_Vector *a, Lgm_Vector *b, double r[3][3] );
void    Lgm_VectorDotTensor( Lgm_Vector *v, double T[3][3], Lgm_Vector *Result );
void    Lgm_TensorDotVector( double T[3][3], Lgm_Vector *v, Lgm_Vector *Result );
void    Lgm_TensorDotTensor( double A[3][3], double B[3][3], double R[3][3] );
double  Lgm_TraceTensor( double T[3][3] );
void    Lgm_TensorTranspose( double T[3][3], double T_transpose[3][3] );
double  Lgm_DoubleDot( double A[3][3], double B[3][3] );
double  Lgm_Burby( Lgm_Vector *q, Lgm_Vector *v, double Gamma, double Mass, double Charge, double *mu0_out, double *mu1_out, double *mu2_out, Lgm_MagModelInfo *mInfo );


/*
 * Misc Routines
 */
double  Lgm_Bradial( double MLT, double mlat, Lgm_MagModelInfo *m );
int     Lgm_FindDipEquator( double MLT, double *mlat, Lgm_MagModelInfo *m );


#endif
