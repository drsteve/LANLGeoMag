#ifndef LGM_RBF
#define LGM_RBF

#include "Lgm_Vec.h"
#include "Lgm_RBF.h"
#include "Lgm_DynamicMemory.h"
#include "uthash.h"
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

#define LGM_CHOLESKY_DECOMP 0
#define LGM_PLU_DECOMP      1    
#define LGM_SVD             2    

#define LGM_RBF_GAUSSIAN         0
#define LGM_RBF_MULTIQUADRIC     1
#define LGM_RBF_INV_MULTIQUADRIC 2
#define LGM_RBF_WENDLAND30       3
#define LGM_RBF_WENDLAND31       4
#define LGM_RBF_WENDLAND32       5
#define LGM_RBF_WENDLAND33       6

/*
 * For  doing div-free interp.
 */
typedef struct _Lgm_DFI_RBF_Info {

    unsigned long int   *LookUpKey;  // keyword for hashing (an array of id numbers)
    int                 RadialBasisFunction;
    int                 n;
    int                 n3;
    double              Bx0;
    double              By0;
    double              Bz0;
    Lgm_Vector          *v;
    Lgm_Vector          *c;
    double              eps;
//NOTE
    UT_hash_handle      hh; // Make structure hashable via uthash

} Lgm_DFI_RBF_Info;




/*
 * For doing vec interp with each comp treated independantly.
 */
typedef struct _Lgm_Vec_RBF_Info {

    unsigned long int   *LookUpKey;  // keyword for hashing (an array of id numbers)
    int                 RadialBasisFunction;
    int                 DoPoly;      // Do a simultaneous linear polynomial fit.
    int                 n;
    double              Bx0;
    double              By0;
    double              Bz0;
    Lgm_Vector          *v;
    double              *cx;
    double              *cy;
    double              *cz;

    double              *cx_new;
    double              *cy_new;
    double              *cz_new;

    double              *gx_new;
    double              *gy_new;
    double              *gz_new;

    double              *eps_x;
    double              *eps_y;
    double              *eps_z;

    double              size;       // Size of memory consumed in MB

    UT_hash_handle      hh;         // Make structure hashable via uthash

} Lgm_Vec_RBF_Info;


Lgm_DFI_RBF_Info *Lgm_DFI_RBF_Init( unsigned long int *I, Lgm_Vector *v, Lgm_Vector *B, int n, double eps, int RadialBasisFunction );
void    Lgm_DFI_RBF_Free( Lgm_DFI_RBF_Info *rbf );
void    Lgm_DFI_RBF_Phi( Lgm_Vector *v, Lgm_Vector *v0, double Phi[3][3], Lgm_DFI_RBF_Info *rbf );
void    Lgm_DFI_RBF_Eval( Lgm_Vector *v, Lgm_Vector *B, Lgm_DFI_RBF_Info *rbf );
void    Lgm_DFI_RBF_dPhi_dx( Lgm_Vector *v, Lgm_Vector *v0, double dPdx[3][3], Lgm_DFI_RBF_Info *rbf  );
void    Lgm_DFI_RBF_dPhi_dy( Lgm_Vector *v, Lgm_Vector *v0, double dPdy[3][3], Lgm_DFI_RBF_Info *rbf  );
void    Lgm_DFI_RBF_dPhi_dz( Lgm_Vector *v, Lgm_Vector *v0, double dPdz[3][3], Lgm_DFI_RBF_Info *rbf  );
void    Lgm_DFI_RBF_Derivs_Eval( Lgm_Vector *v, Lgm_Vector *dBdx, Lgm_Vector *dBdy, Lgm_Vector *dBdz, Lgm_DFI_RBF_Info *rbf );


Lgm_Vec_RBF_Info *Lgm_Vec_RBF_Init( unsigned long int *I_data, Lgm_Vector *v, Lgm_Vector *B, double *eps_x, double *eps_y, double *eps_z, int n, int DoPoly, int RadialBasisFunction );
void    Lgm_Vec_RBF_Free( Lgm_Vec_RBF_Info *rbf );
void    Lgm_Vec_RBF_Psi( Lgm_Vector *v, Lgm_Vector *v0, double eps, double *Psi, int RbfType  );
void    Lgm_Vec_RBF_Psi2( Lgm_Vector *v, Lgm_Vector *v0, double eps_x, double eps_y, double eps_z, double *Psi, int RbfType  );
void    Lgm_Vec_RBF_Derivs( Lgm_Vector *v, Lgm_Vector *v0, double eps, double *dPdx, double *dPdy, double *dPdz, int RbfType  );
void    Lgm_Vec_RBF_Derivs2( Lgm_Vector *v, Lgm_Vector *v0, double eps_x, double eps_y, double eps_z, double *dPdx, double *dPdy, double *dPdz, int RbfType  );
void    Lgm_Vec_RBF_Eval( Lgm_Vector *v, Lgm_Vector *B, Lgm_Vec_RBF_Info *rbf );
void    Lgm_Vec_RBF_Derivs_Eval( Lgm_Vector *v, Lgm_Vector *dBdx, Lgm_Vector *dBdy, Lgm_Vector *dBdz, Lgm_Vec_RBF_Info *rbf );





#endif
