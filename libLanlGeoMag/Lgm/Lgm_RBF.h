#ifndef LGM_RBF
#define LGM_RBF

#include "Lgm/Lgm_Vec.h"
#include "Lgm/Lgm_RBF.h"
#include "Lgm/Lgm_DynamicMemory.h"
#include "Lgm/uthash.h"
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
//#define LGM_RBF_WENDLAND         3

/*
 * For  doing div-free interp.
 */
typedef struct _Lgm_DFI_RBF_Info {

    unsigned long int   *LookUpKey;  // keyword for hashing (an array of id numbers)
    int                 RadialBasisFunction;
    int                 n;
    int                 n3;
    double              eps;
    double              Bx0;
    double              By0;
    double              Bz0;
    Lgm_Vector          *v;
    Lgm_Vector          *c;
    UT_hash_handle      hh; // Make structure hashable via uthash

} Lgm_DFI_RBF_Info;




/*
 * For doing vec interp with each comp treated independantly.
 */
typedef struct _Lgm_Vec_RBF_Info {

    unsigned long int   *LookUpKey;  // keyword for hashing (an array of id numbers)
    int                 RadialBasisFunction;
    int                 n;
    double              eps;
    double              Bx0;
    double              By0;
    double              Bz0;
    Lgm_Vector          *v;
    double              *cx;
    double              *cy;
    double              *cz;
    UT_hash_handle      hh; // Make structure hashable via uthash

} Lgm_Vec_RBF_Info;


Lgm_DFI_RBF_Info *Lgm_DFI_RBF_Init( unsigned long int *I, Lgm_Vector *v, Lgm_Vector *B, int n, double eps, int RadialBasisFunction );
void    Lgm_DFI_RBF_Free( Lgm_DFI_RBF_Info *rbf );
void    Lgm_DFI_RBF_Phi( Lgm_Vector *v, Lgm_Vector *v0, double Phi[3][3], Lgm_DFI_RBF_Info *rbf );
void    Lgm_DFI_RBF_Eval( Lgm_Vector *v, Lgm_Vector *B, Lgm_DFI_RBF_Info *rbf );
void    Lgm_DFI_RBF_dPhi_dx( Lgm_Vector *v, Lgm_Vector *v0, double dPdx[3][3], Lgm_DFI_RBF_Info *rbf  );
void    Lgm_DFI_RBF_dPhi_dy( Lgm_Vector *v, Lgm_Vector *v0, double dPdy[3][3], Lgm_DFI_RBF_Info *rbf  );
void    Lgm_DFI_RBF_dPhi_dz( Lgm_Vector *v, Lgm_Vector *v0, double dPdz[3][3], Lgm_DFI_RBF_Info *rbf  );
void    Lgm_DFI_RBF_Derivs_Eval( Lgm_Vector *v, Lgm_Vector *dBdx, Lgm_Vector *dBdy, Lgm_Vector *dBdz, Lgm_DFI_RBF_Info *rbf );


Lgm_Vec_RBF_Info *Lgm_Vec_RBF_Init( unsigned long int *I_data, Lgm_Vector *v, Lgm_Vector *B, int n, double eps, int RadialBasisFunction );
void    Lgm_Vec_RBF_Free( Lgm_Vec_RBF_Info *rbf );
void    Lgm_Vec_RBF_Psi( Lgm_Vector *v, Lgm_Vector *v0, double *Psi, Lgm_Vec_RBF_Info *rbf  );
void    Lgm_Vec_RBF_Eval( Lgm_Vector *v, Lgm_Vector *B, Lgm_Vec_RBF_Info *rbf );
void    Lgm_Vec_RBF_Derivs( Lgm_Vector *v, Lgm_Vector *v0, double *dPdx, double *dPdy, double *dPdz, Lgm_Vec_RBF_Info *rbf  );
void    Lgm_Vec_RBF_Derivs_Eval( Lgm_Vector *v, Lgm_Vector *dBdx, Lgm_Vector *dBdy, Lgm_Vector *dBdz, Lgm_Vec_RBF_Info *rbf );





#endif
