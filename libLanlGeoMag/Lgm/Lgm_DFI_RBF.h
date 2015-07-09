#ifndef LGM_DFI_RBF
#define LGM_DFI_RBF

#include "Lgm/Lgm_Vec.h"
#include "Lgm/Lgm_DFI_RBF.h"
#include "Lgm/Lgm_DynamicMemory.h"
#include "Lgm/uthash.h"
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

#define LGM_CHOLESKY_DECOMP 0                                                                                                                                                                                    
#define LGM_PLU_DECOMP      1    

#define LGM_RBF_GAUSSIAN         0
#define LGM_RBF_MULTIQUADRIC     1
#define LGM_RBF_INV_MULTIQUADRIC 2
#define LGM_RBF_WENDLAND         3

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


Lgm_DFI_RBF_Info *Lgm_DFI_RBF_Init( unsigned long int *I, Lgm_Vector *v, Lgm_Vector *B, int n, double eps, int RadialBasisFunction );
void              Lgm_DFI_RBF_Phi( Lgm_Vector *v, Lgm_Vector *v0, double Phi[3][3], Lgm_DFI_RBF_Info *rbf );
void              Lgm_DFI_RBF_Free( Lgm_DFI_RBF_Info *rbf );
void              Lgm_DFI_RBF_Eval( Lgm_Vector *v, Lgm_Vector *B, Lgm_DFI_RBF_Info *rbf );

void    Lgm_DFI_RBF_dPhi_dx( Lgm_Vector *v, Lgm_Vector *v0, double dPdx[3][3], Lgm_DFI_RBF_Info *rbf  );
void    Lgm_DFI_RBF_dPhi_dy( Lgm_Vector *v, Lgm_Vector *v0, double dPdy[3][3], Lgm_DFI_RBF_Info *rbf  );
void    Lgm_DFI_RBF_dPhi_dz( Lgm_Vector *v, Lgm_Vector *v0, double dPdz[3][3], Lgm_DFI_RBF_Info *rbf  );


#endif
