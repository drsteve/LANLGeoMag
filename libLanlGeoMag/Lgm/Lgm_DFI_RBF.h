#ifndef LGM_DFI_RBF
#define LGM_DFI_RBF

#include "Lgm_Vec.h"
#include "Lgm/Lgm_DFI_RBF.h"
#include "Lgm/Lgm_DynamicMemory.h"
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

#define LGM_CHOLESKY_DECOMP 0                                                                                                                                                                                    
#define LGM_PLU_DECOMP      1    

#define LGM_RBF_GAUSSIAN        0
#define LGM_RBF_MULTIQUADRIC    1
#define LGM_RBF_WENDLAND        2

typedef struct _Lgm_DFI_RBF_Info {

    int         RadialBasisFunction;
    int         n;
    int         n3;
    double      eps;
    Lgm_Vector  *v;
    Lgm_Vector  *c;

} Lgm_DFI_RBF_Info;


Lgm_DFI_RBF_Info *Lgm_DFI_RBF_Init( Lgm_Vector *v, Lgm_Vector *B, int n, double eps, int RadialBasisFunction );
void              Lgm_DFI_RBF_Phi( Lgm_Vector *v, Lgm_Vector *v0, double Phi[3][3], Lgm_DFI_RBF_Info *rbf );
void              Lgm_DFI_RBF_Free( Lgm_DFI_RBF_Info *rbf );
void              Lgm_DFI_RBF_Eval( Lgm_Vector *v, Lgm_Vector *B, Lgm_DFI_RBF_Info *rbf );


#endif
