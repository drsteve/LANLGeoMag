#ifndef LGM_DFI_RBF
#define LGM_DFI_RBF

#include "Lgm_Vec.h"
#include "Lgm/Lgm_DFI_RBF.h"
#include "Lgm/Lgm_DynamicMemory.h"
#include <stdlib.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>


typedef struct _Lgm_DFI_RBF_Info {

    int                  n;
    int                  n3;
    double              eps;
    Lgm_Vector          *v;
    gsl_vector          *c;

} Lgm_DFI_RBF_Info;


void              Lgm_DFI_RBF_Phi( Lgm_Vector *v, Lgm_Vector *v0, double eps, double **Phi );
Lgm_DFI_RBF_Info *Lgm_DFI_RBF_Init( Lgm_Vector *v, Lgm_Vector *B, int n, double eps );
void              Lgm_DFI_RBF_Free( Lgm_DFI_RBF_Info *rbf );
void              Lgm_DFI_RBF_Eval( Lgm_Vector *v, Lgm_Vector *B, Lgm_DFI_RBF_Info *rbf );


#endif
