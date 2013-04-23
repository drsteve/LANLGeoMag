

#include <stdio.h>     

#include <Lgm_DynamicMemory.h>

double sum_dyn(double *inval, long x, long y, long z) {
  /* given a 3d array input sum over the whole thing using Lgm_DynamicMemory */

  // make a pointer to a pointer to a pointer for each dim of the input
  double ***invaltmp;
  double ans=0.;
  long i, j, k;

  // this macro makes invaltmp indexable via [][][] as inval is not
  // this does allocate memory so be sure a free later
  LGM_ARRAY_FROM_DATA_3D( invaltmp, inval, x, y, z, double );

  for (i=0; i<x; i++) {
    for (j=0; j<x; j++) {
      for (k=0; k<x; k++) {
	ans += invaltmp[i][j][k];
      }
    }
  }

  // free so that we don't leak memory
  LGM_ARRAY_FROM_DATA_3D_FREE(invaltmp);

  return ans;
}

double sum(double *inval, long x, long y, long z) {
  /* given a 3d array input sum over the whole thing using pointer math */


  /* For an array of [ A ][ B ][ C ], we can index it thus: */
  /* (a * B * C) + (b * C) + (c * 1) */
  /* Notice the pattern: for each dimension, its index is multiplied by the product of the cardinally of all dimensions minor to it. */
  /* To help digest that, here's another example. Given array[ A ][ B ][ C ][ D ]: */
  /* (a * B * C * D) + (b * C * D) + (c * D) + (d * 1) */

  double ans=0.;
  long i, j, k;

  for (i=0; i<x; i++) {
    for (j=0; j<x; j++) {
      for (k=0; k<x; k++) {
	ans += inval[i*y*z + j*z + k];
      }
    }
  }
  return ans;
}

int main(void) {
 
  double val[10][5][7];
  long i, j, k;
  
  for (i=0; i<10; i++) {
    for (j=0; j<5; j++) {
      for (k=0; k<7; k++) {
	val[i][j][k] = ((double)(i+j+k))/2.;
      }
    }
  }
  
  printf("Summing via Lgm_DynamicMemory.h gives: %lf\n", sum_dyn( (double*)val, 10, 5, 7));
  printf("Summing via pointer math gives:        %lf\n", sum( (double*)val, 10, 5, 7));



  return(0);
}







