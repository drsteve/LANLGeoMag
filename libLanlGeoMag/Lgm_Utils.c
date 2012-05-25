#include <math.h>
#include "Lgm/Lgm_Utils.h"


/*
 * Fill an array with log spaced numberd form start to stop
 */
void Lgm_LogSpace( double start, double stop, int num, double *array ){
  double delta, logmin, accDelta = 0;
  int i;
  logmin = log10(start);
  delta = (log10(stop) - logmin)/(double)num;
  for (i=0; i<num; i++) {
    array[i] = pow(10, logmin + accDelta);
    accDelta += delta;
  }
}
