#include <math.h>
#include "Lgm/Lgm_Utils.h"


/*
 * Fill an array with log spaced numbers from start to stop
 */
void Lgm_LogSpace( double start, double stop, long num, double *array ){
  double delta, logmin, accDelta = 0;
  long i;
  logmin = log10(start);
  delta = (log10(stop) - logmin)/(double)num;
  for (i=0; i<num; i++) {
    array[i] = pow(10, logmin + accDelta);
    accDelta += delta;
  }
}

/*
 * Fill an array with linear spaced numbers from start to stop, includes the end point
 */
void Lgm_LinSpace(double start, double stop, long num, double *array) {
  double delta, linmin, accDelta = 0;
  long i;
  linmin = start;
  delta = (stop-start)/((double)(num-1));
  for (i=0; i<num; i++) {
    array[i] = linmin + accDelta;
    accDelta += delta;
  }
}


/*
 * given a sorted array, data, find the index where value should be inserted to maintain order
 */
long Lgm_Bisect(double *data, double value, long len) {
  long mid, hi=len, lo = 0;
  while (lo <hi) {
    mid = (lo+hi) / 2;  // integer division
    if (value < data[mid])
      hi = mid;
    else
      lo = mid+1;
  }
  return (lo);
}


