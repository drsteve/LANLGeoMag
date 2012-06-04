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


/*
 * given a sorted array, data, find the index where value should be inserted to maintain order
 */
int Lgm_Bisect(double *data, double value, int len) {
  int mid, hi=len, lo = 0;
  while (lo <hi) {
    mid = (lo+hi) / 2;  // integer division
    if (value < data[mid])
      hi = mid;
    else
      lo = mid+1;
  }
  return (lo);
}


