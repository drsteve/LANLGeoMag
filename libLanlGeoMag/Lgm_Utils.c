#include "Lgm/Lgm_Utils.h"
#include <math.h>
#include <stdint.h>

/*
 * Fill an array with log spaced numbers from start to stop
 */
void Lgm_LogSpace(double start, double stop, long num, double* array) {
    /* logspace is equivalent to
     *   y=linspace(log10(start), log10(stop), num)
     *  then pow(10, y)
     */
    uint_fast64_t i;
    Lgm_LinSpace(log10(start), log10(stop), num, array);
    for (i = 0; i < num; i++)
        array[i] = pow(10., array[i]);
}

/*
 * Fill an array with linear spaced numbers from start to stop, includes the end
 * point
 */
void Lgm_LinSpace(double start, double stop, long num, double* array) {
    double delta, linmin, accDelta = 0;
    uint_fast64_t i;
    linmin = start;
    delta = (stop - start) / ((double)(num - 1));
    for (i = 0; i < num; i++) {
        array[i] = linmin + accDelta;
        accDelta += delta;
    }
}

/*
 * given a sorted array, data, find the index where value should be inserted to
 * maintain order
 */
long Lgm_Bisect(double* data, double value, long len) {
    uint_fast64_t mid, hi = len, lo = 0;
    while (lo < hi) {
        mid = (lo + hi) / 2;  // integer division
        if (value < data[mid])
            hi = mid;
        else
            lo = mid + 1;
    }
    return (lo);
}

/*
 * given an array return the the minimum and maximum of the array
 */
void Lgm_MinMax(double* inval, long len, double* min, double* max) {
    uint_fast64_t i;
    *min = inval[0];
    *max = inval[0];
    for (i = 1; i < len; i++) {
        if (inval[i] > *max)
            *max = inval[i];
        if (inval[i] < *min)
            *min = inval[i];
    }
}
