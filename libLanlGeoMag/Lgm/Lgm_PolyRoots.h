#ifndef LGM_POLY_ROOTS_H
#define LGM_POLY_ROOTS_H
#include <math.h>
#include <complex.h>

int     Lgm_QuadraticRoots( double a, double b, double c, complex double *z1, complex double *z2 );
int     Lgm_CubicRoots( double b, double c, double d, double *z1, complex double *z2, complex double *z3 );
double  Lgm_CubicRealRoot( double b, double c, double d );
int     Lgm_QuarticRoots( double b, double c, double d, double e, double complex *z1, double complex *z2, double complex *z3, double complex *z4 );

#endif
