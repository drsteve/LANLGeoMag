#ifndef LGM_POLY_ROOTS_H
#define LGM_POLY_ROOTS_H
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>

#define LGM_POLYROOTS_MAXITER 500

int     Lgm_QuadraticRoots( double a, double b, double c, complex double *z1, complex double *z2 );
int     Lgm_CubicRoots( double b, double c, double d, double *z1, complex double *z2, complex double *z3 );
double  Lgm_CubicRealRoot( double b, double c, double d );
int     Lgm_QuarticRoots( double b, double c, double d, double e, double complex *z1, double complex *z2, double complex *z3, double complex *z4 );
int     Lgm_QuarticRootsSorted( double b, double c, double d, double e, int *nReal, double *RealRoots, int *nComplex, double complex *ComplexRoots );

int     Lgm_PolyRoots( double *a, int n, double complex *z );
int     Lgm_PolyRoots_Roots( double *a, int n, double *wr, double *wi );
void    Lgm_PolyRootsDeflate( double *a, int n, double *b, double *quad, double *err );
void    Lgm_PolyRoots_FindQuad( double *a, int n, double *b, double *quad, double *err, int *iter );
void    Lgm_PolyRoots_DiffPoly( double *a, int n, double *b );
void    Lgm_PolyRoots_Recurse( double *a, int n, double *b, int m, double *quad, double *err, int *iter );
void    Lgm_PolyRoots_GetQuad( double *a, int n, double *quad, double *x );


#endif
