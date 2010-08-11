#ifndef LGM_QUADPACK_H
#define LGM_QUADPACK_H 1

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/*
 *    QuadPack.h
 */


#define         TRUE    1
#define         FALSE   0
#define         dmax1(a, b)     ( ((a) > (b)) ? (a) : (b) )
#define         dmin1(a, b)     ( ((a) < (b)) ? (a) : (b) )


#ifndef _qpInfo
typedef int _qpInfo;
#endif

double d1mach( int i );

int dqags(double (*f)( double, _qpInfo *), _qpInfo *qpInfo, double a, double b, 
	double epsabs, double epsrel, double *result, double *abserr, int *neval, 
	int *ier, int limit, int lenw, int *last, int *iwork, double *work);

int dqagse(double (*f)( double, _qpInfo *), _qpInfo *qpInfo, double a, double b, 
	double epsabs, double epsrel, int limit, double *result, double *abserr, 
	int *neval, int *ier, double *alist, double *blist, double *rlist, 
	double *elist, int *iord, int *last);


int dqagp(double (*f)( double, _qpInfo *), _qpInfo *qpInfo, double a, double b, 
    int npts2, double   *points, double epsabs, double epsrel, double *result, 
    double *abserr, int *neval, int *ier, int leniw, int lenw, int *last, 
    int *iwork, double *work);

int dqagpe(double (*f)( double, _qpInfo *), _qpInfo *qpInfo, double a, double b, 
    int npts2, double  *points, double epsabs, double epsrel, int limit, 
    double *result, double *abserr, int *neval, int *ier, double *alist, 
    double *blist, double *rlist, double *elist, double *pts, int *iord, 
    int  *level, int  *ndin, int *last);


int dqk21(double (*f)( double, _qpInfo *), _qpInfo *qpInfo, double a, double b, 
	double *result, double *abserr, double *resabs, double *resasc);

int dqelg(int n, double epstab[], double *result, double *abserr, double res3la[], 
	int *nres);

int dqpsrt(int limit, int last, int *maxerr, double *ermax, double elist[], 
	int iord[], int *nrmax);

#endif 
