#ifndef LGM_QUADPACK_H
#define LGM_QUADPACK_H 1

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/**
 * @file QuadPack.h
 * @brief Cified FORTRAN QUADPACK routines
 */

#define         TRUE    1
#define         FALSE   0
#define         dmax1(a, b)     ( ((a) > (b)) ? (a) : (b) )
#define         dmin1(a, b)     ( ((a) < (b)) ? (a) : (b) )


#ifndef _qpInfo
typedef int _qpInfo;
#endif

double d1mach( int i );

/*! \brief FORTRAN QUADPACK routine DQAGS
 *         Definite integral, double-precision, general purpose
 */
int dqags(double (*f)( double, _qpInfo *), //!< integrand function f(x, param)
    _qpInfo *qpInfo, //!< auxiliary information to pass to function
    double a, //!< integral lower limit
    double b, //!< integral upper limit
    double epsabs, //!< absolute accuracy requested
    double epsrel, //!< relative accuracy requested
    double *result, //!< calculated integral of f from a to b
    double *abserr, //!< estimated modulus of absolute error in result
    int *neval, //!< number of function evaluations
    int *ier, //!< flag indicating error if positive
    int limit, //!< positive number determining maximum #subintervals used
    int lenw, //!< must be at least limit*4
    int *last, //!< on return indicates #subintervals produced
    int *iwork, //!< must be allocated with size >= limit
    double *work, //!< must be allocated with size >= lenw
    int verbosity //!< detail level of messaging to stdout
    );

/*! \brief FORTRAN QUADPACK routine DQAGSE
 *         Definite integral, double-precision, general purpose, finer control
 *         and more information returned than DQAGS
 */
int dqagse(double (*f)( double, _qpInfo *), //!< integrand function f(x, param)
  _qpInfo *qpInfo, //!< auxiliary information to pass to function
  double a, //!< integral lower limit
  double b, //!< integral upper limit
  double epsabs, //!< absolute accuracy requested
  double epsrel, //!< relative accuracy requested
  int limit, //!<  positive number determining maximum #subintervals used
  double *result, //!< calculated integral of f from a to b
  double *abserr, //!< estimated modulus of absolute error in result
  int *neval, //!< number of function evaluations
  int *ier, //!< flag indicating error if positive
  double *alist, //!< must be allocated with size >= limit, left endpoints of
                 //!< subintervals in the partition of (a, b)
  double *blist, //!< must be allocated with size >= limit, right endpoints of
                 //!< subintervals in the partition of (a, b)
  double *rlist, //!< must be allocated with size >= limit, approximations of
                 //!< integral on the subintervals
  double *elist, //!< must be allocated with size >= limit, moduli of absolute 
                 //!< error estimates on the subintervals
  int *iord, //!< must be allocated with size>= limit, pointers to largest
             //!< error estimates
  int *last  //!< number of subintervals actually produced in the subdivision
  );

/*! \brief FORTRAN QUADPACK routine DQAGP
 *         Definite integral, double-precision, singularities or discontinuities
 */
int dqagp(double (*f)( double, _qpInfo *), //!< integrand function f(x, param)
  _qpInfo *qpInfo, //!< auxuliary information to pass to function
  double a, //!< integral lower limit
  double b, //!< integral upper limit
  int npts2, //!< number of user-supplied break points + 2
  double *points, //!< user-supplied breakpoints
  double epsabs, //!< absolute accuracy requested
  double epsrel, //!< relative accuracy requested
  double *result, //!< calculated integral of f from a to b
  double *abserr, //!< estimated modulus of absolute error
  int *neval, //!< number of integrand evaluations
  int *ier, //!< flag indicating error if positive
  int leniw, //!< dimensioning parameter for iwork
  int lenw, //!< dimensioning parameter for work
  int *last, //!< on return, #subintervals produced in the subdivision
  int *iwork, //!< must be allocated with size >= leniw, pointers to largest
              //!< error estimates
  double *work, //!< must be allocated with size at least lenw, pointers to
                //!< interval endpoints, integral approximations, corresponding
                //!< error estimates, and integration limits/break points sorted
  int verbosity //!< detail level of messaging to stdout
);


/*! \brief FORTRAN QUADPACK routine DQAGPE
 *         Definite integral, double-precision, singularities or discontinuities
 *         More control and information returned than DQAGP
 */
int dqagpe(double (*f)( double, _qpInfo *), //!< integrand function f(x, param)
  _qpInfo *qpInfo, //!< auxuliary information to pass to function
  double a, //!< integral lower limit
  double b, //!< integral upper limit
  int npts2, //!< number of user-supplied break points + 2
  double *points, //!< user-supplied breakpoints
  double epsabs, //!< absolute accuracy requested
  double epsrel, //!< relative accuracy requested
  int limit, //!< upper bound on number of subintervals of (a, b) to use
  double *result, //!< calculated integral of f from a to b
  double *abserr, //!< estimated modulus of absolute error
  int *neval, //!< number of integrand evaluations
  int *ier, //!< flag indicating error if positive
  double *alist, //!< must be allocated with size at least limit,
                 //!< left endpoints of subintervals
  double *blist, //!< must be allocated with size at least limit,
                 //!< right endpoints of subintervals              
  double *rlist, //!< must be allocated with size at least limit,
                 //!< integral approximations on the subintervals
  double *elist, //!< must be allocated with size at least limit,
                 //!< moduli of absolute error estimates on subintervals
  double *pts, //!< must be allocated with size at least npts2,
               //!< integration limits and breakpoints in ascending sequence
  int *iord, //!< must be allocated with size at least limit,
             //!< pointers to error estimates in decreasing order
  int *level, //!< must be allocated with size at least limit,
              //!< subdivision levels of subintervals
  int *ndin, //!< must be allocated with size at least limit, 
             //!< indicates whether an interval's error estimate was increased
             //!< artificially after first integral step in order to push
             //!< a subdivision
  int *last //!< number of subintervals actually produced in subdivision
  );

/*! \brief FORTRAN QUADPACK routine DQK21
 *         Definite integral using 21-point Gauss-Kronrod rule, with error 
 *         estimate j = integral of |f| over (a, b)
 */
int dqk21(double (*f)( double, _qpInfo *), //!< integrand function f(x, param)
  _qpInfo *qpInfo, //!< auxiliary information to pass to function
  double a, //!< integral lower limit
  double b, //!< integral upper limit
  double *result, //!< calculated integral of f from a to b
  double *abserr, //!< estimate of modulus of absolute error
  double *resabs, //!< approximation to error estimate j
  double *resasc //!< approximation to the integral of |f-i/(b-a)| over (a, b)
  );

/*! \brief FORTRAN QUADPACK routine DQELG
 *         estimates limit of a sequence of approximations and error
 *         Uses P. Wynn's epsilon algorithm
 */
int dqelg(int n, //!< new element in first column of epsilon table
  double epstab[], //!< vector of size 52 containing elements of the two lower
                   //!< diagonals of the triangular epsilon table, numbered
                   //!< starting at the right-hand corner of the triangle
  double *result, //!< resulting approximation
  double *abserr, //!< estimate of absolute error from RESULT and the 3 previous
  double res3la[], //!< vector of dimension 3 containing the last 3 results
  int *nres //!< number of calls to the routine
  );

/**
 * \brief FORTRAN QUADPACK routine DQPSRT
 *        Maintains the descending ordering in the list of the local error
 *        estimated resulting from the interval subdivision process. At each
 *        call two error estimates are inserted using the sequential search
 *        method, top-down for the largest error estimate and bottom-up for the
 *        smallest error estimate.
 */
int dqpsrt(int limit, //!< maximum number of error estimates the list can have
  int last, //!< number of error estimates currently in the list
  int *maxerr, //!< points to NRMAX-th largest error estimate currently in list
  double *ermax, //!< NRMAX-th largest error estimate
  double elist[], //!< dimension LAST vector containing error estimates
  int iord[], //!< dimension LAST vector whose first K elements contain 
              //!< pointers to error estimates in decreasing sequence
  int *nrmax //!< MAXERR = IORD(NRMAX), in accordance with the prophecy
  );

/**
 * \brief PrintQuadpackError Prints plain-English description of ier flag output
 *        Will print "Unknown error" if given something that isn't a known
 *        error, even if that something is "no error happened"
 */
void PrintQuadpackError( int ier //!< error flag given by QUADPACK routine
    );

#endif 


