#include "Lgm/Lgm_QuadPack.h"

/*
 *      QUADPACK DQAGP Routine converted to C 
 */
int dqagp(f, qpInfo, a, b, npts2, points, epsabs, epsrel, result, abserr, neval, ier, leniw, lenw, last, iwork, work)
double  (*f)( double, _qpInfo *); /*  The integrand function -- I.e. the function to integrate    	 */
_qpInfo *qpInfo;        	  /*  Auxilliary information to pass to function (to avoid making globals) */
double   a;             	  /*  Lower Limit of integration.                                 	 */
double   b;             	  /*  Upper limit of integration.                                 	 */
int      npts2;
double   *points;
double   epsabs;        	  /*  Absolute accuracy requested.                                	 */
double   epsrel;        	  /*  Relative accuracy requested.                                	 */
double  *result;        	  /*  The desired result. I.e. integral of f() from a to b        	 */
double  *abserr;        	  /*  Estimate of the modulus of the absolute error in the result 	 */
int     *neval;         	  /*  The number of integrand evaluations performed.              	 */
int     *ier;           	  /*  Error flag. An error occurred if ier > 0. See below.        	 */
int	     leniw;
int	     lenw;
int	    *last;
int	    *iwork;
double	*work;
{


    /*    EGIN PROLOGUE  DQAGP
	 *    URPOSE  The routine calculates an approximation result to a given
	 *            definite integral I = Integral of F over (A,B),
	 *            hopefully satisfying following claim for accuracy
	 *            break points of the integration interval, where local
	 *            difficulties of the integrand may occur (e.g.
	 *            SINGULARITIES, DISCONTINUITIES), are provided by the user.
	 *    IBRARY   SLATEC (QUADPACK)
	 *    ATEGORY  H2A2A1
	 *    YPE      DOUBLE PRECISION (QAGP-S, DQAGP-D)
	 *    EYWORDS  AUTOMATIC INTEGRATOR, EXTRAPOLATION, GENERAL-PURPOSE,
	 *             GLOBALLY ADAPTIVE, QUADPACK, QUADRATURE,
	 *             SINGULARITIES AT USER SPECIFIED POINTS
	 *    UTHOR  Piessens, Robert
	 *             Applied Mathematics and Programming Division
	 *             K. U. Leuven
	 *           de Doncker, Elise
	 *             Applied Mathematics and Programming Division
	 *             K. U. Leuven
	 *    ESCRIPTION
	 *
	 *        Computation of a definite integral
	 *        Standard fortran subroutine
	 *        Double precision version
	 *
	 *        PARAMETERS
	 *         ON ENTRY
	 *            F      - Double precision
	 *                     Function subprogram defining the integrand
	 *                     Function F(X). The actual name for F needs to be
	 *                     declared E X T E R N A L in the driver program.
	 *
	 *            A      - Double precision
	 *                     Lower limit of integration
	 *
	 *            B      - Double precision
	 *                     Upper limit of integration
	 *
	 *            NPTS2  - Integer
	 *                     Number equal to two more than the number of
	 *                     user-supplied break points within the integration
	 *                     range, NPTS.GE.2.
	 *                     If NPTS2.LT.2, The routine will end with IER = 6.
	 *
	 *            POINTS - Double precision
	 *                     Vector of dimension NPTS2, the first (NPTS2-2)
	 *                     elements of which are the user provided break
	 *                     points. If these points do not constitute an
	 *                     ascending sequence there will be an automatic
	 *                     sorting.
	 *
	 *            EPSABS - Double precision
	 *                     Absolute accuracy requested
	 *            EPSREL - Double precision
	 *                     Relative accuracy requested
	 *                     If  EPSABS.LE.0
	 *                     And EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28),
	 *                     The routine will end with IER = 6.
	 *
	 *         ON RETURN
	 *            RESULT - Double precision
	 *                     Approximation to the integral
	 *
	 *            ABSERR - Double precision
	 *                     Estimate of the modulus of the absolute error,
	 *                     which should equal or exceed ABS(I-RESULT)
	 *
	 *            NEVAL  - Integer
	 *                     Number of integrand evaluations
	 *
	 *            IER    - Integer
	 *                     IER = 0 Normal and reliable termination of the
	 *                             routine. It is assumed that the requested
	 *                             accuracy has been achieved.
	 *                     IER.GT.0 Abnormal termination of the routine.
	 *                             The estimates for integral and error are
	 *                             less reliable. it is assumed that the
	 *                             requested accuracy has not been achieved.
	 *            ERROR MESSAGES
	 *                     IER = 1 Maximum number of subdivisions allowed
	 *                             has been achieved. one can allow more
	 *                             subdivisions by increasing the value of
	 *                             LIMIT (and taking the according dimension
	 *                             adjustments into account). However, if
	 *                             this yields no improvement it is advised
	 *                             to analyze the integrand in order to
	 *                             determine the integration difficulties. If
	 *                             the position of a local difficulty can be
	 *                             determined (i.e. SINGULARITY,
	 *                             DISCONTINUITY within the interval), it
	 *                             should be supplied to the routine as an
	 *                             element of the vector points. If necessary
	 *                             an appropriate special-purpose integrator
	 *                             must be used, which is designed for
	 *                             handling the type of difficulty involved.
	 *                         = 2 The occurrence of roundoff error is
	 *                             detected, which prevents the requested
	 *                             tolerance from being achieved.
	 *                             The error may be under-estimated.
	 *                         = 3 Extremely bad integrand behaviour occurs
	 *                             at some points of the integration
	 *                             interval.
	 *                         = 4 The algorithm does not converge.
	 *                             roundoff error is detected in the
	 *                             extrapolation table.
	 *                             It is presumed that the requested
	 *                             tolerance cannot be achieved, and that
	 *                             the returned RESULT is the best which
	 *                             can be obtained.
	 *                         = 5 The integral is probably divergent, or
	 *                             slowly convergent. it must be noted that
	 *                             divergence can occur with any other value
	 *                             of IER.GT.0.
	 *                         = 6 The input is invalid because
	 *                             NPTS2.LT.2 or
	 *                             break points are specified outside
	 *                             the integration range or
	 *                             (EPSABS.LE.0 and
	 *                              EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28))
	 *                             RESULT, ABSERR, NEVAL, LAST are set to
	 *                             zero.  Except when LENIW or LENW or NPTS2
	 *                             is invalid, IWORK(1), IWORK(LIMIT+1),
	 *                             WORK(LIMIT*2+1) and WORK(LIMIT*3+1)
	 *                             are set to zero.
	 *                             WORK(1) is set to A and WORK(LIMIT+1)
	 *                             to B (where LIMIT = (LENIW-NPTS2)/2).
	 *
	 *         DIMENSIONING PARAMETERS
	 *            LENIW - Integer
	 *                    Dimensioning parameter for IWORK
	 *                    LENIW determines LIMIT = (LENIW-NPTS2)/2,
	 *                    which is the maximum number of subintervals in the
	 *                    partition of the given integration interval (A,B),
	 *                    LENIW.GE.(3*NPTS2-2).
	 *                    If LENIW.LT.(3*NPTS2-2), the routine will end with
	 *                    IER = 6.
	 *
	 *            LENW  - Integer
	 *                    Dimensioning parameter for WORK
	 *                    LENW must be at least LENIW*2-NPTS2.
	 *                    If LENW.LT.LENIW*2-NPTS2, the routine will end
	 *                    with IER = 6.
	 *
	 *            LAST  - Integer
	 *                    On return, LAST equals the number of subintervals
	 *                    produced in the subdivision process, which
	 *                    determines the number of significant elements
	 *                    actually in the WORK ARRAYS.
	 *
	 *         WORK ARRAYS
	 *            IWORK - Integer
	 *                    Vector of dimension at least LENIW. on return,
	 *                    the first K elements of which contain
	 *                    pointers to the error estimates over the
	 *                    subintervals, such that WORK(LIMIT*3+IWORK(1)),...,
	 *                    WORK(LIMIT*3+IWORK(K)) form a decreasing
	 *                    sequence, with K = LAST if LAST.LE.(LIMIT/2+2), and
	 *                    K = LIMIT+1-LAST otherwise
	 *                    IWORK(LIMIT+1), ...,IWORK(LIMIT+LAST) Contain the
	 *                     subdivision levels of the subintervals, i.e.
	 *                     if (AA,BB) is a subinterval of (P1,P2)
	 *                     where P1 as well as P2 is a user-provided
	 *                     break point or integration LIMIT, then (AA,BB) has
	 *                     level L if ABS(BB-AA) = ABS(P2-P1)*2**(-L),
	 *                    IWORK(LIMIT*2+1), ..., IWORK(LIMIT*2+NPTS2) have
	 *                     no significance for the user,
	 *                    note that LIMIT = (LENIW-NPTS2)/2.
	 *
	 *            WORK  - Double precision
	 *                    Vector of dimension at least LENW
	 *                    on return
	 *                    WORK(1), ..., WORK(LAST) contain the left
	 *                     end points of the subintervals in the
	 *                     partition of (A,B),
	 *                    WORK(LIMIT+1), ..., WORK(LIMIT+LAST) contain
	 *                     the right end points,
	 *                    WORK(LIMIT*2+1), ..., WORK(LIMIT*2+LAST) contain
	 *                     the integral approximations over the subintervals,
	 *                    WORK(LIMIT*3+1), ..., WORK(LIMIT*3+LAST)
	 *                     contain the corresponding error estimates,
	 *                    WORK(LIMIT*4+1), ..., WORK(LIMIT*4+NPTS2)
	 *                     contain the integration limits and the
	 *                     break points sorted in an ascending sequence.
	 *                    note that LIMIT = (LENIW-NPTS2)/2.
	 *
	 *   REFERENCES  (NONE)
	 *   ROUTINES CALLED  DQAGPE, XERMSG
	 *   REVISION HISTORY  (YYMMDD)
	 *   800101  DATE WRITTEN
	 *   890831  Modified array declarations.  (WRB)
	 *   890831  REVISION DATE from Version 3.2
	 *   891214  Prologue converted to Version 4.0 format.  (BAB)
	 *   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
     *   20071127 Converted to C by M. Henderson
	 *   END PROLOGUE  DQAGP
     */





    int 	lvl=0, l1, l2, l3, l4, limit;


    /* First executable statement  dqagp */


    /*
     *    Check validity of limit and lenw.
     */
    *ier = 6;
    *neval = 0;
    *last = 0;
    *result = 0.0;
    *abserr = 0.0;

    if( (leniw < (3*npts2-2)) || (lenw < (leniw*2-npts2)) || (npts2 < 2) ) {

        if (*ier == 6) lvl = 1;
        if (*ier != 0) fprintf(stderr, "Abnormal return from dqagp, ier = %d,  lvl = %d\n", *ier, lvl);
printf("leniw = %d    lenw = %d    npts2 = %d\n", leniw, lenw, npts2);
        return(-1);

    } else {

        /*
         *    Prepare call for dqagse.
         */
        limit = (leniw-npts2)/2;
        l1    = limit+1;
        l2    = limit+l1;
        l3    = limit+l2;
        l4    = limit+l3;


        dqagpe(f, qpInfo, a, b, npts2, points, epsabs, epsrel, limit, result, abserr, neval, ier, 
				work, work+l1, work+l2, work+l3, work+l4, iwork, iwork+l1, iwork+l2, last);

    }



    /*
     *    Call error handler if necessary.
     */
    lvl = 0;


    if(*ier == 6) lvl = 1;
    if(*ier != 0){ 
	fprintf(stderr, "Abnormal return from dqagp, ier = %d,  lvl = %d\n", *ier, lvl);
        return(-1);
    } else {
        return(1);
    }


}




/*
 *      QUADPACK DQAGPE Routine converted to C
 */
int dqagpe(f, qpInfo, a, b, npts2, points, epsabs, epsrel, limit, result, abserr, neval, ier, alist, blist, rlist, elist, pts, iord, level, ndin, last)
double  (*f)( double, _qpInfo *); /*  The integrand function -- I.e. the function to integrate           */
_qpInfo *qpInfo;                  /*  Auxilliary information to pass to function (to avoid making globals) */
double   a;             	  /*  Lower Limit of integration.                                 */
double   b;             	  /*  Upper limit of integration.                                 */
int      npts2;
double  *points;
double   epsabs;        	  /*  Absolute accuracy requested.                                */
double   epsrel;        	  /*  Relative accuracy requested.                                */
int      limit;
double  *result;        	  /*  The desired result. I.e. integral of f() from a to b        */
double  *abserr;        	  /*  Estimate of the modulus of the absolute error in the result */
int     *neval;         	  /*  The number of integrand evaluations performed.              */
int     *ier;           	  /*  Error flag. An error occurred if ier > 0. See below.        */
double   *alist;
double   *blist;
double   *rlist;
double   *elist;
double   *pts;
int	 *iord;
int	 *level;
int	 *ndin;
int	*last;
{


	/*   BEGIN PROLOGUE  DQAGPE
	 *   PURPOSE  Approximate a given definite integral I = Integral of F
	 *            over (A,B), hopefully satisfying the accuracy claim:
	 *                 ABS(I-RESULT).LE.MAX(EPSABS,EPSREL*ABS(I)).
	 *            Break points of the integration interval, where local
	 *            difficulties of the integrand may occur (e.g. singularities
	 *            or discontinuities) are provided by the user.
	 *   LIBRARY   SLATEC (QUADPACK)
	 *   CATEGORY  H2A2A1
	 *   TYPE      DOUBLE PRECISION (QAGPE-S, DQAGPE-D)
	 *   KEYWORDS  AUTOMATIC INTEGRATOR, EXTRAPOLATION, GENERAL-PURPOSE,
	 *             GLOBALLY ADAPTIVE, QUADPACK, QUADRATURE,
	 *             SINGULARITIES AT USER SPECIFIED POINTS
	 *   AUTHOR  Piessens, Robert
	 *             Applied Mathematics and Programming Division
	 *             K. U. Leuven
	 *           de Doncker, Elise
	 *             Applied Mathematics and Programming Division
	 *             K. U. Leuven
	 *   DESCRIPTION
	 *
	 *        Computation of a definite integral
	 *        Standard fortran subroutine
	 *        Double precision version
	 *
	 *        PARAMETERS
	 *         ON ENTRY
	 *            F      - Double precision
	 *                     Function subprogram defining the integrand
	 *                     function F(X). The actual name for F needs to be
	 *                     declared E X T E R N A L in the driver program.
	 *
	 *            A      - Double precision
	 *                     Lower limit of integration
	 *
	 *            B      - Double precision
	 *                     Upper limit of integration
	 *
	 *            NPTS2  - Integer
	 *                     Number equal to two more than the number of
	 *                     user-supplied break points within the integration
	 *                     range, NPTS2.GE.2.
	 *                     If NPTS2.LT.2, the routine will end with IER = 6.
	 *
	 *            POINTS - Double precision
	 *                     Vector of dimension NPTS2, the first (NPTS2-2)
	 *                     elements of which are the user provided break
	 *                     POINTS. If these POINTS do not constitute an
	 *                     ascending sequence there will be an automatic
	 *                     sorting.
	 *
	 *            EPSABS - Double precision
	 *                     Absolute accuracy requested
	 *            EPSREL - Double precision
	 *                     Relative accuracy requested
	 *                     If  EPSABS.LE.0
	 *                     and EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28),
	 *                     the routine will end with IER = 6.
	 *
	 *            LIMIT  - Integer
	 *                     Gives an upper bound on the number of subintervals
	 *                     in the partition of (A,B), LIMIT.GE.NPTS2
	 *                     If LIMIT.LT.NPTS2, the routine will end with
	 *                     IER = 6.
	 *
	 *         ON RETURN
	 *            RESULT - Double precision
	 *                     Approximation to the integral
	 *
	 *            ABSERR - Double precision
	 *                     Estimate of the modulus of the absolute error,
	 *                     which should equal or exceed ABS(I-RESULT)
	 *
	 *            NEVAL  - Integer
	 *                     Number of integrand evaluations
	 *
	 *            IER    - Integer
	 *                     IER = 0 Normal and reliable termination of the
	 *                             routine. It is assumed that the requested
	 *                             accuracy has been achieved.
	 *                     IER.GT.0 Abnormal termination of the routine.
	 *                             The estimates for integral and error are
	 *                             less reliable. It is assumed that the
	 *                             requested accuracy has not been achieved.
	 *            ERROR MESSAGES
	 *                     IER = 1 Maximum number of subdivisions allowed
	 *                             has been achieved. One can allow more
	 *                             subdivisions by increasing the value of
	 *                             LIMIT (and taking the according dimension
	 *                             adjustments into account). However, if
	 *                             this yields no improvement it is advised
	 *                             to analyze the integrand in order to
	 *                             determine the integration difficulties. If
	 *                             the position of a local difficulty can be
	 *                             determined (i.e. SINGULARITY,
	 *                             DISCONTINUITY within the interval), it
	 *                             should be supplied to the routine as an
	 *                             element of the vector points. If necessary
	 *                             an appropriate special-purpose integrator
	 *                             must be used, which is designed for
	 *                             handling the type of difficulty involved.
	 *                         = 2 The occurrence of roundoff error is
	 *                             detected, which prevents the requested
	 *                             tolerance from being achieved.
	 *                             The error may be under-estimated.
	 *                         = 3 Extremely bad integrand behaviour occurs
	 *                             At some points of the integration
	 *                             interval.
	 *                         = 4 The algorithm does not converge.
	 *                             Roundoff error is detected in the
	 *                             extrapolation table. It is presumed that
	 *                             the requested tolerance cannot be
	 *                             achieved, and that the returned result is
	 *                             the best which can be obtained.
	 *                         = 5 The integral is probably divergent, or
	 *                             slowly convergent. It must be noted that
	 *                             divergence can occur with any other value
	 *                             of IER.GT.0.
	 *                         = 6 The input is invalid because
	 *                             NPTS2.LT.2 or
	 *                             Break points are specified outside
	 *                             the integration range or
	 *                             (EPSABS.LE.0 and
	 *                              EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28))
	 *                             or LIMIT.LT.NPTS2.
	 *                             RESULT, ABSERR, NEVAL, LAST, RLIST(1),
	 *                             and ELIST(1) are set to zero. ALIST(1) and
	 *                             BLIST(1) are set to A and B respectively.
	 *
	 *            ALIST  - Double precision
	 *                     Vector of dimension at least LIMIT, the first
	 *                      LAST  elements of which are the left end points
	 *                     of the subintervals in the partition of the given
	 *                     integration range (A,B)
	 *
	 *            BLIST  - Double precision
	 *                     Vector of dimension at least LIMIT, the first
	 *                      LAST  elements of which are the right end points
	 *                     of the subintervals in the partition of the given
	 *                     integration range (A,B)
	 *
	 *            RLIST  - Double precision
	 *                     Vector of dimension at least LIMIT, the first
	 *                      LAST  elements of which are the integral
	 *                     approximations on the subintervals
	 *
	 *            ELIST  - Double precision
	 *                     Vector of dimension at least LIMIT, the first
	 *                      LAST  elements of which are the moduli of the
	 *                     absolute error estimates on the subintervals
	 *
	 *            PTS    - Double precision
	 *                     Vector of dimension at least NPTS2, containing the
	 *                     integration limits and the break points of the
	 *                     interval in ascending sequence.
	 *
	 *            LEVEL  - Integer
	 *                     Vector of dimension at least LIMIT, containing the
	 *                     subdivision levels of the subinterval, i.e. if
	 *                     (AA,BB) is a subinterval of (P1,P2) where P1 as
	 *                     well as P2 is a user-provided break point or
	 *                     integration limit, then (AA,BB) has level L if
	 *                     ABS(BB-AA) = ABS(P2-P1)*2**(-L).
	 *
	 *            NDIN   - Integer
	 *                     Vector of dimension at least NPTS2, after first
	 *                     integration over the intervals (PTS(I)),PTS(I+1),
	 *                     I = 0,1, ..., NPTS2-2, the error estimates over
	 *                     some of the intervals may have been increased
	 *                     artificially, in order to put their subdivision
	 *                     forward. If this happens for the subinterval
	 *                     numbered K, NDIN(K) is put to 1, otherwise
	 *                     NDIN(K) = 0.
	 *
	 *            IORD   - Integer
	 *                     Vector of dimension at least LIMIT, the first K
	 *                     elements of which are pointers to the
	 *                     error estimates over the subintervals,
	 *                     such that ELIST(IORD(1)), ..., ELIST(IORD(K))
	 *                     form a decreasing sequence, with K = LAST
	 *                     If LAST.LE.(LIMIT/2+2), and K = LIMIT+1-LAST
	 *                     otherwise
	 *
	 *            LAST   - Integer
	 *                     Number of subintervals actually produced in the
	 *                     subdivisions process
	 *
	 *   REFERENCES  (NONE)
	 *   ROUTINES CALLED  D1MACH, DQELG, DQK21, DQPSRT
	 *   REVISION HISTORY  (YYMMDD)
	 *   800101  DATE WRITTEN
	 *   890531  Changed all specific intrinsics to generic.  (WRB)
	 *   890831  Modified array declarations.  (WRB)
	 *   890831  REVISION DATE from Version 3.2
	 *   891214  Prologue converted to Version 4.0 format.  (BAB)
     *   20071127 Converted to C by M. Henderson
	 *   END PROLOGUE  DQAGPE
     */





    double 	area, abseps, area1, area12, area2, a1;
    double  a2, b1, b2, correc=0.0, defabs, defab1, defab2, d1mach();
    double  dres, epmach, erlarg=0.0, erlast, errbnd, errmax;
    double  error1, error2, erro12, errsum, ertest=0.0, oflow, resabs, reseps; 
    double  res3la[4], rlist2[53], uflow;
    int 	id, ierro, iroff1, iroff2, iroff3, jupbnd, k=1, ksgn;
    int		ktmin, maxerr, nres, nrmax, numrl2;
    int 	extrap, noext;
    int		crap = 0;
    int     i, j, nint, nintp1, ip1, jlow, ind1, ind2, levmax, levcur;


	/*            THE DIMENSION OF RLIST2 IS DETERMINED BY THE VALUE OF
	 *            LIMEXP IN SUBROUTINE EPSALG (RLIST2 SHOULD BE OF DIMENSION
	 *            (LIMEXP+2) AT LEAST).
	 *
	 *
	 *            LIST OF MAJOR VARIABLES
	 *            -----------------------
	 *
	 *           ALIST     - LIST OF LEFT END POINTS OF ALL SUBINTERVALS
	 *                       CONSIDERED UP TO NOW
	 *           BLIST     - LIST OF RIGHT END POINTS OF ALL SUBINTERVALS
	 *                       CONSIDERED UP TO NOW
	 *           RLIST(I)  - APPROXIMATION TO THE INTEGRAL OVER
	 *                       (ALIST(I),BLIST(I))
	 *           RLIST2    - ARRAY OF DIMENSION AT LEAST LIMEXP+2
	 *                       CONTAINING THE PART OF THE EPSILON TABLE WHICH
	 *                       IS STILL NEEDED FOR FURTHER COMPUTATIONS
	 *           ELIST(I)  - ERROR ESTIMATE APPLYING TO RLIST(I)
	 *           MAXERR    - POINTER TO THE INTERVAL WITH LARGEST ERROR
	 *                       ESTIMATE
	 *           ERRMAX    - ELIST(MAXERR)
	 *           ERLAST    - ERROR ON THE INTERVAL CURRENTLY SUBDIVIDED
	 *                       (BEFORE THAT SUBDIVISION HAS TAKEN PLACE)
	 *           AREA      - SUM OF THE INTEGRALS OVER THE SUBINTERVALS
	 *           ERRSUM    - SUM OF THE ERRORS OVER THE SUBINTERVALS
	 *           ERRBND    - REQUESTED ACCURACY MAX(EPSABS,EPSREL*
	 *                       ABS(RESULT))
	 *           *****1    - VARIABLE FOR THE LEFT SUBINTERVAL
	 *           *****2    - VARIABLE FOR THE RIGHT SUBINTERVAL
	 *           LAST      - INDEX FOR SUBDIVISION
	 *           NRES      - NUMBER OF CALLS TO THE EXTRAPOLATION ROUTINE
	 *           NUMRL2    - NUMBER OF ELEMENTS IN RLIST2. IF AN APPROPRIATE
	 *                       APPROXIMATION TO THE COMPOUNDED INTEGRAL HAS
	 *                       BEEN OBTAINED, IT IS PUT IN RLIST2(NUMRL2) AFTER
	 *                       NUMRL2 HAS BEEN INCREASED BY ONE.
	 *           ERLARG    - SUM OF THE ERRORS OVER THE INTERVALS LARGER
	 *                       THAN THE SMALLEST INTERVAL CONSIDERED UP TO NOW
	 *           EXTRAP    - LOGICAL VARIABLE DENOTING THAT THE ROUTINE
	 *                       IS ATTEMPTING TO PERFORM EXTRAPOLATION. I.E.
	 *                       BEFORE SUBDIVIDING THE SMALLEST INTERVAL WE
	 *                       TRY TO DECREASE THE VALUE OF ERLARG.
	 *           NOEXT     - LOGICAL VARIABLE DENOTING THAT EXTRAPOLATION IS
	 *                       NO LONGER ALLOWED (TRUE-VALUE)
	 *
	 *            MACHINE DEPENDENT CONSTANTS
	 *            ---------------------------
	 *
	 *           EPMACH IS THE LARGEST RELATIVE SPACING.
	 *           UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
	 *           OFLOW IS THE LARGEST POSITIVE MAGNITUDE.
     */

    int     npts, sign;
    double  temp, resa;






    /*  first executable statement  dqagpe */
    epmach = d1mach(4);




    /*
     *             test on validity of parameters
     */
    *ier = 0;
    *neval = 0;
    *last = 0;
    *result = 0.0;
    *abserr = 0.0;
    alist[1] = a;
    blist[1] = b;
    rlist[1] = 0.0;
    elist[1] = 0.0;
    iord[1]  = 0;
    level[1] = 0;
    npts = npts2 - 2;


    if ( (npts2 < 2) || (limit <= npts) || ( (epsabs <= 0.0) && (epsrel <  dmax1(50.0*epmach, 0.5e-28))) ) {
	    *ier = 6;
	    return(-1);
    }


    /*
     *            IF ANY BREAK POINTS ARE PROVIDED, SORT THEM INTO AN
     *            ASCENDING SEQUENCE.
     */
    sign   = ( a > b ) ? -1.0 : 1.0;
    pts[1] = dmin1(a, b);
    if ( npts != 0 ){
        for (i=1; i<=npts; i++) pts[i+1] = points[i];
    }
    

    pts[npts+2] = dmax1(a,b);
    nint = npts+1;
    a1   = pts[1];

    if ( npts != 0 ) {
        nintp1 = nint+1;
        for (i=1; i<=nint; i++){
            ip1 = i+1;
            for (j=ip1; j<=nintp1; j++){
                if ( pts[i] > pts[j] ) {
                    temp   = pts[i];
                    pts[i] = pts[j];
                    pts[j] = temp;
                }

            }
        }

        if ( (pts[1] != dmin1(a,b)) || ( pts[nintp1] != dmax1(a,b) ) ) {
            *ier = 6;
            return(-1);
        }

    }




    /*
     *           Compute first integral and error approximations.
     *           ------------------------------------------------
     */
    resabs = 0.0;
    for (i=1; i<=nint; i++){
        b1 = pts[i+1];
        dqk21(f, qpInfo, a1, b1, &area1, &error1, &defabs, &resa);
        *abserr += error1;
        *result += area1;
        ndin[i] = 0;
        if ( (error1 == resa) && (error1 != 0.0) ) ndin[i] = 1;
        resabs += defabs;
        level[i] = 0;
        elist[i] = error1;
        alist[i] = a1;
        blist[i] = b1;
        rlist[i] = area1;
        iord[i] = i;
        a1 = b1;
    }
    errsum = 0.0;
    for (i=1; i<=nint; i++){
        if ( ndin[i] == 1) elist[i] = *abserr;
        errsum += elist[i];
    }




    /*
     *           test on accuracy.
     */
    *last = nint;
    *neval = 21*nint;
    dres = fabs(*result);
    errbnd = dmax1(epsabs, epsrel*dres);

    if (  (*abserr <= 100.0*epmach*defabs) && (*abserr > errbnd)  ) *ier = 2;
    if ( nint != 1) {

        for (i=1; i<=npts; i++){

            jlow = i+1;
            ind1 = iord[i];
            for ( j=jlow; j<=nint; j++){
                ind2 = iord[j];
                if ( elist[ind1] <= elist[ind2]) {
                    ind1 = ind2;
                    k    = j;
                }
            }

            if ( ind1 != iord[i] ) {
                iord[k] = iord[i];
                iord[i] = ind1;
            }
        }

        if ( limit < npts2 ) *ier = 1;

    } else if ( (*ier != 0) || (*abserr <= errbnd) ) {
        return(-1);
    }



    /*
     *           initialization
     *           --------------
     */
    rlist2[1] = *result;
    maxerr = iord[1];
    errmax = elist[maxerr];
    area   = *result;
    nrmax  = 1;
    nres   = 0;
    numrl2 = 1;
    ktmin  = 0;
    extrap = FALSE;
    noext  = FALSE;
    erlarg = errsum;
    ertest = errbnd;
    levmax = 1;
    iroff1 = 0;
    iroff2 = 0;
    iroff3 = 0;
    ierro  = 0;
    uflow  = d1mach(1);
    oflow  = d1mach(2);
    *abserr = oflow;
    ksgn   = -1;
    if(dres >= (1.0 - 50.0*epmach)*defabs) ksgn = 1;





    /*
     *           main do-loop
     *           ------------
     */
    for (*last = npts2; *last <= limit; ++(*last)) {


        /*
         *           bisect the subinterval with the nrmax-th largest error
         *           estimate.
         */
        levcur = level[maxerr]+1;
        a1 = alist[maxerr];
        b1 = 0.5*(alist[maxerr] + blist[maxerr]);
        a2 = b1;
        b2 = blist[maxerr];
        erlast = errmax;
        dqk21(f, qpInfo, a1, b1, &area1, &error1, &resa, &defab1);
        dqk21(f, qpInfo, a2, b2, &area2, &error2, &resa, &defab2);



        /*
         *           improve previous approximations to integral
         *           and error and test for accuracy.
         */
        neval  += 42;
        area12 = area1+area2;
        erro12 = error1+error2;
        errsum += erro12-errmax;
        area   += area12-rlist[maxerr];


        if (  (defab1 != error1) && (defab2 != error2)  ) {

            if (  (fabs(rlist[maxerr]-area12) <= 0.1e-4*fabs(area12)) && (erro12 >= 0.99*errmax)  ) { 
                if (extrap) ++iroff2;
                if (!extrap) ++iroff1;
	        }
	        if ( (*last > 10) && (erro12 > errmax) ) ++iroff3;

	    }


        level[maxerr] = levcur;
        level[*last]   = levcur;
	    rlist[maxerr] = area1;
        rlist[*last]  = area2;
        errbnd = dmax1( epsabs, epsrel*fabs(area) );




        /*
         *           test for roundoff error and eventually set error flag.
         */
        if ( (iroff1+iroff2 >= 10) || (iroff3 >= 20) ) *ier = 2;
        if (iroff2 >= 5) ierro = 3;



        /*
         *           set error flag in the case that the number of subintervals
         *           equals limit.
         */
        if (*last == limit) *ier = 1;



        /*
         *           set error flag in the case of bad integrand behaviour
         *           at a point of the integration range.
         */
        if ( dmax1(fabs(a1), fabs(b2)) <= ((1.0+100.0*epmach)*(fabs(a2)+1000.0*uflow)) ) *ier = 4;



        /*
         *           append the newly-created intervals to the list.
         */
        if(error2 > error1) { 

            alist[maxerr] = a2;
            alist[*last]  = a1;
            blist[*last]  = b1;
            rlist[maxerr] = area2;
            rlist[*last]  = area1;
            elist[maxerr] = error2;
            elist[*last]  = error1;

	    } else {

            alist[*last]  = a2;
            blist[maxerr] = b1;
            blist[*last]  = b2;
            elist[maxerr] = error1;
            elist[*last]  = error2;

	    }






        /*
         *           call subroutine dqpsrt to maintain the descending ordering
         *           in the list of error estimates and select the subinterval
         *           with nrmax-th largest error estimate (to be bisected next).
         */
        dqpsrt(limit, *last, &maxerr, &errmax, elist, iord, &nrmax);



        /*  
         *  Jump out of do-loop
         */
        if (errsum <= errbnd) goto L190;


        /*  
         *  Jump out of do-loop
         */
        if (*ier != 0) break;

        if ( !noext ) {

            erlarg -= erlast;
            if ( (levcur+1) <= levmax) erlarg += erro12;

            if ( !extrap ) {
    	        /*
	             *           test whether the interval to be bisected next is the
	             *           smallest interval.
	             */
                if ( (level[maxerr]+1) <= levmax) goto L160;
                extrap = TRUE;
                nrmax  = 2;
	        }



   	        if( (ierro != 3) && (erlarg > ertest) ) {

                /*
                 *           the smallest interval has the largest error.
                 *           before bisecting decrease the sum of the errors over the
                 *           larger intervals (erlarg) and perform extrapolation.
                 */
                id = nrmax;
                jupbnd = *last;
                if(*last > (2+limit/2)) jupbnd = limit+3- (*last);

		        for (k = id; k<=jupbnd; ++k) {

		            maxerr = iord[nrmax];
		            errmax = elist[maxerr];

                    /*  
                     *  Jump out of do-loop
                     */
                    if( (level[maxerr]+1) <= levmax) goto L160;
                    ++nrmax;

		        }

	        }




            /*
             *           perform extrapolation.
             */
            ++numrl2;
            rlist2[numrl2] = area;
            if ( numrl2 <= 2 ) goto L155;

            dqelg(numrl2, rlist2, &reseps, &abseps, res3la, &nres);
            ++ktmin;
            if( (ktmin > 5) && (*abserr < 0.1e-2*errsum) ) *ier = 5;
	        if( abseps >= *abserr ) goto L150;
            ktmin = 0;
		    *abserr = abseps;
		    *result = reseps;
		    correc = erlarg;
		    ertest = dmax1( epsabs, epsrel*fabs(reseps) );

		    /*  
		     *  Jump out of do-loop
		     */
		    if( *abserr <= ertest ) goto L170;




            /*
             *           prepare bisection of the smallest interval.
             */
L150:
            if ( numrl2 == 1 ) noext = TRUE;
            if ( *ier == 5 ) goto L170;
L155:
            maxerr = iord[1];
            errmax = elist[maxerr];
            nrmax = 1;
            extrap = FALSE;
            ++levmax;
            erlarg = errsum;

	    }
L160:
        crap = 0;

	}




L170:




    /*
     *           set final result and error estimate.
     *           ------------------------------------
     *   I dont have time to pretty up all these ugly goto's right now! -- MGH
     */
    if (*abserr == oflow) goto L190;
    if (*ier+ierro == 0) goto L180;
    if (ierro == 3) *abserr += correc;
    if (*ier == 0) *ier = 3;
    if ( (*result != 0.0) && (area != 0.0) ) goto L175;
    if (*abserr > errsum) goto L190;
    if (area == 0.0) goto L210;
    goto L180;

L175:
    if(*abserr/fabs(*result) > errsum/fabs(area)) goto L190;



    /*
     *           test on divergence.
     */
L180:
    if(ksgn == (-1) && dmax1(fabs(*result), fabs(area)) <=  defabs*0.01) goto L210;
    if( (0.01 > (*result/area)) || ((*result/area) > 100.0) || (errsum > fabs(area))  ) *ier = 6;
    goto L210;



    /*
     *           compute global integral sum.
     */
L190:
    *result = 0.0;
    for (k=1; k<= (*last); ++k) *result += rlist[k];
    *abserr = errsum;

L210:
    if (*ier > 2) --(*ier);
    *result *= sign;

    return(0);


}
