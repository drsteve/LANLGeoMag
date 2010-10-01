#include "Lgm/Lgm_QuadPack.h"

/*
 *      QUADPACK DQAGI Routine converted to C 
 */
int dqagi(f, qpInfo, bound, inf, epsabs, epsrel, result, abserr, neval, ier, limit, lenw, last, iwork, work)
double  (*f)( double, _qpInfo *); /*  The integrand function -- I.e. the function to integrate    	 */
_qpInfo *qpInfo;        	  /*  Auxilliary information to pass to function (to avoid making globals) */
double   bound;               /*  One limit of integration.                                 	 */
int      inf;             	  /*  see below                                                      */
double   epsabs;        	  /*  Absolute accuracy requested.                                	 */
double   epsrel;        	  /*  Relative accuracy requested.                                	 */
double  *result;        	  /*  The desired result. I.e. integral of f() from a to b        	 */
double  *abserr;        	  /*  Estimate of the modulus of the absolute error in the result 	 */
int     *neval;         	  /*  The number of integrand evaluations performed.              	 */
int     *ier;           	  /*  Error flag. An error occurred if ier > 0. See below.        	 */
int	    limit;
int	    lenw;
int	   *last;
int	   *iwork;
double *work;
{

    /* ***BEGIN PROLOGUE  DQAGI
     * ***PURPOSE  The routine calculates an approximation result to a given
     *             INTEGRAL   I = Integral of F over (BOUND,+INFINITY)
     *             OR I = Integral of F over (-INFINITY,BOUND)
     *             OR I = Integral of F over (-INFINITY,+INFINITY)
     *             Hopefully satisfying following claim for accuracy
     *             ABS(I-RESULT).LE.MAX(EPSABS,EPSREL*ABS(I)).
     * ***LIBRARY   SLATEC (QUADPACK)
     * ***CATEGORY  H2A3A1, H2A4A1
     * ***TYPE      DOUBLE PRECISION (QAGI-S, DQAGI-D)
     * ***KEYWORDS  AUTOMATIC INTEGRATOR, EXTRAPOLATION, GENERAL-PURPOSE,
     *              GLOBALLY ADAPTIVE, INFINITE INTERVALS, QUADPACK,
     *              QUADRATURE, TRANSFORMATION
     * ***AUTHOR  Piessens, Robert
     *              Applied Mathematics and Programming Division
     *              K. U. Leuven
     *            de Doncker, Elise
     *              Applied Mathematics and Programming Division
     *              K. U. Leuven
     * ***DESCRIPTION
     * 
     *         Integration over infinite intervals
     *         Standard fortran subroutine
     * 
     *         PARAMETERS
     *          ON ENTRY
     *             F      - Double precision
     *                      Function subprogram defining the integrand
     *                      function F(X). The actual name for F needs to be
     *                      declared E X T E R N A L in the driver program.
     * 
     *             BOUND  - Double precision
     *                      Finite bound of integration range
     *                      (has no meaning if interval is doubly-infinite)
     * 
     *             INF    - Integer
     *                      indicating the kind of integration range involved
     *                      INF = 1 corresponds to  (BOUND,+INFINITY),
     *                      INF = -1            to  (-INFINITY,BOUND),
     *                      INF = 2             to (-INFINITY,+INFINITY).
     * 
     *             EPSABS - Double precision
     *                      Absolute accuracy requested
     *             EPSREL - Double precision
     *                      Relative accuracy requested
     *                      If  EPSABS.LE.0
     *                      and EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28),
     *                      the routine will end with IER = 6.
     * 
     * 
     *          ON RETURN
     *             RESULT - Double precision
     *                      Approximation to the integral
     * 
     *             ABSERR - Double precision
     *                      Estimate of the modulus of the absolute error,
     *                      which should equal or exceed ABS(I-RESULT)
     * 
     *             NEVAL  - Integer
     *                      Number of integrand evaluations
     * 
     *             IER    - Integer
     *                      IER = 0 normal and reliable termination of the
     *                              routine. It is assumed that the requested
     *                              accuracy has been achieved.
     *                    - IER.GT.0 abnormal termination of the routine. The
     *                              estimates for result and error are less
     *                              reliable. It is assumed that the requested
     *                              accuracy has not been achieved.
     *             ERROR MESSAGES
     *                      IER = 1 Maximum number of subdivisions allowed
     *                              has been achieved. One can allow more
     *                              subdivisions by increasing the value of
     *                              LIMIT (and taking the according dimension
     *                              adjustments into account). However, if
     *                              this yields no improvement it is advised
     *                              to analyze the integrand in order to
     *                              determine the integration difficulties. If
     *                              the position of a local difficulty can be
     *                              determined (e.g. SINGULARITY,
     *                              DISCONTINUITY within the interval) one
     *                              will probably gain from splitting up the
     *                              interval at this point and calling the
     *                              integrator on the subranges. If possible,
     *                              an appropriate special-purpose integrator
     *                              should be used, which is designed for
     *                              handling the type of difficulty involved.
     *                          = 2 The occurrence of roundoff error is
     *                              detected, which prevents the requested
     *                              tolerance from being achieved.
     *                              The error may be under-estimated.
     *                          = 3 Extremely bad integrand behaviour occurs
     *                              at some points of the integration
     *                              interval.
     *                          = 4 The algorithm does not converge.
     *                              Roundoff error is detected in the
     *                              extrapolation table.
     *                              It is assumed that the requested tolerance
     *                              cannot be achieved, and that the returned
     *                              RESULT is the best which can be obtained.
     *                          = 5 The integral is probably divergent, or
     *                              slowly convergent. It must be noted that
     *                              divergence can occur with any other value
     *                              of IER.
     *                          = 6 The input is invalid, because
     *                              (EPSABS.LE.0 and
     *                               EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28))
     *                               or LIMIT.LT.1 or LENIW.LT.LIMIT*4.
     *                              RESULT, ABSERR, NEVAL, LAST are set to
     *                              zero.  Except when LIMIT or LENIW is
     *                              invalid, IWORK(1), WORK(LIMIT*2+1) and
     *                              WORK(LIMIT*3+1) are set to ZERO, WORK(1)
     *                              is set to A and WORK(LIMIT+1) to B.
     * 
     *          DIMENSIONING PARAMETERS
     *             LIMIT - Integer
     *                     Dimensioning parameter for IWORK
     *                     LIMIT determines the maximum number of subintervals
     *                     in the partition of the given integration interval
     *                     (A,B), LIMIT.GE.1.
     *                     If LIMIT.LT.1, the routine will end with IER = 6.
     * 
     *             LENW  - Integer
     *                     Dimensioning parameter for WORK
     *                     LENW must be at least LIMIT*4.
     *                     If LENW.LT.LIMIT*4, the routine will end
     *                     with IER = 6.
     * 
     *             LAST  - Integer
     *                     On return, LAST equals the number of subintervals
     *                     produced in the subdivision process, which
     *                     determines the number of significant elements
     *                     actually in the WORK ARRAYS.
     * 
     *          WORK ARRAYS
     *             IWORK - Integer
     *                     Vector of dimension at least LIMIT, the first
     *                     K elements of which contain pointers
     *                     to the error estimates over the subintervals,
     *                     such that WORK(LIMIT*3+IWORK(1)),... ,
     *                     WORK(LIMIT*3+IWORK(K)) form a decreasing
     *                     sequence, with K = LAST if LAST.LE.(LIMIT/2+2), and
     *                     K = LIMIT+1-LAST otherwise
     * 
     *             WORK  - Double precision
     *                     Vector of dimension at least LENW
     *                     on return
     *                     WORK(1), ..., WORK(LAST) contain the left
     *                      end points of the subintervals in the
     *                      partition of (A,B),
     *                     WORK(LIMIT+1), ..., WORK(LIMIT+LAST) Contain
     *                      the right end points,
     *                     WORK(LIMIT*2+1), ...,WORK(LIMIT*2+LAST) contain the
     *                      integral approximations over the subintervals,
     *                     WORK(LIMIT*3+1), ..., WORK(LIMIT*3)
     *                      contain the error estimates.
     * 
     * ***REFERENCES  (NONE)
     * ***ROUTINES CALLED  DQAGIE, XERMSG
     * ***REVISION HISTORY  (YYMMDD)
     *    800101  DATE WRITTEN
     *    890831  Modified array declarations.  (WRB)
     *    890831  REVISION DATE from Version 3.2
     *    891214  Prologue converted to Version 4.0 format.  (BAB)
     *    900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
     * ***END PROLOGUE  DQAGI
     *    20071127 Converted to C  M. Henderson
     */





    int 	lvl=0, l1, l2, l3;

/*
double TMPwork[600], TMPwork2[600], TMPwork3[600], TMPwork4[600];
int TMPiwork[600];
*/


    /* First executable statement  dqags */


    /*
     *    Check validity of limit and lenw.
     */
    *ier = 6;
    *neval = 0;
    *last = 0;
    *result = 0.0;
    *abserr = 0.0;

    if( (limit < 1) || (lenw < limit*4) ) {

        if (*ier == 6) lvl = 1;
        if (*ier != 0) fprintf(stderr, "Abnormal return from dqagi, ier = %d,  lvl = %d\n", *ier, lvl);
        return(-1);

    } else {

        /*
         *    Prepare call for dqagse.
         */
        l1 = limit+1;
        l2 = limit+l1;
        l3 = limit+l2;


        dqagie(f, qpInfo, bound, inf, epsabs, epsrel, limit, result, abserr, neval, ier, 
				work, work+l1, work+l2, work+l3, iwork, last);
    }



    /*
     *    Call error handler if necessary.
     */
    lvl = 0;


    if(*ier == 6) lvl = 1;
    if(*ier != 0){ 
	fprintf(stderr, "Abnormal return from dqags, ier = %d,  lvl = %d\n", *ier, lvl);
        return(-1);
    } else {
        return(1);
    }


}




/*
 *      QUADPACK DQAGIE Routine converted to C
 */
int dqagse(f, qpInfo, bound, inf, epsabs, epsrel, limit, result, abserr, neval, ier, alist, blist, rlist, elist, iord, last)
double  (*f)( double, _qpInfo *); /*  The integrand function -- I.e. the function to integrate           */
_qpInfo *qpInfo;                  /*  Auxilliary information to pass to function (to avoid making globals) */
double   a;             	  /*  Lower Limit of integration.                                 */
double   b;             	  /*  Upper limit of integration.                                 */
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
int	 *iord;
int	*last;
{

    /*   Begin prologue:  dqagie
     *   Purpose: The routine calculates an approximation result to a given
     *            integral   I = Integral of F over (BOUND,+INFINITY)
     *            or I = Integral of F over (-INFINITY,BOUND)
     *            or I = Integral of F over (-INFINITY,+INFINITY),
     *            hopefully satisfying following claim for accuracy
     *            ABS(I-RESULT).LE.MAX(EPSABS,EPSREL*ABS(I))
     *
     *   Category  H2A3A1, H2A4A1
     *   Type      double precision (qagie-s, dqagie-d)
     *   Keywords  automatic integrator, extrapolation, general-purpose,
     *             globally adaptive, infinite intervals, quadpack,
     *             quadrature, transformatioN
     *   Author  Piessens, Robert
     *             Applied Mathematics and Programming Division
     *             K. U. Leuven
     *           de Doncker, Elise
     *             Applied Mathematics and Programming Division
     *             K. U. Leuven
     *   Description
     *
     * Integration over infinite intervals
     * Standard fortran subroutine
     *
     *            f      - Double precision
     *                     Function subprogram defining the integrand
     *                     function F(X). The actual name for F needs to be
     *                     declared E X T E R N A L in the driver program.
     *
     *            bound  - Double precision
     *                     Finite bound of integration range
     *                     (has no meaning if interval is doubly-infinite)
     *
     *            inf    - Double precision
     *                     Indicating the kind of integration range involved
     *                     INF = 1 corresponds to  (BOUND,+INFINITY),
     *                     INF = -1            to  (-INFINITY,BOUND),
     *                     INF = 2             to (-INFINITY,+INFINITY).
     *
     *            epsabs - Double precision
     *                     Absolute accuracy requested
     *            epsrel - Double precision
     *                     Relative accuracy requested
     *                     If  EPSABS.LE.0
     *                     and EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28),
     *                     the routine will end with IER = 6.
     *
     *            limit  - Integer
     *                     Gives an upper bound on the number of subintervals
     *                     in the partition of (A,B), LIMIT.GE.1
     *
     *         on return
     *            result - Double precision
     *                     Approximation to the integral
     *
     *            abserr - Double precision
     *                     Estimate of the modulus of the absolute error,
     *                     which should equal or exceed ABS(I-RESULT)
     *
     *            neval  - Integer
     *                     Number of integrand evaluations
     *
     *            ier    - Integer
     *                     ier == 0 Normal and reliable termination of the
     *                              routine. It is assumed that the requested
     *                              accuracy has been achieved.
     *                   - ier > 0  Abnormal termination of the routine. The
     *                              estimates for result and error are less
     *                              reliable. It is assumed that the requested
     *                              accuracy has not been achieved.
     *            Error messages
     *                     ier = 1 Maximum number of subdivisions allowed
     *                             has been achieved. One can allow more
     *                             subdivisions by increasing the value of
     *                             LIMIT (and taking the according dimension
     *                             adjustments into account).  However, if
     *                             this yields no improvement it is advised
     *                             to analyze the integrand in order to
     *                             determine the integration difficulties.
     *                             If the position of a local difficulty can
     *                             be determined (e.g. SINGULARITY,
     *                             DISCONTINUITY within the interval) one
     *                             will probably gain from splitting up the
     *                             interval at this point and calling the
     *                             integrator on the subranges. If possible,
     *                             an appropriate special-purpose integrator
     *                             should be used, which is designed for
     *                             handling the type of difficulty involved.
     *                         = 2 The occurrence of roundoff error is
     *                             detected, which prevents the requested
     *                             tolerance from being achieved.
     *                             The error may be under-estimated.
     *                         = 3 Extremely bad integrand behaviour occurs
     *                             at some points of the integration
     *                             interval.
     *                         = 4 The algorithm does not converge.
     *                             Roundoff error is detected in the
     *                             extrapolation table.
     *                             It is assumed that the requested tolerance
     *                             cannot be achieved, and that the returned
     *                             result is the best which can be obtained.
     *                         = 5 The integral is probably divergent, or
     *                             slowly convergent. It must be noted that
     *                             divergence can occur with any other value
     *                             of IER.
     *                         = 6 The input is invalid, because
     *                             (EPSABS.LE.0 and
     *                              EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28),
     *                             RESULT, ABSERR, NEVAL, LAST, RLIST(1),
     *                             ELIST(1) and IORD(1) are set to zero.
     *                             ALIST(1) and BLIST(1) are set to 0
     *                             and 1 respectively.
     *
     *
     *            alist  - Double precision
     *                     Vector of dimension at least LIMIT, the first
     *                      LAST  elements of which are the left
     *                     end points of the subintervals in the partition
     *                     of the transformed integration range (0,1).
     *
     *            blist  - Double precision
     *                     Vector of dimension at least LIMIT, the first
     *                      LAST  elements of which are the right
     *                     end points of the subintervals in the partition
     *                     of the transformed integration range (0,1).
     *
     *            rlist  - Double precision
     *                     Vector of dimension at least LIMIT, the first
     *                      LAST  elements of which are the integral
     *                     approximations on the subintervals
     *
     *            elist  - Double precision
     *                     Vector of dimension at least LIMIT,  the first
     *                     LAST elements of which are the moduli of the
     *                     absolute error estimates on the subintervals
     *
     *            iord   - Integer
     *                     Vector of dimension LIMIT, the first K
     *                     elements of which are pointers to the
     *                     error estimates over the subintervals,
     *                     such that ELIST(IORD(1)), ..., ELIST(IORD(K))
     *                     form a decreasing sequence, with K = LAST
     *                     If LAST.LE.(LIMIT/2+2), and K = LIMIT+1-LAST
     *                     otherwise
     *
     *            last   - Integer
     *                     Number of subintervals actually produced
     *                     in the subdivision process
     *
     *   References  (none)
     *   Routines called  d1mach, dqelg, dqk15i, dqpsrt
     *   Revision history  (yymmdd)
     *   800101  DATE WRITTEN
     *   890531  Changed all specific intrinsics to generic.  (WRB)
     *   890831  Modified array declarations.  (WRB)
     *   890831  REVISION DATE from Version 3.2
     *   891214  Prologue converted to Version 4.0 format.  (BAB)
     *   20071127 Converted to C by M. Henderson
     *   end prologue:  dqagie
     */






    double 	area, abseps, area1, area12, area2, a1;
    double     	a2, b1, b2, correc=0.0, defabs, defab1, defab2, d1mach();
    double     	dres, epmach, erlarg=0.0, erlast, errbnd, errmax;
    double     	error1, error2, erro12, errsum, ertest=0.0, oflow, resabs, reseps; 
    double     	res3la[4], rlist2[53], small=0.0, uflow;
    int 	id, ierro, iroff1, iroff2, iroff3, jupbnd, k, ksgn;
    int		ktmin, maxerr, nres, nrmax, numrl2;
    int 	extrap, noext;
    int		crap = 0;


    /*  first executable statement  dqagie */
    epmach = d1mach(4);




    /*
     *             test on validity of parameters
     */
    *ier = 0;
    *neval = 0;
    *last = 0;
    *result = 0.0;
    *abserr = 0.0;
    alist[1] = 0.0;
    blist[1] = 1.0;
    rlist[1] = 0.0;
    elist[1] = 0.0;
    if( (epsabs <= 0.0) && (epsrel < dmax1(50.0*epmach, 0.5e-28))  ) {
	    *ier = 6;
	    return(-1);
    }




    /*
     *           first approximation to the integral
     *           -----------------------------------
     *           Determine the interval to be mapped onto (0,1).
     *           If inf = 2 the integral is computed as i = i1+i2, where
     *           i1 = integral of f over (-infinity,0),
     *           i2 = integral of f over (0,+infinity).
     */
    boun = bound;
    if (inf == 2) boun = 0.0;
    dqk15i( f, qpInfo, boun, inf, 0.0, 1.0, result, abserr, defabs, resabs );


    /*
     *           test on accuracy.
     */
    *last    = 1;
    rlist[1] = *result;
    elist[1] = *abserr;
    iord[1]  = 1;
    dres = fabs(*result);
    errbnd = dmax1(epsabs, epsrel*dres);
    if (  (*abserr <= 100.0*epmach*defabs) && (*abserr > errbnd)  ) *ier = 2;
    if (limit == 1) *ier = 1;
    if (  (*ier != 0) || ( (*abserr <= errbnd) && (*abserr != resabs)) ||  (*abserr == 0.0)  ) {
	    *neval = 30 * *last - 15;
        if (inf == 2) NEVAL *= 2;
        if (ier >  2) --ier;
	    return(0);
    }



    /*
     *           initialization
     *           --------------
     */
    uflow = d1mach(1);
    oflow = d1mach(2);
    rlist2[1] = *result;
    errmax = *abserr;
    maxerr = 1;
    area = *result;
    errsum = *abserr;
    *abserr = oflow;
    nrmax = 1;
    nres = 0;
    ktmin = 0;
    numrl2 = 2;
    extrap = FALSE;
    noext = FALSE;
    ierro = 0;
    iroff1 = 0;
    iroff2 = 0;
    iroff3 = 0;
    ksgn = -1;
    if(dres >= (1.0 - 50.0*epmach)*defabs) ksgn = 1;






    /*
     *           main do-loop
     *           ------------
     */
    for (*last = 2; *last <= limit; ++(*last)) {


        /*
         *           bisect the subinterval with the nrmax-th largest error
         *           estimate.
         */
        a1 = alist[maxerr];
        b1 = 0.5*(alist[maxerr] + blist[maxerr]);
        a2 = b1;
        b2 = blist[maxerr];
        erlast = errmax;
        dqk15i(f, qpInfo, boun, inf, a1, b1, &area1, &error1, &resabs, &defab1);
        dqk15i(f, qpInfo, boun, inf, a2, b2, &area2, &error2, &resabs, &defab2);



        /*
         *           improve previous approximations to integral
         *           and error and test for accuracy.
         */
        area12 = area1+area2;
        erro12 = error1+error2;
        errsum += erro12-errmax;
        area += area12-rlist[maxerr];


        if (  (defab1 != error1) && (defab2 != error2)  ) {

            if (  (fabs(rlist[maxerr]-area12) <= 0.1e-4*fabs(area12)) && (erro12 >= 0.99*errmax)  ) { 
                if (extrap) ++iroff2;
                if (!extrap) ++iroff1;
	        }
	        if ( (*last > 10) && (erro12 > errmax) ) ++iroff3;

	    }


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
            alist[*last] = a1;
            blist[*last] = b1;
            rlist[maxerr] = area2;
            rlist[*last] = area1;
            elist[maxerr] = error2;
            elist[*last] = error1;

	    } else {

            alist[*last] = a2;
            blist[maxerr] = b1;
            blist[*last] = b2;
            elist[maxerr] = error1;
            elist[*last] = error2;

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
        if (errsum <= errbnd) goto L115;


	    /*  
	     *  Jump out of do-loop
	     */
        if (*ier != 0) break;

        if (*last == 2) {
	        small = 0.375;
	        erlarg = errsum;
	        ertest = errbnd;
	        rlist2[2] = area;
	
	        goto L90;
	    }




        if ( !noext ) {


            erlarg -= erlast;
            if( fabs(b1-a1) > small ) erlarg += erro12;


            if ( !extrap ) {
    	        /*
	             *           test whether the interval to be bisected next is the
	             *           smallest interval.
	             */
                if ( fabs(blist[maxerr]-alist[maxerr]) > small ) goto L90;
                extrap = TRUE;
                nrmax = 2;
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
                    if ( fabs(blist[maxerr]-alist[maxerr]) > small ) goto L90;
                    ++nrmax;

                }
            }




            /*
             *           perform extrapolation.
             */
            ++numrl2;
            rlist2[numrl2] = area;
            dqelg(numrl2, rlist2, &reseps, &abseps, res3la, &nres);
            ++ktmin;
            if( (ktmin > 5) && (*abserr < 0.1e-2*errsum) ) *ier = 5;




            if( abseps < *abserr ) {
                ktmin = 0;
                *abserr = abseps;
                *result = reseps;
                correc = erlarg;
                ertest = dmax1( epsabs, epsrel*fabs(reseps) );
                /*  
                 *  Jump out of do-loop
                 */
                if( *abserr <= ertest ) break;
            }



            /*
             *           prepare bisection of the smallest interval.
             */
            if ( numrl2 == 1 ) noext = TRUE;
            if ( *ier == 5 ) break;
            maxerr = iord[1];
            errmax = elist[maxerr];
            nrmax = 1;
            extrap = FALSE;
            small *= 0.5;
            erlarg = errsum;

	    }


L90:
	    crap = 0;

    }










    /*
     *           set final result and error estimate.
     *           ------------------------------------
     *   I dont have time to pretty up all these ugly goto's right now! -- MGH
     */
    if (*abserr == oflow) goto L115;
    if (*ier+ierro == 0) goto L110;
    if (ierro == 3) *abserr += correc;
    if (*ier == 0) *ier = 3;
    if ( (*result != 0.0) && (area != 0.0) ) goto L105;
    if (*abserr > errsum) goto L115;
    if (area == 0.0) goto L130;
    goto L110;

L105:
    if(*abserr/fabs(*result) > errsum/fabs(area)) goto L115;



    /*
     *           test on divergence.
     */
L110:
    if(ksgn == (-1) && dmax1(fabs(*result), fabs(area)) <=  defabs*0.01) goto L130;
    if( (0.01 > (*result/area)) || ((*result/area) > 100.0) || (errsum > fabs(area))  ) *ier = 6;
    goto L130;



    /*
     *           compute global integral sum.
     */
L115:
    *result = 0.0;
    for (k=1; k<= (*last); ++k) *result += rlist[k];
    *abserr = errsum;

L130:
    *neval = 30*(*last) - 15;
    if (inf == 2) neval *= 2;
    if (*ier > 2) --(*ier);
//L140:

    return(0);


}


/*
 *   $Id$
 */

