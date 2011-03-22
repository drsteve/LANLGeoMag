#include "Lgm/Lgm_QuadPack.h"

/*
 *      QUADPACK DQAGS Routine converted to C 
 */
int dqags(f, qpInfo, a, b, epsabs, epsrel, result, abserr, neval, ier, limit, lenw, last, iwork, work)
double  (*f)( double, _qpInfo *); /*  The integrand function -- I.e. the function to integrate    	 */
_qpInfo *qpInfo;        	  /*  Auxilliary information to pass to function (to avoid making globals) */
double   a;             	  /*  Lower Limit of integration.                                 	 */
double   b;             	  /*  Upper limit of integration.                                 	 */
double   epsabs;        	  /*  Absolute accuracy requested.                                	 */
double   epsrel;        	  /*  Relative accuracy requested.                                	 */
double  *result;        	  /*  The desired result. I.e. integral of f() from a to b        	 */
double  *abserr;        	  /*  Estimate of the modulus of the absolute error in the result 	 */
int     *neval;         	  /*  The number of integrand evaluations performed.              	 */
int     *ier;           	  /*  Error flag. An error occurred if ier > 0. See below.        	 */
int	 limit;
int	 lenw;
int	*last;
int	*iwork;
double	*work;
{



    /*
     *
     *    Begin prologue:   	dqags
     *    Date written:   	800101   (yymmdd)
     *    Revision date:  	830518   (yymmdd)
     *    Category no.:  	h2a1a1
     *
     *    Keywords:  		automatic integrator, general-purpose,
     *              		(end-point) singularities, extrapolation,
     *              		globally adaptive
     *
     *    Author:  		piessens,robert,appl. math. & progr. div. - k.u.leuven
     *            		de doncker,elise,appl. math. & prog. div. - k.u.leuven
     *
     *    Purpose:  		the routine calculates an approximation result to a given
     *             		definite integral  i = integral of f over (a,b),
     *             		hopefully satisfying following claim for accuracy
     *             		abs(i-result).le.max(epsabs,epsrel*abs(i)).
     *
     *    Description: 		Computation of a definite integral
     *        			standard fortran subroutine
     *        			double precision version
     *
     *
     *    Parameters:
     *         on entry
     *            f      - double precision
     *                     function subprogram defining the integrand
     *                     function f(x). the actual name for f needs to be
     *                     declared e x t e r n a l in the driver program.
     *
     *            a      - double precision
     *                     lower limit of integration
     *
     *            b      - double precision
     *                     upper limit of integration
     *
     *            epsabs - double precision
     *                     absolute accuracy requested
     *            epsrel - double precision
     *                     relative accuracy requested
     *                     if  epsabs.le.0
     *                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
     *                     the routine will end with ier = 6.
     *
     *         on return
     *            result - double precision
     *                     approximation to the integral
     *
     *            abserr - double precision
     *                     estimate of the modulus of the absolute error,
     *                     which should equal or exceed abs(i-result)
     *
     *            neval  - integer
     *                     number of integrand evaluations
     *
     *            ier    - integer
     *                     ier = 0 normal and reliable termination of the
     *                             routine. it is assumed that the requested
     *                             accuracy has been achieved.
     *                     ier.gt.0 abnormal termination of the routine
     *                             the estimates for integral and error are
     *                             less reliable. it is assumed that the
     *                             requested accuracy has not been achieved.
     *            error messages
     *                     ier = 1 maximum number of subdivisions allowed
     *                             has been achieved. one can allow more sub-
     *                             divisions by increasing the value of limit
     *                             (and taking the according dimension
     *                             adjustments into account. however, if
     *                             this yields no improvement it is advised
     *                             to analyze the integrand in order to
     *                             determine the integration difficulties. if
     *                             the position of a local difficulty can be
     *                             determined (e.g. singularity,
     *                             discontinuity within the interval) one
     *                             will probably gain from splitting up the
     *                             interval at this point and calling the
     *                             integrator on the subranges. if possible,
     *                             an appropriate special-purpose integrator
     *                             should be used, which is designed for
     *                             handling the type of difficulty involved.
     *                         = 2 the occurrence of roundoff error is detec-
     *                             ted, which prevents the requested
     *                             tolerance from being achieved.
     *                             the error may be under-estimated.
     *                         = 3 extremely bad integrand behaviour
     *                             occurs at some points of the integration
     *                             interval.
     *                         = 4 the algorithm does not converge.
     *                             roundoff error is detected in the
     *                             extrapolation table. it is presumed that
     *                             the requested tolerance cannot be
     *                             achieved, and that the returned result is
     *                             the best which can be obtained.
     *                         = 5 the integral is probably divergent, or
     *                             slowly convergent. it must be noted that
     *                             divergence can occur with any other value
     *                             of ier.
     *                         = 6 the input is invalid, because
     *                             (epsabs.le.0 and
     *                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28)
     *                             or limit.lt.1 or lenw.lt.limit*4.
     *                             result, abserr, neval, last are set to
     *                             zero.except when limit or lenw is invalid,
     *                             iwork(1), work(limit*2+1) and
     *                             work(limit*3+1) are set to zero, work(1)
     *                             is set to a and work(limit+1) to b.
     *
     *         dimensioning parameters
     *            limit - integer
     *                    dimensioning parameter for iwork
     *                    limit determines the maximum number of subintervals
     *                    in the partition of the given integration interval
     *                    (a,b), limit.ge.1.
     *                    if limit.lt.1, the routine will end with ier = 6.
     *
     *            lenw  - integer
     *                    dimensioning parameter for work
     *                    lenw must be at least limit*4.
     *                    if lenw.lt.limit*4, the routine will end
     *                    with ier = 6.
     *
     *            last  - integer
     *                    on return, last equals the number of subintervals
     *                    produced in the subdivision process, detemines the
     *                    number of significant elements actually in the work
     *                    arrays.
     *
     *         work arrays
     *            iwork - integer
     *                    vector of dimension at least limit, the first k
     *                    elements of which contain pointers
     *                    to the error estimates over the subintervals
     *                    such that work(limit*3+iwork(1)),... ,
     *                    work(limit*3+iwork(k)) form a decreasing
     *                    sequence, with k = last if last.le.(limit/2+2),
     *                    and k = limit+1-last otherwise
     *
     *            work  - double precision
     *                    vector of dimension at least lenw
     *                    on return
     *                    work(1), ..., work(last) contain the left
     *                     end-points of the subintervals in the
     *                     partition of (a,b),
     *                    work(limit+1), ..., work(limit+last) contain
     *                     the right end-points,
     *                    work(limit*2+1), ..., work(limit*2+last) contain
     *                     the integral approximations over the subintervals,
     *                    work(limit*3+1), ..., work(limit*3+last)
     *                     contain the error estimates.
     *
     *     references  (none)
     *     routines called  dqagse,xerror
     *     end prologue  dqags
     *
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
        if (*ier != 0) fprintf(stderr, "dqags: Abnormal return from dqags, ier = %d,  lvl = %d\n", *ier, lvl);
        return(-1);

    } else {

        /*
         *    Prepare call for dqagse.
         */
        l1 = limit+1;
        l2 = limit+l1;
        l3 = limit+l2;


        dqagse(f, qpInfo, a, b, epsabs, epsrel, limit, result, abserr, neval, ier, 
				work, work+l1, work+l2, work+l3, iwork, last);

    }



    /*
     *    Call error handler if necessary.
     */
    lvl = 0;


    if(*ier == 6) lvl = 1;
    if(*ier != 0){ 
	fprintf(stderr, "dqags: Abnormal return from dqags, ier = %d,  lvl = %d\n", *ier, lvl);
        return(-1);
    } else {
        return(1);
    }


}




/*
 *      QUADPACK DQAGSE Routine converted to C
 */
int dqagse(f, qpInfo, a, b, epsabs, epsrel, limit, result, abserr, neval, ier, alist, blist, rlist, elist, iord, last)
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




    /*
     *
     *    Begin prologue:       dqagse
     *    date written   800101   (yymmdd)
     *    revision date  830518   (yymmdd)
     *    category no.  h2a1a1
     *    keywords  automatic integrator, general-purpose,
     *              (end point) singularities, extrapolation,
     *              globally adaptive
     *    author  piessens,robert,appl. math. & progr. div. - k.u.leuven
     *            de doncker,elise,appl. math. & progr. div. - k.u.leuven
     *    purpose  the routine calculates an approximation result to a given
     *             definite integral i = integral of f over (a,b),
     *             hopefully satisfying following claim for accuracy
     *             abs(i-result).le.max(epsabs,epsrel*abs(i)).
     *    description
     *
     *        computation of a definite integral
     *        standard fortran subroutine
     *        double precision version
     *
     *        parameters
     *         on entry
     *            f      - double precision
     *                     function subprogram defining the integrand
     *                     function f(x). the actual name for f needs to be
     *                     declared e x t e r n a l in the driver program.
     *
     *            a      - double precision
     *                     lower limit of integration
     *
     *            b      - double precision
     *                     upper limit of integration
     *
     *            epsabs - double precision
     *                     absolute accuracy requested
     *            epsrel - double precision
     *                     relative accuracy requested
     *                     if  epsabs.le.0
     *                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
     *                     the routine will end with ier = 6.
     *
     *            limit  - integer
     *                     gives an upperbound on the number of subintervals
     *                     in the partition of (a,b)
     *
     *         on return
     *            result - double precision
     *                     approximation to the integral
     *
     *            abserr - double precision
     *                     estimate of the modulus of the absolute error,
     *                     which should equal or exceed abs(i-result)
     *
     *            neval  - integer
     *                     number of integrand evaluations
     *
     *            ier    - integer
     *                     ier = 0 normal and reliable termination of the
     *                             routine. it is assumed that the requested
     *                             accuracy has been achieved.
     *                     ier.gt.0 abnormal termination of the routine
     *                             the estimates for integral and error are
     *                             less reliable. it is assumed that the
     *                             requested accuracy has not been achieved.
     *            error messages
     *                         = 1 maximum number of subdivisions allowed
     *                             has been achieved. one can allow more sub-
     *                             divisions by increasing the value of limit
     *                             (and taking the according dimension
     *                             adjustments into account). however, if
     *                             this yields no improvement it is advised
     *                             to analyze the integrand in order to
     *                             determine the integration difficulties. if
     *                             the position of a local difficulty can be
     *                             determined (e.g. singularity,
     *                             discontinuity within the interval) one
     *                             will probably gain from splitting up the
     *                             interval at this point and calling the
     *                             integrator on the subranges. if possible,
     *                             an appropriate special-purpose integrator
     *                             should be used, which is designed for
     *                             handling the type of difficulty involved.
     *                         = 2 the occurrence of roundoff error is detec-
     *                             ted, which prevents the requested
     *                             tolerance from being achieved.
     *                             the error may be under-estimated.
     *                         = 3 extremely bad integrand behaviour
     *                             occurs at some points of the integration
     *                             interval.
     *                         = 4 the algorithm does not converge.
     *                             roundoff error is detected in the
     *                             extrapolation table.
     *                             it is presumed that the requested
     *                             tolerance cannot be achieved, and that the
     *                             returned result is the best which can be
     *                             obtained.
     *                         = 5 the integral is probably divergent, or
     *                             slowly convergent. it must be noted that
     *                             divergence can occur with any other value
     *                             of ier.
     *                         = 6 the input is invalid, because
     *                             epsabs.le.0 and
     *                             epsrel.lt.max(50*rel.mach.acc.,0.5d-28).
     *                             result, abserr, neval, last, rlist(1),
     *                             iord(1) and elist(1) are set to zero.
     *                             alist(1) and blist(1) are set to a and b
     *                             respectively.
     *
     *            alist  - double precision
     *                     vector of dimension at least limit, the first
     *                      last  elements of which are the left end points
     *                     of the subintervals in the partition of the
     *                     given integration range (a,b)
     *
     *            blist  - double precision
     *                     vector of dimension at least limit, the first
     *                      last  elements of which are the right end points
     *                     of the subintervals in the partition of the given
     *                     integration range (a,b)
     *
     *            rlist  - double precision
     *                     vector of dimension at least limit, the first
     *                      last  elements of which are the integral
     *                     approximations on the subintervals
     *
     *            elist  - double precision
     *                     vector of dimension at least limit, the first
     *                      last  elements of which are the moduli of the
     *                     absolute error estimates on the subintervals
     *
     *            iord   - integer
     *                     vector of dimension at least limit, the first k
     *                     elements of which are pointers to the
     *                     error estimates over the subintervals,
     *                     such that elist(iord(1)), ..., elist(iord(k))
     *                     form a decreasing sequence, with k = last
     *                     if last.le.(limit/2+2), and k = limit+1-last
     *                     otherwise
     *
     *            last   - integer
     *                     number of subintervals actually produced in the
     *                     subdivision process
     *
     *   references  (none)
     *   routines called  d1mach,dqelg,dqk21,dqpsrt
     *   end prologue  dqagse
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




    /*
     *            the dimension of rlist2 is determined by the value of
     *            limexp in subroutine dqelg (rlist2 should be of dimension
     *            (limexp+2) at least).
     *
     *            list of major variables
     *            -----------------------
     *
     *           alist     - list of left end points of all subintervals
     *                       considered up to now
     *           blist     - list of right end points of all subintervals
     *                       considered up to now
     *           rlist(i)  - approximation to the integral over
     *                       (alist(i),blist(i))
     *           rlist2    - array of dimension at least limexp+2 containing
     *                       the part of the epsilon table which is still
     *                       needed for further computations
     *           elist(i)  - error estimate applying to rlist(i)
     *           maxerr    - pointer to the interval with largest error
     *                       estimate
     *           errmax    - elist(maxerr)
     *           erlast    - error on the interval currently subdivided
     *                       (before that subdivision has taken place)
     *           area      - sum of the integrals over the subintervals
     *           errsum    - sum of the errors over the subintervals
     *           errbnd    - requested accuracy max(epsabs,epsrel*
     *                       abs(result))
     *           *****1    - variable for the left interval
     *           *****2    - variable for the right interval
     *           last      - index for subdivision
     *           nres      - number of calls to the extrapolation routine
     *           numrl2    - number of elements currently in rlist2. if an
     *                       appropriate approximation to the compounded
     *                       integral has been obtained it is put in
     *                       rlist2(numrl2) after numrl2 has been increased
     *                       by one.
     *           small     - length of the smallest interval considered up
     *                       to now, multiplied by 1.5
     *           erlarg    - sum of the errors over the intervals larger
     *                       than the smallest interval considered up to now
     *           extrap    - logical variable denoting that the routine is
     *                       attempting to perform extrapolation i.e. before
     *                       subdividing the smallest interval we try to
     *                       decrease the value of erlarg.
     *           noext     - logical variable denoting that extrapolation
     *                       is no longer allowed (true value)
     *
     *            machine dependent constants
     *            ---------------------------
     *
     *           epmach is the largest relative spacing.
     *           uflow is the smallest positive magnitude.
     *           oflow is the largest positive magnitude.
     */




    /*  first executable statement  dqagse */
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
    if( (epsabs <= 0.0) && (epsrel < dmax1(50.0*epmach, 0.5e-28))  ) {
	*ier = 6;
	return(-1);
    }




    /*
     *           first approximation to the integral
     *           -----------------------------------
     */
    uflow = d1mach(1);
    oflow = d1mach(2);
    ierro = 0;
    dqk21(f, qpInfo, a, b, result, abserr, &defabs, &resabs);



    /*
     *           test on accuracy.
     */
    dres = fabs(*result);
    errbnd = dmax1(epsabs, epsrel*dres);
    *last = 1;
    rlist[1] = *result;
    elist[1] = *abserr;
    iord[1] = 1;
    if (  (*abserr <= 100.0*epmach*defabs) && (*abserr > errbnd)  ) *ier = 2;
    if (limit == 1) *ier = 1;
    if (  (*ier != 0) || ( (*abserr <= errbnd) && (*abserr != resabs)) ||  (*abserr == 0.0)  ) {
	*neval = 42 * *last - 21;
	return(0);
    }



    /*
     *           initialization
     *           --------------
     */
    rlist2[1] = *result;
    errmax = *abserr;
    maxerr = 1;
    area = *result;
    errsum = *abserr;
    *abserr = oflow;
    nrmax = 1;
    nres = 0;
    numrl2 = 2;
    ktmin = 0;
    extrap = FALSE;
    noext = FALSE;
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
        dqk21(f, qpInfo, a1, b1, &area1, &error1, &resabs, &defab1);
        dqk21(f, qpInfo, a2, b2, &area2, &error2, &resabs, &defab2);



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
        rlist[*last] = area2;
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
	    small = fabs(b-a)*0.375;
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
    if (*ier > 2) --(*ier);
//L140:
    *neval = 42*(*last) - 21;

    return(0);


}





/*
 *      QUADPACK DQELG Routine converted to C
 */
int dqelg(int n, double epstab[], double *result, double *abserr, double res3la[], int *nres) {


    /*
     *
     *    Begin prologue:       dqpelg
     *    refer to  dqagie,dqagoe,dqagpe,dqagse
     *    routines called  d1mach
     *    revision date  830518   (yymmdd)
     *    keywords  epsilon algorithm, convergence acceleration,
     *              extrapolation
     *    author  piessens,robert,appl. math. & progr. div. - k.u.leuven
     *           de doncker,elise,appl. math & progr. div. - k.u.leuven
     *    purpose  the routine determines the limit of a given sequence of
     *            approximations, by means of the epsilon algorithm of
     *            p.wynn. an estimate of the absolute error is also given.
     *            the condensed epsilon table is computed. only those
     *            elements needed for the computation of the next diagonal
     *            are preserved.
     *    *    description
     *
     *           epsilon algorithm
     *           standard fortran subroutine
     *           double precision version
     *
     *           parameters
     *              n      - integer
     *                       epstab(n) contains the new element in the
     *                       first column of the epsilon table.
     *
     *              epstab - double precision
     *                       vector of dimension 52 containing the elements
     *                       of the two lower diagonals of the triangular
     *                       epsilon table. the elements are numbered
     *                       starting at the right-hand corner of the
     *                       triangle.
     *
     *              result - double precision
     *                       resulting approximation to the integral
     *
     *              abserr - double precision
     *                       estimate of the absolute error computed from
     *                       result and the 3 previous results
     *
     *              res3la - double precision
     *                       vector of dimension 3 containing the last 3
     *                       results
     *
     *              nres   - integer
     *                       number of calls to the routine
     *                       (should be zero at first call)
     *
     *    end prologue  dqelg
     */




    double 	delta1, delta2, delta3, d1mach();
    double	epmach, epsinf, error, err1, err2, err3, e0, e1, e1abs, e2, e3;
    double	oflow, res, ss, tol1, tol2, tol3;
    int 	i, ib, ib2, ie, indx, k1, k2, k3, limexp, newelm, num;



    /*
     *           list of major variables
     *           -----------------------
     *
     *           e0     - the 4 elements on which the computation of a new
     *           e1       element in the epsilon table is based
     *           e2
     *           e3                 e0
     *                        e3    e1    new
     *                              e2
     *           newelm - number of elements to be computed in the new
     *                    diagonal
     *           error  - error = abs(e1-e0)+abs(e2-e1)+abs(new-e2)
     *           result - the element in the new diagonal with least value
     *                    of error
     *
     *           machine dependent constants
     *           ---------------------------
     *
     *           epmach is the largest relative spacing.
     *           oflow is the largest positive magnitude.
     *           limexp is the maximum number of elements the epsilon
     *           table can contain. if this number is reached, the upper
     *           diagonal of the epsilon table is deleted.
     */





    /*
     *    first executable statement  dqelg
     */
    epmach = d1mach(4);
    oflow = d1mach(2);
    ++(*nres);
    *abserr = oflow;
    *result = epstab[n];
    if ( n < 3 ) {
	*abserr = dmax1(*abserr, 5.0*epmach*fabs( *result ));
	return(0);
    }
    limexp = 50;
    epstab[n+2] = epstab[n];
    newelm = (n-1)/2;
    epstab[n] = oflow;
    num = n;
    k1 = n;

    for (i = 1; i <= newelm; ++i) {
        k2 = k1-1;
        k3 = k1-2;
        res = epstab[k1+2];
        e0 = epstab[k3];
        e1 = epstab[k2];
        e2 = res;
        e1abs = fabs(e1);
        delta2 = e2-e1;
        err2 = fabs(delta2);
        tol2 = dmax1(fabs(e2), e1abs)*epmach;
        delta3 = e1-e0;
        err3 = fabs(delta3);
        tol3 = dmax1(e1abs, fabs(e0))*epmach;

        if( (err2 <= tol2) && (err3 <= tol3)) {

	    /*
	     *           if e0, e1 and e2 are equal to within machine
	     *           accuracy, convergence is assumed.
	     *           result = e2
	     *           abserr = abs(e1-e0)+abs(e2-e1)
	     */
	    *result = res;
	    *abserr = err2+err3;
	    /*
	     *   jump out of do-loop
	     */
	    *abserr = dmax1(*abserr, 5.0*epmach*fabs( *result ));
	    return(0);
	}


	e3 = epstab[k1];
        epstab[k1] = e1;
        delta1 = e1-e3;
        err1 = fabs(delta1);
        tol1 = dmax1(e1abs, fabs(e3))*epmach;


	/*
	 *           if two elements are very close to each other, omit
	 *           a part of the table by adjusting the value of n
	 */
        if( (err1 <= tol1) || (err2 <= tol2) || (err3 <= tol3) ) {
	    n = i+i-1;
	    /*
	     *    jump out of do-loop
	     */
	    break;

	} 
        ss = 1.0/delta1+1.0/delta2-1.0/delta3;
        epsinf = fabs(ss*e1);




	/*
	 *           test to detect irregular behaviour in the table, and
	 *           eventually omit a part of the table adjusting the value
	 *           of n.
	 */
        if( epsinf <= 0.1e-3 ) {

	    n = i+i-1;
	    /*
	     *    jump out of do-loop
	     */
	    break;

	} else {

	    /*
             *           compute a new element and eventually adjust
             *           the value of result.
             */
	    res = e1+1.0/ss;
	    epstab[k1] = res;
	    k1 = k1-2;
	    error = err2+fabs(res-e2)+err3;
	    if (error <= *abserr) {
		*abserr = error;
		*result = res;
	    }

	}

    }   




    /*
     *           shift the table.
     */
    if ( n == limexp ) n = 2*(limexp/2)-1;

    ib = 1;
    if ( (num/2)*2 == num ) ib = 2;
    ie = newelm+1;

    for (i=1; i<=ie; ++i) {
        ib2 = ib+2;
        epstab[ib] = epstab[ib2];
        ib = ib2;
    }

    if ( num != n ) {
        indx = num-n+1;
        for (i = 1; i<=n; ++i) {
            epstab[i]= epstab[indx];
            ++indx;
        }
    }


    if ( *nres < 4 ) {

        res3la[*nres] = *result;
        *abserr = oflow;

    } else {

        /*
         *           compute error estimate
         */
	*abserr = fabs( *result-res3la[3] ) + fabs( *result-res3la[2] ) + fabs( *result-res3la[1] );
	res3la[1] = res3la[2];
	res3la[2] = res3la[3];
	res3la[3] = *result;

    }

    *abserr = dmax1(*abserr, 5.0*epmach*fabs( *result ));

    return( 0 );


}





/*
 *      QUADPACK DQK21 Routine converted to C
 */
int dqk21(f, qpInfo, a, b, result, abserr, resabs, resasc)
double  (*f)( double, _qpInfo *); /*  The integrand function -- I.e. the function to integrate           */
_qpInfo *qpInfo;                  /*  Auxilliary information to pass to function (to avoid making globals) */
double   a;             	  /*  Lower Limit of integration.                                 */
double   b;             	  /*  Upper limit of integration.                                 */
double  *result;        	  /*  The desired result. I.e. integral of f() from a to b        */
double  *abserr;        	  /*  Estimate of the modulus of the absolute error in the result */
double  *resabs;        	  /*  */
double  *resasc;        	  /*  */
{




    /*
     *
     *    Begin prologue:       dqk21
     *    date written   800101   (yymmdd)
     *    revision date  830518   (yymmdd)
     *    category no.  h2a1a2
     *    keywords  21-point gauss-kronrod rules
     *    author  piessens,robert,appl. math. & progr. div. - k.u.leuven
     *           de doncker,elise,appl. math. & progr. div. - k.u.leuven
     *    purpose  to compute i = integral of f over (a,b), with error
     *                           estimate
     *                       j = integral of abs(f) over (a,b)
     *    description
     *
     *           integration rules
     *           standard fortran subroutine
     *           double precision version
     *
     *           parameters
     *            on entry
     *              f      - double precision
     *                       function subprogram defining the integrand
     *                       function f(x). the actual name for f needs to be
     *                       declared e x t e r n a l in the driver program.
     *
     *              a      - double precision
     *                       lower limit of integration
     *
     *              b      - double precision
     *                       upper limit of integration
     *
     *            on return
     *              result - double precision
     *                       approximation to the integral i
     *                       result is computed by applying the 21-point
     *                       kronrod rule (resk) obtained by optimal addition
     *                       of abscissae to the 10-point gauss rule (resg).
     *
     *              abserr - double precision
     *                       estimate of the modulus of the absolute error,
     *                       which should not exceed abs(i-result)
     *
     *              resabs - double precision
     *                       approximation to the integral j
     *
     *              resasc - double precision
     *                       approximation to the integral of abs(f-i/(b-a))
     *                       over (a,b)
     *
     *      references  (none)
     *      routines called  d1mach
     *      end prologue  dqk21
     */




    double 	absc, centr, dhlgth, d1mach();
    double	epmach, fc, fsum, fval1, fval2, fv1[11], fv2[11], hlgth;
    double 	resg, resk, reskh, uflow;
    int 	j, jtw, jtwm1;




 
    /*
     *           the abscissae and weights are given for the interval (-1,1).
     *           because of symmetry only the positive abscissae and their
     *           corresponding weights are given.
     *
     *           xgk    - abscissae of the 21-point kronrod rule
     *                    xgk(2), xgk(4), ...  abscissae of the 10-point
     *                    gauss rule
     *                    xgk(1), xgk(3), ...  abscissae which are optimally
     *                    added to the 10-point gauss rule
     *
     *           wgk    - weights of the 21-point kronrod rule
     *
     *           wg     - weights of the 10-point gauss rule
     *
     *
     * gauss quadrature weights and kronron quadrature abscissae and weights
     * as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
     * bell labs, nov. 1981.
     */



    double wg[] = { 	0.0,
		0.066671344308688137593568809893332,
		0.149451349150580593145776339657697,
		0.219086362515982043995534934228163,
		0.269266719309996355091226921569469,
		0.295524224714752870173892994651338 };

    double xgk[] = {	0.0,
		0.995657163025808080735527280689003,
		0.973906528517171720077964012084452,
		0.930157491355708226001207180059508,
		0.865063366688984510732096688423493,
		0.780817726586416897063717578345042,
		0.679409568299024406234327365114874,
		0.562757134668604683339000099272694,
		0.433395394129247190799265943165784,
		0.294392862701460198131126603103866,
		0.148874338981631210884826001129720,
		0.000000000000000000000000000000000 };
 
    double wgk[] = {	0.0,
		0.011694638867371874278064396062192,
		0.032558162307964727478818972459390,
		0.054755896574351996031381300244580,
		0.075039674810919952767043140916190,
		0.093125454583697605535065465083366,
		0.109387158802297641899210590325805,
		0.123491976262065851077958109831074,
		0.134709217311473325928054001771707,
		0.142775938577060080797094273138717,
		0.147739104901338491374841515972068,
		0.149445554002916905664936468389821 };



    /*
     *
     *           list of major variables
     *           -----------------------
     *
     *           centr  - mid point of the interval
     *           hlgth  - half-length of the interval
     *           absc   - abscissa
     *           fval*  - function value
     *           resg   - result of the 10-point gauss formula
     *           resk   - result of the 21-point kronrod formula
     *           reskh  - approximation to the mean value of f over (a,b),
     *                    i.e. to i/(b-a)
     *
     *
     *           machine dependent constants
     *           ---------------------------
     *
     *           epmach is the largest relative spacing.
     *           uflow is the smallest positive magnitude.
     */




    /*
     *   first executable statement  dqk21
     */
    epmach = d1mach(4);
    uflow = d1mach(1);
 
    centr = 0.5*(a+b);
    hlgth = 0.5*(b-a);
    dhlgth = fabs(hlgth);


    /*
     *           compute the 21-point kronrod approximation to
     *           the integral, and estimate the absolute error.
     */
    resg = 0.0;
    fc = (*f)(centr, qpInfo);
    resk = wgk[11]*fc;
    *resabs = fabs(resk);

    for (j=1; j<=5; ++j) {
        jtw = 2*j;
        absc = hlgth*xgk[jtw];
        fval1 = (*f)(centr-absc, qpInfo);
        fval2 = (*f)(centr+absc, qpInfo);
        fv1[jtw] = fval1;
        fv2[jtw] = fval2;
        fsum = fval1+fval2;
        resg += wg[j]*fsum;
        resk += wgk[jtw]*fsum;
        *resabs += wgk[jtw]*(fabs(fval1)+fabs(fval2));
    }



    for (j = 1; j<=5; ++j) {
        jtwm1 = 2*j-1;
        absc = hlgth*xgk[jtwm1];
        fval1 = (*f)(centr-absc, qpInfo);
        fval2 = (*f)(centr+absc, qpInfo);
        fv1[jtwm1] = fval1;
        fv2[jtwm1] = fval2;
        fsum = fval1+fval2;
        resk += wgk[jtwm1]*fsum;
        *resabs += wgk[jtwm1]*(fabs(fval1)+fabs(fval2));
    }



    reskh = resk*0.5;
    *resasc = wgk[11]*fabs(fc-reskh);

    for (j=1; j<=10; ++j) *resasc += wgk[j]*(fabs(fv1[j]-reskh)+fabs(fv2[j]-reskh));

    *result = resk*hlgth;
    *resabs *= dhlgth;
    *resasc *= dhlgth;
    *abserr = fabs((resk-resg)*hlgth);

    if( (*resasc != 0.0) && (*abserr != 0.0)) *abserr = *resasc*dmin1( 1.0, pow(200.0*(*abserr)/(*resasc), 1.5) );

    if ( *resabs > uflow/(50.0*epmach) ) *abserr = dmax1( (epmach*50.0)*(*resabs), *abserr );




    return(1);




}







/*
 *      QUADPACK DQPSRT Routine converted to C
 */
int dqpsrt(int limit, int last, int *maxerr, double *ermax, double elist[], int iord[], int *nrmax) {



    /*
     *
     *    Begin prologue:       dqpsrt
     *    refer to  dqage,dqagie,dqagpe,dqawse
     *    routines called  (none)
     *    revision date  810101   (yymmdd)
     *    keywords  sequential sorting
     *    author  piessens,robert,appl. math. & progr. div. - k.u.leuven
     *           de doncker,elise,appl. math. & progr. div. - k.u.leuven
     *    purpose  this routine maintains the descending ordering in the
     *            list of the local error estimated resulting from the
     *            interval subdivision process. at each call two error
     *            estimates are inserted using the sequential search
     *            method, top-down for the largest error estimate and
     *            bottom-up for the smallest error estimate.
     *    description
     *
     *           ordering routine
     *           standard fortran subroutine
     *           double precision version
     *
     *           parameters (meaning at output)
     *              limit  - integer
     *                       maximum number of error estimates the list
     *                       can contain
     *
     *              last   - integer
     *                       number of error estimates currently in the list
     *
     *              maxerr - integer
     *                       maxerr points to the nrmax-th largest error
     *                       estimate currently in the list
     *
     *              ermax  - double precision
     *                       nrmax-th largest error estimate
     *                       ermax = elist(maxerr)
     *
     *              elist  - double precision
     *                       vector of dimension last containing
     *                       the error estimates
     *
     *              iord   - integer
     *                       vector of dimension last, the first k elements
     *                       of which contain pointers to the error
     *                       estimates, such that
     *                       elist(iord(1)),...,  elist(iord(k))
     *                       form a decreasing sequence, with
     *                       k = last if last.le.(limit/2+2), and
     *                       k = limit+1-last otherwise
     *
     *              nrmax  - integer
     *                       maxerr = iord(nrmax)
     *
     *     end prologue  dqpsrt
     */



    double 	errmax, errmin;
    int 	i, ibeg, ido, isucc, j, jbnd, jupbn, k;




    /*
     *           check whether the list contains more than
     *           two error estimates.
     */



    /*
     *    first executable statement  dqpsrt
     */
    if ( last <= 2 ) {
        iord[1] = 1;
        iord[2] = 2;
    } else {



        /*
         *           this part of the routine is only executed if, due to a
         *           difficult integrand, subdivision increased the error
         *           estimate. in the normal case the insert procedure should
         *           start after the nrmax-th largest error estimate.
         */
	errmax = elist[*maxerr];
	if(*nrmax != 1) {
	    ido = *nrmax-1;
	    for (i = 1; i<=ido; ++i) {
	        isucc = iord[*nrmax-1];
	        /*
	         *    jump out of do-loop
	         */
	        if(errmax <= elist[isucc]) break;
	        iord[*nrmax] = isucc;
	        --(*nrmax);
	    }
	}



        /*
         *           compute the number of elements in the list to be maintained
         *           in descending order. this number depends on the number of
         *           subdivisions still allowed.
         */
	jupbn = last;
	if(last > (limit/2+2)) jupbn = limit+3-last;
	errmin = elist[last];



	/*
	 *           insert errmax by traversing the list top-down,
	 *           starting comparison from the element elist(iord(nrmax+1)).
	 */
	jbnd = jupbn-1;
	ibeg = *nrmax+1;
	if(ibeg > jbnd) goto L50;

	for (i=ibeg; i<=jbnd; ++i) {
            isucc = iord[i];
	    /*
	     *   jump out of do-loop
	     */
            if(errmax >= elist[isucc]) goto L60;
            iord[i-1] = isucc;
	}


L50:
	iord[jbnd] = *maxerr;
	iord[jupbn] = last;
	goto L90;



	/*
	 *           insert errmin by traversing the list bottom-up.
	 */
L60:
	iord[i-1] = *maxerr;
	k = jbnd;

	for (j=i; j<=jbnd; ++j) {
	    isucc = iord[k];
	    /*
	     *   jump out of do-loop
	     */
	    if(errmin < elist[isucc]) goto L80;
	    iord[k+1] = isucc;
	    --k;
	}

	iord[i] = last;
	goto L90;


L80:
	iord[k+1] = last;




    }

    /*
     *           set maxerr and ermax.
     */
L90:
    *maxerr = iord[*nrmax];
    *ermax = elist[*maxerr];
    return(1);


}




    /*
     *   FIX, FIX !!! This routine needs to find these parameters automeatically!!!!
     *   d1mach( 1) = B^EMIN,   	    the smallest positive number.
     *   d1mach( 2) = (1-B^(-T-1)) B^EMAX,   the largest positive number.
     *   d1mach( 4) = B^(-T),   smallest number that when added to 1.0 gives something different from 1.0
     *   d1mach( 3) = B^(-T-1), smallest number that when subtracted from 1.0 gives something different from 1.0
     *   d1mach( 5) = log10(B)
     */

/*
 *  For IEEE-compliant machines we should have (for double precision):
 *
 *  	B    = 2
 *	EMIN = -1022
 *  	EMAX = 1024
 *	T    = 52       (-T is the smallest power of B that added to 1.0 gives something different from 1.0)
 *  
 *   So this gives:
 *	
 *	B^EMIN   = 2.225073698e-308 
 *	B^EMAX   = 1.797693301e308      (B^EMAX - B^(-T-1)*B^EMAX is the same to within precision here )
 *	B^(-T-1) = 1.110223025e-16	(smallest number that when subtracted from 1.0 gives something different from 1.0)
 *	B^(-T)   = 2.220446049e-16	(smallest number that when added to 1.0 gives something different from 1.0)
 *
 */


double d1mach( int i ) {

    switch (i) {

    	case 1:
		    return( 2.2250739e-308 );  	/* assumes IEEE-compliant double precision machine */
		    //return( 1.17549e-38 );
		    break;
    	case 2:
		    return( 1.797693e308 );  	/* assumes IEEE-compliant double precision machine */
		    //return( 3.40282e+38 );
		    break;
    	case 3:
		    return( 2.220446049e-16 );	/* assumes IEEE-compliant double precision machine */
		    //return( 5.96046e-08 );	
		    break;
    	case 4:
		    return( 1.110223025e-16 );	/* assumes IEEE-compliant double precision machine */
		    //return( 1.19209e-07 );
		    break;
    	case 5:
		    return( 0.301029995);	/* assumes IEEE-compliant machine */
		    //return( 0.301029995);
		    break;
    }

    return(-999.9);

}

