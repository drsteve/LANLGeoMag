#include "Lgm/Lgm_MagModelInfo.h"
#define JUMP_METHOD	1







/*
 *   These routines evaluates the "Sb integral" from mirror point to mirror
 *   point. 
 *
 *   There are two different approaches here. Both use Quadpack routines to
 *   evaluate the integral. And Quadpack requires that we be able to evaluate the
 *   integrand at arbitrary "s" vlaues (s is distance along the field line).
 *   The two approaches to doing this are as follows:
 *
 *      (1) SbIntegral() -- this routine gets "s" values (distance along the
 *          FL) from Quadpack. Then we have to trace along the field line until we get to
 *          that s value.  This was the original routine which was written long ago when
 *          B-field models were quite cheap to evaluate. The advantage of this approach is
 *          that it can yield more accurate results since we are tracing along mthe field
 *          line rather than interpolating mbetween pre-traced points. The disadvantage is
 *          that it can be very slow and round \-off error could start to be a problem if
 *          too many evaluationsm are needed.
 *
 *      (2) SbIntegral_intered() -- these routines (the ones with "_interped"
 *          in the name) are currently the preferred method. They relay on use of the
 *          routine BofS() to get the magnitude of B at any value of s. This routine
 *          require a pre-trace of the field line and then an initialization of the spline
 *          function used in BofS().
 *
 *   The Sb integral is as follows:
 *
 *
 *                          / sm_north
 *                         |                 [                ] (-1/2)
 *                         |	             [	      B(s)    ]
 *                 Sb  =   |                 [ 1 -  --------  ]  	   ds
 *                         |                 [         Bm     ]
 *                         |                 [                ]
 *                        / sm_south
 *
 *
 *  I think these are all now reentrant.
 *
 */


/*
 * This version does tracing along the FL to get to each requested s value. Can
 * be slow, no longer the preferred method.
 */
double SbIntegral( Lgm_MagModelInfo *fInfo ) {


    double	a, b;
    double	epsabs, epsrel, result, abserr;
    int		npts, key, neval, ier, leniw, lenw, last, iwork[501];
    double	work[2002], points[20];
    _qpInfo	*qpInfo;


    /*
     *  Type-cast our data structure to a generic type.
     *  The structure holds auzilliary info we need down
     *  in the function calls.
     */
    qpInfo = (_qpInfo *)fInfo;



    /*
     *  Limits of integration.
     */
    a = 0.0;
    b = fInfo->Sm_North - fInfo->Sm_South;



    /*
     *   set tolerances.
     */
    //epsabs = fInfo->epsabs;
    //epsrel = fInfo->epsrel;
    epsabs = fInfo->Lgm_Sb_Integrator_epsabs;
    epsrel = fInfo->Lgm_Sb_Integrator_epsrel;

    /*
     *  Init some vars used in Sb_integrand() (these are not declared static in
     *  Sb_integrand() in order to avoid making it non-reentrant).
     */
    fInfo->Lgm_Sb_integrand_S        = 0.0;
    fInfo->Lgm_Sb_integrand_FirstCall = TRUE;



    leniw = 500;
    lenw  = 4*leniw;
    key   = 6;
    //iwork  = (int *) calloc( limit+1, sizeof(int) );
    //work   = (double *) calloc( lenw+1, sizeof(double) );
    //
    points[1] = a;
    points[2] = b;
    npts = 4;
    dqagp(Sb_integrand, qpInfo, a, b, npts, points, epsabs, epsrel, &result, &abserr, &neval, &ier, leniw, lenw, &last, iwork, work, fInfo->VerbosityLevel );
    //free( iwork );
    //free( work );

    //dqk21(Sb_integrand, qpInfo, a, b, &result, &abserr, &resabs, &resasc);


    return( result );

}


double Sb_integrand( double s, _qpInfo *qpInfo ) {


    //static double 	S=0.0;
    //static Lgm_Vector P, u_scale;
    Lgm_Vector 	        Bvec;
    int			        reset=1, done, Count;
    double		        f, g, Hdone, H, Htry, Hdid, Hnext, sgn=1.0, sdid, B, dS;
    Lgm_MagModelInfo	*fInfo;



    /*
     *  Get pointer to our auxilliary data structure.
     */
    fInfo = (Lgm_MagModelInfo *)qpInfo;



    /* 
     * JUMP_METHOD controls how we get to the s-value.
     */
    if ( JUMP_METHOD == 0 ) {

        /*
         *  Set starting point.
         *  If this is the first call for this integral, set point
         *  to the lower limit.
         */
        if ( fInfo->Lgm_Sb_integrand_FirstCall  == TRUE ) {
            fInfo->Lgm_Sb_integrand_FirstCall = FALSE;
	        fInfo->Lgm_Sb_integrand_u_scale.x = fInfo->Lgm_Sb_integrand_u_scale.y = fInfo->Lgm_Sb_integrand_u_scale.z = 1.0;
            fInfo->Lgm_Sb_integrand_P = fInfo->Pm_South;
	        dS = s;
	        fInfo->Lgm_Sb_integrand_S = 0.0;
        } else {
	        dS = s - fInfo->Lgm_Sb_integrand_S;
        }

        if (dS < 0.0) {
            H = -dS; sgn = -1.0;
        } else {
            H = dS; sgn = 1.0;
        }
        //printf("1. dS, s, H S = %g %g %g %g\n", dS, s, H, fInfo->Lgm_Sb_integrand_S);

    } else {

        /*
         *  this strategy starts at start each time.
         *  Seems to limit roundoff error from accumulating?
         *  But may require more tracing over-all.
         */
        if ( fInfo->Lgm_Sb_integrand_FirstCall  == TRUE ) {
    	    fInfo->Lgm_Sb_integrand_FirstCall = FALSE;
        }
	    fInfo->Lgm_Sb_integrand_u_scale.x = fInfo->Lgm_Sb_integrand_u_scale.y = fInfo->Lgm_Sb_integrand_u_scale.z = 1.0;
    	fInfo->Lgm_Sb_integrand_P = fInfo->Pm_South;
    	H = s; sgn = 1.0;

    }


    /*
     *  Get B-field at the given value of s.
     *  Need to advance along field line by an amount Htry. May need to
     *  do more than one call to get there...
     */
    done = FALSE; Count = 0; Htry = H; Hdone = 0.0; reset = 1;
    if (Htry < 1e-12) done = TRUE;
    while ( !done ) {

        //Lgm_MagStep( &fInfo->Lgm_Sb_integrand_P, &fInfo->Lgm_Sb_integrand_u_scale, Htry, &Hdid, &Hnext, 1.0e-7, sgn, &sdid, &reset, fInfo->Bfield, fInfo );
        Lgm_MagStep( &fInfo->Lgm_Sb_integrand_P, &fInfo->Lgm_Sb_integrand_u_scale, Htry, &Hdid, &Hnext, sgn, &sdid, &reset, fInfo->Bfield, fInfo );
        fInfo->Lgm_Sb_integrand_S += sgn*Hdid;
    	Hdone += Hdid;
    	if ( (Htry < 1e-12) || fabs(Hdone-H) < 1e-12 ){
    	    done = TRUE;
    	} else {
    	    Htry = H - Hdone;
    	}
    	++Count;

    }
    //printf("H = %g , Htry = %g Hdid = %g  Hnext = %g  Hdone = %g Count = %d P = %g %g %g\n", H, Htry, Hdid, Hnext, Hdone, Count, P.x, P.y, P.z);
    //printf("2. dS, s, Hdone S = %g %g %g %g    Count = %d\n", dS, s, Hdone, fInfo->Lgm_Sb_integrand_S, Count);

    fInfo->Bfield( &fInfo->Lgm_Sb_integrand_P, &Bvec, fInfo );
    B = Lgm_Magnitude( &Bvec );

    g = 1.0 - B/fInfo->Bm;
    f = (g > 0.0) ? sqrt( g ) : 0.0;

    ++fInfo->Lgm_n_Sb_integrand_Calls;

//printf("Sb_integrand: s = %g  f = %g   g = %g  B = %g    Info->Pm_South = %g %g %g   P = %g %g %g\n", s, f, g, B, fInfo->Pm_South.x, fInfo->Pm_South.y, fInfo->Pm_South.z, fInfo->Lgm_Sb_integrand_P.x, fInfo->Lgm_Sb_integrand_P.y, fInfo->Lgm_Sb_integrand_P.z);
    if ( f == 0.0 ) {
        return( 0.0 );
    } else {
        return( 1.0/f );
    }


}




/*
 * This uses spline interpolation to get the B(s) values needed.  Note that the
 * fInfo structure must have be pre-initialized with a number of this for this
 * to work. You can do it by hand, or look at SbIntegral_interped_Setup() and
 * SbIntegral_interped_Teardown() to see what is involved.
 */
double SbIntegral_interped( Lgm_MagModelInfo *fInfo ) {


    double	a, b;
    double	epsabs, epsrel, result, abserr;
    int		npts, key, neval, ier, leniw, lenw, last, iwork[501];
    double	work[2002], points[20];
    _qpInfo	*qpInfo;


    /*
     *  Type-cast our data structure to a generic type.
     *  The structure holds auzilliary info we need down
     *  in the function calls.
     */
    qpInfo = (_qpInfo *)fInfo;



    /*
     *  Limits of integration.
     */
    a = fInfo->Sm_South;
    b = fInfo->Sm_North;


    //printf("\n\n\nSbIntegral_interped: a, b = %g %g\n", a, b);



    /*
     *   set tolerances.
     */
    //epsabs = fInfo->epsabs;
    //epsrel = fInfo->epsrel;
    epsabs = fInfo->Lgm_Sb_Integrator_epsabs;
    epsrel = fInfo->Lgm_Sb_Integrator_epsrel;


    leniw = 500;
    lenw  = 4*leniw;
    key   = 6;

    /*
     *  Perform integrations. The points[] array contains a list of points
     *  in the integrand where integrable singularities occur. dqagp() uses
     *  these to break up the integral at these points (dqagp() is an
     *  extrapolation algorithm so it handles singularities very well.)
     *  Note: the value of npts2 must be 2+ the number of points added.
     *  This is because it adds the endpoints internally (so dont add those
     *  here).
     */
/*
*/
//    points[1] = a;
//    points[2] = b;
    npts = 2;
    dqagp(Sb_integrand_interped, qpInfo, a, b, npts, points, epsabs, epsrel, &result, &abserr, &neval, &ier, leniw, lenw, &last, iwork, work, fInfo->VerbosityLevel );
/*
printf("here i am\n");
double s;
for(s=a; s<=b; s+=(b-a)/1000.0){
printf("%g %g\n", s, Sb_integrand_interped(s, qpInfo));
}
exit(0);
    dqags(Sb_integrand_interped, qpInfo, a, b, epsabs, epsrel, &result, &abserr, &neval, &ier, leniw, lenw, &last, iwork, work);
*/

    return( result );

}






/*
 *  This computes the Sb integral (using a spine interpolation of the B(s)
 *  function) but with user-supplied limits. 
 */
double SbIntegral_interped2( Lgm_MagModelInfo *fInfo, double a, double b ) {


//    double	a, b;
    double	epsabs, epsrel, result, abserr;
    int		npts, key, neval, ier, leniw, lenw, last, iwork[501];
    double	work[2002], points[20];
    _qpInfo	*qpInfo;


    /*
     *  Type-cast our data structure to a generic type.
     *  The structure holds auzilliary info we need down
     *  in the function calls.
     */
    qpInfo = (_qpInfo *)fInfo;



    /*
     *  Limits of integration.
     */
// Take user inputs instead.
//    a = fInfo->Sm_South;
//    b = fInfo->Sm_North;
    //printf("\n\n\nSbIntegral_interped: a, b = %g %g\n", a, b);



    /*
     *   set tolerances.
     */
    //epsabs = fInfo->epsabs;
    //epsrel = fInfo->epsrel;
    epsabs = fInfo->Lgm_Sb_Integrator_epsabs;
    epsrel = fInfo->Lgm_Sb_Integrator_epsrel;


    leniw = 500;
    lenw  = 4*leniw;
    key   = 6;

    /*
     *  Perform integrations. The points[] array contains a list of points
     *  in the integrand where integrable singularities occur. dqagp() uses
     *  these to break up the integral at these points (dqagp() is an
     *  extrapolation algorithm so it handles singularities very well.)
     *  Note: the value of npts2 must be 2+ the number of points added.
     *  This is because it adds the endpoints internally (so dont add those
     *  here).
     */
/*
*/
//    points[1] = a;
//    points[2] = b;
    npts = 2;
    dqagp(Sb_integrand_interped, qpInfo, a, b, npts, points, epsabs, epsrel, &result, &abserr, &neval, &ier, leniw, lenw, &last, iwork, work, fInfo->VerbosityLevel );
/*
printf("here i am\n");
double s;
for(s=a; s<=b; s+=(b-a)/1000.0){
printf("%g %g\n", s, Sb_integrand_interped(s, qpInfo));
}
exit(0);
    dqags(Sb_integrand_interped, qpInfo, a, b, epsabs, epsrel, &result, &abserr, &neval, &ier, leniw, lenw, &last, iwork, work);
*/

    return( result );

}





/*
 * This Sb integrand uses the spline interpolation in BofS(). The fInfo
 * structure that is passed through Sb_integrand() and Sb_integrand_interped()
 * must be propoerly initialized.
 */
double Sb_integrand_interped( double s, _qpInfo *qpInfo ) {

    double              B, g, f;
    Lgm_MagModelInfo	*fInfo;
    /*
     *  Get pointer to our auxilliary data structure.
     */
    fInfo = (Lgm_MagModelInfo *)qpInfo;

    B = BofS( s, fInfo );
    g = 1.0 - B/fInfo->Bm;
    f = (g > 0.0) ? sqrt( g ) : 0.0;
    ++fInfo->Lgm_n_Sb_integrand_Calls;

    return( (f==0.0) ? 0.0 : 1.0/f );

}






/*
 * Example routine to setup the info needed in the MagModelInfo structure, before
 * calling SbIntegral_interped(). 
 *
 * In this version, we:
 *
 *  (1) take an initial point in GSM coords and trace it with Lgm_Trace(). This
 *      will give us a variety of critical points like footpoints, Bmin point, FL
 *      length, etc, etc...
 *
 *  (2) From the equatorial pitch agle, we compute a Bmirror value. Bmirror =
 *      Beq/sin^2(alpga_eq).
 *
 *  (3) The we retrace the field line to initialize the splinbe used in BofS().
 *
 *  (4) Finally, we compute what the limits of the integral should be. This
 *  ammounts to figuring out what s-vlues the mirror points are at. Most of the
 *  code here is involved in doing this -- is very straight forward bracketing
 *  and bisection search. There is probably a routine in Lgm or gsl somewhere
 *  that'll do the search for you in 1 or two lines, of code. But this givebsm
 *  you an idea of what is involved.
 *
 *  For optimized applications, you may not want to use something like this.
 *  For example, lets say you are working on 1 FL, but have many pitch angles
 *  to do.  Then you could just do the first part (seting up the spline), and
 *  leave it alone. Do the integrals for the different pitch angles just
 *  requires the bottom part (finding and setting the limits and Bmirror
 *  values.) 
 *
 *  Also, note that this will not give the correct value for a_eq = 90deg. The
 *  reason is that the limits will be the same.  Quadpack will tell you the
 *  integral is zero. However, the correct answer is that its not zero. It
 *  asymtotes to a constant at 90deg. See Roederer, page 37, eqn 2.13a and b.
 *  Elsewhere in Lgm (e.g. in Lstar and Drift shell calcs), this condition is
 *  detected and 2.13a is used. 
 *
 *  The convenience rotuine Lgm_Trace() has an option for lcomputing Sb0. To
 *  turn this on, set mInfo->ComputeSb0 = TRUE in the calling routines (result
 *  will be in mInfo->Sb0).
 */
int SbIntegral_interped_Setup( Lgm_Vector *u, double a_eq, Lgm_MagModelInfo *mInfo ) {

    double      atol_save, rtol_save, g;
    double      s, slo, shi, b, blo, bhi, s_bmin;
    double      Bmin, Stotal, sa, sa2, B_mirror, snorth, ssouth;
    int         Flag, i, ilo, ihi, i_bmin, done;
    Lgm_Vector  v1, v2, v3;

    /*
     * Lets increase tolerance on the tracing stuff to ensure Lgm_Trace() and TraceLine3() give same vals.
     */
    //atol_save = mInfo->Lgm_MagStep_BS_atol;
    //rtol_save = mInfo->Lgm_MagStep_BS_rtol;
    //mInfo->Lgm_MagStep_BS_atol = 1e-10;
    //mInfo->Lgm_MagStep_BS_rtol = 0.0;


    /*
     * Use the Lgm_Trace() convenience rotuinem first.
     */
    Flag = Lgm_Trace( u, &v1, &v2, &v3, 120.0, 1e-7, 1e-7, mInfo );

    /*
     * FL should be closed (technically an open FL could work if there is a
     * Bmin between two vals of Bm, but we'll ignore this case.)
     */
    if ( Flag != LGM_CLOSED ) {
        //mInfo->Lgm_MagStep_BS_atol = atol_save;
        //mInfo->Lgm_MagStep_BS_rtol = rtol_save;
        return( -1 );
    }

    // this is the Bmin value.
    Bmin = mInfo->Bmin;

    // this is the total FL length bwetween footpoints.
    Stotal = mInfo->Stotal;

    // Compute B_mirror
    sa = sin( a_eq *RadPerDeg ); sa2 = sa*sa;
    B_mirror = Bmin/sa2;
    
    /*
     * Now, re-trace the FL with regular steps. Start at the southern
     * footpoint.  Use 1000 steps and trace a totsal distance of Stotal. This
     * will fill arrays in the mInfo structure with the FL points (and B vlas,
     * etc.) The second to last argument is 0 or 1 (1 means add the Bmin point
     * the arrays, 0 means dont do that.) Here we will add it, but be careful!
     * If your atol vals arent very strict, you can get different mresults
     * between Lgm_Trace() and Lgm_TraceLine3(). Then adding the Bmin (which
     * was got from Lgm_Trace(), may cause a glitchy discontinuity.)
     */

    Lgm_TraceLine3( &v1, Stotal, 1000 , 1.0, 1e-7, 1, mInfo );
    if ( !InitSpline( mInfo ) ) {
        printf( "SbIntegral_interped_Setup: Failed to initialize spline.\n");
        //mInfo->Lgm_MagStep_BS_atol = atol_save;
        //mInfo->Lgm_MagStep_BS_rtol = rtol_save;
        return(-1);
    }
    

    /*
     * At this point, we now have an initialized spline curve in the mInfo
     * structure. It can be called using the BofS() routine. I.e. B( s, mInfo )
     * will returnm the magnitude of B at the given s value. Note that s will
     * range from 0 to Stotal. To verify this you could look at mInfo->s[0] and
     * mInfo->s[mInfo->nPnts-1], which are the first and last s-values in the
     * arrays.
     * 
     *  We are not done yet though. For the Sb integral, we need to integrate
     *  from mirror point to mirror point.  This ammounts to figuring out wat
     *  values of s correspond to the Bmirror values.  There are a number of
     *  ways to get these. We could use the Lgm_TraceToMirrorPoint() routine or
     *  we could just find it with bisection on thej BofS() function we just
     *  initialized. Lets try the later approach here....
     */

    // locate Bmin point in the arrays.
    g = 1e99;
    i_bmin = 0;

    for (i=0; i<mInfo->nPnts; i++ ) {
        if ( mInfo->Bmag[i] < g ) {
            g = mInfo->Bmag[i];
            i_bmin = i;
        }
    }
    s_bmin = mInfo->s[i_bmin];



    // find bracket to the north
    i = i_bmin;
    done = FALSE;
    while ( !done ) {
        if ( i == mInfo->nPnts-1 ) {
            // we are already at the last point. Doesnt look like there is any B small enough!
            // should bail...
            //mInfo->Lgm_MagStep_BS_atol = atol_save;
            //mInfo->Lgm_MagStep_BS_rtol = rtol_save;
            if ( mInfo->AllocedSplines ) FreeSpline( mInfo );
            return(-1);
        }
        ++i;
        if ( mInfo->Bmag[i] >= B_mirror ) {
            // we have a bracket, such that B[ilo] <= Bmirror <= B[ihi]
            ilo = i-1;
            ihi = i;
            done = TRUE;
        }
    }

    // Use bisection to get exact s value in the north.
    slo = mInfo->s[ilo]; blo = BofS( slo, mInfo ); // Use function values to avoid inconsistencies
    shi = mInfo->s[ihi]; bhi = BofS( shi, mInfo ); 
    done = FALSE;
    while ( !done ) {
        // test for convergence
        if ( fabs(slo - shi) < 1e-8 ) {
            done = TRUE; // we have converged.
        } else {
            s = 0.5*( slo + shi ); // bisect current bracket.
            b = BofS( s, mInfo );  // evaluate at this s.
            if ( b < B_mirror ) {  // test result
                slo = s; blo = b;  // reset the lower bound on the bracket.
            } else {
                shi = s; bhi = b;  // reset the upper bound on the bracket.
            }
        }
    }
    snorth = slo; // rather than averaging the samll bracket, take lower value to ensure we are still in range.
    
    



    // find bracket to the south
    i = i_bmin;
    done = FALSE;
    while ( !done ) {
        if ( i == 0 ) {
            // we are already at the first point. Doesnt look like there is any B small enough!
            // should bail...
            //mInfo->Lgm_MagStep_BS_atol = atol_save;
            //mInfo->Lgm_MagStep_BS_rtol = rtol_save;
            if ( mInfo->AllocedSplines ) FreeSpline( mInfo );
            return(-1);
        }
        --i;
        if ( mInfo->Bmag[i] >= B_mirror ) {
            // we have a bracket, such that B[ilo] <= Bmirror <= B[ihi]
            ilo = i+1;
            ihi = i;
            done = TRUE;
        }
    }

    // Use bisection to get exact s value in the south.
    slo = mInfo->s[ilo]; blo = BofS( slo, mInfo ); // Use function values to avoid inconsistencies
    shi = mInfo->s[ihi]; bhi = BofS( shi, mInfo ); 
    done = FALSE;
    while ( !done ) {
        // test for convergence
        if ( fabs(slo - shi) < 1e-8 ) {
            done = TRUE; // we have converged.
        } else {
            s = 0.5*( slo + shi ); // bisect current bracket.
            b = BofS( s, mInfo );  // evaluate at this s.
            if ( b < B_mirror ) {  // test result
                slo = s; blo = b;  // reset the lower bound on the bracket.
            } else {
                shi = s; bhi = b;  // reset the upper bound on the bracket.
            }
        }
    }
    ssouth = shi; // rather than averaging the samll bracket, take upper value to ensure we are still in range.


    /*
     * Set the parameters for the integral ([a,b] limits and Bmirror)
     */
    mInfo->Sm_North = snorth;
    mInfo->Sm_South = ssouth;
    mInfo->Bm       = B_mirror;
    //printf("mInfo->Sm_North, mInfo->Sm_South, Bmin, mInfo->Bm = %g %g %g %g\n", mInfo->Sm_North, mInfo->Sm_South, Bmin, mInfo->Bm);
    
    
    //mInfo->Lgm_MagStep_BS_atol = atol_save;
    //mInfo->Lgm_MagStep_BS_rtol = rtol_save;


    return(1);


}

void SbIntegral_interped_Teardown( Lgm_MagModelInfo *mInfo ) {

    if ( mInfo->AllocedSplines ) FreeSpline( mInfo );
    return;

}
