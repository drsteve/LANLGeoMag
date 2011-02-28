#include "Lgm/Lgm_MagModelInfo.h"
//#include "MagStep.h"

#define JUMP_METHOD	    0

#define USE_SIX_POINT   0
#define USE_FOUR_POINT  1
#define USE_TWO_POINT   2

//#define DIFF_SCHEME     USE_TWO_POINT
#define DIFF_SCHEME     USE_FOUR_POINT





/*
 *   This routine evaluates the "integral invariant, I" from mirror point to
 *   mirror point. Instead of tracing the whole field line, this version allows
 *   the Quadpack integration routine to evaluate B(s). The hope is that this
 *   will decrease the amount of tracing needed to evaluate the integral. It may
 *   very well be slower than tracing the whole line and then doing the integral,
 *   but we'll try this for comparison ... The other obvious way to do it is to
 *   pre-trace the FL and then interpolate the points in some manner. I'm hoping
 *   that combining dqags and Bulirsch-Stoer (which are both extrapolation methods)
 *   that we'll get a speedup.
 *   
 *   The integral is as follows:
 *
 *
 *                          / sm_north       
 *                         |                 [                ] (1/2)
 *                         |	             [	      B(s)    ]
 *                 I  =    |                 [ 1 -  --------  ]  	ds
 *                         |                 [         Bm     ]
 *                         |                 [                ]
 *                        / sm_south
 *                       
 *
 *  I think these are now all reentrant.
 *
 */
double Iinv( Lgm_MagModelInfo *mInfo ) {


    double	a, b;
    double	epsabs, epsrel, result, abserr;
    double  resabs, resasc;
    int		key, neval, ier, limit, lenw, last, iwork[501];
    double	work[2001];
    _qpInfo	*qpInfo;


    /*
     *  Type-cast our data structure to a generic type.
     *  The structure holds auzilliary info we need down 
     *  in the function calls.
     */
    qpInfo = (_qpInfo *)mInfo;



    /*
     *  Limits of integration.
     */
    a = mInfo->Sm_South;
    b = mInfo->Sm_North;


    /*
     *   set tolerances. 
     */
    //epsabs = mInfo->epsabs;
    //epsrel = mInfo->epsrel;
    epsabs = mInfo->Lgm_I_Integrator_epsabs;
    epsrel = mInfo->Lgm_I_Integrator_epsrel;


    /*
     *  Init some vars used in I_integrand() (these are not declared static in
     *  I_integrand() in order to avoid making it non-reentrant).
     */
    mInfo->Lgm_I_integrand_S         = 0.0;
    mInfo->Lgm_I_integrand_FirstCall = TRUE;
    mInfo->Lgm_n_I_integrand_Calls    = 0;

                                                                                                                                                                             


    if ( mInfo->Lgm_I_Integrator == DQAGS ) {

        /*
         *  Use DQAGS
         */
        limit = 500; lenw = 4*limit; key = 6;
        dqags(I_integrand, qpInfo, a, b, epsabs, epsrel, &result, &abserr, &neval, &ier, limit, lenw, &last, iwork, work);

    } else if ( mInfo->Lgm_I_Integrator == DQK21 ) {

        /*
         *  Use DQK21
         */
        dqk21(I_integrand, qpInfo, a, b, &result, &abserr, &resabs, &resasc);

    } else {

        /*
         *  Unknown Integrator
         */
        printf("Iinv: Unknown integrator. Lgm_Inv_Integrator = %d\n", mInfo->Lgm_n_I_integrand_Calls );
        result = -9e99;

    }


    return( result );

}

double Iinv_interped( Lgm_MagModelInfo *mInfo ) {

    double	a, b;
    double	epsabs, epsrel, result, abserr, resabs, resasc;
    int		key, limit, lenw;
    int     iwork[501], last, ier, neval;
    double	work[2001];
    _qpInfo	*qpInfo;


    /*
     *  Type-cast our data structure to a generic type.
     *  The structure holds auzilliary info we need down 
     *  in the function calls.
     */
    qpInfo = (_qpInfo *)mInfo;

    /*
     *  Limits of integration.
     */
    a = mInfo->Sm_South;
    b = mInfo->Sm_North;


    /*
     *   set tolerances for QuadPack routines. 
     */
    //epsabs = mInfo->epsabs;
    //epsrel = mInfo->epsrel;
    epsabs = mInfo->Lgm_I_Integrator_epsabs;
    epsrel = mInfo->Lgm_I_Integrator_epsrel;


    /*
     *  Init some vars used in I_integrand_interped() (these are not declared static in
     *  I_integrand() in order to avoid making it non-reentrant).
     */
    mInfo->Lgm_n_I_integrand_Calls  = 0;



    if ( mInfo->Lgm_I_Integrator == DQAGS ) {

        /*
         *  Use DQAGS
         */
        limit = 500; lenw = 4*limit; key = 6;
        dqags(I_integrand_interped, qpInfo, a, b, epsabs, epsrel, &result, &abserr, &neval, &ier, limit, lenw, &last, iwork, work);

    } else if ( mInfo->Lgm_I_Integrator == DQK21 ) {

        /*
         *  Use DQAGS
         */
        dqk21(I_integrand_interped, qpInfo, a, b, &result, &abserr, &resabs, &resasc);

    } else {

        /*
         *  Unknown Integrator
         */
        printf("Iinv_interped: Unknown integrator. Lgm_I_Integrator = %d\n", mInfo->Lgm_n_I_integrand_Calls );
        result = -9e99;

    }

    return( result );

}


double I_integrand_interped( double s, _qpInfo *qpInfo ) {

    double              B, g, f;
    Lgm_MagModelInfo    *mInfo;

    /*
     *  Get pointer to our auxilliary data structure.
     */
    mInfo = (Lgm_MagModelInfo *)qpInfo;

    B = BofS( s, mInfo );
    g = 1.0 - B/mInfo->Bm;
    f = (g > 0.0) ? sqrt( g ) : 0.0;
    ++mInfo->Lgm_n_I_integrand_Calls;

    return( f );

}




/*
 *  Note: The variables that save our state in this (and other fundtions) are
 *  not local static variables, but are contained in the Lgm_MagModelInfo
 *  structure that the user passes to it. This allows the function to be
 *  thread-safe and re-entrant. 
 */
double I_integrand( double s, _qpInfo *qpInfo ) {


    Lgm_Vector 		    Bvec;
    int			        reset=1, done, Count;
    double		        f, g, Hremaining, Hdone, H, Htry, Hdid, Hnext, sgn=1.0, sdid, B, dS;
    Lgm_MagModelInfo    *mInfo;



    /*
     *  Get pointer to our auxilliary data structure.
     */
    mInfo = (Lgm_MagModelInfo *)qpInfo;

	mInfo->Lgm_I_integrand_u_scale.x =  10.0;  mInfo->Lgm_I_integrand_u_scale.y = 1.0; mInfo->Lgm_I_integrand_u_scale.z = 10.0;



    if ( mInfo->Lgm_I_integrand_JumpMethod == LGM_RELATIVE_JUMP_METHOD ) {

        /*
         *  Set starting point.  If this is the first call for this integral,
         *  set point to the lower limit. 
         */
        if ( mInfo->Lgm_I_integrand_FirstCall  == TRUE ) {
	        mInfo->Lgm_I_integrand_FirstCall = FALSE;
            mInfo->Lgm_I_integrand_P = mInfo->Pm_South;
	        mInfo->Lgm_I_integrand_S = 0.0;
	        dS = s;
        } else {
	        dS = s - mInfo->Lgm_I_integrand_S;
        }

        if (dS < 0.0) {
            H = -dS; sgn = -1.0; // H is a positive qnty
        } else {
            H = dS; sgn = 1.0;
        }

    } else if ( mInfo->Lgm_I_integrand_JumpMethod == LGM_ABSOLUTE_JUMP_METHOD ) {

	    /*
         *  This strategy starts at start each time.  Slower(?), but seems to
         *  limit roundoff error from accumulating? The speed increase of the
         *  relative method really depends on how far apart the s points are
         *  that the integrator is picking. If they are bouncing all over the
         *  place the relative method still does alot of tracing.
	     */
        if ( mInfo->Lgm_I_integrand_FirstCall == TRUE ) {
	        mInfo->Lgm_I_integrand_FirstCall = FALSE;
        } 
	    mInfo->Lgm_I_integrand_P = mInfo->Pm_South;
        mInfo->Lgm_I_integrand_S = 0.0;
	    H = s; sgn = 1.0; // H is a positive qnty
 
    } else {

        printf("I_integrand: Error, unknown Jump Method ( Lgm_I_integrand_JumpMethod = %d\n", mInfo->Lgm_I_integrand_JumpMethod );
        exit(-1);

    }


    /*
     *  Get B-field at the given value of s.  Need to advance along field line
     *  by an amount H. May need to do more than one call to get there...
     *  Setting Htry to be H is an attempt to make the full jump in one call to
     *  MagStep. For very precise calculations, the number of jumps we do here
     *  seems to be fairly critical. Making two jumps (i.e. Htry = H/2) and
     *  limiting Hmax to 0.5 Re seems to work pretty well (perhaps its because
     *  symetry between mirror points?). Could still play with the max step
     *  size of 0.5 -- maybe something a bit higher would work just as well and
     *  be more efficient?
     */
    done = FALSE; Count = 0; Htry = 0.5*H; Hdone = 0.0; reset = TRUE;
    if (Htry > 0.5) Htry = 0.5;
    if ( fabs(mInfo->Lgm_I_integrand_S-s) < 1e-12 ) done = TRUE;

    while ( !done ) {

        Lgm_MagStep( &mInfo->Lgm_I_integrand_P, &mInfo->Lgm_I_integrand_u_scale, Htry, &Hdid, &Hnext, 1.0e-7, sgn, &sdid, &reset, mInfo->Bfield, mInfo );
        Hdone += Hdid; // positive qnty

        mInfo->Lgm_I_integrand_S += sgn*Hdid;
        Hremaining = H - Hdone;
	    if ( Count > 1000 ) {
            printf("I_integrand: Warning, early return because Count > 1000. Ill-conditioned integrand?\n");
	        done = TRUE;
        } else if ( fabs(mInfo->Lgm_I_integrand_S-s) < 1e-12 ){
	        done = TRUE;
	    } else {
            if ( Htry > Hremaining ) Htry = Hremaining;
	    }
	    ++Count;
    }

    mInfo->Bfield( &mInfo->Lgm_I_integrand_P, &Bvec, mInfo );
    B = Lgm_Magnitude( &Bvec );

    /*
     * Compute integrand, ( 1 - B/Bm ) ^ (1/2) Make sure 1-B/Bm is not
     * negative. If it is, assume its zero. (Note that sometimes (due to
     * round-off error) it could be very slightly negative).
     */
    g = 1.0 - B/mInfo->Bm;
    f = (g > 0.0) ? sqrt( g ) : 0.0;

    ++mInfo->Lgm_n_I_integrand_Calls;

    return( f );

}



/*
 *  Compute Gradient of I
 */
int Lgm_Grad_I( Lgm_Vector *v0, Lgm_Vector *GradI, Lgm_MagModelInfo *mInfo ) {

    Lgm_Vector  u, Pa, Pb;
    double  H, h, a, b, Sa, Sb, I, f[6], r;
    int     i, N;
    


    /*
     * We want to compute the gradient of I at the point v0.
     * This requires 3 derivatives: one each in x, y and z directions.
     * We will try a fairly acurate difference scehme:
     *
     *      f_0^(1) = 1/(60h)  ( f_3 - 9f_2 + 45f_1  - 45f_-1 + 9f_-2 - f_-3 )
     * See page 450 of CRC standard Math tables 28th edition.
     */


    /*
     *  Set h to a smallish value
     */
    h = 5e-2;
    h = 0.1;
    h = 0.2;


    switch ( DIFF_SCHEME ) {
        case USE_SIX_POINT:
            N = 3;
            break;
        case USE_FOUR_POINT:
            N = 2;
            break;
        case USE_TWO_POINT:
            N = 1;
            break;
    }


    // User should set these?
    mInfo->Lgm_I_Integrator_epsabs = 0.0;
    mInfo->Lgm_I_Integrator_epsrel = 1e-7;
//            mInfo2->FirstCall = TRUE;

    
    /* X-component */
    mInfo->UseInterpRoutines = 1;
    if (mInfo->VerbosityLevel > 0) printf("\t\tComputing dIdx: h = %g\n", h);
printf("\t\tComputing dIdx: h = %g\n", h);
    for (i=-N; i<=N; ++i){

        if (i!=0) { // dont need the center value in our difference scheme
    
            u = *v0; H = (double)i*h; u.x += H;

            /*
             * Trace to northern mirror point
             */
            //printf("A. Lgm_TraceToMirrorPoint\n");
            if ( Lgm_TraceToMirrorPoint( &u, &Pb, &Sb, mInfo->Bm, 1.0, mInfo->Lgm_TraceToMirrorPoint_Tol, mInfo ) > 0 ) {

                /*
                 * Trace to southern mirror point
                 */
                //printf("B. Lgm_TraceToMirrorPoint\n");
                if ( Lgm_TraceToMirrorPoint( &u, &Pa, &Sa, mInfo->Bm, -1.0, mInfo->Lgm_TraceToMirrorPoint_Tol, mInfo ) > 0 ) {

                    r  = Lgm_Magnitude( &Pb );
//                    Lgm_TraceLine2( &Pa, &Pb, (r-1.0)*Re, 0.5*Sa-mInfo->Hmax, 1.0, 1e-7, FALSE, mInfo );
Lgm_TraceLine2( &Pa, &Pb, (r-1.0)*Re, Sa/200.0, 1.0, 1e-7, FALSE, mInfo );
                    mInfo->Sm_South = 0.0;
                    mInfo->Sm_North = Sa;
                    ReplaceFirstPoint( 0.0, mInfo->Bm, &Pa, mInfo );
                    AddNewPoint( Sa,  mInfo->Bm, &Pb, mInfo );
                    mInfo->FirstCall = TRUE;
                    mInfo->Lgm_n_I_integrand_Calls = 0;
                    InitSpline( mInfo );
mInfo->Lgm_I_integrand_S         = 0.0;                                                                                            
mInfo->Lgm_I_integrand_FirstCall = TRUE;                                                                                           
mInfo->Lgm_n_I_integrand_Calls    = 0;
                    I = Iinv_interped( mInfo  );
                    if (mInfo->VerbosityLevel > 2) printf("I = %g Lgm_n_I_integrand_Calls = %d\n", I, mInfo->Lgm_n_I_integrand_Calls );
                    printf("I = %g Lgm_n_I_integrand_Calls = %d\n", I, mInfo->Lgm_n_I_integrand_Calls );
//                    if (VerbosityLevel > 1) printf("\t\t%s  Integral Invariant, I (interped):      %15.8g    I-I0:    %15.8g    mlat:   %12.8lf  (nCalls = %d)%s\n",  PreStr, I, I-I0, b, mInfo->Lgm_n_I_integrand_Calls, PostStr );
                    FreeSpline( mInfo );

//                    a = 0.0; b = Sa+Sb;
//                    mInfo->Pm_South = Pa; mInfo->Sm_South = a;
//                    mInfo->Pm_North = Pb; mInfo->Sm_North = b;
//                    if (mInfo->VerbosityLevel > 2) printf("\t\tComputing Integral Invariant. Limits of integration are: [a,b] = [%g,%g]\n", a, b);
//
//                    mInfo->FirstCall = TRUE;
//                    mInfo->Lgm_n_I_integrand_Calls = 0;
//                    I = Iinv( mInfo );
                    if (mInfo->VerbosityLevel > 2) printf("\t\tI = %g   ", I);
//                    if (mInfo->VerbosityLevel > 2) printf("Lgm_n_I_integrand_Calls = %d\n", mInfo->Lgm_n_I_integrand_Calls );


                } else {
                    printf("\t\tMirror point below %g km in Southern Hemisphere\n", mInfo->Lgm_LossConeHeight);
                    exit(0);
                }

            } else {
                printf("\t\tMirror point below %g km in Northern Hemisphere\n", mInfo->Lgm_LossConeHeight);
                exit(0);
            }


            f[i+N] = I;

        }

    }
    if (DIFF_SCHEME == USE_SIX_POINT){
        GradI->x = (f[6] - 9.0*f[5] + 45.0*f[4] - 45.0*f[2] + 9.0*f[1] - f[0])/(60.0*h);
    } else if (DIFF_SCHEME == USE_FOUR_POINT){
        GradI->x = (-f[4] + 8.0*f[3] - 8.0*f[1] + f[0])/(12.0*h);
    } else if (DIFF_SCHEME == USE_TWO_POINT){
        GradI->x = (f[2] - f[0])/(2.0*h);
    }
    
    

    /* Y-component */
    if (mInfo->VerbosityLevel > 0) printf("\t\tComputing dIdy: h = %g\n", h);
printf("\t\tComputing dIdy: h = %g\n", h);
    for (i=-N; i<=N; ++i){

        if (i!=0) { // dont need the center value in our difference scheme
    
            u = *v0; H = (double)i*h; u.y += H;

            /*
             * Trace to northern mirror point
             */
            //printf("C. Lgm_TraceToMirrorPoint\n");
            if ( Lgm_TraceToMirrorPoint( &u, &Pb, &Sb, mInfo->Bm, 1.0, mInfo->Lgm_TraceToMirrorPoint_Tol, mInfo ) > 0 ) {

                /*
                 * Trace to southern mirror point
                 */
                //printf("D. Lgm_TraceToMirrorPoint\n");
                if ( Lgm_TraceToMirrorPoint( &u, &Pa, &Sa, mInfo->Bm, -1.0, mInfo->Lgm_TraceToMirrorPoint_Tol, mInfo ) > 0 ) {

                    r  = Lgm_Magnitude( &mInfo->Pm_North );
//                    Lgm_TraceLine2( &Pa, &Pb, (r-1.0)*Re, 0.5*Sa-mInfo->Hmax, 1.0, 1e-7, FALSE, mInfo );
Lgm_TraceLine2( &Pa, &Pb, (r-1.0)*Re, Sa/200.0, 1.0, 1e-7, FALSE, mInfo );
                    mInfo->Sm_South = 0.0;
                    mInfo->Sm_North = Sa;
                    ReplaceFirstPoint( 0.0, mInfo->Bm, &Pa, mInfo );
                    AddNewPoint( Sa,  mInfo->Bm, &Pb, mInfo );
                    InitSpline( mInfo );
                    mInfo->Lgm_n_I_integrand_Calls = 0;
                    InitSpline( mInfo );
mInfo->Lgm_I_integrand_S         = 0.0;                                                                                            
mInfo->Lgm_I_integrand_FirstCall = TRUE;                                                                                           
mInfo->Lgm_n_I_integrand_Calls    = 0;
                    I = Iinv_interped( mInfo  );
                    if (mInfo->VerbosityLevel > 2) printf("I = %g Lgm_n_I_integrand_Calls = %d\n", I, mInfo->Lgm_n_I_integrand_Calls );
printf("I = %g Lgm_n_I_integrand_Calls = %d\n", I, mInfo->Lgm_n_I_integrand_Calls );
//                    if (VerbosityLevel > 1) printf("\t\t%s  Integral Invariant, I (interped):      %15.8g    I-I0:    %15.8g    mlat:   %12.8lf  (nCalls = %d)%s\n",  PreStr, I, I-I0, b, mInfo->Lgm_n_I_integrand_Calls, PostStr );
                    FreeSpline( mInfo );
//                    a = 0.0; b = Sa+Sb;
//                    mInfo->Pm_South = Pa; mInfo->Sm_South = a;
//                    mInfo->Pm_North = Pb; mInfo->Sm_North = b;
//                    if (mInfo->VerbosityLevel > 2) printf("\t\tComputing Integral Invariant. Limits of integration are: [a,b] = [%g,%g]\n", a, b);
//
//                    mInfo->FirstCall = TRUE;
//                    mInfo->Lgm_n_I_integrand_Calls = 0;
//                    I = Iinv( mInfo );
//                    if (mInfo->VerbosityLevel > 2) printf("\t\tI = %g   ", I);
//                    if (mInfo->VerbosityLevel > 2) printf("Lgm_n_I_integrand_Calls = %d\n", mInfo->Lgm_n_I_integrand_Calls );


                } else {
                    printf("\t\tMirror point below %g km in Southern Hemisphere\n", mInfo->Lgm_LossConeHeight);
                    exit(0);
                }

            } else {
                printf("\t\tMirror point below %g km in Northern Hemisphere\n", mInfo->Lgm_LossConeHeight);
                exit(0);
            }


            f[i+N] = I;

        }

    }
    if (DIFF_SCHEME == USE_SIX_POINT){
        GradI->y = (f[6] - 9.0*f[5] + 45.0*f[4] - 45.0*f[2] + 9.0*f[1] - f[0])/(60.0*h);
    } else if (DIFF_SCHEME == USE_FOUR_POINT){
        GradI->y = (-f[4] + 8.0*f[3] - 8.0*f[1] + f[0])/(12.0*h);
    } else if (DIFF_SCHEME == USE_TWO_POINT){
        GradI->y = (f[2] - f[0])/(2.0*h);
    }


    /* Z-component */
    if (mInfo->VerbosityLevel > 0) printf("\t\tComputing dIdz: h = %g\n", h);
printf("\t\tComputing dIdz: h = %g\n", h);
    for (i=-N; i<=N; ++i){

        if (i!=0) { // dont need the center value in our difference scheme
    
            u = *v0; H = (double)i*h; u.z += H;

            /*
             * Trace to northern mirror point
             */
            //printf("E. Lgm_TraceToMirrorPoint\n");
            if ( Lgm_TraceToMirrorPoint( &u, &Pb, &Sb, mInfo->Bm, 1.0, mInfo->Lgm_TraceToMirrorPoint_Tol, mInfo ) > 0 ) {

                /*
                 * Trace to southern mirror point
                 */
                //printf("F. Lgm_TraceToMirrorPoint\n");
                if ( Lgm_TraceToMirrorPoint( &u, &Pa, &Sa, mInfo->Bm, -1.0, mInfo->Lgm_TraceToMirrorPoint_Tol, mInfo ) > 0 ) {

                    r  = Lgm_Magnitude( &mInfo->Pm_North );
//                    Lgm_TraceLine2( &Pa, &Pb, (r-1.0)*Re, 0.5*Sa-mInfo->Hmax, 1.0, 1e-7, FALSE, mInfo );
Lgm_TraceLine2( &Pa, &Pb, (r-1.0)*Re, Sa/200.0, 1.0, 1e-7, FALSE, mInfo );
                    mInfo->Sm_South = 0.0;
                    mInfo->Sm_North = Sa;
                    ReplaceFirstPoint( 0.0, mInfo->Bm, &Pa, mInfo );
                    AddNewPoint( Sa,  mInfo->Bm, &Pb, mInfo );
                    InitSpline( mInfo );
                    mInfo->Lgm_n_I_integrand_Calls = 0;
                    InitSpline( mInfo );
mInfo->Lgm_I_integrand_S         = 0.0;                                                                                            
mInfo->Lgm_I_integrand_FirstCall = TRUE;                                                                                           
mInfo->Lgm_n_I_integrand_Calls    = 0;
                    I = Iinv_interped( mInfo  );
                    if (mInfo->VerbosityLevel > 2) printf("I = %g Lgm_n_I_integrand_Calls = %d\n", I, mInfo->Lgm_n_I_integrand_Calls );
printf("I = %g Lgm_n_I_integrand_Calls = %d\n", I, mInfo->Lgm_n_I_integrand_Calls );
//                    if (VerbosityLevel > 1) printf("\t\t%s  Integral Invariant, I (interped):      %15.8g    I-I0:    %15.8g    mlat:   %12.8lf  (nCalls = %d)%s\n",  PreStr, I, I-I0, b, mInfo->Lgm_n_I_integrand_Calls, PostStr );
                    FreeSpline( mInfo );
//                    a = 0.0; b = Sa+Sb;
//                    mInfo->Pm_South = Pa; mInfo->Sm_South = a;
//                    mInfo->Pm_North = Pb; mInfo->Sm_North = b;
//                    if (mInfo->VerbosityLevel > 2) printf("\t\tComputing Integral Invariant. Limits of integration are: [a,b] = [%g,%g]\n", a, b);
//
//                    mInfo->FirstCall = TRUE;
//                    mInfo->Lgm_n_I_integrand_Calls = 0;
//                    I = Iinv( mInfo );
//                    if (mInfo->VerbosityLevel > 2) printf("\t\tI = %g   ", I);
//                    if (mInfo->VerbosityLevel > 2) printf("Lgm_n_I_integrand_Calls = %d\n", mInfo->Lgm_n_I_integrand_Calls );


                } else {
                    printf("\t\tMirror point below %g km in Southern Hemisphere\n", mInfo->Lgm_LossConeHeight);
                    exit(0);
                }

            } else {
                printf("\t\tMirror point below %g km in Northern Hemisphere\n", mInfo->Lgm_LossConeHeight);
                exit(0);
            }


            f[i+N] = I;

        }

    }
    if (DIFF_SCHEME == USE_SIX_POINT){
        GradI->z = (f[6] - 9.0*f[5] + 45.0*f[4] - 45.0*f[2] + 9.0*f[1] - f[0])/(60.0*h);
    } else if (DIFF_SCHEME == USE_FOUR_POINT){
        GradI->z = (-f[4] + 8.0*f[3] - 8.0*f[1] + f[0])/(12.0*h);
    } else if (DIFF_SCHEME == USE_TWO_POINT){
        GradI->z = (f[2] - f[0])/(2.0*h);
    }

    //printf("\t\tGrad_I = %g %g %g\n", GradI->x, GradI->y, GradI->z);

    return(0);

}







/*
 *    $Id: IntegralInvariant.c 141 2011-02-01 16:09:15Z mgh $
 */

