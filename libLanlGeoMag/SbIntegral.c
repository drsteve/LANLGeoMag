#include "Lgm/Lgm_MagModelInfo.h"
#define JUMP_METHOD	1







/*
 *   This routine evaluates the "Sb integral" from mirror point to
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
double SbIntegral( Lgm_MagModelInfo *fInfo ) {


    double	a, b;
    double	epsabs, epsrel, result, abserr;
    int		npts, key, neval, ier, leniw, lenw, last, iwork[501];
    double	work[2001], points[20];
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
    dqagp(Sb_integrand, qpInfo, a, b, npts, points, epsabs, epsrel, &result, &abserr, &neval, &ier, leniw, lenw, &last, iwork, work);
    //free( iwork );
    //free( work );

    //dqk21(Sb_integrand, qpInfo, a, b, &result, &abserr, &resabs, &resasc);


    return( result );

}


double SbIntegral_interped( Lgm_MagModelInfo *fInfo ) {


    double	a, b;
    double	epsabs, epsrel, result, abserr;
    int		npts, key, neval, ier, leniw, lenw, last, iwork[501];
    double	work[2001], points[20];
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


    points[1] = a;
    points[2] = b;
    npts = 4;
    dqagp(Sb_integrand_interped, qpInfo, a, b, npts, points, epsabs, epsrel, &result, &abserr, &neval, &ier, leniw, lenw, &last, iwork, work);

    return( result );

}



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



    if ( JUMP_METHOD == 0 ) {

        /*
         *  Set starting point.
         *  If this is the first call for this integral, set point
         *  to the lower limit.
         */
        if ( fInfo->Lgm_Sb_integrand_FirstCall  == TRUE ) {
            fInfo->Lgm_Sb_integrand_FirstCall = FALSE;
	        fInfo->Lgm_Sb_integrand_u_scale.x =  10.0;  fInfo->Lgm_Sb_integrand_u_scale.y = 1.0; fInfo->Lgm_Sb_integrand_u_scale.z = 10.0;
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
         */
        if ( fInfo->Lgm_Sb_integrand_FirstCall  == TRUE ) {
    	    fInfo->Lgm_Sb_integrand_FirstCall = FALSE;
        }
	    fInfo->Lgm_Sb_integrand_u_scale.x =  10.0;  fInfo->Lgm_Sb_integrand_u_scale.y = 1.0; fInfo->Lgm_Sb_integrand_u_scale.z = 10.0;
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

        Lgm_MagStep( &fInfo->Lgm_Sb_integrand_P, &fInfo->Lgm_Sb_integrand_u_scale, Htry, &Hdid, &Hnext, 1.0e-7, sgn, &sdid, &reset, fInfo->Bfield, fInfo );
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



