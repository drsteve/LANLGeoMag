/*
 * Examples showing how to use the Quadpack routines to do simple integrals.
 */
#include <stdio.h>
#include <math.h>
#include <Lgm_QuadPack.h>

typedef struct MyInfo {
    int a;
    int b;
} MyInfo;


/*
 *   Integrand for the following integral.  This has an integrable singularity at 0.
 *
 *      / 
 *      |
 *  I = |  ln( |x| )/ sqrt( |x| ) dx
 *      |
 *      / 
 *   
 *
 * Note that the definite integral from 0 to 1 has an endpoint singularity with solution;
 *
 *      / 1
 *      |
 *  I = |  ln( |x| )/ sqrt( |x| ) dx = -4
 *      |
 *      / 0
 *
 * Note that the definite integral from -1 to 1 has an interior singularity with solution;
 *
 *      / 1
 *      |
 *  I = |  ln( |x| )/ sqrt( |x| ) dx = -8
 *      |
 *      / -1
 *
 *
 */
double Integrand1( double x, _qpInfo *qpInfo ) {

    double  f;

    f = log(fabs(x))/sqrt(fabs(x));

    return( f );

}

/*
 *   Integrand for the following integral. 
 *
 *      / 
 *      |
 *  I = |  sin(ax) cos(bx) dx  = 2a/(a^2-b^2) 
 *      |
 *      / 
 *
 * Note that the definite integral from 0 to pi has the solution;
 *
 *      / pi
 *      |
 *  I = |  sin(ax) cos(bx) dx  = 2a/(a^2-b^2) if a-b is odd, or 0 if a-b is even
 *      |
 *      / 0
 *   
 *   
 */
double Integrand2( double x, _qpInfo *qpInfo ) {

    double  f;
    MyInfo *mi;

    /*
     * Get pointer to parameters for integrand.
     */
    mi = (MyInfo *)qpInfo;

    f = sin( mi->a * x ) * cos( mi->b * x );

    return( f );

}

int main(){

    int         limit=500, lenw=4*limit, iwork[502], last, ier, neval, VerbosityLevel, npts;
    double      epsabs, epsrel, abserr, work[2002], points[10];
    double      a, b, I, Iexact;
    MyInfo mi;


    /*
     * Integrable singularity (x=0) at one endpoint a=0 -- use dqags
     */
    mi.a = 2.0;     // Parameter needed in our integrand
    a =  0.0;       // lower limit of integration
    b =  1.0;       // upper limit of integration
    epsabs = 0.1;   // absolute error
    epsrel = 0.0;  // reletive error
    VerbosityLevel = 3;
    dqags( Integrand1, NULL, a, b, epsabs, epsrel, &I, &abserr, &neval, &ier, limit, lenw, &last, iwork, work, VerbosityLevel );
    Iexact = -4.0;
    printf("I = %.15lf    neval = %d (Exact answer = %.15lf)\n", I, neval, Iexact);




    /*
     * Internal integrable singularity (x=0)  -- use dqagp and tell it wehere the singularity is.
     */
    a = -1.0;       // lower limit of integration
    b =  1.0;       // upper limit of integration
    epsabs = 0.0;   // absolute error
    epsrel = 1e-1;  // reletive error
    VerbosityLevel = 3;
    npts = 3; points[0] = 0.0; // location (in x) of singularity
    dqagp( Integrand1, NULL, a, b, npts, points, epsabs, epsrel, &I, &abserr, &neval, &ier, limit, lenw, &last, iwork, work, VerbosityLevel );
    Iexact = -8.0;
    printf("I = %.15lf    neval = %d (Exact answer = %.15lf)\n", I, neval, Iexact);



    /*
     * Demonstration of how to get parmater/info into the integrand function.
     */
    a = 0.0;        // lower limit of integration
    b = M_PI;       // upper limit of integration
    epsabs = 0.0;   // absolute error
    epsrel = 1e-1;  // reletive error
    VerbosityLevel = 3;
    mi.a = 3;
    mi.b = 2;
    dqags( Integrand2, (_qpInfo *)&mi, a, b, epsabs, epsrel, &I, &abserr, &neval, &ier, limit, lenw, &last, iwork, work, VerbosityLevel );
    Iexact = (mi.a-mi.b)%2 ? 2.0*mi.a/((double)(mi.a*mi.a-mi.b*mi.b))  : 0.0;
    printf("I = %.15lf    neval = %d (Exact answer = %.15lf)\n", I, neval, Iexact);



    return(0);

}

