#include <math.h>
#include <complex.h>
#include <stdio.h>
#include <stdlib.h>

/*
 *
 *  Returns the three roots of the cubic equation;
 *
 *       z^3 + b z^2 + c z + d = 0
 *
 *  where b, c, d are real coefficients. There will always be at least one real
 *  root (and as many as 3).
 *
 *
 */
int  Lgm_CubicRoots( double b, double c, double d, double *z1, complex double *z2, complex double *z3 ){

    int     nReal;
    double  b2, Q, R, D2, D, RpD, RmD, S, T, Theta, SqrtNQ, f, g, h;

    b2 = b*b;

    Q  = (3.0*c - b2)/9.0;
    R  = (9.0*b*c - 27.0*d - 2.0*b2*b)/54.0;
    D2 = Q*Q*Q + R*R;

    if (D2 == 0){
        // All 3 roots are real and at least two are equal.
        if ( R > 0.0) {
            f = pow( R, 1.0/3.0 );
            *z1 = 2.0*f - b/3.0;
            *z2 = *z3 = -f - b/3.0;
        } else {
            f = -pow( -R, 1.0/3.0 );
            *z1 = 2.0*f - b/3.0;
            *z2 = *z3 = -f - b/3.0;
        }
    } else if ( D2 > 0 ) {
        // only one real root exists. The other two are complex conjugates.
        D = sqrt(D2); RpD = R+D; RmD = R-D;
        S = (RpD > 0.0) ? pow( RpD, 1.0/3.0 ) : -pow( -RpD, 1.0/3.0 );
        T = (RmD > 0.0) ? pow( RmD, 1.0/3.0 ) : -pow( -RmD, 1.0/3.0);
        f = S+T;
        g = sqrt(3.0)/2.0*(S-T);
        h = b/3.0;
        *z1 = f - h;
        *z2 = -0.5*f - h + g*I;
        *z3 = -0.5*f - h - g*I;

    } else if ( D2 < 0 ) {
        // All three roots are real and unequal.
        // This uses the simplest of the alternative trig forms.
        SqrtNQ = sqrt(-Q);
        Theta = acos( R/(SqrtNQ*fabs(Q)) );
        *z1 = 2.0*SqrtNQ*cos(Theta/3.0) - b/3.0; // real (returned as real double z1)
        *z2 = 2.0*SqrtNQ*cos( (Theta + 2.0*M_PI)/3.0 ) - b/3.0; // real (but returned as complex double z2)
        *z3 = 2.0*SqrtNQ*cos( (Theta + 4.0*M_PI)/3.0 ) - b/3.0; // real (but returned as complex double z3)

    } 

    nReal = 1; // there will always be at least one real root.
    if ( fabs(cimag(*z2)) < 1e-10 ) ++nReal;
    if ( fabs(cimag(*z3)) < 1e-10 ) ++nReal;

    return( nReal );

}

/*
 *
 *  Return a real root for the cubic equation;
 *
 *       z^3 + b z^2 + c z + d = 0
 *
 *  Note that there is always at least one real root of any cubic.  This
 *  routine chooses a single (arbitrary) real root. A real root is needed to
 *  solve the "resolvent cubic" in the exact solution of the quartic. If you
 *  need all three roots, use Lgm_CubicRoots() instead.
 *
 */
double  Lgm_CubicRealRoot( double b, double c, double d ){

    double  b2, Q, R, D2, D, RpD, RmD, S, T, x1, Theta, SqrtNQ;

    b2 = b*b;

    Q  = (3.0*c - b2)/9.0;
    R  = (9.0*b*c - 27.0*d - 2.0*b2*b)/54.0;
    D2 = Q*Q*Q + R*R;

    if (D2 == 0){
        // All 3 roots are real and at least two are equal.
        if ( R > 0.0) {
            x1 = 2.0*pow( R, 1.0/3.0 ) - b/3.0;
        } else {
            x1 = -2.0*pow( -R, 1.0/3.0 ) - b/3.0;
        }
    } else if ( D2 > 0 ) {
        // only one real root exists. The other two are complex conjugates.
        D = sqrt(D2); RpD = R+D; RmD = R-D;
        S = (RpD > 0.0) ? pow( RpD, 1.0/3.0 ) : -pow( -RpD, 1.0/3.0 );
        T = (RmD > 0.0) ? pow( RmD, 1.0/3.0 ) : -pow( -RmD, 1.0/3.0);
        x1 = S + T - b/3.0;
    } else if ( D2 < 0 ) {
        // All three roots are real and unequal.
        // This uses the simplest of the alternative trig forms.
        SqrtNQ = sqrt(-Q);
        Theta = acos( R/(SqrtNQ*fabs(Q)) );
        x1 = 2.0*SqrtNQ*cos(Theta/3.0) - b/3.0;
    } 

    return( x1 );

}







/*
 *
 *  Solve the general quartic equation;
 *
 *       z^4 + b z^3 + c z^2 + d z + e = 0
 *  
 *  where b, c, d, e are real coefficients.
 *
 *  This rotuine uses Ferrari's' method. The idea is to recast the quartic in
 *  terms of a quadratic. In the process of doing this, a real root of a cubic
 *  equation needs to be obtained (hence the reason for Lgm_RealCubicRoot()).
 *
 */
int Lgm_QuarticRoots( double b, double c, double d, double e, double complex *z1, double complex *z2, double complex *z3, double complex *z4 ){

    int             nReal;
    double          p, q, r, y1, v, R2, b2, b3, f;
    double complex  D, E, g, R, s, t;


    /*
     * Obtain a real root from the "resolvent cubic".
     *    y^3 + p y^2 + q y + r = 0
     */
    p  = -c;
    q  = (b*d - 4.0*e);
    r  = 4.0*c*e - b*b*e - d*d;
    y1 = Lgm_CubicRealRoot( p, q, r );

    
    /*
     * Construct the quadratic and solve. (Must use complex math of course.)
     */
    b2 = b*b; 
    b3 = b*b2;
    R2 = b2/4.0 - c + y1; 
    R  = csqrt( R2 ); // R will in general be complex

    /*
     *  If R is small do the first one. Its fairly sensitive here...
     */
    if ( fabs(R2) < 1e-7 ) {
        v = y1*y1 - 4.0*e; 
        f = 0.75*b2 - 2.0*c;
        g = 2.0*csqrt( v );
        D = csqrt( f + g );
        E = csqrt( f - g );
        s = t = -0.5*b;
    } else {
        f = 0.75*b2 - R2 - 2.0*c;
        g = 0.25*(4.0*b*c - 8.0*d - b3)/R;
        D = csqrt( f + g );
        E = csqrt( f - g );
        s = -0.5*b + R; 
        t = -0.5*b - R;
    }

    *z1 = 0.5*( s + D);
    *z2 = 0.5*( s - D);
    *z3 = 0.5*( t + E);
    *z4 = 0.5*( t - E);

    nReal = 0;
    if ( fabs(cimag(*z1)) < 1e-10 ) ++nReal;
    if ( fabs(cimag(*z2)) < 1e-10 ) ++nReal;
    if ( fabs(cimag(*z3)) < 1e-10 ) ++nReal;
    if ( fabs(cimag(*z4)) < 1e-10 ) ++nReal;

    return( nReal );
    
}





int main( int argc, char *argv[] ){
    
    int             nReal=0, nErrors=0;
    double          b, c, d, e, x, f;
    double complex  z1, z2, z3, z4, cf1, cf2, cf3, cf4;

    for (b = -10.0; b< 10.0; b += 2.0 ){
        for (c = -10.0; c< 10.0; c += 2.0 ){
            for (d = -10.0; d< 10.0; d += 2.0 ){
                x = Lgm_CubicRealRoot( b, c, d );
                f = x*x*x + b*x*x + c*x + d;
                if ( fabs(f) > 1e-10 ) {
                    ++nErrors;
                    printf("b, c, d = %g %g %g   x = %lf   f = %g\n", b, c, d, x, f );
                }
            }
        }
    }

    for (b = -10.0; b< 10.0; b += 2.0 ){
        for (c = -10.0; c< 10.0; c += 2.0 ){
            for (d = -10.0; d< 10.0; d += 2.0 ){
                nReal = Lgm_CubicRoots( b, c, d, &x, &z2, &z3 );
                  f = x*x*x + b*x*x + c*x + d;
                cf2 = z2*z2*z2 + b*z2*z2 + c*z2 + d;
                cf3 = z3*z3*z3 + b*z3*z3 + c*z3 + d;
                if ( (fabs(f) > 1e-10) || (cabs(cf2) > 1e-10) || (cabs(cf3) > 1e-10) ) {
                    ++nErrors;
                    printf("nReal = %d\n", nReal );
                    printf("b, c, d = %g %g %g   x = %g         f = %g\n", b, c, d, x, f );
                    printf("b, c, d = %g %g %g  z2 = %g+%gi   cf2 = %g+%gi |cf1| = %g\n", b, c, d, creal(z2), cimag(z2), creal(cf2), cimag(cf2), cabs(cf2)  );
                    printf("b, c, d = %g %g %g  z3 = %g+%gi   cf3 = %g+%gi |cf1| = %g\n", b, c, d, creal(z3), cimag(z3), creal(cf3), cimag(cf3), cabs(cf3)  );
                }
            }
        }
    }

    for (b = -10.0; b< 10.0; b += 1.0 ){
        for (c = -10.0; c< 10.0; c += 1.0 ){
            for (d = -10.0; d< 10.0; d += 1.0 ){
                for (e = -10.0; e< 10.0; e += 1.0 ){
                    nReal = Lgm_QuarticRoots( b, c, d, e, &z1, &z2, &z3, &z4 );
                    cf1 = z1*z1*z1*z1 + b*z1*z1*z1 + c*z1*z1 + d*z1 + e;
                    cf2 = z2*z2*z2*z2 + b*z2*z2*z2 + c*z2*z2 + d*z2 + e;
                    cf3 = z3*z3*z3*z3 + b*z3*z3*z3 + c*z3*z3 + d*z3 + e;
                    cf4 = z4*z4*z4*z4 + b*z4*z4*z4 + c*z4*z4 + d*z4 + e;
                    if ( (cabs(cf1) > 1e-6) || (cabs(cf2) > 1e-6) || (cabs(cf3) > 1e-6) || (cabs(cf4) > 1e-6) ) {
                        printf("p = -c = %g\nq = (b*d - 4.0*e) = %g\nr = 4.0*c*e - b*b*e - d*d = %g\n", -c, (b*d - 4.0*e), 4.0*c*e - b*b*e - d*d);;
                        printf("FAC = %g\n", 4.0*b*c - 8.0*d - b*b*b);
                        ++nErrors;
                        printf("nReal = %d\n", nReal );
                        printf("b, c, d, e = %g %g %g %g  z1 = %g+%gi   cf1 = %g+%gi |cf1| = %g\n", b, c, d, e, creal(z1), cimag(z1), creal(cf1), cimag(cf1), cabs(cf1)  );
                        printf("b, c, d, e = %g %g %g %g  z2 = %g+%gi   cf2 = %g+%gi |cf1| = %g\n", b, c, d, e, creal(z2), cimag(z2), creal(cf2), cimag(cf2), cabs(cf2)  );
                        printf("b, c, d, e = %g %g %g %g  z3 = %g+%gi   cf3 = %g+%gi |cf1| = %g\n", b, c, d, e, creal(z3), cimag(z3), creal(cf3), cimag(cf3), cabs(cf3)  );
                        printf("b, c, d, e = %g %g %g %g  z4 = %g+%gi   cf4 = %g+%gi |cf1| = %g\n\n", b, c, d, e, creal(z4), cimag(z4), creal(cf4), cimag(cf4), cabs(cf4)  );
                    }
                }
            }
        }
    }
    
    printf("nErrors = %d\n", nErrors);


    return(0);

}
