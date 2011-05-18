#include "Lgm/Lgm_PolyRoots.h"

/*
 *  Routines to compute roots of quadtratic, cubic and quartic equations exactly.
 *  See CRC standard math tables (Beyer) for details.
 */


/*
 *  \brief
 *      Returns the two roots of the quadratic equation.
 *
 *  \details
 *      Returns the two roots of the quadratic equation.
 *
 *              \f[a z^2 + b z + c = 0\f],
 *
 *      where \f$a, b, c\f$ are real coefficients.
 *
 *      \param[in]      a   Coefficient of the \f$z^2\f$ term.
 *      \param[in]      b   Coefficient of the \f$z\f$ term.
 *      \param[in]      c   Constant term.
 *      \param[out]     z1  First (possibly complex) root.
 *      \param[out]     z2  Second (possibly complex) root.
 *
 *      \return         The number of real roots
 *
 *      \author         Mike Henderson
 *      \date           2011
 *
 */
int  Lgm_QuadraticRoots( double a, double b, double c, complex double *z1, complex double *z2 ){

    int             nReal;
    double          x, D2 = b*b - 4.0*a*c;
    complex double  D;

    if ( fabs(D2) < 1e-10 ) {
        // roots are real and equal
        nReal = 2;
        x = -0.5*b/a;
        *z1 = *z2 = x;
    } else if ( D2 > 0.0 ) {
        // roots are real and unequal
        nReal = 2;
        D = sqrt(D2);
        x = (-b + D)/(2.0*a); *z1 = x;
        x = (-b - D)/(2.0*a); *z2 = x;
    } else {
        // roots are imaginary and unequal
        nReal = 0;
        D = csqrt(D2);
        *z1 = (-b + D)/(2.0*a);
        *z2 = (-b - D)/(2.0*a);
    }

    return(nReal);

}




/*
 *  \brief 
 *      Returns the three roots of the cubic equation with real coefficients.
 *
 *  \details 
 *      Returns the three roots of the cubic equation.
 *
 *              \f[z^3 + b z^2 + c z + d = 0\f],
 *
 *      where \f$b, c, d\f$ are real coefficients, and the coefficient on the
 *      \f$x^3\f$ term is assumed to be 1.
 *
 *      \param[in]      b   Coefficient of the \f$z^2\f$ term.
 *      \param[in]      c   Coefficient of the \f$z\f$ term.
 *      \param[in]      d   Constant term.
 *      \param[out]     z1  First real root.
 *      \param[out]     z2  Second (possibly complex) root.
 *      \param[out]     z3  Third (possibly complex) root.
 *
 *      \return         The number of real roots
 *
 *      \author         Mike Henderson
 *      \date           2011
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
        nReal = 3;
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
        nReal = 1;

    } else if ( D2 < 0 ) {
        // All three roots are real and unequal.
        // This uses the simplest of the alternative trig forms.
        SqrtNQ = sqrt(-Q);
        Theta = acos( R/(SqrtNQ*fabs(Q)) );
        *z1 = 2.0*SqrtNQ*cos(Theta/3.0) - b/3.0; // real (returned as real double z1)
        *z2 = 2.0*SqrtNQ*cos( (Theta + 2.0*M_PI)/3.0 ) - b/3.0; // real (but returned as complex double z2)
        *z3 = 2.0*SqrtNQ*cos( (Theta + 4.0*M_PI)/3.0 ) - b/3.0; // real (but returned as complex double z3)
        nReal = 3;

    }

    return( nReal );

}

/*
 *
 *  Return a real root for the cubic equation;
 *
 *       z^3 + b z^2 + c z + d = 0
 *
 *
 */

/*
 *  \brief 
 *      Returns a real root of the cubic equation with real coefficients.
 *
 *  \details 
 *      Returns a real root of the cubic equation.
 *
 *              \f[z^3 + b z^2 + c z + d = 0\f],
 *
 *      where \f$b, c, d\f$ are real coefficients, and the coefficient on the
 *      \f$4^3\f$ term is assumed to be 1.
 *
 *      Note that there is always at least one real root of any cubic.  This
 *      routine chooses a single (arbitrary) real root. A real root is needed
 *      to solve the "resolvent cubic" in the exact solution of the quartic
 *      (see Lgm_QuarticRoots()). If you need all three roots, use
 *      Lgm_CubicRoots() instead.
 *
 *      \param[in]      b   Coefficient of the \f$z^2\f$ term.
 *      \param[in]      c   Coefficient of the \f$z\f$ term.
 *      \param[in]      d   Constant term.
 *
 *      \return         A single real root of the cubic equation.
 *
 *      \author         Mike Henderson
 *      \date           2011
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

/*
 *  \brief 
 *      Returns the four roots of the quartic equation with real coefficients.
 *
 *  \details 
 *      Returns the four roots of the quartic equation;
 *
 *              \f[z^4 + b z^3 + c z^2 + d z + e = 0\f],
 *
 *      where \f$b, c, d, e\f$ are real coefficients, and the coefficient on the
 *      \f$z^4\f$ term is assumed to be 1.
 *
 *      This routine
 *
 *      \param[in]      b   Coefficient of the \f$z^3\f$ term.
 *      \param[in]      c   Coefficient of the \f$z^2\f$ term.
 *      \param[in]      d   Coefficient of the \f$z\f$ term.
 *      \param[in]      e   Constant term.
 *      \param[out]     z1  First (possibly complex) root.
 *      \param[out]     z2  First (possibly complex) root.
 *      \param[out]     z3  First (possibly complex) root.
 *      \param[out]     z4  First (possibly complex) root.
 *
 *      \return         The number of real roots found.
 *
 *      \author         Mike Henderson
 *      \date           2011
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
