#include "Lgm/Lgm_PolyRoots.h"

/*
 *  Routines to compute roots of quadtratic, cubic and quartic equations exactly.
 *  See CRC standard math tables (Beyer) for details.
 */


/**
 *  \brief
 *      Returns the two roots of the quadratic equation.
 *
 *  \details
 *      Returns the two roots of the quadratic equation;
 *
 *              \f[a z^2 + b z + c = 0,\f]
 *
 *      where \f$a, b, c\f$ are real coefficients.
 *
 *      \param[in]      a   Real coefficient of the \f$z^2\f$ term.
 *      \param[in]      b   Real coefficient of the \f$z\f$ term.
 *      \param[in]      c   Real constant term.
 *      \param[out]     z1  First (possibly complex) root.
 *      \param[out]     z2  Second (possibly complex) root.
 *
 *      \return         The number of real roots
 *
 *      \author         Mike Henderson
 *      \date           2011
 *
 */
int  Lgm_QuadraticRoots( double a, double b, double c, double complex *z1, double complex *z2 ){

    int             nReal;
    double          x, D2 = b*b - 4.0*a*c;
    double complex  D;

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




/**
 *  \brief
 *      Returns the three roots of the cubic equation with real coefficients.
 *
 *  \details
 *      Returns the three roots of the cubic equation;
 *
 *              \f[z^3 + b z^2 + c z + d = 0,\f]
 *
 *      where \f$b, c, d\f$ are real coefficients, and the coefficient on the
 *      \f$x^3\f$ term is assumed to be 1.
 *
 *      \param[in]      b   Real coefficient of the \f$z^2\f$ term.
 *      \param[in]      c   Real coefficient of the \f$z\f$ term.
 *      \param[in]      d   Real constant term.
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
int  Lgm_CubicRoots( double b, double c, double d, double *z1, double complex *z2, double complex *z3 ){

    int     nReal;
    double  b2, Q, R, D2, D, RpD, RmD, S, T, Theta, SqrtNQ, f, g, h;

    b2 = b*b;

    Q  = (3.0*c - b2)/9.0;
    R  = (9.0*b*c - 27.0*d - 2.0*b2*b)/54.0;
    D2 = Q*Q*Q + R*R;

    if (D2 == 0){
        // All 3 roots are real and at least two are equal.
        if ( R > 0.0) {
            f = cbrt( R );
            *z1 = 2.0*f - b/3.0;
            *z2 = *z3 = -f - b/3.0;
        } else {
            f = -cbrt( -R );
            *z1 = 2.0*f - b/3.0;
            *z2 = *z3 = -f - b/3.0;
        }
        nReal = 3;
    } else if ( D2 > 0 ) {
        // only one real root exists. The other two are complex conjugates.
        D = sqrt(D2); RpD = R+D; RmD = R-D;
        S = (RpD > 0.0) ? cbrt( RpD ) : -cbrt( -RpD );
        T = (RmD > 0.0) ? cbrt( RmD ) : -cbrt( -RmD);
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
        *z2 = 2.0*SqrtNQ*cos( (Theta + 2.0*M_PI)/3.0 ) - b/3.0; // real (but returned as double complex z2)
        *z3 = 2.0*SqrtNQ*cos( (Theta + 4.0*M_PI)/3.0 ) - b/3.0; // real (but returned as double complex z3)
        nReal = 3;

    }

    return( nReal );

}

/**
 *  \brief
 *      Returns a real root of the cubic equation with real coefficients.
 *
 *  \details
 *      Returns a real root of the cubic equation;
 *
 *              \f[z^3 + b z^2 + c z + d = 0,\f]
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
 *      \param[in]      b   Real coefficient of the \f$z^2\f$ term.
 *      \param[in]      c   Real coefficient of the \f$z\f$ term.
 *      \param[in]      d   Real constant term.
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
            x1 = 2.0*cbrt( R ) - b/3.0;
        } else {
            x1 = -2.0*cbrt( -R ) - b/3.0;
        }
    } else if ( D2 > 0 ) {
        // only one real root exists. The other two are complex conjugates.
        D = sqrt(D2); RpD = R+D; RmD = R-D;
        S = (RpD > 0.0) ? cbrt( RpD ) : -cbrt( -RpD );
        T = (RmD > 0.0) ? cbrt( RmD ) : -cbrt( -RmD);
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







/**
 *  \brief
 *      Returns the four roots of the quartic equation with real coefficients.
 *
 *  \details
 *      Returns the four roots of the quartic equation;
 *
 *              \f[z^4 + b z^3 + c z^2 + d z + e = 0,\f]
 *
 *      where \f$b, c, d, e\f$ are real coefficients, and the coefficient on the
 *      \f$z^4\f$ term is assumed to be 1.
 *
 *      This rotuine uses Ferrari's' method. The idea is to recast the quartic
 *      in terms of a quadratic. In the process of doing this, a real root of a
 *      cubic equation needs to be obtained (hence the reason for
 *      Lgm_RealCubicRoot()).
 *
 *      \param[in]      b   Real coefficient of the \f$z^3\f$ term.
 *      \param[in]      c   Real coefficient of the \f$z^2\f$ term.
 *      \param[in]      d   Real coefficient of the \f$z\f$ term.
 *      \param[in]      e   Real constant term.
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

/**
 *  \brief
 *      Returns the four roots of the quartic equation with real coefficients. Results are sorted.
 *
 *  \details
 *      Returns the four roots of the quartic equation;
 *
 *              \f[z^4 + b z^3 + c z^2 + d z + e = 0,\f]
 *
 *      where \f$b, c, d, e\f$ are real coefficients, and the coefficient on the
 *      \f$z^4\f$ term is assumed to be 1.
 *
 *      This routine calls Lgm_QuarticRoots() and splits the results into Real and Complex roots.
 *
 *      \param[in]      b            Real coefficient of the \f$z^3\f$ term.
 *      \param[in]      c            Real coefficient of the \f$z^2\f$ term.
 *      \param[in]      d            Real coefficient of the \f$z\f$ term.
 *      \param[in]      e            Real constant term.
 *      \param[out]     nReal        Number of real roots found.
 *      \param[out]     RealRoots    Array of real roots.
 *      \param[out]     nComplex     Number of complex roots found.
 *      \param[out]     ComplexRoots Array of complex roots.
 *
 *      The RealRoots and ComplexRoots arrays need to be suitably allocated by
 *      the user and should have at least 4 elements each.
 *
 *      \return         The number of real roots found.
 *
 *      \author         Mike Henderson
 *      \date           2011
 *
 */
int Lgm_QuarticRootsSorted( double b, double c, double d, double e, int *nReal, double *RealRoots, int *nComplex, double complex *ComplexRoots ){

    int             nr=0, nc = 0;
    double complex  z1, z2, z3, z4;

    Lgm_QuarticRoots( b, c, d, e, &z1, &z2, &z3, &z4 );

    if ( cimag( z1 ) < 1e-16 ) { RealRoots[nr++] = creal(z1); } else { ComplexRoots[nc++] = z1; }
    if ( cimag( z2 ) < 1e-16 ) { RealRoots[nr++] = creal(z2); } else { ComplexRoots[nc++] = z2; }
    if ( cimag( z3 ) < 1e-16 ) { RealRoots[nr++] = creal(z3); } else { ComplexRoots[nc++] = z3; }
    if ( cimag( z4 ) < 1e-16 ) { RealRoots[nr++] = creal(z4); } else { ComplexRoots[nc++] = z4; }

    *nReal    = nr;
    *nComplex = nc;

    return( nr );

}



/*
 *  An arbitrary polynomial root solver written by C. Bond (see
 *  http://www.crbond.com ).
 *
 *  Finds all roots of polynomial by first finding quadratic factors using
 *  Bairstow's method, then extracting roots from quadratics. Implements new
 *  algorithm for managing multiple roots.  (C) 2002, 2003, C. Bond. All rights
 *  reserved.
 *
 */

#define LGM_POLYROOTS_MAXITER 500

/*
 *  Extract individual real or complex roots from list of quadratic factors.
 */
int Lgm_PolyRoots_Roots( double *a, int n, double *wr, double *wi ) {

    double  sq, b2, c, disc;
    int     i, m, numroots;

    m = n;
    numroots = 0;

    while (m > 1) {

        b2   = -0.5*a[m-2];
        c    = a[m-1];
        disc = b2*b2-c;
        if (fabs(disc)/(fabs(b2*b2)+fabs(c)) <= DBL_EPSILON) disc = 0.0;

        if (disc < 0.0) {
            sq = sqrt(-disc);
            wr[m-2] = b2;
            wi[m-2] = sq;
            wr[m-1] = b2;
            wi[m-1] = -sq;
            numroots+=2;
        } else {
            sq = sqrt(disc);
            wr[m-2] = fabs(b2)+sq;
            if (b2 < 0.0) wr[m-2] = -wr[m-2];
            if (wr[m-2] == 0) {
                wr[m-1] == 0;
            } else {
                wr[m-1] = c/wr[m-2];
                numroots+=2;
            }
            wi[m-2] = 0.0;
            wi[m-1] = 0.0;
        }

        m -= 2;

    }

    if (m == 1) {
       wr[0] = -a[0];
       wi[0] = 0.0;
       numroots++;
    }

    return( numroots );

}




/*
 *  Deflate polynomial 'a' by division of 'quad'. Return quotient
 *  polynomial in 'b' and error (metric based on remainder) in 'err'.
 */
void Lgm_PolyRootsDeflate( double *a, int n, double *b, double *quad, double *err ) {

    int     i;
    double  r, s;
    double  *c = (double *)calloc( n+1, sizeof(double) );

    r = quad[1];
    s = quad[0];

    b[1] = a[1] - r;
    c[1] = b[1] - r;

    for (i=2;i<=n;i++){
        b[i] = a[i] - r*b[i-1] - s*b[i-2];
        c[i] = b[i] - r*c[i-1] - s*c[i-2];
    }
    *err = fabs(b[n]) + fabs(b[n-1]);
    free( c );

}



/*
 *  Find quadratic factor using Bairstow's method (quadratic Newton method).
 *  A number of ad hoc safeguards are incorporated to prevent stalls due
 *  to common difficulties, such as zero slope at iteration point, and
 *  convergence problems.
 */
void Lgm_PolyRoots_FindQuad( double *a, int n, double *b, double *quad, double *err, int *iter ) {

    int     i;
    double  dn, dr, ds, drn, dsn, eps, r, s;
    double  *c = (double *)calloc( n+1, sizeof(double) );

    c[0]  = 1.0;
    r     = quad[1];
    s     = quad[0];
	dr    = 1.0;
	ds    = 0;
    eps   = 1e-15;
    *iter = 1;

	while ((fabs(dr)+fabs(ds)) > eps) {

        if (*iter > LGM_POLYROOTS_MAXITER) break;

        if ( ((*iter) % 200) == 0 ) eps*=10.0;

		b[1] = a[1] - r;
		c[1] = b[1] - r;

		for ( i=2; i<=n; i++ ){
			b[i] = a[i] - r*b[i-1] - s*b[i-2];
			c[i] = b[i] - r*c[i-1] - s*c[i-2];
		}

		dn  =c[n-1] * c[n-3] - c[n-2] * c[n-2];
		drn =b[n]   * c[n-3] - b[n-1] * c[n-2];
		dsn =b[n-1] * c[n-1] - b[n]   * c[n-2];

        if (fabs(dn) < 1e-10) {
            if (dn < 0.0) dn = -1e-8;
            else dn = 1e-8;
        }

        dr = drn / dn;
        ds = dsn / dn;

		r += dr;
		s += ds;

        (*iter)++;

	}

    quad[0] = s;
    quad[1] = r;
    *err = fabs(ds)+fabs(dr);

    free(c);

}



/*
 *  Differentiate polynomial 'a' returning result in 'b'.
 */
void Lgm_PolyRoots_DiffPoly( double *a, int n, double *b) {

    int     i;
    double  coef;

    coef = (double)n;
    b[0] = 1.0;

    for ( i=1; i<n; i++) b[i] = a[i]*((double)(n-i))/coef;

}



/*
 * Attempt to find a reliable estimate of a quadratic factor using modified
 * Bairstow's method with provisions for 'digging out' factors associated
 * with multiple roots.
 *
 * This resursive routine operates on the principal that differentiation of
 * a polynomial reduces the order of all multiple roots by one, and has no
 * other roots in common with it. If a root of the differentiated polynomial
 * is a root of the original polynomial, there must be multiple roots at
 * that location. The differentiated polynomial, however, has lower order
 * and is easier to solve.
 *
 * When the original polynomial exhibits convergence problems in the
 * neighborhood of some potential root, a best guess is obtained and tried
 * on the differentiated polynomial. The new best guess is applied
 * recursively on continually differentiated polynomials until failure
 * occurs. At this point, the previous polynomial is accepted as that with
 * the least number of roots at this location, and its estimate is
 * accepted as the root.
 */
void Lgm_PolyRoots_Recurse( double *a, int n, double *b, int m, double *quad, double *err, int *iter) {

    double *c, *x, rs[2], tst, e1, e2;

    if (fabs(b[m]) < 1e-16) m--;    // this bypasses roots at zero

    if (m == 2) {
        quad[0] = b[2];
        quad[1] = b[1];
        *err  = 0;
        *iter = 0;
        return;
    }

    c = (double *)calloc( m+1, sizeof(double) );
    x = (double *)calloc( n+1, sizeof(double) );
    c[0]  = x[0] = 1.0;
    rs[0] = quad[0];
    rs[1] = quad[1];
    *iter = 0;

    Lgm_PolyRoots_FindQuad( b, m, c, rs, err, iter );
    tst = fabs( rs[0]-quad[0] ) + fabs( rs[1]-quad[1] );

    if (*err < 1e-12) {
        quad[0] = rs[0];
        quad[1] = rs[1];
    }

    // tst will be 'large' if we converge to wrong root
    if ( ((*iter > 5) && (tst < 1e-4)) || ((*iter > 20) && (tst < 1e-1)) ) {
        Lgm_PolyRoots_DiffPoly( b, m, c );
        Lgm_PolyRoots_Recurse( a, n, c, m-1, rs, err, iter );
        quad[0] = rs[0];
        quad[1] = rs[1];
    }

    free(c);
    free(x);

}





/*
 *  Top level routine to manage the determination of all roots of the given
 *  polynomial 'a', returning the quadratic factors (and possibly one linear
 *  factor) in 'x'.
 */
void Lgm_PolyRoots_GetQuads( double *a, int n, double *quad, double *x ) {

    int     iter, i, m;
    double  *b, *z, err, tmp;
    double  xr, xs;

    if ( (tmp = a[0]) != 1.0 ) {
        a[0] = 1.0;
        for ( i=1; i<=n; i++ ) a[i] /= tmp;
    }

    if (n == 2) {

        x[0] = a[1];
        x[1] = a[2];
        return;

    } else if (n == 1) {

        x[0] = a[1];
        return;

    }


    m = n;
    b = (double *)calloc( n+1, sizeof(double) );
    z = (double *)calloc( n+1, sizeof(double) );
    b[0] = 1.0;

    for ( i=0; i<=n; i++ ) {
        z[i] = a[i];
        x[i] = 0.0;
    }


    do {

        if (n > m) {
            quad[0] = 3.14159e-1;
            quad[1] = 2.78127e-1;
        }

        do {

            Lgm_PolyRoots_FindQuad( z, m, b, quad, &err, &iter );
            if ( (err > 1e-7) || (iter > LGM_POLYROOTS_MAXITER) ) {
                Lgm_PolyRoots_DiffPoly( z, m, b );
                iter = 0;
                Lgm_PolyRoots_Recurse( z, m, b, m-1, quad, &err, &iter );
            }
            Lgm_PolyRootsDeflate( z, m, b, quad, &err );

        } while ( err > 1 );

        x[m-2] = quad[1];
        x[m-1] = quad[0];
        m -= 2;
        for (i=0;i<=m;i++) {
            z[i] = b[i];
        }

    } while (m > 2);


    if (m == 2) {
        x[0] = b[1];
        x[1] = b[2];
    } else {
        x[0] = b[1];
    }

    free(z);
    free(b);

}


/**
 *  \brief
 *      Find roots of arbitrary polynomials (with real coefficients).
 *
 *      Solves for the roots of the polynomial;
 *
 *          \f[ a_0 x^n + a_1 x^{n-1} + a_2 x^{n-2} + \cdots + a_{n-1} x + a_n = 0 \f]
 *
 *
 *  \details
 *
 *
 *      \param[in]      a   The polynomial coefficients. There should be n+1 of these.
 *      \param[in]      n   The order of the polynomial.
 *      \param[in]      z   The complex roots found.
 *
 *      \return         The number of roots found.
 *
 *
 */
int  Lgm_PolyRoots( double *a, int n, double complex *z ){

    int     numr, i;
    double  *wr, *wi, *x, quad[2];

    if (a[0] == 0) {
        printf("Lgm_PolyRoots(): Error! Highest coefficient cannot be 0.\n");
        return( 0 );
    }
    if (a[n] == 0) {
        printf("Lgm_PolyRoots(): Error! Lowest coefficient (constant term) cannot be 0.\n");
        return( 0 );
    }

    // initialize estimate for 1st root pair
    quad[0] = 2.71828e-1;
    quad[1] = 3.14159e-1;

    x  = (double *)calloc( n+1, sizeof(double) );
    wr = (double *)calloc( n+1, sizeof(double) );
    wi = (double *)calloc( n+1, sizeof(double) );

    // get roots
    Lgm_PolyRoots_GetQuads( a, n, quad, x );
    numr = Lgm_PolyRoots_Roots( x, n, wr, wi );

    for (i=0; i<numr; i++){
        z[i] = wr[i] + wi[i]*I;
    }

    return( numr );

}

