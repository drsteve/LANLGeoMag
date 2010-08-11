#include <stdio.h>
#include <math.h>

/*
 *  Hilton's approximation to compute McIlwain's L-shell parameter:
 *   
 *   		L^3 Bm / M = F( I^3 Bm / M )
 *
 *  This routine uses the F() function as described in 
 *  Hilton, JGR 76, No 28, 6952-6954, October 1, 1971.
 *
 */
double	LFromIBmM_Hilton( double I, double Bm, double M ) {

    double	p, a1, a2, a3, X, X13, X23, F;


    p = 1.0/3.0;
    X   = I*I*I*Bm/M;
    X13 = pow( X, p );
    X23 = X13*X13;

    a1 = 1.35047;
    a2 = 0.465376;
    a3 = 0.0475455;

    F = 1.0 + a1*X13 + a2*X23 + a3*X;

    return( pow(F*M/Bm, p)  );

}



/*
 *  Iterative inversion of Hilton
 */
double	IFromLBmM_Hilton( double L, double Bm, double M ) {

    int		done;
    double	LHS, RHS, p, a1, a2, a3, X13, X23, X, a, b, c, d, I;

    /*
     *  we have;
     *   
     *		L^3Bm/M - 1.0 = ax^1/3 + bx^2/3 + cx 
     *
     *	or;
     * 
     * 		LHS = RHS
     * 
     *  where;
     *   
     *   	LHS = L^3Bm/M - 1.0 
     *   	RHS = ax^1/3 + bx^2/3 + cx
     *
     *  Need to find X such that RHS = LHS
     *
     */

    LHS = L*L*L*Bm/M - 1.0;
    p = 1.0/3.0;
    a1 = 1.35047;
    a2 = 0.465376;
    a3 = 0.0475455;


    /* 
     *  Set bounds. Cant go any lower than 0.0 for X
     *  step until we have an upper bound
     */
    a = 0.0;
    done = 0; c = a;
    while (!done) {
	X = c; X13 = pow( X, p ); X23 = X13*X13;
	RHS = a1*X13 + a2*X23 + a3*X;

	if ( RHS - LHS > 0.0 ) done = 1;
	else if (c == 0.0 )    	    c = 1.5;
	else 		       	    c *= 10.0; 
    }

    if ( fabs(a-c) < 1e-14 ) return( 0.0 );


    /* 
     *  We have a bracket. Zero in with bisection
     */
    done = 0; 
    while (!done) {
	d = c-a;
	b = a + 0.5*d;

	X = b; X13 = pow( X, p ); X23 = X13*X13;
	RHS = a1*X13 + a2*X23 + a3*X;
	
	if ( fabs((RHS-LHS)/LHS) < 1e-8 ) { done = 1; }
	else if ( RHS - LHS > 0.0 ) { c = b; }
	else { a = b; }
    }

    /*
     *  Take average of endpoints as final answer
     */
    X = 0.5*(a+c); I = pow( X*M/Bm, p );
    return( I );

}








/*
 *  McIlwain's original approximation to compute his L-shell parameter:
 *   
 *   		L^3 Bm / M = F( I^3 Bm / M )
 *
 *  This routine uses the F() function as described in 
 *  Roederer's 1970 book, Appendix VI, page 155.
 *
 */
double	LFromIBmM_McIlwain( double I, double Bm, double M ) {

    double	X, X2, X3, X4, X5, X6, X7, X8, X9, Y, G;

    X  = log(I*I*I*Bm/M);

    if        ( X < -22.0 ) { 	/* X < -22.0 */

	Y = 3.0062102e-1 + 3.333389e-1*X;

    } else if ( X < -3.0 ) { 	/* -22.0 <= X < -3.0 */

        X2 = X*X; X3 = X*X2; X4 = X2*X2; X5 = X2*X3; X6 = X3*X3; X7 = X3*X4; X8 = X4*X4; X9 = X4*X5;
	Y = 6.2337691e-1 + 4.3432642e-1*X + 1.5017245e-2*X2 + 1.3714667e-3*X3 + 8.2711096e-5*X4
		+ 3.2916354e-6*X5 + 8.1048663e-8*X6 + 1.0066362e-9*X7 + 8.3232531e-13*X8 - 8.1537735e-14*X9;

    } else if ( X <  3.0 ) { 	/* -3.0 <= X < 3.0 */

        X2 = X*X; X3 = X*X2; X4 = X2*X2; X5 = X2*X3; X6 = X3*X3; X7 = X3*X4; X8 = X4*X4; X9 = X4*X5;
	Y = 6.228644e-1 + 4.3352788e-1*X + 1.4492441e-2*X2 + 1.1784234e-3*X3 + 3.8379917e-5*X4
		- 3.3408822e-6*X5 - 5.3977642e-7*X6 - 2.1997983e-8*X7 + 2.3028767e-9*X8 + 2.6047023e-10*X9;

    } else if ( X < 11.7 ) { 	/* 3.0 <= X < 11.7 */

        X2 = X*X; X3 = X*X2; X4 = X2*X2; X5 = X2*X3; X6 = X3*X3; X7 = X3*X4; X8 = X4*X4; X9 = X4*X5;
	Y = 6.222355e-1 + 4.3510529e-1*X + 1.2817956e-2*X2 + 2.1680398e-3*X3 - 3.2077032e-4*X4 
		+ 7.9451313e-5*X5 - 1.2531932e-5*X6 + 9.9766148e-7*X7 - 3.958306e-8*X8 + 6.3271665e-10*X9;

    } else if ( X < 23.0 ) { 	/* 11.7 <= X < 23.0 */

        X2 = X*X; X3 = X*X2; X4 = X2*X2; X5 = X2*X3; X6 = X3*X3; 
	Y = 2.0007187 - 1.8461796e-1*X + 1.2038224e-1*X2 - 6.7310339e-3*X3 + 2.170224e-4*X4
		- 3.8049276e-6*X5 + 2.8212095e-8*X6;

    } else { 			/* X >= 23.0 */

	Y = -3.0460681 + X;

    }

    G = exp( Y ) + 1.0;
    return( pow( G*M/Bm  , 1.0/3.0 ) );


}

/*
 *  Inverse of above
 */
double	IFromLBmM_McIlwain( double L, double Bm, double M ) {

    double	Y, Y2, Y3, Y4, Y5, Y6, Y7, Y8, Y9, X, H, f;

    f = L*L*L*Bm/M - 1.0;
    if ( f < 0.0 ) {
	return( 0.0 );
    } else {
        Y  = log(L*L*L*Bm/M - 1.0);
    }

    if        ( Y < -7.0 ) { 	/* Y < -7.0 */

	X = -9.0176314e-1 + 2.999966*Y;

    } else if ( Y < -0.7 ) { 	/* -7.0 <= Y < -0.7 */

        Y2 = Y*Y; Y3 = Y*Y2; Y4 = Y2*Y2; Y5 = Y2*Y3; Y6 = Y3*Y3; Y7 = Y3*Y4; Y8 = Y4*Y4; Y9 = Y4*Y5;
	X = -1.5035469 + 2.516886*Y - 1.5473094e-1*Y2 - 1.482533e-2*Y3 + 5.7665822e-3*Y4
		+ 2.5165553e-3*Y5 + 4.6975588e-4*Y6 + 4.9336212e-5*Y7 + 2.8332675e-6*Y8 + 6.9447496e-8*Y9;

    } else if ( Y <  2.0 ) { 	/* -0.7 <= Y < 2.0 */

        Y2 = Y*Y; Y3 = Y*Y2; Y4 = Y2*Y2; Y5 = Y2*Y3; Y6 = Y3*Y3; Y7 = Y3*Y4; Y8 = Y4*Y4; Y9 = Y4*Y5;
	X = -1.5035665 + 2.5166432*Y - 1.5579012e-1*Y2 - 1.6991733e-2*Y3 + 3.3726437e-3*Y4
		+ 9.4374468e-4*Y5 - 1.6048742e-4*Y6 - 5.7606054e-5*Y7 + 1.6653982e-5*Y8 - 1.11444463e-6*Y9;

    } else if ( Y < 8.5 ) { 	/* 2.0 <= Y < 8.5 */

        Y2 = Y*Y; Y3 = Y*Y2; Y4 = Y2*Y2; Y5 = Y2*Y3; Y6 = Y3*Y3; Y7 = Y3*Y4; Y8 = Y4*Y4; Y9 = Y4*Y5;
	X = -1.4921674 + 2.4799708*Y - 1.0499454e-1*Y2 - 5.6314931e-2*Y3 + 2.16412e-2*Y4
		- 3.9980282e-3*Y5 + 4.5340102e-4*Y6 - 3.2040569e-5*Y7 + 1.3002586e-6*Y8 - 2.3206968e-8*Y9;

    } else if ( Y < 20.0 ) { 	/* 8.5 <= Y < 20.0 */

        Y2 = Y*Y; Y3 = Y*Y2; Y4 = Y2*Y2; Y5 = Y2*Y3; Y6 = Y3*Y3; Y7 = Y3*Y4; Y8 = Y4*Y4; Y9 = Y4*Y5;
	X = -1.3415203 + 2.5733781*Y - 2.747596e-1*Y2 + 2.9961512e-2*Y3 - 2.1911161e-3*Y4
		+ 1.0870528e-4*Y5 - 3.5717807e-6*Y6 + 7.2957464e-8*Y7 - 8.0754426e-10*Y8 + 3.3797513e-12*Y9;

    } else { 			/* Y >= 20.0 */

	X = 3.0460681 + Y;

    }

    H = exp( X );
    return( pow( H*M/Bm  , 1.0/3.0 ) );


}


