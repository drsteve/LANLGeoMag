#include "Lgm/Lgm_MagModelInfo.h"
#define ITMAX 100
#define CGOLD 0.3819660
#define ZEPS 1.0e-16
#define EPS  1.0e-16
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

/*
 *  Sa, Sb, Sc are distances along the field line that bracket the Bmin value
 *  Pa, Pb, Pc are the 3D vector points that correspond to the Sa, Sb, Sc vals.
 *
 *  Smin, Pmin are the final result
 */
int Lgm_BrentP(double Sa, double Sb, double Sc, double Bb, Lgm_Vector Pa, Lgm_Vector Pb, Lgm_Vector Pc, BrentFuncInfoP *f, double tol, double *Smin, double *Bmin, Lgm_Vector *Pmin ) {

    Lgm_Vector  P, Px, Pw, Pv, Pu, Btmp;
    int         iter;
    double      a, b, d, etemp, fu, fv, fw, fx, p, q, r, tol1, tol2, u, v, w, x, xm;
    double      e = 0.0;
    double      Htry, Hdid, Hnext, du, s, B;

    a = (Sa < Sc ? Sa : Sc);
    b = (Sa > Sc ? Sa : Sc);

    x  = w = v = Sb;
    Px = Pw = Pv = Pb;
    fx = fw = fv = Bb;

    for ( iter=1; iter<=ITMAX; iter++ ) {

        xm = 0.5*(a+b);

        tol1 = tol*fabs(x);
        tol2 = 2.0*(tol1+ZEPS);

        // Quit if we are done
        if ( fabs(x-xm) <= (tol2-0.5*(b-a)) ) {
            *Smin = x;
            *Bmin = fx;
            *Pmin = Px;
            return(1);
        }

        if ( fabs(e) > tol1 ) {
            r = (x-w)*(fx-fv);
            q = (x-v)*(fx-fw);
            p = (x-v)*q-(x-w)*r;
            q = 2.0*(q-r);
            if (q > 0.0) p = -p;
            q = fabs(q);
            etemp = e;
            e = d;
            if ( fabs(p) >= fabs(0.5*q*etemp ) || (p <= q*(a-x)) || (p >= q*(b-x)) ) {
                e = (x >= xm) ? a-x : b-x;
                d = CGOLD*e;
            } else {
                d = p/q;
                u = x+d;
                if ( (u-a < tol2) || (b-u < tol2) ) d = SIGN( tol1, xm-x );
            }
        } else {
            e = x >= xm ? a-x : b-x;
            d = CGOLD*e;
        }


        /*
         * Evaluate Func at u
         */
        du = ( fabs(d) >= tol1 ) ? d : SIGN( tol1, d );
        u  = x + du;
        // x to u. (From point Px to Pu.)
        P    = Px;
        Htry = du;
        if ( Htry > 1e-16 ) {
            Lgm_MagStep( &P, &f->u_scale, Htry, &Hdid, &Hnext, f->sgn, &s, &f->reset, f->Info->Bfield, f->Info );
        }
        f->Info->Bfield( &P, &Btmp, f->Info );
        B = Lgm_Magnitude( &Btmp );
        Pu = P;
        fu = B;


        if ( fu <= fx ) {

            if ( u >= x ) { a = x; } else { b = x; }
            SHFT(v,w,x,u);
            SHFT(Pv,Pw,Px,Pu);
            SHFT(fv,fw,fx,fu);

        } else {

            if ( u < x ) { a = u; } else { b = u; }
            if ( (fu <= fw) || (w == x) ) {
                v  = w;
                w  = u;
                Pv = Pw;
                Pw = Pu;
                fv = fw;
                fw = fu;
            } else if ( (fu <= fv) ||  (v == x) ||  (v == w) ) {
                v  = u;
                Pv = Pu;
                fv = fu;
            }

        }

    }
    printf("Lgm_BrentP(): Too many iterations in brent.\n");
    *Smin = x;
    *Bmin = fx;
    *Pmin = Px;
    return(0);

}





int Lgm_zBrentP(double S1, double S2, double F1, double F2, Lgm_Vector P1, Lgm_Vector P2, BrentFuncInfoP *f, double tol, double *Sz, double *Fz, Lgm_Vector *Pz ) {

	int         iter, Count;
	double      a, b, c, d, e, min1, min2;
	double      fa, fb, fc, p, q, r, s, tol1, xm;
    double      Htry, Hdid, Hnext, dd, sgn, htry, hdid;
    Lgm_Vector  Pa, Pb, Pc, P;

    a = S1; Pa = P1; fa = F1;
    b = S2; Pb = P2; fb = F2;
    c = S2; Pc = P2; fc = F2;

    sgn = f->sgn;

    if ( fabs(fa) < tol ) {
        *Sz = S1;
        *Fz = F1;
        *Pz = P1;
    } else if ( fabs(fb) < tol ) {
        *Sz = S2;
        *Fz = F2;
        *Pz = P2;
    }


	if ( ((fa > 0.0) && (fb > 0.0)) || ((fa < 0.0) && (fb < 0.0)) ) {
		fprintf(stderr, "Root not bracketed in Lgm_zBrent, fa, fb = %g %g tol = %g\n", fa, fb, tol);
        return(0);
    }

	for ( iter=1; iter<=ITMAX; iter++ ) {

		if ( ((fb > 0.0) && (fc > 0.0)) || ((fb < 0.0) && (fc < 0.0)) ) {
			c  = a;
            Pc = Pa;
			fc = fa;
			e = d = b-a;
		}

		if ( fabs(fc) < fabs(fb) ) {
			a  = b;  b  = c;  c  = a;
            Pa = Pb; Pb = Pc; Pc = Pa;
			fa = fb; fb = fc; fc = fa;
		}

        // Check if we are done
		tol1 = 2.0*EPS*fabs(b)+0.5*tol;
		xm = 0.5*(c-b);
		if ( (fabs(xm) <= tol1) || (fb == 0.0) ) {
		//if ( (fabs(fb) <= tol) || (fb == 0.0) ) {
//printf("xm = %g    fb = %g\n", xm, fb);
//printf("c, b   = %lf %lf\n", c, b);
//printf("fc, fb = %g %g\n", fc, fb);
            *Sz = b;
            *Fz = fb;
            *Pz = Pb;
            return(1);
        }

		if ( (fabs(e) >= tol1) && (fabs(fa) > fabs(fb)) ) {
			s = fb/fa;
			if ( a == c ) {
				p = 2.0*xm*s;
				q = 1.0-s;
			} else {
				q = fa/fc;
				r = fb/fc;
				p = s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
				q = (q-1.0)*(r-1.0)*(s-1.0);
			}
			if (p > 0.0) q = -q;
			p = fabs(p);
			min1 = 3.0*xm*q-fabs(tol1*q);
			min2 = fabs(e*q);
			if ( 2.0*p < ( (min1 < min2) ? min1 : min2) ) {
				e = d;
				d = p/q;
			} else {
				d = xm;
				e = d;
			}
		} else {
			d = xm;
			e = d;
		}

		a  = b;
        Pa = Pb;
		fa = fb;

		if (fabs(d) > tol1) {
            dd = d;
		} else {
            dd = SIGN(tol1,xm);
        }
        b += dd;

		//fb = (*func)(b);
//        if ( dd < 0.0 ) sgn *= -1.0;
//        Htry = fabs(dd);
        Htry = dd;


        /*
         * We want to make sure that we actually do a step of Htry, so keep trying until we get there.
         */
//printf("\n\n\n########################\nIn brent: Htry = %g\n", Htry);
        htry = Htry, Hdid = 0.0; Count = 0;
        while ( (fabs(htry) >= 0.5*tol1 ) && (fabs(xm) >= tol1) && (Count<100) ) {
            if ( Lgm_MagStep( &Pb, &f->u_scale, htry, &hdid, &Hnext, sgn, &s, &f->reset, f->Info->Bfield, f->Info ) < 0 ) { printf("BAILING 4\n");return(-1); }
//printf("In brent: Count=%d htry, hdid = %g %g  Htry, Hdid = %g %g    fabs(Hdid-Htry) = %g\n", Count, htry, hdid, Htry, Hdid, fabs(Hdid-Htry) );
            Hdid += hdid;
            htry = Htry - Hdid;
//printf("In brent: Count=%d htry, hdid = %g %g  Htry, Hdid = %g %g    fabs(Hdid-Htry) = %g\n\n", Count, htry, hdid, Htry, Hdid, fabs(Hdid-Htry) );
            ++Count;
        }


if (Htry != Hdid) printf("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA  Htry = %g  Hdid = %g  AAAAAAAa\n", Htry, Hdid );
//printf("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA  Htry = %g  Hdid = %g  AAAAAAAa\n", Htry, Hdid );
        fb = f->func( &Pb, f->Val, f->Info );
        //printf("fa, fb, fc = %g %g %g\n", fa, fb, fc);

	}
	printf("Lgm_zBrent(): Maximum number of iterations exceeded, iter=%d. \n", iter);
    *Sz = b;
    *Fz = fb;
    *Pz = Pb;
	return(0);

}


int Lgm_zBrent(double S1, double S2, double F1, double F2, BrentFuncInfo *f, double tol, double *Sz, double *Fz ) {

	int         iter;
	double      a, b, c, d, e, min1, min2;
	double      fa, fb, fc, p, q, r, s, tol1, xm;
    double      Htry, Hdid, Hnext, dd, sgn;
    Lgm_Vector  Pa, Pb, Pc, P;

    a = S1; fa = F1;
    b = S2; fb = F2;
    c = S2; fc = F2;


    if ( fabs(fa) < tol ) {
        *Sz = S1;
        *Fz = F1;
    } else if ( fabs(fb) < tol ) {
        *Sz = S2;
        *Fz = F2;
    }


	if ( ((fa > 0.0) && (fb > 0.0)) || ((fa < 0.0) && (fb < 0.0)) ) {
//printf("fa, fb = %g %g\n", fa, fb);
		fprintf(stderr, "Root not bracketed in Lgm_zBrent, tol = %g\n", tol);
        return(0);
    }

    fc = fb;
	for ( iter=1; iter<=ITMAX; iter++ ) {

		if ( ((fb > 0.0) && (fc > 0.0)) || ((fb < 0.0) && (fc < 0.0)) ) {
			c  = a;
			fc = fa;
			e = d = b-a;
		}

		if ( fabs(fc) < fabs(fb) ) {
			a  = b;  b  = c;  c  = a;
			fa = fb; fb = fc; fc = fa;
		}

        // Check if we are done
		tol1 = 2.0*EPS*fabs(b)+0.5*tol;
		xm = 0.5*(c-b);
		if ( (fabs(xm) <= tol1) || (fb == 0.0) ) {
//printf("c, b   = %lf %lf\n", c, b);
//printf("fc, fb = %g %g\n", fc, fb);
            *Sz = b;
            *Fz = fb;
            return(1);
        }

		if ( (fabs(e) >= tol1) && (fabs(fa) > fabs(fb)) ) {
			s = fb/fa;
			if ( a == c ) {
				p = 2.0*xm*s;
				q = 1.0-s;
			} else {
				q = fa/fc;
				r = fb/fc;
				p = s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
				q = (q-1.0)*(r-1.0)*(s-1.0);
			}
			if (p > 0.0) q = -q;
			p = fabs(p);
			min1 = 3.0*xm*q-fabs(tol1*q);
			min2 = fabs(e*q);
			if ( 2.0*p < ( (min1 < min2) ? min1 : min2) ) {
				e = d;
				d = p/q;
			} else {
				d = xm;
				e = d;
			}
		} else {
			d = xm;
			e = d;
		}

		a  = b;
		fa = fb;

		if (fabs(d) > tol1) {
            dd = d;
		} else {
            dd = SIGN(tol1,xm);
        }
        b += dd;

        fb = f->func( b, f->Val, f->Info );
        //printf("a, b, c, fa, fb, fc = %g %g %g   %g %g %g\n", a, b, c, fa, fb, fc);

	}
	printf("Lgm_zBrent(): Maximum number of iterations exceeded\n");
    *Sz = b;
    *Fz = fb;
	return(0);

}



int Lgm_Brent(double xa, double xb, double xc, BrentFuncInfo *fInfo, double tol, double *xmin, double *fmin ) {

    int         iter;
    double      a, b, d, etemp, fu, fv, fw, fx, p, q, r, tol1, tol2, u, v, w, x, xm;
    double      e = 0.0;
    double      Htry, Hdid, Hnext, du, s, B;

    a = (xa < xc ? xa : xc);
    b = (xa > xc ? xa : xc);
    x = w = v = xb;
    fw = fv = fx = fInfo->func( xb, 0.0, fInfo->Info );

    for ( iter=1; iter<=ITMAX; iter++ ) {

        xm = 0.5*(a+b);

        tol1 = tol*fabs(x);
        tol2 = 2.0*(tol1+ZEPS);

        // Quit if we are done
        if ( fabs(x-xm) <= (tol2-0.5*(b-a)) ) {
            *xmin = x;
            *fmin = fx;
            return(1);
        }

        if ( fabs(e) > tol1 ) {
            r = (x-w)*(fx-fv);
            q = (x-v)*(fx-fw);
            p = (x-v)*q-(x-w)*r;
            q = 2.0*(q-r);
            if (q > 0.0) p = -p;
            q = fabs(q);
            etemp = e;
            e = d;
            if ( fabs(p) >= fabs(0.5*q*etemp ) || (p <= q*(a-x)) || (p >= q*(b-x)) ) {
                e = (x >= xm) ? a-x : b-x;
                d = CGOLD*e;
            } else {
                d = p/q;
                u = x+d;
                if ( (u-a < tol2) || (b-u < tol2) ) d = SIGN( tol1, xm-x );
            }
        } else {
            e = x >= xm ? a-x : b-x;
            d = CGOLD*e;
        }


        /*
         * Evaluate Func at u
         */
        u = ( fabs(d) >= tol1 ) ? x+d : x+SIGN( tol1, d );
        fu = fInfo->func( u, 0.0, fInfo->Info );

        if ( fu <= fx ) {

            if ( u >= x ) { a = x; } else { b = x; }
            SHFT(v,w,x,u);
            SHFT(fv,fw,fx,fu);

        } else {

            if ( u < x ) { a = u; } else { b = u; }
            if ( (fu <= fw) || (w == x) ) {
                v  = w;
                w  = u;
                fv = fw;
                fw = fu;
            } else if ( (fu <= fv) ||  (v == x) ||  (v == w) ) {
                v  = u;
                fv = fu;
            }

        }

    }
    printf("Lgm_Brent(): Too many iterations in brent.\n");
    *xmin = x;
    *fmin = fx;
    return(0);

}
