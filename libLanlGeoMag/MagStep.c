/* Copyright (c) 1999 Michael G. Henderson <mghenderson@lanl.gov>
 *
 *     Routines to perform Bulirsch-Stoer step.
 *
 * Permission to use, copy, modify, distribute, and sell this software and its
 * documentation for any purpose is hereby granted without fee, provided that
 * the above copyright notice appear in all copies and that both that
 * copyright notice and this permission notice appear in supporting
 * documentation.  No representations are made about the suitability of this
 * software for any purpose.  It is provided "as is" without express or
 * implied warranty.
 *
 */
#include "Lgm/Lgm_MagModelInfo.h"


void Lgm_ModMid( Lgm_Vector *u, Lgm_Vector *v, double H, int n, double sgn, 
         int (*Mag)(Lgm_Vector *, Lgm_Vector *, Lgm_MagModelInfo *), Lgm_MagModelInfo *Info ) {

    int            m;
    double        h2, h;
    Lgm_Vector    z0, z1, z2, B;


    /*
     *  Set stepsize, h and 2*h.
     */
    h  = sgn * H/(double)n;
    h2 = 2.0 * h;


    /*
     *  Set initial point, z0.
     */
    z0 = *u;


    /*
     *  Do initial Euler step to get z1.
     */
    (*Mag)(&z0, &B, Info); Lgm_NormalizeVector(&B);
    z1.x = z0.x + h*B.x; 
    z1.y = z0.y + h*B.y; 
    z1.z = z0.z + h*B.z;
    

    /*
     *  Do general step to get z2 -> zn. This is midpoint formula.
     */
    for ( m = 1; m < n; ++m ) {

        (*Mag)(&z1, &B, Info); Lgm_NormalizeVector(&B);
        z2.x = z0.x + h2*B.x; 
        z2.y = z0.y + h2*B.y; 
        z2.z = z0.z + h2*B.z;
        z0 = z1;
        z1 = z2;

    }


    /*
     *  Do final Euler step to get z(n+1).
     */
    (*Mag)(&z1, &B, Info); Lgm_NormalizeVector(&B);
    z2.x = z1.x + h*B.x; 
    z2.y = z1.y + h*B.y; 
    z2.z = z1.z + h*B.z;



    /*
     *   The final answer for zn is the average of z(n-1) and z(n+1). 
     */
    v->x = 0.5 * (z0.x + z2.x); 
    v->y = 0.5 * (z0.y + z2.y); 
    v->z = 0.5 * (z0.z + z2.z);
    


}


// NOT thread safe?
void Lgm_RatFunExt( int k, double x_k, Lgm_Vector *u_k, Lgm_Vector *w, Lgm_Vector *dw, Lgm_MagModelInfo *Info ) {

    int             i, j;
    double             yy, v, ddy=0.0, c, b1, b, fx[LGM_MAGSTEP_JMAX];
    double            y_k[3], y[3], dy[3];
//    static double   d[LGM_MAGSTEP_JMAX][LGM_MAGSTEP_JMAX], x[LGM_MAGSTEP_JMAX];

    /*
     *  (x_k, y_k) is the kth "data point" in the sequence.
     */
    y_k[0] = u_k->x; y_k[1] = u_k->y; y_k[2] = u_k->z;

    Info->Lgm_MagStep_x[k] = x_k;
    if (k == 0) {
        for (i=0; i<3; ++i) y[i] = dy[i] = Info->Lgm_MagStep_d[i][0] = y_k[i];
    } else {
        for (j=0; j<k; ++j) fx[j+1] = Info->Lgm_MagStep_x[k-j-1] / x_k;
        for (i=0; i<3; ++i) {
            v = Info->Lgm_MagStep_d[i][0];
            c = yy = Info->Lgm_MagStep_d[i][0] = y_k[i];
            for (j=1; j<=k; ++j) {
            b1 = fx[j] * v;
            b  = b1 - c;
            if (b != 0.0) {
                b    = (c - v)/b;
                ddy  = c * b;
                c    = b1 * b;
            } else{
                ddy = v;
            }
            if (j != k) v = Info->Lgm_MagStep_d[i][j];
            Info->Lgm_MagStep_d[i][j]  = ddy;
            yy      += ddy;
            }
            dy[i] = ddy;
            y[i]  = yy;
        }
    }

    /*
     *  (x_k, u_k) are the input "data points". The rational function should go through
     *  all of these points. w is the rational function evaluated at x_k = 0.0 -- i.e.
     *  an estimate of what the sequence of modified midpoint estimates would converge to
     *  for infinitely small step sizes. dw is the error in this estimate.
     *
     */
    w->x  = y[0];    w->y  = y[1];    w->z  = y[2];
    dw->x = dy[0];   dw->y = dy[1];   dw->z = dy[2];
}




/*
 *
 *
 *  This routine is basically the one in Num. Rec.  The bulk of the code here is
 *  trying to figure out what the optimal value should be for the next step.  It
 *  uses a method developed by Deuflhard, 1985 based on information theory.
 *  However, Shampine, 1987 points out that a much simpler heuristic often does
 *  alot better. The problem with Deuflhard's method is apparently twofold: first
 *  Shampine claims that the information theory is not (necessarily) valid for any
 *  given problem (and when it comes down to it we really do have only one problem
 *  to solve). And second, even if it was applicable, Shampine claims that it is
 *  "the answer" to the wrong question -- i.e. it doesnt address the real
 *  problem.  Well....  We'll stick with this one for a while, but Shampine's
 *  approach sounds better, and the stats he shows indicate that significantly
 *  fewer function evals are needed...
 *
 *
 */
// NOT thread safe?
int Lgm_MagStep( Lgm_Vector *u, Lgm_Vector *u_scale, 
          double Htry, double *Hdid, double *Hnext, 
          double eps, double sgn, double *s, int *reset,
          int (*Mag)(Lgm_Vector *, Lgm_Vector *, Lgm_MagModelInfo *), Lgm_MagModelInfo *Info ){


    Lgm_Vector        u_save, v, uerr, e;
    int                q, k, kk, km=0, n;
    int                reduction, done;
    double            h2, sss, n2, f, H, err[LGM_MAGSTEP_KMAX+1];
    double            eps1, max_error=0.0, fact, red=1.0, scale=1.0, work, workmin;
//static int        Info->Lgm_MagStep_FirstTimeThrough=TRUE, Info->Lgm_MagStep_kmax, Info->Lgm_MagStep_kopt;
//static double   Info->Lgm_MagStep_eps_old=-1.0, Info->Lgm_MagStep_snew;
//static double   A[JMAX+1], Info->Lgm_MagStep_alpha[IMAX+1][IMAX+1];
    /*
     *  This sequence seems to be pretty good in general. The "1" seems to be
     *  really essential to speed things up (on some test runsa, its almost a factor of
     *  (2 slower without it!).
     */
    static int     Seq[] = { 0, 1, 2, 4, 6, 8, 12, 18, 24, 32, 48, 64 };



    if ( ( eps != Info->Lgm_MagStep_eps_old ) || ( *reset ) ){

        *Hnext  = -1.0;
        eps1    = LGM_MAGSTEP_SAFE1*eps;
        Info->Lgm_MagStep_eps_old = eps;

        /*
         *  Compute work to get to column k of the tableau.
         *  Here we assume work is just the number of function
         *  evals needed to get to column k.
         */
        Info->Lgm_MagStep_A[1] = Seq[1] + 1;
        for (k=1; k <= LGM_MAGSTEP_KMAX; ++k) Info->Lgm_MagStep_A[k+1] = Info->Lgm_MagStep_A[k] + Seq[k+1];


        /*
         *  Compute Deuflhard's correction factors.
         */
        for (q=2; q <= LGM_MAGSTEP_KMAX; ++q) {
            for (k=1; k < q; ++k) {
                Info->Lgm_MagStep_alpha[k][q] = pow( eps1, (Info->Lgm_MagStep_A[k+1] - Info->Lgm_MagStep_A[q+1]) 
                / ( (2.0*k+1.0)*(Info->Lgm_MagStep_A[q+1] - Info->Lgm_MagStep_A[1] + 1.0)) );
            }
        }

        /*
         *  Compute optimal row for convergence.
         */
        for (Info->Lgm_MagStep_kopt=2; Info->Lgm_MagStep_kopt < LGM_MAGSTEP_KMAX; Info->Lgm_MagStep_kopt++) 
            if (Info->Lgm_MagStep_A[Info->Lgm_MagStep_kopt+1] > Info->Lgm_MagStep_alpha[Info->Lgm_MagStep_kopt-1][Info->Lgm_MagStep_kopt] * Info->Lgm_MagStep_A[Info->Lgm_MagStep_kopt]) break;
            Info->Lgm_MagStep_kmax = Info->Lgm_MagStep_kopt;

    }




    H = Htry;
    u_save = *u;
    if ( (*s != Info->Lgm_MagStep_snew) || (H != *Hnext) || ( *reset ) ) {
        Info->Lgm_MagStep_snew = 0.0;
        *s   = 0.0;
        Info->Lgm_MagStep_kopt = Info->Lgm_MagStep_kmax;
        Info->Lgm_MagStep_FirstTimeThrough = TRUE;
    }
    reduction = FALSE;
    done      = FALSE;


    while ( !done ) {

        for (k=1; k<=Info->Lgm_MagStep_kmax; ++k) {

            Info->Lgm_MagStep_snew = *s + H;
            if (Info->Lgm_MagStep_snew == (*s)) { 
                printf("H = %g\n", H);
                fprintf(stderr, "step size underflow\n");
                fprintf(stderr, "Htry, Hdid, Hnext = %g %g %g\n", Htry, *Hdid, *Hnext);
                Info->Lgm_MagStep_FirstTimeThrough=TRUE;
                Info->Lgm_MagStep_eps_old = -1.0;
printf("HOW DID I GET HERE? P = \n");
//exit(0);
                return(-1);
            }

            n  = Seq[k];
            n2 = (double)(n*n);
            h2 = H*H;
            Lgm_ModMid( &u_save, &v, H, n, sgn, Mag, Info );
            sss = h2/n2;
            Lgm_RatFunExt( k-1, sss, &v, u, &uerr, Info );
            
            if (k !=  1){

                e.x = fabs( uerr.x/u_scale->x );
                e.y = fabs( uerr.y/u_scale->y );
                e.z = fabs( uerr.z/u_scale->z );
                max_error = 1.0e-30;
                if (e.x > max_error) max_error = e.x;
                if (e.y > max_error) max_error = e.y;
                if (e.z > max_error) max_error = e.z;
                max_error /= eps;
                km = k-1;
                err[km] = pow( max_error/LGM_MAGSTEP_SAFE1, 1.0/(2.0*km + 1.0) );
            }


            if ( (k != 1) && ((k >= Info->Lgm_MagStep_kopt-1) || Info->Lgm_MagStep_FirstTimeThrough )) {

                if ( max_error < 1.0 ) {

                    /*
                     *  We've converged! Bailout and go home...
                     */
                    done = TRUE;

                } else if ( (k == Info->Lgm_MagStep_kmax)||(k == Info->Lgm_MagStep_kopt+1) ) {

                    red = LGM_MAGSTEP_SAFE2/err[km];

                } else if ( (k == Info->Lgm_MagStep_kopt) && (Info->Lgm_MagStep_alpha[Info->Lgm_MagStep_kopt-1][Info->Lgm_MagStep_kopt] < err[km]) ) {

                    red = 1.0/err[km];

                } else if ( (k == Info->Lgm_MagStep_kmax) && (Info->Lgm_MagStep_alpha[km][Info->Lgm_MagStep_kmax-1] < err[km]) ) {

                    red = LGM_MAGSTEP_SAFE2 * Info->Lgm_MagStep_alpha[km][Info->Lgm_MagStep_kmax-1]/err[km];

                } else if ( Info->Lgm_MagStep_alpha[km][Info->Lgm_MagStep_kopt] < err[km] ) {

                    red = Info->Lgm_MagStep_alpha[km][Info->Lgm_MagStep_kopt-1]/err[km];

                }

            }


        }

        if (!done) {
            red = (red < LGM_MAGSTEP_REDMIN) ? red : LGM_MAGSTEP_REDMIN;
            red = (red > LGM_MAGSTEP_REDMIN) ? red : LGM_MAGSTEP_REDMIN;
            H *= red;
            reduction = TRUE;
        }

    }

    *s     = Info->Lgm_MagStep_snew;
    *Hdid  = H;
    *reset = FALSE;
    Info->Lgm_MagStep_FirstTimeThrough = FALSE;
    workmin = 1e99;
    for (kk=1; kk<=km; ++kk) {
        fact = (err[kk] > LGM_MAGSTEP_SCLMAX) ? err[kk] : LGM_MAGSTEP_SCLMAX;

        work = fact*Info->Lgm_MagStep_A[kk+1];
        if (work < workmin){
            scale   = fact;
            workmin = work;
            Info->Lgm_MagStep_kopt    = kk+1;
        }

    }

    *Hnext = H/scale;
    if ( (Info->Lgm_MagStep_kopt >= k) && (Info->Lgm_MagStep_kopt != Info->Lgm_MagStep_kmax) && !reduction ) {
        f = scale/Info->Lgm_MagStep_alpha[Info->Lgm_MagStep_kopt-1][Info->Lgm_MagStep_kopt];
        fact = (f > LGM_MAGSTEP_SCLMAX) ? f : LGM_MAGSTEP_SCLMAX;
        if ( Info->Lgm_MagStep_A[Info->Lgm_MagStep_kopt+1]*fact <= workmin){
            *Hnext = H/fact;
            Info->Lgm_MagStep_kopt++;
        }
    }



    return(1);

}

