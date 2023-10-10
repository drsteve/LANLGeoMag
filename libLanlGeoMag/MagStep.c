/* Michael G. Henderson <mghenderson@lanl.gov>
 *
 *     Routines to perform Bulirsch-Stoer step.
 *
 * BS is based on NR 3rd edition routines.
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
#include <math.h>
#include "Lgm/Lgm_MagModelInfo.h"

#define FMAX(a,b)  (((a)>(b))?(a):(b))
#define FMIN(a,b)  (((a)<(b))?(a):(b))

/* 
 * NOTE: Two very influential parameters are the tolerances atol, rtol used in
 *       these routines.  rtol seems less useful than atol, so it is typically
 *       set to zero. Long ago, atol was set to a fairly small value (e.g. 1e-9
 *       or 1e-10, These days the default seems to be 1e-5 -- which isnt that
 *       accurate.)
 */
int Lgm_MagStep( Lgm_Vector *u, Lgm_Vector *u_scale, double Htry, double *Hdid,
double *Hnext, double sgn, double *s, int *reset, int (*Mag)(Lgm_Vector *,
Lgm_Vector *, Lgm_MagModelInfo *), Lgm_MagModelInfo *Info ){

    Lgm_Vector u0;
    double  eps;
    int status = 1;

    u0 = *u;
    if (        Info->Lgm_MagStep_Integrator == LGM_MAGSTEP_ODE_BS ) {

        eps = Info->Lgm_MagStep_BS_Eps;
        status = Lgm_MagStep_BS( u, u_scale, Htry, Hdid, Hnext, eps, sgn, s, reset, Mag, Info );

    } else if ( Info->Lgm_MagStep_Integrator == LGM_MAGSTEP_ODE_RK5 ) {

        eps = Info->Lgm_MagStep_RK5_Eps;
        status = Lgm_MagStep_RK5( u, u_scale, Htry, Hdid, Hnext, eps, sgn, s, reset, Mag, Info );

    } else {

        printf("Lgm_MagStep: Error. Unknown ODE solver. Info->Lgm_MagStep_Integrator must be either LGM_MAGSTEP_ODE_BS or LGM_MAGSTEP_ODE_RK5\n");
        return(-1);

    }


    return(status);


}





int Lgm_ModMid( Lgm_Vector *u, Lgm_Vector *b0, Lgm_Vector *v, double H, int n, double sgn,
         int (*Mag)(Lgm_Vector *, Lgm_Vector *, Lgm_MagModelInfo *), Lgm_MagModelInfo *Info ) {

    int            m;
    double        h2, h, Bmag;
    Lgm_Vector    z0, z1, z2, B;


    /*
     *  Set stepsize, h and 2*h.
     */
    h  = sgn * H/(double)n;
    h2 = 2.0 * h;
//printf("h=%g\n", h);


    /*
     *  Set initial point, z0.
     */
    z0 = *u;


    /*
     *  Do initial Euler step to get z1. We get b0 from the arg list (its computed once).
     */
    z1.x = z0.x + h*b0->x;
    z1.y = z0.y + h*b0->y;
    z1.z = z0.z + h*b0->z;


    /*
     *  Do general step to get z2 -> zn. This is midpoint formula.
     */
    for ( m = 1; m < n; ++m ) {

        if ( (*Mag)(&z1, &B, Info) == 0 ) {
            // bail if B-field eval had issues.
            printf("Lgm_ModMid(): B-field evaluation during midpoint phase (m = %d and z1 = %g %g %g) returned with errors (returning with 0)\n", m, z1.x, z1.y, z1.z );
            return(0);
        }
        ++(Info->Lgm_nMagEvals);
        Bmag = Lgm_NormalizeVector(&B);
        if ( Bmag < 1e-16 ) {
            // bail if B-field magnitude is too small
            printf("Lgm_ModMid(): Bmag too small during midpoint phase (m = %d and z1 = %g %g %g Bmag = %g) is too small (returning with 0).\n", m, z1.x, z1.y, z1.z, Bmag );
            return(0);
        }
        z2.x = z0.x + h2*B.x;
        z2.y = z0.y + h2*B.y;
        z2.z = z0.z + h2*B.z;
        z0 = z1;
        z1 = z2;

    }


    /*
     *  Do final Euler step to get z(n+1).
     */
    if ( (*Mag)(&z1, &B, Info) == 0 ) {
        // bail if B-field eval had issues.
        printf("Lgm_ModMid(): B-field evaluation during final Euler step (z1 = %g %g %g) returned with errors (returning with 0)\n", z1.x, z1.y, z1.z );
        return(0);
    }
    ++(Info->Lgm_nMagEvals);
    Bmag = Lgm_NormalizeVector(&B);
    if ( Bmag < 1e-16 ) {
        // bail if B-field magnitude is too small
        printf("Lgm_ModMid(): Bmag too small during final Euler step (z1 = %g %g %g B = %g %g %g Bmag = %g) is too small (returning with 0).\n", z1.x, z1.y, z1.z, B.x, B.y, B.z, Bmag );
//printf("z0, z1 = %g %g %g    %g %g %g   z2 = %g %g %g\n", z0.x, z0.y, z0.z, z1.x, z1.y, z1.z, z2.x, z2.y, z2.z );

        return(0);
    }
    z2.x = z1.x + h*B.x;
    z2.y = z1.y + h*B.y;
    z2.z = z1.z + h*B.z;



    /*
     *   The final answer for zn is the average of z(n-1) and z(n+1).
     */
    v->x = 0.5 * (z0.x + z2.x);
    v->y = 0.5 * (z0.y + z2.y);
    v->z = 0.5 * (z0.z + z2.z);
//printf("z0, z1 = %g %g %g    %g %g %g   v = %g %g %g\n", z0.x, z0.y, z0.z, z2.x, z2.y, z2.z, v->x, v->y, v->z );


    return(1);

}


// NOT thread safe? I think it is....
void Lgm_RatFunExt( int k, double x_k, Lgm_Vector *u_k, Lgm_Vector *w, Lgm_Vector *dw, Lgm_MagModelInfo *Info ) {

    int             i, j;
    double             yy, v, ddy=0.0, c, b1, b, fx[LGM_MAGSTEP_JMAX];
    double            y_k[3], y[3], dy[3];
//    static double   d[LGM_MAGSTEP_JMAX][LGM_MAGSTEP_JMAX], x[LGM_MAGSTEP_JMAX];

    /*
     *  (x_k, y_k) is the kth "data point" in the sequence.
     */
    y_k[0] = u_k->x; y_k[1] = u_k->y; y_k[2] = u_k->z;

    Info->Lgm_MagStep_BS_x[k] = x_k;
    // next line is not needed -- just debugginh...
    Info->Lgm_MagStep_BS_u[k] = *u_k;
    if (k == 0) {
        for (i=0; i<3; ++i) y[i] = dy[i] = Info->Lgm_MagStep_BS_d[i][0] = y_k[i];
    } else {
        for (j=0; j<k; ++j) fx[j+1] = Info->Lgm_MagStep_BS_x[k-j-1] / x_k;
        for (i=0; i<3; ++i) {
            v = Info->Lgm_MagStep_BS_d[i][0];
            c = yy = Info->Lgm_MagStep_BS_d[i][0] = y_k[i];
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
                if (j != k) v = Info->Lgm_MagStep_BS_d[i][j];
                Info->Lgm_MagStep_BS_d[i][j]  = ddy;
                yy      += ddy;
            }
            dy[i] = ddy;
            y[i]  = yy;
        }
    }

//printf("Lgm_RatFunExt:\n");
//for (j=0; j<=k; j++){
//printf("k = %d  x,y = %.15lf    %.15lf %.15lf %.15lf\n", k, Info->Lgm_MagStep_BS_x[j], Info->Lgm_MagStep_BS_u[j].x, Info->Lgm_MagStep_BS_u[j].y, Info->Lgm_MagStep_BS_u[j].z);
//}

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



void Lgm_PolFunExt( int k, double x_k, Lgm_Vector *u_k, Lgm_Vector *w, Lgm_Vector *dw, Lgm_MagModelInfo *Info ) {

    int             k1, j;
    double          q, f1, f2, c[LGM_MAGSTEP_JMAX];
    double          y_k[3], y[3], dy[3], delta;

    /*
     *  (x_k, y_k) is the kth "data point" in the sequence.
     */
    y_k[0] = u_k->x; y_k[1] = u_k->y; y_k[2] = u_k->z;

    Info->Lgm_MagStep_BS_x[k] = x_k;
    Info->Lgm_MagStep_BS_u[k] = *u_k; // line not needed -- just debugginh...

    for (j=0; j<3; j++) y[j] = dy[j] = y_k[j];

    if (k == 0) {
        for (j=0; j<3; j++) Info->Lgm_MagStep_BS_d[j][0] = y_k[j];
    } else {
        for (j=0; j<3; j++) c[j] = y_k[j];

        for (k1=1; k1<k; k1++) {
            delta = 1.0/(Info->Lgm_MagStep_BS_x[k-k1] - x_k);
//printf("delta = %g\n", delta);
            f1 = x_k*delta;
            f2 = Info->Lgm_MagStep_BS_x[k-k1]*delta;
            for ( j=0; j<3; j++){
                q = Info->Lgm_MagStep_BS_d[j][k1];
                Info->Lgm_MagStep_BS_d[j][k1] =  dy[j];
                delta = c[j]-q;
                dy[j] = f1*delta;
                c[j] = f2*delta;
                y[j] += dy[j];
            }
        }
        for (j=0; j<3; j++) Info->Lgm_MagStep_BS_d[j][k] = dy[j];
    }


//printf("Lgm_PolFunExt:\n");
//for (j=0; j<=k; j++){
//printf("k = %d  x,y = %.15lf    %.15lf %.15lf %.15lf\n", k, Info->Lgm_MagStep_BS_x[j], Info->Lgm_MagStep_BS_u[j].x, Info->Lgm_MagStep_BS_u[j].y, Info->Lgm_MagStep_BS_u[j].z);
//}

    /*
     *  (x_k, u_k) are the input "data points". The polynomial function should go through
     *  all of these points. w is the polynomial function evaluated at x_k = 0.0 -- i.e.
     *  an estimate of what the sequence of modified midpoint estimates would converge to
     *  for infinitely small step sizes. dw is the error in this estimate.
     *
     */
    w->x  = y[0];    w->y  = y[1];    w->z  = y[2];
//printf("k = %d  w = %.15lf %.15lf %.15lf\n", w->x, w->y, w->z );
    dw->x = dy[0];   dw->y = dy[1];   dw->z = dy[2];
}


void polyextr( int k, double table[][3], double *last, Lgm_MagModelInfo *Info ) {

    int i, j;

    for ( j=k-1; j>0; j-- ) {       // Update the current row using the Neville recursive formula.
        for ( i=0; i<3; i++ ) {
            table[j-1][i] = table[j][i] + Info->Lgm_MagStep_BS_coeff[k][j]*(table[j][i]-table[j-1][i]);
        }
    }
    for ( i=0; i<3; i++) {       // Update the last element.
        last[i] = table[0][i] + Info->Lgm_MagStep_BS_coeff[k][0]*(table[0][i]-last[i]);
    }

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

int Lgm_MagStep_BS( Lgm_Vector *u, Lgm_Vector *u_scale,
          double htry, double *hdid, double *hnext,
          double eps, double sgn, double *s, int *reset,
          int (*Mag)(Lgm_Vector *, Lgm_Vector *, Lgm_MagModelInfo *), Lgm_MagModelInfo *Info ){


    Lgm_Vector  u0, b0, v;
    int         k;
    int         ModMidSuccessfull;
    double      err, Bmag;
    double      hopt[LGM_MAGSTEP_IMAX], work[LGM_MAGSTEP_IMAX];
    double      logfact, ratio, expo, facmin, fac, g;
    double      y[3], ysav[3], yseq[3], scale[3];
    double      rtol, atol;
    double      h, hnew, table[LGM_MAGSTEP_KMAX][3];
    int         i, l, forward;
    int         firstk, kopt;
    /*
     *  For years, we ran with this seq;
     *
     *      static int     Seq[] = { 0, 1, 2, 4, 6, 8, 12, 18, 24, 32, 48, 64 };
     *
     *  but, some tests show this is more efficient;
     *
     *      static int     Seq[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 18, 24, 32, 48, 64, 128, 256 };
     *
     *  Ratio of required func calls (for the test I did) was about 1.6
     */
    static double   STEPFAC1=0.65, STEPFAC2=0.94, STEPFAC3=0.02, STEPFAC4=4.0, KFAC1=0.8, KFAC2=0.9;
    //static int      Seq[] = { 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38 };
    static int      Seq[] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 18, 24, 32, 48, 64, 128, 256 };




//printf("u = %g %g %g\n", u->x, u->y, u->z);

    scale[0] = u_scale->x; scale[1] = u_scale->y; scale[2] = u_scale->z;

    if ( fabs(htry) < 1e-16 ) {
        printf("%s(): Requested stepsize is very small (returning with -1). htry = %g\n", __func__, htry );
        return(-1);
    }

    atol = Info->Lgm_MagStep_BS_atol;
    rtol = Info->Lgm_MagStep_BS_rtol;


    if ( ( eps != Info->Lgm_MagStep_BS_eps_old ) || ( *reset ) ) {

        Info->Lgm_nMagEvals = 0;
        Info->Lgm_MagStep_BS_eps_old = eps;
        *s = 0.0;

        Info->Lgm_MagStep_BS_cost[0] = Seq[0] + 1;
        for ( k=0; k<LGM_MAGSTEP_KMAX; k++ ) Info->Lgm_MagStep_BS_cost[k+1] = Info->Lgm_MagStep_BS_cost[k] + Seq[k+1];
        logfact = -0.6*log10( FMAX(1.0e-12,rtol) ) + 0.5;
        *hnext  = -1e99;
        Info->Lgm_MagStep_BS_k_targ  = FMAX(1,FMIN(LGM_MAGSTEP_KMAX-1,(int)logfact));

        for ( k = 0; k<LGM_MAGSTEP_IMAX; k++ ) {
            for ( l=0; l<k; l++ ) {
                ratio       = (double)Seq[k]/(double)Seq[l];
                Info->Lgm_MagStep_BS_coeff[k][l] = 1.0/(ratio*ratio-1.0);
            }
        }

    }




    /*
     *  Evaluate at u0
     */
    u0 = *u;
    if ( (*Mag)(&u0, &b0, Info) == 0 ) {
        // bail if B-field eval had issues.
        printf("%s(): B-field evaluation at u0 = %g %g %g returned with errors (returning with -1)\n", __func__, u0.x, u0.y, u0.z );
        return(-1);
    }
    ++(Info->Lgm_nMagEvals);
    Bmag = Lgm_NormalizeVector(&b0);
    if ( Bmag < 1e-16 ) {
        // bail if B-field magnitude is too small
        printf("%s(): Bmag too small at u0 = %g %g %g (Bmag = %g; Bx,y,z = %g, %g, %g) (returning with -1).\n", __func__, u0.x, u0.y, u0.z, Bmag, b0.x, b0.y, b0.z );
//printf("Line %d\n", __LINE__ );
//exit(-1);
        return(-1);
    }


    // convert vector into array
    y[0] = u0.x;
    y[1] = u0.y;
    y[2] = u0.z;



    work[0] = 0;
    h = htry;
    forward = ( h>0 ) ? TRUE : FALSE;
    ysav[0] = u0.x; ysav[1] = u0.y; ysav[2] = u0.z;                 // Save the starting values.

    if ( (h != *hnext) && !(Info->Lgm_MagStep_BS_first_step) ) {    // h gets reset in Odeint for the last step. ????
        Info->Lgm_MagStep_BS_last_step = TRUE;
    }
    if ( Info->Lgm_MagStep_BS_reject ) {                            // Previous step was rejected.
        Info->Lgm_MagStep_BS_prev_reject = TRUE;
        Info->Lgm_MagStep_BS_last_step   = FALSE;
    }
    Info->Lgm_MagStep_BS_reject = FALSE;
    firstk = TRUE;
    hnew   = fabs(h);

Info->Lgm_MagStep_BS_last_step = FALSE;

    while ( firstk || Info->Lgm_MagStep_BS_reject ) {                        // Loop until step accepted.
        h = forward ? hnew : -hnew;
        firstk = FALSE;
        Info->Lgm_MagStep_BS_reject = FALSE;
        //if ( fabs(h) <= fabs(*s)*eps ) printf("step size underflow in StepperBS");

        for ( k=0; k<=Info->Lgm_MagStep_BS_k_targ+1; k++ ) {                 // Evaluate the sequence of modified midpoint integrations.

            /*
             * Do a modified midpoint call
             */
            //dy( ysav, h, k, yseq, ipt, derivs );
            ModMidSuccessfull = Lgm_ModMid( &u0, &b0, &v, h, Seq[k], sgn, Mag, Info );
            if ( !ModMidSuccessfull ) return(-1); // bail if Lgm_ModMid() had issues.
            yseq[0] = v.x; yseq[1] = v.y; yseq[2] = v.z;

            if (k == 0) {
                 y[0] = yseq[0];
                 y[1] = yseq[1];
                 y[2] = yseq[2];
            } else {
                // Store result in tableau.
                for (i=0;i<3;i++) table[k-1][i] = yseq[i];
            }

            if ( k != 0 ) {

                // Perform extrapolation.
                polyextr( k, &table[0], y, Info );

                // Compute normalized error estimate err_k.
                err = 0.0;
                for ( i=0; i<3; i++ ) {
                    scale[i]  = atol + rtol*FMAX(fabs(ysav[i]),fabs(y[i]));
                    g         = (y[i]-table[0][i])/scale[i];
                    err      += g*g;
                }
                err    = sqrt( err/3.0 );

                // Compute optimal stepsize for this order.
                expo   = 1.0/(double)(2*k+1);
                facmin = pow( STEPFAC3, expo );
                if (err == 0.0) {
                    fac = 1.0/facmin;
                } else {
                    fac = STEPFAC2/pow(err/STEPFAC1,expo);
                    fac = FMAX(facmin/STEPFAC4,FMIN(1.0/facmin,fac));
                }
                hopt[k] = fabs(h*fac);
                work[k] = Info->Lgm_MagStep_BS_cost[k]/hopt[k];  // Work per unit step

                if ( (Info->Lgm_MagStep_BS_first_step || Info->Lgm_MagStep_BS_last_step) && err <= 1.0 ) break;

                if ( (k == Info->Lgm_MagStep_BS_k_targ-1) && !Info->Lgm_MagStep_BS_prev_reject && !Info->Lgm_MagStep_BS_first_step && !Info->Lgm_MagStep_BS_last_step) {
                    g = Seq[Info->Lgm_MagStep_BS_k_targ]*Seq[Info->Lgm_MagStep_BS_k_targ+1]/(double)(Seq[0]*Seq[0]);
                    
                    if (err <= 1.0) {
                        // Converged within order window.
                        break;
                    } else if ( err>g*g ) {
                        // Criterion (17.3.17) predicts step will fail.
                        Info->Lgm_MagStep_BS_reject = TRUE;
                        Info->Lgm_MagStep_BS_k_targ = k;
                        if ( (Info->Lgm_MagStep_BS_k_targ>1) && (work[k-1]<KFAC1*work[k]) ) Info->Lgm_MagStep_BS_k_targ--;
                        hnew = hopt[Info->Lgm_MagStep_BS_k_targ];
                        break;
                    }
                }
                if ( k == Info->Lgm_MagStep_BS_k_targ ) {
                    g = (double)Seq[k+1]/(double)Seq[0];
                    if ( err <= 1.0 ) {
                        // Converged within order window.
                        break;
                    } else if ( err>g*g ) {
                        // Criterion (17.3.20) predicts step will fail.
                        Info->Lgm_MagStep_BS_reject = TRUE;
                        if ( (Info->Lgm_MagStep_BS_k_targ>1) && (work[k-1]<KFAC1*work[k]) ) Info->Lgm_MagStep_BS_k_targ--;
                        hnew = hopt[Info->Lgm_MagStep_BS_k_targ];
                        break;
                    }
                }
                if ( k == Info->Lgm_MagStep_BS_k_targ+1 ) {
                    if ( err > 1.0 ) {
                        Info->Lgm_MagStep_BS_reject = TRUE;
                        if ( (Info->Lgm_MagStep_BS_k_targ>1) && (work[Info->Lgm_MagStep_BS_k_targ-1]<KFAC1*work[Info->Lgm_MagStep_BS_k_targ]) ) Info->Lgm_MagStep_BS_k_targ--;
                        hnew = hopt[Info->Lgm_MagStep_BS_k_targ];
                    }
                    break;
                }
            }
        } // Go back and try next k.

        // Arrive here from any break in for loop.
        if ( Info->Lgm_MagStep_BS_reject ) Info->Lgm_MagStep_BS_prev_reject = TRUE;

    }   // Go back if step was rejected.


    u->x = y[0]; u->y = y[1]; u->z = y[2];

    Info->Lgm_MagStep_BS_snew = *s + h;
    *s     = Info->Lgm_MagStep_BS_snew;
    *hdid  = h;

    Info->Lgm_MagStep_BS_first_step = FALSE;


    


    // Determine optimal order for next step.
    if ( k == 1 ) {
        kopt = 2;
    } else if (k <= Info->Lgm_MagStep_BS_k_targ) {
        kopt = k;
        if ( work[k-1] < KFAC1*work[k] ) {
            kopt = k-1;
        } else if ( work[k] < KFAC2*work[k-1] ) {
            kopt = FMIN(k+1,LGM_MAGSTEP_KMAX-1);
        }
    } else {
        kopt = k-1;
        if ( (k > 2) && (work[k-2] < KFAC1*work[k-1]) ) kopt = k-2;
        if (work[k] < KFAC2*work[kopt]) kopt = FMIN(k,LGM_MAGSTEP_KMAX-1);
    }

    if ( Info->Lgm_MagStep_BS_prev_reject ) {
        // After a rejected step neither order nor stepsize should increase.
        Info->Lgm_MagStep_BS_k_targ = FMIN(kopt,k);
        hnew   = FMIN(fabs(h),hopt[Info->Lgm_MagStep_BS_k_targ]);
        Info->Lgm_MagStep_BS_prev_reject = FALSE;
    } else {
        // Stepsize control for next step.
        if (kopt <= k) {
            hnew=hopt[kopt];
        } else {
            if ( (k<Info->Lgm_MagStep_BS_k_targ) && (work[k]<KFAC2*work[k-1]) ){
                hnew = hopt[k]*Info->Lgm_MagStep_BS_cost[kopt+1]/Info->Lgm_MagStep_BS_cost[k];
            } else {
                hnew = hopt[k]*Info->Lgm_MagStep_BS_cost[kopt]/Info->Lgm_MagStep_BS_cost[k];
            }
        }
        Info->Lgm_MagStep_BS_k_targ = kopt;
    }
    *hnext = ( forward ) ? hnew : -hnew;


    *reset = FALSE;
    
    return(1);

}


/*
 *
 *
 *  This routine is based on rkqs in in Num. Rec. 5th order Runge-Kutta.
 *
 *
 */

int Lgm_MagStep_RK5( Lgm_Vector *u, Lgm_Vector *u_scale,
          double Htry, double *Hdid, double *Hnext,
          double eps, double sgn, double *s, int *reset,
          int (*Mag)(Lgm_Vector *, Lgm_Vector *, Lgm_MagModelInfo *), Lgm_MagModelInfo *Info ){

    int         Done, Count, i;
    double      h, ErrMax, htemp, yscal[3], yerr[3], snew, Bmag;
    Lgm_Vector  u0, v, b0;


    if (  *reset  ) {
        Info->Lgm_nMagEvals = 0;
        Info->Lgm_MagStep_RK5_snew = 0.0;
        *s   = 0.0;
        Info->Lgm_MagStep_RK5_FirstTimeThrough = TRUE;

    }

    Count = 0;
    Done  = FALSE;


    /*
     * Htry can be positive or negative. The sign is independant of the sgn
     * variable which can also be +/-.
     *
     */
    h  = Htry;
    u0 = *u;
//printf("RK5: u = %g %g %g;   u_scale = %g %g %g\n", u->x, u->y, u->z, u_scale->x, u_scale->y, u_scale->z);
    yerr[0] = yerr[1] = yerr[2] = 0.0;
    yscal[0] = u_scale->x; yscal[1] = u_scale->y; yscal[2] = u_scale->z;
//yscal[0] = 100.0;
//yscal[1] = 100.0;
//yscal[2] = 100.0;

    // get b0
    if ( (*Mag)(&u0, &b0, Info) == 0 ) {
        // bail if B-field eval had issues.
        printf("Lgm_RK5(): B-field evaluation (u = %g %g %g) returned with errors (returning with 0)\n", u0.x, u0.y, u0.z );
        return(0);
    }
    ++(Info->Lgm_nMagEvals);
    Bmag = Lgm_NormalizeVector(&b0);
    if ( Bmag < 1e-16 ) {
        // bail if B-field magnitude is too small
        printf("Lgm_RK5(): Bmag too small (u = %g %g %g Bmag = %g) is too small (returning with 0).\n", u0.x, u0.y, u0.z, Bmag );
        return(0);
    }



    while ( !Done && ( Count < Info->Lgm_MagStep_RK5_MaxCount ) ) {

        /*
         * Take a Cash-Karp Runge-Kutta step (take point from u0->v )
         */
//printf("RK5 C-K: h = %g\n", h);
        Lgm_RKCK( &u0, &b0, &v, h, sgn, yerr, Mag, Info );


        /*
         * Test accuracy.
         */
double dum = 0.0;
        ErrMax = 0.0;
        for ( i=0; i<3; i++ ) {
//            ErrMax = FMAX( ErrMax, fabs( yerr[i]/yscal[i]) );
dum    = fabs( yerr[i]/yscal[i]);
ErrMax += dum*dum;
        }
ErrMax = sqrt(ErrMax); //gets to RMS of scaled error

        ErrMax /= eps;

        if ( ErrMax > 1.0 ) {

            htemp = Info->Lgm_MagStep_RK5_Safety*h*pow( ErrMax, Info->Lgm_MagStep_RK5_pShrnk );
            h = ( (h>=0.0) ? FMAX(htemp, 0.1*h) : FMIN(htemp, 0.1*h) );
            //Info->Lgm_MagStep_RK5_snew = *s + h;
            snew= *s + h;
//printf("Count = %d eps = %g yerr = %g %g %g   yscal = %g %g %g Errmax = %g   snew = %g    *s = %g  htemp = %g    h = %g\n", Count, eps, yerr[0], yerr[1], yerr[2], yscal[0], yscal[1], yscal[2], ErrMax, snew, *s, htemp, h );
            if ( snew == *s ) {
                printf( "Lgm_MagStep2: Stepsize underflow in rkqs. h = %g\n", h );
                return( -1 );
            }

        } else {

            if ( ErrMax > Info->Lgm_MagStep_RK5_ErrCon ) {
                *Hnext = Info->Lgm_MagStep_RK5_Safety*h*pow( ErrMax, Info->Lgm_MagStep_RK5_pGrow );
            } else {
                *Hnext = 5.0*h;
            }
            *Hdid  = h;
            *s    += h;
            *u     = v;
//printf("Hdid = %g Hnext = %g START, FINAL, |DIFF| = %g %g %g   %g %g %g    %g\n", *Hdid, *Hnext, u0.x, u0.y, u0.z, u->x, u->y, u->z, Lgm_VecDiffMag( u, &u0 ) );
            Done   = TRUE;

        }

        ++Count;

    }

    return(1);

}

int Lgm_RKCK( Lgm_Vector *u0, Lgm_Vector *b0, Lgm_Vector *v, double h, double sgn, double *yerr,
            int (*Mag)(Lgm_Vector *, Lgm_Vector *, Lgm_MagModelInfo *), Lgm_MagModelInfo *Info ) {


    int i;
    double  a2=0.2, a3=0.3, a4=0.6, a5=1.0, a6=0.875;

    double  b21 = 0.2;

    // b31=3.0/40.0, b32=9.0/40.0;
    double  b31 = 0.075;
    double  b32 = 0.225;

    // b41=3.0/10.0, b42=-9.0/10.0, b43=6.0/5.0;
    double  b41 =  0.3;
    double  b42 = -0.9;
    double  b43 =  1.2;

    //double  b51=-11.0/54.0, b52=2.5, b53=-70.0/27.0, b54=35.0/27.0;
    double  b51 = -0.20370370370370370370;
    double  b52 =  2.5;
    double  b53 = -2.59259259259259259259;
    double  b54 =  1.29629629629629629629;

    //double  b61=1631.0/55296.0, b62=175.0/512.0, b63=575.0/13824.0, b64=44275.0/110592.0, b65=253.0/4096.0;
    double  b61 = .02949580439814814814;
    double  b62 = .341796875;
    double  b63 = .04159432870370370370;
    double  b64 = .40034541377314814814;
    double  b65 = .061767578125;

    //double  c1=37.0/378.0, c3=250.0/621.0, c4=125.0/594.0, c6=512.0/1771.0;
    double  c1  = 0.09788359788359788359;
    double  c3  = 0.40257648953301127214;
    double  c4  = 0.21043771043771043771;
    double  c6  = 0.28910220214568040654;

    //double  dc1 = c1-2825.0/27648.0, dc3 = c3-18575.0/48384.0, dc4 = c4-13525.0/55296.0, dc5 = -277.00/14336.0, dc6 = c6-0.25;
    double  dc1 = -0.00429377480158730159;
    double  dc3 =  0.01866858609385783299; 
    double  dc4 = -0.03415502683080808080;
    double  dc5 = -0.01932198660714285714;
    double  dc6 =  0.03910220214568040654;

    double  ak1[3], ak2[3], ak3[3], ak4[3], ak5[3], ak6[3], y[3], yout[3], ytemp[3];

    Lgm_Vector  u, B;
    double      Bmag, H;


    y[0]    = u0->x;
    y[1]    = u0->y;
    y[2]    = u0->z;


    H = sgn*h;

//printf("    In cash-karp step: H = %g\n", H);

    // 1st step
    ak1[0] = b0->x; ak1[1] = b0->y; ak1[2] = b0->z;
    for ( i=0; i<3; i++ ) {
        ytemp[i] = y[i] + b21*H*ak1[i]; //WHY DOESN'T ytemp GET USED IN NEXT CALL TO Mag?????
    }



    // 2nd step
    u.x = u0->x + a2*H;
    u.y = u0->y + a2*H;
    u.z = u0->z + a2*H;
    if ( (*Mag)(&u, &B, Info) == 0 ) {
        // bail if B-field eval had issues.
        printf("Lgm_RKCK(): B-field evaluation during cash-karp step (u = %g %g %g) returned with errors (returning with 0)\n", u.x, u.y, u.z );
        return(0);
    }
    ++(Info->Lgm_nMagEvals);
    Bmag = Lgm_NormalizeVector(&B);
    if ( Bmag < 1e-16 ) {
        // bail if B-field magnitude is too small
        printf("Lgm_RKCK(): Bmag too small during cash-karp phase (u = %g %g %g Bmag = %g) is too small (returning with 0).\n", u.x, u.y, u.z, Bmag );
        return(0);
    }
    ak2[0] = B.x; ak2[1] = B.y; ak2[2] = B.z;
    for ( i=0; i<3; i++ ) {
        ytemp[i] = y[i] + H*(b31*ak1[i] + b32*ak2[i]);
    }




    // 3rd step
    u.x = u0->x + a3*H;
    u.y = u0->y + a3*H;
    u.z = u0->z + a3*H;
    if ( (*Mag)(&u, &B, Info) == 0 ) {
        // bail if B-field eval had issues.
        printf("Lgm_RKCK(): B-field evaluation during cash-karp step (u = %g %g %g) returned with errors (returning with 0)\n", u.x, u.y, u.z );
        return(0);
    }
    ++(Info->Lgm_nMagEvals);
    Bmag = Lgm_NormalizeVector(&B);
    if ( Bmag < 1e-16 ) {
        // bail if B-field magnitude is too small
        printf("Lgm_RKCK(): Bmag too small during cash-karp phase (u = %g %g %g Bmag = %g) is too small (returning with 0).\n", u.x, u.y, u.z, Bmag );
        return(0);
    }
    ak3[0] = B.x; ak3[1] = B.y; ak3[2] = B.z;
    for ( i=0; i<3; i++ ) {
        ytemp[i] = y[i] + H*(b41*ak1[i] + b42*ak2[i] + b43*ak3[i]);
    }




    // 4th step
    u.x = u0->x + a4*H;
    u.y = u0->y + a4*H;
    u.z = u0->z + a4*H;
    if ( (*Mag)(&u, &B, Info) == 0 ) {
        // bail if B-field eval had issues.
        printf("Lgm_RKCK(): B-field evaluation during cash-karp step (u = %g %g %g) returned with errors (returning with 0)\n", u.x, u.y, u.z );
        return(0);
    }
    ++(Info->Lgm_nMagEvals);
    Bmag = Lgm_NormalizeVector(&B);
    if ( Bmag < 1e-16 ) {
        // bail if B-field magnitude is too small
        printf("Lgm_RKCK(): Bmag too small during cash-karp phase (u = %g %g %g Bmag = %g) is too small (returning with 0).\n", u.x, u.y, u.z, Bmag );
        return(0);
    }
    ak4[0] = B.x; ak4[1] = B.y; ak4[2] = B.z;
    for ( i=0; i<3; i++ ) {
        ytemp[i] = y[i] + H*(b51*ak1[i] + b52*ak2[i] + b53*ak3[i] + b54*ak4[i]);
    }





    // 5th step
    u.x = u0->x + a5*H;
    u.y = u0->y + a5*H;
    u.z = u0->z + a5*H;
    if ( (*Mag)(&u, &B, Info) == 0 ) {
        // bail if B-field eval had issues.
        printf("Lgm_RKCK(): B-field evaluation during cash-karp step (u = %g %g %g) returned with errors (returning with 0)\n", u.x, u.y, u.z );
        return(0);
    }
    ++(Info->Lgm_nMagEvals);
    Bmag = Lgm_NormalizeVector(&B);
    if ( Bmag < 1e-16 ) {
        // bail if B-field magnitude is too small
        printf("Lgm_RKCK(): Bmag too small during cash-karp phase (u = %g %g %g Bmag = %g) is too small (returning with 0).\n", u.x, u.y, u.z, Bmag );
        return(0);
    }
    ak5[0] = B.x; ak5[1] = B.y; ak5[2] = B.z;
    for ( i=0; i<3; i++ ) {
        ytemp[i] = y[i] + H*(b61*ak1[i] + b62*ak2[i] + b63*ak3[i] + b64*ak4[i] + b65*ak5[i]);
    }





    // 6th step
    u.x = u0->x + a6*H;
    u.y = u0->y + a6*H;
    u.z = u0->z + a6*H;
    if ( (*Mag)(&u, &B, Info) == 0 ) {
        // bail if B-field eval had issues.
        printf("Lgm_RKCK(): B-field evaluation during cash-karp step (u = %g %g %g) returned with errors (returning with 0)\n", u.x, u.y, u.z );
        return(0);
    }
    ++(Info->Lgm_nMagEvals);
    Bmag = Lgm_NormalizeVector(&B);
    if ( Bmag < 1e-16 ) {
        // bail if B-field magnitude is too small
        printf("Lgm_RKCK(): Bmag too small during cash-karp phase (u = %g %g %g Bmag = %g) is too small (returning with 0).\n", u.x, u.y, u.z, Bmag );
        return(0);
    }
    ak6[0] = B.x; ak6[1] = B.y; ak6[2] = B.z;
    for ( i=0; i<3; i++ ) {
        yout[i] = y[i] + H*(c1*ak1[i] + c3*ak3[i] + c4*ak4[i] + c6*ak6[i]);
    }
    v->x = yout[0];
    v->y = yout[1];
    v->z = yout[2];





    // compute yerr
    for ( i=0; i<3; i++ ) {
        yerr[i] = H*(dc1*ak1[i] + dc3*ak3[i] + dc4*ak4[i] + dc5*ak5[i] + dc6*ak6[i]);
    }



    return(1);

}

