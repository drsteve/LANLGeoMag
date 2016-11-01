#include "Lgm/Lgm_TabularBessel.h"



void Lgm_TabularBessel_Init( int M, int N, double xmin, double xmax, Lgm_TabularBessel *tb ){

    int    n, i;
    double x;
    

    LGM_ARRAY_1D( tb->JnTabular_x,       M, double );
    LGM_ARRAY_2D( tb->JnTabular,    N+3, M, double );
    LGM_ARRAY_2D( tb->JDnTabular,   N+3, M, double );
    LGM_ARRAY_2D( tb->JDDnTabular,  N+3, M, double );

    tb->M    = M;       // Number of points in table
    tb->N    = N;       // Max order for Bessel Functions
    tb->xmin = xmin;
    tb->xmax = xmax;
    tb->s    = (xmax-xmin)/(double)M;

    /*
     * Compute Jn(x) n: [0, N+2]
     */
    for ( n=0; n<=(N+2); n++ ) {
        for ( i=0; i<M; i++ ) {
            x = tb->s*i + tb->xmin;
            tb->JnTabular_x[i]   = x;
            tb->JnTabular[n][i]  = gsl_sf_bessel_Jn( n, x );
        }
    }

    /*
     * Compute J'n(x) n: [0, N+1]
     */
    for ( i=0; i<M; i++ ) {
        x = tb->JnTabular_x[i];
        tb->JDnTabular[0][i] = -tb->JnTabular[1][i];
        for ( n=1; n<=(N+1); n++ ) {
            tb->JDnTabular[n][i] = 0.5*(tb->JnTabular[n-1][i] - tb->JnTabular[n+1][i]);
        }
    }

    /*
     * Compute J''n(x) n: [0, N]
     */
    for ( i=0; i<M; i++ ) {
        x = tb->JnTabular_x[i];
        tb->JDDnTabular[0][i] = -tb->JDnTabular[1][i];
        for ( n=1; n<=N; n++ ) {
            tb->JDDnTabular[n][i] = 0.5*(tb->JDnTabular[n-1][i] - tb->JDnTabular[n+1][i]);
        }
    }

}

void Lgm_TabularBessel_Free( Lgm_TabularBessel *tb ){

    LGM_ARRAY_1D_FREE( tb->JnTabular_x );
    LGM_ARRAY_2D_FREE( tb->JnTabular );
    LGM_ARRAY_2D_FREE( tb->JDnTabular );
    LGM_ARRAY_2D_FREE( tb->JDDnTabular );

    free( tb );
}


/*
 * This uses simple linear interp. Need to use large number of points in table.
 */
void Lgm_TabularBessel_Eval( double x, int Q, double *Jn, double *JDn, Lgm_TabularBessel *tb ) {

    double  x0, x1, y0, y1, s;
    int     n, ii; 

    ii = (x-tb->xmin)/tb->s + 1;

    x0 = tb->JnTabular_x[ii];
    x1 = tb->JnTabular_x[ii+1];

    for ( n=0; n<=14; n++ ) {

        y0 = tb->JnTabular[n][ii];
        y1 = tb->JnTabular[n][ii+1];
        s  = (y1-y0)/(x1-x0);

        Jn[n]  = s*(x-x0) + y0;
        //JDn[n] =  etc..

    }

}





/*
 * This uses hermite interp. May be slower than using bess function in the first place !
 */
void Lgm_TabularBessel_Eval2( double x, int Q, double *Jn, double *JDn, Lgm_TabularBessel *tb ) {

    double  za1[12], dd1[12];
    double  za2[12], dd2[12];
    double  xa[6], ya[6], yda[6];
    int     n, i, ii, ipq2, iipi;

    ii = (x-tb->xmin)/tb->s + 1;

    for (i=-Q/2; i<Q/2; i++){
        ipq2 = i+Q/2;
        iipi = ii+i;
        xa[ipq2]  = tb->JnTabular_x[iipi];
    }

    for ( n=0; n<=tb->N; n++ ) {

        for (i=-Q/2; i<Q/2; i++){

            ipq2 = i+Q/2;
            iipi = ii+i;

            ya[ipq2]  = tb->JnTabular[n][iipi];
            yda[ipq2] = tb->JDnTabular[n][iipi];
        }
        gsl_poly_dd_hermite_init( dd1, za1,    xa, ya, yda, Q );

        for (i=-Q/2; i<Q/2; i++){

            ipq2 = i+Q/2;
            iipi = ii+i;

            ya[ipq2]  = tb->JDnTabular[n][iipi];
            yda[ipq2] = tb->JDDnTabular[n][iipi];

        }
        gsl_poly_dd_hermite_init( dd2, za2,    xa, ya, yda, Q );

        Jn[n]  = gsl_poly_dd_eval( dd1, za1, 2*Q, x );
        JDn[n] = gsl_poly_dd_eval( dd2, za2, 2*Q, x );

    }



}







