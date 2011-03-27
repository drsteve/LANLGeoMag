#include <stdio.h>
#include <math.h>
#include <Lgm_MagModelInfo.h>
#include <Lgm_FluxToPsd.h>
int main( ) {

    Lgm_CTrans     *c = Lgm_init_ctrans( 0 );
    Lgm_DateTime   *d;
    Lgm_Vector      u;
    Lgm_FluxToPsd  *f2p = Lgm_F2P_CreateFluxToPsd(1);

    double          **J, *E, *A;
    int             nE, nA, i, j;

    double          *Mu, *K;
    int             nMu, nK;

    double          f, p2c2, df, sa;


    /*
     * Set date/time and position
     */
    d = Lgm_DateTime_Create( 2001, 1, 1, 0.0, LGM_TIME_SYS_UTC, c );
    u.x = -6.6; u.y =  0.0; u.z =  0.0;
    Lgm_F2P_SetDateTimeAndPos( d, &u, f2p );
    Lgm_DateTime_Destroy( d );


    /*
     *  Make up some fake flux data... Lets take a two-component relativistic
     *  maxwellian as described by Cayton et al. [1989].
     *
     *  Maxwellian 1; n= 1e-4 cm^-3; T = 200 keV
     *  Maxwellian 2; n= 5e-3 cm^-3; T =  25 keV
     *
     */
    nE = 100;
    nA = 90;
    LGM_ARRAY_1D( E, nE, double );
    LGM_ARRAY_1D( A, nA, double );
    LGM_ARRAY_2D( J, nE, nA, double );
    for (j=0; j<nA; j++ ) {A[j] = 1.0+(double)j; printf("A[%d] = %g\n", j, A[j]); }
    for (i=0; i<nE; i++ ) E[i] = 1.0+(double)i/20.0;
//    E[0] = 0.1; E[1] = 0.2; E[2] = 0.3; E[3] = 0.5; E[4] = 0.8; E[5] = 1.0; E[6] = 3.0; 
//    A[0] = 10.0; A[1] = 20.0; A[2] = 30.0; A[3] = 40.0; A[4] = 50.0; A[5] = 60.0; A[6] = 70.0; A[7] = 80.0; A[8] = 90.0;
    for ( i=0; i<nE; i++ ) {
        f    = Lgm_MaxJut( 5e-3, 25.0, E[i], LGM_Ee0 ) + Lgm_MaxJut( 1e-4, 200.0, E[i], LGM_Ee0 ); // PSD c^3/cm^3/MeV^3
//        f    = Lgm_MaxJut( 1e-4, 200.0, E[i], LGM_Ee0 ); // PSD c^3/cm^3/MeV^3
        p2c2 = Lgm_p2c2( E[i], LGM_Ee0 ); // p^2c^2 MeV
        df   = Lgm_PsdToDiffFlux( f, p2c2 ); // differential flux #/cm^3/sr/s/MeV
        for ( j=0; j<nA; j++ ) {
            sa = sin(A[j]*RadPerDeg); // sin(PitchAngle)
printf("sa*sa = %g\n", sa*sa);
            J[i][j] = pow(sa, 10.0)*df;       // sin^2() modulation
if (j==20) 
            J[i][j] = .0001*df;       // sin^2() modulation
        }
    }

    // Add diff. flux data/info to f2p structure.
    Lgm_F2P_SetFlux( J, E, nE, A, nA, f2p ); 



    /*
     *  Choose a set of Mu's and K's 
     *  and compute Phase Space Density at these constant Mu's and K's
     */
    nMu = 18;
    nK  = 18;
    LGM_ARRAY_1D( Mu, nMu, double );
    LGM_ARRAY_1D( K, nK, double );
    for (i=0; i<nMu; i++) Mu[i] = 100.0 + 100.0*i;
//    K[0] = 0.01; for (j=1; j<nK; j++) K[j] = K[j-1]*2.0;
for (j=0; j<nK; j++) K[j] = 1.0 + (double)j;
    Lgm_F2P_GetPsdAtConstMusAndKs( Mu, nMu, K, nK, f2p );
    
    

    LGM_ARRAY_1D_FREE( E );
    LGM_ARRAY_1D_FREE( A );
    LGM_ARRAY_2D_FREE( J );
    LGM_ARRAY_1D_FREE( Mu );
    LGM_ARRAY_1D_FREE( K );


    Lgm_free_ctrans( c );
    Lgm_F2P_FreeFluxToPsd( f2p );

    return(0);
}

