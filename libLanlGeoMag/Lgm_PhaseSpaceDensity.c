#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Lgm/Lgm_PhaseSpaceDensity.h"
#include "Lgm/Lgm_CTrans.h"
#include "Lgm/Lgm_DynamicMemory.h"


/*
 *  Routines for converting Flux <-> PSD.
 */




/*
 * Returns First adiabatic invariant, Mu, given: Particle's kinetic energy,
 * Pitch Angle and the local B-field strength. (Mu is Eperp/B.)
 *  
 *  Inputs:
 *          E -- Kinetic Energy of Particle.    ( MeV )
 *          a -- Pitch Angle of Particle.       ( Degrees )
 *          B -- Local magnetic field strength. ( nT )
 *
 *  Returns:
 *          First adiabatic invariant, Mu. ( MeV/nT )
 */
double  Lgm_Energy_to_Mu( double E, double a, double B ) {

    double  sa, sa2;

    if ( B <= 1e-12 ) return( 9e99 );

    sa = sin( a*RadPerDeg );    // sin(Alpha)
    sa2 = sa*sa;                // sin^2(Alpha)
    
    return( E*sa2/B );          // Mu = E*sin^2(Alpha)/B
    
}



/*
 * Returns Particle's Kinetic Energy, E, given: Particle's first invariant, Mu,
 * Pitch Angle and the local B-field strength.
 *  
 *  Inputs:
 *          Mu -- Kinetic Energy of Particle.    ( MeV/nT )
 *          a  -- Pitch Angle of Particle.       ( Degrees )
 *          B  -- Local magnetic field strength. ( nT )
 *
 *  Returns:
 *          First adiabatic invariant, Mu. ( MeV )
 */
double  Lgm_Mu_to_Energy( double Mu, double a, double B ) {

    double  sa, sa2;

    
    sa = sin( a*RadPerDeg );    // sin(Alpha)
    sa2 = sa*sa;                // sin^2(Alpha)

    if ( sa2 < 1e-12 ) return( 9e99 );
    else return( Mu*B/sa2 );    // Mu = E*sin^2(Alpha)/B


}


/*
 * Compute p^2c^2 given a particle's kinetic energy and rest energy.
 *
 *  Some relativistic equations:
 *  
 *  With,
 *      m     = Particle mass.
 *      m0    = Particle rest mass.
 *      v     = Particle speed.
 *      c     = Speed of Light.
 *      gamma = (1-v^2/c^2)^(-1/2).
 *
 *      m = m0 gamma = m0 (1-v^2/c^2)^(-1/2)
 *
 *  Rearrange this to get,
 *      (mc^2)^2 = (m0c^2)^ + p^2c^2
 *
 *  With,
 *      E  = mc^2  = Total Energy of particle
 *      E0 = m0c^2 = Rest Energy of particle
 *      p  = relativistic momentum of particle
 *
 *  this is,
 *      E^2 = E0^2 + (pc)^2
 *  
 *  so, 
 *      (pc)^2  = E^2 - E0^2
 *  
 *  Let E = Ek+E0 (kinetic energy + rest energy). Then,
 *      p^2c^2  = (Ek+E0)^2 - E0^2
 *              = Ek (Ek+2E0)
 *  
 *
 *  Inputs:
 *          E  -- Kinetic Energy of Particle.    ( MeV )
 *          E0 -- Rest energy of Particle.       ( MeV )
 *
 *  Returns:
 *          p^2c^2 = Ek(Ek+2E0) ( MeV^2 )
 */
double  Lgm_p2c2( double Ek, double E0 ) {
    return( Ek*(Ek+2.0*E0) );    // p^2c^2 in units of MeV^2
}


/*
 * returns (v/c)^2
 *
 *  (v/c)^2 = m^2v^2/(m^2c^2) 
 *          = m^2v^2c^2/(m^2c^4) 
 *          = p^2c^2/(m^2c^4) 
 *          = p^2c^2/E^2
 *          = Ek(Ek+2E0)/(Ek+E0)^2
 */
double  Lgm_v2overc2( double Ek, double E0 ) {
    double  E = Ek + E0;
    return( Ek*(Ek+2.0*E0)/(E*E) );    // dimensionless
}


/*
 * returns relativistic factor gamma = [ 1 - (v/c)^2 ]^(-1/2)
 */
double  Lgm_gamma( double Ek, double E0 ) {
    double  E = Ek + E0;
    return( 1.0/sqrt( 1.0 - Ek*(Ek+2.0*E0)/(E*E) ) );    // dimensionless
}


/*
 * Convert differential flux to phase space density.
 *
 *  The basic relationship is;
 *      f = j/p^2
 *
 *  Multiply top and bottom by c^2 gives,
 *      f = (j c^2)/(p^2c^2)
 *      f = j/c * c^3/(p^2c^2)
 *  
 *  Reason for making it c^3 is that the final units become,
 *  c^3 cm^-3 MeV^-3 or (c/cm/MeV)^3
 *
 * -------------------------------------------------------------
 *  Inputs:
 *          j    -- Differential Flux in units of #/cm^2/s/sr/MeV
 *          p2c2 -- (pc)^2 in units of Mev^2
 *  
 *  Output:
 *          Phase space density in units of (c/cm/MeV)^3
 *  
 */
double Lgm_DiffFluxToPsd( double j, double p2c2 ){
    return( j/(p2c2*2.9979e10) ); // f in units of (c/cm/MeV)^3
}




/*
 * Convert phase space density to differential flux.
 *
 *  The basic relationship is;
 *      f = j/p^2
 *
 *  Multiply top and bottom by c^2 gives,
 *      f = (j c^2)/(p^2c^2)
 *      f = j/c * c^3/(p^2c^2)
 *
 *  Reason for making it c^3 is that the final units become,
 *  c^3 cm^-3 MeV^-3 or (c/cm/MeV)^3
 *
 * -------------------------------------------------------------
 *  Inputs:
 *          f    -- Phase space density in units of (c/cm/MeV)^3
 *          p2c2 -- (pc)^2 in units of Mev^2
 *
 *  Output:
 *          Differential Flux in units of #/cm^2/s/sr/MeV
 *
 */
double Lgm_PsdToDiffFlux( double f, double p2c2 ){
    return( f*2.9979e10/p2c2 ); // j in units of #/cm^2/s/sr/MeV
}





/*
 * This allocates memory for and initializes a Lgm_PhaseSpaceDensity structure.
 * It returns a pointer to the alloced structure.
 *
 *
 *  Inputs:
 *                 Flux: 2D array containing the differential flux values as a function of energy and pitch angle.
 *                    E: 1D array containing the energy values implied by the first index of Flux[][] array.
 *                    A: 1D array containing the pitch angles values implied by the second index of Flux[][] array.
 *                   nE: number of energies.
 *                   nA: number of pitch angles.
 *      DumpDiagnostics: Flag to switch on/off diagnostic output.
 *
 */
Lgm_PhaseSpaceDensity *Lgm_InitPhaseSpaceDensity( double **J, double *E, double *A, int nE, int nA, int DumpDiagnostics ) {

    
    int     i, j;
    double  flux, p2c2, fp, Min, Max;
    double  **PsdArray_LoRes, **PsdArray_HiRes;

    Lgm_PhaseSpaceDensity *p;
    

    /*
     * Allocate memory for a Lgm_PhaseSpaceDensity structure.
     */
    p = (Lgm_PhaseSpaceDensity *) calloc( 1, sizeof(*p) );

    /*
     * Set DumpDiagnostics flag to what we got here. This can be changed later as well.
     */
    p->DumpDiagnostics = DumpDiagnostics;


    /*
     * Add Flux array info to p structure. Alloc arrays appropriately.
     */
    p->nE1 = nE; 
    p->nA1 = nA; 
    LGM_ARRAY_1D( p->E1, p->nE1, double );
    LGM_ARRAY_1D( p->A1, p->nA1, double );
    LGM_ARRAY_2D( p->FLUX_EA1, p->nE1, p->nA1, double );
    for (i=0; i<p->nE1; i++) p->E1[i] = E[i];
    for (i=0; i<p->nA1; i++) p->A1[i] = A[i];
    for (i=0; i<p->nE1; i++) {
        for (j=0; j<p->nA1; j++) {
            p->FLUX_EA1[i][j] = J[i][j];
        }
    }


    /*
     * Alloc mem for the PSD array.
     * Convert Flux array into to PSD array. Values are stored as log10(f).
     */
    LGM_ARRAY_2D( p->PSD_EA1, p->nE1, p->nA1, double );
    for (j=0; j<p->nA1; j++) {
        for (i=0; i<p->nE1; i++) {

            flux   = p->FLUX_EA1[i][j];
            p2c2   = Lgm_p2c2( p->E1[i], LGM_Ee0 );
            fp     = Lgm_DiffFluxToPsd( flux, p2c2 );
            fp     = (fp > 0.0) ? log10(fp) : LGM_FILL_VALUE;
            p->PSD_EA1[i][j] = fp;
            
        }
    }
    if ( p->DumpDiagnostics ) {
//        DumpGif( "PSD_Versus_E_and_A_LoRes.gif", p->nE1, p->nA1, p->PSD_EA1 );
    }




    /*
     * Set desired size of HiRes array.
     */
    p->nE2 = 400; // Energy 
    p->nA2 = 400; // Pitch Angles
    LGM_ARRAY_2D( p->PSD_EA2, p->nE2, p->nA2, double );         // (Energy, Pitch Angle) -> (Row/Col)
    UpSizeImage( p->PSD_EA1, p->E1, p->A1, p->nE1, p->nA1, 
                    p->PSD_EA2, p->E2, p->A2, p->nE2, p->nA2 ); // returns p->PSD_EA2, p->E2, p->A2
    if ( p->DumpDiagnostics ) {
//        DumpGif( "PSD_Versus_E_and_A_HiRes.gif", p->nE2, p->nA2, p->PSD_EA2 );
    }
   



    return p;

}


void Lgm_FreePhaseSpaceDensity( Lgm_PhaseSpaceDensity *p ) {

    LGM_ARRAY_1D_FREE( p->E1 );
    LGM_ARRAY_1D_FREE( p->A1 );
    LGM_ARRAY_2D_FREE( p->FLUX_EA1 );
    LGM_ARRAY_2D_FREE( p->PSD_EA1 );
    LGM_ARRAY_2D_FREE( p->PSD_EA2 );

    free( p );

    return;
}




/*
 *  The routine Lgm_ComputePsdVersusEAndAlpha() gives us a Hi-Res array of f(E,
 *  alpha).  BUT, what we really need in the end is f( mu, K ). Although mu is
 *  easy to compute, it is dependant on both E and alpha. K is only dependant
 *  upon alpha, but on the other hand K is not so easy to compute and we dont
 *  want to have to compute lots of K's if we dont have to. So we will use the
 *  following strategy instead,
 * 
 *  Note that f( E, a ) is the same as f( E(mu, a(K)), a(K) ). Thus, for a
 *  given mu and K, we can figure out what E and a they correspond to and then
 *  we can just look up the f value in our HiRes array. The steps are;
 *
 *      1. For each K, compute a(K). We already have this routine ( AlphaOfK() ).
 *
 *      2. Then we compute E from a and the given mu value.
 *
 *      3. Then just look up f(E,a) from the array (interp or fit or whatever).
 *
 *  Inputs:
 *          nMu -- Number of Mu values
 *           Mu -- 1-D array of Mu values
 *           nK -- Number of K values
 *            K -- 1-D array of K values
 *  Outputs:
 *          PSD -- 2-D array of f(Mu, K). The user must alloc memory for this array.
 * 
 *  Usage:
 *          If m->UseInterpRoutines is TRUE, then the user must have pre-traced
 *          a FL with Lgm_TraceLine(). This initializes field line dependant
 *          information that is needed for Lgm_AlphaOfK() to work in the
 *          interpolated mode. Basically Lgm_AlphaOfK() uses bisection to solve
 *          for a(K) and it much faster to use the pre-tracing strategy. E.g.;
 *
 *
 *              m->UseInterpRoutines = TRUE;
 *              Lgm_TraceLine( &u_in, &u_out, m->Lgm_LossConeHeight, -1.0, 1e-8, FALSE, m );
 *              LGM_ARRAY_2D( PSD, nMU, nK, double );
 *              Lgm_PsdAtConstMuAndK( PSD, nMu, MU, nK, K, m );
 *                  ... do stuff with PSD ...
 *              LGM_ARRAY_2D_FREE( PSD );
 * 
 * 
 */
//void Lgm_PsdAtConstMuAndK( double **PSD, int nMu, double *Mu, int nK, double *K, Lgm_MagModelInfo *m ) {
//
//    int     i, k;
//    double  *a, E;
//
//    
//    /*
//     * Compute the alpha's -- 1d array
//     * parallelize this loop?
//     */
//    // #pragma etc..
//    LGM_ARRAY_1D( a, nK, double );
//    for ( k=0; k<nK; k++ ){
////        a[k] = Lgm_AlphaOfK( K[k], m );
//    }
//
//
//    /*
//     * Compute the PSD's -- 2d array
//     */
//    for ( i=0; i<nMu; i++ ){
//        for ( k=0; k<nK; k++ ){
//            E = Lgm_Mu_to_Energy( Mu[i], a[k], B );
////            PSD[i][k] = Lgm_Psd( E, a, INFO THAT CONTAINS THE HIRES array );
//        }
//    }
//
//    LGM_ARRAY_1D_FREE( a );
//
//}






/*
 *  Reverse direction....
 *
 *  Note that f( mu, K ) is the same as f( mu( E, a), K(a) ). Thus, for a given
 *  E and a, we can figure out what mu and K they correspond to and then we
 *  can just look up the f value in our HiRes array. The steps would be;
 *
 *      1. For each a, compute K(a). We already have this routine somewhere.
 *
 *      2. Also compute mu from the given E and a.
 *
 *      3. Then just look f( mu, K ) up from the array (interp or fit or whatever).
 * 
 */












/*
 * Routine to increase the size of an image smoothly. Uses an area-weighting
 * interpolation scheme.
 */
void UpSizeImage( double **Orig, double *X, double *Y, int M, int N, double **New, double *U, double *V, int K, int J ) {

    double  DX, DY, DXO2, DYO2, A;
    double  x, y, xe, ye, xx, yy, xp, yp;
    double  f[4], w0, w1, w2, w3;
    double  dxm, dxp, dym, dyp;
    int     j, k, np, mp, ne, me;
    int     n0, n1, n2, n3;
    int     m0, m1, m2, m3;



    /*
     * compute size and area of a pixel in orig image.
     */
    DX = 1.0/(double)N;
    DY = 1.0/(double)M;
    DXO2 = 0.5*DX;
    DYO2 = 0.5*DY;
    A    = DX*DY;
    printf("A = %g\n", A);


    for (j=0; j<J; j++ ){
        for (k=0; k<K; k++ ){


            /*
             *  Compute the normalized coords (x,y) for the center of this
             *  pixel.
             */
            x = ((double)j+0.5)/(double)J;
            y = ((double)k+0.5)/(double)K;



            /*
             *  Determine which pixel of the original image we are in.
             */
            np = (int)(x*N);
            mp = (int)(y*M);
            xp = ((double)np+0.5)/(double)N;
            yp = ((double)mp+0.5)/(double)M;



            /*
             *  Find the edges in the original image that this pixel is closest to.
             */
            ne = (int)(x*N + 0.5);
            me = (int)(y*M + 0.5);
            xe = (double)ne/(double)N;
            ye = (double)me/(double)M;




            /*
             * Compute coords of the four closest pixels
             */
            //if        ( ( np < ne ) && ( mp < me ) ){  // we are in upper left  Pix #0
            if        ( ( x > xp ) && ( y > yp ) ){  // we are in upper left  Pix #0
                n0 = np;   m0 = mp;
                n1 = np+1; m1 = mp;
                n2 = np+1; m2 = mp+1;
                n3 = np;   m3 = mp+1;
            //} else if ( ( np > ne ) && ( mp < me ) ){  // we are in upper right Pix #1
            } else if ( ( x < xp ) && ( y > yp ) ){  // we are in upper right Pix #1
                n0 = np-1; m0 = mp;
                n1 = np;   m1 = mp;
                n2 = np;   m2 = mp+1;
                n3 = np-1; m3 = mp+1;
            //} else if ( ( np > ne ) && ( mp > me ) ){  // we are in lower right Pix #2
            } else if ( ( x < xp ) && ( y < yp ) ){  // we are in lower right Pix #2
                n0 = np-1; m0 = mp-1;
                n1 = np;   m1 = mp-1;
                n2 = np;   m2 = mp;
                n3 = np-1; m3 = mp;
            } else {                                   // we are in lower left  Pix #3
                n0 = np;   m0 = mp-1;
                n1 = np+1; m1 = mp-1;
                n2 = np+1; m2 = mp;
                n3 = np;   m3 = mp;
            }

  xe = 0.5*(n0+n1+1.0)/(double)N;
  ye = 0.5*(m1+m2+1.0)/(double)M;



            /*
             *  Check for edge problems
             */
            if (n0 < 0)   n0 = 0;
            if (n3 < 0)   n3 = 0;
            if (n1 > N-1) n1 = N-1;
            if (n2 > N-1) n2 = N-1;

            if (m0 < 0) m0 = 0;
            if (m1 < 0) m1 = 0;
            if (m2 > M-1)   m2 = M-1;
            if (m3 > M-1)   m3 = M-1;

            

            /*
             * Extract the pixel vals.
             */
            f[0] = Orig[ m0 ][ n0 ];
            f[1] = Orig[ m1 ][ n1 ];
            f[2] = Orig[ m2 ][ n2 ];
            f[3] = Orig[ m3 ][ n3 ];



            /*
             *  Compute (x',y') (i.e. the coords relative to center of 4 pixels.)
             */
            xx = x - xe;
            yy = y - ye;



            /*
             *  Compute weights (i.e. overlapping areas)
             */
            dxm = DXO2 - xx;
            dxp = DXO2 + xx;
            dym = DYO2 - yy;
            dyp = DYO2 + yy;
//if (dxm < 0.0) { printf("dxm = %g\n", dxm); exit(0); }
//if (dxp < 0.0) { printf("dxp = %g\n", dxp); exit(0); }
//if (dym < 0.0) { printf("dym = %g    DYO2 = %g  x, y = %g %g   xx, yy = %g %g\n", dym, DYO2, x, y, xx, yy); exit(0); }
//if (dyp < 0.0) { printf("dyp = %g\n", dyp); exit(0); }

            w0 = dxm*dym;
            w1 = dxp*dym;
            w2 = dxp*dyp;
            w3 = dxm*dyp;
            /*
            w0=0.1;
            w1=0.2;
            w2=0.3;
            w3=0.4;
            */

if ((j==189)&&(k==118)) {
                printf("np,mp = %d %d    pixels:  %d %d   %d %d   %d %d   %d %d\n", np, mp, n0, m0, n1, m1, n2, m2, n3, m3);
                printf("j,k = %d %d     f = %g %g %g %g     w = %g %g %g %g     xp, yp = %g %g   x,y = %g %g    xx, yy = %g %g\n", j, k, f[0], f[1], f[2], f[3], w0, w1, w2, w3, xp, yp, x, y, xx, yy);
}

            /*
             *  Finally compute weighted average and assign to current pixel of new image.
             *  A is the total area of an original pixel. w's are the partial overlap areas.
             */
            New[k][j] = (f[0]*w0 + f[1]*w1 + f[2]*w2 + f[3]*w3)/A;
            //New[k][j] = (f[0]*w0 + f[1]*w1 )/A;
            //New[k][j] = (f[0]*w0 )/A;
            //New[k][j] = (double)k;


        }
    }






}








