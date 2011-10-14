#include <stdio.h>
#include <stdlib.h>
#include <Lgm_SummersDiffCoeff.h>
#include <Lgm_DynamicMemory.h>
//#include <Lgm_FastPowPoly.h>



void DumpGif( char *FilenameBase, int W, int H, double **Image );


/*
 * Define a function that returns Bw(Latitude)
 */
typedef struct MyBwFuncInfo {

    double  whatever;

} MyBwFuncInfo;

double MyBwFunc( double Lat, void *Data ) {

    double          a[4], b[4], c[4];
    double          Bw, Bw2, y;
    MyBwFuncInfo    *MyInfo;

    // typecast generic data structure type to one we expect
    MyInfo = (MyBwFuncInfo *)Data;
    y = 0.75 + 0.04*Lat*DegPerRad;

    Bw = pow( 10.0, y );
//    Bw = (double)Lgm_FastPowPoly( 10.0, y );
//    a[0] = 10.0; a[1] = 10.0; a[2] = 10.0; a[3] = 10.0;
//    b[0] = y; b[1] = y+0.1; b[2] = y+0.2; b[3] = y+0.3;
//    Lgm_FastPowPoly_v( a, b, c );
    
//   Bw = c[0];
//    Bw = (double)powFastLookup( y, MyInfo->fp );

    return( Bw/1000.0 );
}




int main( ) {

    int             WaveMode, Species, i, j;
    double          aStar, Alpha, Daa_ba, Dap_ba, Dpp_ba, Sig;
    double          Ek, logEk, L, Beq, dB, Omega_e, wm, dw, w1, w2, MaxWaveLat;
    double          Alpha0, Alpha1, dAlpha, logEk0, logEk1, dlogEk;
    double          **ImageDaa, **ImageDap_neg, **ImageDap_pos, **ImageDpp;
    double          n1, n2, n3;
    MyBwFuncInfo    *MyInfo;


    MyInfo = (MyBwFuncInfo *)calloc( 1, sizeof(*MyInfo));


    Ek    = 1.000;  // Kinetic energy in MeV.
    L     = 4.5;    // L-shell parameter (dimensionless);.
    aStar = .16;    // Summer's cold plasma parameter
    aStar = 1.0/4.6;


    // Get Omega_e
    Beq = M_CDIP/(L*L*L);   // nT
    Omega_e = Lgm_GyroFreq( -LGM_e, Beq, LGM_ELECTRON_MASS );     // Omega_e, Hz


    // Set up waves
    dB       = 1.0;             // mean wave amplitude in nT.
    WaveMode = LGM_R_MODE_WAVE; // Wave-mode type (LGM_R_MODE_WAVE or LGM_L_MODE_WAVE).
    Species  = LGM_ELECTRONS;   // Species (LGM_ELECTRONS or LGM_PROTONS).

    






    Sig = 2.0;
    wm  = 0.2*Omega_e/M_2PI;     // Hz
    dw  = 0.1*Omega_e/M_2PI;     // Hz
    w1  = wm - Sig*dw/M_2PI;      // Hz
    w2  = wm + Sig*dw/M_2PI;      // Hz
w1 = 0.1*Omega_e/M_2PI;
w2 = 0.3*Omega_e/M_2PI;
    //MaxWaveLat = 35.0;      // Degrees
    MaxWaveLat = 90.0;      // Degrees

    n1 = 0.7; n2 = 0.2; n3 = 0.1;

    n1 = 1.0; n2 = 0.0; n3 = 0.0;

    dB       = 1.0;             // mean wave amplitude in nT.
    WaveMode = LGM_L_MODE_WAVE; // Wave-mode type (LGM_R_MODE_WAVE or LGM_L_MODE_WAVE).
    Species  = LGM_PROTONS;   // Species (LGM_ELECTRONS or LGM_PROTONS).
    wm = 3.7*Lgm_GyroFreq( -LGM_e, Beq, LGM_OXYGEN_MASS )/M_2PI;
    dw = .25*Lgm_GyroFreq( -LGM_e, Beq, LGM_OXYGEN_MASS )/M_2PI;
    w1 = 3.45*Lgm_GyroFreq( -LGM_e, Beq, LGM_OXYGEN_MASS )/M_2PI;;
    w2 = 3.95*Lgm_GyroFreq( -LGM_e, Beq, LGM_OXYGEN_MASS )/M_2PI;;

    aStar = 1.0/sqrt(15.0);













    int nAlpha  = 50; Alpha0 = 0.0; Alpha1 = 90.0; dAlpha = (Alpha1-Alpha0)/((double)(nAlpha-1));
    int nEnergy = 50; logEk0 = -1.0; logEk1 = 1.0; dlogEk = (logEk1-logEk0)/((double)(nEnergy-1));
    LGM_ARRAY_2D( ImageDaa,     nEnergy, nAlpha, double );
    LGM_ARRAY_2D( ImageDap_neg, nEnergy, nAlpha, double );
    LGM_ARRAY_2D( ImageDap_pos, nEnergy, nAlpha, double );
    LGM_ARRAY_2D( ImageDpp,     nEnergy, nAlpha, double );

    { /***** Start Parallel Execution ****/
//        #pragma omp parallel private(logEk, Ek, i, j, Alpha, Daa_ba, Dap_ba, Dpp_ba)
//        #pragma omp for schedule(dynamic, 8)
        for (i=0; i<nEnergy; i++ ){
            logEk = logEk0 + i*dlogEk;
            Ek = pow( 10.0, logEk );
            for (j=0; j<nAlpha; j++ ){
                Alpha = Alpha0 + j*dAlpha;

//if ( (j==10)&&(i==10)){
//if ( (j>=50)&&(j<=241)&&(i<=500-83)&&(i>=500-150)){
//if ( (j>=241)&&(j<=241)&&(i<=500-83)&&(i>=500-83)){
                Lgm_SummersDxxBounceAvg( LGM_SUMMERS_2007, Alpha, Ek, L, (void *)MyInfo, MyBwFunc, n1, n2, n3, aStar, w1, w2, wm,
                                         dw, WaveMode, Species, MaxWaveLat, &Daa_ba, &Dap_ba, &Dpp_ba );
                ImageDaa[i][j]     = Daa_ba;
                printf("i, j, Alpha = %g Dpp_ba = %d %d %g\n", i, j, Alpha, Dpp_ba);
                if ( Dap_ba < 0.0 ) {
                    ImageDap_neg[i][j] = fabs(Dap_ba);
                } else {
                    ImageDap_pos[i][j] = Dap_ba;
                }
                ImageDpp[i][j]     = Dpp_ba;
//}
            }
        }
    } /***** End Parallel Execution ****/

    DumpGif( "Daa_E_versus_Alpha_2007", nAlpha, nEnergy, ImageDaa );
    DumpGif( "Dap_neg_E_versus_Alpha_2007", nAlpha, nEnergy, ImageDap_neg );
    DumpGif( "Dap_pos_E_versus_Alpha_2007", nAlpha, nEnergy, ImageDap_pos );
    DumpGif( "Dpp_E_versus_Alpha_2007", nAlpha, nEnergy, ImageDpp );


    LGM_ARRAY_2D_FREE( ImageDaa );
    LGM_ARRAY_2D_FREE( ImageDap_neg );
    LGM_ARRAY_2D_FREE( ImageDap_pos );
    LGM_ARRAY_2D_FREE( ImageDpp );
    free( MyInfo );

    return(0);
}

