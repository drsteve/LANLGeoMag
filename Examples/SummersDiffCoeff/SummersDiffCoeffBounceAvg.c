#include <stdio.h>
#include <stdlib.h>
#include <Lgm_SummersDiffCoeff.h>
#include <Lgm_DynamicMemory.h>


void DumpGif( char *FilenameBase, int W, int H, double **Image );


/*
 * Define a function that returns Bw(Latitude)
 */
typedef struct MyBwFuncInfo {

    double  whatever;

} MyBwFuncInfo;

double MyBwFunc( double Lat, void *Data ) {

    double          Bw;
    MyBwFuncInfo    *MyInfo;

    // typecast generic data structure type to one we expect
    MyInfo = (MyBwFuncInfo *)Data;


    Bw = pow( 10.0, 0.75 + 0.04*Lat*DegPerRad );

    return( Bw );

}




int main( ) {

    int             WaveMode, Species, i, j;
    double          aStar, Alpha, Daa_ba, Dap_ba, Dpp_ba, Sig;
    double          Ek, logEk, L, Beq, dB, Omega_e, wm, dw, w1, w2, MaxWaveLat;
    double          Alpha0, Alpha1, dAlpha, logEk0, logEk1, dlogEk;
    double          **ImageDaa, **ImageDap_neg, **ImageDap_pos, **ImageDpp;
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
    dB       = 0.05;             // mean wave amplitude in nT.
    WaveMode = LGM_R_MODE_WAVE; // Wave-mode type (LGM_R_MODE_WAVE or LGM_L_MODE_WAVE).
    Species  = LGM_ELECTRONS;   // Species (LGM_ELECTRONS or LGM_PROTONS).
    Sig = 2.0;
    wm  = 0.2*Omega_e;     // Hz
    dw  = 0.1*Omega_e;     // Hz
    w1  = wm - Sig*dw;      // Hz
    w2  = wm + Sig*dw;      // Hz
w1 = 0.1*Omega_e;
w2 = 0.3*Omega_e;
    MaxWaveLat = 35.0;      // Degrees
    

    int nAlpha  = 100; Alpha0 = 0.0; Alpha1 = 90.0; dAlpha = (Alpha1-Alpha0)/((double)(nAlpha-1));
    int nEnergy = 100; logEk0 = -1.0; logEk1 = 1.0; dlogEk = (logEk1-logEk0)/((double)(nEnergy-1));
    LGM_ARRAY_2D( ImageDaa,     nEnergy, nAlpha, double );
    LGM_ARRAY_2D( ImageDap_neg, nEnergy, nAlpha, double );
    LGM_ARRAY_2D( ImageDap_pos, nEnergy, nAlpha, double );
    LGM_ARRAY_2D( ImageDpp,     nEnergy, nAlpha, double );

    {
        #pragma omp parallel private(logEk, Ek, j, Alpha, Daa_ba, Dap_ba, Dpp_ba)                                                                                                                                          
        #pragma omp for schedule(dynamic, 1)
        for (i=0; i<nEnergy; i++ ){
            logEk = logEk0 + i*dlogEk;
            Ek = pow( 10.0, logEk );
            for (j=0; j<nAlpha; j++ ){
                Alpha = Alpha0 + j*dAlpha;
                Lgm_SummersDxxBounceAvg( Alpha, Ek, L, (void *)MyInfo, MyBwFunc, aStar, w1, w2, wm, dw, WaveMode, Species, MaxWaveLat, &Daa_ba, &Dap_ba, &Dpp_ba );
                ImageDaa[i][j]     = Daa_ba;
                if ( Dap_ba < 0.0 ) {
                    ImageDap_neg[i][j] = fabs(Dap_ba);
                } else {
                    ImageDap_pos[i][j] = Dap_ba;
                }
                ImageDpp[i][j]     = Dpp_ba;
            }
        }
    }


    DumpGif( "Daa_E_versus_Alpha", nAlpha, nEnergy, ImageDaa );
    DumpGif( "Dap_neg_E_versus_Alpha", nAlpha, nEnergy, ImageDap_neg );
    DumpGif( "Dap_pos_E_versus_Alpha", nAlpha, nEnergy, ImageDap_pos );
    DumpGif( "Dpp_E_versus_Alpha", nAlpha, nEnergy, ImageDpp );

    LGM_ARRAY_2D_FREE( ImageDaa );
    LGM_ARRAY_2D_FREE( ImageDap_neg );
    LGM_ARRAY_2D_FREE( ImageDap_pos );
    LGM_ARRAY_2D_FREE( ImageDpp );

    free( MyInfo );

    return(0);
}

