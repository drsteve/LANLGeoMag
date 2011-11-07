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
    y = 0.75 + 0.04*fabs(Lat)*DegPerRad;
    Bw = pow( 10.0, y );
//Bw = 5000.0;

if (Lat*DegPerRad > 35.0) return(0.0);
else return( Bw/1000.0 );
}




int main( ) {

    int             WaveMode, Species, i, j;
    double          aStar, Alpha, Daa_ba, Dap_ba, Dpp_ba, Sig;
    double          Ek, logEk, L, Beq, dB, Omega_e, wm, dw, w1, w2, MaxWaveLat;
    double          Alpha0, Alpha1, dAlpha, logEk0, logEk1, dlogEk;
    double          **ImageDaa, **ImageDap_neg, **ImageDap_pos, **ImageDpp;
    double          **ImageDaa_Li, **ImageDpp_Li;
    double          **ImageDaa_Diff, **ImageDpp_Diff;
    double          n1, n2, n3;
    MyBwFuncInfo    *MyInfo;
    FILE            *fp1, *fp2;


    MyInfo = (MyBwFuncInfo *)calloc( 1, sizeof(*MyInfo));


    Ek    = 1.000;  // Kinetic energy in MeV.
    L     = 4.5;    // L-shell parameter (dimensionless);.
    aStar = .16;    // Summer's cold plasma parameter
    aStar = 1.0/(3.8*3.8);


    // Get Omega_e
    Beq = M_CDIP/(L*L*L);   // nT
    Omega_e = Lgm_GyroFreq( LGM_e, Beq, LGM_ELECTRON_MASS );     // Omega_e, Hz


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
    MaxWaveLat = 15.0;      // Degrees









    /*
     *   Li et al. Figure 2c
     */
    WaveMode = LGM_R_MODE_WAVE; // Wave-mode type (LGM_R_MODE_WAVE or LGM_L_MODE_WAVE).
    Species  = LGM_ELECTRONS;   // Species (LGM_ELECTRONS or LGM_PROTONS).
    wm = 0.35*Lgm_GyroFreq( LGM_e, Beq, LGM_ELECTRON_MASS )/M_2PI;
    dw = .15*Lgm_GyroFreq( LGM_e, Beq, LGM_ELECTRON_MASS )/M_2PI;
    w1 = 0.05*Lgm_GyroFreq( LGM_e, Beq, LGM_ELECTRON_MASS )/M_2PI;;
    w2 = 0.65*Lgm_GyroFreq( LGM_e, Beq, LGM_ELECTRON_MASS )/M_2PI;;
    MaxWaveLat = 15.0;      // Degrees
    aStar = 1.0/(3.8*3.8);
    n1 = 0.7; n2 = 0.2; n3 = 0.1;


    /*
     *   Li et al. Figure 4
     */
    WaveMode = LGM_L_MODE_WAVE; // Wave-mode type (LGM_R_MODE_WAVE or LGM_L_MODE_WAVE).
    Species  = LGM_ELECTRONS;   // Species (LGM_ELECTRONS or LGM_PROTONS).
    wm = 3.5*Lgm_GyroFreq( LGM_e, Beq, LGM_OXYGEN_MASS )/M_2PI;
    dw = .25*Lgm_GyroFreq( LGM_e, Beq, LGM_OXYGEN_MASS )/M_2PI;
    w1 = 3.25*Lgm_GyroFreq( LGM_e, Beq, LGM_OXYGEN_MASS )/M_2PI;;
    w2 = 3.75*Lgm_GyroFreq( LGM_e, Beq, LGM_OXYGEN_MASS )/M_2PI;;
    MaxWaveLat = 15.0;      // Degrees
    aStar = 1.0/(4.0*4.0);
    n1 = 0.7; n2 = 0.2; n3 = 0.1;

    /*
     *   Li et al. Figure 2a
     */
    WaveMode = LGM_R_MODE_WAVE; // Wave-mode type (LGM_R_MODE_WAVE or LGM_L_MODE_WAVE).
    Species  = LGM_ELECTRONS;   // Species (LGM_ELECTRONS or LGM_PROTONS).
    wm = 0.2*Lgm_GyroFreq( LGM_e, Beq, LGM_ELECTRON_MASS )/M_2PI;
    dw = .1*Lgm_GyroFreq( LGM_e, Beq, LGM_ELECTRON_MASS )/M_2PI;
    w1 = 0.1*Lgm_GyroFreq( LGM_e, Beq, LGM_ELECTRON_MASS )/M_2PI;;
    w2 = 0.3*Lgm_GyroFreq( LGM_e, Beq, LGM_ELECTRON_MASS )/M_2PI;;
    MaxWaveLat = 35.0;      // Degrees
    aStar = 1.0/(4.6*4.6);
    n1 = 0.7; n2 = 0.2; n3 = 0.1;



    int nAlpha  = 90; Alpha0 = 1.0; Alpha1 = 90.0; dAlpha = (Alpha1-Alpha0)/((double)(nAlpha-1));
    int nEnergy = 100; logEk0 = -1.0; logEk1 = 1.0; dlogEk = (logEk1-logEk0)/((double)(nEnergy-1));
    LGM_ARRAY_2D( ImageDaa,     nEnergy, nAlpha, double );
    LGM_ARRAY_2D( ImageDap_neg, nEnergy, nAlpha, double );
    LGM_ARRAY_2D( ImageDap_pos, nEnergy, nAlpha, double );
    LGM_ARRAY_2D( ImageDpp,     nEnergy, nAlpha, double );

    { /***** Start Parallel Execution ****/
        #pragma omp parallel private(logEk, Ek, i, j, Alpha, Daa_ba, Dap_ba, Dpp_ba)
        #pragma omp for schedule(dynamic, 1)
        for (i=0; i<nEnergy; i++ ){
        //for (i=0; i<nEnergy; i+=10 ){
            logEk = logEk0 + i*dlogEk;
            Ek = pow( 10.0, logEk );
            for (j=0; j<nAlpha; j++ ){
            //for (j=0; j<nAlpha; j+=2 ){
                Alpha = Alpha0 + j*dAlpha;

//if ( (i==20)){
//if ( (i==20)&&(j==86)){
//if ( (j>=50)&&(j<=241)&&(i<=500-83)&&(i>=500-150)){
//if ( (j>=241)&&(j<=241)&&(i<=500-83)&&(i>=500-83)){
                Lgm_SummersDxxBounceAvg( LGM_SUMMERS_2007, Alpha, Ek, L, (void *)MyInfo, MyBwFunc, n1, n2, n3, aStar, LGM_FRWD_BKWD, w1, w2, wm,
                                         dw, WaveMode, Species, MaxWaveLat, &Daa_ba, &Dap_ba, &Dpp_ba );
                ImageDaa[i][j]     = Daa_ba*86400.0;
                printf("i = %d , j = %d , E = %g  Alpha = %g Daa_ba = %g\n", i, j, Ek, Alpha, Daa_ba);
                if ( Dap_ba < 0.0 ) {
                    ImageDap_neg[i][j] = fabs(Dap_ba*86400.0);
                } else {
                    ImageDap_pos[i][j] = Dap_ba*86400.0;
                }
                ImageDpp[i][j]     = Dpp_ba*86400.0*4.0;
//}
            }
        }
    } /***** End Parallel Execution ****/

    DumpGif2( "Daa_E_versus_Alpha_2007_3", -6.0, 2.0, nAlpha, nEnergy, ImageDaa );
    DumpGif2( "Dap_neg_E_versus_Alpha_2007_3", -6.0, 2.0, nAlpha, nEnergy, ImageDap_neg );
    DumpGif2( "Dap_pos_E_versus_Alpha_2007_3", -6.0, 2.0, nAlpha, nEnergy, ImageDap_pos );
    DumpGif2( "Dpp_E_versus_Alpha_2007_3", -6.0, 2.0, nAlpha, nEnergy, ImageDpp );
    /*
    DumpGif( "Daa_E_versus_Alpha_2007_3", nAlpha, nEnergy, ImageDaa );
    DumpGif( "Dap_neg_E_versus_Alpha_2007_3", nAlpha, nEnergy, ImageDap_neg );
    DumpGif( "Dap_pos_E_versus_Alpha_2007_3", nAlpha, nEnergy, ImageDap_pos );
    DumpGif( "Dpp_E_versus_Alpha_2007_3", nAlpha, nEnergy, ImageDpp );
*/





    /*
     * Do a diff on Li's results
     */
    LGM_ARRAY_2D( ImageDaa_Li,     nEnergy, nAlpha, double );
    LGM_ARRAY_2D( ImageDpp_Li,     nEnergy, nAlpha, double );
    LGM_ARRAY_2D( ImageDaa_Diff,   nEnergy, nAlpha, double );
    LGM_ARRAY_2D( ImageDpp_Diff,   nEnergy, nAlpha, double );


    fp1 = fopen("Li/diffusion coef/chorus_Daa_2d_L45_lam35_day.dat","r");
    fp2 = fopen("Li/diffusion coef/chorus_Dpp_2d_L45_lam35_day.dat","r");
double d;
    for (i=0; i<nEnergy; i++ ){
        for (j=0; j<nAlpha; j++ ){
            fscanf( fp1, "%lf", &ImageDaa_Li[i][j]);
            fscanf( fp2, "%lf", &ImageDpp_Li[i][j]);
            
        }
    }
    fclose(fp1);
    fclose(fp2);

    for (i=0; i<nEnergy; i++ ){
        for (j=0; j<nAlpha; j++ ){
            d = ImageDaa_Li[i][j]/ImageDaa[i][j];
            ImageDaa_Diff[i][j] = (d>=0.0) ? d : 0.0;
            d = ImageDpp_Li[i][j]/ImageDpp[i][j];
            ImageDpp_Diff[i][j] = (d>=0.0) ? d : 0.0;
        }
    }

    DumpGif2( "Daa_E_versus_Alpha_Li", -6.0, 2.0, nAlpha, nEnergy, ImageDaa_Li );
    DumpGif2( "Dpp_E_versus_Alpha_Li", -6.0, 2.0, nAlpha, nEnergy, ImageDpp_Li );
    
    DumpGif2( "Daa_E_versus_Alpha_Diff", -1.0, 1.0, nAlpha, nEnergy, ImageDaa_Diff );
    DumpGif2( "Dpp_E_versus_Alpha_Diff", -1.0, 1.0, nAlpha, nEnergy, ImageDpp_Diff );
    









    LGM_ARRAY_2D_FREE( ImageDaa );
    LGM_ARRAY_2D_FREE( ImageDap_neg );
    LGM_ARRAY_2D_FREE( ImageDap_pos );
    LGM_ARRAY_2D_FREE( ImageDpp );
    LGM_ARRAY_2D_FREE( ImageDaa_Li );
    LGM_ARRAY_2D_FREE( ImageDpp_Li );
    LGM_ARRAY_2D_FREE( ImageDaa_Diff );
    LGM_ARRAY_2D_FREE( ImageDpp_Diff );
    free( MyInfo );

    return(0);
}

