#include <stdio.h>
#include <Lgm_SummersDiffCoeff.h>
int main( ) {

    int     WaveMode, Species;
    double  SinAlpha, SinAlpha2, E, dBoverB2, BoverBeq, Rho, xm, dx, Lambda, s, aStar, Alpha, Daa, Sig;
    double  Ek, L, Beq, dB, Omega_e, Omega_p, wm, dw, w1, w2;

    Ek       = 0.1100;   // Kinetic energy in MeV.
    L        = 4.0;   // L-shell parameter (dimensionless);.
    dB       = 1.0;   // mean wave amplitude in nT.
    WaveMode = LGM_R_MODE_WAVE; // Wave-mode type (LGM_R_MODE_WAVE or LGM_L_MODE_WAVE).
    Species  = LGM_ELECTRONS;   // Species (LGM_ELECTRONS or LGM_PROTONS).
Lambda = LGM_EPS;
s = -1.0;
        aStar = 4.6e-3;



    Beq = M_CDIP/(L*L*L);

    // Assume we are at the equator
    BoverBeq = 1.0;

    // Get Omega_e and Omega_p
    Omega_e = Lgm_GyroFreq( -LGM_e, Beq, LGM_ELECTRON_MASS );     // Omega_e
    Omega_p = Lgm_GyroFreq(  LGM_e, Beq, LGM_PROTON_MASS );       // Omega_p

    E = Ek/LGM_Ep0;
    dBoverB2 = dB*dB/(Beq*Beq);


    Sig = 4.0/3.0;
    wm  = 0.15*Omega_p;
    dw  = 0.0375*Omega_p;
    w1  = wm - Sig*dw;
    w2  = wm + Sig*dw;
    Rho = ( erf( (wm-w1)/dw ) + erf( (w2-wm)/dw ) )/M_2_SQRTPI;
    xm = wm/fabs(Omega_e);
    dx = dw/fabs(Omega_e);


    
    for (Alpha=0.0; Alpha<=90.0; Alpha += 1.0 ) {
    
        // local PA
        SinAlpha = sin(Alpha*M_PI/180.0);
        SinAlpha2 = SinAlpha*SinAlpha;

        Daa = Lgm_SummersDaaLocal( SinAlpha2, E, dBoverB2, BoverBeq, Omega_e, Omega_p, Rho, Sig, xm, dx, Lambda, s, aStar );

        printf("%g %g\n", Alpha, Daa );
        
    }


    return(0);
}

