#include <stdio.h>
#include <Lgm_SummersDiffCoeff.h>
int main( ) {

    int     WaveMode, Species;
    double  SinAlpha, SinAlpha2, E, dBoverB2, BoverBeq, Rho, x1, x2, xm, dx, Lambda, s, aStar, Alpha, Daa, Sig;
    double  Ek, L, Beq, dB, Omega_e, Omega_p, Omega_SigEq, wm, dw, w1, w2, E0;
    double  Lat, CosLat, CosLat2, CosLat3, CosLat6, v, Alpha0, B, Omega_Sig;

    Ek       = 1.0000;   // Kinetic energy in MeV.
    L        = 7.0;   // L-shell parameter (dimensionless);.
    dB       = 0.1;   // mean wave amplitude in nT.
    WaveMode = LGM_R_MODE_WAVE; // Wave-mode type (LGM_R_MODE_WAVE or LGM_L_MODE_WAVE).
    Species  = LGM_ELECTRONS;   // Species (LGM_ELECTRONS or LGM_PROTONS).
    aStar    = 0.16;


    Beq = M_CDIP/(L*L*L);
    if ( Species == LGM_ELECTRONS ) {
        Lambda      = -1.0;
        E0          = LGM_Ee0;      // set rest energy to proton rest energy (MeV)
        Omega_SigEq = Lgm_GyroFreq( -LGM_e, Beq, LGM_ELECTRON_MASS );
    } else if ( Species == LGM_PROTONS ) {
        Lambda      = LGM_EPS;
        E0          = LGM_Ep0;       // set rest energy to electron rest energy (MeV)
        Omega_SigEq = Lgm_GyroFreq( LGM_e, Beq, LGM_PROTON_MASS );
    } else {
        printf("Unknown species = %d\n", Species );
        exit(0);
    }

    if ( WaveMode == LGM_R_MODE_WAVE ) {
        s = 1.0;
    } else if ( WaveMode == LGM_L_MODE_WAVE ) {
        s = -1.0;
    } else {
        printf("Unknown WabveMode = %d\n", WaveMode );
        exit(0);
    }








    // For a given latitude, compute BoverBeq
    Lat = 35.0*RadPerDeg;
    CosLat  = cos( Lat );     CosLat2 = CosLat*CosLat;
    CosLat3 = CosLat2*CosLat; CosLat6 = CosLat3*CosLat3;
    v = sqrt(4.0 - 3.0*CosLat2);
    BoverBeq  = v/CosLat6;              // B/Beq
    B = BoverBeq*M_CDIP/(L*L*L);

    // Get Omega_e and Omega_p
    Omega_e = Lgm_GyroFreq( -LGM_e, B, LGM_ELECTRON_MASS );     // Omega_e
    Omega_p = Lgm_GyroFreq(  LGM_e, B, LGM_PROTON_MASS );       // Omega_p
    Omega_Sig = Omega_SigEq*BoverBeq*BoverBeq;

    E = Ek/E0;
    dBoverB2 = dB*dB/(B*B);


    Sig = 2.0;
    wm  = 0.35*Omega_SigEq;
    dw  = 0.15*Omega_SigEq;
    w1  = wm - Sig*dw;
    w2  = wm + Sig*dw;
    Rho = ( erf( (wm-w1)/dw ) + erf( (w2-wm)/dw ) )/M_2_SQRTPI;
    x1 = w1/fabs(Omega_e);
    x2 = w2/fabs(Omega_e);
    xm = wm/fabs(Omega_e);
    dx = dw/fabs(Omega_e);



    for (Alpha=0.0; Alpha<=90.0; Alpha += 1.0 ) {

        // local PA
        SinAlpha = sin(Alpha*M_PI/180.0);
        SinAlpha2 = SinAlpha*SinAlpha;

        Daa = Lgm_SummersDaaLocal( SinAlpha2, E, dBoverB2, BoverBeq, Omega_e, Omega_Sig, Rho, Sig, x1, x2, xm, dx, Lambda, s, aStar );

        Alpha0 = DegPerRad*asin( sqrt( SinAlpha2/BoverBeq ) );

        printf("%g %g\n", Alpha0, Daa );

    }


    return(0);
}

