#include <Lgm_MagModelInfo.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

int main(){
    long int            Date, n;
    double              L, I, Bm, M, a, UTC, B, glat, glon, alt;
    Lgm_Vector          u, v, Bvec;
    Lgm_MagModelInfo    *mInfo = Lgm_InitMagInfo();
    Lgm_CTrans          *c = Lgm_init_ctrans( 0 );

    mInfo->Lgm_LossConeHeight = 100.0;
    mInfo->Bfield = Lgm_B_JensenCain1960;
    mInfo->Bfield = Lgm_B_igrf;
    mInfo->VerbosityLevel = 0;
    mInfo->MaxDiv = 0.1;

    mInfo->Lgm_I_Integrator = DQAGS;
    mInfo->nDivs  = 200.0;
//    mInfo->Lgm_TraceToBmin_Tol = 1e-2;
//    mInfo->Lgm_TraceToEarth_Tol = 1e-1;
//    mInfo->Lgm_TraceToMirrorPoint_Tol = 1e-5;

    Date = 20110101;
    UTC = 6.0;
    Lgm_Set_Coord_Transforms( Date, UTC, c );

    n = 0;
    for ( glon=0.0; glon<360; glon += 60.0 ){
        for ( glat=-60.0; glat <=60.0; glat += 20.0 ){
            for ( alt=200.0; alt<=1000.0; alt += 200.0 ){
    
                Lgm_GEOD_to_WGS84( glat, glon, alt, &u );
                Lgm_Convert_Coords( &u, &v, GEO_TO_GSM, c );
        

                a = 90.0;
                L = Lgm_McIlwain_L( Date, UTC, &v, a, 1, &I, &Bm, &M, mInfo );
                mInfo->Bfield( &v, &Bvec, mInfo );
                B = Lgm_Magnitude( &Bvec );
                printf("Pitch Angle: %g    Geodetic Lat/Lon/Rad: %3g %3g %3g    McIlwain L  = %.7g   ( I, Bm, M = %g %g %g )  mInfo->nPnts = %d mInfo->Hmax = %g\n", a, glat, glon, alt, L, I, Bm, M, mInfo->nPnts, mInfo->Hmax );
                ++n;
            }
        }
    }
    printf("n = %ld\n", n);

    return(0);
}


