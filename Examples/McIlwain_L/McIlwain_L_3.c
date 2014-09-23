#include <Lgm_MagModelInfo.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

static double  glat[] = { 0.0, 20.0, 40.0, 60.0, 80.0 };
static double  glon[] = { 0.0, 135.0, 225.0 };
static double  alt[] = { 500.0 };


int main(){
    int                 i, j, k;
    long int            Date;
    double              L, I, Bm, M, a, UTC, B;
    Lgm_Vector          u, v, Bvec;
    Lgm_MagModelInfo    *mInfo = Lgm_InitMagInfo();
    Lgm_CTrans          *c = Lgm_init_ctrans( 0 );

    mInfo->Lgm_LossConeHeight = 100.0;
    mInfo->Bfield = Lgm_B_igrf;
    mInfo->Bfield = Lgm_B_JensenCain1960;
    mInfo->VerbosityLevel = 0;


    for ( j=0; j<3; j++ ){
        for ( i=0; i<5; i++ ){
            for ( k=0; k<1; k++ ){
    
                Lgm_GEOD_to_WGS84( glat[i], glon[j], alt[k], &u );
//Lgm_ScaleVector( &u, WGS84_A/6371.2 );
//printf("6371.2/WGS84_A = %g\n", 6371.2/WGS84_A);
    /*
                r = 1.0+alt[k]/WGS84_A;
                u.x = r*cos(glon[j]*RadPerDeg)*cos(glat[i]*RadPerDeg);
                u.y = r*sin(glon[j]*RadPerDeg)*cos(glat[i]*RadPerDeg);
                u.z = r*sin(glat[j]*RadPerDeg);
    */
                

        
                Date = 19600101;
                UTC = 0.0;
                Lgm_Set_Coord_Transforms( Date, UTC, c );
                Lgm_Convert_Coords( &u, &v, GEO_TO_GSM, c );

                
                
            
                a = 90.0;
                L = Lgm_McIlwain_L( Date, UTC, &v, a, 0, &I, &Bm, &M, mInfo );
                mInfo->Bfield( &v, &Bvec, mInfo );
                B = Lgm_Magnitude( &Bvec );
                printf("Pitch Angle: %g    Geodetic Lat/Lon/Rad: %3g %3g %3g    McIlwain L  = %.7g   ( I, Bm, M = %g %g %g )    |Pmin| = %g   B = %g\n", a, glat[i], glon[j], alt[k], L, I, Bm, M, Lgm_Magnitude( &mInfo->Pmin ), B );
            }
        }
    }

    Lgm_FreeMagInfo( mInfo );
    Lgm_free_ctrans( c );

    return(0);
}


