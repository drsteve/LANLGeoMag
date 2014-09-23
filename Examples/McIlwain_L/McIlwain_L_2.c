#include <Lgm_MagModelInfo.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

static double  glat[] = { 17.0, 27.0, 38.0, 48.0, 57.0 };
static double  glon[] = { 192.0, 193.0, 195.0, 198.0, 201.0 };
static double  alt[] = { 100.0, 150.0, 200.0, 400.0, 1000.0 };


int main(){
    int                 i, j;
    long int            Date;
    double              L, I, Bm, M, a, UTC;
    Lgm_Vector          u, v;
    Lgm_MagModelInfo    *mInfo = Lgm_InitMagInfo();
    Lgm_CTrans          *c = Lgm_init_ctrans( 0 );

    mInfo->Lgm_LossConeHeight = 80.0;
    mInfo->Bfield = Lgm_B_JensenCain1960;
    mInfo->Bfield = Lgm_B_igrf;


    for ( i=0; i<5; i++ ){
        for ( j=0; j<5; j++ ){
    
            Lgm_GEOD_to_WGS84( glat[i], glon[i], alt[j], &u );
/*
            r = 1.0+alt[j]/WGS84_A;
            u.x = r*cos(glon[i]*RadPerDeg)*cos(glat[i]*RadPerDeg);
            u.y = r*sin(glon[i]*RadPerDeg)*cos(glat[i]*RadPerDeg);
            u.z = r*sin(glat[i]*RadPerDeg);
*/
            

    
            Date = 19850101;
            UTC = 0.0;
            Lgm_Set_Coord_Transforms( Date, UTC, c );
            Lgm_Convert_Coords( &u, &v, GEO_TO_GSM, c );
            
            
        
            a = 90.0;
            L = Lgm_McIlwain_L( Date, UTC, &v, a, 0, &I, &Bm, &M, mInfo );
            printf("Pitch Angle: %g    Geodetic LAt/Lon/Rad: %g %g %g    McIlwain L  = %.3g   ( I, Bm, M = %g %g %g )    |Pmin| = %g\n", a, glat[i], glon[i], alt[j], L, I, Bm, M, Lgm_Magnitude( &mInfo->Pmin ) );
        }
    }

    Lgm_FreeMagInfo(mInfo);
    Lgm_free_ctrans(c);

    return(0);
}


