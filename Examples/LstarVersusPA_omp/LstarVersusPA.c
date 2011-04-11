#define MAIN
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <Lgm_CTrans.h>
#include <Lgm_MagEphemInfo.h>

#define KP_DEFAULT 0

void LstarVersusPA( long int Date, double UTC, Lgm_Vector *u, int nAlpha, double *Alpha, int Quality, Lgm_MagEphemInfo *MagEphemInfo );

int main( int argc, char *argv[] ){

    double           UTC, Alpha[1000], a;
    long int         Date;
    int              nAlpha, Kp, i;
    char             Filename2[1024];
    Lgm_Vector       Psm, P;
    Lgm_CTrans       *c = Lgm_init_ctrans(0);
    Lgm_MagEphemInfo *MagEphemInfo = Lgm_InitMagEphemInfo(0, 100);
    FILE             *fp;


    // Date and UTC
    Date       = 20010425;
    UTC        = 16.0 + 54.0/60.0 + 9.595436/3600.0;
    Lgm_Set_Coord_Transforms( Date, UTC, c );

    // Position in SM
//    Psm.x = -6.6; Psm.y = 0.0; Psm.z = 0.0;
//    Lgm_Convert_Coords( &Psm, &P, SM_TO_GSM, c );
P.x = -4.847796682966851;
P.y = -1.9408472800988594;
P.z = -1.7773158425650286;

    // Create array of Pitch Angles to compute
    for (nAlpha=0,a=5.0; a<=90.0; a+=5.0,++nAlpha) {
        Alpha[nAlpha] = a ;
        //printf("Alpha[%d] = %g\n", nAlpha, Alpha[nAlpha]);
    }
nAlpha = 1;
Alpha[0] = 74.0890238603;





//USER INPUT STUFF
    MagEphemInfo->LstarQuality   = 2;
    MagEphemInfo->SaveShellLines = TRUE;
    MagEphemInfo->LstarInfo->VerbosityLevel = 0;
    MagEphemInfo->LstarInfo->mInfo->VerbosityLevel = 0;

    Kp = 5;
    MagEphemInfo->LstarInfo->mInfo->Bfield        = Lgm_B_cdip;
    MagEphemInfo->LstarInfo->mInfo->Bfield        = Lgm_B_T89;
    MagEphemInfo->LstarInfo->mInfo->InternalModel = LGM_CDIP;
    MagEphemInfo->LstarInfo->mInfo->InternalModel = LGM_IGRF;
    MagEphemInfo->LstarInfo->mInfo->Kp = ( Kp >= 0 ) ? Kp : KP_DEFAULT;
    if ( MagEphemInfo->LstarInfo->mInfo->Kp > 5 ) MagEphemInfo->LstarInfo->mInfo->Kp = 5;

MagEphemInfo->LstarInfo->mInfo->Bfield        = Lgm_B_OP77;


    /*
     * Compute L*s, Is, Bms, Footprints, etc...
     * These quantities are stored in the MagEphemInfo Structure
     */
    ComputeLstarVersusPA( Date, UTC, &P, nAlpha, Alpha, MagEphemInfo->LstarQuality, MagEphemInfo );

    fp = fopen("Lstar.dat", "w");
    for ( i=0; i<nAlpha; ++i ) {
        fprintf( fp, "%g %.10lf\n", Alpha[i], MagEphemInfo->Lstar[i] );
    }
    fclose(fp);

    Lgm_free_ctrans( c );
    Lgm_FreeMagEphemInfo( MagEphemInfo );

    return(0);

}





