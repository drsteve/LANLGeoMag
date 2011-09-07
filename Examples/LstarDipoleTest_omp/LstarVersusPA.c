#define MAIN
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <Lgm_CTrans.h>
#include <Lgm_MagEphemInfo.h>

#define KP_DEFAULT 0

//void LstarVersusPA( long int Date, double UTC, Lgm_Vector *u, int nAlpha, double *Alpha, int Quality, Lgm_MagEphemInfo *MagEphemInfo );

int main( int argc, char *argv[] ){

    double           UTC, Alpha[1000], a;
    long int         Date;
    int              nAlpha, Kp, Colorize, i;
    char             Filename[1024];
    Lgm_Vector       Psm, P;
    Lgm_CTrans       *c = Lgm_init_ctrans(0);
    Lgm_MagEphemInfo *MagEphemInfo;
    FILE             *fpout;

    // Date and UTC
    Date       = 19800625;
    UTC        = 19.0;
    Lgm_Set_Coord_Transforms( Date, UTC, c );

    // Position in SM
//    Psm.x = -1.5; Psm.y = 0.0; Psm.z = 0.0;
//    Psm.x = -1.25; Psm.y = 0.0; Psm.z = 0.0;
//    Psm.x = -6.6; Psm.y = 0.0; Psm.z = 0.0;
//    Psm.x = -1.05; Psm.y = 0.0; Psm.z = 0.0;
    Psm.x = -3.0; Psm.y = 0.0; Psm.z = 0.0;
    Lgm_Convert_Coords( &Psm, &P, SM_TO_GSM, c );

    // Create array of Pitch Angles to compute
    for (nAlpha=0,a=1.0; a<=90.0; a+=0.1,++nAlpha) {
        Alpha[nAlpha] = a ;
        printf("Alpha[%d] = %g\n", nAlpha, Alpha[nAlpha]);
    }
//nAlpha = 1;
//Alpha[0] = 3.20;

    if ( nAlpha > 0 ){
        MagEphemInfo = Lgm_InitMagEphemInfo(0, nAlpha);
    } else {
        // doesnt seem to like allocating zero size...
        MagEphemInfo = Lgm_InitMagEphemInfo(0, 1);
    }



//USER INPUT STUFF
    MagEphemInfo->LstarQuality   = 8;
    MagEphemInfo->SaveShellLines = TRUE;
    MagEphemInfo->LstarInfo->VerbosityLevel = 0;
    MagEphemInfo->LstarInfo->mInfo->VerbosityLevel = 0;

    Kp = 1;
    MagEphemInfo->LstarInfo->mInfo->Bfield        = Lgm_B_T89;
    MagEphemInfo->LstarInfo->mInfo->Bfield        = Lgm_B_cdip;
    MagEphemInfo->LstarInfo->mInfo->InternalModel = LGM_IGRF;
    MagEphemInfo->LstarInfo->mInfo->InternalModel = LGM_CDIP;
    MagEphemInfo->LstarInfo->mInfo->Kp = ( Kp >= 0 ) ? Kp : KP_DEFAULT;
    if ( MagEphemInfo->LstarInfo->mInfo->Kp > 5 ) MagEphemInfo->LstarInfo->mInfo->Kp = 5;



    /*
     * Compute L*s, Is, Bms, Footprints, etc...
     * These quantities are stored in the MagEphemInfo Structure
     */
    Colorize = TRUE;
    Lgm_ComputeLstarVersusPA( Date, UTC, &P, nAlpha, Alpha, MagEphemInfo->LstarQuality, Colorize, MagEphemInfo );

    /*
     * Dump results
     */
    sprintf( Filename, "DipoleTest_3.0/results_%.0e.dat", MagEphemInfo->LstarInfo->mInfo->Lgm_FindShellLine_I_Tol );
    fpout = fopen(Filename, "w");
    for (i=0; i<nAlpha; i++ ){
        if ( MagEphemInfo->Lstar[i] > 0.0 ) {
            fprintf(fpout, "%.15lf %.15lf\n", MagEphemInfo->Alpha[i], (3.00-MagEphemInfo->Lstar[i])/MagEphemInfo->LstarInfo->mInfo->Lgm_FindShellLine_I_Tol );
        }
    }
    fclose(fpout);

    Lgm_free_ctrans( c );
    Lgm_FreeMagEphemInfo( MagEphemInfo );

    return(0);

}





