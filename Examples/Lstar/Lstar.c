#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <Lgm_CTrans.h>
#include <Lgm_LstarInfo.h>


int main( int argc, char *argv[] ){

    double           UTC, alpha, a;
    long int         Date;
    int              nAlpha, Kp;
    int              ZeroOneTwo = 1, ThreeFourFive=1, SixSevenEight=1;
    Lgm_Vector       Psm, P, v1, v2, v3;
    Lgm_LstarInfo *LstarInfo = InitLstarInfo(0);

    // Date and UTC
    Date       = 19991122;
    UTC        = 19.0;
    Lgm_Set_Coord_Transforms( Date, UTC, LstarInfo->mInfo->c );

    // Position in SM
    Psm.x = 0.0; Psm.y = 6.6; Psm.z = 0.0;
    Lgm_Convert_Coords( &Psm, &P, SM_TO_GSM, LstarInfo->mInfo->c );

    // pitch angles to compute
    LstarInfo->PitchAngle = 87.5;


    Lgm_SetLstarTolerances( 3, 24, LstarInfo );

    LstarInfo->VerbosityLevel = 1;
    LstarInfo->mInfo->VerbosityLevel = 0;
    LstarInfo->mInfo->Bfield        = Lgm_B_T89c;
    LstarInfo->mInfo->InternalModel = LGM_IGRF;
    LstarInfo->mInfo->Kp = 4;

    /*
     * Compute L*s, Is, Bms, Footprints, etc...
     * These quantities are stored in the MagEphemInfo Structure
     */
    if (ZeroOneTwo) {
        LstarInfo->ShabanskyHandling = LGM_SHABANSKY_IGNORE;
        Lstar( &P, LstarInfo );
        printf("L* = %g\n", LstarInfo->LS);

        LstarInfo->ShabanskyHandling = LGM_SHABANSKY_HALVE_I;
        Lstar( &P, LstarInfo );
        printf("L* = %g\n", LstarInfo->LS);

        LstarInfo->ShabanskyHandling = LGM_SHABANSKY_REJECT;
        Lstar( &P, LstarInfo );
        printf("L* = %g\n", LstarInfo->LS);
    }


    Date       = 20130601;
    UTC        = 3.0 + 29.0/60.0;
    Lgm_Set_Coord_Transforms( Date, UTC, LstarInfo->mInfo->c );

    // Position in SM
    Psm.x = -4.286861; Psm.y = 3.207236; Psm.z = 1.710963;
    Lgm_Convert_Coords( &Psm, &P, SM_TO_GSM, LstarInfo->mInfo->c );

    if (ThreeFourFive) {
        // Set of Shabansky orbits that don't trigger a modified I (or if they do, don't notify)
        LstarInfo->PitchAngle = 25.0;
        LstarInfo->mInfo->Kp = 7;
        LstarInfo->mInfo->Bm = 0; //reset Bmirror in LstarInfo
    
        /*
         * Compute L*s, Is, Bms, Footprints, etc...
         * These quantities are stored in the MagEphemInfo Structure
         */
        LstarInfo->ShabanskyHandling = LGM_SHABANSKY_IGNORE;
        Lstar( &P, LstarInfo );
        printf("L* = %g\n", LstarInfo->LS);
    
        LstarInfo->ShabanskyHandling = LGM_SHABANSKY_HALVE_I;
        Lstar( &P, LstarInfo );
        printf("L* = %g\n", LstarInfo->LS);
    
        LstarInfo->ShabanskyHandling = LGM_SHABANSKY_REJECT;
        Lstar( &P, LstarInfo );
        printf("L* = %g\n", LstarInfo->LS);
    }

    if (SixSevenEight) {
        // Set of open (?) Shabansky orbits that that aren't correctly labeled
        LstarInfo->PitchAngle = 30.0;
        LstarInfo->mInfo->Kp = 7;
        LstarInfo->mInfo->Bm = 0; //reset Bmirror in LstarInfo
    
        /*
         * Compute L*s, Is, Bms, Footprints, etc...
         * These quantities are stored in the MagEphemInfo Structure
         */
        LstarInfo->ShabanskyHandling = LGM_SHABANSKY_IGNORE;
        Lstar( &P, LstarInfo );
        printf("L* = %g\n", LstarInfo->LS);
    
        LstarInfo->ShabanskyHandling = LGM_SHABANSKY_HALVE_I;
        Lstar( &P, LstarInfo );
        printf("L* = %g\n", LstarInfo->LS);
    
        LstarInfo->ShabanskyHandling = LGM_SHABANSKY_REJECT;
        Lstar( &P, LstarInfo );
        printf("L* = %g\n", LstarInfo->LS);
    }


    FreeLstarInfo( LstarInfo );

    return(0);

}

