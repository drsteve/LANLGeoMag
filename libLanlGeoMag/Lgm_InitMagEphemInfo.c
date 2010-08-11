#include "Lgm/Lgm_MagEphemInfo.h"


/*
 *  A convenience routine to set some default values
 */

Lgm_MagEphemInfo *Lgm_InitMagEphemInfo( int Verbosity ) {

    Lgm_MagEphemInfo  *MagEphemInfo = (Lgm_MagEphemInfo *) calloc (1, sizeof(*MagEphemInfo));

    MagEphemInfo->LstarInfo = InitLstarInfo( Verbosity );

    return MagEphemInfo;

}



void Lgm_FreeMagEphemInfo( Lgm_MagEphemInfo  *Info ) {

    FreeLstarInfo( Info->LstarInfo );
    free( Info );

}
