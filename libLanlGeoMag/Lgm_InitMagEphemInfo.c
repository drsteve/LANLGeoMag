#include "Lgm/Lgm_MagEphemInfo.h"
#include "DynamicMemory.h"


/*
 *  A convenience routine to set some default values
 */

Lgm_MagEphemInfo *Lgm_InitMagEphemInfo( int Verbosity, int MaxPitchAngles ) {

    Lgm_MagEphemInfo  *MagEphemInfo = (Lgm_MagEphemInfo *) calloc (1, sizeof(*MagEphemInfo));
    MagEphemInfo->LstarInfo = InitLstarInfo( Verbosity );


    /*
     * Allocate Arrays that depend on # of pitch angles.
     */
    ARRAY_1D( MagEphemInfo->Alpha,        MaxPitchAngles, double );
    ARRAY_1D( MagEphemInfo->Pmn_gsm,      MaxPitchAngles, Lgm_Vector );
    ARRAY_1D( MagEphemInfo->Pms_gsm,      MaxPitchAngles, Lgm_Vector );
    ARRAY_1D( MagEphemInfo->Bm,           MaxPitchAngles, double );
    ARRAY_1D( MagEphemInfo->I,            MaxPitchAngles, double );
    ARRAY_1D( MagEphemInfo->Sb,           MaxPitchAngles, double );
    ARRAY_1D( MagEphemInfo->Tb,           MaxPitchAngles, double );
    ARRAY_1D( MagEphemInfo->K,            MaxPitchAngles, double );
    ARRAY_1D( MagEphemInfo->nShellPoints, MaxPitchAngles, int );

    ARRAY_2D( MagEphemInfo->ShellSphericalFootprint_Pn, MaxPitchAngles, 100, Lgm_Vector )
    ARRAY_2D( MagEphemInfo->ShellSphericalFootprint_Sn, MaxPitchAngles, 100, double );
    ARRAY_2D( MagEphemInfo->ShellSphericalFootprint_Bn, MaxPitchAngles, 100, double );
    ARRAY_2D( MagEphemInfo->ShellSphericalFootprint_Ps, MaxPitchAngles, 100, Lgm_Vector )
    ARRAY_2D( MagEphemInfo->ShellSphericalFootprint_Ss, MaxPitchAngles, 100, double );
    ARRAY_2D( MagEphemInfo->ShellSphericalFootprint_Bs, MaxPitchAngles, 100, double );
    ARRAY_2D( MagEphemInfo->ShellEllipsoidFootprint_Ps, MaxPitchAngles, 100, Lgm_Vector )
    ARRAY_2D( MagEphemInfo->ShellEllipsoidFootprint_Ss, MaxPitchAngles, 100, double );
    ARRAY_2D( MagEphemInfo->ShellEllipsoidFootprint_Bs, MaxPitchAngles, 100, double );
    ARRAY_2D( MagEphemInfo->ShellEllipsoidFootprint_Pn, MaxPitchAngles, 100, Lgm_Vector )
    ARRAY_2D( MagEphemInfo->ShellEllipsoidFootprint_Sn, MaxPitchAngles, 100, double );
    ARRAY_2D( MagEphemInfo->ShellEllipsoidFootprint_Bn, MaxPitchAngles, 100, double );
    ARRAY_2D( MagEphemInfo->ShellMirror_Pn,             MaxPitchAngles, 100, Lgm_Vector );
    ARRAY_2D( MagEphemInfo->ShellMirror_Sn,             MaxPitchAngles, 100, double );
    ARRAY_2D( MagEphemInfo->ShellMirror_Ps,             MaxPitchAngles, 100, Lgm_Vector );
    ARRAY_2D( MagEphemInfo->ShellMirror_Ss,             MaxPitchAngles, 100, double );
    ARRAY_2D( MagEphemInfo->ShellI,                     MaxPitchAngles, 100, double );

    ARRAY_2D( MagEphemInfo->nFieldPnts, MaxPitchAngles, 48, int );
    ARRAY_3D( MagEphemInfo->s_gsm,      MaxPitchAngles, 48, 1000, double );
    ARRAY_3D( MagEphemInfo->Bmag,       MaxPitchAngles, 48, 1000, double );
    ARRAY_3D( MagEphemInfo->x_gsm,      MaxPitchAngles, 48, 1000, double );
    ARRAY_3D( MagEphemInfo->y_gsm,      MaxPitchAngles, 48, 1000, double );
    ARRAY_3D( MagEphemInfo->z_gsm,      MaxPitchAngles, 48, 1000, double );

    ARRAY_1D( MagEphemInfo->LHilton,   MaxPitchAngles, double );
    ARRAY_1D( MagEphemInfo->LMcIlwain, MaxPitchAngles, double );
    ARRAY_1D( MagEphemInfo->Lstar,     MaxPitchAngles, double );



    



    return MagEphemInfo;

}



void Lgm_FreeMagEphemInfo( Lgm_MagEphemInfo  *MagEphemInfo ) {

    /*
     * De-allocate Arrays that depend on # of pitch angles.
     */
    ARRAY_1D_FREE( MagEphemInfo->Alpha );
    ARRAY_1D_FREE( MagEphemInfo->Pmn_gsm );
    ARRAY_1D_FREE( MagEphemInfo->Pms_gsm );
    ARRAY_1D_FREE( MagEphemInfo->Bm );
    ARRAY_1D_FREE( MagEphemInfo->I );
    ARRAY_1D_FREE( MagEphemInfo->Sb );
    ARRAY_1D_FREE( MagEphemInfo->Tb );
    ARRAY_1D_FREE( MagEphemInfo->K );
    ARRAY_1D_FREE( MagEphemInfo->nShellPoints );

    ARRAY_2D_FREE( MagEphemInfo->ShellSphericalFootprint_Pn );
    ARRAY_2D_FREE( MagEphemInfo->ShellSphericalFootprint_Sn );
    ARRAY_2D_FREE( MagEphemInfo->ShellSphericalFootprint_Bn );
    ARRAY_2D_FREE( MagEphemInfo->ShellSphericalFootprint_Ps );
    ARRAY_2D_FREE( MagEphemInfo->ShellSphericalFootprint_Ss );
    ARRAY_2D_FREE( MagEphemInfo->ShellSphericalFootprint_Bs );
    ARRAY_2D_FREE( MagEphemInfo->ShellEllipsoidFootprint_Ps );
    ARRAY_2D_FREE( MagEphemInfo->ShellEllipsoidFootprint_Ss );
    ARRAY_2D_FREE( MagEphemInfo->ShellEllipsoidFootprint_Bs );
    ARRAY_2D_FREE( MagEphemInfo->ShellEllipsoidFootprint_Pn );
    ARRAY_2D_FREE( MagEphemInfo->ShellEllipsoidFootprint_Sn );
    ARRAY_2D_FREE( MagEphemInfo->ShellEllipsoidFootprint_Bn );
    ARRAY_2D_FREE( MagEphemInfo->ShellMirror_Pn );
    ARRAY_2D_FREE( MagEphemInfo->ShellMirror_Sn );
    ARRAY_2D_FREE( MagEphemInfo->ShellMirror_Ps );
    ARRAY_2D_FREE( MagEphemInfo->ShellMirror_Ss );
    ARRAY_2D_FREE( MagEphemInfo->ShellI );

    ARRAY_2D_FREE( MagEphemInfo->nFieldPnts );
    ARRAY_3D_FREE( MagEphemInfo->s_gsm );
    ARRAY_3D_FREE( MagEphemInfo->Bmag );
    ARRAY_3D_FREE( MagEphemInfo->x_gsm );
    ARRAY_3D_FREE( MagEphemInfo->y_gsm );
    ARRAY_3D_FREE( MagEphemInfo->z_gsm );

    ARRAY_1D_FREE( MagEphemInfo->LHilton );
    ARRAY_1D_FREE( MagEphemInfo->LMcIlwain );
    ARRAY_1D_FREE( MagEphemInfo->Lstar );

    FreeLstarInfo( MagEphemInfo->LstarInfo );
    free( MagEphemInfo );

}
