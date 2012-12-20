#include "Lgm/Lgm_MagEphemInfo.h"
#include "Lgm/Lgm_DynamicMemory.h"


/*
 *  A convenience routine to set some default values
 */

Lgm_MagEphemInfo *Lgm_InitMagEphemInfo( int Verbosity, int MaxPitchAngles ) {
    Lgm_MagEphemInfo  *MagEphemInfo = (Lgm_MagEphemInfo *) calloc (1, sizeof(*MagEphemInfo));
    Lgm_InitMagEphemInfoDefaults(MagEphemInfo, MaxPitchAngles, Verbosity );
    //printf("MagEphemInfo = %p\n", MagEphemInfo);
    return MagEphemInfo;
}

void Lgm_InitMagEphemInfoDefaults(Lgm_MagEphemInfo *MagEphemInfo, int MaxPitchAngles, int Verbosity) {

    int i;

    MagEphemInfo->LstarInfo = InitLstarInfo( Verbosity );

    MagEphemInfo->LstarInfo->SaveShellLines = TRUE;
    MagEphemInfo->SaveShellLines = TRUE;


    /*
     * Allocate Arrays that depend on # of pitch angles.
     */
    LGM_ARRAY_1D( MagEphemInfo->Alpha,          MaxPitchAngles, double );
    LGM_ARRAY_1D( MagEphemInfo->Pmn_gsm,        MaxPitchAngles, Lgm_Vector );
    LGM_ARRAY_1D( MagEphemInfo->Pms_gsm,        MaxPitchAngles, Lgm_Vector );
    LGM_ARRAY_1D( MagEphemInfo->Bm,             MaxPitchAngles, double );
    LGM_ARRAY_1D( MagEphemInfo->I,              MaxPitchAngles, double );
    LGM_ARRAY_1D( MagEphemInfo->Sb,             MaxPitchAngles, double );
    LGM_ARRAY_1D( MagEphemInfo->Tb,             MaxPitchAngles, double );
    LGM_ARRAY_1D( MagEphemInfo->K,              MaxPitchAngles, double );
    LGM_ARRAY_1D( MagEphemInfo->nShellPoints,   MaxPitchAngles, int );
    LGM_ARRAY_1D( MagEphemInfo->LHilton,        MaxPitchAngles, double );
    LGM_ARRAY_1D( MagEphemInfo->LMcIlwain,      MaxPitchAngles, double );
    LGM_ARRAY_1D( MagEphemInfo->Lstar,          MaxPitchAngles, double );
    LGM_ARRAY_1D( MagEphemInfo->DriftOrbitType, MaxPitchAngles, int );


    for ( i=0; i<MaxPitchAngles; ++i ){
        MagEphemInfo->Pmn_gsm[i].x      = MagEphemInfo->Pmn_gsm[i].y = MagEphemInfo->Pmn_gsm[i].z = LGM_FILL_VALUE;
        MagEphemInfo->Pms_gsm[i].x      = MagEphemInfo->Pms_gsm[i].y = MagEphemInfo->Pms_gsm[i].z = LGM_FILL_VALUE;
        MagEphemInfo->Alpha[i]          = LGM_FILL_VALUE;
        MagEphemInfo->Bm[i]             = LGM_FILL_VALUE;
        MagEphemInfo->I[i]              = LGM_FILL_VALUE;
        MagEphemInfo->Sb[i]             = LGM_FILL_VALUE;
        MagEphemInfo->Tb[i]             = LGM_FILL_VALUE;
        MagEphemInfo->K[i]              = LGM_FILL_VALUE;
        MagEphemInfo->nShellPoints[i]   = 0;
        MagEphemInfo->LHilton[i]        = LGM_FILL_VALUE;
        MagEphemInfo->LMcIlwain[i]      = LGM_FILL_VALUE;
        MagEphemInfo->Lstar[i]          = LGM_FILL_VALUE;
        MagEphemInfo->DriftOrbitType[i] = LGM_DRIFT_ORBIT_OPEN;
    }

    LGM_ARRAY_2D( MagEphemInfo->ShellSphericalFootprint_Pn, MaxPitchAngles, 100, Lgm_Vector );
    LGM_ARRAY_2D( MagEphemInfo->ShellSphericalFootprint_Sn, MaxPitchAngles, 100, double );
    LGM_ARRAY_2D( MagEphemInfo->ShellSphericalFootprint_Bn, MaxPitchAngles, 100, double );
    LGM_ARRAY_2D( MagEphemInfo->ShellSphericalFootprint_Ps, MaxPitchAngles, 100, Lgm_Vector );
    LGM_ARRAY_2D( MagEphemInfo->ShellSphericalFootprint_Ss, MaxPitchAngles, 100, double );
    LGM_ARRAY_2D( MagEphemInfo->ShellSphericalFootprint_Bs, MaxPitchAngles, 100, double );
    LGM_ARRAY_2D( MagEphemInfo->ShellEllipsoidFootprint_Ps, MaxPitchAngles, 100, Lgm_Vector );
    LGM_ARRAY_2D( MagEphemInfo->ShellEllipsoidFootprint_Ss, MaxPitchAngles, 100, double );
    LGM_ARRAY_2D( MagEphemInfo->ShellEllipsoidFootprint_Bs, MaxPitchAngles, 100, double );
    LGM_ARRAY_2D( MagEphemInfo->ShellEllipsoidFootprint_Pn, MaxPitchAngles, 100, Lgm_Vector );
    LGM_ARRAY_2D( MagEphemInfo->ShellEllipsoidFootprint_Sn, MaxPitchAngles, 100, double );
    LGM_ARRAY_2D( MagEphemInfo->ShellEllipsoidFootprint_Bn, MaxPitchAngles, 100, double );
    LGM_ARRAY_2D( MagEphemInfo->ShellMirror_Pn,             MaxPitchAngles, 100, Lgm_Vector );
    LGM_ARRAY_2D( MagEphemInfo->ShellMirror_Sn,             MaxPitchAngles, 100, double );
    LGM_ARRAY_2D( MagEphemInfo->ShellMirror_Ps,             MaxPitchAngles, 100, Lgm_Vector );
    LGM_ARRAY_2D( MagEphemInfo->ShellMirror_Ss,             MaxPitchAngles, 100, double );
    LGM_ARRAY_2D( MagEphemInfo->ShellI,                     MaxPitchAngles, 100, double );
    LGM_ARRAY_2D( MagEphemInfo->Shell_Bmin,                 MaxPitchAngles, 100, Lgm_Vector );
    LGM_ARRAY_2D( MagEphemInfo->Shell_Pmin,                 MaxPitchAngles, 100, Lgm_Vector );
    LGM_ARRAY_2D( MagEphemInfo->Shell_GradI,                MaxPitchAngles, 100, Lgm_Vector );
    LGM_ARRAY_2D( MagEphemInfo->Shell_Vgc,                  MaxPitchAngles, 100, Lgm_Vector );
    LGM_ARRAY_2D( MagEphemInfo->nMinima,                    MaxPitchAngles, 100, int );
    LGM_ARRAY_2D( MagEphemInfo->nMaxima,                    MaxPitchAngles, 100, int );

    LGM_ARRAY_2D( MagEphemInfo->nFieldPnts, MaxPitchAngles, 48, int );
    LGM_ARRAY_3D( MagEphemInfo->s_gsm,      MaxPitchAngles, 48, 1000, double );
    LGM_ARRAY_3D( MagEphemInfo->Bmag,       MaxPitchAngles, 48, 1000, double );
    LGM_ARRAY_3D( MagEphemInfo->x_gsm,      MaxPitchAngles, 48, 1000, double );
    LGM_ARRAY_3D( MagEphemInfo->y_gsm,      MaxPitchAngles, 48, 1000, double );
    LGM_ARRAY_3D( MagEphemInfo->z_gsm,      MaxPitchAngles, 48, 1000, double );


}


void Lgm_FreeMagEphemInfo_Children( Lgm_MagEphemInfo  *MagEphemInfo ) {
    LGM_ARRAY_1D_FREE( MagEphemInfo->Alpha );
    LGM_ARRAY_1D_FREE( MagEphemInfo->Pmn_gsm );
    LGM_ARRAY_1D_FREE( MagEphemInfo->Pms_gsm );
    LGM_ARRAY_1D_FREE( MagEphemInfo->Bm );
    LGM_ARRAY_1D_FREE( MagEphemInfo->I );
    LGM_ARRAY_1D_FREE( MagEphemInfo->Sb );
    LGM_ARRAY_1D_FREE( MagEphemInfo->Tb );
    LGM_ARRAY_1D_FREE( MagEphemInfo->K );
    LGM_ARRAY_1D_FREE( MagEphemInfo->nShellPoints );

    LGM_ARRAY_2D_FREE( MagEphemInfo->ShellSphericalFootprint_Pn );
    LGM_ARRAY_2D_FREE( MagEphemInfo->ShellSphericalFootprint_Sn );
    LGM_ARRAY_2D_FREE( MagEphemInfo->ShellSphericalFootprint_Bn );
    LGM_ARRAY_2D_FREE( MagEphemInfo->ShellSphericalFootprint_Ps );
    LGM_ARRAY_2D_FREE( MagEphemInfo->ShellSphericalFootprint_Ss );
    LGM_ARRAY_2D_FREE( MagEphemInfo->ShellSphericalFootprint_Bs );
    LGM_ARRAY_2D_FREE( MagEphemInfo->ShellEllipsoidFootprint_Ps );
    LGM_ARRAY_2D_FREE( MagEphemInfo->ShellEllipsoidFootprint_Ss );
    LGM_ARRAY_2D_FREE( MagEphemInfo->ShellEllipsoidFootprint_Bs );
    LGM_ARRAY_2D_FREE( MagEphemInfo->ShellEllipsoidFootprint_Pn );
    LGM_ARRAY_2D_FREE( MagEphemInfo->ShellEllipsoidFootprint_Sn );
    LGM_ARRAY_2D_FREE( MagEphemInfo->ShellEllipsoidFootprint_Bn );
    LGM_ARRAY_2D_FREE( MagEphemInfo->ShellMirror_Pn );
    LGM_ARRAY_2D_FREE( MagEphemInfo->ShellMirror_Sn );
    LGM_ARRAY_2D_FREE( MagEphemInfo->ShellMirror_Ps );
    LGM_ARRAY_2D_FREE( MagEphemInfo->ShellMirror_Ss );
    LGM_ARRAY_2D_FREE( MagEphemInfo->ShellI );
    LGM_ARRAY_2D_FREE( MagEphemInfo->Shell_Bmin );
    LGM_ARRAY_2D_FREE( MagEphemInfo->Shell_Pmin );
    LGM_ARRAY_2D_FREE( MagEphemInfo->Shell_GradI );
    LGM_ARRAY_2D_FREE( MagEphemInfo->Shell_Vgc );
    LGM_ARRAY_2D_FREE( MagEphemInfo->nMinima );
    LGM_ARRAY_2D_FREE( MagEphemInfo->nMaxima );

    LGM_ARRAY_2D_FREE( MagEphemInfo->nFieldPnts );
    LGM_ARRAY_3D_FREE( MagEphemInfo->s_gsm );
    LGM_ARRAY_3D_FREE( MagEphemInfo->Bmag );
    LGM_ARRAY_3D_FREE( MagEphemInfo->x_gsm );
    LGM_ARRAY_3D_FREE( MagEphemInfo->y_gsm );
    LGM_ARRAY_3D_FREE( MagEphemInfo->z_gsm );

    LGM_ARRAY_1D_FREE( MagEphemInfo->LHilton );
    LGM_ARRAY_1D_FREE( MagEphemInfo->LMcIlwain );
    LGM_ARRAY_1D_FREE( MagEphemInfo->Lstar );
    LGM_ARRAY_1D_FREE( MagEphemInfo->DriftOrbitType );

    FreeLstarInfo( MagEphemInfo->LstarInfo );
}

void Lgm_FreeMagEphemInfo( Lgm_MagEphemInfo  *MagEphemInfo ) {

    /*
     * De-allocate Arrays that depend on # of pitch angles.
     */
    Lgm_FreeMagEphemInfo_Children(MagEphemInfo);
    free( MagEphemInfo );

}



void WriteMagEphemInfoStruct( char *Filename, int nPitchAngles, Lgm_MagEphemInfo *MagEphemInfo ) {

    long int fd;

    fd = open( Filename, O_CREAT|O_WRONLY, 0664 );
    write( fd, &nPitchAngles, sizeof(nPitchAngles) );
    write( fd, MagEphemInfo, sizeof(*MagEphemInfo) );

    write( fd, MagEphemInfo->Alpha,        nPitchAngles*sizeof( double ) );
    write( fd, MagEphemInfo->Pmn_gsm,      nPitchAngles*sizeof( Lgm_Vector ) );
    write( fd, MagEphemInfo->Pms_gsm,      nPitchAngles*sizeof( Lgm_Vector ) );
    write( fd, MagEphemInfo->Bm,           nPitchAngles*sizeof( double ) );
    write( fd, MagEphemInfo->I,            nPitchAngles*sizeof( double ) );
    write( fd, MagEphemInfo->Sb,           nPitchAngles*sizeof( double ) );
    write( fd, MagEphemInfo->Tb,           nPitchAngles*sizeof( double ) );
    write( fd, MagEphemInfo->K,            nPitchAngles*sizeof( double ) );
    write( fd, MagEphemInfo->nShellPoints, nPitchAngles*sizeof( int ) );

    write( fd, &MagEphemInfo->ShellSphericalFootprint_Pn[0][0], nPitchAngles*100*sizeof( Lgm_Vector ) );

    write( fd, &MagEphemInfo->ShellSphericalFootprint_Sn[0][0], nPitchAngles*100*sizeof( double ) );
    write( fd, &MagEphemInfo->ShellSphericalFootprint_Bn[0][0], nPitchAngles*100*sizeof( double ) );
    write( fd, &MagEphemInfo->ShellSphericalFootprint_Ps[0][0], nPitchAngles*100*sizeof( Lgm_Vector ) );
    write( fd, &MagEphemInfo->ShellSphericalFootprint_Ss[0][0], nPitchAngles*100*sizeof( double ) );
    write( fd, &MagEphemInfo->ShellSphericalFootprint_Bs[0][0], nPitchAngles*100*sizeof( double ) );
    write( fd, &MagEphemInfo->ShellEllipsoidFootprint_Ps[0][0], nPitchAngles*100*sizeof( Lgm_Vector ) );
    write( fd, &MagEphemInfo->ShellEllipsoidFootprint_Ss[0][0], nPitchAngles*100*sizeof( double ) );
    write( fd, &MagEphemInfo->ShellEllipsoidFootprint_Bs[0][0], nPitchAngles*100*sizeof( double ) );
    write( fd, &MagEphemInfo->ShellEllipsoidFootprint_Pn[0][0], nPitchAngles*100*sizeof( Lgm_Vector ) );
    write( fd, &MagEphemInfo->ShellEllipsoidFootprint_Sn[0][0], nPitchAngles*100*sizeof( double ) );
    write( fd, &MagEphemInfo->ShellEllipsoidFootprint_Bn[0][0], nPitchAngles*100*sizeof( double ) );
    write( fd, &MagEphemInfo->ShellMirror_Pn[0][0],             nPitchAngles*100*sizeof( Lgm_Vector ) );
    write( fd, &MagEphemInfo->ShellMirror_Sn[0][0],             nPitchAngles*100*sizeof( double ) );
    write( fd, &MagEphemInfo->ShellMirror_Ps[0][0],             nPitchAngles*100*sizeof( Lgm_Vector ) );
    write( fd, &MagEphemInfo->ShellMirror_Ss[0][0],             nPitchAngles*100*sizeof( double ) );
    write( fd, &MagEphemInfo->ShellI[0][0],                     nPitchAngles*100*sizeof( double ) );


    write( fd, &MagEphemInfo->nFieldPnts[0][0], nPitchAngles*48*sizeof( int ) );
    write( fd, &MagEphemInfo->s_gsm[0][0][0],      nPitchAngles*48*1000*sizeof( double ) );
    write( fd, &MagEphemInfo->Bmag[0][0][0],       nPitchAngles*48*1000*sizeof( double ) );
    write( fd, &MagEphemInfo->x_gsm[0][0][0],      nPitchAngles*48*1000*sizeof( double ) );
    write( fd, &MagEphemInfo->y_gsm[0][0][0],      nPitchAngles*48*1000*sizeof( double ) );
    write( fd, &MagEphemInfo->z_gsm[0][0][0],      nPitchAngles*48*1000*sizeof( double ) );

    write( fd, &MagEphemInfo->LHilton[0],   nPitchAngles*sizeof( double ) );
    write( fd, &MagEphemInfo->LMcIlwain[0], nPitchAngles*sizeof( double ) );
    write( fd, &MagEphemInfo->Lstar[0],     nPitchAngles*sizeof( double ) );

    close(fd);


}


void ReadMagEphemInfoStruct( char *Filename, int *nPitchAngles, Lgm_MagEphemInfo *MagEphemInfo ) {

    double      *ddata;
    Lgm_Vector  *vdata;
    int         *idata;
    int      n;
    long int fd;

    fd = open( Filename, O_RDONLY );

    read( fd, &n, sizeof( n ) );
    *nPitchAngles = n;
    read( fd, MagEphemInfo, sizeof(*MagEphemInfo) );

    ddata = (double *)calloc( n, sizeof(double) );
    read( fd, ddata, n*sizeof( double ) );
    LGM_ARRAY_FROM_DATA_1D( MagEphemInfo->Alpha, ddata, n, double );

    vdata = (Lgm_Vector *)calloc( n, sizeof(Lgm_Vector) );
    read( fd, vdata, n*sizeof( Lgm_Vector ) );
    LGM_ARRAY_FROM_DATA_1D( MagEphemInfo->Pmn_gsm, vdata, n, Lgm_Vector );

    vdata = (Lgm_Vector *)calloc( n, sizeof(Lgm_Vector) );
    read( fd, vdata, n*sizeof( Lgm_Vector ) );
    LGM_ARRAY_FROM_DATA_1D( MagEphemInfo->Pms_gsm, vdata, n, Lgm_Vector );

    ddata = (double *)calloc( n, sizeof(double) );
    read( fd, ddata, n*sizeof( double ) );
    LGM_ARRAY_FROM_DATA_1D( MagEphemInfo->Bm, ddata, n, double );

    ddata = (double *)calloc( n, sizeof(double) );
    read( fd, ddata, n*sizeof( double ) );
    LGM_ARRAY_FROM_DATA_1D( MagEphemInfo->I, ddata, n, double );

    ddata = (double *)calloc( n, sizeof(double) );
    read( fd, ddata, n*sizeof( double ) );
    LGM_ARRAY_FROM_DATA_1D( MagEphemInfo->Sb, ddata, n, double );

    ddata = (double *)calloc( n, sizeof(double) );
    read( fd, ddata, n*sizeof( double ) );
    LGM_ARRAY_FROM_DATA_1D( MagEphemInfo->Tb, ddata, n, double );

    ddata = (double *)calloc( n, sizeof(double) );
    read( fd, ddata, n*sizeof( double ) );
    LGM_ARRAY_FROM_DATA_1D( MagEphemInfo->K, ddata, n, double );

    idata = (int *)calloc( n, sizeof(int) );
    read( fd, idata, n*sizeof( int ) );
    LGM_ARRAY_FROM_DATA_1D( MagEphemInfo->nShellPoints, idata, n, int );


    vdata = (Lgm_Vector *)calloc( n*100, sizeof(Lgm_Vector) );
    read( fd, vdata, n*100*sizeof( Lgm_Vector ) );
    LGM_ARRAY_FROM_DATA_2D( MagEphemInfo->ShellSphericalFootprint_Pn, vdata, n, 100, Lgm_Vector );

    ddata = (double *)calloc( n*100, sizeof(double) );
    read( fd, ddata, n*100*sizeof( double ) );
    LGM_ARRAY_FROM_DATA_2D( MagEphemInfo->ShellSphericalFootprint_Sn, ddata, n, 100, double );

    ddata = (double *)calloc( n*100, sizeof(double) );
    read( fd, ddata, n*100*sizeof( double ) );
    LGM_ARRAY_FROM_DATA_2D( MagEphemInfo->ShellSphericalFootprint_Bn, ddata, n, 100, double );

    vdata = (Lgm_Vector *)calloc( n*100, sizeof(Lgm_Vector) );
    read( fd, vdata, n*100*sizeof( Lgm_Vector ) );
    LGM_ARRAY_FROM_DATA_2D( MagEphemInfo->ShellSphericalFootprint_Ps, vdata, n, 100, Lgm_Vector );

    ddata = (double *)calloc( n*100, sizeof(double) );
    read( fd, ddata, n*100*sizeof( double ) );
    LGM_ARRAY_FROM_DATA_2D( MagEphemInfo->ShellSphericalFootprint_Ss, ddata, n, 100, double );

    ddata = (double *)calloc( n*100, sizeof(double) );
    read( fd, ddata, n*100*sizeof( double ) );
    LGM_ARRAY_FROM_DATA_2D( MagEphemInfo->ShellSphericalFootprint_Bs, ddata, n, 100, double );

    vdata = (Lgm_Vector *)calloc( n*100, sizeof(Lgm_Vector) );
    read( fd, vdata, n*100*sizeof( Lgm_Vector ) );
    LGM_ARRAY_FROM_DATA_2D( MagEphemInfo->ShellEllipsoidFootprint_Ps, vdata, n, 100, Lgm_Vector );

    ddata = (double *)calloc( n*100, sizeof(double) );
    read( fd, ddata, n*100*sizeof( double ) );
    LGM_ARRAY_FROM_DATA_2D( MagEphemInfo->ShellEllipsoidFootprint_Ss, ddata, n, 100, double );

    ddata = (double *)calloc( n*100, sizeof(double) );
    read( fd, ddata, n*100*sizeof( double ) );
    LGM_ARRAY_FROM_DATA_2D( MagEphemInfo->ShellEllipsoidFootprint_Bs, ddata, n, 100, double );

    vdata = (Lgm_Vector *)calloc( n*100, sizeof(Lgm_Vector) );
    read( fd, vdata, n*100*sizeof( Lgm_Vector ) );
    LGM_ARRAY_FROM_DATA_2D( MagEphemInfo->ShellEllipsoidFootprint_Pn, vdata, n, 100, Lgm_Vector );

    ddata = (double *)calloc( n*100, sizeof(double) );
    read( fd, ddata, n*100*sizeof( double ) );
    LGM_ARRAY_FROM_DATA_2D( MagEphemInfo->ShellEllipsoidFootprint_Sn, ddata, n, 100, double );

    ddata = (double *)calloc( n*100, sizeof(double) );
    read( fd, ddata, n*100*sizeof( double ) );
    LGM_ARRAY_FROM_DATA_2D( MagEphemInfo->ShellEllipsoidFootprint_Bn, ddata, n, 100, double );

    vdata = (Lgm_Vector *)calloc( n*100, sizeof(Lgm_Vector) );
    read( fd, vdata, n*100*sizeof( Lgm_Vector ) );
    LGM_ARRAY_FROM_DATA_2D( MagEphemInfo->ShellMirror_Pn, vdata, n, 100, Lgm_Vector );

    ddata = (double *)calloc( n*100, sizeof(double) );
    read( fd, ddata, n*100*sizeof( double ) );
    LGM_ARRAY_FROM_DATA_2D( MagEphemInfo->ShellMirror_Sn, ddata, n, 100, double );

    vdata = (Lgm_Vector *)calloc( n*100, sizeof(Lgm_Vector) );
    read( fd, vdata, n*100*sizeof( Lgm_Vector ) );
    LGM_ARRAY_FROM_DATA_2D( MagEphemInfo->ShellMirror_Ps, vdata, n, 100, Lgm_Vector );

    ddata = (double *)calloc( n*100, sizeof(double) );
    read( fd, ddata, n*100*sizeof( double ) );
    LGM_ARRAY_FROM_DATA_2D( MagEphemInfo->ShellMirror_Ss, ddata, n, 100, double );

    ddata = (double *)calloc( n*100, sizeof(double) );
    read( fd, ddata, n*100*sizeof( double ) );
    LGM_ARRAY_FROM_DATA_2D( MagEphemInfo->ShellI, ddata, n, 100, double );



    idata = (int *)calloc( n*48, sizeof(int) );
    read( fd, idata, n*48*sizeof( int ) );
    LGM_ARRAY_FROM_DATA_2D( MagEphemInfo->nFieldPnts, idata, n, 48, int );

    ddata = (double *)calloc( n*48*1000, sizeof(double) );
    read( fd, ddata, n*48*1000*sizeof( double ) );
    LGM_ARRAY_FROM_DATA_3D( MagEphemInfo->s_gsm, ddata, n, 48, 1000, double );

    ddata = (double *)calloc( n*48*1000, sizeof(double) );
    read( fd, ddata, n*48*1000*sizeof( double ) );
    LGM_ARRAY_FROM_DATA_3D( MagEphemInfo->Bmag, ddata, n, 48, 1000, double );

    ddata = (double *)calloc( n*48*1000, sizeof(double) );
    read( fd, ddata, n*48*1000*sizeof( double ) );
    LGM_ARRAY_FROM_DATA_3D( MagEphemInfo->x_gsm, ddata, n, 48, 1000, double );

    ddata = (double *)calloc( n*48*1000, sizeof(double) );
    read( fd, ddata, n*48*1000*sizeof( double ) );
    LGM_ARRAY_FROM_DATA_3D( MagEphemInfo->y_gsm, ddata, n, 48, 1000, double );

    ddata = (double *)calloc( n*48*1000, sizeof(double) );
    read( fd, ddata, n*48*1000*sizeof( double ) );
    LGM_ARRAY_FROM_DATA_3D( MagEphemInfo->z_gsm, ddata, n, 48, 1000, double );



    ddata = (double *)calloc( n, sizeof(double) );
    read( fd, ddata, n*sizeof( double ) );
    LGM_ARRAY_FROM_DATA_1D( MagEphemInfo->LHilton, ddata, n, double );

    ddata = (double *)calloc( n, sizeof(double) );
    read( fd, ddata, n*sizeof( double ) );
    LGM_ARRAY_FROM_DATA_1D( MagEphemInfo->LMcIlwain, ddata, n, double );

    ddata = (double *)calloc( n, sizeof(double) );
    read( fd, ddata,     n*sizeof( double ) );
    LGM_ARRAY_FROM_DATA_1D( MagEphemInfo->Lstar, ddata, n, double );

    close(fd);


}



/*
 * This allocates/initializes memory for a MagEphemData structure.
 */
Lgm_MagEphemData *Lgm_InitMagEphemData( int nRows, int nPA ) {


    Lgm_MagEphemData  *MagEphemData = (Lgm_MagEphemData *) calloc (1, sizeof(*MagEphemData));

    MagEphemData->H5_nPerigee = 0;
    MagEphemData->H5_nApogee  = 0;
    MagEphemData->H5_nAscend  = 0;

    LGM_ARRAY_2D( MagEphemData->H5_Perigee_IsoTimes,  nRows, 80,      char   );
    LGM_ARRAY_2D( MagEphemData->H5_Apogee_IsoTimes,   nRows, 80,      char   );
    LGM_ARRAY_2D( MagEphemData->H5_Ascend_IsoTimes,   nRows, 80,      char   );
    LGM_ARRAY_2D( MagEphemData->H5_Perigee_Geod,      nRows, 3,       double );
    LGM_ARRAY_2D( MagEphemData->H5_Apogee_Geod,       nRows, 3,       double );
    LGM_ARRAY_2D( MagEphemData->H5_Ascend_Geod,       nRows, 3,       double );
                                             
    LGM_ARRAY_1D( MagEphemData->H5_Alpha,             nPA,            double );

    LGM_ARRAY_2D( MagEphemData->H5_IsoTimes,          nRows, 80,      char   );
    LGM_ARRAY_2D( MagEphemData->H5_FieldLineType,     nRows, 80,      char   );
    LGM_ARRAY_2D( MagEphemData->H5_IntModel,          nRows, 80,      char   );
    LGM_ARRAY_2D( MagEphemData->H5_ExtModel,          nRows, 80,      char   );

    LGM_ARRAY_1D( MagEphemData->H5_Date,              nRows,          long int );
    LGM_ARRAY_1D( MagEphemData->H5_Doy,               nRows,          int    );
    LGM_ARRAY_1D( MagEphemData->H5_UTC,               nRows,          double );
    LGM_ARRAY_1D( MagEphemData->H5_JD,                nRows,          double );
    LGM_ARRAY_1D( MagEphemData->H5_GpsTime,           nRows,          double );
    LGM_ARRAY_1D( MagEphemData->H5_TiltAngle,         nRows,          double );
    LGM_ARRAY_1D( MagEphemData->H5_InOut,             nRows,          int    );
    LGM_ARRAY_1D( MagEphemData->H5_OrbitNumber,       nRows,          int    );


    LGM_ARRAY_2D( MagEphemData->H5_Rgeo,              nRows, 3,       double );
    LGM_ARRAY_2D( MagEphemData->H5_Rgeod,             nRows, 3,       double );
    LGM_ARRAY_2D( MagEphemData->H5_Rgeod_LatLon,      nRows, 2,       double );
    LGM_ARRAY_1D( MagEphemData->H5_Rgeod_Height,      nRows,          double );
    LGM_ARRAY_2D( MagEphemData->H5_Rgsm,              nRows, 3,       double );
    LGM_ARRAY_2D( MagEphemData->H5_Rsm,               nRows, 3,       double );
    LGM_ARRAY_2D( MagEphemData->H5_Rgei,              nRows, 3,       double );
    LGM_ARRAY_2D( MagEphemData->H5_Rgse,              nRows, 3,       double );

    LGM_ARRAY_1D( MagEphemData->H5_CDMAG_MLAT,        nRows,          double );
    LGM_ARRAY_1D( MagEphemData->H5_CDMAG_MLON,        nRows,          double );
    LGM_ARRAY_1D( MagEphemData->H5_CDMAG_MLT,         nRows,          double );
    LGM_ARRAY_1D( MagEphemData->H5_CDMAG_R,           nRows,          double );

    LGM_ARRAY_1D( MagEphemData->H5_EDMAG_MLAT,        nRows,          double );
    LGM_ARRAY_1D( MagEphemData->H5_EDMAG_MLON,        nRows,          double );
    LGM_ARRAY_1D( MagEphemData->H5_EDMAG_MLT,         nRows,          double );
    LGM_ARRAY_1D( MagEphemData->H5_EDMAG_R,           nRows,          double );

    LGM_ARRAY_1D( MagEphemData->H5_Kp,                nRows,          double );
    LGM_ARRAY_1D( MagEphemData->H5_Dst,               nRows,          double );

    LGM_ARRAY_2D( MagEphemData->H5_Bsc_gsm,           nRows, 4,       double );

    LGM_ARRAY_1D( MagEphemData->H5_S_sc_to_pfn,       nRows,          double );
    LGM_ARRAY_1D( MagEphemData->H5_S_sc_to_pfs,       nRows,          double );
    LGM_ARRAY_1D( MagEphemData->H5_S_pfs_to_Bmin,     nRows,          double );
    LGM_ARRAY_1D( MagEphemData->H5_S_Bmin_to_sc,      nRows,          double );
    LGM_ARRAY_1D( MagEphemData->H5_S_total,           nRows,          double );
    LGM_ARRAY_1D( MagEphemData->H5_d2B_ds2,           nRows,          double );
    LGM_ARRAY_1D( MagEphemData->H5_Sb0,               nRows,          double );
    LGM_ARRAY_1D( MagEphemData->H5_RadiusOfCurv,      nRows,          double );


    LGM_ARRAY_2D( MagEphemData->H5_Pfn_geo,           nRows, 3,       double );
    LGM_ARRAY_2D( MagEphemData->H5_Pfn_gsm,           nRows, 3,       double );
    LGM_ARRAY_2D( MagEphemData->H5_Pfn_geod,          nRows, 3,       double );
    LGM_ARRAY_2D( MagEphemData->H5_Pfn_geod_LatLon,   nRows, 2,       double );
    LGM_ARRAY_1D( MagEphemData->H5_Pfn_geod_Height,   nRows,          double );
    LGM_ARRAY_2D( MagEphemData->H5_Pfn_cdmag,         nRows, 3,       double );
    LGM_ARRAY_1D( MagEphemData->H5_Pfn_CD_MLAT,       nRows,          double );
    LGM_ARRAY_1D( MagEphemData->H5_Pfn_CD_MLON,       nRows,          double );
    LGM_ARRAY_1D( MagEphemData->H5_Pfn_CD_MLT,        nRows,          double );
    LGM_ARRAY_2D( MagEphemData->H5_Pfn_edmag,         nRows, 3,       double );
    LGM_ARRAY_1D( MagEphemData->H5_Pfn_ED_MLAT,       nRows,          double );
    LGM_ARRAY_1D( MagEphemData->H5_Pfn_ED_MLON,       nRows,          double );
    LGM_ARRAY_1D( MagEphemData->H5_Pfn_ED_MLT,        nRows,          double );
    LGM_ARRAY_2D( MagEphemData->H5_Bfn_geo,           nRows, 4,       double );
    LGM_ARRAY_2D( MagEphemData->H5_Bfn_gsm,           nRows, 4,       double );
    LGM_ARRAY_1D( MagEphemData->H5_LossConeAngleN,    nRows,          double );

    LGM_ARRAY_2D( MagEphemData->H5_Pfs_geo,           nRows, 3,       double );
    LGM_ARRAY_2D( MagEphemData->H5_Pfs_gsm,           nRows, 3,       double );
    LGM_ARRAY_2D( MagEphemData->H5_Pfs_geod,          nRows, 3,       double );
    LGM_ARRAY_2D( MagEphemData->H5_Pfs_geod_LatLon,   nRows, 2,       double );
    LGM_ARRAY_1D( MagEphemData->H5_Pfs_geod_Height,   nRows,          double );
    LGM_ARRAY_2D( MagEphemData->H5_Pfs_cdmag,         nRows, 3,       double );
    LGM_ARRAY_1D( MagEphemData->H5_Pfs_CD_MLAT,       nRows,          double );
    LGM_ARRAY_1D( MagEphemData->H5_Pfs_CD_MLON,       nRows,          double );
    LGM_ARRAY_1D( MagEphemData->H5_Pfs_CD_MLT,        nRows,          double );
    LGM_ARRAY_2D( MagEphemData->H5_Pfs_edmag,         nRows, 3,       double );
    LGM_ARRAY_1D( MagEphemData->H5_Pfs_ED_MLAT,       nRows,          double );
    LGM_ARRAY_1D( MagEphemData->H5_Pfs_ED_MLON,       nRows,          double );
    LGM_ARRAY_1D( MagEphemData->H5_Pfs_ED_MLT,        nRows,          double );
    LGM_ARRAY_2D( MagEphemData->H5_Bfs_geo,           nRows, 4,       double );
    LGM_ARRAY_2D( MagEphemData->H5_Bfs_gsm,           nRows, 4,       double );
    LGM_ARRAY_1D( MagEphemData->H5_LossConeAngleS,    nRows,          double );


    LGM_ARRAY_2D( MagEphemData->H5_Pmin_gsm,          nRows, 3,       double );
    LGM_ARRAY_2D( MagEphemData->H5_Bmin_gsm,          nRows, 4,       double );

    LGM_ARRAY_1D( MagEphemData->H5_Lsimple,           nRows,          double );
    LGM_ARRAY_1D( MagEphemData->H5_InvLat,            nRows,          double );
    LGM_ARRAY_1D( MagEphemData->H5_Lm_eq,             nRows,          double );
    LGM_ARRAY_1D( MagEphemData->H5_InvLat_eq,         nRows,          double );
    LGM_ARRAY_1D( MagEphemData->H5_BoverBeq,          nRows,          double );
    LGM_ARRAY_1D( MagEphemData->H5_MlatFromBoverBeq,  nRows,          double );

    LGM_ARRAY_1D( MagEphemData->H5_M_used,            nRows,          double );
    LGM_ARRAY_1D( MagEphemData->H5_M_ref,             nRows,          double );
    LGM_ARRAY_1D( MagEphemData->H5_M_igrf,            nRows,          double );

    LGM_ARRAY_2D( MagEphemData->H5_Lstar,             nRows, nPA,     double );
    LGM_ARRAY_2D( MagEphemData->H5_Sb,                nRows, nPA,     double );
    LGM_ARRAY_2D( MagEphemData->H5_Tb,                nRows, nPA,     double );
    LGM_ARRAY_2D( MagEphemData->H5_Kappa,             nRows, nPA,     double );
    LGM_ARRAY_2D( MagEphemData->H5_DriftShellType,    nRows, nPA,     int    );
    LGM_ARRAY_2D( MagEphemData->H5_L,                 nRows, nPA,     double );
    LGM_ARRAY_2D( MagEphemData->H5_Bm,                nRows, nPA,     double );
    LGM_ARRAY_2D( MagEphemData->H5_I,                 nRows, nPA,     double );
    LGM_ARRAY_2D( MagEphemData->H5_K,                 nRows, nPA,     double );

    return( MagEphemData );

}



void Lgm_FreeMagEphemData( Lgm_MagEphemData *MagEphemData ) {

    LGM_ARRAY_2D_FREE( MagEphemData->H5_Perigee_IsoTimes );
    LGM_ARRAY_2D_FREE( MagEphemData->H5_Apogee_IsoTimes );
    LGM_ARRAY_2D_FREE( MagEphemData->H5_Ascend_IsoTimes );
    LGM_ARRAY_2D_FREE( MagEphemData->H5_Perigee_Geod );
    LGM_ARRAY_2D_FREE( MagEphemData->H5_Apogee_Geod );
    LGM_ARRAY_2D_FREE( MagEphemData->H5_Ascend_Geod );

    LGM_ARRAY_1D_FREE( MagEphemData->H5_Alpha );

    LGM_ARRAY_2D_FREE( MagEphemData->H5_IsoTimes );
    LGM_ARRAY_2D_FREE( MagEphemData->H5_FieldLineType );
    LGM_ARRAY_2D_FREE( MagEphemData->H5_IntModel );
    LGM_ARRAY_2D_FREE( MagEphemData->H5_ExtModel );

    LGM_ARRAY_1D_FREE( MagEphemData->H5_Date );
    LGM_ARRAY_1D_FREE( MagEphemData->H5_Doy );
    LGM_ARRAY_1D_FREE( MagEphemData->H5_UTC );
    LGM_ARRAY_1D_FREE( MagEphemData->H5_JD );
    LGM_ARRAY_1D_FREE( MagEphemData->H5_GpsTime );
    LGM_ARRAY_1D_FREE( MagEphemData->H5_TiltAngle );
    LGM_ARRAY_1D_FREE( MagEphemData->H5_InOut );
    LGM_ARRAY_1D_FREE( MagEphemData->H5_OrbitNumber );


    LGM_ARRAY_2D_FREE( MagEphemData->H5_Rgeo );
    LGM_ARRAY_2D_FREE( MagEphemData->H5_Rgeod );
    LGM_ARRAY_2D_FREE( MagEphemData->H5_Rgeod_LatLon );
    LGM_ARRAY_1D_FREE( MagEphemData->H5_Rgeod_Height );
    LGM_ARRAY_2D_FREE( MagEphemData->H5_Rgsm );
    LGM_ARRAY_2D_FREE( MagEphemData->H5_Rsm );
    LGM_ARRAY_2D_FREE( MagEphemData->H5_Rgei );
    LGM_ARRAY_2D_FREE( MagEphemData->H5_Rgse );

    LGM_ARRAY_1D_FREE( MagEphemData->H5_CDMAG_MLAT );
    LGM_ARRAY_1D_FREE( MagEphemData->H5_CDMAG_MLON );
    LGM_ARRAY_1D_FREE( MagEphemData->H5_CDMAG_MLT );
    LGM_ARRAY_1D_FREE( MagEphemData->H5_CDMAG_R );

    LGM_ARRAY_1D_FREE( MagEphemData->H5_EDMAG_MLAT );
    LGM_ARRAY_1D_FREE( MagEphemData->H5_EDMAG_MLON );
    LGM_ARRAY_1D_FREE( MagEphemData->H5_EDMAG_MLT );
    LGM_ARRAY_1D_FREE( MagEphemData->H5_EDMAG_R );

    LGM_ARRAY_1D_FREE( MagEphemData->H5_Kp );
    LGM_ARRAY_1D_FREE( MagEphemData->H5_Dst );

    LGM_ARRAY_2D_FREE( MagEphemData->H5_Bsc_gsm );

    LGM_ARRAY_1D_FREE( MagEphemData->H5_S_sc_to_pfn );
    LGM_ARRAY_1D_FREE( MagEphemData->H5_S_sc_to_pfs );
    LGM_ARRAY_1D_FREE( MagEphemData->H5_S_pfs_to_Bmin );
    LGM_ARRAY_1D_FREE( MagEphemData->H5_S_Bmin_to_sc );
    LGM_ARRAY_1D_FREE( MagEphemData->H5_S_total );
    LGM_ARRAY_1D_FREE( MagEphemData->H5_d2B_ds2 );
    LGM_ARRAY_1D_FREE( MagEphemData->H5_Sb0 );
    LGM_ARRAY_1D_FREE( MagEphemData->H5_RadiusOfCurv );


    LGM_ARRAY_2D_FREE( MagEphemData->H5_Pfn_geo );
    LGM_ARRAY_2D_FREE( MagEphemData->H5_Pfn_gsm );
    LGM_ARRAY_2D_FREE( MagEphemData->H5_Pfn_geod );
    LGM_ARRAY_2D_FREE( MagEphemData->H5_Pfn_geod_LatLon );
    LGM_ARRAY_1D_FREE( MagEphemData->H5_Pfn_geod_Height );
    LGM_ARRAY_2D_FREE( MagEphemData->H5_Pfn_cdmag );
    LGM_ARRAY_1D_FREE( MagEphemData->H5_Pfn_CD_MLAT );
    LGM_ARRAY_1D_FREE( MagEphemData->H5_Pfn_CD_MLON );
    LGM_ARRAY_1D_FREE( MagEphemData->H5_Pfn_CD_MLT );
    LGM_ARRAY_2D_FREE( MagEphemData->H5_Pfn_edmag );
    LGM_ARRAY_1D_FREE( MagEphemData->H5_Pfn_ED_MLAT );
    LGM_ARRAY_1D_FREE( MagEphemData->H5_Pfn_ED_MLON );
    LGM_ARRAY_1D_FREE( MagEphemData->H5_Pfn_ED_MLT );
    LGM_ARRAY_2D_FREE( MagEphemData->H5_Bfn_geo );
    LGM_ARRAY_2D_FREE( MagEphemData->H5_Bfn_gsm );
    LGM_ARRAY_1D_FREE( MagEphemData->H5_LossConeAngleN );

    LGM_ARRAY_2D_FREE( MagEphemData->H5_Pfs_geo );
    LGM_ARRAY_2D_FREE( MagEphemData->H5_Pfs_gsm );
    LGM_ARRAY_2D_FREE( MagEphemData->H5_Pfs_geod );
    LGM_ARRAY_2D_FREE( MagEphemData->H5_Pfs_geod_LatLon );
    LGM_ARRAY_1D_FREE( MagEphemData->H5_Pfs_geod_Height );
    LGM_ARRAY_2D_FREE( MagEphemData->H5_Pfs_cdmag );
    LGM_ARRAY_1D_FREE( MagEphemData->H5_Pfs_CD_MLAT );
    LGM_ARRAY_1D_FREE( MagEphemData->H5_Pfs_CD_MLON );
    LGM_ARRAY_1D_FREE( MagEphemData->H5_Pfs_CD_MLT );
    LGM_ARRAY_2D_FREE( MagEphemData->H5_Pfs_edmag );
    LGM_ARRAY_1D_FREE( MagEphemData->H5_Pfs_ED_MLAT );
    LGM_ARRAY_1D_FREE( MagEphemData->H5_Pfs_ED_MLON );
    LGM_ARRAY_1D_FREE( MagEphemData->H5_Pfs_ED_MLT );
    LGM_ARRAY_2D_FREE( MagEphemData->H5_Bfs_geo );
    LGM_ARRAY_2D_FREE( MagEphemData->H5_Bfs_gsm );
    LGM_ARRAY_1D_FREE( MagEphemData->H5_LossConeAngleS );


    LGM_ARRAY_2D_FREE( MagEphemData->H5_Pmin_gsm );
    LGM_ARRAY_2D_FREE( MagEphemData->H5_Bmin_gsm );

    LGM_ARRAY_1D_FREE( MagEphemData->H5_Lsimple );
    LGM_ARRAY_1D_FREE( MagEphemData->H5_InvLat );
    LGM_ARRAY_1D_FREE( MagEphemData->H5_Lm_eq );
    LGM_ARRAY_1D_FREE( MagEphemData->H5_InvLat_eq );
    LGM_ARRAY_1D_FREE( MagEphemData->H5_BoverBeq );
    LGM_ARRAY_1D_FREE( MagEphemData->H5_MlatFromBoverBeq );

    LGM_ARRAY_1D_FREE( MagEphemData->H5_M_used );
    LGM_ARRAY_1D_FREE( MagEphemData->H5_M_ref );
    LGM_ARRAY_1D_FREE( MagEphemData->H5_M_igrf );

    LGM_ARRAY_2D_FREE( MagEphemData->H5_Lstar );
    LGM_ARRAY_2D_FREE( MagEphemData->H5_Sb );
    LGM_ARRAY_2D_FREE( MagEphemData->H5_Tb );
    LGM_ARRAY_2D_FREE( MagEphemData->H5_Kappa );
    LGM_ARRAY_2D_FREE( MagEphemData->H5_DriftShellType );
    LGM_ARRAY_2D_FREE( MagEphemData->H5_L );
    LGM_ARRAY_2D_FREE( MagEphemData->H5_Bm );
    LGM_ARRAY_2D_FREE( MagEphemData->H5_I );
    LGM_ARRAY_2D_FREE( MagEphemData->H5_K );

    free( MagEphemData );

    return;

}



