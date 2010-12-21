#include "Lgm/Lgm_MagEphemInfo.h"
#include "Lgm/Lgm_DynamicMemory.h"


/*
 *  A convenience routine to set some default values
 */

Lgm_MagEphemInfo *Lgm_InitMagEphemInfo( int Verbosity, int MaxPitchAngles ) {

    Lgm_MagEphemInfo  *MagEphemInfo = (Lgm_MagEphemInfo *) calloc (1, sizeof(*MagEphemInfo));
    MagEphemInfo->LstarInfo = InitLstarInfo( Verbosity );


    /*
     * Allocate Arrays that depend on # of pitch angles.
     */
    LGM_ARRAY_1D( MagEphemInfo->Alpha,        MaxPitchAngles, double );
    LGM_ARRAY_1D( MagEphemInfo->Pmn_gsm,      MaxPitchAngles, Lgm_Vector );
    LGM_ARRAY_1D( MagEphemInfo->Pms_gsm,      MaxPitchAngles, Lgm_Vector );
    LGM_ARRAY_1D( MagEphemInfo->Bm,           MaxPitchAngles, double );
    LGM_ARRAY_1D( MagEphemInfo->I,            MaxPitchAngles, double );
    LGM_ARRAY_1D( MagEphemInfo->Sb,           MaxPitchAngles, double );
    LGM_ARRAY_1D( MagEphemInfo->Tb,           MaxPitchAngles, double );
    LGM_ARRAY_1D( MagEphemInfo->K,            MaxPitchAngles, double );
    LGM_ARRAY_1D( MagEphemInfo->nShellPoints, MaxPitchAngles, int );

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

    LGM_ARRAY_2D( MagEphemInfo->nFieldPnts, MaxPitchAngles, 48, int );
    LGM_ARRAY_3D( MagEphemInfo->s_gsm,      MaxPitchAngles, 48, 1000, double );
    LGM_ARRAY_3D( MagEphemInfo->Bmag,       MaxPitchAngles, 48, 1000, double );
    LGM_ARRAY_3D( MagEphemInfo->x_gsm,      MaxPitchAngles, 48, 1000, double );
    LGM_ARRAY_3D( MagEphemInfo->y_gsm,      MaxPitchAngles, 48, 1000, double );
    LGM_ARRAY_3D( MagEphemInfo->z_gsm,      MaxPitchAngles, 48, 1000, double );

    LGM_ARRAY_1D( MagEphemInfo->LHilton,   MaxPitchAngles, double );
    LGM_ARRAY_1D( MagEphemInfo->LMcIlwain, MaxPitchAngles, double );
    LGM_ARRAY_1D( MagEphemInfo->Lstar,     MaxPitchAngles, double );



    



    return MagEphemInfo;

}



void Lgm_FreeMagEphemInfo( Lgm_MagEphemInfo  *MagEphemInfo ) {

    /*
     * De-allocate Arrays that depend on # of pitch angles.
     */
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

    LGM_ARRAY_2D_FREE( MagEphemInfo->nFieldPnts );
    LGM_ARRAY_3D_FREE( MagEphemInfo->s_gsm );
    LGM_ARRAY_3D_FREE( MagEphemInfo->Bmag );
    LGM_ARRAY_3D_FREE( MagEphemInfo->x_gsm );
    LGM_ARRAY_3D_FREE( MagEphemInfo->y_gsm );
    LGM_ARRAY_3D_FREE( MagEphemInfo->z_gsm );

    LGM_ARRAY_1D_FREE( MagEphemInfo->LHilton );
    LGM_ARRAY_1D_FREE( MagEphemInfo->LMcIlwain );
    LGM_ARRAY_1D_FREE( MagEphemInfo->Lstar );

    FreeLstarInfo( MagEphemInfo->LstarInfo );
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
 *   $Id$
 */

