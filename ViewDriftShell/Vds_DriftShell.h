#ifndef VDS_DRIFT_SHELL_H
#define VDS_DRIFT_SHELL_H

#include <Lgm_MagEphemInfo.h>

#define GSL_INTERP  gsl_interp_akima


typedef struct Vds_ObjectInfo {


    Lgm_MagEphemInfo    *MagEphemInfo;          //!< Structure containing magnetic ephemeris information.
    int                  MagEphemAlloced;       //!< Flag indicating whether memory for MagEphemInfo has been allocated or not.


    int                  nPitchAngles;          //!< Number of pitch angles defined in MagEphemInfo structure.


    int                  nPnts[90][24];         //!< Number of points in the full (foot-to-foot) field line arrays.
    double               s_gsm[90][24][1000];   //!< Distance along FL for full (foot-to-foot) field line arrays.
    double               x_gsm[90][24][1000];   //!< X-position along FL for full (foot-to-foot) field line arrays.
    double               y_gsm[90][24][1000];   //!< Y-Position along FL for full (foot-to-foot) field line arrays.
    double               z_gsm[90][24][1000];   //!< Z-Position along FL for full (foot-to-foot) field line arrays.



    int                  nPnts2[90][24];        //!< Number of points in the partial (mirror-to-mirror) field line arrays.
    double               s2_gsm[90][24][1000];  //!< Distance along FL for full (mirror-to-mirror) field line arrays.
    double               x2_gsm[90][24][1000];  //!< X-position along FL for full (mirror-to-mirror) field line arrays.
    double               y2_gsm[90][24][1000];  //!< Y-Position along FL for full (mirror-to-mirror) field line arrays.
    double               z2_gsm[90][24][1000];  //!< Z-Position along FL for full (mirror-to-mirror) field line arrays.


    int                  nFieldPoints[90];
    double               x3_gsm[90][24][200];   //!< ???? Vertices?
    double               y3_gsm[90][24][200];   //!< ???? Vertices?
    double               z3_gsm[90][24][200];   //!< ???? Vertices?
    double               nx3_gsm[90][24][200];  //!< ???? Normals?
    double               ny3_gsm[90][24][200];  //!< ???? Normals?
    double               nz3_gsm[90][24][200];  //!< ???? Normals?

    int                  nShellPoints4;
    double               x4_gsm[90][500][200];   //!< ???? Vertices?
    double               y4_gsm[90][500][200];   //!< ???? Vertices?
    double               z4_gsm[90][500][200];   //!< ???? Vertices?
    double               nx4_gsm[90][500][200];  //!< ???? Normals?
    double               ny4_gsm[90][500][200];  //!< ???? Normals?
    double               nz4_gsm[90][500][200];  //!< ???? Normals?


    GLuint               DriftShellList2;        //!< OpenGL Display List for Drift Shell Field Lines (Full)
    GLuint               DriftShellList3;        //!< OpenGL Display List for Drift Shell Field Lines (Mirror-to-Miroor)
    GLuint               DriftShellList4;        //!< OpenGL Display List for Drift Shell Surfaces


} Vds_ObjectInfo;



Vds_ObjectInfo *Vds_InitObjectInfo();
void Vds_FreeObjectInfo( Vds_ObjectInfo *ObjInfo );
void CreateFieldLinesAndDriftShells( char *Filename, Vds_ObjectInfo *ObjInfo );
void MakeFieldLines( int nNewPnts, Vds_ObjectInfo *ObjInfo );
void MakeDriftShellMesh( Vds_ObjectInfo *ObjInfo );
void InterpFieldLine( double *x, double *y, double *z, double *s, double Ss, double Sn, int n1, double *xout, double *yout, double *zout, int n2  );

void ReadMagEphemData( char *Filename, long int Date, double UT, double Lat, double Lon, double Rad, int Mode );

void GenerateDriftShellLists( Vds_ObjectInfo *ObjInfo );
void ReGenerateDriftShellLists( Vds_ObjectInfo *ObjInfo );
void GenerateFieldLineLists( Vds_ObjectInfo *ObjInfo );
void ReGenerateFieldLineLists( Vds_ObjectInfo *ObjInfo );




#endif
