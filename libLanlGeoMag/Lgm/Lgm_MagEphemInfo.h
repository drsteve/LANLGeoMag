#ifndef MAGEPHEMINFO_H
#define MAGEPHEMINFO_H

#include "Lgm/Lgm_MagModelInfo.h"
#include "Lgm/Lgm_LstarInfo.h"
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include "Lgm/Lgm_HDF5.h"


#define MAX_PITCH_ANGLES 90



/*! \struct Lgm_MagEphemInfo
 * Structure to hold all relevant magnetic ephemeris quantities.  I.e.,
 * everything from physical location of S/C to magnetic coordinates adiabatic
 * invariants etc...
 */
typedef struct Lgm_MagEphemInfo {

    Lgm_LstarInfo   *LstarInfo;
    int             LstarQuality;
    int             SaveShellLines;

    long int        Date;           //!< Date in YYYYMMDD format
    double          UTC;            //!< UTC in decimal hours

    int             InOut;          //!< Flag indicating whether we are inbound (-1) or outbound (+1)
    int             OrbitNumber;    //!< Orbit Number

    double          Lat;            //!< Geographic Latitude
    double          Lon;            //!< Geographic Longitude
    double          Rad;            //!< Geographic Radius

    Lgm_Vector      P;              //!< S/C position in GSM
    double          S;              //!< Distance along FL from southern footpoint to S/C location in Re.
    double          B;              //!< Local (model) B-field magnitude (i.e. at S/C position)
    double          d2B_ds2;        //!< Value of second deriviative of B wrt s at eq (used to compute Sb0).
    double          RofC;           //!< Field line radius of curvature at Bmin point.
    double          Sb0;            //!< Value of Sb integral for eq mirroring particles (its not generally zero).

    Lgm_Vector      Pmin;           //!< position of minimum |B| in GSM
    double          Bmin;           //!< Value of |Bmin|
    double          Smin;           //!< Distance from southern footpoint to Pmin along FL.
    double          Snorth;         //!< Distance from initial point to northern footpoint along FL.
    double          Ssouth;         //!< Distance from initial point to southern footpoint along FL.

    Lgm_Vector      Spherical_Footprint_Pn;   //!< position of northern footpoint (at 120km above spherical Earth)
    double          Spherical_Footprint_Sn;   //!< Distance along FL from southern foorpoint in Re
    double          Spherical_Footprint_Bn;   //!< Value of |B| at Spherical_Footprint_Pn

    Lgm_Vector      Spherical_Footprint_Ps;   //!< position of southern footpoint (at 120km above spherical Earth)
    double          Spherical_Footprint_Ss;   //!< Distance along FL from southern foorpoint in Re (i.e. this one is zero by definition)
    double          Spherical_Footprint_Bs;   //!< Value of |B| at Spherical_Footprint_Ps

    Lgm_Vector      Ellipsoid_Footprint_Pn;   //!< position of northern footpoint (at 120km above WGS84 ellipsoid)
    double          Ellipsoid_Footprint_Sn;   //!< Distance along FL from southern foorpoint in Re
    double          Ellipsoid_Footprint_Bn;   //!< Value of |B| at Ellipsoid_Footprint_Pn

    Lgm_Vector      Ellipsoid_Footprint_Ps;   //!< position of southern footpoint (at 120km above WGS84 ellipsoid)
    double          Ellipsoid_Footprint_Ss;   //!< Distance along FL from southern foorpoint in Re (i.e. this one is zero by definition)
    double          Ellipsoid_Footprint_Bs;   //!< Value of |B| at Ellipsoid_Footprint_Ps

    int             FieldLineType;  //!< Field line type. (I.e., LGM_OPEN_IMF, LGM_CLOSED, LGM_OPEN_N_LOBE, LGM_OPEN_S_LOBE, LGM_INSIDE_EARTH, LGM_TARGET_HEIGHT_UNREACHABLE)

    int         UseInterpRoutines;  //!< Whether to use the interped integrands or not.
    int         ComputeVgc;         //!< Whether to compute the Grad-I vals at Pmin for each shell FL

    /*
     * Array of PAs we need to compute things for...
     */
    int         nAlpha;     //!< Number of Pitch Angles
    double      *Alpha;     //!< Pitch Angles (Degrees).
                            //!< 1D array (dynamically allocated). Access as follows: Alpha[ PitchAngleIndex ]
    Lgm_Vector  *Pmn_gsm;   //!< 1D array (dynamically allocated). Access as follows: Pmn_gsm[ PitchAngleIndex ]  position of northern |Bmirror|
    Lgm_Vector  *Pms_gsm;   //!< 1D array (dynamically allocated). Access as follows: Pms_gsm[ PitchAngleIndex ]  position of southern |Bmirror|
    double      *Bm;        //!< 1D array (dynamically allocated). Access as follows: Bm[ PitchAngleIndex ]       Value of |Bm| (i.e. Bmirror)
    double      *I;         //!< 1D array (dynamically allocated). Access as follows: I[ PitchAngleIndex ]
    double      *Sb;        //!< 1D array (dynamically allocated). Access as follows: Sb[ PitchAngleIndex ]
    double      *Tb;        //!< 1D array (dynamically allocated). Access as follows: Tb[ PitchAngleIndex ]       we calculate this for a 1Mev electron but could do more...
    double      *K;         //!< 1D array (dynamically allocated). Access as follows: K[ PitchAngleIndex ]

    int         *nShellPoints; //!<  1D array. Access as follows: nShellPoints[ PitchAngleIndex ] # of point (i.e. FLs) in a shell calculation)

    /*
     * Footpoints XXXkm above spherical Earth defined as sphere with
     * radius==WGS84_A (i.e. equatorial radius of WGS84 model)
     */
    Lgm_Vector  **ShellSphericalFootprint_Pn; //!< [ PitchAngleIndex ][ FieldLineIndex ]  north footprints of shell lines
    double      **ShellSphericalFootprint_Sn; //!< [ PitchAngleIndex ][ FieldLineIndex ]  north mirror locations (dist along FL)
    double      **ShellSphericalFootprint_Bn; //!< [ PitchAngleIndex ][ FieldLineIndex ]  north mirror locations (dist along FL)

    Lgm_Vector  **ShellSphericalFootprint_Ps; //!< [ PitchAngleIndex ][ FieldLineIndex ]  north footprints of shell lines
    double      **ShellSphericalFootprint_Ss; //!< [ PitchAngleIndex ][ FieldLineIndex ]  north mirror locations (dist along FL)
    double      **ShellSphericalFootprint_Bs; //!< [ PitchAngleIndex ][ FieldLineIndex ]  north mirror locations (dist along FL)

    /*
     * footpoints at XXXkm above surface of ellipsoid.
     */
    Lgm_Vector  **ShellEllipsoidFootprint_Ps; //!< [ PitchAngleIndex ][ FieldLineIndex ]  north footprints of shell lines
    double      **ShellEllipsoidFootprint_Ss; //!< [ PitchAngleIndex ][ FieldLineIndex ]  north mirror locations (dist along FL)
    double      **ShellEllipsoidFootprint_Bs; //!< [ PitchAngleIndex ][ FieldLineIndex ]  north mirror locations (dist along FL)

    Lgm_Vector  **ShellEllipsoidFootprint_Pn; //!< [ PitchAngleIndex ][ FieldLineIndex ]  north footprints of shell lines
    double      **ShellEllipsoidFootprint_Sn; //!< [ PitchAngleIndex ][ FieldLineIndex ]  north mirror locations (dist along FL)
    double      **ShellEllipsoidFootprint_Bn; //!< [ PitchAngleIndex ][ FieldLineIndex ]  north mirror locations (dist along FL)


    Lgm_Vector  **ShellMirror_Pn; //!< [ PitchAngleIndex ][ FieldLineIndex ]   north mirror locations
    double      **ShellMirror_Sn; //!< [ PitchAngleIndex ][ FieldLineIndex ]   north mirror locations (dist along FL)

    Lgm_Vector  **ShellMirror_Ps; //!< [ PitchAngleIndex ][ FieldLineIndex ]   south mirror locations
    double      **ShellMirror_Ss; //!< [ PitchAngleIndex ][ FieldLineIndex ]   south mirror locations (dist along FL)

    Lgm_Vector  **Shell_Bmin;     //!< [ PitchAngleIndex ][ FieldLineIndex ]   value of Bmin
    Lgm_Vector  **Shell_Pmin;     //!< [ PitchAngleIndex ][ FieldLineIndex ]   Bmin points
    Lgm_Vector  **Shell_GradI;    //!< [ PitchAngleIndex ][ FieldLineIndex ]   Gradient and Curvature Drift velocity at Pmin
    Lgm_Vector  **Shell_Vgc;      //!< [ PitchAngleIndex ][ FieldLineIndex ]   Gradient of I at Pmin

    double      **ShellI;         //!< [ PitchAngleIndex ][ FieldLineIndex ]   Individual I values computed for each FL

    /*
     * these are like the variables in the LstarInfo structure. Except they have an extra
     * dimension to hold pitch angle as well.
     */
    int         **nFieldPnts;   //!< [ PitchAngleIndex ][ FieldLineIndex ]
    double      ***s_gsm;       //!< [ PitchAngleIndex ][ FieldLineIndex ][ FieldLinePointIndex ]
    double      ***Bmag;        //!< [ PitchAngleIndex ][ FieldLineIndex ][ FieldLinePointIndex ]
    double      ***x_gsm;       //!< [ PitchAngleIndex ][ FieldLineIndex ][ FieldLinePointIndex ]
    double      ***y_gsm;       //!< [ PitchAngleIndex ][ FieldLineIndex ][ FieldLinePointIndex ]
    double      ***z_gsm;       //!< [ PitchAngleIndex ][ FieldLineIndex ][ FieldLinePointIndex ]

    double      Mcurr;
    double      Mref;
    double      Mused;
    double      *LHilton;    //!< LHilton[ PitchAngleIndex ]
    double      *LMcIlwain;  //!< LMcIlwain[ PitchAngleIndex ]
    double      *Lstar;      //!< Lstar[ PitchAngleIndex ]

    int         *DriftOrbitType;     // e.g. Open, Closed, Shabansky
    int         **nMinima;           // # of minima on FL
    int         **nMaxima;           // # of maxima on FL (not including endpoints


} Lgm_MagEphemInfo;

/*! \struct Lgm_MagEphemData
 * Structure to hold all relevant magnetic ephemeris quantities.  I.e.,
 * everything from physical location of S/C to magnetic coordinates adiabatic
 * invariants etc...
 * Similar to Lgm_MagEphemInfo, but organized differently. This is easier to write to HDF/CDF files etc...
 */
typedef struct Lgm_MagEphemData {

    int         H5_nApogee;
    int         H5_nPerigee;
    int         H5_nAscend;
    char        **H5_Perigee_IsoTimes;
    char        **H5_Apogee_IsoTimes;
    char        **H5_Ascend_IsoTimes;
    double      **H5_Perigee_Geod;
    double      **H5_Apogee_Geod;
    double      **H5_Ascend_Geod;

    int         H5_nT;
    int         H5_nAlpha;

    double      *H5_Alpha;

    char        **H5_IsoTimes;
    char        **H5_FieldLineType;
    char        **H5_IntModel;
    char        **H5_ExtModel;
    
    long int    *H5_Date;
    int         *H5_Doy;
    double      *H5_UTC;
    double      *H5_JD;
    double      *H5_GpsTime;
    double      *H5_TiltAngle;
    int         *H5_InOut;
    int         *H5_OrbitNumber;

    double      **H5_Rgeo;
    double      **H5_Rgeod;
    double      **H5_Rgeod_LatLon;
    double      *H5_Rgeod_Height;
    double      **H5_Rgsm;
    double      **H5_Rsm;
    double      **H5_Rgei;
    double      **H5_Rgse;

    double      *H5_CDMAG_MLAT;
    double      *H5_CDMAG_MLON;
    double      *H5_CDMAG_MLT;
    double      *H5_CDMAG_R;

    double      *H5_EDMAG_MLAT;
    double      *H5_EDMAG_MLON;
    double      *H5_EDMAG_MLT;
    double      *H5_EDMAG_R;

    double      *H5_Kp;
    double      *H5_Dst;

    double      **H5_Bsc_gsm;
    double      *H5_S_sc_to_pfn;
    double      *H5_S_sc_to_pfs;
    double      *H5_S_pfs_to_Bmin;
    double      *H5_S_Bmin_to_sc;
    double      *H5_S_total;
    double      *H5_d2B_ds2;
    double      *H5_Sb0;
    double      *H5_RadiusOfCurv;

    double      **H5_Pfn_geo;
    double      **H5_Pfn_gsm;
    double      **H5_Pfn_geod;
    double      **H5_Pfn_geod_LatLon;
    double      *H5_Pfn_geod_Height;
    double      **H5_Pfn_cdmag;
    double      *H5_Pfn_CD_MLAT;
    double      *H5_Pfn_CD_MLON;
    double      *H5_Pfn_CD_MLT;
    double      **H5_Pfn_edmag;
    double      *H5_Pfn_ED_MLAT;
    double      *H5_Pfn_ED_MLON;
    double      *H5_Pfn_ED_MLT;
    double      **H5_Bfn_geo;
    double      **H5_Bfn_gsm;
    double      *H5_LossConeAngleN;

    double      **H5_Pfs_geo;
    double      **H5_Pfs_gsm;
    double      **H5_Pfs_geod;
    double      **H5_Pfs_geod_LatLon;
    double      *H5_Pfs_geod_Height;
    double      **H5_Pfs_cdmag;
    double      *H5_Pfs_CD_MLAT;
    double      *H5_Pfs_CD_MLON;
    double      *H5_Pfs_CD_MLT;
    double      **H5_Pfs_edmag;
    double      *H5_Pfs_ED_MLAT;
    double      *H5_Pfs_ED_MLON;
    double      *H5_Pfs_ED_MLT;
    double      **H5_Bfs_geo;
    double      **H5_Bfs_gsm;
    double      *H5_LossConeAngleS;

    double      **H5_Pmin_gsm;
    double      **H5_Bmin_gsm;

    double      *H5_Lsimple;
    double      *H5_InvLat;
    double      *H5_Lm_eq;
    double      *H5_InvLat_eq;
    double      *H5_BoverBeq;
    double      *H5_MlatFromBoverBeq;
    double      *H5_M_used;
    double      *H5_M_ref;
    double      *H5_M_igrf;

    double      **H5_Lstar;
    double      **H5_Sb;
    double      **H5_Tb;
    double      **H5_Kappa;
    int         **H5_DriftShellType;
    double      **H5_L;
    double      **H5_Bm;
    double      **H5_I;
    double      **H5_K;

} Lgm_MagEphemData;


/*
 * Function Prototypes
 */
Lgm_MagEphemInfo *Lgm_InitMagEphemInfo( int Verbosity, int MaxPitchAngles );
void Lgm_InitMagEphemInfoDefaults(Lgm_MagEphemInfo *MagEphemInfo, int MaxPitchAngles, int Verbosity);
void Lgm_FreeMagEphemInfo( Lgm_MagEphemInfo  *Info );
void Lgm_FreeMagEphemInfo_Children( Lgm_MagEphemInfo  *MagEphemInfo );

double  j_to_fp_1( double j, double Ek );
double  j_to_fp_2( double j, double Ek0, double Ek1 );
double  k_to_mu_1( double Ek, double alpha, double B );
double  k_to_mu_2( double Ek0, double Ek1, double alpha, double B );
double  Ek_to_v( double Ek, int Species );
void    Lgm_ComputeLstarVersusPA( long int Date, double UTC, Lgm_Vector *u, int nAlpha, double *Alpha, int Quality, int Colorize, Lgm_MagEphemInfo *MagEphemInfo );

void    ReadMagEphemInfoStruct( char *Filename, int *nPitchAngles, Lgm_MagEphemInfo *MagEphemInfo );
void    WriteMagEphemInfoStruct( char *Filename, int nPitchAngles, Lgm_MagEphemInfo *MagEphemInfo );
void    Lgm_WriteMagEphemHeader( FILE *fp, char *CodeVersion, char *ExtModel, int SpiceBody,  char *Spacecraft, int IdNumber, char *IntDesig, char *CmdLine, int nAscend, Lgm_DateTime *Ascend_UTC, Lgm_Vector *Ascend_U, int nPerigee, Lgm_DateTime *Perigee_UTC, Lgm_Vector *Perigee_U, int nApogee, Lgm_DateTime *Apogee_UTC, Lgm_Vector *Apogee_U, Lgm_MagEphemInfo *m );
void    Lgm_WriteMagEphemHeaderHdf( hid_t file, char *argp_program_version, char *ExtModel, int SpiceBody,  char *Spacecraft, int IdNumber, char *IntDesig, char *CmdLine, int nAscend, Lgm_DateTime *Ascend_UTC, Lgm_Vector *Ascend_U, int nPerigee, Lgm_DateTime *Perigee_UTC, Lgm_Vector *Perigee_U, int nApogee, Lgm_DateTime *Apogee_UTC, Lgm_Vector *Apogee_U, Lgm_MagEphemInfo *m, Lgm_MagEphemData *med  );
void    Lgm_WriteMagEphemData( FILE *fp, char *IntModel, char *ExtModel, double Kp, double Dst, Lgm_MagEphemInfo *m );
void    Lgm_WriteMagEphemDataHdf( hid_t file, int iii, Lgm_MagEphemData *m );


Lgm_MagEphemData *Lgm_InitMagEphemData( int nRows, int nPA );                                                                                                                                                                              
void              Lgm_FreeMagEphemData( Lgm_MagEphemData *MagEphemData );



#endif
