#ifndef MAGEPHEMINFO_H
#define MAGEPHEMINFO_H

#include <Lgm/Lgm_MagModelInfo.h>
#include <Lgm/Lgm_LstarInfo.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>


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


} Lgm_MagEphemInfo;


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
void    Lgm_WriteMagEphemHeader( FILE *fp, char *Spacecraft, int IdNumber, char *IntDesig, char *CmdLine, int nPerigee, Lgm_DateTime *Perigee_UTC, Lgm_Vector *Perigee_U, int nApogee, Lgm_DateTime *Apogee_UTC, Lgm_Vector *Apogee_U, Lgm_MagEphemInfo *m );
void    Lgm_WriteMagEphemData( FILE *fp, char *IntModel, char *ExtModel, double Kp, double Dst, Lgm_MagEphemInfo *m );




#endif
