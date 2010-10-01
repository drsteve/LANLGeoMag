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


/**
 * Structure to hold all relevant magnetic ephemeris quantities.  I.e.,
 * everything from physical location of S/C to magnetic coordinates adiabatic
 * invariants etc...
 */
typedef struct Lgm_MagEphemInfo {

    Lgm_LstarInfo   *LstarInfo;
    int             LstarQuality;
    int             SaveShellLines;
    
    long int        Date;           //< Date in YYYYMMDD format
    double          UTC;            //< UTC in decimal hours

    double          Lat;            //< Geographic Latitude
    double          Lon;            //< Geographic Longitude
    double          Rad;            //< Geographic Radius

    Lgm_Vector      P;              //< S/C position in GSM
    double          S;              //< Distance along FL from southern footpoint to S/C location in Re.
    double          B;              //< Local (model) B-field magnitude (i.e. at S/C position)

    Lgm_Vector      Pmin;           //< position of minimum |B| in GSM
    double          Bmin;           //< Value of |Bmin|
    double          Smin;           //< Distance from southern footpoint to Pmin along FL.

    Lgm_Vector      Spherical_Footprint_Pn;   //< position of northern footpoint (at 120km above spherical Earth)
    double          Spherical_Footprint_Sn;   //< Distance along FL from southern foorpoint in Re
    double          Spherical_Footprint_Bn;   //< Value of |B| at Spherical_Footprint_Pn

    Lgm_Vector      Spherical_Footprint_Ps;   //< position of southern footpoint (at 120km above spherical Earth)
    double          Spherical_Footprint_Ss;   //< Distance along FL from southern foorpoint in Re (i.e. this one is zero by definition)
    double          Spherical_Footprint_Bs;   //< Value of |B| at Spherical_Footprint_Ps

    Lgm_Vector      Ellipsoid_Footprint_Pn;   //< position of northern footpoint (at 120km above WGS84 ellipsoid)
    double          Ellipsoid_Footprint_Sn;   //< Distance along FL from southern foorpoint in Re
    double          Ellipsoid_Footprint_Bn;   //< Value of |B| at Ellipsoid_Footprint_Pn

    Lgm_Vector      Ellipsoid_Footprint_Ps;   //< position of southern footpoint (at 120km above WGS84 ellipsoid)
    double          Ellipsoid_Footprint_Ss;   //< Distance along FL from southern foorpoint in Re (i.e. this one is zero by definition)
    double          Ellipsoid_Footprint_Bs;   //< Value of |B| at Ellipsoid_Footprint_Ps

    int             FieldLineType;  //< Field line type. (I.e., LGM_OPEN_IMF, LGM_CLOSED, LGM_OPEN_N_LOBE, LGM_OPEN_S_LOBE, LGM_INSIDE_EARTH, LGM_TARGET_HEIGHT_UNREACHABLE)






    int         UseInterpRoutines;


    /*
     * Array of PAs we need to compute things for...
     */
    int         nAlpha;     // Number of PAs
    double      *Alpha;     // Alpha[MAX_PITCH_ANGLES]    Pitch Angles (Degrees)
    Lgm_Vector  *Pmn_gsm;   // Pmn_gsm[MAX_PITCH_ANGLES]  position of northern |Bmirror|
    Lgm_Vector  *Pms_gsm;   // Pms_gsm[MAX_PITCH_ANGLES]  position of southern |Bmirror|
    double      *Bm;        // Bm[MAX_PITCH_ANGLES]       Value of |Bm| (i.e. Bmirror)
    double      *I;         // I[MAX_PITCH_ANGLES]
    double      *Sb;        // Sb[MAX_PITCH_ANGLES]
    double      *Tb;        // Tb[MAX_PITCH_ANGLES]       we calculate this for a 1Mev electron but could do more...
    double      *K;         // K[MAX_PITCH_ANGLES]

    int         *nShellPoints; //  nShellPoints[MAX_PITCH_ANGLES] # of point (i.e. FLs) in a shell calculation)


    /*
     * Footpoints XXXkm above spherical Earth defined as sphere with
     * radius==WGS84_A (i.e. equatorial radius of WGS84 model)
     */
    Lgm_Vector  **ShellSphericalFootprint_Pn; // [MAX_PITCH_ANGLES][100]  north footprints of shell lines
    double      **ShellSphericalFootprint_Sn; // [MAX_PITCH_ANGLES][100]  north mirror locations (dist along FL)
    double      **ShellSphericalFootprint_Bn; // [MAX_PITCH_ANGLES][100]  north mirror locations (dist along FL)

    Lgm_Vector  **ShellSphericalFootprint_Ps; // [MAX_PITCH_ANGLES][100]  north footprints of shell lines
    double      **ShellSphericalFootprint_Ss; // [MAX_PITCH_ANGLES][100]  north mirror locations (dist along FL)
    double      **ShellSphericalFootprint_Bs; // [MAX_PITCH_ANGLES][100]  north mirror locations (dist along FL)

    /*
     * footpoints at XXXkm above surface of ellipsoid.
     */
    Lgm_Vector  **ShellEllipsoidFootprint_Ps; // [MAX_PITCH_ANGLES][100]  north footprints of shell lines
    double      **ShellEllipsoidFootprint_Ss; // [MAX_PITCH_ANGLES][100]  north mirror locations (dist along FL)
    double      **ShellEllipsoidFootprint_Bs; // [MAX_PITCH_ANGLES][100]  north mirror locations (dist along FL)

    Lgm_Vector  **ShellEllipsoidFootprint_Pn; // [MAX_PITCH_ANGLES][100]  north footprints of shell lines
    double      **ShellEllipsoidFootprint_Sn; // [MAX_PITCH_ANGLES][100]  north mirror locations (dist along FL)
    double      **ShellEllipsoidFootprint_Bn; // [MAX_PITCH_ANGLES][100]  north mirror locations (dist along FL)


    Lgm_Vector  **ShellMirror_Pn; // [MAX_PITCH_ANGLES][100]   north mirror locations
    double      **ShellMirror_Sn; // [MAX_PITCH_ANGLES][100]   north mirror locations (dist along FL)

    Lgm_Vector  **ShellMirror_Ps; // [MAX_PITCH_ANGLES][100]   south mirror locations
    double      **ShellMirror_Ss; // [MAX_PITCH_ANGLES][100]   south mirror locations (dist along FL)


    double      **ShellI;         // [MAX_PITCH_ANGLES][100]   Individual I values computed for each FL

    /*
     * these are like the variables in the LstarInfo structure. Except they have an extra
     * dimension to hold pitch angle as well.
     */
    int         **nFieldPnts;   // [MAX_PITCH_ANGLES][48]     
    double      ***s_gsm;       // [MAX_PITCH_ANGLES][48][1000]
    double      ***Bmag;        // [MAX_PITCH_ANGLES][48][1000]
    double      ***x_gsm;       // [MAX_PITCH_ANGLES][48][1000]
    double      ***y_gsm;       // [MAX_PITCH_ANGLES][48][1000]
    double      ***z_gsm;       // [MAX_PITCH_ANGLES][48][1000]

    double      Mcurr;
    double      Mref;
    double      Mused;
    double      *LHilton;    // LHilton[MAX_PITCH_ANGLES]
    double      *LMcIlwain;  // LMcIlwain[MAX_PITCH_ANGLES]
    double      *Lstar;      // Lstar[MAX_PITCH_ANGLES]


} Lgm_MagEphemInfo;


/*
 * Protos
 */
double  j_to_fp_1( double j, double Ek );
double  j_to_fp_2( double j, double Ek0, double Ek1 );
double Ek_to_mu_1( double Ek, double alpha, double B );
double Ek_to_mu_2( double Ek0, double Ek1, double alpha, double B );
double  Ek_to_v( double Ek, int Species );
//void ComputeFieldLineQuantities( long int Date, double UTC, Lgm_Vector *u, int nAlpha, double *Alpha, int Quality, Lgm_MagEphemInfo *MagEphemInfo );

Lgm_MagEphemInfo *Lgm_InitMagEphemInfo( int Verbosity, int MaxPitchAngles );
void Lgm_FreeMagEphemInfo( Lgm_MagEphemInfo  *Info );

void ReadMagEphemInfoStruct( char *Filename, int *nPitchAngles, Lgm_MagEphemInfo *MagEphemInfo );
void WriteMagEphemInfoStruct( char *Filename, int nPitchAngles, Lgm_MagEphemInfo *MagEphemInfo );


#endif

/*
 *    $Id$
 */

