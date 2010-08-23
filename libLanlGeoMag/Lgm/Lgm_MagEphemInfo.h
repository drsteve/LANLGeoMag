#ifndef MAGEPHEMINFO_H
#define MAGEPHEMINFO_H

#include <Lgm/Lgm_MagModelInfo.h>
#include <Lgm/Lgm_LstarInfo.h>

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

    Lgm_Vector      P_gsm;          //< S/C position in GSM
    double          B;              //< Local (model) B-field magnitude (i.e. at S/C position)

    Lgm_Vector      Pmin_gsm;       //< position of minimum |B| in GSM
    double          Bmin;           //< Value of |Bmin|
    double          smin;           //< Distance from southern footpoint to Pmin along FL.

    Lgm_Vector      Footprint_Pn;   //< position of northern footpoint (at 120km)
    double          Footprint_Bn;   //< Value of |B| at Footprint_Pn

    Lgm_Vector      Footprint_Ps;   //< position of southern footpoint (at 120km)
    double          Footprint_Bs;   //< Value of |B| at Footprint_Ps

    int             FieldLineType;  //< Field line type. (I.e., LGM_OPEN_IMF, LGM_CLOSED, LGM_OPEN_N_LOBE, LGM_OPEN_S_LOBE, LGM_INSIDE_EARTH, LGM_TARGET_HEIGHT_UNREACHABLE)






    int         UseInterpRoutines;


    /*
     * Array of PAs we need to compute things for...
     */
    int         nAlpha;         // Number of PAs
    double      Alpha[MAX_PITCH_ANGLES];      // Pitch Angles (Degrees)
    Lgm_Vector  Pmn_gsm[MAX_PITCH_ANGLES];    // position of northern |Bmirror|
    Lgm_Vector  Pms_gsm[MAX_PITCH_ANGLES];    // position of southern |Bmirror|
    double      Bm[MAX_PITCH_ANGLES];         // Value of |Bm| (i.e. Bmirror)
    double      I[MAX_PITCH_ANGLES];
    double      Sb[MAX_PITCH_ANGLES];
    double      Tb[MAX_PITCH_ANGLES];         // we calculate this for a 1Mev electron but could do more...
    double      K[MAX_PITCH_ANGLES];

    int         nShellPoints[MAX_PITCH_ANGLES];               // # of point (i.e. FLs) in a shell calculation)
    Lgm_Vector  ShellFootprint_Pn[MAX_PITCH_ANGLES][100];    // north footprints of shell lines
    Lgm_Vector  ShellFootprint_Ps[MAX_PITCH_ANGLES][100];    // south footprints of shell lines
    Lgm_Vector  ShellMirror_Pn[MAX_PITCH_ANGLES][100];       // north mirror locations
    Lgm_Vector  ShellMirror_Ps[MAX_PITCH_ANGLES][100];       // south mirror locations
    double      ShellMirror_Sn[MAX_PITCH_ANGLES][100];       // north mirror locations (dist along FL)
    double      ShellMirror_Ss[MAX_PITCH_ANGLES][100];       // south mirror locations (dist along FL)
    double      ShellI[MAX_PITCH_ANGLES][100];               // Individual I values computed for each FL

    /*
     * these are like the variables in the LstarInfo structure. Except they have an extra
     * dimension to hold pitch angle as well.
     */
    int         nFieldPnts[MAX_PITCH_ANGLES][48];     
    double      s_gsm[MAX_PITCH_ANGLES][48][1000];
    double      Bmag[MAX_PITCH_ANGLES][48][1000];
    double      x_gsm[MAX_PITCH_ANGLES][48][1000];
    double      y_gsm[MAX_PITCH_ANGLES][48][1000];
    double      z_gsm[MAX_PITCH_ANGLES][48][1000];

    double      Mcurr;
    double      Mref;
    double      Mused;
    double      LHilton[MAX_PITCH_ANGLES];
    double      LMcIlwain[MAX_PITCH_ANGLES];
    double      Lstar[MAX_PITCH_ANGLES];


} Lgm_MagEphemInfo;


/*
 * Protos
 */
double  j_to_fp_1( double j, double Ek );
double  j_to_fp_2( double j, double Ek0, double Ek1 );
double Ek_to_mu_1( double Ek, double alpha, double B );
double Ek_to_mu_2( double Ek0, double Ek1, double alpha, double B );
double  Ek_to_v( double Ek, int Species );
void ComputeFieldLineQuantities( long int Date, double UTC, Lgm_Vector *u, int nAlpha, double *Alpha, int Quality, Lgm_MagEphemInfo *MagEphemInfo );

Lgm_MagEphemInfo *Lgm_InitMagEphemInfo( int Verbosity );
void Lgm_FreeMagEphemInfo( Lgm_MagEphemInfo  *Info );


#endif
