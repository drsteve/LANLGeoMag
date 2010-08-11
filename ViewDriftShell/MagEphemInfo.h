#ifndef MAGEPHEMINFO_H
#define MAGEPHEMINFO_H

#include <Lgm/Lgm_MagModelInfo.h>
#include <Lgm/Lgm_LstarInfo.h>



/*
 * Structure to hold all relevant magnetic ephemeris quantities.  I.e.,
 * everything from physical location of S/C to magnetic coordinates adiabatic
 * invariants etc...
 */
typedef struct _MagEphemInfo {
    
    long int    Date;
    double      UT;
    double      Lat;
    double      Lon;
    double      Rad;

    Lgm_Vector  P_gsm;      // S/C position in GSM
    double      B;          // Local B-field magnitude (i.e. at S/C position)

    Lgm_Vector  Pmin_gsm;   // position of minimum |B| in GSM
    double      Bmin;       // Value of |Bmin|

    int         UseInterpRoutines;

    /*
     * Array of PAs we need to compute things for...
     */
    int         nAlpha;         // Number of PAs
    double      Alpha[90];      // Pitch Angles (Degrees)
    Lgm_Vector  Pmn_gsm[90];    // position of northern |Bmirror|
    Lgm_Vector  Pms_gsm[90];    // position of southern |Bmirror|
    double      Bm[90];         // Value of |Bm| (i.e. Bmirror)
    double      I[90];
    double      Sb[90];
    double      Tb[90];         // we calculate this for a 1Mev electron but could do more...
    double      K[90];

    int         nShellPoints[90];               // # of point (i.e. FLs) in a shell calculation)
    Lgm_Vector  ShellFootprint_Pn[90][100];    // north footprints of shell lines
    Lgm_Vector  ShellFootprint_Ps[90][100];    // south footprints of shell lines
    Lgm_Vector  ShellMirror_Pn[90][100];       // north mirror locations
    Lgm_Vector  ShellMirror_Ps[90][100];       // south mirror locations
    double      ShellMirror_Sn[90][100];       // north mirror locations (dist along FL)
    double      ShellMirror_Ss[90][100];       // south mirror locations (dist along FL)

    /*
     * these are like the variables in the LstarInfo structure. Except they have an extra
     * dimension to hold pitch angle as well.
     */
    int         nFieldPnts[90][24];     
    double      s_gsm[90][24][1000];
    double      Bmag[90][24][1000];
    double      x_gsm[90][24][1000];
    double      y_gsm[90][24][1000];
    double      z_gsm[90][24][1000];

    double      Mcurr;
    double      Mref;
    double      Mused;
    double      LHilton[90];
    double      LMcIlwain[90];
    double      Lstar[90];


} _MagEphemInfo;


/*
 * Protos
 */
double  j_to_fp_1( double j, double Ek );
double  j_to_fp_2( double j, double Ek0, double Ek1 );
double Ek_to_mu_1( double Ek, double alpha, double B );
double Ek_to_mu_2( double Ek0, double Ek1, double alpha, double B );
double  Ek_to_v( double Ek, int Species );

/*
 * Global defs
 */
#ifdef MAIN
/*
 *  Make MagEphemInfo Global. We could get around this,
 *  but it adds complexity.
 */
_MagEphemInfo   *MagEphemInfo;
#else
extern _MagEphemInfo   *MagEphemInfo;
#endif



#endif
