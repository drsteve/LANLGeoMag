#ifndef LGM_SGP_H
#define LGM_SGP_H
#include <math.h>

/*
 * This structure contains raw lines for Two Line Elements (TLEs)
 * and also the decoded values. Here is how the TLEs are formated;
 *
 *  Line 0:
 *      Columns     Example                         Description
 *      -------     -----------         ----------------------------------------------------------
 *      1-24        ISS (ZARYA)         The common name for the object based on information from 
                                        the SatCat. 
 *
 *  Line 1:
 *      Columns     Example                         Description
 *      -------     -----------         ----------------------------------------------------------
 *      1           1                   Line Number
 *      3-7         25544               Object Identification Number Search by object ID
 *      8           U                   Elset Classification
 *      10-17       98067A              International Designator Search by International Designator
 *      19-32       04236.56031392      Element Set Epoch (UTC)
 *      34-43       .00020137           1st Derivative of the Mean Motion with respect to Time
 *      45-52(46-53??)       12345-7             2nd Derivative of the Mean Motion with respect to Time 
 *                                      NOTE: decimal point assumed.  This type of format is to be 
 *                                      read as: 0.12345e-7
 *      54-61       16538-3             B* Drag Term
 *      63  0       Element             Set Type
 *      65-68       513                 Element Number
 *      69          5                   Checksum
 *
 *  Line 2:
 *      Columns     Example                         Description
 *      -------     -----------         ----------------------------------------------------------
 *      1           2                   Line Number
 *      3-7         25544               Object Identification Number
 *      9-16        51.6335             Orbit Inclination (degrees)
 *      18-25       344.7760            Right Ascension of Ascending Node (degrees)
 *      27-33       0007976             Eccentricity (decimal point assumed)
 *      35-42       126.2523            Argument of Perigee (degrees)
 *      44-51       325.9359            Mean Anomaly (degrees)
 *      53-63       15.70406856         Mean Motion (revolutions/day)
 *      64-68       32890               Revolution Number at Epoch
 *      69          3                   Checksum
 *
 */
typedef struct _SgpTLE {

    // Raw lines
    char    Line0[80];  // Contains the common name for the object
    char    Line1[80];  // Contains identifiers derivatives of mean motion, drag term etc..
    char    Line2[80];  // Contains orbital elements (RA of Asc Node, mean anomaly, etc..)


    // Decoded Line0
    char    Name[80];           // Common object name


    // Decoded Line1
    int     IdNumber;           // Object identification number
    char    ElsetClass;         // Elset Classification
    char    IntDesig[20];       // International designator
    double  ElementSetEpoch;    // Element Set Epoch Time
    double  dMMdT1;             // 1st Derivative of the Mean Motion with respect to Time
    double  dMMdT2;             // 2nd Derivative of the Mean Motion with respect to Time
    double  BstarDrag;          // B* Drag Term
    int     ElementSetType;     // Element Set Type
    int     ElementSetNum;      // Element Set Number
    int     Line1CheckSum;      // Line1 Checksum


    // Decoded Line2
    double    Inclination;        // Orbit Inclination (degrees)
    double    RAofAscNode;        // Right Ascension of Ascending Node (degrees)
    double    Eccentricity;       // Eccentricity
    double    ArgOfPerigee;       // Argument of Perigee (degrees)
    double    MeanAnomaly;        // Mean Anomaly (degrees)
    double    MeanMotion;         // Mean Motion (revolutions/day)
    int       RevNumAtEpoch;      // Revolution Number at Epoch
    int       Line2CheckSum;      // Line2 Checksum

    // Extra (derived) stuff
    long int  Date;
    double    UT;
    int       Year;
    int       Month;
    int       Day;
    int       Doy;
    char      Dow[5];
    double    JD;             // Julian Date of Epoch
    double    Period;         // in minutes
    char      IntDesig2[20];  // e.g. 1989-046A
    char      ObjectType[20]; // SAT, R/B or DEB
    char      EpochStr[20];   // A more readabl version of the epoch time...
    double    YYYYDDDdFRAC;   // YYYYDDD.FRAC
    


} _SgpTLE;





typedef struct {

    /*
     * THESE WERE FOR THGE STR3 codes
     */
    int     IFLAG;
    double  XMO;
    double  XNODEO;
    double  OMEGAO;
    double  EO;
    double  XINCL;
    double  XNO;
    double  XNDT2O;
    double  XNDD6O;
    double  BSTAR;
//    double  X;
//    double  Y;
//    double  Z;
    double  XDOT;
    double  YDOT;
    double  ZDOT;
    double  EPOCH;
    double  DS50;

    /*
     * New vars
     */
    double  argpdot, argpo, atime, aycof, bstar, cc1, cc4, cc5, con41, d2, d2201, d2211;
    double  d3, d3210, d3222, d4, d4410, d4422, d5220, d5232, d5421, d5433, dedt, del1;
    double  del2, del3, delmo, didt, dmdt, dnodt, domdt, e3, ecco, ee2, error, eta, gsto;
    double  inclo, mdot, mo, no, nodecf, nodedot, nodeo, omgcof, peo;
    double  pgho, pho, pinco, plo, se2, se3, sgh2, sgh3, sgh4, sh2, sh3, si2, si3, sinmao;
    double  sl2, sl3, sl4, t, t2cof, t3cof, t4cof, t5cof, x1mth2, x7thm1, xfact, xgh2;
    double  xgh3, xgh4, xh2, xh3, xi2, xi3, xl2, xl3, xl4, xlamo, xlcof, xli, xmcof, xni;
    double  zmol, zmos;

//    double  a, altp, alta, epochdays, jdsatepoch, nddot, ndot, bstar, rcse, inclo, nodeo; 
//    double  ecco, argpo, mo, no;
    
    int     GravConst, irez;
    char    init, method, isimp;

    double  X;
    double  Y;
    double  Z;

    double  VX;
    double  VY;
    double  VZ;




} _SgpInfo;

#ifndef TRUE
#define TRUE  1
#endif
#ifndef FALSE
#define FALSE 0
#endif

/*
 * Constants used by SGP. These are certainly not as acurate as they could be. But DONT
 * change them. The Two Line Elements were created with these same values -- so we have 
 * to use them here too.
 */
#define SGP_CK2         5.413080e-4 
#define SGP_CK4         0.62098875e-6 
#define SGP_E6A         1.0e-6 
#define SGP_QOMS2T      1.88027916e-9 
#define SGP_S           1.01222928 
#define SGP_TOTHRD      0.66666667          // 2/3
#define SGP_XJ3         (-0.253881e-5)
#define SGP_XKE         0.743669161e-1 
#define SGP_XKMPER      6378.135            // kilometers/earth radii
#define SGP_XMNPDA      1440.0              // time units/day
#define SGP_AE          1.0                 // distance units/earth radii
#define SGP_DE2RA       0.174532925e-1      // radians/degree
#define SGP_PI          3.14159265          // pi
#define SGP_PIO2        1.57079633          // pi/2
#define SGP_TWOPI       6.2831853           // 2*pi
#define SGP_X3PIO2      4.71238898          // 3*pi/2



// Grav constants
#define SGP_wgs72old    0
#define SGP_wgs72       1
#define SGP_wgs84       2


#ifndef M_PI
#define M_PI            3.141592653589793238462643      // M_PI
#endif

#ifndef M_2PI
#define M_2PI           6.283185307179586476925287      // 2*M_PI
#endif





int     LgmSgp_TleChecksum( char *Line );
void    Lgm_SgpDecodeTle( char *Line0, char *Line1, char *Line2, _SgpTLE *TLE, int Verbosity );
int     LgmSgp_ReadTlesFromFile( char *Filename, int *nTLEs, _SgpTLE *TLEs, int Verbosity );
int     LgmSgp_ReadTlesFromStrings( char *Line0, char *Line1, char *Line2, int *nTLEs, _SgpTLE *TLEs, int Verbosity );
int     LgmSgp_FindTLEforGivenTime( int nTLEs, _SgpTLE *TLEs, int SortOrder, double JD, int Verbosity );



/*
 * OLD SGP stuff - from "Spacetrack Report #3"
 */
void LgmSgp_InitElements( _SgpInfo *s, _SgpTLE *t );
int LgmSgp_SGP_STR3( double TSINCE, _SgpInfo *s );
int LgmSgp_SGP4_STR3( double TSINCE, _SgpInfo *s );
int LgmSgp_SDP4_STR3( double TSINCE, _SgpInfo *s );
int LgmSgp_SGP8_STR3( double TSINCE, _SgpInfo *s );
int LgmSgp_SDP8_STR3( double TSINCE, _SgpInfo *s );

/*
 * NEW SGP stuff - from "Revisiting Spacetrack Report #3"
 */
int LgmSgp_SGP4( double TSINCE, _SgpInfo *s );

int LgmSgp_SGP4_Init( _SgpInfo *s, _SgpTLE *t );

void LgmSgp_GetGravConst( int whichconst, double *tumin, double *radiusearthkm,
                    double *xke, double *j2, double *j3, double *j4, double *j3oj2 );

void LgmSgp_dpper( double inclo, char init, double *ep, double *inclp, 
                    double *nodep, double *argpp, double *mp, _SgpInfo *s );

void LgmSgp_dspace( double tc, double *atime, double *em, double *argpm, double *inclm, double *xli, 
                    double *mm, double *xni, double *nodem, double *dndt, double *nm, _SgpInfo *s );

void LgmSgp_initl( int satn, int whichconst, double ecco, double epoch, double inclo,
                    double *no, char *method, double *ainv, double *ao, double *con41, double *con42,
                    double *cosio, double *cosio2, double *eccsq, double *omeosq, double *posq,
                    double *rp, double *rteosq, double *sinio, double *gsto);

void LgmSgp_dscom( double epoch, double ep, double argpp, double tc, double inclp, double nodep, double np,
                    double *snodm, double *cnodm, double *sinim, double *cosim, double *sinomm,
                    double *cosomm,double *day, double *e3, double *ee2, double *em,
                    double *emsq, double *gam, double *peo, double *pgho, double *pho,
                    double *pinco, double *plo, double *rtemsq, double *se2, double *se3,
                    double *sgh2, double *sgh3, double *sgh4, double *sh2, double *sh3,
                    double *si2, double *si3, double *sl2, double *sl3, double *sl4,
                    double *s1, double *s2, double *s3, double *s4, double *s5,
                    double *s6, double *s7, double *ss1, double *ss2, double *ss3,
                    double *ss4, double *ss5, double *ss6, double *ss7, double *sz1,
                    double *sz2, double *sz3, double *sz11, double *sz12, double *sz13,
                    double *sz21, double *sz22, double *sz23, double *sz31, double *sz32,
                    double *sz33, double *xgh2, double *xgh3, double *xgh4, double *xh2,
                    double *xh3, double *xi2, double *xi3, double *xl2, double *xl3,
                    double *xl4, double *nm, double *z1, double *z2, double *z3,
                    double *z11, double *z12, double *z13, double *z21, double *z22,
                    double *z23, double *z31, double *z32, double *z33, double *zmol, double *zmos);

void LgmSgp_dsinit( int whichconst, double cosim, double emsq, double argpo,
                    double s1, double s2, double s3, double s4, double s5, double sinim, double ss1,
                    double ss2, double ss3, double ss4, double ss5, double sz1, double sz3, double sz11,
                    double sz13, double sz21, double sz23, double sz31, double sz33, double t, double tc,
                    double gsto, double mo, double mdot, double no, double nodeo, double nodedot,
                    double xpidot, double z1, double z3, double z11, double z13, double z21, double z23,
                    double z31, double z33, double ecco, double eccsq,
                    double *em, double *argpm, double *inclm, double *mm, double *nm, double *nodem,
                    int *irez, double *atime, double *d2201, double *d2211, double *d3210, double *d3222,
                    double *d4410, double *d4422, double *d5220, double *d5232, double *d5421, double *d5433,
                    double *dedt, double *didt, double *dmdt, double *dndt, double *dnodt, double *domdt,
                    double *del1, double *del2, double *del3, double *xfact, double *xlamo, double *xli, double *xni);

double LgmSgp_gstime( double jdut1 );

#endif
