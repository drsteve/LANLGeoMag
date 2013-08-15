#ifndef LGM_NRLMSISE00_H
#define LGM_NRLMSISE00_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define AMOD(A,P) ( (A-((int)((double)A/(double)P)*P)) )
#define AMIN1(A,B) ( ((A<B)?(A):(B)) )
#define AMAX1(A,B) ( ((A>B)?(A):(B)) )

typedef struct Lgm_Msis00Info {

    char    NAME[80];
    char    ISDATE[80];
    char    ISTIME[80];

    // PARM7 COMMON BLOCK VARS
    double   *PT, **PD, *PS, **PDL, **PTL, **PMA, *SAM;
    double   **PD_rc, **PTL_rc, **PMA_rc;

    // MAVG7 COMMON BLOCK VAR
    double   *PAVGM;

    // LOWER7 COMMON BLOCK VARS
    double   *PTM, **PDM;

    // GTS3C COMMON BLOCK VARS
    double  TLB, S, DB04, DB16, DB28, DB32, DB40, DB48, DB01, ZA, T0, Z0, G0, RL, DD, DB14, TR12;

    // MESO7 COMMON BLOCK VARS
    double  TN1[6], TN2[5], TN3[6], TGN1[3], TGN2[3], TGN3[3];

    // PARMB COMMON BLOCK VARS
    double  GSURF, RE;

    // CSW COMMON BLOCK VARS
    double *SW, ISW, *SWC;

    // LPOLY COMMON BLOCK VARS
    int     IYR;
    double  PLG[10][5], CTLOC, STLOC, C2TLOC, S2TLOC, C3TLOC, S3TLOC, DAY, DF, DFA, APD, APDF, APT[5], LONG;

    // METSEL COMMON BLOCK VAR
    int IMR;

    // DMIX COMMON BLOCK VARS 
    double  DM04, DM16, DM28, DM32, DM40, DM01, DM14;

    // TTEST COMMON BLOCK VARS (only T[] used ? )
    double TINF, GB, ROUT, T[16]; 

    // Cached Vars
    double  DAYL, XL, TLL, SV[26], SAV[26];

    double  IYDL[3], SECL[3], GLATL[3], GLL[3], STLL[3], FAL[3], FL[3], APL[8][3], SWL[26][3], SWCL[26][3]; 
    
    // Some variables in GTD7() that need to be saved.
    double  GTD7_ALAST, GTD7_MSSL, GTS7_ALAST;

    // Some variables in GLOBE7() that need to be saved.
    double  GLOBE7_P14, GLOBE7_P18, GLOBE7_P32, GLOBE7_P39;

    // Some variables in GLOB7S() that need to be saved.
    double  GLOB7S_P14, GLOB7S_P18, GLOB7S_P32, GLOB7S_P39;


} Lgm_Msis00Info;


// Function prototypes
void   GTD7( int IYD, double SEC, double ALT, double GLAT, double GLONG, double STL, double F107A, double F107, double *AP, double MASS, double *D, double *T, Lgm_Msis00Info *p );
void   GTD7D( int IYD, double SEC, double ALT, double GLAT, double GLONG, double STL, double F107A, double F107, double *AP, double MASS, double *D, double *T, Lgm_Msis00Info *p );
void   GHP7( int IYD, double SEC, double ALT, double GLAT, double GLONG, double STL, double F107A, double F107, double *AP, double *D, double *T, double PRESS, Lgm_Msis00Info *p );
void   GLATF( double LAT, double *GV, double *REFF );
int    VTST7( double IYD, double SEC, double GLAT, double GLONG, double STL, double F107A, double F107, double *AP, int IC, Lgm_Msis00Info *p );
void   GTS7( double IYD, double SEC, double ALT, double GLAT, double GLONG, double STL, double F107A, double F107, double *AP, double MASS, double *D, double *T, Lgm_Msis00Info *p );
void   METERS( int METER, Lgm_Msis00Info *p );
double SCALH( double ALT, double XM, double TEMP, Lgm_Msis00Info *p );
double SG0( double EX, double *AP, double *P );
double GLOBE7( double YRD, double SEC, double LAT, double LONG, double TLOC, double F107A, double F107, double *AP, double *P, Lgm_Msis00Info *p );
void   TSELEC( double *SV, Lgm_Msis00Info *p );
void   TRETRV( double *SVV, Lgm_Msis00Info *p );
double GLOB7S( double *P, Lgm_Msis00Info *p );
double DENSU( double ALT, double DLB, double TINF, double TLB, double XM, double ALPHA, double *TZ, double ZLB, double S2, int MN1, double *ZN1, double *TN1, double *TGN1, Lgm_Msis00Info *p );
double DENSM( double ALT, double D0, double XM, double *TZ, int MN3, double *ZN3, double *TN3, double *TGN3, int MN2, double *ZN2, double *TN2, double *TGN2, Lgm_Msis00Info *p );
void   SPLINE( double *X, double *Y, int N, double YP1, double YPN, double *Y2 );
void   SPLINT( double *XA, double *YA, double *Y2A, int N, double X, double *Y );
void   SPLINI( double *XA, double *YA, double *Y2A, int N, double X, double *YI );
double DNET( double DD, double DM, double ZHM, double XMM, double XM );
double CCOR( double ALT, double R, double H1, double ZH );
double CCOR2( double ALT, double R, double H1, double ZH, double H2 );
Lgm_Msis00Info *InitMsis00( );
void   Lgm_FreeMsis00( Lgm_Msis00Info *p );


#endif
