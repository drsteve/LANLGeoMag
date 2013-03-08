#ifndef LGM_TSYG1996_H
#define LGM_TSYG1996_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/*
 *    Lgm_Tsyg1996.h
 */


#ifndef TRUE
#define TRUE    1
#endif
#ifndef FALSE
#define FALSE   0
#endif



/*
 * Define some structures to accomodate the FORTRAN common blocks
 */
typedef struct _CB_T96_DX1        { double  DX, SCALEIN, SCALEOUT;                      } _CB_T96_DX1;
typedef struct _CB_T96_LOOPDIP1   { double  TILT, XCENTRE[3], RADIUS[3], DIPX, DIPY;    } _CB_T96_LOOPDIP1;
typedef struct _CB_T96_RHDR       { double  RH, DR;                                     } _CB_T96_RHDR;
typedef struct _CB_T96_WARP       { double  CPSS, SPSS, DPSRR, RPS, WARP, D, XS, ZS, 
                                            DXSX, DXSY, DXSZ, DZSX, DZSY, DZSZ, DZETAS, 
                                            DDZETADX, DDZETADY, DDZETADZ, ZSWW;         } _CB_T96_WARP;



/*
 * Define a structure to hold all of the info needed in TS04
 */
typedef struct LgmTsyg1996_Info {

    /*
     * These are versions of the original fortran common blocks
     */
    _CB_T96_DX1         CB_T96_DX1;
    _CB_T96_LOOPDIP1    CB_T96_LOOPDIP1;
    _CB_T96_RHDR        CB_T96_RHDR;
    _CB_T96_WARP        CB_T96_WARP;


    /*
     *  SIN(PS) and COS(PS) are calculated 6 times each. This is very wasteful
     *  especially since we have it already. Lets pass sin_psi_op and
     *  cos_psi_op as globals....
     */
    double    sin_psi, cos_psi;


    double      OLD_PS;
    double      OLD_X;
    double      OLD_Y;
    double      OLD_Z;
    double      OLD_PDYN;


    int         INTERCON_M_FLAG;
    double      P[4], R[4], RP[4], RR[4], SQPR[4][4];




} LgmTsyg1996_Info;







/*
 *  Function declarations
 */
void    Lgm_Init_T96( LgmTsyg1996_Info *t );
void    Tsyg_T96( int IOPT, double *PARMOD, double PS, double SINPS, double COSPS, double X, double Y, double Z, double *BX, double *BY, double *BZ, LgmTsyg1996_Info *tInfo ) ;
void    DIPSHLD_T96( double PS, double X, double Y, double Z, double *BX, double *BY, double *BZ, LgmTsyg1996_Info *tInfo ) ;
void    CYLHARM_T96( double A[], double X, double Y, double Z, double *BX, double *BY, double *BZ ) ;
void    CYLHAR1_T96( double A[], double X, double Y, double Z, double *BX, double *BY, double *BZ ) ;
void    INTERCON_T96( double X, double Y, double Z, double *BX, double *BY, double *BZ, LgmTsyg1996_Info *tInfo );
void    TAILRC96_T96( double SPS, double X, double Y, double Z, double *BXRC, double *BYRC, double *BZRC, double *BXT2, double *BYT2, double *BZT2, double *BXT3, double *BYT3, double *BZT3, LgmTsyg1996_Info *t ) ;
void    RINGCURR96_T96( double X, double Y, double Z, double *BX, double *BY, double *BZ, LgmTsyg1996_Info *t ) ;
void    TAILDISK_T96( double X, double Y, double Z, double *BX, double *BY, double *BZ, LgmTsyg1996_Info *t ) ;
void    TAIL87_T96( double X, double Z, double *BX, double *BZ, LgmTsyg1996_Info *t ) ;
void    SHLCAR3X3_T96( double A[] , double X, double Y, double Z, double SPS, double *HX, double *HY, double *HZ ) ;
void    BIRK1TOT_02_T96( double PS, double X, double Y, double Z, double *BX, double *BY, double *BZ, LgmTsyg1996_Info *t ) ;
void    BIRK2TOT_02_T96( double PS, double X, double Y, double Z,double *BX, double *BY,double *BZ, LgmTsyg1996_Info *tInfo ) ;
void    DIPLOOP1_T96( double XI[5], double D[4][27], double *XX, double *YY, LgmTsyg1996_Info *t ) ;
void    CIRCLE_T96( double X, double Y, double Z, double RL, double *BX, double *BY, double *BZ ) ;
void    CROSSLP_T96( double X, double Y, double Z, double *BX, double *BY, double *BZ, double XC, double RL, double AL ) ;
void    DIPXYZ_T96( double X, double Y, double Z, double *BXX, double *BYX, double *BZX, double *BXY, double *BYY, double *BZY, double *BXZ, double *BYZ, double *BZZ ) ;
void    BIRK1SHLD_T96( double PS, double X, double Y, double Z, double *BX, double *BY, double *BZ ) ;
void    BIRK2SHL_T96( double X, double Y, double Z, double PS, double *HX, double *HY, double *HZ ) ;
void    R2_BIRK_T96( double X, double Y, double Z, double PS, double *BX, double *BY, double *BZ ) ;
void    R2INNER_T96( double X, double Y, double Z, double *BX, double *BY, double *BZ ) ;
void    BCONIC_T96( double X, double Y, double Z, double CBX[], double CBY[], double CBZ[], int NMAX ) ;
void    DIPDISTR_T96( double X, double Y, double Z, double *BX, double *BY, double *BZ, int MODE ) ;
void    R2OUTER_T96( double X, double Y, double Z, double *BX, double *BY, double *BZ ) ;
void    LOOPS4_T96(double X, double Y, double Z, double *BX, double *BY, double *BZ, double XC, double YC, double ZC, double R, double THETA, double PHI ) ;
void    R2SHEET_T96( double X, double Y, double Z, double *BX, double *BY, double *BZ ) ;
void    DIPOLE_T96( double PS, double X, double Y, double Z, double *BX, double *BY, double *BZ, LgmTsyg1996_Info *tInfo ) ;
void    CONDIP1_T96( double XI[5], double D[4][80], double *XX, double *YY, double *ZZ, LgmTsyg1996_Info *t );
double  BES0_T96( double X );
double  BES1_T96( double X );
double  BES_T96( double X, int K );
double  XKSI_T96( double X, double Y, double Z );
double  FEXP_T96( double S, double A );
double  FEXP1_T96( double S, double A );
double  TKSI_T96( double XKSI, double XKS0, double DXKSI ) ;




#endif
