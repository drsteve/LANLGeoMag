#ifndef LGM_TSYG2004_H
#define LGM_TSYG2004_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/*
 *    Lgm_Tsyg2004.h
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
typedef struct _CB_TAIL        { double DXSHIFT1, DXSHIFT2, D, DELTADY; } _CB_TAIL;
typedef struct _CB_BIRKPAR     { double XKAPPA1, XKAPPA2;               } _CB_BIRKPAR;
typedef struct _CB_RCPAR       { double SC_SY, SC_AS, PHI;              } _CB_RCPAR;
typedef struct _CB_G           { double G;                              } _CB_G;
typedef struct _CB_RH0         { double RH0;                            } _CB_RH0;
typedef struct _CB_DPHI_B_RHO0 { double DPHI, B, RHO_0, XKAPPA;         } _CB_DPHI_B_RHO0;
typedef struct _CB_MODENUM     { int    M;                              } _CB_MODENUM;
typedef struct _CB_DTHETA      { double DTHETA;                         } _CB_DTHETA;


/*
 * Define a structure to hold all of the info needed in TS04
 */
typedef struct LgmTsyg2004_Info {

    /*
     * These are versions of the original fortran common blocks
     */
    _CB_TAIL         CB_TAIL;
    _CB_BIRKPAR      CB_BIRKPAR;
    _CB_RCPAR        CB_RCPAR;
    _CB_G            CB_G;
    _CB_RH0          CB_RH0;
    _CB_DPHI_B_RHO0  CB_DPHI_B_RHO0;
    _CB_MODENUM      CB_MODENUM;
    _CB_DTHETA       CB_DTHETA;


    /*
     *  SIN(PS) and COS(PS) are calculated 6 times each. This is very wasteful
     *  especially since we have it already. Lets pass sin_psi_op and
     *  cos_psi_op as globals....
     */
    double    sin_psi_op, cos_psi_op;


    double      OLD_PS;
    double      OLD_X;
    double      OLD_Y;
    double      OLD_Z;


    // cache vars used in BIRK_SHL()
    int         DoneJ[4];
    double      CPS, SPS, S3PS, PST1[4], PST2[4], ST1[4], CT1[4], ST2[4], CT2[4], X1[4], Z1[4], X2[4], Z2[4];
    double      P[4][4];
    double      ooP[4][4], ooP2[4][4];
    double      Q[4][4], ooQ[4][4], ooQ2[4][4];
    double      R[4][4], ooR[4][4], ooR2[4][4];
    double      S[4][4], ooS[4][4], ooS2[4][4];
    double      YooP[4][4], YooQ[4][4], CYPI[4][4], CYQI[4][4], SYPI[4][4], SYQI[4][4];
    double      Z1ooR[4][4], Z2ooS[4][4], SZRK[4][4], CZSK[4][4], CZRK[4][4], SZSK[4][4];
    double      SQPR[4][4][4], SQQS[4][4][4], EPR[4][4][4], EQS[4][4][4];


















} LgmTsyg2004_Info;







/*
 *  Function declarations
 */
void Lgm_Init_TS04( LgmTsyg2004_Info *t );
void TS04_EXTERN( int IOPGEN, int IOPT, int IOPB, int IOPR, double *A, int NTOT, double PDYN, double DST, double BXIMF, double BYIMF,
                double BZIMF, double W1, double W2, double W3, double W4, double W5, double W6, double PS,
                double X, double Y, double Z, double *BXCF, double *BYCF, double *BZCF, double *BXT1, double *BYT1,
                double *BZT1, double *BXT2, double *BYT2, double *BZT2, double *BXSRC, double *BYSRC, double *BZSRC,
                double *BXPRC, double *BYPRC, double *BZPRC,  double *BXR11, double *BYR11, double *BZR11, double *BXR12,
                double *BYR12, double *BZR12, double *BXR21, double *BYR21, double *BZR21, double *BXR22, double *BYR22,
                double *BZR22, double *HXIMF, double *HYIMF, double *HZIMF, double *BBX, double *BBY, double *BBZ, LgmTsyg2004_Info *tInfo );
void    SHLCAR3X3( double X, double Y, double Z, double PS, double *BX, double *BY, double *BZ, LgmTsyg2004_Info *tInfo );
void    DEFORMED( int IOPT, double PS, double X, double Y, double Z,
                double *BX1, double *BY1, double *BZ1, double *BX2, double *BY2, double *BZ2, LgmTsyg2004_Info *tInfo );
void    WARPED( int IOPT, double PS, double X, double Y, double Z,
                double *BX1, double *BY1, double *BZ1, double *BX2, double *BY2, double *BZ2, LgmTsyg2004_Info *tInfo );
void    UNWARPED( int IOPT, double X, double Y, double Z,
                double *BX1, double *BY1, double *BZ1, double *BX2, double *BY2, double *BZ2, LgmTsyg2004_Info *tInfo );
void    TAILDISK( double D0, double DELTADX, double DELTADY, double X, double Y, double Z, double *BX, double *BY, double *BZ);
void    SHLCAR5X5( double *A, double X, double Y, double Z, double DSHIFT, double *HX, double *HY, double *HZ );
void    BIRK_TOT( int IOPB, double PS, double X, double Y, double Z,
                double *BX11, double *BY11, double *BZ11, double *BX12, double *BY12, double *BZ12,
                double *BX21, double *BY21, double *BZ21, double *BX22, double *BY22, double *BZ22, LgmTsyg2004_Info *tInfo );
void    BIRK_1N2( int NUMB, int MODE, double PS, double X, double Y, double Z, double *BX, double *BY, double *BZ, LgmTsyg2004_Info *tInfo ) ;
void    TWOCONES( double *A, double X, double Y, double Z, double *BX, double *BY, double *BZ, LgmTsyg2004_Info *tInfo );
void    ONE_CONE( double *A, double X, double Y, double Z, double *BX, double *BY, double *BZ, LgmTsyg2004_Info *tInfo );
double  THETA_S( double *A, double R, double Theta, double SinTheta, double Sin2Theta, double Sin3Theta);
double  R_S( double *A, double R, double CosTheta, double Cos2Theta );
void    FIALCOS( double R, double THETA, double SIN_THETA, double COS_THETA, double SIN_PHI, double COS_PHI, double *BTHETA, double *BPHI, int N, double THETA0, double DT);
void    BIRK_SHL( int, int, int, int, int, double *A, double PS, double X_SC, double X, double Y, double Z, double *BX, double *BY, double *BZ, LgmTsyg2004_Info *tInfo );
void    FULL_RC( int IOPR, double PS, double X, double Y, double Z,
                double *BXSRC, double *BYSRC, double *BZSRC, double *BXPRC, double *BYPRC,  double *BZPRC, LgmTsyg2004_Info *tInfo );
void    SRC_PRC( int IOPR, double SC_SY, double SC_PR, double PHI, double PS, double X, double Y, double Z,
                double *BXSRC, double *BYSRC, double *BZSRC, double *BXPRC, double *BYPRC, double *BZPRC, LgmTsyg2004_Info *tInfo );
void    RC_SYMM( double X, double Y, double Z, double *BX, double *BY, double *BZ);
double  AP( double R, double SINT, double COST);
void    PRC_SYMM( double X,double Y,double Z, double *BX, double *BY, double *BZ );
double  APPRC( double R, double SINT, double COST);
void    PRC_QUAD( double X, double Y, double Z, double *BX, double *BY, double *BZ);
double  BR_PRC_Q( double R, double SINT, double COST );
double  BT_PRC_Q( double R, double SINT, double COST);
void    FFS( double A, double A0, double DA, double *F, double *FA, double *FS );
void    RC_SHIELD( double *A, double PS, double X_SC, double X, double Y, double Z, double *BX, double *BY, double *BZ, LgmTsyg2004_Info *tInfo );
void    DIPOLE( double PS, double X, double Y, double Z, double *BX, double *BY, double *BZ, LgmTsyg2004_Info *tInfo );











#endif

