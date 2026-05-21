#ifndef LGM_TSYG2001_H
#define LGM_TSYG2001_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/*
 *    Lgm_Tsyg2001.h
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
typedef struct _CB_T01S_TAIL        { double DXSHIFT1, DXSHIFT2, D, DELTADY; } _CB_T01S_TAIL;
typedef struct _CB_T01S_BIRKPAR     { double XKAPPA1, XKAPPA2;               } _CB_T01S_BIRKPAR;
typedef struct _CB_T01S_RCPAR       { double SC_SY, SC_AS, PHI;              } _CB_T01S_RCPAR;
typedef struct _CB_T01S_G           { double G;                              } _CB_T01S_G;
typedef struct _CB_T01S_RH0         { double RH0;                            } _CB_T01S_RH0;
typedef struct _CB_T01S_DPHI_B_RHO0 { double DPHI, B, RHO_0, XKAPPA;         } _CB_T01S_DPHI_B_RHO0;
typedef struct _CB_T01S_MODENUM     { int    M;                              } _CB_T01S_MODENUM;
typedef struct _CB_T01S_DTHETA      { double DTHETA;                         } _CB_T01S_DTHETA;



/*
 * Define a structure to hold all of the info needed in TS04
 */
typedef struct LgmTsyg2001_Info {

    /*
     * These are versions of the original fortran common blocks
     */
    _CB_T01S_TAIL             CB_TAIL;
    _CB_T01S_BIRKPAR          CB_BIRKPAR;
    _CB_T01S_RCPAR            CB_RCPAR;
    _CB_T01S_G                CB_G;
    _CB_T01S_RH0              CB_RH0;
    _CB_T01S_DPHI_B_RHO0      CB_DPHI_B_RHO0;
    _CB_T01S_MODENUM          CB_MODENUM;
    _CB_T01S_DTHETA           CB_DTHETA;


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


    int         DoneJ[4];
    double      CPS, SPS, S3PS, PST1[4], PST2[4], ST1[4], CT1[4], ST2[4], CT2[4], X1[4], Z1[4], X2[4], Z2[4];
    double      P[4][4], ooP[4][4], ooP2[4][4];
    double      Q[4][4], ooQ[4][4], ooQ2[4][4];
    double      R[4][4], ooR[4][4], ooR2[4][4];
    double      S[4][4], ooS[4][4], ooS2[4][4];
    double      YooP[4][4], YooQ[4][4], CYPI[4][4], CYQI[4][4], SYPI[4][4], SYQI[4][4];
    double      Z1ooR[4][4], Z2ooS[4][4], SZRK[4][4], CZSK[4][4], CZRK[4][4], SZSK[4][4];
    double      SQPR[4][4][4], SQQS[4][4][4], EPR[4][4][4], EQS[4][4][4];




} LgmTsyg2001_Info;







/*
 *  Function declarations for T01S (aka TSK03)
 */
void T01S_EXTALL( int IOPGEN, int IOPT, int IOPB, int IOPR, double *A, int NTOT, double PDYN, double DST, double BYIMF,
                    double BZIMF, double G1, double G2, double G3, double PS, double X, double Y, double Z,
                    double *BXCF, double *BYCF, double *BZCF, double *BXT1, double *BYT1, double *BZT1,
                    double *BXT2, double *BYT2, double *BZT2, double *BXSRC, double *BYSRC, double *BZSRC,
                    double *BXPRC, double *BYPRC, double *BZPRC,  double *BXR11, double *BYR11, double *BZR11,
                    double *BXR12, double *BYR12, double *BZR12, double *BXR21, double *BYR21, double *BZR21,
                    double *BXR22, double *BYR22, double *BZR22, double *HXIMF, double *HYIMF, double *HZIMF,
                    double *BX, double *BY, double *BZ, LgmTsyg2001_Info *t );
void    T01S_SHLCAR3X3( double X, double Y, double Z, double PS, double *BX, double *BY, double *BZ );
void    T01S_DEFORMED( int IOPT, double PS, double X, double Y, double Z,
                double *BX1, double *BY1, double *BZ1, double *BX2, double *BY2, double *BZ2, LgmTsyg2001_Info *t );
void    T01S_WARPED( int IOPT, double PS, double X, double Y, double Z,
                double *BX1, double *BY1, double *BZ1, double *BX2, double *BY2, double *BZ2, LgmTsyg2001_Info *t );
void    T01S_UNWARPED( int IOPT, double X, double Y, double Z,
                double *BX1, double *BY1, double *BZ1, double *BX2, double *BY2, double *BZ2, LgmTsyg2001_Info *t );
void    T01S_TAILDISK( double D0, double DELTADX, double DELTADY, double X, double Y, double Z, double *BX, double *BY, double *BZ);
void    T01S_SHLCAR5X5( double *A, double X, double Y, double Z, double DSHIFT, double *HX, double *HY, double *HZ );
void    T01S_BIRK_TOT( int IOPB, double PS, double X, double Y, double Z,
                double *BX11, double *BY11, double *BZ11, double *BX12, double *BY12, double *BZ12,
                double *BX21, double *BY21, double *BZ21, double *BX22, double *BY22, double *BZ22, LgmTsyg2001_Info *t );
void    T01S_BIRK_1N2( int NUMB, int MODE, double PS, double X, double Y, double Z, double *BX, double *BY, double *BZ, LgmTsyg2001_Info *t  ) ;
void    T01S_TWOCONES( double *A, double X, double Y, double Z, double *BX, double *BY, double *BZ, LgmTsyg2001_Info *t );
void    T01S_ONE_CONE( double *A, double X, double Y, double Z, double *BX, double *BY, double *BZ, LgmTsyg2001_Info *t );
double  T01S_THETA_S( double *A, double R, double THETA);
double  T01S_R_S( double *A, double R, double THETA);
void    T01S_FIALCOS( double R, double THETA, double PHI, double *BTHETA, double *BPHI, int N, double THETA0, double DT);
void    T01S_BIRK_SHL( int, int, int, int, int, double *A, double PS, double X_SC, double X, double Y, double Z, double *BX, double *BY, double *BZ, LgmTsyg2001_Info *t );
void    T01S_FULL_RC( int IOPR, double PS, double X, double Y, double Z,
                double *BXSRC, double *BYSRC, double *BZSRC, double *BXPRC, double *BYPRC,  double *BZPRC, LgmTsyg2001_Info *t);
void    T01S_SRC_PRC( int IOPR, double SC_SY, double SC_PR, double PHI, double PS, double X, double Y, double Z,
                double *BXSRC, double *BYSRC, double *BZSRC, double *BXPRC, double *BYPRC, double *BZPRC, LgmTsyg2001_Info *t);
void    T01S_RC_SYMM( double X, double Y, double Z, double *BX, double *BY, double *BZ);
double  T01S_AP( double R, double SINT, double COST);
void    PT01S_RC_SYMM( double X,double Y,double Z, double *BX, double *BY, double *BZ );
double  T01S_T01S_APPRC( double R, double SINT, double COST);
void    T01S_PRC_QUAD( double X, double Y, double Z, double *BX, double *BY, double *BZ);
double  T01S_BR_PRC_Q( double R, double SINT, double COST );
double  T01S_BT_PRC_Q( double R, double SINT, double COST);
void    T01S_FFS( double A, double A0, double DA, double *F, double *FA, double *FS );
void    T01S_RC_SHIELD( double *A, double PS, double X_SC, double X, double Y, double Z, double *BX, double *BY, double *BZ, LgmTsyg2001_Info *t );
void    T01S_DIPOLE( double PS, double X, double Y, double Z, double *BX, double *BY, double *BZ, LgmTsyg2001_Info *t );



/*
 *  Function declarations for T02 (aka T01_01)
 */
void T02_EXTALL( int IOPGEN, int IOPT, int IOPB, int IOPR, double *A, int NTOT, double PDYN, double DST, double BYIMF,
                    double BZIMF, double G1, double G2, double PS, double X, double Y, double Z,
                    double *BXCF, double *BYCF, double *BZCF, double *BXT1, double *BYT1, double *BZT1,
                    double *BXT2, double *BYT2, double *BZT2, double *BXSRC, double *BYSRC, double *BZSRC,
                    double *BXPRC, double *BYPRC, double *BZPRC,  double *BXR11, double *BYR11, double *BZR11,
                    double *BXR12, double *BYR12, double *BZR12, double *BXR21, double *BYR21, double *BZR21,
                    double *BXR22, double *BYR22, double *BZR22, double *HXIMF, double *HYIMF, double *HZIMF,
                    double *BX, double *BY, double *BZ, LgmTsyg2001_Info *t );
void    T02_SHLCAR3X3( double X, double Y, double Z, double PS, double *BX, double *BY, double *BZ );
void    T02_DEFORMED( int IOPT, double PS, double X, double Y, double Z,
                double *BX1, double *BY1, double *BZ1, double *BX2, double *BY2, double *BZ2, LgmTsyg2001_Info *t );
void    T02_WARPED( int IOPT, double PS, double X, double Y, double Z,
                double *BX1, double *BY1, double *BZ1, double *BX2, double *BY2, double *BZ2, LgmTsyg2001_Info *t );
void    T02_UNWARPED( int IOPT, double X, double Y, double Z,
                double *BX1, double *BY1, double *BZ1, double *BX2, double *BY2, double *BZ2, LgmTsyg2001_Info *t );
void    T02_TAILDISK( double D0, double DELTADX, double DELTADY, double X, double Y, double Z, double *BX, double *BY, double *BZ);
void    T02_SHLCAR5X5( double *A, double X, double Y, double Z, double DSHIFT, double *HX, double *HY, double *HZ );
void    T02_BIRK_TOT( int IOPB, double PS, double X, double Y, double Z,
                double *BX11, double *BY11, double *BZ11, double *BX12, double *BY12, double *BZ12,
                double *BX21, double *BY21, double *BZ21, double *BX22, double *BY22, double *BZ22, LgmTsyg2001_Info *t );
void    T02_BIRK_1N2( int NUMB, int MODE, double PS, double X, double Y, double Z, double *BX, double *BY, double *BZ, LgmTsyg2001_Info *t  ) ;
void    T02_TWOCONES( double *A, double X, double Y, double Z, double *BX, double *BY, double *BZ, LgmTsyg2001_Info *t );
void    T02_ONE_CONE( double *A, double X, double Y, double Z, double *BX, double *BY, double *BZ, LgmTsyg2001_Info *t );
double  T02_THETA_S( double *A, double R, double THETA);
double  T02_R_S( double *A, double R, double THETA);
void    T02_FIALCOS( double R, double THETA, double PHI, double *BTHETA, double *BPHI, int N, double THETA0, double DT);
void    T02_BIRK_SHL( int, int, int, int, int, double *A, double PS, double X_SC, double X, double Y, double Z, double *BX, double *BY, double *BZ, LgmTsyg2001_Info *t );
void    T02_FULL_RC( int IOPR, double PS, double X, double Y, double Z,
                double *BXSRC, double *BYSRC, double *BZSRC, double *BXPRC, double *BYPRC,  double *BZPRC, LgmTsyg2001_Info *t);
void    T02_SRC_PRC( int IOPR, double SC_SY, double SC_PR, double PHI, double PS, double X, double Y, double Z,
                double *BXSRC, double *BYSRC, double *BZSRC, double *BXPRC, double *BYPRC, double *BZPRC, LgmTsyg2001_Info *t);
void    T02_RC_SYMM( double X, double Y, double Z, double *BX, double *BY, double *BZ);
double  T02_AP( double R, double SINT, double COST);
void    PT02_RC_SYMM( double X,double Y,double Z, double *BX, double *BY, double *BZ );
double  T02_APPRC( double R, double SINT, double COST);
void    T02_PRC_QUAD( double X, double Y, double Z, double *BX, double *BY, double *BZ);
double  T02_BR_PRC_Q( double R, double SINT, double COST );
double  T02_BT_PRC_Q( double R, double SINT, double COST);
void    T02_FFS( double A, double A0, double DA, double *F, double *FA, double *FS );
void    T02_RC_SHIELD( double *A, double PS, double X_SC, double X, double Y, double Z, double *BX, double *BY, double *BZ, LgmTsyg2001_Info *t );
void    T02_DIPOLE( double PS, double X, double Y, double Z, double *BX, double *BY, double *BZ, LgmTsyg2001_Info *t );




#endif
