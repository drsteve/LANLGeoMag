#ifndef LGM_TSYG2007_H
#define LGM_TSYG2007_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/*
 *    Lgm_Tsyg2007.h
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
typedef struct _CB7_TAIL        { double DXSHIFT1, DXSHIFT2, D, DELTADY; } _CB7_TAIL;
typedef struct _CB7_BIRKPAR     { double XKAPPA1, XKAPPA2;               } _CB7_BIRKPAR;
typedef struct _CB7_G           { double G; double TW;                   } _CB7_G;
typedef struct _CB7_RH0         { double RH0;                            } _CB7_RH0;
typedef struct _CB7_DPHI_B_RHO0 { double DPHI, B, RHO_0, XKAPPA;         } _CB7_DPHI_B_RHO0;
typedef struct _CB7_MODENUM     { int    M;                              } _CB7_MODENUM;
typedef struct _CB7_DTHETA      { double DTHETA;                         } _CB7_DTHETA;


/*
 * Define a structure to hold all of the info needed in TS04
 */
typedef struct LgmTsyg2007_Info {

    /*
     * These are versions of the original fortran common blocks (TSS, TSO, TSE just are vars in the Info structure...)
     */
    _CB7_TAIL         CB_TAIL;
    _CB7_BIRKPAR      CB_BIRKPAR;
    _CB7_G            CB_G;
    _CB7_RH0          CB_RH0;
    _CB7_DPHI_B_RHO0  CB_DPHI_B_RHO0;
    _CB7_MODENUM      CB_MODENUM;
    _CB7_DTHETA       CB_DTHETA;


    int     ArraysAlloced;
    double  *A;
    double  **TSS;    //[81][6];
    double  ***TSO;   //[81][6][5];
    double  ***TSE;   //[81][6][5];
    double  Pdyn;



    double      OLD_PS;
    double      OLD_X;
    double      OLD_Y;
    double      OLD_Z;
    double      OLD_PDYN;


    // cache vars used in TS07D_BIRK_SHL()
    int         DoneJ[4];
    double      sin_psi_op, cos_psi_op;
    double      CPS, SPS, S3PS, PST1[4], PST2[4], ST1[4], CT1[4], ST2[4], CT2[4], X1[4], Z1[4], X2[4], Z2[4];
    double      P[4][4];
    double      ooP[4][4], ooP2[4][4];
    double      Q[4][4], ooQ[4][4], ooQ2[4][4];
    double      R[4][4], ooR[4][4], ooR2[4][4];
    double      S[4][4], ooS[4][4], ooS2[4][4];
    double      YooP[4][4], YooQ[4][4], CYPI[4][4], CYQI[4][4], SYPI[4][4], SYQI[4][4];
    double      Z1ooR[4][4], Z2ooS[4][4], SZRK[4][4], CZSK[4][4], CZRK[4][4], SZSK[4][4];
    double      SQPR[4][4][4], SQQS[4][4][4], EPR[4][4][4], EQS[4][4][4];
    double      XAPPA;

    // cache vars used in TS07D_BIRSH_SY()
    int         S_DoneJ[4];
    double      S_sin_psi_op, S_cos_psi_op;
    double      S_CPS, S_SPS, S_S3PS, S_PST1[4], S_PST2[4], S_ST1[4], S_CT1[4], S_ST2[4], S_CT2[4], S_X1[4], S_Z1[4], S_X2[4], S_Z2[4];
    double      S_P[4][4];
    double      S_ooP[4][4], S_ooP2[4][4];
    double      S_Q[4][4], S_ooQ[4][4], S_ooQ2[4][4];
    double      S_R[4][4], S_ooR[4][4], S_ooR2[4][4];
    double      S_S[4][4], S_ooS[4][4], S_ooS2[4][4];
    double      S_YooP[4][4], S_YooQ[4][4], S_CYPI[4][4], S_CYQI[4][4], S_SYPI[4][4], S_SYQI[4][4];
    double      S_Z1ooR[4][4], S_Z2ooS[4][4], S_SZRK[4][4], S_CZSK[4][4], S_CZRK[4][4], S_SZSK[4][4];
    double      S_SQPR[4][4][4], S_SQQS[4][4][4], S_EPR[4][4][4], S_EQS[4][4][4];
    double      S_XAPPA;



} LgmTsyg2007_Info;







/*
 *  Function declarations
 */
void Lgm_Init_TS07( LgmTsyg2007_Info *t );
void Lgm_SetCoeffs_TS07( long int Date, double UTC, LgmTsyg2007_Info *t );

void Tsyg_TS07( int IOPT, double *PARMOD, double PS, double SINPS, double COSPS, double X, double Y, double Z,
                double *BX, double *BY, double *BZ, LgmTsyg2007_Info *tInfo );

void TS07D_EXTERN( int IOPGEN, double *A, int NTOT, double PS, double PDYN, double X, double Y, double Z, double *BXCF, double *BYCF,                                     
        double *BZCF, double BXTS[], double BYTS[], double BZTS[], double BXTO[][5], double BYTO[][5], double BZTO[][5], double BXTE[][5], double BYTE[][5],              
        double BZTE[][5], double *BXR11, double *BYR11, double *BZR11, double *BXR12, double *BYR12, double *BZR12, double *BXR21a, double *BYR21a,                       
        double *BZR21a, double *BXR21s, double *BYR21s, double *BZR21s, double *BX, double *BY, double *BZ,  LgmTsyg2007_Info *tInfo );

void    T07D_SHLCAR3X3( double X, double Y, double Z, double PS, double *BX, double *BY, double *BZ, LgmTsyg2007_Info *tInfo );

void    TS07D_DEFORMED( double PS, double X, double Y, double Z,                                                                                                          
        double BXS[], double BYS[], double BZS[], double BXO[][5], double BYO[][5], double BZO[][5],                                                                      
        double BXE[][5], double BYE[][5], double BZE[][5], LgmTsyg2007_Info *tInfo );

void TS07D_WARPED( double PS, double X, double Y, double Z,                                                                                                               
    double BXS[], double BYS[], double BZS[], double BXO[][5], double BYO[][5], double BZO[][5],                                                                          
    double BXE[][5], double BYE[][5], double BZE[][5], LgmTsyg2007_Info *tInfo );

void    TS07D_UNWARPED( double X, double Y, double Z, double BXS[], double BYS[], double BZS[],
        double BXO[][5], double BYO[][5], double BZO[][5], double BXE[][5], double BYE[][5],
        double BZE[][5], LgmTsyg2007_Info *tInfo );

void TS07D_TAILSHT_S( int M, double X, double Y, double Z, double *BX, double *BY, double *BZ, LgmTsyg2007_Info *tInfo );
void TS07D_SHTBNORM_S( int K, double X, double Y, double Z, double *FX, double *FY, double *FZ, LgmTsyg2007_Info *tInfo );
void TS07D_TAILSHT_OE( int IEVO, int MK, int M, double X, double Y, double Z, double *BX, double *BY, double *BZ, LgmTsyg2007_Info *tInfo );
void TS07D_SHTBNORM_O( int K, int L, double X, double Y, double Z, double *FX, double *FY, double *FZ, LgmTsyg2007_Info *tInfo );
void TS07D_SHTBNORM_E( int K, int L, double X, double Y, double Z, double *FX, double *FY, double *FZ, LgmTsyg2007_Info *tInfo );
double  bessj0( double x );
double  bessj1( double x );
double  bessj( int n, double x );
void     TS07D_BIRK_TOT( double PS, double X, double Y, double Z,
        double *BX11, double *BY11, double *BZ11, double *BX12, double *BY12, double *BZ12,
        double *BX21, double *BY21, double *BZ21, double *BX22, double *BY22, double *BZ22, LgmTsyg2007_Info *tInfo );
void    TS07D_BIRK_1N2( int NUMB, int MODE, double PS, double X, double Y, double Z, double *BX, double *BY, double *BZ, LgmTsyg2007_Info *tInfo );
void     TS07D_TWOCONES( double *A, double X, double Y, double Z, double *BX, double *BY, double *BZ, LgmTsyg2007_Info *tInfo );
void    TS07D_ONE_CONE( double *A, double X, double Y, double Z, double *BX, double *BY, double *BZ, LgmTsyg2007_Info *tInfo );
double TS07D_R_S( double *A, double R, double CosTheta, double Cos2Theta );
double TS07D_THETA_S( double *A, double R, double Theta, double SinTheta, double Sin2Theta, double Sin3Theta);
void     TS07D_FIALCOS( double R, double THETA, double SIN_THETA, double COS_THETA, double SIN_PHI, double COS_PHI, double *BTHETA, double *BPHI, int N, double THETA0, double DT);
void    TS07D_BIRK_SHL( int J, int PSChanged, int XChanged, int YChanged, int ZChanged, double *A, double PS, double X_SC,
                                double X, double Y, double Z, double *BX, double *BY, double *BZ, LgmTsyg2007_Info *tInfo );
void     TS07D_BIRTOTSY( double PS, double X, double Y, double Z,
        double *BX11, double *BY11, double *BZ11, double *BX12, double *BY12, double *BZ12,
        double *BX21, double *BY21, double *BZ21, double *BX22, double *BY22, double *BZ22, LgmTsyg2007_Info *tInfo );
void    TS07D_BIR1N2SY( int NUMB, int MODE, double PS, double X, double Y, double Z, double *BX, double *BY, double *BZ, LgmTsyg2007_Info *tInfo );
void     TS07D_TWOCONSS( double *A, double X, double Y, double Z, double *BX, double *BY, double *BZ, LgmTsyg2007_Info *tInfo );
void    TS07D_BIRSH_SY( int J, int PSChanged, int XChanged, int YChanged, int ZChanged, double *A, double PS, double X_SC,
                                double X, double Y, double Z, double *BX, double *BY, double *BZ, LgmTsyg2007_Info *tInfo );





#endif

