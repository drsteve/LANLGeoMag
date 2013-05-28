#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Lgm/Lgm_Tsyg2001.h"

/*
 *
 *    Converted to C by Michael G. Henderson (mghenderson@lanl.gov) Sept 14, 2010.
 *
 *     RELEASE DATE OF THIS VERSION:   AUGUST 8, 2001.
 *
 *    LATEST MODIFICATIONS/BUGS REMOVED:  JUNE 24, 2006:  REPLACED COEFFICIENTS IN:
 *        (i) DATA statement in FUNCTION T01S_AP,
 *        (ii)  DATA C_SY statement in SUBROUTINE T01S_FULL_RC, and
 *        (iii) DATA A statement in SUBROUTINE T01_01.
 *    This correction was needed because of a bug found in the symmetric ring current module.
 *    Its impact is a minor (a few percent) change of the model field in the inner magnetosphere.
 *
 *  --------------------------------------------------------------------
 *   A DATA-BASED MODEL OF THE EXTERNAL (I.E., WITHOUT EARTH'S CONTRIBUTION) PART OF THE
 *   MAGNETOSPHERIC MAGNETIC FIELD, CALIBRATED BY
 *    (1) SOLAR WIND PRESSURE PDYN (NANOPASCALS),
 *    (2) DST (NANOTESLA),
 *    (3) BYIMF,
 *    (4) BZIMF (NANOTESLA)
 *    (5) G1-INDEX
 *    (6) G2-INDEX  (SEE TSYGANENKO [2001] FOR AN EXACT DEFINITION OF THESE TWO INDICES)
 *
 *   THESE INPUT PARAMETERS SHOULD BE PLACED IN THE FIRST 6 ELEMENTS
 *   OF THE ARRAY PARMOD(10).
 *
 *   THE REST OF THE INPUT VARIABLES ARE: THE GEOT01S_DIPOLE TILT ANGLE PS (RADIANS),
 *     AND   X,Y,Z -  GSM POSITION (RE)
 *
 *   IOPT  IS JUST A DUMMY INPUT PARAMETER, NECESSARY TO MAKE THIS SUBROUTINE
 *   COMPATIBLE WITH THE TRACING SOFTWARE PACKAGE (GEOPACK). IN THIS MODEL
 *   IT DOES NOT AFFECT THE OUTPUT FIELD.
 *
 *  ******************************************************************************************
 *  * ATTENTION:  THE MODEL IS BASED ON DATA TAKEN SUNWARD FROM X=-15Re, AND HENCE BECOMES   *
 *  *              INVALID AT LARGER TAILWARD DISTANCES //////                                  *
 *  ******************************************************************************************
 *
 *   OUTPUT:  GSM COMPONENTS OF THE EXTERNAL MAGNETIC FIELD (BX,BY,BZ, nanotesla)
 *            COMPUTED AS A SUM OF CONTRIBUTIONS FROM PRINCIPAL FIELD SOURCES
 *
 *  (C) Copr. 2001, Nikolai A. Tsyganenko, USRA, Code 690.2, NASA GSFC
 *      Greenbelt, MD 20771, USA
 *
 *                            REFERENCE:
 *
 *    N. A. Tsyganenko, A new data-based model of the near magnetosphere magnetic field:
 *       1. Mathematical structure.
 *       2. Parameterization and fitting to observations.
 *
 *             (submitted to JGR, July 2001)
 *
 *
 *  ----------------------------------------------------------------------
 *
 *
 *
 */

void Lgm_Init_T01S( LgmTsyg2001_Info *t ){

    int                 i, j;

    // Init some params
    t->OLD_PS = -9e99;
    t->OLD_X  = -9e99;
    t->OLD_Y  = -9e99;
    t->OLD_Z  = -9e99;
    for (i=0; i<4; i++ ){
        t->DoneJ[i] = 0;
        for (j=0; j<4; j++ ){
            t->P[i][j] = -9e99;
        }
    }

    return;

}



void Tsyg_T01S( int IOPT, double *PARMOD, double PS, double SINPS, double COSPS, double X, double Y, double Z, double *BX, double *BY, double *BZ, LgmTsyg2001_Info *t ) {

    double       PDYN, DST_AST, BYIMF=0.0, BZIMF=0.0, G1, G2;
//    double       BXIMF=0.0;
    double       PSS, XX, YY, ZZ, BXCF, BYCF, BZCF, BXT1, BYT1, BZT1, BXT2, BYT2, BZT2;
    double       BXSRC, BYSRC, BZSRC, BXPRC, BYPRC, BZPRC,  BXR11, BYR11, BZR11;
    double       BXR12, BYR12, BZR12, BXR21, BYR21, BZR21, BXR22, BYR22, BZR22, HXIMF;
    double       HYIMF, HZIMF, BBX, BBY, BBZ;
//    int           IOPGEN=0, IOPTT=0, IOPB=0, IOPR=0;
    static double A[] = { -9e99, 1.00000,2.47341,0.40791,0.30429,-0.10637,-0.89108,3.29350,
                          -0.05413,-0.00696,1.07869,-0.02314,-0.66173,-0.68018,-0.03246,
                           0.02681,0.28062,0.16535,-0.02939,0.02639,-0.24891,-0.08063,
                           0.08900,-0.02475,0.05887,0.57691,0.65256,-0.03230,2.24733,
                           4.10546,1.13665,0.05506,0.97669,0.21164,0.64594,1.12556,0.01389,
                           1.02978,0.02968,0.15821,9.00519,28.17582,1.35285,0.42279};

    t->sin_psi = SINPS;
    t->cos_psi = COSPS;

    PDYN    = PARMOD[1];
    DST_AST = PARMOD[2]*0.8 - 13.0*sqrt( PDYN );
    BYIMF   = PARMOD[3];
    BZIMF   = PARMOD[4];

    G1   = PARMOD[5];
    G2   = PARMOD[6];
    PSS  = PS;
    XX   = X;
    YY   = Y;
    ZZ   = Z;

    T01S_EXTALL( 0, 0, 0, 0, A, 43, PDYN, DST_AST, BYIMF, BZIMF, G1, G2,
         PSS, XX, YY, ZZ, &BXCF, &BYCF, &BZCF, &BXT1, &BYT1, &BZT1, &BXT2, &BYT2, &BZT2, 
         &BXSRC, &BYSRC, &BZSRC, &BXPRC, &BYPRC, &BZPRC,  &BXR11, &BYR11, &BZR11, 
         &BXR12, &BYR12, &BZR12, &BXR21, &BYR21, &BZR21, &BXR22, &BYR22, &BZR22, &HXIMF, 
         &HYIMF, &HZIMF, &BBX, &BBY, &BBZ, t );

//printf("C: BXCF, BYCF, BZCF = %le %le %le\n", BXCF, BYCF, BZCF );
//printf("C: BXT1, BYT1, BZT1 = %e %e %e\n", BXT1, BYT1, BZT1 );
//printf("C: BXT2, BYT2, BZT2 = %e %e %e\n", BXT2, BYT2, BZT2 );
//printf("C: BXSRC, BYSRC, BZSRC = %e %e %e\n", BXSRC, BYSRC, BZSRC );
//printf("C: BXPRC, BYPRC, BZPRC = %e %e %e\n", BXPRC, BYPRC, BZPRC );
//printf("C: BXR11, BYR11, BZR11 = %e %e %e\n", BXR11, BYR11, BZR11 );
//printf("C: BXR21, BYR21, BZR21 = %e %e %e\n", BXR21, BYR21, BZR21 );
//printf("C: BXR22, BYR22, BZR22 = %e %e %e\n\n", BXR22, BYR22, BZR22 );

    *BX = BBX;
    *BY = BBY;
    *BZ = BBZ;

    return;

}




/*
 * 
 *    IOPGEN - GENERAL OPTION FLAG:  IOPGEN=0 - CALCULATE TOTAL FIELD
 *                                   IOPGEN=1 - T01S_DIPOLE SHIELDING ONLY
 *                                   IOPGEN=2 - TAIL FIELD ONLY
 *                                   IOPGEN=3 - BIRKELAND FIELD ONLY
 *                                   IOPGEN=4 - RING CURRENT FIELD ONLY
 *                                   IOPGEN=5 - INTERCONNECTION FIELD ONLY
 * 
 *    IOPT -  TAIL FIELD FLAG:       IOPT=0  -  BOTH MODES
 *                                   IOPT=1  -  MODE 1 ONLY
 *                                   IOPT=2  -  MODE 2 ONLY
 * 
 *    IOPB -  BIRKELAND FIELD FLAG:  IOPB=0  -  ALL 4 TERMS
 *                                   IOPB=1  -  REGION 1, MODES 1 AND 2
 *                                   IOPB=2  -  REGION 2, MODES 1 AND 2
 * 
 *    IOPR -  RING CURRENT FLAG:     IOPR=0  -  BOTH SRC AND PRC
 *                                   IOPR=1  -  SRC ONLY
 *                                   IOPR=2  -  PRC ONLY
 * 
 */
void T01S_EXTALL( int IOPGEN, int IOPT, int IOPB, int IOPR, double *A, int NTOT, double PDYN, double DST, double BYIMF, 
                    double BZIMF, double VBIMF1, double VBIMF2, double PS, double X, double Y, double Z, 
                    double *BXCF, double *BYCF, double *BZCF, double *BXT1, double *BYT1, double *BZT1, 
                    double *BXT2, double *BYT2, double *BZT2, double *BXSRC, double *BYSRC, double *BZSRC, 
                    double *BXPRC, double *BYPRC, double *BZPRC,  double *BXR11, double *BYR11, double *BZR11, 
                    double *BXR12, double *BYR12, double *BZR12, double *BXR21, double *BYR21, double *BZR21, 
                    double *BXR22, double *BYR22, double *BZR22, double *HXIMF, double *HYIMF, double *HZIMF, 
                    double *BX, double *BY, double *BZ, LgmTsyg2001_Info *t ){

    int       done;
    double    BPERP, THETA, CT, ST, YS, ZS, st, STHETAH, a, aa, b, A_R12, A_R22;
    double    XAPPA, XAPPA3, SPS, X0, AM, S0, FACTIMF, OIMFX, OIMFY, OIMFZ, R, XSS, ZSS;
    double    XSOLD, ZSOLD, RH, SINPSAS, COSPSAS, DD, RHO2, ASQ;
    double    XMXM, AXX0, ARO, SIGMA, CFX, CFY, CFZ, ZNAM;
    double    DLP1, DLP2, TAMP1, TAMP2, A_SRC, A_PRC, A_R11, XX, YY, ZZ;
    double    A_R21, QX, QY, QZ, FINT, FEXT, BBX, BBY, BBZ;



    double    A0_A=34.586, A0_S0=1.1960, A0_X0=3.4397;     // SHUE ET AL. PARAMETERS
    double    DSIG=0.003, RH2=-5.2;


    XAPPA  = pow( 0.5*PDYN, A[39] );   //  NOW THIS IS A VARIABLE PARAMETER
    //t->CB_RH0.RH0  = 8.0; // always overwritten...?
    t->CB_RH0.RH0  = A[40]; // TAIL HINGING DISTANCE
    t->CB_G.G      = A[41]; // TAIL WARPING PARAMETER

    XAPPA3 = XAPPA*XAPPA*XAPPA;

    XX = X*XAPPA;
    YY = Y*XAPPA;
    ZZ = Z*XAPPA;

    SPS = sin( PS );

    X0 = A0_X0/XAPPA;
    AM = A0_A/XAPPA;
    S0 = A0_S0;

    BPERP = sqrt( BYIMF*BYIMF + BZIMF*BZIMF );


    /*
     * CALCULATE THE IMF CLOCK ANGLE:
     */
    if ( (BYIMF==0.0)&&(BZIMF==0.0)) {
        THETA = 0.0;
    } else {
        THETA = atan2( BYIMF, BZIMF );
        if ( THETA <= 0.0 ) THETA += 2.0*M_PI;
    }


    CT = cos(THETA);
    ST = sin(THETA);
    YS = Y*CT - Z*ST;
    ZS = Z*CT + Y*ST;

    st = sin( 0.5*THETA );
    STHETAH = st*st;

    /*
     *  CALCULATE "IMF" COMPONENTS OUTSIDE THE MAGNETOPAUSE LAYER (HENCE BEGIN WITH "O")
     *  THEY ARE NEEDED ONLY IF THE POINT (X,Y,Z) IS WITHIN THE TRANSITION MAGNETOPAUSE LAYER
     *  OR OUTSIDE THE MAGNETOSPHERE:
     *
     */

    FACTIMF = A[24] + A[25]*STHETAH;

    OIMFX = 0.0;
    OIMFY = BYIMF*FACTIMF;
    OIMFZ = BZIMF*FACTIMF;

    R   = sqrt( X*X + Y*Y + Z*Z );
    XSS = X;
    ZSS = Z;

    /*
     * BEGIN ITERATIVE SEARCH OF T01S_UNWARPED COORDS (TO FIND SIGMA)
     */
    done = FALSE;
    while ( !done ) {
        XSOLD = XSS;
        ZSOLD = ZSS;

        a       = ZSS/R;
        RH      = t->CB_RH0.RH0 + RH2*a*a;
        b       = R/RH;
        SINPSAS = SPS/pow(1.0 + b*b*b, 0.33333333333333);
        COSPSAS = sqrt( 1.0 - SINPSAS*SINPSAS);
        ZSS     = X*SINPSAS + Z*COSPSAS;
        XSS     = X*COSPSAS - Z*SINPSAS;
        DD      = fabs(XSS-XSOLD) + fabs(ZSS-ZSOLD);
        if ( DD <= 1e-6 ) done = TRUE;
    }


    RHO2 = Y*Y + ZSS*ZSS;
    ASQ  = AM*AM;
    XMXM = AM + XSS - X0;
    if ( XMXM < 0.0 ) XMXM = 0.0; // THE BOUNDARY IS A CYLINDER TAILWARD OF X=X0-AM
    AXX0  = XMXM*XMXM;
    ARO   = ASQ + RHO2;
    aa    = ARO + AXX0;
    SIGMA = sqrt( (ARO + AXX0 + sqrt( aa*aa-4.0*ASQ*AXX0 ))/(2.0*ASQ) );



    /*
     *   NOW, THERE ARE THREE POSSIBLE CASES:
     *    (1) INSIDE THE MAGNETOSPHERE   (SIGMA
     *    (2) IN THE BOUNDARY LAYER
     *    (3) OUTSIDE THE MAGNETOSPHERE AND B.LAYER
     *        FIRST OF ALL, CONSIDER THE CASES (1) AND (2):
     *
     * 
     */
    if ( SIGMA < S0+DSIG ) {  //CASES (1) OR (2); CALCULATE THE MODEL FIELD (WITH THE POTENTIAL "PENETRATED" INTERCONNECTION FIELD):

        if ( IOPGEN <= 1 ) {
            T01S_SHLCAR3X3( XX, YY, ZZ, PS, &CFX, &CFY, &CFZ );         //  T01S_DIPOLE SHIELDING FIELD
//printf("PS = %g   XAPPA3 = %g\n", PS, XAPPA3);

            *BXCF = CFX*XAPPA3;
            *BYCF = CFY*XAPPA3;
            *BZCF = CFZ*XAPPA3;
        } else {
            *BXCF = 0.0;
            *BYCF = 0.0;
            *BZCF = 0.0;
        }

        if ( (IOPGEN == 0) || (IOPGEN == 2) ) {
            t->CB_TAIL.DXSHIFT1 = A[26] + A[27]*VBIMF2;
            t->CB_TAIL.DXSHIFT2 = 0.0;
            t->CB_TAIL.D        = A[28];
            t->CB_TAIL.DELTADY  = A[29];
            T01S_DEFORMED( IOPT, PS, XX, YY, ZZ, BXT1, BYT1, BZT1, BXT2, BYT2, BZT2, t ); //  TAIL FIELD (THREE MODES)
        } else {
            *BXT1 = 0.0;
            *BYT1 = 0.0;
            *BZT1 = 0.0;
            *BXT2 = 0.0;
            *BYT2 = 0.0;
            *BZT2 = 0.0;
        }

        if ( (IOPGEN == 0) || (IOPGEN == 3) ) {
            t->CB_BIRKPAR.XKAPPA1 = A[35] + A[36]*VBIMF2;
            t->CB_BIRKPAR.XKAPPA2 = A[37] + A[38]*VBIMF2;
            T01S_BIRK_TOT( IOPB, PS, XX, YY, ZZ, BXR11, BYR11, BZR11, BXR12, BYR12, 
                                  BZR12, BXR21, BYR21, BZR21, BXR22, BYR22, BZR22, t  );    //   BIRKELAND FIELD (TWO MODES FOR R1 AND TWO MODES FOR R2)
        } else {
            *BXR11 = 0.0;
            *BYR11 = 0.0;
            *BZR11 = 0.0;
            *BXR12 = 0.0;
            *BYR12 = 0.0;
            *BZR12 = 0.0;
            *BXR21 = 0.0;
            *BYR21 = 0.0;
            *BZR21 = 0.0;
            *BXR22 = 0.0;
            *BYR22 = 0.0;
            *BZR22 = 0.0;
        }

        if ( (IOPGEN == 0) || (IOPGEN == 4) ) {
            ZNAM = fabs(DST);
            t->CB_RCPAR.PHI  = 1.5707963*tanh( ZNAM/A[34] );
            if ( ZNAM < 20.0 ) ZNAM = 20.0;
            a     = 20.0/ZNAM;
            t->CB_RCPAR.SC_SY = A[30]*pow(a, A[31])*XAPPA;
            t->CB_RCPAR.SC_AS = A[32]*pow(a, A[33])*XAPPA;
            T01S_FULL_RC(IOPR, PS, XX, YY, ZZ, BXSRC, BYSRC, BZSRC, BXPRC, BYPRC, BZPRC, t );  //  SHIELDED RING CURRENT (SRC AND PRC)
        } else {
            *BXSRC = 0.0;
            *BYSRC = 0.0;
            *BZSRC = 0.0;
            *BXPRC = 0.0;
            *BYPRC = 0.0;
            *BZPRC = 0.0;
        }

        if ( (IOPGEN == 0) || (IOPGEN == 5) ) {
            *HXIMF = 0.0;
            *HYIMF = BYIMF;
            *HZIMF = BZIMF;  /* THESE ARE COMPONENTS OF THE PENETRATED FIELD PER
                             * UNIT OF THE PENETRATION COEFFICIENT.  IN OTHER
                             * WORDS, THESE ARE DERIVATIVES OF THE PENETRATION
                             * FIELD COMPONENTS WITH RESPECT TO THE PENETRATION
                             * COEFFICIENT. WE ASSUME THAT ONLY THE TRANSVERSE
                             * COMPONENT OF THE FIELD PENETRATES INSIDE.
                             */

        } else { 
            *HXIMF = 0.0; 
            *HYIMF = 0.0; 
            *HZIMF = 0.0; 
        }


        /*
         *  -------------------------------------
         *  NOW, ADD UP ALL THE COMPONENTS:
         *
         */

        a    = 0.5*PDYN;
        DLP1 = pow( a, A[42] );
        DLP2 = pow( a, A[43] );

        a     = sqrt(PDYN);
printf("DST = %g\n", DST);
        TAMP1 = A[2]  + A[3]*DLP1 + A[4]*VBIMF1 + A[5]*DST;
        TAMP2 = A[6]  + A[7]*DLP2 + A[8]*VBIMF1 + A[9]*DST;
        A_SRC = A[10] + A[11]*DST + A[12]*a;
        A_PRC = A[13] + A[14]*DST + A[15]*a;
        A_R11 = A[16] + A[17]*VBIMF2;
        A_R12 = A[18] + A[19]*VBIMF2;
        A_R21 = A[20] + A[21]*VBIMF2;
        A_R22 = A[22] + A[23]*VBIMF2;

        BBX = A[1]* *BXCF + TAMP1* *BXT1 + TAMP2* *BXT2 + A_SRC* *BXSRC + A_PRC* *BXPRC
                 + A_R11* *BXR11 + A_R12* *BXR12 + A_R21* *BXR21 + A_R22* *BXR22
                 + A[24]* *HXIMF + A[25]* *HXIMF*STHETAH;

        BBY = A[1]* *BYCF + TAMP1* *BYT1 + TAMP2* *BYT2 + A_SRC* *BYSRC + A_PRC* *BYPRC
                 + A_R11* *BYR11 + A_R12* *BYR12 + A_R21* *BYR21 + A_R22* *BYR22
                 + A[24]* *HYIMF + A[25]* *HYIMF*STHETAH;

        BBZ = A[1]* *BZCF + TAMP1* *BZT1 + TAMP2* *BZT2 + A_SRC* *BZSRC + A_PRC* *BZPRC
                 + A_R11* *BZR11 + A_R12* *BZR12 + A_R21* *BZR21 + A_R22* *BZR22
                 + A[24]* *HZIMF + A[25]* *HZIMF*STHETAH;


        /*
         *   AND WE HAVE THE TOTAL EXTERNAL FIELD.
         *   NOW, LET US CHECK WHETHER WE HAVE THE CASE (1). IF YES - WE ARE DONE:
         *
         */
        if ( SIGMA < S0-DSIG ) {    //  (X,Y,Z) IS INSIDE THE MAGNETOSPHERE
            *BX = BBX;
            *BY = BBY;
            *BZ = BBZ;
        } else {                    //  THIS IS THE MOST COMPLEX CASE: WE ARE INSIDE THE INTERPOLATION REGION
            a = 0.5*(SIGMA-S0)/DSIG;
            FINT = 0.5 - a;
            FEXT = 0.5 + a;
    
            T01S_DIPOLE( PS, X, Y, Z, &QX, &QY, &QZ, t );
            *BX = (BBX+QX)*FINT + OIMFX*FEXT - QX;
            *BY = (BBY+QY)*FINT + OIMFY*FEXT - QY;
            *BZ = (BBZ+QZ)*FINT + OIMFZ*FEXT - QZ;
        } // THE CASES (1) AND (2) ARE EXHAUSTED; THE ONLY REMAINING POSSIBILITY IS NOW THE CASE (3):

    } else {
        T01S_DIPOLE( PS, X, Y, Z, &QX, &QY, &QZ, t );
        *BX = OIMFX - QX;
        *BY = OIMFY - QY;
        *BZ = OIMFZ - QZ;
    }

}



/*
 *$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
 *
 *       SUBROUTINE  T01S_SHLCAR3X3(X,Y,Z,PS,BX,BY,BZ)
 *
 *   THIS S/R RETURNS THE SHIELDING FIELD FOR THE EARTH'S T01S_DIPOLE,
 *   REPRESENTED BY  2x3x3=18 "CARTESIAN" HARMONICS, tilted with respect
 *   to the z=0 plane  (NB#4, p.74)
 *
 * - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 *  The 36 coefficients enter in pairs in the amplitudes of the "cartesian"
 *    harmonics (A(1)-A(36).
 *  The 14 nonlinear parameters (A(37)-A(50) are the scales Pi,Ri,Qi,and Si
 *   entering the arguments of exponents, sines, and cosines in each of the
 *   18 "Cartesian" harmonics  PLUS TWO TILT ANGLES FOR THE CARTESIAN HARMONICS
 *       (ONE FOR THE PSI=0 MODE AND ANOTHER FOR THE PSI=90 MODE)
 * - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 */
void T01S_SHLCAR3X3( double X, double Y, double Z, double PS, double *BX, double *BY, double *BZ ) {

    static double A[] = { -9e99, -901.2327248,895.8011176,817.6208321,-845.5880889,
                          -83.73539535,86.58542841,336.8781402,-329.3619944,-311.2947120,
                          308.6011161,31.94469304,-31.30824526,125.8739681,-372.3384278,
                          -235.4720434,286.7594095,21.86305585,-27.42344605,-150.4874688,
                          2.669338538,1.395023949,-.5540427503,-56.85224007,3.681827033,
                          -43.48705106,5.103131905,1.073551279,-.6673083508,12.21404266,
                          4.177465543,5.799964188,-.3977802319,-1.044652977,.5703560010,
                          3.536082962,-3.222069852,9.620648151,6.082014949,27.75216226,
                          12.44199571,5.122226936,6.982039615,20.12149582,6.150973118,
                          4.663639687,15.73319647,2.303504968,5.840511214,.8385953499e-01, 
                          .3477844929 };

    double  P1, P2, P3, ooP1, ooP2, ooP3, ooP1_2, ooP2_2, ooP3_2;
    double  R1, R2, R3, ooR1, ooR2, ooR3, ooR1_2, ooR2_2, ooR3_2;
    double  Q1, Q2, Q3, ooQ1, ooQ2, ooQ3, ooQ1_2, ooQ2_2, ooQ3_2;
    double  S1, S2, S3, ooS1, ooS2, ooS3, ooS1_2, ooS2_2, ooS3_2;
    double  T1, T2;
    double  CYP, SYP, CZR, SZR, CYQ, SYQ, CZS, SZS;
    double  sqrtP1R1, sqrtP1R2, sqrtP1R3, sqrtP2R1, sqrtP2R2, sqrtP2R3, sqrtP3R1, sqrtP3R2, sqrtP3R3;
    double  sqrtQ1S1, sqrtQ1S2, sqrtQ1S3, sqrtQ2S1, sqrtQ2S2, sqrtQ2S3, sqrtQ3S1, sqrtQ3S2, sqrtQ3S3;
    double  YoP1, YoP2, YoP3, Z1oR1, Z1oR2, Z1oR3;
    double  YoQ1, YoQ2, YoQ3, Z2oS1, Z2oS2, Z2oS3;
    double  COSZ1oR1, COSZ1oR2, COSZ1oR3, SINZ1oR1, SINZ1oR2, SINZ1oR3;
    double  COSZ2oS1, COSZ2oS2, COSZ2oS3, SINZ2oS1, SINZ2oS2, SINZ2oS3;
    double  COSYoP1, COSYoP2, COSYoP3, SINYoP1, SINYoP2, SINYoP3;
    double  COSYoQ1, COSYoQ2, COSYoQ3, SINYoQ1, SINYoQ2, SINYoQ3;
    double  X1, X2, Z1, Z2, SQPR, EXPR, SQQS, EXQS, CPS, SPS, S2PS;
    double  ST1, ST2, CT1, CT2;
    double  FX1, HY1, FZ1, HX1, HZ1, FX2, HY2, FZ2, HX2, HZ2, FX3, HY3, FZ3, HX3, HZ3;
    double  FX4, HY4, FZ4, HX4, HZ4, FX5, HY5, FZ5, HX5, HZ5, FX6, HY6, FZ6, HX6, HZ6;
    double  FX7, HY7, FZ7, HX7, HZ7, FX8, HY8, FZ8, HX8, HZ8, FX9, HY9, FZ9, HX9, HZ9;
    double  A1, A2, A3, A4, A5, A6, A7, A8, A9;


    P1 = A[37]; P2 = A[38]; P3 = A[39];
    R1 = A[40]; R2 = A[41]; R3 = A[42];
    Q1 = A[43]; Q2 = A[44]; Q3 = A[45];
    S1 = A[46]; S2 = A[47]; S3 = A[48];
    T1 = A[49]; T2 = A[50];


    CPS  = cos(PS); SPS  = sin(PS); S2PS = 2.0*CPS; 

    ST1 = sin( PS*T1 ); CT1 = cos( PS*T1 );
    ST2 = sin( PS*T2 ); CT2 = cos( PS*T2 );

    X1 = X*CT1 - Z*ST1; Z1 = X*ST1 + Z*CT1;
    X2 = X*CT2 - Z*ST2; Z2 = X*ST2 + Z*CT2;

    ooP1 = 1.0/P1; ooP1_2 = ooP1*ooP1;
    ooP2 = 1.0/P2; ooP2_2 = ooP2*ooP2;
    ooP3 = 1.0/P3; ooP3_2 = ooP3*ooP3;
    ooR1 = 1.0/R1; ooR1_2 = ooR1*ooR1;
    ooR2 = 1.0/R2; ooR2_2 = ooR2*ooR2;
    ooR3 = 1.0/R3; ooR3_2 = ooR3*ooR3;

    sqrtP1R1 = sqrt(ooP1_2 + ooR1_2);
    sqrtP1R2 = sqrt(ooP1_2 + ooR2_2);
    sqrtP1R3 = sqrt(ooP1_2 + ooR3_2);

    sqrtP2R1 = sqrt(ooP2_2 + ooR1_2);
    sqrtP2R2 = sqrt(ooP2_2 + ooR2_2);
    sqrtP2R3 = sqrt(ooP2_2 + ooR3_2);

    sqrtP3R1 = sqrt(ooP3_2 + ooR1_2);
    sqrtP3R2 = sqrt(ooP3_2 + ooR2_2);
    sqrtP3R3 = sqrt(ooP3_2 + ooR3_2);

    YoP1 = Y*ooP1;
    YoP2 = Y*ooP2;
    YoP3 = Y*ooP3;

    Z1oR1 = Z1*ooR1;
    Z1oR2 = Z1*ooR2;
    Z1oR3 = Z1*ooR3;

    COSZ1oR1 = cos(Z1oR1); SINZ1oR1 = sin(Z1oR1);
    COSZ1oR2 = cos(Z1oR2); SINZ1oR2 = sin(Z1oR2);
    COSZ1oR3 = cos(Z1oR3); SINZ1oR3 = sin(Z1oR3);

    COSYoP1 = cos(YoP1); SINYoP1 = sin(YoP1);
    COSYoP2 = cos(YoP2); SINYoP2 = sin(YoP2);
    COSYoP3 = cos(YoP3); SINYoP3 = sin(YoP3);


    /*
     *
     *  MAKE THE TERMS IN THE 1ST SUM ("PERPENDICULAR" SYMMETRY):
     *
     *       I=1:
     */
    SQPR = sqrtP1R1;
    CYP  = COSYoP1;
    SYP  = SINYoP1;
    CZR  = COSZ1oR1;
    SZR  = SINZ1oR1;
    EXPR = exp(SQPR*X1);
    FX1  = -SQPR*EXPR*CYP*SZR;
    HY1  = EXPR/P1*SYP*SZR;
    FZ1  = -EXPR*CYP/R1*CZR;
    HX1  = FX1*CT1 + FZ1*ST1;
    HZ1  = -FX1*ST1 + FZ1*CT1;

    SQPR = sqrtP1R2;
    CYP  = COSYoP1;
    SYP  = SINYoP1;
    CZR  = COSZ1oR2;
    SZR  = SINZ1oR2;
    EXPR = exp(SQPR*X1);
    FX2  = -SQPR*EXPR*CYP*SZR;
    HY2  = EXPR/P1*SYP*SZR;
    FZ2  = -EXPR*CYP/R2*CZR;
    HX2  = FX2*CT1+FZ2*ST1;
    HZ2  = -FX2*ST1+FZ2*CT1;

    SQPR = sqrtP1R3;
    CYP  = COSYoP1;
    SYP  = SINYoP1;
    CZR  = COSZ1oR3;
    SZR  = SINZ1oR3;
    EXPR = exp(SQPR*X1);
    FX3  = -EXPR*CYP*(SQPR*Z1*CZR+SZR/R3*(X1+1.0/SQPR));
    HY3  = EXPR/P1*SYP*(Z1*CZR+X1/R3*SZR/SQPR);
    FZ3  = -EXPR*CYP*(CZR*(1.0+X1/(R3*R3)/SQPR)-Z1/R3*SZR);
    HX3  = FX3*CT1+FZ3*ST1;
    HZ3  = -FX3*ST1+FZ3*CT1;


    /*
     *       I=2:
     */
    SQPR =  sqrtP2R1;
    CYP  =  COSYoP2;
    SYP  =  SINYoP2;
    CZR  =  COSZ1oR1;
    SZR  =  SINZ1oR1;
    EXPR =  exp(SQPR*X1);
    FX4  = -SQPR*EXPR*CYP*SZR;
    HY4  =  EXPR/P2*SYP*SZR;
    FZ4  = -EXPR*CYP/R1*CZR;
    HX4  =  FX4*CT1+FZ4*ST1;
    HZ4  = -FX4*ST1+FZ4*CT1;

    SQPR =  sqrtP2R2;
    CYP  =  COSYoP2;
    SYP  =  SINYoP2;
    CZR  =  COSZ1oR2;
    SZR  =  SINZ1oR2;
    EXPR =  exp(SQPR*X1);
    FX5  = -SQPR*EXPR*CYP*SZR;
    HY5  =  EXPR/P2*SYP*SZR;
    FZ5  = -EXPR*CYP/R2*CZR;
    HX5  =  FX5*CT1+FZ5*ST1;
    HZ5  = -FX5*ST1+FZ5*CT1;

    SQPR =  sqrtP2R3;
    CYP  =  COSYoP2;
    SYP  =  SINYoP2;
    CZR  =  COSZ1oR3;
    SZR  =  SINZ1oR3;
    EXPR =  exp(SQPR*X1);
    FX6  = -EXPR*CYP*(SQPR*Z1*CZR+SZR/R3*(X1+1.0/SQPR));
    HY6  =  EXPR/P2*SYP*(Z1*CZR+X1/R3*SZR/SQPR);
    FZ6  = -EXPR*CYP*(CZR*(1.0+X1/(R3*R3)/SQPR)-Z1/R3*SZR);
    HX6  =  FX6*CT1+FZ6*ST1;
    HZ6  = -FX6*ST1+FZ6*CT1;



    /*
     *      I=3:
     */
    SQPR =  sqrtP3R1;
    CYP  =  COSYoP3;
    SYP  =  SINYoP3;
    CZR  =  COSZ1oR1;
    SZR  =  SINZ1oR1;
    EXPR =  exp(SQPR*X1);
    FX7  = -SQPR*EXPR*CYP*SZR;
    HY7  =  EXPR/P3*SYP*SZR;
    FZ7  = -EXPR*CYP/R1*CZR;
    HX7  =  FX7*CT1+FZ7*ST1;
    HZ7  = -FX7*ST1+FZ7*CT1;

    SQPR =  sqrtP3R2;
    CYP  =  COSYoP3;
    SYP  =  SINYoP3;
    CZR  =  COSZ1oR2;
    SZR  =  SINZ1oR2;
    EXPR =  exp(SQPR*X1);
    FX8  = -SQPR*EXPR*CYP*SZR;
    HY8  =  EXPR/P3*SYP*SZR;
    FZ8  = -EXPR*CYP/R2*CZR;
    HX8  =  FX8*CT1+FZ8*ST1;
    HZ8  = -FX8*ST1+FZ8*CT1;

    SQPR =  sqrtP3R3;
    CYP  =  COSYoP3;
    SYP  =  SINYoP3;
    CZR  =  COSZ1oR3;
    SZR  =  SINZ1oR3;
    EXPR =  exp(SQPR*X1);
    FX9  = -EXPR*CYP*(SQPR*Z1*CZR+SZR/R3*(X1+1.0/SQPR));
    HY9  =  EXPR/P3*SYP*(Z1*CZR+X1/R3*SZR/SQPR);
    FZ9  = -EXPR*CYP*(CZR*(1.0+X1/(R3*R3)/SQPR)-Z1/R3*SZR);
    HX9  =  FX9*CT1+FZ9*ST1;
    HZ9  = -FX9*ST1+FZ9*CT1;

    A1 = A[1]  + A[2]*CPS;
    A2 = A[3]  + A[4]*CPS;
    A3 = A[5]  + A[6]*CPS;
    A4 = A[7]  + A[8]*CPS;
    A5 = A[9]  + A[10]*CPS;
    A6 = A[11] + A[12]*CPS;
    A7 = A[13] + A[14]*CPS;
    A8 = A[15] + A[16]*CPS;
    A9 = A[17] + A[18]*CPS;

    *BX = A1*HX1 + A2*HX2 + A3*HX3 + A4*HX4 + A5*HX5 + A6*HX6 + A7*HX7 + A8*HX8 + A9*HX9;
    *BY = A1*HY1 + A2*HY2 + A3*HY3 + A4*HY4 + A5*HY5 + A6*HY6 + A7*HY7 + A8*HY8 + A9*HY9;
    *BZ = A1*HZ1 + A2*HZ2 + A3*HZ3 + A4*HZ4 + A5*HZ5 + A6*HZ6 + A7*HZ7 + A8*HZ8 + A9*HZ9;


    /*  
     *  MAKE THE TERMS IN THE 2ND SUM ("PARALLEL" SYMMETRY):
     */
    ooQ1 = 1.0/Q1; ooQ1_2 = ooQ1*ooQ1;
    ooQ2 = 1.0/Q2; ooQ2_2 = ooQ2*ooQ2;
    ooQ3 = 1.0/Q3; ooQ3_2 = ooQ3*ooQ3;
    ooS1 = 1.0/S1; ooS1_2 = ooS1*ooS1;
    ooS2 = 1.0/S2; ooS2_2 = ooS2*ooS2;
    ooS3 = 1.0/S3; ooS3_2 = ooS3*ooS3;

    sqrtQ1S1 = sqrt(ooQ1_2 + ooS1_2);
    sqrtQ1S2 = sqrt(ooQ1_2 + ooS2_2);
    sqrtQ1S3 = sqrt(ooQ1_2 + ooS3_2);

    sqrtQ2S1 = sqrt(ooQ2_2 + ooS1_2);
    sqrtQ2S2 = sqrt(ooQ2_2 + ooS2_2);
    sqrtQ2S3 = sqrt(ooQ2_2 + ooS3_2);

    sqrtQ3S1 = sqrt(ooQ3_2 + ooS1_2);
    sqrtQ3S2 = sqrt(ooQ3_2 + ooS2_2);
    sqrtQ3S3 = sqrt(ooQ3_2 + ooS3_2);


    YoQ1 = Y*ooQ1;
    YoQ2 = Y*ooQ2;
    YoQ3 = Y*ooQ3;

    Z2oS1 = Z2*ooS1;
    Z2oS2 = Z2*ooS2;
    Z2oS3 = Z2*ooS3;

    COSZ2oS1 = cos(Z2oS1); SINZ2oS1 = sin(Z2oS1);
    COSZ2oS2 = cos(Z2oS2); SINZ2oS2 = sin(Z2oS2);
    COSZ2oS3 = cos(Z2oS3); SINZ2oS3 = sin(Z2oS3);

    COSYoQ1 = cos(YoQ1); SINYoQ1 = sin(YoQ1);
    COSYoQ2 = cos(YoQ2); SINYoQ2 = sin(YoQ2);
    COSYoQ3 = cos(YoQ3); SINYoQ3 = sin(YoQ3);

    // I = 1
    SQQS =  sqrtQ1S1; 
    CYQ  =  COSYoQ1;
    SYQ  =  SINYoQ1;
    CZS  =  COSZ2oS1;
    SZS  =  SINZ2oS1;
    EXQS =  exp(SQQS*X2);
    FX1  = -SQQS*EXQS*CYQ*CZS*SPS;
    HY1  =  EXQS/Q1*SYQ*CZS*SPS;
    FZ1  =  EXQS*CYQ/S1*SZS*SPS;
    HX1  =  FX1*CT2+FZ1*ST2;
    HZ1  = -FX1*ST2+FZ1*CT2;

    SQQS =  sqrtQ1S2;
    CYQ  =  COSYoQ1;
    SYQ  =  SINYoQ1;
    CZS  =  COSZ2oS2;
    SZS  =  SINZ2oS2;
    EXQS =  exp(SQQS*X2);
    FX2  = -SQQS*EXQS*CYQ*CZS*SPS;
    HY2  =  EXQS/Q1*SYQ*CZS*SPS;
    FZ2  =  EXQS*CYQ/S2*SZS*SPS;
    HX2  =  FX2*CT2+FZ2*ST2;
    HZ2  = -FX2*ST2+FZ2*CT2;

    SQQS =  sqrtQ1S3;
    CYQ  =  COSYoQ1;
    SYQ  =  SINYoQ1;
    CZS  =  COSZ2oS3;
    SZS  =  SINZ2oS3;
    EXQS =  exp(SQQS*X2);
    FX3  = -SQQS*EXQS*CYQ*CZS*SPS;
    HY3  =  EXQS/Q1*SYQ*CZS*SPS;
    FZ3  =  EXQS*CYQ/S3*SZS*SPS;
    HX3  =  FX3*CT2+FZ3*ST2;
    HZ3  = -FX3*ST2+FZ3*CT2;



    /*
     *       I=2:
     */
    SQQS =  sqrtQ2S1;
    CYQ  =  COSYoQ2;
    SYQ  =  SINYoQ2;
    CZS  =  COSZ2oS1;
    SZS  =  SINZ2oS1;
    EXQS =  exp(SQQS*X2);
    FX4  = -SQQS*EXQS*CYQ*CZS *SPS;
    HY4  =  EXQS/Q2*SYQ*CZS   *SPS;
    FZ4  =  EXQS*CYQ/S1*SZS   *SPS;
    HX4  =  FX4*CT2+FZ4*ST2;
    HZ4  = -FX4*ST2+FZ4*CT2;

    SQQS =  sqrtQ2S2;
    CYQ  =  COSYoQ2;
    SYQ  =  SINYoQ2;
    CZS  =  COSZ2oS2;
    SZS  =  SINZ2oS2;
    EXQS =  exp(SQQS*X2);
    FX5  = -SQQS*EXQS*CYQ*CZS *SPS;
    HY5  =  EXQS/Q2*SYQ*CZS   *SPS;
    FZ5  =  EXQS*CYQ/S2*SZS   *SPS;
    HX5  =  FX5*CT2+FZ5*ST2;
    HZ5  = -FX5*ST2+FZ5*CT2;

    SQQS =  sqrtQ2S3;
    CYQ  =  COSYoQ2;
    SYQ  =  SINYoQ2;
    CZS  =  COSZ2oS3;
    SZS  =  SINZ2oS3;
    EXQS =  exp(SQQS*X2);
    FX6  = -SQQS*EXQS*CYQ*CZS *SPS;
    HY6  =  EXQS/Q2*SYQ*CZS   *SPS;
    FZ6  =  EXQS*CYQ/S3*SZS   *SPS;
    HX6  =  FX6*CT2+FZ6*ST2;
    HZ6  = -FX6*ST2+FZ6*CT2;



    /*
     *       I=3:
     */
    SQQS =  sqrtQ3S1;
    CYQ  =  COSYoQ3;
    SYQ  =  SINYoQ3;
    CZS  =  COSZ2oS1;
    SZS  =  SINZ2oS1;
    EXQS =  exp(SQQS*X2);
    FX7  = -SQQS*EXQS*CYQ*CZS*SPS;
    HY7  =  EXQS/Q3*SYQ*CZS*SPS;
    FZ7  =  EXQS*CYQ/S1*SZS*SPS;
    HX7  =  FX7*CT2+FZ7*ST2;
    HZ7  = -FX7*ST2+FZ7*CT2;

    SQQS =  sqrtQ3S2;
    CYQ  =  COSYoQ3; 
    SYQ  =  SINYoQ3; 
    CZS  =  COSZ2oS2;
    SZS  =  SINZ2oS2;
    EXQS =  exp(SQQS*X2);
    FX8  = -SQQS*EXQS*CYQ*CZS*SPS;
    HY8  =  EXQS/Q3*SYQ*CZS*SPS;
    FZ8  =  EXQS*CYQ/S2*SZS*SPS;
    HX8  =  FX8*CT2+FZ8*ST2;
    HZ8  = -FX8*ST2+FZ8*CT2;

    SQQS =  sqrtQ3S3;
    CYQ  =  COSYoQ3; 
    SYQ  =  SINYoQ3; 
    CZS  =  COSZ2oS3;
    SZS  =  SINZ2oS3;
    EXQS =  exp(SQQS*X2);
    FX9  = -SQQS*EXQS*CYQ*CZS*SPS;
    HY9  =  EXQS/Q3*SYQ*CZS*SPS;
    FZ9  =  EXQS*CYQ/S3*SZS*SPS;
    HX9  =  FX9*CT2+FZ9*ST2;
    HZ9  = -FX9*ST2+FZ9*CT2;

    A1 = A[19] + A[20]*S2PS;
    A2 = A[21] + A[22]*S2PS;
    A3 = A[23] + A[24]*S2PS;
    A4 = A[25] + A[26]*S2PS;
    A5 = A[27] + A[28]*S2PS;
    A6 = A[29] + A[30]*S2PS;
    A7 = A[31] + A[32]*S2PS;
    A8 = A[33] + A[34]*S2PS;
    A9 = A[35] + A[36]*S2PS;

    *BX += A1*HX1 + A2*HX2 + A3*HX3 + A4*HX4 + A5*HX5 + A6*HX6 + A7*HX7 + A8*HX8 + A9*HX9;
    *BY += A1*HY1 + A2*HY2 + A3*HY3 + A4*HY4 + A5*HY5 + A6*HY6 + A7*HY7 + A8*HY8 + A9*HY9;
    *BZ += A1*HZ1 + A2*HZ2 + A3*HZ3 + A4*HZ4 + A5*HZ5 + A6*HZ6 + A7*HZ7 + A8*HZ8 + A9*HZ9;

    return;

}






void    T01S_DEFORMED( int IOPT, double PS, double X, double Y, double Z,
        double *BX1, double *BY1, double *BZ1, double *BX2, double *BY2, double *BZ2, LgmTsyg2001_Info *t ) {

    /*
     *  RH0,RH1,RH2, AND IEPS CONTROL THE TILT-RELATED DEFORMATION OF THE TAIL FIELD
     */
    //static int      IEPS = 3; replaced explicitly below
    static double   RH2  = -5.2;
    double      SPS, CPS, R2, R, ZR, RH, DRHDR, DRHDZ, RRH, F, DFDR, DFDRH;
    double      SPSAS, CPSAS, XAS, ZAS, FACPS, PSASX, PSASY, PSASZ;
    double      DXASDX, DXASDY, DXASDZ, DZASDX, DZASDY, DZASDZ, FAC1, FAC2, FAC3;
    double      ooR, ZR2, RRH2, RRH3, F2, F4, ooCPSAS;
    double      BXAS1, BYAS1, BZAS1, BXAS2, BYAS2, BZAS2;


    /*
     *
     *   IOPT - TAIL FIELD MODE FLAG:   IOPT=0 - THE TWO TAIL MODES ARE ADDED UP
     *                                  IOPT=1 - MODE 1 ONLY
     *                                  IOPT=2 - MODE 2 ONLY
     *
     *   CALCULATES GSM COMPONENTS OF TWO UNIT-AMPLITUDE TAIL FIELD MODES,
     *    TAKING INTO ACCOUNT BOTH EFFECTS OF T01S_DIPOLE TILT:
     *    WARPING IN Y-Z (DONE BY THE S/R T01S_WARPED) AND BENDING IN X-Z (DONE BY THIS SUBROUTINE)
     */



    SPS = t->sin_psi; CPS = t->cos_psi;
    R2    =  X*X + Y*Y + Z*Z;
    R     =  sqrt( R2 ); ooR = 1.0/R;
    ZR    =  Z*ooR; ZR2 = ZR*ZR; 
    RH    =  t->CB_RH0.RH0 + RH2*ZR2;
    DRHDZ =  ZR*ooR*2.0*RH2;
    DRHDR = -DRHDZ*ZR;

    RRH   =  R/RH; RRH2 = RRH*RRH;  RRH3 = RRH2*RRH;

    F     =  1.0/pow(1.0+RRH3, 1.0/3.0); F2 = F*F; F4 = F2*F2;
    DFDR  = -RRH2*F4/RH;
    DFDRH = -RRH*DFDR;

    SPSAS = SPS*F; 
    CPSAS = sqrt(1.0-SPSAS*SPSAS);

    XAS   = X*CPSAS - Z*SPSAS;
    ZAS   = X*SPSAS + Z*CPSAS;

    ooCPSAS = 1.0/CPSAS;
    FACPS = SPS*ooCPSAS*(DFDR+DFDRH*DRHDR)*ooR;
    PSASX = FACPS*X;
    PSASY = FACPS*Y;
    PSASZ = FACPS*Z+SPS*ooCPSAS*DFDRH*DRHDZ;

    DXASDX = CPSAS-ZAS*PSASX;
    DXASDY = -ZAS*PSASY;
    DXASDZ = -SPSAS-ZAS*PSASZ;
    DZASDX = SPSAS+XAS*PSASX;
    DZASDY = XAS*PSASY;
    DZASDZ = CPSAS+XAS*PSASZ;
    FAC1 = DXASDZ*DZASDY-DXASDY*DZASDZ;
    FAC2 = DXASDX*DZASDZ-DXASDZ*DZASDX;
    FAC3 = DZASDX*DXASDY-DXASDX*DZASDY;
//printf("FAC1, FAC2, FAC3 = %g %g %g\n", FAC1, FAC2, FAC3);

    /*
     *      DEFORM:
     */
    T01S_WARPED( IOPT, PS, XAS, Y, ZAS, &BXAS1, &BYAS1, &BZAS1, &BXAS2, &BYAS2, &BZAS2, t );
//printf("IOPT, PS, XAS, Y, ZAS = %d %g %g %g %g\n", IOPT, PS, XAS, Y, ZAS);

    *BX1 = BXAS1*DZASDZ - BZAS1*DXASDZ + BYAS1*FAC1;
    *BY1 = BYAS1*FAC2;
    *BZ1 = BZAS1*DXASDX - BXAS1*DZASDX + BYAS1*FAC3;
//printf("BX1, BY1, BZ1 = %g %g %g\n", *BX1, *BY1, *BZ1);

    *BX2 = BXAS2*DZASDZ - BZAS2*DXASDZ + BYAS2*FAC1;
    *BY2 = BYAS2*FAC2;
    *BZ2 = BZAS2*DXASDX - BXAS2*DZASDX + BYAS2*FAC3;
//printf("BX2, BY2, BZ2 = %g %g %g\n", *BX2, *BY2, *BZ2);

    return;

}




void T01S_WARPED( int IOPT, double PS, double X, double Y, double Z,
    double *BX1, double *BY1, double *BZ1, double *BX2, double *BY2, double *BZ2, LgmTsyg2001_Info *t ) {

    /*
     *   CALCULATES GSM COMPONENTS OF THE T01S_WARPED FIELD FOR TWO TAIL UNIT MODES.
     *   THE WARPING DEFORMATION IS IMPOSED ON THE T01S_UNWARPED FIELD, COMPUTED
     *   BY THE S/R "T01S_UNWARPED".  THE WARPING PARAMETER G WAS OBTAINED BY LEAST
     *   SQUARES FITTING TO THE ENTIRE DATASET.
     *
     *   IOPT - TAIL FIELD MODE FLAG:   IOPT=0 - THE TWO TAIL MODES ARE ADDED UP
     *                                  IOPT=1 - MODE 1 ONLY
     *                                  IOPT=2 - MODE 2 ONLY
     */

    double  DGDX, XL, DXLDX, SPS, RHO2, RHO, PHI, CPHI, SPHI, XL2, XL4, XL3, RHO4, RR4L4;
    double  F, G, DFDPHI, RR4L42, DFDRHO, DFDX, CF, SF, YAS, ZAS, BX_AS1, BY_AS1, BZ_AS1;
    double  BX_AS2, BY_AS2, BZ_AS2, BRHO_AS, BPHI_AS, BRHO_S, BPHI_S;


    G = t->CB_G.G;

    DGDX  = 0.0;
    XL    = 20.0;
    DXLDX = 0.0;

    SPS  = sin( PS );
    RHO2 = Y*Y+Z*Z;
    RHO  = sqrt( RHO2 );

    if ( (Y == 0.0) && (Z == 0.0) ) {
        PHI  = 0.0;
        CPHI = 1.0;
        SPHI = 0.0;
    } else {
        PHI  = atan2( Z, Y );
        CPHI = Y/RHO;
        SPHI = Z/RHO;
    }

    XL2 = XL*XL; XL4 = XL2*XL2; XL3 = XL2*XL; RHO4 = RHO2*RHO2;
    RR4L4 = RHO/(RHO4+XL4);


    F      = PHI+G*RHO2*RR4L4*CPHI*SPS;
    DFDPHI = 1.0-G*RHO2*RR4L4*SPHI*SPS;
    RR4L42 = RR4L4*RR4L4;
    DFDRHO = G*RR4L42*(3.0*XL4-RHO4)*CPHI*SPS;
    DFDX   = RR4L4*CPHI*SPS*(DGDX*RHO2-G*RHO*RR4L4*4.0*XL3*DXLDX);

    CF  = cos(F);
    SF  = sin(F);
    YAS = RHO*CF;
    ZAS = RHO*SF;

    T01S_UNWARPED( IOPT, X, YAS, ZAS, &BX_AS1, &BY_AS1, &BZ_AS1, &BX_AS2, &BY_AS2, &BZ_AS2, t );
//printf("IOPT, X, YAS, ZAS = %d %g %g %g\n", IOPT, X, YAS, ZAS);
//printf("BX_AS1, BY_AS1, BZ_AS1, BX_AS2, BY_AS2, BZ_AS2 = %g %g %g %g %g %g\n", BX_AS1, BY_AS1, BZ_AS1, BX_AS2, BY_AS2, BZ_AS2);

    BRHO_AS  =   BY_AS1*CF+BZ_AS1*SF;       //   DEFORM THE 1ST MODE
    BPHI_AS  =  -BY_AS1*SF+BZ_AS1*CF;

    BRHO_S  =  BRHO_AS*DFDPHI;
    BPHI_S  =  BPHI_AS-RHO*(BX_AS1*DFDX+BRHO_AS*DFDRHO);
    *BX1     =  BX_AS1*DFDPHI;

    *BY1     =  BRHO_S*CPHI-BPHI_S*SPHI;
    *BZ1     =  BRHO_S*SPHI+BPHI_S*CPHI;        //   DONE

    BRHO_AS =   BY_AS2*CF+BZ_AS2*SF;        //   DEFORM THE 2ND MODE
    BPHI_AS =  -BY_AS2*SF+BZ_AS2*CF;

    BRHO_S  =  BRHO_AS*DFDPHI;
    BPHI_S  =  BPHI_AS-RHO*(BX_AS2*DFDX+BRHO_AS*DFDRHO);
    *BX2     =  BX_AS2*DFDPHI;

    *BY2     = BRHO_S*CPHI-BPHI_S*SPHI;
    *BZ2     = BRHO_S*SPHI+BPHI_S*CPHI;     //   DONE

    return;

}

void    T01S_UNWARPED( int IOPT, double X, double Y, double Z,
        double *BX1, double *BY1, double *BZ1, double *BX2, double *BY2, double *BZ2, LgmTsyg2001_Info *t ) {

    /*   IOPT - TAIL FIELD MODE FLAG:   IOPT=0 - THE TWO TAIL MODES ARE ADDED UP
     *                                  IOPT=1 - MODE 1 ONLY
     *                                  IOPT=2 - MODE 2 ONLY
     *
     *    CALCULATES GSM COMPONENTS OF THE SHIELDED FIELD OF TWO TAIL MODES WITH UNIT
     *    AMPLITUDES,  WITHOUT ANY WARPING OR BENDING.  NONLINEAR PARAMETERS OF THE MODES
     *    ARE FORWARDED HERE VIA A COMMON BLOCK /TAIL/.
     */

    // TAIL SHIELDING FIELD PARAMETERS FOR THE MODES #1 & #2
    double    DELTADX1=1.0, ALPHA1=1.1,  XSHIFT1=6.0;
    static double A1[] = { -9e99, -25.45869857,57.35899080,317.5501869,-2.626756717,
                     -93.38053698,-199.6467926,-858.8129729,34.09192395,845.4214929,
                     -29.07463068,47.10678547,-128.9797943,-781.7512093,6.165038619,
                     167.8905046,492.0680410,1654.724031,-46.77337920,-1635.922669,
                     40.86186772,-.1349775602,-.9661991179E-01,-.1662302354,
                     .002810467517,.2487355077,.1025565237,-14.41750229,-.8185333989,
                     11.07693629,.7569503173,-9.655264745,112.2446542,777.5948964,
                     -5.745008536,-83.03921993,-490.2278695,-1155.004209,39.08023320,
                     1172.780574,-39.44349797,-14.07211198,-40.41201127,-313.2277343,
                     2.203920979,8.232835341,197.7065115,391.2733948,-18.57424451,
                     -437.2779053,23.04976898,11.75673963,13.60497313,4.691927060,
                     18.20923547,27.59044809,6.677425469,1.398283308,2.839005878,
                     31.24817706,24.53577264};

    double    DELTADX2=0.0, ALPHA2=0.25, XSHIFT2=4.0;
    static double A2[] = { -9e99, -287187.1962,4970.499233,410490.1952,-1347.839052,
                     -386370.3240,3317.983750,-143462.3895,5706.513767,171176.2904,
                     250.8882750,-506570.8891,5733.592632,397975.5842,9771.762168,
                     -941834.2436,7990.975260,54313.10318,447.5388060,528046.3449,
                     12751.04453,-21920.98301,-21.05075617,31971.07875,3012.641612,
                     -301822.9103,-3601.107387,1797.577552,-6.315855803,142578.8406,
                     13161.93640,804184.8410,-14168.99698,-851926.6360,-1890.885671,
                     972475.6869,-8571.862853,26432.49197,-2554.752298,-482308.3431,
                     -4391.473324,105155.9160,-1134.622050,-74353.53091,-5382.670711,
                     695055.0788,-916.3365144,-12111.06667,67.20923358,-367200.9285,
                     -21414.14421,14.75567902,20.75638190,59.78601609,16.86431444,
                     32.58482365,23.69472951,17.24977936,13.64902647,68.40989058,
                     11.67828167};

    double  XM1=-12.0, XM2=-12.0;
    double  XSC1, YSC1, ZSC1, D0SC1, FX1, FY1, FZ1, HX1, HY1, HZ1;
    double  XSC2, YSC2, ZSC2, D0SC2, FX2, FY2, FZ2, HX2, HY2, HZ2;

    if (IOPT != 2) {

        XSC1  = (X-XSHIFT1-t->CB_TAIL.DXSHIFT1)*ALPHA1-XM1*(ALPHA1-1.0);
        YSC1  = Y*ALPHA1;
        ZSC1  = Z*ALPHA1;
        D0SC1 = t->CB_TAIL.D*ALPHA1;       // HERE WE USE A SINGLE VALUE D0 OF THE THICKNESS FOR BOTH MODES

        T01S_TAILDISK( D0SC1, DELTADX1, t->CB_TAIL.DELTADY, XSC1, YSC1, ZSC1, &FX1, &FY1, &FZ1 );
        T01S_SHLCAR5X5( A1, X, Y, Z, t->CB_TAIL.DXSHIFT1, &HX1, &HY1, &HZ1 );
//printf("X, Y, Z, t->CB_TAIL.DXSHIFT1 = %g %g %g %g\n", X, Y, Z, t->CB_TAIL.DXSHIFT1);
//printf("FX1, FY1, FZ1 = %g %g %g\n", FX1, FY1, FZ1);
//printf("HX1, HY1, HZ1 = %g %g %g\n", HX1, HY1, HZ1);

        *BX1 = FX1 + HX1;
        *BY1 = FY1 + HY1;
        *BZ1 = FZ1 + HZ1;

        if (IOPT == 1) {
            *BX2 = 0.0;
            *BY2 = 0.0;
            *BZ2 = 0.0;
            return;
        }

    }

    XSC2  = (X-XSHIFT2-t->CB_TAIL.DXSHIFT2)*ALPHA2-XM2*(ALPHA2-1.0);
    YSC2  = Y*ALPHA2;
    ZSC2  = Z*ALPHA2;
    D0SC2 = t->CB_TAIL.D*ALPHA2;           // HERE WE USE A SINGLE VALUE D0 OF THE THICKNESS FOR BOTH MODES


    T01S_TAILDISK( D0SC2, DELTADX2, t->CB_TAIL.DELTADY, XSC2, YSC2, ZSC2, &FX2, &FY2, &FZ2);
    T01S_SHLCAR5X5( A2, X, Y, Z, t->CB_TAIL.DXSHIFT2, &HX2, &HY2, &HZ2);
//printf("X, Y, Z, t->CB_TAIL.DXSHIFT2 = %g %g %g %g\n", X, Y, Z, t->CB_TAIL.DXSHIFT2);
//printf("FX2, FY2, FZ2 = %g %g %g\n", FX2, FY2, FZ2);
//printf("HX2, HY2, HZ2 = %g %g %g\n", HX2, HY2, HZ2);

    *BX2 = FX2 + HX2;
    *BY2 = FY2 + HY2;
    *BZ2 = FZ2 + HZ2;

    if (IOPT == 2) {
    *BX1 = 0.0;
    *BY1 = 0.0;
    *BZ1 = 0.0;
    }

    return;

}

void    T01S_TAILDISK( double D0, double DELTADX, double DELTADY, double X, double Y, double Z, double *BX, double *BY, double *BZ) {

    /*
     *  THIS SUBROUTINE COMPUTES THE COMPONENTS OF THE TAIL CURRENT
     *  FIELD, SIMILAR TO THAT DESCRIBED BY TSYGANENKO AND PEREDO (1994).
     *  THE DIFFERENCE IS THAT NOW WE USE SPACEWARPING, AS DESCRIBED IN
     *  OUR PT01S_APER ON MODELING BIRKELAND CURRENTS (TSYGANENKO AND STERN,
     *  1996) INSTEAD OF SHEARING IT IN THE SPIRIT OF T89 TAIL MODEL.
     */

    static double   F[] = { -9e99, -71.09346626, -1014.308601, -1272.939359, -3224.935936, -44546.86232 };
    static double   B[] = { -9e99, 10.90101242, 12.68393898, 13.51791954, 14.86775017, 15.12306404 };
    static double   C[] = { -9e99, .7954069972, .6716601849, 1.174866319, 2.565249920, 10.01986790 };
    double      Y2, RHO, DRHODX, DRHODY, ooSeven, DEX, D, DDDY, DDDX, DZETA, ooDZETA, DDZETADX, DDZETADY, DDZETADZ;
    double      DBX, CI, RHOpBI, BI, RHOpBI2, RHOmBI, RHOmBI2, DZETApCI, DZETApCI2, S1, S2, ooS1, ooS2;
    double      DS1DRHO, DS2DRHO, DS1DDZ, DS2DDZ, DS1DX, DS1DY, DS1DZ, DS2DX, DS2DY, DS2DZ, S1TS2, S1PS2, S1PS2SQ;
    double      BI2, FAC1, AS, DASDS1, DASDS2, DASDX, DASDY, DASDZ, DBY, DBZ;
    int         I;


    Y2 = Y*Y;
    RHO    = sqrt(X*X+Y2);
    DRHODX = X/RHO;
    DRHODY = Y/RHO;

    ooSeven = 1.0/7.0;

    DEX  = exp(X*ooSeven);
    D    = D0+DELTADY*0.0025*Y2  + DELTADX*DEX;     // THE LAST TERM (INTRODUCED 10/11/2000) MAKES THE SHEET
    DDDY = DELTADY*Y*0.005;             // THICKEN SUNWARD, TO AVOID PROBLEMS IN THE SUBSOLAR REGION
    DDDX = DELTADX*ooSeven*DEX;

    DZETA    = sqrt(Z*Z+D*D);               // THIS IS THE SAME SIMPLE WAY TO SPREAD OUT THE SHEET, AS THAT USED IN T89
    ooDZETA  = 1.0/DZETA;
    DDZETADX = D*DDDX*ooDZETA;
    DDZETADY = D*DDDY*ooDZETA;
    DDZETADZ = Z*ooDZETA;

    DBX = 0.0;
    DBY = 0.0;
    DBZ = 0.0;
    for ( I=1; I<=5; I++ ){

        BI = B[I];
        CI = C[I];

        RHOpBI = RHO+BI; RHOpBI2 = RHOpBI*RHOpBI;
        RHOmBI = RHO-BI; RHOmBI2 = RHOmBI*RHOmBI;
        DZETApCI = DZETA+CI; DZETApCI2 = DZETApCI*DZETApCI;
        S1 = sqrt(RHOpBI2+DZETApCI2);
        S2 = sqrt(RHOmBI2+DZETApCI2);

        ooS1 = 1.0/S1; ooS2 = 1.0/S2;
        DS1DRHO = RHOpBI*ooS1;
        DS2DRHO = RHOmBI*ooS2;
        DS1DDZ  = DZETApCI*ooS1;
        DS2DDZ  = DZETApCI*ooS2;

        DS1DX = DS1DRHO*DRHODX  + DS1DDZ*DDZETADX;
        DS1DY = DS1DRHO*DRHODY  + DS1DDZ*DDZETADY;
        DS1DZ =                   DS1DDZ*DDZETADZ;

        DS2DX = DS2DRHO*DRHODX  + DS2DDZ*DDZETADX;
        DS2DY = DS2DRHO*DRHODY  + DS2DDZ*DDZETADY;
        DS2DZ =                   DS2DDZ*DDZETADZ;

        S1TS2   = S1*S2;
        S1PS2   = S1+S2;
        S1PS2SQ = S1PS2*S1PS2;

        BI2 = BI*BI;
        FAC1   = sqrt(S1PS2SQ-4.0*BI2);
        AS     = FAC1/(S1TS2*S1PS2SQ);
        DASDS1 = (1.0/(FAC1*S2)-AS/S1PS2*(S2*S2+S1*(3.0*S1+4.0*S2))) /(S1*S1PS2);
        DASDS2 = (1.0/(FAC1*S1)-AS/S1PS2*(S1*S1+S2*(3.0*S2+4.0*S1))) /(S2*S1PS2);

        DASDX = DASDS1*DS1DX + DASDS2*DS2DX;
        DASDY = DASDS1*DS1DY + DASDS2*DS2DY;
        DASDZ = DASDS1*DS1DZ + DASDS2*DS2DZ;

        DBX = DBX - F[I]*X*DASDZ;
        DBY = DBY - F[I]*Y*DASDZ;
        DBZ = DBZ + F[I]*(2.0*AS+X*DASDX+Y*DASDY);

    }

    *BX = DBX;
    *BY = DBY;
    *BZ = DBZ;

    return;


}

void    T01S_SHLCAR5X5( double *A, double X, double Y, double Z, double DSHIFT, double *HX, double *HY, double *HZ ) {

    /*
     *
     *  THIS CODE RETURNS THE SHIELDING FIELD REPRESENTED BY  5x5=25 "CARTESIAN"
     *     HARMONICS
     *
     *
     *  - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     *   The NLIN coefficients are the amplitudes of the "cartesian"
     *     harmonics (A(1)-A(NLIN).
     *   The NNP nonlinear parameters (A(NLIN+1)-A(NTOT) are the scales Pi and Ri
     *    entering the arguments of exponents, sines, and cosines in each of the
     *    NLIN "Cartesian" harmonics
     *  - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     */

    int     I, K, L;
    double  DHX, DHY, DHZ, RP, RP2, YRP, CYPI, SYPI, RR, RR2;
    double  SZRK, CZRK, SQPR, EPR, DBX, DBY, DBZ, COEF;

    DHX = 0.0;
    DHY = 0.0;
    DHZ = 0.0;

    L = 0;

    for ( I=1; I<=5; I++ ) {

        RP   = 1.0/A[50+I];   // we could get rid of these divides by hand computing the inverse cooeffs...
        RP2 = RP*RP;
        YRP = Y*RP;
        CYPI = cos(YRP);
        SYPI = sin(YRP);

        for ( K=1; K<=5; K++ ) {

            RR   = 1.0/A[55+K];   // we could get rid of these divides by hand computing the inverse cooeffs...
            RR2  = RR*RR;
            SZRK = sin(Z*RR);
            CZRK = cos(Z*RR);
            SQPR = sqrt(RP2+RR2);
            EPR  = exp(X*SQPR);

            DBX = -SQPR*EPR*CYPI*SZRK;
            DBY =  RP*EPR*SYPI*SZRK;
            DBZ = -RR*EPR*CYPI*CZRK;

            L += 2;
            COEF = A[L-1]+A[L]*DSHIFT;

            DHX += COEF*DBX;
            DHY += COEF*DBY;
            DHZ += COEF*DBZ;

        }
    }


    *HX = DHX;
    *HY = DHY;
    *HZ = DHZ;


    return;

}








void    T01S_BIRK_TOT( int IOPB, double PS, double X, double Y, double Z,
        double *BX11, double *BY11, double *BZ11, double *BX12, double *BY12, double *BZ12,
        double *BX21, double *BY21, double *BZ21, double *BX22, double *BY22, double *BZ22, LgmTsyg2001_Info *t ){

    /*
     *
     *
     *        IOPB -  BIRKELAND FIELD MODE FLAG:
     *           IOPB=0 - ALL COMPONENTS
     *           IOPB=1 - REGION 1, MODES 1 & 2
     *           IOPB=2 - REGION 2, MODES 1 & 2
     */

    static double SH11[] = { -9e99, 46488.84663,-15541.95244,-23210.09824,-32625.03856,
                         -109894.4551,-71415.32808,58168.94612,55564.87578,-22890.60626,
                         -6056.763968,5091.368100,239.7001538,-13899.49253,4648.016991,
                         6971.310672,9699.351891,32633.34599,21028.48811,-17395.96190,
                         -16461.11037,7447.621471,2528.844345,-1934.094784,-588.3108359,
                         -32588.88216,10894.11453,16238.25044,22925.60557,77251.11274,
                         50375.97787,-40763.78048,-39088.60660,15546.53559,3559.617561,
                         -3187.730438,309.1487975,88.22153914,-243.0721938,-63.63543051,
                         191.1109142,69.94451996,-187.9539415,-49.89923833,104.0902848,
                         -120.2459738,253.5572433,89.25456949,-205.6516252,-44.93654156,
                         124.7026309,32.53005523,-98.85321751,-36.51904756,98.88241690,
                         24.88493459,-55.04058524,61.14493565,-128.4224895,-45.35023460,
                         105.0548704,-43.66748755,119.3284161,31.38442798,-92.87946767,
                         -33.52716686,89.98992001,25.87341323,-48.86305045,59.69362881,
                         -126.5353789,-44.39474251,101.5196856,59.41537992,41.18892281,
                         80.86101200,3.066809418,7.893523804,30.56212082,10.36861082,
                         8.222335945,19.97575641,2.050148531,4.992657093,2.300564232,
                         .2256245602,-.05841594319};

    static double SH12[] = { -9e99, 210260.4816,-1443587.401,-1468919.281,281939.2993,
                         -1131124.839,729331.7943,2573541.307,304616.7457,468887.5847,
                         181554.7517,-1300722.650,-257012.8601,645888.8041,-2048126.412,
                         -2529093.041,571093.7972,-2115508.353,1122035.951,4489168.802,
                         75234.22743,823905.6909,147926.6121,-2276322.876,-155528.5992,
                         -858076.2979,3474422.388,3986279.931,-834613.9747,3250625.781,
                         -1818680.377,-7040468.986,-414359.6073,-1295117.666,-346320.6487,
                         3565527.409,430091.9496,-.1565573462,7.377619826,.4115646037,
                         -6.146078880,3.808028815,-.5232034932,1.454841807,-12.32274869,
                         -4.466974237,-2.941184626,-.6172620658,12.64613490,1.494922012,
                         -21.35489898,-1.652256960,16.81799898,-1.404079922,-24.09369677,
                         -10.99900839,45.94237820,2.248579894,31.91234041,7.575026816,
                         -45.80833339,-1.507664976,14.60016998,1.348516288,-11.05980247,
                         -5.402866968,31.69094514,12.28261196,-37.55354174,4.155626879,
                         -33.70159657,-8.437907434,36.22672602,145.0262164,70.73187036,
                         85.51110098,21.47490989,24.34554406,31.34405345,4.655207476,
                         5.747889264,7.802304187,1.844169801,4.867254550,2.941393119,
                         .1379899178,.06607020029};

    static double SH21[] = { -9e99, 162294.6224,503885.1125,-27057.67122,-531450.1339,
                         84747.05678,-237142.1712,84133.61490,259530.0402,69196.05160,
                         -189093.5264,-19278.55134,195724.5034,-263082.6367,-818899.6923,
                         43061.10073,863506.6932,-139707.9428,389984.8850,-135167.5555,
                         -426286.9206,-109504.0387,295258.3531,30415.07087,-305502.9405,
                         100785.3400,315010.9567,-15999.50673,-332052.2548,54964.34639,
                         -152808.3750,51024.67566,166720.0603,40389.67945,-106257.7272,
                         -11126.14442,109876.2047,2.978695024,558.6019011,2.685592939,
                         -338.0004730,-81.99724090,-444.1102659,89.44617716,212.0849592,
                         -32.58562625,-982.7336105,-35.10860935,567.8931751,-1.917212423,
                         -260.2023543,-1.023821735,157.5533477,23.00200055,232.0603673,
                         -36.79100036,-111.9110936,18.05429984,447.0481000,15.10187415,
                         -258.7297813,-1.032340149,-298.6402478,-1.676201415,180.5856487,
                         64.52313024,209.0160857,-53.85574010,-98.52164290,14.35891214,
                         536.7666279,20.09318806,-309.7349530,58.54144539,67.45226850,
                         97.92374406,4.752449760,10.46824379,32.91856110,12.05124381,
                         9.962933904,15.91258637,1.804233877,6.578149088,2.515223491,
                         .1930034238,-.02261109942};

    static double SH22[] = { -9e99, -131287.8986,-631927.6885,-318797.4173,616785.8782,
                         -50027.36189,863099.9833,47680.20240,-1053367.944,-501120.3811,
                         -174400.9476,222328.6873,333551.7374,-389338.7841,-1995527.467,
                         -982971.3024,1960434.268,297239.7137,2676525.168,-147113.4775,
                         -3358059.979,-2106979.191,-462827.1322,1017607.960,1039018.475,
                         520266.9296,2627427.473,1301981.763,-2577171.706,-238071.9956,
                         -3539781.111,94628.16420,4411304.724,2598205.733,637504.9351,
                         -1234794.298,-1372562.403,-2.646186796,-31.10055575,2.295799273,
                         19.20203279,30.01931202,-302.1028550,-14.78310655,162.1561899,
                         .4943938056,176.8089129,-.2444921680,-100.6148929,9.172262228,
                         137.4303440,-8.451613443,-84.20684224,-167.3354083,1321.830393,
                         76.89928813,-705.7586223,18.28186732,-770.1665162,-9.084224422,
                         436.3368157,-6.374255638,-107.2730177,6.080451222,65.53843753,
                         143.2872994,-1028.009017,-64.22739330,547.8536586,-20.58928632,
                         597.3893669,10.17964133,-337.7800252,159.3532209,76.34445954,
                         84.74398828,12.76722651,27.63870691,32.69873634,5.145153451,
                         6.310949163,6.996159733,1.971629939,4.436299219,2.904964304,
                         .1486276863,.06859991529};


//    double        XKAPPA;
    double      X_SC, FX11, FY11, FZ11, HX11, HY11, HZ11, FX12, FY12, FZ12, HX12, HY12, HZ12;
    double      FX21, FY21, FZ21, HX21, HY21, HZ21, FX22, FY22, FZ22, HX22, HY22, HZ22;
    int         PSChanged, XChanged, YChanged, ZChanged;

    PSChanged = ( fabs(PS-t->OLD_PS) > 1e-10 ) ? TRUE : FALSE;
    XChanged = ( fabs(X-t->OLD_X) > 1e-10 ) ? TRUE : FALSE;
    YChanged = ( fabs(Y-t->OLD_Y) > 1e-10 ) ? TRUE : FALSE;
    ZChanged = ( fabs(Z-t->OLD_Z) > 1e-10 ) ? TRUE : FALSE;

    t->CB_DPHI_B_RHO0.XKAPPA = t->CB_BIRKPAR.XKAPPA1;     // FORWARDED IN T01S_BIRK_1N2
    X_SC   = t->CB_BIRKPAR.XKAPPA1-1.1;    // FORWARDED IN T01S_BIRK_SHL

    if ( (IOPB == 0) || (IOPB == 1) ) {

        T01S_BIRK_1N2( 1, 1, PS, X, Y, Z, &FX11, &FY11, &FZ11, t );     // REGION 1,  MODE 1
        T01S_BIRK_SHL( 0, PSChanged, XChanged, YChanged, ZChanged, SH11, PS, X_SC, X, Y, Z, &HX11, &HY11, &HZ11, t );
        *BX11 = FX11 + HX11;
        *BY11 = FY11 + HY11;
        *BZ11 = FZ11 + HZ11;

        T01S_BIRK_1N2( 1, 2, PS, X, Y, Z, &FX12, &FY12, &FZ12, t );     // REGION 1,  MODE 2
        T01S_BIRK_SHL( 1, PSChanged, XChanged, YChanged, ZChanged, SH12, PS, X_SC, X, Y, Z, &HX12, &HY12, &HZ12, t );
        *BX12 = FX12 + HX12;
        *BY12 = FY12 + HY12;
        *BZ12 = FZ12 + HZ12;

    }

    t->CB_DPHI_B_RHO0.XKAPPA = t->CB_BIRKPAR.XKAPPA2;     // FORWARDED IN T01S_BIRK_1N2
    X_SC   = t->CB_BIRKPAR.XKAPPA2-1.0;    // FORWARDED IN T01S_BIRK_SHL

    if ( (IOPB == 0)  || (IOPB == 2) ) {

        T01S_BIRK_1N2( 2, 1, PS, X, Y, Z, &FX21, &FY21, &FZ21, t );     // REGION 2,  MODE 1
        T01S_BIRK_SHL( 2, PSChanged, XChanged, YChanged, ZChanged, SH21, PS, X_SC, X, Y, Z, &HX21, &HY21, &HZ21, t );
        *BX21 = FX21 + HX21;
        *BY21 = FY21 + HY21;
        *BZ21 = FZ21 + HZ21;

        T01S_BIRK_1N2( 2, 2, PS, X, Y, Z, &FX22, &FY22, &FZ22, t );     // REGION 2,  MODE 2
        T01S_BIRK_SHL( 3, PSChanged, XChanged, YChanged, ZChanged, SH22, PS, X_SC, X, Y, Z, &HX22, &HY22, &HZ22, t );
        *BX22 = FX22 + HX22;
        *BY22 = FY22 + HY22;
        *BZ22 = FZ22 + HZ22;

    }

    t->OLD_PS = PS;
    t->OLD_X = X;
    t->OLD_Y = Y;
    t->OLD_Z = Z;

    return;

}



void T01S_BIRK_1N2( int NUMB, int MODE, double PS, double X, double Y, double Z, double *BX, double *BY, double *BZ, LgmTsyg2001_Info *t ) {

    /*
     *
     *
     *
     *    CALCULATES COMPONENTS  OF REGION 1/2 FIELD IN SPHERICAL COORDS.  DERIVED FROM THE S/R DIPDEF2C (WHICH
     *      DOES THE SAME JOB, BUT INPUT/OUTPUT THERE WAS IN SPHERICAL COORDS, WHILE HERE WE USE CARTESIAN ONES)
     *
     *     INPUT:  NUMB=1 (2) FOR REGION 1 (2) CURRENTS
     *             MODE=1 YIELDS SIMPLE SINUSOIDAL MLT VARIATION, WITH MAXIMUM CURRENT AT DAWN/DUSK MERIDIAN
     *       WHILE MODE=2 YIELDS THE SECOND HARMONIC.
     *
     *  (1) DPHI:   HALF-DIFFERENCE (IN RADIANS) BETWEEN DAY AND NIGHT LATITUDE OF FAC OVAL AT IONOSPHERIC ALTITUDE;
     *              TYPICAL VALUE: 0.06
     *  (2) B:      AN ASYMMETRY FACTOR AT HIGH-ALTITUDES;  FOR B=0, THE ONLY ASYMMETRY IS THAT FROM DPHI
     *              TYPICAL VALUES: 0.35-0.70
     *  (3) RHO_0:  A FIXED PARAMETER, DEFINING THE DISTANCE RHO, AT WHICH THE LATITUDE SHIFT GRADUALLY SATURATES AND
     *              STOPS INCREASING
     *              ITS VALUE WAS ASSUMED FIXED, EQUAL TO 7.0.
     *  (4) XKAPPA: AN OVERALL SCALING FACTOR, WHICH CAN BE USED FOR CHANGING THE SIZE OF THE F.A.C. OVAL
     *
     */

    static double A11[] = { -9e99, .1618068350,-.1797957553,2.999642482,-.9322708978,
                             -.6811059760,.2099057262,-8.358815746,-14.86033550,.3838362986,
                             -16.30945494,4.537022847,2.685836007,27.97833029,6.330871059,
                             1.876532361,18.95619213,.9651528100,.4217195118,-.08957770020,
                             -1.823555887,.7457045438,-.5785916524,-1.010200918,.01112389357,
                             .09572927448,-.3599292276,8.713700514,.9763932955,3.834602998,
                             2.492118385,.7113544659};
    static double A12[] = { -9e99, .7058026940,-.2845938535,5.715471266,-2.472820880,
                             -.7738802408,.3478293930,-11.37653694,-38.64768867,.6932927651,
                             -212.4017288,4.944204937,3.071270411,33.05882281,7.387533799,
                             2.366769108,79.22572682,.6154290178,.5592050551,-.1796585105,
                             -1.654932210,.7309108776,-.4926292779,-1.130266095,-.009613974555,
                             .1484586169,-.2215347198,7.883592948,.02768251655,2.950280953,
                             1.212634762,.5567714182};
    static double A21[] = { -9e99, .1278764024,-.2320034273,1.805623266,-32.37241440,
                             -.9931490648,.3175085630,-2.492465814,-16.21600096,.2695393416,
                             -6.752691265,3.971794901,14.54477563,41.10158386,7.912889730,
                             1.258297372,9.583547721,1.014141963,.5104134759,-.1790430468,
                             -1.756358428,.7561986717,-.6775248254,-.04014016420,.01446794851,
                             .1200521731,-.2203584559,4.508963850,.8221623576,1.779933730,
                             1.102649543,.8867880020};
    static double A22[] = { -9e99, .4036015198,-.3302974212,2.827730930,-45.44405830,
                             -1.611103927,.4927112073,-.003258457559,-49.59014949,.3796217108,
                             -233.7884098,4.312666980,18.05051709,28.95320323,11.09948019,
                             .7471649558,67.10246193,.5667096597,.6468519751,-.1560665317,
                             -1.460805289,.7719653528,-.6658988668,.2515179349E-05,
                             .02426021891,.1195003324,-.2625739255,4.377172556,.2421190547,
                             2.503482679,1.071587299,.7247997430};

    double  BETA=0.9, RH=10.0;  // parameters of the tilt-dependent deformation of the untilted F.A.C. field
//    double    EPS=3.0;            // parameters of the tilt-dependent deformation of the untilted F.A.C. field
    double  RHO2, Xsc, Xsc2, Ysc, Zsc, Zsc2, RHO, RHOSQ;
    double  Rsc, Ysc2, PHI, SPHIC, CPHIC, BRACK, R1RH, R1RH2, R1RH3, PSIAS, PHIS, DPHISPHI;
    double  RHO2pRHOSQ, RHO2pRHOSQ2, DPHISRHO, DPHISDY, SPHICS, CPHICS, XS, ZS, BXS, BYAS;
    double  BZS, BRHOAS, BPHIAS, BRHO_S, BPHI_S, BY_S;

    t->CB_DPHI_B_RHO0.B     = 0.5;
    t->CB_DPHI_B_RHO0.RHO_0 = 7.0; RHO2 = 49.0;

    t->CB_MODENUM.M = MODE;
    if (NUMB == 1) {
        t->CB_DPHI_B_RHO0.DPHI   = 0.055;
        t->CB_DTHETA.DTHETA      = 0.06;
    } else if (NUMB == 2) {
        t->CB_DPHI_B_RHO0.DPHI   = 0.030;
        t->CB_DTHETA.DTHETA      = 0.09;
    }

    Xsc = X*t->CB_DPHI_B_RHO0.XKAPPA; Xsc2 = Xsc*Xsc;
    Ysc = Y*t->CB_DPHI_B_RHO0.XKAPPA; Ysc2 = Ysc*Ysc;
    Zsc = Z*t->CB_DPHI_B_RHO0.XKAPPA; Zsc2 = Zsc*Zsc;
    RHO = sqrt(Xsc2+Zsc2);
    RHOSQ = Xsc2+Zsc2;

    Rsc  = sqrt(Xsc2+Ysc2+Zsc2);        // SCALED

    if ( (Xsc == 0.0) && (Zsc == 0.0) ) {
        PHI = 0.0;
    } else {
        PHI = atan2(-Zsc, Xsc);         // FROM CARTESIAN TO CYLINDRICAL (RHO,PHI,Y)
    }

    SPHIC = sin(PHI);
    CPHIC = cos(PHI);               // "C" means "CYLINDRICAL", TO DISTINGUISH FROM SPHERICAL PHI

    BRACK = t->CB_DPHI_B_RHO0.DPHI+t->CB_DPHI_B_RHO0.B*RHO2/(RHO2+1.0)*(RHOSQ-1.0)/(RHO2+RHOSQ);
    R1RH  = (Rsc-1.0)/RH;
    R1RH2 = R1RH*R1RH; R1RH3 = R1RH2*R1RH;
    PSIAS = BETA*PS/pow( 1.0+R1RH3, 1.0/3.0 );

    PHIS = PHI-BRACK*SPHIC - PSIAS;
    DPHISPHI = 1.0-BRACK*CPHIC;
    RHO2pRHOSQ = RHO2+RHOSQ; RHO2pRHOSQ2 = RHO2pRHOSQ*RHO2pRHOSQ;
    DPHISRHO = -2.0*t->CB_DPHI_B_RHO0.B*RHO2*RHO/RHO2pRHOSQ2*SPHIC +BETA*PS*R1RH2*RHO/(RH*Rsc*pow( 1.0+R1RH3, 4.0/3.0 ));
    DPHISDY= BETA*PS*R1RH2*Ysc/(RH*Rsc*pow( 1.0+R1RH3, 4.0/3.0 ));

    SPHICS = sin(PHIS);
    CPHICS = cos(PHIS);

    XS =  RHO*CPHICS;
    ZS = -RHO*SPHICS;

    if (NUMB == 1) {
        if (MODE == 1) T01S_TWOCONES( A11, XS, Ysc, ZS, &BXS, &BYAS, &BZS, t );
        if (MODE == 2) T01S_TWOCONES( A12, XS, Ysc, ZS, &BXS, &BYAS, &BZS, t );
    } else {
        if (MODE == 1) T01S_TWOCONES( A21, XS, Ysc, ZS, &BXS, &BYAS, &BZS, t );
        if (MODE == 2) T01S_TWOCONES( A22, XS, Ysc, ZS, &BXS, &BYAS, &BZS, t );
    }

    BRHOAS = BXS*CPHICS-BZS*SPHICS;
    BPHIAS = -BXS*SPHICS-BZS*CPHICS;

    BRHO_S = BRHOAS*DPHISPHI                             *t->CB_DPHI_B_RHO0.XKAPPA;    // SCALING
    BPHI_S = (BPHIAS-RHO*(BYAS*DPHISDY+BRHOAS*DPHISRHO)) *t->CB_DPHI_B_RHO0.XKAPPA;
    BY_S   = BYAS*DPHISPHI                               *t->CB_DPHI_B_RHO0.XKAPPA;

    *BX = BRHO_S*CPHIC-BPHI_S*SPHIC;
    *BY = BY_S;
    *BZ = -BRHO_S*SPHIC-BPHI_S*CPHIC;

    return;

}


void    T01S_TWOCONES( double *A, double X, double Y, double Z, double *BX, double *BY, double *BZ, LgmTsyg2001_Info *t ){

    /*
     * ADDS FIELDS FROM TWO CONES (NORTHERN AND SOUTHERN), WITH A PROPER
     * SYMMETRY OF THE CURRENT AND FIELD, CORRESPONDING TO THE REGION 1
     * BIRKELAND CURRENTS. (NB #6, P.58).
     */
    double  BXN, BYN, BZN, BXS, BYS, BZS;

    T01S_ONE_CONE( A, X, Y, Z, &BXN, &BYN, &BZN, t );
    T01S_ONE_CONE( A, X, -Y, -Z, &BXS, &BYS, &BZS, t );
    *BX = BXN - BXS;
    *BY = BYN + BYS;
    *BZ = BZN + BZS;

    return;

}


void    T01S_ONE_CONE( double *A, double X, double Y, double Z, double *BX, double *BY, double *BZ, LgmTsyg2001_Info *t ) {

    /*
     *
     *
     *
     *  RETURNS FIELD COMPONENTS FOR A T01S_DEFORMED CONICAL CURRENT SYSTEM, FITTED TO A BIOSAVART FIELD
     *    BY SIM_14.FOR.  HERE ONLY THE NORTHERN CONE IS TAKEN INTO ACCOUNT.
     *
     */

    double  DR=1e-6, DT=1e-6;   // JUST FOR NUMERICAL DIFFERENTIATION
    double  THETA0, RHO2, RHO, R, THETA, PHI, RS, THETAS, PHIS, BTAST, BFAST;
    double  DRSDR, DRSDT, DTSDR, DTSDT, STSST, RSR, BR, BTHETA, BPHI;
    double  S, C, SF, CF, BE;


    THETA0 = A[31];

    RHO2  = X*X+Y*Y;
    RHO   = sqrt(RHO2);
    R     = sqrt(RHO2+Z*Z);
    THETA = atan2(RHO,Z);
    PHI   = atan2(Y,X);


    /*
     *   MAKE THE DEFORMATION OF COORDINATES:
     */
    RS     = T01S_R_S( A, R, THETA );
    THETAS = T01S_THETA_S( A, R, THETA );
    PHIS   = PHI;


    /*
     *   CALCULATE FIELD COMPONENTS AT THE NEW POSITION (ASTERISKED):
     */
    T01S_FIALCOS( RS, THETAS, PHIS, &BTAST, &BFAST, t->CB_MODENUM.M, THETA0, t->CB_DTHETA.DTHETA); // MODE #M


    /*
     *   NOW TRANSFORM B{R,T,F}_AST BY THE DEFORMATION TENSOR:
     *
     *      FIRST OF ALL, FIND THE DERIVATIVES:
     */
    DRSDR = (T01S_R_S(A,R+DR,THETA)-T01S_R_S(A,R-DR,THETA))/(2.0*DR);
    DRSDT = (T01S_R_S(A,R,THETA+DT)-T01S_R_S(A,R,THETA-DT))/(2.0*DT);
    DTSDR = (T01S_THETA_S(A,R+DR,THETA)-T01S_THETA_S(A,R-DR,THETA))/(2.0*DR);
    DTSDT = (T01S_THETA_S(A,R,THETA+DT)-T01S_THETA_S(A,R,THETA-DT))/(2.0*DT);

    STSST = sin(THETAS)/sin(THETA);
    RSR  =  RS/R;

    BR      =  -RSR/R*STSST*BTAST*DRSDT;
    BTHETA  =   RSR*STSST*BTAST*DRSDR;
    BPHI    =  RSR*BFAST*(DRSDR*DTSDT-DRSDT*DTSDR);

    S  = RHO/R;
    C  = Z/R;
    SF = Y/RHO;
    CF = X/RHO;

    BE = BR*S+BTHETA*C;

    *BX = A[1]*(BE*CF-BPHI*SF);
    *BY = A[1]*(BE*SF+BPHI*CF);
    *BZ = A[1]*(BR*C-BTHETA*S);

    return;

}





double T01S_R_S( double *A, double R, double THETA) {

    double  RS, R2, A11_2, A12_2, A13_2, A14_2, A15_2, A16_2, R2pA16_2, R2pA16_2_2;

    R2 = R*R;
    A11_2 = A[11]*A[11];
    A12_2 = A[12]*A[12];
    A13_2 = A[13]*A[13];
    A14_2 = A[14]*A[14];
    A15_2 = A[15]*A[15];
    A16_2 = A[16]*A[16];
    R2pA16_2 = R2+A16_2; R2pA16_2_2 = R2pA16_2*R2pA16_2;

    RS = R+A[2]/R +A[3]*R/sqrt(R2+A11_2)+A[4]*R/(R2+A12_2)
           +(A[5]+A[6]/R+A[7]*R/sqrt(R2+A13_2)+A[8]*R/(R2+A14_2))*cos(THETA)
           +(A[9]*R/sqrt(R2+A15_2)+A[10]*R/R2pA16_2_2)*cos(2.0*THETA);

      return( RS );

}
double T01S_THETA_S( double *A, double R, double THETA) {

    double  TS, R2, A27_2, A28_2, A29_2, A30_2;

    R2 = R*R;
    A27_2 = A[27]*A[27];
    A28_2 = A[28]*A[28];
    A29_2 = A[29]*A[29];
    A30_2 = A[30]*A[30];

    TS = THETA +(A[17]+A[18]/R+A[19]/R2
            +A[20]*R/sqrt(R2+A27_2))*sin(THETA)
            +(A[21]+A[22]*R/sqrt(R2+A28_2)
            +A[23]*R/(R2+A29_2))*sin(2.0*THETA)
            +(A[24]+A[25]/R+A[26]*R/(R2+A30_2))*sin(3.0*THETA);

      return( TS );
}


void    T01S_FIALCOS( double R, double THETA, double PHI, double *BTHETA, double *BPHI, int N, double THETA0, double DT) {

    /*
     *
     *
     *    CONICAL MODEL OF BIRKELAND CURRENT FIELD; BASED ON THE OLD S/R FIALCO (OF 1990-91)
     *
     *    BTN, AND BPN ARE THE ARRAYS OF BTHETA AND BPHI (BTN(i), BPN(i) CORRESPOND TO i-th MODE).
     *     ONLY FIRST  N  MODE AMPLITUDES ARE COMPUTED (N<=10).
     *      THETA0 IS THE ANGULAR HALF-WIDTH OF THE CONE, DT IS THE ANGULAR H.-W. OF THE CURRENT LAYER
     *
     *     NOTE:  BR=0  (BECAUSE ONLY RADIAL CURRENTS ARE PRESENT IN THIS MODEL)
     */

    int     M;
    double  SINTE, RO, COSTE, SINFI, COSFI, TG, CTG, TETANP, TETANM, TGP=0.0, TGM=0.0, TGM2=0.0, TGP2=0.0;
    double  COSM1, SINM1, TM, TGM2M, TGP2M, T, FC, FC1, TGM2M1, TG21, DTT, DTT0;
    double  BTN[10], BPN[10], CCOS[10], SSIN[10];

    SINTE = sin(THETA);
    RO    = R*SINTE;
    COSTE = cos(THETA);
    SINFI = sin(PHI);
    COSFI = cos(PHI);
    TG    = SINTE/(1.0+COSTE);      //  TAN(THETA/2)
    CTG   = SINTE/(1.0-COSTE);      //  CTG(THETA/2)


    TETANP = THETA0+DT;
    TETANM = THETA0-DT;

    if (THETA >= TETANM) {
        TGP  = tan(TETANP*0.5);
        TGM  = tan(TETANM*0.5);
        TGM2 = TGM*TGM;
        TGP2 = TGP*TGP;
    }

    COSM1 = 1.0;
    SINM1 = 0.0;
    TM    = 1.0;
    TGM2M = 1.0;
    TGP2M = 1.0;

    for ( M=1; M<= N; M++ ) {

    TM *= TG;
    CCOS[M] = COSM1*COSFI-SINM1*SINFI;
    SSIN[M] = SINM1*COSFI+COSM1*SINFI;
    COSM1   = CCOS[M];
    SINM1   = SSIN[M];
    if (THETA < TETANM) {
        T    = TM;
        DTT  = 0.5*M*TM*(TG+CTG);
        DTT0 = 0.0;
    } else if (THETA < TETANP) {
        TGM2M  = TGM2M*TGM2;
        FC     = 1.0/(TGP-TGM);
        FC1    = 1.0/(2*M+1);
        TGM2M1 = TGM2M*TGM;
        TG21   = 1.0+TG*TG;
        T      = FC*(TM*(TGP-TG)+FC1*(TM*TG-TGM2M1/TM));
        DTT    = 0.5*M*FC*TG21*(TM/TG*(TGP-TG)-FC1*(TM-TGM2M1/(TM*TG)));
        DTT0   = 0.5*FC*((TGP+TGM)*(TM*TG-FC1*(TM*TG-TGM2M1/TM))+
                TM*(1.0-TGP*TGM)-(1.0+TGM2)*TGM2M/TM);
    } else {
        TGP2M = TGP2M*TGP2;
        TGM2M = TGM2M*TGM2;
        FC    = 1.0/(TGP-TGM);
        FC1   = 1.0/(2*M+1);
        T     = FC*FC1*(TGP2M*TGP-TGM2M*TGM)/TM;
        DTT   = -T*M*0.50*(TG+CTG);
    }

    BTN[M] = M*T*CCOS[M]/RO;
    BPN[M] = -DTT*SSIN[M]/R;

    }

    *BTHETA = BTN[N]*800.0;
    *BPHI   = BPN[N]*800.0;

}

void    T01S_BIRK_SHL( int J, int PSChanged, int XChanged, int YChanged, int ZChanged, double *A, double PS, double X_SC, double X, double Y, double Z, double *BX, double *BY, double *BZ, LgmTsyg2001_Info *t ) {

    int        L, M, I, N, NN, K;
    double    GX, GY, GZ, FX, FY, FZ, HX, HY, HZ, HXR, HZR;


    /*
     *  There are 4 different sets of coefficients that come to us. We added an argument to the
     *  calling routine so that we know which set we are dealing with. Then many things dont need to
     *  be recomputed. Probably dont need to cache all quantities here -- i.e. could prob optimize memory
     *  usage as well...
     */
    for ( I=1;  I<=3; I++ ) {
        if ( !t->DoneJ[J] ) {

            t->P[J][I] = A[72+I]; t->ooP[J][I] = 1.0/t->P[J][I]; t->ooP2[J][I] = t->ooP[J][I]*t->ooP[J][I];
            t->Q[J][I] = A[78+I]; t->ooQ[J][I] = 1.0/t->Q[J][I]; t->ooQ2[J][I] = t->ooQ[J][I]*t->ooQ[J][I];

            t->R[J][I] = A[75+I]; t->ooR[J][I] = 1.0/t->R[J][I]; t->ooR2[J][I] = t->ooR[J][I]*t->ooR[J][I];
            t->S[J][I] = A[81+I]; t->ooS[J][I] = 1.0/t->S[J][I]; t->ooS2[J][I] = t->ooS[J][I]*t->ooS[J][I];
        }
    }

    /*
     *  if PS has changed we need to recompute things
     */
    if ( PSChanged || !t->DoneJ[J] ) {
        t->CPS  = t->cos_psi;
        t->SPS  = t->sin_psi;
        t->S3PS = 2.0*t->CPS;
        t->PST1[J] = PS*A[85];
        t->PST2[J] = PS*A[86];
        t->ST1[J] = sin(t->PST1[J]);
        t->CT1[J] = cos(t->PST1[J]);
        t->ST2[J] = sin(t->PST2[J]);
        t->CT2[J] = cos(t->PST2[J]);
    }


    /*
     *  if X or Z or PS have changed we need to recompute things
     */
    if ( XChanged || ZChanged || PSChanged || !t->DoneJ[J] ) {
        t->X1[J] = X*t->CT1[J] - Z*t->ST1[J];
        t->Z1[J] = X*t->ST1[J] + Z*t->CT1[J];
        t->X2[J] = X*t->CT2[J] - Z*t->ST2[J];
        t->Z2[J] = X*t->ST2[J] + Z*t->CT2[J];
    }


    /*
     *  cache the trig stuff. This doesnt depend on PS
     */
    if ( YChanged || !t->DoneJ[J] ) {
        for (I=1; I<=3; I++ ){
            t->YooP[J][I] = Y*t->ooP[J][I];
            t->YooQ[J][I] = Y*t->ooQ[J][I];
            t->CYPI[J][I] = cos(t->YooP[J][I]);
            t->CYQI[J][I] = cos(t->YooQ[J][I]);
            t->SYPI[J][I] = sin(t->YooP[J][I]);
            t->SYQI[J][I] = sin(t->YooQ[J][I]);
        }
    }

    if ( XChanged || ZChanged || PSChanged || !t->DoneJ[J] ) {
        for (K=1; K<=3; K++ ){
            t->Z1ooR[J][K] = t->Z1[J]*t->ooR[J][K];
            t->Z2ooS[J][K] = t->Z2[J]*t->ooS[J][K];
            t->SZRK[J][K] = sin(t->Z1ooR[J][K]);
            t->CZSK[J][K] = cos(t->Z2ooS[J][K]);
            t->CZRK[J][K] = cos(t->Z1ooR[J][K]);
            t->SZSK[J][K] = sin(t->Z2ooS[J][K]);
        }
    }

    if ( XChanged || YChanged || ZChanged || PSChanged || !t->DoneJ[J] ) {
        for (I=1; I<=3; I++ ){
            for (K=1; K<=3; K++ ){
                t->SQPR[J][I][K] = sqrt(t->ooP2[J][I] + t->ooR2[J][K]);
                t->SQQS[J][I][K] = sqrt(t->ooQ2[J][I] + t->ooS2[J][K]);
                t->EPR[J][I][K]  = exp(t->X1[J]*t->SQPR[J][I][K]);
                t->EQS[J][I][K]  = exp(t->X2[J]*t->SQQS[J][I][K]);
            }
        }
    }


    /*
     *  Flag that we've been through once already with this value of J
     */
    t->DoneJ[J] = TRUE;
t->DoneJ[J] = FALSE;


    L  = 0;
    GX = 0.0;
    GY = 0.0;
    GZ = 0.0;


    for ( M=1;  M<=2; M++ ) {    //  M=1 IS FOR THE 1ST SUM ("PERP." SYMMETRY)
                //  AND M=2 IS FOR THE SECOND SUM ("PARALL." SYMMETRY)

        for ( I=1;  I<=3; I++ ) {

            for ( K=1; K<= 3; K++ ) {

                for (N=1; N<=2; N++ ) {        // N=1 IS FOR THE FIRST PART OF EACH COEFFICIENT
                                // AND N=2 IS FOR THE SECOND ONE

                    for (NN=1; NN<=2; NN++ ) {        // NN = 1,2 FURTHER SPLITS THE COEFFICIENTS INTO 2 PARTS,
                                    // TO TAKE INTO ACCOUNT THE SCALE FACTOR DEPENDENCE
                        if (M == 1) {

                            FX = -t->SQPR[J][I][K] * t->EPR[J][I][K] * t->CYPI[J][I] * t->SZRK[J][K];
                            FY =  t->EPR[J][I][K]  * t->SYPI[J][I]   * t->SZRK[J][K] * t->ooP[J][I];
                            FZ = -t->EPR[J][I][K]  * t->CYPI[J][I]   * t->CZRK[J][K] * t->ooR[J][K];

                            if (N == 1) {
                            if (NN == 1) {
                                HX = FX;
                                HY = FY;
                                HZ = FZ;
                            } else {
                                HX = FX*X_SC;
                                HY = FY*X_SC;
                                HZ = FZ*X_SC;
                            }
                            } else {
                            if (NN == 1) {
                                HX = FX*t->CPS;
                                HY = FY*t->CPS;
                                HZ = FZ*t->CPS;
                            } else {
                                HX = FX*t->CPS*X_SC;
                                HY = FY*t->CPS*X_SC;
                                HZ = FZ*t->CPS*X_SC;
                            }
                            }

                        } else {    //   M.EQ.2

                            FX = -t->SPS * t->SQQS[J][I][K] * t->EQS[J][I][K] * t->CYQI[J][I] * t->CZSK[J][K];
                            FY =  t->SPS * t->ooQ[J][I]     * t->EQS[J][I][K] * t->SYQI[J][I] * t->CZSK[J][K];
                            FZ =  t->SPS * t->ooS[J][K]     * t->EQS[J][I][K] * t->CYQI[J][I] * t->SZSK[J][K];

                            if (N == 1) {
                            if (NN == 1) {
                                HX = FX;
                                HY = FY;
                                HZ = FZ;
                            } else {
                                HX = FX*X_SC;
                                HY = FY*X_SC;
                                HZ = FZ*X_SC;
                            }
                            } else {
                            if (NN == 1) {
                                HX = FX*t->S3PS;
                                HY = FY*t->S3PS;
                                HZ = FZ*t->S3PS;
                            } else {
                                HX = FX*t->S3PS*X_SC;
                                HY = FY*t->S3PS*X_SC;
                                HZ = FZ*t->S3PS*X_SC;
                            }
                            }

                        }

                        ++L;

                        if (M == 1) {
                            HXR =  HX*t->CT1[J] + HZ*t->ST1[J];
                            HZR = -HX*t->ST1[J] + HZ*t->CT1[J];
                        } else {
                            HXR =  HX*t->CT2[J] + HZ*t->ST2[J];
                            HZR = -HX*t->ST2[J] + HZ*t->CT2[J];
                        }

                        GX = GX + HXR*A[L];
                        GY = GY + HY *A[L];
                        GZ = GZ + HZR*A[L];


                    }
                }
            }
        }
    }

    *BX = GX;
    *BY = GY;
    *BZ = GZ;


    return;

}



void     T01S_FULL_RC( int IOPR, double PS, double X, double Y, double Z,
    double *BXSRC, double *BYSRC, double *BZSRC, double *BXPRC, double *BYPRC,  double *BZPRC, LgmTsyg2001_Info *t ) {

    /*
     *
     *   CALCULATES GSM FIELD COMPONENTS OF THE SYMMETRIC (SRC) AND PARTIAL (PRC) COMPONENTS OF THE RING CURRENT
     *   SRC  PROVIDES A DEPRESSION OF -28 nT AT EARTH
     *   PRC  CORRESPONDS TO THE PRESSURE DIFFERENCE OF 2 nPa BETWEEN MIDNIGHT AND NOON RING CURRENT
     *             PARTICLE PRESSURE AND YIELDS A DEPRESSION OF -17 nT AT X=-6Re
     *
     *   SC_SY AND SC_PR ARE SCALING FACTORS FOR THE SYMMETRIC AND PARTIAL COMPONENTS:
     *          VALUES LARGER THAN 1 RESULT IN SPATIALLY LARGER CURRENTS
     *
     *   PHI IS THE ROTATION ANGLE IN RADIANS OF THE PARTIAL RING CURRENT (MEASURED FROM MIDNIGHT TOWARD DUSK)
     *
     *     IOPR -  A RING CURRENT CALCULATION FLAG (FOR LEAST-SQUARES FITTING ONLY):
     *             IOPR=0 - BOTH SRC AND PRC FIELDS ARE CALCULATED
     *             IOPR=1 - SRC ONLY
     *             IOPR=2 - PRC ONLY
     */


    static double C_SY[] = { -9e99, -957.2534900,-817.5450246,583.2991249,758.8568270,     //   CORRECTED VALUES (AS OF MAY 2006)
                         13.17029064,68.94173502,-15.29764089,-53.43151590,27.34311724,
                         149.5252826,-11.00696044,-179.7031814,953.0914774,817.2340042,
                         -581.0791366,-757.5387665,-13.10602697,-68.58155678,15.22447386,
                         53.15535633,-27.07982637,-149.1413391,10.91433279,179.3251739,
                         -6.028703251,1.303196101,-1.345909343,-1.138296330,-0.06642634348,
                         -0.3795246458,.07487833559,.2891156371,-.5506314391,-.4443105812,
                         0.2273682152,0.01086886655,-9.130025352,1.118684840,1.110838825,
                         .1219761512,-.06263009645,-.1896093743,.03434321042,.01523060688,
                         -.4913171541,-.2264814165,-.04791374574,.1981955976,-68.32678140,
                         -48.72036263,14.03247808,16.56233733,2.369921099,6.200577111,
                         -1.415841250,-0.8184867835,-3.401307527,-8.490692287,3.217860767,
                         -9.037752107,66.09298105,48.23198578,-13.67277141,-16.27028909,
                         -2.309299411,-6.016572391,1.381468849,0.7935312553,3.436934845,
                          8.260038635,-3.136213782,8.833214943,8.041075485,8.024818618,
                          35.54861873,12.55415215,1.738167799,3.721685353,23.06768025,
                          6.871230562,6.806229878,21.35990364,1.687412298,3.500885177,
                          0.3498952546,0.6595919814};

    static double C_PR[] = { -9e99, -64820.58481,-63965.62048,66267.93413,135049.7504,
                         -36.56316878,124.6614669,56.75637955,-87.56841077,5848.631425,
                         4981.097722,-6233.712207,-10986.40188,68716.52057,65682.69473,
                         -69673.32198,-138829.3568,43.45817708,-117.9565488,-62.14836263,
                         79.83651604,-6211.451069,-5151.633113,6544.481271,11353.03491,
                         23.72352603,-256.4846331,25.77629189,145.2377187,-4.472639098,
                         -3.554312754,2.936973114,2.682302576,2.728979958,26.43396781,
                         -9.312348296,-29.65427726,-247.5855336,-206.9111326,74.25277664,
                         106.4069993,15.45391072,16.35943569,-5.965177750,-6.079451700,
                         115.6748385,-35.27377307,-32.28763497,-32.53122151,93.74409310,
                         84.25677504,-29.23010465,-43.79485175,-6.434679514,-6.620247951,
                         2.443524317,2.266538956,-43.82903825,6.904117876,12.24289401,
                         17.62014361,152.3078796,124.5505289,-44.58690290,-63.02382410,
                         -8.999368955,-9.693774119,3.510930306,3.770949738,-77.96705716,
                         22.07730961,20.46491655,18.67728847,9.451290614,9.313661792,
                         644.7620970,418.2515954,7.183754387,35.62128817,19.43180682,
                         39.57218411,15.69384715,7.123215241,2.300635346,21.90881131,
                         -.01775839370,.3996346710};

    double    HXSRC, HYSRC, HZSRC, HXPRC, HYPRC, HZPRC, X_SC, FSX, FSY, FSZ;
    double    FPX, FPY, FPZ;

    T01S_SRC_PRC( IOPR, t->CB_RCPAR.SC_SY, t->CB_RCPAR.SC_AS, t->CB_RCPAR.PHI, PS, X, Y, Z, &HXSRC, &HYSRC, &HZSRC, &HXPRC, &HYPRC, &HZPRC, t );
//printf("HXSRC, HYSRC, HZSRC, HXPRC, HYPRC, HZPRC = %g %g %g %g %g %g \n", HXSRC, HYSRC, HZSRC, HXPRC, HYPRC, HZPRC);

    X_SC = t->CB_RCPAR.SC_SY-1.0;
    if ( (IOPR == 0) || (IOPR == 1) ) {
        T01S_RC_SHIELD( C_SY, PS, X_SC, X, Y, Z, &FSX, &FSY, &FSZ, t );
    } else {
        FSX = 0.0;
        FSY = 0.0;
        FSZ = 0.0;
    }

    X_SC = t->CB_RCPAR.SC_AS-1.0;
    if ( (IOPR == 0) || (IOPR == 2) ) {
        T01S_RC_SHIELD( C_PR, PS, X_SC, X, Y, Z, &FPX, &FPY, &FPZ, t );
    } else {
        FPX = 0.0;
        FPY = 0.0;
        FPZ = 0.0;
    }

    *BXSRC = HXSRC + FSX;
    *BYSRC = HYSRC + FSY;
    *BZSRC = HZSRC + FSZ;

    *BXPRC = HXPRC + FPX;
    *BYPRC = HYPRC + FPY;
    *BZPRC = HZPRC + FPZ;

    return;
}


void     T01S_SRC_PRC( int IOPR, double SC_SY, double SC_PR, double PHI, double PS, double X, double Y, double Z,
        double *BXSRC, double *BYSRC, double *BZSRC, double *BXPRC, double *BYPRC, double *BZPRC, LgmTsyg2001_Info *t ) {
    /*
     *   RETURNS FIELD COMPONENTS FROM A MODEL RING CURRENT, INCLUDING ITS SYMMETRIC PART
     *     AND A PARTIAL RING CURRENT, CLOSED VIA BIRKELAND CURRENTS. BASED ON RESULTS, DESCRIBED
     *     IN A PT01S_APER "MODELING THE INNER MAGNETOSPHERE: ASYMMETRIC RING CURRENT AND REGION 2
     *     BIRKELAND CURRENTS REVISITED" (JGR, DEC.2000).
     *
     *     IOPR -  A RING CURRENT CALCULATION FLAG (FOR LEAST-SQUARES FITTING ONLY):
     *             IOPR=0 - BOTH SRC AND PRC FIELDS ARE CALCULATED
     *             IOPR=1 - SRC ONLY
     *             IOPR=2 - PRC ONLY
     *
     *     SC_SY &  SC_PR ARE SCALE FACTORS FOR THE ABOVE COMPONENTS;  TAKING SC<1 OR SC>1 MAKES THE CURRENTS
     *                      SHRINK OR EXPAND, RESPECTIVELY.
     *
     *   PHI IS THE ROTATION ANGLE (RADIANS) OF THE PARTIAL RING CURRENT (MEASURED FROM MIDNIGHT TOWARD DUSK)
     *
     *
     *   1.  TRANSFORM TO TILTED COORDINATES (i.e., SM coordinates):
     */

    double    CPS, SPS, XT, XTS, ooSC_SY, YTS, ZTS, ZT, ooSC_PR, XTA, YTA, ZTA, BXS, BYS, BZS;
    double    BXA_S, BYA_S, BZA_S, CP, SP, XR, YR, BXA_QR, BYA_QR, BZA_Q, BXA_Q, BYA_Q, BXP, BYP, BZP;


    //CPS = cos(PS); SPS = sin(PS);
    CPS = t->cos_psi; SPS = t->sin_psi;

    XT = X*CPS - Z*SPS;
    ZT = Z*CPS + X*SPS;


    /*
     *     2.  SCALE THE COORDINATES FOR THE SYMMETRIC AND PARTIAL RC COMPONENTS:
     */
    ooSC_SY = 1.0/SC_SY;
    XTS = XT*ooSC_SY;    //  SYMMETRIC
    YTS = Y *ooSC_SY;
    ZTS = ZT*ooSC_SY;

    ooSC_PR = 1.0/SC_PR;
    XTA = XT*ooSC_PR;    //  PARTIAL
    YTA = Y *ooSC_PR;
    ZTA = ZT*ooSC_PR;


    /*
     *     3.  CALCULATE COMPONENTS OF THE TOTAL FIELD IN THE TILTED (SOLAR-MAGNETIC) COORDINATE SYSTEM:
     *
     *
     *      3a. SYMMETRIC FIELD:
     */
    if (IOPR <= 1)                      T01S_RC_SYMM( XTS, YTS, ZTS, &BXS, &BYS, &BZS);
    if ( (IOPR == 0) || (IOPR == 2) )   PT01S_RC_SYMM( XTA, YTA, ZTA, &BXA_S, &BYA_S, &BZA_S);

    /*     3b. ROTATE THE SCALED SM COORDINATES BY PHI AROUND ZSM AXIS AND CALCULATE QUADRUPOLE PRC FIELD
     *         IN THOSE COORDS:
     */
    CP = cos(PHI);
    SP = sin(PHI);
    XR = XTA*CP - YTA*SP;
    YR = XTA*SP + YTA*CP;
//printf("PHI, CP, SP, XR, YR = %g %g %g %g %g\n", PHI, CP, SP, XR, YR);

    if ( (IOPR == 0) || (IOPR == 2) )   T01S_PRC_QUAD( XR, YR, ZTA, &BXA_QR, &BYA_QR, &BZA_Q );

    /*     3c. TRANSFORM THE QUADRUPOLE FIELD COMPONENTS BACK TO THE SM COORDS:
     */
    BXA_Q =  BXA_QR*CP + BYA_QR*SP;
    BYA_Q = -BXA_QR*SP + BYA_QR*CP;

    /*     3d. FIND THE TOTAL FIELD OF PRC (SYMM.+QUADR.) IN THE SM COORDS:
     */
    BXP = BXA_S + BXA_Q;
    BYP = BYA_S + BYA_Q;
    BZP = BZA_S + BZA_Q;



    /*
     *     4.  TRANSFORM THE FIELDS OF BOTH PARTS OF THE RING CURRENT BACK TO THE GSM SYSTEM:
     */
    *BXSRC = BXS*CPS+BZS*SPS;   //     SYMMETRIC RC
    *BYSRC = BYS;
    *BZSRC = BZS*CPS-BXS*SPS;

    *BXPRC = BXP*CPS+BZP*SPS;   //     PARTIAL RC
    *BYPRC = BYP;
    *BZPRC = BZP*CPS-BXP*SPS;

    return;

}






void    T01S_RC_SYMM( double X, double Y, double Z, double *BX, double *BY, double *BZ){

    double    DS=1e-2, DC=0.99994999875, D=1e-4, DRD=5e3;     //  DS=SIN(THETA) AT THE BOUNDARY OF THE LINEARITY
                                //  REGION; DC = SQRT(1-DS**2);  DRD=1/(2*D)

    double    RHO2, R2, R, ooR, RP, RM, SINT, COST, A, DARDR, FXY, THETA;
    double    TP, TM, SINTP, SINTM, COSTP, COSTM, BR, BT;

    RHO2 = X*X+Y*Y;
    R2   = RHO2+Z*Z;
    R    = sqrt(R2);
    ooR  = 1.0/R;
    RP   = R+D;
    RM   = R-D;
    SINT = sqrt(RHO2)*ooR;
    COST = Z*ooR;

    if ( SINT < DS ) {    // TOO CLOSE TO THE Z-AXIS; USING A LINEAR T01S_APPROXIMATION A_PHI~SINT,
                          // TO AVOID THE SINGULARITY PROBLEM
        A = T01S_AP(R,DS,DC)/DS;
        DARDR = (RP*T01S_AP(RP,DS,DC)-RM*T01S_AP(RM,DS,DC))*DRD;
        FXY   = Z*(2.0*A-DARDR)/(R*R2);
        *BX    = FXY*X;
        *BY    = FXY*Y;
        *BZ    = (2.0*A*COST*COST+DARDR*SINT*SINT)*ooR;

    } else {

        THETA = atan2(SINT,COST);
        TP    = THETA+D;
        TM    = THETA-D;
        SINTP = sin(TP);
        SINTM = sin(TM);
        COSTP = cos(TP);
        COSTM = cos(TM);
        BR    = (SINTP*T01S_AP(R,SINTP,COSTP)-SINTM*T01S_AP(R,SINTM,COSTM))/(R*SINT)*DRD;
        BT    = (RM*T01S_AP(RM,SINT,COST)-RP*T01S_AP(RP,SINT,COST))/R*DRD;
        FXY   = (BR+BT*COST/SINT)/R;
        *BX    = FXY*X;
        *BY    = FXY*Y;
        *BZ    = BR*COST-BT*SINT;

    }

    return;

}




double    T01S_AP( double R, double SINT, double COST) {

    /*
     *  Calculates azimuthal component of the vector potential of the symmetric
     *  part of the model ring current.
     */
    int        PROX;    // INDICATES WHETHER WE ARE TOO CLOSE TO THE AXIS OF SYMMETRY, WHERE THE INVERSION
                        // OF DIPOLAR COORDINATES BECOMES INACCURATE

    // CORRECTED VALUES (UPDATED 04/20/06 (SEE NB#9, P.37))
    double      A1=-456.5289941, A2=375.9055332, RRC1=4.274684950,   DD1=2.439528329;
    double    RRC2=3.367557287, DD2=3.146382545,   P1=-0.2291904607, R1=3.746064740,  DR1=1.508802177;
    double    DLA1=0.5873525737, P2=0.1556236119,  R2=4.993638842,  DR2=3.324180497, DLA2=0.4368407663;
    double      P3=0.1855957207, R3=2.969226745,  DR3=2.243367377;
    double    ooR, ooR2, SINT1, COST1, ALPHA, GAMMA, RR1oDR1, RR1oDR1_2, RR2oDR2, RR2oDR2_2;
    double    RR3oDR3, RR3oDR3_2, COST1oDLA1, COST1oDLA1_2, COST1oDLA2, COST1oDLA2_2, ARG1;
    double    ARG2, ARG3, DEXP1, DEXP2, DEXP3, ALPHA_S, ALPHA_S2, GAMMA_S, GAMMAS2, ALSQH2;
    double    ALSQH, F, Q, C, G, RS, COSTS, SINTS, RHOS, RHOS2, ZS, RRC1pRHOS, RRC1pRHOS2;
    double    P, XK2, XK, XKRHO12, XK2S, DL, ELK, ELE, T01S_APHI1, RRC2pRHOS, RRC2pRHOS2, T01S_APHI2, Result;

    ooR = 1.0/R;
    ooR2 = ooR*ooR;

    PROX  = FALSE;
    SINT1 = SINT;
    COST1 = COST;
    if (SINT1 < 1e-2) {            // TOO CLOSE TO Z-AXIS;  USE LINEAR INTERPOLATION BETWEEN SINT=0 & SINT=0.01
    SINT1 = 1e-2;
    COST1 = .99994999875;
    PROX  = TRUE;
    }

    ALPHA = SINT1*SINT1*ooR;        // R,THETA -> ALPHA,GAMMA
    GAMMA = COST1*ooR2;


    RR1oDR1 = ((R-R1)/DR1); RR1oDR1_2 = RR1oDR1*RR1oDR1;
    RR2oDR2 = ((R-R2)/DR2); RR2oDR2_2 = RR2oDR2*RR2oDR2;
    RR3oDR3 = ((R-R3)/DR3); RR3oDR3_2 = RR3oDR3*RR3oDR3;
    COST1oDLA1 = COST1/DLA1; COST1oDLA1_2 = COST1oDLA1*COST1oDLA1;
    COST1oDLA2 = COST1/DLA2; COST1oDLA2_2 = COST1oDLA2*COST1oDLA2;
    ARG1 = -RR1oDR1_2-COST1oDLA1_2;
    ARG2 = -RR2oDR2_2-COST1oDLA2_2;
    ARG3 = -RR3oDR3_2;

    DEXP1 = (ARG1 < -500.0) ? 0.0 : exp(ARG1); // TO PREVENT "FLOATING UNDERFLOW" CRASHES
    DEXP2 = (ARG2 < -500.0) ? 0.0 : exp(ARG2);
    DEXP3 = (ARG3 < -500.0) ? 0.0 : exp(ARG3);

    ALPHA_S = ALPHA*(1.0+P1*DEXP1+P2*DEXP2+P3*DEXP3);    // ALPHA -> ALPHA_S  (T01S_DEFORMED)
    ALPHA_S2 = ALPHA_S*ALPHA_S;


    GAMMA_S = GAMMA;
    GAMMAS2 = GAMMA_S*GAMMA_S;





    ALSQH = 0.5*ALPHA_S2;        //  ALPHA_S,GAMMA_S -> RS,SINTS,COSTS
    ALSQH2 = ALSQH*ALSQH;
    F = 64.0/27.0*GAMMAS2 + ALSQH2;
    Q = pow( sqrt(F)+ALSQH, 1.0/3.0 );
    C = Q-4.0*pow( GAMMAS2, 1.0/3.0 )/(3.0*Q);



    if (C < 0.0) C = 0.0;

    G  = sqrt(C*C+4.0*pow( GAMMAS2, 1.0/3.0 ));
    RS = 4.0/((sqrt(2.0*G-C)+sqrt(C))*(G+C));
    COSTS = GAMMA_S*RS*RS;
    SINTS = sqrt(1.0-COSTS*COSTS);
    RHOS  = RS*SINTS;
    RHOS2 = RHOS*RHOS;
    ZS    = RS*COSTS;


    /*
     *    1st loop:
     */

    RRC1pRHOS = RRC1+RHOS; RRC1pRHOS2 = RRC1pRHOS*RRC1pRHOS;
    P   = RRC1pRHOS2+ZS*ZS+DD1*DD1;
    XK2 = 4.0*RRC1*RHOS/P;
    XK  = sqrt(XK2);
    XKRHO12 = XK*sqrt(RHOS);            //   SEE NB#4, P.3

    XK2S = 1.0-XK2;
    DL   = log(1.0/XK2S);
    ELK  = 1.38629436112+XK2S*(0.09666344259+XK2S*(0.03590092383+
        XK2S*(0.03742563713+XK2S*0.01451196212))) +DL*
        (0.5+XK2S*(0.12498593597+XK2S*(0.06880248576+
        XK2S*(0.03328355346+XK2S*0.00441787012))));
    ELE  = 1.0+XK2S*(0.44325141463+XK2S*(0.0626060122+XK2S*
        (0.04757383546+XK2S*0.01736506451))) +DL*
        XK2S*(0.2499836831+XK2S*(0.09200180037+XK2S*
        (0.04069697526+XK2S*0.00526449639)));

    T01S_APHI1 = ((1.0-XK2*0.5)*ELK-ELE)/XKRHO12;



    /*
     *    2nd loop:
     */

    RRC2pRHOS = RRC2+RHOS; RRC2pRHOS2 = RRC2pRHOS*RRC2pRHOS;
    P   = RRC2pRHOS2+ZS*ZS+DD2*DD2;
    XK2 = 4.0*RRC2*RHOS/P;
    XK  = sqrt(XK2);
    XKRHO12 = XK*sqrt(RHOS);        //   SEE NB#4, P.3

    XK2S = 1.0-XK2;
    DL   = log(1.0/XK2S);
    ELK  = 1.38629436112+XK2S*(0.09666344259+XK2S*(0.03590092383+
        XK2S*(0.03742563713+XK2S*0.01451196212))) +DL*
        (0.5+XK2S*(0.12498593597+XK2S*(0.06880248576+
        XK2S*(0.03328355346+XK2S*0.00441787012))));
    ELE  = 1.0+XK2S*(0.44325141463+XK2S*(0.0626060122+XK2S*
        (0.04757383546+XK2S*0.01736506451))) +DL*
        XK2S*(0.2499836831+XK2S*(0.09200180037+XK2S*
        (0.04069697526+XK2S*0.00526449639)));

    T01S_APHI2 = ((1.0-XK2*0.5)*ELK-ELE)/XKRHO12;

    Result = A1*T01S_APHI1+A2*T01S_APHI2;
    if (PROX) Result *= SINT/SINT1;        // LINEAR INTERPOLATION, IF TOO CLOSE TO THE Z-AXIS


    return( Result ) ;


}









void    PT01S_RC_SYMM( double X,double Y,double Z, double *BX, double *BY, double *BZ ) {

    double    DS=1e-2, DC=0.99994999875, D=1e-4, DRD=5e3;   // DS=SIN(THETA) AT THE BOUNDARY OF THE LINEARITY
                                                            // REGION; DC=SQRT(1-DS**2);  DRD=1/(2*D)
    double    RHO2, R2, R, ooR, RP, RM, SINT, COST, A, DARDR, FXY, THETA;
    double    TP, TM, SINTP, SINTM, COSTP, COSTM, BR, BT;

    RHO2 = X*X+Y*Y;
    R2   = RHO2+Z*Z;
    R    = sqrt(R2); ooR = 1.0/R;
    RP   = R+D;
    RM   = R-D;
    SINT = sqrt(RHO2)*ooR;
    COST = Z*ooR;

    if (SINT < DS) {    //  TOO CLOSE TO THE Z-AXIS; USING A LINEAR T01S_APPROXIMATION A_PHI~SINT,
                        //  TO AVOID THE SINGULARITY PROBLEM
        A     = T01S_T01S_APPRC(R,DS,DC)/DS;
        DARDR = (RP*T01S_T01S_APPRC(RP,DS,DC)-RM*T01S_T01S_APPRC(RM,DS,DC))*DRD;
        FXY   = Z*(2.0*A-DARDR)/(R*R2);
        *BX    = FXY*X;
        *BY    = FXY*Y;
        *BZ    = (2.0*A*COST*COST+DARDR*SINT*SINT)*ooR;

    } else {

        THETA = atan2(SINT,COST);
        TP    = THETA+D;
        TM    = THETA-D;
        SINTP = sin(TP);
        SINTM = sin(TM);
        COSTP = cos(TP);
        COSTM = cos(TM);
        BR    = (SINTP*T01S_T01S_APPRC(R,SINTP,COSTP)-SINTM*T01S_T01S_APPRC(R,SINTM,COSTM))/(R*SINT)*DRD;
        BT    = (RM*T01S_T01S_APPRC(RM,SINT,COST)-RP*T01S_T01S_APPRC(RP,SINT,COST))*ooR*DRD;
        FXY   = (BR+BT*COST/SINT)*ooR;
        *BX    = FXY*X;
        *BY    = FXY*Y;
        *BZ    = BR*COST-BT*SINT;

    }

    return;

}

double    T01S_T01S_APPRC( double R, double SINT, double COST) {

    /*
     *        Calculates azimuthal component of the vector potential of the symmetric
     *    part of the model PARTIAL ring current.
     */
    int        PROX;



    double A1=-80.11202281, A2=12.58246758, RRC1=6.560486035, DD1=1.930711037, RRC2=3.827208119;
    double DD2=.7789990504, P1=.3058309043, ALPHA1=.1817139853, DAL1=.1257532909, BETA1=3.422509402;
    double DG1=.04742939676, P2=-4.800458958, ALPHA2=-.02845643596, DAL2=.2188114228, BETA2=2.545944574;
    double DG2=.00813272793, BETA3=.35868244, P3=103.1601001, ALPHA3=-.00764731187, DAL3=.1046487459;
    double BETA4=2.958863546, DG3=.01172314188, BETA5=.4382872938, Q0=.01134908150, Q1=14.51339943;
    double ALPHA4=.2647095287, DAL4=.07091230197, DG4=.01512963586, Q2=6.861329631, ALPHA5=.1677400816;
    double DAL5=.04433648846, DG5=.05553741389, BETA6=.7665599464, BETA7=.7277854652;



    double SINT1, COST1, ALPHA, GAMMA, GAMMAoDG1, GAMMAoDG12, ARG1;
    double ALPHAmALPHA4oDAL4, ALPHAmALPHA4oDAL42, GAMMAoDG4, GAMMAoDG42;
    double ARG2, DEXP1, DEXP2, ALPHAmALPHA1oDAL1, ALPHAmALPHA1oDAL12, ALPHAmALPHA2oDAL2;
    double ALPHAmALPHA2oDAL22, ALPHAmALPHA3oDAL3, ALPHAmALPHA3oDAL32, ALPHAmALPHA3, ALPHAmALPHA32;
    double GAMMAoDG2, GAMMAoDG22, GAMMAoDG3, GAMMAoDG32, ALPHA_S, GAMMAoDG5, GAMMAoDG52;
    double ALPHApALPHA5oDAL5, ALPHApALPHA5oDAL52, GAMMA_S, GAMMAS2, ALSQH, F, Q, C, G, RS;
    double COSTS, SINTS, RHOS, RHOS2, ZS, RRC1pRHOS, RRC1pRHOS2, P, XK2, XK, XKRHO12, XK2S, DL;
    double ELK, ELE, T01S_APHI1, RRC2pRHOS, RRC2pRHOS2, T01S_APHI2, Result;



    PROX  = FALSE;
    SINT1 = SINT;
    COST1 = COST;
    if (SINT1 < 1e-2) {        // TOO CLOSE TO Z-AXIS;  USE LINEAR INTERPOLATION BETWEEN SINT=0 & SINT=0.01
    SINT1 = 1e-2;
    COST1 = .99994999875;
    PROX  = TRUE;
    }

    ALPHA = SINT1*SINT1/R;    // R,THETA -> ALPHA,GAMMA
    GAMMA = COST1/(R*R);

    GAMMAoDG1 = GAMMA/DG1; GAMMAoDG12 = GAMMAoDG1*GAMMAoDG1;
    ARG1 = -GAMMAoDG12;
    ALPHAmALPHA4oDAL4 = (ALPHA-ALPHA4)/DAL4; ALPHAmALPHA4oDAL42 = ALPHAmALPHA4oDAL4*ALPHAmALPHA4oDAL4;
    GAMMAoDG4 = GAMMA/DG4; GAMMAoDG42 = GAMMAoDG4*GAMMAoDG4;
    ARG2 = -ALPHAmALPHA4oDAL42-GAMMAoDG42;


    DEXP1 = (ARG1 < -500.0) ? 0.0 : exp(ARG1);    // TO PREVENT "FLOATING UNDERFLOW" CRASHES
    DEXP2 = (ARG2 < -500.0) ? 0.0 : exp(ARG2);    // TO PREVENT "FLOATING UNDERFLOW" CRASHES



    ALPHAmALPHA1oDAL1 = (ALPHA-ALPHA1)/DAL1; ALPHAmALPHA1oDAL12 = ALPHAmALPHA1oDAL1*ALPHAmALPHA1oDAL1;
    ALPHAmALPHA2oDAL2 = (ALPHA-ALPHA2)/DAL2; ALPHAmALPHA2oDAL22 = ALPHAmALPHA2oDAL2*ALPHAmALPHA2oDAL2;
    ALPHAmALPHA3oDAL3 = (ALPHA-ALPHA3)/DAL3; ALPHAmALPHA3oDAL32 = ALPHAmALPHA3oDAL3*ALPHAmALPHA3oDAL3;
    ALPHAmALPHA3 = ALPHA-ALPHA3; ALPHAmALPHA32 = ALPHAmALPHA3*ALPHAmALPHA3;
    GAMMAoDG2 = GAMMA/DG2; GAMMAoDG22 = GAMMAoDG2*GAMMAoDG2;
    GAMMAoDG3 = GAMMA/DG3; GAMMAoDG32 = GAMMAoDG3*GAMMAoDG3;

    ALPHA_S = ALPHA*(1.+P1/pow(1.+ALPHAmALPHA1oDAL12, BETA1)
        *DEXP1+P2*(ALPHA-ALPHA2)/pow(1.+ALPHAmALPHA2oDAL22, BETA2)
        /pow(1.+GAMMAoDG22, BETA3)
        +P3*ALPHAmALPHA32/pow(1.+ALPHAmALPHA3oDAL32, BETA4)
        /pow(1.+GAMMAoDG32, BETA5));        // ALPHA -> ALPHA_S  (T01S_DEFORMED)



    GAMMAoDG5 = GAMMA/DG5; GAMMAoDG52 = GAMMAoDG5*GAMMAoDG5;
    ALPHApALPHA5oDAL5 = (ALPHA-ALPHA5)/DAL5; ALPHApALPHA5oDAL52 = ALPHApALPHA5oDAL5*ALPHApALPHA5oDAL5;
    GAMMA_S = GAMMA*(1.+Q0+Q1*(ALPHA-ALPHA4)*DEXP2    // GAMMA -> GAMMA_  (T01S_DEFORMED)
        +Q2*(ALPHA-ALPHA5)/pow( 1.+ALPHApALPHA5oDAL52, BETA6 )
        /pow( 1.+GAMMAoDG52, BETA7 ));

    GAMMAS2 = GAMMA_S*GAMMA_S;

    ALSQH = 0.5*ALPHA_S*ALPHA_S;        // ALPHA_S,GAMMA_S -> RS,SINTS,COSTS
    F = 64./27.*GAMMAS2+ALSQH*ALSQH;
    Q = pow( sqrt(F)+ALSQH, 1./3. );
    C = Q-4.*pow( GAMMAS2, 1./3. )/(3.*Q);


    if (C < 0.) C = 0.;

    G  = sqrt(C*C+4.*pow( GAMMAS2, 1./3. ));
    RS = 4./((sqrt(2.*G-C)+sqrt(C))*(G+C));
    COSTS = GAMMA_S*RS*RS;
    SINTS = sqrt(1.-COSTS*COSTS);
    RHOS  = RS*SINTS;
    RHOS2 = RHOS*RHOS;
    ZS = RS*COSTS;


    /*
     *  1st loop:
     */
    RRC1pRHOS = RRC1+RHOS; RRC1pRHOS2 = RRC1pRHOS*RRC1pRHOS;
    P   = RRC1pRHOS2+ZS*ZS+DD1*DD1;
    XK2 = 4.*RRC1*RHOS/P;
    XK  = sqrt(XK2);
    XKRHO12 = XK*sqrt(RHOS);

    XK2S = 1.-XK2;
    DL   = log(1./XK2S);
    ELK  = 1.38629436112+XK2S*(0.09666344259+XK2S*(0.03590092383+
        XK2S*(0.03742563713+XK2S*0.01451196212))) +DL*
        (0.5+XK2S*(0.12498593597+XK2S*(0.06880248576+
        XK2S*(0.03328355346+XK2S*0.00441787012))));
    ELE  = 1.+XK2S*(0.44325141463+XK2S*(0.0626060122+XK2S*
        (0.04757383546+XK2S*0.01736506451))) +DL*
        XK2S*(0.2499836831+XK2S*(0.09200180037+XK2S*
        (0.04069697526+XK2S*0.00526449639)));

    T01S_APHI1 = ((1.-XK2*0.5)*ELK-ELE)/XKRHO12;



    /*
     *  2nd loop:
     */

    RRC2pRHOS = RRC2+RHOS; RRC2pRHOS2 = RRC2pRHOS*RRC2pRHOS;
    P   = RRC2pRHOS2+ZS*ZS+DD2*DD2;
    XK2 = 4.*RRC2*RHOS/P;
    XK  = sqrt(XK2);
    XKRHO12 = XK*sqrt(RHOS);

    XK2S = 1.-XK2;
    DL   = log(1./XK2S);
    ELK  = 1.38629436112+XK2S*(0.09666344259+XK2S*(0.03590092383+
          XK2S*(0.03742563713+XK2S*0.01451196212))) +DL*
          (0.5+XK2S*(0.12498593597+XK2S*(0.06880248576+
          XK2S*(0.03328355346+XK2S*0.00441787012))));
    ELE  = 1.+XK2S*(0.44325141463+XK2S*(0.0626060122+XK2S*
           (0.04757383546+XK2S*0.01736506451))) +DL*
          XK2S*(0.2499836831+XK2S*(0.09200180037+XK2S*
            (0.04069697526+XK2S*0.00526449639)));

    T01S_APHI2 = ((1.-XK2*0.5)*ELK-ELE)/XKRHO12;

    Result = A1*T01S_APHI1+A2*T01S_APHI2;
    if (PROX) Result *= SINT/SINT1;        // LINEAR INTERPOLATION, IF TOO CLOSE TO THE Z-AXIS


    return( Result );


}

void    T01S_PRC_QUAD( double X, double Y, double Z, double *BX, double *BY, double *BZ) {

//    double    DD=2e-4;
    double    D=1e-4, DS=1e-2, DC=0.99994999875, ooDD=5000.0;


    double    X2, Y2, RHO2, R, ooR, RHO, SINT, COST, RP, RM, ooRHO, CPHI, SPHI;
    double    BR, BT, DBRR, THETA, TP, TM, SINTP, COSTP, SINTM, COSTM, DBTT, ST;
    double    CT, FCXY, RtST, RtST2, ooRtST2;


    X2 = X*X; Y2 = Y*Y;
    RHO2 = X2+Y2;
    R    = sqrt(RHO2+Z*Z); ooR = 1.0/R;
    RHO  = sqrt(RHO2);
    SINT = RHO*ooR;
    COST = Z*ooR;
    RP = R+D;
    RM = R-D;

    if (SINT > DS) {
        ooRHO = 1.0/RHO;
        CPHI  = X*ooRHO;
        SPHI  = Y*ooRHO;
        BR = T01S_BR_PRC_Q(R,SINT,COST);
        BT = T01S_BT_PRC_Q(R,SINT,COST);
        DBRR  = (T01S_BR_PRC_Q(RP,SINT,COST)-T01S_BR_PRC_Q(RM,SINT,COST))*ooDD;
        THETA = atan2(SINT,COST);
        TP = THETA+D;
        TM = THETA-D;
        SINTP = sin(TP);
        COSTP = cos(TP);
        SINTM = sin(TM);
        COSTM = cos(TM);
        DBTT  = (T01S_BT_PRC_Q(R,SINTP,COSTP)-T01S_BT_PRC_Q(R,SINTM,COSTM))*ooDD;
        *BX = SINT*(BR+(BR+R*DBRR+DBTT)*SPHI*SPHI)+COST*BT;
        *BY = -SINT*SPHI*CPHI*(BR+R*DBRR+DBTT);
        *BZ = (BR*COST-BT*SINT)*CPHI;
    } else {
        ST = DS;
        CT = DC;
        if (Z < 0.0) CT = -DC;
        THETA = atan2(ST,CT);
        TP = THETA+D;
        TM = THETA-D;
        SINTP = sin(TP);
        COSTP = cos(TP);
        SINTM = sin(TM);
        COSTM = cos(TM);
        BR = T01S_BR_PRC_Q(R,ST,CT);
        BT = T01S_BT_PRC_Q(R,ST,CT);
        DBRR = (T01S_BR_PRC_Q(RP,ST,CT)-T01S_BR_PRC_Q(RM,ST,CT))*ooDD;
        DBTT = (T01S_BT_PRC_Q(R,SINTP,COSTP)-T01S_BT_PRC_Q(R,SINTM,COSTM))*ooDD;
        FCXY = R*DBRR+DBTT;
        RtST = R*ST; RtST2 = RtST*RtST; ooRtST2 = 1.0/RtST2;
        *BX = (BR*(X2+2.0*Y2)+FCXY*Y2)*ooRtST2+BT*COST;
        *BY = -(BR+FCXY)*X*Y*ooRtST2;
        *BZ = (BR*COST/ST-BT)*X/R;
    }

    return;

}


double    T01S_BR_PRC_Q( double R, double SINT, double COST ) {

    /*
     *   Calculates the radial component of the "quadrupole" part of the model partial ring current.
     */

    double A1=-21.2666329;
    double A2=32.24527521, A3=-6.062894078, A4=7.515660734, A5=233.7341288, A6=-227.1195714;        // ALL LINEAR PARAMETERS HERE
    double A7=8.483233889, A8=16.80642754, A9=-24.63534184, A10=9.067120578, A11=-1.052686913;        // WERE MULTIPLIED BY 0.1,
    double A12=-12.08384538, A13=18.61969572, A14=-12.71686069, A15=47017.35679, A16=-50646.71204;    // SO THAT THEY CORRESPOND TO P_0=1 nPa,
    double A17=7746.058231,  A18=1.531069371, XK1=2.318824273, AL1=.1417519429, DAL1=.6388013110E-02;    // RATHER THAN THE ORIGINAL VALUE OF 10 nPa
    double B1=5.303934488, BE1=4.213397467, XK2=.7955534018, AL2=.1401142771, DAL2=.2306094179E-01;    // ASSUMED IN THE BIOT-SAVART INTEGRAL
    double B2=3.462235072, BE2=2.568743010, XK3=3.477425908, XK4=1.922155110, AL3=.1485233485;
    double DAL3=.2319676273E-01, B3=7.830223587,  BE3=8.492933868, AL4=.1295221828, DAL4=.01753008801;
    double DG1=.01125504083, AL5=.1811846095, DAL5=.04841237481, DG2=.01981805097, C1=6.557801891;
    double C2=6.348576071, C3=5.744436687, AL6=.2265212965, DAL6=.1301957209, DRM=.5654023158;

    double ooR, ooR2, SINT2, COST2, SC, ALPHA, GAMMA, F, FA, FS;
    double D1, D2, D3, D4, D5, D6, ALPHAmAL4oDAL4, ALPHAmAL4oDAL42;
    double GAMMAoDG1, GAMMAoDG12, ARGA, ooARGA, ARGG, D7, D8, D9, D10;
    double ALPHAmAL5oDAL5, ALPHAmAL5oDAL52, GAMMAoDG2, GAMMAoDG22, D11;
    double D12, D13, D14, R2, R4, C12, C14, C22, C24, C32, C34, D15, D16;
    double D17, RmNUMoDRM, RmNUMoDRM2, D18, Result;

    ooR = 1.0/R; ooR2 = ooR*ooR;
    SINT2 = SINT*SINT;
    COST2 = COST*COST;
    SC    = SINT*COST;
    ALPHA = SINT2*ooR;
    GAMMA = COST*ooR2;

    T01S_FFS(ALPHA,AL1,DAL1,&F,&FA,&FS);
    D1 = SC*pow(F, XK1)/(pow( R/B1, BE1) +1.0);
    D2 = D1*COST2;

    T01S_FFS(ALPHA,AL2,DAL2,&F,&FA,&FS);
    D3 = SC*pow(FS, XK2)/(pow( R/B2, BE2) +1.0);
    D4 = D3*COST2;

    T01S_FFS(ALPHA,AL3,DAL3,&F,&FA,&FS);
    D5 = SC*pow(ALPHA, XK3)*pow(FS, XK4)/(pow(R/B3, BE3)+1.0);
    D6 = D5*COST2;

    ALPHAmAL4oDAL4 = (ALPHA-AL4)/DAL4; ALPHAmAL4oDAL42 = ALPHAmAL4oDAL4*ALPHAmAL4oDAL4;
    GAMMAoDG1 = GAMMA/DG1; GAMMAoDG12 = GAMMAoDG1*GAMMAoDG1;
    ARGA = ALPHAmAL4oDAL42+1.0;    ooARGA = 1.0/ARGA;
    ARGG = 1.0+GAMMAoDG12;

    D7  = SC*ooARGA/ARGG;
    D8  = D7*ooARGA;
    D9  = D8*ooARGA;
    D10 = D9*ooARGA;

    ALPHAmAL5oDAL5 = (ALPHA-AL5)/DAL5; ALPHAmAL5oDAL52 = ALPHAmAL5oDAL5*ALPHAmAL5oDAL5;
    GAMMAoDG2 = GAMMA/DG2; GAMMAoDG22 = GAMMAoDG2*GAMMAoDG2;
    ARGA = ALPHAmAL5oDAL52+1.0; ooARGA = 1.0/ARGA;
    ARGG = 1.0+GAMMAoDG22;

    D11 = SC*ooARGA/ARGG;
    D12 = D11*ooARGA;
    D13 = D12*ooARGA;
    D14 = D13*ooARGA;


    R2 = R*R; R4 = R2*R2;
    C12 = C1*C1; C14 = C12*C12;
    C22 = C2*C2; C24 = C22*C22;
    C32 = C3*C3; C34 = C32*C32;

    D15 = SC/(R4+C14);
    D16 = SC/(R4+C24)*COST2;
    D17 = SC/(R4+C34)*COST2*COST2;

    T01S_FFS(ALPHA,AL6,DAL6,&F,&FA,&FS);

    RmNUMoDRM = (R-1.20)/DRM; RmNUMoDRM2 = RmNUMoDRM*RmNUMoDRM;
    D18 = SC*FS/(1.0+RmNUMoDRM2);

    Result = A1*D1+A2*D2+A3*D3+A4*D4+A5*D5+A6*D6+A7*D7+A8*D8+A9*D9
        + A10*D10+A11*D11+A12*D12+A13*D13+A14*D14+A15*D15+A16*D16+A17*D17
        + A18*D18;

    return( Result );
}



double    T01S_BT_PRC_Q( double R, double SINT, double COST) {


    /*
     *  Calculates the Theta component of the "quadrupole" part of the model partial ring current.
     */

    double A1=12.74640393, A2=-7.516393516;
    double A3=-5.476233865, A4=3.212704645, A5=-59.10926169, A6=46.62198189, A7=-.01644280062;        // ALL LINEAR PARAMETERS HERE
    double A8=.1234229112, A9=-.08579198697, A10=.01321366966, A11=.8970494003, A12=9.136186247;    // WERE MULTIPLIED BY 0.1,
    double A13=-38.19301215, A14=21.73775846, A15=-410.0783424, A16=-69.90832690, A17=-848.8543440;    // SO THAT THEY CORRESPOND TO P_0=1 nPa,
    double XK1=1.243288286, AL1=.2071721360, DAL1=.05030555417, B1=7.471332374, BE1=3.180533613;    // RATHER THAN THE ORIGINAL VALUE OF 10 nPa
    double XK2=1.376743507, AL2=.1568504222, DAL2=.02092910682, BE2=1.985148197, XK3=.3157139940;    // ASSUMED IN THE BIOT-SAVART INTEGRAL
    double XK4=1.056309517, AL3=.1701395257, DAL3=.1019870070, B3=6.293740981, BE3=5.671824276;
    double AL4=.1280772299, DAL4=.02189060799, DG1=.01040696080, AL5=.1648265607, DAL5=.04701592613;
    double DG2=.01526400086, C1=12.88384229, C2=3.361775101, C3=23.44173897;

    double ooR, ooR2, SINT2, COST2, SC, ALPHA, GAMMA, F, FA, FS, D1, D2;
    double D3, D4, D5, D6, ALPHAmAL4oDAL42, ALPHAmAL4oDAL4, FCC, ooFCC, D7;
    double ALPHAmAL5oDAL52, ALPHAmAL5oDAL5, GAMMAoDG2, GAMMAoDG22, ARG, ooARG;
    double D8, D9, D10, D11, D12, D13, D14, R2, R4, C12, C22, C32, D15, Result;
    double D16, D17;


    ooR = 1.0/R; ooR2 = ooR*ooR;
    SINT2 = SINT*SINT;
    COST2 = COST*COST;
    SC    = SINT*COST;
    ALPHA = SINT2*ooR;
    GAMMA = COST*ooR2;

    T01S_FFS(ALPHA,AL1,DAL1,&F,&FA,&FS);
    D1 = pow(F, XK1)/(pow(R/B1, BE1)+1.0);
    D2 = D1*COST2;

    T01S_FFS(ALPHA,AL2,DAL2,&F,&FA,&FS);
    D3 = pow(FA, XK2)/pow(R, BE2);
    D4 = D3*COST2;

    T01S_FFS(ALPHA,AL3,DAL3,&F,&FA,&FS);
    D5 = pow(FS, XK3)*pow(ALPHA, XK4)/(pow(R/B3, BE3)+1.0);
    D6 = D5*COST2;

    T01S_FFS(GAMMA,0.0,DG1,&F,&FA,&FS);
    ALPHAmAL4oDAL4 = (ALPHA-AL4)/DAL4; ALPHAmAL4oDAL42 = ALPHAmAL4oDAL4*ALPHAmAL4oDAL4;
    FCC = (1.0+ALPHAmAL4oDAL42); ooFCC = 1.0/FCC;
    D7  = ooFCC*FS;
    D8  = D7*ooFCC;
    D9  = D8*ooFCC;
    D10 = D9*ooFCC;

    ALPHAmAL5oDAL5 = (ALPHA-AL5)/DAL5; ALPHAmAL5oDAL52 = ALPHAmAL5oDAL5*ALPHAmAL5oDAL5;
    GAMMAoDG2 = GAMMA/DG2; GAMMAoDG22 = GAMMAoDG2*GAMMAoDG2;
    ARG = 1.0+ALPHAmAL5oDAL52; ooARG = 1.0/ARG;
    D11 = ooARG/(1.0+GAMMAoDG22);
    D12 = D11*ooARG;
    D13 = D12*ooARG;
    D14 = D13*ooARG;

    R2 = R*R; R4 = R2*R2;
    C12 = C1*C1;
    C22 = C2*C2;
    C32 = C3*C3;
    D15 = 1.0/(R4+C12);
    D16 = COST2/(R4+C22);
    D17 = COST2*COST2/(R4+C32);

    Result = A1*D1+A2*D2+A3*D3+A4*D4+A5*D5+A6*D6+A7*D7+A8*D8+A9*D9
        + A10*D10+A11*D11+A12*D12+A13*D13+A14*D14+A15*D15+A16*D16+A17*D17;

    return( Result );

}

void    T01S_FFS( double A, double A0, double DA, double *F, double *FA, double *FS ) {

    double    DA2, ApA0, ApA02, AmA0, AmA02, SQ1, SQ2;

    DA2 = DA*DA;
    ApA0 = A+A0; ApA02 = ApA0*ApA0;
    AmA0 = A-A0; AmA02 = AmA0*AmA0;

    SQ1 = sqrt(ApA02+DA2);
    SQ2 = sqrt(AmA02+DA2);
    *FA  = 2.0/(SQ1+SQ2);
    *F   = *FA*A;
    *FS  = 0.50*(SQ1+SQ2)/(SQ1*SQ2)*(1.0-(*F)*(*F));

    return;

}


void    T01S_RC_SHIELD( double *A, double PS, double X_SC, double X, double Y, double Z, double *BX, double *BY, double *BZ, LgmTsyg2001_Info *t ) {


    int    L, M, I, K, N, NN ;
    double X_SCp1, X_SCp12, X_SCp13, FAC_SC, CPS, SPS, S3PS, PST1;
    double PST2, ST1, CT1, ST2, CT2, X1, Z1, X2, Z2, GX, GY, GZ, P;
    double ooP, ooP2, Q, ooQ, ooQ2, YoP, YoQ, CYPI, CYQI, SYPI, SYQI;
    double R, ooR, ooR2, S, ooS, ooS2, Z1oR, Z2oS, SZRK, CZSK, CZRK;
    double SZSK, SQPR, SQQS, EPR, EQS, FX, FY, FZ, HX, HY, HZ, HXR, HZR;


    X_SCp1 = X_SC+1.0; X_SCp12 = X_SCp1*X_SCp1; X_SCp13 = X_SCp1*X_SCp12;

    FAC_SC = X_SCp13;

    //CPS = cos(PS); SPS = sin(PS);
    CPS = t->cos_psi; SPS = t->sin_psi;

    S3PS = 2.0*CPS;

    PST1 = PS*A[85];
    PST2 = PS*A[86];

    ST1 = sin(PST1);
    CT1 = cos(PST1);
    ST2 = sin(PST2);
    CT2 = cos(PST2);

    X1 = X*CT1-Z*ST1;
    Z1 = X*ST1+Z*CT1;
    X2 = X*CT2-Z*ST2;
    Z2 = X*ST2+Z*CT2;

    L  = 0;
    GX = 0.0;
    GY = 0.0;
    GZ = 0.0;

    for (M=1; M<=2; M++ ){      // M=1 IS FOR THE 1ST SUM ("PERP." SYMMETRY)
                                // AND M=2 IS FOR THE SECOND SUM ("PARALL." SYMMETRY)

    for (I=1; I<=3; I++ ){

        P    = A[72+I]; ooP = 1.0/P; ooP2 = ooP*ooP;
        Q    = A[78+I]; ooQ = 1.0/Q; ooQ2 = ooQ*ooQ;
        YoP = Y*ooP; YoQ = Y*ooQ;
        CYPI = cos(YoP);
        CYQI = cos(YoQ);
        SYPI = sin(YoP);
        SYQI = sin(YoQ);

        for (K=1; K<=3; K++ ){

        R    = A[75+K]; ooR = 1.0/R; ooR2 = ooR*ooR;
        S    = A[81+K]; ooS = 1.0/S; ooS2 = ooS*ooS;
        Z1oR = Z1*ooR; Z2oS = Z2*ooS;
        SZRK = sin(Z1oR);
        CZSK = cos(Z2oS);
        CZRK = cos(Z1oR);
        SZSK = sin(Z2oS);
        SQPR = sqrt(ooP2+ooR2);
        SQQS = sqrt(ooQ2+ooS2);
        EPR  = exp(X1*SQPR);
        EQS  = exp(X2*SQQS);

        for (N=1; N<=2; N++ ){      // N=1 IS FOR THE FIRST PART OF EACH COEFFICIENT
                                    // AND N=2 IS FOR THE SECOND ONE

            for (NN=1; NN<=2; NN++ ){   // NN = 1,2 FURTHER SPLITS THE COEFFICIENTS INTO 2 PARTS,
                                        // TO TAKE INTO ACCOUNT THE SCALE FACTOR DEPENDENCE

                if (M == 1) {

                    FX = -SQPR*EPR*CYPI*SZRK  *FAC_SC;
                    FY = EPR*SYPI*SZRK*ooP   *FAC_SC;
                    FZ = -EPR*CYPI*CZRK*ooR  *FAC_SC;

                    if (N == 1) {
                        if (NN == 1) {
                            HX = FX;
                            HY = FY;
                            HZ = FZ;
                        } else {
                            HX = FX*X_SC;
                            HY = FY*X_SC;
                            HZ = FZ*X_SC;
                        }
                    } else {
                        if (NN == 1) {
                            HX = FX*CPS;
                            HY = FY*CPS;
                            HZ = FZ*CPS;
                        } else {
                            HX = FX*CPS*X_SC;
                            HY = FY*CPS*X_SC;
                            HZ = FZ*CPS*X_SC;
                        }
                    }

            } else {        // M.EQ.2

                FX = -SPS*SQQS*EQS*CYQI*CZSK  *FAC_SC;
                FY = SPS*ooQ*EQS*SYQI*CZSK   *FAC_SC;
                FZ = SPS*ooS*EQS*CYQI*SZSK   *FAC_SC;

                if (N == 1) {

                    if (NN == 1) {
                        HX = FX;
                        HY = FY;
                        HZ = FZ;
                    } else {
                        HX = FX*X_SC;
                        HY = FY*X_SC;
                        HZ = FZ*X_SC;
                    }
                } else {

                    if (NN == 1) {
                        HX = FX*S3PS;
                        HY = FY*S3PS;
                        HZ = FZ*S3PS;
                    } else {
                        HX = FX*S3PS*X_SC;
                        HY = FY*S3PS*X_SC;
                        HZ = FZ*S3PS*X_SC;
                    }
                }
            }

            ++L;

            if (M == 1) {
                HXR = HX*CT1+HZ*ST1;
                HZR = -HX*ST1+HZ*CT1;
            } else {
                HXR = HX*CT2+HZ*ST2;
                HZR = -HX*ST2+HZ*CT2;
            }

            GX = GX + HXR*A[L];
            GY = GY + HY *A[L];
            GZ = GZ + HZR*A[L];

            }
        }
        }
    }
    }

    *BX = GX;
    *BY = GY;
    *BZ = GZ;

    return;

}

void    T01S_DIPOLE( double PS, double X, double Y, double Z, double *BX, double *BY, double *BZ, LgmTsyg2001_Info *t ) {

    /*
     *        A DOUBLE PRECISION ROUTINE
     *
     *    CALCULATES GSM COMPONENTS OF A GEOT01S_DIPOLE FIELD WITH THE T01S_DIPOLE MOMENT
     *    CORRESPONDING TO THE EPOCH OF 2000.
     *
     *  ----INPUT PARAMETERS:
     *       PS - GEOT01S_DIPOLE TILT ANGLE IN RADIANS,
     *       X,Y,Z - GSM COORDINATES IN RE (1 RE = 6371.2 km)
     *
     *  ----OUTPUT PARAMETERS:
     *       BX,BY,BZ - FIELD COMPONENTS IN GSM SYSTEM, IN NANOTESLA.
     */
    double    SPS, CPS, P, U, V, T, s, s2, s3, s5, Q;

    //SPS = sin(PS); CPS = cos(PS);
    SPS = t->sin_psi; CPS = t->cos_psi;
    P = X*X;
    U = Z*Z;
    V = 3.0*Z*X;
    T = Y*Y;
    s = sqrt(P+T+U); s2 = s*s; s3 = s2*s; s5 = s2*s3;
    Q = 30115.0/s5;

    *BX = Q*((T+U-2.0*P)*SPS-V*CPS);
    *BY = -3.0*Y*Q*(X*SPS+Z*CPS);
    *BZ = Q*((P+T-2.0*U)*CPS-V*SPS);

    return;

}


/*
 *   $Id$
 */

