/*
 *
 *  The T89 Magnetic field model. See following reference (note: there are some
 *  errors in the paper);
 *
 *      Tsyganenko, N. A., A magnetospheric magnetic field model with a warped
 *      tail current sheet, Planetary and Space Science Volume 37, Issue 1,
 *      January 1989, Pages 5-20.
 *
 *      N. A. Tsyganenko, Correction to A magnetospheric magnetic field model
 *      with a warped tail current sheet. Planet. Space Sci. 37, S-20 (1989).
 *
 *
 *  Note that this version also has 7 Kp levels.
 *
 */

#include "Lgm/Lgm_MagModelInfo.h"

/*
 *  The 30 constant parameters for the 7 T89c Kp models -- there
 *  is one set of 39 for each Kp level. Specifically;
 *
 *                a[0][] -> Kp =     0, 0+
 *                a[1][] -> Kp = 1-, 1, 1+
 *                a[2][] -> Kp = 2-, 2, 2+
 *                a[3][] -> Kp = 3-, 3, 3+
 *                a[4][] -> Kp = 4-, 4, 4+
 *                a[5][] -> Kp = 5-, 5, 5+
 *                a[6][] -> Kp = >= 6-
 *
 *  The organization of these coefficients is different from the T89.c version.
 *
 */
static double Lgm_T89c_a[7][31] = {
    {  -1e31,   -116.53,    -10719.0,   42.375,     59.753,     -11363.0,   1.7844,
                30.268,     -0.035372,  -0.066832,  0.016456,   -1.3024,    0.0016529,
                0.0020293,  20.289,     -0.025203,  224.91,     -9234.8,    22.788,       /* Kp =     0, 0+ */
                7.8813,     1.8362,     -0.27228,   8.8184,     2.8714,     14.468,
                32.177,     0.01,       0.0,        7.0459,     4.0,        20.0 },

    {  -1e31,   -55.553,    -13198.0,   60.647,     61.072,     -16064.0,   2.2534,
                34.407,     -0.038887,  -0.094571,  0.027154,   -1.3901,    0.0013460,
                0.0013238,  23.005,     -0.030565,  55.047,     -3875.7,    20.178,     /* Kp = 1-, 1, 1+ */
                7.9693,     1.4575,     0.89471,    9.4039,     3.5215,     14.474,
                36.555,     0.01,       0.0,        7.0787,     4.0,        20.0},

    {  -1e31,   -101.34,    -13480.0,   111.35,     12.386,     -24699.0,   2.6459,
                38.948,     -0.034080,  -0.12404,   0.029702,   -1.4052,    0.0012103,
                0.0016381,  24.49,      -0.037705,  -298.32,    4400.9,     18.692,     /* Kp = 2-, 2, 2+ */
                7.9064,     1.3047,     2.4541,     9.7012,     7.1624,     14.288,
                33.822,     0.01,       0.0,        6.7442,     4.0,        20.0},

    {  -1e31,   -181.69,    -12320.0,   173.79,     -96.664,    -39051.0,   3.2633,
                44.968,     -0.046377,  -0.16686,   0.048298,   -1.5473,    0.0010277,
                0.0031632,  27.341,     -0.050655,  -514.10,    12482.0,    16.257,     /* Kp = 3-, 3, 3+ */
                8.5834,     1.0194,     3.6148,     8.6042,     5.5057,     13.778,
                32.373,     0.01,       0.0,        7.3195,     4.0,        20.0},

    {  -1e31,   -436.54,    -9001.0,    323.66,     -410.08,    -50340.0,   3.9932,
                58.524,     -0.038519,  -0.26822,   0.074528,   -1.4268,    -0.0010985,
                0.0096613,  27.557,     -0.056522,  -867.03,    20652.0,    14.101,     /* Kp = 4-, 4, 4+ */
                8.3501,     0.72996,    3.8149,     9.2908,     6.4674,     13.729,
                28.353,     0.01,       0.0,        7.4237,     4.0,        20.0},

    {  -1e31,   -707.77,    -4471.9,    432.81,     -435.51,    -60400.0,   4.6229,
                68.178,     -0.088245,  -0.21002,   0.11846,    -2.6711,    0.0022305,
                0.010910,   27.547,     -0.054080,  -424.23,    1100.2,     13.954,     /* Kp = 5-, 5, 5+ */
                7.5337,     0.89714,    3.7813,     8.2945,     5.174,      14.213,
                25.237,     0.01,       0.0,        7.0037,     4.0,        20.0},

    {  -1e31,   -1190.4,    2749.9,     742.56,     -1110.3,    -77193.0,   7.6727,
                102.05,     -0.096015,  -0.74507,   0.11214,    -1.3614,    0.0015157,
                0.022283,   23.164,     -0.074146,  -2219.1,    48253.0,    12.714,     /* Kp >= 6-      */
                7.6777,     0.57138,    2.9633,     9.3909,     9.7263,     11.123,
                21.558,     0.01,       0.0,        4.4518,     4.0,        20.0}
};


void T89c( int IOPT, double *PARMOD, double PS, double X, double Y, double Z, double *BX, double *BY, double *BZ, Lgm_MagModelInfo *Info ) {

    /*
     *
     *
     *   COMPUTES GSM COMPONENTS OF THE MAGNETIC FIELD PRODUCED BY EXTRA-
     *   TERRESTRIAL CURRENT SYSTEMS IN THE GEOMAGNETOSPHERE. THE MODEL IS
     *   VALID UP TO GEOCENTRIC DISTANCES OF 70 RE AND IS BASED ON THE MER-
     *   GED IMP-A,C,D,E,F,G,H,I,J (1966-1974), HEOS-1 AND -2 (1969-1974),
     *   AND ISEE-1 AND -2  SPACECRAFT DATA SET.
     *
     *   THIS IS A MODIFIED VERSION (T89c), WHICH REPLACED THE ORIGINAL ONE
     *     IN 1992 AND DIFFERS FROM IT IN THE FOLLOWING:
     *
     *   (1)  ISEE-1,2 DATA WERE ADDED TO THE ORIGINAL IMP-HEOS DATASET
     *   (2)  TWO TERMS WERE ADDED TO THE ORIGINAL TAIL FIELD MODES, ALLOWING
     *          A MODULATION OF THE CURRENT BY THE GEODIPOLE TILT ANGLE
     *
     *
     *  REFERENCE FOR THE ORIGINAL MODEL: N.A. TSYGANENKO, A MAGNETOSPHERIC MAGNETIC
     *       FIELD MODEL WITH A WARPED TAIL CURRENT SHEET: PLANET.SPACE SCI., V.37,
     *         PP.5-20, 1989.
     *
     *----INPUT PARAMETERS: IOPT - SPECIFIES THE GROUND DISTURBANCE LEVEL:
     *
     *   IOPT= 1       2        3        4        5        6      7
     *                  CORRESPOND TO:
     *    KP= 0,0+  1-,1,1+  2-,2,2+  3-,3,3+  4-,4,4+  5-,5,5+  > =6-
     *
     *    PS - GEODIPOLE TILT ANGLE IN RADIANS
     *    X, Y, Z  - GSM COORDINATES OF THE POINT IN EARTH RADII
     *
     *----OUTPUT PARAMETERS: BX,BY,BZ - GSM COMPONENTS OF THE MODEL MAGNETIC
     *                        FIELD IN NANOTESLAS
     *
     *   THE PARAMETER PARMOD(10) IS A DUMMY ARRAY.  IT IS NOT USED IN THIS
     *        SUBROUTINE AND IS PROVIDED JUST FOR MAKING IT COMPATIBLE WITH THE
     *           NEW VERSION (4/16/96) OF THE GEOPACK SOFTWARE.
     *
     *   THIS RELEASE OF T89C IS DATED  FEB 12, 1996;
     *--------------------------------------------------------------------------
     *
     *
     *              AUTHOR:     NIKOLAI A. TSYGANENKO
     *                          HSTX CORP./NASA GSFC
     *
     *
     *
     *  Converted to C by M. G. Henderson, LANL.
     *
     */
    double  *A;
    double  XI[5], F[4];
    int     indx;


    indx = Info->Kp;
    if (indx < 0) indx = 0;
    if (indx > 6) indx = 6;

    A = Lgm_T89c_a[ indx ];

    XI[1] = X;
    XI[2] = Y;
    XI[3] = Z;
    XI[4] = Info->c->psi;

    T89c_Field( 1, A, XI, F, Info );

    *BX = F[1];
    *BY = F[2];
    *BZ = F[3];

}






void T89c_Field( int ID, double *A, double *XI, double *F, Lgm_MagModelInfo *Info ) {

    /*
     *-------------------------------------------------------------------
     *
     *
     *        ***  N.A. Tsyganenko ***  8-10.12.1991  ***
     *
     *      Calculates dependent model variables and their deriva-
     *  tives for given independent variables and model parame-
     *  ters.  Specifies model functions with free parameters which
     *  must be determined by means of least squares fits (RMS
     *  minimization procedure).
     *
     *      Description of parameters:
     *
     *  ID  - number of the data point in a set (initial assignments are performed
     *        only for ID=1, saving thus CPU time)
     *  A   - input vector containing model parameters;
     *  XI  - input vector containing independent variables;
     *  F   - output double precision vector containing
     *        calculated values of dependent variables;
     *  DER   - output double precision vector containing
     *        calculated values for derivatives of dependent
     *        variables with respect to model parameters;
     *
     * - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     *
     *      T89 represents external magnetospheric magnetic field
     *  in Cartesian SOLAR MAGNETOSPHERIC coordinates (Tsyganenko N.A.,
     *  Planet. Space Sci., 1989, v.37, p.5-20; the "T89 model" with the warped
     *  tail current sheet) + A MODIFICATION ADDED IN APRIL 1992 (SEE BELOW)
     *
     *      Model formulas for the magnetic field components contain in total
     *  30 free parameters (17 linear and 13 nonlinear parameters).
     *      First 2 independent linear parameters A(1)-A(2) correspond to contribu-
     *  tion from the tail current system, then follow A(3) and A(4) which are the
     *  amplitudes of symmetric and antisymmetric terms in the contribution from
     *  the closure currents; A(5) is the ring current amplitude. Then follow the
     * coefficients A(6)-A(15) which define Chapman-Ferraro+Birkeland current field.
     *    The coefficients c16-c19  (see Formula 20 in the original paper),
     *   due to DivB=0 condition, are expressed through A(6)-A(15) and hence are not
     *    independent ones.
     *  A(16) AND A(17) CORRESPOND TO THE TERMS WHICH YIELD THE TILT ANGLE DEPEN-
     *    DENCE OF THE TAIL CURRENT INTENSITY (ADDED ON APRIL 9, 1992)
     *
     *      Nonlinear parameters:
     *
     *    A(18) : DX - Characteristic scale of the Chapman-Ferraro field along the
     *        X-axis
     *    A(19) : ADR (aRC) - Characteristic radius of the ring current
     *    A(20) : D0 - Basic half-thickness of the tail current sheet
     *    A(21) : DD (GamRC)- defines rate of thickening of the ring current, as
     *             we go from night- to dayside
     *    A(22) : Rc - an analog of "hinging distance" entering formula (11)
     *    A(23) : G - amplitude of tail current warping in the Y-direction
     *    A(24) : aT - Characteristic radius of the tail current
     *    A(25) : Dy - characteristic scale distance in the Y direction entering
     *                 in W(x,y) in (13)
     *    A(26) : Delta - defines the rate of thickening of the tail current sheet
     *                 in the Y-direction (in T89 it was fixed at 0.01)
     *    A(27) : Q - this parameter was fixed at 0 in the final version of T89;
     *              initially it was introduced for making Dy to depend on X
     *    A(28) : Sx (Xo) - enters in W(x,y) ; see (13)
     *    A(29) : Gam (GamT) - enters in DT in (13) and defines rate of tail sheet
     *              thickening on going from night to dayside; in T89 fixed at 4.0
     *    A(30) : Dyc - the Dy parameter for closure current system; in T89 fixed
     *               at 20.0
     *  - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     *
     */

    double  DER[4][31];
    double  A02, XLW2, YN, RPI, RT, XD, XLD2, SXC, XLWC2, DXL;
    double  DYC, DYC2, DX, HA02, RDX2M, RDX2, RDYC2, HLWC2M, DRDYC2;
    double  DRDYC3, HXLW2M, ADR, D0, DD, RC, G, AT, DT, DEL, P, Q, SX, GAM;
    double  HXLD2M, ADSL, XGHS, H, HS, GAMH, W1, DBLDEL, W2, W3, W4, W5, W6, AK1, AK2, AK3;
    double  AK4, AK5, AK6, AK7, AK8, AK9, AK10, AK11, AK12, AK13, AK14, AK15, AK16, AK17, SXA;
    double  SYA, SZA, AK610, AK711, AK812, AK913, RDXL, HRDXL, A6H, A9T, YNP, YND, X, Y, Z;
    double  TILT, TLT2, SPS, CPS, X2, Y2, Z2, TPS, HTP, GSP, XSM, ZSM, XRC, XRC16, SXRC;
    double  Y4, Y410, GSY4, SY4, ZS1, DZSX, ZS, D2ZSGY, DZSY, XSM2, DSQT, DSQT3, FA0, DDR;
    double  DFA0, ZR, TR, RTR, RO2, ADRT, ADRT2, FK, DSFC, FC, FACXY, XZR, YZR, DBXDP;
    double  XZYZ, FAQ, DBZDP, DELY2, D, XXD, RQD, RQDS, D2, T, XSMX, RDSQ2, RDSQ, V, DVX;
    double  OM, OMS, RDY, OMSV, RDY2, FY, W, YFY1, FYPR, FYDY, DWX, YDWY, DDY, ATT, S1, F5;
    double  F7, F1, F3, F9, FS, XDWX, RTT, WT, BRRZ1, BRRZ2, DBXC1, DBXC2, WTFS, DBZC1;
    double  DBZC2, ZPL, ZMN, ROGSM2, SPL, SMN, XSXC, RQC2, RQC, FYC, WC, DWCX, DWCY, SZRP;
    double  SZRM, XYWC, WCSP, WCSM, FXYP, FXYM, FXPL, FXMN, FYPL, FYMN, FZPL, FZMN, EX, EC;
    double  ES, ECZ, ESZ, ESZY2, ESZZ2, ECZ2, ESY, SX1, SY1, SZ1, BXCL, BYCL, BZCL, BXT;
    double  BYT, BZT;
    int     i, l;



    /*
     * The last four quantities define variation of tail sheet thickness along X
     */
    A02  = 25.0;
    XLW2 = 170.0;
    YN   = 30.0;
    RPI  = 0.31830989;
    RT   = 30.0;
    XD   = 0.0;
    XLD2 = 40.0;



    /*
     * The two quantities belong to the function WC which confines tail closure
     * current in X- and Y- direction
     */
    SXC = 4.0;
    XLWC2 = 50.0;


    DXL = 20.0;


//    if ( ID == 1 ) {

        for (i=1; i<=30; i++ ) {
            for (l=1; l<=3; l++ ) {
                DER[l][i] = 0.0;
            }
        }


        DYC     = A[30];
        DYC2    = DYC*DYC;
        DX      = A[18];
        HA02    = 0.5*A02;
        RDX2M   = -1.0/(DX*DX);
        RDX2    = -RDX2M;
        RDYC2   = 1.0/DYC2;
        HLWC2M  = -0.50*XLWC2;
        DRDYC2  = -2.0*RDYC2;
        DRDYC3  = 2.0*RDYC2*sqrt(RDYC2);
        HXLW2M  = -0.50*XLW2;

        ADR  = A[19];
        D0   = A[20];
        DD   = A[21];
        RC   = A[22];
        G    = A[23];
        AT   = A[24];
        DT   = D0;
        DEL  = A[26];
        P    = A[25];
        Q    = A[27];
        SX   = A[28];
        GAM  = A[29];

        HXLD2M  = -0.5*XLD2;
        ADSL    = 0.0;
        XGHS    = 0.0;
        H       = 0.0;
        HS      = 0.0;
        GAMH    = 0.0;
        W1      = -0.5/DX;
        DBLDEL  = 2.0*DEL;
        W2      = W1*2.0;
        W4      = -1.0/3.0;
        W3      = W4/DX;
        W5      = -0.5;
        W6      = -3.0;

        AK1  = A[1];
        AK2  = A[2];
        AK3  = A[3];
        AK4  = A[4];
        AK5  = A[5];
        AK6  = A[6];
        AK7  = A[7];
        AK8  = A[8];
        AK9  = A[9];
        AK10 = A[10];
        AK11 = A[11];
        AK12 = A[12];
        AK13 = A[13];
        AK14 = A[14];
        AK15 = A[15];
        AK16 = A[16];
        AK17 = A[17];

        SXA = 0.0;
        SYA = 0.0;
        SZA = 0.0;

        AK610 = AK6*W1 + AK10*W5;
        AK711 = AK7*W2 - AK11;
        AK812 = AK8*W2 + AK12*W6;
        AK913 = AK9*W3 + AK13*W4;
        RDXL  = 1.0/DXL;
        HRDXL = 0.5*RDXL;
        A6H = AK6*0.5;
        A9T = AK9/3.0;
        YNP = RPI/YN*0.50;
        YND = 2.0*YN;

//    }


    X    = XI[1];
    Y    = XI[2];
    Z    = XI[3];
    //TILT = XI[4];
    TILT = Info->c->psi;
    TLT2 = TILT*TILT;
    //SPS  = sin(TILT);
    //CPS  = sqrt( 1.0 - SPS*SPS );
    SPS  = Info->c->sin_psi;
    CPS  = Info->c->cos_psi;


    X2 = X*X;
    Y2 = Y*Y;
    Z2 = Z*Z;

    //TPS = SPS/CPS;
    TPS = Info->c->tan_psi;
    HTP = TPS*0.5;
    GSP = G*SPS;
    XSM = X*CPS - Z*SPS;
    ZSM = X*SPS + Z*CPS;


    /*
     *  CALCULATE THE FUNCTION ZS DEFINING THE SHAPE OF THE TAIL CURRENT SHEET
     *  AND ITS SPATIAL DERIVATIVES:
     */
       XRC    = XSM + RC;
       XRC16  = XRC*XRC + 16.0;
       SXRC   = sqrt( XRC16 );
       Y4     = Y2*Y2;
       Y410   = Y4 + 1.0e4;
       SY4    = SPS/Y410;
       GSY4   = G*SY4;
       ZS1    = HTP*(XRC-SXRC);
       DZSX   = -ZS1/SXRC;
       ZS     = ZS1 - GSY4*Y4;
       D2ZSGY = -SY4/Y410*4.0e4*Y2*Y;
       DZSY   = G*D2ZSGY;



    /*
     *  CALCULATE THE COMPONENTS OF THE RING CURRENT CONTRIBUTION:
     */
    XSM2  = XSM*XSM;
    DSQT  = sqrt( XSM2 + A02 );
    DSQT3 = DSQT*DSQT*DSQT;
    FA0   = 0.5*(1.0 + XSM/DSQT);
    DDR   = D0 + DD*FA0;
    DFA0  = HA02/DSQT3;
    ZR    = ZSM - ZS;
    TR    = sqrt( ZR*ZR + DDR*DDR );
    RTR   = 1.0/TR;
    RO2   = XSM2 + Y2;
    ADRT  = ADR + TR;
    ADRT2 = ADRT*ADRT;
    FK    = 1.0/(ADRT2 + RO2);
    DSFC  = sqrt( FK );
    FC    = FK*FK*DSFC;
    FACXY = 3.0*ADRT*FC*RTR;
    XZR   = XSM*ZR;
    YZR   = Y*ZR;
    DBXDP = FACXY*XZR;

    DER[2][5] = FACXY*YZR;
    XZYZ     = XSM*DZSX + Y*DZSY;
    FAQ      = ZR*XZYZ - DDR*DD*DFA0*XSM;
    DBZDP    = FC*(2.0*ADRT2 - RO2) + FACXY*FAQ;
    DER[1][5] = DBXDP*CPS + DBZDP*SPS;
    DER[3][5] = DBZDP*CPS - DBXDP*SPS;


    /*
     *  CALCULATE THE TAIL CURRENT SHEET CONTRIBUTION:
     */
    DELY2 = DEL*Y2;
    D     = DT + DELY2;
    if ( fabs(GAM) >= 1e-6 ) {

        XXD  =  XSM-XD;
        RQD  =  1.0/(XXD*XXD + XLD2);
        RQDS =  sqrt( RQD );
        H    =  0.5*( 1.0 +XXD*RQDS );
        HS   =  -HXLD2M*RQD*RQDS;
        GAMH =  GAM*H;
        D    += GAMH;
        XGHS =  XSM*GAM*HS;
        ADSL = -D*XGHS;

    }

    D2 = D*D;

    T     = sqrt( ZR*ZR + D2 );
    XSMX  = XSM-SX;
    RDSQ2 = 1.0/( XSMX*XSMX + XLW2 );
    RDSQ  = sqrt(RDSQ2);
    V     = 0.5*(1.0 - XSMX*RDSQ);
    DVX   = HXLW2M*RDSQ*RDSQ2;
    OM    = sqrt( sqrt( XSM2 + 16.0 ) - XSM );
    OMS   = -OM/(OM*OM + XSM)*0.5;
    RDY   = 1.0/(P + Q*OM);
    OMSV  = OMS*V;
    RDY2  = RDY*RDY;
    FY    = 1.0/(1.0 + Y2*RDY2);
    W     = V*FY;

    YFY1 = 2.0*FY*Y2*RDY2;
    FYPR = YFY1*RDY;
    FYDY = FYPR*FY;
    DWX  = DVX*FY + FYDY*Q*OMSV;
    YDWY = -V*YFY1*FY;
    DDY  = DBLDEL*Y;
    ATT  = AT + T;
    S1   = sqrt( ATT*ATT + RO2 );

    F5   = 1.0/S1;
    F7   = 1.0/(S1 + ATT);
    F1   = F5*F7;
    F3   = F5*F5*F5;
    F9   = ATT*F3;
    FS   = ZR*XZYZ - D*Y*DDY + ADSL;
    XDWX = XSM*DWX + YDWY;
    RTT  = 1.0/T;
    WT   = W*RTT;

    BRRZ1 = WT*F1;
    BRRZ2 = WT*F3;
    DBXC1 = BRRZ1*XZR;
    DBXC2 = BRRZ2*XZR;

    DER[2][1]  = BRRZ1*YZR;
    DER[2][2]  = BRRZ2*YZR;
    DER[2][16] = DER[2][1]*TLT2;
    DER[2][17] = DER[2][2]*TLT2;

    WTFS  = WT*FS;
    DBZC1 = W*F5 + XDWX*F7 + WTFS*F1;
    DBZC2 = W*F9 + XDWX*F1 + WTFS*F3;

    DER[1][1]  = DBXC1*CPS + DBZC1*SPS;
    DER[1][2]  = DBXC2*CPS + DBZC2*SPS;
    DER[3][1]  = DBZC1*CPS - DBXC1*SPS;
    DER[3][2]  = DBZC2*CPS - DBXC2*SPS;
    DER[1][16] = DER[1][1]*TLT2;
    DER[1][17] = DER[1][2]*TLT2;
    DER[3][16] = DER[3][1]*TLT2;
    DER[3][17] = DER[3][2]*TLT2;


    /*
     *  CALCULATE CONTRIBUTION FROM THE CLOSURE CURRENTS
     */
    ZPL    = Z + RT;
    ZMN    = Z - RT;
    ROGSM2 = X2 + Y2;
    SPL    = sqrt( ZPL*ZPL + ROGSM2 );
    SMN    = sqrt( ZMN*ZMN + ROGSM2 );

    XSXC = X-SXC;
    RQC2 = 1.0/(XSXC*XSXC + XLWC2);
    RQC  = sqrt( RQC2 );
    FYC  = 1.0/(1.0 +Y2*RDYC2);
    WC   = 0.5*(1.0 - XSXC*RQC)*FYC;

    DWCX = HLWC2M*RQC2*RQC*FYC;
    DWCY = DRDYC2*WC*FYC*Y;
    SZRP = 1.0/(SPL+ZPL);
    SZRM = 1.0/(SMN-ZMN);
    XYWC = X*DWCX + Y*DWCY ;

    WCSP = WC/SPL;
    WCSM = WC/SMN;
    FXYP = WCSP*SZRP;
    FXYM = WCSM*SZRM;
    FXPL = X*FXYP;
    FXMN = -X*FXYM;
    FYPL = Y*FXYP;
    FYMN = -Y*FXYM;
    FZPL = WCSP + XYWC*SZRP;
    FZMN = WCSM + XYWC*SZRM;

    DER[1][3] =  FXPL + FXMN;
    DER[1][4] = (FXPL - FXMN)*SPS;
    DER[2][3] =  FYPL + FYMN;
    DER[2][4] = (FYPL - FYMN)*SPS;
    DER[3][3] =  FZPL + FZMN;
    DER[3][4] = (FZPL - FZMN)*SPS;


    /*
     *  NOW CALCULATE CONTRIBUTION FROM CHAPMAN-FERRARO SOURCES + ALL OTHER
     */
    EX    = exp(X/DX);
    EC    = EX*CPS;
    ES    = EX*SPS;
    ECZ   = EC*Z;
    ESZ   = ES*Z;
    ESZY2 = ESZ*Y2;
    ESZZ2 = ESZ*Z2;
    ECZ2  = ECZ*Z;
    ESY   = ES*Y;


    DER[1][6]  = ECZ;
    DER[1][7]  = ES;
    DER[1][8]  = ESY*Y;
    DER[1][9]  = ESZ*Z;
    DER[2][10] = ECZ*Y;
    DER[2][11] = ESY;
    DER[2][12] = ESY*Y2;
    DER[2][13] = ESY*Z2;
    DER[3][14] = EC;
    DER[3][15] = EC*Y2;
    DER[3][6]  = ECZ2*W1;
    DER[3][10] = ECZ2*W5;
    DER[3][7]  = ESZ*W2;
    DER[3][11] = -ESZ;
    DER[3][8]  = ESZY2*W2;
    DER[3][12] = ESZY2*W6;
    DER[3][9]  = ESZZ2*W3;
    DER[3][13] = ESZZ2*W4;

    /*
     *  FINALLY, CALCULATE NET EXTERNAL MAGNETIC FIELD COMPONENTS,
     *  BUT FIRST OF ALL THOSE FOR C.-F. FIELD:
     */
    SX1  = AK6*DER[1][6]   + AK7*DER[1][7]   + AK8*DER[1][8]   + AK9*DER[1][9];
    SY1  = AK10*DER[2][10] + AK11*DER[2][11] + AK12*DER[2][12] + AK13*DER[2][13];
    SZ1  = AK14*DER[3][14] + AK15*DER[3][15] + AK610*ECZ2 + AK711*ESZ + AK812*ESZY2 + AK913*ESZZ2;
    BXCL = AK3*DER[1][3]   + AK4*DER[1][4];
    BYCL = AK3*DER[2][3]   + AK4*DER[2][4];
    BZCL = AK3*DER[3][3]   + AK4*DER[3][4];
    BXT  = AK1*DER[1][1]   + AK2*DER[1][2] + BXCL + AK16*DER[1][16] + AK17*DER[1][17];
    BYT  = AK1*DER[2][1]   + AK2*DER[2][2] + BYCL + AK16*DER[2][16] + AK17*DER[2][17];
    BZT  = AK1*DER[3][1]   + AK2*DER[3][2] + BZCL + AK16*DER[3][16] + AK17*DER[3][17];
    F[1] = BXT + AK5*DER[1][5] + SX1 + SXA;
    F[2] = BYT + AK5*DER[2][5] + SY1 + SYA;
    F[3] = BZT + AK5*DER[3][5] + SZ1 + SZA;


}

int Lgm_B_T89c( Lgm_Vector *v, Lgm_Vector *B, Lgm_MagModelInfo *Info ) {

    Lgm_Vector      B1, B5;
    double          PARMOD[11];


    T89c( Info->Kp, PARMOD, Info->c->psi, v->x, v->y, v->z, &B1.x, &B1.y, &B1.z, Info );
    switch ( Info->InternalModel ){

        case LGM_CDIP:
                        Lgm_B_cdip( v, &B5, Info );
                        break;
        case LGM_EDIP:
                        Lgm_B_edip( v, &B5, Info );
                        break;
        case LGM_IGRF:
                        Lgm_B_igrf( v, &B5, Info );
                        break;
        default:
                        fprintf(stderr, "Lgm_B_T89c: Unknown internal model (%d)\n", Info->InternalModel );
                        break;

    }





    B->x = B1.x + B5.x;
    B->y = B1.y + B5.y;
    B->z = B1.z + B5.z;

    ++Info->nFunc;

    return(1);

}
