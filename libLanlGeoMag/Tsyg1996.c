#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Lgm/Lgm_Tsyg1996.h"


/*
 *  This is the T96 model converted to C by M Henderson.
 *  For problems with the C version contact mghenderson@lanl.gov
 *
 *
 *
 *
 *
 * The FORTRAN source code T96_01.FOR is the last version (first release June 22,
 * 1996, amended in July 2010 for compatibility with Intel Fortran compilers) of
 * a data-based model of the geomagnetospheric magnetic field with an explicitly
 * defined realistic magnetopause, large-scale Region 1 and 2 Birkeland current
 * systems, and the IMF penetration  across the boundary.
 * 
 * The file T96_01.FOR contains a set of 33 subroutines and functions. The first
 * subroutine, named T96_01, is the primary one, accepting the input values of the
 * solar wind pressure, Dst-index, By- and Bz-components of the interplanetary
 * magnetic field, the geodipole tilt angle, and GSM position of the observation
 * point (X,Y,Z). The subroutine returns GSM components of the external field
 * (i.e., the total vector of B minus the Earth's contribution).  The remaining
 * 32 subroutines are invoked by T96_01.
 * 
 * 
 * 
 * The source code T96_01.FOR differs in several aspects from the previous version
 * T95_06.FOR (released in November 1995).
 * 
 * (1) The T96_01 does not include the AE-index as an input parameter.
 * 
 * (2) The number of terms in the ring current, tail modes, and in the shielding
 * fields was reduced to a necessary minimum, in order to increase the speed of the
 * code. The tail field is now a sum of two modes; the first one is responsible
 * for the near-Earth tail field and is similar to that in the previous version.
 * The second mode provides the asymptotic field in the far magnetotail.
 * 
 * (3) The way of representing the effects of the dipole tilt on the tail/ring
 * current field was revised:  instead of a "shear" transformation (introduced in
 * the T89 model), a radially dependent "space-warping" is used in this model,
 * which decreases tilt-induced spurious currents.
 * 
 * (4) The representation for the Region 2 Birkeland current field was completely
 * revised:  in the present version, a smooth approximation was developed for
 * the field inside the current layer. As a result, unphysical kinks in the Bz
 * profile on the nightside were eliminated.
 * 
 * 
 * 
 *             *******************************************
 *             | Users should be aware of the following. |
 *             *******************************************
 * 
 * 
 *   (1) A simple linear dependence of the amplitudes of the field sources on the
 * SQRT(Pdyn), Dst, and the IMF-related parameter EPS=SQRT(N)*V*Bt*sin(theta/2)
 * was employed.  Hence, the best results should be expected near the most probable
 * values of the input parameters,  corresponding to the regions in the Pdyn-Dst-
 * ByIMF-BzIMF space with the highest density of the spacecraft measurements. For
 * the same reason, caution is needed in modeling situations with unusually low or
 * high values of these parameters: extrapolating the model too far beyond the
 * range of reliable approximation can lead to unrealistic results.  As a rough
 * estimate, the parameter values should remain within the intervals:
 * Pdyn:  between 0.5 and 10 nPa,
 * Dst:  between -100 and +20,
 * ByIMF and BzIMF: between -10 and +10 nT.
 * 
 *   (2) The only parameter which controls the size of the model magnetopause is
 * the solar wind ram pressure Pdyn. No IMF dependence has been introduced so far
 * in the magnetopause shape/size.  This is planned to be done in future versions
 * of the model.
 *      To facilitate using the model, we provide users with two supplementary
 * FORTRAN subroutines, named LOCATE and CROSSING.  The first one finds the point
 * on the model magnetopause which is closest to a given point of space, for any
 * values of the solar wind density and velocity (or, optionally, solar wind
 * pressure).  The second subroutine estimates the current value of the solar wind
 * ram pressure, based on the observed GSM position of the magnetopause at any
 * local time or latitude, sunward from XGSM=-60Re.
 * 
 *   (3) In its present form, the subroutine T96_01 is compatible with new version
 * (April 16, 1996) of the software package GEOPACK for coordinate transformation
 * and line tracing, which replaces the old version and is available from the same
 * WWW site.
 * 
 *   (4) This is not a "final version":  the model is supposed to be further
 * improved and upgraded in the future. In this regard, any kind of feedback from
 * the users is very important for us and will be greatly appreciated. In
 * particular, it is very important that any problems encountered in adapting and
 * using the code be reported to us as soon as possible.  Please send your
 * questions and comments to the address given in the end of this file.
 * 
 *   (5) More details on the approach used in devising this model can be found in
 * the following publications:
 * 
 * 
 *       Tsyganenko, N.A. and M. Peredo, Analytical models of the magnetic field
 *         of disk-shaped current sheets,  J.Geophys.Res., v.99, pp.199-205, 1994.
 * 
 *       Tsyganenko, N.A., Modeling the Earth's magnetospheric magnetic field
 *         confined within a realistic magnetopause, J.Geophys.Res., v.100,
 *         pp.5599-5612, 1995.
 * 
 *       Fairfield, D.H., N.A. Tsyganenko, A.V. Usmanov, and M.V. Malkov, A large
 *         magnetosphere magnetic field database, J.Geophys.Res., v.99,
 *         pp.11319-11326, 1994.
 * 
 *       Tsyganenko, N.A. and D.P. Stern, Modeling the global magnetic field
 *         the large-scale Birkeland current systems, J.Geophys.Res., v.101,
 *         p.27187-27198, 1996.
 * 
 *       Tsyganenko, N.A., Effects of the solar wind conditions on the global
 *          magnetospheric configuration as deduced from data-based field
 *          models, in:  Proc.of 3rd International Conference on Substorms
 *          (ICS-3), Versailles, France, 12-17 May 1996, ESA SP-389, p.181-185,
 *          1996.
 * 
 *        (A PostScript file of the last paper, named versail.ps, can be ftp-ed
 *        from anonymous ftp-area at:   www-spof.gsfc.nasa.gov;   /pub/kolya)
 * 
 * 
 *    Please send your questions, comments, and requests to:
 * 
 *    Nikolai Tsyganenko
 * 
 *    email:   nikolai.tsyganenko@gmail.com
 */

void Lgm_Init_T96( LgmTsyg1996_Info *t ){

    int                 i, j;

    // Init some params
    t->OLD_PS = -9e99;
    t->OLD_X  = -9e99;
    t->OLD_Y  = -9e99;
    t->OLD_Z  = -9e99;

    t->INTERCON_M_FLAG = 0;


    return;

}


void Tsyg_T96( int IOPT, double *PARMOD, double PS, double SINPS, double COSPS, double X, double Y, double Z, double *BX, double *BY, double *BZ, LgmTsyg1996_Info *t) {

    /*
     *
     *     RELEASE DATE OF THIS (FORTRAN) VERSION:   JUNE 22, 1996.
     *     LAST UPDATE: MAY 01, 2006:  IN THE S/R DIPOLE_T96, SPS AND CPS WERE ADDED IN THE SAVE STATEMENT
     *     ----------------------------------------------------------------------
     *     
     *       WITH TWO CORRECTIONS, SUGGESTED BY T.SOTIRELIS' COMMENTS (APR.7, 1997)
     *     
     *       (1) A "STRAY "  CLOSING PARENTHESIS WAS REMOVED IN THE S/R   R2_BIRK_T96
     *       (2) A 0/0 PROBLEM ON THE Z-AXIS WAS SIDESTEPPED (LINES 44-46 OF THE
     *            DOUBLE PRECISION FUNCTION XKSI_T96)
     *     --------------------------------------------------------------------
     *      DATA-BASED MODEL CALIBRATED BY (1) SOLAR WIND PRESSURE PDYN (NANOPASCALS),
     *                (2) DST (NANOTESLA),  (3) BYIMF, AND (4) BZIMF (NANOTESLA).
     *      THESE INPUT PARAMETERS SHOULD BE PLACED IN THE FIRST 4 ELEMENTS
     *      OF THE ARRAY PARMOD(10).
     *     
     *        THE REST OF THE INPUT VARIABLES ARE: THE GEODIPOLE_T96 TILT ANGLE PS (RADIANS),
     *      AND   X,Y,Z -  GSM POSITION (RE)
     *     
     *        IOPT  IS JUST A DUMMY INPUT PARAMETER, NECESSARY TO MAKE THIS SUBROUTINE
     *      COMPATIBLE WITH THE NEW RELEASE (APRIL 1996) OF THE TRACING SOFTWARE
     *      PACKAGE (GEOPACK). IOPT VALUE DOES NOT AFFECT THE OUTPUT FIELD.
     *     
     *     
     *      OUTPUT:  GSM COMPONENTS OF THE EXTERNAL MAGNETIC FIELD (BX,BY,BZ, nanotesla)
     *                 COMPUTED AS A SUM OF CONTRIBUTIONS FROM PRINCIPAL FIELD SOURCES
     *     
     *       (C) Copr. 1995, 1996, Nikolai A. Tsyganenko, Hughes STX, Code 695, NASA GSFC
     *           Greenbelt, MD 20771, USA
     *     
     *                                 REFERENCES:
     *     
     *                    (1) N.A. TSYGANENKO AND D.P. STERN, A NEW-GENERATION GLOBAL
     *                MAGNETOSPHERE FIELD MODEL  , BASED ON SPACECRAFT MAGNETOMETER DATA,
     *                ISTP NEWSLETTER, V.6, NO.1, P.21, FEB.1996.
     *     
     *                   (2) N.A.TSYGANENKO,  MODELING THE EARTH'S MAGNETOSPHERIC
     *                MAGNETIC FIELD CONFINED WITHIN A REALISTIC MAGNETOPAUSE,
     *                J.GEOPHYS.RES., V.100, P. 5599, 1995.
     *     
     *                   (3) N.A. TSYGANENKO AND M.PEREDO, ANALYTICAL MODELS OF THE
     *                MAGNETIC FIELD OF DISK-SHAPED CURRENT SHEETS, J.GEOPHYS.RES.,
     *                V.99, P. 199, 1994.
     *     
     */

    double  CFX, CFY, CFZ, BXRC, BYRC, BZRC, BXT2, BYT2, BZT2, BXT3, BYT3, BZT3;
    double  R1X, R1Y, R1Z, R2X, R2Y, R2Z, RIMFX, RIMFYS, RIMFZS, QX, QY, QZ;

    double  A[]    = { -9e99, 1.162, 22.344, 18.50, 2.602, 6.903, 5.287, 0.5790, 0.4462, 0.7850 };
    double PDYN0   = 2.0;
    double EPS10   = 3630.7;
    double AM0     = 70.0;
    double S0      = 1.08;
    double X00     = 5.48;
    double DSIG    = 0.005;
    double DELIMFX = 20.0;
    double DELIMFY = 10.0;

    double  PDYN, DST, BYIMF, BZIMF, PPS, SqrtPDYN, DEPR, Bt, THETA, ST, CT, EPS, FACTEPS;
    double  FACTPD, RCAMPL, TAMPL2, TAMPL3, B1AMPL, B2AMPL, RECONN, XAPPA, XAPPA3, YS, ZS;
    double  g, g2, FACTIMF, OIMFX, OIMFY, OIMFZ, RIMFAMPL, XX, YY, ZZ, X0, AM, RHO2, ASQ;
    double  XMXM, AXX0, ARO, SIGMA, SPS, RIMFY, RIMFZ, FX, FY, FZ, FINT, FEXT;

 
 
    PDYN  = PARMOD[1];
    DST   = PARMOD[2];
    BYIMF = PARMOD[3];
    BZIMF = PARMOD[4];
 
    PPS = PS;
//    SPS = SINPS;
    t->sin_psi = SINPS;
 
    SqrtPDYN = sqrt( PDYN );
    DEPR = 0.8*DST - 13.0*SqrtPDYN;  // DEPR is an estimate of total near-Earth depression, based on DST and Pdyn (usually, DEPR < 0 )
 

    /*
     * CALCULATE THE IMF-RELATED QUANTITIES:
     */
    Bt = sqrt( BYIMF*BYIMF + BZIMF*BZIMF );

    if ( (BYIMF == 0.0) && (BZIMF == 0.0) ) {
        THETA = 0.0;
    } else {
        THETA = atan2( BYIMF, BZIMF );
        if ( THETA <= 0.0) THETA += 6.2831853; // MGH - precision for 2pi is pretty low....
    }

    ST  = sin( THETA );
    CT  = cos( THETA );
    EPS = 718.5*SqrtPDYN*Bt*sin( 0.5*THETA );
 
    FACTEPS = EPS/EPS10 - 1.0;
    FACTPD  = sqrt( PDYN/PDYN0 ) - 1.0;
 
    RCAMPL = -A[1]*DEPR; //   RCAMPL is the amplitude of the ring current (positive and equal to abs.value of RC depression at origin)
 
    TAMPL2 = A[2] + A[3]*FACTPD + A[4]*FACTEPS;
    TAMPL3 = A[5] + A[6]*FACTPD;
    B1AMPL = A[7] + A[8]*FACTEPS;
    B2AMPL = 20.0*B1AMPL;  // IT IS EQUIVALENT TO ASSUMING THAT THE TOTAL CURRENT IN THE REGION 2 SYSTEM IS 40% OF THAT IN REGION 1
    RECONN = A[9];
 
    XAPPA  = pow( PDYN/PDYN0, 0.14 );
    XAPPA3 = XAPPA*XAPPA*XAPPA;
    YS = Y*CT - Z*ST;
    ZS = Z*CT + Y*ST;
     
    g = YS/DELIMFY; g2 = g*g;
    FACTIMF = exp( X/DELIMFX - g2 );



    /*
     * CALCULATE THE "IMF" COMPONENTS OUTSIDE THE LAYER  (HENCE BEGIN WITH "O")
     */
    g = RECONN*FACTIMF;
    OIMFX = 0.0;
    OIMFY = g*BYIMF;
    OIMFZ = g*BZIMF;
 
    RIMFAMPL = RECONN*Bt;
 
    PPS = PS;
    XX  = X*XAPPA;
    YY  = Y*XAPPA;
    ZZ  = Z*XAPPA;



    /*
     *
     *  SCALE AND CALCULATE THE MAGNETOPAUSE PARAMETERS FOR THE INTERPOLATION ACROSS
     *   THE BOUNDARY LAYER (THE COORDINATES XX,YY,ZZ  ARE ALREADY SCALED)
     */
    X0   = X00/XAPPA;
    AM   = AM0/XAPPA;
    RHO2 = Y*Y + Z*Z;
    ASQ  = AM*AM;
    XMXM = AM + X - X0;

    if (XMXM < 0.0) XMXM = 0.0; // THE BOUNDARY IS A CYLINDER TAILWARD OF X=X0-AM
    AXX0  = XMXM*XMXM;
    ARO   = ASQ + RHO2;

    g = ARO + AXX0; g2 = g*g;
    SIGMA = sqrt( ( g + sqrt( g2 - 4.0*ASQ*AXX0) )/(2.0*ASQ) );





    /*
     *
     *   NOW, THERE ARE THREE POSSIBLE CASES:
     *    (1) INSIDE THE MAGNETOSPHERE
     *    (2) IN THE BOUNDARY LAYER
     *    (3) OUTSIDE THE MAGNETOSPHERE AND B.LAYER
     *       FIRST OF ALL, CONSIDER THE CASES (1) AND (2):
     */
    if ( SIGMA < (S0+DSIG) ) { // CALCULATE THE T95_06 FIELD (WITH THE POTENTIAL "PENETRATED" INTERCON_T96NECTION FIELD):


        DIPSHLD_T96( PPS, XX, YY, ZZ, &CFX, &CFY, &CFZ, t );
        TAILRC96_T96( SPS, XX, YY, ZZ, &BXRC, &BYRC, &BZRC, &BXT2, &BYT2, &BZT2, &BXT3, &BYT3, &BZT3, t );
        BIRK1TOT_02_T96( PPS, XX, YY, ZZ, &R1X, &R1Y, &R1Z, t );
        BIRK2TOT_02_T96( PPS, XX, YY, ZZ, &R2X, &R2Y, &R2Z, t );
        INTERCON_T96( XX, YS*XAPPA, ZS*XAPPA, &RIMFX, &RIMFYS, &RIMFZS, t );

        RIMFY = RIMFYS*CT + RIMFZS*ST;
        RIMFZ = RIMFZS*CT - RIMFYS*ST;

        FX = CFX*XAPPA3 + RCAMPL*BXRC + TAMPL2*BXT2 + TAMPL3*BXT3 + B1AMPL*R1X + B2AMPL*R2X + RIMFAMPL*RIMFX;
        FY = CFY*XAPPA3 + RCAMPL*BYRC + TAMPL2*BYT2 + TAMPL3*BYT3 + B1AMPL*R1Y + B2AMPL*R2Y + RIMFAMPL*RIMFY;
        FZ = CFZ*XAPPA3 + RCAMPL*BZRC + TAMPL2*BZT2 + TAMPL3*BZT3 + B1AMPL*R1Z + B2AMPL*R2Z + RIMFAMPL*RIMFZ;


        /*
         *  NOW, LET US CHECK WHETHER WE HAVE THE CASE (1). IF YES - WE ARE DONE:
         */
        if ( SIGMA < (S0-DSIG) ) {
             *BX = FX;
             *BY = FY;
             *BZ = FZ;
        } else { //  THIS IS THE MOST COMPLEX CASE: WE ARE INSIDE THE INTERPOLATION REGION
            g = (SIGMA-S0)/DSIG;
            FINT = 0.5*(1.0-g);
            FEXT = 0.5*(1.0+g);
                 
            DIPOLE_T96( PS, X, Y, Z, &QX, &QY, &QZ, t );
            *BX = (FX+QX)*FINT + OIMFX*FEXT - QX;
            *BY = (FY+QY)*FINT + OIMFY*FEXT - QY;
            *BZ = (FZ+QZ)*FINT + OIMFZ*FEXT - QZ;
        } //   THE CASES (1) AND (2) ARE EXHAUSTED; THE ONLY REMAINING POSSIBILITY IS NOW THE CASE (3):

    } else {
        DIPOLE_T96( PS, X, Y, Z, &QX, &QY, &QZ, t );
        *BX = OIMFX - QX;
        *BY = OIMFY - QY;
        *BZ = OIMFZ - QZ;
    }

    return;

}


/*
 *   CALCULATES GSM COMPONENTS OF THE EXTERNAL MAGNETIC FIELD DUE TO
 *    SHIELDING OF THE EARTH'S DIPOLE_T96 ONLY
 */
void DIPSHLD_T96( double PS, double X, double Y, double Z, double *BX, double *BY, double *BZ, LgmTsyg1996_Info *tInfo ) {

    double HX, HY, HZ;
    double FX, FY, FZ;
    double A1[] = {-9e99,  0.24777, -27.003, -0.46815,  7.0637, -1.5918, -0.90317e-01, 57.522, 13.757, 2.0100, 10.458, 4.5798, 2.1695 };
    double A2[] = {-9e99, -0.65385, -18.061, -0.40457, -5.0995, 1.2846,   0.78231e-01, 39.592, 13.291, 1.9970, 10.062, 4.5140, 2.1558 };

    CYLHARM_T96( A1, X, Y, Z, &HX, &HY, &HZ );
    CYLHAR1_T96( A2, X, Y, Z, &FX, &FY, &FZ );

    *BX = HX*tInfo->cos_psi + FX*tInfo->sin_psi;
    *BY = HY*tInfo->cos_psi + FY*tInfo->sin_psi;
    *BZ = HZ*tInfo->cos_psi + FZ*tInfo->sin_psi;

    return;
}


/*
*  THIS CODE YIELDS THE SHIELDING FIELD FOR THE PERPENDICULAR DIPOLE_T96
*
*
*
*   ***  N.A. Tsyganenko ***  Sept. 14-18, 1993; revised March 16, 1994 ***
*
*   An approximation for the Chapman-Ferraro field by a sum of 6 cylin-
*   drical harmonics (see pp. 97-113 in the brown GSFC notebook #1)
*
*      Description of parameters:
*
*  A   - input vector containing model parameters;
*  X,Y,Z   -  input GSM coordinates
*  BX,BY,BZ - output GSM components of the shielding field
* - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*  The 6 linear parameters A(1)-A(6) are amplitudes of the cylindrical harmonic
*       terms.
*  The 6 nonlinear parameters A(7)-A(12) are the corresponding scale lengths
*       for each term (see GSFC brown notebook).
*
* - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*
*/
void CYLHARM_T96( double A[], double X, double Y, double Z, double *BX, double *BY, double *BZ ) {

    int     I;
    double  RHO, SINFI, COSFI, SINFI2, COSFI2, SI2CO2;
    double  DZETA, XJ0, XJ1, XEXP, XKSI_T96, BRHO, BPHI;

    RHO = sqrt( Y*Y + Z*Z );
    if ( RHO < 1.0e-8 ) {
        SINFI = 1.0;
        COSFI = 0.0;
        RHO   = 1.0e-8;
    } else {
        SINFI = Z/RHO;
        COSFI = Y/RHO;
    }



    SINFI2 = SINFI*SINFI;
    COSFI2 = COSFI*COSFI;
    SI2CO2 = SINFI2 - COSFI2;
    *BX = *BY = *BZ = 0.0;

    for ( I=1; I<=3; I++ ) {
        DZETA = RHO/A[I+6];
        XJ0   = BES_T96( DZETA, 0 );
        XJ1   = BES_T96( DZETA, 1 );
        XEXP  = exp( X/A[I+6] );

        *BX = *BX - A[I]*XJ1*XEXP*SINFI;
        *BY = *BY + A[I]*( 2.0*XJ1/DZETA - XJ0 ) * XEXP*SINFI*COSFI;
        *BZ = *BZ + A[I]*( XJ1/DZETA*SI2CO2 - XJ0*SINFI2 )*XEXP;
    }

    for ( I=4; I<=6; I++ ) {
        DZETA = RHO/A[I+6];
        XKSI_T96  = X/A[I+6];
        XJ0   = BES_T96( DZETA, 0 );
        XJ1   = BES_T96( DZETA, 1 );
        XEXP  = exp( XKSI_T96 );
        BRHO  = ( XKSI_T96*XJ0-(DZETA*DZETA + XKSI_T96 - 1.0 )*XJ1/DZETA ) * XEXP*SINFI;
        BPHI  = ( XJ0 + XJ1/DZETA*(XKSI_T96-1.0) ) * XEXP*COSFI;

        *BX += A[I]*( DZETA*XJ0 + XKSI_T96*XJ1 )*XEXP*SINFI;
        *BY += A[I]*( BRHO*COSFI - BPHI*SINFI );
        *BZ += A[I]*( BRHO*SINFI + BPHI*COSFI );
    }


    return;
}



/*
 *  THIS CODE YIELDS THE SHIELDING FIELD FOR THE PARALLEL DIPOLE_T96
 *
 *
 *
 *   ***  N.A. Tsyganenko ***  Sept. 14-18, 1993; revised March 16, 1994 ***
 *
 *   An approximation of the Chapman-Ferraro field by a sum of 6 cylin-
 *   drical harmonics (see pages 97-113 in the brown GSFC notebook #1)
 *
 *      Description of parameters:
 *
 *  A   - input vector containing model parameters;
 *  X,Y,Z - input GSM coordinates,
 *  BX,BY,BZ - output GSM components of the shielding field
 *
 * - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 *
 *      The 6 linear parameters A(1)-A(6) are amplitudes of the cylindrical
 *  harmonic terms.
 *      The 6 nonlinear parameters A(7)-A(12) are the corresponding scale
 *  lengths for each term (see GSFC brown notebook).
 *
 * - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 *
 */
void CYLHAR1_T96( double A[], double X, double Y, double Z, double *BX, double *BY, double *BZ ) {

    int     I;
    double  RHO, SINFI, COSFI, DZETA, XJ0, XJ1, XEXP, BRHO, XKSI_T96;

    RHO = sqrt( Y*Y + Z*Z );

    if ( RHO < 1.0e-10 ) {
        SINFI = 1.0;
        COSFI = 0.0;
    } else {
        SINFI = Z/RHO;
        COSFI = Y/RHO;
    }

    *BX = *BY = *BZ = 0.0;

    for ( I=1; I<=3; I++ ) {
        DZETA = RHO/A[I+6];
        XKSI_T96  = X/A[I+6];
        XJ0   = BES_T96( DZETA, 0 );
        XJ1   = BES_T96( DZETA, 1 );
        XEXP  = exp(XKSI_T96);
        BRHO  = XJ1*XEXP;
        *BX   = *BX - A[I]*XJ0*XEXP;
        *BY   = *BY + A[I]*BRHO*COSFI;
        *BZ   = *BZ + A[I]*BRHO*SINFI;
    }
    for ( I=4; I<=6; I++ ) {
        DZETA = RHO/A[I+6];
        XKSI_T96  = X/A[I+6];
        XJ0   = BES_T96(DZETA,0);
        XJ1   = BES_T96(DZETA,1);
        XEXP  = exp(XKSI_T96);
        BRHO  = (DZETA*XJ0+XKSI_T96*XJ1)*XEXP;
        *BX   += A[I]*( DZETA*XJ1 - XJ0*(XKSI_T96+1.0) ) * XEXP;
        *BY   += A[I]*BRHO*COSFI;
        *BZ   += A[I]*BRHO*SINFI;
    }


    return;
}


/*
 * Bessel functions of x of the first kind of order 0
 * Its much easier to just use the built-in math function j0()....
 */
double BES_T960_T96( double X ) {

    double X32, R, XD3, F0, T0;

    if ( fabs(X) < 3.0 ) {
        X32 = X*X/9.0;
        R = 1.0 - X32*( 2.2499997 - X32*(1.2656208 - X32* (0.3163866-X32*(0.0444479-X32*(0.0039444 -X32*0.00021)))));
    } else {
        XD3 = 3.0/X;
        F0 = 0.79788456-XD3*(0.00000077+XD3*(0.00552740+XD3* (0.00009512-XD3*(0.00137237-XD3*(0.00072805 -XD3*0.00014476)))));
        T0 = X - 0.78539816-XD3*(0.04166397+XD3*(0.00003954-XD3* (0.00262573-XD3*(0.00054125+XD3*(0.00029333 -XD3*0.00013558)))));
        R = F0/sqrt(X)*cos(T0);
    }

    return( R );

}

/*
 * Bessel functions of x of the first kind of order 1
 * Its much easier to just use the built-in math function j1()....
 */
double BES_T961_T96( double X ) {

    double X32, R, BES_T961_T96XM1, XD3, F1, T1;

    if ( fabs(X) < 3.0 ) {
        X32 = X*X/9.0;
        BES_T961_T96XM1 =0.5-X32*(0.56249985-X32*(0.21093573-X32* (0.03954289-X32*(0.00443319-X32*(0.00031761 -X32*0.00001109)))));
        R = BES_T961_T96XM1*X;
    } else {
        XD3  = 3.0/X;
        F1   =  0.79788456+XD3*(0.00000156+XD3*(0.01659667+XD3* (0.00017105-XD3*(0.00249511-XD3*(0.00113653 -XD3*0.00020033)))));
        T1   =  X-2.35619449+XD3*(0.12499612+XD3*(0.0000565-XD3* (0.00637879-XD3*(0.00074348+XD3*(0.00079824 -XD3*0.00029166)))));
        R = F1/sqrt(X)*cos(T1);
    }

    return( R );
}

/*
 * Bessel functions of x of the first kind of order K
 * Its much easier to just use the built-in math function jn(K,X)....
 */
double BES_T96( double X, int K ) {

    int N;
    double  R, G, XJN, XJNM1, XJNP1, SUM;

    if ( K == 0 ) return( BES_T960_T96( X ) );
    if ( K == 1 ) return( BES_T961_T96( X ) );
    if ( X == 0.0 ) return( 0.0 );

    G = 2.0/X;

    if ( X > (double)K ) {

        N = 1;
        XJN   = BES_T961_T96(X);
        XJNM1 = BES_T960_T96(X);

        for(;;){

            XJNP1 = G*N*XJN - XJNM1;
            ++N;

            if ( N >= K ) return( XJNP1 );

            XJNM1 = XJN;
            XJN   = XJNP1;
        }

    }

    N     = 24;
    XJN   = 1.0;
    XJNP1 = 0.0;

    for(;;) {
        if ( N%2 ==  0) SUM += XJN;
        XJNM1 = G*N*XJN - XJNP1;
        --N;
        XJNP1 = XJN;
        XJN   = XJNM1;
        if (N == K) R = XJN;

       if ( fabs(XJN) >  1.0e5 ) {
            XJNP1 = XJNP1*1.e-5;
            XJN   = XJN*1.0e-5;
            SUM   = SUM*1.0e-5;
            if ( N <= K )  R *= 1.e-5;
        }

        if ( N == 0 ) {
            SUM = XJN + 2.0*SUM;
            R /= SUM;
            return( R );
        }
    }

    return( R );

}







/*------------------------------------------------------------
 *
 *      Calculates the potential interconnection field inside the magnetosphere,
 *  corresponding to  DELTA_X = 20Re and DELTA_Y = 10Re (NB#3, p.90, 6/6/1996).
 *  The position (X,Y,Z) and field components BX,BY,BZ are given in the rotated
 *   coordinate system, in which the Z-axis is always directed along the BzIMF
 *   (i.e. rotated by the IMF clock angle Theta)
 *   It is also assumed that the IMF Bt=1, so that the components should be
 *     (i) multiplied by the actual Bt, and
 *     (ii) transformed to standard GSM coords by rotating back around X axis
 *              by the angle -Theta.
 *
 *      Description of parameters:
 *
 *     X,Y,Z -   GSM POSITION
 *      BX,BY,BZ - INTERCON_T96NECTION FIELD COMPONENTS INSIDE THE MAGNETOSPHERE
 *        OF A STANDARD SIZE (TO TAKE INTO ACCOUNT EFFECTS OF PRESSURE CHANGES,
 *         APPLY THE SCALING TRANSFORMATION)
 *
 * - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 *
 *     The 9 linear parameters are amplitudes of the "cartesian" harmonics
 *     The 6 nonlinear parameters are the scales Pi and Ri entering
 *    the arguments of exponents, sines, and cosines in the 9 "Cartesian"
 *       harmonics (3+3)
 * - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 */
void INTERCON_T96( double X, double Y, double Z, double *BX, double *BY, double *BZ, LgmTsyg1996_Info *tInfo ){
    int     L, I, K;
    double  CYPI, SYPI, SZRK, CZRK, SQPR, EPR, HX, HY, HZ;
    double  RP[4], RR[4], P[4], R[4];
    double  A[] = { -9e99, -8.411078731, 5932254.951, -9073284.93, -11.68794634, 6027598.824, -9218378.368, -6.508798398,
                           -11824.42793, 18015.66212, 7.99754043,    13.9669886, 90.24475036,  16.75728834,  1015.645781, 
                           1553.493216 };

    if ( tInfo->INTERCON_M_FLAG == 0 ) {
        tInfo->P[1] = A[10];
        tInfo->P[2] = A[11];
        tInfo->P[3] = A[12];
        tInfo->R[1] = A[13];
        tInfo->R[2] = A[14];
        tInfo->R[3] = A[15];
        for ( I=1; I<=3; I++ ) {
            tInfo->RP[I] = 1.0/tInfo->P[I];
            tInfo->RR[I] = 1.0/tInfo->R[I];
            for ( K=1; K<=3; K++ ) {
                tInfo->SQPR[I][K] = sqrt( RP[I]*RP[I] + RR[K]*RR[K] );
            }
        }
        tInfo->INTERCON_M_FLAG = 1; // flag that we've been here.
    }


    /*
     *        "PERPENDICULAR" KIND OF SYMMETRY ONLY
     */
    *BX = *BY = *BZ = 0.0;
    for ( L=0, I=1; I<=3; I++ ) {
        CYPI = cos( Y*RP[I] );
        SYPI = sin( Y*RP[I] );
 
        for ( K=1; K<=3; K++ ) {
            SZRK = sin( Z*RR[K] );
            CZRK = cos( Z*RR[K] );
            //SQPR = sqrt( RP[I]*RP[I] + RR[K]*RR[K] );
            SQPR = tInfo->SQPR[I][K];
            EPR  = exp( X*SQPR );
 
            HX = -SQPR*EPR*CYPI*SZRK;
            HY =  RP[I]*EPR*SYPI*SZRK;
            HZ = -RR[K]*EPR*CYPI*CZRK;
            ++L;
 
            *BX += A[L]*HX;
            *BY += A[L]*HY;
            *BZ += A[L]*HZ;
        }
    }

    return;

}


/*
 *  COMPUTES THE COMPONENTS OF THE FIELD OF THE MODEL RING CURRENT AND THREE
 *                   TAIL MODES WITH UNIT AMPLITUDES
 *      (FOR THE RING CURRENT, IT MEANS THE DISTURBANCE OF Bz=-1nT AT ORIGIN,
 *   AND FOR THE TAIL MODES IT MEANS MAXIMAL BX JUST ABOVE THE SHEET EQUAL 1 nT.
 */
void TAILRC96_T96( double SPS, double X, double Y, double Z, double *BXRC, double *BYRC, double *BZRC, double *BXT2, double *BYT2, double *BZT2, double *BXT3, double *BYT3, double *BZT3, LgmTsyg1996_Info *t ) {


    double  ARC[] = { -9e99,    -3.087699646,     3.516259114,        18.81380577,       -13.95772338,   -5.497076303,     0.1712890838,   
                                 2.392629189,    -2.728020808,       -14.79349936,        11.08738083,    4.388174084,     0.2492163197e-01,   
                                 0.7030375685,   -0.7966023165,       -3.835041334,        2.642228681,  -0.2405352424,   -0.7297705678, 
                                -0.3680255045,    0.1333685557,        2.795140897,       -1.078379954,   0.8014028630,    0.1245825565,   
                                 0.6149982835,   -0.2207267314,       -4.424578723,        1.730471572,  -1.716313926,    -0.2306302941,  
                                -0.2450342688,    0.8617173961e-01,    1.54697858,        -0.6569391113, -0.6537525353,    0.2079417515, 
                                12.75434981,     11.37659788,        636.4346279,          1.752483754,   3.604231143,    12.83078674, 
                                 7.412066636,     9.434625736,       676.7557193,          1.701162737,   3.580307144,    14.64298662         };

    double  ATAIL2[] = { -9e99,  .8747515218, -.9116821411, 2.209365387, -2.159059518, 
                               -7.059828867, 5.924671028, -1.916935691, 1.996707344, -3.877101873, 
                               3.947666061, 11.38715899, -8.343210833, 1.194109867, -1.244316975, 
                               3.73895491, -4.406522465, -20.66884863, 3.020952989, .2189908481, 
                               -.09942543549, -.927225562, .1555224669, .6994137909, -.08111721003, 
                               -.7565493881, .4686588792, 4.266058082, -.3717470262, -3.920787807, 
                               .02298569870, .7039506341, -.5498352719, -6.675140817, .8279283559, 
                               -2.234773608, -1.622656137, 5.187666221, 6.802472048, 39.13543412, 
                                2.784722096, 6.979576616, 25.71716760, 4.495005873, 8.068408272, 
                               93.47887103, 4.158030104, 9.313492566, 57.18240483 };

    double  ATAIL3[] = { -9e99, -19091.95061, -3011.613928, 20582.16203, 4242.918430, 
                               -2377.091102, -1504.820043, 19884.04650, 2725.150544, -21389.04845, 
                               -3990.475093, 2401.610097, 1548.171792, -946.5493963, 490.1528941, 
                               986.9156625, -489.3265930, -67.99278499, 8.711175710, -45.15734260, 
                               -10.76106500, 210.7927312, 11.41764141, -178.0262808, .7558830028, 
                                339.3806753, 9.904695974, 69.50583193, -118.0271581, 22.85935896, 
                               45.91014857, -425.6607164, 15.47250738, 118.2988915, 65.58594397, 
                               -201.4478068, -14.57062940, 19.69877970, 20.30095680, 86.45407420, 
                               22.50403727, 23.41617329, 48.48140573, 24.61031329, 123.5395974, 
                               223.5367692, 39.50824342, 65.83385762, 266.2948657 };

    double  WX, WY, WZ, HX, HY, HZ, RH, DR;
    double  DR2, g2, g, C11, C12, SPSC1, RPS, R, gp, gp2, gm, gm2, SQ2, C, CS;
    double  SPSS, CPSS, h, DPSRR, Y2, Y3, Y4, WFAC, W, WS, WARP, XS, ZSWW, ZS;
    double  DXSX, DXSY, DXSZ, DZSX, DZSY, DZSZ, D, DDDY, DZETAS, DDZETADX, DDZETADY, DDZETADZ;
    double  G, D0, DELTADY, C1, SQ1, h2;

    CPSS     = t->CB_T96_WARP.CPSS;
    SPSS     = t->CB_T96_WARP.SPSS;
    DPSRR    = t->CB_T96_WARP.DPSRR;
    RPS      = t->CB_T96_WARP.RPS;
    WARP     = t->CB_T96_WARP.WARP;
    D        = t->CB_T96_WARP.D;
    XS       = t->CB_T96_WARP.XS;
    ZS       = t->CB_T96_WARP.ZS;
    DXSX     = t->CB_T96_WARP.DXSX; 
    DXSY     = t->CB_T96_WARP.DXSY;
    DXSZ     = t->CB_T96_WARP.DXSZ;
    DZSX     = t->CB_T96_WARP.DZSX;
    DZSY     = t->CB_T96_WARP.DZSY;
    DZSZ     = t->CB_T96_WARP.DZSZ;
    DZETAS   = t->CB_T96_WARP.DZETAS;
    DDZETADX = t->CB_T96_WARP.DDZETADX;
    DDZETADY = t->CB_T96_WARP.DDZETADY;
    DDZETADZ = t->CB_T96_WARP.DDZETADZ;
    ZSWW     = t->CB_T96_WARP.ZSWW;


    
    RH      = t->CB_T96_RHDR.RH = 9.0; // Set in CB also
    DR      = t->CB_T96_RHDR.DR = 4.0; // Set in CB also
    G       =  10.0;
    D0      =  2.0;
    DELTADY =  10.0;


    /*
     *   TO ECONOMIZE THE CODE, WE FIRST CALCULATE COMMON VARIABLES, WHICH ARE
     *      THE SAME FOR ALL MODES, AND PUT THEM IN THE COMMON-BLOCK /WARP/
         COMMON /WARP/ CPSS,SPSS,DPSRR,RPS,WARP,D,XS,ZS,DXSX,DXSY,DXSZ, DZSX,DZSY,DZSZ,DZETAS,DDZETADX,DDZETADY,DDZETADZ,ZSWW
     */
    DR2   = DR*DR;
    g     = 1.0+RH; g2 = g*g;
    C11   = sqrt( g2 + DR2 );
    g     = 1.0-RH; g2 = g*g;
    C12   = sqrt( g2 + DR2 );
    C1    = C11 - C12;
    SPSC1 = SPS/C1;
    RPS   = 0.5*(C11+C12)*SPS;      //  THIS IS THE SHIFT OF OF THE SHEET WITH RESPECT TO GSM EQ.PLANE FOR THE 3RD (ASYMPTOTIC) TAIL MODE

    R     = sqrt( X*X + Y*Y + Z*Z );
    gp = R+RH; gp2 = g*g;
    SQ1   = sqrt( gp2 + DR2 );
    gm = R-RH; gm2 = gm*gm;
    SQ2   = sqrt( gm2 + DR2 );
    C     = SQ1-SQ2;
    CS    = gp/SQ1 - gm/SQ2;
    SPSS  = SPSC1/R * C;
    CPSS  = sqrt( 1.0 - SPSS*SPSS );
    g     = R*C1; g2 = g*g;
    h     = C*SPS; h2 = h*h;
    DPSRR = SPS/(R*R)*(CS*R-C) / sqrt( g2 - h2 );

    Y2    = Y*Y; Y3 = Y2*Y; Y4 = Y2*Y2;
    WFAC  = Y/(Y4 + 1.0e4);             //  WARPING
    W     = WFAC*Y3;
    WS    = 4.0e4*Y*WFAC*WFAC;
    WARP  = G*SPS*W;
    XS    = X*CPSS - Z*SPSS;
    ZSWW  = Z*CPSS + X*SPSS;            // "WW" MEANS "WITHOUT Y-Z WARPING" (IN X-Z ONLY)
    ZS    = ZSWW + WARP;

    DXSX  = CPSS-X*ZSWW*DPSRR;
    DXSY  = -Y*ZSWW*DPSRR;
    DXSZ  = -SPSS-Z*ZSWW*DPSRR;
    DZSX  = SPSS+X*XS*DPSRR;
    DZSY  = XS*Y*DPSRR  +G*SPS*WS;      //  THE LAST TERM IS FOR THE Y-Z WARP
    DZSZ  = CPSS+XS*Z*DPSRR;            //      (TAIL MODES ONLY)

    g     = Y/20.0; g2 = g*g;
    D     = D0 + DELTADY*g2;            //  SHEET HALF-THICKNESS FOR THE TAIL MODES
    DDDY  = DELTADY*Y*0.005;            //  (THICKENS TO FLANKS, BUT NO VARIATION ALONG X, IN CONTRAST TO RING CURRENT)

    DZETAS = sqrt( ZS*ZS + D*D );       //  THIS IS THE SAME SIMPLE WAY TO SPREAD OUT THE SHEET, AS THAT USED IN T89

    DDZETADX = ZS*DZSX/DZETAS;
    DDZETADY = (ZS*DZSY + D*DDDY)/DZETAS;
    DDZETADZ = ZS*DZSZ/DZETAS;


    SHLCAR3X3_T96( ARC, X, Y, Z, SPS, &WX, &WY, &WZ );
    RINGCURR96_T96( X, Y, Z, &HX, &HY, &HZ, t );
    *BXRC = WX + HX;
    *BYRC = WY + HY;
    *BZRC = WZ + HZ;

    SHLCAR3X3_T96( ATAIL2, X, Y, Z, SPS, &WX, &WY, &WZ );
    TAILDISK_T96( X, Y, Z, &HX, &HY, &HZ, t );
    *BXT2 = WX + HX;
    *BYT2 = WY + HY;
    *BZT2 = WZ + HZ;

    SHLCAR3X3_T96( ATAIL3, X, Y, Z, SPS, &WX, &WY, &WZ );
    TAIL87_T96( X, Z, &HX, &HZ, t );
    *BXT3 = WX + HX;
    *BYT3 = WY;
    *BZT3 = WZ + HZ;


    t->CB_T96_WARP.CPSS      = CPSS;
    t->CB_T96_WARP.SPSS      = SPSS;
    t->CB_T96_WARP.DPSRR     = DPSRR;
    t->CB_T96_WARP.RPS       = RPS;
    t->CB_T96_WARP.WARP      = WARP;
    t->CB_T96_WARP.D         = D;
    t->CB_T96_WARP.XS        = XS;
    t->CB_T96_WARP.ZS        = ZS;
    t->CB_T96_WARP.DXSX      = DXSX;
    t->CB_T96_WARP.DXSY      = DXSY;
    t->CB_T96_WARP.DXSZ      = DXSZ;
    t->CB_T96_WARP.DZSX      = DZSX;
    t->CB_T96_WARP.DZSY      = DZSY;
    t->CB_T96_WARP.DZSZ      = DZSZ;
    t->CB_T96_WARP.DZETAS    = DZETAS;
    t->CB_T96_WARP.DDZETADX  = DDZETADX;
    t->CB_T96_WARP.DDZETADY  = DDZETADY;
    t->CB_T96_WARP.DDZETADZ  = DDZETADZ;
    t->CB_T96_WARP.ZSWW      = ZSWW;

    return;

}



/*
 *
 *       THIS SUBROUTINE COMPUTES THE COMPONENTS OF THE RING CURRENT FIELD,
 *       SIMILAR TO THAT DESCRIBED BY TSYGANENKO AND PEREDO (1994).  THE
 *       DIFFERENCE IS THAT NOW WE USE SPACEWARPING, AS DESCRIBED IN THE
 *       PAPER ON MODELING BIRKELAND CURRENTS (TSYGANENKO AND STERN, 1996),
 *       INSTEAD OF SHEARING IT IN THE SPIRIT OF THE T89 TAIL MODEL.
 *
 *       IN  ADDITION, INSTEAD OF 7 TERMS FOR THE RING CURRENT MODEL, WE USE
 *       NOW ONLY 2 TERMS;  THIS SIMPLIFICATION ALSO GIVES RISE TO AN
 *       EASTWARD RING CURRENT LOCATED EARTHWARD FROM THE MAIN ONE,
 *       IN LINE WITH WHAT IS ACTUALLY OBSERVED
 *
 *       FOR DETAILS, SEE NB #3, PAGES 70-73
 */
void RINGCURR96_T96( double X, double Y, double Z, double *BX, double *BY, double *BZ, LgmTsyg1996_Info *t ) {

    int     I;
    double  BI, g, g2, hp, hp2, hm, hm2,S1, S2, DS1DDZ, DS2DDZ, DS1DRHOS, DS2DRHOS;
    double  DS1DX, DS1DY, DS1DZ, DS2DX, DS2DY, DS2DZ, S1TS2, S1PS2, S1PS2SQ, FAC1;
    double  AS, TERM1, FAC2, DASDS1, DASDS2, DASDX, DASDY, DASDZ, DZSY, XXD, h, h2;
    double  q,FDX, DDDX, D, DZETAS, RHOS, DDZETADX, DDZETADY, DDZETADZ, DRHOSDX, DRHOSDY, DRHOSDZ;
    double  CPSS, SPSS, DPSRR, XS, ZSWARPED, DXSX, DXSY, DXSZ, DZSX, DZSYWARPED, DZSZ, ZS;
    double  F[3], BETA[3], D0, DELTADX, XD, XLDX;


    CPSS        = t->CB_T96_WARP.CPSS;
    SPSS        = t->CB_T96_WARP.SPSS;
    DPSRR       = t->CB_T96_WARP.DPSRR;
    XS          = t->CB_T96_WARP.XS;
    ZSWARPED    = t->CB_T96_WARP.ZS;
    DXSX        = t->CB_T96_WARP.DXSX;
    DXSY        = t->CB_T96_WARP.DXSY;
    DXSZ        = t->CB_T96_WARP.DXSZ;
    DZSX        = t->CB_T96_WARP.DZSX;
    DZSYWARPED  = t->CB_T96_WARP.DZSY;
    DZSZ        = t->CB_T96_WARP.DZSZ;
    ZS          = t->CB_T96_WARP.ZSWW; // ZS HERE IS WITHOUT Y-Z WARP







    D0      = 2.0;
    DELTADX = 0.0; // ACHTUNG !!  THE RC IS NOW COMPLETELY SYMMETRIC (DELTADX=0)
    XD      = 0.0;
    XLDX    = 4.0;
 
    F[1]     = 569.895366;
    F[2]     = -1603.386993;
    BETA[1]  = 2.722188;
    BETA[2]  = 3.766875;
 

    /*
     *  THE ORIGINAL VALUES OF F(I) WERE MULTIPLIED BY BETA(I) (TO REDUCE THE
     *     NUMBER OF MULTIPLICATIONS BELOW)  AND BY THE FACTOR -0.43, NORMALIZING
     *      THE DISTURBANCE AT ORIGIN  TO  B=-1nT
     */
    DZSY = XS*Y*DPSRR;  // NO WARPING IN THE Y-Z PLANE (ALONG X ONLY), AND THIS IS WHY WE DO NOT USE  DZSY FROM THE COMMON-BLOCK
    XXD  = X-XD;
    g = XXD; g2 = g*g;
    h = XLDX; h2 = h*h;
    q = sqrt( g2 + h2 );
    FDX  = 0.5*(1.0+XXD/q);
    DDDX = DELTADX*0.5*h2 / (q*q*q);
    D    = D0 + DELTADX*FDX;


    DZETAS   = sqrt(ZS*ZS + D*D);  //  THIS IS THE SAME SIMPLE WAY TO SPREAD OUT THE SHEET, AS THAT USED IN T89
    RHOS     = sqrt(XS*XS + Y*Y);
    DDZETADX = (ZS*DZSX + D*DDDX)/DZETAS;
    g = ZS/DZETAS;
    DDZETADY = g*DZSY;
    DDZETADZ = g*DZSZ;


    if (RHOS < 1.0e-5) {
        DRHOSDX = 0.0;
        //DRHOSDY = DSIGN(1.0,Y);
        DRHOSDY = ( (Y>=0.0) ? +1.0 : -1.0);
        DRHOSDZ = 0.0;
    } else {
        g = XS/RHOS;
        DRHOSDX = g*DXSX;
        DRHOSDY = (XS*DXSY+Y)/RHOS;
        DRHOSDZ = g*DXSZ;
    }



    *BX = *BY = *BZ = 0.0;
    for ( I=1; I<=2; I++ ) {

        BI = BETA[I];
 
        g = DZETAS + BI; g2 = g*g;
        hp = RHOS + BI; hp2 = hp*hp;
        hm = RHOS - BI; hm2 = hm*hm;;
        S1       = sqrt( g2 + hp2 );
        S2       = sqrt( g2 + hm2 );
        DS1DDZ   = g/S1;
        DS2DDZ   = g/S2;
        DS1DRHOS = hp/S1;
        DS2DRHOS = hm/S2;
 
        DS1DX = DS1DDZ*DDZETADX + DS1DRHOS*DRHOSDX;
        DS1DY = DS1DDZ*DDZETADY + DS1DRHOS*DRHOSDY;
        DS1DZ = DS1DDZ*DDZETADZ + DS1DRHOS*DRHOSDZ;
 
        DS2DX = DS2DDZ*DDZETADX + DS2DRHOS*DRHOSDX;
        DS2DY = DS2DDZ*DDZETADY + DS2DRHOS*DRHOSDY;
        DS2DZ = DS2DDZ*DDZETADZ + DS2DRHOS*DRHOSDZ;
 
        S1TS2   = S1*S2;
        S1PS2   = S1 + S2;
        S1PS2SQ = S1PS2*S1PS2;
        FAC1    = sqrt( S1PS2SQ - 4.0*BI*BI );
        AS      = FAC1/(S1TS2*S1PS2SQ);
        TERM1   = 1.0/(S1TS2*S1PS2*FAC1);
        FAC2    = AS/S1PS2SQ;
        DASDS1  = TERM1 - FAC2/S1*(S2*S2 + S1*(3.0*S1 + 4.0*S2));
        DASDS2  = TERM1 - FAC2/S2*(S1*S1 + S2*(3.0*S2 + 4.0*S1));
 
        DASDX = DASDS1*DS1DX + DASDS2*DS2DX;
        DASDY = DASDS1*DS1DY + DASDS2*DS2DY;
        DASDZ = DASDS1*DS1DZ + DASDS2*DS2DZ;
 
        *BX = *BX + F[I]*((2.0*AS + Y*DASDY)*SPSS - XS*DASDZ + AS*DPSRR * (Y*Y*CPSS + Z*ZS) );
        *BY = *BY - F[I]*Y*(AS*DPSRR*XS + DASDZ*CPSS + DASDX*SPSS);
        *BZ = *BZ + F[I]*((2.0*AS + Y*DASDY)*CPSS + XS*DASDX - AS*DPSRR * (X*ZS + Y*Y*SPSS));


    }

    t->CB_T96_WARP.CPSS  = CPSS;
    t->CB_T96_WARP.SPSS  = SPSS;
    t->CB_T96_WARP.DPSRR = DPSRR;
    t->CB_T96_WARP.XS    = XS;
    t->CB_T96_WARP.ZS    = ZSWARPED;
    t->CB_T96_WARP.DXSX  = DXSX;
    t->CB_T96_WARP.DXSY  = DXSY;
    t->CB_T96_WARP.DXSZ  = DXSZ;
    t->CB_T96_WARP.DZSX  = DZSX;
    t->CB_T96_WARP.DZSY  = DZSYWARPED;
    t->CB_T96_WARP.DZSZ  = DZSZ;
    t->CB_T96_WARP.ZSWW  = ZS; // ZS HERE IS WITHOUT Y-Z WARP

    return;

}








/*
 *       THIS SUBROUTINE COMPUTES THE COMPONENTS OF THE TAIL CURRENT FIELD,
 *        SIMILAR TO THAT DESCRIBED BY TSYGANENKO AND PEREDO (1994).  THE
 *          DIFFERENCE IS THAT NOW WE USE SPACEWARPING, AS DESCRIBED IN OUR
 *           PAPER ON MODELING BIRKELAND CURRENTS (TSYGANENKO AND STERN, 1996)
 *            INSTEAD OF SHEARING IT IN THE SPIRIT OF T89 TAIL MODEL.
 *
 *          IN  ADDITION, INSTEAD OF 8 TERMS FOR THE TAIL CURRENT MODEL, WE USE
 *           NOW ONLY 4 TERMS
 *
 *             FOR DETAILS, SEE NB #3, PAGES 74-
 */
void  TAILDISK_T96( double X, double Y, double Z, double *BX, double *BY, double *BZ, LgmTsyg1996_Info *t ) {

    static double F[]    = { -9e99, -745796.7338, 1176470.141, -444610.529, -57508.01028 };
    static double BETA[] = { -9e99, 7.9250000, 8.0850000, 8.4712500, 27.89500 };
    static double XSHIFT = 4.5;

    int     I;
    double  CPSS, SPSS, DPSRR, XS, ZS, DXSX, DXSY, DXSZ, DZETAS, DDZETADX, DDZETADY, DDZETADZ, ZSWW;
    double  g, g2, RHOS, DRHOSDX, DRHOSDY, DRHOSDZ, hp, hp2, hm, hm2, S1, S2, DS1DDZ, DS2DDZ;
    double  DS1DRHOS, DS2DRHOS, DS1DX, DS1DY, DS1DZ, DS2DX, DS2DY, DS2DZ, S1TS2, S1PS2, S1PS2SQ, FAC1;
    double  BI, AS, TERM1, FAC2, DASDS1, DASDS2, DASDX, DASDY, DASDZ;

    CPSS     = t->CB_T96_WARP.CPSS;
    SPSS     = t->CB_T96_WARP.SPSS;
    DPSRR    = t->CB_T96_WARP.DPSRR;
    XS       = t->CB_T96_WARP.XS;
    ZS       = t->CB_T96_WARP.ZS;
    DXSX     = t->CB_T96_WARP.DXSX;
    DXSY     = t->CB_T96_WARP.DXSY;
    DXSZ     = t->CB_T96_WARP.DXSZ;
    DZETAS   = t->CB_T96_WARP.DZETAS;
    DDZETADX = t->CB_T96_WARP.DDZETADX;
    DDZETADY = t->CB_T96_WARP.DDZETADY;
    DDZETADZ = t->CB_T96_WARP.DDZETADZ;
    ZSWW     = t->CB_T96_WARP.ZSWW;



    /*
     *  here original F(I) are multiplied by BETA(I), to economize
     *    calculations
     */
    g  = XS-XSHIFT; g2 = g*g;
    RHOS = sqrt( g2 + Y*Y );
    if ( RHOS < 1.e-5 ) {
        DRHOSDX = 0.0;
        //DRHOSDY = DSIGN( 1.0, Y );
        DRHOSDY = ( (Y>=0.0) ? +1.0 : -1.0);
        DRHOSDZ = 0.0;
    } else {
        DRHOSDX = (XS-XSHIFT)*DXSX/RHOS;
        DRHOSDY = ((XS-XSHIFT)*DXSY+Y)/RHOS;
        DRHOSDZ = (XS-XSHIFT)*DXSZ/RHOS;
    }

    *BX = *BY = *BZ = 0.0;

    for ( I=1; I<=4; I++ ) {

        BI = BETA[I];

        g = DZETAS+BI; g2 = g*g;
        hp = RHOS+BI; hp2 = hp*hp;
        hm = RHOS-BI; hm2 = hm*hm;
        S1       = sqrt( g2 + hp2 );
        S2       = sqrt( g2 + hm2 );
        DS1DDZ   = g/S1;
        DS2DDZ   = g/S2;
        DS1DRHOS = hp/S1;
        DS2DRHOS = hm/S2;

        DS1DX = DS1DDZ*DDZETADX + DS1DRHOS*DRHOSDX;
        DS1DY = DS1DDZ*DDZETADY + DS1DRHOS*DRHOSDY;
        DS1DZ = DS1DDZ*DDZETADZ + DS1DRHOS*DRHOSDZ;
 
        DS2DX = DS2DDZ*DDZETADX + DS2DRHOS*DRHOSDX;
        DS2DY = DS2DDZ*DDZETADY + DS2DRHOS*DRHOSDY;
        DS2DZ = DS2DDZ*DDZETADZ + DS2DRHOS*DRHOSDZ;
 
        S1TS2   = S1*S2;
        S1PS2   = S1 + S2;
        S1PS2SQ = S1PS2*S1PS2;
        FAC1    = sqrt( S1PS2SQ - 4.0*BI*BI );
        AS      = FAC1/(S1TS2*S1PS2SQ);
        TERM1   = 1.0/(S1TS2*S1PS2*FAC1);
        FAC2    = AS/S1PS2SQ;
        DASDS1  = TERM1 - FAC2/S1*(S2*S2 + S1*(3.0*S1 + 4.0*S2));
        DASDS2  = TERM1 - FAC2/S2*(S1*S1 + S2*(3.0*S2 + 4.0*S1));

        DASDX = DASDS1*DS1DX + DASDS2*DS2DX;
        DASDY = DASDS1*DS1DY + DASDS2*DS2DY;
        DASDZ = DASDS1*DS1DZ + DASDS2*DS2DZ;

        *BX = *BX + F[I]*((2.0*AS+Y*DASDY)*SPSS-(XS-XSHIFT)*DASDZ + AS*DPSRR*(Y*Y*CPSS + Z*ZSWW));
        *BY = *BY - F[I]*Y*(AS*DPSRR*XS + DASDZ*CPSS + DASDX*SPSS);
        *BZ = *BZ + F[I]*((2.0*AS+Y*DASDY)*CPSS+(XS-XSHIFT)*DASDX - AS*DPSRR*(X*ZSWW + Y*Y*SPSS));
    }

    t->CB_T96_WARP.CPSS      = CPSS;
    t->CB_T96_WARP.SPSS      = SPSS;
    t->CB_T96_WARP.DPSRR     = DPSRR;
    t->CB_T96_WARP.XS        = XS;
    t->CB_T96_WARP.ZS        = ZS;
    t->CB_T96_WARP.DXSX      = DXSX;
    t->CB_T96_WARP.DXSY      = DXSY;
    t->CB_T96_WARP.DXSZ      = DXSZ;
    t->CB_T96_WARP.DZETAS    = DZETAS;
    t->CB_T96_WARP.DDZETADX  = DDZETADX;
    t->CB_T96_WARP.DDZETADY  = DDZETADY;
    t->CB_T96_WARP.DDZETADZ  = DDZETADZ;
    t->CB_T96_WARP.ZSWW      = ZSWW;

    return;

}





/*
 *      'LONG' VERSION OF THE 1987 TAIL MAGNETIC FIELD MODEL
 *              (N.A.TSYGANENKO, PLANET. SPACE SCI., V.35, P.1347, 1987)
 *
 *      D   IS THE Y-DEPENDENT SHEET HALF-THICKNESS (INCREASING TOWARDS FLANKS)
 *      RPS  IS THE TILT-DEPENDENT SHIFT OF THE SHEET IN THE Z-DIRECTION,
 *           CORRESPONDING TO THE ASYMPTOTIC HINGING DISTANCE, DEFINED IN THE
 *           MAIN SUBROUTINE (TAILRC96_T96) FROM THE PARAMETERS RH AND DR OF THE
 *           T96-TYPE MODULE, AND
 *      WARP  IS THE BENDING OF THE SHEET FLANKS IN THE Z-DIRECTION, DIRECTED
 *           OPPOSITE TO RPS, AND INCREASING WITH DIPOLE_T96 TILT AND |Y|
 */
void  TAIL87_T96( double X, double Z, double *BX, double *BZ, LgmTsyg1996_Info *t ) {

    double  ZS, ZP, ZM, XNX, XNX2, XC1, XC2, XC22, XR2, XC12, D2, B20, B2P, B2M;
    double  B, BP, BM, XA1, XAP1, XAM1, XA2, XAP2, XAM2, XNA, XNAP, XNAM, F, FP, FM;
    double  XLN1, XLNP1, XLNM1, XLN2, XLNP2, XLNM2, ALN, S0, S0P, S0M, S1, S1P, S1M;
    double  S2, S2P, S2M, G1, G1P, G1M, G2, G2P, G2M, RPS, WARP, D, DD, HPI, RT, XN;
    double  X1, X2, XNR, XN21, ADLN, B0, B1, B2;

    RPS      = t->CB_T96_WARP.RPS;
    WARP     = t->CB_T96_WARP.WARP;
    D        = t->CB_T96_WARP.D;





    DD   = 3.0;
    HPI  = 1.5707963; // MGH - half of Pi, really?
    RT   = 40.0;
    XN   = -10.0;
    X1   = -1.261;
    X2   = -0.663;
    B0   = 0.391734;
    B1   = 5.89715;
    B2   = 24.6833;
    XN21 = 76.37;
    XNR  = -0.1071;
    ADLN = 0.13238005;

    /*                !!!   THESE ARE NEW VALUES OF  X1, X2, B0, B1, B2,
     *                       CORRESPONDING TO TSCALE=1, INSTEAD OF TSCALE=0.6
     *
     *  THE ABOVE QUANTITIES WERE DEFINED AS FOLLOWS:------------------------
     *       HPI=PI/2
     *       RT=40.      !  Z-POSITION OF UPPER AND LOWER ADDITIONAL SHEETS
     *       XN=-10.     !  INNER EDGE POSITION
     *
     *       TSCALE=1  !  SCALING FACTOR, DEFINING THE RATE OF INCREASE OF THE
     *                       CURRENT DENSITY TAILWARDS
     *
     *  ATTENTION !  NOW I HAVE CHANGED TSCALE TO:  TSCALE=1.0, INSTEAD OF 0.6
     *                  OF THE PREVIOUS VERSION
     *
     *       B0=0.391734
     *       B1=5.89715 *TSCALE
     *       B2=24.6833 *TSCALE**2
     *
     *    HERE ORIGINAL VALUES OF THE MODE AMPLITUDES (P.77, NB#3) WERE NORMALIZED
     *      SO THAT ASYMPTOTIC  BX=1  AT X=-200RE
     *
     *      X1=(4.589  -5.85) *TSCALE -(TSCALE-1.)*XN ! NONLINEAR PARAMETERS OF THE
     *                                                         CURRENT FUNCTION
     *      X2=(5.187  -5.85) *TSCALE -(TSCALE-1.)*XN
     *
     *
     *      XN21=(XN-X1)**2
     *      XNR=1./(XN-X2)
     *      ADLN=-DLOG(XNR**2*XN21)
     *
     */

    ZS = Z - RPS + WARP;
    ZP = Z - RT;
    ZM = Z + RT;

    XNX  = XN-X;
    XNX2 = XNX*XNX;
    XC1  = X-X1;
    XC2  = X-X2;
    XC22 = XC2*XC2;
    XR2  = XC2*XNR;

    XC12 = XC1*XC1;
    D2   = DD*DD;    //  SQUARE OF THE TOTAL HALFTHICKNESS (DD=3Re for this mode)
    B20  = ZS*ZS + D2;
    B2P  = ZP*ZP + D2;
    B2M  = ZM*ZM + D2;

    B   = sqrt(B20);
    BP  = sqrt(B2P);
    BM  = sqrt(B2M);

    XA1  = XC12+B20;
    XAP1 = XC12+B2P;
    XAM1 = XC12+B2M;

    XA2  = 1.0/(XC22+B20);
    XAP2 = 1.0/(XC22+B2P);
    XAM2 = 1.0/(XC22+B2M);

    XNA  = XNX2+B20;
    XNAP = XNX2+B2P;
    XNAM = XNX2+B2M;

    F  = B20-XC22;
    FP = B2P-XC22;
    FM = B2M-XC22;

    XLN1  = log( XN21/XNA );
    XLNP1 = log( XN21/XNAP );
    XLNM1 = log( XN21/XNAM );

    XLN2  = XLN1+ADLN;
    XLNP2 = XLNP1+ADLN;
    XLNM2 = XLNM1+ADLN;

    ALN = 0.25*( XLNP1 + XLNM1 - 2.0*XLN1 );

    S0  = ( atan(XNX/B) + HPI )/B;
    S0P = ( atan(XNX/BP) + HPI )/BP;
    S0M = ( atan(XNX/BM) + HPI )/BM;

    S1  = (XLN1*0.5  + XC1*S0) /XA1;
    S1P = (XLNP1*0.5 + XC1*S0P)/XAP1;
    S1M = (XLNM1*0.5 + XC1*S0M)/XAM1;

    S2  = (XC2*XA2*XLN2 - XNR - F*XA2*S0)*XA2;
    S2P = (XC2*XAP2*XLNP2 - XNR - FP*XAP2*S0P)*XAP2;
    S2M = (XC2*XAM2*XLNM2 - XNR - FM*XAM2*S0M)*XAM2;

    G1  = (B20*S0 -0.5*XC1*XLN1)/XA1;
    G1P = (B2P*S0P -0.5*XC1*XLNP1)/XAP1;
    G1M = (B2M*S0M -0.5*XC1*XLNM1)/XAM1;

    G2  = ((0.5*F*XLN2 + 2.0*S0*B20*XC2)*XA2+XR2)*XA2;
    G2P = ((0.5*FP*XLNP2 + 2.0*S0P*B2P*XC2)*XAP2+XR2)*XAP2;
    G2M = ((0.5*FM*XLNM2 + 2.0*S0M*B2M*XC2)*XAM2+XR2)*XAM2;

    *BX = B0*(ZS*S0-0.5*(ZP*S0P+ZM*S0M)) + B1*(ZS*S1-0.5*(ZP*S1P+ZM*S1M)) + B2*(ZS*S2-0.5*(ZP*S2P+ZM*S2M));
    *BZ = B0*ALN+B1*(G1-0.5*(G1P+G1M)) + B2*(G2-0.5*(G2P+G2M));

    t->CB_T96_WARP.RPS  = RPS;;
    t->CB_T96_WARP.WARP = WARP;;
    t->CB_T96_WARP.D    = D;;

    return;

}





/*
 * THIS CODE RETURNS THE SHIELDING FIELD REPRESENTED BY  2x3x3=18 "CARTESIAN"
 *    HARMONICS
 *
 *
 * - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 *  The 36 coefficients enter in pairs in the amplitudes of the "cartesian"
 *    harmonics (A(1)-A(36).
 *  The 12 nonlinear parameters (A(37)-A(48) are the scales Pi,Ri,Qi,and Si
 *   entering the arguments of exponents, sines, and cosines in each of the
 *   18 "Cartesian" harmonics
 * - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 */
void SHLCAR3X3_T96( double A[] , double X, double Y, double Z, double SPS, double *HX, double *HY, double *HZ ) {

    int     L, M, K, N, I;
    double  CPS, S3PS, P, Q, CYPI, CYQI, SYPI, SYQI, R, S, SZRK, CZSK, CZRK, SZSK;
    double  SQPR, SQQS, EPR, EQS, DX, DY, DZ;
    

    CPS  = sqrt( 1.0-SPS*SPS );
    S3PS = 4.0*CPS*CPS - 1.0;   //  THIS IS SIN(3*PS)/SIN(PS)
 
    *HX = *HY = *HZ = 0.0;

    L = 0;
    for ( M=1; M<=2; M++ ) {     //    M=1 IS FOR THE 1ST SUM ("PERP." SYMMETRY) AND M=2 IS FOR THE SECOND SUM ("PARALL." SYMMETRY)
        for ( I=1; I<=3; I++ ) {

            P = A[36+I];
            Q = A[42+I];
            CYPI = cos(Y/P);
            CYQI = cos(Y/Q);
            SYPI = sin(Y/P);
            SYQI = sin(Y/Q);
 
            for ( K=1; K<=3; K++ ) {

                R = A[39+K];
                S = A[45+K];
                SZRK = sin(Z/R);
                CZSK = cos(Z/S);
                CZRK = cos(Z/R);
                SZSK = sin(Z/S);
                SQPR = sqrt(1.0/(P*P) + 1.0/(R*R));
                SQQS = sqrt(1.0/(Q*Q) + 1.0/(S*S));
                EPR  = exp(X*SQPR);
                EQS  = exp(X*SQQS);

                for ( N=1; N<=2; N++ ) { //  N=1 IS FOR THE FIRST PART OF EACH COEFFICIENT AND N=2 IS FOR THE SECOND ONE
                    ++L;
                    if ( M == 1 ) {
                        if ( N == 1 ) {
                            DX = -SQPR*EPR*CYPI*SZRK;
                            DY =  EPR/P*SYPI*SZRK;
                            DZ = -EPR/R*CYPI*CZRK;
                            *HX += A[L]*DX;
                            *HY += A[L]*DY;
                            *HZ += A[L]*DZ;
                        } else {
                            DX = DX*CPS;
                            DY = DY*CPS;
                            DZ = DZ*CPS;
                            *HX += A[L]*DX;
                            *HY += A[L]*DY;
                            *HZ += A[L]*DZ;
                        }
                    } else {
                        if ( N == 1 ) {
                            DX = -SPS*SQQS*EQS*CYQI*CZSK;
                            DY =  SPS*EQS/Q*SYQI*CZSK;
                            DZ =  SPS*EQS/S*CYQI*SZSK;
                            *HX += A[L]*DX;
                            *HY += A[L]*DY;
                            *HZ += A[L]*DZ;
                        } else {
                            DX = DX*S3PS;
                            DY = DY*S3PS;
                            DZ = DZ*S3PS;
                            *HX += A[L]*DX;
                            *HY += A[L]*DY;
                            *HZ += A[L]*DZ;
                        }
                    }

                }
            }
        }
    }

    return;

}





/*
*  THIS IS THE SECOND VERSION OF THE ANALYTICAL MODEL OF THE REGION 1 FIELD
*   BASED ON A SEPARATE REPRESENTATION OF THE POTENTIAL FIELD IN THE INNER AND
*   OUTER SPACE, MAPPED BY MEANS OF A SPHERO-DIPOLAR COORDINATE SYSTEM (NB #3,
*   P.91).   THE DIFFERENCE FROM THE FIRST ONE IS THAT INSTEAD OF OCTAGONAL
*   CURRENT LOOPS, CIRCULAR ONES ARE USED IN THIS VERSION FOR APPROXIMATING THE
*   FIELD IN THE OUTER REGION, WHICH IS FASTER.
*/
void  BIRK1TOT_02_T96( double PS, double X, double Y, double Z, double *BX, double *BY, double *BZ, LgmTsyg1996_Info *t ) {





    static double   C1[] = { -9e99, -0.911582E-03,-0.376654E-02,-0.727423E-02,-0.270084E-02,
                            -0.123899E-02,-0.154387E-02,-0.340040E-02,-0.191858E-01,
                            -0.518979E-01,0.635061E-01,0.440680,-0.396570,0.561238E-02,
                             0.160938E-02,-0.451229E-02,-0.251810E-02,-0.151599E-02,
                            -0.133665E-02,-0.962089E-03,-0.272085E-01,-0.524319E-01,
                             0.717024E-01,0.523439,-0.405015,-89.5587,23.2806 };

    static double   C2[] = { -9e99, 6.04133,.305415,.606066E-02,.128379E-03,-.179406E-04,
                            1.41714,-27.2586,-4.28833,-1.30675,35.5607,8.95792,.961617E-03,
                           -.801477E-03,-.782795E-03,-1.65242,-16.5242,-5.33798,.424878E-03,
                           .331787E-03,-.704305E-03,.844342E-03,.953682E-04,.886271E-03,
                           25.1120,20.9299,5.14569,-44.1670,-51.0672,-1.87725,20.2998,
                           48.7505,-2.97415,3.35184,-54.2921,-.838712,-10.5123,70.7594,
                           -4.94104,.106166E-03,.465791E-03,-.193719E-03,10.8439,-29.7968,
                            8.08068,.463507E-03,-.224475E-04,.177035E-03,-.317581E-03,
                           -.264487E-03,.102075E-03,7.71390,10.1915,-4.99797,-23.1114,
                          -29.2043,12.2928,10.9542,33.6671,-9.3851,.174615E-03,-.789777E-06,
                           .686047E-03,.460104E-04,-.345216E-02,.221871E-02,.110078E-01,
                           -.661373E-02,.249201E-02,.343978E-01,-.193145E-05,.493963E-05,
                           -.535748E-04,.191833E-04,-.100496E-03,-.210103E-03,-.232195E-02,
                           .315335E-02,-.134320E-01,-.263222E-01 };

    static double   XX1[] = { -9e99, -11.0, -7.0, -7.0, -3.0, -3.0,  1.0, 1.0, 1.0, 5.0, 5.0, 9.0, 9.0 };
    static double   YY1[] = { -9e99,   2.0,  0.0,  4.0,  2.0,  6.0,  0.0, 4.0, 8.0, 2.0, 6.0, 0.0, 4.0 } ;

    static double   XX2[] = { -9e99,   -10.0, -7.0, -4.0, -4.0, 0.0, 4.0,  4.0,  7.0, 10.0, 0.0, 0.0, 0.0, 0.0,  0.0 };
    static double   YY2[] = { -9e99,     3.0,  6.0,  3.0,  9.0, 6.0, 3.0,  9.0,  6.0,  3.0, 0.0, 0.0, 0.0, 0.0,  0.0 };
    static double   ZZ2[] = { -9e99,    20.0, 20.0,  4.0, 20.0, 4.0, 4.0, 20.0, 20.0, 20.0, 2.0, 3.0, 4.5, 7.0, 10.0 };


    static double   DTET0   = 0.034906;  //  THIS IS THE LATITUDINAL HALF-THICKNESS OF THE * R-1 OVAL (THE INTERPOLATION REGION BETWEEN * THE HIGH-LAT. AND THE PLASMA SHEET)
    static double   XLTDAY  = 78.0;      //  THESE ARE LATITUDES OF THE R-1 OVAL AT NOON AND AT MIDNIGHT
    static double   XLTNGHT = 70.0; 
    double          RH, DR;

    double D1[4][27], D2[3][79], XI[4];

    int     LOC, I;
    double  TNOONN, TNOONS, DTETDN, DR2, SPS, R, R2, R3, RMRH, RPRH, SQM, SQP, C, g, g2, h, h2, Q;
    double  SPSAS, CPSAS, XAS, ZAS, PAS, TAS, STAS, g3, g6, F, TET0, DTET, TETR1N, TETR1S;
    double  T01, T02, SQR, ST01AS, ST02AS, CT01AS, CT02AS, XAS1, Y1, ZAS1, X1, Z1, BX1, BY1, BZ1;
    double  XAS2, Y2, ZAS2, X2, Z2, BX2, BY2, BZ2, q, q2, SS, DS, FRAC, BSX, BSY, BSZ;
    


    t->CB_T96_DX1.DX       = -0.16;
    t->CB_T96_DX1.SCALEIN  =  0.08;
    t->CB_T96_DX1.SCALEOUT =  0.4;


    t->CB_T96_LOOPDIP1.TILT       = 1.00891;
    t->CB_T96_LOOPDIP1.XCENTRE[1] = 2.28397;
    t->CB_T96_LOOPDIP1.XCENTRE[2] = -5.60831;
    t->CB_T96_LOOPDIP1.RADIUS[1]  = 1.86106;
    t->CB_T96_LOOPDIP1.RADIUS[2]  = 7.83281;
    t->CB_T96_LOOPDIP1.DIPX       = 1.12541;
    t->CB_T96_LOOPDIP1.DIPY       = 0.945719;



    // RH IS THE "HINGING DISTANCE" AND DR IS THE * TRANSITION SCALE LENGTH, DEFINING THE * CURVATURE  OF THE WARPING (SEE P.89, NB #2)
    RH      = t->CB_T96_RHDR.RH = 9.0; // Set in CB also
    DR      = t->CB_T96_RHDR.DR = 4.0; // Set in CB also



    TNOONN = (90.0-XLTDAY)*0.01745329;
    TNOONS = 3.141592654-TNOONN; //HERE WE ASSUME THAT THE POSITIONS OF THE NORTHERN AND SOUTHERN R-1 OVALS ARE SYMMETRIC IN THE SM-COORDINATES
    DTETDN = (XLTDAY-XLTNGHT)*0.017453290;
    DR2    = DR*DR;

    SPS = sin(PS);
    R2  = X*X + Y*Y + Z*Z;
    R   = sqrt( R2 );
    R3  = R*R2;

    RMRH = R-RH;
    RPRH = R+RH;
    SQM  = sqrt( RMRH*RMRH + DR2 );
    SQP  = sqrt( RPRH*RPRH + DR2 );

    C = SQP-SQM;
    g  = RH+1.0; g2 = g*g;
    h  = RH-1.0; h2 = h*h;
    Q = sqrt( g2 + DR2 ) - sqrt( h2 + DR2 );

    SPSAS = SPS/R * C/Q;
    CPSAS = sqrt( 1.0 - SPSAS*SPSAS );
    XAS   = X*CPSAS - Z*SPSAS;
    ZAS   = X*SPSAS + Z*CPSAS;

    if ( (XAS != 0.0) || (Y != 0.0) ) {
        PAS = atan2( Y, XAS );
    } else {
        PAS = 0.0;
    }

    TAS  = atan2( sqrt(XAS*XAS + Y*Y), ZAS );
    STAS = sin(TAS);
    g = STAS; g2 = g*g; g3 = g2*g; g6 = g3*g3;
    F    = pow ( STAS/( g6*(1.0-R3) + R3 ), 0.1666666667 );

    TET0 = asin(F);
    if ( TAS > 1.5707963 ) TET0 = 3.141592654-TET0;
    
    g = sin(PAS*0.5); g2 = g*g; // cant we use a half-angle rel here?
    DTET   = DTETDN*g2;
    TETR1N = TNOONN+DTET;
    TETR1S = TNOONS-DTET;



    /*
     * NOW LET'S DEFINE WHICH OF THE FOUR REGIONS (HIGH-LAT., NORTHERN PSBL,
     *   PLASMA SHEET, SOUTHERN PSBL) DOES THE POINT (X,Y,Z) BELONG TO:
     */
    if ( ( TET0 <  (TETR1N-DTET0) ) || ( TET0 >  (TETR1S+DTET0) ) ) LOC = 1; // HIGH-LAT.
    if ( ( TET0 >  (TETR1N+DTET0) ) && ( TET0 <  (TETR1S-DTET0) ) ) LOC = 2; // PL.SHEET
    if ( ( TET0 >= (TETR1N-DTET0) ) && ( TET0 <= (TETR1N+DTET0) ) ) LOC = 3; // NORTH PSBL
    if ( ( TET0 >= (TETR1S-DTET0) ) && ( TET0 <= (TETR1S+DTET0) ) ) LOC = 4; // SOUTH PSBL


    if ( LOC == 1 ) {   // IN THE HIGH-LAT. REGION USE THE SUBROUTINE DIPOCT
        XI[1] = X;
        XI[2] = Y;
        XI[3] = Z;
        XI[4] = PS;
        DIPLOOP1_T96( XI, D1, XX1, YY1, t );
        *BX = *BY = *BZ = 0.0;
        for ( I=1; I<=26; I++ ) {
            *BX += C1[I]*D1[1][I];
            *BY += C1[I]*D1[2][I];
            *BZ += C1[I]*D1[3][I];
        }
    }

    if ( LOC == 2 ) {
        XI[1] = X;
        XI[2] = Y;
        XI[3] = Z;
        XI[4] = PS;
        CONDIP1( XI, D2, XX2, YY2, ZZ2 );
        *BX = *BY = *BZ = 0.0;
        for ( I=1; I<=79; I++ ) {
            *BX += C2[I]*D2[1][I];
            *BY += C2[I]*D2[2][I];
            *BZ += C2[I]*D2[3][I];
        }
    } //   END OF THE CASE 2


    if ( LOC == 3 ) {
        T01 = TETR1N-DTET0;
        T02 = TETR1N+DTET0;
        SQR = sqrt(R);

        g = sin(T01); g2  =g*g; g3 = g*g2; g6 = g3*g3;
        ST01AS = pow( SQR/(R3+1.0/g6-1.0), 0.1666666667 );
        g = sin(T02); g2  =g*g; g3 = g*g2; g6 = g3*g3;
        ST02AS = pow( SQR/(R3+1.0/g6-1.0), 0.1666666667 );
        CT01AS = sqrt(1.0-ST01AS*ST01AS);
        CT02AS = sqrt(1.0-ST02AS*ST02AS);
        XAS1   = R*ST01AS*cos(PAS);
        Y1     = R*ST01AS*sin(PAS);
        ZAS1   = R*CT01AS;
        X1     = XAS1*CPSAS+ZAS1*SPSAS;
        Z1     = -XAS1*SPSAS+ZAS1*CPSAS; // X1,Y1,Z1 ARE COORDS OF THE NORTHERN BOUNDARY POINT

        XI[1] = X1;
        XI[2] = Y1;
        XI[3] = Z1;
        XI[4] = PS;
        DIPLOOP1_T96( XI, D1, XX1, YY1, t );
        BX1 = BY1 = BZ1 = 0.0;
        for ( I=1; I<=26; I++ ) {
            BX1 += C1[I]*D1[1][I];  //   BX1,BY1,BZ1  ARE FIELD COMPONENTS
            BY1 += C1[I]*D1[2][I];  //  IN THE NORTHERN BOUNDARY POINT
            BZ1 += C1[I]*D1[3][I];  //
        }

        XAS2 = R*ST02AS*cos(PAS);
        Y2   = R*ST02AS*sin(PAS);
        ZAS2 = R*CT02AS;
        X2   = XAS2*CPSAS+ZAS2*SPSAS;
        Z2   = -XAS2*SPSAS+ZAS2*CPSAS; // X2,Y2,Z2 ARE COORDS OF THE SOUTHERN BOUNDARY POINT

        XI[1] = X2;
        XI[2] = Y2;
        XI[3] = Z2;
        XI[4] = PS;
        CONDIP1( XI, D2, XX2, YY2, ZZ2 );
        BX2 = BY2 = BZ2 = 0.0;
        for ( I=1; I<=79; I++ ) {
            BX2 += C2[I]*D2[1][I]; // BX2,BY2,BZ2  ARE FIELD COMPONENTS IN THE SOUTHERN BOUNDARY POINT
            BY2 += C2[I]*D2[2][I];
            BZ2 += C2[I]*D2[3][I];
        }

        /*
         *  NOW INTERPOLATE:
         */
        g = X2-X1; g2 = g*g; h = Y2-Y1; h2 = h*h; q = Z2-Z1; q2 = q*q;
        SS = sqrt( g2 + h2 + q2 );
        g = X-X1; g2 = g*g; h = Y-Y1; h2 = h*h; q = Z-Z1; q2 = q*q;
        DS = sqrt( g2 + h2 + q2 );
        FRAC = DS/SS;
        g = 1.0-FRAC;
        *BX = BX1*g + BX2*FRAC;
        *BY = BY1*g + BY2*FRAC;
        *BZ = BZ1*g + BZ2*FRAC;

    } // END OF THE CASE 3

    if ( LOC == 4 ) {
        T01 = TETR1S-DTET0;
        T02 = TETR1S+DTET0;
        SQR = sqrt(R);

        g = sin(T01); g2 = g*g; g3 = g2*g; g6 = g3*g3;
        ST01AS = pow( SQR/(R3+1.0/g6-1.0), 0.1666666667 );
        g = sin(T02); g2 = g*g; g3 = g2*g; g6 = g3*g3;
        ST02AS = pow( SQR/(R3+1.0/g6-1.0), 0.1666666667 );
        CT01AS = -sqrt(1.0-ST01AS*ST01AS);
        CT02AS = -sqrt(1.0-ST02AS*ST02AS);
        XAS1   = R*ST01AS*cos(PAS);
        Y1     = R*ST01AS*sin(PAS);
        ZAS1   = R*CT01AS;
        X1     = XAS1*CPSAS+ZAS1*SPSAS;
        Z1     = -XAS1*SPSAS+ZAS1*CPSAS; // X1,Y1,Z1 ARE COORDS OF THE NORTHERN BOUNDARY POINT

        XI[1] = X1;
        XI[2] = Y1;
        XI[3] = Z1;
        XI[4] = PS;
        CONDIP1( XI, D2, XX2, YY2, ZZ2 );
        BX1 = BY1 = BZ1 = 0.0;
        for ( I=1; I<=79; I++ ) {
            BX1 += C2[I]*D2[1][I]; // !  BX1,BY1,BZ1  ARE FIELD COMPONENTS IN THE NORTHERN BOUNDARY POINT
            BY1 += C2[I]*D2[2][I];
            BZ1 += C2[I]*D2[3][I];
        }

        XAS2 = R*ST02AS*cos(PAS);
        Y2   = R*ST02AS*sin(PAS);
        ZAS2 = R*CT02AS;
        X2   = XAS2*CPSAS+ZAS2*SPSAS;
        Z2   = -XAS2*SPSAS+ZAS2*CPSAS; // X2,Y2,Z2 ARE COORDS OF THE SOUTHERN BOUNDARY POINT
        XI[1] = X2;
        XI[2] = Y2;
        XI[3] = Z2;
        XI[4] = PS;
        DIPLOOP1_T96( XI, D1, XX1, YY1, t );
        BX2 = BY2 = BZ2 = 0.0;
        for ( I=1; I<=26; I++ ) {
            BX2 += C1[I]*D1[1][I]; //  BX2,BY2,BZ2  ARE FIELD COMPONENTS IN THE SOUTHERN BOUNDARY POINT
            BY2 += C1[I]*D1[2][I];
            BZ2 += C1[I]*D1[3][I];
        }


        /*
         *  NOW INTERPOLATE:
         */
        g = X2-X1; g2 = g*g; h = Y2-Y1; h2 = h*h; q = Z2-Z1; q2 = q*q;
        SS   = sqrt( g2 + h2 + q2 );
        g = X-X1; g2 = g*g; h = Y-Y1; h2 = h*h; q = Z-Z1; q2 = q*q;
        DS   = sqrt( g2 + h2 + q2 );
        FRAC = DS/SS;
        g = 1.0 - FRAC;
        *BX = BX1*g + BX2*FRAC;
        *BY = BY1*g + BY2*FRAC;
        *BZ = BZ1*g + BZ2*FRAC;

    } // END OF THE CASE 4



    /*
     *   NOW, LET US ADD THE SHIELDING FIELD
     */
    BIRK1SHLD_T96( PS, X, Y, Z, &BSX, &BSY, &BSZ );
    *BX += BSX;
    *BY += BSY;
    *BZ += BSZ;

    return;

}



/*
 *
 *  Calculates dependent model variables and their deriva-
 *  tives for given independent variables and model parame-
 *  ters.  Specifies model functions with free parameters which
 *  must be determined by means of least squares fits (RMS
 *  minimization procedure).
 *
 *      Description of parameters:
 *
 *  XI  - input vector containing independent variables;
 *  D   - output double precision vector containing
 *        calculated values for derivatives of dependent
 *        variables with respect to LINEAR model parameters;
 *
 * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 *
 *  The  26 coefficients are moments (Z- and X-components) of 12 dipoles placed
 *    inside the  R1-shell,  PLUS amplitudes of two octagonal double loops.
 *     The dipoles with nonzero  Yi appear in pairs with equal moments.
 *                  (see the notebook #2, pp.102-103, for details)
 * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 */

void DIPLOOP1_T96( double XI[5], double D[4][27], double *XX, double *YY, LgmTsyg1996_Info *t ) {

    int     I;
    double  X, Y, Z, PS, SPS, RH, DR;
    double  g, g2, R2, R, RMRH, RPRH, DR2, SQM, SQP, C, h, h2, Q, SPSAS, CPSAS, XD, YD, ZD;
    double  BX1X, BY1X, BZ1X, BX1Y, BY1Y, BZ1Y, BX1Z, BY1Z, BZ1Z;
    double  BX2X, BY2X, BZ2X, BX2Y, BY2Y, BZ2Y, BX2Z, BY2Z, BZ2Z;
    double  XOCT1, YOCT1, ZOCT1, XOCT2, YOCT2, ZOCT2, BX, BY, BZ, BXOCT1, BYOCT1, BZOCT1;
    double  DIPX, DIPY;

    RH = t->CB_T96_RHDR.RH;
    DR = t->CB_T96_RHDR.DR;

    DIPX = t->CB_T96_LOOPDIP1.DIPX;
    DIPY = t->CB_T96_LOOPDIP1.DIPY;


    X   = XI[1];
    Y   = XI[2];
    Z   = XI[3];
    PS  = XI[4];
    SPS = sin(PS);

    for ( I=1; I<=12; I++ ) {
        g = XX[I]*DIPX; g2 = g*g; h = YY[I]*DIPY;
        R2 =  g2 + h2;
        R     = sqrt( R2 );
        RMRH  = R-RH;
        RPRH  = R+RH;
        DR2   = DR*DR;;
        SQM   = sqrt( RMRH*RMRH + DR2 );
        SQP   = sqrt( RPRH*RPRH + DR2 );
        C     = SQP-SQM;
        g = RH+1.0; g2 = g*g;
        h = RH-1.0; h2 = h*h;
        Q     = sqrt( g2 + DR2 ) - sqrt( h2 + DR2 );
        SPSAS = SPS/R * C/Q;
        CPSAS = sqrt( 1.0-SPSAS*SPSAS );
        XD =  ( XX[I]*DIPX )*CPSAS;
        YD =  ( YY[I]*DIPY );
        ZD = -( XX[I]*DIPX )*SPSAS;
        DIPXYZ_T96( X-XD, Y-YD, Z-ZD, &BX1X, &BY1X, &BZ1X, &BX1Y, &BY1Y, &BZ1Y, &BX1Z, &BY1Z, &BZ1Z);
        if ( fabs(YD) > 1.0e-10 ) {
            DIPXYZ_T96( X-XD, Y+YD, Z-ZD, &BX2X, &BY2X, &BZ2X, &BX2Y, &BY2Y, &BZ2Y, &BX2Z, &BY2Z, &BZ2Z );
        } else {
            BX2X = BY2X = BZ2X = 0.0;
            BX2Z = BY2Z = BZ2Z = 0.0;
        }

        D[1][I]    = BX1Z+BX2Z;
        D[2][I]    = BY1Z+BY2Z;
        D[3][I]    = BZ1Z+BZ2Z;
        D[1][I+12] = (BX1X+BX2X)*SPS;
        D[2][I+12] = (BY1X+BY2X)*SPS;
        D[3][I+12] = (BZ1X+BZ2X)*SPS;
    }

    g = t->CB_T96_LOOPDIP1.XCENTRE[1] + t->CB_T96_LOOPDIP1.RADIUS[1];
    R2 = g*g;
    R     = sqrt(R2);
    RMRH  = R-RH;
    RPRH  = R+RH;
    DR2   = DR*DR;
    SQM   = sqrt( RMRH*RMRH + DR2 );
    SQP   = sqrt( RPRH*RPRH + DR2 );
    C     = SQP-SQM;
    g = RH+1.0; g2 = g*g; h = RH-1.0; h2 = h*h;
    Q     = sqrt( g2 + DR2 ) - sqrt( h2 + DR2 );
    SPSAS = SPS/R * C/Q;
    CPSAS = sqrt( 1.0 - SPSAS*SPSAS );
    XOCT1 = X*CPSAS - Z*SPSAS;
    YOCT1 = Y;
    ZOCT1 = X*SPSAS + Z*CPSAS;


    CROSSLP_T96( XOCT1, YOCT1, ZOCT1, &BXOCT1, &BYOCT1, &BZOCT1, t->CB_T96_LOOPDIP1.XCENTRE[1], t->CB_T96_LOOPDIP1.RADIUS[1], t->CB_T96_LOOPDIP1.TILT );
    D[1][25] =  BXOCT1*CPSAS + BZOCT1*SPSAS;
    D[2][25] =  BYOCT1;
    D[3][25] = -BXOCT1*SPSAS + BZOCT1*CPSAS;

    g = t->CB_T96_LOOPDIP1.RADIUS[2] - t->CB_T96_LOOPDIP1.XCENTRE[2];
    R2 = g*g;
    R  = sqrt(R2);
    RMRH = R-RH;
    RPRH = R+RH;
    DR2  = DR*DR;;
    SQM = sqrt( RMRH*RMRH + DR2 );
    SQP = sqrt( RPRH*RPRH + DR2 );
    C   = SQP-SQM;
    g = RH+1.0; g2 = g*g; h = RH-1.0; h2 = h*h;
    Q   = sqrt( g2 + DR2) - sqrt( h2 + DR2 );
    SPSAS = SPS/R * C/Q;
    CPSAS = sqrt(1.0-SPSAS*SPSAS);
    XOCT2 = X*CPSAS - Z*SPSAS - t->CB_T96_LOOPDIP1.XCENTRE[2];
    YOCT2 = Y;
    ZOCT2 = X*SPSAS + Z*CPSAS;
    CIRCLE_T96( XOCT2, YOCT2, ZOCT2, t->CB_T96_LOOPDIP1.RADIUS[2], &BX, &BY, &BZ );

    D[1][26] =  BX*CPSAS+BZ*SPSAS;
    D[2][26] =  BY;
    D[3][26] = -BX*SPSAS+BZ*CPSAS;

    return;

}






/*
 *  RETURNS COMPONENTS OF THE FIELD FROM A CIRCULAR CURRENT LOOP OF RADIUS RL
 *  USES THE SECOND (MORE ACCURATE) APPROXIMATION GIVEN IN ABRAMOWITZ AND STEGUN
 */
void CIRCLE_T96( double X, double Y, double Z, double RL, double *BX, double *BY, double *BZ ) {

    double  PI, RHO2, RHO, g, g2, R22, R2, R12, R32, XK2, XK2S, DL, K, E, BRHO;

    PI = 3.141592654;
 
    RHO2 = X*X + Y*Y;
    RHO  = sqrt( RHO2 );
    g    = RHO+RL; g2 = g*g;
    R22  = Z*Z + g2;
    R2   = sqrt(R22);
    R12  = R22 - 4.0*RHO*RL;
    R32  = 0.5*(R12+R22);
    XK2  = 1.0-R12/R22;
    XK2S = 1.0-XK2;
    DL   = log(1.0/XK2S);
    K    = 1.38629436112 
            + XK2S*(0.09666344259+XK2S*(0.03590092383+ XK2S*(0.03742563713+XK2S*0.01451196212))) 
            + DL*(0.5+XK2S*(0.12498593597+XK2S*(0.06880248576+ XK2S*(0.03328355346+XK2S*0.00441787012))));
    E    = 1.0 
            + XK2S*(0.44325141463+XK2S*(0.0626060122+XK2S* (0.04757383546+XK2S*0.01736506451))) 
            + DL*XK2S*(0.2499836831+XK2S*(0.09200180037+XK2S* (0.04069697526+XK2S*0.00526449639)));

    if ( RHO > 1.0e-6 ) {
        BRHO = Z/(RHO2*R2)*(R32/R12*E-K); // THIS IS NOT EXACTLY THE B-RHO COMPONENT - NOTE THE ADDITIONAL
    } else {           
        BRHO=PI*RL/R2*(RL-RHO)/R12*Z/(R32-RHO2);  // DIVISION BY RHO
    }

    *BX = BRHO*X;
    *BY = BRHO*Y;
    *BZ = (K-E*(R32-2.0*RL*RL)/R12)/R2;

    return;
}



/*
 *   RETURNS FIELD COMPONENTS OF A PAIR OF LOOPS WITH A COMMON CENTER AND
 *    DIAMETER,  COINCIDING WITH THE X AXIS. THE LOOPS ARE INCLINED TO THE
 *    EQUATORIAL PLANE BY THE ANGLE AL (RADIANS) AND SHIFTED IN THE POSITIVE
 *     X-DIRECTION BY THE DISTANCE  XC.
 */
void CROSSLP_T96( double X, double Y, double Z, double *BX, double *BY, double *BZ, double XC, double RL, double AL ) {

    double  CAL, SAL, Y1, Z1, Y2, Z2, BX1, BY1, BZ1, BX2, BY2, BZ2;

    CAL = cos(AL);
    SAL = sin(AL);

    Y1 =  Y*CAL - Z*SAL;
    Z1 =  Y*SAL + Z*CAL;
    Y2 =  Y*CAL + Z*SAL;
    Z2 = -Y*SAL + Z*CAL;
    CIRCLE_T96( X-XC, Y1, Z1, RL, &BX1, &BY1, &BZ1 );
    CIRCLE_T96( X-XC, Y2, Z2, RL, &BX2, &BY2, &BZ2 );
    *BX =   BX1+BX2;
    *BY =  (BY1+BY2)*CAL + (BZ1-BZ2)*SAL;
    *BZ = -(BY1-BY2)*SAL + (BZ1+BZ2)*CAL;

    return;

}



/*
 *       RETURNS THE FIELD COMPONENTS PRODUCED BY THREE DIPOLE_T96S, EACH
 *        HAVING M=Me AND ORIENTED PARALLEL TO X,Y, and Z AXIS, RESP.
 */
void DIPXYZ_T96( double X, double Y, double Z, double *BXX, double *BYX, double *BZX, double *BXY, double *BYY, double *BZY, double *BXZ, double *BYZ, double *BZZ ) {

    double  X2, Y2, Z2, R2, XMR5, XMR53;

    X2 = X*X;
    Y2 = Y*Y;
    Z2 = Z*Z;
    R2 = X2 + Y2 + Z2;

    XMR5  = 30574.0/(R2*R2*sqrt(R2));
    XMR53 = 3.0*XMR5;
    *BXX   = XMR5*(3.0*X2-R2);
    *BYX   = XMR53*X*Y;
    *BZX   = XMR53*X*Z;

    *BXY = *BYX;
    *BYY = XMR5*(3.0*Y2-R2);
    *BZY = XMR53*Y*Z;

    *BXZ = *BZX;
    *BYZ = *BZY;
    *BZZ = XMR5*(3.0*Z2 - R2);

    return;

}




/*
 *      Calculates dependent model variables and their derivatives for given
 *  independent variables and model parameters.  Specifies model functions with
 *  free parameters which must be determined by means of least squares fits
 *  (RMS minimization procedure).
 *
 *      Description of parameters:
 *
 *  XI  - input vector containing independent variables;
 *  D   - output double precision vector containing
 *        calculated values for derivatives of dependent
 *        variables with respect to LINEAR model parameters;
 *
 * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 *
 *  The  79 coefficients are (1) 5 amplitudes of the conical harmonics, plus
 *                           (2) (9x3+5x2)x2=74 components of the dipole moments
 *              (see the notebook #2, pp.113-..., for details)
 * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 */
CONDIP1( double XI[5], double D[4][80], double *XX, double *YY, double *ZZ, LgmTsyg1996_Info *t ) {


    int     M, I, IX, IY, IZ;
    double  a, b, g, h, CH2, SH2;
    double  BX1X, BY1X, BZ1X, BX1Y, BY1Y, BZ1Y, BX1Z, BY1Z, BZ1Z;
    double  BX2X, BY2X, BZ2X, BX2Y, BY2Y, BZ2Y, BX2Z, BY2Z, BZ2Z;
    double  BX3X, BY3X, BZ3X, BX3Y, BY3Y, BZ3Y, BX3Z, BY3Z, BZ3Z;
    double  BX4X, BY4X, BZ4X, BX4Y, BY4Y, BZ4Y, BX4Z, BY4Z, BZ4Z;
    double  X, Y, Z, PS, SPS, CPS, XSM, ZSM, RO2, RO, R2;
    double  R, C, S, CH, SH, TNH, CNH, BT, BF, BXSM, BY, BZSM, XD;
    double  YD, ZD; 
    double  CF[6], SF[6];


 
    X   = XI[1];
    Y   = XI[2];
    Z   = XI[3];
    PS  = XI[4];
    SPS = sin(PS);
    CPS = cos(PS);

    XSM = X*CPS - Z*SPS  - t->CB_T96_DX1.DX;
    ZSM = Z*CPS + X*SPS;
    RO2 = XSM*XSM + Y*Y;
    RO  = sqrt(RO2);

    CF[1] = XSM/RO;
    SF[1] = Y/RO;

    CF[2] = CF[1]*CF[1] - SF[1]*SF[1];
    SF[2] = 2.0*SF[1]*CF[1];
    CF[3] = CF[2]*CF[1] - SF[2]*SF[1];
    SF[3] = SF[2]*CF[1] + CF[2]*SF[1];
    CF[4] = CF[3]*CF[1] - SF[3]*SF[1];
    SF[4] = SF[3]*CF[1] + CF[3]*SF[1];
    CF[5] = CF[4]*CF[1] - SF[4]*SF[1];
    SF[5] = SF[4]*CF[1] + CF[4]*SF[1];

    R2  = RO2 + ZSM*ZSM;
    R   = sqrt(R2);
    C   = ZSM/R;
    S   = RO/R;
    CH  = sqrt(0.5*(1.0+C)); CH2 = CH*CH;
    SH  = sqrt(0.5*(1.0-C)); SH2 = SH*SH;
    TNH = SH/CH;
    CNH = 1.0/TNH;

    g = TNH; h = CNH;
    a = 1.0; b = 1.0;


    for ( M=1; M<=5; M++ ) {


        //BT      = M*CF[M]/(R*S) * (TNH**M+CNH**M);
        //BF      = -0.5*M*SF[M]/R* (TNH**(M-1) /CH2 - CNH**(M-1)/SH2 );
        BT      = M*CF[M]/(R*S) * (g + h);
        BF      = -0.5*M*SF[M]/R* (a/CH2 - b/SH2 );

        BXSM    = BT*C*CF[1]-BF*SF[1];
        BY      = BT*C*SF[1]+BF*CF[1];
        BZSM    = -BT*S;
        D[1][M] =  BXSM*CPS + BZSM*SPS;
        D[2][M] =  BY;
        D[3][M] = -BXSM*SPS + BZSM*CPS;

        g *= TNH;
        h *= CNH;
        a *= TNH;
        b *= CNH;
    }


    XSM = X*CPS - Z*SPS;
    ZSM = Z*CPS + X*SPS;

    for ( I=1; I<=9; I++ ) {
        if ( (I==3) || (I==5) || (I==6) ) {
            XD = XX[I]*t->CB_T96_DX1.SCALEIN;
            YD = YY[I]*t->CB_T96_DX1.SCALEIN;
        } else {
            XD = XX[I]*t->CB_T96_DX1.SCALEOUT;
            YD = YY[I]*t->CB_T96_DX1.SCALEOUT;
        }

        ZD = ZZ[I];

        DIPXYZ_T96( XSM-XD, Y-YD, ZSM-ZD, &BX1X, &BY1X, &BZ1X, &BX1Y, &BY1Y, &BZ1Y, &BX1Z, &BY1Z, &BZ1Z );
        DIPXYZ_T96( XSM-XD, Y+YD, ZSM-ZD, &BX2X, &BY2X, &BZ2X, &BX2Y, &BY2Y, &BZ2Y, &BX2Z, &BY2Z, &BZ2Z );
        DIPXYZ_T96( XSM-XD, Y-YD, ZSM+ZD, &BX3X, &BY3X, &BZ3X, &BX3Y, &BY3Y, &BZ3Y, &BX3Z, &BY3Z, &BZ3Z );
        DIPXYZ_T96( XSM-XD, Y+YD, ZSM+ZD, &BX4X, &BY4X, &BZ4X, &BX4Y, &BY4Y, &BZ4Y, &BX4Z, &BY4Z, &BZ4Z );

        IX = I*3 + 3;
        IY = IX + 1;
        IZ = IY + 1;
     
        D[1][IX] = (BX1X+BX2X-BX3X-BX4X)*CPS + (BZ1X+BZ2X-BZ3X-BZ4X)*SPS;
        D[2][IX] =  BY1X+BY2X-BY3X-BY4X;
        D[3][IX] = (BZ1X+BZ2X-BZ3X-BZ4X)*CPS - (BX1X+BX2X-BX3X-BX4X)*SPS;
      
        D[1][IY] = (BX1Y-BX2Y-BX3Y+BX4Y)*CPS + (BZ1Y-BZ2Y-BZ3Y+BZ4Y)*SPS;
        D[2][IY] =  BY1Y-BY2Y-BY3Y+BY4Y;
        D[3][IY] = (BZ1Y-BZ2Y-BZ3Y+BZ4Y)*CPS - (BX1Y-BX2Y-BX3Y+BX4Y)*SPS;
      
        D[1][IZ] = (BX1Z+BX2Z+BX3Z+BX4Z)*CPS + (BZ1Z+BZ2Z+BZ3Z+BZ4Z)*SPS;
        D[2][IZ] =  BY1Z+BY2Z+BY3Z+BY4Z;
        D[3][IZ] = (BZ1Z+BZ2Z+BZ3Z+BZ4Z)*CPS - (BX1Z+BX2Z+BX3Z+BX4Z)*SPS;
     
        IX += 27;
        IY += 27;
        IZ += 27;
     
        D[1][IX] = SPS*((BX1X+BX2X+BX3X+BX4X)*CPS + (BZ1X+BZ2X+BZ3X+BZ4X)*SPS);
        D[2][IX] = SPS*(BY1X+BY2X+BY3X+BY4X);
        D[3][IX] = SPS*((BZ1X+BZ2X+BZ3X+BZ4X)*CPS - (BX1X+BX2X+BX3X+BX4X)*SPS);
     
        D[1][IY] = SPS*((BX1Y-BX2Y+BX3Y-BX4Y)*CPS + (BZ1Y-BZ2Y+BZ3Y-BZ4Y)*SPS);
        D[2][IY] = SPS*(BY1Y-BY2Y+BY3Y-BY4Y);
        D[3][IY] = SPS*((BZ1Y-BZ2Y+BZ3Y-BZ4Y)*CPS - (BX1Y-BX2Y+BX3Y-BX4Y)*SPS);
     
        D[1][IZ] = SPS*((BX1Z+BX2Z-BX3Z-BX4Z)*CPS + (BZ1Z+BZ2Z-BZ3Z-BZ4Z)*SPS);
        D[2][IZ] = SPS*(BY1Z+BY2Z-BY3Z-BY4Z);
        D[3][IZ] = SPS*((BZ1Z+BZ2Z-BZ3Z-BZ4Z)*CPS - (BX1Z+BX2Z-BX3Z-BX4Z)*SPS);
    }

    for ( I=1; I<=5; I++ ) {

        ZD = ZZ[I+9];
        DIPXYZ_T96( XSM, Y, ZSM-ZD, &BX1X, &BY1X, &BZ1X, &BX1Y, &BY1Y, &BZ1Y, &BX1Z, &BY1Z, &BZ1Z );
        DIPXYZ_T96( XSM, Y, ZSM+ZD, &BX2X, &BY2X, &BZ2X, &BX2Y, &BY2Y, &BZ2Y, &BX2Z, &BY2Z, &BZ2Z );
        IX = 58 + I*2;
        IZ = IX + 1;
        D[1][IX] = (BX1X-BX2X)*CPS + (BZ1X-BZ2X)*SPS;
        D[2][IX] = BY1X-BY2X;
        D[3][IX] = (BZ1X-BZ2X)*CPS - (BX1X-BX2X)*SPS;
 
        D[1][IZ] = (BX1Z+BX2Z)*CPS + (BZ1Z+BZ2Z)*SPS;
        D[2][IZ] = BY1Z+BY2Z;
        D[3][IZ] = (BZ1Z+BZ2Z)*CPS - (BX1Z+BX2Z)*SPS;
 
        IX += 10;
        IZ += 10;
        D[1][IX] = SPS*((BX1X+BX2X)*CPS + (BZ1X+BZ2X)*SPS);
        D[2][IX] = SPS*(BY1X+BY2X);
        D[3][IX] = SPS*((BZ1X+BZ2X)*CPS - (BX1X+BX2X)*SPS);
 
        D[1][IZ] = SPS*((BX1Z-BX2Z)*CPS + (BZ1Z-BZ2Z)*SPS);
        D[2][IZ] = SPS*(BY1Z-BY2Z);
        D[3][IZ] = SPS*((BZ1Z-BZ2Z)*CPS - (BX1Z-BX2Z)*SPS);
    }

    return;
}



/*
 *  - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 *
 *  The 64 linear parameters are amplitudes of the "box" harmonics.
 *  The 16 nonlinear parameters are the scales Pi, and Qk entering the arguments
 *  of sines/cosines and exponents in each of  32 cartesian harmonics
 *  N.A. Tsyganenko, Spring 1994, adjusted for the Birkeland field Aug.22, 1995
 *    Revised  June 12, 1996.
 * - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 */
void BIRK1SHLD_T96( double PS, double X, double Y, double Z, double *BX, double *BY, double *BZ ) {

 
    static double A[] = { -9e99, 1.174198045,-1.463820502,4.840161537,-3.674506864,
                                   82.18368896,-94.94071588,-4122.331796,4670.278676,-21.54975037,
                                   26.72661293,-72.81365728,44.09887902,40.08073706,-51.23563510,
                                   1955.348537,-1940.971550,794.0496433,-982.2441344,1889.837171,
                                   -558.9779727,-1260.543238,1260.063802,-293.5942373,344.7250789,
                                   -773.7002492,957.0094135,-1824.143669,520.7994379,1192.484774,
                                   -1192.184565,89.15537624,-98.52042999,-0.8168777675E-01,
                                   0.4255969908E-01,0.3155237661,-0.3841755213,2.494553332,
                                   -0.6571440817E-01,-2.765661310,0.4331001908,0.1099181537,
                                   -0.6154126980E-01,-0.3258649260,0.6698439193,-5.542735524,
                                   0.1604203535,5.854456934,-0.8323632049,3.732608869,-3.130002153,
                                   107.0972607,-32.28483411,-115.2389298,54.45064360,-0.5826853320,
                                   -3.582482231,-4.046544561,3.311978102,-104.0839563,30.26401293,
                                   97.29109008,-50.62370872,-296.3734955,127.7872523,5.303648988,
                                   10.40368955,69.65230348,466.5099509,1.645049286,3.825838190,
                                   11.66675599,558.9781177,1.826531343,2.066018073,25.40971369,
                                   990.2795225,2.319489258,4.555148484,9.691185703,591.8280358 };

    double  P1[5], R1[5], Q1[5], S1[5], RP[5], RR[5], RQ[5], RS[5];

    int     I, L, M, K, N;
    double  CPS, SPS, S3PS, CYPI, CYQI, SYPI, SYQI, SZRK, CZSK, CZRK, SZSK, SQPR, SQQS, EPR, EQS, HX, HY, HZ;


    P1[1] = A[65];
    R1[1] = A[69];
    Q1[1] = A[73];
    S1[1] = A[77];

    *BX = *BY = *BZ = 0.0;
    CPS  = cos(PS);
    SPS  = sin(PS);
    S3PS = 4.0*CPS*CPS-1.0;

    for ( I=1; I<=4; I++ ) {
        RP[I] = 1.0/P1[I];
        RR[I] = 1.0/R1[I];
        RQ[I] = 1.0/Q1[I];
        RS[I] = 1.0/S1[I];
    }

    L = 0;
    for ( M=1; M<=2; M++ ) { // M=1 IS FOR THE 1ST SUM ("PERP." SYMMETRY) AND M=2 IS FOR THE SECOND SUM ("PARALL." SYMMETRY)

        for ( I=1; I<=4; I++ ) {

            CYPI = cos( Y*RP[I] );
            CYQI = cos( Y*RQ[I] );
            SYPI = sin( Y*RP[I] );
            SYQI = sin( Y*RQ[I] );

            for ( K=1; K<=4; K++ ) {

                SZRK = sin( Z*RR[K] );
                CZSK = cos( Z*RS[K] );
                CZRK = cos( Z*RR[K] );
                SZSK = sin( Z*RS[K] );
                SQPR = sqrt( RP[I]*RP[I] + RR[K]*RR[K] );
                SQQS = sqrt( RQ[I]*RQ[I] + RS[K]*RS[K] );
                EPR  = exp( X*SQPR );
                EQS  = exp( X*SQQS );

                for ( N=1; N<=2; N++ ) { // N=1 IS FOR THE FIRST PART OF EACH COEFFICIENT AND N=2 IS FOR THE SECOND ONE

                    if ( M == 1 ) {
                        if (N == 1) {
                            HX = -SQPR*EPR*CYPI*SZRK;
                            HY =  RP[I]*EPR*SYPI*SZRK;
                            HZ = -RR[K]*EPR*CYPI*CZRK;
                        } else {
                            HX *= CPS;
                            HY *= CPS;
                            HZ *= CPS;
                        }
                    } else {
                        if (N == 1) {
                            HX = -SPS*SQQS*EQS*CYQI*CZSK;
                            HY =  SPS*RQ[I]*EQS*SYQI*CZSK;
                            HZ =  SPS*RS[K]*EQS*CYQI*SZSK;
                        } else {
                            HX *= S3PS;
                            HY *= S3PS;
                            HZ *= S3PS;
                       }
                    }

                    ++L;
                    *BX += A[L]*HX;
                    *BY += A[L]*HY;
                    *BZ += A[L]*HZ;

                } // end N loop
            } // end K loop
        } // end I loop
    } // end M loop

    return;

}




void BIRK2TOT_02_T96( double PS, double X, double Y, double Z,double *BX, double *BY,double *BZ, LgmTsyg1996_Info *tInfo ) {

    double  WX, WY, WZ, HX, HY, HZ;

    BIRK2SHL_T96( X, Y, Z, PS, &WX, &WY, &WZ );
    R2_BIRK_T96( X, Y, Z, PS, &HX, &HY, &HZ );

    *BX = WX + HX;
    *BY = WY + HY;
    *BZ = WZ + HZ;

    return;

}





/*
 * THIS CODE IS FOR THE FIELD FROM  2x2x2=8 "CARTESIAN" HARMONICS
 *
 *
 * - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 *    The model parameters are provided to this module via common-block /A/.
 *  The 16 linear parameters enter in pairs in the amplitudes of the
 *       "cartesian" harmonics.
 *    The 8 nonlinear parameters are the scales Pi,Ri,Qi,and Si entering the
 *  arguments of exponents, sines, and cosines in each of the 8 "Cartesian"
 *   harmonics
 * - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 */
void BIRK2SHL_T96( double X, double Y, double Z, double PS, double *HX, double *HY, double *HZ ) {

    double  P[3], R[3], Q[3], S[3];

    static double A[] = { -9e99, -111.6371348,124.5402702,110.3735178,-122.0095905,
                               111.9448247,-129.1957743,-110.7586562,126.5649012,-0.7865034384,
                               -0.2483462721,0.8026023894,0.2531397188,10.72890902,0.8483902118,
                               -10.96884315,-0.8583297219,13.85650567,14.90554500,10.21914434,
                               10.09021632,6.340382460,14.40432686,12.71023437,12.83966657 };

    int     L, M, I, K, N;
    double  SPS, CPS, S3PS, SYPI, CYPI, SYQI, CYQI, SZRK, CZSK, SZSK;
    double  CZRK, SQPR, SQQS, EPR, EQS, DX, DY, DZ;

    P[1] = A[17];
    R[1] = A[19];
    Q[1] = A[21];
    S[1] = A[23];

    SPS  = sin(PS);
    CPS  = cos(PS);
    S3PS = 4.0*CPS*CPS-1.0; //   THIS IS SIN(3*PS)/SIN(PS)

    *HX = *HY= *HZ = 0.0; 
    L = 0;

    for ( M=1; M<=2; M++ ) { // M=1 IS FOR THE 1ST SUM ("PERP." SYMMETRY) AND M=2 IS FOR THE SECOND SUM ("PARALL." SYMMETRY)

        for ( I=1; I<=2; I++ ) {

            SYPI = sin( Y/P[I] );
            CYPI = cos( Y/P[I] );

            SYQI = sin( Y/Q[I] );
            CYQI = cos( Y/Q[I] );

            for ( K=1; K<=2; K++ ) {

                SZRK = sin( Z/R[K] );
                CZSK = cos( Z/S[K] );

                SZSK = sin( Z/S[K] );
                CZRK = cos( Z/R[K] );

                SQPR = sqrt( 1.0/(P[I]*P[I]) + 1.0/(R[K]*R[K]) );
                SQQS = sqrt( 1.0/(Q[I]*Q[I]) + 1.0/(S[K]*S[K]) );
                EPR  = exp( X*SQPR );
                EQS  = exp( X*SQQS );

                for ( N=1; N<=2; N++ ) { // N=1 IS FOR THE FIRST PART OF EACH COEFFICIENT AND N=2 IS FOR THE SECOND ONE

                    ++L;
                    if (M == 1) {
                        if (N == 1) {
                            DX = -SQPR*EPR*CYPI*SZRK;
                            DY =  EPR/P[I]*SYPI*SZRK;
                            DZ = -EPR/R[K]*CYPI*CZRK;
                            *HX += A[L]*DX;
                            *HY += A[L]*DY;
                            *HZ += A[L]*DZ;
                         } else {
                            DX *= CPS;
                            DY *= CPS;
                            DZ *= CPS;
                            *HX += A[L]*DX;
                            *HY += A[L]*DY;
                            *HZ += A[L]*DZ;
                        }
                     } else {
                        if (N == 1) {
                            DX = -SPS*SQQS*EQS*CYQI*CZSK;
                            DY =  SPS*EQS/Q[I]*SYQI*CZSK;
                            DZ =  SPS*EQS/S[K]*CYQI*SZSK;
                            *HX += A[L]*DX;
                            *HY += A[L]*DY;
                            *HZ += A[L]*DZ;
                       } else {
                            DX *= S3PS;
                            DY *= S3PS;
                            DZ *= S3PS;
                            *HX += A[L]*DX;
                            *HY += A[L]*DY;
                            *HZ += A[L]*DZ;
    
                        }
                    }

                } // end N loop
            } // end K loop
        } // end I loop
    } // end M loop

    return;

}






/*
 *  RETURNS THE MODEL FIELD FOR THE REGION 2 BIRKELAND CURRENT/PARTIAL RC
 *    (WITHOUT SHIELDING FIELD)
 */
void R2_BIRK_T96( double X, double Y, double Z, double PS, double *BX, double *BY, double *BZ ) {
//!!!!SAVE PSI,CPS,SPS
    static double   DELARG  = 0.030;
    static double   DELARG1 = 0.015;
    static double   PSI     = 10.0;
    double          CPS, SPS, XSM, ZSM, BXSM, BZSM, F2, F1, BXSM1, BY1, BZSM1, BXSM2, BY2, BZSM2, XKS;

    if ( fabs(PSI-PS) > 1.e-10) {
         PSI = PS;
         CPS = cos(PS);
         SPS = sin(PS);
    }

    XSM = X*CPS - Z*SPS;
    ZSM = Z*CPS + X*SPS;

    XKS = XKSI_T96( XSM, Y, ZSM );
    if ( XKS < -(DELARG+DELARG1) ) {
        R2OUTER( XSM, Y, ZSM, &BXSM, BY, &BZSM );
        BXSM = -BXSM*0.02; //  ALL COMPONENTS ARE MULTIPLIED BY THE
        *BY   = -(*BY)*0.02;   //  FACTOR -0.02, IN ORDER TO NORMALIZE THE
        BZSM = -BZSM*0.02; //  FIELD (SO THAT Bz=-1 nT at X=-5.3 RE, Y=Z=0)
    }


    if ( (XKS >= -(DELARG+DELARG1)) && (XKS < (-DELARG+DELARG1)) ) {
        R2OUTER( XSM, Y, ZSM, &BXSM1, &BY1, &BZSM1 );
        R2SHEET_T96( XSM, Y, ZSM, &BXSM2, &BY2, &BZSM2 );
        F2 = -0.02*TKSI_T96( XKS, -DELARG, DELARG1 );
        F1 = -0.02 - F2;
        BXSM = BXSM1*F1 + BXSM2*F2;
        *BY   = BY1*F1   + BY2*F2;
        BZSM = BZSM1*F1 + BZSM2*F2;
    }

    if ( (XKS >= (-DELARG+DELARG1)) && (XKS < (DELARG-DELARG1)) ) {
        R2SHEET_T96( XSM, Y, ZSM, &BXSM, BY, &BZSM );
        BXSM = -BXSM*0.02;
        *BY   = -(*BY)*0.02;
        BZSM = -BZSM*0.02;
    }

    if ( (XKS >= (DELARG-DELARG1)) && (XKS < (DELARG+DELARG1)) ) {
        R2INNER_T96( XSM, Y, ZSM, &BXSM1, &BY1, &BZSM1 );
        R2SHEET_T96( XSM, Y, ZSM, &BXSM2, &BY2, &BZSM2 );
        F1 = -0.02*TKSI_T96( XKS, DELARG, DELARG1 );
        F2 = -0.02 - F1;
        BXSM = BXSM1*F1 + BXSM2*F2;
        *BY   = BY1*F1   + BY2*F2;
        BZSM = BZSM1*F1 + BZSM2*F2;
    }

    if ( XKS >= (DELARG+DELARG1) ) {
        R2INNER_T96( XSM, Y, ZSM, &BXSM, BY, &BZSM );
        BXSM = -BXSM*0.02;
        *BY   = -(*BY)*0.02;
        BZSM = -BZSM*0.02;
    }

    *BX = BXSM*CPS + BZSM*SPS;
    *BZ = BZSM*CPS - BXSM*SPS;

    return;
}



void  R2INNER_T96( double X, double Y, double Z, double *BX, double *BY, double *BZ ) {

    static double PL[] = { -9e99, 154.185, -2.12446, 0.601735E-01, -0.153954E-02, 0.355077E-04, 29.9996, 262.886, 99.9132 };
    static double PN[] = { -9e99, -8.1902, 6.5239, 5.504, 7.7815, 0.8573, 3.0986, 0.0774, -0.038 };
    double  CBX[6], CBY[6], CBZ[6];
    double  DBX8, DBY8, DBZ8, DBX6, DBY6, DBZ6, DBX7, DBY7, DBZ7;

    BCONIC_T96( X, Y, Z, CBX, CBY, CBZ, 5 );

    /*
     *   NOW INTRODUCE  ONE  4-LOOP SYSTEM:
     */
    LOOPS4_T96( X, Y, Z, &DBX8, &DBY8, &DBZ8, PN[1], PN[2], PN[3], PN[4], PN[5], PN[6] );
    DIPDISTR_T96( X-PN[7], Y, Z, &DBX6, &DBY6, &DBZ6, 0);
    DIPDISTR_T96( X-PN[8], Y, Z, &DBX7, &DBY7, &DBZ7, 1);

    // NOW COMPUTE THE FIELD COMPONENTS:
    *BX = PL[1]*CBX[1] + PL[2]*CBX[2] + PL[3]*CBX[3] + PL[4]*CBX[4] + PL[5]*CBX[5] + PL[6]*DBX6 + PL[7]*DBX7 + PL[8]*DBX8;
    *BY = PL[1]*CBY[1] + PL[2]*CBY[2] + PL[3]*CBY[3] + PL[4]*CBY[4] + PL[5]*CBY[5] + PL[6]*DBY6 + PL[7]*DBY7 + PL[8]*DBY8;
    *BZ = PL[1]*CBZ[1] + PL[2]*CBZ[2] + PL[3]*CBZ[3] + PL[4]*CBZ[4] + PL[5]*CBZ[5] + PL[6]*DBZ6 + PL[7]*DBZ7 + PL[8]*DBZ8;
    
    return;
}





/*
 *   "CONICAL" HARMONICS
 */
void  BCONIC_T96( double X, double Y, double Z, double CBX[], double CBY[], double CBZ[], int NMAX ) {

    int     M;
    double RO2, RO, CF, SF, CFM1, SFM1, R2, R, C, S, CH, SH, TNHM1, CNHM1, TNH;
    double CNH, CFM, SFM, TNHM, CNHM, BT, BF;


    RO2  = X*X + Y*Y;
    RO   = sqrt(RO2);
 
    CF   = X/RO;
    SF   = Y/RO;
    CFM1 = 1.0;
    SFM1 = 0.0;
 
    R2    = RO2 + Z*Z;
    R     = sqrt(R2);
    C     = Z/R;
    S     = RO/R;
    CH    = sqrt(0.5*(1.0+C));
    SH    = sqrt(0.5*(1.0-C));
    TNHM1 = 1.0;
    CNHM1 = 1.0;
    TNH   = SH/CH;
    CNH   = 1.0/TNH;
 
    for ( M=1; M<=NMAX; M++ ){
        CFM    = CFM1*CF - SFM1*SF;
        SFM    = CFM1*SF + SFM1*CF;
        CFM1   = CFM;
        SFM1   = SFM;
        TNHM   = TNHM1*TNH;
        CNHM   = CNHM1*CNH;
        BT     = M*CFM/(R*S)*(TNHM+CNHM);
        BF     = -0.5*M*SFM/R*(TNHM1/(CH*CH) - CNHM1/(SH*SH));
        TNHM1  = TNHM;
        CNHM1  = CNHM;
        CBX[M] =  BT*C*CF - BF*SF;
        CBY[M] =  BT*C*SF + BF*CF;
        CBZ[M] = -BT*S;
    }

    return;
}




/*
 *   RETURNS FIELD COMPONENTS FROM A LINEAR DISTRIBUTION OF DIPOLAR SOURCES
 *     ON THE Z-AXIS.  THE PARAMETER MODE DEFINES HOW THE DIPOLE_T96 STRENGTH
 *     VARIES ALONG THE Z-AXIS:  MODE=0 IS FOR A STEP-FUNCTION (Mx=const > 0
 *         FOR Z > 0, AND Mx=-const < 0 FOR Z < 0)
 *      WHILE MODE=1 IS FOR A LINEAR VARIATION OF THE DIPOLE_T96 MOMENT DENSITY
 *       SEE NB#3, PAGE 53 FOR DETAILS.
 *
 *
 * INPUT: X,Y,Z OF A POINT OF SPACE, AND MODE
 */
void  DIPDISTR_T96( double X, double Y, double Z, double *BX, double *BY, double *BZ, int MODE ) {

    double  X2, Y2, RHO2, RHO4, R2, R3;

    X2   = X*X; Y2 = Y*Y; 
    RHO2 = X2 + Y2;
    RHO4 = RHO2*RHO2;
    R2   = RHO2 + Z*Z;
    R3   = R2*sqrt(R2);

    if ( MODE == 0 ) {
        *BX = Z/RHO4*(R2*(Y2-X2)-RHO2*X2)/R3;
        *BY = -X*Y*Z/RHO4*(2.0*R2+RHO2)/R3;
        *BZ = X/R3;
    } else {
        *BX = Z/RHO4*(Y2-X2);
        *BY = -2.0*X*Y*Z/RHO4;
        *BZ = X/RHO2;
    }
    return;
}




void R2OUTER_T9_T966 ( double X, double Y, double Z, double *BX, double *BY, double *BZ ) {

    double DBX1, DBY1, DBZ1, DBX2, DBY2, DBZ2, DBX3, DBY3, DBZ3, DBX4, DBY4, DBZ4, DBX5, DBY5, DBZ5;
    
    static double PL[] = { -9e99, -34.105, -2.00019, 628.639, 73.4847, 12.5162 };
    static double PN[] = { -9e99, 0.55, 0.694, 0.0031, 1.55, 2.8, 0.1375, -0.7, 0.2, 0.9625, -2.994, 2.925, -1.775, 4.3, -0.275, 2.7, 0.4312, 1.55 };

    /* 
     * THREE PAIRS OF CROSSED LOOPS:
     */
    CROSSLP_T96( X, Y, Z, &DBX1, &DBY1, &DBZ1, PN[1], PN[2], PN[3] );
    CROSSLP_T96( X, Y, Z, &DBX2, &DBY2, &DBZ2, PN[4], PN[5], PN[6] );
    CROSSLP_T96( X, Y, Z, &DBX3, &DBY3, &DBZ3, PN[7], PN[8], PN[9] );

    /*
     *  NOW AN EQUATORIAL LOOP ON THE NIGHTSIDE
     */
    CIRCLE_T96( X-PN[10], Y, Z, PN[11], &DBX4, &DBY4, &DBZ4 );

    /*
     *  NOW A 4-LOOP SYSTEM ON THE NIGHTSIDE
     */
    LOOPS4_T96( X, Y, Z, &DBX5, &DBY5, &DBZ5, PN[12], PN[13], PN[14], PN[15], PN[16], PN[17] );

    /*
     *  NOW COMPUTE THE FIELD COMPONENTS:
     */
    *BX = PL[1]*DBX1 + PL[2]*DBX2 + PL[3]*DBX3 + PL[4]*DBX4 + PL[5]*DBX5;
    *BY = PL[1]*DBY1 + PL[2]*DBY2 + PL[3]*DBY3 + PL[4]*DBY4 + PL[5]*DBY5;
    *BZ = PL[1]*DBZ1 + PL[2]*DBZ2 + PL[3]*DBZ3 + PL[4]*DBZ4 + PL[5]*DBZ5;

    return;

}



/*
 *   RETURNS FIELD COMPONENTS FROM A SYSTEM OF 4 CURRENT LOOPS, POSITIONED
 *     SYMMETRICALLY WITH RESPECT TO NOON-MIDNIGHT MERIDIAN AND EQUATORIAL
 *      PLANES.
 *  INPUT: X,Y,Z OF A POINT OF SPACE
 *        XC,YC,ZC (YC > 0 AND ZC > 0) - POSITION OF THE CENTER OF THE
 *                                         1ST-QUADRANT LOOP
 *        R - LOOP RADIUS (THE SAME FOR ALL FOUR)
 *        THETA, PHI  -  SPECIFY THE ORIENTATION OF THE NORMAL OF THE 1ST LOOP
 */
void LOOPS4_T96( double X, double Y, double Z, double *BX, double *BY, double *BZ, double XC, double YC, double ZC, double R, double THETA, double PHI ) {

    double  CT, ST, CP, SP, XS, YSS, ZS, XSS, ZSS, BXSS, BYS, BZSS, BXS;
    double  BZ1, BX1, BY1, BZ2, BX2, BY2, BZ3, BX3, BY3, BZ4, BY4, BX4;

    CT = cos(THETA); ST = sin(THETA);
    CP = cos(PHI);   SP = sin(PHI);

    //------------------------------------1ST QUADRANT:
    XS  = (X-XC)*CP + (Y-YC)*SP;
    YSS = (Y-YC)*CP - (X-XC)*SP;
    ZS  = Z-ZC;
    XSS = XS*CT - ZS*ST;
    ZSS = ZS*CT + XS*ST;

    CIRCLE_T96( XSS, YSS, ZSS, R, &BXSS, &BYS, &BZSS );
    BXS = BXSS*CT + BZSS*ST;
    BZ1 = BZSS*CT - BXSS*ST;
    BX1 = BXS*CP - BYS*SP;
    BY1 = BXS*SP + BYS*CP;

    //-------------------------------------2nd QUADRANT:
    XS  = (X-XC)*CP - (Y+YC)*SP;
    YSS = (Y+YC)*CP + (X-XC)*SP;
    ZS  = Z-ZC;
    XSS = XS*CT - ZS*ST;
    ZSS = ZS*CT + XS*ST;

    CIRCLE_T96( XSS, YSS, ZSS, R, &BXSS, &BYS, &BZSS );
    BXS = BXSS*CT + BZSS*ST;
    BZ2 = BZSS*CT - BXSS*ST;
    BX2 = BXS*CP  + BYS*SP;
    BY2 = -BXS*SP + BYS*CP;

    //-------------------------------------3RD QUADRANT:
    XS  = -(X-XC)*CP + (Y+YC)*SP;
    YSS = -(Y+YC)*CP - (X-XC)*SP;
    ZS  = Z+ZC;
    XSS = XS*CT - ZS*ST;
    ZSS = ZS*CT + XS*ST;

    CIRCLE_T96( XSS, YSS, ZSS, R, &BXSS, &BYS, &BZSS );
    BXS = BXSS*CT + BZSS*ST;
    BZ3 = BZSS*CT - BXSS*ST;
    BX3 = -BXS*CP - BYS*SP;
    BY3 = BXS*SP - BYS*CP;

    //-------------------------------------4TH QUADRANT:
    XS  = -(X-XC)*CP - (Y-YC)*SP;
    YSS = -(Y-YC)*CP + (X-XC)*SP;
    ZS  = Z+ZC;
    XSS = XS*CT - ZS*ST;
    ZSS = ZS*CT + XS*ST;

    CIRCLE_T96( XSS, YSS, ZSS, R, &BXSS, &BYS, &BZSS );
    BXS = BXSS*CT + BZSS*ST;
    BZ4 = BZSS*CT - BXSS*ST;
    BX4 = -BXS*CP + BYS*SP;
    BY4 = -BXS*SP - BYS*CP;

    *BX = BX1 + BX2 + BX3 + BX4;
    *BY = BY1 + BY2 + BY3 + BY4;
    *BZ = BZ1 + BZ2 + BZ3 + BZ4;

    return;

}










void R2SHEET_T96( double X, double Y, double Z, double *BX, double *BY, double *BZ ) {

    static double PNONX[] = { -9e99, -19.0969, -9.28828, -0.129687, 5.58594, 22.5055, 0.483750e-01, 0.396953e-01, 0.579023e-01};
    static double PNONY[] = { -9e99, -13.6750, -6.70625, 2.31875, 11.4062, 20.4562, 0.478750e-01, 0.363750e-01, 0.567500e-01 };
    static double PNONZ[] = { -9e99, -16.7125, -16.4625, -0.1625, 5.1, 23.7125, 0.355625e-01, 0.318750e-01, 0.538750e-01 };

    static double A[] = { -9e99, 8.07190, -7.39582, -7.62341, 0.684671, -13.5672, 11.6681, 
                                 13.1154, -0.890217, 7.78726, -5.38346, -8.08738, 0.609385, 
                                 -2.70410,  3.53741, 3.15549, -1.11069, -8.47555, 0.278122, 
                                 2.73514, 4.55625, 13.1134, 1.15848, -3.52648, -8.24698, 
                                 -6.85710, -2.81369,  2.03795,  4.64383, 2.49309, -1.22041, 
                                 -1.67432, -0.422526, -5.39796, 7.10326, 5.53730, -13.1918, 
                                 4.67853, -7.60329, -2.53066,  7.76338,  5.60165, 5.34816, 
                                 -4.56441, 7.05976, -2.62723, -0.529078, 1.42019, -2.93919, 
                                 55.6338, -1.55181, 39.8311, -80.6561, -46.9655, 32.8925, 
                                 -6.32296, 19.7841, 124.731, 10.4347, -30.7581, 102.680, 
                                 -47.4037, -3.31278, 9.37141, -50.0268, -533.319, 110.426, 
                                 1000.20, -1051.40,  1619.48, 589.855, -1462.73, 1087.10, 
                                 -1994.73, -1654.12, 1263.33, -260.210, 1424.84, 1255.71, 
                                 -956.733,  219.946 } ;

    static double B[] = { -9e99, -9.08427, 10.6777, 10.3288, -0.969987, 6.45257, -8.42508, 
                                  -7.97464, 1.41996, -1.92490, 3.93575, 2.83283, -1.48621, 
                                  0.244033, -0.757941, -0.386557, 0.344566, 9.56674, -2.5365, 
                                  -3.32916, -5.86712, -6.19625, 1.83879, 2.52772, 4.34417, 
                                  1.87268, -2.13213, -1.69134, -.176379, -.261359, .566419, 
                                  0.3138, -0.134699, -3.83086, -8.4154, 4.77005, -9.31479, 
                                  37.5715, 19.3992, -17.9582, 36.4604, -14.9993, -3.1442, 
                                  6.17409, -15.5519, 2.28621, -0.891549e-2, -.462912, 2.47314, 
                                  41.7555, 208.614, -45.7861, -77.8687, 239.357, -67.9226, 
                                  66.8743, 238.534, -112.136, 16.2069, -40.4706, -134.328, 
                                  21.56, -0.201725, 2.21, 32.5855, -108.217, -1005.98, 
                                  585.753, 323.668, -817.056, 235.750, -560.965, -576.892, 
                                  684.193, 85.0275, 168.394, 477.776, -289.253, -123.216, 
                                  75.6501, -178.605 };
    
    static double C[] = { -9e99,  1167.61, -917.782, -1253.2, -274.128, -1538.75, 1257.62, 
                                  1745.07, 113.479, 393.326, -426.858, -641.1, 190.833, 
                                  -29.9435, -1.04881, 117.125, -25.7663, -1168.16, 910.247, 
                                  1239.31, 289.515, 1540.56, -1248.29, -1727.61, -131.785, 
                                  -394.577, 426.163, 637.422, -187.965, 30.0348, 0.221898, 
                                  -116.68, 26.0291, 12.6804, 4.84091, 1.18166, -2.75946, 
                                  -17.9822, -6.80357, -1.47134, 3.02266, 4.79648, 0.665255, 
                                  -0.256229, -0.857282e-1, -0.588997, 0.634812e-1, 0.164303, 
                                  -0.15285, 22.2524, -22.4376, -3.85595, 6.07625, -105.959, 
                                  -41.6698, 0.378615, 1.55958, 44.3981, 18.8521, 3.19466, 
                                   5.89142, -8.63227, -2.36418, -1.027, -2.31515, 1035.38, 
                                   2040.66, -131.881, -744.533, -3274.93, -4845.61, 482.438, 
                                  1567.43, 1354.02, 2040.47, -151.653, -845.012, -111.723, 
                                  -265.343, -26.1171, 216.632 };
    double  XKS, q2, T1X, g, g2 , g3, h, h3, T2X, PNONX8, g4, h2, h5, T3X, T2Y, T3Y, T3Z, T2Z;
    double  RHO2, R, RHO, C1P, S1P, S2P, C2P, S3P, C3P, S4P, CT, ST, S1, S2, S3, S4, S5, T1Y, T1Z;

 
    XKS = XKSI_T96( X, Y, Z );    //  variation across the current sheet
    q2 = XKS*XKS;

    T1X = XKS/sqrt( q2 + PNONX[6]*PNONX[6] );
    g = PNONX[7]; g2 = g*g; g3 = g*g2;
    h = sqrt( q2 + g2 ); h3 = h*h*h;
    T2X = g3/h3;
    g = PNONX[8]; g2 = g*g; g4 = g2*g2;
    h = sqrt( q2 + g2 ); h2 = h*h; h3 = h2*h; h5 = h3*h2;
    T3X = XKS/h5 * 3.4938560*g4;
 
    T1Y = XKS/sqrt( q2 + PNONY[6]*PNONY[6] );
    g = PNONY[7]; g2 = g*g; g3 = g2*g;
    h = sqrt( q2 + g2 ); h3 = h*h*h;
    T2Y = g3/h3;
    g = PNONY[8]; g2 = g*g; g4=g2*g2;
    h = sqrt( q2 + g2 ); h2 = h*h; h3 = h*h2; h5 = h3*h2;
    T3Y = XKS/h5 * 3.4938560*g4;
 
    
    T1Z = XKS/sqrt( q2 + PNONZ[6]*PNONZ[6]);
    g = PNONZ[7]; g2 = g*g; g3 = g2*g;
    h = sqrt( q2 + g2 ); h3 = h*h*h;
    T2Z = g3/h3;
    g = PNONZ[8]; g2 = g*g; g4=g2*g2;
    h = sqrt( q2 + g2 ); h2 = h*h; h3 = h*h2; h5 = h3*h2;
    T3Z = XKS/h5 * 3.4938560*g4;
 
    RHO2 = X*X + Y*Y;
    R    = sqrt( RHO2 + Z*Z );
    RHO  = sqrt( RHO2 );
 
    C1P = X/RHO;
    S1P = Y/RHO;
    S2P = 2.0*S1P*C1P;
    C2P = C1P*C1P-S1P*S1P;
    S3P = S2P*C1P+C2P*S1P;
    C3P = C2P*C1P-S2P*S1P;
    S4P = S3P*C1P+C3P*S1P;
    CT  = Z/R;
    ST  = RHO/R;


    /*
     *                   NOW COMPUTE THE GSM FIELD COMPONENTS:
     *
     */
    S1 = FEXP_T96( CT, PNONX[1] ); S2 = FEXP_T96( CT, PNONX[2] ); S3 = FEXP_T96( CT, PNONX[3] ); S4 = FEXP_T96( CT, PNONX[4] ); S5 = FEXP_T96( CT, PNONX[5] );
    *BX =   S1*((A[1]+A[2]*T1X+A[3]*T2X+A[4]*T3X)     + C1P*(A[5]+A[6]*T1X+A[7]*T2X+A[8]*T3X)     + C2P*(A[9]+A[10]*T1X+A[11]*T2X+A[12]*T3X)  + C3P*(A[13]+A[14]*T1X+A[15]*T2X+A[16]*T3X))
          + S2*((A[17]+A[18]*T1X+A[19]*T2X+A[20]*T3X) + C1P*(A[21]+A[22]*T1X+A[23]*T2X+A[24]*T3X) + C2P*(A[25]+A[26]*T1X+A[27]*T2X+A[28]*T3X) + C3P*(A[29]+A[30]*T1X+A[31]*T2X+A[32]*T3X))
          + S3*((A[33]+A[34]*T1X+A[35]*T2X+A[36]*T3X) + C1P*(A[37]+A[38]*T1X+A[39]*T2X+A[40]*T3X) + C2P*(A[41]+A[42]*T1X+A[43]*T2X+A[44]*T3X) + C3P*(A[45]+A[46]*T1X+A[47]*T2X+A[48]*T3X))
          + S4*((A[49]+A[50]*T1X+A[51]*T2X+A[52]*T3X) + C1P*(A[53]+A[54]*T1X+A[55]*T2X+A[56]*T3X) + C2P*(A[57]+A[58]*T1X+A[59]*T2X+A[60]*T3X) + C3P*(A[61]+A[62]*T1X+A[63]*T2X+A[64]*T3X))
          + S5*((A[65]+A[66]*T1X+A[67]*T2X+A[68]*T3X) + C1P*(A[69]+A[70]*T1X+A[71]*T2X+A[72]*T3X) + C2P*(A[73]+A[74]*T1X+A[75]*T2X+A[76]*T3X) + C3P*(A[77]+A[78]*T1X+A[79]*T2X+A[80]*T3X));

    S1 = FEXP_T96( CT, PNONY[1] ); S2 = FEXP_T96( CT, PNONY[2] ); S3 = FEXP_T96( CT, PNONY[3] ); S4 = FEXP_T96( CT, PNONY[4] ); S5 = FEXP_T96( CT, PNONY[5] );
    *BY =    S1*(S1P*(B[1] +B[2]*T1Y+ B[3]*T2Y+ B[4]*T3Y)  + S2P*(B[5]+ B[6]*T1Y +B[7]*T2Y +B[8]*T3Y)  + S3P*(B[9] +B[10]*T1Y+B[11]*T2Y+B[12]*T3Y) + S4P*(B[13]+B[14]*T1Y+B[15]*T2Y+B[16]*T3Y))
           + S2*(S1P*(B[17]+B[18]*T1Y+B[19]*T2Y+B[20]*T3Y) + S2P*(B[21]+B[22]*T1Y+B[23]*T2Y+B[24]*T3Y) + S3P*(B[25]+B[26]*T1Y+B[27]*T2Y+B[28]*T3Y) + S4P*(B[29]+B[30]*T1Y+B[31]*T2Y+B[32]*T3Y))
           + S3*(S1P*(B[33]+B[34]*T1Y+B[35]*T2Y+B[36]*T3Y) + S2P*(B[37]+B[38]*T1Y+B[39]*T2Y+B[40]*T3Y) + S3P*(B[41]+B[42]*T1Y+B[43]*T2Y+B[44]*T3Y) + S4P*(B[45]+B[46]*T1Y+B[47]*T2Y+B[48]*T3Y))
           + S4*(S1P*(B[49]+B[50]*T1Y+B[51]*T2Y+B[52]*T3Y) + S2P*(B[53]+B[54]*T1Y+B[55]*T2Y+B[56]*T3Y) + S3P*(B[57]+B[58]*T1Y+B[59]*T2Y+B[60]*T3Y) + S4P*(B[61]+B[62]*T1Y+B[63]*T2Y+B[64]*T3Y))
           + S5*(S1P*(B[65]+B[66]*T1Y+B[67]*T2Y+B[68]*T3Y) + S2P*(B[69]+B[70]*T1Y+B[71]*T2Y+B[72]*T3Y) + S3P*(B[73]+B[74]*T1Y+B[75]*T2Y+B[76]*T3Y) + S4P*(B[77]+B[78]*T1Y+B[79]*T2Y+B[80]*T3Y));

    S1 = FEXP1_T96( CT, PNONZ[1] ); S2 = FEXP1_T96( CT, PNONZ[2] ); S3 = FEXP1_T96( CT, PNONZ[3] ); S4 = FEXP1_T96( CT, PNONZ[4] ); S5 = FEXP1_T96( CT, PNONZ[5] );
    *BZ =   S1*((C[1] +C[2]*T1Z +C[3]*T2Z +C[4]*T3Z)  + C1P*(C[5] +C[6]*T1Z +C[7]*T2Z +C[8]*T3Z)  + C2P*(C[9] +C[10]*T1Z+C[11]*T2Z+C[12]*T3Z) + C3P*(C[13]+C[14]*T1Z+C[15]*T2Z+C[16]*T3Z))
          + S2*((C[17]+C[18]*T1Z+C[19]*T2Z+C[20]*T3Z) + C1P*(C[21]+C[22]*T1Z+C[23]*T2Z+C[24]*T3Z) + C2P*(C[25]+C[26]*T1Z+C[27]*T2Z+C[28]*T3Z) + C3P*(C[29]+C[30]*T1Z+C[31]*T2Z+C[32]*T3Z))
          + S3*((C[33]+C[34]*T1Z+C[35]*T2Z+C[36]*T3Z) + C1P*(C[37]+C[38]*T1Z+C[39]*T2Z+C[40]*T3Z) + C2P*(C[41]+C[42]*T1Z+C[43]*T2Z+C[44]*T3Z) + C3P*(C[45]+C[46]*T1Z+C[47]*T2Z+C[48]*T3Z))
          + S4*((C[49]+C[50]*T1Z+C[51]*T2Z+C[52]*T3Z) + C1P*(C[53]+C[54]*T1Z+C[55]*T2Z+C[56]*T3Z) + C2P*(C[57]+C[58]*T1Z+C[59]*T2Z+C[60]*T3Z) + C3P*(C[61]+C[62]*T1Z+C[63]*T2Z+C[64]*T3Z))
          + S5*((C[65]+C[66]*T1Z+C[67]*T2Z+C[68]*T3Z) + C1P*(C[69]+C[70]*T1Z+C[71]*T2Z+C[72]*T3Z) + C2P*(C[73]+C[74]*T1Z+C[75]*T2Z+C[76]*T3Z) + C3P*(C[77]+C[78]*T1Z+C[79]*T2Z+C[80]*T3Z));

    return;

}





double XKSI_T96( double X, double Y, double Z ) {

    double  A11A12, A21A22, A41A42, A51A52, A61A62, B11B12, B21B22, C61C62, C71C72, R0, DR, TNOON, DTETA;
    double  DR2, X2, Y2, Z2, XY, XYZ, R2, R, R3, R4, XR, YR, ZR;
    double  PR, F, G, H, G2, FGH, FGH32, FCHSG2, SQFCHSG2, ALPHA, THETA, PHI, g, g2;


    /*
     *   A11 - C72, R0, and DR below  ARE STRETCH PARAMETERS (P.26-27, NB# 3),
     */
    A11A12 = 0.305662;
    A21A22 = -0.383593;
    A41A42 = 0.2677733;
    A51A52 = -0.097656;
    A61A62 = -0.636034;
    B11B12 = -0.359862;
    B21B22 = 0.424706;
    C61C62 = -0.126366;
    C71C72 = 0.292578;
    R0     = 1.21563;
    DR     = 7.50937;

    TNOON = 0.3665191;  // Correspond to noon and midnight latitudes 69 and 63.5 degs, resp.
    DTETA = 0.09599309; // Correspond to noon and midnight latitudes 69 and 63.5 degs, resp.

    DR2 = DR*DR;

    X2  = X*X;
    Y2  = Y*Y;
    Z2  = Z*Z;
    XY  = X*Y;
    XYZ = XY*Z;
    R2  = X2 + Y2 + Z2;
    R   = sqrt(R2);
    R3  = R2*R;
    R4  = R2*R2;
    XR  = X/R;
    YR  = Y/R;
    ZR  = Z/R;


    g = R-R0; g2 = g*g;
    PR = ( R < R0 ) ? 0.0 : sqrt( g2 + DR2 ) - DR;


    F  = X + PR*(A11A12 + A21A22*XR + A41A42*XR*XR + A51A52*YR*YR + A61A62*ZR*ZR);
    G  = Y + PR*(B11B12*YR + B21B22*XR*YR);
    H  = Z + PR*(C61C62*ZR + C71C72*XR*ZR);
    G2 = G*G;

    FGH    = F*F + G2 + H*H;
    g = sqrt( FGH );
    FGH32  = g*g*g;
    FCHSG2 = F*F + G2;


    if ( FCHSG2 < 1.0e-5 ) return( -1.0 ); // THIS IS JUST FOR ELIMINATING PROBLEMS ON THE Z-AXIS

    SQFCHSG2 = sqrt(FCHSG2);
    ALPHA    = FCHSG2/FGH32;
    THETA    = TNOON + 0.5*DTETA*(1.0 - F/SQFCHSG2);
    g = sin(THETA);
    PHI      = g*g;
 
    return( ALPHA-PHI );

}



double FEXP_T96( double S, double A ) {
    if ( A <  0.0 ) return( sqrt(-2.0*A*M_E)*S*exp(A*S*S) );
    if ( A >= 0.0 ) return( S*exp(A*(S*S-1.0)) );
}


double  FEXP_T961_T96( double S, double A ) {
    if ( A <= 0.0 ) return( exp( A*S*S ) );
    if ( A > 0.0 )  return( exp( A*(S*S-1.0) ) );
}



double  TKSI_T96( double XKSI_T96, double XKS0, double DXKSI_T96 ) {

    double  R, TDZ3, g, BR3;

    TDZ3 = 2.0*DXKSI_T96*DXKSI_T96*DXKSI_T96;

    if ( (XKSI_T96-XKS0) < -DXKSI_T96 ) R = 0.0;
    if ( (XKSI_T96-XKS0) >= DXKSI_T96 ) R = 1.0;

    if ( (XKSI_T96 >= (XKS0-DXKSI_T96)) && (XKSI_T96 < XKS0) ) {
        g = XKSI_T96-XKS0+DXKSI_T96;
        BR3   = g*g*g;
        R = 1.5*BR3/(TDZ3+BR3);
    }

    if ( (XKSI_T96 >= XKS0) && (XKSI_T96 < (XKS0+DXKSI_T96)) ) {
        g = XKSI_T96-XKS0-DXKSI_T96;
        BR3 = g*g*g;
        R   = 1.0 + 1.5*BR3/(TDZ3-BR3);
    }

    return( R );

}


void    DIPOLE_T96( double PS, double X, double Y, double Z, double *BX, double *BY, double *BZ, LgmTsyg1996_Info *tInfo ) {

    /*
     *        A DOUBLE PRECISION ROUTINE
     *
     *    CALCULATES GSM COMPONENTS OF A GEODIPOLE_T96 FIELD WITH THE DIPOLE_T96 MOMENT
     *    CORRESPONDING TO THE EPOCH OF 1980.
     *
     *  ----INPUT PARAMETERS:
     *       PS - GEODIPOLE_T96 TILT ANGLE IN RADIANS,
     *       X,Y,Z - GSM COORDINATES IN RE (1 RE = 6371.2 km)
     *
     *  ----OUTPUT PARAMETERS:
     *       BX,BY,BZ - FIELD COMPONENTS IN GSM SYSTEM, IN NANOTESLA.
     */
    double    SPS, CPS, P, U, V, T, s, s2, s3, s5, Q;

    //SPS = sin(PS); CPS = cos(PS);
    SPS = tInfo->sin_psi; CPS = tInfo->cos_psi;
    P = X*X;
    U = Z*Z;
    V = 3.0*Z*X;
    T = Y*Y;
    s = sqrt(P+T+U); s2 = s*s; s3 = s2*s; s5 = s2*s3;
    Q = 30574.0/s5;

    *BX = Q*((T+U-2.0*P)*SPS-V*CPS);
    *BY = -3.0*Y*Q*(X*SPS+Z*CPS);
    *BZ = Q*((P+T-2.0*U)*CPS-V*SPS);

    return;

}


