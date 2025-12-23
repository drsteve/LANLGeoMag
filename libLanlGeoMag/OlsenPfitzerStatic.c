#include <stdio.h>
#include <math.h>
#include "Lgm/Lgm_MagModelInfo.h"
/**
 *  Ported to C by Michael G. Henderson Spetember 17, 2010.
 *
 *  VERSION 11/01/76
 *
 *     PURPOSE
 *       TO CALCULATE THE CONTRIBUTION TO THE EARTHS MAGNETIC FIELD BY
 *       SOURCES EXTERNAL TO THE EARTH.  NO INTERNAL FIELD IS INCLUDED
 *       IN THIS ROUTINE.
 *
 *     METHOD
 *        THE ROUTINE INCLUDES THE FIELD CONTRIBUTIONS FROM THE
 *        MAGNETOPAUSE CURRENTS,  AND CURRENTS DISTRIBUTED THROUGHOUT
 *        THE MAGNETOSPHERE (THE TAIL AND RING CURRENTS).  IT IS VALID
 *        FOR ALL TILTS OF THE EARTHS DIPOLE AXIS AND IS VALID DURING
 *        QUIET MAGNETIC CONDITIONS.
 *        A GENERALIZED ORTHONORMAL LEAST SQUARES PROGRAM WAS USED
 *        TO FIT THE COEFFICIENTS OF A POWER SERIES (INCLUDING
 *        EXPONENTIAL TERMS) THROUGH FOURTH ORDER IN SPACE AND
 *        THIRD ORDER IN TILT.  THIS EXPANSION HAS BEEN OPTIMIZED
 *        FOR THE NEAR EARTH REGION AND IS VALID TO 15 EARTH RADII.
 *        OUTSIDE OF THIS REGION THE FIELD DIVERGES RAPIDLY AND A
 *        TEMPLATE SETS THE FIELD TO ZERO.  IN ORDER TO IMPROVE
 *        COMPUTATIONAL SPEED THE FIELD IS SET TO ZERO BELOW 2 EARTH
 *        RADII.  (IN THIS REGION THE EARTHS INTERNAL FIELD DOMINATES
 *        AND THE VARIATIONS EXCPRESSED BY THIS EXPANSION IS NOT
 *        SUFFICIENTLY ACCURATE THE PREDICT VARIATIONS ON THE EARTHS
 *        SURFACE)
 *
 *        THE POWER SERIES REPRESENTING THE MAGNETIC FIELD IS
 *        BX = SUM OVER I, J, K OF ( A(I, J, K)*X**(I-1)*Y**(2*J-2)*Z**(K-1)
 *               + B(I, J, K)*X**(I-1)*Y**(2*J-2)*Z**(K-1)*EXP(-.06*R**2))
 *            I GOES FROM 1 TO 5, J FROM 1 TO 3,  K FROM 1 TO 5
 *            THE SUM OF I + 2*J + K IS LESS THAN OR EQUAL TO 9
 *        BY = SUM OVER I, J, K OF ( C(I, J, K)*X**(I-1)*Y**(2*J-1)*Z**(K-1)
 *               + D(I, J, K)*X**(I-1)*Y**(2*J-1)*Z**(K-1)*EXP(-.06*R**2))
 *            I GOES FROM 1 TO 5,  J FROM 1 TO 3,  K FROM 1 TO 5
 *            THE SUM OF I + 2*J+1 + K IS LESS THAN OR EQUAL TO 9
 *        BZ = SUM OVER I, J, K OF ( E(I, J, K)*X**(I-1)*Y**(2*J-2)*Z**(K-1)
 *               + F(I, J, K)*X**(I-1)*Y**(2*J-2)*Z**(K-1)*EXP(-.06*R**2))
 *            I GOES FROM 1 TO 5,  J FROM 1 TO 3,  K FROM 1 TO 5
 *            THE SUM OF I + 2*J + K IS LESS THAN OR EQUAL TO 9
 *        THE COEFFICIENTS A-F ARE DEPENDENT ONLY ON POSITION AND
 *        ARE RECALCULATED EACH TIME THE TILT OF THE DIPOLE IS CHANGED.
 *        THE COEEFICIENTS A-F ARE DETERMINED FROM THE TILT DEPENDENT
 *        CONSTANTS AA-FF BY THE FOLLOWING EXPRESSIONS
 *        A(I, J, K) = AA(I, J, K, 1)*TILT**(K-1-(K-1)/2*2)
 *                 +AA(I, J, K, 2)*TILT**(K+1-(K-1)/2*2)
 *        B(I, J, K) = BB......
 *        C(I, J, K) = CC......
 *        D(I, J, K) = DD......
 *        E(I, J, K) = EE(I, J, K, 1)*TILT**(K-(K)/2*2)
 *                 +EE(I, J, K, 2)*TILT**(K+2-(K)/2*2)
 *        F(I, J, K) = FF......
 *
 *     INPUT -- CALLING SEQUENCE
 *        XX     A REAL ARRAY GIVING THE POSITION WHERE THE MAGNETIC
 *               FIELD IS TO BE EVALUATED.  XX(1),  XX(2),  XX(3) ARE
 *               RESPECTIVELY THE X,  Y,  AND Z SOLAR MAGNETIC
 *               COORDINATES IN EARTH RADII. Z IS ALONG THE EARTHS
 *               NORTH DIPOLE AXIS,  X IS PERPENDICULAR TO Z AND IN THE
 *               PLANE CONTAINING THE Z AXIS AND THE SUN-EARTH LINE, 
 *               Y IS PERPENDICULAR TO X AND Z FORMING A RIGHT HANDED
 *               COORDINATE SYSTEM. X IS POSITIVE IN THE SOLAR DIRECTION.
 *        TILT   IS THE TILT OF THE DIPOLE AXIS IN DEGREES.  IT IS
 *               THE COMPLEMENT OF THE ANGLE BETWEEEN THE NORTH DIPOLE
 *               AXIS AND THE SOLAR DIRECTION (PSI). TILT = 90-PSI.
 *
 *     OUTPUT -- CALLING SEQUENCE
 *        BF     A REAL ARRAY CONTAING THE X,  Y,  AND Z COMPONENTS OF
 *               THE MAGNETOSPHERIC MAGNETIC FIELD IN GAMMA. BF(1), 
 *               BF(2) AND BF(3) ARE THE BX,  BY,  BZ COMPONENTS.
 *
 *     CONSTANTS
 *        AA, BB, CC, DD, EE, FF ARE REAL ARRAYS CONTAING THE TILT DEPENDED
 *               COEFFICIENTS.
 *               AA(I, J, K, L) ARE STORED SUCH THAT L VARIES MOST RAPIDLY
 *               FOLLOWED IN ORDER BY K,  J AND I.  I VARIES THE SLOWEST.
 *               THE ARRAY IS CLOSE PACKED AND ALL COEFFICIENTS THAT
 *               ARE ZERO BECAUSE OF SYMMETRY OR BECAUSE THE CROSS TERM
 *               POWER IS TOO LARGE ARE DELETED.
 *
 *     VARIABLES
 *        A, B, C, D, E, F  THE TILT INDEPENDENT COEFFICIENTS.  THEIR USE
 *               IS DESCRIBED UNDER METHOD.
 *        ITA    A REAL ARRAY WHICH CONTAINS THE SYMMETRY OF THE TILT
 *               DEPENDENCE FOR EACH OF THE A AND B COEFFICIENTS
 *               ITA(1) HAS THE SYMMETRY INFORMATION FOR A(1, 1, 1, 1)
 *                      AND A(1, 1, 1, 2)
 *               ITA(2) HAS THE SYMMETRY INFORMATION FOR A(1, 1, 2, 1)
 *                      AND A(1, 1, 2, 2)  ETC.
 *               IF ITA  =  1 TILT SYMMETRY IS EVEN WITH RESPECT TO Z SYM.
 *               IF ITA  =  2 TILT SYMMETRY IS ODD WITH RESPECT TO Z SYM.
 *        ITB    SYMMETRY POINTER FOR C AND D ARRAYS
 *        ITC    SYMMETRY POINTER FOR E AND F ARRAYS
 *        X      X COMPONENT OF POSITION
 *        Y      Y COMPONENT OF POSITION
 *        Z      Z COMPONENT OF POSITION
 *        Y2     Y**2
 *        Z2     Z**2
 *        R2     X**2 + Y**2 + Z**2
 *        R      SQRT(R2)
 *        I      DO LOOP VARIABLE.  IN THE FIELD EXPANSION LOOP IT
 *               REPRESENTS THE POWER TO WHICH X IS CURRENTLY RAISED
 *               I.E. X**(I-1)
 *        J      DO LOOP VARIABLE.  ALSO Y**(2*J-2)
 *        K      DO LOOP VARIABLE.  ALSO Z**(K-1)
 *        XB     X**(I-1)
 *        YEXB   X**(I-1)*Y**(2*J-2)
 *        ZEYEXB X**(I-1)*Y**(2*J-2)*Z**(K-1)
 *        IJK    I + 2*J + K
 *        II     POINTS TO THE ARRAY LOCATION WHERE THE CURRENT POWER
 *               SERIES COEFFICIENT FOR BX IS LOCATED
 *        JJ     BY COEFFICIENT LOCATION POINTER
 *        KK     BZ COEFFICIENT LOCATION POINTER
 *        BX, BY, BZ ARE USED TO CONSTRUCT THE MAGNETIC FIELD WITHIN THE
 *               POWER SERIES LOOP.
 *        EXPR   EXP(-.06*R2)
 *        TILTL  HOLDS THE LAST VALUE OF THE TILT FOR WHICH THE TILT
 *               INDEPENDENT COEFFICIENTS A-F WERE CALCULATED
 *        TT     A REAL ARRAY HOLDING THE POWERS OF THE TILT.
 *               TT(1) = TILT**0,  TT(2) = TILT**1,  ETC.
 *        CON     = 0 FOR R LESS THAN 2
 *                = 1 FOR R GREATER THAN 2.5
 *               GOES FROM 0 TO 1 IN THE REGION 2 TO 2.5
 *
 *     FOR MORE INFORMATION CALL OR WRITE K. A. PFITZER OR W. P. OLSON
 *     AT MCDONNEL DOUGLAS ASTRONAUTICS CO. 5301 BOLSA AVE,  HUNTINGTON
 *     CALIF.,  PHONE (714) 896-3231.
 *
 *
 *  \warning    The original code artificially bails out with zero field when
 *              r<2 and when r>15 Re. This may cause discontinuities that can cause
 *              trouble. This version tries to smoothly transition to dipole at large
 *              distances.
 *
 */

void OlsenPfitzerStatic( double XX[], double BF[], double TILT, Lgm_MagModelInfo *m ) {

    int     IDX, J, K, Jp1, II, JJ, KK, IJK;
    double  X, Y, Z, X2, Y2, Z2, R2, MYFAC;
    double  BX, BY, BZ, CON, EXPR, XB, YEXB, ZEYEXB;

    static int ITA[]  =  { -99,  2, 1, 2, 1, 2, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 2, 1, 2, 2, 2, 1, 2, 1, 2, 1, 2, 1, 2, 2, 2, 1};
    static int ITB[]  =  { -99,  2, 1, 2, 1, 2, 2, 1, 2, 2, 2, 1, 2, 1, 2, 1, 2, 1, 2, 2, 2, 1, 2};
    static int ITC[]  =  { -99,  1, 2, 1, 2, 1, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 1, 2, 1, 1, 1, 2, 1, 2, 1, 2, 1, 2, 1, 1, 1, 2};

    static double AA[]  =  { -9e99,  -2.26836E-02, -1.01863E-04,  3.42986E+00, 
      -3.12195E-04,  9.50629E-03, -2.91512E-06, -1.57317E-03,  8.62856E-08, 
      -4.26478E-05,  1.62924E-08, -1.27549E-04,  1.90732E-06, -1.65983E-02, 
       8.46680E-09, -5.55850E-05,  1.37404E-08,  9.91815E-05,  1.59296E-08, 
       4.52864E-07, -7.17669E-09,  4.98627E-05,  3.33662E-10, -5.97620E-02, 
       1.60669E-05, -2.29457E-01, -1.43777E-04,  1.09403E-03, -9.15606E-07, 
       1.60658E-03, -4.01198E-07, -3.15064E-06,  2.03125E-09,  4.92887E-04, 
      -1.80676E-07, -1.12022E-03,  5.98568E-07, -5.90009E-06,  5.16504E-09, 
      -1.48737E-06,  4.83477E-10, -7.44379E-04,  3.82472E-06,  7.41737E-04, 
      -1.31468E-05, -1.24729E-04,  1.92930E-08, -1.91764E-04, -5.30371E-08, 
       1.38186E-05, -2.81594E-08,  7.46386E-06,  2.64404E-08,  2.45049E-04, 
      -1.81802E-07, -1.00278E-03,  1.98742E-06, -1.16425E-05,  1.17556E-08, 
      -2.46079E-06, -3.45831E-10,  1.02440E-05, -1.90716E-08, -4.00855E-05, 
       1.25818E-07};

    static double BB[]  =  { -9e99,  9.47753E-02,  1.45981E-04, -1.82933E+00, 
       5.54882E-04,  5.03665E-03, -2.07698E-06,  1.10959E-01, -3.45837E-05, 
      -4.40075E-05,  5.06464E-07, -1.20112E-03,  3.64911E-06,  1.49849E-01, 
      -7.44929E-05,  2.46382E-04,  9.65870E-07, -9.54881E-04,  2.43647E-07, 
       3.06520E-04,  3.07836E-07,  6.48301E-03,  1.26251E-06, -7.09548E-03, 
      -1.55596E-05,  3.06465E+00, -7.84893E-05,  4.95145E-03,  3.71921E-06, 
      -1.52002E-01,  6.81988E-06, -8.55686E-05, -9.01230E-08, -3.71458E-04, 
       1.30476E-07, -1.82971E-01,  1.51390E-05, -1.45912E-04, -2.22778E-07, 
       6.49278E-05, -3.72758E-08, -1.59932E-03,  8.04921E-06,  5.38012E-01, 
      -1.43182E-04,  1.50000E-04,  5.88020E-07, -1.59000E-02,  1.60744E-06, 
       3.17837E-04,  1.78959E-07, -8.93794E-03,  6.37549E-06,  1.27887E-03, 
      -2.45878E-07, -1.93210E-01,  6.91233E-06, -2.80637E-04, -2.57073E-07, 
       5.78343E-05,  4.52128E-10,  1.89621E-04, -4.84911E-08, -1.50058E-02, 
       6.21772E-06};

    static double CC[]  =  { -9e99,  -1.88177e-02, -1.92493e-06, -2.89064e-01, 
      -8.49439E-05, -4.76380E-04, -4.52998E-08,  1.61086E-03,  3.18728E-07, 
       1.29159E-06,  5.52259E-10,  3.95543E-05,  5.61209E-08,  1.38287E-03, 
       5.74237E-07,  1.86489E-06,  7.10175E-10,  1.45243E-07, -2.97591E-10, 
      -2.43029E-03, -6.70000E-07, -2.30624E-02, -6.22193E-06, -2.40815E-05, 
       2.01689E-08,  1.76721E-04,  3.78689E-08,  9.88496E-06,  7.33820E-09, 
       7.32126E-05,  8.43986E-08,  8.82449E-06, -6.11708E-08,  1.78881E-04, 
       8.62589E-07,  3.43724E-06,  2.53783E-09, -2.04239E-07,  8.16641E-10, 
       1.68075E-05,  7.62815E-09,  2.26026E-04,  3.66341E-08,  3.44637E-07, 
       2.25531E-10};

    static double DD[]  =  { -9e99,  2.50143E-03,  1.01200E-06,  3.23821E+00, 
       1.08589E-05, -3.39199E-05, -5.27052E-07, -9.46161E-02, -1.95413E-09, 
      -4.23614E-06,  1.43153E-08, -2.62948E-04,  1.05138E-07, -2.15784E-01, 
      -2.20717E-07, -2.65687E-05,  1.26370E-08,  5.88917E-07, -1.13658E-08, 
       1.64385E-03,  1.44263E-06, -1.66045E-01, -1.46096E-05,  1.22811E-04, 
       3.43922E-08,  9.66760E-05, -6.32150E-07, -4.97400E-05, -2.78578E-08, 
       1.77366E-02,  2.05401E-07, -1.91756E-03, -9.49392E-07, -1.99488E-01, 
      -2.07170E-06, -5.40443E-05,  1.59289E-08,  7.30914E-05,  3.38786E-08, 
      -1.59537E-04, -1.65504E-07,  1.90940E-02,  2.03238E-06,  1.01148E-04, 
       5.20815E-08};

    static double EE[]  =  { -9e99,  -2.77924E+01, -1.01457E-03,  9.21436E-02, 
      -8.52177E-06,  5.19106E-01,  8.28881E-05, -5.59651E-04,  1.16736E-07, 
      -2.11206E-03, -5.35469E-07,  4.41990E-01, -1.33679E-05, -7.18642E-04, 
       6.17358E-08, -3.51990E-03, -5.29070E-07,  1.88443E-06, -6.60696E-10, 
      -1.34708E-03,  1.02160E-07,  1.58219E-06,  2.05040E-10,  1.18039E+00, 
       1.58903E-04,  1.86944E-02, -4.46477E-06,  5.49869E-02,  4.94690E-06, 
      -1.18335E-04,  6.95684E-09, -2.73839E-04, -9.17883E-08,  2.79126E-02, 
      -1.02567E-05, -1.25427E-04,  3.07143E-08, -5.31826E-04, -2.98476E-08, 
      -4.89899E-05,  4.91480E-08,  3.85563E-01,  4.16966E-05,  6.74744E-04, 
      -2.08736E-07, -3.42654E-03, -3.13957E-06, -6.31361E-06, -2.92981E-09, 
      -2.63883E-03, -1.32235E-07, -6.19406E-06,  3.54334E-09,  6.65986E-03, 
      -5.81949E-06, -1.88809E-04,  3.62055E-08, -4.64380E-04, -2.21159E-07, 
      -1.77496E-04,  4.95560E-08, -3.18867E-04, -3.17697E-07, -1.05815E-05, 
       2.22220E-09};

    static double FF[]  =  { -9e99,  -5.07092E+00,  4.71960E-03, -3.79851E-03, 
      -3.67309E-06, -6.02439E-01,  1.08490E-04,  5.09287E-04,  5.62210E-07, 
       7.05718E-02,  5.13160E-06, -2.85571E+00, -4.31728E-05,  1.03185E-03, 
       1.05332E-07,  1.04106E-02,  1.60749E-05,  4.18031E-05,  3.32759E-08, 
       1.20113E-01,  1.40486E-05, -3.37993E-05,  5.48340E-09,  9.10815E-02, 
      -4.00608E-04,  3.75393E-03, -4.69939E-07, -2.48561E-02,  1.31836E-04, 
      -2.67755E-04, -7.60285E-08,  3.04443E-03, -3.28956E-06,  5.82367E-01, 
       5.39496E-06, -6.15261E-04,  4.05316E-08,  1.13546E-02, -4.26493E-06, 
      -2.72007E-02,  5.72523E-08, -2.98576E+00,  3.07325E-05,  1.51645E-03, 
       1.25098E-06,  4.07213E-02,  1.05964E-05,  1.04232E-04,  1.77381E-08, 
       1.92781E-01,  2.15734E-05, -1.65741E-05, -1.88683E-09,  2.44803E-01, 
       1.51316E-05, -3.01157E-04,  8.47006E-08,  1.86971E-02, -6.94074E-06, 
       9.13198E-03, -2.38052E-07,  1.28552E-01,  6.92595E-06, -8.36792E-05, 
      -6.10021E-08};


    /*
     * SET UP SOME OF THE INITIAL POSITION VARIABLES
     */
    X  = XX[1]; X2 = X*X;
    Y  = XX[2]; Y2 = Y*Y;
    Z  = XX[3]; Z2 = Z*Z;
    R2 = X2 + Y2 + Z2;


    /*
     * SET MAGNETIC FIELD VARIABLES TO ZERO
     */
    BX = BY = BZ = 0.0;


    /*
     * CHECK TO SEE IF POSITION IS WITHIN REGION OF VALIDITY
     */
    CON = 1.0;


    /*
     * IF DISTANCE TOO LARGE TAKE ERROR EXIT
     */
    if ( R2 > 225.0 ) {
        MYFAC = R2*R2*R2/11390625.0; // Kludge to avoid power series from blowing up at large R2
//        printf("OlsenPfitzerStatic: ( X, Y, Z ) = ( %10.3g, %10.3g, %10.3g ) BX, BY, BZ = %g %g %g IS OUTSIDE THE VALID REGION--POWER SERIES DIVERGES BFIELD IS SET TO ZERO\n", X, Y, Z, BX, BY, BZ );
//        BF[1] = 0.0;
//        BF[2] = 0.0;
//        BF[3] = 0.0;
//        return;
    }
/*
*/


    /*
     * IF DISTANCE TOO SMALL SET FIELD TO ZERO AND EXIT
     */
    if ( R2 < 4.0 ){
        BF[1] = 0.0;
        BF[2] = 0.0;
        BF[3] = 0.0;
        return;
    }

    if ( R2 < 6.25 ) CON *= (R2-4.0)/2.25;





    /*
     * IF TILT HAS NOT CHANGED,  GO DRECTLY TO FIELD CALCULATION
     */
    if ( m->OP77_TILTL != TILT ) {

        /*
         * SET UP POWERS OF TILT
         */
        m->OP77_TILTL = TILT;
        m->OP77_TT[1] = 1.0;
        m->OP77_TT[2] = m->OP77_TILTL;
        m->OP77_TT[3] = m->OP77_TILTL*m->OP77_TILTL;
        m->OP77_TT[4] = m->OP77_TILTL*m->OP77_TT[3];

        /*
         * SET UP THE X AND Z COMPONENT TILT INDEPENDENT COEFFICIENTS
         */
        for ( IDX=1; IDX<=32; IDX++ ) {
            J = 2*IDX-1;
            Jp1 = J+1;
            K = ITA[IDX];
            m->OP77_A[IDX] = AA[J]*m->OP77_TT[K] + AA[Jp1]*m->OP77_TT[K+2];
            m->OP77_B[IDX] = BB[J]*m->OP77_TT[K] + BB[Jp1]*m->OP77_TT[K+2];
            K = ITC[IDX];
            m->OP77_E[IDX] = EE[J]*m->OP77_TT[K] + EE[Jp1]*m->OP77_TT[K+2];
            m->OP77_F[IDX] = FF[J]*m->OP77_TT[K] + FF[Jp1]*m->OP77_TT[K+2];
        }

        /*
         * SET UP THE Y COMPONENT TILT INDEPENDENT COEFFICIENTS
         */
        for ( IDX=1; IDX<=22; IDX++ ) {
            J = 2*IDX-1;
            K = ITB[IDX];
            m->OP77_C[IDX] = CC[J]*m->OP77_TT[K] + CC[J+1]*m->OP77_TT[K+2];
            m->OP77_D[IDX] = DD[J]*m->OP77_TT[K] + DD[J+1]*m->OP77_TT[K+2];
        }

    }

    EXPR = exp( -0.06*R2 );






    /*
     *  INITIALIZE THE POINTERS
     */
    II = 1; JJ = 1; KK = 1; XB = 1.0;


    /* 
     * BEGIN SUM OVER X
     */
    for ( IDX=1; IDX<=5; IDX++ ) {

        YEXB = XB;

        /* 
         * BEGIN SUM OVER X
         */
        for ( J=1; J<=3; J++ ) {
    
            if (IDX+2*J >  8) break;
            ZEYEXB = YEXB;
            IJK = IDX + 2*J + 1;
            K = 1;

            /* 
             * Z LOOP STARTS HERE
             */
            do {
                BZ += (m->OP77_E[KK]+m->OP77_F[KK]*EXPR)*ZEYEXB;
                ++KK;
                BX += (m->OP77_A[II]+m->OP77_B[II]*EXPR )*ZEYEXB;
                ++II;
                if ( IJK > 8 ) break;
                BY += (m->OP77_C[JJ]+m->OP77_D[JJ]*EXPR )*ZEYEXB*Y;
                ++JJ;
                ZEYEXB *= Z;
                ++IJK;
                ++K;
            } while ( (IJK <= 9) && (K <= 5) ) ;


            YEXB *= Y2;

        }

        XB *= X;

    }


    /*
     * SET UP THE OUTPUT ARRAY,  MULTIPLY BY CON. CON IS NORMALY ONE
     * BUT INSIDE OF R = 2.5 IT GOES TO ZERO. INSIDE R = 2 IT IS ZERO.
     */
    BF[1] = BX*CON;
    BF[2] = BY*CON;
    BF[3] = BZ*CON;

    if ( R2 > 225.0 ) {
        // Kludge to avoid diverging power series. MYFAC has an R^6 in the demon.
        BF[1] = BX/MYFAC;
        BF[2] = BY/MYFAC;
        BF[3] = BZ/MYFAC;
        //printf("OlsenPfitzerStatic: ( X, Y, Z ) = ( %10.3g, %10.3g, %10.3g ) BX, BY, BZ = %g %g %g IS OUTSIDE THE VALID REGION--POWER SERIES DIVERGES BFIELD IS SET TO ZERO\n", X, Y, Z, BX, BY, BZ );
    }


    return;



}



/*
 *   $Id$
 */

