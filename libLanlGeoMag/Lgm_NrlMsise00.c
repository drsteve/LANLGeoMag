#include "Lgm/Lgm_NrlMsise00.h"

static const double RGAS = 831.4;

#define ZETA(ZZ,ZL,p) ( (ZZ-ZL)*(p->RE+ZL)/(p->RE+ZZ) )

/*----------------------------------------------------------------------
 *    SUBROUTINE GTD7(IYD,SEC,ALT,GLAT,GLONG,STL,F107A,F107,AP,MASS,D,T)
 *
 *    NRLMSISE-00
 *    -----------
 *       Neutral Atmosphere Empirical Model from the surface to lower
 *       exosphere
 *
 *       NEW FEATURES:
 *         *Extensive satellite drag database used in model generation
 *         *Revised O2 (and O) in lower thermosphere
 *         *Additional nonlinear solar activity term
 *         *"ANOMALOUS OXYGEN" NUMBER DENSITY, OUTPUT D(9)
 *          At high altitudes (> 500 km), hot atomic oxygen or ionized
 *          oxygen can become appreciable for some ranges of subroutine
 *          inputs, thereby affecting drag on satellites and debris. We
 *          group these species under the term "anomalous oxygen," since
 *          their individual variations are not presently separable with
 *          the drag data used to define this model component.
 *
 *       SUBROUTINES FOR SPECIAL OUTPUTS:
 *       
 *       HIGH ALTITUDE DRAG: EFFECTIVE TOTAL MASS DENSITY 
 *       (SUBROUTINE GTD7D, OUTPUT D(6))
 *          For atmospheric drag calculations at altitudes above 500 km,
 *          call SUBROUTINE GTD7D to compute the "effective total mass
 *          density" by including contributions from "anomalous oxygen."
 *          See "NOTES ON OUTPUT VARIABLES" below on D(6).
 *
 *       PRESSURE GRID (SUBROUTINE GHP7)
 *         See subroutine GHP7 to specify outputs at a pressure level
 *         rather than at an altitude.
 *
 *       OUTPUT IN M-3 and KG/M3:   CALL METERS(.TRUE.)
 *
 *    INPUT VARIABLES:
 *       IYD - YEAR AND DAY AS YYDDD (day of year from 1 to 365 (or 366))
 *             (Year ignored in current model)
 *       SEC - UT(SEC)
 *       ALT - ALTITUDE(KM)
 *       GLAT - GEODETIC LATITUDE(DEG)
 *       GLONG - GEODETIC LONGITUDE(DEG)
 *       STL - LOCAL APPARENT SOLAR TIME(HRS; see Note below)
 *       F107A - 81 day AVERAGE OF F10.7 FLUX (centered on day DDD)
 *       F107 - DAILY F10.7 FLUX FOR PREVIOUS DAY
 *       AP - MAGNETIC INDEX(DAILY) OR WHEN SW(9)=-1. :
 *          - ARRAY CONTAINING:
 *            (1) DAILY AP
 *            (2) 3 HR AP INDEX FOR CURRENT TIME
 *            (3) 3 HR AP INDEX FOR 3 HRS BEFORE CURRENT TIME
 *            (4) 3 HR AP INDEX FOR 6 HRS BEFORE CURRENT TIME
 *            (5) 3 HR AP INDEX FOR 9 HRS BEFORE CURRENT TIME
 *            (6) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 12 TO 33 HRS PRIOR
 *                   TO CURRENT TIME
 *            (7) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 36 TO 57 HRS PRIOR
 *                   TO CURRENT TIME
 *       MASS - MASS NUMBER (ONLY DENSITY FOR SELECTED GAS IS
 *                CALCULATED.  MASS 0 IS TEMPERATURE.  MASS 48 FOR ALL.
 *                MASS 17 IS Anomalous O ONLY.)
 *
 *    NOTES ON INPUT VARIABLES: 
 *       UT, Local Time, and Longitude are used independently in the
 *       model and are not of equal importance for every situation.  
 *       For the most physically realistic calculation these three
 *       variables should be consistent (STL=SEC/3600+GLONG/15).
 *       The Equation of Time departures from the above formula
 *       for apparent local time can be included if available but
 *       are of minor importance.
 *
 *       F107 and F107A values used to generate the model correspond
 *       to the 10.7 cm radio flux at the actual distance of the Earth
 *       from the Sun rather than the radio flux at 1 AU. The following
 *       site provides both classes of values:
 *       ftp://ftp.ngdc.noaa.gov/STP/SOLAR_DATA/SOLAR_RADIO/FLUX/
 *
 *       F107, F107A, and AP effects are neither large nor well
 *       established below 80 km and these parameters should be set to
 *       150., 150., and 4. respectively.
 *
 *    OUTPUT VARIABLES:
 *       D(1) - HE NUMBER DENSITY(CM-3)
 *       D(2) - O NUMBER DENSITY(CM-3)
 *       D(3) - N2 NUMBER DENSITY(CM-3)
 *       D(4) - O2 NUMBER DENSITY(CM-3)
 *       D(5) - AR NUMBER DENSITY(CM-3)                       
 *       D(6) - TOTAL MASS DENSITY(GM/CM3)
 *       D(7) - H NUMBER DENSITY(CM-3)
 *       D(8) - N NUMBER DENSITY(CM-3)
 *       D(9) - Anomalous oxygen NUMBER DENSITY(CM-3)
 *       T(1) - EXOSPHERIC TEMPERATURE
 *       T(2) - TEMPERATURE AT ALT
 *
 *    NOTES ON OUTPUT VARIABLES:
 *       TO GET OUTPUT IN M-3 and KG/M3:   CALL METERS(.TRUE.) 
 *
 *       O, H, and N are set to zero below 72.5 km
 *
 *       T(1), Exospheric temperature, is set to global average for
 *       altitudes below 120 km. The 120 km gradient is left at global
 *       average value for altitudes below 72 km.
 *
 *       D(6), TOTAL MASS DENSITY, is NOT the same for subroutines GTD7 
 *       and GTD7D
 *
 *         SUBROUTINE GTD7 -- D(6) is the sum of the mass densities of the
 *         species labeled by indices 1-5 and 7-8 in output variable D.
 *         This includes He, O, N2, O2, Ar, H, and N but does NOT include
 *         anomalous oxygen (species index 9).
 *
 *         SUBROUTINE GTD7D -- D(6) is the "effective total mass density
 *         for drag" and is the sum of the mass densities of all species
 *         in this model, INCLUDING anomalous oxygen.
 *       
 *    SWITCHES: The following is for test and special purposes:
 *         
 *       TO TURN ON AND OFF PARTICULAR VARIATIONS CALL TSELEC(SW),
 *       WHERE SW IS A 25 ELEMENT ARRAY CONTAINING 0. FOR OFF, 1. 
 *       FOR ON, OR 2. FOR MAIN EFFECTS OFF BUT CROSS TERMS ON
 *       FOR THE FOLLOWING VARIATIONS
 *              1 - F10.7 EFFECT ON MEAN  2 - TIME INDEPENDENT
 *              3 - SYMMETRICAL ANNUAL    4 - SYMMETRICAL SEMIANNUAL
 *              5 - ASYMMETRICAL ANNUAL   6 - ASYMMETRICAL SEMIANNUAL
 *              7 - DIURNAL               8 - SEMIDIURNAL
 *              9 - DAILY AP             10 - ALL UT/LONG EFFECTS
 *             11 - LONGITUDINAL         12 - UT AND MIXED UT/LONG
 *             13 - MIXED AP/UT/LONG     14 - TERDIURNAL
 *             15 - DEPARTURES FROM DIFFUSIVE EQUILIBRIUM
 *             16 - ALL TINF VAR         17 - ALL TLB VAR
 *             18 - ALL TN1 VAR           19 - ALL S VAR
 *             20 - ALL TN2 VAR           21 - ALL NLB VAR
 *             22 - ALL TN3 VAR           23 - TURBO SCALE HEIGHT VAR
 *
 *       To get current values of SW: CALL TRETRV(SW)
 */
void GTD7( int IYD, double SEC, double ALT, double GLAT, double GLONG, double STL, double F107A, double F107, double *AP, double MASS, double *D, double *T, Lgm_Msis00Info *p ) {

    int     J, V1;
//    double  DS[10], TS[3];
    double  XLAT, ALTT, MSS, DM28M, g, g2, DMC, DZ28, DMR, TZ, XMM;


/*
    static double MN3   = 5;
    static double ZN3[] = { -9e99, 32.5, 20.0, 15.0, 10.0, 0.0 };
    static double MN2   = 4;
    static double ZN2[] = { -9e99, 72.5, 55.0, 45.0, 32.5 };
    static const double ZMIX  = 62.5;
*/
    double MN3   = 5;
    double ZN3[] = { -9e99, 32.5, 20.0, 15.0, 10.0, 0.0 };
    double MN2   = 4;
    double ZN2[] = { -9e99, 72.5, 55.0, 45.0, 32.5 };
    const double ZMIX  = 62.5;
//    static double SV[]  = {  /25*1./ };



    if ( p->ISW != 64999 ) TSELEC( p->SV, p );

    /*
     * Test for changed input
     */
    V1 = VTST7( IYD, SEC, GLAT, GLONG, STL, F107A, F107, AP, 1, p );
V1 =1;

    /*
     * Latitude variation of gravity (none for SW(2)=0)
     */
    XLAT = GLAT;
    if ( p->SW[2] == 0 ) XLAT = 45.0;
    GLATF( XLAT, &p->GSURF, &p->RE );

    XMM = p->PDM[5][3];


    /*
     * THERMOSPHERE/MESOSPHERE (above ZN2(1))
     */
    ALTT = AMAX1( ALT, ZN2[1] );
    MSS  = MASS;


    /*
     * Only calculate N2 in thermosphere if alt in mixed region
     */
    if( (ALT < ZMIX) && (MASS > 0) ) MSS = 28;


    /*
     * Only calculate thermosphere if input parameters changed
     * or altitude above ZN2(1) in mesosphere
     */
    if ( (V1 == 1) || (ALT > ZN2[1]) || (p->GTD7_ALAST > ZN2[1]) || (MSS != p->GTD7_MSSL) ) { 
        GTS7( IYD, SEC, ALTT, GLAT, GLONG, STL, F107A, F107, AP, MSS, p->DS, p->TS, p );
        DM28M = p->DM28;
        if (p->IMR == 1 ) DM28M = p->DM28*1.0e6; // metric adjustment
        p->GTD7_MSSL = MSS;
    }

    T[1] = p->TS[1];
    T[2] = p->TS[2];

    if ( ALT >= ZN2[1] ) {
        for ( J=1; J<=9; J++ ) D[J] = p->DS[J];
        p->GTD7_ALAST = ALT;
        return;
    }

    /*
     *  LOWER MESOSPHERE/UPPER STRATOSPHERE [between ZN3[1] and ZN2[1]]
     *     Temperature at nodes and gradients at end nodes
     *     Inverse temperature a linear function of spherical harmonics
     *     Only calculate nodes if input changed
     */
    if ( (V1 == 1) || (p->GTD7_ALAST >= ZN2[1]) ) {
        p->TGN2[1] = p->TGN1[2];
        p->TN2[1]  = p->TN1[5];
        p->TN2[2]  = p->PMA[1][1]*p->PAVGM[1]/(1.-p->SW[20]*GLOB7S(p->PMA_rc[1], p));
        p->TN2[3]  = p->PMA[1][2]*p->PAVGM[2]/(1.-p->SW[20]*GLOB7S(p->PMA_rc[2], p));
        p->TN2[4]  = p->PMA[1][3]*p->PAVGM[3]/(1.-p->SW[20]*p->SW[22]* GLOB7S( p->PMA_rc[3], p) );
        g = p->PMA[1][3]*p->PAVGM[3]; g2 = g*g;
        p->TGN2[2] = p->PAVGM[9]*p->PMA[1][10]*(1.0+p->SW[20]*p->SW[22]* GLOB7S( p->PMA_rc[10], p ) )*p->TN2[4]*p->TN2[4]/g2;
        p->TN3[1]  = p->TN2[4];
    }


    if ( ALT < ZN3[1] ) {

        /*
         *
         *  LOWER STRATOSPHERE AND TROPOSPHERE [below ZN3(1)]
         *     Temperature at nodes and gradients at end nodes
         *     Inverse temperature a linear function of spherical harmonics
         *     Only calculate nodes if input changed
         */
        if ( (V1 == 1) || (p->GTD7_ALAST >= ZN3[1]) ) {
            p->TGN3[1] = p->TGN2[2];
            p->TN3[2]  = p->PMA[1][4]*p->PAVGM[4]/(1.-p->SW[22]* GLOB7S( p->PMA_rc[4], p ) );
            p->TN3[3]  = p->PMA[1][5]*p->PAVGM[5]/(1.-p->SW[22]* GLOB7S( p->PMA_rc[5], p ) );
            p->TN3[4]  = p->PMA[1][6]*p->PAVGM[6]/(1.-p->SW[22]* GLOB7S( p->PMA_rc[6], p ) );
            p->TN3[5]  = p->PMA[1][7]*p->PAVGM[7]/(1.-p->SW[22]* GLOB7S( p->PMA_rc[7], p ) );
            g = p->PMA[1][7]*p->PAVGM[7]; g2 = g*g;
            p->TGN3[2] = p->PMA[1][8]*p->PAVGM[8]*(1.+p->SW[22]* GLOB7S( p->PMA_rc[1], p ) ) *p->TN3[5]*p->TN3[5]/g2;
        }

    }



    if ( MASS != 0 ) {

        // LINEAR TRANSITION TO FULL MIXING BELOW ZN2(1)
        DMC  = 0;
        if ( ALT > ZMIX ) DMC = 1.-(ZN2[1]-ALT)/(ZN2[1]-ZMIX);
        DZ28 = p->DS[3];

        // ***** N2 DENSITY ****
        DMR  = p->DS[3]/DM28M-1.0;
        D[3] = DENSM(ALT,DM28M,XMM,&TZ,MN3,ZN3,p->TN3,p->TGN3,MN2,ZN2,p->TN2,p->TGN2, p);
        D[3] = D[3]*(1.+DMR*DMC);


        // ***** HE DENSITY ****
        D[1] = 0;
        if ( (MASS == 4) || (MASS == 48) ) {
            DMR  = p->DS[1]/( DZ28*p->PDM[2][1] ) - 1.0;
            D[1] = D[3]*p->PDM[2][1]*(1.0 + DMR*DMC);
        }


        // **** O DENSITY ****
        D[2] = 0;
        D[9] = 0;



        // ***** O2 DENSITY ****
        D[4] = 0;
        if ( MASS == 32 || MASS == 48 ) {
          DMR  = p->DS[4]/(DZ28*p->PDM[2][4] ) - 1.0;
          D[4] = D[3]*p->PDM[2][4]*(1.0 + DMR*DMC);
        }


        // ***** AR DENSITY ****
        D[5] = 0;
        if ( (MASS == 40) || (MASS == 48) ) {
            DMR  = p->DS[5]/(DZ28*p->PDM[2][5]) - 1.0;
            D[5] = D[3]*p->PDM[2][5]*(1.0 + DMR*DMC);
        }


        // ***** HYDROGEN DENSITY ****
        D[7] = 0;


        // ***** ATOMIC NITROGEN DENSITY ****
        D[8] = 0;


        /*
         *
         *  TOTAL MASS DENSITY
         */
        if ( MASS == 48 ) {
            D[6] = 1.66E-24*( 4.0*D[1] + 16.0*D[2] + 28.0*D[3] + 32.0*D[4] + 40.0*D[5] + D[7] + 14.0*D[8] )  ;
            if ( p->IMR == 1 ) D[6] /= 1000.0;
        }
        T[2] = TZ;

    } else {

        p->DD   = DENSM( ALT, 1., 0, &TZ, MN3, ZN3, p->TN3, p->TGN3, MN2, ZN2, p->TN2, p->TGN2, p );
        T[2] = TZ;

    }


    p->GTD7_ALAST = ALT;
    return;

} // end void function GDT7





/*
 *    NRLMSISE-00
 *    -----------
 *       This subroutine provides Effective Total Mass Density for
 *       output D(6) which includes contributions from "anomalous
 *       oxygen" which can affect satellite drag above 500 km.  This
 *       subroutine is part of the distribution package for the 
 *       Neutral Atmosphere Empirical Model from the surface to lower
 *       exosphere.  See subroutine GTD7 for more extensive comments.
 *
 *    INPUT VARIABLES:
 *       IYD - YEAR AND DAY AS YYDDD (day of year from 1 to 365 (or 366))
 *             (Year ignored in current model)
 *       SEC - UT(SEC)
 *       ALT - ALTITUDE(KM)
 *       GLAT - GEODETIC LATITUDE(DEG)
 *       GLONG - GEODETIC LONGITUDE(DEG)
 *       STL - LOCAL APPARENT SOLAR TIME(HRS; see Note below)
 *       F107A - 81 day AVERAGE OF F10.7 FLUX (centered on day DDD)
 *       F107 - DAILY F10.7 FLUX FOR PREVIOUS DAY
 *       AP - MAGNETIC INDEX(DAILY) OR WHEN SW(9)=-1. :
 *          - ARRAY CONTAINING:
 *            (1) DAILY AP
 *            (2) 3 HR AP INDEX FOR CURRENT TIME
 *            (3) 3 HR AP INDEX FOR 3 HRS BEFORE CURRENT TIME
 *            (4) 3 HR AP INDEX FOR 6 HRS BEFORE CURRENT TIME
 *            (5) 3 HR AP INDEX FOR 9 HRS BEFORE CURRENT TIME
 *            (6) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 12 TO 33 HRS PRIOR
 *                   TO CURRENT TIME
 *            (7) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 36 TO 57 HRS PRIOR
 *                   TO CURRENT TIME
 *       MASS - MASS NUMBER (ONLY DENSITY FOR SELECTED GAS IS
 *                CALCULATED.  MASS 0 IS TEMPERATURE.  MASS 48 FOR ALL.
 *                MASS 17 IS Anomalous O ONLY.)
 *
 *    NOTES ON INPUT VARIABLES: 
 *       UT, Local Time, and Longitude are used independently in the
 *       model and are not of equal importance for every situation.  
 *       For the most physically realistic calculation these three
 *       variables should be consistent (STL=SEC/3600+GLONG/15).
 *       The Equation of Time departures from the above formula
 *       for apparent local time can be included if available but
 *       are of minor importance.
 *
 *       F107 and F107A values used to generate the model correspond
 *       to the 10.7 cm radio flux at the actual distance of the Earth
 *       from the Sun rather than the radio flux at 1 AU.
 *
 *    OUTPUT VARIABLES:
 *       D(1) - HE NUMBER DENSITY(CM-3)
 *       D(2) - O NUMBER DENSITY(CM-3)
 *       D(3) - N2 NUMBER DENSITY(CM-3)
 *       D(4) - O2 NUMBER DENSITY(CM-3)
 *       D(5) - AR NUMBER DENSITY(CM-3)                       
 *       D(6) - TOTAL MASS DENSITY(GM/CM3) [includes anomalous oxygen]
 *       D(7) - H NUMBER DENSITY(CM-3)
 *       D(8) - N NUMBER DENSITY(CM-3)
 *       D(9) - Anomalous oxygen NUMBER DENSITY(CM-3)
 *       T(1) - EXOSPHERIC TEMPERATURE
 *       T(2) - TEMPERATURE AT ALT
 */

void GTD7D( int IYD, double SEC, double ALT, double GLAT, double GLONG, double STL, double F107A, double F107, double *AP, double MASS, double *D, double *T, Lgm_Msis00Info *p ) {



    GTD7( IYD, SEC, ALT, GLAT, GLONG, STL, F107A, F107, AP, MASS, D, T, p );

    /*
     * TOTAL MASS DENSITY
     */
    if ( MASS == 48 ) {
        D[6] = 1.66E-24*( 4.0*D[1] + 16.0*D[2] + 28.0*D[3] + 32.0*D[4] + 40.0*D[5] + D[7] + 14.0*D[8] + 16.0*D[9] )  ;
        if( p->IMR == 1 ) D[6] /= 1000.0;
    }

    return;
}






/*
 *      FIND ALTITUDE OF PRESSURE SURFACE (PRESS) FROM GTD7
 *    INPUT:
 *       IYD - YEAR AND DAY AS YYDDD
 *       SEC - UT(SEC)
 *       GLAT - GEODETIC LATITUDE(DEG)
 *       GLONG - GEODETIC LONGITUDE(DEG)
 *       STL - LOCAL APPARENT SOLAR TIME(HRS)
 *       F107A - 3 MONTH AVERAGE OF F10.7 FLUX
 *       F107 - DAILY F10.7 FLUX FOR PREVIOUS DAY
 *       AP - MAGNETIC INDEX(DAILY) OR WHEN SW(9)=-1. :
 *          - ARRAY CONTAINING:
 *            (1) DAILY AP
 *            (2) 3 HR AP INDEX FOR CURRENT TIME
 *            (3) 3 HR AP INDEX FOR 3 HRS BEFORE CURRENT TIME
 *            (4) 3 HR AP INDEX FOR 6 HRS BEFORE CURRENT TIME
 *            (5) 3 HR AP INDEX FOR 9 HRS BEFORE CURRENT TIME
 *            (6) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 12 TO 33 HRS PRIOR
 *                   TO CURRENT TIME
 *            (7) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 36 TO 59 HRS PRIOR
 *                   TO CURRENT TIME
 *       PRESS - PRESSURE LEVEL(MB)
 *    OUTPUT:
 *       ALT - ALTITUDE(KM) 
 *       D(1) - HE NUMBER DENSITY(CM-3)
 *       D(2) - O NUMBER DENSITY(CM-3)
 *       D(3) - N2 NUMBER DENSITY(CM-3)
 *       D(4) - O2 NUMBER DENSITY(CM-3)
 *       D(5) - AR NUMBER DENSITY(CM-3)
 *       D(6) - TOTAL MASS DENSITY(GM/CM3)
 *       D(7) - H NUMBER DENSITY(CM-3)
 *       D(8) - N NUMBER DENSITY(CM-3)
 *       D(9) - HOT O NUMBER DENSITY(CM-3)
 *       T(1) - EXOSPHERIC TEMPERATURE
 *       T(2) - TEMPERATURE AT ALT
 */
void GHP7( int IYD, double SEC, double ALT, double GLAT, double GLONG, double STL, double F107A, double F107, double *AP, double *D, double *T, double PRESS, Lgm_Msis00Info *p ) {


    int     L, done, IDAY;
    double  PL, ZI, CL, CL2, CD, CA, Z, g, g2, XN, P, DIFF, XM, G, SH;

/*
    static double BM    = 1.3806E-19;
    static double TEST  = 0.00043;
    static double LTEST = 12;
*/

    double BM    = 1.3806E-19;
    double TEST  = 0.00043;
    double LTEST = 12;

    PL = log10( PRESS );

    // Initial altitude estimate
    if ( PL >= -5.0) {


        if (  PL > 2.5 )                     { ZI = 18.06*(3.00-PL); }
        else if ( (PL >  0.75) && (PL <=  2.5)  ) { ZI = 14.98*(3.08-PL); }
        else if ( (PL > -1.0)  && (PL <=  0.75) ) { ZI = 17.8*(2.72-PL); }
        else if ( (PL > -2.0)  && (PL <= -1.0)  ) { ZI = 14.28*(3.64-PL); }
        else if ( (PL > -4.0)  && (PL <= -2.0)  ) { ZI = 12.72*(4.32-PL); }
        //else if (  PL <= -4.0 )                   { ZI = 25.3*(0.11-PL); }
        else                                      { ZI = 25.3*(0.11-PL); }

        //IDAY = MOD(IYD,1000)
        IDAY = AMOD(IYD,1000);
        CL  = GLAT/90.0;
        CL2 = CL*CL;
        if ( IDAY <  182 ) CD = 1.0 - IDAY/91.25;
        if ( IDAY >= 182 ) CD = IDAY/91.25 - 3.0;

        CA = 0;
        if (PL > -1.11 && (PL <= -0.23)) CA = 1.0;
        if (PL > -0.23)                  CA = (2.79-PL)/(2.79+.23);
        if (PL <= -1.11 && (PL > -3.0) ) CA = (-2.93-PL)/(-2.93+1.11);

        Z = ZI - 4.87*CL*CD*CA - 1.64*CL2*CA + .31*CA*CL;

    } else {

        g = PL+4.0; g2 = g*g;
        Z = 22.0*g2 + 110.0;

    }

    //  ITERATION LOOP
    L = 0;
    
    done = 0;
    while ( !done ) {
        ++L;
        GTD7( IYD, SEC, Z, GLAT, GLONG, STL, F107A, F107, AP, 48, D, T, p );
        XN = D[1] + D[2] + D[3] + D[4] + D[5] + D[7] + D[8];
        P = BM*XN*T[2];
        if ( p->IMR == 1) P *= 1.e-6;
        DIFF = PL - log10( P );
        if ( (fabs(DIFF) < TEST) || (L == LTEST) ) break;
        XM = D[6]/(XN*1.66E-24);
        if ( p->IMR == 1 ) XM *= 1.0e3;
        g  = 1.0 + Z/p->RE; g2 = g*g;
        G  = p->GSURF/g2;
        SH = RGAS*T[2]/(XM*G);
        // New altitude estimate using scale height
        Z -= ( L < 6 )  ? SH*DIFF*2.302 : SH*DIFF;
    }


    if ( L == LTEST ) printf("GHP7 NOT CONVERGING FOR PRESS %12.2lf %12.2lf\n", PRESS, DIFF );
    ALT=Z;
    
    return;
}






















/*
 *     CALCULATE LATITUDE VARIABLE GRAVITY (GV) AND EFFECTIVE
 *     RADIUS (REFF)
 */
void GLATF( double LAT, double *GV, double *REFF ) {

    double          C2;
    //static double   DGTR = 1.74533e-2;
    double   DGTR = 1.74533e-2;

    C2 = cos( 2.0*DGTR*LAT );
    *GV = 980.616*(1.0 - 0.0026373*C2);
    *REFF = 2.0*(*GV)/(3.085462e-6 + 2.27e-9*C2)*1.0e-5;
    return;
}



/*
 *      Test if geophysical variables or switches changed and save
 *      Return 0 if unchanged and 1 if changed
 */
int VTST7( double IYD, double SEC, double GLAT, double GLONG, double STL, double F107A, double F107, double *AP, int IC, Lgm_Msis00Info *p ) {

    int fVTST7, Flag, I;
//return(1);

    fVTST7 = 0;

    Flag = 0;
    if      ( IYD   != p->IYDL[IC]  ) Flag = 1;
    else if ( SEC   != p->SECL[IC]  ) Flag = 1;
    else if ( GLAT  != p->GLATL[IC] ) Flag = 1;
    else if ( GLONG != p->GLL[IC]   ) Flag = 1;
    else if ( STL   != p->STLL[IC]  ) Flag = 1;
    else if ( F107A != p->FAL[IC]   ) Flag = 1;
    else if ( F107  != p->FL[IC]    ) Flag = 1;
    else {
        for ( I=1; I<=7; I++ ) {
            if ( AP[I]  != p->APL[I][IC]  ) {
                Flag = 1;
                break;
            }
        }
        for ( I=1; I<=25; I++ ) {
            if ( p->SW[I]  != p->SWL[I][IC]  ) {
                Flag = 1;
                break;
            }
            if ( p->SWC[I] != p->SWCL[I][IC] ) {
                Flag = 1;
                break;
            }
        }
    }

    if ( Flag == 1 ) {
        fVTST7    = 1;
        p->IYDL[IC]  = IYD;
        p->SECL[IC]  = SEC;
        p->GLATL[IC] = GLAT;
        p->GLL[IC]   = GLONG;
        p->STLL[IC]  = STL;
        p->FAL[IC]   = F107A;
        p->FL[IC]    = F107;
        for ( I=1; I<=7; I++ ) {
            p->APL[I][IC] = AP[I];
        }
        for ( I=1; I<=25; I++ ) {
            p->SWL[I][IC]  = p->SW[I];
            p->SWCL[I][IC] = p->SWC[I];
        }
    }

    return( fVTST7 );
}



/*
 *
 *    Thermospheric portion of NRLMSISE-00
 *    See GTD7 for more extensive comments
 *
 *       OUTPUT IN M-3 and KG/M3:   CALL METERS(.TRUE.)
 *
 *    INPUT VARIABLES:
 *       IYD - YEAR AND DAY AS YYDDD (day of year from 1 to 365 (or 366))
 *             (Year ignored in current model)
 *       SEC - UT(SEC)
 *       ALT - ALTITUDE(KM) (>72.5 km)
 *       GLAT - GEODETIC LATITUDE(DEG)
 *       GLONG - GEODETIC LONGITUDE(DEG)
 *       STL - LOCAL APPARENT SOLAR TIME(HRS; see Note below)
 *       F107A - 81 day AVERAGE OF F10.7 FLUX (centered on day DDD)
 *       F107 - DAILY F10.7 FLUX FOR PREVIOUS DAY
 *       AP - MAGNETIC INDEX(DAILY) OR WHEN SW(9)=-1. :
 *          - ARRAY CONTAINING:
 *            (1) DAILY AP
 *            (2) 3 HR AP INDEX FOR CURRENT TIME
 *            (3) 3 HR AP INDEX FOR 3 HRS BEFORE CURRENT TIME
 *            (4) 3 HR AP INDEX FOR 6 HRS BEFORE CURRENT TIME
 *            (5) 3 HR AP INDEX FOR 9 HRS BEFORE CURRENT TIME
 *            (6) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 12 TO 33 HRS PRIOR
 *                   TO CURRENT TIME
 *            (7) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 36 TO 57 HRS PRIOR
 *                   TO CURRENT TIME
 *       MASS - MASS NUMBER (ONLY DENSITY FOR SELECTED GAS IS
 *                CALCULATED.  MASS 0 IS TEMPERATURE.  MASS 48 FOR ALL.
 *                MASS 17 IS Anomalous O ONLY.)
 *
 *    NOTES ON INPUT VARIABLES: 
 *       UT, Local Time, and Longitude are used independently in the
 *       model and are not of equal importance for every situation.  
 *       For the most physically realistic calculation these three
 *       variables should be consistent (STL=SEC/3600+GLONG/15).
 *       The Equation of Time departures from the above formula
 *       for apparent local time can be included if available but
 *       are of minor importance.
 *
 *       F107 and F107A values used to generate the model correspond
 *       to the 10.7 cm radio flux at the actual distance of the Earth
 *       from the Sun rather than the radio flux at 1 AU. The following
 *       site provides both classes of values:
 *       ftp://ftp.ngdc.noaa.gov/STP/SOLAR_DATA/SOLAR_RADIO/FLUX/
 *
 *       F107, F107A, and AP effects are neither large nor well
 *       established below 80 km and these parameters should be set to
 *       150., 150., and 4. respectively.
 *
 *    OUTPUT VARIABLES:
 *       D(1) - HE NUMBER DENSITY(CM-3)
 *       D(2) - O NUMBER DENSITY(CM-3)
 *       D(3) - N2 NUMBER DENSITY(CM-3)
 *       D(4) - O2 NUMBER DENSITY(CM-3)
 *       D(5) - AR NUMBER DENSITY(CM-3)                       
 *       D(6) - TOTAL MASS DENSITY(GM/CM3) [Anomalous O NOT included]
 *       D(7) - H NUMBER DENSITY(CM-3)
 *       D(8) - N NUMBER DENSITY(CM-3)
 *       D(9) - Anomalous oxygen NUMBER DENSITY(CM-3)
 *       T(1) - EXOSPHERIC TEMPERATURE
 *       T(2) - TEMPERATURE AT ALT
 */
void GTS7( double IYD, double SEC, double ALT, double GLAT, double GLONG, double STL, double F107A, double F107, double *AP, double MASS, double *D, double *T, Lgm_Msis00Info *p ) {

/*
    static int      MN1     = 5;
    static double   DGTR    = 1.74533e-2;
    static double   DR      = 1.72142e-2;
    static int      MT[]    = { -9999, 48, 0, 4, 16, 28, 32, 40, 1, 49, 14, 17 };
    static double   ALTL[]  = { -1e99, 200.0, 300.0, 160.0, 250.0, 240.0, 450.0, 320.0, 450.0 };
    static double   ZN1[]   = { -1e99, 120.0, 110.0, 100.0, 90.0, 72.5 };
    static double   ALPHA[] = { -1e99, -0.38, 0.0, 0.0, 0.0, 0.17, 0.0, -0.38, 0.0, 0.0 };
*/

    int      MN1     = 5;
    double   DGTR    = 1.74533e-2;
    double   DR      = 1.72142e-2;
    int      MT[]    = { -9999, 48, 0, 4, 16, 28, 32, 40, 1, 49, 14, 17 };
    double   ALTL[]  = { -1e99, 200.0, 300.0, 160.0, 250.0, 240.0, 450.0, 320.0, 450.0 };
    double   ZN1[]   = { -1e99, 120.0, 110.0, 100.0, 90.0, 72.5 };
    double   ALPHA[] = { -1e99, -0.38, 0.0, 0.0, 0.0, 0.17, 0.0, -0.38, 0.0, 0.0 };

    int     I, J, ValidMass;
double TINF;
double  V2, YRD, S, g, g2, G28, DAY, ZHF, XMM, Z;
double  ZH28, ZHM28, XMD, B28, TZ, G4, ZH04, B04, ZHM04, ZC04, HC04;
double  G16, ZH16, B16, ZHM16, HC16, ZC16, HC216, HCC16, ZCC16, RC16, G32, ZH32;
double  B32, ZHM32, HC32, ZC32, HCC32, HCC232, ZCC32, RC32, G40, ZH40, B40, ZHM40, HC40, ZC40;
double  G1, ZH01, B01, ZHM01, HC01, ZC01, HCC01, ZCC01, RC01, G14, ZH14, B14, ZHM14, HC14, ZC14, HCC14, ZCC14, RC14;
double  G16H, DB16H, THO, T2, ZSHT, ZMHO, ZSHO;
//    double  DDUM;





    /*
     * Test for changed input
     */
    V2 = VTST7( IYD, SEC, GLAT, GLONG, STL, F107A, F107, AP, 2, p );
V2 =1;

    YRD    = IYD;
    p->ZA     = p->PDL[16][2];
    ZN1[1] = p->ZA;
    for ( J=1; J<=9; J++ ) D[J] = 0.0;

    /*
     * TINF VARIATIONS NOT IMPORTANT BELOW ZA OR ZN1(1)
     */
    if ( ALT > ZN1[1] ) {
        if( (V2 == 1) || (p->GTS7_ALAST <= ZN1[1]) ) TINF = p->PTM[1]*p->PT[1]*( 1.0 + p->SW[16]*GLOBE7( YRD, SEC, GLAT, GLONG, STL, F107A, F107, AP, p->PT, p ) );
    } else {
        TINF = p->PTM[1]*p->PT[1];
    }
    T[1] = TINF;

    /*
     * GRADIENT VARIATIONS NOT IMPORTANT BELOW ZN1(5)
     */
    if ( ALT > ZN1[5] ) {
        if( (V2 == 1) || (p->GTS7_ALAST <= ZN1[5]) ) p->G0 = p->PTM[4]*p->PS[1]*(1.0 + p->SW[19]*GLOBE7( YRD, SEC, GLAT, GLONG, STL, F107A, F107, AP, p->PS, p ));
    } else {
        p->G0 = p->PTM[4]*p->PS[1];
    }


    /*
     * Calculate these temperatures only if input changed
     */
    if ( (V2 == 1) || (ALT < 300.0) ) {
        p->TLB = p->PTM[2]*(1.0 + p->SW[17]*GLOBE7( YRD, SEC, GLAT, GLONG, STL, F107A, F107, AP, p->PD_rc[4], p ))*p->PD[1][4];
        S = p->G0/(TINF-p->TLB);
    }


    /*
     *   Lower thermosphere temp variations not significant for
     *   density above 300 km
     */
    if ( ALT < 300.0 ) {
        if ( (V2 == 1.0) || (p->GTS7_ALAST >= 300.0) ) {
            p->TN1[2]  = p->PTM[7]*p->PTL[1][1]/(1.0 - p->SW[18]*GLOB7S( p->PTL_rc[1], p ));
            p->TN1[3]  = p->PTM[3]*p->PTL[1][2]/(1.0 - p->SW[18]*GLOB7S( p->PTL_rc[2], p ));
            p->TN1[4]  = p->PTM[8]*p->PTL[1][3]/(1.0 - p->SW[18]*GLOB7S( p->PTL_rc[3], p ));
            p->TN1[5]  = p->PTM[5]*p->PTL[1][4]/(1.0 - p->SW[18]*p->SW[20]*GLOB7S( p->PTL_rc[4], p ));
            g = p->PTM[5]*p->PTL[1][4]; g2 = g*g;
            p->TGN1[2] = p->PTM[9]*p->PMA[1][9]*(1.0 + p->SW[18]*p->SW[20]*GLOB7S( p->PMA_rc[9], p ))*p->TN1[5]*p->TN1[5]/g2;
        }
    } else {
        p->TN1[2] = p->PTM[7]*p->PTL[1][1];
        p->TN1[3] = p->PTM[3]*p->PTL[1][2];
        p->TN1[4] = p->PTM[8]*p->PTL[1][3];
        p->TN1[5] = p->PTM[5]*p->PTL[1][4];
        g = p->PTM[5]*p->PTL[1][4]; g2 = g*g;
        p->TGN1[2]=p->PTM[9]*p->PMA[1][9]*p->TN1[5]*p->TN1[5]/g2;
    }



    // Are these used at all?
    p->Z0   = ZN1[4];
    p->T0   = p->TN1[4];
    p->TR12 = 1.0;

    if ( MASS != 0 ) {

        // N2 variation factor at Zlb
        G28 = p->SW[21]*GLOBE7( YRD, SEC, GLAT, GLONG, STL, F107A, F107, AP, p->PD_rc[3], p );
        DAY = AMOD(YRD,1000);

        // VARIATION OF TURBOPAUSE HEIGHT
        ZHF  = p->PDL[25][2]*(1.0 + p->SW[5]*p->PDL[25][1]*sin(DGTR*GLAT)*cos(DR*(DAY-p->PT[14])));
        YRD  = IYD;
        T[1] = TINF;
        XMM  = p->PDM[5][3];
        Z    = ALT;

        ValidMass = 0;
        for ( J=1; J<=11; J++ ) {
            if ( MASS == MT[J] ) {
                ValidMass = 1;
                break;
            }
        }

        if ( !ValidMass ) { // Invalid Mass

            printf( "MASS %g NOT VALID\n", MASS );

        } else {            // Valid Mass

            if ( (Z <= ALTL[6]) || (MASS == 28) || (MASS == 48) ) {

                /*
                 * **** N2 DENSITY ****
                 */

                // Diffusive density at Zlb
                p->DB28 = p->PDM[1][3]*exp(G28)*p->PD[1][3];

                // Diffusive density at Alt


                D[3] = DENSU( Z, p->DB28, TINF, p->TLB, 28.0, ALPHA[3], &T[2], p->PTM[6], S, MN1, ZN1, p->TN1,p->TGN1, p);
                p->DD   = D[3];

                // Turbopause
                ZH28  = p->PDM[3][3]*ZHF;
                ZHM28 = p->PDM[4][3]*p->PDL[6][2];
                XMD   = 28.0 - XMM;

                // Mixed density at Zlb
                B28 = DENSU( ZH28, p->DB28, TINF, p->TLB, XMD, ALPHA[3]-1.0, &TZ, p->PTM[6], S, MN1, ZN1, p->TN1, p->TGN1, p );
                if ( (Z <= ALTL[3]) && (p->SW[15] != 0.0) ) {

                    // Mixed density at Alt
                    p->DM28 = DENSU( Z, B28, TINF, p->TLB, XMM, ALPHA[3], &TZ, p->PTM[6], S, MN1, ZN1, p->TN1, p->TGN1, p );

                    // Net density at Alt
                    D[3] = DNET( D[3], p->DM28, ZHM28, XMM, 28.0 );

                }

            }


            switch ( J ) {

                case 1:
                case 3:

                        /*
                         *      **** HE DENSITY ****
                         */

                        // Density variation factor at Zlb
                        G4 = p->SW[21]*GLOBE7( YRD, SEC, GLAT, GLONG, STL, F107A, F107, AP, p->PD_rc[1], p );

                        // Diffusive density at Zlb
                        p->DB04 = p->PDM[1][1]*exp(G4)*p->PD[1][1];

                        // Diffusive density at Alt
                        D[1] = DENSU( Z, p->DB04, TINF, p->TLB, 4.0, ALPHA[1], &T[2], p->PTM[6], S, MN1, ZN1, p->TN1, p->TGN1, p );
                        p->DD   = D[1];

                        if ( (Z <= ALTL[1]) && (p->SW[15] != 0.0) ) {

                            // Turbopause
                            ZH04 = p->PDM[3][1];

                            // Mixed density at Zlb
                            B04 = DENSU( ZH04, p->DB04, TINF, p->TLB, 4.0-XMM , ALPHA[1]-1.0, &T[2], p->PTM[6], S, MN1, ZN1, p->TN1, p->TGN1, p );

                            // Mixed density at Alt
                            p->DM04  = DENSU( Z, B04, TINF, p->TLB, XMM, 0.0, &T[2], p->PTM[6], S, MN1, ZN1, p->TN1, p->TGN1, p );
                            ZHM04 = ZHM28;

                            // Net density at Alt
                            D[1] = DNET( D[1], p->DM04, ZHM04, XMM, 4.0 );

                            // Correction to specified mixing ratio at ground
                            p->RL   = log( B28*p->PDM[2][1]/B04 );
                            ZC04 = p->PDM[5][1]*p->PDL[1][2];
                            HC04 = p->PDM[6][1]*p->PDL[2][2];

                            // Net density corrected at Alt
                            D[1] *= CCOR( Z, p->RL, HC04, ZC04 );

                        }
                        if ( MASS != 48) break;

                case 4:
                case 9:
                        /*
                         *    **** O DENSITY ****
                         */

                        // Density variation factor at Zlb
                        G16 = p->SW[21]*GLOBE7( YRD, SEC, GLAT, GLONG, STL, F107A, F107, AP, p->PD_rc[2], p );

                        // Diffusive density at Zlb
                        p->DB16 =  p->PDM[1][2]*exp(G16)*p->PD[1][2];

                        // Diffusive density at Alt
                        D[2] = DENSU( Z, p->DB16, TINF, p->TLB, 16.0, ALPHA[2], &T[2], p->PTM[6], S, MN1, ZN1, p->TN1, p->TGN1, p );
                        p->DD   = D[2];

                        if ( (Z <= ALTL[2]) && (p->SW[15] != 0.0) ) {

                            //  Corrected from PDM(3,1) to PDM(3,2)  12/2/85
                            //  Turbopause
                            ZH16 = p->PDM[3][2];

                            //  Mixed density at Zlb
                            B16 = DENSU( ZH16, p->DB16, TINF, p->TLB, 16-XMM, ALPHA[2]-1.0, &T[2], p->PTM[6], S, MN1, ZN1, p->TN1, p->TGN1, p );

                            //  Mixed density at Alt
                            p->DM16  = DENSU( Z, B16, TINF, p->TLB, XMM, 0.0, &T[2], p->PTM[6], S, MN1, ZN1, p->TN1, p->TGN1, p );
                            ZHM16 = ZHM28;

                            //  Net density at Alt
                            D[2] = DNET( D[2], p->DM16, ZHM16, XMM, 16.0 );

                            //  3/16/99 Change form to match O2 departure from diff equil near 150
                            //  km and add dependence on F10.7
                            //  RL=ALOG(B28*PDM(2,2)*fabs(PDL(17,2))/B16)
                            p->RL = p->PDM[2][2]*p->PDL[17][2]*( 1.0 + p->SW[1]*p->PDL[24][1]*(F107A-150.0) );
                            HC16  = p->PDM[6][2]*p->PDL[4][2];
                            ZC16  = p->PDM[5][2]*p->PDL[3][2];
                            HC216 = p->PDM[6][2]*p->PDL[5][2];
                            D[2]  *= CCOR2( Z, p->RL, HC16, ZC16, HC216 );

                            //  Chemistry correction
                            HCC16 = p->PDM[8][2]*p->PDL[14][2];
                            ZCC16 = p->PDM[7][2]*p->PDL[13][2];
                            RC16  = p->PDM[4][2]*p->PDL[15][2];

                            //  Net density corrected at Alt
                            D[2] *= CCOR( Z, RC16, HCC16, ZCC16 );

                        }
                        if ( (MASS != 48) && (MASS != 49) ) break;

                case 6:
                        /*
                         *      **** O2 DENSITY ****
                         */

                        //  Density variation factor at Zlb
                        G32= p->SW[21]*GLOBE7( YRD, SEC, GLAT, GLONG, STL, F107A, F107, AP, p->PD_rc[5], p );

                        //  Diffusive density at Zlb
                        p->DB32 = p->PDM[1][4]*exp(G32)*p->PD[1][5];

                        //  Diffusive density at Alt
                        D[4] = DENSU( Z, p->DB32, TINF, p->TLB, 32.0, ALPHA[4], &T[2], p->PTM[6], S, MN1, ZN1, p->TN1, p->TGN1, p );
                        p->DD = ( MASS == 49 ) ? p->DD + 2.0*D[4] : D[4];

                        if ( p->SW[15] != 0.0 ) {

                            if ( Z <= ALTL[4] ) {
                                //  Turbopause
                                ZH32 = p->PDM[3][4];

                                //  Mixed density at Zlb
                                B32 = DENSU( ZH32, p->DB32, TINF, p->TLB, 32.0-XMM, ALPHA[4]-1.0, &T[2], p->PTM[6], S, MN1, ZN1, p->TN1, p->TGN1, p );

                                //  Mixed density at Alt
                                p->DM32  = DENSU( Z, B32, TINF, p->TLB, XMM, 0.0, &T[2], p->PTM[6], S, MN1, ZN1, p->TN1, p->TGN1, p );
                                ZHM32 = ZHM28;

                                //  Net density at Alt
                                D[4] = DNET( D[4], p->DM32, ZHM32, XMM, 32.0 );

                                //  Correction to specified mixing ratio at ground
                                p->RL   = log( B28*p->PDM[2][4]/B32 );
                                HC32 = p->PDM[6][4]*p->PDL[8][2];
                                ZC32 = p->PDM[5][4]*p->PDL[7][2];
                                D[4] *= CCOR( Z, p->RL, HC32, ZC32 );
                            }

                            //  Correction for general departure from diffusive equilibrium above Zlb
                            HCC32  = p->PDM[8][4]*p->PDL[23][2];
                            HCC232 = p->PDM[8][4]*p->PDL[23][1];
                            ZCC32  = p->PDM[7][4]*p->PDL[22][2];
                            RC32   = p->PDM[4][4]*p->PDL[24][2]*( 1.0 + p->SW[1]*p->PDL[24][1]*(F107A-150.0) );

                            //  Net density corrected at Alt
                            D[4] *= CCOR2( Z, RC32, HCC32, ZCC32, HCC232 );

                        }
                        if ( MASS != 48 ) break;

                case 7:
                        /*
                         *      **** AR DENSITY ****
                         */

                        // Density variation factor at Zlb
                        G40= p->SW[21]*GLOBE7( YRD, SEC, GLAT, GLONG, STL, F107A, F107, AP, p->PD_rc[6], p );

                        // Diffusive density at Zlb
                        p->DB40 = p->PDM[1][5]*exp(G40)*p->PD[1][6];

                        // Diffusive density at Alt
                        D[5] = DENSU( Z, p->DB40, TINF, p->TLB, 40.0 , ALPHA[5], &T[2], p->PTM[6], S, MN1, ZN1, p->TN1, p->TGN1, p );
                        p->DD   = D[5];

                        if ( (Z <= ALTL[5]) && (p->SW[15] != 0.0) ) {

                            // Turbopause
                            ZH40 = p->PDM[3][5];

                            // Mixed density at Zlb
                            B40 = DENSU( ZH40, p->DB40, TINF, p->TLB, 40.0-XMM, ALPHA[5]-1., &T[2], p->PTM[6], S, MN1, ZN1, p->TN1, p->TGN1, p );

                            // Mixed density at Alt
                            p->DM40  = DENSU( Z, B40, TINF, p->TLB, XMM, 0.0, &T[2], p->PTM[6], S, MN1, ZN1, p->TN1, p->TGN1, p );
                            ZHM40 = ZHM28;

                            // Net density at Alt
                            D[5] = DNET( D[5], p->DM40, ZHM40, XMM, 40.0 );

                            // Correction to specified mixing ratio at ground
                            p->RL   = log( B28*p->PDM[2][5]/B40 );
                            HC40 = p->PDM[6][5]*p->PDL[10][2];
                            ZC40 = p->PDM[5][5]*p->PDL[9][2];

                            // Net density corrected at Alt
                            D[5] *= CCOR( Z, p->RL, HC40, ZC40 );

                        }
                        if ( MASS != 48 ) break;

                case 8:
                        /*
                         *       **** HYDROGEN DENSITY ****
                         */

                        //  Density variation factor at Zlb
                        G1 = p->SW[21]*GLOBE7( YRD, SEC, GLAT, GLONG, STL, F107A, F107, AP, p->PD_rc[7], p );

                        //  Diffusive density at Zlb
                        p->DB01 = p->PDM[1][6]*exp(G1)*p->PD[1][7];

                        //  Diffusive density at Alt
                        D[7] = DENSU(Z, p->DB01, TINF, p->TLB, 1.0, ALPHA[7], &T[2], p->PTM[6], S, MN1, ZN1, p->TN1, p->TGN1, p);
                        p->DD   = D[7];

                        if ( (Z <= ALTL[7]) && (p->SW[15] != 0.0) ) {

                            // Turbopause
                            ZH01 = p->PDM[3][6];

                            // Mixed density at Zlb
                            B01 = DENSU( ZH01, p->DB01, TINF, p->TLB, 1.0-XMM, ALPHA[7]-1., &T[2], p->PTM[6], S, MN1, ZN1, p->TN1, p->TGN1, p);

                            // Mixed density at Alt
                            p->DM01  = DENSU(Z, B01, TINF, p->TLB, XMM, 0., &T[2], p->PTM[6], S, MN1, ZN1, p->TN1, p->TGN1, p);
                            ZHM01 = ZHM28;

                            // Net density at Alt
                            D[7] = DNET( D[7], p->DM01, ZHM01, XMM, 1.0 );

                            // Correction to specified mixing ratio at ground
                            p->RL   =  log( B28*p->PDM[2][6]*fabs(p->PDL[18][2])/B01 );
                            HC01 =  p->PDM[6][6]*p->PDL[12][2];
                            ZC01 =  p->PDM[5][6]*p->PDL[11][2];
                            D[7] *= CCOR( Z, p->RL, HC01, ZC01 );

                            // Chemistry correction
                            HCC01 = p->PDM[8][6]*p->PDL[20][2];
                            ZCC01 = p->PDM[7][6]*p->PDL[19][2];
                            RC01  = p->PDM[4][6]*p->PDL[21][2];

                            // Net density corrected at Alt
                            D[7] *= CCOR( Z, RC01, HCC01, ZCC01 );

                        }
                        if ( MASS != 48 ) break;

                case 10:
                        /*
                         *  **** ATOMIC NITROGEN DENSITY ****
                         */

                        //  Density variation factor at Zlb
                        G14 = p->SW[21]*GLOBE7( YRD, SEC, GLAT, GLONG, STL, F107A, F107, AP, p->PD_rc[8], p );

                        //  Diffusive density at Zlb
                        p->DB14 = p->PDM[1][7]*exp(G14)*p->PD[1][8];

                        //  Diffusive density at Alt
                        D[8] = DENSU( Z, p->DB14, TINF, p->TLB, 14., ALPHA[8], &T[2], p->PTM[6], S, MN1, ZN1, p->TN1, p->TGN1, p );
                        p->DD   = D[8];

                        if ( (Z <= ALTL[8]) && (p->SW[15] != 0.0) ) {

                            // Turbopause
                            ZH14 = p->PDM[3][7];

                            // Mixed density at Zlb
                            B14 = DENSU( ZH14, p->DB14, TINF, p->TLB, 14.0-XMM, ALPHA[8]-1.0, &T[2], p->PTM[6], S, MN1, ZN1, p->TN1, p->TGN1, p );

                            // Mixed density at Alt
                            p->DM14  = DENSU( Z, B14, TINF, p->TLB, XMM, 0.0, &T[2], p->PTM[6], S, MN1, ZN1, p->TN1, p->TGN1, p );
                            ZHM14 = ZHM28;

                            // Net density at Alt
                            D[8] = DNET( D[8], p->DM14, ZHM14, XMM, 14.0 );

                            // Correction to specified mixing ratio at ground
                            p->RL   = log( B28*p->PDM[2][7]*fabs(p->PDL[3][1])/B14 );
                            HC14 = p->PDM[6][7]*p->PDL[2][1];
                            ZC14 = p->PDM[5][7]*p->PDL[1][1];
                            D[8] *= CCOR( Z, p->RL, HC14, ZC14 );

                            // Chemistry correction
                            HCC14 = p->PDM[8][7]*p->PDL[5][1];
                            ZCC14 = p->PDM[7][7]*p->PDL[4][1];
                            RC14  = p->PDM[4][7]*p->PDL[6][1];

                            // Net density corrected at Alt
                            D[8] *= CCOR( Z, RC14, HCC14, ZCC14 );

                        }
                        if ( MASS != 48 ) break;


                case 11:
                        /* 
                         *  **** Anomalous OXYGEN DENSITY ****
                         */
                        G16H  = p->SW[21]*GLOBE7( YRD, SEC, GLAT, GLONG, STL, F107A, F107, AP, p->PD_rc[9], p );
                        DB16H = p->PDM[1][8]*exp(G16H)*p->PD[1][9];
                        THO   = p->PDM[10][8]*p->PDL[7][1];
                        p->DD    = DENSU( Z, DB16H, THO, THO, 16.0, ALPHA[9], &T2, p->PTM[6], S, MN1, ZN1, p->TN1, p->TGN1, p );
                        ZSHT  = p->PDM[6][8];
                        ZMHO  = p->PDM[5][8];
                        ZSHO  = SCALH( ZMHO, 16.0, THO, p );
                        D[9]  = p->DD*exp( -ZSHT/ZSHO*(exp(-(Z-ZMHO)/ZSHT)-1.0) );
                        if ( MASS != 48 ) break;

                default:
                        /*
                         *      TOTAL MASS DENSITY
                         */
                        D[6] = 1.66E-24*( 4.0*D[1] + 16.0*D[2] + 28.0*D[3] + 32.0*D[4] + 40.0*D[5] + D[7] + 14.0*D[8] );
                        p->DB48 = 1.66E-24*( 4.0*p->DB04 + 16.0*p->DB16 + 28.0*p->DB28 + 32.0*p->DB32 + 40.0*p->DB40 + p->DB01 + 14.0*p->DB14 );
//NOT USED ANYWHERE????
                        break;

            } // end switch

        } // end valid mass

    } else {  // case for MASS = 0

        // TEMPERATURE AT ALTITUDE
        Z = fabs(ALT);
//        DDUM  = DENSU( Z, 1.0, TINF, p->TLB, 0.0, 0.0, &T[2], p->PTM[6], S, MN1, ZN1, p->TN1, p->TGN1, p );

    } 


    // ADJUST DENSITIES FROM CGS TO KGM
    if ( p->IMR == 1 ) {
        for ( I=1; I<=9; I++ ) {
            D[I] *= 1.0e6;
        }
        D[6] /= 1000.0;
    }
    p->GTS7_ALAST = ALT;

    return;

}

















/*
 * Convert outputs to Kg & Meters if METER true
 */
void METERS( int METER, Lgm_Msis00Info *p ) {
    p->IMR = ( METER ) ? 1 : 0;
    return;
}




/*
 * Calculate scale height (km)
 */
double SCALH( double ALT, double XM, double TEMP, Lgm_Msis00Info *p ) {
    double  q, q2, G;
//    static double   RGAS = 831.4;
    q = 1.0+ALT/p->RE; q2 = q*q;
    G = p->GSURF/q2;
    return( RGAS*TEMP/(G*XM) );
}



#define G0(A) ( (A-4.0+(P[26]-1.0)*(A-4.0+(exp(-fabs(P[25])*(A-4.0))-1.0)/fabs(P[25]))) )
double SG0( double EX, double *AP, double *P ) {
    double  F, EX2, EX3, EX4, EX8, EX12, EX16, EX19, SUMEX;
    EX2 = EX*EX; EX3 = EX2*EX; EX4 = EX2*EX2; EX8 = EX4*EX4; EX12 = EX8*EX4; EX16 = EX8*EX8; EX19 = EX16*EX3;
    SUMEX = 1.0 + (1.0 - EX19 )/(1.0-EX) * sqrt(EX);
    F = (G0(AP[2]) + (G0(AP[3])*EX + G0(AP[4])*EX2 + G0(AP[5])*EX3 + ( G0(AP[6])*EX4 + G0(AP[7])*EX12)*(1.-EX8)/(1.-EX)) )/SUMEX;
    return( F );
}


/*
 *      CALCULATE G(L) FUNCTION 
 *      Upper Thermosphere Parameters
 */
double GLOBE7( double YRD, double SEC, double LAT, double LONG, double TLOC, double F107A, double F107, double *AP, double *P, Lgm_Msis00Info *p ) {

/*
    static int    NSW   =  14;
    static double DGTR  =  1.74533e-2;
    static double DR    =  1.72142e-2;
    static double SW9   =  1.0;


    static double P14   = -1000.0;
    static double P18   = -1000.0;
    static double P32   = -1000.0;
    static double HR    =  0.2618;
    static double SR    =  7.2722e-5;
    static double P39   = -1000.0;
*/
    int    NSW   =  14;
    double DGTR  =  1.74533e-2;
    double DR    =  1.72142e-2;
    double SW9   =  1.0;


    double HR    =  0.2618;
    double SR    =  7.2722e-5;





    int     I, J;
    double  C, S, C2, C4, S2, F1, F2;
    double  T71, T72, T81, T82, P44, P45, EXP1;


    /*
     * 3hr Magnetic activity functions
     */
    // Eq. A24d  Originally a fortran statemen function -- changed to macro
    //G0(A) = (A-4.0+(P[26]-1.0)*(A-4.0+(exp(-fabs(P[25])*(A-4.0))-1.0)/fabs(P[25])));


    // next two functions, originally a fortran statemen function -- changed to functions
    // Eq. A24c Originally a fortran statemen function -- changed to functions
    //EX2 = EX*EX; EX3 = EX2*EX; EX4 = EX2*EX2; EX8 = EX4*EX4; EX12 = EX8*EX4; EX16 = EX8*EX8; EX19 = EX16*EX3;
    //SUMEX(EX) = 1.0 + (1.0 - EX19 )/(1.0-EX) * sqrt(EX);

    // Eq. A24a
    //SG0(EX)= ((G0(AP[2]) + (G0(AP[3])*EX + G0(AP[4])*EX2 + G0(AP[5])*EX3 + ( G0(AP[6])*EX4 + G0(AP[7])*EX12)*(1.-EX8)/(1.-EX)) )/SUMEX(EX));







    if ( p->ISW != 64999 ) TSELEC( p->SV, p );
    for ( J=1; J<=14; J++ ) p->T[J] = 0;

    if ( p->SW[9] > 0) SW9 =  1.0;
    if ( p->SW[9] < 0) SW9 = -1.0;
    p->IYR   = YRD/1000.0;
    p->DAY   = YRD - p->IYR*1000.0;
    p->LONG = LONG;

    // Eq. A22 (remainder of code)
    if ( p->XL != LAT ) {
        // CALCULATE LEGENDRE POLYNOMIALS
        C  = sin(LAT*DGTR);
        S  = cos(LAT*DGTR);
        C2 = C*C;
        C4 = C2*C2;
        S2 = S*S;
        p->PLG[2][1] = C;
        p->PLG[3][1] = 0.5*(3.*C2 -1.);
        p->PLG[4][1] = 0.5*(5.*C*C2-3.*C);
        p->PLG[5][1] = (35.*C4 - 30.*C2 + 3.)/8.;
        p->PLG[6][1] = (63.*C2*C2*C - 70.*C2*C + 15.*C)/8.;
        p->PLG[7][1] = (11.*C*p->PLG[6][1] - 5.*p->PLG[5][1])/6.;
        //p->PLG[8][1] = (13.*C*p->PLG[7][1] - 6.*p->PLG[6][1])/7.;
        p->PLG[2][2] = S;
        p->PLG[3][2] = 3.*C*S;
        p->PLG[4][2] = 1.5*(5.*C2-1.)*S;
        p->PLG[5][2] = 2.5*(7.*C2*C-3.*C)*S;
        p->PLG[6][2] = 1.875*(21.*C4 - 14.*C2 +1.)*S;
        p->PLG[7][2] = (11.*C*p->PLG[6][2]-6.*p->PLG[5][2])/5.;
        //p->PLG[8][2] = (13.*C*p->PLG[7][2]-7.*p->PLG[6][2])/6.;
        //p->PLG[9][2] = (15.*C*p->PLG[8][2]-8.*p->PLG[7][2])/7.;
        p->PLG[3][3] = 3.*S2;
        p->PLG[4][3] = 15.*S2*C;
        p->PLG[5][3] = 7.5*(7.*C2 -1.)*S2;
        p->PLG[6][3] = 3.*C*p->PLG[5][3]-2.*p->PLG[4][3];
        p->PLG[7][3] = (11.*C*p->PLG[6][3]-7.*p->PLG[5][3])/4.;
        p->PLG[8][3] = (13.*C*p->PLG[7][3]-8.*p->PLG[6][3])/5.;
        p->PLG[4][4] = 15.*S2*S;
        p->PLG[5][4] = 105.*S2*S*C ;
        p->PLG[6][4] = (9.*C*p->PLG[5][4]-7.*p->PLG[4][4])/2.;
        p->PLG[7][4] = (11.*C*p->PLG[6][4]-8.*p->PLG[5][4])/3.;
        p->XL = LAT;
    }



    if ( p->TLL != TLOC) {  
        if ( (p->SW[7] != 0) || (p->SW[8] != 0) || (p->SW[14] != 0) ) {
            p->STLOC  = sin( HR*TLOC );
            p->CTLOC  = cos( HR*TLOC );
            p->S2TLOC = sin( 2.*HR*TLOC );
            p->C2TLOC = cos( 2.*HR*TLOC );
            p->S3TLOC = sin( 3.*HR*TLOC );
            p->C3TLOC = cos( 3.*HR*TLOC );
            p->TLL = TLOC;
        }
    }



    if ( (p->DAY != p->DAYL) || (P[14] != p->GLOBE7_P14) ) p->GLOBE7_CD14 = cos(     DR*(p->DAY-P[14]) );
    if ( (p->DAY != p->DAYL) || (P[18] != p->GLOBE7_P18) ) p->GLOBE7_CD18 = cos( 2.0*DR*(p->DAY-P[18]) );
    if ( (p->DAY != p->DAYL) || (P[32] != p->GLOBE7_P32) ) p->GLOBE7_CD32 = cos(     DR*(p->DAY-P[32]) );
    if ( (p->DAY != p->DAYL) || (P[39] != p->GLOBE7_P39) ) p->GLOBE7_CD39 = cos( 2.0*DR*(p->DAY-P[39]) );
    p->DAYL        = p->DAY;
    p->GLOBE7_P14  = P[14];
    p->GLOBE7_P18  = P[18];
    p->GLOBE7_P32  = P[32];
    p->GLOBE7_P39  = P[39];

    // F10.7 EFFECT
    p->DF   = F107 - F107A;
    p->DFA  = F107A - 150.0;
    p->T[1] = P[20]*p->DF*(1.0+P[60]*p->DFA) + P[21]*p->DF*p->DF + P[22]*p->DFA + P[30]*p->DFA*p->DFA;
    F1   = 1.0 + (P[48]*p->DFA + P[20]*p->DF + P[21]*p->DF*p->DF )*p->SWC[1];
    F2   = 1.0 + (P[50]*p->DFA + P[20]*p->DF + P[21]*p->DF*p->DF )*p->SWC[1];


    // TIME INDEPENDENT
    p->T[2] = P[2]*p->PLG[3][1] + P[3]*p->PLG[5][1] + P[23]*p->PLG[7][1] 
                + (P[15]*p->PLG[3][1])*p->DFA*p->SWC[1] 
                    + P[27]*p->PLG[2][1];


    // SYMMETRICAL ANNUAL
    p->T[3] = P[19]*p->GLOBE7_CD32;


    // SYMMETRICAL SEMIANNUAL
    p->T[4] = ( P[16] + P[17]*p->PLG[3][1] )*p->GLOBE7_CD18;

    // ASYMMETRICAL ANNUAL
    p->T[5] = F1*( P[10]*p->PLG[2][1] + P[11]*p->PLG[4][1] )*p->GLOBE7_CD14;

    // ASYMMETRICAL SEMIANNUAL
    p->T[6] = P[38]*p->PLG[2][1]*p->GLOBE7_CD39;

    // DIURNAL
    if ( p->SW[7] != 0 ) {
      T71 = (P[12]*p->PLG[3][2])*p->GLOBE7_CD14*p->SWC[5];
      T72 = (P[13]*p->PLG[3][2])*p->GLOBE7_CD14*p->SWC[5];
      p->T[7] = F2* (  (P[4]*p->PLG[2][2] + P[5]*p->PLG[4][2] + P[28]*p->PLG[6][2] + T71)*p->CTLOC
                  + (P[7]*p->PLG[2][2] + P[8]*p->PLG[4][2] + P[29]*p->PLG[6][2] + T72)*p->STLOC );
    }


    // SEMIDIURNAL
    if ( p->SW[8] != 0 ) {
        T81  = ( P[24]*p->PLG[4][3]+P[36]*p->PLG[6][3] )*p->GLOBE7_CD14*p->SWC[5];
        T82  = ( P[34]*p->PLG[4][3]+P[37]*p->PLG[6][3] )*p->GLOBE7_CD14*p->SWC[5];
        p->T[8] = F2* (  (P[6]*p->PLG[3][3] + P[42]*p->PLG[5][3] + T81)*p->C2TLOC 
                    + (P[9]*p->PLG[3][3] + P[43]*p->PLG[5][3] + T82)*p->S2TLOC );
    }


    // TERDIURNAL
    if ( p->SW[14] != 0 ) {
        p->T[14] = F2*(  ( P[40]*p->PLG[4][4] + (P[94]*p->PLG[5][4] + P[47]*p->PLG[7][4])*p->GLOBE7_CD14*p->SWC[5] )*p->S3TLOC
                    + ( P[41]*p->PLG[4][4] + (P[95]*p->PLG[5][4] + P[49]*p->PLG[7][4])*p->GLOBE7_CD14*p->SWC[5] )*p->C3TLOC );
    }

    // MAGNETIC ACTIVITY BASED ON DAILY AP
    if ( SW9 != -1.0 ) {
        p->APD = AP[1]-4.0;
        P44 = P[44];
        P45 = P[45];
        if ( P44 < 0) P44 = 1.E-5;
        p->APDF = p->APD + (P45-1.0)*( p->APD + (exp(-P44*p->APD)-1.0)/P44 );
        if( p->SW[9] != 0) {
            p->T[9] = p->APDF*( P[33] + P[46]*p->PLG[3][1] + P[35]*p->PLG[5][1]
                        + (P[101]*p->PLG[2][1]+P[102]*p->PLG[4][1]+P[103]*p->PLG[6][1])*p->GLOBE7_CD14*p->SWC[5]
                        + (P[122]*p->PLG[2][2]+P[123]*p->PLG[4][2]+P[124]*p->PLG[6][2])*p->SWC[7]*cos(HR*(TLOC-P[125])) );
        }
    } else if ( P[52] != 0) {
        EXP1 = exp( -10800.0*fabs(P[52]) / (1.0+P[139]*(45.0-fabs(LAT))) );
        if ( EXP1 > 0.99999) EXP1 = 0.99999;
        if ( P[25] < 1.0e-4) P[25] = 1.e-4;
        p->APT[1] = SG0( EXP1, AP, P );
        //p->APT[2] = SG2(EXP1);
        //p->APT[3] = SG0( EXP2, AP, P );
        //p->APT[4] = SG2(EXP2);
        if ( p->SW[9] != 0) {
            p->T[9] = p->APT[1]*( P[51] + P[97]*p->PLG[3][1] + P[55]*p->PLG[5][1]
                            + (P[126]*p->PLG[2][1]+P[127]*p->PLG[4][1]+P[128]*p->PLG[6][1])*p->GLOBE7_CD14*p->SWC[5]
                            + (P[129]*p->PLG[2][2]+P[130]*p->PLG[4][2]+P[131]*p->PLG[6][2])*p->SWC[7]*cos(HR*(TLOC-P[132])) );
        }
    }





    if ( (p->SW[10] != 0) && (LONG > -1000.0) ) {

        // LONGITUDINAL
        if ( p->SW[11] != 0 ) {
            p->T[11] = ( 1.0 + P[81]*p->DFA*p->SWC[1] )*
                      ( (P[65]*p->PLG[3][2]+P[66]*p->PLG[5][2]+P[67]*p->PLG[7][2] +P[104]*p->PLG[2][2]+P[105]*p->PLG[4][2]+P[106]*p->PLG[6][2] +p->SWC[5]*(P[110]*p->PLG[2][2]+P[111]*p->PLG[4][2]+P[112]*p->PLG[6][2])*p->GLOBE7_CD14)*cos(DGTR*LONG)
                      + (P[91]*p->PLG[3][2]+P[92]*p->PLG[5][2]+P[93]*p->PLG[7][2] +P[107]*p->PLG[2][2]+P[108]*p->PLG[4][2]+P[109]*p->PLG[6][2] +p->SWC[5]*(P[113]*p->PLG[2][2]+P[114]*p->PLG[4][2]+P[115]*p->PLG[6][2])*p->GLOBE7_CD14)*sin(DGTR*LONG) );
        }

        // UT AND MIXED UT,LONGITUDE
        if ( p->SW[12] != 0) {
            p->T[12] = (1.0+P[96]*p->PLG[2][1]) * (1.0+P[82]*p->DFA*p->SWC[1]) * (1.0+P[120]*p->PLG[2][1]*p->SWC[5]*p->GLOBE7_CD14) * ((P[69]*p->PLG[2][1]+P[70]*p->PLG[4][1]+P[71]*p->PLG[6][1])*cos(SR*(SEC-P[72])));
            p->T[12] += p->SWC[11]*(P[77]*p->PLG[4][3]+P[78]*p->PLG[6][3]+P[79]*p->PLG[8][3]) * cos(SR*(SEC-P[80])+2.*DGTR*LONG)*(1.0+P[138]*p->DFA*p->SWC[1]);
        }

        // UT,LONGITUDE MAGNETIC ACTIVITY
        if ( (p->SW[13] != 0) && (SW9 != -1.0) ) {
            p->T[13] = p->APDF*p->SWC[11]*(1.+P[121]*p->PLG[2][1])* ((P[ 61]*p->PLG[3][2]+P[ 62]*p->PLG[5][2]+P[ 63]*p->PLG[7][2])* cos(DGTR*(LONG-P[ 64])))
                        + p->APDF*p->SWC[11]*p->SWC[5]* (P[116]*p->PLG[2][2]+P[117]*p->PLG[4][2]+P[118]*p->PLG[6][2])* p->GLOBE7_CD14*cos(DGTR*(LONG-P[119]))
                        + p->APDF*p->SWC[12]* (P[ 84]*p->PLG[2][1]+P[ 85]*p->PLG[4][1]+P[ 86]*p->PLG[6][1])* cos(SR*(SEC-P[ 76]));
        } else if ( P[52] != 0 ) {
            p->T[13] = p->APT[1]*p->SWC[11]*(1.+P[133]*p->PLG[2][1])* ((P[53]*p->PLG[3][2]+P[99]*p->PLG[5][2]+P[68]*p->PLG[7][2])* cos(DGTR*(LONG-P[98])))
                        +p->APT[1]*p->SWC[11]*p->SWC[5]* (P[134]*p->PLG[2][2]+P[135]*p->PLG[4][2]+P[136]*p->PLG[6][2])* p->GLOBE7_CD14*cos(DGTR*(LONG-P[137]))
                        +p->APT[1]*p->SWC[12]* (P[56]*p->PLG[2][1]+P[57]*p->PLG[4][1]+P[58]*p->PLG[6][1])* cos(SR*(SEC-P[59]));
        }
    }


    // PARMS NOT USED: 83, 90,100,140-150


//printf("GLOBE7: p->TINF = %g\n", p->TINF);
    //p->TINF = P[31];
double TINF;
    TINF = P[31];
    for ( I=1; I<=NSW; I++ ) {
        //p->TINF += fabs( p->SW[I] ) * p->T[I];
//printf("    A. GLOBE7: TINF, T[%d] = %g %g \n", I, TINF, p->T[I] );
        TINF += fabs( p->SW[I] ) * p->T[I];
    }

    //return( p->TINF );
    return( TINF );

}


/* 
 *       SET SWITCHES
 *       Output in  COMMON/CSW/SW(25),ISW,SWC(25)
 *       SW FOR MAIN TERMS, SWC FOR CROSS TERMS
 * 
 *       TO TURN ON AND OFF PARTICULAR VARIATIONS CALL TSELEC(SV),
 *       WHERE SV IS A 25 ELEMENT ARRAY CONTAINING 0. FOR OFF, 1. 
 *       FOR ON, OR 2. FOR MAIN EFFECTS OFF BUT CROSS TERMS ON
 *
 *       To get current values of SW: CALL TRETRV(SW)
 */
void TSELEC( double *SV, Lgm_Msis00Info *p ) {
    int I;
    for ( I=1; I<=25; I++ ) {
        p->SAV[I] = SV[I];
        p->SW[I]  = AMOD( SV[I], 2.0 );
        p->SWC[I] = ( (fabs(SV[I]) == 1) || (fabs(SV[I]) == 2.0) ) ? 1.0 : 0.0;
    }
    p->ISW = 64999;
    return;

}

void TRETRV( double *SVV, Lgm_Msis00Info *p ) {
    int I;
    for ( I=1; I<=25; I++ ) SVV[I] = p->SAV[I];
    return;
}





/*
 *     VERSION OF GLOBE FOR LOWER ATMOSPHERE 10/26/99
 */
double GLOB7S( double *P, Lgm_Msis00Info *p ) {

    int     I, J;
    double  TT, T[15];
    double  T71, T72, T81, T82;

/*
    static double DR   = 1.72142E-2;
    static double DGTR =  1.74533E-2;
    static double PSET =  2.0;
    static double P32  = -1000.0;
    static double P18  = -1000.0;
    static double P14  = -1000.0;
    static double P39  = -1000.0;
*/
    double DR   = 1.72142E-2;
    double DGTR =  1.74533E-2;
    double PSET =  2.0;


    // CONFIRM PARAMETER SET
    if ( P[100] == 0 ) P[100] = PSET;
    if ( P[100] != PSET) {
        printf("WRONG PARAMETER SET FOR GLOB7S %g %g\n", PSET, P[100] );
        exit( 1 );
    }

    for ( J=1; J<=14; J++ ) T[J] = 0.0;

    if ( (p->DAY != p->DAYL) || (p->GLOB7S_P14 != P[14]) ) p->GLOB7S_CD14 = cos(     DR*(p->DAY-P[14]) );
    if ( (p->DAY != p->DAYL) || (p->GLOB7S_P18 != P[18]) ) p->GLOB7S_CD18 = cos( 2.0*DR*(p->DAY-P[18]) )       ;
    if ( (p->DAY != p->DAYL) || (p->GLOB7S_P32 != P[32]) ) p->GLOB7S_CD32 = cos(     DR*(p->DAY-P[32]) );
    if ( (p->DAY != p->DAYL) || (p->GLOB7S_P39 != P[39]) ) p->GLOB7S_CD39 = cos( 2.0*DR*(p->DAY-P[39]) );
    p->DAYL = p->DAY;
    p->GLOB7S_P32  = P[32];
    p->GLOB7S_P18  = P[18];
    p->GLOB7S_P14  = P[14];
    p->GLOB7S_P39  = P[39];

    // F10.7
    T[1] = P[22]*p->DFA;

    // TIME INDEPENDENT
    T[2] = P[2]*p->PLG[3][1]+P[3]*p->PLG[5][1]+P[23]*p->PLG[7][1] +P[27]*p->PLG[2][1]+P[15]*p->PLG[4][1]+P[60]*p->PLG[6][1];

    // SYMMETRICAL ANNUAL
    T[3] = (P[19]+P[48]*p->PLG[3][1]+P[30]*p->PLG[5][1])*p->GLOB7S_CD32;

    // SYMMETRICAL SEMIANNUAL
    T[4] = (P[16]+P[17]*p->PLG[3][1]+P[31]*p->PLG[5][1])*p->GLOB7S_CD18;

    // ASYMMETRICAL ANNUAL
    T[5] = (P[10]*p->PLG[2][1]+P[11]*p->PLG[4][1]+P[21]*p->PLG[6][1])*p->GLOB7S_CD14;

    // ASYMMETRICAL SEMIANNUAL
    T[6] = (P[38]*p->PLG[2][1])*p->GLOB7S_CD39;

    // DIURNAL
    if ( p->SW[7] != 0 ) {
        T71  = P[12]*p->PLG[3][2]*p->GLOB7S_CD14*p->SWC[5];
        T72  = P[13]*p->PLG[3][2]*p->GLOB7S_CD14*p->SWC[5];
        T[7] = ((P[4]*p->PLG[2][2] + P[5]*p->PLG[4][2] + T71)*p->CTLOC + (P[7]*p->PLG[2][2] + P[8]*p->PLG[4][2] + T72)*p->STLOC);
    }

    // SEMIDIURNAL
    if ( p->SW[8] != 0 ) {
        T81  = (P[24]*p->PLG[4][3]+P[36]*p->PLG[6][3])*p->GLOB7S_CD14*p->SWC[5] ;
        T82  = (P[34]*p->PLG[4][3]+P[37]*p->PLG[6][3])*p->GLOB7S_CD14*p->SWC[5];
        T[8] = ((P[6]*p->PLG[3][3] + P[42]*p->PLG[5][3] + T81)*p->C2TLOC +(P[9]*p->PLG[3][3] + P[43]*p->PLG[5][3] + T82)*p->S2TLOC);
    }


    // TERDIURNAL
    if ( p->SW[14] != 0) {
        T[14] = P[40]*p->PLG[4][4]*p->S3TLOC +P[41]*p->PLG[4][4]*p->C3TLOC;
    }

    // MAGNETIC ACTIVITY
    if ( p->SW[9] != 0 ) {
        if ( p->SW[9] ==  1) T[9] = p->APDF*(P[33]+P[46]*p->PLG[3][1]*p->SWC[2]);
        if ( p->SW[9] == -1) T[9] = (P[51]*p->APT[1]+P[97]*p->PLG[3][1]*p->APT[1]*p->SWC[2]);
    }

    if ( (p->SW[10] != 0) && (p->SW[11] != 0) && (p->LONG > -1000.0) ) {

        // LONGITUDINAL
        T[11] = ( 1.0 + p->PLG[2][1]*(P[81]*p->SWC[5]*cos(DR*(p->DAY-P[82])) + P[86]*p->SWC[6]*cos(2.*DR*(p->DAY-P[87]))) 
                    + P[84]*p->SWC[3]*cos(DR*(p->DAY-P[85])) + P[88]*p->SWC[4]*cos(2.*DR*(p->DAY-P[89])) )
                *(  (P[65]*p->PLG[3][2]+P[66]*p->PLG[5][2]+P[67]*p->PLG[7][2] +P[75]*p->PLG[2][2]+P[76]*p->PLG[4][2]+P[77]*p->PLG[6][2])*cos(DGTR*p->LONG)
                   +(P[91]*p->PLG[3][2]+P[92]*p->PLG[5][2]+P[93]*p->PLG[7][2] +P[78]*p->PLG[2][2]+P[79]*p->PLG[4][2]+P[80]*p->PLG[6][2])*sin(DGTR*p->LONG) );
    }

    for ( TT=0.0, I=1; I<=14; I++ ) TT += fabs( p->SW[I] )*T[I];
    return( TT );

}















/* 
 *      Calculate Temperature and Density Profiles for MSIS models
 *      New lower thermo polynomial 10/30/89
 */
double DENSU( double ALT, double DLB, double TINF, double TLB, double XM, double ALPHA, double *TZ, double ZLB, double S2, int MN1, double *ZN1, double *TN1, double *TGN1, Lgm_Msis00Info *p ) {

//    static double RGAS = 831.4;

//    int     MN, K;
//    double  XS[6], YS[6], Y2OUT[6];
//    double  fDENSU, Z, ZG2, TT, TA, g, g2, DTA, Z1, Z2, T1, T2, ZG, ZGDIF;
//    double  YD1, YD2, X, Y, GLB, GAMMA, EXPL, DENSA, GAMM, YI;

/*
static     int     MN, K;
static     double  XS[6], YS[6], Y2OUT[6];
static     double  fDENSU, Z, ZG2, TT, TA, g, g2, DTA, Z1, Z2, T1, T2, ZG, ZGDIF;
static     double  YD1, YD2, X, Y, GLB, GAMMA, EXPL, DENSA, GAMM, YI;
*/
int     MN, K;
double  XS[6], YS[6], Y2OUT[6];
double  fDENSU, Z, ZG2, TT, TA, g, g2, DTA, Z1, Z2, T1, T2, ZG, ZGDIF;
double  YD1, YD2, X, Y, GLB, GAMMA, EXPL, DENSA, GAMM, YI;


    //ZETA(ZZ, ZL) = (ZZ-ZL)*(p->RE+ZL)/(p->RE+ZZ);
    //////printf("DB %g %g %g %g %g %g %g %g %g %g %g\n", ALT, DLB, TINF, TLB, XM, ALPHA, ZLB, S2, MN1, ZN1, TN1 );
    fDENSU = 1.0;

    // Joining altitude of Bates and spline
    p->ZA = ZN1[1];
    Z  = AMAX1( ALT, p->ZA );

    // Geopotential altitude difference from ZLB
    ZG2 = ZETA( Z, ZLB, p );

    // Bates temperature
    TT    = TINF-(TINF-TLB)*exp(-S2*ZG2);
    TA    = TT;
    *TZ    = TT;
    fDENSU = *TZ;


    if ( ALT < p->ZA ) {
        // CALCULATE TEMPERATURE BELOW ZA
        // Temperature gradient at ZA from Bates profile
        g       = (p->RE+ZLB)/(p->RE+p->ZA); g2 = g*g;
        DTA     = (TINF-TA)*S2*g2;
        TGN1[1] = DTA ;
        TN1[1]  = TA;
        Z  = AMAX1( ALT, ZN1[MN1] );
        MN = MN1;
        Z1 = ZN1[1];
        Z2 = ZN1[MN];
        T1 = TN1[1];
        T2 = TN1[MN];

        // Geopotental difference from Z1
        ZG    = ZETA( Z, Z1, p );
        ZGDIF = ZETA( Z2, Z1, p );

        // Set up spline nodes
        for ( K=1; K<=MN; K++ ) {
            XS[K] = ZETA( ZN1[K], Z1, p )/ZGDIF;
            YS[K] = 1.0/TN1[K];
        }

        // End node derivatives
        YD1 = -TGN1[1]/(T1*T1)*ZGDIF;
        g   = (p->RE+Z2)/(p->RE+Z1); g2 = g*g;
        YD2 = -TGN1[2]/(T2*T2)*ZGDIF*g2;

        // Calculate spline coefficients
        SPLINE( XS, YS, MN, YD1, YD2, Y2OUT );
        X = ZG/ZGDIF;
        SPLINT( XS, YS, Y2OUT, MN, X, &Y );
    
        // temperature at altitude
        *TZ    = 1.0/Y;
        fDENSU = *TZ;

    }


    if ( XM != 0.0 ) {

        // CALCULATE DENSITY ABOVE ZA
        g     = 1.0 + ZLB/p->RE; g2 = g*g;
        GLB   = p->GSURF/g2;
        GAMMA = XM*GLB/(S2*RGAS*TINF);
        EXPL  = exp(-S2*GAMMA*ZG2);
        if ( (EXPL > 50.0) || (TT <= 0.0) ) EXPL = 50.0;

        // Density at altitude
        DENSA = DLB*pow(TLB/TT, 1.0+ALPHA+GAMMA)*EXPL;
        fDENSU = DENSA;

        if ( ALT < p->ZA ) {

            // CALCULATE DENSITY BELOW ZA
            g    = 1.0 + Z1/p->RE; g2 = g*g;
            GLB  = p->GSURF/g2;
            GAMM = XM*GLB*ZGDIF/RGAS;

            // integrate spline temperatures
            SPLINI( XS, YS, Y2OUT, MN, X, &YI );
            EXPL = GAMM*YI;
            if ( (EXPL > 50.0) || (*TZ <= 0.0) ) EXPL = 50.0;

            // Density at altitude
            fDENSU *= pow(T1/(*TZ), 1.+ALPHA)*exp(-EXPL);

        }

    }


    return( fDENSU );

}





/*
 *      Calculate Temperature and Density Profiles for lower atmos.
 */
double DENSM( double ALT, double D0, double XM, double *TZ, int MN3, double *ZN3, double *TN3, double *TGN3, int MN2, double *ZN2, double *TN2, double *TGN2, Lgm_Msis00Info *p ) {

//    static double RGAS = 831.4;

    int     MN, K;
    double  fDENSM, Z, Z1, Z2, T1, T2, ZG, ZGDIF, XS[11], YS[11], YD1, g, g2;
    double  YD2, Y2OUT[11], X, Y, GLB, GAMM, YI, EXPL;

    // changed this fortran statement function to a macro
    //ZETA(ZZ,ZL) = (ZZ-ZL)*(p->RE+ZL)/(p->RE+ZZ);
    fDENSM = D0;

    if ( ALT <= ZN2[1] ) {

        // STRATOSPHERE/MESOSPHERE TEMPERATURE
        Z     = AMAX1(ALT,ZN2[MN2]);
        MN    = MN2;
        Z1    = ZN2[1];
        Z2    = ZN2[MN];
        T1    = TN2[1];
        T2    = TN2[MN];
        ZG    = ZETA( Z, Z1, p );
        ZGDIF = ZETA( Z2, Z1, p );

        // Set up spline nodes
        for ( K=1; K<=MN; K++ ) {
            XS[K] = ZETA( ZN2[K], Z1, p )/ZGDIF;
            YS[K] = 1.0/TN2[K];
        }


        YD1 = -TGN2[1]/(T1*T1)*ZGDIF;
        g = (p->RE+Z2)/(p->RE+Z1); g2 = g*g;
        YD2 = -TGN2[2]/(T2*T2)*ZGDIF*g2;

        // Calculate spline coefficients
        SPLINE( XS, YS, MN, YD1, YD2, Y2OUT );
        X = ZG/ZGDIF;
        SPLINT( XS, YS, Y2OUT, MN, X, &Y );

        // Temperature at altitude
        *TZ =1.0/Y;

        if ( XM != 0.0 ) {
            // CALCULATE STRATOSPHERE/MESOSPHERE DENSITY 
            g    = 1.+Z1/p->RE; g2 = g*g;
            GLB  = p->GSURF/g2;
            GAMM = XM*GLB*ZGDIF/RGAS;

            // Integrate temperature profile
            SPLINI( XS, YS, Y2OUT, MN, X, &YI );
            EXPL = GAMM*YI;
            if ( EXPL > 50.0 ) EXPL=50.0;

            // Density at altitude
            fDENSM *= (T1/(*TZ))*exp(-EXPL);
        }

        if ( ALT <= ZN3[1] ) {
            // TROPOSPHERE/STRATOSPHERE TEMPERATURE
            Z     = ALT;
            MN    = MN3;
            Z1    = ZN3[1];
            Z2    = ZN3[MN];
            T1    = TN3[1];
            T2    = TN3[MN];
            ZG    = ZETA( Z, Z1, p );
            ZGDIF = ZETA( Z2, Z1, p );

            // Set up spline nodes
            for ( K=1; K<= MN; K++ ) {
                XS[K] = ZETA( ZN3[K], Z1, p )/ZGDIF;
                YS[K] = 1.0/TN3[K];
            }

            YD1 = -TGN3[1]/(T1*T1)*ZGDIF;
            g   = (p->RE+Z2)/(p->RE+Z1); g2 = g*g;
            YD2 = -TGN3[2]/(T2*T2)*ZGDIF*g2;

            // Calculate spline coefficients
            SPLINE( XS, YS, MN, YD1, YD2, Y2OUT );
            X = ZG/ZGDIF;
            SPLINT( XS, YS, Y2OUT, MN, X, &Y );

            // temperature at altitude
            *TZ = 1.0/Y;
            if( XM != 0.0 ) {
                // CALCULATE TROPOSPHERIC/STRATOSPHERE DENSITY 
                g    = 1.0 + Z1/p->RE; g2 = g*g;
                GLB  = p->GSURF/g2;
                GAMM = XM*GLB*ZGDIF/RGAS;

                // Integrate temperature profile
                SPLINI( XS, YS, Y2OUT, MN, X, &YI );
                EXPL = GAMM*YI;
                if ( EXPL > 50.0 ) EXPL = 50.0;

                // Density at altitude
                fDENSM *= (T1/(*TZ))*exp(-EXPL);
            }
        }
    }

    if ( XM == 0.0 ) fDENSM = *TZ;

    return( fDENSM );

}




/*
 *       CALCULATE 2ND DERIVATIVES OF CUBIC SPLINE INTERP FUNCTION
 *       ADAPTED FROM NUMERICAL RECIPES BY PRESS ET AL
 *       X,Y: ARRAYS OF TABULATED FUNCTION IN ASCENDING ORDER BY X
 *       N: SIZE OF ARRAYS X,Y
 *       YP1,YPN: SPECIFIED DERIVATIVES AT X(1) AND X(N); VALUES
 *                >= 1E30 SIGNAL SIGNAL SECOND DERIVATIVE ZERO
 *       Y2: OUTPUT ARRAY OF SECOND DERIVATIVES
*/
void SPLINE( double *X, double *Y, int N, double YP1, double YPN, double *Y2 ) {

//    static int  NMAX = 100;
    int         I, K;
    double      SIG, P, QN, UN, U[100];

    if ( YP1 > .99E30 ) {
        Y2[1] = 0.0;
        U[1]  = 0.0;
    } else {
        Y2[1] = -0.5;
        U[1]  = (3./(X[2]-X[1]))*((Y[2]-Y[1])/(X[2]-X[1])-YP1);
    }

    for ( I=2; I<=N-1; I++ ) {
        SIG   = (X[I]-X[I-1])/(X[I+1]-X[I-1]);
        P     = SIG*Y2[I-1]+2.;
        Y2[I] = (SIG-1.)/P;
        U[I]  = (6.*((Y[I+1]-Y[I])/(X[I+1]-X[I])-(Y[I]-Y[I-1])/(X[I]-X[I-1]))/(X[I+1]-X[I-1])-SIG*U[I-1])/P;
    }

    if ( YPN > 0.99e30 ) {
        QN = 0.0;
        UN = 0.0;
    } else {
        QN = 0.5;
        UN = (3./(X[N]-X[N-1]))*(YPN-(Y[N]-Y[N-1])/(X[N]-X[N-1]));
    }

    Y2[N] = (UN-QN*U[N-1])/(QN*Y2[N-1]+1.);
    for ( K=N-1; K>=1; --K ) {
        Y2[K] = Y2[K]*Y2[K+1]+U[K];
    }

    return;

}


/*
 *       CALCULATE CUBIC SPLINE INTERP VALUE
 *       ADAPTED FROM NUMERICAL RECIPES BY PRESS ET AL.
 *       XA,YA: ARRAYS OF TABULATED FUNCTION IN ASCENDING ORDER BY X
 *       Y2A: ARRAY OF SECOND DERIVATIVES
 *       N: SIZE OF ARRAYS XA,YA,Y2A
 *       X: ABSCISSA FOR INTERPOLATION
 *       Y: OUTPUT VALUE
 */
void SPLINT( double *XA, double *YA, double *Y2A, int N, double X, double *Y ) {
    double  H, A, B;
    int     KLO, KHI, K;

    KLO = 1; KHI = N;
    while ( (KHI-KLO) > 1 ) {
        K = (KHI+KLO)/2;
        if ( XA[K] > X) {
            KHI = K;
        } else {
            KLO = K;
        }
    }

    H = XA[KHI]-XA[KLO];
    if ( H == 0) printf("BAD XA INPUT TO SPLINT, H = %g", H );
    A = (XA[KHI]-X)/H;
    B = (X-XA[KLO])/H;
    *Y = A*YA[KLO] + B*YA[KHI] + ((A*A*A-A)*Y2A[KLO] + (B*B*B-B)*Y2A[KHI])*H*H/6.0;

    return;

}


/*
 *      INTEGRATE CUBIC SPLINE FUNCTION FROM XA(1) TO X
 *       XA,YA: ARRAYS OF TABULATED FUNCTION IN ASCENDING ORDER BY X
 *       Y2A: ARRAY OF SECOND DERIVATIVES
 *       N: SIZE OF ARRAYS XA,YA,Y2A
 *       X: ABSCISSA ENDPOINT FOR INTEGRATION
 *       YI: OUTPUT VALUE
 */
void SPLINI( double *XA, double *YA, double *Y2A, int N, double X, double *YI ) {
    //static int     KLO, KHI;
    //static double  XX, H, A, B, A2, B2;
    int     KLO, KHI;
    double  XX, H, A, B, A2, B2;
    *YI  = 0.0; KLO = 1; KHI = 2;
    while ( (X > XA[KLO]) && (KHI <= N) ){
        XX = X;
        if ( KHI < N ) XX = AMIN1( X, XA[KHI] );
        H   = XA[KHI]-XA[KLO];
        A   = (XA[KHI]-XX)/H;
        B   = (XX-XA[KLO])/H;
        A2  = A*A;
        B2  = B*B;
        *YI  += ((1.-A2)*YA[KLO]/2.+B2*YA[KHI]/2.+ ((-(1.+A2*A2)/4.+A2/2.)*Y2A[KLO] + (B2*B2/4.-B2/2.)*Y2A[KHI])*H*H/6.)*H;
        ++KLO;
        ++KHI;
    }
    return;

}


/*
 *      TURBOPAUSE CORRECTION FOR MSIS MODELS
 *        Root mean density
 *      8/20/80
 *         DD - diffusive density
 *         DM - full mixed density
 *         ZHM - transition scale length
 *         XMM - full mixed molecular weight
 *         XM  - species molecular weight
 *         DNET - combined density
 */
double DNET( double DD, double DM, double ZHM, double XMM, double XM ) {

    double  A, YLOG;

    A = ZHM/(XMM-XM);
    if ( (DM <= 0.0) && (DD <= 0.0) ) {
        printf("DNET LOG ERROR %g %g %g\n", DM, DD, XM );
        if ( (DD == 0.0) && (DM == 0.0) ) DD = 1.0;
        if ( DM == 0.0 ) return( DD );
        if ( DD == 0.0 ) return( DM );
    }

    YLOG = A*log( DM/DD );
    if ( YLOG < -10.0 ) return( DD );
    if ( YLOG >  10.0 ) return( DM );

    return( DD*pow( 1.0 + exp(YLOG), 1.0/A ) );
}




/*
 *       CHEMISTRY/DISSOCIATION CORRECTION FOR MSIS MODELS
 *       ALT - altitude
 *       R - target ratio
 *       H1 - transition scale length
 *       ZH - altitude of 1/2 R
 */
double CCOR( double ALT, double R, double H1, double ZH ) {
    double  E;
    E = (ALT-ZH)/H1;
    if ( E >  70.0 ) return( 1.0 );
    if ( E < -70.0 ) return( exp( R ) );
    return( exp( R/( 1.0 + exp(E) ) ) );
}

/*
 *      O&O2 CHEMISTRY/DISSOCIATION CORRECTION FOR MSIS MODELS
 */
double CCOR2( double ALT, double R, double H1, double ZH, double H2 ) {
    double  E1, E2;
    E1 = (ALT-ZH)/H1;
    E2 = (ALT-ZH)/H2;
    if ( (E1 >  70.0) || (E2 >  70.0)) return( 1.0 );
    if ( (E1 < -70.0) && (E2 < -70.0)) return( exp( R ) );
    return( exp( R/( 1.0 + 0.5*( exp(E1) + exp(E2) ) ) ) );
}


