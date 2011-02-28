/*
 *  Lgm_AE8_AP8.c - 2010 M. G. Henderson
 *
 *  This is the AE8/AP8 model converted to C. The original fortran is
 *  an extremely awefull example of spaghetti progamming...
 *
 */




// Original FORTRAN COMMENTS RETAINED...
// **********************************************************************
// * 8 of June 2006                                                     *
// * Ph.D. Alexey Nikolaevich Petrov                                    *
// * Email: gluk@srd.sinp.msu.ru , alex.petrov@mail.ru                  *
// * Space Physics Research Department,                                 *
// * Scobeltsin Institute of Nuclear Research,                          *
// * Leninskie Gory, 1                                                  *
// * Moscow State Univercity, Moscow, Russia                            *
// * Succesfully compiled by WATCOM FORTRAN/77 Ver.1.4 on WindowsXP SP2 *
// * based on                                                           *
// * RADBELT.FOR   SEPTEMBER 88 by Dieter Bilitza                       *
// *                                                                    *
//***********************************************************************
//*                                                                     *
//*                    TRAPPED RADIATION MODELS                         *
//*                         Program RADBELT                             *
//*                                                                     *
//***  Dieter Bilitza  ******************************  25 March 1988  ***
//***********************************************************************
//***************  Goddard Space Flight Center, code 633  ***************
//**  National Space Science Data Center, Greenbelt, MD 20771, U.S.A.  **
//***********************************************************************
//**********************  NSSDC-ID  PT-14A  *****************************
//***********************************************************************
//***   This program  gives an example of how to use NSSDC's Trapped  ***
//***   Radiation Models. It determines omnidirectional integral and  ***
//***   differential electron or proton fluxes for given energies,    ***
//***   L-values, magnetic field strengths and map-type.              ***
//***   The program will ask you for:                                 ***
//***      NE     number of energies you want displayed               ***
//***      E(1),... E(NE)   energies   or   EBEGIN,EEND,ESTEP  begin, *** 
//***             end and stepsize of energy range                    ***
//***      NL     number of L-values                                  ***
//***      L(1),... L(NL)   L-values   or   LBEGIN,LEND,LSTEP         ***
//***      NB     number of B/B0-values                               ***
//***      B/B0(1),... B/B0(NL)   B/B0-values   or  BBEGIN,BEND,BSTEP ***
//***      MTYPE  map type:  1 AP8MAX   2 AP8MIN  3 AE4MAX  4 AE4MIN  ***
//***                        5 AEI7HI   6 AEI7LO  7 AE8MAX  8 AE8MIN  ***
//***      JTAB   output options: integral or differential fluxes     ***
//***             versus L or B/B0                                    ***
//***   The program interpolates the NSSDC model map in B/B0, L       ***
//***   (function TRARA2), and energy (subroutine TRARA1).            ***
//***   The program opens the map as binary data file (e.g.           ***
//***   AE8MAX.BIN) on unit 15. Make sure that the map you want to    ***
//***   use is available under the correct name in your account       ***
//***   ( see MTYPE for available models and names ).                 ***
//***********************************************************************
//***********************************************************************
//**************************  NOMENCLATURE  *****************************
//********  omnidirectional flux is flux per unit time and      *********
//********      unit sphere surface.                            *********
//********  integral flux for energy E is flux per unit time    *********
//********      and unit surface of particles with energies     *********
//********      greater than E.                                 *********
//********  differential flux for energy E and energy range DE  *********
//********      is average flux per unit time, unit surface and *********
//********      unit energy of particles with energies between  *********
//********      E-DE/2 and E+DE/2.                              *********
//***********************************************************************
//******************************  UNITS  ********************************
//********                 energies: MeV                        *********
//********          integral fluxes: particles/(cm*cm*sec)      ********* 
//********      differential fluxes: particles/(cm*cm*sec*MeV)  *********
//********  magnetic field strength: Gauss                      *********
//***********************************************************************
//***********************************************************************
//*************  DESCRIPTION OF MODEL DATA FILE FORMAT  *****************
//***********************************************************************
//***  THE FILE CONSISTS OF A HEADER ARRAY (IHEAD(8)) AND A MODEL MAP ***
//***  ARRAY (MAP(...)). ALL ELEMENTS ARE INTEGER.                    ***
//***                                                                 ***
//***  IHEAD(1)   MODEL MAP TYPE (SEE ABOVE)                          ***
//***       (2)   INCREMENTS PER DECADE OF LOGARITHMIC FLUX           ***
//***       (3)   EPOCH OF MODEL                                      ***
//***       (4)   SCALE FACTOR FOR ENERGY; E/MEV=E(MAP)/IHEAD(4)      ***
//***                                      =6400 (AE-8),  =100 (AP-8) ***
//***       (5)   SCALE FACTOR FOR L-VALUE =2100 (AE-8), =2048 (AP-8) ***
//***       (6)   SCALE FACTOR FOR B/B0    =1024 (AE-8), =2048 (AP-8) ***
//***       (7)   SCALE FACTOR FOR LOGARITHM OF FLUXES =1024 (AE,AP-8)***
//***       (8)   NUMBER OF ELEMENTS IN MAP =13548 (AE8MAX),          ***
//***              =13168 (AE8MIN), =6509 (AP8MAX), =6688 (AP8MIN)    ***
//***                                                                 ***
//***  LAYOUT OF MAP:                                                 ***
//***      MAP CONSISTS OF SEVERAL VARIABLE-LENGTH SUB-MAPS, EACH     ***
//***      FOR A DIFFERENT ENERGY. EACH SUB-MAP CONSISTS OF SEVERAL   ***
//***      VARIABLE-LENGTH SUB-SUB-MAPS EACH FOR A DIFFERENT L-VALUE. ***
//***      EACH SUB-SUB-MAP CONTAINS THE CURVE LOG(F) [DECADIC        ***
//***      LOGARITHM OF OMNIDIRECTIONAL INTEGRAL PARTICLE FLUX]       ***
//***      VERSUS B/B0 [MAGNETIC FIELD STRENGTH NORMALIZED TO THE     ***
//***      EQUATORIAL VALUE]. THE CURVE IS PARAMETERIZED BY USING     ***
//***      EQUAL INCREMENTS IN LOG(F); THE NUMBER OF INCREMENTS       ***
//***      PER DECADE IS LISTED IN THE HEADER ARRAY [IHEAD(2)]:       ***
//***                                                                 ***
//***         I     B(I)/B(0)   (B(I)-B(I-1))/B(0)   LOG(F(I))        ***
//***       ----------------------------------------------------      ***
//***         0        1                 -              Y             ***
//***         1     B(1)/B(0)   (B(1)-B(0))/B(0)      Y-1/IHEAD(2)    ***
//***         2     B(2)/B(0)   (B(2)-B(1))/B(0)      Y-2/IHEAD(2)    ***
//***         .       ....            .....            ....           ***
//***                                                                 ***
//***      THE SUB-SUB-MAP CONTAINS THE EQUATORIAL FLUX LOGARITHM Y   ***
//***      AND THE B/B0-INCREMENTS (THIRD COLUMN) MULTIPLIED BY       ***
//***      THEIR CORRESPONDING SCALE VALUES ( IHEAD(7) AND (8) ).     ***
//***                                                                 ***
//***  MAP(1)  NUMBER OF ELEMENTS IN SUB-MAP                          ***
//***  MAP(2)  ENERGY FOR THIS SUB-MAP; MAP(2)=E/MEV*IHEAD(4)         ***
//***    MAP(3)  NUMBER OF ELEMENTS IN SUB-SUB-MAP                    ***
//***    MAP(4)  L-VALUE FOR THIS SUB-SUB-MAP; MAP(4)=L*IHEAD(5)      ***
//***    MAP(5)  LOGARITHM OF FLUX AT EQUATOR; MAP(5)=LOG(F0)*IHEAD(7)***
//***      MAP(6)  =(B1-B0)/B0; B1 IS THE MAGNETIC FIELD STRENGTH     ***
//***              THAT CORRESPONDS TO LOG(F1)=LOG(F0)-1/IHEAD(2)     ***
//***      MAP(7)  =(B2-B1)/B0; LOG(F2)=LOG(F1)-1/IHEAD(2)            ***
//***       ...              ....                                     ***
//***      MAP(L)  LAST ELEMENT IN SUB-SUB-MAP; L=MAP(3)+2            ***
//***    MAP(I)  NUMBER OF ELEMENTS IN NEXT SUB-SUB-MAP; I=L+1        ***
//***       ...              ....                                     ***
//***     ...                ....                                     ***
//***    MAP( )  NUMBER OF ELEMENTS IN LAST SUB-SUB-MAP               ***
//***       ...              ....                                     ***
//***      MAP(K)  LAST ELEMENT IN SUB-MAP; K=MAP(1)                  ***
//***  MAP(J)  NUMBER OF ELEMENTS IN NEXT SUB-MAP; J=MAP(1)+1         ***
//***     ...                ....                                     ***
//***       ...              ....                                     ***
//***   ...                  ....                                     ***   
//***  MAP( )  NUMBER OF ELEMENTS IN LAST SUB-MAP                     ***
//***     ...                ....                                     ***
//***       ...              ....                                     ***
//***      MAP(M)  LAST ELEMENT OF MAP; M=IHEAD(8)                    ***
//***********************************************************************
//**                        ENERGY/MEV GRID                           ***
//**  AE-8:    0.04  0.1   0.25  0.5   0.75  1.0  1.5  2.0  2.5  3.0  ***
//**           3.5   4.0   4.5   5.0   5.5   6.0  6.5  7.0 (18 GRID P)***
//**  AE-5,6:  0.04  0.1   0.25  0.5   0.75  1.0  1.5  2.0  2.5  3.0  ***
//**           4.0   4.5   (12 GRID POINTS)                           ***
//**  AE-4:    0.04  0.1   0.3   0.5   1.0   2.0  2.5  3.0  3.5  4.0  ***
//**           4.1   4.25  4.35  4.5   4.65  4.85 (16 GRID POINTS)    ***
//**                                                                  ***
//**                         L-VALUE GRID                             ***
//**           BEGIN  STEP  END   STEP  END   STEP  END   STEP  END   ***
//**  AE-8:     1.2   0.05  1.5    0.1  2.0    0.2  2.4    0.1  3.0   ***
//**            3.0   0.2   3.4    0.1  3.6    0.2  4.4    0.1  4.6   ***
//**            4.6   0.2   5.0    0.5* 8.0    1.0 12.0 (43 GRID P.)  ***
//**                    * 6.6 INSTEAD OF 6.5                          ***
//**  AE-5,6:   1.2   0.05  1.5    0.1  2.0    0.2  2.8 (15 GRID P.)  ***
//**  AE-4:     2.8   0.2   4.0    0.5  6.0    0.6  6.6    0.4  7.0   ***
//**            7.0   1.0  11.0   (16 GRID POINTS)                    ***
//***********************************************************************
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Lgm_AE8_AP8.h"

#define TRUE 1
#define FALSE 0
#define AMIN1(a,b)  ((a)<(b))?(a):(b)
#define AMAX1(a,b)  ((a)>(b))?(a):(b)







double  AE8_AP8_Flux( double L, double BB0, int MODEL, int FLUXTYPE, double E1, double E2 ) {

    int     i, j, k, l;
    int     *MAP, *IHEAD;
    int     NE, NL, NB;
    int     ITEST, IEI, ITT, JTAB, IBL, N, NN, IBLTAB, INDE, IL, IB;
    int     EI, MODEL, FLUXTYPE;
    double  XX, EDIFF, BLV, BB0, FL;
    double  L_FROMFILE, BB0_FROMFILE, E1, E2;
    double  XL[3][11], E[8], FLUX[8], AF[8][11][11], DF[6][11][11], EDA[4], OUTPUT_ARR[8];
    char    ITE[6][11][11];
    char    BLTEX[5], LBTEX[5];
    char    NAME[7], MNAME[9];
    char    fName[11], InFilename[12], OutFilename[12];
    FILE    *fp_in, *fp_out;

      
    EDA[1] = 0.05;
    EDA[2] = 0.10;
    EDA[3] = 0.20;

    for (i=1; i<=2; i++) {
        for (k=1; k<=10; k++) XL[i][k] = 0.0;
    }

    for (i=1; i<=7; i++) E[i] = 0.0;
    for (j=1; j<=10; j++){
        for (k=1; k<=10; k++){
            for (i=1; i<=5; i++){
                DF[i][j][k]  = 0.0;
                ITE[i][j][k] = ' ';
            }
            for (i=1; i<=7; i++){
                AF[i][j][k] = 0.0;
            }
        }
    }

    if ( FLUXTYPE == 1 ) {
        NE   = 1;
        E[1] = E1;
    } else if ( FLUXTYPE == 2 ) {
        NE   = 2;
        E[1] = E1;
        E[2] = E2;
    }

    switch ( MODEL ) {
        case 1:
            MAP   = AP8MAX_MAP;
            IHEAD = AP8MAX_IHEAD;
            break;
        case 2:
            MAP   = AP8MIN_MAP;
            IHEAD = AP8MIN_IHEAD;
            break;
        case 7:
            MAP   = AE8MAX_MAP;
            IHEAD = AE8MAX_IHEAD;
            break;
        case 8:
            MAP   = AE8MIN_MAP;
            IHEAD = AE8MIN_IHEAD;
            break;
        default:
            printf( "Lgm_AE8: Unknown model. MODEL = %d\n", MODEL);
    }
    


    /*
     * L values
     */
    NL = 1;                 // Number of L values
    XL[1][1] = L_FROMFILE;


    /*
     * B/B0 values
     */
    NB = 1;                 // Number of B/B0 values
    XL[2][1] = BB0_FROMFILE;




    // THE L-VALUE LOOP (MGH LOOP?? l is only ever 1)
    l  = 1;
    FL = XL[1][l];

    // THE B LOOP (MGH LOOP?? i is only ever 1)
    i = 1;
    BB0 = XL[2][i];
    TRARA1( IHEAD, MAP, FL, BB0, E, FLUX, NE );


    // THE ENERGY LOOP
    for ( k=1; k<= NE; k++ ){
        AF[k][l][i] = 0.0;
        if ( FLUX[k] > 0.0 ) AF[k][l][i] = pow( 10.0, FLUX[k] );
        if ( k > 1 ) {
            DF[k-1][l][i] = fabs( AF[k][l][i]-AF[k-1][l][i] )/(E[k]-E[k-1]);
            if ( AF[k][l][i] <= 0.0 ) DF[k-1][l][i] = 0.0;
        }
    }


    /* 
     * TESTS VALIDITY OF DIFFERENTIAL FLUX
     */
    for ( k=1; k<NE; k++ ) {

        EI    = E[k+1];
        ITEST = 0;
        EDIFF = EI-E[k];

        // IS ENERGY INTERVALL LARGE ENOUGH?
        IEI = 1;
        if ( EI > 0.10 ) IEI = 2;
        if ( EI > 0.25 ) IEI = 3;
        if ( EDIFF < EDA[IEI] ) ITEST = 1;

        for ( l=1; l<=NL; l++ ) {
            ITT = ( XL[1][l] < 1.2 ) ? 1 : 0;
            for ( i=1; i<=NB; i++ ) {
                ITE[k][l][i] = ( (ITEST+ITT) != 0 ) ? '?' : ' ';
            }
        }

    }


        
    if (FLUXTYPE == 1) {
        JTAB = 1;    // integral flux: FLUXTYPE=1
    } else if ( FLUXTYPE == 2 ) {
        JTAB = 3;    // differential flux: FLUXTYPE=2
    }
            
    IBL    = 2;
    N      = NL;
    NN     = NB;
    IBLTAB = 1;
    if (IBL == 1) IBLTAB=2;
    INDE = 1;
   
    BLV = XL[IBL][INDE];
    i = 1;
    IL = INDE;
    IB = INDE;
    if ( IBL == 1 ) {
        IB = i;
        XX = XL[2][i];
    } else {
        IL = i;
        XX = XL[1][i];
    }
    if (JTAB < 3 ) {
            Flux = AF[1][IL][IB];
    } else {
            Flux = DF[1][IL][IB];
    }

    return( Flux );


}



/* M. G. Henderson - Converted to C, 2010. This is some of the worst
 * spaghetti-coding Ive ever seen!  Just aweful..... It is a complete and
 * totally incomprehensible mess.
 *
 * TRMFUN.FOR	1987
 *
 *********************************************************************
 **************** SUBROUTINES, FUNCTIONS *****************************
 *********************************************************************
 ******************** TRARA1, TRARA2 *********************************
 *********************************************************************
 *
 *    SUBROUTINE TRARA1(DESCR,MAP,FL,BB0,E,F,N)                         
 ************************************************************************
 **** TRARA1 FINDS PARTICLE FLUXES FOR GIVEN ENERGIES, MAGNETIC FIELD *** 
 **** STRENGTH AND L-VALUE. FUNCTION TRARA2 IS USED TO INTERPOLATE IN ***
 **** B-L-SPACE.                                                      ***
 ****   INPUT: DESCR(8)   HEADER OF SPECIFIED TRAPPED RADITION MODEL  ***
 ****          MAP(...)   MAP OF TRAPPED RADITION MODEL               ***
 ****                     (DESCR AND MAP ARE EXPLAINED AT THE BEGIN   ***
 ****                     OF THE MAIN PROGRAM MODEL)                  ***
 ****          N          NUMBER OF ENERGIES                          ***
 ****          E(N)       ARRAY OF ENERGIES IN MEV                    ***
 ****          FL         L-VALUE                                     ***
 ****          BB0        =B/B0  MAGNETIC FIELD STRENGTH NORMALIZED   ***
 ****                     TO FIELD STRENGTH AT MAGNETIC EQUATOR       ***
 ****  OUTPUT: F(N)       DECADIC LOGARITHM OF INTEGRAL FLUXES IN     ***
 ****                     PARTICLES/(CM*CM*SEC)                       ***
 ************************************************************************
 */


void TRARA1( int DESCR[], int MAP[], double FL, double BB0, double E[], double F[], int N ) {

    int     S0, S1, S2;
    int     NL, NB, I0, I1, I2, I3, L3, IE;
    double  FISTEP, ESCALE, FSCALE, XNL, E0, E1, E2, F0;
    double  F1 = 1.001;
    double  F2 = 1.002;

    FISTEP = (double)DESCR[7]/(double)DESCR[2];
    ESCALE = (double)DESCR[4];
    FSCALE = (double)DESCR[7];
    XNL    = AMIN1( 15.6, fabs(FL) );
    NL     = (int)(XNL*DESCR[5]);

    if ( BB0 < 1.0 ) BB0 = 1.0;
    NB = (int)((BB0-1.0)*DESCR[6]);

    /*
     *  I2 IS THE NUMBER OF ELEMENTS IN THE FLUX MAP FOR THE FIRST ENERGY.  
     *  I3 IS THE INDEX OF THE LAST ELEMENT OF THE SECOND ENERGY MAP.       
     *  L3 IS THE LENGTH OF THE MAP FOR THE THIRD ENERGY.                   
     *  E1 IS THE ENERGY OF THE FIRST ENERGY MAP (UNSCALED)                 
     *  E2 IS THE ENERGY OF THE SECOND ENERGY MAP (UNSCALED)                
     */ 
    I1 = 0;
    I2 = MAP[1];
    I3 = I2+MAP[I2+1];
    L3 = MAP[I3+1];
    E1 = MAP[I1+2]/ESCALE;
    E2 = MAP[I2+2]/ESCALE;

    /* 
     *  ENERGY LOOP
     */
    for ( IE=1; IE<=N; ++IE ) {

        /*
         *
         * FOR EACH ENERGY E(I) FIND THE SUCCESSIVE ENERGIES E0,E1,E2 IN MODEL
         * MAP, WHICH OBEY  E0 < E1 < E(I) < E2 . 
         */
        while ( (E[IE] > E2) && (L3 != 0) ) {
            I0 = I1; I1 = I2; I2 = I3;
            I3 += L3;
            L3 = MAP[I3+1];
            E0 = E1; E1 = E2;
            E2 = MAP[I2+2]/ESCALE;
            F0 = F1; F1 = F2;
        }

        /*
         * CALL TRARA2 TO INTERPOLATE THE FLUX-MAPS FOR E1,E2 IN L-B/B0- SPACE
         * TO FIND FLUXES F1,F2 [IF THEY HAVE NOT ALREADY BEEN CALCULATED FOR A
         * PREVIOUS E(I)].
         *
         */
        F1 = TRARA2( MAP + I1+3-1, NL, NB, FISTEP )/FSCALE;
        F2 = TRARA2( MAP + I2+3-1, NL, NB, FISTEP )/FSCALE;

        // FINALLY, INTERPOLATE IN ENERGY.
        F[IE] = F1 + (F2-F1)*(E[IE]-E1)/(E2-E1);
        if ( (F2 <= 0.0) && (I1 != 0) ) {

            /*                                                                      
             *  ---------------- SPECIAL INTERPOLATION -------------------------
             *  IF THE FLUX FOR THE SECOND ENERGY CANNOT BE FOUND (I.E. F2=0.0),
             *  AND THE ZEROTH ENERGY MAP HAS BEEN DEFINED (I.E. I1 NOT EQUAL 0),
             *  THEN INTERPOLATE USING THE FLUX MAPS FOR THE ZEROTH AND FIRST
             *  ENERGY AND CHOOSE THE MINIMUM OF THIS INTERPOLATIONS AND THE
             *  INTERPOLATION THAT WAS DONE WITH F2=0. 
             *                                                                         
             */
            F0 = TRARA2( MAP + I0+3-1, NL, NB, FISTEP )/FSCALE; 
            F[IE] = AMIN1( F[IE], F0 + (F1-F0)*(E[IE]-E0)/(E1-E0) );

        }



        /*
         * THE LOGARITHMIC FLUX IS ALWAYS KEPT GREATER OR EQUAL ZERO.
         */
        F[IE] = AMAX1( F[IE], 0.0 );
    }

    return;
}




// FUNCTION TRARA2(MAP,IL,IB)                                        
//*****************************************************************
//***  TRARA2 INTERPOLATES LINEARLY IN L-B/B0-MAP TO OBTAIN     ***
//***  THE LOGARITHM OF INTEGRAL FLUX AT GIVEN L AND B/B0.      ***
//***    INPUT: MAP(..) IS SUB-MAP (FOR SPECIFIC ENERGY) OF     ***
//***                   TRAPPED RADIATION MODEL MAP             ***
//***           IL      SCALED L-VALUE                          ***
//***           IB      SCALED B/B0-1                           ***
//***   OUTPUT: TRARA2  SCALED LOGARITHM OF PARTICLE FLUX       ***
//*****************************************************************
//***  SEE MAIN PROGRAM 'MODEL' FOR EXPLANATION OF MAP FORMAT   ***
//***  SCALING FACTORS.                                         ***
//***  THE STEPSIZE FOR THE PARAMETERIZATION OF THE LOGARITHM   ***
//***  OF FLUX IS OBTAINED FROM 'COMMON/TRA2/'.                 ***
//*****************************************************************
// Whoever coded up the original FORTRAN should be shot!
double  TRARA2( int MAP[], int IL, int IB, double FISTEP ) {

    
    int     done, flag, flag2, flag3, flag4, flag5;
    int     ITIME, I1, I2, J1, J2, KT, L1, L2;
    double  FINCR1, FKBJ1, FLOGM, FLOG, FKBM, FKB, FKBJ2, Result;
    double  FNL, FNB, FLL1, FLL2, DFL, FLOG1, FLOG2, FKB1, FKB2, FINCR2, SL2, SL1;

    FNL   = (double)IL;
    FNB   = (double)IB;
    ITIME = 0;



    /*
     *  FIND CONSECUTIVE SUB-SUB-MAPS FOR SCALED L-VALUES LS1,LS2, WITH IL LESS
     *  OR EQUAL LS2.  L1,L2 ARE LENGTHS OF SUB-SUB-MAPS.  I1,I2 ARE INDICES OF
     *  FIRST ELEMENTS MINUS 1.
     */
    I2    = 0;
    done = FALSE;
    while ( !done ) {
        L2 = MAP[I2+1];
        if ( MAP[I2+2] > IL) {
            done = TRUE;
        } else {
	        I1 = I2;
    	    L1 = L2;
      	    I2 += L2;
        }
    }





    /*  
     * IF SUB-SUB-MAPS ARE EMPTY, I. E. LENGTH LESS 4, THAN TRARA2=0
     */       
    if ( (L1<4) && (L2<4)) {
        return( 0.0 );                            
     }

    flag  = 0;
    flag2 = 0;
    flag3 = 0;
    flag4 = 0;
    flag5 = 0;
    done  = 0;
    while (!done) {

        if ( (flag==1) || (MAP[I2+3]<=MAP[I1+3]) ) {
            KT = I1; I1 = I2; I2 = KT; KT = L1; L1 = L2; L2 = KT;
        }

        //
        // DETERMINE INTERPOLATE IN SCALED L-VALUE
        // 
        FLL1  = MAP[I1+2]; FLL2 = MAP[I2+2];
        DFL   = (FNL-FLL1)/(FLL2-FLL1);
        FLOG1 = MAP[I1+3]; FLOG2 = MAP[I2+3];
        FKB1  = 0.0; FKB2  = 0.0;

        if (L1>=4) {

            /*
             * B/B0 LOOP
             */                                               
            for (J2=4; J2<=L2; J2++){
                FINCR2 = MAP[I2+J2];
                if ( (FKB2+FINCR2) > FNB) { flag5 = 1; break; }
                FKB2  = FKB2+FINCR2; FLOG2 = FLOG2-FISTEP;
            }

            if ( flag5 == 0 ){
                ++ITIME;
                if (ITIME == 1) { flag = 1; } 
                else { 
                    return(0.0); 
                }
            }

            if ( ITIME != 1 ) {
                if ( J2 == 4 ) {
                    done  = 1;
                    flag3 = 1;
                } else {
                    SL2 = FLOG2/FKB2;
                    for ( J1=4; J1<=L1; J1++ ) {
                        FINCR1 = MAP[I1+J1];
                        FKB1   = FKB1+FINCR1;
                        FLOG1  = FLOG1-FISTEP;
                        FKBJ1  = ((FLOG1/FISTEP)*FINCR1+FKB1)/((FINCR1/FISTEP)*SL2+1.0);
                        if ( FKBJ1 <= FKB1 ) break;
                    }

                    if ( (FKBJ1 > FKB1) && (FKBJ1 <= FKB2) ) {
                        return(0.0);                                   
                    }

                    if ( FKBJ1 <= FKB2 ) {
                        FKBM  = FKBJ1+(FKB2-FKBJ1)*DFL;
                        FLOGM = FKBM*SL2;
                        FLOG2 = FLOG2-FISTEP;
                        FKB2  = FKB2+FINCR2;
                        SL1   = FLOG1/FKB1;
                        SL2   = FLOG2/FKB2;

                        done  = 1;
                        flag4 = 1;
                    } else {
                        FKB1 = 0.0;
                    }
                }
            }
            if ( (flag == 0) && (flag3 == 0) && (flag4 == 0) ) FKB2 = 0.0;
        } else {
            done = 1;
        }
    }

    if (flag4 == 0 ){
        if ( flag3 == 0 ){
            J2 = 4;
            FINCR2 = MAP[I2+J2];
            FLOG2  = MAP[I2+3];
            FLOG1  = MAP[I1+3];
        }



        FLOGM = FLOG1+(FLOG2-FLOG1)*DFL;
        FKBM  = 0.0;
        FKB2  = FKB2+FINCR2;
        FLOG2 = FLOG2-FISTEP;
        SL2   = FLOG2/FKB2;
        if ( L1 < 4 ) {
            FINCR1 = 0.0;
            SL1    = -900000.0;
            FKBJ1  = ((FLOG1/FISTEP)*FINCR1+FKB1)/((FINCR1/FISTEP)*SL2+1.0);
            FKB    = FKBJ1+(FKB2-FKBJ1)*DFL;
            FLOG   = FKB*SL2;
            if ( FKB >= FNB ) {
                if ( FKB < (FKBM+1e-10) ) {
                    return(0.0);                               
                } else {
                    Result = FLOGM+(FLOG-FLOGM)*((FNB-FKBM)/(FKB-FKBM));
                    Result = AMAX1(Result,0.0);
                    return( Result );
                }
            }
            FKBM  = FKB;
            FLOGM = FLOG;
            if ( J2 >= L2 ) {
                return(0.0);                                        
            }
            ++J2;
            FINCR2 = MAP[I2+J2];
            FLOG2  = FLOG2-FISTEP;
            FKB2   = FKB2+FINCR2;
            SL2    = FLOG2/FKB2;
        } else {
            J1 = 4;
            FINCR1 = MAP[I1+J1];
            FKB1   = FKB1+FINCR1;
            FLOG1  = FLOG1-FISTEP;
            SL1    = FLOG1/FKB1;
        }
    }


    while (1){
        if ( SL1 < SL2 ) {
            FKBJ1 = ((FLOG1/FISTEP)*FINCR1+FKB1)/((FINCR1/FISTEP)*SL2+1.0);
            FKB   = FKBJ1+(FKB2-FKBJ1)*DFL;
            FLOG  = FKB*SL2;
            if ( FKB >= FNB ) {
                if ( FKB < (FKBM+1e-10) ) {
                    return(0.0);                               
                } else {
                    Result = FLOGM+(FLOG-FLOGM)*((FNB-FKBM)/(FKB-FKBM));
                    Result = AMAX1(Result,0.0);
                    return( Result );
                }
            }
            FKBM  = FKB;
            FLOGM = FLOG;
            if ( J2 >= L2 ) {
                return(0.0);                                        
            }
            ++J2;
            FINCR2 = MAP[I2+J2];
            FLOG2  = FLOG2-FISTEP;
            FKB2   = FKB2+FINCR2;
            SL2    = FLOG2/FKB2;
        } else {
            FKBJ2 = ((FLOG2/FISTEP)*FINCR2+FKB2)/((FINCR2/FISTEP)*SL1+1.0);
            FKB   = FKB1+(FKBJ2-FKB1)*DFL;
            FLOG  = FKB*SL1;
            if ( FKB >= FNB) {
                if ( FKB < (FKBM+1e-10) ) {
                    return(0.0);                               
                } else {
                    Result = FLOGM+(FLOG-FLOGM)*((FNB-FKBM)/(FKB-FKBM));
                    Result = AMAX1(Result,0.0);
                    return( Result );
                }
            }
            FKBM  = FKB;
            FLOGM = FLOG;
            if ( J1 >= L1 ) {
                return(0.0);                                       
            }
            ++J1;
            FINCR1 = MAP[I1+J1];
            FLOG1  = FLOG1-FISTEP;
            FKB1   = FKB1+FINCR1;
            SL1    = FLOG1/FKB1;
        }
    }


    if ( FKB < (FKBM+1e-10) ) {
        return(0.0);                               
    } else {
        Result = FLOGM+(FLOG-FLOGM)*((FNB-FKBM)/(FKB-FKBM));
        Result = AMAX1(Result,0.0);
        return( Result );                                                         
    }


}
