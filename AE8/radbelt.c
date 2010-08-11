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
#include "Lgm_AE8.h"
static char *ModelName[] = { " ", "AP8MAX", "AP8MIN", "AE4MAX", "AE4MIN", "AEI7HI", "AEI7LO", "AE8MAX", "AE8MIN" };
static char *FluxType[] = { " ", "Integral", "Differential" };

void TRARA1( int DESCR[], int MAP[], double FL, double BB0, double E[], double F[], int N );


int main (){

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


    /*
     * Open Input File for Reading
     */
    strcpy( InFilename, "input" );
    if ( (fp_in = fopen( InFilename, "r" )) == NULL ){
        printf("Couldn't open file for reading: %s\n", InFilename);
        exit(1);
    }

    /*
     * Open Output File for Writing
     */
    strcpy( OutFilename, "mike.dat" );
    if ( (fp_out = fopen( OutFilename, "w" )) == NULL ){
        printf("Couldn't open file for writing: %s\n", OutFilename);
        exit(1);
    }


    /*
     *  Read Input File
     */
    printf( "     L         B/B0     MODEL    FLUXTYPE      E1         E2        FLUX1\n" );
    while ( fscanf( fp_in, "%lf %lf %d %d %lf %lf", &L_FROMFILE, &BB0_FROMFILE, &MODEL, &FLUXTYPE, &E1, &E2 ) != EOF ) {


        if ( FLUXTYPE == 1 ) {
            NE   = 1;
            E[1] = E1;
        } else if ( FLUXTYPE == 2 ) {
            NE   = 2;
            E[1] = E1;
            E[2] = E2;
        }


        /*
         * Set Model Type, and read the datas file
         */
//        MTYPE = MODEL;
//        strcpy( NAME, MNAME[MTYPE] );
//
////---------------4. window: model type-------------------------------
//        MTYPE = MODEL
//        NAME=MNAME(MTYPE)
//        WRITE(FNAME,1129) NAME
//
//// ASCII
//// Using the ASCII coefficient files instead of the binary 
//1129  FORMAT(A6,'.ASC')
//
//      OPEN(IUAEAP,FILE=FNAME,STATUS='OLD',FORM='FORMATTED')
//      READ(IUAEAP,1301) IHEAD
//      WRITE(*, *) IHEAD
//      NMAP=IHEAD(8)
//      READ(IUAEAP,1301) (MAP(I),I=1,NMAP)
//1301  FORMAT(1X,12I6)
//      CLOSE(IUAEAP)


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


        /*
         * TABLE OUTPUT
         */
        
        if (FLUXTYPE == 1) {
            JTAB = 1;    // integral flux: FLUXTYPE=1
        } else if ( FLUXTYPE == 2 ) {
            JTAB = 3;    // differential flux: FLUXTYPE=2
        }
                
        IBL    = 2;
        sprintf( BLTEX, "B/B0" );
        sprintf( LBTEX, " L  " );
        N      = NL;
        NN     = NB;
        IBLTAB = 1;
        if (IBL == 1) IBLTAB=2;
        INDE = 1;
       
        OUTPUT_ARR[1] = XL[1][1];   // L
        OUTPUT_ARR[2] = XL[2][1];   // B/B0
        OUTPUT_ARR[3] = MODEL;      // MODEL NUMBER       
        OUTPUT_ARR[4] = FLUXTYPE;   // FLUXTYPE                     
        OUTPUT_ARR[5] = E1;         // E1                     
        OUTPUT_ARR[6] = E2;         // E2                     
           
        BLV = XL[IBL][INDE];
        for ( i=1; i<=N; i++){
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
                OUTPUT_ARR[7] = AF[1][IL][IB];
            } else {
                OUTPUT_ARR[7] = DF[1][IL][IB];
            }

            printf("%-10.6lf %-10.6lf %-12s %-12s %-10.6lf %-10.6lf %-10.6lf\n", OUTPUT_ARR[1], OUTPUT_ARR[2], ModelName[(int)OUTPUT_ARR[3]], FluxType[(int)OUTPUT_ARR[4]], OUTPUT_ARR[5], OUTPUT_ARR[6], OUTPUT_ARR[7]);
            //WRITE(OUTFILE, *) OUTPUT_ARR
            //WRITE(MONITO, *) OUTPUT_ARR
        }


    }

    return(0);

}
