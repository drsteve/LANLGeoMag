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
#include <stdio.h>
#include <math.h>
#define TRUE 1
#define FALSE 0
#define AMIN1(a,b)  ((a)<(b))?(a):(b)
#define AMAX1(a,b)  ((a)>(b))?(a):(b)


double  TRARA2( int MAP[], int IL, int IB, double FISTEP );

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



// STATIC VARIABLES NOT THREAD SAFE!
    /*
     *  S0, S1, S2 ARE LOGICAL VARIABLES WHICH INDICATE WHETHER THE FLUX FOR A
     *  PARTICULAR E, B, L POINT HAS ALREADY BEEN FOUND IN A PREVIOUS CALL  TO
     *  FUNCTION TRARA2. IF NOT, S.. =.TRUE.
     */
//    S1 = TRUE;
//    S2 = TRUE;


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
//            S0 = S1; S1 = S2;
//            S2 = TRUE;
            F0 = F1; F1 = F2;
        }

        /*
         * CALL TRARA2 TO INTERPOLATE THE FLUX-MAPS FOR E1,E2 IN L-B/B0- SPACE
         * TO FIND FLUXES F1,F2 [IF THEY HAVE NOT ALREADY BEEN CALCULATED FOR A
         * PREVIOUS E(I)].
         *
         */
//        if ( S1 ) 
        F1 = TRARA2( MAP + I1+3-1, NL, NB, FISTEP )/FSCALE;
//printf("F1 = %g (%d %d %d %g)\n", F1, *(MAP + I1+3), NL, NB, FISTEP);

//        if ( S2 ) 
        F2 = TRARA2( MAP + I2+3-1, NL, NB, FISTEP )/FSCALE;
//printf("F2 = %g (%d %d %d %g)\n", F2, *(MAP + I1+3), NL, NB, FISTEP);
//        S1 = FALSE;
//        S2 = FALSE;

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
            //        if ( S0 ) 
            F0 = TRARA2( MAP + I0+3-1, NL, NB, FISTEP )/FSCALE; 
            //        S0 = FALSE; 
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
//printf("$$$$$$$$$ I2,J2,FKB2, FINCR2, FNB, FLOG2, FISTEP = %d %d %g %g %g %g %g\n", I2, J2, FKB2, FINCR2, FNB, FLOG2, FISTEP);
            }

            if ( flag5 == 0 ){
//printf("DFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFS ITIME = %d\n", ITIME);
                ++ITIME;
                if (ITIME == 1) { flag = 1; } 
                else { 
                    return(0.0); 
                }
            }

            if ( ITIME != 1 ) {

//printf("HEREEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEe\n");
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
//printf("00000000000000000000 SL1, SL2 = %g %g\n", SL1, SL2);

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
//printf("11111111111111111111 SL2 = %g\n", SL2);
        if ( L1 < 4 ) {
            FINCR1 = 0.0;
            SL1    = -900000.0;
            FKBJ1  = ((FLOG1/FISTEP)*FINCR1+FKB1)/((FINCR1/FISTEP)*SL2+1.0);
            FKB    = FKBJ1+(FKB2-FKBJ1)*DFL;
            FLOG   = FKB*SL2;
//printf("1. ####################################   FKB, FKB1, FKBJ2, FKB1, DFL, FNB = %g %g %g %g %g %g\n", FKB, FKB1, FKBJ2, FKB1, DFL, FNB);
            if ( FKB >= FNB ) {
                if ( FKB < (FKBM+1e-10) ) {
                    return(0.0);                               
                } else {
                    Result = FLOGM+(FLOG-FLOGM)*((FNB-FKBM)/(FKB-FKBM));
                    Result = AMAX1(Result,0.0);
//printf("1. Result,FLOGM,FLOG,FNB,FKBM=%g %g %g %g %g\n", Result,FLOGM,FLOG,FNB,FKBM);
                    return( Result );
                }
            }
//printf("DFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFS\n");
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

//printf("FKBM = %g\n", FKBM);
//printf("99999999999999999999 SL2 = %g\n", SL2);

    while (1){
        if ( SL1 < SL2 ) {
            FKBJ1 = ((FLOG1/FISTEP)*FINCR1+FKB1)/((FINCR1/FISTEP)*SL2+1.0);
            FKB   = FKBJ1+(FKB2-FKBJ1)*DFL;
            FLOG  = FKB*SL2;
//printf("2. ####################################   FKB, FKB1, FKBJ2, FKB1, DFL, FNB, SL1, SL2 = %g %g %g %g %g %g %g %g\n", FKB, FKB1, FKBJ2, FKB1, DFL, FNB, SL1, SL2);
            if ( FKB >= FNB ) {
                if ( FKB < (FKBM+1e-10) ) {
                    return(0.0);                               
                } else {
                    Result = FLOGM+(FLOG-FLOGM)*((FNB-FKBM)/(FKB-FKBM));
                    Result = AMAX1(Result,0.0);
//printf("2. Result,FLOGM,FLOG,FNB,FKBM=%g %g %g %g %g\n", Result,FLOGM,FLOG,FNB,FKBM);
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
//printf("3. ####################################   FKB, FKB1, FKBJ2, FKB1, DFL, FNB, SL1, SL2  = %g %g %g %g %g %g %g %g\n", FKB, FKB1, FKBJ2, FKB1, DFL, FNB, SL1, SL2 );
            if ( FKB >= FNB) {
                if ( FKB < (FKBM+1e-10) ) {
                    return(0.0);                               
                } else {
                    Result = FLOGM+(FLOG-FLOGM)*((FNB-FKBM)/(FKB-FKBM));
                    Result = AMAX1(Result,0.0);
//printf("3. Result,FLOGM,FLOG,FNB,FKBM=%g %g %g %g %g\n", Result,FLOGM,FLOG,FNB,FKBM);
                    return( Result );
                }
            }
//printf("DFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFS\n");
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
