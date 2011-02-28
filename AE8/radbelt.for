C **********************************************************************
C * 8 of June 2006                                                     *
C * Ph.D. Alexey Nikolaevich Petrov                                    *
C * Email: gluk@srd.sinp.msu.ru , alex.petrov@mail.ru                  *
C * Space Physics Research Department,                                 *
C * Scobeltsin Institute of Nuclear Research,                          *
C * Leninskie Gory, 1                                                  *
C * Moscow State Univercity, Moscow, Russia                            *
C * Succesfully compiled by WATCOM FORTRAN/77 Ver.1.4 on WindowsXP SP2 *
C * based on                                                           *
C * RADBELT.FOR   SEPTEMBER 88 by Dieter Bilitza                       *
C *                                                                    *
C***********************************************************************
C*                                                                     *
C*                    TRAPPED RADIATION MODELS                         *
C*                         Program RADBELT                             *
C*                                                                     *
C***  Dieter Bilitza  ******************************  25 March 1988  ***
C***********************************************************************
C***************  Goddard Space Flight Center, code 633  ***************
C**  National Space Science Data Center, Greenbelt, MD 20771, U.S.A.  **
C***********************************************************************
C**********************  NSSDC-ID  PT-14A  *****************************
C***********************************************************************
C***   This program  gives an example of how to use NSSDC's Trapped  ***
C***   Radiation Models. It determines omnidirectional integral and  ***
C***   differential electron or proton fluxes for given energies,    ***
C***   L-values, magnetic field strengths and map-type.              ***
C***   The program will ask you for:                                 ***
C***      NE     number of energies you want displayed               ***
C***      E(1),... E(NE)   energies   or   EBEGIN,EEND,ESTEP  begin, *** 
C***             end and stepsize of energy range                    ***
C***      NL     number of L-values                                  ***
C***      L(1),... L(NL)   L-values   or   LBEGIN,LEND,LSTEP         ***
C***      NB     number of B/B0-values                               ***
C***      B/B0(1),... B/B0(NL)   B/B0-values   or  BBEGIN,BEND,BSTEP ***
C***      MTYPE  map type:  1 AP8MAX   2 AP8MIN  3 AE4MAX  4 AE4MIN  ***
C***                        5 AEI7HI   6 AEI7LO  7 AE8MAX  8 AE8MIN  ***
C***      JTAB   output options: integral or differential fluxes     ***
C***             versus L or B/B0                                    ***
C***   The program interpolates the NSSDC model map in B/B0, L       ***
C***   (function TRARA2), and energy (subroutine TRARA1).            ***
C***   The program opens the map as binary data file (e.g.           ***
C***   AE8MAX.BIN) on unit 15. Make sure that the map you want to    ***
C***   use is available under the correct name in your account       ***
C***   ( see MTYPE for available models and names ).                 ***
C***********************************************************************
C***********************************************************************
C**************************  NOMENCLATURE  *****************************
C********  omnidirectional flux is flux per unit time and      *********
C********      unit sphere surface.                            *********
C********  integral flux for energy E is flux per unit time    *********
C********      and unit surface of particles with energies     *********
C********      greater than E.                                 *********
C********  differential flux for energy E and energy range DE  *********
C********      is average flux per unit time, unit surface and *********
C********      unit energy of particles with energies between  *********
C********      E-DE/2 and E+DE/2.                              *********
C***********************************************************************
C******************************  UNITS  ********************************
C********                 energies: MeV                        *********
C********          integral fluxes: particles/(cm*cm*sec)      ********* 
C********      differential fluxes: particles/(cm*cm*sec*MeV)  *********
C********  magnetic field strength: Gauss                      *********
C***********************************************************************
C***********************************************************************
C*************  DESCRIPTION OF MODEL DATA FILE FORMAT  *****************
C***********************************************************************
C***  THE FILE CONSISTS OF A HEADER ARRAY (IHEAD(8)) AND A MODEL MAP ***
C***  ARRAY (MAP(...)). ALL ELEMENTS ARE INTEGER.                    ***
C***                                                                 ***
C***  IHEAD(1)   MODEL MAP TYPE (SEE ABOVE)                          ***
C***       (2)   INCREMENTS PER DECADE OF LOGARITHMIC FLUX           ***
C***       (3)   EPOCH OF MODEL                                      ***
C***       (4)   SCALE FACTOR FOR ENERGY; E/MEV=E(MAP)/IHEAD(4)      ***
C***                                      =6400 (AE-8),  =100 (AP-8) ***
C***       (5)   SCALE FACTOR FOR L-VALUE =2100 (AE-8), =2048 (AP-8) ***
C***       (6)   SCALE FACTOR FOR B/B0    =1024 (AE-8), =2048 (AP-8) ***
C***       (7)   SCALE FACTOR FOR LOGARITHM OF FLUXES =1024 (AE,AP-8)***
C***       (8)   NUMBER OF ELEMENTS IN MAP =13548 (AE8MAX),          ***
C***              =13168 (AE8MIN), =6509 (AP8MAX), =6688 (AP8MIN)    ***
C***                                                                 ***
C***  LAYOUT OF MAP:                                                 ***
C***      MAP CONSISTS OF SEVERAL VARIABLE-LENGTH SUB-MAPS, EACH     ***
C***      FOR A DIFFERENT ENERGY. EACH SUB-MAP CONSISTS OF SEVERAL   ***
C***      VARIABLE-LENGTH SUB-SUB-MAPS EACH FOR A DIFFERENT L-VALUE. ***
C***      EACH SUB-SUB-MAP CONTAINS THE CURVE LOG(F) [DECADIC        ***
C***      LOGARITHM OF OMNIDIRECTIONAL INTEGRAL PARTICLE FLUX]       ***
C***      VERSUS B/B0 [MAGNETIC FIELD STRENGTH NORMALIZED TO THE     ***
C***      EQUATORIAL VALUE]. THE CURVE IS PARAMETERIZED BY USING     ***
C***      EQUAL INCREMENTS IN LOG(F); THE NUMBER OF INCREMENTS       ***
C***      PER DECADE IS LISTED IN THE HEADER ARRAY [IHEAD(2)]:       ***
C***                                                                 ***
C***         I     B(I)/B(0)   (B(I)-B(I-1))/B(0)   LOG(F(I))        ***
C***       ----------------------------------------------------      ***
C***         0        1                 -              Y             ***
C***         1     B(1)/B(0)   (B(1)-B(0))/B(0)      Y-1/IHEAD(2)    ***
C***         2     B(2)/B(0)   (B(2)-B(1))/B(0)      Y-2/IHEAD(2)    ***
C***         .       ....            .....            ....           ***
C***                                                                 ***
C***      THE SUB-SUB-MAP CONTAINS THE EQUATORIAL FLUX LOGARITHM Y   ***
C***      AND THE B/B0-INCREMENTS (THIRD COLUMN) MULTIPLIED BY       ***
C***      THEIR CORRESPONDING SCALE VALUES ( IHEAD(7) AND (8) ).     ***
C***                                                                 ***
C***  MAP(1)  NUMBER OF ELEMENTS IN SUB-MAP                          ***
C***  MAP(2)  ENERGY FOR THIS SUB-MAP; MAP(2)=E/MEV*IHEAD(4)         ***
C***    MAP(3)  NUMBER OF ELEMENTS IN SUB-SUB-MAP                    ***
C***    MAP(4)  L-VALUE FOR THIS SUB-SUB-MAP; MAP(4)=L*IHEAD(5)      ***
C***    MAP(5)  LOGARITHM OF FLUX AT EQUATOR; MAP(5)=LOG(F0)*IHEAD(7)***
C***      MAP(6)  =(B1-B0)/B0; B1 IS THE MAGNETIC FIELD STRENGTH     ***
C***              THAT CORRESPONDS TO LOG(F1)=LOG(F0)-1/IHEAD(2)     ***
C***      MAP(7)  =(B2-B1)/B0; LOG(F2)=LOG(F1)-1/IHEAD(2)            ***
C***       ...              ....                                     ***
C***      MAP(L)  LAST ELEMENT IN SUB-SUB-MAP; L=MAP(3)+2            ***
C***    MAP(I)  NUMBER OF ELEMENTS IN NEXT SUB-SUB-MAP; I=L+1        ***
C***       ...              ....                                     ***
C***     ...                ....                                     ***
C***    MAP( )  NUMBER OF ELEMENTS IN LAST SUB-SUB-MAP               ***
C***       ...              ....                                     ***
C***      MAP(K)  LAST ELEMENT IN SUB-MAP; K=MAP(1)                  ***
C***  MAP(J)  NUMBER OF ELEMENTS IN NEXT SUB-MAP; J=MAP(1)+1         ***
C***     ...                ....                                     ***
C***       ...              ....                                     ***
C***   ...                  ....                                     ***   
C***  MAP( )  NUMBER OF ELEMENTS IN LAST SUB-MAP                     ***
C***     ...                ....                                     ***
C***       ...              ....                                     ***
C***      MAP(M)  LAST ELEMENT OF MAP; M=IHEAD(8)                    ***
C***********************************************************************
C**                        ENERGY/MEV GRID                           ***
C**  AE-8:    0.04  0.1   0.25  0.5   0.75  1.0  1.5  2.0  2.5  3.0  ***
C**           3.5   4.0   4.5   5.0   5.5   6.0  6.5  7.0 (18 GRID P)***
C**  AE-5,6:  0.04  0.1   0.25  0.5   0.75  1.0  1.5  2.0  2.5  3.0  ***
C**           4.0   4.5   (12 GRID POINTS)                           ***
C**  AE-4:    0.04  0.1   0.3   0.5   1.0   2.0  2.5  3.0  3.5  4.0  ***
C**           4.1   4.25  4.35  4.5   4.65  4.85 (16 GRID POINTS)    ***
C**                                                                  ***
C**                         L-VALUE GRID                             ***
C**           BEGIN  STEP  END   STEP  END   STEP  END   STEP  END   ***
C**  AE-8:     1.2   0.05  1.5    0.1  2.0    0.2  2.4    0.1  3.0   ***
C**            3.0   0.2   3.4    0.1  3.6    0.2  4.4    0.1  4.6   ***
C**            4.6   0.2   5.0    0.5* 8.0    1.0 12.0 (43 GRID P.)  ***
C**                    * 6.6 INSTEAD OF 6.5                          ***
C**  AE-5,6:   1.2   0.05  1.5    0.1  2.0    0.2  2.8 (15 GRID P.)  ***
C**  AE-4:     2.8   0.2   4.0    0.5  6.0    0.6  6.6    0.4  7.0   ***
C**            7.0   1.0  11.0   (16 GRID POINTS)                    ***
C***********************************************************************
      integer           EI
      INTEGER           INFILE,OUTFILE   
      INTEGER           MODEL, FLUXTYPE
      REAL               L_FROMFILE, BB0_FROMFILE, E1, E2   
      DIMENSION         XL(2,10),E(7),FLUX(7),AF(7,10,10),
     &                  IHEAD(8),MAP(20000),DF(5,10,10),EDA(3)
      DIMENSION         OUTPUT_ARR(7)
      CHARACTER*1       ITE(5,10,10)
      CHARACTER*4       BLTEX,LBTEX
      CHARACTER*6       NAME,MNAME(8)
      CHARACTER*10      FNAME
      CHARACTER*11      OUTFILENAME
      CHARACTER*11      INFILENAME      

              
      
      DATA MNAME        /'AP8MAX','AP8MIN','AE4MAX','AE4MIN','AEI7HI',
     &                  'AEI7LO','AE8MAX','AE8MIN'/

      DATA E,XL,AF,DF   /1227*0.0/,EDA/.05,.1,.2/
C
        DO 1567 I=1,5
        DO 1567 K=1,10
        DO 1567 L=1,10
1567    ITE(I,K,L)=' '
C
C               I/O UNIT NUMBERS
C
        EGNR=5                  ! INPUT
        MONITO=6                ! MONITOR
        IUAEAP=15               ! MODEL COEFFICIENTS INPUT
        OUTFILE=16
        INFILE=17
C
C               INPUT OF PARAMETERS
C
      WRITE(MONITO, *) 'TYPE INPUT FILE NAME: '
      READ '(13A)', INFILENAME
      
      OPEN(UNIT=INFILE,FILE=INFILENAME,STATUS='OLD',
     &     FORM='FORMATTED',ERR=8000)     
      
      WRITE(MONITO, *) 'TYPE RESULT FILE NAME: '
      READ '(13A)', OUTFILENAME      
      
      OPEN(UNIT=OUTFILE,FILE=OUTFILENAME,
     &     STATUS='NEW',FORM='FORMATTED') 
     
               
      WRITE(OUTFILE,3910)       
3910  FORMAT('L B/B0 MODEL FLUXTYPE E1 E2 FLUX1')  
      
1000  READ(INFILE, *, END=9999) L_FROMFILE, BB0_FROMFILE, MODEL, 
     &  FLUXTYPE, E1, E2  
     
      IF (FLUXTYPE.EQ.1) THEN
          NE = 1  
          E(1)=E1  
      endif
      if (FLUXTYPE.EQ.2) THEN          
          NE = 2
          E(1)=E1
          E(2)=E2
      ENDIF     

c---------------4. window: model type-------------------------------
        MTYPE = MODEL
        NAME=MNAME(MTYPE)
        WRITE(FNAME,1129) NAME

C ASCII
C Using the ASCII coefficient files instead of the binary 
1129  FORMAT(A6,'.ASC')

      OPEN(IUAEAP,FILE=FNAME,STATUS='OLD',FORM='FORMATTED')
      READ(IUAEAP,1301) IHEAD
      NMAP=IHEAD(8)
      READ(IUAEAP,1301) (MAP(I),I=1,NMAP)
1301  FORMAT(1X,12I6)

      CLOSE(IUAEAP)

C-------------------7. window: number of L-values--------------------
7012  NL=1
C----------------8. window (b): L-values-------------------------
2798  XL(1,1) = L_FROMFILE
C-----------------9. window: number of B/B0 values------------------
7013  NB=1
C------------10. window (a): B/B0 grid--------------------------------
3812  XL(2,1) = BB0_FROMFILE
C----------------THE L-VALUE LOOP-------------------------------
      L=1 
      FL=XL(1,L) 
C----------------THE B LOOP------------------------------------- 
      I=1
      BB0 = XL(2,I)
      CALL TRARA1(IHEAD,MAP,FL,BB0,E,FLUX,NE)
C----------------THE ENERGY LOOP--------------------------------
      DO 1110 K=1,NE      
          AF(K,L,I)=0.0
          IF(FLUX(K).GT.0.0) AF(K,L,I)=10.**FLUX(K) 
          IF(K.EQ.1) GOTO 1110
          DF(K-1,L,I)=ABS(AF(K,L,I)-AF(K-1,L,I))/
     & (E(K)-E(K-1))
          IF(AF(K,L,I).LE.0.0) DF(K-1,L,I)=0.0
1110  CONTINUE      
C-----------------TESTS VALIDITY OF DIFFERENTIAL FLUX------------
       DO 7788 K=1,NE-1
                EI=E(K+1)
                ITEST=0
                EDIFF=EI-E(K)     
C-----------------IS ENERGY INTERVALL LARGE ENOUGH ?-------------
                IEI=1      
                IF(EI.GT.0.10) IEI=2
                IF(EI.GT.0.25) IEI=3
                IF(EDIFF.LT.EDA(IEI)) ITEST=1
                DO 7788 L=1,NL
                        ITT=0
                        IF(XL(1,L).LT.1.2) ITT=1
                        DO 7788 I=1,NB
                                ITE(K,L,I)=' '
                                IF(ITEST+ITT.NE.0) ITE(K,L,I)='?'
7788                            CONTINUE
C------------------TABLE OUTPUT------------------------------------
C        WRITE(MONITO,9101)
9101    FORMAT(1X/)

C      integral flux: FLUXTYPE=1
       if(FLUXTYPE.EQ.1) then
           JTAB = 1
       end if

C      differential flux: FLUXTYPE=2
       if(FLUXTYPE.EQ.2) then
           JTAB = 3
       end if
            
       IBL=2
       BLTEX='B/B0'
       LBTEX=' L  '
       N=NL
       NN=NB
       IBLTAB=1
       IF(IBL.EQ.1) IBLTAB=2
       INDE = 1
       
       OUTPUT_ARR(1) = XL(1,1)  !L
       OUTPUT_ARR(2) = XL(2,1)  !B/B0
       OUTPUT_ARR(3) = MODEL    !MODEL NUMBER       
       OUTPUT_ARR(4) = FLUXTYPE !FLUXTUPE                     
       OUTPUT_ARR(5) = E1 !E1                     
       OUTPUT_ARR(6) = E2 !E2                     
       
       BLV=XL(IBL,INDE)
       DO 1094 I=1,N
        IL=INDE
        IB=INDE          
        IF(IBL.EQ.1) THEN
                IB=I
                XX=XL(2,I)
        ELSE
                IL=I
                XX=XL(1,I)
        ENDIF
        IF(JTAB.LT.3) THEN
            OUTPUT_ARR(7) = AF(1,IL,IB)
        ELSE
            OUTPUT_ARR(7) = DF(1,IL,IB)
        
        ENDIF
        WRITE(OUTFILE, *) OUTPUT_ARR
        WRITE(MONITO, *) OUTPUT_ARR
 
1094    CONTINUE
        goto 1000
8000    write(MONITO, *) 'input file not found' 
9999    STOP     
        END
