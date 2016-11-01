#ifdef HAVE_CONFIG_H
// MGH
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Lgm/Lgm_Tsyg2007.h"
#include "Lgm/Lgm_DynamicMemory.h"
#include "Lgm/Lgm_CTrans.h"
#include <gsl/gsl_sf_bessel.h>
#include <float.h>

#define USE_CACHEING 0


void mysincos(double val, double *sin_val, double *cos_val);

void Lgm_SetTabulatedBessel_TS07( int Flag, LgmTsyg2007_Info *t ){

    t->UseTabulatedBessel = ( Flag <= 0 ) ? 0 : 1;

    if ( !(t->BesselTableAlloced) && t->UseTabulatedBessel ) {

        t->BesselTable = (Lgm_TabularBessel *)calloc( 1, sizeof(Lgm_TabularBessel) );
        Lgm_TabularBessel_Init( 200000, 14, -1.0, 14.5, t->BesselTable );
        t->BesselTableAlloced = 1;

    }

}

int J_N_Arr( int n, double x, double *JnArr, LgmTsyg2007_Info *tInfo  ) {

    int i;
    double xa, JDnArr[20];


    if ( tInfo->UseTabulatedBessel ) {

        xa = fabs( x );
        if ( xa < 14.0 ) {
// Check sign!!!!!
            Lgm_TabularBessel_Eval( xa, 4, JnArr, JDnArr, tInfo->BesselTable );
            return( 1 );
        }

    }


    // use Lau's routine
    //bessjlau( n, x, JnArr );
    //return( 1 );

/*
if (x > max_x ) {
    max_x = x;
    printf("min_x, max_x = %g %g\n", min_x, max_x);
}
if (x < min_x ) {
    min_x = x;
    printf("min_x, max_x = %g %g\n", min_x, max_x);
}
*/

    // use Jay Albert's NR mod
    bessjj( n, x, JnArr );
    return( 1 );

    // use gsl's array func
    //gsl_sf_bessel_Jn_array( 0, n, x, JnArr );
    //return( 1 );


    // use num recip.
    if ( n<0 ) {
        printf("Bad value of n: n = %d (must be >=0)\n", n);
        return( -1 );
    }

    JnArr[0] = bessj0( x );
    if ( n == 0 ) return( 1 );

    JnArr[1] = bessj1( x );
    if ( n == 1 ) return( 1 );

    for ( i=2; i<=n; i++ ) {
        // This is an unstable recurrence -- dont use!
        //JnArr[i] = 2.0*(i-1)*JnArr[i-1]/x - JnArr[i-2];

        JnArr[i] = bessj( i, x );
    }

    return( 1 );

}

int SinN_CosN_Arr( int n, double phi, double *SinArr, double *CosArr ) {

    double  g;
    int     i;


    if ( n<0 ) {
        printf("Bad value of n: n = %d (must be >=0)\n", n);
        return( -1 );
    }

    SinArr[0] = 0.0;
    CosArr[0] = 1.0;
    if ( n == 0 ) return( 1 );

    sincos( phi, &SinArr[1], &CosArr[1] );
    if ( n == 1 ) return( 1 );

    g = 2.0*CosArr[1];
    for ( i=2; i<=n; i++ ) {
        SinArr[i] = g*SinArr[i-1] - SinArr[i-2];
        CosArr[i] = g*CosArr[i-1] - CosArr[i-2];
    }

    return( 1 );

}


/*
 *  Converted to C by Michael G. Henderson (mghenderson@lanl.gov) Aug 10, 2012.
 */

void Lgm_SetCoeffs_TS07( long int Date, double UTC, LgmTsyg2007_Info *t ){

    int     k, year, month, day, doy, hour, minute, min5, isLeap;
    int     foundP=FALSE;
    double  fpart;
    char    Filename[1024], tmpstr[512];
    char    *p_str="Pdyn";
    FILE    *fp;
    const char* TS07_DATA_PATH = getenv("TS07_DATA_PATH");
    if (TS07_DATA_PATH==NULL) {
        TS07_DATA_PATH = LGM_TS07_DATA_DIR;
    }

    //get time and round to nearest 5 minutes... 
    Lgm_Doy(Date, &year, &month, &day, &doy);
    //TODO: should read two files and interpolate coeffs
    hour = (int)UTC;
    fpart = UTC - (int)UTC;
    minute = (int)(fpart*60.0);
    min5 = ((minute + 5/2) / 5) * 5;
    //at this point we can have minute 60, so we need to make sure the hours/minutes are right
    //solution must be robust to day and year boundaries
    if (min5>=60) { // if minute is 60 (or more) then recycle and increment hour
        min5 = min5 % 60;
        hour++;
    }
    if (hour>=24) {
        hour = hour % 24;
        doy++;
    }
    //to increment year, check if it's a leap year...
    isLeap = Lgm_LeapYear(year);
    if ((isLeap) && (doy>366)) {
        doy = 1;
        year++;
    }
    if ((!isLeap) && (doy>365)) {
        doy = 1;
        year++;
    }

    /*
     *  Read in coeffs
     */
    sprintf( Filename, "%s/Coeffs/%d_%03d/%d_%03d_%02d_%02d.par", TS07_DATA_PATH, year, doy, year, doy, hour, min5 );

    if ( (fp = fopen( Filename, "r" )) != NULL ) {

        for ( k=1; k<=101; k++ ) {
            fgets( &tmpstr, 512, fp);
            sscanf( &tmpstr, "%lf", &t->A[k] );
            //fscanf( fp, "%lf%*[\n]\n", &t->A[k] );
            //printf("t->A[%d] = %g\n", k, t->A[k]);
        }
        while ((!foundP) && (!feof(fp))) {
            fgets( &tmpstr, 512, fp);
            if ( strstr( &tmpstr, p_str) != NULL ) { //check line for Pdyn, if present read value
                sscanf( &tmpstr, "%*s %lf", &t->Pdyn);
                foundP = TRUE;
            }
        }
        //printf("t->Pdyn = %g\n", t->Pdyn);
        
        fclose(fp);

    } else {

        printf("Lgm_SetCoeffs_TS07(): Line %d in file %s. Could not open file %s\n", __LINE__, __FILE__, Filename );
	    //perror("Error");
        exit(-1);

    }

    //printf("\n********\nLgm_SetCoeffs_TS07(): Loaded %s for (date, time) = %d, %g\n********\n\n", Filename, Date, UTC );
    //printf("1. t->P[1][1], t->S_P[1][1] = %g, %g\n", t->P[1][1], t->S_P[1][1]);
}



/*
 * Copy from source ('s') to taget ('t')
 */
int Lgm_Copy_TS07_Info( LgmTsyg2007_Info *t, LgmTsyg2007_Info *s ){

    int i, j, k, N, M;

    if ( s == NULL) {
        printf("Lgm_Copy_TS07_Info: Error, source structure is NULL\n");
        return(-1);
    }


    /*
     * Do memcpy. Note that for things that are dynamically allocated in
     * LgmTsyg2007_Info structure, this will copy pointers to to memory that belong
     * to the source.  Need to take care of these things.
     */
    memcpy( t, s, sizeof(LgmTsyg2007_Info) );


    /*
     * Allocate memory for parameter arrays.
     */
    LGM_ARRAY_1D( t->A,   102, double );
    LGM_ARRAY_2D( t->TSS, 81, 6, double );
    LGM_ARRAY_3D( t->TSO, 81, 6, 5, double );
    LGM_ARRAY_3D( t->TSE, 81, 6, 5, double );
    t->ArraysAlloced = TRUE;

    for ( i=0; i<102; i++ ) { t->A[i] = s->A[i]; }
    for ( i=0; i<81; i++ ) { 
        for ( j=0; j<6; j++ ) { 
            t->TSS[i][j] = s->TSS[i][j]; 
            for ( k=0; k<5; k++ ) { 
                t->TSO[i][j][k] = s->TSO[i][j][k]; 
                t->TSE[i][j][k] = s->TSE[i][j][k]; 
            }
        }
    }

// This should work, but Ive had problems with it!?
//    memcpy( t->A,   s->A,      102*sizeof(double) );
//    memcpy( t->TSS, s->TSS,   81*6*sizeof(double) );
//    memcpy( t->TSO, s->TSO, 81*6*5*sizeof(double) );
//    memcpy( t->TSE, s->TSE, 81*6*5*sizeof(double) );



    if ( t->BesselTableAlloced ) {
        M = s->BesselTable->M;
        N = s->BesselTable->N;
        t->BesselTable = (Lgm_TabularBessel *)calloc( 1, sizeof(Lgm_TabularBessel) );
        t->BesselTable->M    = M;
        t->BesselTable->N    = N;
        t->BesselTable->xmin = s->BesselTable->xmin;
        t->BesselTable->xmax = s->BesselTable->xmax;
        t->BesselTable->s    = s->BesselTable->s;
        LGM_ARRAY_1D( t->BesselTable->JnTabular_x,       M, double );
        LGM_ARRAY_2D( t->BesselTable->JnTabular,    N+3, M, double );
        LGM_ARRAY_2D( t->BesselTable->JDnTabular,   N+3, M, double );
        LGM_ARRAY_2D( t->BesselTable->JDDnTabular,  N+3, M, double );
        for ( i=0; i<M; i++ ) {
            t->BesselTable->JnTabular_x[i] = s->BesselTable->JnTabular_x[i];
            for ( j=0; j<=(N+2); j++ ) {
                t->BesselTable->JnTabular[j][i]   = s->BesselTable->JnTabular[j][i];
                t->BesselTable->JDnTabular[j][i]  = s->BesselTable->JDnTabular[j][i];
                t->BesselTable->JDDnTabular[j][i] = s->BesselTable->JDDnTabular[j][i];
            }
        }
        //memcpy( t->BesselTable->JnTabular_x, s->BesselTable->JnTabular_x,    M*sizeof(double) );
        //memcpy( t->BesselTable->JnTabular,   s->BesselTable->JnTabular,   17*M*sizeof(double) );
        //memcpy( t->BesselTable->JDnTabular,  s->BesselTable->JDnTabular,  17*M*sizeof(double) );
        //memcpy( t->BesselTable->JDDnTabular, s->BesselTable->JDDnTabular, 17*M*sizeof(double) );
    }




    return( 1 );

}


void Lgm_Init_TS07( LgmTsyg2007_Info *t ){

    int     i, j, k;
    char    Filename[1024];
    FILE    *fp;
    const char* TS07_DATA_PATH = getenv("TS07_DATA_PATH");
//printf("[TS07] Got path: %s\n", TS07_DATA_PATH);
    if (TS07_DATA_PATH==NULL) {
        TS07_DATA_PATH = LGM_TS07_DATA_DIR;
    }


    // Init some params
    t->OLD_PS = -9e99;
    t->OLD_X  = -9e99;
    t->OLD_Y  = -9e99;
    t->OLD_Z  = -9e99;
    for (i=0; i<4; i++ ){
        t->DoneJ[i] = 0;
        t->S_DoneJ[i] = 0;
        for (j=0; j<4; j++ ){
            t->P[i][j] = -9e99;
            t->Q[i][j] = -9e99;
            t->R[i][j] = -9e99;
            t->S[i][j] = -9e99;
            t->S_P[i][j] = -9e99;
            t->S_Q[i][j] = -9e99;
            t->S_R[i][j] = -9e99;
            t->S_S[i][j] = -9e99;
        }
    }


    /*
     * Allocate memory for parameter arrays.
     */
    LGM_ARRAY_1D( t->A,   102, double );
    LGM_ARRAY_2D( t->TSS, 81, 6, double );
    LGM_ARRAY_3D( t->TSO, 81, 6, 5, double );
    LGM_ARRAY_3D( t->TSE, 81, 6, 5, double );
    t->ArraysAlloced = TRUE;




    /*
     *  Read in the TSS, TSO, and TSE .par files
     */
    for ( i=1; i<=5; i++ ) {

        sprintf( Filename, "%s/TAIL_PAR/tailamebhr%1d.par", TS07_DATA_PATH, i );
        if ( (fp = fopen( Filename, "r" )) != NULL ) {

            for ( k=1; k<=80; k++ ) fscanf( fp, "%lf", &t->TSS[k][i] );
	        fclose(fp);

        } else {

            printf("Lgm_Init_TS07(): Line %d in file %s. Could not open file %s\n", __LINE__, __FILE__, Filename );
	        perror("Error");
            exit(-1);

        }

    }

    for ( i=1; i<=5; i++ ) {
        for ( j=1; j<=4; j++ ) {

            sprintf( Filename, "%s/TAIL_PAR/tailamhr_o_%1d%1d.par", TS07_DATA_PATH, i, j );
            if ( (fp = fopen( Filename, "r" )) != NULL ) {
                for ( k=1; k<=80; k++ ) fscanf( fp, "%lf", &t->TSO[k][i][j] );
		        fclose(fp);

            } else {

                printf("Lgm_Init_TS07(): Line %d in file %s. Could not open file %s\n", __LINE__, __FILE__, Filename );
		        perror("Error");
                exit(-1);

            }

        }
    }

    for ( i=1; i<=5; i++ ) {
        for ( j=1; j<=4; j++ ) {

            sprintf( Filename, "%s/TAIL_PAR/tailamhr_e_%1d%1d.par", TS07_DATA_PATH, i, j );
            if ( (fp = fopen( Filename, "r" )) != NULL ) {

                for ( k=1; k<=80; k++ ) fscanf( fp, "%lf", &t->TSE[k][i][j] );
		        fclose(fp);

            } else {

                printf("Lgm_Init_TS07(): Line %d in file %s. Could not open file %s\n", __LINE__, __FILE__, Filename );
		        perror("Error");
                exit(-1);

            }

        }
    }


    /*
     * Init the cache "flags for SHTBNORM func
     */
    t->SHTBNORM_S_RHO_LAST = LGM_FILL_VALUE;
    t->SHTBNORM_E_RHO_LAST = LGM_FILL_VALUE;
    t->SHTBNORM_O_RHO_LAST = LGM_FILL_VALUE;


    t->BesselTableAlloced = 0;


    return;

}

void Lgm_DeAllocate_TS07( LgmTsyg2007_Info *t ){

    if ( t->ArraysAlloced == TRUE ) {
        LGM_ARRAY_1D_FREE( t->A );
        LGM_ARRAY_2D_FREE( t->TSS );
        LGM_ARRAY_3D_FREE( t->TSO );
        LGM_ARRAY_3D_FREE( t->TSE );
        t->ArraysAlloced = FALSE;
    }

    if ( t->BesselTableAlloced ) {
        Lgm_TabularBessel_Free( t->BesselTable );
        t->BesselTableAlloced = FALSE;
    }

}


/*
 *  CALCULATES GSM COMPONENTS OF THE SHIELDING FIELD FOR THE EARTH'S DIPOLE
 *    WITHIN A MAGNETOPAUSE WITH VARIABLE SHAPE AND SIZE
 */
void Tsyg_TS07( int IOPT, double *PARMOD, double PS, double SINPS, double COSPS, double X, double Y, double Z, double *BX, double *BY, double *BZ, LgmTsyg2007_Info *tInfo ) {

    double  BXTS[6], BXTO[6][5], BXTE[6][5];
    double  BYTS[6], BYTO[6][5], BYTE[6][5];
    double  BZTS[6], BZTO[6][5], BZTE[6][5];
    double  BXCF, BYCF, BZCF;
    double  BXR11, BYR11, BZR11, BXR12, BYR12, BZR12;
    double  BXR21a, BYR21a, BZR21a, BXR21s, BYR21s, BZR21s;
    int     NTOT=101;
    double  PDYN, *A;


    PDYN = PARMOD[1]; // following other versions...
    A = tInfo->A;

    tInfo->sin_psi_op = SINPS;
    tInfo->cos_psi_op = COSPS;

    TS07D_EXTERN( 0, A, NTOT, PS, PDYN, X, Y, Z, &BXCF, &BYCF, &BZCF, BXTS, BYTS, BZTS,
            BXTO, BYTO, BZTO, BXTE, BYTE, BZTE, &BXR11, &BYR11, &BZR11, &BXR12, &BYR12, &BZR12,
            &BXR21a, &BYR21a, &BZR21a, &BXR21s, &BYR21s, &BZR21s, BX, BY, BZ, tInfo );

    //printf("BXCF,   BYCF,   BZCF   = %.8lf %.8lf %.8lf\n", BXCF, BYCF, BZCF );
    //printf("BXR11,  BYR11,  BZR11  = %.8lf %.8lf %.8lf\n", BXR11, BYR11, BZR11 );
    //printf("BXR12,  BYR12,  BZR12  = %.8lf %.8lf %.8lf\n", BXR12, BYR12, BZR12 );
    //printf("BXR21a, BYR21a, BZR21a = %.8lf %.8lf %.8lf\n", BXR21a, BYR21a, BZR21a );
    //printf("BXR21s, BYR21s, BZR21s = %.8lf %.8lf %.8lf\n", BXR21s, BYR21s, BZR21s );
    //printf("BX,     BY,     BZ     = %.8lf %.8lf %.8lf\n", *BX, *BY, *BZ);



    return;

}






/*
 *  IOPGEN - GENERAL OPTION FLAG:  IOPGEN=0 - CALCULATE TOTAL FIELD
 *                                 IOPGEN=1 - DIPOLE SHIELDING ONLY
 *                                 IOPGEN=2 - TAIL FIELD ONLY
 *                                 IOPGEN=3 - BIRKELAND FIELD ONLY
 *                                 IOPGEN=4 - RING CURRENT FIELD ONLY
 */
void TS07D_EXTERN( int IOPGEN, double *A, int NTOT, double PS, double PDYN, double X, double Y, double Z, double *BXCF, double *BYCF,
        double *BZCF, double BXTS[], double BYTS[], double BZTS[], double BXTO[][5], double BYTO[][5], double BZTO[][5], double BXTE[][5], double BYTE[][5],
        double BZTE[][5], double *BXR11, double *BYR11, double *BZR11, double *BXR12, double *BYR12, double *BZR12, double *BXR21a, double *BYR21a,
        double *BZR21a, double *BXR21s, double *BYR21s, double *BZR21s, double *BX, double *BY, double *BZ,  LgmTsyg2007_Info *tInfo ) {

    int       done, K, L, IND;
    double    XAPPA, XAPPA2, XAPPA3, SPS, X0, AM, S0, FACTIMF, OIMFX, OIMFY, OIMFZ, R, XSS, ZSS;
    double    XSOLD, ZSOLD, ZSSoR, ZSSoR2, RH, RoRH, RoRH2, RoRH3, SINPSAS, SINPSAS2, COSPSAS, DD, RHO2, ASQ;
    double    XMXM, AXX0, ARO, AROpAXX0, AROpAXX02, SIGMA, CFX, CFY, CFZ, DSTT, ZNAM;
    double    DLP1, DLP2, TAMP1, TAMP2, A_SRC, A_PRC, A_R11, XX, YY, ZZ;
    double    A_R21, A_R12, A_R21a, A_R21s, QX, QY, QZ, FINT, FEXT, ZNAM05, ooZNAM20, BBX, BBY, BBZ;
    double    BXR22a, BYR22a, BZR22a, BXR22s, BYR22s, BZR22s;
    double    BXR11s, BYR11s, BZR11s, BXR12s, BYR12s, BZR12s;
    double    TX, TY, TZ, PDYN_0, P_FACTOR;


    double    A0_A=34.586, A0_S0=1.1960, A0_X0=3.4397;     // SHUE ET AL. PARAMETERS
    double    DSIG=0.005, RH2=-5.2;



    //if ( PDYN != tInfo->OLD_PDYN ) {
        //XAPPA = mypow( 0.5*PDYN, 0.155 );   //  0.155 is the value obtained in TS05
        XAPPA = pow( 0.5*PDYN, 0.155 );   //  0.155 is the value obtained in TS05
        tInfo->XAPPA = XAPPA;
    //} else {
    //    XAPPA = tInfo->XAPPA;
    //}

    XAPPA2 = XAPPA*XAPPA;
    XAPPA3 = XAPPA*XAPPA2;


    tInfo->CB_TAIL.D          = A[96];      // FORWARDS TAIL SHEET THICKNESS
    tInfo->CB_RH0.RH0         = A[97];
    tInfo->CB_G.G             = A[98];
    tInfo->CB_BIRKPAR.XKAPPA1 = A[99];      // SCALING FACTORS FOR BIRKELAND CURRENTS
    tInfo->CB_BIRKPAR.XKAPPA2 = A[100];     // SCALING FACTORS FOR BIRKELAND CURRENTS
    tInfo->CB_G.TW            = A[101];     //  THIS PARAMETER CONTROLS THE IMF-INDUCED TWISTING (ADDED 04/21/06)

    XX = X*XAPPA;        //  pressure scaling has been reinstated here
    YY = Y*XAPPA;
    ZZ = Z*XAPPA;


    //printf("%g %g\n", XAPPA,PDYN);

    //SPS = sin( PS );
    SPS = tInfo->sin_psi_op;

    X0 = A0_X0/XAPPA;   // pressure scaling has been reinstated, even though these parameters are not used in this code
    AM = A0_A/XAPPA;    // pressure scaling has been reinstated, even though these parameters are not used in this code
    S0 = A0_S0;


    if ( IOPGEN <= 1 ) {
        TS07D_SHLCAR3X3( XX, YY, ZZ, PS, &CFX, &CFY, &CFZ, tInfo ); //  DIPOLE SHIELDING FIELD

        *BXCF = CFX*XAPPA3;
        *BYCF = CFY*XAPPA3;
        *BZCF = CFZ*XAPPA3;
    } else {
        *BXCF = 0.0;
        *BYCF = 0.0;
        *BZCF = 0.0;
    }



    if ( (IOPGEN == 0) || (IOPGEN == 2) ) {

      TS07D_DEFORMED( PS, XX, YY, ZZ, BXTS, BYTS, BZTS, BXTO, BYTO, BZTO, BXTE, BYTE, BZTE, tInfo ); // TAIL FIELD (THREE MODES)
//printf("BXTS[2], BYTS[2], BZTS[2] = %g %g %g\n", BXTS[2], BYTS[2], BZTS[2]);

    } else {

        for ( K=1; K<=5; K++ ) {
            BXTS[K] = 0.0;
            BYTS[K] = 0.0;
            BZTS[K] = 0.0;
        }

        for ( K=1; K<=5; K++ ) {
            for ( L=1; L<=4; L++ ) {

                BXTO[K][L] = BYTO[K][L] = BZTO[K][L] = 0.0;
                BXTE[K][L] = BYTE[K][L] = BZTE[K][L] = 0.0;

            }
        }
    }




    if ( (IOPGEN == 0) || (IOPGEN == 3) ) {

         TS07D_BIRK_TOT( PS,XX,YY,ZZ,BXR11,BYR11,BZR11,BXR12,BYR12,
                   BZR12,BXR21a,BYR21a,BZR21a,&BXR22a,&BYR22a,&BZR22a, tInfo );          // BIRKELAND FIELD (TWO MODES FOR R1 AND TWO MODES FOR R2)

         TS07D_BIRTOTSY( PS,XX,YY,ZZ,&BXR11s,&BYR11s,&BZR11s, &BXR12s,
                   &BYR12s,&BZR12s,BXR21s,BYR21s,BZR21s,&BXR22s,&BYR22s,&BZR22s, tInfo );  //   "SYMMETRIC" BIRKELAND FIELD
                                                                               // (TWO MODES FOR R1s AND TWO MODES FOR R2s)
                                                                               // (but we actually use from here only R2s modes)
    } else {

         *BXR11  = *BYR11  = *BZR11  = 0.0;
         *BXR12  = *BYR12  = *BZR12  = 0.0;
         *BXR21a = *BYR21a = *BZR21a = 0.0;
         *BXR21s = *BYR21s = *BZR21s = 0.0;

    }


    /*
     *  NOW, ADD UP ALL THE COMPONENTS:
     */
    A_R11  = A[92];
    A_R12  = A[93];
    A_R21a = A[94];
    A_R21s = A[95];

    TX = TY = TZ = 0.0;


    // --- New tail structure -------------
    PDYN_0 = 2.0;   //   AVERAGE PRESSURE USED FOR NORMALIZATION
    P_FACTOR = sqrt( PDYN/PDYN_0 ) - 1.0;
    IND = 1;


    for ( K=1; K<=5; K++ ){
        ++IND;
        TX += ( A[IND] + A[IND+45]*P_FACTOR )*BXTS[K];  //   2 - 6  &  47 - 51
        TY += ( A[IND] + A[IND+45]*P_FACTOR )*BYTS[K];
        TZ += ( A[IND] + A[IND+45]*P_FACTOR )*BZTS[K];
//printf("TX, A[IND], BXTS[K]  = %g %g %g\n", TX, A[IND], BXTS[K] );
    }


    for ( K=1; K<=5; K++ ){
        for ( L=1; L<=4; L++ ){

            ++IND;

            TX += ( A[IND] + A[IND+45]*P_FACTOR )*BXTO[K][L];    //   7 -26  &  52 - 71
            TY += ( A[IND] + A[IND+45]*P_FACTOR )*BYTO[K][L];
            TZ += ( A[IND] + A[IND+45]*P_FACTOR )*BZTO[K][L];

            TX += ( A[IND+20] + A[IND+65]*P_FACTOR )*BXTE[K][L]; //   27 -46  &  72 - 91
            TY += ( A[IND+20] + A[IND+65]*P_FACTOR )*BYTE[K][L];
            TZ += ( A[IND+20] + A[IND+65]*P_FACTOR )*BZTE[K][L];

        }
    }


    BBX = A[1]*(*BXCF) + TX + A_R11*(*BXR11) + A_R12*(*BXR12) + A_R21a*(*BXR21a) + A_R21s*(*BXR21s);
    BBY = A[1]*(*BYCF) + TY + A_R11*(*BYR11) + A_R12*(*BYR12) + A_R21a*(*BYR21a) + A_R21s*(*BYR21s);
    BBZ = A[1]*(*BZCF) + TZ + A_R11*(*BZR11) + A_R12*(*BZR12) + A_R21a*(*BZR21a) + A_R21s*(*BZR21s);

//printf("TX, P_FACTOR = %g %g\n", TX, P_FACTOR );


    *BX = BBX;
    *BY = BBY;
    *BZ = BBZ;


    return;

}


/*
 *
 *    THIS S/R RETURNS THE SHIELDING FIELD FOR THE EARTH'S DIPOLE,
 *    REPRESENTED BY  2x3x3=18 "CARTESIAN" HARMONICS, tilted with respect
 *    to the z=0 plane (see NB#4, p.74-74)
 *
 *  - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 *   The 36 coefficients enter in pairs in the amplitudes of the "cartesian"
 *     harmonics (A(1)-A(36).
 *   The 14 nonlinear parameters (A(37)-A(50) are the scales Pi,Ri,Qi,and Si
 *    entering the arguments of exponents, sines, and cosines in each of the
 *    18 "Cartesian" harmonics  PLUS TWO TILT ANGLES FOR THE CARTESIAN HARMONICS
 *        (ONE FOR THE PSI=0 MODE AND ANOTHER FOR THE PSI=90 MODE)
 *  - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 */
void    TS07D_SHLCAR3X3( double X, double Y, double Z, double PS, double *BX, double *BY, double *BZ, LgmTsyg2007_Info *tInfo ) {

    static double    A[] = { -9e99, -901.2327248,895.8011176,817.6208321,-845.5880889,
                      -83.73539535,86.58542841,336.8781402,-329.3619944,-311.2947120,
                      308.6011161,31.94469304,-31.30824526,125.8739681,-372.3384278,
                      -235.4720434,286.7594095,21.86305585,-27.42344605,-150.4874688,
                      2.669338538,1.395023949,-.5540427503,-56.85224007,3.681827033,
                      -43.48705106,5.103131905,1.073551279,-.6673083508,12.21404266,
                      4.177465543,5.799964188,-.3977802319,-1.044652977,.5703560010,
                      3.536082962,-3.222069852,9.620648151,6.082014949,27.75216226,
                      12.44199571,5.122226936,6.982039615,20.12149582,6.150973118,
                      4.663639687,15.73319647,2.303504968,5.840511214,.8385953499E-01,
                      .3477844929 };

    double    P1, P2, P3, ooP1, ooP2, ooP3, ooP1_2, ooP2_2, ooP3_2;
    double    R1, R2, R3, ooR1, ooR2, ooR3, ooR1_2, ooR2_2, ooR3_2;
    double    Q1, Q2, Q3, ooQ1, ooQ2, ooQ3, ooQ1_2, ooQ2_2, ooQ3_2;
    double    S1, S2, S3, ooS1, ooS2, ooS3, ooS1_2, ooS2_2, ooS3_2;
    double    T1, T2;
    double    sqrtP1R1, sqrtP1R2, sqrtP1R3, sqrtP2R1, sqrtP2R2, sqrtP2R3, sqrtP3R1, sqrtP3R2, sqrtP3R3;
    double    sqrtQ1S1, sqrtQ1S2, sqrtQ1S3, sqrtQ2S1, sqrtQ2S2, sqrtQ2S3, sqrtQ3S1, sqrtQ3S2, sqrtQ3S3;
    double    YoP1, YoP2, YoP3, Z1oR1, Z1oR2, Z1oR3;
    double    YoQ1, YoQ2, YoQ3, Z2oS1, Z2oS2, Z2oS3;
    double    COSZ1oR1, COSZ1oR2, COSZ1oR3, SINZ1oR1, SINZ1oR2, SINZ1oR3;
    double    COSZ2oS1, COSZ2oS2, COSZ2oS3, SINZ2oS1, SINZ2oS2, SINZ2oS3;
    double    COSYoP1, COSYoP2, COSYoP3, SINYoP1, SINYoP2, SINYoP3;
    double    COSYoQ1, COSYoQ2, COSYoQ3, SINYoQ1, SINYoQ2, SINYoQ3;
    double    X1, X2, Z1, Z2, SQPR, EXPR, SQQS, EXQS, CPS, SPS, S2PS, ooSQPR;
    double    ST1, ST2, CT1, CT2;
    double    FX1, HY1, FZ1, HX1, HZ1, FX2, HY2, FZ2, HX2, HZ2, FX3, HY3, FZ3, HX3, HZ3;
    double    FX4, HY4, FZ4, HX4, HZ4, FX5, HY5, FZ5, HX5, HZ5, FX6, HY6, FZ6, HX6, HZ6;
    double    FX7, HY7, FZ7, HX7, HZ7, FX8, HY8, FZ8, HX8, HZ8, FX9, HY9, FZ9, HX9, HZ9;
    double    A1, A2, A3, A4, A5, A6, A7, A8, A9;

    P1 = A[37]; P2 = A[38]; P3=A[39];
    R1 = A[40]; R2 = A[41]; R3=A[42];
    Q1 = A[43]; Q2 = A[44]; Q3=A[45];
    S1 = A[46]; S2 = A[47]; S3=A[48];
    T1 = A[49]; T2 = A[50];

    //CPS = cos(PS); SPS = sin(PS); S2PS = 2.0*CPS;
    CPS = tInfo->cos_psi_op; SPS = tInfo->sin_psi_op; S2PS = 2.0*CPS; // MODIFIED HERE (INSTEAD OF SIN(3*PS) I TRY SIN(2*PS)


    //ST1 = sin( PS*T1 ); CT1 = cos( PS*T1 );
    sincos( PS*T1, &ST1, &CT1 );
    //ST2 = sin( PS*T2 ); CT2 = cos( PS*T2 );
    sincos( PS*T2, &ST2, &CT2 );

    X1 = X*CT1 - Z*ST1; Z1 = X*ST1 + Z*CT1;
    X2 = X*CT2 - Z*ST2; Z2 = X*ST2 + Z*CT2;


    /*
     *   MAKE THE TERMS IN THE 1ST SUM ("PERPENDICULAR" SYMMETRY):
      *
     *   MGH - the original fortran code repeats many calculations.
     *         In addition, there are far more divides than necessary.
     *           Very expensive. (Is there a recurrence relation for any of this?)
     *         sqrt()'s andm trigs are very expensive so dont do more than we need.
     */
// do these change? Isnt this wasteful if they dont?
// especially the sqrt() below...
    ooP1 = 1.0/P1; ooP1_2 = ooP1*ooP1;
    ooP2 = 1.0/P2; ooP2_2 = ooP2*ooP2;
    ooP3 = 1.0/P3; ooP3_2 = ooP3*ooP3;
    ooR1 = 1.0/R1; ooR1_2 = ooR1*ooR1;
    ooR2 = 1.0/R2; ooR2_2 = ooR2*ooR2;
    ooR3 = 1.0/R3; ooR3_2 = ooR3*ooR3;

// are these constants?!
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


    //COSZ1oR1 = cos(Z1oR1); SINZ1oR1 = sin(Z1oR1);
    sincos( Z1oR1, &SINZ1oR1, &COSZ1oR1 );
    //COSZ1oR2 = cos(Z1oR2); SINZ1oR2 = sin(Z1oR2);
    sincos( Z1oR2, &SINZ1oR2, &COSZ1oR2 );
    //COSZ1oR3 = cos(Z1oR3); SINZ1oR3 = sin(Z1oR3);
    sincos( Z1oR3, &SINZ1oR3, &COSZ1oR3 );

    //COSYoP1 = cos(YoP1); SINYoP1 = sin(YoP1);
    sincos( YoP1, &SINYoP1, &COSYoP1 );
    //COSYoP2 = cos(YoP2); SINYoP2 = sin(YoP2);
    sincos( YoP2, &SINYoP2, &COSYoP2 );
    //COSYoP3 = cos(YoP3); SINYoP3 = sin(YoP3);
    sincos( YoP3, &SINYoP3, &COSYoP3 );


    // I = 1
    SQPR =  sqrtP1R1; EXPR =  exp(SQPR*X1);
    FX1  = -SQPR*EXPR*COSYoP1*SINZ1oR1;
    HY1  =  EXPR*ooP1*SINYoP1*SINZ1oR1;
    FZ1  = -EXPR*COSYoP1*ooR1*COSZ1oR1;
    HX1  =  FX1*CT1 + FZ1*ST1;
    HZ1  = -FX1*ST1 + FZ1*CT1;

    SQPR =  sqrtP1R2; EXPR =  exp(SQPR*X1);
    FX2  = -SQPR*EXPR*COSYoP1*SINZ1oR2;
    HY2  =  EXPR*ooP1*SINYoP1*SINZ1oR2;
    FZ2  = -EXPR*COSYoP1*ooR2*COSZ1oR2;
    HX2  =  FX2*CT1+FZ2*ST1;
    HZ2  = -FX2*ST1+FZ2*CT1;

    SQPR =  sqrtP1R3; ooSQPR = 1.0/SQPR; EXPR =  exp(SQPR*X1);
    FX3  = -EXPR*COSYoP1*(SQPR*Z1*COSZ1oR3 + SINZ1oR3*ooR3*(X1+ooSQPR));
    HY3  =  EXPR*ooP1*SINYoP1*(Z1*COSZ1oR3 + X1*ooR3*SINZ1oR3*ooSQPR);
    FZ3  = -EXPR*COSYoP1*( COSZ1oR3*(1.0 + X1*ooR3_2*ooSQPR) - Z1*ooR3*SINZ1oR3);
    HX3  =  FX3*CT1 + FZ3*ST1;
    HZ3  = -FX3*ST1 + FZ3*CT1;

    // I = 2
    SQPR =  sqrtP2R1; EXPR =  exp(SQPR*X1);
    FX4  = -SQPR*EXPR*COSYoP2*SINZ1oR1;
    HY4  =  EXPR*ooP2*SINYoP2*SINZ1oR1;
    FZ4  = -EXPR*COSYoP2*ooR1*COSZ1oR1;
    HX4  =  FX4*CT1+FZ4*ST1;
    HZ4  = -FX4*ST1+FZ4*CT1;

    SQPR =  sqrtP2R2; EXPR =  exp(SQPR*X1);
    FX5  = -SQPR*EXPR*COSYoP2*SINZ1oR2;
    HY5  =  EXPR*ooP2*SINYoP2*SINZ1oR2;
    FZ5  = -EXPR*COSYoP2*ooR2*COSZ1oR2;
    HX5  =  FX5*CT1+FZ5*ST1;
    HZ5  = -FX5*ST1+FZ5*CT1;

    SQPR =  sqrtP2R3; ooSQPR = 1.0/SQPR; EXPR =  exp(SQPR*X1);
    FX6  = -EXPR*COSYoP2*(SQPR*Z1*COSZ1oR3+SINZ1oR3/R3*(X1+ooSQPR));
    HY6  =  EXPR*ooP2*SINYoP2*(Z1*COSZ1oR3+X1*ooR3*SINZ1oR3*ooSQPR);
    FZ6  = -EXPR*COSYoP2*(COSZ1oR3*(1.0+X1*ooR3_2*ooSQPR)-Z1*ooR3*SINZ1oR3);
    HX6  =  FX6*CT1+FZ6*ST1;
    HZ6  = -FX6*ST1+FZ6*CT1;


    // I = 3
    SQPR =  sqrtP3R1; EXPR =  exp(SQPR*X1);
    FX7  = -SQPR*EXPR*COSYoP3*SINZ1oR1;
    HY7  =  EXPR*ooP3*SINYoP3*SINZ1oR1;
    FZ7  = -EXPR*COSYoP3*ooR1*COSZ1oR1;
    HX7  =  FX7*CT1+FZ7*ST1;
    HZ7  = -FX7*ST1+FZ7*CT1;

    SQPR =  sqrtP3R2; EXPR =  exp(SQPR*X1);
    FX8  = -SQPR*EXPR*COSYoP3*SINZ1oR2;
    HY8  =  EXPR*ooP3*SINYoP3*SINZ1oR2;
    FZ8  = -EXPR*COSYoP3*ooR2*COSZ1oR2;
    HX8  =  FX8*CT1+FZ8*ST1;
    HZ8  = -FX8*ST1+FZ8*CT1;

    SQPR =  sqrtP3R3; ooSQPR = 1.0/SQPR; EXPR =  exp(SQPR*X1);
    FX9  = -EXPR*COSYoP3*(SQPR*Z1*COSZ1oR3+SINZ1oR3*ooR3*(X1+ooSQPR));
    HY9  =  EXPR*ooP3*SINYoP3*(Z1*COSZ1oR3+X1*ooR3*SINZ1oR3*ooSQPR);
    FZ9  = -EXPR*COSYoP3*(COSZ1oR3*(1.0+X1*ooR3_2*ooSQPR)-Z1*ooR3*SINZ1oR3);
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

    //COSZ2oS1 = cos(Z2oS1); SINZ2oS1 = sin(Z2oS1);
    sincos( Z2oS1, &SINZ2oS1, &COSZ2oS1 );
    //COSZ2oS2 = cos(Z2oS2); SINZ2oS2 = sin(Z2oS2);
    sincos( Z2oS2, &SINZ2oS2, &COSZ2oS2 );
    //COSZ2oS3 = cos(Z2oS3); SINZ2oS3 = sin(Z2oS3);
    sincos( Z2oS3, &SINZ2oS3, &COSZ2oS3 );

    //COSYoQ1 = cos(YoQ1); SINYoQ1 = sin(YoQ1);
    sincos( YoQ1, &SINYoQ1, &COSYoQ1 );
    //COSYoQ2 = cos(YoQ2); SINYoQ2 = sin(YoQ2);
    sincos( YoQ2, &SINYoQ2, &COSYoQ2 );
    //COSYoQ3 = cos(YoQ3); SINYoQ3 = sin(YoQ3);
    sincos( YoQ3, &SINYoQ3, &COSYoQ3 );

    // I = 1
    SQQS =  sqrtQ1S1; EXQS =  exp(SQQS*X2);
    FX1  = -SQQS*EXQS*COSYoQ1*COSZ2oS1 *SPS;
    HY1  =  EXQS*ooQ1*SINYoQ1*COSZ2oS1   *SPS;
    FZ1  =  EXQS*COSYoQ1*ooS1*SINZ2oS1   *SPS;
    HX1  =  FX1*CT2+FZ1*ST2;
    HZ1  = -FX1*ST2+FZ1*CT2;

    SQQS =  sqrtQ1S2; EXQS =  exp(SQQS*X2);
    FX2  = -SQQS*EXQS*COSYoQ1*COSZ2oS2 *SPS;
    HY2  =  EXQS*ooQ1*SINYoQ1*COSZ2oS2   *SPS;
    FZ2  =  EXQS*COSYoQ1*ooS2*SINZ2oS2   *SPS;
    HX2  =  FX2*CT2+FZ2*ST2;
    HZ2  = -FX2*ST2+FZ2*CT2;

    SQQS =  sqrtQ1S3; EXQS =  exp(SQQS*X2);
    FX3  = -SQQS*EXQS*COSYoQ1*COSZ2oS3 *SPS;
    HY3  =  EXQS*ooQ1*SINYoQ1*COSZ2oS3   *SPS;
    FZ3  =  EXQS*COSYoQ1*ooS3*SINZ2oS3   *SPS;
    HX3  =  FX3*CT2+FZ3*ST2;
    HZ3  = -FX3*ST2+FZ3*CT2;


    // I = 2
    SQQS =  sqrtQ2S1; EXQS =  exp(SQQS*X2);
    FX4  = -SQQS*EXQS*COSYoQ2*COSZ2oS1 *SPS;
    HY4  =  EXQS*ooQ2*SINYoQ2*COSZ2oS1   *SPS;
    FZ4  =  EXQS*COSYoQ2*ooS1*SINZ2oS1   *SPS;
    HX4  =  FX4*CT2+FZ4*ST2;
    HZ4  = -FX4*ST2+FZ4*CT2;

    SQQS =  sqrtQ2S2; EXQS =  exp(SQQS*X2);
    FX5  = -SQQS*EXQS*COSYoQ2*COSZ2oS2 *SPS;
    HY5  =  EXQS*ooQ2*SINYoQ2*COSZ2oS2   *SPS;
    FZ5  =  EXQS*COSYoQ2*ooS2*SINZ2oS2   *SPS;
    HX5  =  FX5*CT2+FZ5*ST2;
    HZ5  = -FX5*ST2+FZ5*CT2;

    SQQS =  sqrtQ2S3; EXQS =  exp(SQQS*X2);
    FX6  = -SQQS*EXQS*COSYoQ2*COSZ2oS3 *SPS;
    HY6  =  EXQS*ooQ2*SINYoQ2*COSZ2oS3   *SPS;
    FZ6  =  EXQS*COSYoQ2*ooS3*SINZ2oS3   *SPS;
    HX6  =  FX6*CT2+FZ6*ST2;
    HZ6  = -FX6*ST2+FZ6*CT2;

    // I = 3
    SQQS =  sqrtQ3S1; EXQS =  exp(SQQS*X2);
    FX7  = -SQQS*EXQS*COSYoQ3*COSZ2oS1 *SPS;
    HY7  =  EXQS*ooQ3*SINYoQ3*COSZ2oS1   *SPS;
    FZ7  =  EXQS*COSYoQ3*ooS1*SINZ2oS1   *SPS;
    HX7  =  FX7*CT2+FZ7*ST2;
    HZ7  = -FX7*ST2+FZ7*CT2;

    SQQS =  sqrtQ3S2; EXQS =  exp(SQQS*X2);
    FX8  = -SQQS*EXQS*COSYoQ3*COSZ2oS2 *SPS;
    HY8  =  EXQS*ooQ3*SINYoQ3*COSZ2oS2   *SPS;
    FZ8  =  EXQS*COSYoQ3*ooS2*SINZ2oS2   *SPS;
    HX8  =  FX8*CT2+FZ8*ST2;
    HZ8  = -FX8*ST2+FZ8*CT2;

    SQQS =  sqrtQ3S3; EXQS =  exp(SQQS*X2);
    FX9  = -SQQS*EXQS*COSYoQ3*COSZ2oS3 *SPS;
    HY9  =  EXQS*ooQ3*SINYoQ3*COSZ2oS3   *SPS;
    FZ9  =  EXQS*COSYoQ3*ooS3*SINZ2oS3   *SPS;
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



/*
 *
 *
 *   CALCULATES GSM COMPONENTS OF 104 UNIT-AMPLITUDE TAIL FIELD MODES,
 *   TAKING INTO ACCOUNT BOTH EFFECTS OF DIPOLE TILT:
 *   WARPING IN Y-Z (DONE BY THE S/R WARPED) AND BENDING IN X-Z (DONE BY THIS SUBROUTINE)
 */
void    TS07D_DEFORMED( double PS, double X, double Y, double Z,
        double BXS[], double BYS[], double BZS[], double BXO[][5], double BYO[][5], double BZO[][5],
        double BXE[][5], double BYE[][5], double BZE[][5], LgmTsyg2007_Info *tInfo ){


    double BXASS[6],BXASO[6][5],BXASE[6][5];
    double BYASS[6],BYASO[6][5],BYASE[6][5];
    double BZASS[6],BZASO[6][5],BZASE[6][5];

    static double RH2=-5.2;
    static int    IEPS=3;

    int     K, L;
    double  RH0, SPS, CPS, X2, Y2, Z2, R, ooR, ZR, ZR2, RH, DRHDR, DRHDZ, RRH, RRH2, RRH3;
    double  F, F2, F4, DFDR, DFDRH, SPSAS, SPSAS2, CPSAS, XAS, ZAS, ooCPSAS, FACPS, PSASX;
    double  PSASY, PSASZ, DXASDX, DXASDY, DXASDZ, DZASDX, DZASDY, DZASDZ, FAC1, FAC2, FAC3;
    double  R2;

    RH0 = tInfo->CB_RH0.RH0;
//printf("RH0 = %g\n", RH0);

    /*
     *  RH0,RH1,RH2, AND IEPS CONTROL THE TILT-RELATED DEFORMATION OF THE TAIL FIELD
     */
   //SPS = sin(PS); SPS2 = SPS*SPS;
    SPS = tInfo->sin_psi_op; CPS = tInfo->cos_psi_op;
    X2 = X*X; Y2 = Y*Y; Z2 = Z*Z; R2 = X2+Y2+Z2;
    R = sqrt(R2); ooR = 1.0/R;
    ZR = Z*ooR; ZR2 = ZR*ZR;
    RH = tInfo->CB_RH0.RH0+RH2*ZR2;
    DRHDR = -ZR*ooR*2.0*RH2*ZR;
    DRHDZ =  2.0*RH2*ZR*ooR;

    RRH = R/RH; RRH2 = RRH*RRH;  RRH3 = RRH2*RRH;

//could create a LUT for this pow.
    F = cbrt( 1.0/(1.0+RRH3) ); F2 = F*F; F4 = F2*F2;
    DFDR = -RRH2*F4/RH;
    DFDRH = -RRH*DFDR;

    SPSAS = SPS*F; SPSAS2 = SPSAS*SPSAS;
    CPSAS = sqrt(1.0-SPSAS2);

    XAS = X*CPSAS-Z*SPSAS;
    ZAS = X*SPSAS+Z*CPSAS;

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
     *     DEFORM:
     */
    TS07D_WARPED( PS, XAS, Y, ZAS, BXASS, BYASS, BZASS, BXASO, BYASO, BZASO, BXASE, BYASE, BZASE, tInfo );
//printf("BXASS[2], BYASS[2], BZASS[2] = %g %g %g\n", BXASS[2], BYASS[2], BZASS[2]);



    /*
     *   --- New tail structure -------------
     */
    for ( K=1; K<=5; K++ ) {
        BXS[K] = BXASS[K]*DZASDZ - BZASS[K]*DXASDZ + BYASS[K]*FAC1;
        BYS[K] = BYASS[K]*FAC2;
        BZS[K] = BZASS[K]*DXASDX - BXASS[K]*DZASDX + BYASS[K]*FAC3;
    }

    for ( K=1; K<=5; K++ ) {
        for ( L=1; L<=4; L++ ) {

          BXO[K][L] = BXASO[K][L]*DZASDZ - BZASO[K][L]*DXASDZ + BYASO[K][L]*FAC1;
          BYO[K][L] = BYASO[K][L]*FAC2;
          BZO[K][L] = BZASO[K][L]*DXASDX - BXASO[K][L]*DZASDX + BYASO[K][L]*FAC3;

          BXE[K][L] = BXASE[K][L]*DZASDZ - BZASE[K][L]*DXASDZ + BYASE[K][L]*FAC1;
          BYE[K][L] = BYASE[K][L]*FAC2;
          BZE[K][L] = BZASE[K][L]*DXASDX - BXASE[K][L]*DZASDX + BYASE[K][L]*FAC3;
        }
    } // ------------------------------------

    return;

}




void TS07D_WARPED( double PS, double X, double Y, double Z,
    double BXS[], double BYS[], double BZS[], double BXO[][5], double BYO[][5], double BZO[][5],
    double BXE[][5], double BYE[][5], double BZE[][5], LgmTsyg2007_Info *tInfo ) {

    /*
     *    CALCULATES GSM COMPONENTS OF THE WARPED FIELD FOR TWO TAIL UNIT MODES.
     *    THE WARPING DEFORMATION IS IMPOSED ON THE UNWARPED FIELD, COMPUTED
     *    BY THE S/R "UNWARPED".  THE WARPING PARAMETERS WERE TAKEN FROM THE
     *    RESULTS OF GEOTAIL OBSERVATIONS (TSYGANENKO ET AL. [1998]).
     *    NB # 6, P.106, OCT 12, 2000.
     *
     *    IOPT - TAIL FIELD MODE FLAG:   IOPT=0 - THE TWO TAIL MODES ARE ADDED UP
     *                                   IOPT=1 - MODE 1 ONLY
     *                                   IOPT=2 - MODE 2 ONLY
     */
    int       K, L;
    double    DGDX, XL, DXLDX, SPS, RHO2, RHO, PHI, CPHI, SPHI, XL2, XL4, XL3, RHO4, RR4L4;
    double    F, G, TW, DFDPHI, RR4L4_2, DFDRHO, DFDX, CF, SF, YAS, ZAS, BX_AS1, BY_AS1, BZ_AS1;
    double    BX_AS2, BY_AS2, BZ_AS2, BRHO_AS, BPHI_AS, BRHO_S, BPHI_S;
    double    BX_ASS[6],BX_ASO[6][5],BX_ASE[6][5];
    double    BY_ASS[6],BY_ASO[6][5],BY_ASE[6][5];
    double    BZ_ASS[6],BZ_ASO[6][5],BZ_ASE[6][5];

    G  = tInfo->CB_G.G;
    TW = tInfo->CB_G.TW;

    DGDX  = 0.0;
    XL    = 20.0;
    DXLDX = 0.0;

    SPS  = tInfo->sin_psi_op;
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
    RR4L4_2 = RR4L4*RR4L4;

    F = PHI+G*RHO2*RR4L4*CPHI*SPS +TW*(X/10.0);
    DFDPHI=1.0-G*RHO2*RR4L4*SPHI*SPS;
    DFDRHO=G* RR4L4_2 *(3.0*XL4-RHO4)*CPHI*SPS;
    DFDX=RR4L4*CPHI*SPS*(DGDX*RHO2-G*RHO*RR4L4*4.0*XL3*DXLDX)+TW/10.0;  //  THE LAST TERM DESCRIBES THE IMF-INDUCED TWISTING (ADDED 04/21/06)
//printf("F, DFDPHI, DFDRHO, DFDX = %g %g %g %g\n", F, DFDPHI, DFDRHO, DFDX);


    //CF  = cos(F); SF  = sin(F);
    mysincos( F, &SF, &CF );
    YAS = RHO*CF;
    ZAS = RHO*SF;
//printf("YAS, ZAS = %g %g\n", YAS, ZAS );


    TS07D_UNWARPED( X, YAS, ZAS, BX_ASS, BY_ASS, BZ_ASS, BX_ASO, BY_ASO, BZ_ASO, BX_ASE, BY_ASE, BZ_ASE, tInfo );
//printf("BX_ASS[2], BY_ASS[2], BZ_ASS[2] = %g %g %g\n", BX_ASS[2], BY_ASS[2], BZ_ASS[2] );



    for ( K=1; K<=5; K++ ){

        /*
         * ------------------------------------------- Deforming symmetric modules
         */
        BRHO_AS =  BY_ASS[K]*CF + BZ_ASS[K]*SF;
        BPHI_AS = -BY_ASS[K]*SF + BZ_ASS[K]*CF;

        BRHO_S  = BRHO_AS*DFDPHI;
        BPHI_S  = BPHI_AS - RHO*( BX_ASS[K]*DFDX + BRHO_AS*DFDRHO );

        BXS[K]  = BX_ASS[K]*DFDPHI;
        BYS[K]  = BRHO_S*CPHI - BPHI_S*SPHI;
        BZS[K]  = BRHO_S*SPHI + BPHI_S*CPHI;

    }


    for ( K=1; K<=5; K++ ){
        for ( L=1; L<=4; L++ ){

            /*
             * -------------------------------------------- Deforming odd modules
             */
            BRHO_AS =  BY_ASO[K][L]*CF + BZ_ASO[K][L]*SF;
            BPHI_AS = -BY_ASO[K][L]*SF + BZ_ASO[K][L]*CF;

            BRHO_S  = BRHO_AS*DFDPHI;
            BPHI_S  = BPHI_AS - RHO*( BX_ASO[K][L]*DFDX + BRHO_AS*DFDRHO );

            BXO[K][L] = BX_ASO[K][L]*DFDPHI;
            BYO[K][L] = BRHO_S*CPHI - BPHI_S*SPHI;
            BZO[K][L] = BRHO_S*SPHI + BPHI_S*CPHI;


            /*
             * ------------------------------------------- Deforming even modules
             */
            BRHO_AS =  BY_ASE[K][L]*CF + BZ_ASE[K][L]*SF;
            BPHI_AS = -BY_ASE[K][L]*SF + BZ_ASE[K][L]*CF;

            BRHO_S  = BRHO_AS*DFDPHI;
            BPHI_S  = BPHI_AS - RHO*( BX_ASE[K][L]*DFDX + BRHO_AS*DFDRHO );

            BXE[K][L] = BX_ASE[K][L]*DFDPHI;
            BYE[K][L] = BRHO_S*CPHI - BPHI_S*SPHI;
            BZE[K][L] = BRHO_S*SPHI + BPHI_S*CPHI;

        }
    }

    return;

}




void    TS07D_UNWARPED( double X, double Y, double Z, double BXS[6], double BYS[6], double BZS[6],
        double BXO[6][5], double BYO[6][5], double BZO[6][5], double BXE[6][5], double BYE[6][5],
        double BZE[6][5], LgmTsyg2007_Info *tInfo ) {

    int     K, L;
    double  **TSS, ***TSO, ***TSE;
    double  D0, BXSK, BYSK, BZSK, BXOKL, BYOKL, BZOKL, BXEKL, BYEKL, BZEKL;
    double  HXSK, HYSK, HZSK, HXOKL, HYOKL, HZOKL, HXEKL, HYEKL, HZEKL;
    int     XChanged, YChanged, ZChanged;

    /*
     *    CALCULATES GSM COMPONENTS OF THE SHIELDED FIELD OF 45 TAIL MODES WITH UNIT
     *    AMPLITUDES,  WITHOUT ANY WARPING OR BENDING.  NONLINEAR PARAMETERS OF THE MODES
     *    ARE FORWARDED HERE VIA A COMMON BLOCK /TAIL/.
     *
     */
    TSS = tInfo->TSS;
    TSO = tInfo->TSO;
    TSE = tInfo->TSE;
    D0  = tInfo->CB_TAIL.D;



    // MGH
    tInfo->phi   = atan2( Y, X );
    tInfo->RHO   = sqrt( X*X + Y*Y );
    tInfo->RHOI  = (  tInfo->RHO < 1e-8 ) ? 1e8 : 1.0/tInfo->RHO;
    tInfo->RHOI2 = tInfo->RHOI*tInfo->RHOI;
    tInfo->DPDX  = -Y*tInfo->RHOI2;
    tInfo->DPDY  =  X*tInfo->RHOI2;
    SinN_CosN_Arr( 14, tInfo->phi, tInfo->SMP, tInfo->CMP );
    tInfo->ZD    = sqrt( Z*Z + tInfo->CB_TAIL.D*tInfo->CB_TAIL.D );



    /*
     *
     * --- New tail structure -------------
     *
     *
     */
    for ( K=1; K<=5; K++ ){

        TS07D_TAILSHT_S( K, X, Y, Z, &BXSK, &BYSK, &BZSK, tInfo );
        TS07D_SHTBNORM_S( K, X, Y, Z, &HXSK, &HYSK, &HZSK, tInfo );
        //printf("C: HXSK,HYSK,HZSK = %lf %lf %lf\n", HXSK,HYSK,HZSK);

        BXS[K] = BXSK + HXSK;
        BYS[K] = BYSK + HYSK;
        BZS[K] = BZSK + HZSK;

    }


    for ( K=1; K<=5; K++ ){
        for ( L=1; L<=4; L++ ){

            TS07D_TAILSHT_OE( 1, K, L, X, Y, Z, &BXOKL, &BYOKL, &BZOKL, tInfo );
            TS07D_SHTBNORM_O( K, L, X, Y, Z, &HXOKL, &HYOKL, &HZOKL, tInfo );


            BXO[K][L] = BXOKL + HXOKL;
            BYO[K][L] = BYOKL + HYOKL;
            BZO[K][L] = BZOKL + HZOKL;

            TS07D_TAILSHT_OE( 0, K, L, X, Y, Z, &BXEKL, &BYEKL, &BZEKL, tInfo );
            TS07D_SHTBNORM_E( K, L, X, Y, Z, &HXEKL, &HYEKL, &HZEKL, tInfo );

            BXE[K][L] = BXEKL + HXEKL;
            BYE[K][L] = BYEKL + HYEKL;
            BZE[K][L] = BZEKL + HZEKL;
//            printf("C: BXE[K][L], BYE[K][L], BZE[K][L] = %lf %lf %lf\n", BXE[K][L], BYE[K][L], BZE[K][L]);

        }
    }

    return;

}








void TS07D_TAILSHT_S( int M, double X, double Y, double Z, double *BX, double *BY, double *BZ, LgmTsyg2007_Info *tInfo ) {

    double  D, RNOT, DLTK, RHO, CSPHI, SNPHI, DKM, RKM;
    double  RKMZ, RKMR, ZD, RJ0, RJ1, REX;

    D = tInfo->CB_TAIL.D; // THE COMMON BLOCKS FORWARDS TAIL SHEET THICKNESS
    RNOT = 20.0;          // This can be replaced by introducing them
    DLTK = 1.0;           // through the above common block


    //RHO   = sqrt( X*X + Y*Y );
    RHO   = tInfo->RHO;
    CSPHI = X/RHO;
    SNPHI = Y/RHO;

    DKM = 1.0 + (double)(M-1)*DLTK;
    RKM = DKM/RNOT;

    RKMZ = RKM*Z;
    RKMR = RKM*RHO;

    //ZD = sqrt( Z*Z + D*D );
    ZD = tInfo->ZD;

// MGH
// for a given M,X,Y,Z, these dont change.
// can these be cached?
    RJ0 = bessj0( RKMR );
    RJ1 = bessj1( RKMR );
    REX = exp( RKM*ZD );


    *BX = RKMZ*RJ1*CSPHI/ZD/REX;
    *BY = RKMZ*RJ1*SNPHI/ZD/REX;
    *BZ = RKM*RJ0/REX;

    /*
     * CALCULATION OF THE MAGNETOTAIL CURRENT CONTRIBUTION IS FINISHED
     */

    return;

}











/*
 * MGH - these routines are very wasteful.
 *       AJM  - is the bessel function of order m, J_m( AKNR )
 *       AJMD - looks to be its derivative, J'_m( AKNR ) which is given by J'_m( AKNR ) = J_m-1( AKNR ) - m*J_m( AKNR )/AKNR
 *       There are at least two problems here;
 *          1) Many calcs are repeated over and over.
 *          2) AKN only depends on the AK[n] parameters -- which never change (for a gievn set of parameters)!
 *       Two fixes: 1) compute all at once, and 2) cache results and only recompute if necessary
 *       
 *       Similar problem for other costly functions.
 *       cos( m*phi ) and sin( m*phi ) 
 *       Fix: Switch to recurrence relation and only compite if X or Y changes.
 */
void TS07D_SHTBNORM_S( int K, double X, double Y, double Z, double *FX, double *FY, double *FZ, LgmTsyg2007_Info *tInfo ) {

    int     k, m, mm, n;
    double  AK[6], AKN, AKNR, CHZ, SHZ;
    double  AKNRI, AJM[15], AJM1, AJMD[15];
    double  HX1, HX2, HX, HY1, HY2, HY, HZ;
    double  **TSS;

    TSS = tInfo->TSS;

    //AK[1] = TSS[76][K];
    //AK[2] = TSS[77][K];
    //AK[3] = TSS[78][K];
    //AK[4] = TSS[79][K];
    //AK[5] = TSS[80][K];


    /* fed in through tInfo
    phi   = atan2( Y, X );
    RHO   = sqrt( X*X + Y*Y );
    RHOI  = ( RHO < 1e-8 ) ? 1.0e8 : 1.0/RHO;
    RHOI2 = RHOI*RHOI;
    DPDX  = -Y*RHOI2;
    DPDY  =  X*RHOI2;

    SinN_CosN_Arr( 14, phi, SMP, CMP );
    */

    if ( tInfo->RHO != tInfo->SHTBNORM_S_RHO_LAST ) {
        // do all K's here to simplfy book-keeping
        for ( k=1; k<=5; k++ ) {
            for ( n=1; n<=5; n++ ) {
                AKN   = fabs( TSS[75+n][k] );
                AKNR  = AKN*tInfo->RHO;
                AKNRI = ( AKNR < 1e-8 ) ? 1.0e8 : 1.0/AKNR;
                // cache results
                tInfo->SHTBNORM_S_AKN[k][n]  = AKN;
                tInfo->SHTBNORM_S_AKNR[k][n] = AKNR;

                J_N_Arr( 14, AKNR, AJM, tInfo );
                // cache results
                for ( m=0; m<=14; m++ ) tInfo->SHTBNORM_S_AJM[k][n][m] = AJM[m];
                tInfo->SHTBNORM_S_AJMD[k][n][0] = -tInfo->SHTBNORM_S_AJM[k][n][1]; 
                for ( m=1; m<=14; m++ ) {   
                    tInfo->SHTBNORM_S_AJMD[k][n][m] = tInfo->SHTBNORM_S_AJM[k][n][m-1] 
                                                    - m*tInfo->SHTBNORM_S_AJM[k][n][m]*AKNRI;
                }
            }
        }
        // flag that we've cached them for this RHO value.
        tInfo->SHTBNORM_S_RHO_LAST = tInfo->RHO;
    }





    *FX = *FY = *FZ = 0.0;


    for ( n=1; n<=5; n++ ) {

        //AKN   = fabs( AK[n] );
        //AKNR  = AKN*tInfo->RHO;
        //AKNRI = ( AKNR < 1e-8 ) ? 1.0e8 : 1.0/AKNR;
        AKN = tInfo->SHTBNORM_S_AKN[K][n];

        CHZ = cosh( Z*AKN ); // only changes if Z changes
        SHZ = sinh( Z*AKN ); // only changes if Z changes

        //J_N_Arr( 14, AKNR, AJM, tInfo );
        // using gsl bessel funcs
        //gsl_sf_bessel_Jn_array( 0, 14, AKNR, AJM );

        // using original bessel funcs
        //AJM[0] = bessj0( AKNR );
        //AJM[1] = bessj1( AKNR );
        //for ( mm=2; mm<=14; mm++ ) AJM[mm] = bessj( mm, AKNR );

        //AJMD[0] = -AJM[1]; for ( mm=1; mm<=14; mm++ ) AJMD[mm] = AJM[mm-1] - mm*AJM[mm]*AKNRI;
        

        for ( m=0; m<=14; m++ ) {

            HX1 =  m*tInfo->DPDX*tInfo->SMP[m]*SHZ * tInfo->SHTBNORM_S_AJM[K][n][m];
            HX2 = -AKN*X*tInfo->RHOI*tInfo->CMP[m]*SHZ * tInfo->SHTBNORM_S_AJMD[K][n][m];

            HX = HX1 + HX2;

            HY1 =  m*tInfo->DPDY*tInfo->SMP[m]*SHZ * tInfo->SHTBNORM_S_AJM[K][n][m];
            HY2 = -AKN*Y*tInfo->RHOI*tInfo->CMP[m]*SHZ * tInfo->SHTBNORM_S_AJMD[K][n][m];

            HY = HY1 + HY2;

            HZ = -AKN*tInfo->CMP[m]*CHZ * tInfo->SHTBNORM_S_AJM[K][n][m];


            *FX += HX*TSS[n+5*m][K];
            *FY += HY*TSS[n+5*m][K];
            *FZ += HZ*TSS[n+5*m][K];

        }
    }

    return;

}



void TS07D_TAILSHT_OE( int IEVO, int MK, int M, double X, double Y, double Z, double *BX, double *BY, double *BZ, LgmTsyg2007_Info *tInfo ) {

    double  D0, RNOT, DLTK, RHO, CSPHI, SNPHI, phi, CSMPHI, SNMPHI, DKM, RKM, RKMZ, RKMR;
    double  ZD, REX, AJM, AJMD, AJM1, BRO, BPHI ;


    D0 = tInfo->CB_TAIL.D;      // THE COMMON BLOCKS FORWARDS TAIL SHEET THICKNESS
    RNOT = 20.0;                // Rho_0 - scale parameter along the tail axis
    DLTK = 1.0;                 // step in Km


    //RHO = sqrt( X*X + Y*Y );
    RHO   = tInfo->RHO;
    CSPHI = X/RHO;
    SNPHI = Y/RHO;

    //phi    = atan2( Y, X );
    //CSMPHI = cos( M*phi );
    //SNMPHI = sin( M*phi );
    CSMPHI = tInfo->CMP[M];
    SNMPHI = tInfo->SMP[M];

    DKM = 1.0 + (MK-1)*DLTK;
    RKM = DKM/RNOT;

    RKMZ = RKM*Z;
    RKMR = RKM*RHO;

    //ZD = sqrt( Z*Z + D0*D0 );
    ZD = tInfo->ZD;

// MGH
// for a given M,X,Y,Z, these dont change.
// can these be cached?
    REX = exp( RKM*ZD );


    /*
     * ---- calculating Jm and its derivatives ------
     */
    if ( M > 2 ) {
        AJM  = bessj( M, RKMR );
        AJM1 = bessj( M-1, RKMR );
        AJMD = AJM1 - M*AJM/RKMR;
    } else {
        if ( M == 2 ) {
            AJM  = bessj( 2, RKMR );
            AJM1 = bessj1( RKMR );
            AJMD = AJM1 - M*AJM/RKMR;
        } else {
            AJM  = bessj1( RKMR );
            AJM1 = bessj0( RKMR );
            AJMD = AJM1 - AJM/RKMR;
        }
    }


    if ( IEVO == 0 ) {

        /*-------------------------------------------
            calculating symmetric modes
         --------------------------------------------*/
        BRO  = -M*SNMPHI*Z*AJMD/ZD/REX;
        BPHI = -M*M*CSMPHI*Z*AJM/RKMR/ZD/REX;
        *BZ   = M*SNMPHI*AJM/REX;

    } else {
        /*-------------------------------------------
            calculating asymmetric modes
         --------------------------------------------*/
        BRO  = M*CSMPHI*Z*AJMD/ZD/REX;
        BPHI = -M*M*SNMPHI*Z*AJM/RKMR/ZD/REX;
        *BZ   = -M*CSMPHI*AJM/REX;
    }


    /*
     * --- transformation from cylindrical ccordinates to GSM ---
     */
    *BX = BRO*CSPHI - BPHI*SNPHI;
    *BY = BRO*SNPHI + BPHI*CSPHI;


    //  CALCULATION OF THE MAGNETOTAIL CURRENT CONTRIBUTION IS FINISHED
    return;

}




void TS07D_SHTBNORM_O( int K, int L, double X, double Y, double Z, double *FX, double *FY, double *FZ, LgmTsyg2007_Info *tInfo ) {

    int     m, mm, n, k, l;
    double  AK[6], AKN;
    double  AKNR, CHZ, SHZ, AKNRI, AJM[15], AJM1, AJMD[15];
    double  HX1, HX2,HY1, HY2;
    double  HX, HY, HZ;
    double  ***TSO;

    TSO = tInfo->TSO;

    /*
    AK[1] = TSO[76][K][L];
    AK[2] = TSO[77][K][L];
    AK[3] = TSO[78][K][L];
    AK[4] = TSO[79][K][L];
    AK[5] = TSO[80][K][L];
    */



    /* These values all fed in throuygh tInfo
    phi   = atan2( Y, X );
    RHO   = sqrt( X*X + Y*Y );
    RHOI  = (  RHO < 1e-8 ) ? 1e8 : 1.0/RHO;
    RHOI2 = RHOI*RHOI;
    DPDX  = -Y*RHOI2;
    DPDY  =  X*RHOI2;
    SinN_CosN_Arr( 14, phi, SMP, CMP );
    */

    if ( tInfo->RHO != tInfo->SHTBNORM_O_RHO_LAST ) {
        // do all K's and L's here to simplfy book-keeping
        for ( k=1; k<=5; k++ ) {
            for ( l=1; l<=4; l++ ) {
                for ( n=1; n<=5; n++ ) {
                    AKN   = fabs( TSO[75+n][k][l] );
                    AKNR  = AKN*tInfo->RHO;
                    AKNRI = ( AKNR < 1e-8 ) ? 1.0e8 : 1.0/AKNR;
                    // cache results
                    tInfo->SHTBNORM_O_AKN[k][l][n]  = AKN;
                    tInfo->SHTBNORM_O_AKNR[k][l][n] = AKNR;

                    J_N_Arr( 14, AKNR, AJM, tInfo );
                    // cache results
                    for ( m=0; m<=14; m++ ) tInfo->SHTBNORM_O_AJM[k][l][n][m] = AJM[m];
                    tInfo->SHTBNORM_O_AJMD[k][l][n][0] = -tInfo->SHTBNORM_O_AJM[k][l][n][1]; 
                    for ( m=1; m<=14; m++ ) {   
                        tInfo->SHTBNORM_O_AJMD[k][l][n][m] = tInfo->SHTBNORM_O_AJM[k][l][n][m-1] 
                                                        - m*tInfo->SHTBNORM_O_AJM[k][l][n][m]*AKNRI;
                    }
                }
            }
        }
        // flag that we've cached them for this RHO value.
        tInfo->SHTBNORM_O_RHO_LAST = tInfo->RHO;
    }
    

    *FX = *FY = *FZ = 0.0;


    for ( n=1; n<=5; n++ ){

        //AKN   = fabs( AK[n] );
        //AKNR  = AKN*tInfo->RHO;
        //AKNRI = ( AKNR < 1e-8 ) ? 1e8 : 1.0/AKNR;
        AKN = tInfo->SHTBNORM_O_AKN[K][L][n];

        CHZ = cosh( Z*AKN );
        SHZ = sinh( Z*AKN );

        //J_N_Arr( 14, AKNR, AJM, tInfo );

        // using gsl bessel funcs
        //gsl_sf_bessel_Jn_array( 0, 14, AKNR, AJM );

        // using original bessel funcs
        //AJM[0] = bessj0( AKNR );
        //AJM[1] = bessj1( AKNR );
        //for ( mm=2; mm<=14; mm++ ) AJM[mm] = bessj( mm, AKNR );

        //AJMD[0] = -AJM[1]; for ( mm=1; mm<=14; mm++ ) AJMD[mm] = AJM[mm-1] - mm*AJM[mm]*AKNRI;


        for ( m=0; m<=14; m++ ){

            HX1 =  m*tInfo->DPDX*tInfo->SMP[m]*SHZ*tInfo->SHTBNORM_O_AJM[K][L][n][m];
            HX2 =  -AKN*X*tInfo->RHOI*tInfo->CMP[m]*SHZ*tInfo->SHTBNORM_O_AJMD[K][L][n][m];
            HX  = HX1+HX2;

            HY1 =  m*tInfo->DPDY*tInfo->SMP[m]*SHZ*tInfo->SHTBNORM_O_AJM[K][L][n][m];
            HY2 = -AKN*Y*tInfo->RHOI*tInfo->CMP[m]*SHZ*tInfo->SHTBNORM_O_AJMD[K][L][n][m];
            HY  = HY1 + HY2;

            HZ = -AKN*tInfo->CMP[m]*CHZ*tInfo->SHTBNORM_O_AJM[K][L][n][m];


            *FX += HX*TSO[n+5*m][K][L];
            *FY += HY*TSO[n+5*m][K][L];
            *FZ += HZ*TSO[n+5*m][K][L];

        }
    }

    return;
}





void TS07D_SHTBNORM_E( int K, int L, double X, double Y, double Z, double *FX, double *FY, double *FZ, LgmTsyg2007_Info *tInfo ) {

    int     m, mm, n, k, l;
    double  AK[6], AKN;
    double  AKNR, CHZ, SHZ, AKNRI, RHOI, RHOI2, AJM[15], AJM1, JMD;
    double  DPDX, DPDY, HX1, HX2, HY1, HY2, AJMD[15], HX, HY, HZ;
    double  ***TSE;

    TSE = tInfo->TSE;

    /*
    AK[1] = TSE[76][K][L];
    AK[2] = TSE[77][K][L];
    AK[3] = TSE[78][K][L];
    AK[4] = TSE[79][K][L];
    AK[5] = TSE[80][K][L];
    */


    /* fed in from tInfo
    phi   = atan2( Y, X );
    RHO   = sqrt( X*X + Y*Y );
    RHOI  = (  RHO < 1e-8 ) ? 1e8 : 1.0/RHO;
    RHOI2 = RHOI*RHOI;
    DPDX  = -Y*RHOI2;
    DPDY  = X*RHOI2;

    SinN_CosN_Arr( 14, phi, SMP, CMP );
    */

    if ( tInfo->RHO != tInfo->SHTBNORM_E_RHO_LAST ) {
        // do all K's and L's here to simplfy book-keeping
        for ( k=1; k<=5; k++ ) {
            for ( l=1; l<=4; l++ ) {
                for ( n=1; n<=5; n++ ) {
                    AKN   = fabs( TSE[75+n][k][l] );
                    AKNR  = AKN*tInfo->RHO;
                    AKNRI = ( AKNR < 1e-8 ) ? 1.0e8 : 1.0/AKNR;
                    // cache results
                    tInfo->SHTBNORM_E_AKN[k][l][n]  = AKN;
                    tInfo->SHTBNORM_E_AKNR[k][l][n] = AKNR;

                    J_N_Arr( 14, AKNR, AJM, tInfo );
                    // cache results
                    for ( m=0; m<=14; m++ ) tInfo->SHTBNORM_E_AJM[k][l][n][m] = AJM[m];
                    tInfo->SHTBNORM_E_AJMD[k][l][n][0] = -tInfo->SHTBNORM_E_AJM[k][l][n][1]; 
                    for ( m=1; m<=14; m++ ) {   
                        tInfo->SHTBNORM_E_AJMD[k][l][n][m] = tInfo->SHTBNORM_E_AJM[k][l][n][m-1] 
                                                        - m*tInfo->SHTBNORM_E_AJM[k][l][n][m]*AKNRI;
                    }
                }
            }
        }
        // flag that we've cached them for this RHO value.
        tInfo->SHTBNORM_E_RHO_LAST = tInfo->RHO;
    }



    *FX = *FY = *FZ = 0.0;



    for ( n=1; n<=5; n++ ) {

        //AKN   = fabs( AK[n] );
        //AKNR  = AKN*tInfo->RHO;
        //AKNRI = ( AKNR < 1e-8 ) ? 1e8 : 1.0/AKNR;
        AKN = tInfo->SHTBNORM_E_AKN[K][L][n];

        CHZ = cosh( Z*AKN );
        SHZ = sinh( Z*AKN );


        //J_N_Arr( 14, AKNR, AJM, tInfo );

        // using gsl bessel funcs
        //gsl_sf_bessel_Jn_array( 0, 14, AKNR, AJM );

        // using original bessel funcs
        //AJM[0] = bessj0( AKNR );
        //AJM[1] = bessj1( AKNR );
        //for ( mm=2; mm<=14; mm++ ) AJM[mm] = bessj( mm, AKNR );

        //AJMD[0] = -AJM[1]; for ( mm=1; mm<=14; mm++ ) AJMD[mm] = AJM[mm-1] - mm*AJM[mm]*AKNRI;

        for ( m=0; m<=14; m++ ) {

            HX1 = -m*tInfo->DPDX*tInfo->CMP[m]*SHZ*tInfo->SHTBNORM_E_AJM[K][L][n][m];
            HX2 = -AKN*X*tInfo->RHOI*tInfo->SMP[m]*SHZ*tInfo->SHTBNORM_E_AJMD[K][L][n][m];
            HX  = HX1 + HX2;

            HY1 = -m*tInfo->DPDY*tInfo->CMP[m]*SHZ*tInfo->SHTBNORM_E_AJM[K][L][n][m];
            HY2 = -AKN*Y*tInfo->RHOI*tInfo->SMP[m]*SHZ*tInfo->SHTBNORM_E_AJMD[K][L][n][m];
            HY  = HY1 + HY2;

            HZ = -AKN*tInfo->SMP[m]*CHZ*tInfo->SHTBNORM_E_AJM[K][L][n][m];


            *FX += HX*TSE[n+5*m][K][L];
            *FY += HY*TSE[n+5*m][K][L];
            *FZ += HZ*TSE[n+5*m][K][L];

        }
    }

    return;

}


double  bessj0( double x ) {

    /*
     * (C) Copr. 1986-92 Numerical Recipes Software .)+1YX39'=K+1.
     */

    static double  p1 =  1.0;
    static double  p2 = -0.1098628627e-2;
    static double  p3 =  0.2734510407e-4;
    static double  p4 = -0.2073370639e-5;
    static double  p5 =  0.2093887211e-6;

    static double  q1 = -0.1562499995e-1;
    static double  q2 =  0.1430488765e-3;
    static double  q3 = -0.6911147651e-5;
    static double  q4 =  0.7621095161e-6;
    static double  q5 = -0.934945152e-7;

    static double  r1 =  57568490574.;
    static double  r2 = -13362590354.;
    static double  r3 =  651619640.7;
    static double  r4 = -11214424.18;
    static double  r5 =  77392.33017;
    static double  r6 = -184.9052456;

    static double  s1 = 57568490411.;
    static double  s2 = 1029532985.;
    static double  s3 = 9494680.718;
    static double  s4 = 59272.64853;
    static double  s5 = 267.8532712;
    static double  s6 = 1.0;

    double  y, result, ax, z, xx;

    if ( fabs(x) < 8.0 ) {
        y = x*x;
        result = (r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))/(s1+y*(s2+y*(s3+y*(s4+y*(s5+y*s6)))));
    } else {
        ax = fabs(x);
        z  = 8.0/ax;
        y  = z*z;
        xx = ax - 0.785398164;
        result = sqrt( 0.636619772/ax )*(cos(xx)*(p1+y*(p2+y*(p3+y*(p4+y*p5))))-z*sin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))));
    }

    return( result );

}




double  bessj1( double x ) {

    /*
     * (C) Copr. 1986-92 Numerical Recipes Software .)+1YX39'=K+1.
     */


    static double  r1 =  72362614232.;
    static double  r2 = -7895059235.;
    static double  r3 =  242396853.1;
    static double  r4 = -2972611.439;
    static double  r5 =  15704.48260;
    static double  r6 = -30.16036606;

    static double  s1 = 144725228442.;
    static double  s2 = 2300535178.;
    static double  s3 = 18583304.74;
    static double  s4 = 99447.43394;
    static double  s5 = 376.9991397;
    static double  s6 = 1.0;

    static double  p1 =  1.0;
    static double  p2 =  0.183105e-2;
    static double  p3 = -0.3516396496e-4;
    static double  p4 =  0.2457520174e-5;
    static double  p5 = -0.240337019e-6;

    static double  q1 =  0.04687499995;
    static double  q2 = -0.2002690873e-3;
    static double  q3 =  0.8449199096e-5;
    static double  q4 = -0.88228987e-6;
    static double  q5 =  0.105787412e-6;

    double  y, result, ax, z, xx;

    if ( fabs(x) < 8.0 ) {
        y = x*x;
        result = x*(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))/(s1+y*(s2+y*(s3+y*(s4+y*(s5+y*s6)))));
    } else {
        ax = fabs(x);
        z  = 8.0/ax;
        y  = z*z;
        xx = ax - 2.356194491;
        if ( x >= 0.0 ) {
            result =  sqrt(0.636619772/ax)*(cos(xx)*(p1+y*(p2+y*(p3+y*(p4+y*p5))))-z*sin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))));
        } else {
            result = -sqrt(0.636619772/ax)*(cos(xx)*(p1+y*(p2+y*(p3+y*(p4+y*p5))))-z*sin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))));
        }
    }

    return( result );

}




double  bessj( int n, double x ) {

    /*
     * (C) Copr. 1986-92 Numerical Recipes Software .)+1YX39'=K+1.
     */
    int     IACC=40, j;
    double  BIGNO=1e10;
    double  BIGNI=1e-10;
    int     m, jsum;
    double  ax, result, tox, bjm, bjp, bj, sum;


    if ( n < 2 ) {
        printf("bessj(): bad argument n in bessj\n");
        exit(-1);
    }

    ax = fabs(x);

    if( ax == 0.0 ) {

        result = 0.0;

    } else if ( ax > (double)n ) {

        tox = 2.0/ax;
        bjm = bessj0( ax );
        bj  = bessj1( ax );

        for ( j=1; j<n; j++ ) {
            bjp = j*tox*bj - bjm;
            bjm = bj;
            bj  = bjp;
        }

        result = bj;

    } else {

        tox = 2.0/ax;
        m = 2*( ( n + (int)(sqrt(IACC*n)) )/2 ); // CHECK THIS

        jsum = 0;
        result = sum = 0.0;
        bjp = 0.0;
        bj  = 1.0;

        for ( j=m; j>=1; j-- ) {

            bjm = j*tox*bj - bjp;
            bjp = bj;
            bj  = bjm;

            if ( fabs(bj) > BIGNO ) {
                bj     *= BIGNI;
                bjp    *= BIGNI;
                result *= BIGNI;
                sum    *= BIGNI;
            }

            if ( jsum != 0 ) sum += bj;
            jsum = 1 - jsum;
            if ( j == n ) result = bjp;

        }

        sum   = 2.0*sum - bj;
        result /= sum;
    }

    if( (x < 0.0)&&( n%2 == 1) ) result = -result; // CHECK

    return( result );
}


// mod due to Jay Albert (its the NR code with mods to save intermediate vals)...
void bessjj( int n, double x, double *Jarr ) {

    /*
     * (C) Copr. 1986-92 Numerical Recipes Software .)+1YX39'=K+1.
     */
    int     IACC=40;
    double  BIGNO=1e10;
    double  BIGNI=1e-10;
    int     m, jsum, i, j, k;
    double  ax, result, tox, bjm, bjp, bj, sum, dum;
    int     IEXP = DBL_MAX_EXP/2;

    if ( n < 0 ) {
        printf("bessjj(): bad argument n in bessjj: n = %d\n", n);
        exit(-1);
    }

    if ( n == 0 ) {
        Jarr[0] = bessj0( x );
        return;
    }

    if ( n == 1 ) {
        Jarr[0] = bessj0( x );
        Jarr[1] = bessj1( x );
        return;
    }


    ax = fabs(x);

    if( ax == 0.0 ) {

        Jarr[0] = 0.0;

    } else if ( ax > (double)n ) {

        tox = 2.0/ax;
        Jarr[0] = bjm = bessj0( ax );   // J_0(x)
        Jarr[1] = bj  = bessj1( ax );   // J_1(x)

        for ( j=1; j<n; j++ ) {
            bjp = j*tox*bj - bjm;
            bjm = bj;
            bj  = bjp;
            Jarr[j+1] = bjp;            // J_j+1(x)
        }


    } else {

        tox = 2.0/ax;
        m = 2*( ( n + (int)(sqrt(IACC*n)) )/2 );

        jsum = 0;
        sum = 0.0;
        bjp = 0.0;
        bj  = 1.0;
        for ( i=1; i<=n; i++ ) Jarr[i] = 0.0;

        for ( j=m; j>=1; j-- ) {

            bjm = j*tox*bj - bjp;
            bjp = bj;
            bj  = bjm;

            if ( fabs(bj) > BIGNO ) {
                bj     *= BIGNI;
                bjp    *= BIGNI;
                sum    *= BIGNI;
                for ( i=j+1; i<=n; i++ ) *(Jarr+i) *= BIGNI;
            }
            //dum = frexp( bj, &k );
            //if ( k > IEXP ) {
            //    bj  = ldexp( bj,  -IEXP );
            //    bjp = ldexp( bjp, -IEXP );
            //    sum = ldexp( sum, -IEXP );
            //    for ( i=j+1; i<=n; i++ ) Jarr[i] = ldexp( Jarr[i], -IEXP );
            //}

            if ( jsum != 0 ) sum += bj;
            jsum = 1 - jsum;
            if ( j <= n ) Jarr[j] = bjp; // J_j(x)

        }

        sum   = 2.0*sum - bj;
        for ( i=1; i<=n; i++ ) Jarr[i] /= sum;
        Jarr[0] = bj/sum;
    }

    if( x < 0.0 ) {
        for ( i=1; i<=n; i+=2 ) Jarr[i] = -Jarr[i];
    }

    return;
}

int startlau( int n, double x, int t ) {

    int s;
    double  p, q, r, y;
    
    s = 2*t-1;
    p = 36.0/x-t;
    r=n/x;
    if ( (r>1.0) || (t==1) ){
        q = sqrt(r*r+s);
        r = r*log(q+r)-q;
    } else {
        r = 0.0;
    }
    q = 18.0/x+r;
    r = (p>q) ? p : q;
    p = sqrt(2.0*(t+r));
    p = x*((1.0+r)+p)/(1.0+p);
    y = 0.0;
    q = y;
    do{
        y=p;
        p/=x;
        q=sqrt(p*p+s);
        p = x*(r+q)/log(p +q);
        q = y;
    } while ( (p>q) || (p < (q-1.0)) );
    return( (t==1) ? floor(p+1.0) : -floor(-p/2.0)*2);
    

}


void bessjlau( int n, double x, double *Jarr ) {

    int l, m, nu, signx;
    double x2, r, s;
    

    if ( x == 0.0 ) {

        Jarr[0] = 1.0;
        for (; n>=1; n--) Jarr[n] = 0.0;

    } else {

        signx = ( x> 0.0) ? 1 : -1;
        x = fabs(x);
        r = s = 0.0;
        x2 = 2.0/x;
        l = 0;
        nu = startlau(n, x, 0);
        for (m=nu; m>=1; m--){
            r = 1.0/(x2*m-r);
            l = 2-l;
            s = r*(l+s);
            if ( m<=n) Jarr[m] = r;
        }   
        Jarr[0] = r=1.0/(1.0+s);

        for (m=1; m<=n; m++) r = Jarr[m] *= r;
        if (signx < 0.0) 
            for (m=1; m<=n; m +=2) Jarr[n] = Jarr[m];
        
        

    }

    
        
}



void     TS07D_BIRK_TOT( double PS, double X, double Y, double Z,
        double *BX11, double *BY11, double *BZ11, double *BX12, double *BY12, double *BZ12,
        double *BX21, double *BY21, double *BZ21, double *BX22, double *BY22, double *BZ22, LgmTsyg2007_Info *tInfo ){


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
                        .2256245602,-.05841594319 };

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
                        .1379899178,.06607020029 };

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
                        .1930034238,-.02261109942 };

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
                        .1486276863,.06859991529 };

    double      X_SC, FX11, FY11, FZ11, HX11, HY11, HZ11, FX12, FY12, FZ12, HX12, HY12, HZ12;
    double      FX21, FY21, FZ21, HX21, HY21, HZ21, FX22, FY22, FZ22, HX22, HY22, HZ22, XKAPPA;
    int         J, PSChanged, XChanged, YChanged, ZChanged;

    PSChanged = ( fabs(PS-tInfo->OLD_PS) > 1e-10 ) ? TRUE : FALSE;
    XChanged = ( fabs(X-tInfo->OLD_X) > 1e-10 ) ? TRUE : FALSE;
    YChanged = ( fabs(Y-tInfo->OLD_Y) > 1e-10 ) ? TRUE : FALSE;
    ZChanged = ( fabs(Z-tInfo->OLD_Z) > 1e-10 ) ? TRUE : FALSE;

    tInfo->CB_DPHI_B_RHO0.XKAPPA = tInfo->CB_BIRKPAR.XKAPPA1;               // FORWARDED IN BIRK_1N2
    X_SC                         = tInfo->CB_BIRKPAR.XKAPPA1 - 1.1;         // FORWARDED IN BIRK_SHL

    TS07D_BIRK_1N2( 1, 1, PS, X, Y, Z, &FX11, &FY11, &FZ11, tInfo );        //  REGION 1, MODE 1
// where are J, PSChanged, XChanged, YChanged, ZChanged defined???
    TS07D_BIRK_SHL( 0, PSChanged, XChanged, YChanged, ZChanged, SH11, PS, X_SC, X, Y, Z, &HX11, &HY11, &HZ11, tInfo);
//printf("C 11 >>>>>>>>>>>> FX11,FY11,FZ11 = %g %g %g\n", FX11,FY11,FZ11 );
//printf("C 11 >>>>>>>>>>>> HX11,HY11,HZ11 = %g %g %g\n", HX11,HY11,HZ11 );
    *BX11 = FX11 + HX11;
    *BY11 = FY11 + HY11;
    *BZ11 = FZ11 + HZ11;

    TS07D_BIRK_1N2( 1, 2, PS, X, Y, Z, &FX12, &FY12, &FZ12, tInfo );        //  REGION 1, MODE 2
    TS07D_BIRK_SHL( 1, PSChanged, XChanged, YChanged, ZChanged, SH12, PS, X_SC, X, Y, Z, &HX12, &HY12, &HZ12, tInfo );
//printf("C 12 >>>>>>>>>>>> FX12,FY12,FZ12 = %g %g %g\n", FX12,FY12,FZ12 );
//printf("C 12 >>>>>>>>>>>> HX12,HY12,HZ12 = %g %g %g\n", HX12,HY12,HZ12 );
    *BX12 = FX12 + HX12;
    *BY12 = FY12 + HY12;
    *BZ12 = FZ12 + HZ12;




    tInfo->CB_DPHI_B_RHO0.XKAPPA = tInfo->CB_BIRKPAR.XKAPPA2;               // FORWARDED IN BIRK_1N2
    X_SC                         = tInfo->CB_BIRKPAR.XKAPPA2 - 1.0;         // FORWARDED IN BIRK_SHL

    TS07D_BIRK_1N2( 2, 1, PS, X, Y, Z, &FX21, &FY21, &FZ21, tInfo );        //  REGION 2, MODE 1
    TS07D_BIRK_SHL( 2, PSChanged, XChanged, YChanged, ZChanged, SH21, PS, X_SC, X, Y, Z, &HX21, &HY21, &HZ21, tInfo );
//printf("C 21 >>>>>>>>>>>> FX21,FY21,FZ21 = %g %g %g\n", FX21,FY21,FZ21 );
//printf("C 21 >>>>>>>>>>>> HX21,HY21,HZ21 = %g %g %g\n", HX21,HY21,HZ21 );
    *BX21 = FX21 + HX21;
    *BY21 = FY21 + HY21;
    *BZ21 = FZ21 + HZ21;

    TS07D_BIRK_1N2( 2, 2, PS, X, Y, Z, &FX22, &FY22, &FZ22, tInfo );        //  REGION 2, MODE 2
    TS07D_BIRK_SHL( 3, PSChanged, XChanged, YChanged, ZChanged, SH22, PS, X_SC, X, Y, Z, &HX22, &HY22, &HZ22, tInfo );
//printf("C 22 >>>>>>>>>>>> FX22,FY22,FZ22 = %g %g %g\n", FX22,FY22,FZ22 );
//printf("C 22 >>>>>>>>>>>> HX22,HY22,HZ22 = %g %g %g\n", HX22,HY22,HZ22 );
    *BX22 = FX22 + HX22;
    *BY22 = FY22 + HY22;
    *BZ22 = FZ22 + HZ22;

    return;

}

void    TS07D_BIRK_1N2( int NUMB, int MODE, double PS, double X, double Y, double Z, double *BX, double *BY, double *BZ, LgmTsyg2007_Info *tInfo ) {


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
                   2.492118385,.7113544659 };

    static double A12[] = { -9e99, .7058026940,-.2845938535,5.715471266,-2.472820880,
                   -.7738802408,.3478293930,-11.37653694,-38.64768867,.6932927651,
                   -212.4017288,4.944204937,3.071270411,33.05882281,7.387533799,
                   2.366769108,79.22572682,.6154290178,.5592050551,-.1796585105,
                   -1.654932210,.7309108776,-.4926292779,-1.130266095,-.009613974555,
                   .1484586169,-.2215347198,7.883592948,.02768251655,2.950280953,
                   1.212634762,.5567714182 };

    static double A21[] = { -9e99, .1278764024,-.2320034273,1.805623266,-32.37241440,
                   -.9931490648,.3175085630,-2.492465814,-16.21600096,.2695393416,
                   -6.752691265,3.971794901,14.54477563,41.10158386,7.912889730,
                   1.258297372,9.583547721,1.014141963,.5104134759,-.1790430468,
                   -1.756358428,.7561986717,-.6775248254,-.04014016420,.01446794851,
                   .1200521731,-.2203584559,4.508963850,.8221623576,1.779933730,
                   1.102649543,.8867880020 };

    static double A22[] = { -9e99, .4036015198,-.3302974212,2.827730930,-45.44405830,
                   -1.611103927,.4927112073,-.003258457559,-49.59014949,.3796217108,
                   -233.7884098,4.312666980,18.05051709,28.95320323,11.09948019,
                   .7471649558,67.10246193,.5667096597,.6468519751,-.1560665317,
                   -1.460805289,.7719653528,-.6658988668,.2515179349E-05,
                   .02426021891,.1195003324,-.2625739255,4.377172556,.2421190547,
                   2.503482679,1.071587299,.7247997430 };


    double    BETA=0.9, RH=10.0;    // parameters of the tilt-dependent deformation of the untilted F.A.C. field
//    double    EPS=3.0;            // parameters of the tilt-dependent deformation of the untilted F.A.C. field
    double    RHO2, Xsc, Xsc2, Ysc, Zsc, Zsc2, RHO, RHOSQ;
    double    Rsc, Ysc2, PHI, SPHIC, CPHIC, BRACK, R1RH, R1RH2, R1RH3, PSIAS, PHIS, DPHISPHI;
    double    RHO2pRHOSQ, RHO2pRHOSQ2, DPHISRHO, DPHISDY, SPHICS, CPHICS, XS, ZS, BXS, BYAS;
    double    BZS, BRHOAS, BPHIAS, BRHO_S, BPHI_S, BY_S;


    tInfo->CB_DPHI_B_RHO0.B     = 0.5;
    tInfo->CB_DPHI_B_RHO0.RHO_0 = 7.0; RHO2 = 49.0;

    tInfo->CB_MODENUM.M = MODE;
    if (NUMB == 1) {
        tInfo->CB_DPHI_B_RHO0.DPHI   = 0.055;
        tInfo->CB_DTHETA.DTHETA = 0.06;
    } else if (NUMB == 2) {
        tInfo->CB_DPHI_B_RHO0.DPHI   = 0.030;
        tInfo->CB_DTHETA.DTHETA = 0.09;
    }

    Xsc = X*tInfo->CB_DPHI_B_RHO0.XKAPPA; Xsc2 = Xsc*Xsc;
    Ysc = Y*tInfo->CB_DPHI_B_RHO0.XKAPPA; Ysc2 = Ysc*Ysc;
    Zsc = Z*tInfo->CB_DPHI_B_RHO0.XKAPPA; Zsc2 = Zsc*Zsc;
    RHO = sqrt(Xsc2+Zsc2);
    RHOSQ = Xsc2+Zsc2;

    Rsc  = sqrt(Xsc2+Ysc2+Zsc2);        // SCALED

    if ( (Xsc == 0.0) && (Zsc == 0.0) ) {
        PHI   = 0.0;
        SPHIC = 0.0;
        CPHIC = 1.0;                // "C" means "CYLINDRICAL", TO DISTINGUISH FROM SPHERICAL PHI
    } else {
        PHI   = atan2(-Zsc, Xsc);            // FROM CARTESIAN TO CYLINDRICAL (RHO,PHI,Y)
        //SPHIC = sin(PHI);
        //CPHIC = cos(PHI);                // "C" means "CYLINDRICAL", TO DISTINGUISH FROM SPHERICAL PHI
        SPHIC = -Zsc/RHO;
        CPHIC =  Xsc/RHO;
    }


    BRACK = tInfo->CB_DPHI_B_RHO0.DPHI+tInfo->CB_DPHI_B_RHO0.B*RHO2/(RHO2+1.0)*(RHOSQ-1.0)/(RHO2+RHOSQ);
    R1RH  = (Rsc-1.0)/RH;
    R1RH2 = R1RH*R1RH; R1RH3 = R1RH2*R1RH;
double R1RH3p1 = 1.0+R1RH3;
double  cbrt_R1RH3p1 = cbrt( R1RH3p1 );
    PSIAS = BETA*PS/cbrt_R1RH3p1;

    PHIS = PHI-BRACK*SPHIC - PSIAS;
    DPHISPHI = 1.0-BRACK*CPHIC;
    RHO2pRHOSQ = RHO2+RHOSQ; RHO2pRHOSQ2 = RHO2pRHOSQ*RHO2pRHOSQ;
double  fff = 1.0/(RH*Rsc*R1RH3p1*cbrt_R1RH3p1);
    DPHISRHO = -2.0*tInfo->CB_DPHI_B_RHO0.B*RHO2*RHO/RHO2pRHOSQ2*SPHIC + BETA*PS*R1RH2*RHO*fff;
    DPHISDY= BETA*PS*R1RH2*Ysc*fff;

    //SPHICS = sin(PHIS); CPHICS = cos(PHIS);
    mysincos( PHIS, &SPHICS, &CPHICS );

    XS =  RHO*CPHICS;
    ZS = -RHO*SPHICS;

    if (NUMB == 1) {
        if (MODE == 1) TS07D_TWOCONES( A11, XS, Ysc, ZS, &BXS, &BYAS, &BZS, tInfo );
        if (MODE == 2) TS07D_TWOCONES( A12, XS, Ysc, ZS, &BXS, &BYAS, &BZS, tInfo );
    } else {
        if (MODE == 1) TS07D_TWOCONES( A21, XS, Ysc, ZS, &BXS, &BYAS, &BZS, tInfo );
        if (MODE == 2) TS07D_TWOCONES( A22, XS, Ysc, ZS, &BXS, &BYAS, &BZS, tInfo );
    }

    BRHOAS = BXS*CPHICS-BZS*SPHICS;
    BPHIAS = -BXS*SPHICS-BZS*CPHICS;

    BRHO_S = BRHOAS*DPHISPHI                             * tInfo->CB_DPHI_B_RHO0.XKAPPA;    // SCALING
    BPHI_S = (BPHIAS-RHO*(BYAS*DPHISDY+BRHOAS*DPHISRHO)) * tInfo->CB_DPHI_B_RHO0.XKAPPA;
    BY_S   = BYAS*DPHISPHI                               * tInfo->CB_DPHI_B_RHO0.XKAPPA;

    *BX = BRHO_S*CPHIC-BPHI_S*SPHIC;
    *BY = BY_S;
    *BZ = -BRHO_S*SPHIC-BPHI_S*CPHIC;

    return;

}


void     TS07D_TWOCONES( double *A, double X, double Y, double Z, double *BX, double *BY, double *BZ, LgmTsyg2007_Info *tInfo ){

    /*
     *
     *
     *     ADDS FIELDS FROM TWO CONES (NORTHERN AND SOUTHERN), WITH A PROPER SYMMETRY OF THE CURRENT AND FIELD,
     *       CORRESPONDING TO THE REGION 1 BIRKELAND CURRENTS.
     */
    double    BXN, BYN, BZN, BXS, BYS, BZS;

    TS07D_ONE_CONE( A, X, Y, Z, &BXN, &BYN, &BZN, tInfo );
    TS07D_ONE_CONE( A, X, -Y, -Z, &BXS, &BYS, &BZS, tInfo );
    *BX = BXN - BXS;
    *BY = BYN + BYS;
    *BZ = BZN + BZS;

    return;

}


void    TS07D_ONE_CONE( double *A, double X, double Y, double Z, double *BX, double *BY, double *BZ, LgmTsyg2007_Info *tInfo ) {

    /*
     *
     *
     *
     *  RETURNS FIELD COMPONENTS FOR A DEFORMED CONICAL CURRENT SYSTEM, FITTED TO A BIOSAVART FIELD
     *    BY SIM_14.FOR.  HERE ONLY THE NORTHERN CONE IS TAKEN INTO ACCOUNT.
     *
     */

    double  DR=1e-6, DT=1e-6;    // JUST FOR NUMERICAL DIFFERENTIATION
    double  THETA0, RHO2, RHO, R, THETA, PHI, RS, THETAS, PHIS, BTAST, BFAST;
    double  DRSDR, DRSDT, DTSDR, DTSDT, STSST, RSR, BR, BTHETA, BPHI;
    double  S, C, SF, CF, BE;
    double  SinTheta, Sin2Theta, Sin3Theta, CosTheta, Cos2Theta;
    double  CosACosB, SinASinB, CosApB, CosAmB, Cos2ApB, Cos2AmB;
    double  SinACosB, CosASinB, SinApB, Sin2ApB, Sin3ApB, SinAmB, Sin2AmB, Sin3AmB;
    double  SIN_PHIS, COS_PHIS, SIN_THETAS, COS_THETAS;


    THETA0 = A[31];

// MGH
// Dont we have these in the tInfop structure already?
    RHO2  = X*X+Y*Y;
    RHO   = sqrt(RHO2);
    R     = sqrt(RHO2+Z*Z);
    THETA = atan2(RHO,Z);
    PHI   = atan2(Y,X);


    /*
     *   MAKE THE DEFORMATION OF COORDINATES:
     */
    SinTheta = RHO/R;
    CosTheta = Z/R;

    Sin2Theta = 2.0*SinTheta*CosTheta;                  // Used identity sin(2a) = 2sin(a)cos(a)
    Sin3Theta = SinTheta*(3.0 - 4.0*SinTheta*SinTheta); // Used identity sin(3a) = 3sin(a) - 4sin^3(a)
    Cos2Theta = 2.0*CosTheta*CosTheta - 1.0;            // Used identity cos(2a) = 2cos^2(a) - 1
    RS     = TS07D_R_S( A, R, CosTheta, Cos2Theta );
    THETAS = TS07D_THETA_S( A, R, THETA, SinTheta, Sin2Theta, Sin3Theta );
    mysincos( THETAS, &SIN_THETAS, &COS_THETAS );

    PHIS     = PHI;
    SIN_PHIS = Y/RHO;
    COS_PHIS = X/RHO;

    /*
     * Efficiently compute cos( THETA+DT ), cos( 2(THETA+DT) ), cos( THETA-DT ), cos( 2(THETA-DT) )
     * Use identities cos(a+b) = cos(a)cos(b) - sin(a)sin(b)
     *                cos(a-b) = cos(a)cos(b) + sin(a)sin(b)
     * NOTE: this ASSUMES that b == 1e-6 (Used linux bc to get cos(1e-6) and sin(1e-6) )
     */
    CosACosB = CosTheta*0.99999999999949995555;
    SinASinB = SinTheta*0.00000099999999999983;

    CosApB = CosACosB - SinASinB;
    Cos2ApB = 2.0*CosApB*CosApB - 1.0;

    CosAmB = CosACosB + SinASinB;
    Cos2AmB = 2.0*CosAmB*CosAmB - 1.0;




    SinACosB = SinTheta*0.99999999999949995555;
    CosASinB = CosTheta*0.00000099999999999983;

    SinApB  = SinACosB + CosASinB;
    Sin2ApB = 2.0*SinApB*CosApB;
    Sin3ApB = SinApB*(3.0 - 4.0*SinApB*SinApB);

    SinAmB  = SinACosB - CosASinB;
    Sin2AmB = 2.0*SinAmB*CosAmB;

    Sin3AmB = SinAmB*(3.0 - 4.0*SinAmB*SinAmB);



    /*
     *   CALCULATE FIELD COMPONENTS AT THE NEW POSITION (ASTERISKED):
     */
    TS07D_FIALCOS( RS, THETAS, SIN_THETAS, COS_THETAS, SIN_PHIS, COS_PHIS, &BTAST, &BFAST, tInfo->CB_MODENUM.M, THETA0, tInfo->CB_DTHETA.DTHETA);    // MODE #M
//FIX what gets returned here?


    /*
     *   NOW TRANSFORM B{R,T,F}_AST BY THE DEFORMATION TENSOR:
     *
     *      FIRST OF ALL, FIND THE DERIVATIVES:
     */
    DRSDR = (TS07D_R_S(A, R+DR, CosTheta, Cos2Theta)-TS07D_R_S(A, R-DR, CosTheta, Cos2Theta))/(2.0*DR);
    DRSDT = (TS07D_R_S(A,R, CosApB, Cos2ApB)-TS07D_R_S(A,R, CosAmB, Cos2AmB))/(2.0*DT);

    DTSDR = (TS07D_THETA_S(A,R+DR,THETA, SinTheta, Sin2Theta, Sin3Theta)-TS07D_THETA_S(A,R-DR,THETA, SinTheta, Sin2Theta, Sin3Theta))/(2.0*DR);
    DTSDT = (TS07D_THETA_S(A,R,THETA+DT, SinApB, Sin2ApB, Sin3ApB)-TS07D_THETA_S(A,R,THETA-DT, SinAmB, Sin2AmB, Sin3AmB))/(2.0*DT);

    STSST = SIN_THETAS/SinTheta;
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


double TS07D_R_S( double *A, double R, double CosTheta, double Cos2Theta ) {

    double    RS, R2, A11_2, A12_2, A13_2, A14_2, A15_2, A16_2, R2pA16_2, R2pA16_2_2;

    R2 = R*R;
    A11_2 = A[11]*A[11];
    A12_2 = A[12]*A[12];
    A13_2 = A[13]*A[13];
    A14_2 = A[14]*A[14];
    A15_2 = A[15]*A[15];
    A16_2 = A[16]*A[16];
    R2pA16_2 = R2+A16_2; R2pA16_2_2 = R2pA16_2*R2pA16_2;

    RS = R+A[2]/R +A[3]*R/sqrt(R2+A11_2)+A[4]*R/(R2+A12_2)
           +(A[5]+A[6]/R+A[7]*R/sqrt(R2+A13_2)+A[8]*R/(R2+A14_2))*CosTheta
           +(A[9]*R/sqrt(R2+A15_2)+A[10]*R/R2pA16_2_2)*Cos2Theta;

      return( RS );

}
double TS07D_THETA_S( double *A, double R, double Theta, double SinTheta, double Sin2Theta, double Sin3Theta) {

    double    TS, R2, A27_2, A28_2, A29_2, A30_2;

    R2 = R*R;
    A27_2 = A[27]*A[27];
    A28_2 = A[28]*A[28];
    A29_2 = A[29]*A[29];
    A30_2 = A[30]*A[30];

    TS = Theta +(A[17]+A[18]/R+A[19]/R2
            +A[20]*R/sqrt(R2+A27_2))*SinTheta
            +(A[21]+A[22]*R/sqrt(R2+A28_2)
            +A[23]*R/(R2+A29_2))*Sin2Theta
            +(A[24]+A[25]/R+A[26]*R/(R2+A30_2))*Sin3Theta;

      return( TS );
}



void     TS07D_FIALCOS( double R, double THETA, double SIN_THETA, double COS_THETA, double SIN_PHI, double COS_PHI, double *BTHETA, double *BPHI, int N, double THETA0, double DT) {

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

    int        M;
    double    SINTE, RO, COSTE, SINFI, COSFI, TG, CTG, TETANP, TETANM, TGP=0.0, TGM=0.0, TGM2=0.0, TGP2=0.0;
    double    COSM1, SINM1, TM, TGM2M, TGP2M, T, FC, FC1, TGM2M1, TG21, DTT, DTT0;
    double    BTN[10], BPN[10], CCOS[10], SSIN[10];

    //SINTE = sin(THETA); COSTE = cos(THETA);
    //mysincos( THETA, &SINTE, &COSTE );
    SINTE = SIN_THETA; COSTE = COS_THETA;

    RO    = R*SINTE;

    //SINFI = sin(PHI); COSFI = cos(PHI);
    //mysincos( PHI, &SINFI, &COSFI );
    SINFI = SIN_PHI; COSFI = COS_PHI;

    TG    = SINTE/(1.0+COSTE);        //  TAN(THETA/2)
    CTG   = SINTE/(1.0-COSTE);        //  CTG(THETA/2)


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


// This is one of the costliest routines -- mostly due to sin,cos,exp
// Almost 13% of Lstar calc is done in here.
void    TS07D_BIRK_SHL( int J, int PSChanged, int XChanged, int YChanged, int ZChanged, double *A, double PS, double X_SC,
                                double X, double Y, double Z, double *BX, double *BY, double *BZ, LgmTsyg2007_Info *tInfo ) {

    int        L, M, I, N, NN, K;
    double    GX, GY, GZ, FX, FY, FZ, HX, HY, HZ, HXR, HZR;



    /*
     *  There are 4 different sets of coefficients that come to us. We added an argument to the
     *  calling routine so that we know which set we are dealing with. Then many things dont need to
     *  be recomputed. Probably dont need to cache all quantities here -- i.e. could prob optimize memory
     *  usage as well...
     */
    for ( I=1;  I<=3; I++ ) {
        if ( !tInfo->DoneJ[J] ) {

            tInfo->P[J][I] = A[72+I]; tInfo->ooP[J][I] = 1.0/tInfo->P[J][I]; tInfo->ooP2[J][I] = tInfo->ooP[J][I]*tInfo->ooP[J][I];
            tInfo->Q[J][I] = A[78+I]; tInfo->ooQ[J][I] = 1.0/tInfo->Q[J][I]; tInfo->ooQ2[J][I] = tInfo->ooQ[J][I]*tInfo->ooQ[J][I];

            tInfo->R[J][I] = A[75+I]; tInfo->ooR[J][I] = 1.0/tInfo->R[J][I]; tInfo->ooR2[J][I] = tInfo->ooR[J][I]*tInfo->ooR[J][I];
            tInfo->S[J][I] = A[81+I]; tInfo->ooS[J][I] = 1.0/tInfo->S[J][I]; tInfo->ooS2[J][I] = tInfo->ooS[J][I]*tInfo->ooS[J][I];
        }
    }


    /*
     *  if PS has changed we need to recompute things
     */
    if ( PSChanged || !tInfo->DoneJ[J] ) {
        tInfo->CPS  = tInfo->cos_psi_op;
        tInfo->SPS  = tInfo->sin_psi_op;
        tInfo->S3PS = 2.0*tInfo->CPS;
        tInfo->PST1[J] = PS*A[85];
        tInfo->PST2[J] = PS*A[86];
        tInfo->ST1[J] = sin(tInfo->PST1[J]);
        tInfo->CT1[J] = cos(tInfo->PST1[J]);
        tInfo->ST2[J] = sin(tInfo->PST2[J]);
        tInfo->CT2[J] = cos(tInfo->PST2[J]);
    }


    /*
     *  if X or Z or PS have changed we need to recompute things
     */
    if ( XChanged || ZChanged || PSChanged || !tInfo->DoneJ[J] ) {
        tInfo->X1[J] = X*tInfo->CT1[J] - Z*tInfo->ST1[J];
        tInfo->Z1[J] = X*tInfo->ST1[J] + Z*tInfo->CT1[J];
        tInfo->X2[J] = X*tInfo->CT2[J] - Z*tInfo->ST2[J];
        tInfo->Z2[J] = X*tInfo->ST2[J] + Z*tInfo->CT2[J];
    }


    /*
     *  cache the trig stuff. This doesnt depend on PS
     */
    if ( YChanged || !tInfo->DoneJ[J] ) {
        for (I=1; I<=3; I++ ){
            tInfo->YooP[J][I] = Y*tInfo->ooP[J][I];
            tInfo->YooQ[J][I] = Y*tInfo->ooQ[J][I];
            tInfo->CYPI[J][I] = cos(tInfo->YooP[J][I]);
            tInfo->CYQI[J][I] = cos(tInfo->YooQ[J][I]);
            tInfo->SYPI[J][I] = sin(tInfo->YooP[J][I]);
            tInfo->SYQI[J][I] = sin(tInfo->YooQ[J][I]);
        }
    }

    if ( XChanged || ZChanged || PSChanged || !tInfo->DoneJ[J] ) {
        for (K=1; K<=3; K++ ){
            tInfo->Z1ooR[J][K] = tInfo->Z1[J]*tInfo->ooR[J][K];
            tInfo->Z2ooS[J][K] = tInfo->Z2[J]*tInfo->ooS[J][K];
            tInfo->SZRK[J][K] = sin(tInfo->Z1ooR[J][K]);
            tInfo->CZSK[J][K] = cos(tInfo->Z2ooS[J][K]);
            tInfo->CZRK[J][K] = cos(tInfo->Z1ooR[J][K]);
            tInfo->SZSK[J][K] = sin(tInfo->Z2ooS[J][K]);
        }
    }

    if ( XChanged || YChanged || ZChanged || PSChanged || !tInfo->DoneJ[J] ) {
        for (I=1; I<=3; I++ ){
            for (K=1; K<=3; K++ ){
                tInfo->SQPR[J][I][K] = sqrt(tInfo->ooP2[J][I] + tInfo->ooR2[J][K]);
                tInfo->SQQS[J][I][K] = sqrt(tInfo->ooQ2[J][I] + tInfo->ooS2[J][K]);
                tInfo->EPR[J][I][K]  = exp(tInfo->X1[J]*tInfo->SQPR[J][I][K]);
                tInfo->EQS[J][I][K]  = exp(tInfo->X2[J]*tInfo->SQQS[J][I][K]);
            }
        }
    }


    /*
     *  Flag that we've been through once already with this value of J
     */
    tInfo->DoneJ[J] = TRUE;




//    printf("PS = %.15lf\n", PS);


    L  = 0;
    GX = 0.0;
    GY = 0.0;
    GZ = 0.0;


    for ( M=1;  M<=2; M++ ) {   //  M=1 IS FOR THE 1ST SUM ("PERP." SYMMETRY)
                                //  AND M=2 IS FOR THE SECOND SUM ("PARALL." SYMMETRY)

        for ( I=1;  I<=3; I++ ) {

            for ( K=1; K<= 3; K++ ) {

                for (N=1; N<=2; N++ ) { // N=1 IS FOR THE FIRST PART OF EACH COEFFICIENT
                                        // AND N=2 IS FOR THE SECOND ONE

                    for (NN=1; NN<=2; NN++ ) {  // NN = 1,2 FURTHER SPLITS THE COEFFICIENTS INTO 2 PARTS,
                                                // TO TAKE INTO ACCOUNT THE SCALE FACTOR DEPENDENCE
                        if (M == 1) {

                            FX = -tInfo->SQPR[J][I][K]*tInfo->EPR[J][I][K]*tInfo->CYPI[J][I]*tInfo->SZRK[J][K];
                            FY =  tInfo->EPR[J][I][K]*tInfo->SYPI[J][I]*tInfo->SZRK[J][K]*tInfo->ooP[J][I];
                            FZ = -tInfo->EPR[J][I][K]*tInfo->CYPI[J][I]*tInfo->CZRK[J][K]*tInfo->ooR[J][K];
                            //printf("1. FX, FY, FZ = %.10lf %.10lf %.10lf\n", FX, FY, FZ );

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
                                    HX = FX*tInfo->CPS;
                                    HY = FY*tInfo->CPS;
                                    HZ = FZ*tInfo->CPS;
                                } else {
                                    HX = FX*tInfo->CPS*X_SC;
                                    HY = FY*tInfo->CPS*X_SC;
                                    HZ = FZ*tInfo->CPS*X_SC;
                                }
                            }

                        } else {    //   M.EQ.2

                            FX = -tInfo->SPS*tInfo->SQQS[J][I][K]*tInfo->EQS[J][I][K]*tInfo->CYQI[J][I]*tInfo->CZSK[J][K];
                            FY = tInfo->SPS*tInfo->ooQ[J][I]*tInfo->EQS[J][I][K]*tInfo->SYQI[J][I]*tInfo->CZSK[J][K];
                            FZ = tInfo->SPS*tInfo->ooS[J][K]*tInfo->EQS[J][I][K]*tInfo->CYQI[J][I]*tInfo->SZSK[J][K];
                            //printf("2. FX, FY, FZ = %.10lf %.10lf %.10lf\n", FX, FY, FZ );

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
                                    HX = FX*tInfo->S3PS;
                                    HY = FY*tInfo->S3PS;
                                    HZ = FZ*tInfo->S3PS;
                                } else {
                                    HX = FX*tInfo->S3PS*X_SC;
                                    HY = FY*tInfo->S3PS*X_SC;
                                    HZ = FZ*tInfo->S3PS*X_SC;
                                }
                            }

                        }

                        ++L;

                        if (M == 1) {
                            HXR =  HX*tInfo->CT1[J] + HZ*tInfo->ST1[J];
                            HZR = -HX*tInfo->ST1[J] + HZ*tInfo->CT1[J];
                        } else {
                            HXR =  HX*tInfo->CT2[J] + HZ*tInfo->ST2[J];
                            HZR = -HX*tInfo->ST2[J] + HZ*tInfo->CT2[J];
                        }

                        GX = GX + HXR*A[L];
                        GY = GY + HY *A[L];
                        GZ = GZ + HZR*A[L];
//printf("GX, GY, GZ, A[L] = %.8lf %.8lf %.8lf %.8lf\n", GX, GY, GZ, A[L]);


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





/*
C
C   THIS S/R IS ALMOST IDENTICAL TO BIRK_TOT, BUT IT IS FOR THE SYMMETRIC MODE, IN WHICH
C   J_parallel IS AN EVEN FUNCTION OF Ygsm.
C
C
C    IOPBS -  BIRKELAND FIELD MODE FLAG:
C         IOPBS=0 - ALL COMPONENTS
C         IOPBS=1 - REGION 1, MODES 1 & 2 (SYMMETRIC !)
C         IOPBS=2 - REGION 2, MODES 1 & 2 (SYMMETRIC !)
*/
void     TS07D_BIRTOTSY( double PS, double X, double Y, double Z,
        double *BX11, double *BY11, double *BZ11, double *BX12, double *BY12, double *BZ12,
        double *BX21, double *BY21, double *BZ21, double *BX22, double *BY22, double *BZ22, LgmTsyg2007_Info *tInfo ){


    static double SH11[] = { -9e99, 4956703.683, -26922641.21, -11383659.85, 29604361.65,
                         -38919785.97,70230899.72,34993479.24,-90409215.02,30448713.69,
                         -48360257.19,-35556751.23,57136283.60,-8013815.613,30784907.86,
                          13501620.50,-35121638.52,50297295.45,-84200377.18,-46946852.58,
                          107526898.8,-39003263.47,59465850.17,47264335.10,-68892388.73,
                          3375901.533,-9181255.754,-4494667.217,10812618.51,-17351920.97,
                          27016083.00,18150032.11,-33186882.96,13340198.63,-19779685.30,
                          -17891788.15,21625767.23,16135.32442,133094.0241,-13845.61859,
                         -79159.98442,432.1215298,-85438.10368,1735.386707,41891.71284,
                          18158.14923,-105465.8135,-11685.73823,62297.34252,-10811.08476,
                          -87631.38186,9217.499261,52079.94529,-68.29127454,56023.02269,
                          -1246.029857,-27436.42793,-11972.61726,69607.08725,7702.743803,
                          -41114.36810,12.08269108,-21.30967022,-9.100782462,18.26855933,
                          -7.000685929,26.22390883,6.392164144,-21.99351743,2.294204157,
                          -16.10023369,-1.344314750,9.342121230,148.5493329,99.79912328,
                          70.78093196,35.23177574,47.45346891,58.44877918,139.8135237,
                          91.96485261,6.983488815,9.055554871,19.80484284,2.860045019,
                          .8213262337e-01,-.7962186676e-05 };

    static double SH12[] = { -9e99, -1210748.720,-52324903.95,-14158413.33,19426123.60,
                          6808641.947,-5138390.983,-1118600.499,-4675055.459,2059671.506,
                          -1373488.052,-114704.4353,-1435920.472,1438451.655,61199067.17,
                           16549301.39,-22802423.47,-7814550.995,5986478.728,1299443.190,
                           5352371.724,-2994351.520,1898553.337,203158.3658,2270182.134,
                          -618083.3112,-25950806.16,-7013783.326,9698792.575,3253693.134,
                          -2528478.464,-546323.4095,-2217735.237,1495336.589,-914647.4222,
                          -114374.1054,-1200441.634,-507068.4700,1163189.975,998411.8381,
                          -861919.3631,5252210.872,-11668550.16,-4113899.385,6972900.950,
                          -2546104.076,7704014.310,2273077.192,-5134603.198,256205.7901,
                          -589970.8086,-503821.0170,437612.8956,-2648640.128,5887640.735,
                           2074286.234,-3519291.144,1283847.104,-3885817.147,-1145936.942,
                           2589753.651,-408.7788403,1234.054185,739.8541716,-965.8068853,
                           3691.383679,-8628.635819,-2855.844091,5268.500178,-1774.372703,
                           5515.010707,1556.089289,-3665.434660,204.8672197,110.7748799,
                           87.36036207,5.522491330,31.06364270,73.57632579,281.5331360,
                           140.3461448,17.07537768,6.729732641,4.100970449,2.780422877,
                           .8742978101e-01,-.1028562327e-04 };

    static double SH21[] = { -9e99, -67763516.61,-49565522.84,10123356.08,51805446.10,
                          -51607711.68,164360662.1,-4662006.024,-191297217.6,-7204547.103,
                          30372354.93,-750371.9365,-36564457.17,61114395.65,45702536.50,
                          -9228894.939,-47893708.68,47290934.33,-149155112.0,4226520.638,
                          173588334.5,7998505.443,-33150962.72,832493.2094,39892545.84,
                          -11303915.16,-8901327.398,1751557.110,9382865.820,-9054707.868,
                          27918664.50,-788741.7146,-32481294.42,-2264443.753,9022346.503,
                          -233526.0185,-10856269.53,-244450.8850,1908295.272,185445.1967,
                          -1074202.863,41827.75224,-241553.7626,-20199.12580,123235.6084,
                          199501.4614,-1936498.464,-178857.4074,1044724.507,121044.9917,
                          -946479.9247,-91808.28803,532742.7569,-20742.28628,120633.2193,
                          10018.49534,-61599.11035,-98709.58977,959095.1770,88500.43489,
                          -517471.5287,-81.56122911,816.2472344,55.30711710,-454.5368824,
                          25.74693810,-202.5007350,-7.369350794,104.9429812,58.14049362,
                          -685.5919355,-51.71345683,374.0125033,247.9296982,159.2471769,
                          102.3151816,15.81062488,34.99767599,133.0832773,219.6475201,
                          107.9582783,10.00264684,7.718306072,25.22866153,5.013583103,
                          .8407754233e-01,-.9613356793e-05 };

    static double SH22[] = { -9e99, -43404887.31,8896854.538,-8077731.036,-10247813.65,
                          6346729.086,-9416801.212,-1921670.268,7805483.928,2299301.127,
                          4856980.170,-1253936.462,-4695042.690,54305735.91,-11158768.10,
                          10051771.85,12837129.47,-6380785.836,12387093.50,1687850.192,
                          -10492039.47,-5777044.862,-6916507.424,2855974.911,7027302.490,
                          -26176628.93,5387959.610,-4827069.106,-6193036.589,2511954.143,
                          -6205105.083,-553187.2984,5341386.847,3823736.361,3669209.068,
                          -1841641.700,-3842906.796,281561.7220,-5013124.630,379824.5943,
                          2436137.901,-76337.55394,548518.2676,42134.28632,-281711.3841,
                          -365514.8666,-2583093.138,-232355.8377,1104026.712,-131536.3445,
                           2320169.882,-174967.6603,-1127251.881,35539.82827,-256132.9284,
                          -19620.06116,131598.7965,169033.6708,1194443.500,107320.3699,
                          -510672.0036,1211.177843,-17278.19863,1140.037733,8347.612951,
                          -303.8408243,2405.771304,174.0634046,-1248.722950,-1231.229565,
                          -8666.932647,-754.0488385,3736.878824,227.2102611,115.9154291,
                          94.34364830,3.625357304,64.03192907,109.0743468,241.4844439,
                          107.7583478,22.36222385,6.282634037,27.79399216,2.270602235,
                          .8708605901e-01,-.1256706895e-04 };

    double      X_SC, FX11, FY11, FZ11, HX11, HY11, HZ11, FX12, FY12, FZ12, HX12, HY12, HZ12;
    double      FX21, FY21, FZ21, HX21, HY21, HZ21, FX22, FY22, FZ22, HX22, HY22, HZ22, XKAPPA;
    int         PSChanged, XChanged, YChanged, ZChanged;

    PSChanged = ( fabs(PS-tInfo->OLD_PS) > 1e-10 ) ? TRUE : FALSE;
    XChanged = ( fabs(X-tInfo->OLD_X) > 1e-10 ) ? TRUE : FALSE;
    YChanged = ( fabs(Y-tInfo->OLD_Y) > 1e-10 ) ? TRUE : FALSE;
    ZChanged = ( fabs(Z-tInfo->OLD_Z) > 1e-10 ) ? TRUE : FALSE;


    tInfo->CB_DPHI_B_RHO0.XKAPPA = tInfo->CB_BIRKPAR.XKAPPA1;               // FORWARDED IN BIRK_1N2
    X_SC                         = tInfo->CB_BIRKPAR.XKAPPA1 - 1.1;         // FORWARDED IN BIRK_SHL

    TS07D_BIR1N2SY( 1, 1, PS, X, Y, Z, &FX11, &FY11, &FZ11, tInfo );        //  REGION 1, MODE 1
    TS07D_BIRSH_SY( 0, PSChanged, XChanged, YChanged, ZChanged, SH11,PS, X_SC, X, Y, Z, &HX11, &HY11, &HZ11, tInfo );
    *BX11 = FX11 + HX11;
    *BY11 = FY11 + HY11;
    *BZ11 = FZ11 + HZ11;
//    printf( "FX11, FY11, FZ11, HX11, HY11, HZ11 = %g %g %g %g %g %g\n", FX11, FY11, FZ11, HX11, HY11, HZ11 );

    TS07D_BIR1N2SY( 1, 2, PS, X, Y, Z, &FX12, &FY12, &FZ12, tInfo );        //  REGION 1, MODE 2
    TS07D_BIRSH_SY( 1, PSChanged, XChanged, YChanged, ZChanged, SH12, PS, X_SC, X, Y, Z, &HX12, &HY12, &HZ12, tInfo );
    *BX12 = FX12 + HX12;
    *BY12 = FY12 + HY12;
    *BZ12 = FZ12 + HZ12;
//    printf( "FX12, FY12, FZ12, HX12, HY12, HZ12 = %g %g %g %g %g %g\n", FX12, FY12, FZ12, HX12, HY12, HZ12 );




    tInfo->CB_DPHI_B_RHO0.XKAPPA = tInfo->CB_BIRKPAR.XKAPPA2;               // FORWARDED IN BIRK_1N2
    X_SC                         = tInfo->CB_BIRKPAR.XKAPPA2 - 1.0;         // FORWARDED IN BIRK_SHL

    TS07D_BIR1N2SY( 2, 1, PS, X, Y, Z, &FX21, &FY21, &FZ21, tInfo );        //  REGION 2, MODE 1
    TS07D_BIRSH_SY( 2, PSChanged, XChanged, YChanged, ZChanged, SH21, PS, X_SC, X, Y, Z, &HX21, &HY21, &HZ21, tInfo );
    *BX21 = FX21 + HX21;
    *BY21 = FY21 + HY21;
    *BZ21 = FZ21 + HZ21;
//    printf( "FX21, FY21, FZ21, HX21, HY21, HZ21 = %g %g %g %g %g %g\n", FX21, FY21, FZ21, HX21, HY21, HZ21 );

    TS07D_BIR1N2SY( 2, 2, PS, X, Y, Z, &FX22, &FY22, &FZ22, tInfo );        //  REGION 2, MODE 2
    TS07D_BIRSH_SY( 3, PSChanged, XChanged, YChanged, ZChanged, SH22, PS, X_SC, X, Y, Z, &HX22, &HY22, &HZ22, tInfo );
    *BX22 = FX22 + HX22;
    *BY22 = FY22 + HY22;
    *BZ22 = FZ22 + HZ22;
//    printf( "FX22, FY22, FZ22, HX22, HY22, HZ22 = %g %g %g %g %g %g\n", FX22, FY22, FZ22, HX22, HY22, HZ22 );

    return;

}




















void    TS07D_BIR1N2SY( int NUMB, int MODE, double PS, double X, double Y, double Z, double *BX, double *BY, double *BZ, LgmTsyg2007_Info *tInfo ) {


    /*
     *   SEE NB# 6, P.60 and NB#7, P.35-...
     *
     *   THIS CODE IS VERY SIMILAR TO BIRK_1N2, BUT IT IS FOR THE "SYMMETRICAL" MODE, IN WHICH J_parallel
     *   IS A SYMMETRIC (EVEN) FUNCTION OF Ygsm
     *
     *   CALCULATES COMPONENTS  OF REGION 1/2 FIELD IN SPHERICAL COORDS.  DERIVED FROM THE S/R DIPDEF2C (WHICH
     *   DOES THE SAME JOB, BUT INPUT/OUTPUT THERE WAS IN SPHERICAL COORDS, WHILE HERE WE USE CARTESIAN ONES)
     *
     *   INPUT:  NUMB=1 (2) FOR REGION 1 (2) CURRENTS
     *           MODE=1 YIELDS SIMPLE SINUSOIDAL MLT VARIATION, WITH MAXIMUM CURRENT AT DAWN/DUSK MERIDIAN
     *   WHILE MODE=2 YIELDS THE SECOND HARMONIC.
     *
     *
     *      (1) DPHI:   HALF-DIFFERENCE (IN RADIANS) BETWEEN DAY AND NIGHT LATITUDE OF FAC OVAL AT IONOSPHERIC ALTITUDE;
     *                  TYPICAL VALUE: 0.06
     *      (2) B:      AN ASYMMETRY FACTOR AT HIGH-ALTITUDES;  FOR B=0, THE ONLY ASYMMETRY IS THAT FROM DPHI
     *                  TYPICAL VALUES: 0.35-0.70
     *      (3) RHO_0:  A FIXED PARAMETER, DEFINING THE DISTANCE RHO, AT WHICH THE LATITUDE SHIFT GRADUALLY SATURATES AND
     *                  STOPS INCREASING
     *                  ITS VALUE WAS ASSUMED FIXED, EQUAL TO 7.0.
     *      (4) XKAPPA: AN OVERALL SCALING FACTOR, WHICH CAN BE USED FOR CHANGING THE SIZE OF THE F.A.C. OVAL
     *
     *
     *
     */

    static double A11[] = { -9e99, .1618068350,-.1797957553,2.999642482,-.9322708978,
                   -.6811059760,.2099057262,-8.358815746,-14.86033550,.3838362986,
                   -16.30945494,4.537022847,2.685836007,27.97833029,6.330871059,
                   1.876532361,18.95619213,.9651528100,.4217195118,-.08957770020,
                   -1.823555887,.7457045438,-.5785916524,-1.010200918,.01112389357,
                   .09572927448,-.3599292276,8.713700514,.9763932955,3.834602998,
                   2.492118385,.7113544659 };

    static double A12[] = { -9e99, .7058026940,-.2845938535,5.715471266,-2.472820880,
                   -.7738802408,.3478293930,-11.37653694,-38.64768867,.6932927651,
                   -212.4017288,4.944204937,3.071270411,33.05882281,7.387533799,
                   2.366769108,79.22572682,.6154290178,.5592050551,-.1796585105,
                   -1.654932210,.7309108776,-.4926292779,-1.130266095,-.009613974555,
                   .1484586169,-.2215347198,7.883592948,.02768251655,2.950280953,
                   1.212634762,.5567714182 };

    static double A21[] = { -9e99, .1278764024,-.2320034273,1.805623266,-32.37241440,
                   -.9931490648,.3175085630,-2.492465814,-16.21600096,.2695393416,
                   -6.752691265,3.971794901,14.54477563,41.10158386,7.912889730,
                   1.258297372,9.583547721,1.014141963,.5104134759,-.1790430468,
                   -1.756358428,.7561986717,-.6775248254,-.04014016420,.01446794851,
                   .1200521731,-.2203584559,4.508963850,.8221623576,1.779933730,
                   1.102649543,.8867880020 };

    static double A22[] = { -9e99, .4036015198,-.3302974212,2.827730930,-45.44405830,
                   -1.611103927,.4927112073,-.003258457559,-49.59014949,.3796217108,
                   -233.7884098,4.312666980,18.05051709,28.95320323,11.09948019,
                   .7471649558,67.10246193,.5667096597,.6468519751,-.1560665317,
                   -1.460805289,.7719653528,-.6658988668,.2515179349E-05,
                   .02426021891,.1195003324,-.2625739255,4.377172556,.2421190547,
                   2.503482679,1.071587299,.7247997430 };


    double    BETA=0.9, RH=10.0;    // parameters of the tilt-dependent deformation of the untilted F.A.C. field
//    double    EPS=3.0;              // parameters of the tilt-dependent deformation of the untilted F.A.C. field
    double    RHO2, Xsc, Xsc2, Ysc, Zsc, Zsc2, RHO, RHOSQ;
    double    Rsc, Ysc2, PHI, SPHIC, CPHIC, BRACK, R1RH, R1RH2, R1RH3, PSIAS, PHIS, DPHISPHI;
    double    RHO2pRHOSQ, RHO2pRHOSQ2, DPHISRHO, DPHISDY, SPHICS, CPHICS, XS, ZS, BXS, BYAS;
    double    BZS, BRHOAS, BPHIAS, BRHO_S, BPHI_S, BY_S;


    tInfo->CB_DPHI_B_RHO0.B     = 0.5;
    tInfo->CB_DPHI_B_RHO0.RHO_0 = 7.0; RHO2 = 49.0;

    tInfo->CB_MODENUM.M = MODE;
    if (NUMB == 1) {
        tInfo->CB_DPHI_B_RHO0.DPHI   = 0.055;
        tInfo->CB_DTHETA.DTHETA = 0.06;
    } else if (NUMB == 2) {
        tInfo->CB_DPHI_B_RHO0.DPHI   = 0.030;
        tInfo->CB_DTHETA.DTHETA = 0.09;
    }

    Xsc = X*tInfo->CB_DPHI_B_RHO0.XKAPPA; Xsc2 = Xsc*Xsc;
    Ysc = Y*tInfo->CB_DPHI_B_RHO0.XKAPPA; Ysc2 = Ysc*Ysc;
    Zsc = Z*tInfo->CB_DPHI_B_RHO0.XKAPPA; Zsc2 = Zsc*Zsc;
    RHO = sqrt(Xsc2+Zsc2);
    RHOSQ = Xsc2+Zsc2;

    Rsc  = sqrt(Xsc2+Ysc2+Zsc2);        // SCALED

    if ( (Xsc == 0.0) && (Zsc == 0.0) ) {
        PHI   = 0.0;
        SPHIC = 0.0;
        CPHIC = 1.0;                // "C" means "CYLINDRICAL", TO DISTINGUISH FROM SPHERICAL PHI
    } else {
        PHI   = atan2(-Zsc, Xsc);            // FROM CARTESIAN TO CYLINDRICAL (RHO,PHI,Y)
        //SPHIC = sin(PHI);
        //CPHIC = cos(PHI);                // "C" means "CYLINDRICAL", TO DISTINGUISH FROM SPHERICAL PHI
        SPHIC = -Zsc/RHO;
        CPHIC =  Xsc/RHO;
    }


    BRACK = tInfo->CB_DPHI_B_RHO0.DPHI+tInfo->CB_DPHI_B_RHO0.B*RHO2/(RHO2+1.0)*(RHOSQ-1.0)/(RHO2+RHOSQ);
    R1RH  = (Rsc-1.0)/RH;
    R1RH2 = R1RH*R1RH; R1RH3 = R1RH2*R1RH;
double R1RH3p1 = 1.0+R1RH3;
double  cbrt_R1RH3p1 = cbrt( R1RH3p1 );
    PSIAS = BETA*PS/cbrt_R1RH3p1; // original is PSIAS=BETA*PS/(1.D0+R1RH**EPS)**(1.D0/EPS)

    PHIS = PHI-BRACK*SPHIC - PSIAS;
    DPHISPHI = 1.0-BRACK*CPHIC;
    RHO2pRHOSQ = RHO2+RHOSQ; RHO2pRHOSQ2 = RHO2pRHOSQ*RHO2pRHOSQ;
double  fff = 1.0/(RH*Rsc*R1RH3p1*cbrt_R1RH3p1);
    DPHISRHO = -2.0*tInfo->CB_DPHI_B_RHO0.B*RHO2*RHO/RHO2pRHOSQ2*SPHIC + BETA*PS*R1RH2*RHO*fff;
    DPHISDY= BETA*PS*R1RH2*Ysc*fff;

    //SPHICS = sin(PHIS); CPHICS = cos(PHIS);
    mysincos( PHIS, &SPHICS, &CPHICS );

    XS =  RHO*CPHICS;
    ZS = -RHO*SPHICS;

    if (NUMB == 1) {
        if (MODE == 1) TS07D_TWOCONSS( A11, XS, Ysc, ZS, &BXS, &BYAS, &BZS, tInfo );
        if (MODE == 2) TS07D_TWOCONSS( A12, XS, Ysc, ZS, &BXS, &BYAS, &BZS, tInfo );
    } else {
        if (MODE == 1) TS07D_TWOCONSS( A21, XS, Ysc, ZS, &BXS, &BYAS, &BZS, tInfo );
        if (MODE == 2) TS07D_TWOCONSS( A22, XS, Ysc, ZS, &BXS, &BYAS, &BZS, tInfo );
    }

    BRHOAS =  BXS*CPHICS - BZS*SPHICS;
    BPHIAS = -BXS*SPHICS - BZS*CPHICS;

    BRHO_S = BRHOAS*DPHISPHI                             * tInfo->CB_DPHI_B_RHO0.XKAPPA;    // SCALING
    BPHI_S = (BPHIAS-RHO*(BYAS*DPHISDY+BRHOAS*DPHISRHO)) * tInfo->CB_DPHI_B_RHO0.XKAPPA;
    BY_S   = BYAS*DPHISPHI                               * tInfo->CB_DPHI_B_RHO0.XKAPPA;

    *BX = BRHO_S*CPHIC-BPHI_S*SPHIC;
    *BY = BY_S;
    *BZ = -BRHO_S*SPHIC-BPHI_S*CPHIC;

    return;

}













void     TS07D_TWOCONSS( double *A, double X, double Y, double Z, double *BX, double *BY, double *BZ, LgmTsyg2007_Info *tInfo ){

    /*
     *
     *   DIFFERS FROM TWOCONES:  THIS S/R IS FOR THE "SYMMETRIC" MODE OF BIRKELAND CURRENTS IN THAT
     *                           HERE THE FIELD IS ROTATED BY 90 DEGS FOR M=1 AND BY 45 DEGS FOR M=2
     *
     *   ADDS FIELDS FROM TWO CONES (NORTHERN AND SOUTHERN), WITH A PROPER SYMMETRY OF THE CURRENT AND FIELD,
     *     CORRESPONDING TO THE REGION 1 BIRKELAND CURRENTS. (SEE NB #6, P.58).
     *
     */
    int             M;
    static double   HSQR2 = 0.707106781;
    double          XAS, YAS, BXAS, BYAS;
    double          BXN, BYN, BZN, BXS, BYS, BZS;

    M = tInfo->CB_MODENUM.M;

    if ( M == 1 ) {         // ROTATION BY 90 DEGS
        XAS =  Y;
        YAS = -X;
    } else {
        XAS = (X+Y)*HSQR2;
        YAS = (Y-X)*HSQR2;
    }

    TS07D_ONE_CONE( A, XAS, YAS, Z, &BXN, &BYN, &BZN, tInfo );
    TS07D_ONE_CONE( A, XAS, -YAS, -Z, &BXS, &BYS, &BZS, tInfo );

    BXAS = BXN - BXS;
    BYAS = BYN + BYS;
    *BZ  = BZN + BZS;

    if ( M == 1 ) {         // ROTATION BY 90 DEGS
        *BX = -BYAS;
        *BY =  BXAS;
    } else {
        *BX = (BXAS - BYAS)*HSQR2;
        *BY = (BXAS + BYAS)*HSQR2;
    }


    return;

}
















// This is one of the costliest routines -- mostly due to sin,cos,exp
void    TS07D_BIRSH_SY( int J, int PSChanged, int XChanged, int YChanged, int ZChanged, double *A, double PS, double X_SC,
                                double X, double Y, double Z, double *BX, double *BY, double *BZ, LgmTsyg2007_Info *tInfo ) {

    /*
     *   THIS S/R IS QUITE SIMILAR TO BIRK_SHL, BUT IT IS FOR THE SYMMETRIC MODE OF BIRKELAND CURRENT FIELD
     *     AND FOR THAT REASON THE FIELD COMPONENTS HAVE A DIFFERENT KIND OF SYMMETRY WITH RESPECT TO Y_gsm
     */

    int        L, M, I, N, NN, K;
    double    GX, GY, GZ, FX, FY, FZ, HX, HY, HZ, HXR, HZR;



    /*
     *  There are 4 different sets of coefficients that come to us. We added an argument to the
     *  calling routine so that we know which set we are dealing with. Then many things dont need to
     *  be recomputed. Probably dont need to cache all quantities here -- i.e. could prob optimize memory
     *  usage as well...
     */
    for ( I=1;  I<=3; I++ ) {
        if ( !tInfo->S_DoneJ[J] ) {

            tInfo->S_P[J][I] = A[72+I]; tInfo->S_ooP[J][I] = 1.0/tInfo->S_P[J][I]; tInfo->S_ooP2[J][I] = tInfo->S_ooP[J][I]*tInfo->S_ooP[J][I];
            tInfo->S_Q[J][I] = A[78+I]; tInfo->S_ooQ[J][I] = 1.0/tInfo->S_Q[J][I]; tInfo->S_ooQ2[J][I] = tInfo->S_ooQ[J][I]*tInfo->S_ooQ[J][I];

            tInfo->S_R[J][I] = A[75+I]; tInfo->S_ooR[J][I] = 1.0/tInfo->S_R[J][I]; tInfo->S_ooR2[J][I] = tInfo->S_ooR[J][I]*tInfo->S_ooR[J][I];
            tInfo->S_S[J][I] = A[81+I]; tInfo->S_ooS[J][I] = 1.0/tInfo->S_S[J][I]; tInfo->S_ooS2[J][I] = tInfo->S_ooS[J][I]*tInfo->S_ooS[J][I];
        }
    }


    /*
     *  if PS has changed we need to recompute things
     */
    if ( PSChanged || !tInfo->S_DoneJ[J] ) {
        //tInfo->S_CPS  = tInfo->S_cos_psi_op;
        //tInfo->S_SPS  = tInfo->S_sin_psi_op;
        //tInfo->S_S3PS = 2.0*tInfo->S_CPS;
        tInfo->S_CPS  = tInfo->cos_psi_op;
        tInfo->S_SPS  = tInfo->sin_psi_op;
        tInfo->S_S3PS = 2.0*tInfo->CPS;
        tInfo->S_PST1[J] = PS*A[85];
        tInfo->S_PST2[J] = PS*A[86];
        tInfo->S_ST1[J] = sin(tInfo->S_PST1[J]);
        tInfo->S_CT1[J] = cos(tInfo->S_PST1[J]);
        //sincos( tInfo->S_PST1[J], &(tInfo->S_ST1[J]), &(tInfo->S_CT1[J]) );
        tInfo->S_ST2[J] = sin(tInfo->S_PST2[J]);
        tInfo->S_CT2[J] = cos(tInfo->S_PST2[J]);
        //sincos( tInfo->S_PST2[J], &(tInfo->S_ST2[J]), &(tInfo->S_CT2[J]) );
    }


    /*
     *  if X or Z or PS have changed we need to recompute things
     */
    if ( XChanged || ZChanged || PSChanged || !tInfo->S_DoneJ[J] ) {
        tInfo->S_X1[J] = X*tInfo->S_CT1[J] - Z*tInfo->S_ST1[J];
        tInfo->S_Z1[J] = X*tInfo->S_ST1[J] + Z*tInfo->S_CT1[J];
        tInfo->S_X2[J] = X*tInfo->S_CT2[J] - Z*tInfo->S_ST2[J];
        tInfo->S_Z2[J] = X*tInfo->S_ST2[J] + Z*tInfo->S_CT2[J];
    }


    /*
     *  cache the trig stuff. This doesnt depend on PS
     */
    if ( YChanged || !tInfo->S_DoneJ[J] ) {
        for (I=1; I<=3; I++ ){
            tInfo->S_YooP[J][I] = Y*tInfo->S_ooP[J][I];
            tInfo->S_YooQ[J][I] = Y*tInfo->S_ooQ[J][I];

            //tInfo->S_SYQI[J][I] = sin(tInfo->S_YooQ[J][I]);
            //tInfo->S_CYQI[J][I] = cos(tInfo->S_YooQ[J][I]);
            sincos( tInfo->S_YooQ[J][I], &(tInfo->S_SYQI[J][I]), &(tInfo->S_CYQI[J][I]) );

            //tInfo->S_SYPI[J][I] = sin(tInfo->S_YooP[J][I]);
            //tInfo->S_CYPI[J][I] = cos(tInfo->S_YooP[J][I]);
            sincos( tInfo->S_YooP[J][I], &(tInfo->S_SYPI[J][I]), &(tInfo->S_CYPI[J][I]) );
        }
    }

    if ( XChanged || ZChanged || PSChanged || !tInfo->DoneJ[J] ) {
        for (K=1; K<=3; K++ ){
            tInfo->S_Z1ooR[J][K] = tInfo->S_Z1[J]*tInfo->S_ooR[J][K];
            tInfo->S_Z2ooS[J][K] = tInfo->S_Z2[J]*tInfo->S_ooS[J][K];

            //tInfo->S_SZRK[J][K] = sin(tInfo->S_Z1ooR[J][K]);
            //tInfo->S_CZRK[J][K] = cos(tInfo->S_Z1ooR[J][K]);
            sincos( tInfo->S_Z1ooR[J][K], &(tInfo->S_SZRK[J][K]), &(tInfo->S_CZRK[J][K]) );

            //tInfo->S_SZSK[J][K] = sin(tInfo->S_Z2ooS[J][K]);
            //tInfo->S_CZSK[J][K] = cos(tInfo->S_Z2ooS[J][K]);
            sincos( tInfo->S_Z2ooS[J][K], &(tInfo->S_SZSK[J][K]), &(tInfo->S_CZSK[J][K]) );
        }
    }

    if ( XChanged || YChanged || ZChanged || PSChanged || !tInfo->S_DoneJ[J] ) {
        for (I=1; I<=3; I++ ){
            for (K=1; K<=3; K++ ){
                tInfo->S_SQPR[J][I][K] = sqrt(tInfo->S_ooP2[J][I] + tInfo->S_ooR2[J][K]);
                tInfo->S_SQQS[J][I][K] = sqrt(tInfo->S_ooQ2[J][I] + tInfo->S_ooS2[J][K]);
                tInfo->S_EPR[J][I][K]  = exp(tInfo->S_X1[J]*tInfo->S_SQPR[J][I][K]);
                tInfo->S_EQS[J][I][K]  = exp(tInfo->S_X2[J]*tInfo->S_SQQS[J][I][K]);
            }
        }
    }


    /*
     *  Flag that we've been through once already with this value of J
     */
    tInfo->S_DoneJ[J] = TRUE;






    L  = 0;
    GX = 0.0;
    GY = 0.0;
    GZ = 0.0;


    for ( M=1;  M<=2; M++ ) {   //  M=1 IS FOR THE 1ST SUM ("PERP." SYMMETRY)
                                //  AND M=2 IS FOR THE SECOND SUM ("PARALL." SYMMETRY)

        for ( I=1;  I<=3; I++ ) {

            for ( K=1; K<= 3; K++ ) {

                for (N=1; N<=2; N++ ) { // N=1 IS FOR THE FIRST PART OF EACH COEFFICIENT
                                        // AND N=2 IS FOR THE SECOND ONE

                    for (NN=1; NN<=2; NN++ ) {  // NN = 1,2 FURTHER SPLITS THE COEFFICIENTS INTO 2 PARTS,
                                                // TO TAKE INTO ACCOUNT THE SCALE FACTOR DEPENDENCE
                        if (M == 1) {

                            FX =  tInfo->S_SQPR[J][I][K]*tInfo->S_EPR[J][I][K]*tInfo->S_SYPI[J][I]*tInfo->S_SZRK[J][K];
                            FY =  tInfo->S_EPR[J][I][K]*tInfo->S_CYPI[J][I]*tInfo->S_SZRK[J][K]*tInfo->S_ooP[J][I];
                            FZ =  tInfo->S_EPR[J][I][K]*tInfo->S_SYPI[J][I]*tInfo->S_CZRK[J][K]*tInfo->S_ooR[J][K];

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
                                    HX = FX*tInfo->S_CPS;
                                    HY = FY*tInfo->S_CPS;
                                    HZ = FZ*tInfo->S_CPS;
                                } else {
                                    HX = FX*tInfo->S_CPS*X_SC;
                                    HY = FY*tInfo->S_CPS*X_SC;
                                    HZ = FZ*tInfo->S_CPS*X_SC;
                                }
                            }

                        } else {    //   M.EQ.2

                            FX =  tInfo->S_SPS*tInfo->S_SQQS[J][I][K]*tInfo->S_EQS[J][I][K]*tInfo->S_SYQI[J][I]*tInfo->S_CZSK[J][K];
                            FY =  tInfo->S_SPS*tInfo->S_ooQ[J][I]*tInfo->S_EQS[J][I][K]*tInfo->S_CYQI[J][I]*tInfo->S_CZSK[J][K];
                            FZ = -tInfo->S_SPS*tInfo->S_ooS[J][K]*tInfo->S_EQS[J][I][K]*tInfo->S_SYQI[J][I]*tInfo->S_SZSK[J][K];

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
                                    HX = FX*tInfo->S_S3PS;
                                    HY = FY*tInfo->S_S3PS;
                                    HZ = FZ*tInfo->S_S3PS;
                                } else {
                                    HX = FX*tInfo->S_S3PS*X_SC;
                                    HY = FY*tInfo->S_S3PS*X_SC;
                                    HZ = FZ*tInfo->S_S3PS*X_SC;
                                }
                            }

                        }

                        ++L;

                        if (M == 1) {
                            HXR =  HX*tInfo->S_CT1[J] + HZ*tInfo->S_ST1[J];
                            HZR = -HX*tInfo->S_ST1[J] + HZ*tInfo->S_CT1[J];
                        } else {
                            HXR =  HX*tInfo->S_CT2[J] + HZ*tInfo->S_ST2[J];
                            HZR = -HX*tInfo->S_ST2[J] + HZ*tInfo->S_CT2[J];
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

