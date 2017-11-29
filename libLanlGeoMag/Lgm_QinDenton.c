/*! \file Lgm_QinDenton.c
 *
 *  \brief Routines for reading in and setting the QinDenton magnetic field model input parameters.
 *
 *
 *
 *  \author M.G. Henderson
 *  \date   20??
 *
 *
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <string.h>
#include "Lgm/Lgm_CTrans.h"
#include "Lgm/Lgm_MagModelInfo.h"
#include "Lgm/Lgm_DynamicMemory.h"
#include "Lgm/Lgm_QinDenton.h"


#define NMAX    (3*1440)


#ifndef LGM_INDEX_DATA_DIR
#warning "hard-coding LGM_INDEX_DATA_DIR because it was not in config.h"
#define LGM_INDEX_DATA_DIR    /usr/local/share/LanlGeoMag/Data
#endif

Lgm_QinDenton *Lgm_init_QinDenton( int Verbose ) {
  // call this from C, never from Python
    Lgm_QinDenton      *q;
    q = (Lgm_QinDenton *) calloc( 1, sizeof(*q) );
    Lgm_init_QinDentonDefaults(q, Verbose);
    return q;
}

void Lgm_init_QinDentonDefaults( Lgm_QinDenton *q, int Verbose ) {
  // call this from C or Python
    q->Verbosity = Verbose;
    LGM_ARRAY_1D( q->Date,          NMAX, long int );
    LGM_ARRAY_1D( q->MJD,           NMAX, double );

    LGM_ARRAY_2D( q->IsoTimeStr,    NMAX, 80, char );

    LGM_ARRAY_1D( q->Year,          NMAX, int );
    LGM_ARRAY_1D( q->Month,         NMAX, int );
    LGM_ARRAY_1D( q->Day,           NMAX, int );
    LGM_ARRAY_1D( q->Hour,          NMAX, int );
    LGM_ARRAY_1D( q->Minute,        NMAX, int );
    LGM_ARRAY_1D( q->Second,        NMAX, int );

    LGM_ARRAY_1D( q->ByIMF,         NMAX, double );
    LGM_ARRAY_1D( q->BzIMF,         NMAX, double );
    LGM_ARRAY_1D( q->V_SW,          NMAX, double );
    LGM_ARRAY_1D( q->Den_P,         NMAX, double );
    LGM_ARRAY_1D( q->Pdyn,          NMAX, double );

    LGM_ARRAY_1D( q->G1,            NMAX, double );
    LGM_ARRAY_1D( q->G2,            NMAX, double );
    LGM_ARRAY_1D( q->G3,            NMAX, double );

    LGM_ARRAY_1D( q->ByIMF_status,  NMAX, int );
    LGM_ARRAY_1D( q->BzIMF_status,  NMAX, int );
    LGM_ARRAY_1D( q->V_SW_status,   NMAX, int );
    LGM_ARRAY_1D( q->Den_P_status,  NMAX, int );
    LGM_ARRAY_1D( q->Pdyn_status,   NMAX, int );

    LGM_ARRAY_1D( q->G1_status,     NMAX, int );
    LGM_ARRAY_1D( q->G2_status,     NMAX, int );
    LGM_ARRAY_1D( q->G3_status,     NMAX, int );

    LGM_ARRAY_1D( q->fKp,           NMAX, double );
    LGM_ARRAY_1D( q->akp3,          NMAX, double );
    LGM_ARRAY_1D( q->Dst,           NMAX, double );

    LGM_ARRAY_1D( q->Bz1,           NMAX, double );
    LGM_ARRAY_1D( q->Bz2,           NMAX, double );
    LGM_ARRAY_1D( q->Bz3,           NMAX, double );
    LGM_ARRAY_1D( q->Bz4,           NMAX, double );
    LGM_ARRAY_1D( q->Bz5,           NMAX, double );
    LGM_ARRAY_1D( q->Bz6,           NMAX, double );

    LGM_ARRAY_1D( q->W1,            NMAX, double );
    LGM_ARRAY_1D( q->W2,            NMAX, double );
    LGM_ARRAY_1D( q->W3,            NMAX, double );
    LGM_ARRAY_1D( q->W4,            NMAX, double );
    LGM_ARRAY_1D( q->W5,            NMAX, double );
    LGM_ARRAY_1D( q->W6,            NMAX, double );

    LGM_ARRAY_1D( q->W1_status,     NMAX, int );
    LGM_ARRAY_1D( q->W2_status,     NMAX, int );
    LGM_ARRAY_1D( q->W3_status,     NMAX, int );
    LGM_ARRAY_1D( q->W4_status,     NMAX, int );
    LGM_ARRAY_1D( q->W5_status,     NMAX, int );
    LGM_ARRAY_1D( q->W6_status,     NMAX, int );
}

void  Lgm_destroy_QinDenton_children( Lgm_QinDenton *q ) {
  // call this from Python, probably not C
    // first free arrays inside structure
    LGM_ARRAY_1D_FREE( q->Date );
    LGM_ARRAY_1D_FREE( q->MJD );

    LGM_ARRAY_2D_FREE( q->IsoTimeStr );

    LGM_ARRAY_1D_FREE( q->Year );
    LGM_ARRAY_1D_FREE( q->Month );
    LGM_ARRAY_1D_FREE( q->Day );
    LGM_ARRAY_1D_FREE( q->Hour );
    LGM_ARRAY_1D_FREE( q->Minute );
    LGM_ARRAY_1D_FREE( q->Second );

    LGM_ARRAY_1D_FREE( q->ByIMF );
    LGM_ARRAY_1D_FREE( q->BzIMF );
    LGM_ARRAY_1D_FREE( q->V_SW );
    LGM_ARRAY_1D_FREE( q->Den_P );
    LGM_ARRAY_1D_FREE( q->Pdyn );

    LGM_ARRAY_1D_FREE( q->G1 );
    LGM_ARRAY_1D_FREE( q->G2 );
    LGM_ARRAY_1D_FREE( q->G3 );

    LGM_ARRAY_1D_FREE( q->ByIMF_status );
    LGM_ARRAY_1D_FREE( q->BzIMF_status );
    LGM_ARRAY_1D_FREE( q->V_SW_status );
    LGM_ARRAY_1D_FREE( q->Den_P_status );
    LGM_ARRAY_1D_FREE( q->Pdyn_status );

    LGM_ARRAY_1D_FREE( q->G1_status );
    LGM_ARRAY_1D_FREE( q->G2_status );
    LGM_ARRAY_1D_FREE( q->G3_status );

    LGM_ARRAY_1D_FREE( q->fKp );
    LGM_ARRAY_1D_FREE( q->akp3 );
    LGM_ARRAY_1D_FREE( q->Dst );

    LGM_ARRAY_1D_FREE( q->Bz1 );
    LGM_ARRAY_1D_FREE( q->Bz2 );
    LGM_ARRAY_1D_FREE( q->Bz3 );
    LGM_ARRAY_1D_FREE( q->Bz4 );
    LGM_ARRAY_1D_FREE( q->Bz5 );
    LGM_ARRAY_1D_FREE( q->Bz6 );

    LGM_ARRAY_1D_FREE( q->W1 );
    LGM_ARRAY_1D_FREE( q->W2 );
    LGM_ARRAY_1D_FREE( q->W3 );
    LGM_ARRAY_1D_FREE( q->W4 );
    LGM_ARRAY_1D_FREE( q->W5 );
    LGM_ARRAY_1D_FREE( q->W6 );

    LGM_ARRAY_1D_FREE( q->W1_status );
    LGM_ARRAY_1D_FREE( q->W2_status );
    LGM_ARRAY_1D_FREE( q->W3_status );
    LGM_ARRAY_1D_FREE( q->W4_status );
    LGM_ARRAY_1D_FREE( q->W5_status );
    LGM_ARRAY_1D_FREE( q->W6_status );

    return;

}

void  Lgm_destroy_QinDenton( Lgm_QinDenton *q ) {
  // call this from C, probably not from Python
  Lgm_destroy_QinDenton_children(q);  
    // then free structure itself
    free( q );
    return;
}


void Lgm_read_QinDenton( long int Date, Lgm_QinDenton *q ) {

    FILE        *fp;
    double      Time, MJD, Prev_MJD, Next_MJD, Prev_UT, Next_UT;
    long int    n, Prev_Date, Next_Date;
    int         Year, Month, Day, Doy, j, done;
    int         Prev_Year, Prev_Month, Prev_Day, Next_Year, Next_Month, Next_Day;
    int         success1, success2, success3;
    char        *Line, *Filename;
    static char *ftype[] = {"1min", "1hr" };
    char        *Path, QinDentonPath[2048];
    Lgm_CTrans  *c = Lgm_init_ctrans(0);
    char        IsoTimeStr[80];
    int         nMatches, tYear, tMonth, tDay, tHour, tMinute, tSecond, ByIMF_status, BzIMF_status, V_SW_status, Den_P_status, Pdyn_status, G1_status;
    int         G2_status, G3_status, W1_status, W2_status, W3_status, W4_status, W5_status, W6_status;
    double      ByIMF, BzIMF, V_SW, Den_P, Pdyn, G1, G2, G3, fKp, akp3, Dst, Bz1, Bz2, Bz3, Bz4, Bz5, Bz6, W1, W2, W3, W4, W5, W6, tTime, tMJD;


    Path = getenv( "QIN_DENTON_PATH" );
    if ( Path == NULL ) {
        strcpy( QinDentonPath, LGM_INDEX_DATA_DIR );
        strcat( QinDentonPath, "/QinDenton" );
    } else {

        /*
         * Test for existence
         */
        struct stat sts;
        if ( ( stat( Path, &sts ) ) == -1 ) {
            strcpy( QinDentonPath, LGM_INDEX_DATA_DIR );
            strcat( QinDentonPath, "/QinDenton" );
            printf("Environment variable QIN_DENTON_PATH points to a non-existent directory: %s. Setting QinDentonPath to: %s \n", Path, QinDentonPath );
        } else {
            strcpy( QinDentonPath, Path );
        }
    }


    

    Lgm_Doy( Date, &Year, &Month, &Day, &Doy);
    MJD = Lgm_MJD( Year, Month, Day, 12.0, LGM_TIME_SYS_UTC, c );


    Prev_MJD = MJD-1.0;
    Lgm_mjd_to_ymdh( Prev_MJD, &Prev_Date, &Prev_Year, &Prev_Month, &Prev_Day, &Prev_UT );

    Next_MJD = MJD+1.0;
    Lgm_mjd_to_ymdh( Next_MJD, &Next_Date, &Next_Year, &Next_Month, &Next_Day, &Next_UT );




    Filename = (char *)calloc( 512, sizeof(char) );
    Line = (char *)calloc( 2050, sizeof(char) );
    n = 0;



    // Read in Previous Date -- Try 1min first, 1hr next...
    tTime = -1e30;
    j = 0; done = FALSE; success1 = FALSE;
    while ( !done ){

        sprintf( Filename, "%s/%4d/QinDenton_%8ld_%s.txt", QinDentonPath, Prev_Year, Prev_Date, ftype[j] );
        if ( (fp = fopen( Filename, "r" )) != NULL ) {
            while( fgets( Line, 2048, fp ) != NULL ) {
                if ( Line[0] != '#' ) {
                    nMatches = sscanf( Line, "%s %d %d %d %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %d %d %d %d %d %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d %d %d %d %d",
                                IsoTimeStr, &tYear, &tMonth, &tDay, &tHour, &tMinute, &tSecond, &ByIMF, &BzIMF, &V_SW, &Den_P, &Pdyn, &G1, &G2, &G3,
                                &ByIMF_status, &BzIMF_status, &V_SW_status, &Den_P_status, &Pdyn_status, &G1_status, &G2_status, &G3_status,
                                &fKp, &akp3, &Dst, &Bz1, &Bz2, &Bz3, &Bz4, &Bz5, &Bz6, &W1, &W2, &W3, &W4, &W5, &W6, &W1_status, &W2_status, &W3_status, &W4_status, &W5_status, &W6_status );


                    if ( nMatches == 44 ) {

                        tTime = tHour + tMinute/60.0 + tSecond/3600.0;
                        tMJD = Lgm_MJD( tYear, tMonth, tDay, tTime, LGM_TIME_SYS_UTC, c );

                        if ( (MJD > 33282.0) && ((n==0) || (MJD > q->MJD[n])) ) {  // make sure MJD > Jan 1, 1950 and time is increasing

                            strcpy( q->IsoTimeStr[n], IsoTimeStr );
                            q->Year[n]   = tYear;
                            q->Month[n]  = tMonth;
                            q->Day[n]    = tDay;
                            q->Hour[n]   = tHour;
                            q->Minute[n] = tMinute;
                            q->Second[n] = tSecond;
                            q->ByIMF[n]  = ByIMF;
                            q->BzIMF[n]  = BzIMF;
                            q->V_SW[n]   = V_SW;
                            q->Den_P[n]  = Den_P;
                            q->Pdyn[n]   = Pdyn;
                            q->G1[n]     = G1;
                            q->G2[n]     = G2;
                            q->G3[n]     = G3;
                            q->ByIMF_status[n] = ByIMF_status;
                            q->BzIMF_status[n] = BzIMF_status;
                            q->V_SW_status[n]  = V_SW_status;
                            q->Den_P_status[n] = Den_P_status;
                            q->Pdyn_status[n]  = Pdyn_status;
                            q->G1_status[n]    = G1_status;
                            q->G2_status[n]    = G2_status;
                            q->G3_status[n]    = G3_status;
                            q->fKp[n]  = fKp;
                            q->akp3[n] = akp3;
                            q->Dst[n]  = Dst;
                            q->Bz1[n]  = Bz1;
                            q->Bz2[n]  = Bz2;
                            q->Bz3[n]  = Bz3;
                            q->Bz4[n]  = Bz4;
                            q->Bz5[n]  = Bz5;
                            q->Bz6[n]  = Bz6;
                            q->W1[n]   = W1;
                            q->W2[n]   = W2;
                            q->W3[n]   = W3;
                            q->W4[n]   = W4;
                            q->W5[n]   = W5;
                            q->W6[n]   = W6;
                            q->W1_status[n] = W1_status;
                            q->W2_status[n] = W2_status;
                            q->W3_status[n] = W3_status;
                            q->W4_status[n] = W4_status;
                            q->W5_status[n] = W5_status;
                            q->W6_status[n] = W6_status;

                            q->MJD[n] = tMJD;
                            ++n;

                        } else {
                            printf( "Warning. Times may be corrupted in QinDenton File: %s   Time = %g  MJD = %g\n", Filename, tTime, tMJD );
                        }
                    }
                }
            }
            fclose( fp );
            done = TRUE;
            success1 = TRUE;
        } else {
            // only complain if this is our last try.
            if (j==1) {
                if ( Path == NULL ) { //i.e. QIN_DENTON_PATH environment variable not set.
                    printf( "Cannot open %s file. Try setting QIN_DENTON_PATH environment variable to path containing QinDenton files.\n", Filename );
                } else {
                    printf( "Cannot open %s file.\n", Filename );
                }
            }
        }

        if (j == 1) {
            done = TRUE;
        } else {
            ++j;
        }

    }


    // Read in Current Date
    j = 0; done = FALSE; success2 = FALSE;
    while ( !done ){

        sprintf( Filename, "%s/%4d/QinDenton_%8ld_%s.txt", QinDentonPath, Year, Date, ftype[j] );
        if ( (fp = fopen( Filename, "r" )) != NULL ) {
            while( fgets( Line, 2048, fp ) != NULL ) {
                if ( Line[0] != '#' ) {


                    nMatches = sscanf( Line, "%s %d %d %d %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %d %d %d %d %d %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d %d %d %d %d",
                                IsoTimeStr, &tYear, &tMonth, &tDay, &tHour, &tMinute, &tSecond, &ByIMF, &BzIMF, &V_SW, &Den_P, &Pdyn, &G1, &G2, &G3,
                                &ByIMF_status, &BzIMF_status, &V_SW_status, &Den_P_status, &Pdyn_status, &G1_status, &G2_status, &G3_status,
                                &fKp, &akp3, &Dst, &Bz1, &Bz2, &Bz3, &Bz4, &Bz5, &Bz6, &W1, &W2, &W3, &W4, &W5, &W6, &W1_status, &W2_status, &W3_status, &W4_status, &W5_status, &W6_status );

                    if ( nMatches == 44 ) {

                        tTime = tHour + tMinute/60.0 + tSecond/3600.0;
                        tMJD = Lgm_MJD( tYear, tMonth, tDay, tTime, LGM_TIME_SYS_UTC, c );

                        if ( (MJD > 33282.0) && ((n==0) || (MJD > q->MJD[n])) ) {  // make sure MJD > Jan 1, 1950 and time is increasing
                            strcpy( q->IsoTimeStr[n], IsoTimeStr );
                            q->Year[n]   = tYear;
                            q->Month[n]  = tMonth;
                            q->Day[n]    = tDay;
                            q->Hour[n]   = tHour;
                            q->Minute[n] = tMinute;
                            q->Second[n] = tSecond;
                            q->ByIMF[n]  = ByIMF;
                            q->BzIMF[n]  = BzIMF;
                            q->V_SW[n]   = V_SW;
                            q->Den_P[n]  = Den_P;
                            q->Pdyn[n]   = Pdyn;
                            q->G1[n]     = G1;
                            q->G2[n]     = G2;
                            q->G3[n]     = G3;
                            q->ByIMF_status[n] = ByIMF_status;
                            q->BzIMF_status[n] = BzIMF_status;
                            q->V_SW_status[n]  = V_SW_status;
                            q->Den_P_status[n] = Den_P_status;
                            q->Pdyn_status[n]  = Pdyn_status;
                            q->G1_status[n]    = G1_status;
                            q->G2_status[n]    = G2_status;
                            q->G3_status[n]    = G3_status;
                            q->fKp[n]  = fKp;
                            q->akp3[n] = akp3;
                            q->Dst[n]  = Dst;
                            q->Bz1[n]  = Bz1;
                            q->Bz2[n]  = Bz2;
                            q->Bz3[n]  = Bz3;
                            q->Bz4[n]  = Bz4;
                            q->Bz5[n]  = Bz5;
                            q->Bz6[n]  = Bz6;
                            q->W1[n]   = W1;
                            q->W2[n]   = W2;
                            q->W3[n]   = W3;
                            q->W4[n]   = W4;
                            q->W5[n]   = W5;
                            q->W6[n]   = W6;
                            q->W1_status[n] = W1_status;
                            q->W2_status[n] = W2_status;
                            q->W3_status[n] = W3_status;
                            q->W4_status[n] = W4_status;
                            q->W5_status[n] = W5_status;
                            q->W6_status[n] = W6_status;

                            q->MJD[n] = tMJD;
                            ++n;

                        } else {
                            printf( "Warning. Times may be corrupted in QinDenton File: %s   Time = %g  MJD = %g\n", Filename, tTime, tMJD );
                        }

                    }

                }
            }
            fclose( fp );
            done = TRUE;
            success2 = TRUE;
        } else {
            // only complain if this is our last try.
            if ( Path == NULL ) { //i.e. QIN_DENTON_PATH environment variable not set.
                printf( "Cannot open %s file. Try setting QIN_DENTON_PATH environment variable to path containing QinDenton files.\n", Filename );
            } else {
                printf( "Cannot open %s file.\n", Filename );
            }
        }

        if (j == 1) {
            done = TRUE;
        } else {
            ++j;
        }

    }



    // Read in Next Date
    j = 0; done = FALSE; success3 = FALSE;
    while ( !done ){

        sprintf( Filename, "%s/%4d/QinDenton_%8ld_%s.txt", QinDentonPath, Next_Year, Next_Date, ftype[j] );
        if ( (fp = fopen( Filename, "r" )) != NULL ) {
            while( fgets( Line, 2048, fp ) != NULL ) {
                if ( Line[0] != '#' ) {
                    nMatches = sscanf( Line, "%s %d %d %d %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %d %d %d %d %d %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d %d %d %d %d",
                                IsoTimeStr, &tYear, &tMonth, &tDay, &tHour, &tMinute, &tSecond, &ByIMF, &BzIMF, &V_SW, &Den_P, &Pdyn, &G1, &G2, &G3,
                                &ByIMF_status, &BzIMF_status, &V_SW_status, &Den_P_status, &Pdyn_status, &G1_status, &G2_status, &G3_status,
                                &fKp, &akp3, &Dst, &Bz1, &Bz2, &Bz3, &Bz4, &Bz5, &Bz6, &W1, &W2, &W3, &W4, &W5, &W6, &W1_status, &W2_status, &W3_status, &W4_status, &W5_status, &W6_status );

                    if ( nMatches == 44 ) {

                        tTime = tHour + tMinute/60.0 + tSecond/3600.0;
                        tMJD = Lgm_MJD( tYear, tMonth, tDay, tTime, LGM_TIME_SYS_UTC, c );

                        if ( (MJD > 33282.0) && ((n==0) || (MJD > q->MJD[n])) ) {  // make sure MJD > Jan 1, 1950 and time is increasing
                            strcpy( q->IsoTimeStr[n], IsoTimeStr );
                            q->Year[n]   = tYear;
                            q->Month[n]  = tMonth;
                            q->Day[n]    = tDay;
                            q->Hour[n]   = tHour;
                            q->Minute[n] = tMinute;
                            q->Second[n] = tSecond;
                            q->ByIMF[n]  = ByIMF;
                            q->BzIMF[n]  = BzIMF;
                            q->V_SW[n]   = V_SW;
                            q->Den_P[n]  = Den_P;
                            q->Pdyn[n]   = Pdyn;
                            q->G1[n]     = G1;
                            q->G2[n]     = G2;
                            q->G3[n]     = G3;
                            q->ByIMF_status[n] = ByIMF_status;
                            q->BzIMF_status[n] = BzIMF_status;
                            q->V_SW_status[n]  = V_SW_status;
                            q->Den_P_status[n] = Den_P_status;
                            q->Pdyn_status[n]  = Pdyn_status;
                            q->G1_status[n]    = G1_status;
                            q->G2_status[n]    = G2_status;
                            q->G3_status[n]    = G3_status;
                            q->fKp[n]  = fKp;
                            q->akp3[n] = akp3;
                            q->Dst[n]  = Dst;
                            q->Bz1[n]  = Bz1;
                            q->Bz2[n]  = Bz2;
                            q->Bz3[n]  = Bz3;
                            q->Bz4[n]  = Bz4;
                            q->Bz5[n]  = Bz5;
                            q->Bz6[n]  = Bz6;
                            q->W1[n]   = W1;
                            q->W2[n]   = W2;
                            q->W3[n]   = W3;
                            q->W4[n]   = W4;
                            q->W5[n]   = W5;
                            q->W6[n]   = W6;
                            q->W1_status[n] = W1_status;
                            q->W2_status[n] = W2_status;
                            q->W3_status[n] = W3_status;
                            q->W4_status[n] = W4_status;
                            q->W5_status[n] = W5_status;
                            q->W6_status[n] = W6_status;

                            q->MJD[n] = tMJD;
                            ++n;

                        } else {
                            printf( "Warning. Times may be corrupted in QinDenton File: %s   Time = %g  MJD = %g\n", Filename, tTime, tMJD );
                        }
                    }
                }
            }
            fclose( fp );
            done = TRUE;
            success3 = TRUE;
        } else {
            // only complain if this is our last try.
            if ( Path == NULL ) { //i.e. QIN_DENTON_PATH environment variable not set.
                printf( "Cannot open %s file. Try setting QIN_DENTON_PATH environment variable to path containing QinDenton files.\n", Filename );
            } else {
                printf( "Cannot open %s file.\n", Filename );
            }
        }

        if (j == 1) {
            done = TRUE;
        } else {
            ++j;
        }

    }


    q->nPnts = n;






    free(Filename);
    free(Line);
    Lgm_free_ctrans( c );

}


/** 
 *   \brief
 *      This routine interpolates the loaded Qin-Denton parameters to a given
 *      Julian day and poulates an Lgm_QinDentonOne structure
 *
 *   \details
 *      This routine has a number of return values, which are linear sums of:
 *          1  - No data available at this time; default values returned.
 *          2  - Kp data has been filled with default values.
 *          4  - Dst data has been filled with default values.
 *          8  - IMF data has been filled with default values.
 *          16 - Velocity data has been filled with default values.
 *          32 - Density/Pressure data has been filled with default values.
 *
 *   \param[in]      JD           The date/time represented as a Julian date.
 *   \param[in,out]  p            Pointer to an Lgm_QinDentonOne structure.
 *   \param[in]      verbose      A flag to set the verbosity level.
 *   \param[in]      Persistence  A flag to allow persistence population of empty fields.
 *
 *
 *   \returns        int
 *
 */
int Lgm_get_QinDenton_at_JD( double JD, Lgm_QinDentonOne *p, int Verbose, int Persistence ) {

    int                 t, nq, i, ny, nm, nd, nGood;
    double              *x, *y, MJD, UTC;
    gsl_interp_accel    *acc;
    gsl_spline          *spline;
    Lgm_QinDenton       *q = Lgm_init_QinDenton( Verbose );
    int                 UsePersistence;
    double              Wdefaults[10] = {0.44, 0.42, 0.66, 0.48, 0.49, 0.91}; //Avg values at status flag 2 (Table 3, Qin et al.)
    double              Gdefaults[10] = {6.0, 10.0, 60.0};
    int                 retval = 0;

    UsePersistence = FALSE;
    MJD = JD - 2400000.5;
    p->JD    = JD;
    p->MJD   = MJD;
    p->Date  = Lgm_JD_to_Date( JD, &p->Year, &p->Month, &p->Day, &UTC );
    p->UTC   = UTC;
    Lgm_UT_to_HMS( UTC, &p->Hour, &p->Minute, &p->Second );
    p->Persistence = Persistence;

    Lgm_read_QinDenton( p->Date, q );
    p->nPnts = q->nPnts;


    if (q->nPnts < 2) {
        retval += 1;

        printf("Not enough QinDenton values to interpolate\n");
        p->Dst   =   -5.0;  // nT
        p->fKp   =    2.0;  // dimensionless
        p->V_SW  =  400.0;  // km/s
        p->Den_P =    5.0;  // #/cm^-3
        p->Pdyn  =    p->Den_P * 1e6 * LGM_PROTON_MASS * p->V_SW*p->V_SW*1e6 * 1e9;   // nPa
        p->ByIMF =    5.0; // nT 
        p->BzIMF =   -5.0; // nT 
        p->G1    =    Gdefaults[0]; // units?
        p->G2    =    Gdefaults[1]; // units?
        p->G3    =    Gdefaults[2]; // units?
        p->akp3  =    2.0; // unitless
        p->W1    =    Wdefaults[0]; // units?
        p->W2    =    Wdefaults[1]; // units?
        p->W3    =    Wdefaults[2]; // units?
        p->W4    =    Wdefaults[3]; // units?
        p->W5    =    Wdefaults[4]; // units?
        p->W6    =    Wdefaults[5]; // units?
        p->Bz1   =    0.0; // nT
        p->Bz2   =    0.0; // nT
        p->Bz3   =    0.0; // nT
        p->Bz4   =    0.0; // nT
        p->Bz5   =    0.0; // nT
        p->Bz6   =    0.0; // nT

    } else if ( (MJD < q->MJD[0]) || (MJD > q->MJD[q->nPnts-1]) ){
        retval += 1;

        if (p->Persistence == 1) {
            UsePersistence = TRUE;
            printf("No Qin Denton data in range -- using persistence. Data MJD range: [%lf, %lf], requested MJD: %lf\n", q->MJD[0], q->MJD[q->nPnts-1], MJD);
            p->Dst   =    q->Dst[q->nPnts-1];  // km/s
            p->fKp   =    q->fKp[q->nPnts-1];  // km/s
            p->V_SW  =    q->V_SW[q->nPnts-1];  // km/s
            p->Den_P =    q->Den_P[q->nPnts-1];  // #/cm^-3
            p->Pdyn  =    q->Pdyn[q->nPnts-1];   // nPa
            p->ByIMF =    q->ByIMF[q->nPnts-1]; // nT 
            p->BzIMF =    q->BzIMF[q->nPnts-1]; // nT 
            p->G1    =    q->G1[q->nPnts-1]; // units? what's a good nominal value here?
            p->G2    =    q->G2[q->nPnts-1]; // units? what's a good nominal value here?
            p->G3    =    q->G2[q->nPnts-1]; // units? what's a good nominal value here?
            p->akp3  =    q->akp3[q->nPnts-1]; // unitless
            p->W1    =    q->W1[q->nPnts-1]; // units? what's a good nominal value here?
            p->W2    =    q->W2[q->nPnts-1]; // units? what's a good nominal value here?
            p->W3    =    q->W3[q->nPnts-1]; // units? what's a good nominal value here?
            p->W4    =    q->W4[q->nPnts-1]; // units? what's a good nominal value here?
            p->W5    =    q->W5[q->nPnts-1]; // units? what's a good nominal value here?
            p->W6    =    q->W6[q->nPnts-1]; // units? what's a good nominal value here?
            p->Bz1   =    q->Bz1[q->nPnts-1]; // nT
            p->Bz2   =    q->Bz2[q->nPnts-1]; // nT
            p->Bz3   =    q->Bz3[q->nPnts-1]; // nT
            p->Bz4   =    q->Bz4[q->nPnts-1]; // nT
            p->Bz5   =    q->Bz5[q->nPnts-1]; // nT
            p->Bz6   =    q->Bz6[q->nPnts-1]; // nT

        } else {
            printf("No Qin Denton data in range -- would require extrapolation. Setting defaults. Data MJD range: [%lf, %lf], requested MJD: %lf\n", q->MJD[0], q->MJD[q->nPnts-1], MJD);
            p->Dst = -5.0;
            p->fKp = 2.0;
            p->V_SW  =  400.0;  // km/s
            p->Den_P =    5.0;  // #/cm^-3
            p->Pdyn  =    p->Den_P * 1e6 * LGM_PROTON_MASS * p->V_SW*p->V_SW*1e6 * 1e9;   // nPa
            p->ByIMF =    2.0; // nT 
            p->BzIMF =   -2.0; // nT 
            p->G1    =    Gdefaults[0]; // units?
            p->G2    =    Gdefaults[1]; // units?
            p->G3    =    Gdefaults[2]; // units?
            p->akp3  =    2.0; // unitless
            p->W1    =    Wdefaults[0]; // units?
            p->W2    =    Wdefaults[1]; // units?
            p->W3    =    Wdefaults[2]; // units?
            p->W4    =    Wdefaults[3]; // units?
            p->W5    =    Wdefaults[4]; // units?
            p->W6    =    Wdefaults[5]; // units?
            p->Bz1   =    0.0; // nT
            p->Bz2   =    0.0; // nT
            p->Bz3   =    0.0; // nT
            p->Bz4   =    0.0; // nT
            p->Bz5   =    0.0; // nT
            p->Bz6   =    0.0; // nT
        }

    } else {

        nq = q->nPnts;
        x  = (double *)calloc( nq, sizeof(double) );
        y  = (double *)calloc( nq, sizeof(double) );

        acc    = gsl_interp_accel_alloc( );

        // interpolate ByIMF
        for ( nGood=0, i=0; i<nq; i++ ){
            if ( (q->ByIMF[i] > -100.0) && (q->ByIMF[i] < 100.0) ) {
                x[nGood] = q->MJD[i];
                y[nGood] = q->ByIMF[i];
                ++nGood;
            }
        }
        if ( nGood > 2 ) {
            spline = ( nGood < 5 ) ? gsl_spline_alloc( gsl_interp_linear, nGood ) : gsl_spline_alloc( gsl_interp_akima, nGood );
            gsl_spline_init( spline, x, y, nGood ); 
            p->ByIMF = gsl_spline_eval( spline, MJD, acc );
            gsl_spline_free (spline);
        } else {
            p->ByIMF = 0.0;
            printf("No Good Qin Denton data in range for ByIMF. Setting ByIMF to %g. Data MJD range: [%lf, %lf], requested MJD: %lf\n", p->ByIMF, q->MJD[0], q->MJD[q->nPnts-1], MJD);
        }

        // interpolate BzIMF
        for ( nGood=0, i=0; i<nq; i++ ){
            if ( (q->BzIMF[i] > -100.0) && (q->BzIMF[i] < 100.0) ) {
                x[nGood] = q->MJD[i];
                y[nGood] = q->BzIMF[i];
                ++nGood;
            }
        }
        if ( nGood > 2 ) {
            spline = ( nGood < 5 ) ? gsl_spline_alloc( gsl_interp_linear, nGood ) : gsl_spline_alloc( gsl_interp_akima, nGood );
            gsl_spline_init( spline, x, y, nGood ); 
            p->BzIMF = gsl_spline_eval( spline, MJD, acc );
            gsl_spline_free (spline);
        } else {
            p->BzIMF = 0.0;
            printf("No Good Qin Denton data in range for BzIMF. Setting BzIMF to %g. Data MJD range: [%lf, %lf], requested MJD: %lf\n", p->BzIMF, q->MJD[0], q->MJD[q->nPnts-1], MJD);
            retval += 8;
        }

        // interpolate V_SW
        for ( nGood=0, i=0; i<nq; i++ ){
            if ( (q->V_SW[i] > 100.0) && (q->V_SW[i] < 2000.0) ) {
                x[nGood] = q->MJD[i];
                y[nGood] = q->V_SW[i];
                ++nGood;
            }
        }
        if ( nGood > 2 ) {
            spline = ( nGood < 5 ) ? gsl_spline_alloc( gsl_interp_linear, nGood ) : gsl_spline_alloc( gsl_interp_akima, nGood );
            gsl_spline_init( spline, x, y, nGood ); 
            p->V_SW = gsl_spline_eval( spline, MJD, acc );
            gsl_spline_free (spline);
        } else {
            p->V_SW = 400.0;
            printf("No Good Qin Denton data in range for V_SW. Setting V_SW to %g. Data MJD range: [%lf, %lf], requested MJD: %lf\n", p->V_SW, q->MJD[0], q->MJD[q->nPnts-1], MJD);
            retval += 16;
        }


        /*
         * Interpolate Den_P. We really should expect (or trust) densities less than about 0.2 cm^-3.
         * And values less than 1 cm^-3 should be suspect....
         */
        for ( nGood=0, i=0; i<nq; i++ ){
            if ( q->Den_P[i] > 0.1 ) {
                x[nGood] = q->MJD[i];
                y[nGood] = q->Den_P[i];
                ++nGood;
            }
        }
        if ( nGood >= 2 ) {
            spline = gsl_spline_alloc( gsl_interp_linear, nGood );
            gsl_spline_init( spline, x, y, nGood ); 
            p->Den_P = gsl_spline_eval( spline, MJD, acc );
            gsl_spline_free( spline );
        } else {
            p->Den_P = 1.0;
            printf("No Good Qin Denton data in range for Den_P. Setting Den_P to %g. Data MJD range: [%lf, %lf], requested MJD: %lf\n", p->Den_P, q->MJD[0], q->MJD[q->nPnts-1], MJD);
            retval += 32;
        }


        /*
         * Interpolate Pdyn. Since Pdyn is basned one Den_P, use Den_P as a filter on what's good.
         */
        for ( nGood=0, i=0; i<nq; i++ ){
            if ( (q->Den_P[i] > 0.1) && (q->Pdyn[i] > 0.0) ) {
                x[nGood] = q->MJD[i];
                y[nGood] = q->Pdyn[i];
                ++nGood;
            }
        }
        if ( nGood >= 2 ) {
            spline = gsl_spline_alloc( gsl_interp_linear, nGood );
            gsl_spline_init( spline, x, y, nGood ); 
            p->Pdyn = gsl_spline_eval( spline, MJD, acc );
            gsl_spline_free( spline );
        } else {
            p->Pdyn = 2.3;
            printf("No Good Qin Denton data in range for Pdyn. Setting Pdyn to %g. Data MJD range: [%lf, %lf], requested MJD: %lf\n", p->Pdyn, q->MJD[0], q->MJD[q->nPnts-1], MJD);
        }


        // interpolate G1
        for ( nGood=0, i=0; i<nq; i++ ){
            x[nGood] = q->MJD[i]; 
            y[nGood] = q->G1[i]; 
            ++nGood;
        }
        if ( nGood > 2 ) {
            spline = ( nGood < 5 ) ? gsl_spline_alloc( gsl_interp_linear, nGood ) : gsl_spline_alloc( gsl_interp_akima, nGood );
            gsl_spline_init( spline, x, y, nGood ); 
            p->G1 = gsl_spline_eval( spline, MJD, acc );
            gsl_spline_free (spline);
        } else {
            p->G1 = 0.0;
            printf("No Good Qin Denton data in range for G1. Setting G1 to %g. Data MJD range: [%lf, %lf], requested MJD: %lf\n", p->G1, q->MJD[0], q->MJD[q->nPnts-1], MJD);
        }

        // interpolate G2
        for ( nGood=0, i=0; i<nq; i++ ){
            x[nGood] = q->MJD[i]; 
            y[nGood] = q->G2[i]; 
            ++nGood;
        }
        if ( nGood > 2 ) {
            spline = ( nGood < 5 ) ? gsl_spline_alloc( gsl_interp_linear, nGood ) : gsl_spline_alloc( gsl_interp_akima, nGood );
            gsl_spline_init( spline, x, y, nGood ); 
            p->G2 = gsl_spline_eval( spline, MJD, acc );
            gsl_spline_free (spline);
        } else {
            p->G2 = 0.0;
            printf("No Good Qin Denton data in range for G2. Setting G2 to %g. Data MJD range: [%lf, %lf], requested MJD: %lf\n", p->G2, q->MJD[0], q->MJD[q->nPnts-1], MJD);
        }

        // interpolate G3
        for ( nGood=0, i=0; i<nq; i++ ){
            x[nGood] = q->MJD[i]; 
            y[nGood] = q->G3[i]; 
            ++nGood;
        }
        if ( nGood > 2 ) {
            spline = ( nGood < 5 ) ? gsl_spline_alloc( gsl_interp_linear, nGood ) : gsl_spline_alloc( gsl_interp_akima, nGood );
            gsl_spline_init( spline, x, y, nGood ); 
            p->G3 = gsl_spline_eval( spline, MJD, acc );
            gsl_spline_free (spline);
        } else {
            p->G3 = 0.0;
            printf("No Good Qin Denton data in range for G3. Setting G3 to %g. Data MJD range: [%lf, %lf], requested MJD: %lf\n", p->G3, q->MJD[0], q->MJD[q->nPnts-1], MJD);
        }

        // interpolate Kp
        for ( nGood=0, i=0; i<nq; i++ ){
            if ( (q->fKp[i] >= 0.0) && (q->fKp[i] <= 9.0)){
                x[nGood] = q->MJD[i]; 
                y[nGood] = q->fKp[i]; 
                ++nGood;
            }
        }
        if ( nGood > 2 ) {
            spline = ( nGood < 5 ) ? gsl_spline_alloc( gsl_interp_linear, nGood ) : gsl_spline_alloc( gsl_interp_akima, nGood );
            gsl_spline_init( spline, x, y, nGood ); 
            p->fKp = gsl_spline_eval( spline, MJD, acc );
            gsl_spline_free (spline);
        } else {
            p->fKp = 2.0;
            printf("No Good Qin Denton data in range for Kp. Setting Kp to %g. Data MJD range: [%lf, %lf], requested MJD: %lf\n", p->fKp, q->MJD[0], q->MJD[q->nPnts-1], MJD);
            retval += 2;
        }

        // interpolate akp3
        for ( nGood=0, i=0; i<nq; i++ ){
            if ( (q->akp3[i] >= 0.0) && (q->akp3[i] <= 9.0)){
                x[nGood] = q->MJD[i]; 
                y[nGood] = q->akp3[i]; 
                ++nGood;
            }
        }
        if ( nGood > 2 ) {
            spline = ( nGood < 5 ) ? gsl_spline_alloc( gsl_interp_linear, nGood ) : gsl_spline_alloc( gsl_interp_akima, nGood );
            gsl_spline_init( spline, x, y, nGood );
            p->akp3 = gsl_spline_eval( spline, MJD, acc );
            gsl_spline_free (spline);
        } else {
            p->akp3 = 2.0;
            printf("No Good Qin Denton data in range for akp3. Setting akp3 to %g. Data MJD range: [%lf, %lf], requested MJD: %lf\n", p->akp3, q->MJD[0], q->MJD[q->nPnts-1], MJD);
        }

        // interpolate Dst
        for ( nGood=0, i=0; i<nq; i++ ){
            if ( (q->Dst[i] >= -2500.0) && (q->Dst[i] < 500.0)){
                x[nGood] = q->MJD[i]; 
                y[nGood] = q->Dst[i]; 
                ++nGood;
            }
        }
        if ( nGood > 2 ) {
            spline = ( nGood < 5 ) ? gsl_spline_alloc( gsl_interp_linear, nGood ) : gsl_spline_alloc( gsl_interp_akima, nGood );
            gsl_spline_init( spline, x, y, nGood );
            p->Dst = gsl_spline_eval( spline, MJD, acc );
            gsl_spline_free (spline);
        } else {
            p->Dst = 0.0;
            printf("No Good Qin Denton data in range for Dst. Setting Dst to %g. Data MJD range: [%lf, %lf], requested MJD: %lf\n", p->Dst, q->MJD[0], q->MJD[q->nPnts-1], MJD);
            retval += 4;
        }

        // interpolate Bz1 -- are these even used?
        for ( nGood=0, i=0; i<nq; i++ ){
            x[nGood] = q->MJD[i]; 
            y[nGood] = q->Bz1[i]; 
            ++nGood;
        }
        if ( nGood > 2 ) {
            spline = ( nGood < 5 ) ? gsl_spline_alloc( gsl_interp_linear, nGood ) : gsl_spline_alloc( gsl_interp_akima, nGood );
            gsl_spline_init( spline, x, y, nGood );
            p->Bz1 = gsl_spline_eval( spline, MJD, acc );
            gsl_spline_free (spline);
        } else {
            p->Bz1 = 0.0;
            printf("No Good Qin Denton data in range for Bz1. Setting Bz1 to %g. Data MJD range: [%lf, %lf], requested MJD: %lf\n", p->Bz1, q->MJD[0], q->MJD[q->nPnts-1], MJD);
        }

        // interpolate Bz2 -- are these even used?
        for ( nGood=0, i=0; i<nq; i++ ){
            x[nGood] = q->MJD[i]; 
            y[nGood] = q->Bz2[i]; 
            ++nGood;
        }
        if ( nGood > 2 ) {
            spline = ( nGood < 5 ) ? gsl_spline_alloc( gsl_interp_linear, nGood ) : gsl_spline_alloc( gsl_interp_akima, nGood );
            gsl_spline_init( spline, x, y, nGood );
            p->Bz2 = gsl_spline_eval( spline, MJD, acc );
            gsl_spline_free (spline);
        } else {
            p->Bz2 = 0.0;
            printf("No Good Qin Denton data in range for Bz2. Setting Bz2 to %g. Data MJD range: [%lf, %lf], requested MJD: %lf\n", p->Bz2, q->MJD[0], q->MJD[q->nPnts-1], MJD);
        }

        // interpolate Bz3 -- are these even used?
        for ( nGood=0, i=0; i<nq; i++ ){
            x[nGood] = q->MJD[i]; 
            y[nGood] = q->Bz3[i]; 
            ++nGood;
        }
        if ( nGood > 2 ) {
            spline = ( nGood < 5 ) ? gsl_spline_alloc( gsl_interp_linear, nGood ) : gsl_spline_alloc( gsl_interp_akima, nGood );
            gsl_spline_init( spline, x, y, nGood );
            p->Bz3 = gsl_spline_eval( spline, MJD, acc );
            gsl_spline_free (spline);
        } else {
            p->Bz3 = 0.0;
            printf("No Good Qin Denton data in range for Bz3. Setting Bz3 to %g. Data MJD range: [%lf, %lf], requested MJD: %lf\n", p->Bz3, q->MJD[0], q->MJD[q->nPnts-1], MJD);
        }

        // interpolate Bz4 -- are these even used?
        for ( nGood=0, i=0; i<nq; i++ ){
            x[nGood] = q->MJD[i]; 
            y[nGood] = q->Bz4[i]; 
            ++nGood;
        }
        if ( nGood > 2 ) {
            spline = ( nGood < 5 ) ? gsl_spline_alloc( gsl_interp_linear, nGood ) : gsl_spline_alloc( gsl_interp_akima, nGood );
            gsl_spline_init( spline, x, y, nGood );
            p->Bz4 = gsl_spline_eval( spline, MJD, acc );
            gsl_spline_free (spline);
        } else {
            p->Bz4 = 0.0;
            printf("No Good Qin Denton data in range for Bz4. Setting Bz4 to %g. Data MJD range: [%lf, %lf], requested MJD: %lf\n", p->Bz4, q->MJD[0], q->MJD[q->nPnts-1], MJD);
        }

        // interpolate Bz5 -- are these even used?
        for ( nGood=0, i=0; i<nq; i++ ){
            x[nGood] = q->MJD[i]; 
            y[nGood] = q->Bz5[i]; 
            ++nGood;
        }
        if ( nGood > 2 ) {
            spline = ( nGood < 5 ) ? gsl_spline_alloc( gsl_interp_linear, nGood ) : gsl_spline_alloc( gsl_interp_akima, nGood );
            gsl_spline_init( spline, x, y, nGood );
            p->Bz5 = gsl_spline_eval( spline, MJD, acc );
            gsl_spline_free (spline);
        } else {
            p->Bz5 = 0.0;
            printf("No Good Qin Denton data in range for Bz5. Setting Bz5 to %g. Data MJD range: [%lf, %lf], requested MJD: %lf\n", p->Bz5, q->MJD[0], q->MJD[q->nPnts-1], MJD);
        }

        // interpolate Bz6 -- are these even used?
        for ( nGood=0, i=0; i<nq; i++ ){
            x[nGood] = q->MJD[i]; 
            y[nGood] = q->Bz6[i]; 
            ++nGood;
        }
        if ( nGood > 2 ) {
            spline = ( nGood < 5 ) ? gsl_spline_alloc( gsl_interp_linear, nGood ) : gsl_spline_alloc( gsl_interp_akima, nGood );
            gsl_spline_init( spline, x, y, nGood );
            p->Bz6 = gsl_spline_eval( spline, MJD, acc );
            gsl_spline_free (spline);
        } else {
            p->Bz6 = 0.0;
            printf("No Good Qin Denton data in range for Bz6. Setting Bz6 to %g. Data MJD range: [%lf, %lf], requested MJD: %lf\n", p->Bz6, q->MJD[0], q->MJD[q->nPnts-1], MJD);
        }

        // interpolate W1
        for ( nGood=0, i=0; i<nq; i++ ){
            x[nGood] = q->MJD[i]; 
            y[nGood] = q->W1[i]; 
            ++nGood;
        }
        if ( nGood > 2 ) {
            spline = ( nGood < 5 ) ? gsl_spline_alloc( gsl_interp_linear, nGood ) : gsl_spline_alloc( gsl_interp_akima, nGood );
            gsl_spline_init( spline, x, y, nGood );
            p->W1 = gsl_spline_eval( spline, MJD, acc );
            gsl_spline_free (spline);
        } else {
            p->W1 = 0.0;
            printf("No Good Qin Denton data in range for W1. Setting W1 to %g. Data MJD range: [%lf, %lf], requested MJD: %lf\n", p->W1, q->MJD[0], q->MJD[q->nPnts-1], MJD);
        }

        // interpolate W2
        for ( nGood=0, i=0; i<nq; i++ ){
            x[nGood] = q->MJD[i]; 
            y[nGood] = q->W2[i]; 
            ++nGood;
        }
        if ( nGood > 2 ) {
            spline = ( nGood < 5 ) ? gsl_spline_alloc( gsl_interp_linear, nGood ) : gsl_spline_alloc( gsl_interp_akima, nGood );
            gsl_spline_init( spline, x, y, nGood );
            p->W2 = gsl_spline_eval( spline, MJD, acc );
            gsl_spline_free (spline);
        } else {
            p->W2 = 0.0;
            printf("No Good Qin Denton data in range for W2. Setting W2 to %g. Data MJD range: [%lf, %lf], requested MJD: %lf\n", p->W2, q->MJD[0], q->MJD[q->nPnts-1], MJD);
        }

        // interpolate W3
        for ( nGood=0, i=0; i<nq; i++ ){
            x[nGood] = q->MJD[i]; 
            y[nGood] = q->W3[i]; 
            ++nGood;
        }
        if ( nGood > 2 ) {
            spline = ( nGood < 5 ) ? gsl_spline_alloc( gsl_interp_linear, nGood ) : gsl_spline_alloc( gsl_interp_akima, nGood );
            gsl_spline_init( spline, x, y, nGood );
            p->W3 = gsl_spline_eval( spline, MJD, acc );
            gsl_spline_free (spline);
        } else {
            p->W3 = 0.0;
            printf("No Good Qin Denton data in range for W3. Setting W3 to %g. Data MJD range: [%lf, %lf], requested MJD: %lf\n", p->W3, q->MJD[0], q->MJD[q->nPnts-1], MJD);
        }

        // interpolate W4
        for ( nGood=0, i=0; i<nq; i++ ){
            x[nGood] = q->MJD[i]; 
            y[nGood] = q->W4[i]; 
            ++nGood;
        }
        if ( nGood > 2 ) {
            spline = ( nGood < 5 ) ? gsl_spline_alloc( gsl_interp_linear, nGood ) : gsl_spline_alloc( gsl_interp_akima, nGood );
            gsl_spline_init( spline, x, y, nGood );
            p->W4 = gsl_spline_eval( spline, MJD, acc );
            gsl_spline_free (spline);
        } else {
            p->W4 = 0.0;
            printf("No Good Qin Denton data in range for W4. Setting W4 to %g. Data MJD range: [%lf, %lf], requested MJD: %lf\n", p->W4, q->MJD[0], q->MJD[q->nPnts-1], MJD);
        }

        // interpolate W5
        for ( nGood=0, i=0; i<nq; i++ ){
            x[nGood] = q->MJD[i]; 
            y[nGood] = q->W5[i]; 
            ++nGood;
        }
        if ( nGood > 2 ) {
            spline = ( nGood < 5 ) ? gsl_spline_alloc( gsl_interp_linear, nGood ) : gsl_spline_alloc( gsl_interp_akima, nGood );
            gsl_spline_init( spline, x, y, nGood );
            p->W5 = gsl_spline_eval( spline, MJD, acc );
            gsl_spline_free (spline);
        } else {
            p->W5 = 0.0;
            printf("No Good Qin Denton data in range for W5. Setting W5 to %g. Data MJD range: [%lf, %lf], requested MJD: %lf\n", p->W5, q->MJD[0], q->MJD[q->nPnts-1], MJD);
        }

        // interpolate W6
        for ( nGood=0, i=0; i<nq; i++ ){
            x[nGood] = q->MJD[i]; 
            y[nGood] = q->W6[i]; 
            ++nGood;
        }
        if ( nGood > 2 ) {
            spline = ( nGood < 5 ) ? gsl_spline_alloc( gsl_interp_linear, nGood ) : gsl_spline_alloc( gsl_interp_akima, nGood );
            gsl_spline_init( spline, x, y, nGood );
            p->W6 = gsl_spline_eval( spline, MJD, acc );
            gsl_spline_free (spline);
        } else {
            p->W6 = 0.0;
            printf("No Good Qin Denton data in range for W6. Setting W6 to %g. Data MJD range: [%lf, %lf], requested MJD: %lf\n", p->W6, q->MJD[0], q->MJD[q->nPnts-1], MJD);
        }

        free(x);
        free(y);
        gsl_interp_accel_free (acc);

    }

    if ( q->Verbosity > 0 ) {
        printf("\n");
        if (!UsePersistence) {
            printf("\t\t         QinDenton Parameters\n");
            printf("\t\t    --------------------------------\n");
        } else {
            printf("\t\t         QinDenton Parameters (Persistence)\n");
            printf("\t\t    --------------------------------------------\n");
        }
        printf("\t\t        Date = %8ld\n", p->Date );
        printf("\t\t          JD = %11.7lf\n", p->JD );
        printf("\t\t         MJD = %11.7lf\n", p->MJD );
        printf("\t\t         UTC = %11.7lf  ( ", p->UTC );     Lgm_Print_HMSd( UTC ); printf(" )\n");
        printf("\t\t       ByIMF = %11.7g nT\n", p->ByIMF );
        printf("\t\t       BzIMF = %11.7g nT\n", p->BzIMF );
        printf("\t\t        V_SW = %11.7g km/s\n", p->V_SW );
        printf("\t\t       Den_P = %11.7g #/cm^3\n", p->Den_P );
        printf("\t\t        Pdyn = %11.7g nPa\n", p->Pdyn );
        printf("\t\t          G1 = %11.7g\n", p->G1 );
        printf("\t\t          G2 = %11.7g\n", p->G2 );
        printf("\t\t          G3 = %11.7g\n", p->G3 );
        printf("\t\t          Kp = %11.7g\n", p->fKp );
        printf("\t\t        akp3 = %11.7g\n", p->akp3 );
        printf("\t\t         Dst = %11.7g nT\n", p->Dst );
        printf("\t\t         Bz1 = %11.7g nT\n", p->Bz1 );
        printf("\t\t         Bz2 = %11.7g nT\n", p->Bz2 );
        printf("\t\t         Bz3 = %11.7g nT\n", p->Bz3 );
        printf("\t\t         Bz4 = %11.7g nT\n", p->Bz4 );
        printf("\t\t         Bz5 = %11.7g nT\n", p->Bz5 );
        printf("\t\t         Bz6 = %11.7g nT\n", p->Bz6 );
        printf("\t\t          W1 = %11.7g\n", p->W1 );
        printf("\t\t          W2 = %11.7g\n", p->W2 );
        printf("\t\t          W3 = %11.7g\n", p->W3 );
        printf("\t\t          W4 = %11.7g\n", p->W4 );
        printf("\t\t          W5 = %11.7g\n", p->W5 );
        printf("\t\t          W6 = %11.7g\n", p->W6 );

        printf("\n");
    }


    Lgm_destroy_QinDenton( q );

    return retval;
}

void Lgm_set_QinDenton( Lgm_QinDentonOne *p, Lgm_MagModelInfo *m ) {
    m->Bx   = 0.0;
    m->By   = p->ByIMF;
    m->Bz   = p->BzIMF;
    m->V    = p->V_SW;
    m->Den  = p->Den_P;
    m->P    = p->Pdyn;
    m->G1   = p->G1;
    m->G2   = p->G2;
    m->G3   = p->G3;
    m->fKp  = p->fKp;
    m->Kp   = (int)(p->fKp+0.5); if (m->Kp > 6) m->Kp = 6; if (m->Kp < 0 ) m->Kp = 0;
    m->aKp3 = p->akp3;
    m->Dst  = p->Dst;
    m->W[0] = p->W1;
    m->W[1] = p->W2;
    m->W[2] = p->W3;
    m->W[3] = p->W4;
    m->W[4] = p->W5;
    m->W[5] = p->W6;
}

