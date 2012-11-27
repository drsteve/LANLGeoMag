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


    Path = getenv( "QIN_DENTON_PATH" );
    if ( Path == NULL ) {
        strcpy( QinDentonPath, LGM_INDEX_DATA_DIR );
    } else {

        /*
         * Test for existence
         */
        struct stat sts;
        if ( ( stat( Path, &sts ) ) == -1 ) {
            printf("Environment variable QIN_DENTON_PATH points to a non-existent directory: %s\n", Path );
            strcpy( QinDentonPath, LGM_INDEX_DATA_DIR );
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
    j = 0; done = FALSE; success1 = FALSE;
    while ( !done ){

        sprintf( Filename, "%s/QinDenton/%4d/QinDenton_%8ld_%s.txt", QinDentonPath, Prev_Year, Prev_Date, ftype[j] );
        if ( (fp = fopen( Filename, "r" )) != NULL ) {
            while( fgets( Line, 2048, fp ) != NULL ) {
                if ( Line[0] != '#' ) {
                    sscanf( Line, "%s %d %d %d %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %d %d %d %d %d %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d %d %d %d %d",
                                q->IsoTimeStr[n], &q->Year[n], &q->Month[n], &q->Day[n], &q->Hour[n], &q->Minute[n], &q->Second[n],
                                &q->ByIMF[n], &q->BzIMF[n], &q->V_SW[n], &q->Den_P[n], &q->Pdyn[n],
                                &q->G1[n], &q->G2[n], &q->G3[n],
                                &q->ByIMF_status[n], &q->BzIMF_status[n], &q->V_SW_status[n], &q->Den_P_status[n], &q->Pdyn_status[n],
                                &q->G1_status[n], &q->G2_status[n], &q->G3_status[n],
                                &q->fKp[n], &q->akp3[n], &q->Dst[n],
                                &q->Bz1[n], &q->Bz2[n], &q->Bz3[n], &q->Bz4[n], &q->Bz5[n], &q->Bz6[n],
                                &q->W1[n], &q->W2[n], &q->W3[n], &q->W4[n], &q->W5[n], &q->W6[n],
                                &q->W1_status[n], &q->W2_status[n], &q->W3_status[n], &q->W4_status[n], &q->W5_status[n], &q->W6_status[n] );
                    Time = q->Hour[n] + q->Minute[n]/60.0 + q->Second[n]/3600.0;
                    q->MJD[n] = Lgm_MJD( q->Year[n], q->Month[n], q->Day[n], Time, LGM_TIME_SYS_UTC, c );
                    ++n;
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

        sprintf( Filename, "%s/QinDenton/%4d/QinDenton_%8ld_%s.txt", QinDentonPath, Year, Date, ftype[j] );
        if ( (fp = fopen( Filename, "r" )) != NULL ) {
            while( fgets( Line, 2048, fp ) != NULL ) {
                if ( Line[0] != '#' ) {
                    sscanf( Line, "%s %d %d %d %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %d %d %d %d %d %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d %d %d %d %d",
                                q->IsoTimeStr[n], &q->Year[n], &q->Month[n], &q->Day[n], &q->Hour[n], &q->Minute[n], &q->Second[n],
                                &q->ByIMF[n], &q->BzIMF[n], &q->V_SW[n], &q->Den_P[n], &q->Pdyn[n],
                                &q->G1[n], &q->G2[n], &q->G3[n],
                                &q->ByIMF_status[n], &q->BzIMF_status[n], &q->V_SW_status[n], &q->Den_P_status[n], &q->Pdyn_status[n],
                                &q->G1_status[n], &q->G2_status[n], &q->G3_status[n],
                                &q->fKp[n], &q->akp3[n], &q->Dst[n],
                                &q->Bz1[n], &q->Bz2[n], &q->Bz3[n], &q->Bz4[n], &q->Bz5[n], &q->Bz6[n],
                                &q->W1[n], &q->W2[n], &q->W3[n], &q->W4[n], &q->W5[n], &q->W6[n],
                                &q->W1_status[n], &q->W2_status[n], &q->W3_status[n], &q->W4_status[n], &q->W5_status[n], &q->W6_status[n] );
                    Time = q->Hour[n] + q->Minute[n]/60.0 + q->Second[n]/3600.0;
                    q->MJD[n] = Lgm_MJD( q->Year[n], q->Month[n], q->Day[n], Time, LGM_TIME_SYS_UTC, c );
                    ++n;
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

        sprintf( Filename, "%s/QinDenton/%4d/QinDenton_%8ld_%s.txt", QinDentonPath, Next_Year, Next_Date, ftype[j] );
        if ( (fp = fopen( Filename, "r" )) != NULL ) {
            while( fgets( Line, 2048, fp ) != NULL ) {
                if ( Line[0] != '#' ) {
                    sscanf( Line, "%s %d %d %d %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %d %d %d %d %d %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d %d %d %d %d",
                                q->IsoTimeStr[n], &q->Year[n], &q->Month[n], &q->Day[n], &q->Hour[n], &q->Minute[n], &q->Second[n],
                                &q->ByIMF[n], &q->BzIMF[n], &q->V_SW[n], &q->Den_P[n], &q->Pdyn[n],
                                &q->G1[n], &q->G2[n], &q->G3[n],
                                &q->ByIMF_status[n], &q->BzIMF_status[n], &q->V_SW_status[n], &q->Den_P_status[n], &q->Pdyn_status[n],
                                &q->G1_status[n], &q->G2_status[n], &q->G3_status[n],
                                &q->fKp[n], &q->akp3[n], &q->Dst[n],
                                &q->Bz1[n], &q->Bz2[n], &q->Bz3[n], &q->Bz4[n], &q->Bz5[n], &q->Bz6[n],
                                &q->W1[n], &q->W2[n], &q->W3[n], &q->W4[n], &q->W5[n], &q->W6[n],
                                &q->W1_status[n], &q->W2_status[n], &q->W3_status[n], &q->W4_status[n], &q->W5_status[n], &q->W6_status[n] );
                    Time = q->Hour[n] + q->Minute[n]/60.0 + q->Second[n]/3600.0;
                    q->MJD[n] = Lgm_MJD( q->Year[n], q->Month[n], q->Day[n], Time, LGM_TIME_SYS_UTC, c );
                    ++n;
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


void Lgm_get_QinDenton_at_JD( double JD, Lgm_QinDentonOne *p, int Verbose ) {

    int                 t, nq, i, ny, nm, nd;
    double              *x, *y, MJD, UTC;
    gsl_interp_accel    *acc;
    gsl_spline          *spline;
    Lgm_QinDenton       *q = Lgm_init_QinDenton( Verbose );

    MJD = JD - 2400000.5;
    p->JD   = JD;
    p->MJD  = MJD;
    p->Date = Lgm_JD_to_Date( JD, &p->Year, &p->Month, &p->Day, &UTC );
    p->UTC  = UTC;
    Lgm_UT_to_HMS( UTC, &p->Hour, &p->Minute, &p->Second );

    Lgm_read_QinDenton( p->Date, q );

    if (q->nPnts < 2) {

        printf("Not enough QinDenton values to interpolate\n");
        p->Dst = -5.0;
        p->fKp = 2.0;

    } else if ( (JD < q->MJD[0]) || (JD > q->MJD[q->nPnts-1]) ){
        printf("No Qin Denton data in range -- would require extrapolation.\n");
        p->Dst = -5.0;
        p->fKp = 2.0;
    } else {

        nq = q->nPnts;
        x  = (double *)calloc( nq, sizeof(double) );
        y  = (double *)calloc( nq, sizeof(double) );

        acc    = gsl_interp_accel_alloc( );
        //spline = gsl_spline_alloc( gsl_interp_cspline, nq );
        if ( nq < 5 ) {
            spline = gsl_spline_alloc( gsl_interp_linear, nq );
        } else {
            spline = gsl_spline_alloc( gsl_interp_akima, nq );
        }



        // interpolate ByIMF
        for ( i=0; i<nq; i++ ){ x[i] = q->MJD[i]; y[i] = q->ByIMF[i]; }
        gsl_spline_init( spline, x, y, nq ); p->ByIMF = gsl_spline_eval( spline, MJD, acc );

        // interpolate BzIMF
        for ( i=0; i<nq; i++ ){ x[i] = q->MJD[i]; y[i] = q->BzIMF[i]; }
        gsl_spline_init( spline, x, y, nq ); p->BzIMF = gsl_spline_eval( spline, MJD, acc );

        // interpolate V_SW
        for ( i=0; i<nq; i++ ){ x[i] = q->MJD[i]; y[i] = q->V_SW[i]; }
        gsl_spline_init( spline, x, y, nq ); p->V_SW = gsl_spline_eval( spline, MJD, acc );

        // interpolate Den_P
        for ( i=0; i<nq; i++ ){ x[i] = q->MJD[i]; y[i] = q->Den_P[i]; }
        gsl_spline_init( spline, x, y, nq ); p->Den_P = gsl_spline_eval( spline, MJD, acc );

        // interpolate Pdyn
        for ( i=0; i<nq; i++ ){ x[i] = q->MJD[i]; y[i] = q->Pdyn[i]; }
        gsl_spline_init( spline, x, y, nq ); p->Pdyn = gsl_spline_eval( spline, MJD, acc );

        // interpolate G1
        for ( i=0; i<nq; i++ ){ x[i] = q->MJD[i]; y[i] = q->G1[i]; }
        gsl_spline_init( spline, x, y, nq ); p->G1 = gsl_spline_eval( spline, MJD, acc );

        // interpolate G2
        for ( i=0; i<nq; i++ ){ x[i] = q->MJD[i]; y[i] = q->G2[i]; }
        gsl_spline_init( spline, x, y, nq ); p->G2 = gsl_spline_eval( spline, MJD, acc );

        // interpolate G3
        for ( i=0; i<nq; i++ ){ x[i] = q->MJD[i]; y[i] = q->G3[i]; }
        gsl_spline_init( spline, x, y, nq ); p->G3 = gsl_spline_eval( spline, MJD, acc );

        // interpolate Kp
        for ( i=0; i<nq; i++ ){ x[i] = q->MJD[i]; y[i] = q->fKp[i]; }
        gsl_spline_init( spline, x, y, nq ); p->fKp = gsl_spline_eval( spline, MJD, acc );

        // interpolate akp3
        for ( i=0; i<nq; i++ ){ x[i] = q->MJD[i]; y[i] = q->akp3[i]; }
        gsl_spline_init( spline, x, y, nq ); p->akp3 = gsl_spline_eval( spline, MJD, acc );

        // interpolate Dst
        for ( i=0; i<nq; i++ ){ x[i] = q->MJD[i]; y[i] = q->Dst[i]; }
        gsl_spline_init( spline, x, y, nq ); p->Dst = gsl_spline_eval( spline, MJD, acc );

        // interpolate Bz1
        for ( i=0; i<nq; i++ ){ x[i] = q->MJD[i]; y[i] = q->Bz1[i]; }
        gsl_spline_init( spline, x, y, nq ); p->Bz1 = gsl_spline_eval( spline, MJD, acc );

        // interpolate Bz2
        for ( i=0; i<nq; i++ ){ x[i] = q->MJD[i]; y[i] = q->Bz2[i]; }
        gsl_spline_init( spline, x, y, nq ); p->Bz2 = gsl_spline_eval( spline, MJD, acc );

        // interpolate Bz3
        for ( i=0; i<nq; i++ ){ x[i] = q->MJD[i]; y[i] = q->Bz3[i]; }
        gsl_spline_init( spline, x, y, nq ); p->Bz3 = gsl_spline_eval( spline, MJD, acc );

        // interpolate Bz4
        for ( i=0; i<nq; i++ ){ x[i] = q->MJD[i]; y[i] = q->Bz4[i]; }
        gsl_spline_init( spline, x, y, nq ); p->Bz4 = gsl_spline_eval( spline, MJD, acc );

        // interpolate Bz5
        for ( i=0; i<nq; i++ ){ x[i] = q->MJD[i]; y[i] = q->Bz5[i]; }
        gsl_spline_init( spline, x, y, nq ); p->Bz5 = gsl_spline_eval( spline, MJD, acc );

        // interpolate Bz6
        for ( i=0; i<nq; i++ ){ x[i] = q->MJD[i]; y[i] = q->Bz6[i]; }
        gsl_spline_init( spline, x, y, nq ); p->Bz6 = gsl_spline_eval( spline, MJD, acc );

        // interpolate W1
        for ( i=0; i<nq; i++ ){ x[i] = q->MJD[i]; y[i] = q->W1[i]; }
        gsl_spline_init( spline, x, y, nq ); p->W1 = gsl_spline_eval( spline, MJD, acc );

        // interpolate W2
        for ( i=0; i<nq; i++ ){ x[i] = q->MJD[i]; y[i] = q->W2[i]; }
        gsl_spline_init( spline, x, y, nq ); p->W2 = gsl_spline_eval( spline, MJD, acc );

        // interpolate W3
        for ( i=0; i<nq; i++ ){ x[i] = q->MJD[i]; y[i] = q->W3[i]; }
        gsl_spline_init( spline, x, y, nq ); p->W3 = gsl_spline_eval( spline, MJD, acc );

        // interpolate W4
        for ( i=0; i<nq; i++ ){ x[i] = q->MJD[i]; y[i] = q->W4[i]; }
        gsl_spline_init( spline, x, y, nq ); p->W4 = gsl_spline_eval( spline, MJD, acc );

        // interpolate W5
        for ( i=0; i<nq; i++ ){ x[i] = q->MJD[i]; y[i] = q->W5[i]; }
        gsl_spline_init( spline, x, y, nq ); p->W5 = gsl_spline_eval( spline, MJD, acc );

        // interpolate W6
        for ( i=0; i<nq; i++ ){ x[i] = q->MJD[i]; y[i] = q->W6[i]; }
        gsl_spline_init( spline, x, y, nq ); p->W6 = gsl_spline_eval( spline, MJD, acc );



        free(x);
        free(y);
        gsl_spline_free (spline);
        gsl_interp_accel_free (acc);

    }

    if ( q->Verbosity > 0 ) {
        printf("\n");
        printf("         QinDenton Parameters\n");
        printf("    --------------------------------\n");
        printf("        Date = %8ld\n", p->Date );
        printf("          JD = %11.7lf\n", p->JD );
        printf("         MJD = %11.7lf\n", p->MJD );
        printf("         UTC = %11.7lf  ( ", p->UTC );     Lgm_Print_HMSd( UTC ); printf(" )\n");
        printf("       ByIMF = %11.7g nT\n", p->ByIMF );
        printf("       BzIMF = %11.7g nT\n", p->BzIMF );
        printf("        V_SW = %11.7g km/s\n", p->V_SW );
        printf("       Den_P = %11.7g nPa\n", p->Den_P );
        printf("        Pdyn = %11.7g nPa\n", p->Pdyn );
        printf("          G1 = %11.7g\n", p->G1 );
        printf("          G2 = %11.7g\n", p->G2 );
        printf("          G3 = %11.7g\n", p->G3 );
        printf("          Kp = %11.7g\n", p->fKp );
        printf("        akp3 = %11.7g\n", p->akp3 );
        printf("         Dst = %11.7g nT\n", p->Dst );
        printf("         Bz1 = %11.7g nT\n", p->Bz1 );
        printf("         Bz2 = %11.7g nT\n", p->Bz2 );
        printf("         Bz3 = %11.7g nT\n", p->Bz3 );
        printf("         Bz4 = %11.7g nT\n", p->Bz4 );
        printf("         Bz5 = %11.7g nT\n", p->Bz5 );
        printf("         Bz6 = %11.7g nT\n", p->Bz6 );
        printf("          W1 = %11.7g\n", p->W1 );
        printf("          W2 = %11.7g\n", p->W2 );
        printf("          W3 = %11.7g\n", p->W3 );
        printf("          W4 = %11.7g\n", p->W4 );
        printf("          W5 = %11.7g\n", p->W5 );
        printf("          W6 = %11.7g\n", p->W6 );

        printf("\n");
    }




    Lgm_destroy_QinDenton( q );




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
    m->Kp   = (int)(p->fKp+0.5); if (m->Kp > 5) m->Kp = 5; if (m->Kp < 0 ) m->Kp = 0;
    m->aKp3 = p->akp3;
    m->Dst  = p->Dst;
    m->W[0] = p->W1;
    m->W[1] = p->W2;
    m->W[2] = p->W3;
    m->W[3] = p->W4;
    m->W[4] = p->W5;
    m->W[5] = p->W6;
}

