#ifdef HAVE_CONFIG_H 
#include <config.h> 
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Lgm/Lgm_CTrans.h"
#include "Lgm/Lgm_Eop.h"

#define JD1962 2437665.5


Lgm_Eop *Lgm_init_eop( int Verbose ) {

    Lgm_Eop      *e;
    long int      days;
    double        cJD;
    Lgm_CTrans   *c = Lgm_init_ctrans(0);

    e = (Lgm_Eop *) calloc (1, sizeof(*e));

    e->Verbosity = Verbose;

    /*
     *  Guess roughly how many values we need to calloc. The data goes from
     *  1962->now+1year or so.
     */
    cJD    = Lgm_GetCurrentJD( c );
    days = (long int)(cJD-JD1962+400.5);
    
    e->Size     = days; // size of arrays
    e->nEopVals = 0;    // number of values stored in arrays
    e->Date = (long int *)calloc( days, sizeof(long int) );
    e->MJD  = (double *)calloc( days, sizeof(double) );
    e->xp   = (double *)calloc( days, sizeof(double) );
    e->yp   = (double *)calloc( days, sizeof(double) );
    e->DUT1 = (double *)calloc( days, sizeof(double) );
    e->LOD  = (double *)calloc( days, sizeof(double) );
    e->dPsi = (double *)calloc( days, sizeof(double) );
    e->dEps = (double *)calloc( days, sizeof(double) );
    e->dX   = (double *)calloc( days, sizeof(double) );
    e->dY   = (double *)calloc( days, sizeof(double) );
    e->DAT  = (double *)calloc( days, sizeof(double) );

    Lgm_free_ctrans( c );

    return e;

}

void  Lgm_destroy_eop( Lgm_Eop *e ) {

    // first free arrays inside structure
    free( e->Date );
    free( e->MJD );
    free( e->xp );
    free( e->yp );
    free( e->DUT1 );
    free( e->LOD );
    free( e->dPsi );
    free( e->dEps );
    free( e->dX );
    free( e->dY );
    free( e->DAT );

    // tehn free structure itself
    free( e );

    return;
    
}


void Lgm_read_eop( Lgm_Eop *e ) {

    FILE        *fp;
    long int    n, N;
    char        *Line, *Filename;

    Filename = (char *)calloc( 512, sizeof(char) );
    sprintf( Filename, "%s/LgmEop.dat", LGM_EOP_DATA_DIR );
    if ( (fp = fopen( Filename, "r" )) != NULL ) {

        Line = (char *)calloc( 260, sizeof(char) );
        while( fgets( Line, 256, fp ) != NULL ) {
            if ( Line[0] != '#' ) {
                n = e->nEopVals;
                N = e->Size;
                if ( n == N ) {
                    // Out array sizes need to be increased  (try adding another years worth)
                    N += 365;
                    e->Date = realloc( e->Date, N );
                    e->MJD  = realloc( e->MJD, N );
                    e->xp   = realloc( e->xp, N );
                    e->yp   = realloc( e->yp, N );
                    e->DUT1 = realloc( e->DUT1, N );
                    e->LOD  = realloc( e->LOD, N );
                    e->dPsi = realloc( e->dPsi, N );
                    e->dEps = realloc( e->dEps, N );
                    e->dX   = realloc( e->dX, N );
                    e->dY   = realloc( e->dY, N );
                    e->DAT  = realloc( e->DAT, N );
                    e->Size = N;
                }
                sscanf( Line, "%ld %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &e->Date[n], &e->MJD[n], &e->xp[n],
                    &e->yp[n], &e->DUT1[n], &e->LOD[n], &e->dPsi[n], &e->dEps[n], &e->dX[n], &e->dY[n], &e->DAT[n]);
                ++(e->nEopVals);
            }
        }
        free(Filename);
        free(Line);
        fclose( fp );

    } else {
        printf("Cannot open LgmEop.dat file\n");
    }

}


void Lgm_get_eop_at_JD( double JD, Lgm_EopOne *eop, Lgm_Eop *e ) { 
     
    long int            q, ql, qh, t;
    int                 nq, i, ny, nm, nd;
    double              *x, *y, MJD, UTC;
    gsl_interp_accel    *acc;
    gsl_spline          *spline;

    MJD = JD - 2400000.5;
    eop->JD   = JD;
    eop->MJD  = MJD;
    eop->Date = Lgm_JD_to_Date( JD, &ny, &nm, &nd, &UTC );
    eop->UTC  = UTC;

    if (e->nEopVals < 3) {

        printf("Not enough EOP values to interpolate\n");
        eop->DUT1 = 0.0;
        eop->xp   = 0.0;
        eop->yp   = 0.0;
        eop->dPsi = 0.0;
        eop->dEps = 0.0;
        eop->LOD  = 0.0;
        eop->dX   = 0.0;
        eop->dY   = 0.0;

    } else {

        q = (int)( JD-JD1962 ); // approximate index where JD is found in the LgmEop.dat file
        ql = q-7; // set low index to -7days
        qh = q+7; // set high index to +7days
        if (ql < 0 ) ql = 0;
        if (qh >= e->nEopVals ) qh = e->nEopVals;
        nq = qh-ql+1;
        x = (double *)calloc( nq, sizeof(double) );
        y = (double *)calloc( nq, sizeof(double) );

        acc    = gsl_interp_accel_alloc( );
        spline = gsl_spline_alloc( gsl_interp_cspline, nq );


        

        // interpolate xp
        for ( i=0, t=ql; t<=qh; t++ ){ x[i] = e->MJD[t]; y[i++] = e->xp[t]; }
        gsl_spline_init (spline, x, y, nq); eop->xp = gsl_spline_eval( spline, MJD, acc );

        // interpolate yp
        for ( i=0, t=ql; t<=qh; t++ ){ x[i] = e->MJD[t]; y[i++] = e->yp[t]; }
        gsl_spline_init (spline, x, y, nq); eop->yp = gsl_spline_eval( spline, MJD, acc );

        // interpolate DUT1
        for ( i=0, t=ql; t<=qh; t++ ){ x[i] = e->MJD[t]; y[i++] = e->DUT1[t]; }
        gsl_spline_init (spline, x, y, nq); eop->DUT1 = gsl_spline_eval( spline, MJD, acc );

        // interpolate LOD
        for ( i=0, t=ql; t<=qh; t++ ){ x[i] = e->MJD[t]; y[i++] = e->LOD[t]; }
        gsl_spline_init (spline, x, y, nq); eop->LOD = gsl_spline_eval( spline, MJD, acc );
        
        // interpolate dPsi
        for ( i=0, t=ql; t<=qh; t++ ){ x[i] = e->MJD[t]; y[i++] = e->dPsi[t]; }
        gsl_spline_init (spline, x, y, nq); eop->dPsi = gsl_spline_eval( spline, MJD, acc );
        
        // interpolate dEps
        for ( i=0, t=ql; t<=qh; t++ ){ x[i] = e->MJD[t]; y[i++] = e->dEps[t]; }
        gsl_spline_init (spline, x, y, nq); eop->dEps = gsl_spline_eval( spline, MJD, acc );
        
        // interpolate dX
        for ( i=0, t=ql; t<=qh; t++ ){ x[i] = e->MJD[t]; y[i++] = e->dX[t]; }
        gsl_spline_init (spline, x, y, nq); eop->dX = gsl_spline_eval( spline, MJD, acc );

        // interpolate dY
        for ( i=0, t=ql; t<=qh; t++ ){ x[i] = e->MJD[t]; y[i++] = e->dY[t]; }
        gsl_spline_init (spline, x, y, nq); eop->dY = gsl_spline_eval( spline, MJD, acc );

        // set DAT
        for ( t=ql; t<=qh; t++ ){ 
            if ( MJD <= e->MJD[t] ) {
                eop->DAT = e->DAT[t];
                break;
            }
        }

        // Some of the values may end up being zero if they are in the future
        // and no prediction was given for them. Its probably better to just
        // pin the value at the last real prediction that was made?
        




        free(x);
        free(y);
        gsl_spline_free (spline);
        gsl_interp_accel_free (acc);

    }

    if ( e->Verbosity > 0 ) {
        printf("\n");
        printf("    Earth Orientation Parameters\n");
        printf("    ----------------------------\n");
        printf("      Date   = %8ld\n", eop->Date );
        printf("        JD   = %11.7lf\n", eop->JD );
        printf("       MJD   = %11.7lf\n", eop->MJD );
        printf("       UTC   = %11.7lf  ( ", eop->UTC );     Lgm_Print_HMSd( UTC ); printf(" )\n");
        printf("        xp   = %11.7lf \u2033\n", eop->xp );
        printf("        yp   = %11.7lf \u2033\n", eop->yp );
        printf("        DUT1 = %11.7lf s\n", eop->DUT1 );
        printf("        LOD  = %11.7lf s\n", eop->LOD );
        printf("        dPsi = %11.7lf \u2033\n", eop->dPsi);
        printf("        dEps = %11.7lf \u2033\n", eop->dEps);
        printf("        dX   = %11.7lf \u2033\n", eop->dX );
        printf("        dY   = %11.7lf \u2033\n", eop->dY );
        printf("        DAT  = %11.7lf s\n", eop->DAT );
        printf("\n");
    }



        
    



}

void Lgm_set_eop( Lgm_EopOne *eop, Lgm_CTrans *c ) { 
    c->DUT1 = eop->DUT1;
    c->DAT  = eop->DAT;
    c->xp   = eop->xp;
    c->yp   = eop->yp;
    c->dPsi = eop->dPsi;
    c->dEps = eop->dEps;
    c->LOD  = eop->LOD;
    c->dX   = eop->dX;
    c->dY   = eop->dY;
}

void Lgm_unset_eop( Lgm_EopOne *eop, Lgm_CTrans *c ) { 
    //c->DAT  = eop->DAT; dont unset these ones
    c->DUT1 = 0.0;
    c->xp   = 0.0;
    c->yp   = 0.0;
    c->dPsi = 0.0;
    c->dEps = 0.0;
    c->LOD  = 0.0;
    c->dX   = 0.0;
    c->dY   = 0.0;
}




/*
 * Use the NGA (National Geospatial Intelligence Agency) EP precition
 * polynomials to predict xp, yp, DUT1
 *
 *  JD -- Time to predict quantities for (in JD)
 *
 */
void Lgm_NgaEoppPred( double JD, Lgm_EopOne *eop, Lgm_NgaEopp *e ){


    double  t, tmta, tmtb; 

    t    = JD - 2400000.5; // the t in the model is in MJD
    tmta = t - e->ta;
    tmtb = t - e->tb;


    eop->xp   = e->A + e->B * tmta 
                + e->C1*sin( M_2PI*tmta/e->P1 ) + e->C2*sin( M_2PI*tmta/e->P2 )
                + e->D1*cos( M_2PI*tmta/e->P1 ) + e->D2*cos( M_2PI*tmta/e->P2 );

    eop->yp   = e->E + e->F * tmta 
                + e->G1*sin( M_2PI*tmta/e->Q1 ) + e->G2*sin( M_2PI*tmta/e->Q2 )
                + e->H1*cos( M_2PI*tmta/e->Q1 ) + e->H2*cos( M_2PI*tmta/e->Q2 );

    eop->DUT1 = e->I + e->J * tmtb 
                + e->K1*sin( M_2PI*tmtb/e->R1 ) + e->K2*sin( M_2PI*tmtb/e->R2 )
                + e->K3*sin( M_2PI*tmtb/e->R3 ) + e->K4*sin( M_2PI*tmtb/e->R4 )
                + e->L1*cos( M_2PI*tmtb/e->R1 ) + e->L2*cos( M_2PI*tmtb/e->R2 )
                + e->L3*cos( M_2PI*tmtb/e->R3 ) + e->L4*cos( M_2PI*tmtb/e->R4 );

}


int Lgm_ReadNgaEopp( Lgm_NgaEopp *e, int Verbosity ) {

    int     y, m, d;
    double  UTC;
    char    *Line, Word[11];
    char    *Filename;
    FILE    *fp;

    Filename = (char *)calloc( 512, sizeof(char) );
    sprintf( Filename, "%s/LgmNgaEopp.dat", LGM_EOP_DATA_DIR );
//    strcpy( Filename, "/home/mgh/DREAM/Dream/libLanlGeoMag/EopData/LgmNgaEopp.dat" );

    if ( (fp = fopen( Filename, "r" )) == NULL ) {
        printf( "Coulndt open file %s\n", Filename );
        free( Filename );
        return(0);
    }
    free( Filename );


    Line     = (char *)calloc( 5000, sizeof(char) );

    fgets( Line, 4096, fp );
    strncpy(Word, Line +  0, 10); Word[10] = '\0'; e->ta = atof( Word );
    strncpy(Word, Line + 10, 10); Word[10] = '\0'; e->A  = atof( Word );
    strncpy(Word, Line + 20, 10); Word[10] = '\0'; e->B  = atof( Word );
    strncpy(Word, Line + 30, 10); Word[10] = '\0'; e->C1 = atof( Word );
    strncpy(Word, Line + 40, 10); Word[10] = '\0'; e->C2 = atof( Word );
    strncpy(Word, Line + 50, 10); Word[10] = '\0'; e->D1 = atof( Word );
    strncpy(Word, Line + 60, 10); Word[10] = '\0'; e->D2 = atof( Word );
    strncpy(Word, Line + 70,  6); Word[6]  = '\0'; e->P1 = atof( Word );


    fgets( Line, 4096, fp );
    strncpy(Word, Line +  0,  6); Word[6]  = '\0'; e->P2 = atof( Word );
    strncpy(Word, Line +  6, 10); Word[10] = '\0'; e->E  = atof( Word );
    strncpy(Word, Line + 16, 10); Word[10] = '\0'; e->F  = atof( Word );
    strncpy(Word, Line + 26, 10); Word[10] = '\0'; e->G1 = atof( Word );
    strncpy(Word, Line + 36, 10); Word[10] = '\0'; e->G2 = atof( Word );
    strncpy(Word, Line + 46, 10); Word[10] = '\0'; e->H1 = atof( Word );
    strncpy(Word, Line + 56, 10); Word[10] = '\0'; e->H2 = atof( Word );
    strncpy(Word, Line + 66,  6); Word[6]  = '\0'; e->Q1 = atof( Word );
    strncpy(Word, Line + 72,  6); Word[6]  = '\0'; e->Q2 = atof( Word );


    fgets( Line, 4096, fp );
    strncpy(Word, Line +  0, 10); Word[10] = '\0'; e->tb = atof( Word );
    strncpy(Word, Line + 10, 10); Word[10] = '\0'; e->I  = atof( Word );
    strncpy(Word, Line + 20, 10); Word[10] = '\0'; e->J  = atof( Word );
    strncpy(Word, Line + 30, 10); Word[10] = '\0'; e->K1 = atof( Word );
    strncpy(Word, Line + 40, 10); Word[10] = '\0'; e->K2 = atof( Word );
    strncpy(Word, Line + 50, 10); Word[10] = '\0'; e->K3 = atof( Word );
    strncpy(Word, Line + 60, 10); Word[10] = '\0'; e->K4 = atof( Word );

    
    fgets( Line, 4096, fp );
    strncpy(Word, Line +  0, 10); Word[10] = '\0'; e->L1 = atof( Word );
    strncpy(Word, Line + 10, 10); Word[10] = '\0'; e->L2 = atof( Word );
    strncpy(Word, Line + 20, 10); Word[10] = '\0'; e->L3 = atof( Word );
    strncpy(Word, Line + 30, 10); Word[10] = '\0'; e->L4 = atof( Word );
    strncpy(Word, Line + 40,  9); Word[9]  = '\0'; e->R1 = atof( Word );
    strncpy(Word, Line + 49,  9); Word[9]  = '\0'; e->R2 = atof( Word );
    strncpy(Word, Line + 58,  9); Word[9]  = '\0'; e->R3 = atof( Word );
    strncpy(Word, Line + 67,  9); Word[9]  = '\0'; e->R4 = atof( Word );

    fgets( Line, 4096, fp );
    strncpy(Word, Line +  0, 4); Word[4] = '\0'; e->dat    = atoi( Word );
    strncpy(Word, Line +  4, 5); Word[5] = '\0'; e->EOPPWk = atoi( Word );
    strncpy(Word, Line +  9, 6); Word[6] = '\0'; e->teff   = atoi( Word );

    free( Line );
    fclose( fp );


    /* 
     *  Print out what we read for debugging purposes
     */
    if ( Verbosity > 0 ) {

        printf("    ta:    %10.2lf\n", e->ta);
        Lgm_MJD_to_Date( e->ta, &y, &m, &d, &UTC);
        printf("         %4d %2d %2d %lf\n", y, m, d, UTC);
        printf("     A:    %10.6lf\n", e->A);
        printf("     B:    %10.6lf\n", e->B);
        printf("    C1:    %10.6lf\n", e->C1);
        printf("    C2:    %10.6lf\n", e->C2);
        printf("    D1:    %10.6lf\n", e->D1);
        printf("    D2:    %10.6lf\n", e->D2);
        printf("    P1:    %6.2lf\n",  e->P1);
        printf("    P2:    %6.2lf\n",  e->P2);
        printf("     E:    %10.6lf\n", e->E);
        printf("     F:    %10.6lf\n", e->F);
        printf("    G1:    %10.6lf\n", e->G1);
        printf("    G2:    %10.6lf\n", e->G2);
        printf("    H1:    %10.6lf\n", e->H1);
        printf("    H2:    %10.6lf\n", e->H2);
        printf("    Q1:    %6.2lf\n",  e->Q1);
        printf("    Q2:    %6.2lf\n",  e->Q2);

        printf("    tb:    %10.2lf\n", e->tb);
        Lgm_MJD_to_Date( e->tb, &y, &m, &d, &UTC);
        printf("         %4d %2d %2d %lf\n", y, m, d, UTC);
        printf("     I:    %10.6lf\n", e->I);
        printf("     J:    %10.6lf\n", e->J);
        printf("    K1:    %10.6lf\n", e->K1);
        printf("    K2:    %10.6lf\n", e->K2);
        printf("    K3:    %10.6lf\n", e->K3);
        printf("    K4:    %10.6lf\n", e->K4);
        printf("    L1:    %10.6lf\n", e->L1);
        printf("    L2:    %10.6lf\n", e->L2);
        printf("    L3:    %10.6lf\n", e->L3);
        printf("    L4:    %10.6lf\n", e->L4);
        printf("    R1:    %9.4lf\n",  e->R1);
        printf("    R2:    %9.4lf\n",  e->R2);
        printf("    R3:    %9.4lf\n",  e->R3);
        printf("    R4:    %9.4lf\n",  e->R4);
        printf("   dat:    %4d\n",     e->dat);
        printf("EOPPWk:    %5d\n",     e->EOPPWk);
        printf("  teff:    %6d\n",     e->teff);
        Lgm_MJD_to_Date( (double)(e->teff), &y, &m, &d, &UTC);
        printf("         %4d %2d %2d %lf\n", y, m, d, UTC);

    }

    return(1);

}

/*
 *   $Id$
 */

