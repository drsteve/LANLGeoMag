// THIS CODE REALLY NEEDS CLEANING UP.
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "Lgm/Lgm_QuadPack.h"
#include "Lgm/Lgm_MagModelInfo.h"
#include "Lgm/Lgm_CTrans.h"
#include "Lgm/Lgm_LstarInfo.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <string.h>

#define TRACE_TOL   1e-7

void PredictMlat1( double *MirrorMLT, double *MirrorMlat, int k, double MLT, double *pred_mlat, double *pred_delta_mlat, double *delta );
void PredictMlat2( double *MirrorMLT, double *MirrorMlat, int k, double MLT, double *pred_mlat, double *pred_delta_mlat, double *delta, Lgm_LstarInfo *LstarInfo );
int FitQuadAndFindZero( double *x, double *y, double *dy, int n, double *res );

/*
 * Routine determines type of FL that is present in the LstarInfo->mInfo
 * structure (FL has to have already been traced!)
 *
 * Returns a flag that is set as follows;
 *
 *     -1: Too many minima -- bad or weird field line.
 *      0: Field line has only a single minimum.
 *      1: Field line has multiple minima, but these wont affect particles for the pitch angle implied by LstarInfo->mInfo->Bm.
 *      2: Field line has multiple minima, and these will affect particles for the pitch angle implied by LstarInfo->mInfo->Bm.
 *
 */
int ClassifyFL( int k, Lgm_LstarInfo *LstarInfo ) {

    Lgm_MagModelInfo    *m;
    
    int     i, iMin, Type, iMax, done, Verbosity, nBounceRegions;
    double  Min, Max, Curr, Prev;

    double  Diff, OldDiff, Minima[ LGM_LSTARINFO_MAX_FL ], Maxima[ LGM_LSTARINFO_MAX_FL ];
    int     nMinima, iMinima[ LGM_LSTARINFO_MAX_FL ];
    int     nMaxima, iMaxima[ LGM_LSTARINFO_MAX_FL ];

    Verbosity = LstarInfo->VerbosityLevel;
    m = LstarInfo->mInfo;

    Type = 0;

    Max  = LGM_FILL_VALUE;
    iMax = -1;
    
    /*
     *   Find Global minimum
     */
    Min = 9e99;
    for (i=0; i<m->nPnts; i++){
        if ( m->Bmag[i] < Min ) {
            Min = m->Bmag[i];
            iMin = i;
        }
    }


    /*
     *   Find number of minimums and maximums
     */
    nMinima = 0;
    nMaxima = 0;
    Diff = LGM_FILL_VALUE;
    for ( i=1; i<m->nPnts; i++ ){
        
        OldDiff = Diff;
        Diff    = m->Bmag[i] - m->Bmag[i-1];

        if ( (Diff > 0.0) && (OldDiff < 0.0) && (nMinima <  LGM_LSTARINFO_MAX_MINIMA )  ) {
            iMinima[nMinima] = i-1;
            Minima[nMinima] = m->Bmag[i-1];
            ++nMinima;
        }

        if ( (Diff < 0.0) && (OldDiff > 0.0) && (nMaxima <  LGM_LSTARINFO_MAX_MINIMA ) ) {
            iMaxima[nMaxima] = i-1;
            Maxima[nMaxima] = m->Bmag[i-1];
//printf("m->Bmag[i-1], m->Bmag[i] = %g %g\n", m->Bmag[i-1], m->Bmag[i]);
            ++nMaxima;
        }

    }
//Verbosity = 3;
    if ( Verbosity > 1 ) {

        if ( nMinima > 1 ) {
            if ( nMinima >=  LGM_LSTARINFO_MAX_MINIMA-1 ) {
                printf("\t\tThis Field line has multiple minima (at least %d minima detected). Probably on a bizzare field line (TS04 model?)\n", nMinima );
            } else {
                printf("\t\tThis Field line has multiple minima (%d minima detected). Probably on Shabansky orbit!\n", nMinima );
            }
            if ( Verbosity > 2 ) {
                for ( i=0; i<nMinima; i++ ){
                    printf( "\t\t\tLocal Minima %02d: Bmag[%d] = %g     (Bmirror = %g)\n", i, iMinima[i], Minima[i], m->Bm );
                }
            }

        }

        if ( nMaxima > 1 ) {
            if ( nMinima >= LGM_LSTARINFO_MAX_MINIMA-1 ) {
                printf("\t\tThis Field line has multiple maxima (at least %d maxima detected). Probably on a bizzare field line (TS04 model?)\n", nMaxima );
            } else {
                printf("\t\tThis Field line has multiple maxima (%d minima detected excluding endpoints). Probably on Shabansky orbit!\n", nMaxima );
            }
            if ( Verbosity > 2 ) {
                for ( i=0; i<nMaxima; i++ ){
                    printf( "\t\t\tLocal Maxima %02d: Bmag[%d] = %g     (Bmirror = %g)\n", i, iMaxima[i], Maxima[i], m->Bm );
                }
            }
        }

    }

    

    /*
     * Save some info.
     * Probably want to save more...?
     */
    LstarInfo->nMinima[k] = nMinima;
    LstarInfo->nMaxima[k] = nMaxima;

    if ( nMinima == 0 ) { // can this ever happen???

        Type = -9;

    } else if ( nMinima == 1 ) {

        Type = 0;

    } else if ( (nMinima == 2) && (nMaxima == 1) && (iMinima[0] < iMaxima[0]) && (iMinima[1] > iMaxima[0]) ) { // typical expected case for Shab.

        nBounceRegions = 0;
        if ( Maxima[0] < m->Bm ) {

            // can only bounce in 1 region.
            nBounceRegions = 1;

        } else {

            // can potentially bounce in 2 regions.
            if ( Minima[0] <= m->Bm ) ++nBounceRegions;
            if ( Minima[1] <= m->Bm ) ++nBounceRegions;

        }

        Type = ( nBounceRegions > 0 ) ?  nBounceRegions : -8;

    } else if ( nMinima > LGM_LSTARINFO_MAX_MINIMA-5 ) { // way too many minima -- bad/bizzare FL (TS04 seems to give these?)

        Type = -1;

    } else { // there are perhaps multiple minima/maxima. We dont expect this, but....

        /*
         *  OK, there are more than 2 minima.  Loop over them and determine how
         *  many separate regions on the FL we can mirror on.  If there is only
         *  one, then we should be a Type 1 Shabansky orbit (i.e. no
         *  biffurcation, but still have multiple minima).
         *  Count up how many cases we have where;
         *     1) a maxima exceeds the Bm values.
         *     2) minima that drop below Bm on either side.
         *  The number of mirroring regions should be this number plus one.
         */
        nBounceRegions = 1; // there should be at least one region....
        for ( i=0; i<nMaxima; i++ ) {
            if ( (Minima[i] <= m->Bm) && (Maxima[i] > m->Bm) && (Minima[i+1] <= m->Bm) ) {
                ++nBounceRegions;
            }
        }

Type = nBounceRegions;

    }



//if ( Type == -8 ) {
//    FILE *fp=fopen("TwoMinima.txt", "a");
//    for ( i=0; i<m->nPnts; i++ ){
//        fprintf(fp, "%g %g\n", m->s[i], m->Bmag[i] );
//    }
//    fclose(fp);
//}

    
//if (Type > 0 ) printf("Type = %d\n", Type);

    return( Type );

}














/*
 * There are many tolerances involved in an Lstar calculation.  Here we try and
 * set them qualitatively given a single "quality" value.
 */
void Lgm_SetLstarTolerances( int Quality, int nFLsInDriftShell, Lgm_LstarInfo *s ) {

    if ( ( Quality < 0 ) || ( Quality > 8 ) ) {
        printf("%sLgm_SetLstarTolerances: Quality value (of %d) not in range [0, 8]. Setting to 3.%s\n", s->PreStr, Quality, s->PostStr );
        s->LstarQuality = 3;
    } else {
        s->LstarQuality = Quality;
    }

    if ( ( nFLsInDriftShell < 0 ) || ( nFLsInDriftShell > 240 ) ) {
        printf("%sLgm_SetLstarTolerances: nFLsInDriftShell value (of %d) not in range [0, 240]. Setting to 24.%s\n", s->PreStr, nFLsInDriftShell, s->PostStr );
        s->nFLsInDriftShell = 24;
    } else {
        s->nFLsInDriftShell = nFLsInDriftShell;
    }

    // These tend to be critical to keep things working smoothly.
//    s->mInfo->Lgm_FindBmRadius_Tol       = 1e-10;
//    s->mInfo->Lgm_TraceToMirrorPoint_Tol = 1e-10;
    s->mInfo->Lgm_FindBmRadius_Tol       = 1e-10;
    s->mInfo->Lgm_TraceToMirrorPoint_Tol = 1e-10;


    switch ( s->LstarQuality ) {

        case 8: // Highest Quality -- (although 7 may be better?)

            s->mInfo->Lgm_MagFlux_Integrator_epsabs = 1e-9;
            s->mInfo->Lgm_MagFlux_Integrator_epsrel = 1e-9;

            s->mInfo->Lgm_LambdaIntegral_Integrator_epsabs = 1e-9;
            s->mInfo->Lgm_LambdaIntegral_Integrator_epsrel = 1e-9;

            s->mInfo->Lgm_I_Integrator        = DQAGS;
            s->mInfo->Lgm_I_Integrator_epsrel = 1e-9;
            s->mInfo->Lgm_I_Integrator_epsabs = 1e-9;

            s->mInfo->Lgm_FindShellLine_I_Tol = 1e-9;

            s->mInfo->nDivs = 500;

            s->mInfo->Lgm_MagStep_BS_atol       = 1e-7;
            s->mInfo->Lgm_MagStep_BS_rtol       = 0.0;

            break;

        case 7:

            s->mInfo->Lgm_MagFlux_Integrator_epsabs = 1e-9;
            s->mInfo->Lgm_MagFlux_Integrator_epsrel = 1e-9;

            s->mInfo->Lgm_LambdaIntegral_Integrator_epsabs = 1e-9;
            s->mInfo->Lgm_LambdaIntegral_Integrator_epsrel = 1e-9;

            s->mInfo->Lgm_I_Integrator        = DQAGS;
            s->mInfo->Lgm_I_Integrator_epsrel = 1e-8;
            s->mInfo->Lgm_I_Integrator_epsabs = 1e-8;

            s->mInfo->Lgm_FindShellLine_I_Tol = 1e-9;

            s->mInfo->nDivs = 400;

            s->mInfo->Lgm_MagStep_BS_atol       = 1e-7;
            s->mInfo->Lgm_MagStep_BS_rtol       = 0.0;

            break;

        case 6:

            s->mInfo->Lgm_MagFlux_Integrator_epsabs = 1e-8;
            s->mInfo->Lgm_MagFlux_Integrator_epsrel = 1e-8;

            s->mInfo->Lgm_LambdaIntegral_Integrator_epsabs = 1e-8;
            s->mInfo->Lgm_LambdaIntegral_Integrator_epsrel = 1e-8;

            s->mInfo->Lgm_I_Integrator        = DQAGS;
            s->mInfo->Lgm_I_Integrator_epsrel = 1e-7;
            s->mInfo->Lgm_I_Integrator_epsabs = 1e-7;

            s->mInfo->Lgm_FindShellLine_I_Tol = 1e-6;

            s->mInfo->nDivs = 300;

            s->mInfo->Lgm_MagStep_BS_atol       = 1e-7;
            s->mInfo->Lgm_MagStep_BS_rtol       = 0.0;

            break;

        case 5:

            s->mInfo->Lgm_MagFlux_Integrator_epsabs = 1e-7;
            s->mInfo->Lgm_MagFlux_Integrator_epsrel = 1e-7;

            s->mInfo->Lgm_LambdaIntegral_Integrator_epsabs = 1e-7;
            s->mInfo->Lgm_LambdaIntegral_Integrator_epsrel = 1e-7;

            s->mInfo->Lgm_I_Integrator        = DQAGS;
            s->mInfo->Lgm_I_Integrator_epsrel = 1e-6;
            s->mInfo->Lgm_I_Integrator_epsabs = 1e-6;

            s->mInfo->Lgm_FindShellLine_I_Tol = 1e-5;

            s->mInfo->nDivs = 200;

            s->mInfo->Lgm_MagStep_BS_atol       = 1e-7;
            s->mInfo->Lgm_MagStep_BS_rtol       = 0.0;

            break;

        case 4:

            s->mInfo->Lgm_MagFlux_Integrator_epsabs = 1e-6;
            s->mInfo->Lgm_MagFlux_Integrator_epsrel = 1e-6;

            s->mInfo->Lgm_LambdaIntegral_Integrator_epsabs = 1e-6;
            s->mInfo->Lgm_LambdaIntegral_Integrator_epsrel = 1e-6;

            s->mInfo->Lgm_I_Integrator        = DQAGS;
            s->mInfo->Lgm_I_Integrator_epsrel = 1e-5;
            s->mInfo->Lgm_I_Integrator_epsabs = 1e-5;

            s->mInfo->Lgm_FindShellLine_I_Tol = 1e-4;

            s->mInfo->nDivs = 200;

            s->mInfo->Lgm_MagStep_BS_atol       = 1e-6;
            s->mInfo->Lgm_MagStep_BS_rtol       = 0.0;

            break;

        case 3:

            s->mInfo->Lgm_MagFlux_Integrator_epsabs = 1e-5;
            s->mInfo->Lgm_MagFlux_Integrator_epsrel = 1e-5;

            s->mInfo->Lgm_LambdaIntegral_Integrator_epsabs = 1e-5;
            s->mInfo->Lgm_LambdaIntegral_Integrator_epsrel = 1e-5;

            s->mInfo->Lgm_I_Integrator        = DQAGS;
            s->mInfo->Lgm_I_Integrator_epsabs = 1e-5;
            s->mInfo->Lgm_I_Integrator_epsrel = 1e-5;

            s->mInfo->Lgm_FindShellLine_I_Tol = 1e-3;

            s->mInfo->nDivs = 200;

            s->mInfo->Lgm_MagStep_BS_atol       = 1e-6;
            s->mInfo->Lgm_MagStep_BS_rtol       = 0.0;

            break;

        case 2:

            s->mInfo->Lgm_MagFlux_Integrator_epsabs = 1e-4;
            s->mInfo->Lgm_MagFlux_Integrator_epsrel = 1e-4;

            s->mInfo->Lgm_LambdaIntegral_Integrator_epsabs = 1e-4;
            s->mInfo->Lgm_LambdaIntegral_Integrator_epsrel = 1e-4;

            s->mInfo->Lgm_I_Integrator        = DQAGS;
            s->mInfo->Lgm_I_Integrator_epsrel = 1e-3;
            s->mInfo->Lgm_I_Integrator_epsabs = 1e-3;

            s->mInfo->Lgm_FindShellLine_I_Tol = 5e-3;

            s->mInfo->nDivs = 200;

            s->mInfo->Lgm_MagStep_BS_atol       = 1e-4;
            s->mInfo->Lgm_MagStep_BS_rtol       = 0.0;

            break;

        case 1:

            s->mInfo->Lgm_MagFlux_Integrator_epsabs = 1e-4;
            s->mInfo->Lgm_MagFlux_Integrator_epsrel = 1e-4;

            s->mInfo->Lgm_LambdaIntegral_Integrator_epsabs = 1e-4;
            s->mInfo->Lgm_LambdaIntegral_Integrator_epsrel = 1e-4;

            s->mInfo->Lgm_I_Integrator        = DQAGS;
            s->mInfo->Lgm_I_Integrator_epsrel = 1e-2;
            s->mInfo->Lgm_I_Integrator_epsabs = 1e-2;

            s->mInfo->Lgm_FindShellLine_I_Tol = 1e-1;

            s->mInfo->nDivs = 100;

            s->mInfo->Lgm_MagStep_BS_atol       = 1e-3;
            s->mInfo->Lgm_MagStep_BS_rtol       = 0.0;

            break;

        case 0:

            s->mInfo->Lgm_MagFlux_Integrator_epsabs = 5e-4;
            s->mInfo->Lgm_MagFlux_Integrator_epsrel = 5e-4;

            s->mInfo->Lgm_LambdaIntegral_Integrator_epsabs = 5e-4;
            s->mInfo->Lgm_LambdaIntegral_Integrator_epsrel = 5e-4;

            s->mInfo->Lgm_I_Integrator        = DQK21; // Note - changed to simpler integrator.
            s->mInfo->Lgm_I_Integrator_epsrel = 1e-2;
            s->mInfo->Lgm_I_Integrator_epsabs = 1e-2;

            s->mInfo->Lgm_FindShellLine_I_Tol = 1e-1;

            s->mInfo->nDivs = 50;

            s->mInfo->Lgm_MagStep_BS_atol       = 1e-3;
            s->mInfo->Lgm_MagStep_BS_rtol       = 0.0;

            break;

    }


    return;


}



Lgm_LstarInfo *InitLstarInfo( int VerbosityLevel ) {

    Lgm_LstarInfo	*LstarInfo;


    /*
     *  Allocate memory for LstarInfo structure
     */
    LstarInfo = (Lgm_LstarInfo *) calloc( 1, sizeof( *LstarInfo));




    /*
     *  Allocate memory for mInfo structure inside the LstarInfo structure.
     */
    LstarInfo->mInfo = Lgm_InitMagInfo( );

    Lgm_InitLstarInfoDefaults(LstarInfo);
    return(LstarInfo);
}


void Lgm_InitLstarInfoDefaults( Lgm_LstarInfo	*LstarInfo ) {
    /*
     *  Default Settings
     */
    LstarInfo->VerbosityLevel = 2;
    LstarInfo->LSimpleMax     = 10.0;
    LstarInfo->ISearchMethod  = 1;
    LstarInfo->ShabanskyHandling = LGM_SHABANSKY_IGNORE;
    LstarInfo->LstarMoment    = LGM_LSTAR_MOMENT_CDIP_2010;

    LstarInfo->PreStr[0]  = '\0';
    LstarInfo->PostStr[0] = '\0';
    LstarInfo->ComputeVgc = FALSE;

    LstarInfo->ComputeSbIntegral = TRUE;
    LstarInfo->mInfo->ComputeSb0 = TRUE;

}



void FreeLstarInfo( Lgm_LstarInfo *s ) {

    Lgm_FreeMagInfo( s->mInfo );
    free( s );

}



/*
 *  The Lgm_LstarInfo structure has pointers in it, so simple
 *  asignments (e.g. *t = *s) are dangerous. Here we make sure that
 *  the target gets an independent copy of the structure.
 */
Lgm_LstarInfo *Lgm_CopyLstarInfo( Lgm_LstarInfo *s ) {

    Lgm_LstarInfo       *t;
    Lgm_MagModelInfo *m;
    Lgm_CTrans       *c;

    if ( s == NULL) {
        printf("%sLgm_CopyLstarInfo: Error, source structure is NULL%s\n", s->PreStr, s->PostStr);
        return( NULL );
    }


    /*
     *  Allocate memory for target Lgm_LstarInfo struct
     */
    t = (Lgm_LstarInfo *)calloc(1, sizeof(Lgm_LstarInfo));

    /*
     *  Copy memory
     */
    memcpy( t, s, sizeof(*s) );


    /*
     *  Copy MagInfo struct
     */
    t->mInfo = Lgm_CopyMagInfo( s->mInfo );


    t->mInfo->Lgm_MagStep_RK5_FirstTimeThrough = TRUE;
    t->mInfo->Lgm_MagStep_BS_FirstTimeThrough = TRUE;
    t->mInfo->Lgm_MagStep_BS_eps_old = -1.0;


    return( t );
}








void NewTimeLstarInfo( long int Date, double UT, double PitchAngle, int (*Mag)(Lgm_Vector *, Lgm_Vector *, Lgm_MagModelInfo *), Lgm_LstarInfo *LstarInfo ){


    /*
     *  Set coord transforms for given data and time
     */
    Lgm_Set_Coord_Transforms( Date, UT, LstarInfo->mInfo->c );


    /*
     *  Set mag field model
     */
    LstarInfo->mInfo->Bfield = Mag;


    /*
     *  Set Pitch Angle
     */
    LstarInfo->PitchAngle = PitchAngle;


}




int Lstar( Lgm_Vector *vin, Lgm_LstarInfo *LstarInfo ){


    Lgm_Vector	u, v, w, v1, v2, v3, Bvec, uu;
    int		i, j, k, nk, nLines, koffset, tkk, nfp, nnn;
    int		done2, Count, FoundShellLine, nIts, Type;
    double	rat, B, dSa, dSb, smax, SS, L, Hmax, epsabs, epsrel;
    double	I=-999.9, Ifound, M, MLT0, MLT, DeltaMLT, mlat, r;
    double	Phi, Phi1, Phi2, sl, cl, MirrorMLT[3*LGM_LSTARINFO_MAX_FL], MirrorMlat[3*LGM_LSTARINFO_MAX_FL], pred_mlat, mlat_try, pred_delta_mlat=0.0, mlat0, mlat1, delta;
    double	MirrorMLT_Old[3*LGM_LSTARINFO_MAX_FL], MirrorMlat_Old[3*LGM_LSTARINFO_MAX_FL], res;
    char    *PreStr, *PostStr;
    Lgm_MagModelInfo    *mInfo2;
//    FILE	*fp;

    PreStr = LstarInfo->PreStr;
    PostStr = LstarInfo->PostStr;
    double PredMinusActualMlat = 0.0;

    switch(LstarInfo->LstarMoment) {
        case LGM_LSTAR_MOMENT_CDIP :
            LstarInfo->Mused = LstarInfo->mInfo->c->M_cd;
            break;
        case LGM_LSTAR_MOMENT_CDIP_2010 :
            LstarInfo->Mused = LstarInfo->mInfo->c->M_cd_2010;
            break;
        case LGM_LSTAR_MOMENT_MCILWAIN :
            LstarInfo->Mused = LstarInfo->mInfo->c->M_cd_McIllwain;
            break;
        default :
            LstarInfo->Mused = LstarInfo->mInfo->c->M_cd_2010;
    }


    /*
     * Initialize some values to FILL in case we bail early
     */
    LstarInfo->LS = LGM_FILL_VALUE;


    if (LstarInfo->VerbosityLevel > 1) {
        printf("\n\n\t\t%s        Computing L* for, Date: %ld, UT: %g%s\n", PreStr, LstarInfo->mInfo->c->UTC.Date, LstarInfo->mInfo->c->UTC.Time, PostStr );
        printf(    "\t\t%s=========================================================================%s\n", PreStr, PostStr );
        printf(    "\t\t%sDate (yyyymmdd):                         %ld%s\n", PreStr, LstarInfo->mInfo->c->UTC.Date, PostStr );
        printf(    "\t\t%sUTC (hours):                             %g%s\n", PreStr, LstarInfo->mInfo->c->UTC.Time, PostStr );
        printf(    "\t\t%sPitch Angle (deg.):                      %g%s\n", PreStr, LstarInfo->PitchAngle, PostStr );
        printf(    "\t\t%sInitial Position, vin (Re):              < %g, %g, %g >%s\n", PreStr, vin->x, vin->y, vin->z, PostStr);
        printf(    "\t\t%sMirror Mag. Field Strength, Bm (nT):     %g%s\n", PreStr, LstarInfo->mInfo->Bm, PostStr );
    }
    LstarInfo->nPnts = 0;


    if ((LstarInfo->PitchAngle < 0.0)||(LstarInfo->PitchAngle>90.0)) return(-1);

    for (k=0; k<LGM_LSTARINFO_MAX_FL; k++){
        LstarInfo->I[k] = LGM_FILL_VALUE;
    }




    /*
     *  Do Initial field Line to get Bm and I
     */
    u = *vin;
    if (LstarInfo->VerbosityLevel > 1) {
        printf("\n\t\t%sInitial Position, U_gsm (Re):            < %g, %g, %g >%s\n", PreStr, u.x, u.y, u.z, PostStr);
	    LstarInfo->mInfo->Bfield( &u, &Bvec, LstarInfo->mInfo );
	    B = Lgm_Magnitude( &Bvec );
        printf("\t\t%sMag. Field Strength, B at U_gsm (nT):    %g%s\n", PreStr, B, PostStr);
    }
    if ( Lgm_Trace( &u, &v1, &v2, &v3, LstarInfo->mInfo->Lgm_LossConeHeight, TRACE_TOL, TRACE_TOL, LstarInfo->mInfo ) == LGM_CLOSED ) {

        LstarInfo->Sb0     = LstarInfo->mInfo->Sb0; // Equatorial value of Sb Integral.
        LstarInfo->d2B_ds2 = LstarInfo->mInfo->d2B_ds2; // second derivative of B wrt s at equator.
        LstarInfo->RofC    = LstarInfo->mInfo->d2B_ds2; // radius of curvature at Bmin point.

	    if (LstarInfo->VerbosityLevel > 1) {
            printf("\n\t\t%sMin-B  Point Location, Pmin (Re):      < %g, %g, %g >%s\n", PreStr, LstarInfo->mInfo->Pmin.x, LstarInfo->mInfo->Pmin.y, LstarInfo->mInfo->Pmin.z, PostStr);
	        LstarInfo->mInfo->Bfield( &u, &Bvec, LstarInfo->mInfo );
	        B = Lgm_Magnitude( &Bvec );
            printf("\t\t%sMag. Field Strength, B at Pmin (nT):    %g%s\n", PreStr, B, PostStr);
        }

        /*
         *  If we are here, we know the field line is closed.  
         */

/*  Now, lets classify what type of FL it is. This call determines
*  whether or not there are multiple minima on the FL. 
*
*  If there are multiple minima, then we need to figure out if, for
*  the current pitch angle, we can mirror in multiple regions (i.e.
*  Shabansky drift orbit type bifurcations.)
*
*/
//// DAMN!!!!!!!!!!!!!!!!!!!1
//This wont work!
//you need to trace the FL with TraceLine for ClassifyFL() to work!!!!!!!!!!!!!
//        Type = ClassifyFL( 0, LstarInfo );


        /*
         * Trace from Bmin point up to northern mirror point and down to
         * southern mirror point. dSa and dSb are the distances along the
         * FL from the starting points point. So dSa is from the Bmin value.
         * And dSb is from the Psouth value.
         */
        if ( Lgm_TraceToMirrorPoint( &(LstarInfo->mInfo->Pmin), &(LstarInfo->mInfo->Pm_South), &dSa, LstarInfo->mInfo->Bm, -1.0, LstarInfo->mInfo->Lgm_TraceToMirrorPoint_Tol, LstarInfo->mInfo ) >= 0 ) {

	        if (LstarInfo->VerbosityLevel > 1) {
                printf("\n\t\t%sMirror Point Location, Pm_South (Re):      < %g, %g, %g >%s\n", PreStr, LstarInfo->mInfo->Pm_South.x, LstarInfo->mInfo->Pm_South.y, LstarInfo->mInfo->Pm_South.z, PostStr);
	            LstarInfo->mInfo->Bfield( &LstarInfo->mInfo->Pm_South, &Bvec, LstarInfo->mInfo );
	            B = Lgm_Magnitude( &Bvec );
                printf("\t\t%sMag. Field Strength, Bm at Pm_South (nT):  %g     (LstarInfo->mInfo->Bm = %g    dSa = %g)%s\n", PreStr, B, LstarInfo->mInfo->Bm, dSa, PostStr);
            }

            if ( Lgm_TraceToMirrorPoint( &(LstarInfo->mInfo->Pmin), &(LstarInfo->mInfo->Pm_North), &dSb, LstarInfo->mInfo->Bm,  1.0, LstarInfo->mInfo->Lgm_TraceToMirrorPoint_Tol, LstarInfo->mInfo ) >= 0 ) {


                /*
                 *  We found both mirror points. Dump some diagnostics.
                 */
                if (LstarInfo->VerbosityLevel > 1) {
                    printf("\n\t\t%sMirror Point Location, Pm_North (Re):      < %g, %g, %g >%s\n",
                            PreStr, LstarInfo->mInfo->Pm_North.x, LstarInfo->mInfo->Pm_North.y,
                            LstarInfo->mInfo->Pm_North.z, PostStr);
                    LstarInfo->mInfo->Bfield( &LstarInfo->mInfo->Pm_North, &Bvec, LstarInfo->mInfo );
                    B = Lgm_Magnitude( &Bvec );
                    printf("\t\t%sMag. Field Strength, Bm at Pm_North (nT):  %g     (LstarInfo->mInfo->Bm = %g    dSb = %g)%s\n",
                            PreStr, B, LstarInfo->mInfo->Bm, dSb, PostStr);
                }




                /*
                 *  Set the limits of integration. Define s=0 at the sourthern mirror point. Then, sm_North will just be dSb
                 */
                //LstarInfo->mInfo->Sm_South = LstarInfo->mInfo->smin - dSa;
                //LstarInfo->mInfo->Sm_North = LstarInfo->mInfo->Sm_South + dSb;
                //SS = dSb;
                SS = dSa+dSb;

                //LstarInfo->mInfo->Hmax = SS/200.0;
                LstarInfo->mInfo->Hmax = SS/(double)LstarInfo->mInfo->nDivs;
                r  = Lgm_Magnitude( &LstarInfo->mInfo->Pm_North );
                LstarInfo->mInfo->Sm_South = 0.0;
                LstarInfo->mInfo->Sm_North = SS;

                if ( SS <= 1e-5 ) {
                    // if FL length is small, use an approx expression for I
                    rat = LstarInfo->mInfo->Bmin/LstarInfo->mInfo->Bm;
                    if ((1.0-rat) < 0.0) {
                        I = 0.0;
                    } else {
                        // Eqn 2.66b in Roederer
                        I = SS*sqrt(1.0 - rat);
                    }

                } else if ( LstarInfo->mInfo->UseInterpRoutines ) {

                    //if ( Lgm_TraceLine2( &(LstarInfo->mInfo->Pm_South), &LstarInfo->mInfo->Pm_North, (r-1.0)*Re, 0.5*SS-LstarInfo->mInfo->Hmax, 1.0, 1e-7, FALSE, LstarInfo->mInfo ) < 0 ) return(LGM_FILL_VALUE);
                    if ( Lgm_TraceLine3( &(LstarInfo->mInfo->Pm_South), SS, LstarInfo->mInfo->nDivs, 1.0, 1e-7, FALSE, LstarInfo->mInfo ) < 0 ) return( -1 );
                    // Lgm_TraceLine4() doesnt seem to work -- dont use it....
                    //if ( Lgm_TraceLine4( &(LstarInfo->mInfo->Pm_South), &(LstarInfo->mInfo->Pm_North), dSa, dSb, LstarInfo->mInfo->nDivs, FALSE, LstarInfo->mInfo ) < 0 ) return( -1 );

                    //ReplaceFirstPoint( 0.0, LstarInfo->mInfo->Bm, &LstarInfo->mInfo->Pm_South, LstarInfo->mInfo );
                    //AddNewPoint( SS,  LstarInfo->mInfo->Bm, &LstarInfo->mInfo->Pm_North, LstarInfo->mInfo );
                    //ReplaceLastPoint( SS,  LstarInfo->mInfo->Bm, &LstarInfo->mInfo->Pm_North, LstarInfo->mInfo );
                    if ( InitSpline( LstarInfo->mInfo ) ) {

                        /*
                         *  Do interped I integral.
                         */
//                        epsabs = LstarInfo->mInfo->Lgm_I_Integrator_epsabs;
//                        epsrel = LstarInfo->mInfo->Lgm_I_Integrator_epsrel;
//                        LstarInfo->mInfo->Lgm_I_Integrator_epsabs /= 10.0;
//                        LstarInfo->mInfo->Lgm_I_Integrator_epsrel /= 10.0;
                        I = Iinv_interped( LstarInfo->mInfo  );
                        if (LstarInfo->VerbosityLevel > 1) {
                            printf("\t\t  %sIntegral Invariant, I (interped):      %g%s\n",  PreStr, I, PostStr );
                        }
//                        LstarInfo->mInfo->Lgm_I_Integrator_epsabs = epsabs;
//                        LstarInfo->mInfo->Lgm_I_Integrator_epsrel = epsrel;

                        /*
                         *  Do interped Sb integral if desired. Note this is
                         *  just for the initial FL (which is why we call it
                         *  SbIntegral0). Can also add for each MLT if we want
                         *  to (then we should probably add an array called
                         *  SbIntegral[]).
                         */
                        LstarInfo->SbIntegral0 = LGM_FILL_VALUE;
                        if ( LstarInfo->ComputeSbIntegral ) {
                            LstarInfo->SbIntegral0 = SbIntegral_interped( LstarInfo->mInfo  );
                            if (LstarInfo->VerbosityLevel > 1) {
                                printf("\t\t  %sSb Integral Equatorially Mirroring:      %g%s\n",  PreStr, LstarInfo->Sb0, PostStr );
                                printf("\t\t  %sSb Integral, (interped):      %g%s\n",  PreStr, LstarInfo->SbIntegral0, PostStr );
                            }
                        }



                        FreeSpline( LstarInfo->mInfo );

                    } else {

                        I = LGM_FILL_VALUE;

                    }



                } else {

                    /*
                     *  Do full blown I integral. (Integrand is evaluated by tracing to required s-values.)
                     *  (This strategy isnt used very much anymore - 20120524, MGH)
                     */
                    I = Iinv( LstarInfo->mInfo  );
                    if (LstarInfo->VerbosityLevel > 1) printf("\t\t  %sIntegral Invariant, I (full integral): %g%s\n",  PreStr, I, PostStr );

                }
                LstarInfo->I0 = I; // save initial I in LstarInfo structure.
                Ifound = I;



                if (LstarInfo->VerbosityLevel > 1) {
                    printf("\t\t  %sLgm_n_I_integrand_Calls:               %d%s\n\n", PreStr, LstarInfo->mInfo->Lgm_n_I_integrand_Calls, PostStr );
                    printf("\t\t%sCurrent Dipole Moment, M_cd:             %g%s\n", PreStr, LstarInfo->mInfo->c->M_cd, PostStr);
                    printf("\t\t%sReference Dipole Moment, M_cd_McIllwain: %g%s\n", PreStr, LstarInfo->mInfo->c->M_cd_McIllwain, PostStr);
                    printf("\t\t%sReference Dipole Moment, M_cd_2010:      %g%s\n", PreStr, LstarInfo->mInfo->c->M_cd_2010, PostStr);
                    printf("\t\t%sDipole Moment Used, Mused:               %g%s\n", PreStr, LstarInfo->Mused, PostStr);
                    printf("\t\t%sMcIlwain L (Hilton):                     %.15g%s\n", PreStr, L = LFromIBmM_Hilton( I, LstarInfo->mInfo->Bm, LstarInfo->Mused ), PostStr );
                    printf("\t\t%sMcIlwain L (McIlwain):                   %.15g%s\n", PreStr, L = LFromIBmM_McIlwain( I, LstarInfo->mInfo->Bm, LstarInfo->Mused ), PostStr );
                }

            } else {

	            if (LstarInfo->VerbosityLevel > 0) printf("\t\t%sMirror point below %g km in Northern Hemisphere%s\n", PreStr, LstarInfo->mInfo->Lgm_LossConeHeight, PostStr );
	            return(-1);

	        }

        } else {

	        if (LstarInfo->VerbosityLevel > 0) printf("\t\t%sMirror point below %g km in Southern Hemisphere%s\n", PreStr, LstarInfo->mInfo->Lgm_LossConeHeight, PostStr );
	        return(-2);

        }

    } else {

	        if (LstarInfo->VerbosityLevel > 0) printf("\t\t%sOpen Field Line%s\n", PreStr, PostStr );
	        return(-5);

    }






    /*
     *  If we get here, then the field line going through our starting point is
     *  closed and we could find the required mirror points.
     */




    // initialize DriftOrbitType as open
    //for ( i=0; i<LGM_LSTARINFO_MAX_FL; ++i ) if ( LstarInfo->nMinima[i] > 1 ) LstarInfo->DriftOrbitType = LGM_DRIFT_ORBIT_OPEN;
    //for ( i=0; i<LGM_LSTARINFO_MAX_FL; ++i ) LstarInfo->DriftOrbitType = LGM_DRIFT_ORBIT_OPEN;



    /*
     *  Construct drift shell. Get a good initial estimate for mlat
     */
    if        ( LstarInfo->ISearchMethod  == 1 ) {
        // this method uses the mlat as the mlat of the mirror point.
        Lgm_Convert_Coords( &LstarInfo->mInfo->Pm_North, &u, GSM_TO_SM, LstarInfo->mInfo->c );
   
    } else if ( LstarInfo->ISearchMethod  == 2 ) {

        // this method uses the mlat as the mlat of the footpoint.
        Lgm_Convert_Coords( &v2, &u, GSM_TO_SM, LstarInfo->mInfo->c );

    } else {

        printf("Unrecognized value for LstarInfo->ISearchMethod ( = %d)\n", LstarInfo->ISearchMethod );

    }
    mlat = DegPerRad*asin(u.z/Lgm_Magnitude( &u ));


    MLT0 = atan2( u.y, u.x )*DegPerRad/15.0 + 12.0;
    if (LstarInfo->VerbosityLevel > 2) printf("\t\t%smlat = %g%s\n", PreStr, mlat, PostStr );
    LstarInfo->nPnts = 0;
    pred_delta_mlat  = 0.0;
    for (k=0; k<3*LGM_LSTARINFO_MAX_FL; ++k){
	    MirrorMLT[k]  = 0.0;
	    MirrorMlat[k] = 0.0;
    }



    nLines = LstarInfo->nFLsInDriftShell;
    DeltaMLT = 24.0/((double)nLines);

    delta = 3.0; // default
    for ( k=0, MLT=MLT0; MLT<(MLT0+24.0-1e-10); MLT += DeltaMLT){

	    /*
	     *  Try to predict what the next mlat should be so we can really
	     *  narrow down the bracket.
	     */
	    if (k == 0 ){

	        pred_mlat 	    = mlat;
	        pred_delta_mlat = 0.001;
	        delta           = 0.001;
	        if (LstarInfo->VerbosityLevel > 2) printf("\t\t%sPredicted mlat1 = %g ( %g : %g )%s ", PreStr, pred_mlat, pred_delta_mlat, delta, PostStr ); fflush(stdout);

	    } else if ( k > 0 ) {

	        PredictMlat2( MirrorMLT, MirrorMlat, k, MLT, &pred_mlat, &pred_delta_mlat, &delta, LstarInfo );
	        if (LstarInfo->VerbosityLevel > 2) printf("\t\t%sPredicted mlat2 = %g ( %g : %g )%s ", PreStr, pred_mlat, pred_delta_mlat, delta, PostStr ); fflush(stdout);

	    } else {

	        PredictMlat1( MirrorMLT, MirrorMlat, k, MLT, &pred_mlat, &pred_delta_mlat, &delta );
	        if (LstarInfo->VerbosityLevel > 2) printf("\t\t%sPredicted mlat3 = %g ( %g : %g )%s ", PreStr, pred_mlat, pred_delta_mlat, delta, PostStr ); fflush(stdout);

	    }


        /*
         * Set the range to search over. If this turns out to be too small,
         * it'll get expanded in subsequent attempts.
         */
        if ( k == 0 ){
            delta           = 0.001;
        } else if ( k < 3 ){
            delta = 3.0;
        } else {
            if (nIts > 1) delta = 1.5*fabs( PredMinusActualMlat );
        }



        /*
         * Loop until we are done. We will be done when we have either;
         *
         *    1) found I to the requested tolerance or
         *    2) I cannot be found after several attempts.
         *
         *  Each attempt is kept track of by the 'Count' variable.  The
         *  strategy is to make a reasonable estimate for the search range on
         *  mlat. If this doesnt work, we expand the size of the range.
         *
         */
	    done2 = FALSE; FoundShellLine = FALSE; Count = 0;
        LstarInfo->nImI0 = 0;
	    while ( !done2 && (k > 0) ) {

	        if ( Count == 0 ) {

                /*
                 *  First time through -- lets use the predicted range.
                 */
                //mlat0 = ((pred_mlat-delta) <  0.0) ?  0.0 : (pred_mlat-delta);
                if (delta > 20.0) delta = 20.0;

                mlat0 = pred_mlat-delta;
                mlat_try = pred_mlat;
                mlat1 = pred_mlat+delta;

    	    } else if ( Count == 1 ) {

                /*
                 * Perhaps our bracket was no good -- too narrow. Enlarge it.
                 * Try double what we had.
                 */
                delta *= 2;
                if (delta < 1.0) delta = 1.0;
                mlat0 = pred_mlat-delta;
                mlat_try = pred_mlat;
                mlat1 = pred_mlat+delta;

                if ( LstarInfo->nImI0 > 2 ) {
                    for (i=0; i<LstarInfo->nImI0; i++) LstarInfo->Earr[i] = 1.0;
                    //FitQuadAndFindZero2( LstarInfo->MLATarr, LstarInfo->ImI0arr, LstarInfo->Earr, LstarInfo->nImI0, 4, &res );
                    //FitLineAndFindZero( LstarInfo->MLATarr, LstarInfo->ImI0arr, LstarInfo->Earr, LstarInfo->nImI0, &res );
                    FitQuadAndFindZero( LstarInfo->MLATarr, LstarInfo->ImI0arr, LstarInfo->Earr, LstarInfo->nImI0, &res );
                    if (LstarInfo->VerbosityLevel > 1){
                        printf("\t\t\t1> Fitting to available values. Predicted mlat: %g\n", res );
                    }
                    if ( fabs(res) < 90.0 ){
                        mlat_try = res;
                        mlat0 = res-1.0;
                        mlat1 = res+1.0;
                    }
                }


            } else if ( Count == 2 ) {

	            /*
	             * Hmmm... This ones stuborn! Lets be more conservative now.
	            * Try +1/-5 around the returned mlat
	            */

                //if ( FoundShellLine == -3 ) {
                //    mlat0 = ((mlat-5.0) >  -10.0) ?  -10.0 : (mlat-5.0);
                //    mlat1 = ((mlat+1.0) > 90.0) ? 90.0 : (mlat+1.0);
                //} else {
                //    mlat0 = ((mlat-1.0) >  -10.0) ?  -10.0 : (mlat-1.0);
                //    mlat1 = ((mlat+5.0) > 90.0) ? 90.0 : (mlat+5.0);
               //// }
                mlat_try = pred_mlat;
                mlat0 = mlat_try-5.0;
                mlat1 = mlat_try+1.0;

                if ( LstarInfo->nImI0 > 3 ) {
                    for (i=0; i<LstarInfo->nImI0; i++) LstarInfo->Earr[i] = 1.0;
                    //FitQuadAndFindZero2( LstarInfo->MLATarr, LstarInfo->ImI0arr, LstarInfo->Earr, LstarInfo->nImI0, 4, &res );
                    FitQuadAndFindZero( LstarInfo->MLATarr, LstarInfo->ImI0arr, LstarInfo->Earr, LstarInfo->nImI0, &res );
                    if (LstarInfo->VerbosityLevel > 1){ 
                        printf("\t\t\t2> Fitting to available values. Predicted mlat: %g\n", res ); 
                    }
                    if ( fabs(res) < 90.0 ){
                        mlat_try = res;
                        mlat0 = res-5.0;
                        mlat1 = res+5.0;
                    }
                }


    	    } else {

                /*
                 * OK, there may not be a good line -- i.e. drift shell may not be closed.
                 * To make sure, lets try 0->90 deg.
                 */
                //mlat0 = 0.0;
                mlat_try = pred_mlat;
                mlat0 = -60.0;
//mlat0 = 0.0;
                mlat1 = 90.0;
//mlat1 = 45.0;
//mlat1 = 23.855537644144878;



//printf("?? LstarInfo->nImI0 = %d\n", LstarInfo->nImI0);
                if ( (LstarInfo->nImI0 > 3) && (LstarInfo->nImI0%4 == 0) ){
                    for (i=0; i<LstarInfo->nImI0; i++) LstarInfo->Earr[i] = 1.0;
                    //FitQuadAndFindZero2( LstarInfo->MLATarr, LstarInfo->ImI0arr, LstarInfo->Earr, LstarInfo->nImI0, 4, &res );
                    FitQuadAndFindZero( LstarInfo->MLATarr, LstarInfo->ImI0arr, LstarInfo->Earr, LstarInfo->nImI0, &res );
                    if (LstarInfo->VerbosityLevel > 1){
                        printf("\t\t\t3> Fitting to available values. Predicted mlat: %g\n", res );
                    }
                    if ( fabs(res) < 90.0 ){
                        mlat_try = res;
                        mlat0 = res-45.0;
                        mlat1 = res+45.0;
                    }
                }


            }


            if (LstarInfo->VerbosityLevel > 1) {
                printf("\n\t\t%s________________________________________________________________________________________________________________________________%s\n\n", PreStr, PostStr );
                printf("\t\t%s             Field Line %02d of %02d   MLT: %g  (Predicted mlat: %g %g %g  delta: %g)   Count: %d          %s\n", PreStr, k+1, nLines, MLT, mlat0, pred_mlat, mlat1, delta, Count, PostStr );
                printf("\t\t%s________________________________________________________________________________________________________________________________%s\n", PreStr, PostStr );
            }

            FoundShellLine = FindShellLine( I, &Ifound, LstarInfo->mInfo->Bm, MLT, &mlat, &r, mlat0, mlat_try, mlat1, &nIts, LstarInfo );







v = LstarInfo->mInfo->Pm_North;
LstarInfo->mInfo->Hmax = 0.1;
if ( !Lgm_TraceToSphericalEarth( &v, &w, LstarInfo->mInfo->Lgm_LossConeHeight, 1.0, 1e-9, LstarInfo->mInfo ) ){ return(-4); }
LstarInfo->Spherical_Footprint_Pn[k] = w;
LstarInfo->mInfo->Hmax = 0.05;
//Lgm_Vector puke;
//Lgm_Convert_Coords( &LstarInfo->mInfo->Pm_North, &puke, GSM_TO_SM, LstarInfo->mInfo->c );
//printf("LstarInfo->mInfo->Pm_North (SM Coords) = %g %g %g\n", puke.x, puke.y, puke.z );
//Lgm_Convert_Coords( &LstarInfo->mInfo->Pm_South, &puke, GSM_TO_SM, LstarInfo->mInfo->c );
//printf("LstarInfo->mInfo->Pm_South (SM Coords) = %g %g %g\n", puke.x, puke.y, puke.z );
//Lgm_Convert_Coords( &LstarInfo->Spherical_Footprint_Pn[k], &puke, GSM_TO_SM, LstarInfo->mInfo->c );
//printf("LstarInfo->Spherical_Footprint_Pn[k] (SM Coords) = %g %g %g\n", puke.x, puke.y, puke.z );
if (LstarInfo->VerbosityLevel > 1) {
    printf("\t\t%sTracing Full FL so that we can classify it. Starting at Spherical_Footprint_Pn[%d] = %g %g %g%s\n", PreStr, k, LstarInfo->Spherical_Footprint_Pn[k].x, LstarInfo->Spherical_Footprint_Pn[k].y, LstarInfo->Spherical_Footprint_Pn[k].z, PostStr );
}
Lgm_TraceLine( &LstarInfo->Spherical_Footprint_Pn[k], &v2, LstarInfo->mInfo->Lgm_LossConeHeight, -1.0, 1e-8, FALSE, LstarInfo->mInfo );
if (LstarInfo->VerbosityLevel > 1) {
    printf("\t\t%sTraced Full FL so that we can classify it. Step size along FL: %g. Number of points: %d.%s\n", PreStr, LstarInfo->mInfo->Hmax, LstarInfo->mInfo->nPnts, PostStr );
}
LstarInfo->Spherical_Footprint_Ps[k] = v2;
//Lgm_Convert_Coords( &LstarInfo->Spherical_Footprint_Ps[k], &puke, GSM_TO_SM, LstarInfo->mInfo->c );
//printf("LstarInfo->Spherical_Footprint_Ps[k] (SM Coords) = %g %g %g\n", puke.x, puke.y, puke.z );

            Type = ClassifyFL( k, LstarInfo );
            if (LstarInfo->VerbosityLevel > 1) {
                printf("\t\t%sClassifying FL: Type = %d. %s\n", PreStr, Type, PostStr );
            }

//for ( i=0; i<LGM_LSTARINFO_MAX_FL; ++i ) if ( LstarInfo->nMinima[i] > 2 ) LstarInfo->DriftOrbitType = LGM_DRIFT_ORBIT_OPEN;
            if (LstarInfo->ShabanskyHandling==LGM_SHABANSKY_IGNORE) {

                PredMinusActualMlat = pred_mlat - mlat;
                if (LstarInfo->VerbosityLevel > 1) {
                    printf("\t\t%s________________________________________________________________________________________________________________________________%s\n\n", PreStr, PostStr );
                    printf("\t\t%s  >>  Pred/Actual/Diff mlat:  %g/%g/%g  MLT/MLAT: %g %g  I0: %g I: %g I-I0: %g %s\n", PreStr, pred_mlat, mlat, PredMinusActualMlat, MLT, mlat, I, Ifound, Ifound-I, PostStr );
                    printf("\t\t%s________________________________________________________________________________________________________________________________ %s\n\n\n", PreStr, PostStr );
                }

            } else if ( (Type > 1) && (LstarInfo->ShabanskyHandling==LGM_SHABANSKY_HALVE_I) ){
                if (LstarInfo->VerbosityLevel > 1) {
                    printf("\t\t\t%sShabansky orbit. Re-doing FL. Target I adjusted to: %g . (Original is: %g) %s\n", PreStr, I/2.0, I, PostStr );
                }

                FoundShellLine = FindShellLine( I/2.0, &Ifound, LstarInfo->mInfo->Bm, MLT, &mlat, &r, mlat0, mlat_try, mlat1, &nIts, LstarInfo );
                //TODO: Do we need to test to make sure that the adjusted I is being found on a field line with multiple minima??
                PredMinusActualMlat = pred_mlat - mlat;
                if (LstarInfo->VerbosityLevel > 1) {
                    printf("\t\t%s________________________________________________________________________________________________________________________________%s\n\n", PreStr, PostStr );
                    printf("\t\t%s  >>  Pred/Actual/Diff mlat:  %g/%g/%g  MLT/MLAT: %g %g  I0: %g I: %g I-I0/2: %g (SHABANSKY)%s\n", PreStr, pred_mlat, mlat, PredMinusActualMlat, MLT, mlat, I, Ifound, Ifound-I/2.0, PostStr );
                    printf("\t\t%s________________________________________________________________________________________________________________________________ %s\n\n\n", PreStr, PostStr );
                }

            }
            //else if (LstarInfo->ShabanskyHandling==LGM_SHABANSKY_REJECT) {
            //    //We should exit from the search if we want to treat Shabansky orbits as undefined drift shells
            //    LstarInfo->DriftOrbitType = LGM_DRIFT_ORBIT_OPEN;
            //
            //    if (LstarInfo->VerbosityLevel >1) { printf(" \t%sNo valid I - Drift Shell not closed: L* = undefined  (FoundShellLine = %d)%s\n", PreStr, FoundShellLine, PostStr); fflush(stdout); }
            //    FoundShellLine = 0;
            //    return(-3);
            //    
            //}






            if ( FoundShellLine > 0 ) {
                done2 = TRUE;
            } else if ( Count > 2 ) {
                done2 =  TRUE;
                if (LstarInfo->VerbosityLevel >1) { printf(" \t%sNo valid I - Drift Shell not closed: L* = undefined  (FoundShellLine = %d)%s\n", PreStr, FoundShellLine, PostStr); fflush(stdout); }
                FoundShellLine = 0;
                return(-3);
            } else {
                ++Count;
            }
	    }


        if (LstarInfo->VerbosityLevel > 2) { printf("\t\t%sActual mlat = %g  MLT = %g   r = %g Ifound = %g\t Count = %d%s", PreStr, mlat, MLT, r, Ifound, Count, PostStr ); fflush(stdout); }



        /*
         *  Note that FindShellLine() takes in MLT and returns mlat, rad. The three
         *  values  (MLT, mlat, rad) are notionally supposed to represent the
         *  position of the northern mirror point.  However, if we are close to the
         *  equator, we may have initially confused the north and south mirror
         *  points. We end up sorting the confusion out, but we really need to make
         *  sure that as we proceed, we refer to the correct values.
         *
         *
         *
         */
//        u = LstarInfo->mInfo->Pm_North;
//        Lgm_Convert_Coords( &u, &v, GSM_TO_SM, LstarInfo->mInfo->c );
//        r    = Lgm_Magnitude( &v );
//        mlat = asin( v.z/r )*DegPerRad;




        /*
         * Save individual I values
         */
        LstarInfo->I[k] = Ifound;


        MirrorMLT[k]  = MLT;
        MirrorMlat[k] = mlat;


        /*
         *  convert mirror point to GSM.
         */
//why are we doing this?
// why does this seem to change?
//        Phi = 15.0*(MLT-12.0)*RadPerDeg;
//        cl = cos( mlat * RadPerDeg ); sl = sin( mlat * RadPerDeg );
//        u.x = r*cl*cos(Phi); u.y = r*cl*sin(Phi); u.z = r*sl;
//        Lgm_Convert_Coords( &u, &v, SM_TO_GSM, LstarInfo->mInfo->c );

        /*
         * Save GSM cartesian as well...
         */
        v = LstarInfo->mInfo->Pm_North;
        LstarInfo->Mirror_Pn[k] = v;
        LstarInfo->Mirror_Sn[k] = LstarInfo->mInfo->Sm_North;

        LstarInfo->Mirror_Ps[k] = LstarInfo->mInfo->Pm_South;
        LstarInfo->Mirror_Ss[k] = LstarInfo->mInfo->Sm_South;

        /*
         *  Trace to earth to get the footpoint. Note that we are trying to
         *  compute L* here, and to do that, we eventually need to integrate B
         *  dot dA to get magnetic flux. We could trace down to the ellipsoid,
         *  but that would complicate the integral. Its much easier to do the
         *  integral on a sphere. So, instead of tracing to ellipsoid, we will
         *  trace to the spherical approx to the Earth instead. There is no
         *  loss of generality in doing this when calculating L*, but we need
         *  to take note that the footpoints obtained are not relative to the
         *  elipsoid.
         */
        LstarInfo->mInfo->Hmax = 0.1;
//        if ( !Lgm_TraceToSphericalEarth( &v, &w, LstarInfo->mInfo->Lgm_LossConeHeight, 1.0, 1e-7, LstarInfo->mInfo ) ){ return(-4); }
        if ( !Lgm_TraceToSphericalEarth( &v, &w, LstarInfo->mInfo->Lgm_LossConeHeight, 1.0, 1e-11, LstarInfo->mInfo ) ){ return(-4); }
        LstarInfo->Spherical_Footprint_Pn[k] = w;

        /*
         *  convert footpoint back to SM
         */
        Lgm_Convert_Coords( &w, &v, GSM_TO_SM, LstarInfo->mInfo->c );
        Phi = atan2( v.y, v.x );
        LstarInfo->MLT[k]  = Phi*DegPerRad/15.0 + 12.0;
        LstarInfo->mlat[k] = asin( v.z/Lgm_Magnitude(&v) )*DegPerRad;
        if (LstarInfo->VerbosityLevel > 2)  { printf(" \t\t\t%sMLT_foot, mlat_foot = %g %g%s\n\n", PreStr, LstarInfo->MLT[k], LstarInfo->mlat[k], PostStr); fflush(stdout); }


        /*
         * If SaveShellLines is set true, then retrace the FL and save the
         * whole FL.  Note that since we appear to only have the north
         * footpoint, we start there and trace to south. So lets pack them in
         * the saved arrays backwards so that they go from south to north.
         */
LstarInfo->FindShellPmin = TRUE;
        if ( LstarInfo->FindShellPmin || LstarInfo->ComputeVgc ) {
            //Lgm_TraceToMinBSurf( &LstarInfo->Spherical_Footprint_Pn[k], &v2, 0.1, 1e-8, LstarInfo->mInfo );
            Lgm_TraceToMinBSurf( &LstarInfo->Spherical_Footprint_Pn[k], &v2, 0.1, 1e-8, LstarInfo->mInfo );
            LstarInfo->mInfo->Bfield( &v2, &LstarInfo->Bmin[k], LstarInfo->mInfo );
            LstarInfo->Pmin[k] = v2;
        }

        /*
         * Compute the gradient of I, Sb and Vgc
         */
        if ( LstarInfo->ComputeVgc ) {
            // Lgm_Grad_I() and other routines may modify mInfo in undesirable ways, so give it a copy.
            mInfo2 = Lgm_CopyMagInfo( LstarInfo->mInfo );

            mInfo2->FirstCall = TRUE;
            mInfo2->Lgm_n_Sb_integrand_Calls = 0;
            mInfo2->Lgm_Sb_Integrator_epsabs = 0.0;
            mInfo2->Lgm_Sb_Integrator_epsrel = 1e-3;
            double Sb = SbIntegral( mInfo2 );
Sb = 1.0;

            Lgm_Grad_I( &LstarInfo->Pmin[k], &LstarInfo->GradI[k], mInfo2 );

            LstarInfo->mInfo->Bfield( &LstarInfo->Pmin[k], &Bvec, LstarInfo->mInfo );
            printf("\t\tB = %g %g %g\n", Bvec.x, Bvec.y, Bvec.z);
            B = Lgm_NormalizeVector( &Bvec );   // nT
            B *= 1e-9;                          // T
            printf("\t\tB = %g T\n", B);

            Lgm_FreeMagInfo( mInfo2 );

            double K = LstarInfo->KineticEnergy;       // keV
K = 1000.0; //keV
            K *= 1e3*LGM_e;                        // Joules
            double M = LstarInfo->Mass;                // kg
M = ELECTRON_MASS; // kg


            /*
             * Note: Kinetic energy is difference between total E and rest Energy
             * In non-rel. limit, the factor g=2. So ration of Non Rel to Rel velocity is
             * (2Eta+2)/(Eta+2)
             */
            double Eta = K/(M*CC*CC);
            double g = K*(2.0+Eta)/(1.0+Eta);


            double q = 1.0;            // units of elementary charge
            q *= 1.6021e-19;    // Coulombs
            Sb *= 6371e3;       // m
            double f = g/(q*Sb*B)/1000.0;

            Lgm_Vector Vcg;
            Lgm_CrossProduct( &LstarInfo->GradI[k], &Bvec, &Vcg );
            Vcg.x *= f; Vcg.y *= f; Vcg.z *= f;
            printf("\t\tEta = %g    <Vcg>_non./<vcg>_rel. = %g (%g)\n", Eta, (2.0*Eta+2.0)/(Eta+2.0), (Eta+2.0)/(2.0*Eta+2.0));
            printf("\t\t<Vcg> = %g, %g, %g    km/s\n\n\n", Vcg.x, Vcg.y, Vcg.z);
            LstarInfo->Vgc[k] = Vcg;

        }

        if ( LstarInfo->SaveShellLines ) {

            Lgm_TraceLine( &LstarInfo->Spherical_Footprint_Pn[k], &v2, LstarInfo->mInfo->Lgm_LossConeHeight, -1.0, 1e-8, FALSE, LstarInfo->mInfo );
            LstarInfo->Spherical_Footprint_Ps[k] = v2;
//FILE *fppp;
//double AlphaEq;
//fppp = fopen("FL.txt","a");
//AlphaEq = asin( sqrt( Lgm_Magnitude( &LstarInfo->Bmin[k]) /LstarInfo->mInfo->Bm ) );
//for (i=0; i<LstarInfo->mInfo->nPnts; i++){
////fprintf(fppp, "%g %g\n", LstarInfo->mInfo->s[i], LstarInfo->mInfo->Bmag[i]);
//fprintf(fppp, "%g:     %g %g %g %g\n", LstarInfo->MLT[k], LstarInfo->mInfo->Px[i], LstarInfo->mInfo->Py[i], LstarInfo->mInfo->Pz[i], AlphaEq*DegPerRad );
//}
//fclose(fppp);

//Type = ClassifyFL( k, LstarInfo );
//printf("k, Type = %d %d\n", k, Type);
            LstarInfo->nMinMax = k;

            nnn = LstarInfo->mInfo->nPnts; smax = LstarInfo->mInfo->s[nnn-1];
            for (tkk=0, nfp=nnn-1; nfp>=0; nfp--){
                if ( LstarInfo->mInfo->Bmag[nfp] > 0.0 ) {
                    LstarInfo->s_gsm[k][tkk] = smax - LstarInfo->mInfo->s[nfp];
                    LstarInfo->Bmag[k][tkk]  = LstarInfo->mInfo->Bmag[nfp];
                    LstarInfo->x_gsm[k][tkk] = LstarInfo->mInfo->Px[nfp];
                    LstarInfo->y_gsm[k][tkk] = LstarInfo->mInfo->Py[nfp];
                    LstarInfo->z_gsm[k][tkk] = LstarInfo->mInfo->Pz[nfp];
                    ++tkk;
                }
            }
            LstarInfo->nFieldPnts[k] = tkk;



            LstarInfo->Spherical_Footprint_Ss[k] = LstarInfo->s_gsm[k][0];
            LstarInfo->Spherical_Footprint_Sn[k] = smax;


            /*
             * Find true footpoints (i.e. relative to ellipsoid), we may want to do an additional trace here...
             */
            Hmax = LstarInfo->mInfo->Hmax;
            LstarInfo->mInfo->Hmax = 0.001;
            //if ( Lgm_TraceToEarth( &LstarInfo->Spherical_Footprint_Ps[k], &LstarInfo->Ellipsoid_Footprint_Ps[k], LstarInfo->mInfo->Lgm_LossConeHeight, -1.0, 1e-7, LstarInfo->mInfo ) ) {
            if ( Lgm_TraceToEarth( &LstarInfo->Spherical_Footprint_Ps[k], &LstarInfo->Ellipsoid_Footprint_Ps[k], LstarInfo->mInfo->Lgm_LossConeHeight, -1.0, 1e-7, LstarInfo->mInfo ) ) {

                LstarInfo->Ellipsoid_Footprint_Ss[k] = LstarInfo->Spherical_Footprint_Ss[k] - LstarInfo->mInfo->Trace_s; // should be slightly negative

                LstarInfo->mInfo->Hmax = 0.001;
                if ( Lgm_TraceToEarth( &LstarInfo->Spherical_Footprint_Pn[k], &LstarInfo->Ellipsoid_Footprint_Pn[k], LstarInfo->mInfo->Lgm_LossConeHeight, 1.0, 1e-7, LstarInfo->mInfo ) ) {

                    LstarInfo->Ellipsoid_Footprint_Sn[k] = LstarInfo->Spherical_Footprint_Sn[k] + LstarInfo->mInfo->Trace_s;

                }

            }
            LstarInfo->mInfo->Hmax = Hmax;

            /*
             *  So far, the way we have done this, we have never really needed
             *  to know the distance of the mirror points from the footpoints.
             *  To compute these, trace from Pm_south back to southern
             *  spherical footpoint. Then add this distance to Sm_South and
             *  Sm_North. That should give the distance of both mirror points
             *  rtelativen to the southern spherical footpoint -- which is how
             *  the field line is defined that we are trying to save.
             */
            //if ( Lgm_TraceToSphericalEarth( &LstarInfo->Mirror_Ps[k], &uu, LstarInfo->mInfo->Lgm_LossConeHeight, -1.0, 1e-7, LstarInfo->mInfo ) ){
             if ( Lgm_TraceToSphericalEarth( &LstarInfo->Mirror_Ps[k], &uu, LstarInfo->mInfo->Lgm_LossConeHeight, -1.0, 1e-7, LstarInfo->mInfo ) ) {
                LstarInfo->Mirror_Ss[k] += LstarInfo->mInfo->Trace_s;
                LstarInfo->Mirror_Sn[k] += LstarInfo->mInfo->Trace_s;
            }





        }



        ++k;

    } // end MLT for loop
    LstarInfo->nPnts = k;


    /*
     *  Save drift shell -- it will help us predict the next one.
     */
    for (i=0; i<LstarInfo->nPnts; ++i){
	    MirrorMLT_Old[i]  = MirrorMLT[i];
	    MirrorMlat_Old[i] = MirrorMlat[i];
    }

    /*
     *  To get Lstar all we need to do now is one final integral.
     *  See roederer's eq 3.6 It should be trivial at this point!
     */
//printf("LstarInfo->nPnts = %d\n", LstarInfo->nPnts);
//for(i=0; i<LstarInfo->nPnts; i++){
//printf("HERE 1a: LstarInfo->MLT[%d] = %g  LstarInfo->mlat[%d] = %g    \n", i, LstarInfo->MLT[i], i, LstarInfo->mlat[i] );
//}
//for(i=0; i<LstarInfo->nPnts; i++){
//printf("HERE 1b: MirrorMLT[%d] = %g  MirrorMlat[%d] = %g    \n", i, MirrorMLT[i], i, MirrorMlat[i] );
//}


    /*
     *  sort arrays so they are monotonically increasing in MLT
     *  (needed for spline interp).
     */
    quicksort2(  LstarInfo->nPnts, LstarInfo->MLT-1, LstarInfo->mlat-1 );



    /*
     *  Create new arrays that have the last two points repeated at the start
     *  and the first two points repeated at the end. We have to add/subtract
     *  24 hours at the end/start. The repeated values ensure that the spline
FIX
     *  will be periodic and continous around the globe.
     */
/*
    j = 0;
    for (i=0; i<LstarInfo->nPnts; ++i){ LstarInfo->xa[j] = LstarInfo->MLT[i]-24.0; LstarInfo->ya[j] = LstarInfo->mlat[i]; ++j; }
    for (i=0; i<LstarInfo->nPnts; ++i){ LstarInfo->xa[j] = LstarInfo->MLT[i]     ; LstarInfo->ya[j] = LstarInfo->mlat[i]; ++j; }
    for (i=0; i<LstarInfo->nPnts; ++i){ LstarInfo->xa[j] = LstarInfo->MLT[i]+24.0; LstarInfo->ya[j] = LstarInfo->mlat[i]; ++j; }
    spline( LstarInfo->xa, LstarInfo->ya, LstarInfo->nSplnPnts, 0.0, 0.0, LstarInfo->y2 );
*/

    j = 0;
    for (i=0; i<LstarInfo->nPnts; ++i){ LstarInfo->xa[j] = LstarInfo->MLT[i]-24.0; LstarInfo->ya[j] = LstarInfo->mlat[i]; ++j; }
    for (i=0; i<LstarInfo->nPnts; ++i){ LstarInfo->xa[j] = LstarInfo->MLT[i]     ; LstarInfo->ya[j] = LstarInfo->mlat[i]; ++j; }
    for (i=0; i<LstarInfo->nPnts; ++i){ LstarInfo->xa[j] = LstarInfo->MLT[i]+24.0; LstarInfo->ya[j] = LstarInfo->mlat[i]; ++j; }
    LstarInfo->nSplnPnts = j;


    /*
     *  Use GSL periodic spline (can also try others if needed....)
     */
//    gsl_set_error_handler_off(); // Turn off gsl default error handler
    LstarInfo->acc = gsl_interp_accel_alloc( );
//    LstarInfo->pspline = gsl_interp_alloc( gsl_interp_cspline_periodic, LstarInfo->nSplnPnts );
    LstarInfo->pspline = gsl_interp_alloc( gsl_interp_cspline_periodic, LstarInfo->nSplnPnts );
//    LstarInfo->pspline = gsl_interp_alloc( gsl_interp_linear, LstarInfo->nSplnPnts );
    //printf ("spline uses '%s' interpolation.\n", gsl_interp_name( LstarInfo->pspline ));
    gsl_interp_init( LstarInfo->pspline, LstarInfo->xa, LstarInfo->ya, LstarInfo->nSplnPnts );

//for(i=0; i<LstarInfo->m; i++){
//printf("HERE 2: LstarInfo->xa[%d] = %g  LstarInfo->ya[%d] = %g    \n", i, LstarInfo->xa[i], i, LstarInfo->ya[i] );
//}


    Phi1 = MagFlux( LstarInfo );
    LstarInfo->LS_dip_approx = -2.0*M_PI*LstarInfo->Mused/Phi1;
    Phi2 = MagFlux2( LstarInfo );
    LstarInfo->LS = -2.0*M_PI*LstarInfo->Mused /Phi2;
    LstarInfo->LS_McIlwain_M = -2.0*M_PI*LstarInfo->mInfo->c->M_cd_McIllwain /Phi2;




    /*
     *  Determine the type of the orbit. We will not be here if we bailed early.. So LGM_DRIFT_ORBIT_OPEN_SHABANSKY
     *  cant be determined here ... 
     */
    LstarInfo->DriftOrbitType = LGM_DRIFT_ORBIT_CLOSED;
    for ( i=0; i<LstarInfo->nPnts; ++i ) if ( LstarInfo->nMinima[i] > 1 ) LstarInfo->DriftOrbitType = LGM_DRIFT_ORBIT_CLOSED_SHABANSKY;
    if ( LstarInfo->LS < 0.0 ) LstarInfo->DriftOrbitType = LGM_DRIFT_ORBIT_OPEN;





    if (LstarInfo->VerbosityLevel > 1) {
        printf("\n\t\t%sL*, Dipole Approximation.\n%s", PreStr, PostStr );
        printf("\t\t%s  Magnetic Flux:                         %.15lf%s\n", PreStr, Phi1, PostStr );
        printf("\t\t%s  L*:                                    %.15lf%s\n", PreStr, LstarInfo->LS_dip_approx, PostStr );
        printf("\t\t%s  L* (Using McIllwain M):                %.15lf%s\n", PreStr, -2.0*M_PI*LstarInfo->mInfo->c->M_cd_McIllwain /Phi1, PostStr );
        printf("\n\t\t%sL*, Full Field.%s\n", PreStr, PostStr );
        printf("\t\t%s  Magnetic Flux:                         %.15lf%s\n", PreStr, Phi2, PostStr );
        printf("\t\t%s  L*:                                    %.15lf%s\n", PreStr, LstarInfo->LS, PostStr );
        printf("\t\t%s  L* (Using McIllwain M):                %.15lf%s\n", PreStr, -2.0*M_PI*LstarInfo->mInfo->c->M_cd_McIllwain /Phi2, PostStr );
        if      ( LstarInfo->DriftOrbitType == LGM_DRIFT_ORBIT_CLOSED )          printf("\n\t\t%sDrift Orbit Type: CLOSED%s\n", PreStr, PostStr );
        else if ( LstarInfo->DriftOrbitType == LGM_DRIFT_ORBIT_CLOSED_SHABANSKY) printf("\n\t\t%sDrift Orbit Type: SHABANSKY%s\n", PreStr, PostStr );
        else if ( LstarInfo->DriftOrbitType == LGM_DRIFT_ORBIT_OPEN )            printf("\n\t\t%sDrift Orbit Type: OPEN%s\n", PreStr, PostStr );
        printf("\n\n\n" );
    }

    gsl_interp_free( LstarInfo->pspline );
    gsl_interp_accel_free( LstarInfo->acc );


    return( (LstarInfo->LS < 0.0) ?  -1 : 1 );

}



double MagFluxIntegrand( double Phi, _qpInfo *qpInfo ) {


    double	MLT, mlat, c;
    Lgm_LstarInfo  *LstarInfo;


    MLT = Phi*DegPerRad/15.0;


    /*
     *  Get pointer to our auxilliary data structure.
     */
    LstarInfo = (Lgm_LstarInfo *)qpInfo;

//    gsl_set_error_handler_off(); // Turn off gsl default error handler
    mlat = gsl_interp_eval( LstarInfo->pspline, LstarInfo->xa, LstarInfo->ya, MLT, LstarInfo->acc );
    //splint( LstarInfo->xa, LstarInfo->ya, LstarInfo->y2, LstarInfo->nSplnPnts, MLT, &mlat);

    c = cos( mlat*RadPerDeg );

    return( c*c );

}



double MagFlux( Lgm_LstarInfo *LstarInfo ) {

    double      a, b, r;
    double      epsabs, epsrel, result, abserr;
    int         key, neval, ier, limit, lenw, last, iwork[502];
    double      work[2002];
    _qpInfo     *qpInfo;


    /*
     *  Type-cast our data structure to a generic type.
     *  The structure holds auzilliary info we need down
     *  in the function calls.
     */
    qpInfo = (_qpInfo *)LstarInfo;



    /*
     *  Limits of integration.
     */
    a = 0.0;
    b = 2.0*M_PI;



    /*
     *   set tolerances.
     */
    //epsabs = 0.0;
    //epsrel = 1e-5;
    epsabs = LstarInfo->mInfo->Lgm_MagFlux_Integrator_epsabs;
    epsrel = LstarInfo->mInfo->Lgm_MagFlux_Integrator_epsrel;



    limit = 500; lenw = 4*limit; key = 6;
/*
    iwork  = (int *) calloc( limit+1, sizeof(int) );
    work   = (double *) calloc( lenw+1, sizeof(double) );
*/
    dqags(MagFluxIntegrand, qpInfo, a, b, epsabs, epsrel, &result, &abserr, &neval, &ier, limit, lenw, &last, iwork, work, LstarInfo->mInfo->VerbosityLevel );
/*
    free( iwork );
    free( work );
*/


    r = 1.0 + LstarInfo->mInfo->Lgm_LossConeHeight/WGS84_A;

    return( -result*LstarInfo->Mused/r );

}





/*
 *  Lets integrate the real field. I.e. dont assume a dipole field
 */
double LambdaIntegrand( double Lambda, _qpInfo *qpInfo ) {


    double	cl, sl, st, ct, sp, cp, Br, Bx, By, Bz, f;
    double	MLT, phi;
    Lgm_Vector	u, w, Bvec;
    Lgm_LstarInfo  *LstarInfo;

    /*
     *  Get pointer to our auxilliary data structure.
     */
    LstarInfo = (Lgm_LstarInfo *)qpInfo;

    MLT = LstarInfo->Phi*DegPerRad/15.0;
    phi = 15.0*(MLT-12.0)*RadPerDeg;
    cl = cos( Lambda ); sl = sin( Lambda );
    u.x = cl*cos( phi );
    u.y = cl*sin( phi );
    u.z = sl;
    Lgm_Convert_Coords( &u, &w, SM_TO_GSM, LstarInfo->mInfo->c );


    LstarInfo->mInfo->Bfield( &w, &Bvec, LstarInfo->mInfo );
    Lgm_Convert_Coords( &Bvec, &u, GSM_TO_SM, LstarInfo->mInfo->c );
    Bx = u.x; By = u.y; Bz = u.z;

    st = sin( M_PI/2.0 - Lambda ); ct = cos( M_PI/2.0 - Lambda ); sp = sin( phi ); cp = cos( phi );
    Br = st*cp*Bx + st*sp*By + ct*Bz;

    f = Br*cos( Lambda );

    return( f );

}


double LambdaIntegral( Lgm_LstarInfo *LstarInfo ) {

    double      a, b;
    double      epsabs, epsrel, result, abserr;
    int         key, neval, ier, limit, lenw, last, iwork[502];
    double      work[2002], MLT, mlat;
    _qpInfo     *qpInfo;


    /*
     *  Type-cast our data structure to a generic type.
     *  The structure holds auzilliary info we need down
     *  in the function calls.
     */
    qpInfo = (_qpInfo *)LstarInfo;

    /*
     *  Get Latitude
     */
    MLT = LstarInfo->Phi*DegPerRad/15.0;
//    gsl_set_error_handler_off(); // Turn off gsl default error handler
//int i;
//for(i=0; i<LstarInfo->m; i++){
//printf("LstarInfo->xa[%d] = %g  LstarInfo->ya[%d] = %g \n", i, LstarInfo->xa[i], i, LstarInfo->ya[i] );
//}
    mlat = gsl_interp_eval( LstarInfo->pspline, LstarInfo->xa, LstarInfo->ya, MLT, LstarInfo->acc );
    //splint( LstarInfo->xa, LstarInfo->ya, LstarInfo->y2, LstarInfo->nSplnPnts, MLT, &mlat);
    mlat *= RadPerDeg;


    /*
     *  Limits of integration for lambda integral.
     */
    a = mlat;
    b = M_PI/2.0;


    /*
     *   set tolerances.
     */
    //epsabs = 0.0;
    //epsrel = 1e-3;
    epsabs = LstarInfo->mInfo->Lgm_MagFlux_Integrator_epsabs;
    epsrel = LstarInfo->mInfo->Lgm_MagFlux_Integrator_epsrel;

    limit = 500; lenw = 4*limit; key = 6;
/*
    iwork  = (int *) calloc( limit+1, sizeof(int) );
    work   = (double *) calloc( lenw+1, sizeof(double) );
*/
//printf("a, b = %g %g\n", a, b);
    dqags(LambdaIntegrand, qpInfo, a, b, epsabs, epsrel, &result, &abserr, &neval, &ier, limit, lenw, &last, iwork, work, LstarInfo->mInfo->VerbosityLevel );
/*
    free( iwork );
    free( work );
*/



    return( result );

}

double MagFluxIntegrand2( double Phi, _qpInfo *qpInfo ) {


    double	f;
    Lgm_LstarInfo  *LstarInfo;

    /*
     *  Get pointer to our auxilliary data structure.
     */
    LstarInfo = (Lgm_LstarInfo *)qpInfo;

    LstarInfo->Phi = Phi;
    f = LambdaIntegral( LstarInfo );

    return( f );

}

double MagFlux2( Lgm_LstarInfo *LstarInfo ) {

    double      a, b, r;
    double      epsabs, epsrel, result, abserr;
    int         key, neval, ier, limit, lenw, last, iwork[502];
    double      work[2002];
    _qpInfo     *qpInfo;


    /*
     *  Type-cast our data structure to a generic type.
     *  The structure holds auzilliary info we need down
     *  in the function calls.
     */
    qpInfo = (_qpInfo *)LstarInfo;



    /*
     *  Limits of integration.
     */
    a = 0.0;
    b = 2.0*M_PI;


//int i;
//for(i=0; i<LstarInfo->m; i++){
//printf("MagFlux2: LstarInfo->xa[%d] = %g  LstarInfo->ya[%d] = %g \n", i, LstarInfo->xa[i], i, LstarInfo->ya[i] );
//}

    /*
     *   set tolerances.
     */
    epsabs = LstarInfo->mInfo->Lgm_LambdaIntegral_Integrator_epsabs;
    epsrel = LstarInfo->mInfo->Lgm_LambdaIntegral_Integrator_epsrel;

    limit = 500; lenw = 4*limit; key = 6;
/*
    iwork  = (int *) calloc( limit+1, sizeof(int) );
    work   = (double *) calloc( lenw+1, sizeof(double) );
*/
    dqags(MagFluxIntegrand2, qpInfo, a, b, epsabs, epsrel, &result, &abserr, &neval, &ier, limit, lenw, &last, iwork, work, LstarInfo->mInfo->VerbosityLevel );
/*
    free( iwork );
    free( work );
*/


    r = 1.0 + LstarInfo->mInfo->Lgm_LossConeHeight/WGS84_A;

    return( result/r );

}


/*
 * Given a history of MirrorMLT[] and MirrorMlat[] vals, this routine tries to
 * predict what the next one will be via rational function extrapolation.
 */
void PredictMlat1( double *MirrorMLT, double *MirrorMlat, int k, double MLT, double *pred_mlat, double *pred_delta_mlat, double *delta ) {

    int		nk, koffset;
    double	pred_mlat1, pred_mlat2, pred_delta_mlat1, pred_delta_mlat2;


    if ( k == 0 ){

        /*
         *  Really shouldnt be here. We should realy lonly be called
         *  with some data to predict with
         */
        *pred_mlat  	 = 45.0;
        *pred_delta_mlat = 45.0;
        *delta 		 = 45.0;

    } else if ( k == 1 ) {

        *pred_mlat 	 = MirrorMlat[0];
        *pred_delta_mlat = 5.0;
        *delta 		 = 5.0;

    } else if ( k == 2 ) {

        *pred_mlat 	 = 0.5*(MirrorMlat[0]+MirrorMlat[1]);
        *pred_delta_mlat = 2.5;
        *delta 		 = 2.5;

    } else {

        if ( k > 4 ) {
	        koffset = k-5;
	        nk = 5;
        } else {
	        koffset = 0;
	        nk = k;
        }


        Lgm_PolFunInt( MirrorMLT+koffset-1, MirrorMlat+koffset-1, nk, MLT, &pred_mlat2, &pred_delta_mlat2 );
	    Lgm_RatFunInt( MirrorMLT+koffset-1, MirrorMlat+koffset-1, nk, MLT, &pred_mlat1, &pred_delta_mlat1 );

       	/*
	     *  Make certain we have sane predictions
	     */
        if ( !isnan(pred_mlat1) && !isnan(pred_mlat2) ) {
            *pred_mlat = 0.5*(pred_mlat1 + pred_mlat2);
            *pred_delta_mlat = fabs(pred_delta_mlat1) + fabs(pred_delta_mlat2);
        } else if ( !isnan(pred_mlat1) ) {
            *pred_mlat = pred_mlat1;
            *pred_delta_mlat = fabs(pred_delta_mlat1);
        } else if ( !isnan(pred_mlat2) ) {
            *pred_mlat = pred_mlat2;
            *pred_delta_mlat = fabs(pred_delta_mlat2);
        } else {
            *pred_mlat = MirrorMlat[0];
            *pred_delta_mlat = 2.0;
        }

        if ( ( *pred_mlat < 0.0) || ( *pred_mlat > 90.0 ) ){
            *pred_mlat = MirrorMlat[0];
            *pred_delta_mlat = 2.0;
            }

        }


        if ( (*pred_delta_mlat < 1e-12) || isnan(*pred_delta_mlat) ) {
            *delta = 2.0;
        } else {
            *delta = 1.5*(*pred_delta_mlat);
        }

        if ( *delta < .2 ) *delta = 0.2;


}


/*
 * Given a history of MirrorMLT[] and MirrorMlat[] vals, this routine tries to
 * predict what the next one will be via a periodic akima spine.
 */
void PredictMlat2( double *MirrorMLT, double *MirrorMlat, int k, double MLT, double *pred_mlat, double *pred_delta_mlat, double *delta, Lgm_LstarInfo *LstarInfo ) {

    double	mlat;
    int 	j, i;



    if ( k == 0 ){

	    /*
	     *  Really shouldnt be here. We should realy lonly be called
	     *  with some data to predict with
	     */
	    *pred_mlat  	 = 45.0;
	    *pred_delta_mlat = 45.0;
	    *delta 		 = 45.0;

    } else if ( k == 1 ) {

	    *pred_mlat 	 = MirrorMlat[0];
	    *pred_delta_mlat = 5.0;
	    *delta 		 = 5.0;

    } else if ( k == 2 ) {

	    *pred_mlat 	 = 0.5*(MirrorMlat[0]+MirrorMlat[1]);
	    *pred_delta_mlat = 2.5;
	    *delta 		 = 2.5;

    } else {

        /*
         *  Create new arrays that have the last two points repeated at the start
         *  and the first two points repeated at the end. We have to add/subtract
         *  24 hours at the end/start. The repeated values ensure that the spline
         *  will be periodic and continous around the globe.
         */

	    /*
	     *  sort arrays so they are monotonically increasing in MLT
	     *  (needed for spline interp).
	     */
	    quicksort2(  k, MirrorMLT-1, MirrorMlat-1 );

/*
        j = 0;
        LstarInfo->xma[j] = MirrorMLT[k-2]-24.0; LstarInfo->yma[j] = MirrorMlat[k-2]; ++j;
        LstarInfo->xma[j] = MirrorMLT[k-1]-24.0; LstarInfo->yma[j] = MirrorMlat[k-1]; ++j;
        for (i=0; i<k; ++i){ LstarInfo->xma[j] = MirrorMLT[i]; LstarInfo->yma[j] = MirrorMlat[i]; ++j; }
        LstarInfo->xma[j] = MirrorMLT[0]+24.0; LstarInfo->yma[j] = MirrorMlat[0]; ++j;
        LstarInfo->xma[j] = MirrorMLT[1]+24.0; LstarInfo->yma[j] = MirrorMlat[1]; ++j;
        LstarInfo->m = j;
*/

        j = 0;
        for (i=0; i<k; ++i){ LstarInfo->xma[j] = MirrorMLT[i]-24.0; LstarInfo->yma[j] = MirrorMlat[i]; ++j; }
        for (i=0; i<k; ++i){ LstarInfo->xma[j] = MirrorMLT[i]     ; LstarInfo->yma[j] = MirrorMlat[i]; ++j; }
        for (i=0; i<k; ++i){ LstarInfo->xma[j] = MirrorMLT[i]+24.0; LstarInfo->yma[j] = MirrorMlat[i]; ++j; }
	    LstarInfo->m = j;

//for (i=0; i<j; ++i){
//    printf("%g %g\n", LstarInfo->xma[i], LstarInfo->yma[i] );
//}

        gsl_set_error_handler_off(); // Turn off gsl default error handler
	    gsl_interp_accel *acc = gsl_interp_accel_alloc( );
	    gsl_interp *pspline = gsl_interp_alloc( gsl_interp_akima_periodic, LstarInfo->m );
	    //gsl_interp *pspline = gsl_interp_alloc( gsl_interp_akima, LstarInfo->m );
	    //gsl_interp *pspline = gsl_interp_alloc( gsl_interp_cspline_periodic, LstarInfo->m );
	    gsl_interp_init( pspline, LstarInfo->xma, LstarInfo->yma, LstarInfo->m );
	    mlat = gsl_interp_eval( pspline, LstarInfo->xma, LstarInfo->yma, MLT, acc );
	    gsl_interp_free( pspline );
	    gsl_interp_accel_free( acc );

	    /* old num rec. stuff
        spline( LstarInfo->xma, LstarInfo->yma, LstarInfo->m, 0.0, 0.0, LstarInfo->ym2 );
        splint( LstarInfo->xma, LstarInfo->yma, LstarInfo->ym2, LstarInfo->m, MLT, &mlat);
	    */


        *pred_mlat       = mlat;
        *pred_delta_mlat = 0.75;
        *delta           = 0.75;

    }



}




/*
 *      Input Variables:
 *
 *                      Date: 
 *                       UTC: 
 *                     brac1:  Radial distance for inner edge of search bracket
 *                     brac2:  Radial distance for outer edge of search bracket
 *                     Alpha:  Equatorial pitch angle to compute 
 *                   Quality:  Quality factor (0-8)
 *
 *      Input/OutPut Variables:
 *                         K:  K value for LCDS at given alpha
 *                 LstarInfo:  LstarInfo structure to specify B-field model, store last valid Lstar, etc.
 *  
 */
int Lgm_LCDS( long int Date, double UTC, double brac1, double brac2, double Kin, double LT, double tol, int Quality, int nFLsInDriftShell, double *K, Lgm_LstarInfo *LstarInfo ) {
    Lgm_LstarInfo   *LstarInfo_brac1, *LstarInfo_brac2, *LstarInfo_test;
    Lgm_Vector      v1, v2, v3, LCDSv3, Bvec, Ptest, PtestSM, Pinner, PinnerSM, Pouter, PouterSM;
    Lgm_DateTime    DT_UTC;
    double          Blocal, Xtest, sa, sa2, LCDS, Alpha;
    double          nTtoG = 1.0e-5;
    int             LS_Flag, TFlag, nn, k;
    int             maxIter = 20;

    //Start by creating necessary structures... need an Lstar info for each bracket, plus test point.
    Lgm_Make_UTC( Date, UTC, &DT_UTC, LstarInfo->mInfo->c );
    LstarInfo_brac1 = Lgm_CopyLstarInfo( LstarInfo );
    LstarInfo_brac2 = Lgm_CopyLstarInfo( LstarInfo );
    LstarInfo_test = Lgm_CopyLstarInfo( LstarInfo );

    // Set Tolerances
    Lgm_SetLstarTolerances( Quality, nFLsInDriftShell, LstarInfo_brac1 );
    Lgm_SetLstarTolerances( Quality, nFLsInDriftShell, LstarInfo_brac2 );
    Lgm_SetLstarTolerances( Quality, nFLsInDriftShell, LstarInfo_test );

    // set coord transformation 
    Lgm_Set_Coord_Transforms( Date, UTC, LstarInfo_brac1->mInfo->c );
    Lgm_Set_Coord_Transforms( Date, UTC, LstarInfo_brac2->mInfo->c );
    Lgm_Set_Coord_Transforms( Date, UTC, LstarInfo_test->mInfo->c );

    //Check for closed drift shell at brackets
    LT *= 15.0;
    LT *= RadPerDeg;
    PinnerSM.x = brac1*cos(LT); PinnerSM.y = brac1*sin(LT); PinnerSM.z = 0.0;
    Lgm_Convert_Coords( &PinnerSM, &Pinner, SM_TO_GSM, LstarInfo_brac1->mInfo->c );

    /*
     *  Test inner bracket location.
     */

    //Trace to minimum-B
    if ( Lgm_Trace( &Pinner, &v1, &v2, &v3, 120.0, 0.01, TRACE_TOL, LstarInfo_brac1->mInfo ) == LGM_CLOSED ) {

        if ( Lgm_Setup_AlphaOfK( &DT_UTC, &v3, LstarInfo_brac1->mInfo ) > 0 ) {
            Alpha = Lgm_AlphaOfK( Kin, LstarInfo_brac1->mInfo );
            LstarInfo_brac1->PitchAngle = Alpha;
            Lgm_TearDown_AlphaOfK(LstarInfo_brac1->mInfo);
        } else {
            return(-8);
        }
        if (LstarInfo->VerbosityLevel > 0) printf("[Inner bracket] Alpha (of K) is %g (%g)\n", Alpha, Kin);
        sa = sin( LstarInfo_brac1->PitchAngle*RadPerDeg ); sa2 = sa*sa;
        LstarInfo_brac1->mInfo->Bm = LstarInfo_brac1->mInfo->Bmin/sa2;
        if (LstarInfo->VerbosityLevel > 0) printf("[Inner bracket] Bm, Bmin = is %g, %g\n", LstarInfo_brac1->mInfo->Bm, LstarInfo_brac1->mInfo->Bmin);
        //Only continue if bracket 1 is closed FL
        //Get L*
        LS_Flag = Lstar( &v3, LstarInfo_brac1);

        if (LstarInfo_brac1->LS <= 0) {
            if (LstarInfo->VerbosityLevel > 0) printf("Undefined DS at inner bracket. L* = %g (Flag = %d); P = %g, %g, %g\n",
                                                       LstarInfo_brac1->LS, LS_Flag, Pinner.x, Pinner.y, Pinner.z);
            FreeLstarInfo( LstarInfo_brac1 );
            FreeLstarInfo( LstarInfo_brac2 );
            FreeLstarInfo( LstarInfo_test );
            return(-8); //Inner bracket is bad - bail
        }
        LCDSv3 = v3;
        LCDS = LstarInfo_brac1->LS;
        *K = (LstarInfo_brac1->I0)*sqrt(Lgm_Magnitude(&LstarInfo_brac1->Bmin[0]));
        if (LstarInfo->VerbosityLevel > 1) printf("Found valid inner bracket. Pinner, Pmin_GSM = (%g, %g, %g), (%g, %g, %g)\n", Pinner.x, Pinner.y, Pinner.z, LstarInfo_brac1->mInfo->Pmin.x, LstarInfo_brac1->mInfo->Pmin.y, LstarInfo_brac1->mInfo->Pmin.z);
        if (LstarInfo->VerbosityLevel > 1) printf("Lstar = %g\n", LstarInfo_brac1->LS);

    } else {
        FreeLstarInfo( LstarInfo_brac1 );
        FreeLstarInfo( LstarInfo_brac2 );
        FreeLstarInfo( LstarInfo_test );
        return(-8); //Inner bracket is bad - bail
    }


    /*
     *  Test outer bracket location.
     */
    PouterSM.x = brac2*cos(LT); PouterSM.y = brac2*sin(LT); PouterSM.z = 0.0;
    Lgm_Convert_Coords( &PouterSM, &Pouter, SM_TO_GSM, LstarInfo_brac2->mInfo->c );

    //Trace to minimum-B
    if ( Lgm_Trace( &Pouter, &v1, &v2, &v3, 120.0, 0.01, TRACE_TOL, LstarInfo_brac2->mInfo ) == LGM_CLOSED ) {
        //If bracket 2 is closed FL then check for undefined L*. If L* is defined, we have a problem
        //Get L*
        if ( Lgm_Setup_AlphaOfK( &DT_UTC, &v3, LstarInfo_brac2->mInfo ) > 0 ) {
            Alpha = Lgm_AlphaOfK( Kin, LstarInfo_brac2->mInfo );
            Lgm_TearDown_AlphaOfK(LstarInfo_brac2->mInfo);
        }
        else {
            Alpha = LGM_FILL_VALUE;
        }
        LstarInfo_brac2->PitchAngle = Alpha;
        sa = sin( LstarInfo_brac2->PitchAngle*RadPerDeg ); sa2 = sa*sa;

        LstarInfo_brac2->mInfo->Bm = LstarInfo_brac2->mInfo->Bmin/sa2;
        LS_Flag = Lstar( &v3, LstarInfo_brac2);
        if (LstarInfo_brac2->LS != LGM_FILL_VALUE) {
            //move outer bracket out and try again
            PouterSM.x = 1.7*Lgm_Magnitude(&Pouter)*cos(LT);
            PouterSM.y = 1.7*Lgm_Magnitude(&Pouter)*sin(LT); PouterSM.z = 0.0;
            Lgm_Convert_Coords( &PouterSM, &Pouter, SM_TO_GSM, LstarInfo_brac2->mInfo->c );


            if ( Lgm_Setup_AlphaOfK( &DT_UTC, &v3, LstarInfo_brac2->mInfo ) > 0 ) {
                Alpha = Lgm_AlphaOfK( Kin, LstarInfo_brac2->mInfo );
                Lgm_TearDown_AlphaOfK(LstarInfo_brac2->mInfo);
            }
            else {
                Alpha = LGM_FILL_VALUE;
            }
            LstarInfo_brac2->PitchAngle = Alpha;
            sa = sin( LstarInfo_brac2->PitchAngle*RadPerDeg ); sa2 = sa*sa;
            if ( Lgm_Trace( &Pouter, &v1, &v2, &v3, 120.0, 0.01, TRACE_TOL, LstarInfo_brac2->mInfo ) == LGM_CLOSED ) {
                LstarInfo_brac2->mInfo->Bm = LstarInfo_brac2->mInfo->Bmin/sa2;
                LS_Flag = Lstar( &v3, LstarInfo_brac2);
                if (LstarInfo_brac2->LS != LGM_FILL_VALUE) {
                    FreeLstarInfo( LstarInfo_brac1 );
                    FreeLstarInfo( LstarInfo_brac2 );
                    FreeLstarInfo( LstarInfo_test );
                    return(-9); //Outer bracket is bad - bail
                }
            }
        }
        if (LstarInfo->VerbosityLevel > 0) printf("Found valid outer bracket. Pouter_GSM, Pmin_GSM = (%g, %g, %g), (%g, %g, %g)\n", Pouter.x, Pouter.y, Pouter.z, LstarInfo_brac2->mInfo->Pmin.x, LstarInfo_brac2->mInfo->Pmin.y, LstarInfo_brac2->mInfo->Pmin.z);
    }

    //if brackets are okay, we've moved on without changing anything except setting initial LCDS value as L* at inner bracket
    nn = 0;
    while (Lgm_Magnitude(&Pouter)-Lgm_Magnitude(&Pinner) > tol){
        if (LstarInfo->VerbosityLevel > 2) printf("Current LCDS iteration, bracket width = %d, %g\n", nn, Lgm_Magnitude(&Pouter)-Lgm_Magnitude(&Pinner));
        if (nn > maxIter) {
            printf("********* EXCEEDED MAXITER\n");
            //free structures before returning
            FreeLstarInfo( LstarInfo_brac1 );
            FreeLstarInfo( LstarInfo_brac2 );
            FreeLstarInfo( LstarInfo_test );
            return(-2); //reached max iterations without achieving tolerance - bail
        }

        Xtest = (Lgm_Magnitude(&Pinner) + (Lgm_Magnitude(&Pouter)-Lgm_Magnitude(&Pinner))/2.0);
        PtestSM.x = -1.0*Xtest*cos(LT); PtestSM.y = -1.0*Xtest*sin(LT); PtestSM.z = 0.0;
        Lgm_Convert_Coords( &PtestSM, &Ptest, SM_TO_GSM, LstarInfo_test->mInfo->c );

        //Trace to minimum-B
        TFlag = Lgm_Trace( &Ptest, &v1, &v2, &v3, 120.0, 0.01, TRACE_TOL, LstarInfo_test->mInfo );
        if ( TFlag == LGM_CLOSED ) {
            if ( Lgm_Setup_AlphaOfK( &DT_UTC, &v3, LstarInfo_test->mInfo ) > 0 ) {
                Alpha = Lgm_AlphaOfK( Kin, LstarInfo_test->mInfo );
                Lgm_TearDown_AlphaOfK(LstarInfo_test->mInfo);
            }
            else {
                Alpha = LGM_FILL_VALUE;
            }
            LstarInfo_test->PitchAngle = Alpha;
            sa = sin( LstarInfo_test->PitchAngle*RadPerDeg ); sa2 = sa*sa;

            //If test point is closed FL then check for undefined L*.
            //Get L*
            LstarInfo_test->mInfo->Bm = LstarInfo_test->mInfo->Bmin/sa2;
            LS_Flag = Lstar( &v3, LstarInfo_test);
            if ( (LS_Flag > 0) || (LstarInfo_test->LS != LGM_FILL_VALUE) ){
                //Drift shell defined
                Pinner = Ptest;
                LCDSv3 = v3;
                LCDS = LstarInfo_test->LS;
                *K = (LstarInfo_test->I[0])*sqrt(LstarInfo_test->mInfo->Bm*nTtoG);

                //Determine the type of the orbit
                LstarInfo_test->DriftOrbitType = LGM_DRIFT_ORBIT_CLOSED;
                for ( k=0; k<LstarInfo_test->nMinMax; ++k ) {
                    if ( LstarInfo_test->nMinima[k] > 1 ) LstarInfo_test->DriftOrbitType = LGM_DRIFT_ORBIT_CLOSED_SHABANSKY;
                    if ( LstarInfo_test->nMinima[k] < 1 ) {printf("Less than one minimum defined on field line: Exit due to impossibility."); exit(-1); }
                }

                if (LstarInfo->VerbosityLevel > 0) printf("Current LCDS, K, alpha, PtestSM is %g, %g, %g, (%g, %g, %g)\n", 
                        LCDS, *K, Alpha, PtestSM.x, PtestSM.y, PtestSM.z);
                LstarInfo->mInfo->Bm = LstarInfo_test->mInfo->Bm;
                LstarInfo->DriftOrbitType = LstarInfo_test->DriftOrbitType;
            } else {
                Pouter = Ptest;
                Lgm_Convert_Coords( &Ptest, &PtestSM, GSM_TO_SM, LstarInfo_test->mInfo->c );
                if (LstarInfo->VerbosityLevel > 0) {
                    printf("Failed DS trace for alpha = %g (K=%g) at test point (%g %g %g GSM; %g %g %g SM; TFlag = %d), moving outer bracket to current location\n",
                        LstarInfo_test->PitchAngle, *K, Ptest.x, Ptest.y, Ptest.z, PtestSM.x, PtestSM.y, PtestSM.z, TFlag);
                    }
            }
        } else {
            //FL open -> drift shell open
            Pouter = Ptest;
            if (LstarInfo->VerbosityLevel > 0) {
                printf("Failed FL trace at test point (%g %g %g GSM; %g %g %g SM; TFlag = %d), moving outer bracket to current location\n",
                    Ptest.x, Ptest.y, Ptest.z, PtestSM.x, PtestSM.y, PtestSM.z, TFlag);
                }
        }

    nn++;
    }
    if (LstarInfo->VerbosityLevel > 0) printf("Final LCDS, K is %g, %g. Eq. radius = %g \n", LCDS, *K, Lgm_Magnitude( &LCDSv3 ));
    LstarInfo->LS = LCDS;

    //free structures
    FreeLstarInfo( LstarInfo_brac1 );
    FreeLstarInfo( LstarInfo_brac2 );
    FreeLstarInfo( LstarInfo_test );
    return(0);


}



