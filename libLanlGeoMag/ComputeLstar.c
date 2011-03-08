#include <stdio.h>
#include <math.h>
#include "Lgm/Lgm_QuadPack.h"
#include "Lgm/Lgm_MagModelInfo.h"
#include "Lgm/Lgm_LstarInfo.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <string.h>

#define DeltaMLT 1.0


void PredictMlat1( double *MirrorMLT, double *MirrorMlat, int k, double MLT, double *pred_mlat, double *pred_delta_mlat, double *delta );
void PredictMlat2( double *MirrorMLT, double *MirrorMlat, int k, double MLT, double *pred_mlat, double *pred_delta_mlat, double *delta, Lgm_LstarInfo *LstarInfo );


/*
 * There are many tolerances involved in an Lstar calculation.  Here we try and
 * set them qualitatively given a single "quality" value.
 */
void SetLstarTolerances( int Quality, Lgm_LstarInfo *s ) {

    if ( ( Quality < 0 ) || ( Quality > 8 ) ) {
        printf("%sSetLstarTolerances: Quality value (of %d) not in range [0, 8]. Setting to 5.%s\n", s->PreStr, Quality, s->PostStr );
        Quality = 5;
    }

    // These tend to be critical to keep things working smoothly.
    s->mInfo->Lgm_FindBmRadius_Tol       = 1e-10;
    s->mInfo->Lgm_TraceToMirrorPoint_Tol = 1e-10;


    switch ( Quality ) {

        case 8: // Highest Quality -- (although 7 may be better?)

            s->mInfo->Lgm_MagFlux_Integrator_epsabs = 0.0;
            s->mInfo->Lgm_MagFlux_Integrator_epsrel = 1e-6;

            s->mInfo->Lgm_LambdaIntegral_Integrator_epsabs = 0.0;
            s->mInfo->Lgm_LambdaIntegral_Integrator_epsrel = 1e-6;

            s->mInfo->Lgm_I_Integrator        = DQAGS;
            s->mInfo->Lgm_I_Integrator_epsrel = 0.0;
            s->mInfo->Lgm_I_Integrator_epsabs = 1e-8;

            s->mInfo->Lgm_FindShellLine_I_Tol = 1e-8;

            break;

        case 7:

            s->mInfo->Lgm_MagFlux_Integrator_epsabs = 0.0;
            s->mInfo->Lgm_MagFlux_Integrator_epsrel = 1e-6;

            s->mInfo->Lgm_LambdaIntegral_Integrator_epsabs = 0.0;
            s->mInfo->Lgm_LambdaIntegral_Integrator_epsrel = 1e-6;

            s->mInfo->Lgm_I_Integrator        = DQAGS;
            s->mInfo->Lgm_I_Integrator_epsrel = 0.0;
            s->mInfo->Lgm_I_Integrator_epsabs = 1e-7;

            s->mInfo->Lgm_FindShellLine_I_Tol = 1e-7;

            break;

        case 6:

            s->mInfo->Lgm_MagFlux_Integrator_epsabs = 0.0;
            s->mInfo->Lgm_MagFlux_Integrator_epsrel = 1e-5;

            s->mInfo->Lgm_LambdaIntegral_Integrator_epsabs = 0.0;
            s->mInfo->Lgm_LambdaIntegral_Integrator_epsrel = 1e-5;

            s->mInfo->Lgm_I_Integrator        = DQAGS;
            s->mInfo->Lgm_I_Integrator_epsrel = 0.0;
            s->mInfo->Lgm_I_Integrator_epsabs = 1e-6;

            s->mInfo->Lgm_FindShellLine_I_Tol = 1e-6;

            break;

        case 5:

            s->mInfo->Lgm_MagFlux_Integrator_epsabs = 0.0;
            s->mInfo->Lgm_MagFlux_Integrator_epsrel = 1e-4;

            s->mInfo->Lgm_LambdaIntegral_Integrator_epsabs = 0.0;
            s->mInfo->Lgm_LambdaIntegral_Integrator_epsrel = 1e-4;

            s->mInfo->Lgm_I_Integrator        = DQAGS;
            s->mInfo->Lgm_I_Integrator_epsrel = 0.0;
            s->mInfo->Lgm_I_Integrator_epsabs = 1e-5;

            s->mInfo->Lgm_FindShellLine_I_Tol = 1e-5;

            break;

        case 4:

            s->mInfo->Lgm_MagFlux_Integrator_epsabs = 0.0;
            s->mInfo->Lgm_MagFlux_Integrator_epsrel = 1e-3;

            s->mInfo->Lgm_LambdaIntegral_Integrator_epsabs = 0.0;
            s->mInfo->Lgm_LambdaIntegral_Integrator_epsrel = 1e-3;

            s->mInfo->Lgm_I_Integrator        = DQAGS;
            s->mInfo->Lgm_I_Integrator_epsrel = 0.0;
            s->mInfo->Lgm_I_Integrator_epsabs = 1e-4;

            s->mInfo->Lgm_FindShellLine_I_Tol = 1e-4;

            break;

        case 3:

            s->mInfo->Lgm_MagFlux_Integrator_epsabs = 0.0;
            s->mInfo->Lgm_MagFlux_Integrator_epsrel = 1e-3;

            s->mInfo->Lgm_LambdaIntegral_Integrator_epsabs = 0.0;
            s->mInfo->Lgm_LambdaIntegral_Integrator_epsrel = 1e-3;

            s->mInfo->Lgm_I_Integrator        = DQAGS;
            s->mInfo->Lgm_I_Integrator_epsrel = 0.0;
            s->mInfo->Lgm_I_Integrator_epsabs = 1e-3;

            s->mInfo->Lgm_FindShellLine_I_Tol = 1e-3;

            break;

        case 2:

            s->mInfo->Lgm_MagFlux_Integrator_epsabs = 0.0;
            s->mInfo->Lgm_MagFlux_Integrator_epsrel = 1e-3;

            s->mInfo->Lgm_LambdaIntegral_Integrator_epsabs = 0.0;
            s->mInfo->Lgm_LambdaIntegral_Integrator_epsrel = 1e-3;

            s->mInfo->Lgm_I_Integrator        = DQAGS;
            s->mInfo->Lgm_I_Integrator_epsrel = 0.0;
            s->mInfo->Lgm_I_Integrator_epsabs = 1e-2;

            s->mInfo->Lgm_FindShellLine_I_Tol = 1e-2;

            break;

        case 1:

            s->mInfo->Lgm_MagFlux_Integrator_epsabs = 0.0;
            s->mInfo->Lgm_MagFlux_Integrator_epsrel = 1e-3;

            s->mInfo->Lgm_LambdaIntegral_Integrator_epsabs = 0.0;
            s->mInfo->Lgm_LambdaIntegral_Integrator_epsrel = 1e-3;

            s->mInfo->Lgm_I_Integrator        = DQAGS;
            s->mInfo->Lgm_I_Integrator_epsrel = 0.0;
            s->mInfo->Lgm_I_Integrator_epsabs = 1e-1;

            s->mInfo->Lgm_FindShellLine_I_Tol = 1e-1;

            break;

        case 0:

            s->mInfo->Lgm_MagFlux_Integrator_epsabs = 0.0;
            s->mInfo->Lgm_MagFlux_Integrator_epsrel = 1e-3;

            s->mInfo->Lgm_LambdaIntegral_Integrator_epsabs = 0.0;
            s->mInfo->Lgm_LambdaIntegral_Integrator_epsrel = 1e-3;

            s->mInfo->Lgm_I_Integrator        = DQK21; // Note - changed to simpler integrator.
            s->mInfo->Lgm_I_Integrator_epsrel = 0.0;
            s->mInfo->Lgm_I_Integrator_epsabs = 1e-1;

            s->mInfo->Lgm_FindShellLine_I_Tol = 1e-1;

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

    LstarInfo->PreStr[0]  = '\0';
    LstarInfo->PostStr[0] = '\0';

    return LstarInfo;
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


    t->mInfo->Lgm_MagStep_FirstTimeThrough = TRUE;
    t->mInfo->Lgm_MagStep_eps_old = -1.0;


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
    int		done2, Count, FoundShellLine;
    double	B, dSa, dSb, smax, SS, L, Hmax;
    double	I=-999.9, Ifound, M, MLT0, MLT, mlat, r;
    double	Phi, Phi1, Phi2, sl, cl, MirrorMLT[500], MirrorMlat[500], pred_mlat, pred_delta_mlat=0.0, mlat0, mlat1, delta;
    double	MirrorMLT_Old[500], MirrorMlat_Old[500];
    char    *PreStr, *PostStr;
    Lgm_MagModelInfo    *mInfo2;
//    FILE	*fp;

    PreStr = LstarInfo->PreStr;
    PostStr = LstarInfo->PostStr;



    if (LstarInfo->VerbosityLevel > 0) {
        printf("\n\n\t\t%s        Computing L* for, Date: %ld, UT: %g%s\n", PreStr, LstarInfo->mInfo->c->UTC.Date, LstarInfo->mInfo->c->UTC.Time, PostStr );
        printf(    "\t\t%s=========================================================================%s\n", PreStr, PostStr );
        printf(    "\t\t%sDate (yyyymmdd):                         %ld%s\n", PreStr, LstarInfo->mInfo->c->UTC.Date, PostStr );
        printf(    "\t\t%sUTC (hours):                             %g%s\n", PreStr, LstarInfo->mInfo->c->UTC.Time, PostStr );
        printf(    "\t\t%sPitch Angle (deg.):                      %g%s\n", PreStr, LstarInfo->PitchAngle, PostStr );
        printf(    "\t\t%sMirror Mag. Field Strength, Bm (nT):     %g%s\n", PreStr, LstarInfo->mInfo->Bm, PostStr );
    }
    LstarInfo->nPnts = 0;


    if ((LstarInfo->PitchAngle < 0.0)||(LstarInfo->PitchAngle>90.0)) return(-1);



    /*
     *  Do Initial field Line to get Bm and I
     */
    u = *vin;
	if (LstarInfo->VerbosityLevel > 0) {
        printf("\n\t\t%sInitial Position, U_gsm (Re):            < %g, %g, %g >%s\n", PreStr, u.x, u.y, u.z, PostStr);
	    LstarInfo->mInfo->Bfield( &u, &Bvec, LstarInfo->mInfo );
	    B = Lgm_Magnitude( &Bvec );
        printf("\t\t%sMag. Field Strength, B at U_gsm (nT):    %g%s\n", PreStr, B, PostStr);
    }
    if ( Lgm_Trace( &u, &v1, &v2, &v3, LstarInfo->mInfo->Lgm_LossConeHeight, 1e-7, 1e-7, LstarInfo->mInfo ) == 1 ) {



        /*
         * Trace from Bmin point up to northern mirror point and down to
         * southern mirror point. dSa and dSb are the distances along the
         * FL from the starting points point. So dSa is from the Bmin value.
         * And dSb is from the Psouth value.
         */
        if ( Lgm_TraceToMirrorPoint( &(LstarInfo->mInfo->Pmin), &(LstarInfo->mInfo->Pm_South), &dSa, LstarInfo->mInfo->Bm, -1.0, LstarInfo->mInfo->Lgm_TraceToMirrorPoint_Tol, LstarInfo->mInfo ) >= 0 ) {

	        if (LstarInfo->VerbosityLevel > 0) {
                printf("\n\t\t%sMin-B  Point Location, Pmin (Re):      < %g, %g, %g >%s\n", PreStr, LstarInfo->mInfo->Pmin.x, LstarInfo->mInfo->Pmin.y, LstarInfo->mInfo->Pmin.z, PostStr);
                printf("\n\t\t%sMirror Point Location, Pm_South (Re):      < %g, %g, %g >%s\n", PreStr, LstarInfo->mInfo->Pm_South.x, LstarInfo->mInfo->Pm_South.y, LstarInfo->mInfo->Pm_South.z, PostStr);
	            LstarInfo->mInfo->Bfield( &LstarInfo->mInfo->Pm_South, &Bvec, LstarInfo->mInfo );
	            B = Lgm_Magnitude( &Bvec );
                printf("\t\t%sMag. Field Strength, Bm at Pm_South (nT):  %g     (LstarInfo->mInfo->Bm = %g)%s\n", PreStr, B, LstarInfo->mInfo->Bm, PostStr);
            }

            if ( Lgm_TraceToMirrorPoint( &(LstarInfo->mInfo->Pm_South), &(LstarInfo->mInfo->Pm_North), &dSb, LstarInfo->mInfo->Bm,  1.0, LstarInfo->mInfo->Lgm_TraceToMirrorPoint_Tol, LstarInfo->mInfo ) >= 0 ) {

                if (LstarInfo->VerbosityLevel > 0) {
                    printf("\n\t\t%sMin-B  Point Location, Pmin (Re):      < %g, %g, %g >%s\n", PreStr, LstarInfo->mInfo->Pmin.x, LstarInfo->mInfo->Pmin.y, LstarInfo->mInfo->Pmin.z, PostStr);
                    printf("\n\t\t%sMirror Point Location, Pm_North (Re):      < %g, %g, %g >%s\n", PreStr, LstarInfo->mInfo->Pm_North.x, LstarInfo->mInfo->Pm_North.y, LstarInfo->mInfo->Pm_North.z, PostStr);
                    LstarInfo->mInfo->Bfield( &LstarInfo->mInfo->Pm_North, &Bvec, LstarInfo->mInfo );
                    B = Lgm_Magnitude( &Bvec );
                    printf("\t\t%sMag. Field Strength, Bm at Pm_North (nT):  %g%s\n", PreStr, B, PostStr);
                }


                /*
                 *  Set the limits of integration. Define s=0 at the sourthern mirror point. Then, sm_North will just be dSb
                 */
                //LstarInfo->mInfo->Sm_South = LstarInfo->mInfo->smin - dSa;
                //LstarInfo->mInfo->Sm_North = LstarInfo->mInfo->Sm_South + dSb;
                SS = dSb;
                LstarInfo->mInfo->Hmax = SS/200.0;
                r  = Lgm_Magnitude( &LstarInfo->mInfo->Pm_North );
                LstarInfo->mInfo->Sm_South = 0.0;
                LstarInfo->mInfo->Sm_North = SS;

                if ( LstarInfo->mInfo->UseInterpRoutines ) {

                    Lgm_TraceLine2( &(LstarInfo->mInfo->Pm_South), &LstarInfo->mInfo->Pm_North, (r-1.0)*Re, 0.5*SS-LstarInfo->mInfo->Hmax, 1.0, 1e-7, FALSE, LstarInfo->mInfo );
                    ReplaceFirstPoint( 0.0, LstarInfo->mInfo->Bm, &LstarInfo->mInfo->Pm_South, LstarInfo->mInfo );
                    AddNewPoint( SS,  LstarInfo->mInfo->Bm, &LstarInfo->mInfo->Pm_North, LstarInfo->mInfo );
                    InitSpline( LstarInfo->mInfo );

                    /*
                     *  Do interped I integral.
                     */
                    I = Iinv_interped( LstarInfo->mInfo  );
                    if (LstarInfo->VerbosityLevel > 0) printf("\t\t  %sIntegral Invariant, I (interped):      %g%s\n",  PreStr, I, PostStr );

                    FreeSpline( LstarInfo->mInfo );


                } else {

                    /*
                     *  Do full blown I integral. (Integrand is evaluated by tracing to required s-values.)
                     */
                    I = Iinv( LstarInfo->mInfo  );
                    if (LstarInfo->VerbosityLevel > 0) printf("\t\t  %sIntegral Invariant, I (full integral): %g%s\n",  PreStr, I, PostStr );

                }



                if (LstarInfo->VerbosityLevel > 0) {
                    // sort this out. FIX User should decide what M they want to use.
                    M = LstarInfo->mInfo->c->M_cd_McIllwain;
                    M = LstarInfo->mInfo->c->M_cd;
                    printf("\t\t  %sLgm_n_I_integrand_Calls:               %d%s\n\n", PreStr, LstarInfo->mInfo->Lgm_n_I_integrand_Calls, PostStr );
                    printf("\t\t%sCurrent Dipole Moment, M_cd:             %g%s\n", PreStr, LstarInfo->mInfo->c->M_cd, PostStr);
                    printf("\t\t%sReference Dipole Moment, M_cd_McIllwain: %g%s\n", PreStr, LstarInfo->mInfo->c->M_cd_McIllwain, PostStr);
                    printf("\t\t%sDipole Moment Used, Mused:               %g%s\n", PreStr, M, PostStr);
                    printf("\t\t%sMcIlwain L (Hilton):                     %g%s\n", PreStr, L = LFromIBmM_Hilton( I, LstarInfo->mInfo->Bm, M ), PostStr );
                    printf("\t\t%sMcIlwain L (McIlwain):                   %g%s\n", PreStr, L = LFromIBmM_McIlwain( I, LstarInfo->mInfo->Bm, M ), PostStr );
                }

            } else {

	            printf("\t\t%sMirror point below %g km in Northern Hemisphere%s\n", PreStr, LstarInfo->mInfo->Lgm_LossConeHeight, PostStr );
	            return(-1);

	        }

        } else {

	        printf("\t\t%sMirror point below %g km in Southern Hemisphere%s\n", PreStr, LstarInfo->mInfo->Lgm_LossConeHeight, PostStr );
	        return(-2);

        }

    }






    /*
     *  Construct drift shell. Get a good initial estimate for mlat
     */
    Lgm_Convert_Coords( &LstarInfo->mInfo->Pm_North, &u, GSM_TO_SM, LstarInfo->mInfo->c );
    mlat = DegPerRad*asin(u.z/Lgm_Magnitude( &u ));
    MLT0 = atan2( u.y, u.x )*DegPerRad/15.0 + 12.0;
    if (LstarInfo->VerbosityLevel > 2) printf("\t\t%smlat = %g%s\n", PreStr, mlat, PostStr );
    LstarInfo->nPnts = 0;
    pred_delta_mlat  = 0.0;
    for (k=0; k<200; ++k){
	    MirrorMLT[k]  = 0.0;
	    MirrorMlat[k] = 0.0;
    }



    nLines = (int)(24.0/DeltaMLT + 0.5);
    for ( koffset=nk=k=0, MLT=MLT0; MLT<(MLT0+24.0-1e-10); MLT += DeltaMLT){


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
if ( k < 3 )       delta = 3.0;
else if ( k < 6 )  delta = 2.0;
else if ( k < 9 )  delta = 2.0;
else if ( k < 12 ) delta = 2.0;
else if ( k < 15 ) delta = 2.0;
else if ( k < 18 ) delta = 1.0;
else if ( k < 21 ) delta = 0.5;
else               delta = 0.5;
delta = 5.0;
//fprintf(fp_predict, "%g %g\n", MLT, pred_mlat );



	    /*
	     *  Now set a range to search. It should be pretty narrow to start.
	     *  If its too small the search may fail. In which case we need to
	     *  enlarge the search.
    	 */
	    done2 = FALSE; FoundShellLine = FALSE; Count = 0;
	    while ( !done2 ) {

	        if ( Count == 0 ) {

                /*
                 *  First time through -- lets use the predicted range.
                 */
                //mlat0 = ((pred_mlat-delta) <  0.0) ?  0.0 : (pred_mlat-delta);
                mlat0 = ((pred_mlat-delta) <  0.0) ?  0.0 : (pred_mlat-delta);
                mlat1 = ((pred_mlat+delta) > 90.0) ? 90.0 : (pred_mlat+delta);

    	    } else if ( Count == 1 ) {

                /*
                 * Perhaps our bracket was no good -- too narrow. Enlarge it.
                 * Try double what we had.
                 */
                delta *= 2;
                mlat0 = ((mlat-delta) >  0.0) ?  0.0 : (mlat-delta);
                mlat1 = ((mlat+delta) > 90.0) ? 90.0 : (mlat+delta);

            } else if ( Count == 2 ) {

	            /*
	             * Hmmm... This ones stuborn! Lets be more conservative now.
	            * Try +1/-5 around the returned mlat
	            */

                if ( FoundShellLine == -3 ) {
                    mlat0 = ((mlat-5.0) >  0.0) ?  0.0 : (mlat-5.0);
                    mlat1 = ((mlat+1.0) > 90.0) ? 90.0 : (mlat+1.0);
                } else {
                    mlat0 = ((mlat-1.0) >  0.0) ?  0.0 : (mlat-1.0);
                    mlat1 = ((mlat+5.0) > 90.0) ? 90.0 : (mlat+5.0);
                }

    	    } else {

                /*
                 * OK, there may not be a good line -- i.e. drift shell may not be closed.
                 * To make sure, lets try 0->90 deg.
                 */
                mlat0 = 0.0;
mlat0 = -30.0;
                mlat1 = 90.0;

            }


            if (LstarInfo->VerbosityLevel > 1) printf("\n\t\t%s-------------------  Line %02d of %02d   MLT: %g  (mlat0, mlat1 = %g %g) ---------------------%s\n", PreStr, k, nLines, MLT, mlat0, mlat1, PostStr );
            FoundShellLine = FindShellLine( I, &Ifound, LstarInfo->mInfo->Bm, MLT, &mlat, &r, mlat0, mlat1, LstarInfo );

            if ( FoundShellLine > 0 ) {
                done2 = TRUE;
            } else if ( Count > 2 ) {
                done2 =  TRUE;
                printf(" \t%sNo valid I - Drift Shell not closed: L* = undefined  (FoundShellLine = %d)%s\n", PreStr, FoundShellLine, PostStr); fflush(stdout);
                FoundShellLine = 0;
                return(-3);
            } else {
                ++Count;
            }
	    }

        /*
         * Save individual I values
         */
        LstarInfo->I[k] = Ifound;


        MirrorMLT[k]  = MLT;
        MirrorMlat[k] = mlat;

        if (LstarInfo->VerbosityLevel > 2) printf("\t\t%sActual mlat         = %g  \t Count = %d%s", PreStr, mlat, Count, PostStr ); fflush(stdout);


        /*
         *  convert mirror point to GSM.
         */
        Phi = 15.0*(MLT-12.0)*RadPerDeg;
        cl = cos( mlat * RadPerDeg ); sl = sin( mlat * RadPerDeg );
        u.x = r*cl*cos(Phi); u.y = r*cl*sin(Phi); u.z = r*sl;
        Lgm_Convert_Coords( &u, &v, SM_TO_GSM, LstarInfo->mInfo->c );

        /*
         * Save GSM cartesian as well...
         */
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
        if ( !Lgm_TraceToSphericalEarth( &v, &w, LstarInfo->mInfo->Lgm_LossConeHeight, 1.0, 1e-7, LstarInfo->mInfo ) ){ return(-4); }
        LstarInfo->Spherical_Footprint_Pn[k] = w;





        /*
         *  convert footpoint back to SM
         */
        Lgm_Convert_Coords( &w, &v, GSM_TO_SM, LstarInfo->mInfo->c );
        Phi = atan2( v.y, v.x );
        LstarInfo->MLT[k]  = Phi*DegPerRad/15.0 + 12.0;
        LstarInfo->mlat[k] = asin( v.z/Lgm_Magnitude(&v) )*DegPerRad;
        if (LstarInfo->VerbosityLevel > 2) printf(" \t%sMLT_foot, mlat_foot = %g %g%s\n\n", PreStr, LstarInfo->MLT[k], LstarInfo->mlat[k], PostStr); fflush(stdout);




        /*
         * If SaveShellLines is set true, then retrace the FL and save the
         * whole FL.  Note that since we appear to only have the north
         * footpoint, we start there and trace to south. So lets pack them in
         * the the saved arrays backwards so that they go from south to north.
         */
        if ( LstarInfo->FindShellPmin || LstarInfo->ComputeVgc ) {
            Lgm_TraceToMinBSurf( &LstarInfo->Spherical_Footprint_Pn[k], &v2, 0.1, 1e-8, LstarInfo->mInfo );
            LstarInfo->mInfo->Bfield( &v2, &LstarInfo->Bmin[k], LstarInfo->mInfo );
            LstarInfo->Pmin[k] = v2;
        }

        /*
         * Compute the gradient of I, Sb and Vgc
         */
        if ( LstarInfo->ComputeVgc ) {
            // Lgm_Grad_I() and other rotuines may modify mInfo in undesirable ways, so give it a copy.
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
            K *= 1e3*EE;                        // Joules
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
            Lgm_TraceToSphericalEarth( &LstarInfo->Mirror_Ps[k], &uu, LstarInfo->mInfo->Lgm_LossConeHeight, -1.0, 1e-7, LstarInfo->mInfo );
                LstarInfo->Mirror_Ss[k] += LstarInfo->mInfo->Trace_s;
                LstarInfo->Mirror_Sn[k] += LstarInfo->mInfo->Trace_s;
            //}





        }










        ++k;

    } // end MLT for loop
    LstarInfo->nPnts = k;


    /*
     *  Save drift shell -- it will help us predict the next one.
     */
    for (i=0; i<k; ++i){
	    MirrorMLT_Old[i]  = MirrorMLT[i];
	    MirrorMlat_Old[i] = MirrorMlat[i];
    }

    /*
     *  To get Lstar all we need to do now is one final integral.
     *  See roederer's eq 3.6 It should be trivial at this point!
     */


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
    for (i=0; i<k; ++i){ LstarInfo->xa[j] = LstarInfo->MLT[i]-24.0; LstarInfo->ya[j] = LstarInfo->mlat[i]; ++j; }
    for (i=0; i<k; ++i){ LstarInfo->xa[j] = LstarInfo->MLT[i]     ; LstarInfo->ya[j] = LstarInfo->mlat[i]; ++j; }
    for (i=0; i<k; ++i){ LstarInfo->xa[j] = LstarInfo->MLT[i]+24.0; LstarInfo->ya[j] = LstarInfo->mlat[i]; ++j; }
    spline( LstarInfo->xa, LstarInfo->ya, LstarInfo->nSplnPnts, 0.0, 0.0, LstarInfo->y2 );
*/

    j = 0;
    for (i=0; i<k; ++i){ LstarInfo->xa[j] = LstarInfo->MLT[i]-24.0; LstarInfo->ya[j] = LstarInfo->mlat[i]; ++j; }
    for (i=0; i<k; ++i){ LstarInfo->xa[j] = LstarInfo->MLT[i]     ; LstarInfo->ya[j] = LstarInfo->mlat[i]; ++j; }
    for (i=0; i<k; ++i){ LstarInfo->xa[j] = LstarInfo->MLT[i]+24.0; LstarInfo->ya[j] = LstarInfo->mlat[i]; ++j; }
    LstarInfo->nSplnPnts = j;

    /*
     *  Use GSL periodic spline (can also try others if needed....)
     */
//    gsl_set_error_handler_off(); // Turn off gsl default error handler
    LstarInfo->acc = gsl_interp_accel_alloc( );
    LstarInfo->pspline = gsl_interp_alloc( gsl_interp_cspline_periodic, LstarInfo->nSplnPnts );
    //printf ("spline uses '%s' interpolation.\n", gsl_interp_name( LstarInfo->pspline ));
    gsl_interp_init( LstarInfo->pspline, LstarInfo->xa, LstarInfo->ya, LstarInfo->nSplnPnts );








/*
For debugging...
    fp = fopen("ShellFootTrace.dat", "w");
    for (MLT=-24.0; MLT<=48.0; MLT += 0.01){
	    //splint( LstarInfo->xa, LstarInfo->ya, LstarInfo->y2, LstarInfo->nSplnPnts, MLT, &mlat);
	    mlat = gsl_interp_eval( LstarInfo->pspline, LstarInfo->xa, LstarInfo->ya, MLT, LstarInfo->acc );
	    fprintf(fp, "%g %g\n", MLT, mlat);
    }
    fclose(fp);
*/


    Phi1 = MagFlux( LstarInfo );
    LstarInfo->LS_dip_approx = -2.0*M_PI*LstarInfo->mInfo->c->M_cd /Phi1;
    Phi2 = MagFlux2( LstarInfo );
    LstarInfo->LS = -2.0*M_PI*LstarInfo->mInfo->c->M_cd /Phi2;
    LstarInfo->LS_McIlwain_M = -2.0*M_PI*LstarInfo->mInfo->c->M_cd_McIllwain /Phi2;

    if (LstarInfo->VerbosityLevel > 0) {
        printf("\n\t\t%sL*, Dipole Approximation.\n%s", PreStr, PostStr );
        printf("\t\t%s  Magnetic Flux:                         %.15lf%s\n", PreStr, Phi1, PostStr );
        printf("\t\t%s  L*:                                    %.15lf%s\n", PreStr, LstarInfo->LS_dip_approx, PostStr );
        printf("\t\t%s  L* (Using McIllwain M):                %.15lf%s\n", PreStr, -2.0*M_PI*LstarInfo->mInfo->c->M_cd_McIllwain /Phi1, PostStr );
        printf("\n\t\t%sL*, Full Field.%s\n", PreStr, PostStr );
        printf("\t\t%s  Magnetic Flux:                         %.15lf%s\n", PreStr, Phi2, PostStr );
        printf("\t\t%s  L*:                                    %.15lf%s\n", PreStr, LstarInfo->LS, PostStr );
        printf("\t\t%s  L* (Using McIllwain M):                %.15lf%s\n\n\n\n", PreStr, -2.0*M_PI*LstarInfo->mInfo->c->M_cd_McIllwain /Phi2, PostStr );
    }

    gsl_interp_free( LstarInfo->pspline );
    gsl_interp_accel_free( LstarInfo->acc );


    return(0);

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
    int         key, neval, ier, limit, lenw, last, iwork[501];
    double      work[2001];
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
    epsabs = LstarInfo->mInfo->Lgm_MagFlux_Integrator_epsrel;



    limit = 500; lenw = 4*limit; key = 6;
/*
    iwork  = (int *) calloc( limit+1, sizeof(int) );
    work   = (double *) calloc( lenw+1, sizeof(double) );
*/
    dqags(MagFluxIntegrand, qpInfo, a, b, epsabs, epsrel, &result, &abserr, &neval, &ier, limit, lenw, &last, iwork, work);
/*
    free( iwork );
    free( work );
*/


    r = 1.0 + LstarInfo->mInfo->Lgm_LossConeHeight/WGS84_A;

    return( -result*LstarInfo->mInfo->c->M_cd/r );

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
    int         key, neval, ier, limit, lenw, last, iwork[501];
    double      work[2001], MLT, mlat;
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
    epsabs = LstarInfo->mInfo->Lgm_LambdaIntegral_Integrator_epsabs;
    epsrel = LstarInfo->mInfo->Lgm_LambdaIntegral_Integrator_epsrel;

    limit = 500; lenw = 4*limit; key = 6;
/*
    iwork  = (int *) calloc( limit+1, sizeof(int) );
    work   = (double *) calloc( lenw+1, sizeof(double) );
*/
    dqags(LambdaIntegrand, qpInfo, a, b, epsabs, epsrel, &result, &abserr, &neval, &ier, limit, lenw, &last, iwork, work);
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
    int         key, neval, ier, limit, lenw, last, iwork[501];
    double      work[2001];
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
    epsabs = 0.0;
    epsrel = 1e-5;

    limit = 500; lenw = 4*limit; key = 6;
/*
    iwork  = (int *) calloc( limit+1, sizeof(int) );
    work   = (double *) calloc( lenw+1, sizeof(double) );
*/
    dqags(MagFluxIntegrand2, qpInfo, a, b, epsabs, epsrel, &result, &abserr, &neval, &ier, limit, lenw, &last, iwork, work);
/*
    free( iwork );
    free( work );
*/


    r = 1.0 + LstarInfo->mInfo->Lgm_LossConeHeight/WGS84_A;

    return( result/r );

}

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


	    Lgm_RatFunInt( MirrorMLT+koffset-1, MirrorMlat+koffset-1, nk, MLT, &pred_mlat1, &pred_delta_mlat1 );
        Lgm_PolFunInt( MirrorMLT+koffset-1, MirrorMlat+koffset-1, nk, MLT, &pred_mlat2, &pred_delta_mlat2 );

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


        gsl_set_error_handler_off(); // Turn off gsl default error handler
	    gsl_interp_accel *acc = gsl_interp_accel_alloc( );
	    gsl_interp *pspline = gsl_interp_alloc( gsl_interp_akima_periodic, LstarInfo->m );
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
 *    $Id$
 */
