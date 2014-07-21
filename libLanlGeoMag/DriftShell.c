#include "Lgm/Lgm_MagModelInfo.h"
#include "Lgm/Lgm_LstarInfo.h"
#include "Lgm/qsort.h"
#include <gsl/gsl_multifit.h>
#define MAX_ITS 100

typedef struct BracketType {
    double  a, b, c;
    double  Da, Db, Dc;
    double  mlat_min, Dmin;
    int     FoundZeroBracket;
} BracketType;

typedef struct FitVals {
    double  key;
    double  x, y, dy;
} FitVals;

static void elt_qsort( struct FitVals *arr, unsigned n ) {
    #define elt_lt(a,b) ((a)->key < (b)->key)
    QSORT( struct FitVals, arr, n, elt_lt );
}




int FitLineAndFindZero( double *x, double *y, double *dy, int n, double *res );
int FitQuadAndFindZero( double *x, double *y, double *dy, int n, double *res );
int FitQuadAndFindZero2( double *x, double *y, double *dy, int n, int nmax, double *res );
int BracketZero( double I0, double *Ifound, double Bm, double MLT, double *mlat, double *rad, double mlat0, double mlat1, double mlat2, BracketType *Bracket, Lgm_LstarInfo *LstarInfo );

/*
 *   FindShellLine
 *   -------------
 */
//! Find the field line that has the given I and Bm values at the specified MLT.
/**
 *            \param[in]        I0          the value of the integral invariant, I,  for the desired drift shell.
 *            \param[out]       Ifound      the value of the integral invariant, I,  that was actually found for this shell line.
 *            \param[in]        Bm          the value of the mirror field strength for the desired drift shell.
 *            \param[in]        MLT         magnetic local time meridian in which to search for shell field line.
 *            \param[in,out]    mlat        an initial guess for the magnetic latitude of the mirror pointat which.
 *            \param[in,out]    rad         the geocentric radial distance of the mirror point (at which B=Bm).
 *            \param[in]        mlat0       lower end of mlat bracket to search over.
 *            \param[in]        mlat1       upper end of mlat bracket to search over.
 *            \param[in,out]    LstarInfo   A properly initialized Lgm_LstarInfo structure.
 *
 *            \return 1 - normal exit. Field line found with acceptable accuracy.
 *            \return 0 - no field line could be found
 *
 */
int FindShellLine(  double I0, double *Ifound, double Bm, double MLT, double *mlat, double *rad, double mlat0, double mlat1, double mlat2, int *Iterations, Lgm_LstarInfo *LstarInfo) {
    Lgm_Vector      u, w, Pm_North, Pmirror, v1, v2, v3;
    double          F, F0, F1, rat, a, b, c, d, d0, d1, Da, Db, Dc, De, I, r, Phi, cl, sl;
    double          SS, Sn, Ss, mlat_min=0.0, Dmin=9e99, e, D, D0, D1, D2, Sign, Dbest, mlatbest, res;
    int             done, FoundValidI, FirstHalf, nIts, FoundZeroBracket;
    int             i, Flag, nbFits;
    BracketType     Bracket;
    double          BRACKET_EPS;

    *Iterations = 0;

    /*
     *  Set the bracket for the mlat range. We got mlat0 and mlat2 from the
     *  caller. The hope is that the field line with the I we are looking for
     *  is between these two latitudes.  (If it isnt we will eventually have to
     *  bail and the user has to try a different range.) We also got mlat1
     *  which is the user's best guess at where it ought to be.
     *
     *  Define D to be I-I0.
     *  We need to find the latitude that gives D=0 (or 0 to within tolerance).
     *
     *  If there is a solution (zero-crossing) in the given range
     *  (mlat0->mlat2), we can get the following cases:
     *
     *      1) Of the 3 mlat's we are given, at least one positive and negative
     *      value of D results straight away. This is the best case scenario
     *      because we immediately have a bracket on the zero. However, try to
     *      minimize number of I evals needed to determine which pair to take.
     *
     *      2) None of the 3 mlats lead to a negative (negatives are difficult
     *      to capture when the D(mlat) curve only slightly dips below the zero
     *      line. Typically whn this happens, there will be two roots between
     *      mlat0 and mlat2 (and we werent lucky with mlat1). In this case, we
     *      will:
     *              a) First make sure the 3 lats form a triple bracket for
     *              minimization.
     *
     *              b) Use a bisection minimizer to zoom in on the minimim
     *              value of D in the range. Note that this is not the same as
     *              the zero crossing. There are two zeros -- one on either
     *              side of the minimum.
     *
     *              c) Stop doing minimization when we have located a value of
     *              D that is negative.
     *
     *              d) Switch to bisection root finder using the new bracketed
     *              zero.
     *
     */
    Dmin = 9e99;



    /*
     * Attempt to bracket the zero in I-I0. Note that this may not be easy if
     * the curve doesnt drop very substantially below the zero point. This is
     * often the case for large eq. pitch angles.
     */
    Bracket.Dmin = 9e99;
    Flag = BracketZero( I0, Ifound, Bm, MLT, mlat, rad, mlat0, mlat1, mlat2, &Bracket, LstarInfo );
    if ( Flag < 0 )  {
        // should we be bailing here?
        return( -11 );   // An evaluation in BracketZero() hit on a value
    }
                                     //  of I that was undefined -- bail.
    if ( Flag == 2 ) return( TRUE ); // An evaluation in BracketZero() hit on a value
                                     // of I that is within tolerance -- we're done.

    FoundZeroBracket = ( Flag == 1 ) ? TRUE : FALSE; // flag whether we found a zero bracket or not.
    if (LstarInfo->VerbosityLevel > 1){
        if ( FoundZeroBracket ) {
            printf("\t\t\t> Found zero bracket: (a,b) = %g %g  (Da, Db) = %g %g\n", Bracket.a, Bracket.b, Bracket.Da, Bracket.Db );
        } else {
            printf("\t\t\t> No zero bracket found: (a,b,c) = %g %g %g  (Da,Db,Dc) = %g %g %g\n",
                Bracket.a, Bracket.b, Bracket.c,
                Bracket.Da, Bracket.Db, Bracket.Dc );
        }
    }





    /*
     * Use the results of sampling the (mlat,I-I0) space that we have
     * accumulated so far in order to make a quadtratic fit to the "data". Then
     * use the fit to estimate where the zero point is.
     */
a = Bracket.a;
b = Bracket.b;
c = Bracket.c;
Da = Bracket.Da;
Db = Bracket.Db;
Dc = Bracket.Dc;
mlat_min = Bracket.mlat_min;
Dmin = Bracket.Dmin;

// save best value so far
Sign  = (Da < 0.0) ? -1.0 : 1.0;
Dbest = Sign*Dmin;
mlatbest = mlat_min;





    if ( !FoundZeroBracket ) {

        /*
         * We were unsuccesful in finding a zero bracket. Attempt to fit to available values.
         */
        if ( LstarInfo->nImI0>2 ){

            for (i=0; i<LstarInfo->nImI0; i++) LstarInfo->Earr[i] = 1.0;
            FitQuadAndFindZero( LstarInfo->MLATarr, LstarInfo->ImI0arr, LstarInfo->Earr, LstarInfo->nImI0, &res );
            if ( (res > -90.0) && (res < 90.0) ){
                if (LstarInfo->VerbosityLevel > 1){
                    printf("\t\t\t> Fitting to available values. Predicted mlat: %g\n", res );
                }
                a = res-2.0;
                b = res;
                c = res+2.0;

                I  = ComputeI_FromMltMlat( Bm, MLT, b, &r, I0, LstarInfo );
                if ( fabs(I) < 1e98 ) {
                    D0 = I-I0;
                    LstarInfo->MLATarr[LstarInfo->nImI0] = b;
                    LstarInfo->ImI0arr[LstarInfo->nImI0++] = D0;
                    if (fabs(D0) < Dmin){ Dmin = fabs(D0); mlat_min = b; }
                }
                Db = D0;
                if ( Dmin < LstarInfo->mInfo->Lgm_FindShellLine_I_Tol ) {
                    *rad    = r;
                    *Ifound = I;
                    *mlat   = mlat_min;
                    FoundValidI = TRUE;
                    if (LstarInfo->VerbosityLevel > 1){
                        printf( "\t\t\t> Converged with requested tolerance: |I-I0|=%g < %g\n", Dmin, LstarInfo->mInfo->Lgm_FindShellLine_I_Tol );
                    }
                    return( FoundValidI );
                } else if ( D0*Sign < 0.0 ) {

                    FoundZeroBracket = TRUE;
                    if ( Dbest < 0.0 ) {
                        a  = mlatbest;
                        Da = Dbest;
                        b  = b;
                        Db = D0;
                    } else {
                        a  = b;
                        Da = D0;
                        b  = mlatbest;
                        Db = Dbest;
                    }
                    if (LstarInfo->VerbosityLevel > 1){
                        printf("\t\t\t> 1. Zero bracket detected!: (a,b) = %g %g  (Da, Db) = %g %g\n", a, b, Da, Db );
                    }

                }


            

                if ( !FoundZeroBracket ) {


                    I  = ComputeI_FromMltMlat( Bm, MLT, a, &r, I0, LstarInfo );
                    if ( fabs(I) < 1e98 ) {
                        LstarInfo->MLATarr[LstarInfo->nImI0] = a;
                        LstarInfo->ImI0arr[LstarInfo->nImI0++] = I-I0;
                        D0 = I-I0;
                        if (fabs(D0) < Dmin){ Dmin = fabs(D0); mlat_min = a; }
                    }
                    Da = D0;
                    if ( Dmin < LstarInfo->mInfo->Lgm_FindShellLine_I_Tol ) {
                        *rad    = r;
                        *Ifound = I;
                        *mlat   = mlat_min;
                        FoundValidI = TRUE;
                        if (LstarInfo->VerbosityLevel > 1){
                            printf( "\t\t\t> Converged with requested tolerance: |I-I0|=%g < %g\n", Dmin, LstarInfo->mInfo->Lgm_FindShellLine_I_Tol );
                        }
                        return( FoundValidI );
                    } else if ( D0*Sign < 0.0 ) {

                        FoundZeroBracket = TRUE;
                        if ( Dbest < 0.0 ) {
                            a  = mlatbest;
                            Da = Dbest;
                            b  = a;
                            Db = D0;
                        } else {
                            //a  = a;
                            Da = D0;
                            b  = mlatbest;
                            Db = Dbest;
                        }
                        if (LstarInfo->VerbosityLevel > 1){
                            printf("\t\t\t> 2. Zero bracket detected!: (a,b) = %g %g  (Da, Db) = %g %g\n", a, b, Da, Db );
                        }

                    }


                    if ( !FoundZeroBracket ) {

                        I  = ComputeI_FromMltMlat( Bm, MLT, c, &r, I0, LstarInfo );
                        if ( fabs(I) < 1e98 ) {
                            LstarInfo->MLATarr[LstarInfo->nImI0] = c;
                            LstarInfo->ImI0arr[LstarInfo->nImI0++] = I-I0;
                            D0 = I-I0;
                            if (fabs(D0) < Dmin){ Dmin = fabs(D0); mlat_min = c; }
                        }
                        Dc = D0;
                        if ( Dmin < LstarInfo->mInfo->Lgm_FindShellLine_I_Tol ) {
                            *rad    = r;
                            *Ifound = I;
                            *mlat   = mlat_min;
                            FoundValidI = TRUE;
                            if (LstarInfo->VerbosityLevel > 1){
                                printf( "\t\t\t> Converged with requested tolerance: |I-I0|=%g < %g\n", Dmin, LstarInfo->mInfo->Lgm_FindShellLine_I_Tol );
                            }
                            return( FoundValidI );
                        } else if ( D0*Sign < 0.0 ) {

                            FoundZeroBracket = TRUE;
                            if ( Dbest < 0.0 ) {
                                a  = mlatbest;
                                Da = Dbest;
                                b  = c;
                                Db = D0;
                            } else {
                                a  = c;
                                Da = D0;
                                b  = mlatbest;
                                Db = Dbest;
                            }
                            if (LstarInfo->VerbosityLevel > 1){
                                printf("\t\t\t> 3. Zero bracket detected!: (a,b) = %g %g  (Da, Db) = %g %g\n", a, b, Da, Db );
                            }

                        }

                    } //!FoundZeroBracket


                    if (LstarInfo->VerbosityLevel > 1){
                        printf("\t\t\t> New trial range: a,b,c = %g %g %g  Da,Db,Dc = %g %g %g\n", a, b, c, Da, Db, Dc);
                    }


                } //!FoundZeroBracket
            }
        }
    }






    if ( LstarInfo->mInfo->Lgm_FindShellLine_I_Tol > 1e-4 ) {
        BRACKET_EPS = 1e-6;
    } else {
        BRACKET_EPS = 1e-10;
    }



    /*
     *  Use minimization to try and get to the requested tolerance, or get a negatibe i-I0
     */
    if ( !FoundZeroBracket ){

        /*
         * We did not find a Zero bracket. We will have to try the minimization
         * strategy to obtain a value of D that is negative.
         */
         if (LstarInfo->VerbosityLevel > 1){ printf( "\t\t\t> No zero bracket found. Attempting a minimization strategy to get a value of I-I0 < 0 ...\n"); }


        /*
         *  Test to see if we have a valid minimization bracket. If we dont,
         *  then either the minimum is outside the range we chose or there is
         *  no solution.
         */
        if ( (Db > Da) || (Db > Dc) ) {
            *Ifound = 9e99;
            return(-5);
        }


        F0 = 0.5; F1 = 0.5;
        done = FALSE; FoundValidI = FALSE; nIts = 0;
        while ( !done ) {


            /*
             *  Find a new point that is a fraction, F through the largest interval.
             *  Compute I and the difference between I and the desired I (i.e. I0
             *  that the caller gave us)
             */
            d0 = b-a;
            d1 = c-b;
            if ( d0 > d1 ) {
                FirstHalf = TRUE;
                e = a + F0*d0;
            } else {
                FirstHalf = FALSE;
                e = b + F1*d1;
            }
            I = ComputeI_FromMltMlat( Bm, MLT, e, &r, I0, LstarInfo );
            LstarInfo->MLATarr[LstarInfo->nImI0] = e;
            LstarInfo->ImI0arr[LstarInfo->nImI0++] = I-I0;
            De = I-I0;
            if (fabs(De) < Dmin){ Dmin = fabs(De); mlat_min = e; }
            //printf("Initially:  a, b, c, [e]  = %g %g %g [%g]   Da, Db, Dc, [De] = %g %g %g [%g]   Dmin = %g\n", a, b, c, e, Da, Db, Dc, De, Dmin );



            /*
             * Analyze the new value.
             */
            if ( nIts > MAX_ITS ) {
                /*
                 *  Too many iterations. Bail with errror.
                 */
                done = TRUE;
                FoundValidI = -6;
                I = 9e99;
            } else if ( fabs(De) < LstarInfo->mInfo->Lgm_FindShellLine_I_Tol ) {
                /*
                 * Converged with requested tolerance.
                 */
                if (LstarInfo->VerbosityLevel > 1){ printf("\t\t\t> Converged with requested tol. during minimization phase. Requested tol. on I-I0 is: %g  Value of I-I0 achieved: %g\n", LstarInfo->mInfo->Lgm_FindShellLine_I_Tol, De ); }
                done = TRUE;
                FoundValidI = TRUE;
                *mlat = mlat_min;
            } else if ( De < 0.0 ) {
                done = TRUE;
                /*
                 *  Set bracket for finding the zero.
                 *  Lower end should be this new negative value.
                 *  If FirstHalf==TRUE, then b is already set to the value we want.
                 */
                a = e; Da = De;
                if ( !FirstHalf ) {
                    b = c; Db = Dc;
                }
                FoundZeroBracket = TRUE;
                //printf("Found bracket!!!!\n");
                //printf("[a, b] = %g %g   [Da, Db] = %g %g\n", a, b, Da, Db);

                //exit(0);
            } else if ( fabs(e-mlat0)  < BRACKET_EPS ) {
                /*
                 * Converged  to lower endpoint -- no valid mlat found within
                 * interval - probably have to enlarge initial bracket.
                 */
                done = TRUE;
                FoundValidI = -2;
            } else if ( fabs(e-mlat2) < BRACKET_EPS ) {
                /*
                 * Converged  to upper endpoint -- no valid mlat found within
                 * interval - probably have to enlarge initial bracket.
                 */
                done = TRUE;
                FoundValidI = -3;
            } else if ( fabs(a-c) < BRACKET_EPS ) {
                /*
                 * Converged  to something, but not I?  try to enlarge initial
                 * bracket.
                 */
                done = TRUE;
                if ( fabs(De) < 0.001 ) {
                    FoundValidI = TRUE;     // lets just take what we get here....
                    *mlat = mlat_min;       // Use the best value we got
                } else {
                    FoundValidI = -4;
                }
            } else if ( De < Db ) {
                /*
                 * We found a better value for the center value.
                 */
                b  = e;
                Db = De;
                F0 = 0.9; F1 = 0.1;
            //} else if ( FirstHalf && ( fabs(De) < fabs(Da) ) && (De < 0.0) ) {
            } else if ( FirstHalf && ( De < Da )  ) {
                /*
                 * We found a better value for the lower end of the bracket. (And its negative).
                 */
                a  = e;
                Da = De;
                F0 = 0.9; F1 = 0.1;
            //} else if ( !FirstHalf && ( fabs(De) < fabs(Dc) && (De > 0.0) ) ) {
            } else if ( !FirstHalf && ( De < Dc ) ) {
                /*
                 * We found a better value for the upper end of the bracket. (And its positive).
                 */
                c  = e;
                Dc = De;
                F0 = 0.9; F1 = 0.1;
            } else {
                /*
                 * For some reason, the new value is not within the bracket.
                 * (Pathological function behavior).  Since we will eventually bail
                 * out if we do too many iterations, lets try a random selection
                 * for the new point (but still between the current bracket).
                 */
                F0 = rand()/(double)RAND_MAX; // F is in range [0, 1]
                F1 = 1.0 - F0;
            }
            //printf("Setting To: a, b, c, [e]  = %g %g %g [%g]   Da, Db, Dc, [De] = %g %g %g [%g]\n\n", a, b, c, e, Da, Db, Dc, De );

            ++nIts;

        }

    } 






    /*
     *  If we still do not have a good value, and we managed to get a bracket
     *  out of the minimization step, then Use bisecction to zoom-in on the
     *  solution (i.e. get the mlat that make I-I0 = 0 within tolerance.)
     */
    if ( FoundValidI > 0 ) {
        /*
         * We probably converged in minimization phase above.
         */

    } else if ( FoundZeroBracket ){

        /*
         * We have a bracket on a zero value. Go in for the kill using bisection.
         */
        if (LstarInfo->VerbosityLevel > 1){ printf("\t\t\t> Using bisection. Requested tolerance on I-I0 is: %g\n", LstarInfo->mInfo->Lgm_FindShellLine_I_Tol ); }
        F = 0.5;
        done = FALSE; FoundValidI = FALSE; nIts = 0;
        nbFits = 0;
        while ( !done ) {


            /*
             *  Find a new point to evaluate.
             *
             *  If we have less than 3 points so far, take it as a a fraction,
             *  F through the interval.  Compute I and the difference between I
             *  and the desired I (i.e. I0 that the caller gave us).
             *
             *  If we have 3 or more points lets also try a fit -- it may
             *  radically accelerate convergence.
             *
             */
            d = b-a;
            e = a + F*d;
            //if (( LstarInfo->nImI0 > 2 )&&( nbFits%2 )){
            //if (( LstarInfo->nImI0 > 2 )&&( nbFits>5 )){
            //if ( LstarInfo->nImI0 == 6){
            if ( (LstarInfo->nImI0 > 3) && (LstarInfo->nImI0%4 == 0) ){
                for (i=0; i<LstarInfo->nImI0; i++) LstarInfo->Earr[i] = 1.0;
                //FitQuadAndFindZero2( LstarInfo->MLATarr, LstarInfo->ImI0arr, LstarInfo->Earr, LstarInfo->nImI0, 4, &res );
                FitQuadAndFindZero( LstarInfo->MLATarr, LstarInfo->ImI0arr, LstarInfo->Earr, LstarInfo->nImI0, &res );
                if (LstarInfo->VerbosityLevel > 1){
                    printf("\t\t\t> Fitting to available values. Predicted mlat: %g\n", res );
                }
                if ( (res > a) && ( res < b) ) {
                    e = res;
                    ++nbFits;
                }
            }


            I = ComputeI_FromMltMlat( Bm, MLT, e, &r, I0, LstarInfo );
            LstarInfo->MLATarr[LstarInfo->nImI0] = e;
            LstarInfo->ImI0arr[LstarInfo->nImI0++] = I-I0;
            De = I-I0;
            if (fabs(De) < Dmin){ Dmin = fabs(De); mlat_min = e; }
            //printf("Initially:  a, b, [e]  = %g %g [%g]   Da, Db, [De] = %g %g [%g]   Dmin = %g\n", a, b, e, Da, Db, De, Dmin );


            /*
             * Analyze the new value.
             */
            if ( nIts > MAX_ITS ) {
                /*
                 *  Too many iterations. Bail with errror.
                 */
                done = TRUE;
                FoundValidI = -6;
                I = 9e99;
                if (LstarInfo->VerbosityLevel > 1){
                    printf( "\t\t\t> Maximum number of iterations reached (%d)\n", MAX_ITS );
                }
            } else if ( Dmin < LstarInfo->mInfo->Lgm_FindShellLine_I_Tol ) {
                /*
                 * Converged with requested tolerance.
                 */
                //printf("mlat_min= %g Dmin = %g\n", mlat_min, Dmin);
                done = TRUE;
                FoundValidI = TRUE;
                *mlat = mlat_min;
                if (LstarInfo->VerbosityLevel > 1){
                    printf( "\t\t\t> Converged with requested tolerance: |I-I0|=%g < %g\n", Dmin, LstarInfo->mInfo->Lgm_FindShellLine_I_Tol );
                }
            } else if ( fabs(e-mlat0)  < BRACKET_EPS ) {
                /*
                 * Converged  to lower endpoint -- no valid mlat found within
                 * interval - probably have to enlarge initial bracket.
                 */
                done = TRUE;
                FoundValidI = -2;
                if (LstarInfo->VerbosityLevel > 1){
                    printf( "\t\t\t> Converged to lower endpoint: mlat, mlat0 = %g %g\n", e, mlat0 );
                }
            } else if ( fabs(e-mlat1) < BRACKET_EPS ) {
                /*
                 * Converged  to upper endpoint -- no valid mlat found within
                 * interval - probably have to enlarge initial bracket.
                 */
                done = TRUE;
                FoundValidI = -3;
                if (LstarInfo->VerbosityLevel > 1){
                    printf( "\t\t\t> Converged to upper endpoint: mlat, mlat1 = %g %g\n", e, mlat0 );
                }
            } else if ( fabs(a-b) < BRACKET_EPS ) {
                /*
                 * Converged  to something, but not I?  try to enlarge initial
                 * bracket.
                 */
                done = TRUE;
                if ( fabs(De) < 0.001 ) {
                    FoundValidI = TRUE;     // lets just take what we get here....
                    *mlat = mlat_min;       // Use the best value we got
                } else {
                    FoundValidI = -4;
                }
                if (LstarInfo->VerbosityLevel > 1){
                    printf( "\t\t\t> Converged to something (but not I?): fabs(a-b) = %g, but |I-I0| = %g\n", fabs(a-b), De );
                    printf( "\t\t\t> Bracket: a, b, c = %g %g %g    Da, Db, Dc = %g %g %g\n", a, b, c, Da, Db, Dc );
                }
            } else if ( De > 0.0 ) {
                if ( Da > 0.0 ) {
                    // We found a better value for the lower end of the bracket.
                    a  = e;
                    Da = De;
                } else {
                    // We found a better value for the upper end of the bracket.
                    b  = e;
                    Db = De;
                }
            } else if ( De < 0.0 ) {
                if ( Da < 0.0 ) {
                    // We found a better value for the lower end of the bracket.
                    a  = e;
                    Da = De;
                } else {
                    // We found a better value for the upper end of the bracket.
                    b  = e;
                    Db = De;
                }
            } else {
                // WE SHOULD NEVER GET HERE
                /*
                 * For some reason, the new value is not within the bracket.
                 * (Pathological function behavior).  Since we will eventually bail
                 * out if we do too many iterations, lets try a random selection
                 * for the new point (but still between the current bracket).
                 */
                F = 0.01+0.98*rand()/(double)RAND_MAX; // F is in range [0.01, .99]
                if (LstarInfo->VerbosityLevel > 1){
                    printf( "\t\t\t> How did we get here. Trying random bisection: F = %g\n", F );
                }
            }
            //printf("Setting To: a, b, [e]  = %g %g [%g]   Da, Db, [De] = %g %g [%g]\n\n", a, b, e, Da, Db, De );

            ++nIts;

        }


    } else {

        /*
         *  After all this, its unlikely we can find a solution.
         */
        FoundValidI = -7;
        I = 9e99;
        printf( "\t\t\t> Bisection failed to find a root.\n");
    }










    *Iterations = nIts;




    *rad    = r;
    *Ifound = I;
    *mlat   = mlat_min;
    return( FoundValidI );

}



double RadBmFunc( double r, double Bm, void *data ) {

    Lgm_Vector          u, v, Bvec;
    Lgm_MagModelInfo    *mInfo;

    mInfo = (Lgm_MagModelInfo *)data;

    u.x = r*mInfo->Ptmp.x; u.y = r*mInfo->Ptmp.y; u.z = r*mInfo->Ptmp.z;
    Lgm_Convert_Coords( &u, &v, SM_TO_GSM, mInfo->c );
    mInfo->Bfield( &v, &Bvec, mInfo );
    return( Lgm_Magnitude( &Bvec ) - Bm );

}


/** Find radius of specified Bm for a given MLT and mlat.
 *
 *  This routine is used by FindShellLine(). For the given MLT and mlat, it
 *  locates the vertical
 *  radius where the field strength is equal to the specified mirror field
 *  strength, Bm.
 *
 *  Need to add error handling for cases where we cant find a radius. This
 *  could happen if the actual radius lies outside of the assumed bracket.
 *  E.g., could drop below Loss Cone height.
 *
 *            \param[in]        Bm          the value of the mirror field strength for the desired drift shell.
 *            \param[in]        MLT         magnetic local time.
 *            \param[in]        mlat        magnetic latitude.
 *            \param[out]       rad         the geocentric radial distance of the mirror point (at which B=Bm).
 *            \param[in]        tol         Tolerance.
 *            \param[in,out]    LstarInfo   A properly initialized Lgm_LstarInfo structure.
 *
 *            \return 1 - normal exit. Bm value found.
 *            \return 0 - no Bm value found.
 *
 */
int FindBmRadius( double Bm, double MLT, double mlat, double *r, double tol, Lgm_LstarInfo *LstarInfo ) {

    Lgm_Vector  u, v, Bvec;
    int         done, FoundValidBm, Flag;
    double      Phi, D, a0, c0, Da0, Dc0;
    double      a, b, c, d, B, cl, sl, cp, sp, lat, f, g;

    LstarInfo->mInfo->nFunc = 0;

    Phi = 15.0*(MLT-12.0)*RadPerDeg; lat = mlat * RadPerDeg;
    cl = cos( lat ); sl = sin( lat );
    cp = cos( Phi ); sp = sin( Phi );
    f = cl*cp; g = cl*sp;


    /*
     *  Get bracket on r.
     *  Originally, I thought we could just set this bracket straight away
     *  without actually do a search.  But it would appear that this is not the
     *  case. Lets search for it in 0.5 Re increments. Assume a0, but lets step
     *  out to find c0.
     */
    a0 = 1.0 + LstarInfo->mInfo->Lgm_LossConeHeight/Re;      // 110 km altitude
//a0 = 1.0;

    u.x = a0*f; u.y = a0*g; u.z = a0*sl;
    Lgm_Convert_Coords( &u, &v, SM_TO_GSM, LstarInfo->mInfo->c );
    LstarInfo->mInfo->Bfield( &v, &Bvec, LstarInfo->mInfo );
    B = Lgm_Magnitude( &Bvec );
    Da0 = B - Bm;

    //c0 = 20.0;                // this ought to be big enough!?
    //c0 = 15.0;                // this ought to be big enough!?
    c0 = a0;
    done = FALSE;
    while ( !done ) {

        // Increment c0 and recompute Dc0.
        c0 += 1.0;
        u.x = c0*f; u.y = c0*g; u.z = c0*sl;
        Lgm_Convert_Coords( &u, &v, SM_TO_GSM, LstarInfo->mInfo->c );
        LstarInfo->mInfo->Bfield( &v, &Bvec, LstarInfo->mInfo );
        B = Lgm_Magnitude( &Bvec );
        Dc0 = B - Bm;

        // Check to see if we went beyond Bm
        if ( Dc0 < 0.0 ) {
            done = TRUE;
        } else {
            a0  = c0;
            Da0 = Dc0;
        }


        //printf("c0 = %g     B = %g    Bm = %g\n", c0, B, Bm);
        // We may never be able to get a B value low enough...
        if (c0>20.0) return(FALSE);
    }




    if (LstarInfo->VerbosityLevel > 5) {
        printf( "%sFindBmRadius: a0, c0 = %g %g   Da0, Dc0 = %g %g%s\n", LstarInfo->PreStr, a0, c0, Da0, Dc0, LstarInfo->PostStr  );
    }


    /*
     *  For the search to work, Da0 (B-Bm at a0) must be positive to start with.
     *  If it isnt it means the value of Bm we are looking for is below the lowest
     *  point we are willing to allow (i.e. its in the loss cone). Also, (for whatever
     *  reason) if Dc0 is > 0.0 we wont be able to find the point because its higher than our
     *  initial c0.
     */
    if ( ( Da0 < 0.0 ) || ( Dc0 > 0.0 ) ) return( FALSE ); // bracket not found



    /*
     * We have a bracket. Use brent's method to find root.
     */
    LstarInfo->mInfo->Ptmp.x = f;
    LstarInfo->mInfo->Ptmp.y = g;
    LstarInfo->mInfo->Ptmp.z = sl;
    BrentFuncInfo bfi;
    bfi.Info    = (void *)LstarInfo->mInfo;
    bfi.func    = &RadBmFunc;
    bfi.Val     = Bm;
    Flag = Lgm_zBrent( a0, c0, Da0, Dc0, &bfi, tol, r, &D );
//printf("tol = %g\n", tol);
    FoundValidBm = (Flag) ? TRUE : FALSE;

//printf("GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG     D = %g\n", D);
    if (LstarInfo->VerbosityLevel > 5) {
      printf("%sFindBmRadius: Final r = %.15lf  (B-Bm = %g nFunc = %n)%s\n", LstarInfo->PreStr, *r, D, LstarInfo->mInfo->nFunc, LstarInfo->PostStr );
    }

    return( FoundValidBm );

}

int FitLineAndFindZero( double *x, double *y, double *dy, int n, double *res ) {

    int     i, Flag;
    double  chisq, root, A, B, C, D;

    gsl_multifit_linear_workspace   *work;
    gsl_matrix                      *X, *cov;
    gsl_vector                      *v, *w, *c, *xx;


    X   = gsl_matrix_alloc( n, 2 );
    v   = gsl_vector_alloc( n );
    w   = gsl_vector_alloc( n );
    c   = gsl_vector_alloc(2);
    cov = gsl_matrix_alloc(2, 2);
    xx  = gsl_vector_alloc(2);

    for ( i=0; i<n; i++ ) {

        gsl_matrix_set( X, i, 0, 1.0 );
        gsl_matrix_set( X, i, 1, x[i] );

        gsl_vector_set( v, i, y[i] );
        gsl_vector_set( w, i, 1.0/(dy[i]*dy[i]) );

    }


    work = gsl_multifit_linear_alloc( n, 2 );
    gsl_multifit_wlinear( X, w, v, c, cov, &chisq, work );
    gsl_multifit_linear_free( work );

    // y = b*x + c
    C = gsl_vector_get( c, 0 );
    B = gsl_vector_get( c, 1 );



    Flag = 0;
    if ( fabs(B) > 0.0 ) {
        *res = -C/B;
        Flag = 1;
    } else {
        *res = -1e31;
    }


    gsl_matrix_free( X );
    gsl_vector_free( xx );
    gsl_vector_free( v );
    gsl_vector_free( w );
    gsl_vector_free( c );
    gsl_matrix_free( cov );


    return( Flag );


}


int FitQuadAndFindZero( double *x, double *y, double *dy, int n, double *res ) {

    int     i, Flag;
    double  chisq, root, A, B, C, D;

    gsl_multifit_linear_workspace   *work;
    gsl_matrix                      *X, *cov;
    gsl_vector                      *v, *w, *c, *xx;


    X   = gsl_matrix_alloc( n, 3 );
    v   = gsl_vector_alloc( n );
    w   = gsl_vector_alloc( n );
    c   = gsl_vector_alloc (3);
    cov = gsl_matrix_alloc (3, 3);
    xx  = gsl_vector_alloc (3);

    for ( i=0; i<n; i++ ) {

        gsl_matrix_set( X, i, 0, 1.0 );
        gsl_matrix_set( X, i, 1, x[i] );
        gsl_matrix_set( X, i, 2, x[i]*x[i] );

        gsl_vector_set( v, i, y[i] );
        gsl_vector_set( w, i, 1.0/(dy[i]*dy[i]) );

    }


    work = gsl_multifit_linear_alloc( n, 3 );
    gsl_multifit_wlinear( X, w, v, c, cov, &chisq, work );
    gsl_multifit_linear_free( work );

    C = gsl_vector_get( c, 0 );
    B = gsl_vector_get( c, 1 );
    A = gsl_vector_get( c, 2 );


    Flag = 0;
    D = B*B - 4.0*A*C;
    if ( D >= 0.0 ) {
        *res  = (-B + sqrt( D ))/(2.0*A);
//printf("A, B, C = %g %g %g   res = %g\n", A, B, C, *res);
//        gsl_vector_set( xx, 0, 1.0 );
//        gsl_vector_set( xx, 1, root );
//        gsl_vector_set( xx, 2, root*root );
//        gsl_multifit_linear_est( xx, c, cov, res, res_err );
        Flag = 1;
    } else {
        *res = -1e31;
    }


    gsl_matrix_free( X );
    gsl_vector_free( xx );
    gsl_vector_free( v );
    gsl_vector_free( w );
    gsl_vector_free( c );
    gsl_matrix_free( cov );


    return( Flag );


}


/*
 * same as above, but only use the best 4 vals
 */
int FitQuadAndFindZero2( double *xin, double *yin, double *dyin, int n, int nmax, double *res ) {

    int         i, Flag;
    double      chisq, root, A, B, C, D;
    double      *x, *y, *dy, *fabsy;
    FitVals     *f;


    gsl_multifit_linear_workspace   *work;
    gsl_matrix                      *X, *cov;
    gsl_vector                      *v, *w, *c, *xx;

    LGM_ARRAY_1D( f, n, FitVals );
    for (i=0; i<n; i++){
        f[i].x   = xin[i];
        f[i].y   = yin[i];
        f[i].dy  = dyin[i];
        f[i].key = fabs(yin[i]);
        //printf("|y|, y, x, = %g %g %g\n", fabs(yin[i]), yin[i], xin[i] );
    }
    elt_qsort( f, n );

    n = nmax;
    LGM_ARRAY_1D( x, n, double );
    LGM_ARRAY_1D( y, n, double );
    LGM_ARRAY_1D( dy, n, double );
    for (i=0; i<n; i++){
        x[i] = f[i].x;
        y[i] = f[i].y;
        dy[i] = f[i].dy;
    }


    X   = gsl_matrix_alloc( n, 3 );
    v   = gsl_vector_alloc( n );
    w   = gsl_vector_alloc( n );
    c   = gsl_vector_alloc (3);
    cov = gsl_matrix_alloc (3, 3);
    xx  = gsl_vector_alloc (3);

    for ( i=0; i<n; i++ ) {

        gsl_matrix_set( X, i, 0, 1.0 );
        gsl_matrix_set( X, i, 1, x[i] );
        gsl_matrix_set( X, i, 2, x[i]*x[i] );

        gsl_vector_set( v, i, y[i] );
        gsl_vector_set( w, i, 1.0/(dy[i]*dy[i]) );

    }


    work = gsl_multifit_linear_alloc( n, 3 );
    gsl_multifit_wlinear( X, w, v, c, cov, &chisq, work );
    gsl_multifit_linear_free( work );

    C = gsl_vector_get( c, 0 );
    B = gsl_vector_get( c, 1 );
    A = gsl_vector_get( c, 2 );


    Flag = 0;
    D = B*B - 4.0*A*C;
    if ( D >= 0.0 ) {
        *res  = (-B + sqrt( D ))/(2.0*A);
//        gsl_vector_set( xx, 0, 1.0 );
//        gsl_vector_set( xx, 1, root );
//        gsl_vector_set( xx, 2, root*root );
//        gsl_multifit_linear_est( xx, c, cov, res, res_err );
        Flag = 1;
    } else {
        *res = -1e31;
    }


    gsl_matrix_free( X );
    gsl_vector_free( xx );
    gsl_vector_free( v );
    gsl_vector_free( w );
    gsl_vector_free( c );
    gsl_matrix_free( cov );


    LGM_ARRAY_1D_FREE( f );
    LGM_ARRAY_1D_FREE( x );
    LGM_ARRAY_1D_FREE( y );
    LGM_ARRAY_1D_FREE( dy );

    return( Flag );


}


int BracketZero( double I0, double *Ifound, double Bm, double MLT, double *mlat, double *rad, double mlat0, double mlat1, double mlat2, BracketType *Bracket, Lgm_LstarInfo *LstarInfo ){

    double  I, r, D0, D1, D2, Dmin, mlat_min;
    int     FoundZeroBracket;

    Dmin = Bracket->Dmin;


    /*
     *    I( mlat1 )
     *
     * Compute I-I0 at best-guess mlat. If the caller predicted well, then this
     * may be dead on, so try it first.
     */
    I  = ComputeI_FromMltMlat( Bm, MLT, mlat1, &r, I0, LstarInfo );
    if ( fabs(I) > 1e99 ) return(-5);
    D1 = I-I0;
    LstarInfo->MLATarr[LstarInfo->nImI0]   = mlat1;
    LstarInfo->ImI0arr[LstarInfo->nImI0++] = D1;
    if (fabs(D1) < Dmin){ Dmin = fabs(D1); mlat_min = mlat1; }
    Bracket->Dmin = Dmin; Bracket->mlat_min = mlat_min;

    if ( Dmin < LstarInfo->mInfo->Lgm_FindShellLine_I_Tol ) {
        /*
         *  Already Converged with requested tolerance.
         */
        //printf("mlat_min= %g Dmin = %g\n", mlat_min, Dmin);
        *rad    = r;
        *Ifound = I;
        *mlat   = mlat_min;

        //printf("Dmin = %g   LstarInfo->mInfo->Lgm_FindShellLine_I_Tol = %g\n", Dmin, LstarInfo->mInfo->Lgm_FindShellLine_I_Tol );
        return( 2 );
    }
    //printf("Bracket->Dmin = %g\n", Bracket->Dmin);




    /*
     * Attempt to bracket the zero in I-I0. Note that this may not be easy if
     * the curve doesnt drop very substantially below the zero point. This is
     * often the case for large eq. pitch angles.
     */
    FoundZeroBracket = FALSE;
    if ( D1 > 0.0 ) {
        /*
         *    I( mlat0 )
         *
         *  Then we would like to have the other side of the bracket be < 0.0.
         *  I.e. we need a smaller I which (usually) is a smaller mlat. So try
         *  mlat0 next.
         */
        I  = ComputeI_FromMltMlat( Bm, MLT, mlat0, &r, I0, LstarInfo );
        if ( fabs(I) > 1e99 ) return(-5);
        D0 = I-I0;
        LstarInfo->MLATarr[LstarInfo->nImI0] = mlat0;
        LstarInfo->ImI0arr[LstarInfo->nImI0++] = D0;
        if (fabs(D0) < Dmin){ Dmin = fabs(D0); mlat_min = mlat0; }
        Bracket->Dmin = Dmin; Bracket->mlat_min = mlat_min;

        if ( Dmin < LstarInfo->mInfo->Lgm_FindShellLine_I_Tol ) {
            /*
             * Already Converged with requested tolerance.
             */
            //printf("mlat_min= %g Dmin = %g\n", mlat_min, Dmin);
            *rad    = r;
            *Ifound = I;
            *mlat   = mlat_min;
            return( 2 );
        }

        if ( D0 < 0.0 ) {

            /*
             *  Yay, Found Zero Bracket!!!
             */
            FoundZeroBracket = TRUE;
            Bracket->a = mlat0; Bracket->Da = D0;
            Bracket->b = mlat1; Bracket->Db = D1;

        } else {

            /*
             *    I( mlat2 )
             *
             * Still no bracket. Evaluate the other side.
             */
            I  = ComputeI_FromMltMlat( Bm, MLT, mlat2, &r, I0, LstarInfo );
            if ( fabs(I) > 1e99 ) return(-5);
            D2 = I-I0;
            LstarInfo->MLATarr[LstarInfo->nImI0] = mlat2;
            LstarInfo->ImI0arr[LstarInfo->nImI0++] = D2;
            if (fabs(D2) < Dmin){ Dmin = fabs(D2); mlat_min = mlat2; }
            Bracket->Dmin = Dmin; Bracket->mlat_min = mlat_min;

            if ( Dmin < LstarInfo->mInfo->Lgm_FindShellLine_I_Tol ) {
                /*
                 * Already Converged with requested tolerance.
                 */
                //printf("mlat_min= %g Dmin = %g\n", mlat_min, Dmin);
                *rad    = r;
                *Ifound = I;
                *mlat   = mlat_min;
                return( 2 );
            }

            if ( D2 < 0.0 ) {
                /*
                 *  Yay, Found Zero Bracket!!!
                 */
                FoundZeroBracket = TRUE;
                Bracket->a = mlat1; Bracket->Da = D1;
                Bracket->b = mlat2; Bracket->Db = D2;
            } else {
                /*
                 *  The three point dont provide a zero bracket. We need to
                 *  treat the problem as a minimization instead.  Set up a
                 *  potential bracket for minimization.
                 */
                Bracket->a = mlat0; Bracket->Da = D0;
                Bracket->b = mlat1; Bracket->Db = D1;
                Bracket->c = mlat2; Bracket->Dc = D2;
            }

        }

    } else {  // D1 is < 0.0

        /*
         *  Then we would like to have the other side of the bracket be > 0.0.
         *  I.e. we need a bigger I which (usually) is a bigger mlat.
         */
        I  = ComputeI_FromMltMlat( Bm, MLT, mlat2, &r, I0, LstarInfo );
        if ( fabs(I) > 1e99 ) return(-5);
        D2 = I-I0;
        LstarInfo->MLATarr[LstarInfo->nImI0]   = mlat2;
        LstarInfo->ImI0arr[LstarInfo->nImI0++] = D2;
        if (fabs(D2) < Dmin){ Dmin = fabs(D2); mlat_min = mlat2; }
        Bracket->Dmin = Dmin; Bracket->mlat_min = mlat_min;

        if ( Dmin < LstarInfo->mInfo->Lgm_FindShellLine_I_Tol ) {
            /*
             * Already Converged with requested tolerance.
             */
            //printf("mlat_min= %g Dmin = %g\n", mlat_min, Dmin);
            *rad    = r;
            *Ifound = I;
            *mlat   = mlat_min;
            return( 2 );
        }

        if ( D2 > 0.0 ) {

            /*
             *  Yay, Found Zero Bracket!!!
             */
            FoundZeroBracket = TRUE;
            Bracket->a = mlat1; Bracket->Da = D1;
            Bracket->b = mlat2; Bracket->Db = D2;

        } else {

            /*
             * Still no bracket. Evaluate the other side.
             */
            I  = ComputeI_FromMltMlat( Bm, MLT, mlat0, &r, I0, LstarInfo );
            if ( fabs(I) > 1e99 ) return(-5);
            D0 = I-I0;
            LstarInfo->MLATarr[LstarInfo->nImI0]   = mlat0;
            LstarInfo->ImI0arr[LstarInfo->nImI0++] = D0;
            if (fabs(D0) < Dmin){ Dmin = fabs(D0); mlat_min = mlat0; }
            Bracket->Dmin = Dmin; Bracket->mlat_min = mlat_min;

            if ( Dmin < LstarInfo->mInfo->Lgm_FindShellLine_I_Tol ) {
                /*
                 * Already Converged with requested tolerance.
                 */
                //printf("mlat_min= %g Dmin = %g\n", mlat_min, Dmin);
                *rad    = r;
                *Ifound = I;
                *mlat   = mlat_min;
                return( 2 );
            }

            if ( D0 > 0.0 ) {
                /*
                 *  Yay, Found Zero Bracket!!!
                 */
                FoundZeroBracket = TRUE;
                Bracket->a = mlat0; Bracket->Da = D0;
                Bracket->b = mlat1; Bracket->Db = D1;
            } else {
                /*
                 *  The three points dont provide a zero bracket. We may need
                 *  to treat the problem as a minimization instead.  Set up a
                 *  potential bracket for minimization.
                 */
                Bracket->a = mlat0; Bracket->Da = D0;
                Bracket->b = mlat1; Bracket->Db = D1;
                Bracket->c = mlat2; Bracket->Dc = D2;
            }

        }

    }

    return( FoundZeroBracket );

}
