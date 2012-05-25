#include "Lgm/Lgm_MagModelInfo.h"
#include "Lgm/Lgm_LstarInfo.h"
#include <gsl/gsl_multifit.h>
#define MAX_ITS 100

int FitQuadAndFindZero( double *x, double *y, double *dy, int n, double *res );

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
    Lgm_Vector    u, w, Pm_North, Pmirror, v1, v2, v3;
    double        F, F0, F1, rat, a, b, c, d, d0, d1, Da, Db, Dc, De, I, r, Phi, cl, sl;
    double        SS, Sn, Ss, mlat_min=0.0, Dmin=9e99, e, D, D0, D1, D2;
    int           done, FoundValidI, FirstHalf, nIts, FoundZeroBracket;
    int           i;

    *Iterations = 0;


    /*
     *  Set the bracket for the mlat range. We got mlat0 and mlat2 from the
     *  caller. The assumption is that the field line with the I we are
     *  looking fo is between these two latitudes.  (If it isnt we will
     *  eventually have to bail and the user has to try a different range.) We
     *  also got mlat1 which is the user's best guess at where it ought to be.
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
     * Compute I-I0 at best-guess mlat. If the caller predicted well, then
     * this may be dead on, so try it first.
     */
    I  = ComputeI_FromMltMlat( Bm, MLT, mlat1, &r, I0, LstarInfo );
    if ( fabs(I) > 1e99 ) return(-5);
    LstarInfo->MLATarr[LstarInfo->nImI0] = mlat1;
    LstarInfo->ImI0arr[LstarInfo->nImI0++] = I-I0;
    D1 = I-I0;
    if (fabs(D1) < Dmin){ Dmin = fabs(D1); mlat_min = mlat1; }

    if ( Dmin < LstarInfo->mInfo->Lgm_FindShellLine_I_Tol ) {
        /*
         * Already Converged with requested tolerance.
         */
        //printf("mlat_min= %g Dmin = %g\n", mlat_min, Dmin);
        *rad    = r;
        *Ifound = I;
        *mlat   = mlat_min;
        FoundValidI = TRUE;

        return( FoundValidI );
    }

//printf("Dmin = %g\n", Dmin);

     


    /*
     * Attempt to bracket the zero in I-I0. Note that this may not be easy if
     * the curve doesnt drop very substantially below the zero point. This is
     * often the case for large eq. pitch angles.
     */
    FoundZeroBracket = FALSE;
    if ( D1 > 0.0 ) {
        /*
         *  Then we would like to have the other side of the bracket be < 0.0.
         *  I.e. we need a smaller I which (usually) is a smaller mlat. So try
         *  mlat0 next.
         */
        I  = ComputeI_FromMltMlat( Bm, MLT, mlat0, &r, I0, LstarInfo );
        if ( fabs(I) > 1e99 ) return(-5);
        LstarInfo->MLATarr[LstarInfo->nImI0] = mlat0;
        LstarInfo->ImI0arr[LstarInfo->nImI0++] = I-I0;
        D0 = I-I0;
        if (fabs(D0) < Dmin){ Dmin = fabs(D0); mlat_min = mlat0; }

        if ( Dmin < LstarInfo->mInfo->Lgm_FindShellLine_I_Tol ) {
            /*
             * Already Converged with requested tolerance.
             */
            //printf("mlat_min= %g Dmin = %g\n", mlat_min, Dmin);
            *rad    = r;
            *Ifound = I;
            *mlat   = mlat_min;
            FoundValidI = TRUE;
            return( FoundValidI );
        }

        if ( D0 < 0.0 ) {

            /*
             *  Yay, Found Zero Bracket!!!
             */
            FoundZeroBracket = TRUE;
            a = mlat0; Da = D0;
            b = mlat1; Db = D1;

        } else {

            /*
             * Still no bracket. Evaluate the other side.
             */
            I  = ComputeI_FromMltMlat( Bm, MLT, mlat2, &r, I0, LstarInfo );
            if ( fabs(I) > 1e99 ) return(-5);
            LstarInfo->MLATarr[LstarInfo->nImI0] = mlat2;
            LstarInfo->ImI0arr[LstarInfo->nImI0++] = I-I0;
            D2 = I-I0;
            if (fabs(D2) < Dmin){ Dmin = fabs(D2); mlat_min = mlat2; }

            if ( Dmin < LstarInfo->mInfo->Lgm_FindShellLine_I_Tol ) {
                /*
                 * Already Converged with requested tolerance.
                 */
                //printf("mlat_min= %g Dmin = %g\n", mlat_min, Dmin);
                *rad    = r;
                *Ifound = I;
                *mlat   = mlat_min;
                FoundValidI = TRUE;
                return( FoundValidI );
            }

            if ( D2 < 0.0 ) {
                /*
                 *  Yay, Found Zero Bracket!!!
                 */
                FoundZeroBracket = TRUE;
                a = mlat1; Da = D1;
                b = mlat2; Db = D2;
            } else {
                /*
                 *  The three point dont provide a zero bracket. We need to
                 *  treat the problem as a minimization instead.  Set up a
                 *  potential bracket for minimization.
                 */
                a = mlat0; Da = D0;
                b = mlat1; Db = D1;
                c = mlat2; Dc = D2;
            }

        }

    } else {  // D1 is < 0.0

        /*
         *  Then we would like to have the other side of the bracket be > 0.0.
         *  I.e. we need a bigger I which (usually) is a bigger mlat.
         */
        I  = ComputeI_FromMltMlat( Bm, MLT, mlat2, &r, I0, LstarInfo );
        if ( fabs(I) > 1e99 ) return(-5);
        LstarInfo->MLATarr[LstarInfo->nImI0] = mlat2;
        LstarInfo->ImI0arr[LstarInfo->nImI0++] = I-I0;
        D2 = I-I0;
        if (fabs(D2) < Dmin){ Dmin = fabs(D2); mlat_min = mlat2; }

        if ( Dmin < LstarInfo->mInfo->Lgm_FindShellLine_I_Tol ) {
            /*
             * Already Converged with requested tolerance.
             */
            //printf("mlat_min= %g Dmin = %g\n", mlat_min, Dmin);
            *rad    = r;
            *Ifound = I;
            *mlat   = mlat_min;
            FoundValidI = TRUE;
            return( FoundValidI );
        }

        if ( D2 > 0.0 ) {

            /*
             *  Yay, Found Zero Bracket!!!
             */
            FoundZeroBracket = TRUE;
            a = mlat1; Da = D1;
            b = mlat2; Db = D2;

        } else {

            /*
             * Still no bracket. Evaluate the other side.
             */
            I  = ComputeI_FromMltMlat( Bm, MLT, mlat0, &r, I0, LstarInfo );
            if ( fabs(I) > 1e99 ) return(-5);
            LstarInfo->MLATarr[LstarInfo->nImI0] = mlat0;
            LstarInfo->ImI0arr[LstarInfo->nImI0++] = I-I0;
            D0 = I-I0;
            if (fabs(D0) < Dmin){ Dmin = fabs(D0); mlat_min = mlat0; }

            if ( Dmin < LstarInfo->mInfo->Lgm_FindShellLine_I_Tol ) {
                /*
                 * Already Converged with requested tolerance.
                 */
                //printf("mlat_min= %g Dmin = %g\n", mlat_min, Dmin);
                *rad    = r;
                *Ifound = I;
                *mlat   = mlat_min;
                FoundValidI = TRUE;
                return( FoundValidI );
            }

            if ( D0 > 0.0 ) {
                /*
                 *  Yay, Found Zero Bracket!!!
                 */
                FoundZeroBracket = TRUE;
                a = mlat0; Da = D0;
                b = mlat1; Db = D1;
            } else {
                /*
                 *  The three points dont provide a zero bracket. We need to
                 *  treat the problem as a minimization instead.  Set up a
                 *  potential bracket for minimization.
                 */
                a = mlat0; Da = D0;
                b = mlat1; Db = D1;
                c = mlat2; Dc = D2;
            }

        }

    }

//printf("Dmin = %g\n", Dmin);


for (i=0; i<LstarInfo->nImI0; i++) {
//printf("%g %g\n", LstarInfo->MLATarr[i], LstarInfo->ImI0arr[i] );
LstarInfo->Earr[i] = 1.0;
}

if (1==1){
double res;
if ( LstarInfo->nImI0>2 ){
FitQuadAndFindZero( LstarInfo->MLATarr, LstarInfo->ImI0arr, LstarInfo->Earr, LstarInfo->nImI0, &res );
if ( (res > -90.0) && (res < 90.0) ){
printf( "res = %g\n", res );
a = res-2.0;
b = res;
c = res+2.0;

I  = ComputeI_FromMltMlat( Bm, MLT, b, &r, I0, LstarInfo );
if ( fabs(I) < 1e98 ) {
    LstarInfo->MLATarr[LstarInfo->nImI0] = b;
    LstarInfo->ImI0arr[LstarInfo->nImI0++] = I-I0;
    D0 = I-I0;
    if (fabs(D0) < Dmin){ Dmin = fabs(D0); mlat_min = b; }
}
Db = D0;
if ( Dmin < LstarInfo->mInfo->Lgm_FindShellLine_I_Tol ) {
    *rad    = r;
    *Ifound = I;
    *mlat   = mlat_min;
    FoundValidI = TRUE;
    return( FoundValidI );
}

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
    return( FoundValidI );
}



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
    return( FoundValidI );
}
printf("new trial range: a,b,c = %g %g %g  Da, Db, Dc = %g %g %g\n", a, b, c, Da, Db, Dc);

}

}

}
    

    if ( !FoundZeroBracket ){

//printf("1. HERE\n");
        /*
         * We did not find a Zero bracket. We will have to try the minimization
         * strategy to obtain a value of D that is negative.
         */


        /*
         *  Test to see if we have a valid minimization bracket. If we dont,
         *  then either the minimum is outside the range we chose or there is
         *  no solution.
         */
        if ( (Db > Da) || (Db > Dc) ) {
            *Ifound = 9e99;
            return(-5);
        }

//printf("2. HERE\n");

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
                //printf("mlat_min= %g Dmin = %g\n", mlat_min, Dmin);
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
            } else if ( fabs(e-mlat0)  < 1e-5 ) {
                /*
                 * Converged  to lower endpoint -- no valid mlat found within
                 * interval - probably have to enlarge initial bracket.
                 */
                done = TRUE;
                FoundValidI = -2;
            } else if ( fabs(e-mlat2) < 1e-5 ) {
                /*
                 * Converged  to upper endpoint -- no valid mlat found within
                 * interval - probably have to enlarge initial bracket.
                 */
                done = TRUE;
                FoundValidI = -3;
            } else if ( fabs(a-c) < 1e-8 ) {
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
//printf("3. HERE\n");


/*
FILE *fp;
double bb;
fp = fopen("crap.txt", "w");
for(bb=a; bb<=b; bb+=0.01){
fprintf(fp, "%g %g\n", bb, ComputeI_FromMltMlat( Bm, MLT, bb, &r, I0, LstarInfo ) - I0);
}
fclose(fp);
exit(0);
*/









    if ( FoundZeroBracket ) {

        /*
         * We have a bracket on a zero value. Go in for the kill using bisection.
         */
        F = 0.5;
        done = FALSE; FoundValidI = FALSE; nIts = 0;
        while ( !done ) {
            /*
             *  Find a new point that is a fraction, F through the interval.
             *  Compute I and the difference between I and the desired I (i.e. I0
             *  that the caller gave us)
             */
            d = b-a;
            e = a + F*d;
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
            } else if ( Dmin < LstarInfo->mInfo->Lgm_FindShellLine_I_Tol ) {
                /*
                 * Converged with requested tolerance.
                 */
                //printf("mlat_min= %g Dmin = %g\n", mlat_min, Dmin);
                done = TRUE;
                FoundValidI = TRUE;
                *mlat = mlat_min;
            } else if ( fabs(e-mlat0)  < 1e-5 ) {
                /*
                 * Converged  to lower endpoint -- no valid mlat found within
                 * interval - probably have to enlarge initial bracket.
                 */
                done = TRUE;
                FoundValidI = -2;
            } else if ( fabs(e-mlat1) < 1e-5 ) {
                /*
                 * Converged  to upper endpoint -- no valid mlat found within
                 * interval - probably have to enlarge initial bracket.
                 */
                done = TRUE;
                FoundValidI = -3;
            } else if ( fabs(a-b) < 1e-8 ) {
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
            //} else if ( (De < 0.0) && (De > Da)  ) {
            } else if ( De < 0.0 ) {
                /*
                 * We found a better value for the lower end of the bracket.
                 */
                a  = e;
                Da = De;
            } else if ( De > 0.0 ) {
            //} else if ( (De > 0.0) && (De < Db)  ) {
                /*
                 * We found a better value for the lower end of the bracket.
                 */
                b  = e;
                Db = De;
            } else {
                // WE SHOULD NEVER GET HERE
                /*
                 * For some reason, the new value is not within the bracket.
                 * (Pathological function behavior).  Since we will eventually bail
                 * out if we do too many iterations, lets try a random selection
                 * for the new point (but still between the current bracket).
                 */
                F = 0.01+0.98*rand()/(double)RAND_MAX; // F is in range [0.01, .99]
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
