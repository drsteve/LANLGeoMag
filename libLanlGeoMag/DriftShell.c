#include "Lgm/Lgm_MagModelInfo.h"
#include "Lgm/Lgm_LstarInfo.h"
#define MAX_ITS 100


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
int FindShellLine(  double I0, double *Ifound, double Bm, double MLT, double *mlat, double *rad, double mlat0, double mlat1, Lgm_LstarInfo *LstarInfo) {
    Lgm_Vector    u, w, Pm_North, Pmirror, v1, v2, v3;
    double        F0, F1, rat, a, b, c, d, d0, d1, Da, Db, Dc, De, I, r, Phi, cl, sl, SS, Sn, Ss, mlat_min=0.0, Dmin=9e99, e, D;
    int           done, FoundValidI, FirstHalf, nIts;


    /*
     *   Set the bracket for the mlat range. The safest values here would be 0->90 deg.
     *   since if a shell field line exists for the given I,Bm values, it must have Bm
     *   between 0 and 90. (Is this strictly true for very low L-shells that involve FLs
     *   near the equator?)
     *
     *   A problem with using 0->90 is that long FLs take some time to integrate. The
     *   search algorithm moves away from these FLs in a few steps, but those steps cost
     *   us alot.
     *
     *   A compromise I have adopted here is to allow the user to set the bracket manually.
     *   Be careful because there *will* be problems if the bracket is no good (i.e. its not
     *   a true bracket).
     *
     */
    Dmin = 9e99;


    /*
     * Compute I-I0 at midpoint.
     */
    b  = 0.5*(mlat0+mlat1);
    I  = ComputeI_FromMltMlat( Bm, MLT, b, &r, I0, LstarInfo );
    if ( fabs(I) > 1e99 ) return(-5);
    Db = I-I0;
    if (fabs(Db) < Dmin){ Dmin = fabs(Db); mlat_min = b; }

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




    /*
     * Compute I-I0 at lower side of potential bracket.
     */
    a  = mlat0;
    I  = ComputeI_FromMltMlat( Bm, MLT, a, &r, I0, LstarInfo );
    if ( fabs(I) > 1e99 ) return(-5);
    Da = I-I0;
    if (fabs(Da) < Dmin){ Dmin = fabs(Da); mlat_min = b; }

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



    /*
     * Compute I-I0 at upper side of potential bracket.
     */
    c  = mlat1;
    I  = ComputeI_FromMltMlat( Bm, MLT, c, &r, I0, LstarInfo );
    if ( fabs(I) > 1e99 ) return(-5);
    Dc = I-I0;
    if (fabs(Dc) < Dmin){ Dmin = fabs(Dc); mlat_min = c; }

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







    /*
     *  Test to see if we have a valid bracket. If we dont, then either the
     *  minimum is outside the range we chose or there is no solution.
     */
    if ( (fabs(Db) > fabs(Da)) || (fabs(Db) > fabs(Dc)) ) {
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
        } else if ( fabs(De) < fabs(Db) ) {
            /*
             * We found a better value for the center value.
             */
            b  = e;
            Db = De;
            F0 = 0.9; F1 = 0.1;
        } else if ( FirstHalf && ( fabs(De) < fabs(Da) ) && (De < 0.0) ) {
            /*
             * We found a better value for the lower end of the bracket. (And its negative).
             */
            a  = e;
            Da = De;
            F0 = 0.9; F1 = 0.1;
        } else if ( !FirstHalf && ( fabs(De) < fabs(Dc) && (De > 0.0) ) ) {
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

