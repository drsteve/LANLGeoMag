#include "Lgm/Lgm_MagModelInfo.h"
#include "Lgm/Lgm_LstarInfo.h"


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
    Lgm_Vector    u, w, Pm_North;
    double        rat, a, b, c, d, D, I, r, Phi, cl, sl, SS, mlat_min=0.0, Dmin;
    int           done, FoundValidI;


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
    a = mlat0;
    c = mlat1;
    Dmin = 9e99;
    done = FALSE; FoundValidI = FALSE;
    while ( !done ) {


        d = c-a;        // range of bracket
//        b = a + 0.5*d;    // pick new mlat point halfway through range
        b = a + 0.6*d;    // pick new mlat point halfway through range

//printf("a, b, c = %g %g %g\n", a, b, c );

        /*
         *  For this MLT/mlat value, find the radius at which
         *  B = Bm. We assume MLT/mlat is defined relative to
         *  SM coords. All of my tracing routines assume GSM,
         *  so convert to GSM.
         */

//printf("Bm, MLT, b, r = %g %g %g %g\n", Bm, MLT, b, r );
        if ( !FindBmRadius( Bm, MLT, b, &r, LstarInfo->mInfo->Lgm_FindBmRadius_Tol, LstarInfo ) ) {

            /*
             *  Couldnt get a valid Bm. (The bracket is pretty huge,so
             *  we probably ought to believe there really isnt a valid one.)
             */
            printf("%sNo Bm found: setting I to 9e99%s\n", LstarInfo->PreStr, LstarInfo->PostStr);
            I = 9e99;

        } else {

            /*
             *  We found a (candidate) mirror point.
             *  Compute its GSM coords.
             */
            Phi = 15.0*(MLT-12.0)*RadPerDeg;
            //printf("b = %.15lf\n", b);
            cl = cos( b * RadPerDeg ); sl = sin( b * RadPerDeg );
            w.x = r*cl*cos(Phi); w.y = r*cl*sin(Phi); w.z = r*sl;
            Lgm_Convert_Coords( &w, &u, SM_TO_GSM, LstarInfo->mInfo->c );
            if (LstarInfo->VerbosityLevel > 4) {
                printf("%sResults of FindBmRadius: Bm, MLT, mlat, r = %g %g %g %g%s\n", LstarInfo->PreStr, Bm, MLT, b, r, LstarInfo->PostStr);
                printf("%sResults of FindBmRadius: u_sm  = %g %g %g%s\n", LstarInfo->PreStr, w.x, w.y, w.z, LstarInfo->PostStr);
                printf("%sResults of FindBmRadius: u_gsm = %g %g %g%s\n", LstarInfo->PreStr, u.x, u.y, u.z, LstarInfo->PostStr);
            }


            /*
             *  This point is the northern mirror point.
             */
            LstarInfo->mInfo->Pm_North = u;




            /*
             *  Test to see if the Northern Mirror  Point is already close to the Bmin point or the
             *  B's are almost the same; And the Pitch angle is close to 90.  If
             *  so, use an approximation to I.
             */
SS=Lgm_VecDiffMag( &LstarInfo->mInfo->Pm_North, &LstarInfo->mInfo->Pmin );
if (LstarInfo->VerbosityLevel > 1) printf("SS = %g\n", SS);
//            if ( ( ((SS=Lgm_VecDiffMag( &LstarInfo->mInfo->Pm_North, &LstarInfo->mInfo->Pmin )) < 1e-4) || (fabs( LstarInfo->mInfo->Bm - LstarInfo->mInfo->Bmin) < 1e-2) ) && (fabs(90.0-LstarInfo->PitchAngle) < 1e-2)  ) {
            if (  (SS < 1e-2)  && (fabs(90.0-LstarInfo->PitchAngle) < 1e-2)  ) {

                // if FL length is small, use an approx expression for I
                rat = LstarInfo->mInfo->Bmin/LstarInfo->mInfo->Bm;
                if ((1.0-rat) < 0.0) {
                    I = 0.0;
                } else {
                    // Eqn 2.66b in Roederer
                    I = SS*sqrt(1.0 - rat);
                }

            } else {


                /*
                 *  Compute I
                 *
                 *     Compute I for the field line that goes through this point.
                 *     Before we do the integral, we need to locate the mirror points.
                 *     We only have to find the southern mirror point, because we have
                 *     the northern one already.
                 *
                 *     Trace from Pm_North to Pm_South
                 */
                 I = 9e99;
                 //LstarInfo->mInfo->Hmax = 10.0;
                 //LstarInfo->mInfo->Hmax = 0.1;
                 LstarInfo->mInfo->Hmax = 0.1;
                 if ( Lgm_TraceToMirrorPoint( &(LstarInfo->mInfo->Pm_North), &(LstarInfo->mInfo->Pm_South), &SS, LstarInfo->mInfo->Bm, -1.0, LstarInfo->mInfo->Lgm_TraceToMirrorPoint_Tol, LstarInfo->mInfo ) > 0 ) {
                    if ( SS <= 1e-5 ) {
printf("FUCK\n");
                        // if FL length is small, use an approx expression for I
                        rat = LstarInfo->mInfo->Bmin/LstarInfo->mInfo->Bm;
                        if ((1.0-rat) < 0.0) {
                            I = 0.0;
                        } else {
                            // Eqn 2.66b in Roederer
                            I = SS*sqrt(1.0 - rat);
                        }

                    } else {

                        if ( LstarInfo->mInfo->UseInterpRoutines ) {

                            /*
                             *  Do interped I integral. For this to work, we need to trace out the FL with TraceLine().
                             *  This is additional overhead to start with, but it may be faster in the end.
                             *
                             *  Note we start at Pm_South and trace to Pm_North (which is at an altitude of (r-1.0) Re above the Earth.
                             */
                            LstarInfo->mInfo->Hmax = SS/200.0;




                            // Do not include Bmin here (second to last arg must be FALSE). We dont have a proper Bmin here.
                            if ( Lgm_TraceLine2( &(LstarInfo->mInfo->Pm_South), &Pm_North, (r-1.0)*Re, 0.5*SS-LstarInfo->mInfo->Hmax, 1.0, 1e-7, FALSE, LstarInfo->mInfo ) < 0 ) return( -1 );


                            /*
                             *  Set the limits of integration.
                             */
                            LstarInfo->mInfo->Sm_South = 0.0;
                            LstarInfo->mInfo->Sm_North = SS;

                            /*
                             *  Add the mirror points explicity. Update: Actually the
                             *  first should already be there so dont include it.
                             */
                            //AddNewPoint( 0.0, LstarInfo->mInfo->Bm, &LstarInfo->mInfo->Pm_South, LstarInfo->mInfo );
    //MGH MGH                        ReplaceFirstPoint( 0.0, LstarInfo->mInfo->Bm, &LstarInfo->mInfo->Pm_South, LstarInfo->mInfo );
    //MGH MGH                        AddNewPoint( SS,  LstarInfo->mInfo->Bm, &Pm_North, LstarInfo->mInfo );


                            if ( InitSpline( LstarInfo->mInfo ) ) {

                                /*
                                 *  Do I integral with interped integrand.
                                 */
                                I = Iinv_interped( LstarInfo->mInfo  );
                                if (LstarInfo->VerbosityLevel > 1) printf("\t\t%s  Integral Invariant, I (interped):      %15.8g    I-I0:    %15.8g    mlat:   %12.8lf  (nCalls = %d)%s\n",  LstarInfo->PreStr, I, I-I0, b, LstarInfo->mInfo->Lgm_n_I_integrand_Calls, LstarInfo->PostStr );

                                FreeSpline( LstarInfo->mInfo );

                            } else {

                                I = -9e99;

                            }

                        } else {

                            /*
                             *  Set the limits of integration. Also set tolerances for
                             *  Quadpack routines.
                             */
                            LstarInfo->mInfo->Sm_South = 0.0;
                            LstarInfo->mInfo->Sm_North = SS;

                            /*
                             *  Do full blown I integral.
                             */
                            I = Iinv( LstarInfo->mInfo  );
                            if (LstarInfo->VerbosityLevel > 1) printf("\t\t%s  Integral Invariant, I (full integral): %15.8g    I-I0:    %15.8g    mlat:   %12.8lf  (nCalls = %d)%s\n",  LstarInfo->PreStr, I, I-I0, b, LstarInfo->mInfo->Lgm_n_I_integrand_Calls, LstarInfo->PostStr );

                        }
                    }


                } else {


                    /*
                     * open field line
                     */
                    if (LstarInfo->VerbosityLevel > 2) printf("%sOpen Field line: setting I to 9e99%s\n", LstarInfo->PreStr, LstarInfo->PostStr );

                 }

            }

        }













        /*
         *  Compute difference between I and the desired I (i.e. I0 that the caller gave us)
         */
        D = I-I0;

        /*
         *  Hold on to closest value
         */
        if (fabs(D) < Dmin){
            Dmin = fabs(D);
            mlat_min = b;
        }



        if ( fabs( D ) < LstarInfo->mInfo->Lgm_FindShellLine_I_Tol ) {
            /* converged */
            done = TRUE;
            FoundValidI = TRUE;
//            *mlat = 0.5*(a+c);  // Take average of endpoints as final answer
            *mlat = mlat_min;
        } else if ( fabs(b-mlat0)  < 1e-5 ) {
            /* converged  to lower endpoint -- no valid mlat found within interval - probably have
             * to enlarge initial bracket.
             */
            done = TRUE;
            FoundValidI = -2;
        //} else if ( fabs(b-mlat1) < 1e-5 ) {
        } else if ( fabs(b-mlat1) < 1e-8 ) {
            /* converged  to upper endpoint -- no valid mlat found within interval - probably have
             * to enlarge initial bracket.
             */
            done = TRUE;
            FoundValidI = -3;
        } else if ( fabs(a-c) < 1e-8 ) {
            /* converged  to something, but not I?
             * try to enlarge initial bracket.
             */
            done = TRUE;
            if ( fabs(D) < 0.1 ) {
                FoundValidI = TRUE; // lets just take what we get here....
                *mlat = mlat_min;  // Use the best value we got
            } else {
                FoundValidI = -4;
            }
printf("D = %g\n", D);
FoundValidI = -4;
        } else if ( D > 0.0 ) {
            c = b;
        } else {
            a = b;
        }

    }
//printf("FoundValidI = %d\n", FoundValidI);




    *rad    = r;
    *Ifound = I;


    return( FoundValidI );



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

    Lgm_Vector    u, v, Bvec;
    int            done, FoundValidBm;
    double        Phi, D, a0, c0, Da0, Dc0;
    double        a, b, c, d, B, cl, sl, cp, sp, lat, f, g;


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
        c0 += 0.1;
        u.x = c0*f; u.y = c0*g; u.z = c0*sl;
        Lgm_Convert_Coords( &u, &v, SM_TO_GSM, LstarInfo->mInfo->c );
        LstarInfo->mInfo->Bfield( &v, &Bvec, LstarInfo->mInfo );
        B = Lgm_Magnitude( &Bvec );
        Dc0 = B - Bm;

        // Check to see if we went beyond Bm
        if ( Dc0 < 0.0 ) {
            done = TRUE;
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

    if ( ( Da0 < 0.0 ) || ( Dc0 > 0.0 ) ) return( FALSE );




    /*
     *  Find a bracket by stepping out in small increments.
     */
    done = FALSE;
    b    = a0;
    while ( !done ) {

        u.x = b*f; u.y = b*g; u.z = b*sl;
        Lgm_Convert_Coords( &u, &v, SM_TO_GSM, LstarInfo->mInfo->c );
        LstarInfo->mInfo->Bfield( &v, &Bvec, LstarInfo->mInfo );
        B = Lgm_Magnitude( &Bvec );
        D = B - Bm;

        if ( D < 0.0 ) {
            done = TRUE;
            a0 = b - 0.1;
            c0 = b + 0.1;
            break;
        } else {
            b += 0.1;
        }

    }
    if (LstarInfo->VerbosityLevel > 5) {
        printf("%sFindBmRadius: a0, c0 = %g %g%s\n", LstarInfo->PreStr, a0, c0, LstarInfo->PostStr );
    }


    a = a0, c = c0;

    done = FALSE;
    FoundValidBm = FALSE;
    while ( !done ) {

        d = c-a;
        b = a + 0.5*d;


        u.x = b*f; u.y = b*g; u.z = b*sl;
        Lgm_Convert_Coords( &u, &v, SM_TO_GSM, LstarInfo->mInfo->c );
        LstarInfo->mInfo->Bfield( &v, &Bvec, LstarInfo->mInfo );
        B = Lgm_Magnitude( &Bvec );
        D = B - Bm;


        if ( (fabs((b-a0)/a0) < tol) || (fabs((b-c0)/c0) < tol) ) {
            /* converged to an endpoint -> no Bm in interval! */
            done = TRUE;
            FoundValidBm = FALSE;
        } else if ( fabs(D/Bm) < tol ) {
            /* converged */
            done = TRUE;
            FoundValidBm = TRUE;
        } else if ( fabs(d) < tol ) {
            done = TRUE;
            printf("%sFindBmRadius: Warning, BmRadius converged before Bm (might not have found proper Bm)%s\n", LstarInfo->PreStr, LstarInfo->PostStr );
            FoundValidBm = TRUE; // not really true
        } else if ( D > 0.0 ){
            a = b;
        } else {
            c = b;
        }

        if (LstarInfo->VerbosityLevel > 5) {
            printf("%sFindBmRadius: a, b, c = %g %g %g    D, Bm, fabs(D/Bm) = %g %g %g%s\n", LstarInfo->PreStr, a, b, c, D, Bm, fabs(D/Bm), LstarInfo->PostStr );
        }

    }



    /*
     *  Take average of endpoints as final answer
     */
    *r = 0.5*(a+c);
    if (LstarInfo->VerbosityLevel > 5) {
        printf("%sFindBmRadius: Final r = %.15lf%s\n", LstarInfo->PreStr, *r, LstarInfo->PostStr );
    }
    return( FoundValidBm );

}

