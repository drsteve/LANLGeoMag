#include "Lgm/Lgm_MagModelInfo.h"
#include "Lgm/Lgm_LstarInfo.h"

/*
 *  This is a wrapper function to switch between ways of computing I when doing the L* calcs.
 *
 *  1. The original method (ROEDERER's) seacrhes for Bm along a ray launched from the gievn MLT/MLAT.
 *
 *  2. The new method seacrhes for Bm along a FL launched traced from the surface at the given MLT/MLAT.
 *
 */
double ComputeI_FromMltMlat( double Bm, double MLT, double mlat, double *r, double I0, Lgm_LstarInfo *LstarInfo ) {

    if ( LstarInfo->ISearchMethod == 1 ) {

        return( ComputeI_FromMltMlat1( Bm, MLT, mlat, r, I0, LstarInfo ) );

    } else if ( LstarInfo->ISearchMethod == 2 ) {

        return( ComputeI_FromMltMlat2( Bm, MLT, mlat, r, I0, LstarInfo ) );

    } else {

        printf("Unkown valid for ISearchMethod (must be 1 or 2)\n");
        exit(0);

    }

}



double ComputeI_FromMltMlat1( double Bm, double MLT, double mlat, double *r, double I0, Lgm_LstarInfo *LstarInfo ) {

    int         reset=1, reset2;

    double      I, Phi, cl, sl, rat, SS1, SS2, SS, Sn, Ss, Htry, Hdid, Hnext, Bs, Be, s, sgn;
    Lgm_Vector  w, u, Pmirror1, Pmirror2, v1, v2, v3, Bvec, P, Ps, u_scale, Bvectmp, Ptmp;
    double      stmp, Btmp;


    /*
     *  For this MLT/mlat value, find the radius at which
     *  B = Bm. We assume MLT/mlat is defined relative to
     *  SM coords. All of my tracing routines assume GSM,
     *  so convert to GSM.
     */
//printf("LstarInfo->mInfo->Lgm_FindBmRadius_Tol = %g\n", LstarInfo->mInfo->Lgm_FindBmRadius_Tol);
    if ( !FindBmRadius( Bm, MLT, mlat, r, LstarInfo->mInfo->Lgm_FindBmRadius_Tol, LstarInfo ) ) {

        /*
         *  Couldnt get a valid Bm. (The bracket is pretty huge,so
         *  we probably ought to believe there really isnt a valid one.)
         */
        if (LstarInfo->VerbosityLevel > 1) printf("\t%sNo Bm found: setting I to 9e99%s\n", LstarInfo->PreStr, LstarInfo->PostStr);
        I = 9e99;

    } else {

        /*
         *  We found a (candidate) mirror point.
         *  Compute its GSM coords.
         */
        Phi = 15.0*(MLT-12.0)*RadPerDeg;
        //printf("b = %.15lf\n", b);
        cl = cos( mlat * RadPerDeg ); sl = sin( mlat * RadPerDeg );
        w.x = (*r)*cl*cos(Phi); w.y = (*r)*cl*sin(Phi); w.z = (*r)*sl;
        Lgm_Convert_Coords( &w, &u, SM_TO_GSM, LstarInfo->mInfo->c );
        if (LstarInfo->VerbosityLevel > 4) {
            printf("\t%sResults of FindBmRadius: Bm, MLT, mlat, r = %g %g %g %g%s\n", LstarInfo->PreStr, Bm, MLT, mlat, (*r), LstarInfo->PostStr);
            printf("\t%sResults of FindBmRadius: u_sm  = %g %g %g%s\n", LstarInfo->PreStr, w.x, w.y, w.z, LstarInfo->PostStr);
            printf("\t%sResults of FindBmRadius: u_gsm = %g %g %g%s\n", LstarInfo->PreStr, u.x, u.y, u.z, LstarInfo->PostStr);
        }

        /*
         * Check to see if FL is opne/closed/etc.
         */
        Lgm_Vector Foot_n, Foot_s;
        if ( !Lgm_TraceToSphericalEarth( &u, &Foot_n, 120.0, -1.0, 0.01, LstarInfo->mInfo ) ) return( -1 );
        if ( !Lgm_TraceToSphericalEarth( &u, &Foot_s, 120.0,  1.0, 0.01, LstarInfo->mInfo ) ) return( -1 );


        Pmirror1 = u;
        /*
         *  The point Pmirror1 is a candidate northern mirror point. Note that
         *  there is no guarantee that it actually *is* the northern mirror
         *  point. If we are close to the Bmin surface, then it could be the
         *  southern mirror point (mix-ups could occur due to
         *  tolerance/roundoff issues.).
         *
         *  If LstarInfo->mInfo->Lgm_FindBmRadius_Tol was set to a low enough
         *  value, we will have this mirror point determined to a very high
         *  precision, so |B-Bm| at this point should be lower than we really
         *  need.  However, while it should be close to zero, it could be
         *  positive or negative.
         *
         *  To find the other mirror point, we need to know which way Bmin is.
         *  Lgm_TraceToMirrorPoint() should do the rest is for us.
         */
        Htry = 1e-5; // (used to be 1e-6) we probably dont ever need to split the mirror points to any finer precision than this(?).
        //u_scale.x =  100.0;  u_scale.y = 100.0; u_scale.z = 100.0;
        u_scale.x =  1.0;  u_scale.y = 1.0; u_scale.z = 1.0;
        P = Pmirror1;
        //printf("P = %g %g %g\n", P.x, P.y, P.z);
        LstarInfo->mInfo->Bfield( &P, &Bvec, LstarInfo->mInfo );
        Bs = Lgm_Magnitude( &Bvec );
        //printf("Ps = %g %g %g\n", P.x, P.y, P.z);
//printf("Bs-Bm = %g\n", Bs-Bm);

        if ( Lgm_MagStep( &P, &u_scale, Htry, &Hdid, &Hnext, -1.0, &s, &reset, LstarInfo->mInfo->Bfield, LstarInfo->mInfo ) < 0 ) return( -1.0 );
//printf("Hdid = %g\n", Hdid);

        LstarInfo->mInfo->Bfield( &P, &Bvec, LstarInfo->mInfo );
        Be  = Lgm_Magnitude( &Bvec );

        if ( Be < Bs ) {
            // our assumption that its a northern mirror point is probably correct.
//printf("1. Be-Bs = %g        Be, Bs = %.15g %.15g  Bm = %.15g\n", Be - Bs, Be, Bs, LstarInfo->mInfo->Bm);
            sgn = -1.0;
        } else {

//printf("2a. Be-Bs = %g        Be, Bs = %.15g %.15g Bm = %.15g\n", Be - Bs, Be, Bs, LstarInfo->mInfo->Bm);
            // try the other direction
            P = Pmirror1;
            LstarInfo->mInfo->Bfield( &P, &Bvec, LstarInfo->mInfo );
            Bs = Lgm_Magnitude( &Bvec );

            if ( Lgm_MagStep( &P, &u_scale, Htry, &Hdid, &Hnext, 1.0, &s, &reset, LstarInfo->mInfo->Bfield, LstarInfo->mInfo ) < 0 ) return( -1.0 );
//printf("Hdid = %g\n", Hdid);

            LstarInfo->mInfo->Bfield( &P, &Bvec, LstarInfo->mInfo );
            Be  = Lgm_Magnitude( &Bvec );

//printf("2b. Be-Bs = %g        Be, Bs = %.15g %.15g Bm = %.15g\n", Be - Bs, Be, Bs, LstarInfo->mInfo->Bm);
            if ( Be < Bs ) {
                sgn = 1.0;
            } else {
                // we are probably very close to Pmin. So I=0.
                LstarInfo->mInfo->Pm_North = Pmirror1;
                LstarInfo->mInfo->Pm_South = Pmirror1;
                return( 0.0 );
            }

        }

        SS1 = Hdid;
//printf("Hdid = %g\n", Hdid);
        //SS1 = 0.0;



        SS2 = 0.0;
//        if ( Lgm_TraceToMirrorPoint( &Pmirror1, &Pmirror2, &SS2, LstarInfo->mInfo->Bm,  sgn, LstarInfo->mInfo->Lgm_TraceToMirrorPoint_Tol, LstarInfo->mInfo ) > 0 ) {
        if ( Lgm_TraceToMirrorPoint( &P, &Pmirror2, &SS2, LstarInfo->mInfo->Bm,  sgn, LstarInfo->mInfo->Lgm_TraceToMirrorPoint_Tol, LstarInfo->mInfo ) > 0 )  {
//printf("LstarInfo->mInfo->Lgm_TraceToMirrorPoint_Tol = %g\n", LstarInfo->mInfo->Lgm_TraceToMirrorPoint_Tol);

            SS = SS1 + SS2;
//printf("SS1, SS2, SS = %g %g %g   sgn = %g\n", SS1, SS2, SS, sgn);
//printf("%g %g\n", mlat, SS);

            if ( sgn < 0.0 ) {
                LstarInfo->mInfo->Pm_North = Pmirror1;
                LstarInfo->mInfo->Pm_South = Pmirror2;
            } else {
                LstarInfo->mInfo->Pm_North = Pmirror2;
                LstarInfo->mInfo->Pm_South = Pmirror1;
//printf("GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGg\n");
            }
            SS = fabs( SS );
            if (SS < 1e-7) {
//printf("AHA! SS = %.15g\n", SS);
//printf("Pm_South = %.15g %.15g %.15g\n", LstarInfo->mInfo->Pm_South.x, LstarInfo->mInfo->Pm_South.y, LstarInfo->mInfo->Pm_South.z );
//printf("Pm_North = %.15g %.15g %.15g\n", LstarInfo->mInfo->Pm_North.x, LstarInfo->mInfo->Pm_North.y, LstarInfo->mInfo->Pm_North.z );
                return(0.0);
            }
//printf("Pm_South = %.15g %.15g %.15g\n", LstarInfo->mInfo->Pm_South.x, LstarInfo->mInfo->Pm_South.y, LstarInfo->mInfo->Pm_South.z );
//printf("Pm_North = %.15g %.15g %.15g\n", LstarInfo->mInfo->Pm_North.x, LstarInfo->mInfo->Pm_North.y, LstarInfo->mInfo->Pm_North.z );

//Lgm_Vector vvv;
//Lgm_Convert_Coords( &LstarInfo->mInfo->Pm_South, &vvv, GSM_TO_SM, LstarInfo->mInfo->c );
//LstarInfo->mInfo->Bfield( &LstarInfo->mInfo->Pm_South, &Bvec, LstarInfo->mInfo );
//printf("Pm_South_sm = %.15g %.15g %.15g      B-Bm = %.15g\n", vvv.x, vvv.y, vvv.z, Lgm_Magnitude( &Bvec)-LstarInfo->mInfo->Bm );
//LstarInfo->mInfo->Bfield( &LstarInfo->mInfo->Pm_North, &Bvec, LstarInfo->mInfo );
//Lgm_Convert_Coords( &LstarInfo->mInfo->Pm_North, &vvv, GSM_TO_SM, LstarInfo->mInfo->c );
//printf("Pm_North_sm = %.15g %.15g %.15g      B-Bm = %.15g\n", vvv.x, vvv.y, vvv.z, Lgm_Magnitude( &Bvec)-LstarInfo->mInfo->Bm );
//printf("%.15g %.15g\n", mlat, LstarInfo->mInfo->Pm_North.z);

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

//            printf("0. HERE: Pmin, Pm_North, Pm_South = %g %g %g   %g %g %g   %g %g %g   SS = %g\n",   LstarInfo->mInfo->Pmin.x, LstarInfo->mInfo->Pmin.y, LstarInfo->mInfo->Pmin.z,
//                                                                                                       LstarInfo->mInfo->Pm_North.x,  LstarInfo->mInfo->Pm_North.y,  LstarInfo->mInfo->Pm_North.z,
//                                                                                                       LstarInfo->mInfo->Pm_South.x,  LstarInfo->mInfo->Pm_South.y,  LstarInfo->mInfo->Pm_South.z, SS );

            if ( LstarInfo->mInfo->UseInterpRoutines ) {

                /*
                 *  Do interped I integral. For this to work, we need to trace out the FL with TraceLine().
                 *  This is additional overhead to start with, but it may be faster in the end.
                 *
                 *  Note we start at Pm_South and trace to Pm_North (which is at an altitude of (r-1.0) Re above the Earth.
                 */
                //LstarInfo->mInfo->Hmax = SS/200.0;
                LstarInfo->mInfo->Hmax = SS/(double)LstarInfo->mInfo->nDivs;
//printf("LstarInfo->mInfo->Hmax = %g\n", LstarInfo->mInfo->Hmax);



                // Do not include Bmin here (second to last arg must be FALSE). We dont have a proper Bmin here.
                //if ( Lgm_TraceLine2( &(LstarInfo->mInfo->Pm_South), &(LstarInfo->mInfo->Pm_North), (*r-1.0)*Re, 0.5*SS-LstarInfo->mInfo->Hmax, 1.0, 1e-7, FALSE, LstarInfo->mInfo ) < 0 ) return( LGM_FILL_VALUE );
                //if ( Lgm_TraceLine3( &(LstarInfo->mInfo->Pm_South), SS, 200, 1.0, 1e-7, FALSE, LstarInfo->mInfo ) < 0 ) return( 9e99 );
//printf("MIKE: SS = %g\n", SS);
//printf("MIKE: SS/LstarInfo->mInfo->nDivs = %g\n", SS/LstarInfo->mInfo->nDivs);


                /*
                 *  This little section is attempting to solve an annoying
                 *  precision issue.  If the distance over which we are trying
                 *  to trace is too small, too many sub-divisions (i.e. total
                 *  steps) will lead to a step size that is too small.  We can
                 *  end up with the gridded FL points computed in
                 *  Lgm_TraceLine3() such that the distance, s of the final
                 *  point is very slightly less than we expect.  This tries to
                 *  reduce the number of divisions to avoid this in cases where
                 *  the total distance to trace is very small.
                 */
                int nDivs;
                if ( SS/LstarInfo->mInfo->nDivs < 1e-6 ) {
                    nDivs = SS/1e-6;
                    if (nDivs < 10) nDivs = 10;
                } else {
                    nDivs = LstarInfo->mInfo->nDivs;
                }
                //printf("nDivs = %d\n", nDivs);



                //if ( Lgm_TraceLine3( &(LstarInfo->mInfo->Pm_South), SS, LstarInfo->mInfo->nDivs, 1.0, 1e-7, FALSE, LstarInfo->mInfo ) < 0 ) return( 9e99 );
                if ( Lgm_TraceLine3( &(LstarInfo->mInfo->Pm_South), SS, nDivs, 1.0, 1e-7, FALSE, LstarInfo->mInfo ) < 0 ) return( 9e99 );
//printf("LstarInfo->mInfo->Pm_South = %g %g %g\n", LstarInfo->mInfo->Pm_South.x, LstarInfo->mInfo->Pm_South.y, LstarInfo->mInfo->Pm_South.z );
//printf("P0 = %g %g %g\n", LstarInfo->mInfo->Px[0], LstarInfo->mInfo->Py[0], LstarInfo->mInfo->Pz[0]);
//printf("Plast = %g %g %g\n", LstarInfo->mInfo->Px[LstarInfo->mInfo->nPnts-1], LstarInfo->mInfo->Py[LstarInfo->mInfo->nPnts-1], LstarInfo->mInfo->Pz[LstarInfo->mInfo->nPnts-1]);






                /*
                 *  Set the limits of integration.
                 */
                LstarInfo->mInfo->Sm_South = 0.0;
                LstarInfo->mInfo->Sm_North = SS;
//printf("LstarInfo->mInfo->Sm_South, LstarInfo->mInfo->Sm_North = %g %g\n", LstarInfo->mInfo->Sm_South, LstarInfo->mInfo->Sm_North);

                /*
                 *  Add the mirror points explicity. Update: Actually the
                 *  first should already be there so dont include it.
                 */
                //AddNewPoint( 0.0, LstarInfo->mInfo->Bm, &LstarInfo->mInfo->Pm_South, LstarInfo->mInfo );
//MGH MGH                            ReplaceFirstPoint( 0.0, LstarInfo->mInfo->Bm, &LstarInfo->mInfo->Pm_South, LstarInfo->mInfo );
//MGH MGH                            AddNewPoint( SS,  LstarInfo->mInfo->Bm, &Pm_North, LstarInfo->mInfo );

if (0==1){
                ReplaceFirstPoint( 0.0, LstarInfo->mInfo->Bm, &LstarInfo->mInfo->Pm_South, LstarInfo->mInfo );
                ReplaceLastPoint( SS, LstarInfo->mInfo->Bm, &LstarInfo->mInfo->Pm_North, LstarInfo->mInfo );

                /*
                 * Make sure we have a small margin before and after so
                 * we dont end up trying to extrapolate if s ever gets
                 * slightly out of bounds.
                 */
                 Ptmp = LstarInfo->mInfo->Pm_South; stmp = 0.0; reset2 = FALSE;
                 if ( Lgm_MagStep( &Ptmp, &u_scale, 0.01, &Hdid, &Hnext, -1.0, &stmp, &reset2, LstarInfo->mInfo->Bfield, LstarInfo->mInfo ) < 0 ) { return(-1); }
                 LstarInfo->mInfo->Bfield( &Ptmp, &Bvectmp, LstarInfo->mInfo ); Btmp = Lgm_Magnitude( &Bvectmp );
                 //printf("-stmp, Btmp = %g %g\n", -stmp, Btmp );
                 AddNewPoint( -stmp, Btmp, &Ptmp, LstarInfo->mInfo );

                 Ptmp = LstarInfo->mInfo->Pm_North; stmp = 0.0; reset2 = FALSE;
                 if ( Lgm_MagStep( &Ptmp, &u_scale, 0.01, &Hdid, &Hnext, 1.0, &stmp, &reset2, LstarInfo->mInfo->Bfield, LstarInfo->mInfo ) < 0 ) { return(-1); }
                 LstarInfo->mInfo->Bfield( &Ptmp, &Bvectmp, LstarInfo->mInfo ); Btmp = Lgm_Magnitude( &Bvectmp );
                 //printf("stmp, Btmp = %g %g\n", SS+stmp, Btmp );
                 AddNewPoint( SS+stmp, Btmp, &Ptmp, LstarInfo->mInfo );

//int iiii;
//for (iiii=0; iiii<LstarInfo->mInfo->nPnts; iiii++){
//printf("%g %g\n", LstarInfo->mInfo->s[iiii], LstarInfo->mInfo->Bmag[iiii]);
//}
}


                if ( InitSpline( LstarInfo->mInfo ) ) {

                    /*
                     *  Do I integral with interped integrand.
                     */
//printf("I = %g\n", I);
                    I = Iinv_interped( LstarInfo->mInfo );
//printf("I = %g     Sm_South, Sm_North = %g %g\n", I, LstarInfo->mInfo->Sm_South, LstarInfo->mInfo->Sm_North);
//                    if (LstarInfo->VerbosityLevel > 1) printf("\t\t%s  Integral Invariant, I (interped):      %15.8g    I-I0:    %15.8g    [a,b]: %.15g  %.15g  mlat:   %12.8lf  (nCalls = %d)%s\n",  LstarInfo->PreStr, I, I-I0, LstarInfo->mInfo->Sm_South, LstarInfo->mInfo->Sm_North, mlat, LstarInfo->mInfo->Lgm_n_I_integrand_Calls, LstarInfo->PostStr );
                    if (LstarInfo->VerbosityLevel > 1) {
                        printf("\t\t%s  mlat: %13.6g   I: %13.6g   I0: %13.6g   I-I0: %13.6g    [Sa,Sb]: %.8g  %.8g  (nCalls = %d)%s\n",  LstarInfo->PreStr, mlat, I, I0, I-I0, LstarInfo->mInfo->Sm_South, LstarInfo->mInfo->Sm_North, LstarInfo->mInfo->Lgm_n_I_integrand_Calls, LstarInfo->PostStr );
                    }
                    FreeSpline( LstarInfo->mInfo );

                } else {

                    I = 9e99;

                }

            } else {

                /*
                 *  Set the limits of integration.
                 */
                LstarInfo->mInfo->Sm_South = 0.0;
                LstarInfo->mInfo->Sm_North = SS;

                /*
                 *  Do full blown I integral.
                 */
                I = Iinv( LstarInfo->mInfo  );
                if (LstarInfo->VerbosityLevel > 1) printf("\t\t%s  Integral Invariant, I (full integral): %15.8g    I-I0:    %15.8g    mlat:   %12.8lf  (nCalls = %d)%s\n",  LstarInfo->PreStr, I, I-I0, mlat, LstarInfo->mInfo->Lgm_n_I_integrand_Calls, LstarInfo->PostStr );
            }

        } else {
            I = 9e99;
        }

    }




    return( I );

}
