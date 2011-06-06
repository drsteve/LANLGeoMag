#include "Lgm/Lgm_MagModelInfo.h"
#include "Lgm/Lgm_LstarInfo.h"

double ComputeI_FromMltMlat( double Bm, double MLT, double mlat, double *r, double I0, Lgm_LstarInfo *LstarInfo ) {

    int         reset=1;
    
    double      I, Phi, cl, sl, rat, SS1, SS2, SS, Sn, Ss, Htry, Hdid, Hnext, Bs, Be, s, sgn;
    Lgm_Vector  w, u, Pmirror1, Pmirror2, v1, v2, v3, Bvec, P, Ps, u_scale;


    /*
     *  For this MLT/mlat value, find the radius at which
     *  B = Bm. We assume MLT/mlat is defined relative to
     *  SM coords. All of my tracing routines assume GSM,
     *  so convert to GSM.
     */
    //printf("Bm, MLT, b, r = %g %g %g %g\n", Bm, MLT, b, r );
//printf("LstarInfo->mInfo->Lgm_FindBmRadius_Tol = %g\n", LstarInfo->mInfo->Lgm_FindBmRadius_Tol);
    if ( !FindBmRadius( Bm, MLT, mlat, r, LstarInfo->mInfo->Lgm_FindBmRadius_Tol, LstarInfo ) ) {

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
        cl = cos( mlat * RadPerDeg ); sl = sin( mlat * RadPerDeg );
        w.x = (*r)*cl*cos(Phi); w.y = (*r)*cl*sin(Phi); w.z = (*r)*sl;
        Lgm_Convert_Coords( &w, &u, SM_TO_GSM, LstarInfo->mInfo->c );
        if (LstarInfo->VerbosityLevel > 4) {
            printf("%sResults of FindBmRadius: Bm, MLT, mlat, r = %g %g %g %g%s\n", LstarInfo->PreStr, Bm, MLT, mlat, (*r), LstarInfo->PostStr);
            printf("%sResults of FindBmRadius: u_sm  = %g %g %g%s\n", LstarInfo->PreStr, w.x, w.y, w.z, LstarInfo->PostStr);
            printf("%sResults of FindBmRadius: u_gsm = %g %g %g%s\n", LstarInfo->PreStr, u.x, u.y, u.z, LstarInfo->PostStr);
        }


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
        Htry = 1e-3; // we probably dont ever need to split the mirror points to any finer precision than this(?).
        u_scale.x =  100.0;  u_scale.y = 100.0; u_scale.z = 100.0;
        P = Pmirror1;
        LstarInfo->mInfo->Bfield( &P, &Bvec, LstarInfo->mInfo );
        Bs = Lgm_Magnitude( &Bvec );
//printf("Bs-Bm = %g\n", Bs-Bm);

        if ( Lgm_MagStep( &P, &u_scale, Htry, &Hdid, &Hnext, 1.0e-7, -1.0, &s, &reset, LstarInfo->mInfo->Bfield, LstarInfo->mInfo ) < 0 ) return( LGM_BAD_TRACE );

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

            if ( Lgm_MagStep( &P, &u_scale, Htry, &Hdid, &Hnext, 1.0e-7, 1.0, &s, &reset, LstarInfo->mInfo->Bfield, LstarInfo->mInfo ) < 0 ) return( LGM_BAD_TRACE );

            LstarInfo->mInfo->Bfield( &P, &Bvec, LstarInfo->mInfo );
            Be  = Lgm_Magnitude( &Bvec );

//printf("2b. Be-Bs = %g        Be, Bs = %.15g %.15g Bm = %.15g\n", Be - Bs, Be, Bs, LstarInfo->mInfo->Bm);
            if ( Be < Bs ) {
                sgn = 1.0;
            } else {
                // we are probably very close to Pmin. So I=0.
                return( 0.0 ); 
            }

        }

        SS1 = Hdid;
//SS1 = 0.0;



//        if ( Lgm_TraceToMirrorPoint( &Pmirror1, &Pmirror2, &SS, LstarInfo->mInfo->Bm,  sgn, LstarInfo->mInfo->Lgm_TraceToMirrorPoint_Tol, LstarInfo->mInfo ) > 0 ) {
        SS2 = 0.0;
        if ( Lgm_TraceToMirrorPoint( &P, &Pmirror2, &SS2, LstarInfo->mInfo->Bm,  sgn, LstarInfo->mInfo->Lgm_TraceToMirrorPoint_Tol, LstarInfo->mInfo ) > 0 ) {

            SS = SS1 + SS2;
//printf("SS = %g\n", SS);

            if ( sgn < 0.0 ) {
                LstarInfo->mInfo->Pm_North = Pmirror1;
                LstarInfo->mInfo->Pm_South = Pmirror2;
            } else {
                LstarInfo->mInfo->Pm_North = Pmirror2;
                LstarInfo->mInfo->Pm_South = Pmirror1;
            }
            SS = fabs( SS );
            if (SS < 1e-7) {
//printf("AHA! SS = %.15g\n", SS);
//printf("Pm_South = %.15g %.15g %.15g\n", LstarInfo->mInfo->Pm_South.x, LstarInfo->mInfo->Pm_South.y, LstarInfo->mInfo->Pm_South.z );
//printf("Pm_North = %.15g %.15g %.15g\n", LstarInfo->mInfo->Pm_North.x, LstarInfo->mInfo->Pm_North.y, LstarInfo->mInfo->Pm_North.z );
                return(0.0);
            }


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
                LstarInfo->mInfo->Hmax = SS/200.0;
//printf("LstarInfo->mInfo->Hmax = %g\n", LstarInfo->mInfo->Hmax);



                // Do not include Bmin here (second to last arg must be FALSE). We dont have a proper Bmin here.
                //if ( Lgm_TraceLine2( &(LstarInfo->mInfo->Pm_South), &(LstarInfo->mInfo->Pm_North), (r-1.0)*Re, 0.5*SS-LstarInfo->mInfo->Hmax, 1.0, 1e-7, FALSE, LstarInfo->mInfo ) < 0 ) return( -9e99 );
                //if ( Lgm_TraceLine2( &(LstarInfo->mInfo->Pm_South), &(LstarInfo->mInfo->Pm_North), (r-1.0)*Re, 0.5*SS-LstarInfo->mInfo->Hmax, 1.0, 1e-7, FALSE, LstarInfo->mInfo ) < 0 ) return( -9e99 );
                if ( Lgm_TraceLine3( &(LstarInfo->mInfo->Pm_South), SS, 200, 1.0, 1e-7, FALSE, LstarInfo->mInfo ) < 0 ) return( 9e99 );


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
//MGH MGH                            ReplaceFirstPoint( 0.0, LstarInfo->mInfo->Bm, &LstarInfo->mInfo->Pm_South, LstarInfo->mInfo );
//MGH MGH                            AddNewPoint( SS,  LstarInfo->mInfo->Bm, &Pm_North, LstarInfo->mInfo );


                if ( InitSpline( LstarInfo->mInfo ) ) {

                    /*
                     *  Do I integral with interped integrand.
                     */
                    I = Iinv_interped( LstarInfo->mInfo  );
                    if (LstarInfo->VerbosityLevel > 1) printf("\t\t%s  Integral Invariant, I (interped):      %15.8g    I-I0:    %15.8g    mlat:   %12.8lf  (nCalls = %d)%s\n",  LstarInfo->PreStr, I, I-I0, mlat, LstarInfo->mInfo->Lgm_n_I_integrand_Calls, LstarInfo->PostStr );
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

        }

    }




    return( I );

}
