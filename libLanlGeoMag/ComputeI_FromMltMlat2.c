#include "Lgm/Lgm_MagModelInfo.h"
#include "Lgm/Lgm_LstarInfo.h"

#define   TRACE_TOL  1e-7

/*
 * This version traces FLs from the Earth at the given MLT/mlat instead of trying to find the Bm radially out first.
 */
double ComputeI_FromMltMlat2( double Bm, double MLT, double mlat, double *r, double I0, int *ErrorStatus, Lgm_LstarInfo *LstarInfo ) {

    int         reset=1, reset2, TraceFlag;

    double      Bmin, I1, Phi, cl, sl, rat, SS1, SS2, SS, Sn, Ss, Htry, Hdid, Hnext, Bs, Be, s, sgn;
    Lgm_Vector  w, u, Pmirror1, Pmirror2, v1, v2, v3, Bvec, P, Ps, u_scale, Bvectmp, Ptmp;
    double      stmp, Btmp;


    // Assume its good to start
    *ErrorStatus = 1;



    /*
     *  First do a trace to identify the FL type and some of its critical points.
     */
    *r = 1.0 + 100.0/Re;
    Phi = 15.0*(MLT-12.0)*RadPerDeg;
    cl = cos( mlat * RadPerDeg ); sl = sin( mlat * RadPerDeg );
    w.x = (*r)*cl*cos(Phi); w.y = (*r)*cl*sin(Phi); w.z = (*r)*sl;

    Lgm_Convert_Coords( &w, &u, SM_TO_GSM, LstarInfo->mInfo->c );
    TraceFlag = Lgm_Trace( &u, &v1, &v2, &v3, LstarInfo->mInfo->Lgm_LossConeHeight, TRACE_TOL, TRACE_TOL, LstarInfo->mInfo );
    LstarInfo->mInfo->Pmin = v3;

    LstarInfo->mInfo->Bfield( &v3, &Bvec, LstarInfo->mInfo );
    Bmin = Lgm_Magnitude( &Bvec );
    //printf("Pmin = %g %g %g\n", v3.x, v3.y, v3.z);

    if ( TraceFlag != LGM_CLOSED ) {

        I1 = 9e99;
        if (LstarInfo->VerbosityLevel > 1) {
            printf("\t\t%s  mlat: %13.6g   I1: %13.6g   I0: %13.6g   > Field Line not closed%s\n",  LstarInfo->PreStr, mlat, I1, I0, LstarInfo->PostStr );
        }

        *ErrorStatus = -1; // FL open. I undefined.
        return( I1 );

    } else if ( Bmin <= Bm ) {

        /*
         * From the minimum B point, attempt to trace along the field to get the northern mirror point.
         */
        P = v3;
        SS1 = 0.0;
        if ( Lgm_TraceToMirrorPoint( &P, &Pmirror1, &SS1, LstarInfo->mInfo->Bm,  1.0, LstarInfo->mInfo->Lgm_TraceToMirrorPoint_Tol, LstarInfo->mInfo ) > 0 )  {

            LstarInfo->mInfo->Pm_North = Pmirror1;
            //printf("Pmirror1 = %g %g %g   Bm = %g    LstarInfo->mInfo->Bm = %g\n", Pmirror1.x, Pmirror1.y, Pmirror1.z, Bm, LstarInfo->mInfo->Bm);


            P = v3;
            SS2 = 0.0;
            if ( Lgm_TraceToMirrorPoint( &P, &Pmirror2, &SS2, LstarInfo->mInfo->Bm,  -1.0, LstarInfo->mInfo->Lgm_TraceToMirrorPoint_Tol, LstarInfo->mInfo ) > 0 )  {

                LstarInfo->mInfo->Pm_South = Pmirror2;
                //printf("Pmirror2 = %g %g %g\n", Pmirror2.x, Pmirror2.y, Pmirror2.z);

            } else {
                I1 = 9e99;
                if (LstarInfo->VerbosityLevel > 1) {
                    printf("\t\t%s  mlat: %13.6g   I1: %13.6g   I0: %13.6g   > Unable to find southern mirror point%s\n",  LstarInfo->PreStr, mlat, I1, I0, LstarInfo->PostStr );
                }
                *ErrorStatus = -3; // No valid mirror point in the south
                return( I1 );
            }

            // total distance between mirror point.
            SS = SS1 + SS2;
            
            // If its really small, just return 0.0 for I
            if ( fabs(SS) < 1e-7) {
                I1 = 0.0;
                if (LstarInfo->VerbosityLevel > 1) {
                    printf("\t\t%s  mlat: %13.6g   I1: %13.6g   I0: %13.6g   > Distance between mirror points is < 1e-7Re, Assuming I1=0%s\n",  LstarInfo->PreStr, mlat, I1, I0, LstarInfo->PostStr );
                }
                *ErrorStatus = 2; // Flag that we did this
                return( I1 );
            }

        } else {
            I1 = 9e99;
            if (LstarInfo->VerbosityLevel > 1) {
                printf("\t\t%s  mlat: %13.6g   I1: %13.6g   I0: %13.6g   > Unable to find northern mirror point%s\n",  LstarInfo->PreStr, mlat, I1, I0, LstarInfo->PostStr );
            }
            *ErrorStatus = -2; // No valid mirror point in the north
            return( I1 );
        }

        /*
         * OK, we have both mirror points. Lets compute I
         */
        I1 = 9e99;
        LstarInfo->mInfo->Hmax = 0.1;
        LstarInfo->mInfo->Hmax = SS/(double)LstarInfo->mInfo->nDivs;

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

         
        if ( Lgm_TraceLine3( &(LstarInfo->mInfo->Pm_South), SS, nDivs, 1.0, 1e-7, FALSE, LstarInfo->mInfo ) < 0 ) return( 9e99 );

        /*
         *  Set the limits of integration.
         */
        LstarInfo->mInfo->Sm_South = 0.0;
        LstarInfo->mInfo->Sm_North = SS;

        if ( InitSpline( LstarInfo->mInfo ) ) {

            /*
             *  Do I integral with interped integrand.
             */
             I1 = Iinv_interped( LstarInfo->mInfo );
             if (LstarInfo->VerbosityLevel > 1) {
                printf("\t\t%s  mlat: %13.6g   I1: %13.6g   I0: %13.6g   I1-I0: %13.6g    [Sa,Sb]: %.8g  %.8g  (nCalls = %d)%s\n",  LstarInfo->PreStr, mlat, I1, I0, I1-I0, LstarInfo->mInfo->Sm_South, LstarInfo->mInfo->Sm_North, LstarInfo->mInfo->Lgm_n_I_integrand_Calls, LstarInfo->PostStr );
             }
             FreeSpline( LstarInfo->mInfo );

        } else {
            I1 = 9e99;
            if (LstarInfo->VerbosityLevel > 1) {
                printf("\t\t%s  mlat: %13.6g   I1: %13.6g   I0: %13.6g   Couldnt initialize spline%s\n",  LstarInfo->PreStr, mlat, I1, I0, LstarInfo->PostStr );
            }
        }


    } else {

//This seems problematic.
//I see the brackets being affected by this, but often in the wrong mdirection!

        if (LstarInfo->VerbosityLevel > 1){ printf( "\t\t\t> Field line min-B is greater than Bm. mlat = %g MLT = %g  Bm = %g Bmin = %g (u = %g %g %g   v1 = %g %g %g  v2 = %g %g %g  v3 = %g %g %g)\n", mlat, MLT, Bm, Bmin, u.x, u.y, u.z, v1.x, v1.y, v1.z, v2.x, v2.y, v2.z, v3.x, v3.y, v3.z ); }

        *ErrorStatus = -10; // We didnt even try to compute I because we Field line min-B is greater than Bm.
        if (LstarInfo->VerbosityLevel > 1) {
            printf("\t\t%s  mlat: %13.6g   I1: undefined   I0: %13.6g   I1 undefined. min-B is greater than Bm.%s\n",  LstarInfo->PreStr, mlat, I0, LstarInfo->PostStr );
        }
        return( 9e99 );

    }


    return( I1 );
    

}
