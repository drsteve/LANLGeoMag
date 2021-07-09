#include "Lgm/Lgm_MagModelInfo.h"
//#include "MagStep.h"

#define JUMP_METHOD	    0

#define USE_SIX_POINT   0
#define USE_FOUR_POINT  1
#define USE_TWO_POINT   2

#define DIFF_SCHEME     USE_TWO_POINT
//#define DIFF_SCHEME     USE_FOUR_POINT

/*
 *  Computes the following integral
 *
 *
 *                          / north_footpoint       
 *                         |
 *                         |	               [	 1    ]
 *                 V  =    |                   [ -------- ]  	ds
 *                         |                   [    B(s)  ]
 *                         |
 *                        / south_footpoint
 *                       
 *
 *
 */





/**
 *   This routine evaluates the "flux tube volume, V" from footpoint to footpoint
 *   
 *   The integral is as follows:
 *
 *      \f[
 *          V = \int_{sf_{south}}^{sf_{north}}
 *              \left\{
 *                  1\over B(s)}
 *              \right\} ds
 *      \f]
 *
 *
 *
 *
 *      \param[in,out]  mInfo    A properly initialized Lgm_MagModelInfo structure.
 *
 *      \return         V, The flux tube volume.
 *
 *
 *      \note
 *           - The routine needs the following values set properly in the mInfo structure;
 *              - mInfo->Stotal
 *              - mInfo->Lgm_V_Integrator_epsabs
 *              - mInfo->Lgm_V_Integrator_epsrel
 *              - mInfo->Lgm_V_Integrator
 *              - other things too (model info etc...)
 *      \note
 *          - On exit, the following will be set;
 *              - other things...
 *
 *
 */
double FluxTubeVolume( Lgm_MagModelInfo *mInfo ) {

    double	a, b;
    double	epsabs, epsrel, result, abserr, resabs, resasc;
    int		key, limit, lenw;
    int     iwork[501], last, ier, neval;
    double	work[2002];
    _qpInfo	*qpInfo;


    /*
     *  Type-cast our data structure to a generic type.
     *  The structure holds auzilliary info we need down 
     *  in the function calls.
     */
    qpInfo = (_qpInfo *)mInfo;

    /*
     *  Limits of integration.
     */
    a = 0.0;             // Southern footpoint
    b = mInfo->Stotal;   // Northern footpoint


    /*
     *   set tolerances for QuadPack routines. 
     */
    //epsabs = mInfo->epsabs;
    //epsrel = mInfo->epsrel;
    epsabs = mInfo->Lgm_V_Integrator_epsabs;
    epsrel = mInfo->Lgm_V_Integrator_epsrel;


    /*
     *  Init some vars used in V_integrand() (these are not declared static in
     *  V_integrand() in order to avoid making it non-reentrant).
     */
    mInfo->Lgm_n_V_integrand_Calls  = 0;



    if ( mInfo->Lgm_V_Integrator == DQAGS ) {

        /*
         *  Use DQAGS
         */
        limit = 500; lenw = 4*limit; key = 6;
        dqags(V_integrand, qpInfo, a, b, epsabs, epsrel, &result, &abserr, &neval, &ier, limit, lenw, &last, iwork, work, mInfo->VerbosityLevel );

    } else if ( mInfo->Lgm_V_Integrator == DQK21 ) {

        /*
         *  Use DQK21
         */
        dqk21(V_integrand, qpInfo, a, b, &result, &abserr, &resabs, &resasc);

    } else {

        /*
         *  Unknown Integrator
         */
        printf("FluxTubeVolume: Unknown integrator. Lgm_V_Integrator = %d\n", mInfo->Lgm_n_V_integrand_Calls );
        result = LGM_FILL_VALUE;

    }

    return( result );

}


double V_integrand( double s, _qpInfo *qpInfo ) {

    double              B, g, f;
    Lgm_MagModelInfo    *mInfo;

    /*
     *  Get pointer to our auxilliary data structure.
     */
    mInfo = (Lgm_MagModelInfo *)qpInfo;


    B = BofS( s, mInfo );

    f = 1.0/B;

    ++mInfo->Lgm_n_V_integrand_Calls;

    return(f);

}
