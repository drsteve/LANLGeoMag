/*
I'm guessing mthej _g's mean 'global'???

MODULE mod_global
  
  IMPLICIT NONE

  REAL(8) :: alpha0_g, Ekin_g, L_g, astar_g
  LOGICAL :: sum_res = .FALSE., select_bandwidth = .FALSE.

END MODULE mod_global
*/











/**
 *  \brief
 *      Computes the bounce-averaged Summer's [2005] diffusion coefficients.
 *
 *  \details
 *
 *      \param[in]      Alpha0      equatoria PA in radians.
 *      \param[in]      Ek_in       kinetic energy in MeV.
 *      \param[in]      L           L-shell parameter (dimensionless).
 *      \param[in]      astar       Cold plasma parameter.
 *      \param[out]     Daa_ba      Bounce-averaged value of Daa.
 *      \param[out]     Dap_ba      Bounce-averaged value of Dap.
 *      \param[out]     Dpp_ba      Bounce-averaged value of Dpp.
 *
 *      \return         void
 *
 *      \author         J. Koller, S. Zaharia, M. Henderson
 *      \date           2010-2011
 *
 */
void getDs_ba( double alpha0,  double E_kin,  double L,  double astar,  double *Daa_ba,  double *Dap_ba,  double *Dpp_ba) {

    double  tau, lambda_m;
    double  SinAlpha0 = sin( alpha0 );


    /*
     *  Quadpack integration tolerances.
1e-10 seems exccessive.
     */
    epsabs = 1e-10;
    epsrel = 1e-3;
   

    /*
     *  Calculate bounce period. THIS IS PROBABLY WRONG. THIS JUST LOOKS LIKE THE f() FACTOR IN AN APPROXIMATE TAU FROM ROEDERER.
*  DO THIS CORRECTLY.
     */
    tau = 1.3 - 0.56*SinAlpha0;


    /*
     *  Compute mirror latitude.
ACTUALLY returns cos(lambda_m)
     */
    lambda_m = Lgm_CdipMirrorLat( SinAlpha0 ); // returns radians
     

    // why is this commented out?
    //lambda_m = MIN( lambda_m, 25.*RadPerDeg );  // equatatorial confinement of waves.
   

    /*
     *  Pack integrand parameters into structure
*  MAKE INTO STRUCTURE.
     */
    alpha0_g = alpha0;
    Ekin_g   = E_kin;
    L_g      = L;
    astar_g  = astar;


    /*
     *  Perform integration using QUADPACK (DQAGS - adaptive integration with possible singularities).
*  NO! DQAGS does not handle singularities well. Change this if we really need it.
     */
    dqags( integaa, 0.0, lambda_m, epsabs, epsrel, Daa_ba, abserr, neval, ier, limit, lenw, last, iwork, work );
    dqags( integap, 0.0, lambda_m, epsabs, epsrel, Dap_ba, abserr, neval, ier, limit, lenw, last, iwork, work );
    dqags( integpp, 0.0, lambda_m, epsabs, epsrel, Dpp_ba, abserr, neval, ier, limit, lenw, last, iwork, work );
   

// I HAVE DOUBTS ABOUT THIS BEING A CORRECT BOUNCE-AVG
    Daa_ba /= tau;
    Dap_ba /= tau;
    Dpp_ba /= tau;


    return;
} 






/*
 *  Integrand for computing the bounce-averaged value of Daa
 */
double  SummersIntegrand_Gaa( double Lat, _qpInfo *qpInfo ) {

    double  LocPa, SinLocPa, CosLocPa;
    double  SinLat, CosLat, CosLat2, CosLat3, CosLat7;
    double  SinAlpha0_g, CosAlpha0_g2;
    double  integaa, Daa, Dap, Dpp;
    
    /*
     * Pull parameters for integrand from Lgm_SummersInfo structure.
     */
    SummersInfo  = (Lgm_SummersInfo *)qpInfo;
    SinAlpha0_g  = SummersInfo->SinAlpha0_g;    // pre-computed sin( Alpha0_g )
    CosAlpha0_g2 = SummersInfo->CosAlpha0_g2;   // pre-computed cos( Alpha0_g )*cos( Alpha0_g )

    
    /*
     * sin's and cos's are costly. A lookup table would speed this up radically.
     */
    SinLat  = sin( Lat );
    CosLat  = cos( Lat );

    CosLat2 = CosLat*CosLat;
    CosLat3 = CosLat2*CosLat;
    CosLat7 = CosLat3*CosLat3*CosLat;
    

    // The pow() should be replaced by sqrt's?.
    SinLocPa = SinAlpha0_g * pow( 4.0 - 3.0*CosLat2, 0.25) / CosLat3;
    if (SinLocPA > 1.0) SinLocPA = 1.0;
    CosLocPa = sqrt( 1.0 - SinLocPa*SinLocPa );
    


    // check on this. do we really need LocPa? Or do we only ever need sin or cos of it?
    LocPa = asin( SinlocPa );
    getDs( LocPa, Ekin_g, L_g, lat, astar_g, Daa, Dap, Dpp );

    integaa = Daa * CosLocPa * CosLat7 / CosAlpha0_g2;

    return( integaa );
  
}

/*
 *  Integrand for computing the bounce-averaged value of Dap
 */
double  SummersIntegrand_Gap( double Lat, _qpInfo *qpInfo ) {

    double  LocPa, SinLocPa, CosLocPa;
    double  SinLat, CosLat, CosLat2, CosLat3, CosLat7;
    double  SinAlpha0_g, CosAlpha0_g;
    double  integaa, Daa, Dap, Dpp;
    
    /*
     * Pull parameters for integrand from Lgm_SummersInfo structure.
     */
    SummersInfo  = (Lgm_SummersInfo *)qpInfo;
    SinAlpha0_g  = SummersInfo->SinAlpha0_g;    // pre-computed sin( Alpha0_g )
    CosAlpha0_g  = SummersInfo->CosAlpha0_g;    // pre-computed cos( Alpha0_g )

    
    /*
     * sin's and cos's are costly. A lookup table would speed this up radically.
     */
    SinLat  = sin( Lat );
    CosLat  = cos( Lat );

    CosLat2 = CosLat*CosLat;
    CosLat3 = CosLat2*CosLat;
    CosLat7 = CosLat3*CosLat3*CosLat;
    

    // The pow() should be replaced by sqrt's?.
    SinLocPa = SinAlpha0_g * pow( 4.0 - 3.0*CosLat2, 0.25) / CosLat3;
    if (SinLocPA > 1.0) SinLocPA = 1.0;
    CosLocPa = sqrt( 1.0 - SinLocPa*SinLocPa );
    


    // check on this. do we really need LocPa? Or do we only ever need sin or cos of it?
    LocPa = asin( SinlocPa );
    getDs( LocPa, Ekin_g, L_g, lat, astar_g, Daa, Dap, Dpp );

    integrand = Dap * SinLocPa * CosLat7/ ( SinAlpha0_g*CosAlpha0_g );

    return( integrand );
  
}

/*
 *  Integrand for computing the bounce-averaged value of Dpp
 */
double  SummersIntegrand_Gpp( double Lat, _qpInfo *qpInfo ) {

    double  LocPa, SinLocPa, CosLocPa;
    double  SinLat, CosLat, CosLat2, CosLat3, CosLat7;
    double  SinAlpha0_g, SinAlpha0_g2;
    double  integaa, Daa, Dap, Dpp;
    
    /*
     * Pull parameters for integrand from Lgm_SummersInfo structure.
     */
    SummersInfo  = (Lgm_SummersInfo *)qpInfo;
    SinAlpha0_g  = SummersInfo->SinAlpha0_g;     // pre-computed sin( Alpha0_g )
    SinAlpha0_g2 = SummersInfo->CosAlpha0_g2;    // pre-computed sin( Alpha0_g )*sin( Alpha0_g )

    
    /*
     * sin's and cos's are costly. A lookup table would speed this up radically.
     */
    SinLat  = sin( Lat );
    CosLat  = cos( Lat );

    CosLat2 = CosLat*CosLat;
    CosLat3 = CosLat2*CosLat;
    CosLat7 = CosLat3*CosLat3*CosLat;
    

    // The pow() should be replaced by sqrt's?.
    SinLocPa = SinAlpha0_g * pow( 4.0 - 3.0*CosLat2, 0.25) / CosLat3;
    if (SinLocPA > 1.0) SinLocPA = 1.0;
    CosLocPa = sqrt( 1.0 - SinLocPa*SinLocPa );
    


    // check on this. do we really need LocPa? Or do we only ever need sin or cos of it?
    LocPa = asin( SinlocPa );
    getDs( LocPa, Ekin_g, L_g, lat, astar_g, Daa, Dap, Dpp );

    integrand = Dpp * SinLocPa*SinLocPa * CosLat7 / ( CosLocPa*SinAlpha0_g2 );
    

    return( integrand );
  
}



// FLESH this out.
/**
 *  \brief
 *      Computes the local Summer's [2005] diffusion coefficients.
 *
 *  \details
 *
 *      \param[in]      Alpha0      equatoria PA in radians.
 *      \param[in]      Ek_in       kinetic energy in MeV.
 *      \param[in]      L           L-shell parameter (dimensionless).
 *      \param[in]      astar       Cold plasma parameter.
 *      \param[out]     Daa         Local value of Daa.
 *      \param[out]     Dap         Local value of Dap.
 *      \param[out]     Dpp         Local value of Dpp.
 *
 *      \return         void
 *
 *      \author         J. Koller, S. Zaharia, M. Henderson
 *      \date           2010-2011
 *
 */
void getDs( double alpha, double E_kin, double L, double lat, double astar, double *Daa, double *Dap, double *Dpp ) {

    // test
    *Daa = 1.0;
    *Dap = 1.0;
    *Dpp = 1.0;

    return;

}
