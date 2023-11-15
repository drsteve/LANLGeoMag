/*! \file Lgm_Elcipse.c
 *
 *  \brief Determine if a position is in umbral or penumbral elicpse.
 *
 *  \details
 *      See discussion at Celestrak: http://www.celestrak.com/columns/v03n01/
 *
 *          
 *
 *  \author M.G. Henderson
 *  \date   2015
 *
 *
 *
 */

#include <stdio.h>
#include <math.h>
#include "Lgm/Lgm_CTrans.h"
#include "Lgm/Lgm_Vec.h"


/*
 *   Compute Earth Eclipse type for S/C position given in MOD coords (input Re).
 */
int Lgm_EarthEclipse( Lgm_Vector *u, Lgm_CTrans *c ) {

    Lgm_Vector  Rsun, Psun, Psc, Pe;
    double      Rsun_mag, Psun_mag, Psc_mag;
    double      Theta, ThetaE, ThetaS;
    int         Type;

    Psc.x = Re*u->x; //km
    Psc.y = Re*u->y; //km
    Psc.z = Re*u->z; //km
    Psc_mag = Lgm_Magnitude( &Psc ); // km

    // vector pointing from sc to earth
    Pe.x = -Psc.x;
    Pe.y = -Psc.y;
    Pe.z = -Psc.z;

    /*
     * We need (all in MOD);
     *          Rsun - Earth to Sun Vector.
     *          Psun - Satellite to Sun Vector.
     *          Pe   - Satellite to Earth Vector (=-Psc).
     *
     */
    Rsun_mag = c->earth_sun_dist*Re;    // km
    Rsun     = c->Sun;                  
    Lgm_ScaleVector( &Rsun, Rsun_mag ); // km

    // vector pointing from sc to sun
    Psun.x = Rsun.x - Psc.x;  //km
    Psun.y = Rsun.y - Psc.y;  //km
    Psun.z = Rsun.z - Psc.z;  //km
    Psun_mag = Lgm_Magnitude( &Psun ); // km
    //printf("Rsun_mag, Psun_mag, Psc_mag = %g %g %g\n", Rsun_mag, Psun_mag, Psc_mag );


    /*
     * Compute angle between Psun and Psc
     */
    Theta = acos( Lgm_DotProduct( &Psun, &Pe )/(Psun_mag*Psc_mag) );


    /*
     * Compute anglular radius of Earth and Sun as seen at S/C
     */
    //printf("Re/Psc_mag, Re, Psc_mag = %g %g %g\n", Re/Psc_mag, Re, Psc_mag );
    ThetaE = asin( Re/Psc_mag );
    ThetaS = asin( SOLAR_RADIUS/Psun_mag );
    //printf("Theta, ThetaE, ThetaS = %g %g %g\n", Theta, ThetaE, ThetaS );


    /*
     * To have an umbral eclipse, need to have ThetaE > ThetaS 
     * and Theta < ThetaE - ThetaS.
     *
     *  For a penumbral eclipse need to have |ThetaE - ThetaS| < Theta < ThetaE + ThetaS
     *
     */
    if ( ( ThetaE > ThetaS ) && ( Theta < (ThetaE - ThetaS) ) ){

        Type = LGM_UMBRAL_ECLIPSE;

    } else if ( (Theta < (ThetaE + ThetaS)) && ( Theta > fabs(ThetaE - ThetaS))) {

        Type = LGM_PENUMBRAL_ECLIPSE;

    } else {

        Type = LGM_NO_ECLIPSE;

    }


    return( Type );

}




/*
 *   Compute Moon Eclipse type for S/C position given in MOD coords (input Re).
 */
int Lgm_MoonEclipse( Lgm_Vector *u, Lgm_CTrans *c ) {

    Lgm_Vector  R_earth_to_sun, P_earth_to_sc, Rmoon, R_moon_to_sc;
    Lgm_Vector  R_moon_to_sun, R_sc_to_sun, Psc, Psun, Rsun;
    double      R_earth_to_sun_mag, Psc_mag, Psun_mag, Rsun_mag;
    double      Theta, ThetaE, ThetaS;
    int         Type;

    // Earth to sun vector in km
    R_earth_to_sun_mag = c->earth_sun_dist*Re;    // km
    R_earth_to_sun     = c->Sun;                  // unit vec.
    Lgm_ScaleVector( &R_earth_to_sun, R_earth_to_sun_mag ); // km

    // Earth to S/C vector in km
    P_earth_to_sc.x = Re*u->x; //km
    P_earth_to_sc.y = Re*u->y; //km
    P_earth_to_sc.z = Re*u->z; //km

    // Earth to Moon unit vector in MOD
    Lgm_Radec_to_Cart( c->RA_sun, c->DEC_sun, &Rmoon );
    
    // Earth to Moon Vector in km
    Lgm_ScaleVector( &Rmoon, c->EarthMoonDistance*Re );






    // Moon to S/C Vector in km
    R_moon_to_sc.x = P_earth_to_sc.x - Rmoon.x;
    R_moon_to_sc.y = P_earth_to_sc.y - Rmoon.y;
    R_moon_to_sc.z = P_earth_to_sc.z - Rmoon.z;

    // Moon to Sun Vector in km
    R_moon_to_sun.x = R_earth_to_sun.x - Rmoon.x;
    R_moon_to_sun.y = R_earth_to_sun.y - Rmoon.y;
    R_moon_to_sun.z = R_earth_to_sun.z - Rmoon.z;

    // S/C to Sun Vector in km
    R_sc_to_sun.x = R_moon_to_sun.x - R_moon_to_sc.x;
    R_sc_to_sun.y = R_moon_to_sun.y - R_moon_to_sc.y;
    R_sc_to_sun.z = R_moon_to_sun.z - R_moon_to_sc.z;


    Psc  = R_moon_to_sc;  Psc_mag  = Lgm_Magnitude( &Psc );
    Psun = R_sc_to_sun;   Psun_mag = Lgm_Magnitude( &Psun );
    Rsun = R_moon_to_sun; Rsun_mag = Lgm_Magnitude( &Rsun );

    

    //printf("Moon: Rsun_mag, Psun_mag, Psc_mag = %g %g %g\n", Rsun_mag, Psun_mag, Psc_mag );


    /*
     * Compute angle between Psun and Psc
     */
    Theta = acos( Lgm_DotProduct( &Psun, &Psc )/(Psun_mag*Psc_mag) );


    /*
     * Compute anglular radius of Earth and Sun as seen at S/C
     */
    //printf("Re/Psc_mag, Re, Psc_mag = %g %g %g\n", Re/Psc_mag, Re, Psc_mag );
    ThetaE = asin( Re/Psc_mag );
    ThetaS = asin( LUNAR_RADIUS/Psun_mag );
    //printf("Theta, ThetaE, ThetaS = %g %g %g\n", Theta, ThetaE, ThetaS );

    /*
     * To have an umbral eclipse, need to have ThetaE > ThetaS 
     * and Theta < ThetaE - ThetaS.
     *
     *  For a penumbral eclipse need to have |ThetaE - ThetaS| < Theta < ThetaE + ThetaS
     *
     */
    if ( ( ThetaE > ThetaS ) && ( Theta < (ThetaE - ThetaS) ) ){

        Type = LGM_UMBRAL_ECLIPSE;

    } else if ( (Theta < (ThetaE + ThetaS)) && ( Theta > fabs(ThetaE - ThetaS))) {

        Type = LGM_PENUMBRAL_ECLIPSE;

    } else {

        Type = LGM_NO_ECLIPSE;

    }


    return( Type );

}

