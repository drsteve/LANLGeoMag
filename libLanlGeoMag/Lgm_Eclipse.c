/*! \file Lgm_Elcipse.c
 *
 *  \brief Determine if a position is in unmral or penumbral elicpse.
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
 *   
 */
int Lgm_Eclipse( Lgm_Vector Psc, Lgm_CTrans *c ) {

    Lgm_Vector  Rsun, Psun;
    double      Rsun_mag, Psun_mag, Psc_mag;
    double      Theta, ThetaE, ThetaS;
    int         Type;

    /*
     * We need (all in MOD);
     *          Rsun - Earth to Sun Vector.
     *          Psun - Satellite to Sun Vector.
     *          Psc  - Earth to Satellite Vector.
     *
     */
    Rsun = c->Sun; 
    Rsun_mag = Lgm_Magnitude( &Rsun );

    Psun.x = Rsun.x - Psc.x; 
    Psun.y = Rsun.y - Psc.y; 
    Psun.z = Rsun.z - Psc.z;
    Psun_mag = Lgm_Magnitude( &Psun );

    



    /*
     * Compute angle between Psun and Psc
     */
    Theta = acos( Lgm_DotProduct( &Psun, &Psc )/(Psun_mag*Psc_mag) );


    /*
     * Compute anglular radius of Earth and Sun as seen at S/C
     */
    ThetaE = asin( Re/Psc_mag );
    ThetaS = asin( SOLAR_RADIUS/Psun_mag );


    /*
     * To have an umbral eclipse, need to have ThetaE > ThetaS 
     * and Theta < ThetaE - ThetaS.
     *
     *  For a penumbral eclipse need to have |ThetaE - ThetaS| < Theta < ThetaE + ThetaS
     *
     */
    if ( ( ThetaE > ThetaS ) && ( Theta < (ThetaE - ThetaS) ) ){

        Type = LGM_UMBRAL_ECLIPSE;

    } else if ( (Theta < (ThetaE +  ThetaS)) && ( Theta > fabs(ThetaE +  ThetaS))) {

        Type = LGM_PENUMBRAL_ECLIPSE;

    } else {

        Type = LGM_NO_ECLIPSE;

    }


    return( Type );

}

