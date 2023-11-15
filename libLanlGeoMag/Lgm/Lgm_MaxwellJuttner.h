#ifndef LGM_MAXWELLJUTTNER_H
#define LGM_MAXWELLJUTTNER_H

/*
 *  Lgm_MaxwellJuttner.h
 *
 *  The maxwell-juttner distribution. I.e. the so-called relativistic maxwellian.
 *
 */

#include <gsl/gsl_sf_bessel.h>

double  Lgm_MaxJut( double n, double T, double Ek, double E0 );
void    Lgm_MaxJut_Derivs( double n, double T, double Ek, double E0, double *dfdn, double *dfdT );
double  Lgm_Maxwellian( double n, double T, double Ek, double E0 );


#endif

/*
 *    $Id: Lgm_MaxwellJuttner.h 46 2010-10-01 20:46:59Z mgh $
 */
