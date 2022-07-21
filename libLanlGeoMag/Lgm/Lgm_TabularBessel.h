#ifndef LGM_TABULARBESSEL
#define LGM_TABULARBESSEL

#include "Lgm_DynamicMemory.h"
#include <math.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_poly.h>

/*
 * For doing tabular Bessel Functions of orders 0-14 (needed in TS07 model).
 */
typedef struct _Lgm_TabularBessel {

    int     M;              //!< Number of points in tables
    double  N;              //!< Maximum order to compute
    double  xmin;           //!< Minimum value for tabulating
    double  xmax;           //!< Maximum value for tabulating
    double  s;              //!< (xmax-xmin)/M
    double  *JnTabular_x;   //!< Table of x vals
    double  **JnTabular;    //!< Table of Jn(x) values 
    double  **JDnTabular;   //!< Table of J'n(x) values
    double  **JDDnTabular;  //!< Table of J''n(x) values

} Lgm_TabularBessel;

void    Lgm_TabularBessel_Init( int M, int N, double xmin, double xmax, Lgm_TabularBessel *tb );
void    Lgm_TabularBessel_Eval( double x, int Q, double *Jn, double *Jdn, Lgm_TabularBessel *tb );
void    Lgm_TabularBessel_Free( Lgm_TabularBessel *tb );


#endif
