#ifndef LGM_EOP_H
#define LGM_EOP_H
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>


#ifndef M_2PI
#define M_2PI            6.283185307179586476925286766559       /*  2PI           */
#endif


/*
 * Structures to contain Earth Observing parameters from various sources.
 */


typedef struct Lgm_NgaEopp {

    double  ta, A, B, C1, C2, D1, D2, P1;
    double  P2, E, F, G1, G2, H1, H2, Q1, Q2;
    double  tb, I, J, K1, K2, K3, K4;
    double  L1, L2, L3, L4, R1, R2, R3, R4;
    int     dat, EOPPWk, teff;

} Lgm_NgaEopp;

typedef struct Lgm_Eop {

    long int    Size;     // Size of calloced arrays
    long int    nEopVals; // Number of values
    int         Verbosity;
    long int    *Date;
    double      *MJD;
    double      *xp;
    double      *yp;
    double      *DUT1;
    double      *LOD;
    double      *dPsi;  // actually ddPsi
    double      *dEps;  // actually ddEps
    double      *dX;    // actually ddX?
    double      *dY;    // actually ddY?
    double      *DAT;
    
} Lgm_Eop;

typedef struct Lgm_EopOne {

    long int    Date;
    double      JD;
    double      MJD;
    double      UTC;
    double      xp;
    double      yp;
    double      DUT1;
    double      LOD;
    double      dPsi;  // actually ddPsi
    double      dEps;  // actually ddEps
    double      dX;    // actually ddX?
    double      dY;    // actually ddY?
    double      DAT;
    
} Lgm_EopOne;

Lgm_Eop *Lgm_init_eop( int Verbose );
void    Lgm_destroy_eop( Lgm_Eop *e );
void    Lgm_read_eop( Lgm_Eop *e );
void    Lgm_NgaEoppPred( double JD, Lgm_EopOne *eop, Lgm_NgaEopp *e );
int     Lgm_ReadNgaEopp( Lgm_NgaEopp *e, int Verbosity );
void    Lgm_get_eop_at_JD( double JD, Lgm_EopOne *eop, Lgm_Eop *e );
void    Lgm_set_eop( Lgm_EopOne *eop, Lgm_CTrans *c );
void    Lgm_unset_eop( Lgm_EopOne *eop, Lgm_CTrans *c );




#endif

/*
 *    $Id: Lgm_Eop.h 46 2010-10-01 20:46:59Z mgh $
 */

