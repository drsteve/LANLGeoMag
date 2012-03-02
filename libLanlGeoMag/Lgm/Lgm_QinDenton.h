#ifndef LGM_QINDENTON_H
#define LGM_QINDENTON_H
#include <math.h>
#include <Lgm/Lgm_MagModelInfo.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>


#ifndef M_2PI
#define M_2PI            6.283185307179586476925286766559       /*  2PI           */
#endif


/*
 * Structures to contain Earth Observing parameters from various sources.
 */


typedef struct Lgm_QinDenton {

    int         nPnts;
    int         Verbosity;
    long int    *Date;
    double      *MJD;
    char        **IsoTimeStr;
    int         *Year, *Month, *Day, *Hour, *Minute, *Second;
    double      *ByIMF, *BzIMF, *V_SW, *Den_P, *Pdyn;
    double      *G1, *G2, *G3;
    int         *ByIMF_status, *BzIMF_status, *V_SW_status, *Den_P_status, *Pdyn_status;
    int         *G1_status, *G2_status, *G3_status;
    double      *fKp, *akp3, *Dst;
    double      *Bz1, *Bz2, *Bz3, *Bz4, *Bz5, *Bz6;
    double      *W1, *W2, *W3, *W4, *W5, *W6;
    int         *W1_status, *W2_status, *W3_status, *W4_status, *W5_status, *W6_status;

} Lgm_QinDenton;

typedef struct Lgm_QinDentonOne {

    long int    Date;
    double      JD;
    double      MJD;
    double      UTC;
    char        IsoTimeStr[80];
    int         Year, Month, Day, Hour, Minute, Second;
    double      ByIMF, BzIMF, V_SW, Den_P, Pdyn;
    double      G1, G2, G3;
    int         ByIMF_status, BzIMF_status, V_SW_status, Den_P_status, Pdyn_status;
    int         G1_status, G2_status, G3_status;
    double      fKp, akp3, Dst;
    double      Bz1, Bz2, Bz3, Bz4, Bz5, Bz6;
    double      W1, W2, W3, W4, W5, W6;
    int         W1_status, W2_status, W3_status, W4_status, W5_status, W6_status;

} Lgm_QinDentonOne;

Lgm_QinDenton   *Lgm_init_QinDenton( int Verbose );
void            Lgm_init_QinDentonDefaults( Lgm_QinDenton *q, int Verbose );
void            Lgm_destroy_QinDenton( Lgm_QinDenton *q );
void            Lgm_destroy_QinDenton_children( Lgm_QinDenton *q );
void            Lgm_read_QinDenton( long int Date, Lgm_QinDenton *q );
void            Lgm_get_QinDenton_at_JD( double JD, Lgm_QinDentonOne *p, int Verbose );
void            Lgm_set_QinDenton( Lgm_QinDentonOne *p, Lgm_MagModelInfo *m );




#endif
