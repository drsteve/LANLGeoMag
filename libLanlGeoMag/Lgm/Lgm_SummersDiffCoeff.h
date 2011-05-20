#ifndef LGM_SUMMERS_DIFF_COEFF_H
#define LGM_SUMMERS_DIFF_COEFF_H
#include <math.h>
#include "Lgm/Lgm_CTrans.h"
#include "Lgm/Lgm_QuadPack.h"

typedef struct Lgm_SummersInfo {

    double  Alpha0;
    double  SinAlpha0;
    double  CosAlpha0;
    double  SinAlpha02;
    double  CosAlpha02;
    double  Ek_in;
    double  L;
    double  astar;

} Lgm_SummersInfo;

void    getDs_ba( double Alpha0,  double Ek_in,  double L,  double astar,  double *Daa_ba,  double *Dap_ba,  double *Dpp_ba);
void    getDs( double alpha, double Ek_in, double L, double Lat, double astar, double *Daa, double *Dap, double *Dpp );
double  CdipIntegrand_Sb( double Lat, _qpInfo *qpInfo );
double  SummersIntegrand_Gaa( double Lat, _qpInfo *qpInfo );
double  SummersIntegrand_Gap( double Lat, _qpInfo *qpInfo );
double  SummersIntegrand_Gpp( double Lat, _qpInfo *qpInfo );

#endif
