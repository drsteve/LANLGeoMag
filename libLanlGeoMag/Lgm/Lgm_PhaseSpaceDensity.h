#ifndef LGM_PHASE_SPACE_DENSITY_H
#define LGM_PHASE_SPACE_DENSITY_H

double Lgm_Energy_to_Mu( double E, double a, double B );
double Lgm_Mu_to_Energy( double Mu, double a, double B );
double Lgm_PsdToDiffFlux( double f, double p2c2 );
double Lgm_DiffFluxToPsd( double j, double p2c2 );


#endif
