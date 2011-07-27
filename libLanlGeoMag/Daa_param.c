#include <math.h>
#include "Lgm/Lgm_SummersDiffCoeff.h"
/*
 * Daa_param.c
 */

/**
 *
 *\brief 
 *    This subroutine calculates the bounce-averaged pitch-angle diffusivity
 *    (radians^2/s) using a parameterization  of the Summers (2005) model for:
 *
 *      dB = hiss amplitude (Tesla)
 *      alpha = equatorial pitch angle (Radians)
 *      f_cen = frequency center of hiss (Hz)
 *      E_MeV = electron energy (MeV)
 *      L = L-Shell
 *
 *      and a background plasma density of the form
 *
 *      No = (3500 cm^(-3))*(2/L)^gamma
 *
 *      such that the concentration is fixed at 3500 cm^(-3) at L=2, and
 *      decreases more rapidly with increasing gamma
 *
 *      Valid Parameter Ranges:
 *
 *(1) 0.05 < E_MeV < 10    MeV
 *(2) 550 < f_cen < 16000  Hz
 *(3) 1 < L < 4.5
 *(4) 2.5 < gamma < 5
 *
 *Hiss frequency Spectrum:
 *
 *The hiss frequency spectrum is a Gaussian centered at f_cen (Hz) and has a standard 
 *deviation of 300 Hz
 *
 *Hiss Amplitude:
 *
 *Note that the input parameter dB appears as a simple multiplicative factor, (dB)^2,
 *in the expression for Daa.  Hence, one can calculate Daa using dB=1, and then
 *trivially scale afterwards for different dB.
 *
 *----
 *
 *Version 1.0
 *Christopher Jeffery
 *7/14/2011
 *
 */

double Daa_param(double dB, double alpha, double f_cen, double E_MeV, double L, double gamma)
{

  /*
    Here are the fitting parameters
  */

  static const double a1_L = 1.3 ;
  static const double a1_mu = 2.7 ;
  static const double a1_fcrit = -2.6 ;
  static const double a2_fcrit = -1.4 ;
  static const double a1_asmp = 1.2 ;
  
  static const double a1_fit = 2.5e17 ;
  static const double a1_mu_f = 0.7 ;
  static const double a2_mu_f = 0.35 ;
  static const double a1_const = 1.5 ;
  
  static const double a1_sig = 0.016 ;
  static const double a1_exp_mu = 0.5 ;
  static const double a1_exp_f = 0.26 ;
  static const double a1_exp_E = 0.47 ;

  /*
    Here is the fitting function
  */

  double Daa ;
  double f_0 = 550. ;
  double mu = cos(alpha) ;  
  double fcrit,mu_f_expnt,fnorm,fcrit_expnt,fcrit_sig,E,E_kin ;
  double get_fcrit(double,double,double,double) ;

  E_kin = E_MeV*LGM_e*1e6 ;    // kinetic energy in [joule]
  E = E_kin/(LGM_ELECTRON_MASS*LGM_c*LGM_c) ;  // dimensionless particle energy

  Daa = pow(dB,2.)*a1_fit*pow(L,a1_L)/pow(E+1,2.)*exp(-a1_mu*pow(fabs(mu-0.2),1.5)) ;

  fcrit = get_fcrit(alpha,E_MeV,L,gamma) ;

  mu_f_expnt = a1_mu_f - a2_mu_f*mu ;
  Daa = Daa*pow(f_0/f_cen,mu_f_expnt) ;

  fnorm = f_cen/fcrit ;

  if(f_cen>=fcrit)
    {
      fcrit_expnt = a1_fcrit*exp(a2_fcrit*mu) ;
      Daa = Daa*pow((1 + a1_const*pow(fnorm-1,a1_asmp)),fcrit_expnt/a1_asmp) ;
    }
  else
    {
      fcrit_sig = a1_sig*exp(a1_exp_mu*(1-mu) + a1_exp_f*(f_cen/f_0) + a1_exp_E*(E_MeV-1.)) ;
      Daa = Daa*exp(-pow(fnorm-1,4.)/fcrit_sig) ;
    }

  return(Daa) ;
}

/*
This subroutine calculates the frequency that maximizes the pitch angle diffusivity for
a given alpha, E_MeV, L and gamma
*/

double get_fcrit(double alpha, double E_MeV, double L, double gamma)
{
  /*
    Here are the fitting parameters
  */

  static const double c1_beta = -7.5 ;
  static const double c2_beta = 0.8 ;
  static const double c1_alpha = 400 ;

  /*
    Here is the fitting function
  */

  double fcrit ;
  double f_0 = 550. ;
  double mu = cos(alpha) ;

  fcrit = f_0*pow( (c1_alpha/pow(gamma,2.)/pow(mu,1.6))/E_MeV*
		   pow(L,c1_beta + c2_beta*gamma), 1/0.9) ;
  
  return(fcrit) ;
}
