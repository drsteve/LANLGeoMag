#include <stdio.h>
#include <math.h>

void Lgm_ComputeW( double W[], int i, double Nk[], double Vk[], double Bsk[], int nk ){

    /*
     *   Nk, Vk, Bsk are arrays of SW/IMF inputs. Assume they are 1 minute resolution.
     *
     *   Nk[] is array of solar wind density values from beginning of storm
     *   Vk[] is array of solar wind speed values from beginning of storm
     *   Bsk[] is array of IMF Bs values from beginning of storm
     *   nk is number of values in these arrays
     *
     *	 i is the index of the current time
     *   W[] are the 6 W index parameters defined in T04S
     *
     */
     



    static double	Lambda[6] = { 0.39, 0.46, 0.39, 0.42, 0.41, 1.29 };
    static double	Beta[6]  = { 0.80, 0.18, 2.32, 1.25, 1.6, 2.4 };
    static double	Gamma[6] = { 0.87, 0.67, 1.32, 1.29, 0.69, 0.53 };
    static double	Rate[6]  = { 0.39, 0.7, 0.031, 0.58, 1.15, 0.88 };
    double		    ti, tk, sum, Sk;
    int			    k, q;



    
    /*
     *  Time is measured in minutes. But each step is 5 minutes so we need to multiply by 5 on the indices...
     */
    ti = (double)i*5.0;

    
    for (q=0; q<6; ++q) {

        for (sum=0.0, k=0; k<=i; ++k){

	        tk = (double)k*5;
	        Sk = pow( Nk[k]/5.0, Lambda[q] ) * pow( Vk[k]/400.0, Beta[q] ) * pow( Bsk[k]/5.0, Gamma[q] );
	        sum += Sk*exp(Rate[q]/60.0*(tk-ti));

	    }
	    W[q] = Rate[q]/12.0*sum;

    }


}
/*
 *   $id$
 */
