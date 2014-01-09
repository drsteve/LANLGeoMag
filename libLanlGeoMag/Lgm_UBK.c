#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>


void Lgm_UBK_InitModelInfo( double Kp, double gamma, double a, Lgm_UBK_ModelInfo *uInfo ) {


    Info->gamma       = gamma;
    Info->Kp          = Kp;
    Info->KpModel     = KpModel;
    Info->Be          = 30400.0;


    Info->c           = Info->Be;
    Info->a       = a;



    Info->Bgeo        = 30400.0/pow(6.6, 3.0);
    Info->Side        = 1.0;



    if (KpModel == 0){

        Info->b       = 0.045/pow( 1.0 - 0.159*Kp + 0.0093*Kp*Kp, 3.0 );

    } else if (KpModel == 1){

        Info->l_m     = 67.8 - 2.07*Kp;
        Info->cos_l_m = cos( Info->l_m * M_PI/180.0 );
        Info->b       = a*pow( Info->cos_l_m*Info->cos_l_m/(1.0+1.0/gamma) , gamma+1.0 )/gamma;

    }

    printf("Rs = %f\n", pow( a/(gamma*Info->b), 1.0/(1.0+gamma) ));
    printf("a, b = %f %f\n",  Info->a, Info->b );



}



