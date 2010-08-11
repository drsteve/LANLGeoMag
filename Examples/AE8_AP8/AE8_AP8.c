#include <stdio.h>
#include <math.h>
#include <Lgm_AE8_AP8.h>
int main( ) {



    printf("Lgm_AE8_AP8_Flux( 1.15, 1.0, LGM_AE8MAX, LGM_INTEGRAL_FLUX, 0.5, 0.7 ) = %.13lf\n", Lgm_AE8_AP8_Flux( 1.15, 1.0, LGM_AE8MAX, LGM_INTEGRAL_FLUX, 0.5, 0.7 ) );
    printf("Lgm_AE8_AP8_Flux( 9.00, 13.9, LGM_AE8MAX, LGM_INTEGRAL_FLUX, 0.5, 0.7 ) = %.13lf\n", Lgm_AE8_AP8_Flux( 9.00, 13.9, LGM_AE8MAX, LGM_INTEGRAL_FLUX, 0.5, 0.7 ) );
    

    return(0);
}

