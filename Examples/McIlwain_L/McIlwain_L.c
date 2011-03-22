#include <Lgm_MagModelInfo.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>


int main(){
    double              L, I, Bm, M;
    Lgm_Vector          u;
    Lgm_MagModelInfo    *mInfo = Lgm_InitMagInfo();


    u.x = -6.6; u.y = 0.0; u.z = 0.0;

    
    L = Lgm_McIlwain_L( 20100101, 12.0, &u, 45.0, 0, &I, &Bm, &M, mInfo );
    printf("McIlwain L  = %.10g   ( I, Bm, M = %g %g %g )\n", L, I, Bm, M);

    L = Lgm_McIlwain_L( 20100101, 12.0, &u, 45.0, 1, &I, &Bm, &M, mInfo );
    printf("McIlwain L  = %.10g   ( I, Bm, M = %g %g %g )\n", L, I, Bm, M);


    return(0);
}


