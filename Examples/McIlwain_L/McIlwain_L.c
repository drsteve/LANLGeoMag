#include <Lgm_MagModelInfo.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>


int main(){
    long int            Date;
    double              L, I, Bm, M, a, UTC;
    Lgm_Vector          u, v;
    Lgm_MagModelInfo    *mInfo = Lgm_InitMagInfo();
    Lgm_CTrans          *c = Lgm_init_ctrans( 0 );



    
    
    Date = 20100101;
    UTC = 12.0;
    Lgm_Set_Coord_Transforms( Date, UTC, c );
    u.x = -16.6; u.y = 0.0; u.z = 0.0; // SM
    Lgm_Convert_Coords( &u, &v, SM_TO_GSM, c );
                                                     

    
    
    mInfo->Bfield = Lgm_B_cdip;

    

    
    for (a=10.0; a<=90.0; a+= 5.0){
        L = Lgm_McIlwain_L( Date, UTC, &v, a, 0, &I, &Bm, &M, mInfo );
        printf("Pitch Angle: %g    McIlwain L  = %.10g   ( I, Bm, M = %g %g %g )   Delta = %g\n", a, L, I, Bm, M, 16.6-L);
    }

//    L = Lgm_McIlwain_L( 20100101, 12.0, &u, 45.0, 1, &I, &Bm, &M, mInfo );
//    printf("McIlwain L  = %.10g   ( I, Bm, M = %g %g %g )\n", L, I, Bm, M);


    return(0);
}


