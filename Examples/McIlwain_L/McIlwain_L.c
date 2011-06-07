#include <Lgm_MagModelInfo.h>
#include <Lgm_ElapsedTime.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>


int main(){
    long int             Date;
    double               L, I, Bm, M, a, UTC, Ltarget;
    Lgm_Vector           u, v;
    Lgm_MagModelInfo    *mInfo = Lgm_InitMagInfo();
    Lgm_CTrans          *c = Lgm_init_ctrans( 0 );
    Lgm_ElapsedTimeInfo  tInfo;

    Lgm_ElapsedTimeInit( &tInfo, 255, 95, 0 );

    for (Ltarget=1.5; Ltarget<16.0; Ltarget += 0.1){
        //Ltarget = 5.0;
    
        Date = 20100101;
        UTC = 0.0;
        Lgm_Set_Coord_Transforms( Date, UTC, c );
        u.x = -Ltarget; u.y = 0.0; u.z = 0.0; // SM
        Lgm_Convert_Coords( &u, &v, SM_TO_GSM, c );
        
        
        mInfo->Bfield = Lgm_B_cdip;
        
        for (a=5.0; a<=90.0; a+= 5.0){
            L = Lgm_McIlwain_L( Date, UTC, &v, a, 0, &I, &Bm, &M, mInfo );
            printf("Pitch Angle: %g    McIlwain L  = %.10g   ( I, Bm, M = %g %g %g )   Delta = %g\n", a, L, I, Bm, M, Ltarget-L);
        }
    }



    Lgm_PrintElapsedTime( &tInfo );

    return(0);
}


