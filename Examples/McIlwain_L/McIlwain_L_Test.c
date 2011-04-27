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

    //vals = magcoords.Lvalue([-4, 0, 1], datetime.datetime(2009,1,1))
    Date = 20090101;
    UTC  = 0.0;
    mInfo->Kp = 2;
    Lgm_Set_Coord_Transforms( Date, UTC, c );
    u.x = -4.0; u.y = 0.0; u.z = 1.0;
    Lgm_Convert_Coords( &u, &v, GSM_TO_WGS84, c );
    printf("v_wgs84 = %.15lf %.15lf %.15lf\n", v.x, v.y, v.z);
    a = 90.0;

    L = Lgm_McIlwain_L( Date, UTC, &u, a, 1, &I, &Bm, &M, mInfo );
    printf("Pitch Angle: %g    McIlwain L  = %.15lf   ( I, Bm, M = %.15g %g %g )\n", a, L, I, Bm, M);



    return(0);
}


