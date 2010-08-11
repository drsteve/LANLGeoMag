#define MAIN
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <Lgm_CTrans.h>





int main( int argc, char *argv[] ){

    char         Command[1024], HostName[128], UserName[128], UsedFluxFilename[512];
    int          flag, ii, jj, i, ny, nm, nd, DifferencesExist, mode, FluxesOK;
    double       sJD, eJD, cJD, JD, UT, r, diff, old, val, phi, theta;
    long int     Date;
    Lgm_CTrans   *c = Lgm_init_ctrans(0);
    Lgm_Vector   v, B;
    FILE         *fp2;
    char         Filename2[1024];

c->Lgm_IGRF_FirstCall = TRUE;
    
    fp2 = fopen("results.dat", "w");

    Lgm_Set_Coord_Transforms( 20100625, 19.0, c );
    
for (phi=60.0; phi<=60.0; phi+=30.0){
for (theta=0.0; theta<=55.0; theta+=1.0){
    flag = 1;
    for (r=1.1; r<10.0; r+=0.1){


        v.x = r*cos(RadPerDeg*phi)*cos(RadPerDeg*theta);
        v.y = r*sin(RadPerDeg*phi)*cos(RadPerDeg*theta);
        v.z = r*sin(RadPerDeg*theta);
//v.x = r;
//v.y = theta*RadPerDeg;
//v.z = phi*RadPerDeg;
 //       Lgm_IGRF( &v, &B, c );
        Lgm_B_edip_ctrans( &v, &B, c );
        val = fabs(B.x);
        val = Lgm_Magnitude(&B);

        if ( flag == 0 ) {
            diff = val-old;
            //printf("B = %g %g %g\n", B.x, B.y, B.z);
            //fprintf(fp2, "%g %g\n", r, log10(fabs(B.z)));
            //fprintf(fp2, "%g %g\n", r, log10(fabs(diff)));
            fprintf(fp2, "%g %g\n", r, fabs(val));
        }
        
        flag = 0;
        old = val;

    }
}
}


    free(c);

    

}





