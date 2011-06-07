#include <Lgm_ElapsedTime.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>


int main(){
    int                 t;
    Lgm_ElapsedTimeInfo tInfo;



    /*
     *  The xterm colorizing is not 24-bit -- there are less than 256 colors
     *  accessible via escape sequences Lgm_ElapsedTimeInit() takes RGB vals
     *  between 0-255, but they will map to only 256 colors...
     */
    Lgm_ElapsedTimeInit( &tInfo, 255, 95, 0 );
    Lgm_PrintCurrentTime( &tInfo );
    for ( t=0; t<10; t++) {
        printf(" ... Do Work ...\n");
        usleep( 1000000 );
        Lgm_PrintElapsedTime( &tInfo );
    }


    /*
     * Here is an example cycling through the various colors available.
     */
    int i, j, k;
    int r, g, b;
    for (i=0; i<6; i++ ){
        r = (i>0) ? i*40+55 : 0;
        for (j=0; j<6; j++ ){
            g = (j>0) ? j*40+55 : 0;
            for (k=0; k<6; k++ ){
                b = (k>0) ? k*40+55 : 0;
                printf("    r, g, b = %3d %3d %3d\n", r, g, b );
                Lgm_ElapsedTimeInit( &tInfo, r, g, b );
                Lgm_PrintElapsedTime( &tInfo );
            }
        }
    }
 

    Lgm_PrintCurrentTime( &tInfo );
    return(0);
}


