#include <Lgm_NrlMsise00.h>
#include <Lgm_DynamicMemory.h>
#include <omp.h>
int main() {

    int     i, j, k, Ni, Nj;
    double  **Image;
    double  AP[8], D[10], T[3], SEC, GLAT, GLONG, ALT;
    Lgm_Msis00Info *p2, *p = InitMsis00();
    char    Filename[80];
    FILE    *fp;


    AP[1] = 4.0;
    AP[2] = 4.0;
    AP[3] = 4.0;
    AP[4] = 4.0;
    AP[5] = 4.0;
    AP[6] = 4.0;
    AP[7] = 4.0;
    SEC   = 29000.0;
    GLONG = 120.0;
    ALT   = 75.0;
    GTD7( 99172, SEC, ALT, 60.0, GLONG, 16.0, 150.0, 150.0, AP, 48, D, T, p );
    for (i=1; i<=9; i++ ) printf( "D[%d] = %g\n", i, D[i] );
    for (i=1; i<=2; i++ ) printf( "T[%d] = %lf\n", i, T[i] );
    printf("\n");


    // Altitude profile
    fp = fopen( "D6c.txt", "w" );
    for ( ALT=0.0; ALT<= 150.0; ALT+=1.0 ) {
        GTD7( 99172, SEC, ALT, 60.0, GLONG, 16.0, 150.0, 150.0, AP, 48, D, T, p );
        fprintf( fp, "%lf %g\n", ALT, D[6] );
    }
    fclose( fp );


    // LAT/LON MAPs at different ALTs
    Ni = 360;
    Nj = 180;
    ALT = 400.0; // km 
    LGM_ARRAY_2D( Image, Nj, Ni, double );
    for ( k=0; k<=1200; k += 1 ) {
        ALT = (double)k;
        printf("ALT = %g\n", ALT );

        {   // BEGIN PARALLEL EXECUTION

            #pragma omp parallel private(i,GLONG,j,GLAT,D,T,p2)
            #pragma omp for schedule(dynamic, 1)
            for ( i=0; i<Ni; i++ ) {
                p2 = InitMsis00();
                GLONG = 360.0*i/(double)Ni;
                for ( j=0; j<Nj; j++ ) {
                    GLAT = 180.0*j/(double)Nj - 90.0;
                    GTD7( 99172, SEC, ALT, GLAT, GLONG, 16.0, 150.0, 150.0, AP, 48, D, T, p2 );
                    Image[j][i] = D[6];
                }
                Lgm_FreeMsis00( p2 );
            }

            sprintf(Filename, "MSISE00_%04d", k);
            DumpGif( Filename, Ni, Nj, Image );
            //DumpGif2( Filename, -8.0, -1.0, Ni, Nj, Image );
        }   // END PARALLEL EXECUTION



    }

    LGM_ARRAY_2D_FREE( Image );

    Lgm_FreeMsis00( p );

    

    return(0);

}
