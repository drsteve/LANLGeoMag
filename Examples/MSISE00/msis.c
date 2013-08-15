#include <Lgm_NrlMsise00.h>
int main() {

    int     i;
    double  AP[8], D[10], T[3], SEC, GLONG, ALT;
    Lgm_Msis00Info *p = InitMsis00();
    FILE    *fp = fopen( "D6c.txt", "w" );

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

    GTD7( 99172, SEC, ALT, 60.0, GLONG, 16.0, 150.0, 150.0, AP, 48, D, T, p );
    for (i=1; i<=9; i++ ) printf( "D[%d] = %g\n", i, D[i] );
    for (i=1; i<=2; i++ ) printf( "T[%d] = %lf\n", i, T[i] );

    for ( ALT=0.0; ALT<= 150.0; ALT+=1.0 ) {
        GTD7( 99172, SEC, ALT, 60.0, GLONG, 16.0, 150.0, 150.0, AP, 48, D, T, p );
        fprintf( fp, "%lf %g\n", ALT, D[6] );
    }


    Lgm_FreeMsis00( p );

    
    fclose( fp );

    return(0);

}
