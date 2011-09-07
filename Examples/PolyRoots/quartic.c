#include <stdio.h>
#include <Lgm_PolyRoots.h>
int main( ) {

    int             nReal;
    double          a1, a2, a3, a4;
    double complex  z1, z2, z3, z4;

    a1 =  0.0;
    a2 = -1.0;
    a3 =  0.0;
    a4 =  1.0;

    printf("a1 = %.15g\n", a1);
    printf("a2 = %.15g\n", a2);
    printf("a3 = %.15g\n", a3);
    printf("a4 = %.15g\n", a4);



    nReal = Lgm_QuarticRoots( a1, a2, a3, a4, &z1, &z2, &z3, &z4 );

    printf("nReal = %d\n", nReal );
    printf("z1    = %.16g + %.16g i\n", creal(z1), cimag(z1) );
    printf("z2    = %.16g + %.16g i\n", creal(z2), cimag(z2) );
    printf("z3    = %.16g + %.16g i\n", creal(z3), cimag(z3) );
    printf("z4    = %.16g + %.16g i\n", creal(z4), cimag(z4) );



    return(0);
}

