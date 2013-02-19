#include <stdio.h>
#include <Lgm_Vec.h>


int main(){


    Lgm_Vector      a, b, z;
    Lgm_SlerpInfo   s;
    double          alpha;

    a.x = 1.0; a.y = 1.0; a.z = 1.0;
    b.x = 0.0; b.y = -1.0; b.z = 1.0;

    Lgm_NormalizeVector( &a );
    Lgm_NormalizeVector( &b );

    Lgm_InitSlerp( &a, &b, &s );
    

    alpha = 0.1;
    Lgm_Slerp( &a, &b, &z, alpha, &s );

    printf("Results of slerping:\n");
    printf("a = %g %g %g\n", a.x, a.y, a.z);
    printf("b = %g %g %g\n", b.x, b.y, b.z);
    printf("z = %g %g %g\n", z.x, z.y, z.z);




    z.x = (1.0-alpha)*a.x + alpha*b.x;
    z.y = (1.0-alpha)*a.y + alpha*b.y;
    z.z = (1.0-alpha)*a.z + alpha*b.z;
    Lgm_NormalizeVector( &z );

    printf("Results of interpolating components (and re-normalizing):\n");
    printf("a = %g %g %g\n", a.x, a.y, a.z);
    printf("b = %g %g %g\n", b.x, b.y, b.z);
    printf("z = %g %g %g\n", z.x, z.y, z.z);

    return(0);


}

