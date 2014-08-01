#include <stdio.h>
#include <Lgm_Vec.h>


int main(){


  Lgm_Vector      a, b, z1, z2;
    Lgm_SlerpInfo   s;
    double          alpha;

    a.x = 1.0; a.y = 1.0; a.z = 1.0;
    b.x = 0.0; b.y = -1.0; b.z = 1.0;

    Lgm_NormalizeVector( &a );
    Lgm_NormalizeVector( &b );

    Lgm_InitSlerp( &a, &b, &s );
    

    alpha = 0.1;
    Lgm_Slerp( &a, &b, &z1, alpha, &s );

    printf("Results of slerping:\n");
    Lgm_PrintVector(&a);
    Lgm_PrintVector(&b);
    Lgm_PrintVector(&z1);

    z2.x = (1.0-alpha)*a.x + alpha*b.x;
    z2.y = (1.0-alpha)*a.y + alpha*b.y;
    z2.z = (1.0-alpha)*a.z + alpha*b.z;
    Lgm_NormalizeVector( &z2 );

    printf("Results of interpolating components (and re-normalizing):\n");
    Lgm_PrintVector(&a);
    Lgm_PrintVector(&b);
    Lgm_PrintVector(&z2);

    Lgm_VecSub(&a, &z1, &z2); 
    printf("\nThey are different by:\n");
    Lgm_PrintVector(&a);

    return(0);


}

