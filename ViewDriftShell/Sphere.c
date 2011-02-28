#include <math.h>
#include "Objects.h"

/*
 *  Returns both intersection points (if they exist)
 * and t will be the closest one in front of us.
 */
int SphereIntersect(SphereType *Sphere, RayType *Ray, double *tmin, double *tmax, double *t ){

    double A, B, C, D, Q;
    double t0, t1;


    A = Ray->Direction.x*Ray->Direction.x 
            + Ray->Direction.y*Ray->Direction.y 
                + Ray->Direction.z*Ray->Direction.z;

    B = 2.0*(Ray->Direction.x*Ray->Origin.x
            + Ray->Direction.y*Ray->Origin.y
                + Ray->Direction.z*Ray->Origin.z);

    C = Ray->Origin.x*Ray->Origin.x
            + Ray->Origin.y*Ray->Origin.y
                + Ray->Origin.z*Ray->Origin.z
                    - Sphere->Radius2;

    // compute discriminant
    D = B*B - 4.0*A*C;
    if ( D < 0.0 ) return(0);

    // take sqrt of dicsrim.
    D = sqrt(D); 

    // compute Q
    Q = 0.5*( (B<0.0) ?  -B - D : -B + D );

    t0 = Q/A;
    t1 = C/Q;


    if ( t0 < t1 ){
        *tmin = t0;
        *tmax = t1;
    } else {
        *tmin = t1;
        *tmax = t0;
    }


    if ( t1 < 0.0 ) {
        // both intersection points are behind us.
        return(0); 
    }

    if ( t0 < 0.0 ) {
        // then t0 is behind us and t1 is in front of us.
        *t = *tmax;
        return(1);
    } else {
        // then t0 is the closest in the forward direction
        *t = *tmin;
        return(1);
    }
    
}
