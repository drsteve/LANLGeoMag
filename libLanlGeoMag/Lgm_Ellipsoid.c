#include <math.h>
#include <stdlib.h>
#include "Lgm/Lgm_Objects.h"

/*
 * Returns both intersection points (if they exist)
 * and t will be the closest one in front of us.
 *
 * Ellipsoid has equation:
 *
 *      x^2/r_a^2 + y^2/r_a^2 + z^2/r_b^2 = 1
 *
 * I.e. its flattened in the z direction (E.g., if we're using it for the
 * spheroid of the Earth, r_a is the eq. radius and r_b is the polar radius.)
 *
 */
int Lgm_EllipsoidIntersect( EllipsoidType *Ellipsoid, RayType *Ray, double *tmin, double *tmax, double *t ){

    double A, B, C, D, Q, g;
    double t0, t1;

    g = Ellipsoid->Radius2_a/Ellipsoid->Radius2_b;

    A = Ray->Direction.x*Ray->Direction.x
            + Ray->Direction.y*Ray->Direction.y
                + g*Ray->Direction.z*Ray->Direction.z;

    B = 2.0*(Ray->Direction.x*Ray->Origin.x
            + Ray->Direction.y*Ray->Origin.y
                + g*Ray->Direction.z*Ray->Origin.z);

    C = Ray->Origin.x*Ray->Origin.x
            + Ray->Origin.y*Ray->Origin.y
                + g*Ray->Origin.z*Ray->Origin.z
                    - Ellipsoid->Radius2_a;

    // compute discriminant
    D = B*B - 4.0*A*C;
    if ( D < 0.0 ) return(0);

    // take sqrt of dicsrim.
    D = sqrt(D);

    // compute Q
    //Q = 0.5*( (B<0.0) ?  -B - D : -B + D );
    //
    //t0 = Q/A;
    //t1 = C/Q;

    Q = 0.5/A;
    t0 = (-B - D)*Q;
    t1 = (-B + D)*Q;


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


/*
 * Returns the volume of an Ellipsoid.
 *
 * Ellipsoid has equation:
 *
 *      x^2/r_a^2 + y^2/r_a^2 + z^2/r_b^2 = 1
 *
 * I.e. its flattened in the z direction (E.g., if we're using it for the
 * spheroid of the Earth, r_a is the eq. radius and r_b is the polar radius.)
 *
 *  The Volume is then:
 *
 *      4/3 \pi a^2 b
 *
 *
 */
double Lgm_EllipsoidVolume( EllipsoidType *Ellipsoid){
    double area;

    return((4./3.)*M_PI*Ellipsoid->Radius2_a*Ellipsoid->Radius_b);
}


/*
 * Define an Ellipsoid strcture.
 *
 * Ellipsoid has equation:
 *
 *      x^2/r_a^2 + y^2/r_a^2 + z^2/r_b^2 = 1
 *
 * I.e. its flattened in the z direction (E.g., if we're using it for the
 * spheroid of the Earth, r_a is the eq. radius and r_b is the polar radius.)
 *
 */
EllipsoidType *Lgm_CreateEllipsoid( Lgm_Vector *Origin, double a, double b){
    EllipsoidType *Ellipsoid;
    Ellipsoid = (EllipsoidType *) calloc( 1, sizeof(*Ellipsoid) );
    Ellipsoid->Radius_a = a;
    Ellipsoid->Radius2_a = a*a;
    Ellipsoid->Radius_b = b;
    Ellipsoid->Radius2_b = b*b;

    return( Ellipsoid );
}

/*
 * Free an EllipsoidType
 */
void Lgm_FreeEllipsoid( EllipsoidType *Ellipsoid ) {
    free(Ellipsoid);
}

