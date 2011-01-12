#include <math.h>
#include <stdlib.h>
#include "Lgm/Lgm_Vec.h"

/*
 * Create an Lgm_Vector
 */
Lgm_Vector *Lgm_CreateVector( double x, double y, double z ){
    Lgm_Vector *v;
    v = (Lgm_Vector *) calloc( 1, sizeof(*v) );
    v->x = x;
    v->y = y;
    v->z = z;
    return( v );
}


/*
 * Given vectors a and b, compute the cross product, c.
 */
void Lgm_CrossProduct(Lgm_Vector *a, Lgm_Vector *b, Lgm_Vector *c) {

    c->x = (a->y * b->z) - (a->z * b->y);
    c->y = (a->z * b->x) - (a->x * b->z);
    c->z = (a->x * b->y) - (a->y * b->x);
}


/*
 * Given vectors a and b, compute and return the dot product
 */
double Lgm_DotProduct(Lgm_Vector *a, Lgm_Vector *b) {
    return(a->x * b->x + a->y * b->y + a->z * b->z);
}


/*
 * Normalize the given vector. I.e. make it into a unit vector.
 * plus return magnitude since you already have it! 
 */
double Lgm_NormalizeVector(Lgm_Vector *a) {
    double magnitude, inv;

    magnitude = sqrt((a->x * a->x) + (a->y * a->y) + (a->z * a->z));
    if (magnitude > 0.0) {
        inv = 1.0/magnitude;
        a->x *= inv;
        a->y *= inv;
        a->z *= inv;
        return(magnitude);
    } else {
	    return(-1.0);
    }
}


/* 
 * multiply the given vector by a scalar value
 */
void Lgm_ScaleVector(Lgm_Vector *a, double value) {
    a->x *= value;
    a->y *= value;
    a->z *= value;
}




/* 
 * return the magnitude of a vector
 */
double Lgm_Magnitude( Lgm_Vector *a) {
    return( sqrt((a->x * a->x) + (a->y * a->y) + (a->z * a->z)) );
}


/* 
 *  Subtract two vectors
 *  c = a-b
 */
void Lgm_VecSub(Lgm_Vector *c, Lgm_Vector *a, Lgm_Vector *b ) {
    c->x = a->x - b->x;
    c->y = a->y - b->y;
    c->z = a->z - b->z;
}

/* 
 *  Add two vectors
 *  c = a+b
 */
void Lgm_VecAdd(Lgm_Vector *c, Lgm_Vector *a, Lgm_Vector *b ) {
    c->x = a->x + b->x;
    c->y = a->y + b->y;
    c->z = a->z + b->z;
}

/* 
 *  Find Magnitude of difference between to vectors
 */
double Lgm_VecDiffMag(Lgm_Vector *a, Lgm_Vector *b ) {
    Lgm_Vector c;
    Lgm_VecSub( &c, a, b ); // c = a-b
    return( Lgm_Magnitude( &c ) );
}

/* 
 * makes the given vector have a magnitude of mag
 */
void Lgm_ForceMagnitude(Lgm_Vector *a, double mag) {
    Lgm_NormalizeVector(a);
    Lgm_ScaleVector(a, mag);
}


/*
 *  Routine to compute A dot V
 *    double A[3][3];         -- Trans Matrix
 *    Lgm_Vector *V;              -- Lgm_Vector
 *    Lgm_Vector *Result;         -- the result of the product
 */
void Lgm_MatTimesVec(double A[3][3], Lgm_Vector *V, Lgm_Vector *Result) {
    Result->x = A[0][0]*V->x + A[1][0]*V->y + A[2][0]*V->z;
    Result->y = A[0][1]*V->x + A[1][1]*V->y + A[2][1]*V->z;
    Result->z = A[0][2]*V->x + A[1][2]*V->y + A[2][2]*V->z;
}

/*
 * Return transpose of A in B
 */
void Lgm_Transpose( double A[3][3], double B[3][3] ) {

    int i, j;

    for (i=0; i<3; i++){
        for (j=0; j<3; j++){
            B[i][j] = A[j][i];
        }
    }

}


/*
 *  Routine to compute A x B
 *    double A[3][3];         -- Trans Matrix
 *    double B[3][3];         -- Trans Matrix
 *    double Result[3][3];    -- Trans Matrix
 */
void Lgm_MatTimesMat( double A[3][3], double B[3][3], double R[3][3] ) {

    R[0][0] = A[0][0]*B[0][0] + A[1][0]*B[0][1] + A[2][0]*B[0][2];
    R[0][1] = A[0][1]*B[0][0] + A[1][1]*B[0][1] + A[2][1]*B[0][2];
    R[0][2] = A[0][2]*B[0][0] + A[1][2]*B[0][1] + A[2][2]*B[0][2];

    R[1][0] = A[0][0]*B[1][0] + A[1][0]*B[1][1] + A[2][0]*B[1][2];
    R[1][1] = A[0][1]*B[1][0] + A[1][1]*B[1][1] + A[2][1]*B[1][2];
    R[1][2] = A[0][2]*B[1][0] + A[1][2]*B[1][1] + A[2][2]*B[1][2];

    R[2][0] = A[0][0]*B[2][0] + A[1][0]*B[2][1] + A[2][0]*B[2][2];
    R[2][1] = A[0][1]*B[2][0] + A[1][1]*B[2][1] + A[2][1]*B[2][2];
    R[2][2] = A[0][2]*B[2][0] + A[1][2]*B[2][1] + A[2][2]*B[2][2];

}


/*
 * Convert Spherical-Polar coords to Cartesian
 *
 * Inputs:
 *          Lat: Latitude in Spherical Polar - Degrees
 *          Lon: Longitude in Spherical Polar - Degrees
 *            r: Geocentric Radiusin Spherical Polar
 *
 * Output:
 *            u: Cartesian vector (units of whatever you used for r)
 *
 */
void Lgm_SphToCartCoords( double Lat, double Lon, double r, Lgm_Vector *c ) {

    double  CosLat;

    Lat *= RadPerDeg;
    Lon *= RadPerDeg;
    CosLat = cos(Lat);
    
    u->x = r*CosLat*cos(Lon);
    u->y = r*CosLat*sin(Lon);
    u->z = r*sin(Lat);

    return;
}


/*
 *   $Id$
 */
