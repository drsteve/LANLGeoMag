#include <math.h>
#include "Lgm/LgmVec.h"

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
 *   Multiply two matrices together 
 *          Result = A * B
 */
void Lgm_MatTimeMat( double A[3][3], double B[3][3], double Result[3][3] ) {

    int i, j;

    for (j=0; j<3; ++j){
        for (i=0; i<3; ++i){
            Result[i][j] = A[0][j]*B[i][0] + A[1][j]*B[i][1] + A[2][j]*B[i][2];
        }
    }

}


/*
 *   $Id$
 */
