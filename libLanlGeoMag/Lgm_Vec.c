/*! \file Lgm_Vec.c
 *
 *  \brief Routines for vector operations (e.g. dot-product, cross-product, etc.)
 *
 *
 *
 *  \author M.G. Henderson
 *  \date   1999
 *
 *
 *
 */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "Lgm/Lgm_Vec.h"
#include "Lgm/Lgm_CTrans.h"

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
 * Free an Lgm_Vector
 */
void Lgm_FreeVector( Lgm_Vector *v ) {
  free(v);
}

/*
 * Print an Lgm_Vector
 */
void Lgm_PrintVector(Lgm_Vector *v) {
  printf("x:%lf y:%lf z:%lf mag:%lf\n", v->x, v->y, v->z, Lgm_Magnitude(v));
  return;
}

/*
 * Divide each component of a vector by a scalar 
 */
void Lgm_VecDivideScalar( Lgm_Vector *v, double x) {
  v->x /= x;
  v->y /= x;
  v->z /= x;
}

/*
 * Multiply each component of a vector by a scalar 
 */
void Lgm_VecMultiplyScalar( Lgm_Vector *v, double x) {
  v->x *= x;
  v->y *= x;
  v->z *= x;
}


/*
 * Given vectors a and b, compute the angle between them in degrees
 */
double Lgm_VectorAngle(Lgm_Vector *a, Lgm_Vector *b) {
  double ang;
  Lgm_Vector  u, v;
  u = *a;
  v = *b;
  Lgm_NormalizeVector(&u);
  Lgm_NormalizeVector(&v);
  ang = acos(Lgm_DotProduct(&u, &v));
  return (ang*DegPerRad);
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
 * Given vectors a, b, and c compute the scalar triple product
 *  a \cdot (b \crossprod c)
 */
double Lgm_ScalarTripleProduct(Lgm_Vector *a, Lgm_Vector *b, Lgm_Vector *c) {
    Lgm_Vector tmp;
    Lgm_CrossProduct(b,c,&tmp);
    return( Lgm_DotProduct(a, &tmp ) );
}

/*
 * Given vectors a, b, and c compute the vector triple product
 *  a x (b x c) or (c x b) x a
 */
void Lgm_VectorTripleProduct(Lgm_Vector *a, Lgm_Vector *b, Lgm_Vector *c, Lgm_Vector *d) {
    Lgm_Vector tmp;
    Lgm_CrossProduct(b, c, &tmp);
    Lgm_CrossProduct(a, &tmp, d) ;
}


/*
 * Normalize the given vector. I.e. make it into a unit vector.
 * plus return magnitude since you already have it! 
 */
double Lgm_NormalizeVector(Lgm_Vector *a) {
    double magnitude, inv;

    magnitude = Lgm_Magnitude(a);
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


/**
 * Convert Spherical-Polar coords to Cartesian
 *
 *          \param[in] Lat Latitude in Spherical Polar - Degrees
 *          \param[in] Lon: Longitude in Spherical Polar - Degrees
 *          \param[in] r: Geocentric Radius in Spherical Polar
 *
 *          \param[out] u: Cartesian vector (units of whatever you used for r)
 *
 */
void Lgm_SphToCartCoords( double Lat, double Lon, double r, Lgm_Vector *u ) {

    double  rCosLat;

    Lat *= RadPerDeg;
    Lon *= RadPerDeg;
    rCosLat = r*cos(Lat);
    
    u->x = rCosLat*cos(Lon);
    u->y = rCosLat*sin(Lon);
    u->z = r*sin(Lat);

    return;
}


/**
 * Convert Cartesian to Spherical-Polar coords
 *
 *          \param[in] u: Cartesian vector (units of whatever you used for r)
 *
 *          \param[out] Lat: Latitude in Spherical Polar - Degrees
 *          \param[out] Lon: Longitude in Spherical Polar - Degrees
 *          \param[out] r: Geocentric Radius in Spherical Polar
 *
 */
void Lgm_CartToSphCoords( Lgm_Vector *u, double *Lat, double *Lon, double *r) {

    double sq, x2, y2, z2;
    
    x2 = u->x * u->x;
    y2 = u->y * u->y;
    z2 = u->z * u->z;
    
    *r = sqrt(x2 + y2 + z2);
    //sq = sqrt(x2 + y2);
    
    //take care of the poles
    if ((x2+y2 <= 1e-8)) {
        *Lon = 0.0;
        *Lat = (u->z < 0.0) ? -90.0 : 90.0;
    }
    else {
        *Lon = atan2(u->y,u->x)*DegPerRad;
        //*Lat = 90.0 - atan2(sq, u->z)*DegPerRad;
        *Lat = asin(u->z/(*r))*DegPerRad;
    }

    return;
}

/**
 * Copies elements of a vector to an array
 *
 *          \param[in]  u: Cartesian vector (units of whatever you used for r)
 *
 *          \param[out] A: 3-element array. A[0] set to u.x, A[1] set to u.y, A[2] set to u.z.
 *
 */
void Lgm_VecToArr( Lgm_Vector *u, double *A ) {

    A[0] = u->x;
    A[1] = u->y;
    A[2] = u->z;
    
    return;
}

/**
 * Copies elements of a 3 element array to a vector
 *
 *          \param[in]  A: 3-element array. A[0] set to u.x, A[1] set to u.y, A[2] set to u.z.
 *
 *          \param[out] u: Cartesian vector (units of whatever you used for r)
 *
 */
void Lgm_ArrToVec( double *A, Lgm_Vector *u) {

    u->x = A[0];
    u->y = A[1];
    u->z = A[2];
    
    return;
}



/**
 * Sets elements of a vector to a value
 *
 *          \param[in]  u: Cartesian vector (units of whatever you used for r)
 *
 *          \param[out] f: value to set all vector comonents to.
 *
 */
void Lgm_SetVecVal( Lgm_Vector *u, double f ) {

    u->x = f;
    u->y = f;
    u->z = f;
    
    return;

}

/**
 * Sets individual components of a vector.
 *
 *          \param[in]  u: Cartesian vector (units of whatever you used for r)
 *
 *          \param[out] x: x-component value.
 *          \param[out] y: y-component value.
 *          \param[out] z: z-component value.
 *
 */
void Lgm_SetVecElements( Lgm_Vector *u, double x, double y, double z ) {

    u->x = x;
    u->y = y;
    u->z = z;
    
    return;

}

/**
 * Sets elements of a 2-component array to a value
 *
 *          \param[in]  A: 2-element array. 
 *
 *          \param[out] f: value to set all vector comonents to.
 *
 */
void Lgm_SetArrVal2( double *A, double f ) {

    A[0] = f;
    A[1] = f;
    
    return;

}
/**
 * Sets elements of a 3-component array to a value
 *
 *          \param[in]  A: 3-element array. 
 *
 *          \param[out] f: value to set all vector comonents to.
 *
 */
void Lgm_SetArrVal3( double *A, double f ) {

    A[0] = f;
    A[1] = f;
    A[2] = f;
    
    return;

}
/**
 * Sets elements of a 4-component array to a value
 *
 *          \param[in]  A: 4-element array. 
 *
 *          \param[out] f: value to set all vector comonents to.
 *
 */
void Lgm_SetArrVal4( double *A, double f ) {

    A[0] = f;
    A[1] = f;
    A[2] = f;
    A[3] = f;
    
    return;

}


/**
 * Sets individual components of a vector.
 *
 *          \param[in]  A: 2-element array.
 *
 *          \param[out] x: x-component value.
 *          \param[out] y: y-component value.
 *
 */
void Lgm_SetArrElements2( double *A, double x, double y ) {

    A[0] = x;
    A[1] = y;
    
    return;

}


/**
 * Sets individual components of an array.
 *
 *          \param[in]  A: 3-element array.
 *
 *          \param[out] x: x-component value.
 *          \param[out] y: y-component value.
 *          \param[out] z: z-component value.
 *
 */
void Lgm_SetArrElements3( double *A, double x, double y, double z ) {

    A[0] = x;
    A[1] = y;
    A[2] = z;
    
    return;

}


/**
 * Sets individual components of an array.
 *
 *          \param[in]  A: 4-element array.
 *
 *          \param[out] a: x-component value.
 *          \param[out] b: y-component value.
 *          \param[out] c: z-component value.
 *          \param[out] d: z-component value.
 *
 */
void Lgm_SetArrElements4( double *A, double a, double b, double c, double d ) {

    A[0] = a;
    A[1] = b;
    A[2] = c;
    A[3] = d;
    
    return;

}



/**
 *  \brief
 *      Initialize a Lgm_SlerpInfo structure. 
 *  \details
 *      This routine computes \f$\sin(\phi)\f$, \f$\cos(\phi)\f$, and
 *      \f$\phi\f$. Where \f$\phi\f$ is the angle between the initial and final
 *      unit vectors.  This is done outside the main slerp routine because you
 *      may want to perform many interps with the same endpoint vectors.
 *
 *          \param[in]  a: initial unit vector.
 *          \param[in]  b: final unit vector.
 *
 *          \param[in,out] si: Lgm_SlerpInfo structure.
 *
 */
void Lgm_InitSlerp( Lgm_Vector *a, Lgm_Vector *b, Lgm_SlerpInfo *si ) {

    Lgm_Vector  v, w;

    si->CosPhi = Lgm_DotProduct( a, b );

    w = *a; Lgm_ScaleVector( &w, si->CosPhi );  // w = cos(phi)*a-hat
    Lgm_VecSub( &v, b, &w );                    // v = b-hat - cos(phi)*a-hat
    si->SinPhi = Lgm_Magnitude( &v );           // sin(phi)

    si->Phi = atan2( si->SinPhi, si->CosPhi );  // Phi

}

/**
 *  \brief
 *      Computes Spherical Linear Interpolation (SLERP) between two unit vectors.
 *  \details
 *      Spherical linear interpolation (or slerping) is a method for smoothly
 *      interpolating points on a sphere (its a spherical geometry version of
 *      linear interpolation).  Let \f$\hat{a}\f$ and \f$\hat{b}\f$ be units
 *      vectors and \f$\phi\f$ be the angle between them. We can interpolate to
 *      a fraction alpha through the rotation between them using the formula;
 *
 *          \f[ \hat{z} = { \sin((1-\alpha)\phi)\over\sin(\phi)} \hat{a} + {
 *                                \sin(\alpha\phi)\over\sin(\phi)} \hat{b} \f]
 *
 *      The routine Lgm_InitSlerp() finds \f$\sin(\phi)\f$ and \f$\phi\f$ given
 *      the unbit vectors \f$\hat{a}\f$ and \f$\hat{b}\f$. Note that for small
 *      values of \f$\phi\f$, the routine uses the approximation;
 *
 *          \f[ \hat{z} =  (1-\alpha) \hat{a} + \alpha \hat{b} \f].
 *
 *
 * 
 *          \param[in]  a: initial unit vector.
 *          \param[in]  b: final unit vector.
 *          \param[in]  z: interpolated unit vector.
 *          \param[in]  alpha: fraction of the angle phi to rotate to.
 *          \param[in]  si: Lgm_SlerpInfo structure.
 */
void Lgm_Slerp( Lgm_Vector *a, Lgm_Vector *b, Lgm_Vector *z, double alpha, Lgm_SlerpInfo *si ) {

    Lgm_Vector  w1, w2;

    /*
     * Check to see that Phi is small. If it is, use approximation.
     */
    if ( si->Phi < 1e-9 ) {
        w1 = *a; Lgm_ScaleVector( &w1, (1.0-alpha) );
        w2 = *b; Lgm_ScaleVector( &w2, alpha );
    } else {
        w1 = *a; Lgm_ScaleVector( &w1, sin( (1.0-alpha)*si->Phi )/ si->SinPhi );
        w2 = *b; Lgm_ScaleVector( &w2, sin( alpha*si->Phi )/ si->SinPhi );
    }

    Lgm_VecAdd( z, &w1, &w2 ); 


}




