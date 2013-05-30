/*! \file Lgm_Quat.c
 *
 *  \brief Set of routines for performing operations on quaternions.
 *
 *
 *
 *  \author M.G. Henderson
 *  \date   2013
 *
 *
 *
 */

#include <stdio.h>
#include <math.h>
#include "Lgm/Lgm_Quat.h"

/*
 * Print a quaternian
 */
void Lgm_PrintQuat(double Q[4]) {
  printf("x:%lf y:%lf z:%lf w:%lf mag:%lf\n", Q[0], Q[1], Q[2], Q[3], Lgm_QuatMagnitude(Q));
  return;
}

double Lgm_NormalizeQuat( double Q[4] ) {

    double  mag, mag_inv;

    mag     = sqrt(Q[0]*Q[0] + Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3]);
    if ( mag > 1e-12 ) {

        mag_inv = 1.0/mag;
        Q[0] *= mag_inv;
        Q[1] *= mag_inv;
        Q[2] *= mag_inv;
        Q[3] *= mag_inv;

    } else {
        /*
         *  Normalize to the identity quaternion
         */
        Q[0] = Q[1] = Q[2] = 0.0;
        Q[3] = 1.0;
    }

    return( mag );

}

double Lgm_MatrixTrace( double A[3][3] ) {
    return( A[0][0] + A[1][1] + A[2][2] );
}



/*
 *  Routine to compute a quaternion rotation from a 3x3 rotation matrix.
 *    double A[3][3];         -- Trans Matrix
 *    double Q[4]             -- Quaternion (qx, qy, qz, qw)
 */
void Lgm_MatrixToQuat( double A[3][3], double *Q ) {

    double  T, SqrtT, f, finv;;

    T = 1.0 + Lgm_MatrixTrace( A );

    if ( T >= 0.0 ) {
        SqrtT = sqrt( T );
        Q[3] = 0.5*SqrtT;

        f = 0.25/Q[3];
        Q[0] = (A[2][1]-A[1][2]) * f;
        Q[1] = (A[0][2]-A[2][0]) * f;
        Q[2] = (A[1][0]-A[0][1]) * f;

    } else {

        if ( (A[0][0] > A[1][1]) && (A[0][0] > A[2][2]) ) { 
           f = sqrt( 1.0 + A[0][0] - A[1][1] - A[2][2] ) * 2.0; // f=4*qx 
           finv = 1.0/f;
           Q[0] = 0.25 * f;
           Q[1] = (A[0][1] + A[1][0]) * finv; 
           Q[2] = (A[0][2] + A[2][0]) * finv; 
           Q[3] = (A[1][2] - A[2][1]) * finv;
        } else if (A[1][1] > A[2][2]) { 
           f = sqrt( 1.0 + A[1][1] - A[0][0] - A[2][2] ) * 2.0; // f=4*qy
           finv = 1.0/f;
           Q[0] = (A[0][1] + A[1][0]) * finv; 
           Q[1] = 0.25 * f;
           Q[2] = (A[1][2] + A[2][1]) * finv; 
           Q[3] = (A[0][2] - A[2][0]) * finv;
        } else { 
           f = sqrt( 1.0 + A[2][2] - A[0][0] - A[1][1] ) * 2.0; // f=4*qz
           finv = 1.0/f;
           Q[0] = (A[0][2] + A[2][0]) * finv; 
           Q[1] = (A[1][2] + A[2][1]) * finv; 
           Q[2] = 0.25 * f;
           Q[3] = (A[0][1] - A[1][0]) * finv;
        }

    }


}

void Lgm_Quat_To_Matrix( double Q[4], double A[3][3] ) {
    A[0][0] = 1.0 - 2.0 * (Q[1] * Q[1] + Q[2] * Q[2]);
    A[0][1] = 2.0 * (Q[0] * Q[1] - Q[2] * Q[3]);
    A[0][2] = 2.0 * (Q[2] * Q[0] + Q[1] * Q[3]);

    A[1][0] = 2.0 * (Q[0] * Q[1] + Q[2] * Q[3]);
    A[1][1]= 1.0 - 2.0 * (Q[2] * Q[2] + Q[0] * Q[0]);
    A[1][2] = 2.0 * (Q[1] * Q[2] - Q[0] * Q[3]);

    A[2][0] = 2.0 * (Q[2] * Q[0] - Q[1] * Q[3]);
    A[2][1] = 2.0 * (Q[1] * Q[2] + Q[0] * Q[3]);
    A[2][2] = 1.0 - 2.0 * (Q[1] * Q[1] + Q[0] * Q[0]);
}


/*
 *  This routinem determines the quaternion implied by the rotation angle and axis.
 *
 *  Output:
 *      Q[4]  -- four element array containing Qx, Qy, Qz,Qw
 *
 *  Input:
 *
 *      Angle -- angle (in degrees) of rotation.
 *      u     -- unit vector of rotation axis.
 *
 *
 */
void Lgm_AxisAngleToQuat( Lgm_Vector *u, double Angle, double *Q ) {

    double a, sa;

    // rotate by Angle deg. around u axis
    a = 0.5*Angle*RadPerDeg; sa = sin( a );
    Q[0] = sa*u->x;
    Q[1] = sa*u->y;
    Q[2] = sa*u->z;
    Q[3] = cos( a );

}



/*
 *  This routinem determines the angle and axis of rotation implied
 *  by a (normalized) quaternion.
 *
 *  Input:
 *      Q[4]  -- four element array containing Qx, Qy, Qz, Qw
 *
 *  Output:
 *
 *      Angle -- angle (in degrees) of rotation.
 *      u     -- unit vector of rotation axis.
 *
 *
 */
void Lgm_QuatToAxisAngle( double *Q, double *Angle, Lgm_Vector *u ) {

    double  Aover2, SinAover2, SinAover2_inv;


    /*
     * Force normalization of Q
     */
    Lgm_NormalizeQuat( Q );

    if ( (1.0-fabs(Q[3])) < 1e-8 ) {

        /*
         *  Rotation angle is close to zero
         *  (Its either A/2 = 0 (i.e. A = 0) or
         *  A/2 = 180 (i.e. A = 360))
         *  Rotation axis is arbitrary
         */
        *Angle = 0.0;
        u->x = 0.0;
        u->y = 0.0;
        u->z = 1.0;

    } else {

        /*
         *  Determine Angle/2
         */
        Aover2 = acos( Q[3] );
        *Angle = 2.0*Aover2*DegPerRad;

        /* 
         * Compute 1.0/sin(Angle/2)
         */
        SinAover2 = sqrt( 1.0 - Q[3]*Q[3] );
        SinAover2_inv = 1.0/SinAover2;

        /*
         *  Compute components of rotation axis
         */
        u->x = Q[0]*SinAover2;
        u->y = Q[1]*SinAover2;
        u->z = Q[2]*SinAover2;


    }



}



/*
 *  Given a Quaternion, and a vector, computes the rotated vector.
 */
void Lgm_QuatRotateVector( double *Q, Lgm_Vector *v, Lgm_Vector *vp ) {

    double xx = Q[0]*Q[0];
    double yy = Q[1]*Q[1];
    double zz = Q[2]*Q[2];
    double ww = Q[3]*Q[3];

    double xw = Q[0]*Q[3];
    double yw = Q[1]*Q[3];
    double zw = Q[2]*Q[3];

    double xy = Q[0]*Q[1];
    double xz = Q[0]*Q[2];
    double yz = Q[1]*Q[2];

    vp->x = ww*v->x + 2.0*yw*v->z - 2.0*zw*v->y + xx*v->x + 2*xy*v->y + 2.0*xz*v->z - zz*v->x - yy*v->x;
    vp->y = 2.0*xy*v->x + yy*v->y + 2.0*yz*v->z + 2.0*zw*v->x - zz*v->y + ww*v->y - 2.0*xw*v->z - xx*v->y;
    vp->z = 2.0*xz*v->x + 2.0*yz*v->y + zz*v->z - 2.0*yw*v->x - yy*v->z + 2.0*xw*v->y - xx*v->z + ww*v->z;

    return;


}



/* 
 * return the magnitude of a vector
 */
double Lgm_QuatMagnitude( double *Q ) {
    return( sqrt( Q[0]*Q[0] + Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3]) );
}

double Lgm_QuatVecLength( double *v ) { return ( sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]) ); }
double Lgm_QuatVecDot( double *v1, double *v2 ) { return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]; }
void Lgm_QuatVecZero( double *v ) { v[0] = 0.0; v[1] = 0.0; v[2] = 0.0; }
void Lgm_QuatVecSet( double *v, double x, double y, double z ) { v[0] = x; v[1] = y; v[2] = z; }
void Lgm_QuatVecAdd( double *a, double *b, double *c ) { c[0] = a[0] + b[0]; c[1] = a[1] + b[1]; c[2] = a[2] + b[2]; }
void Lgm_QuatVecSub( double *a, double *b, double *c) { c[0] = a[0] - b[0]; c[1] = a[1] - b[1]; c[2] = a[2] - b[2]; }
void Lgm_QuatVecCopy( double *v1, double *v2 ) { int i; for (i = 0 ; i < 3 ; i++) v2[i] = v1[i]; }
void Lgm_QuatVecScale( double *v, double f ) { v[0] *= f; v[1] *= f; v[2] *= f; }
void Lgm_QuatVecNormalize( double *v ) { Lgm_QuatVecScale( v, 1.0/Lgm_QuatVecLength(v) ); }

void Lgm_QuatVecCross( double *a, double *b, double *c ) {
    c[0] = (a[1]*b[2]) - (a[2]*b[1]);
    c[1] = (a[2]*b[0]) - (a[0]*b[2]);
    c[2] = (a[0]*b[1]) - (a[1]*b[0]);
}

void Lgm_QuatCombineQuats( double Q1[4], double Q2[4], double Q[4] ) {

    double  t1[4], t2[4], t3[4];
    double  tf[4];

    Lgm_QuatVecCopy( Q1,t1 );
    Lgm_QuatVecScale( t1,Q2[3] );

    Lgm_QuatVecCopy( Q2,t2 );
    Lgm_QuatVecScale( t2,Q1[3] );

    Lgm_QuatVecCross( Q2, Q1, t3 );
    Lgm_QuatVecAdd( t1, t2, tf );
    Lgm_QuatVecAdd( t3, tf, tf );
    tf[3] = Q1[3]*Q2[3] - Lgm_QuatVecDot( Q1, Q2 );

    Q[0] = tf[0]; Q[1] = tf[1]; Q[2] = tf[2]; Q[3] = tf[3];

}



