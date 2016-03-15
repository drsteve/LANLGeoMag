/*! \file Lgm_Quat.c
 *
 *  \brief Set of routines for performing operations on quaternions.
 *
 *  \details
 *      See papers by;
 *
 *          1. Dam, E. B., M. Kock, M. Lillholm, Quaternions, Interpolation and
 *             Animation, Technical Report DIKU-TR-98/5, Department of Computer
 *             Science, University of Copenhagen, Universitetsparken 1, DK-2100
 *             Kbh 0, Denmark (July 17, 1998).
 *
 *          2. Eberly, D., Quaternion Algebra and Calculus, Geometric Tools,
 *             LLC,  http://www.geometrictools.com/Documentation/Quaternions.pdf (August 18, 2010).
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

double Lgm_QuatNorm( double Q[4] ) {
    return( Q[0]*Q[0] + Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3] );
}

void Lgm_QuatConjugate( double Qin[4], double Qout[4] ) {
    Qout[0] = -Qin[0];
    Qout[1] = -Qin[1];
    Qout[2] = -Qin[2];
    Qout[3] =  Qin[3];
}

void Lgm_QuatInverse( double Qin[4], double Qout[4] ) {

    Lgm_QuatConjugate( Qin, Qout );
    Lgm_QuatScale( Qout, 1.0/Lgm_QuatNorm( Qin ) );

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
 *  \brief
 *      Routine to compute a quaternion rotation from a 3x3 rotation matrix.
 *  \details
 *    double A[3][3];         -- Trans Matrix
 *    double Q[4]             -- Quaternion (qx, qy, qz, qw)
 *
 *          \param[in]  A: transformation matrix.
 *          \param[out] Q: 
 *
 */
void Lgm_MatrixToQuat( double A[3][3], double *Q ) {

    double  T, SqrtT, f, finv;;

    T = 1.0 + Lgm_MatrixTrace( A );

    if ( T >= 0.0 ) {
        SqrtT = sqrt( T );
        Q[3] = 0.5*SqrtT;

        f = 0.25/Q[3];
        Q[0] = (A[1][2]-A[2][1]) * f;
        Q[1] = (A[2][0]-A[0][2]) * f;
        Q[2] = (A[0][1]-A[1][0]) * f;

    } else {

        if ( (A[0][0] > A[1][1]) && (A[0][0] > A[2][2]) ) { 
           f = sqrt( 1.0 + A[0][0] - A[1][1] - A[2][2] ) * 2.0; // f=4*qx 
           finv = 1.0/f;
           Q[0] = 0.25 * f;
           Q[1] = (A[1][0] + A[0][1]) * finv; 
           Q[2] = (A[2][0] + A[0][2]) * finv; 
           Q[3] = (A[2][1] - A[1][2]) * finv;
        } else if (A[1][1] > A[2][2]) { 
           f = sqrt( 1.0 + A[1][1] - A[0][0] - A[2][2] ) * 2.0; // f=4*qy
           finv = 1.0/f;
           Q[0] = (A[1][0] + A[0][1]) * finv; 
           Q[1] = 0.25 * f;
           Q[2] = (A[2][1] + A[1][2]) * finv; 
           Q[3] = (A[2][0] - A[0][2]) * finv;
        } else { 
           f = sqrt( 1.0 + A[2][2] - A[0][0] - A[1][1] ) * 2.0; // f=4*qz
           finv = 1.0/f;
           Q[0] = (A[2][0] + A[0][2]) * finv; 
           Q[1] = (A[2][1] + A[1][2]) * finv; 
           Q[2] = 0.25 * f;
           Q[3] = (A[1][0] - A[0][1]) * finv;
        }

    }


}

void Lgm_Quat_To_Matrix( double Q[4], double A[3][3] ) {
    A[0][0] = 1.0 - 2.0 * (Q[1] * Q[1] + Q[2] * Q[2]);
    A[1][0] = 2.0 * (Q[0] * Q[1] - Q[2] * Q[3]);
    A[2][0] = 2.0 * (Q[2] * Q[0] + Q[1] * Q[3]);

    A[0][1] = 2.0 * (Q[0] * Q[1] + Q[2] * Q[3]);
    A[1][1]= 1.0 - 2.0 * (Q[2] * Q[2] + Q[0] * Q[0]);
    A[2][1] = 2.0 * (Q[1] * Q[2] - Q[0] * Q[3]);

    A[0][2] = 2.0 * (Q[2] * Q[0] - Q[1] * Q[3]);
    A[1][2] = 2.0 * (Q[1] * Q[2] + Q[0] * Q[3]);
    A[2][2] = 1.0 - 2.0 * (Q[1] * Q[1] + Q[0] * Q[0]);
}


/*
/*
 *  \brief
 *      Routine determines the quaternion implied by the rotation angle and axis.
 *  \details
 *      The Quaternion is of the form;
 *            \f[ Q = ( \sin(\theta/2) \hat{u}, \cos(\theta/2) ) = ( \sin(\theta/2) u_x, \sin(\theta/2) u_y, \sin(\theta/2) u_z, \cos(\theta/2) ) \f].
 *
 *          \param[in]  u: pointer to unit vector rotation axis.
 *          \param[in]  Angle: angle to rotate about rotation axis.
 *          \param[out] Q[4]: resulting quatoernion
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
        Lgm_NormalizeVector( u );


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

void Lgm_QuatScale( double *Q, double t ) { Q[0] *= t; Q[1] *= t; Q[2] *= t; Q[3] *= t; }
void Lgm_QuatAdd( double *Q1, double *Q2, double *Q3 ) { Q3[0] = Q1[0]+Q2[0]; Q3[1] = Q1[1]+Q2[1]; Q3[2] = Q1[2]+Q2[2]; Q3[3] = Q1[3]+Q2[3]; }

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




// Q = Q1*Q2
void Lgm_QuatMultiply( double Q1[4], double Q2[4], double Q[4] ) {

    // vector part
    Q[0] = Q1[3]*Q2[0] + Q1[0]*Q2[3] + Q1[1]*Q2[2] - Q1[2]*Q2[1];
    Q[1] = Q1[3]*Q2[1] - Q1[0]*Q2[2] + Q1[1]*Q2[3] + Q1[2]*Q2[0];
    Q[2] = Q1[3]*Q2[2] + Q1[0]*Q2[1] - Q1[1]*Q2[0] + Q1[2]*Q2[3];

    // real part
    Q[3] = Q1[3]*Q2[3] - Q1[0]*Q2[0] - Q1[1]*Q2[1] - Q1[2]*Q2[2];

}





/**
 *  \brief
 *      Computes q^t .
 *  \details
 *      Exponentiates a quaternion. If,
 *            \f[  q = ( \sin(\theta) \hat{u}, \cos(\theta) ) \f].
 *      then,
 *
 *            \f[  q^t = ( \sin(t \theta) \hat{u}, \cos(t \theta) ) \f].
 *
 */
void Lgm_QuatPow( double Qin[4], double t, double Qout[4] ) {

    double      Angle, Theta, tTheta, s, c;
    Lgm_Vector  u;

    Lgm_QuatToAxisAngle( Qin, &Angle, &u );
    Theta  = Angle/2.0*RadPerDeg;

    tTheta = t*Theta;
    s = sin( tTheta );
    c = cos( tTheta );
    Qout[0] = s*u.x;
    Qout[1] = s*u.y;
    Qout[2] = s*u.z;
    Qout[3] = c;

}



/**
 *  \brief
 *      Computes Spherical Linear Quaternion Interpolation (SLERP).
 *  \details
 *      Spherical linear interpolation (or slerping) is a method for smoothly
 *      interpolating points on a sphere (its a spherical geometry version of
 *      linear interpolation).  
 */
void Lgm_QuatSlerp( double p[4], double q[4], double h, double Qout[4] ){

    double  A[4], B[4], pconj[4];

    /*
    Lgm_QuatPow( p, 1.0-h, A );
    Lgm_QuatPow( q, h, B );
    Lgm_QuatMultiply( A, B, Qout );
    */

    // probably more efficient
    Lgm_QuatConjugate( p, pconj );
    Lgm_QuatMultiply( pconj, q, A );
    Lgm_QuatPow( A, h, B );
    Lgm_QuatMultiply( p, B, Qout );
    
}



/**
 *  \brief
 *      Computes Spherical Spline Quaternion Interpolation (SQUAD).
 *  \details
 *      Spherical spline interpolation is a method for smoothly
 *      interpolating points on a sphere.
 */
void Lgm_QuatSquad( double Q0[4], double Q1[4], double S0[4], double S1[4], double h, double Qout[4] ){

    double  A[4], B[4], hh;

    //Lgm_QuatSlerp( Q0, Q1, h, Qout );
    //return;

    Lgm_QuatSlerp( Q0, Q1, h, A );
    Lgm_QuatSlerp( S0, S1, h, B );
    hh = 2.0*h*(1.0-h);

    Lgm_QuatSlerp( A, B, hh, Qout );

}




/*
 * Take log() of a unit quaternion. It will be purely 'imaginary' -- real part
 * is zero.  Result is not necessarily a unit quaternion.
 */
void Lgm_QuatLog( double Q[4], double logQ[4] ) {

    double     Angle, Theta;
    Lgm_Vector u;

    Lgm_QuatToAxisAngle( Q, &Angle, &u );

    // the Angle is actually 2*Theta
    Theta = 0.5*Angle*RadPerDeg;

    logQ[0] = Theta*u.x;
    logQ[1] = Theta*u.y;
    logQ[2] = Theta*u.z;
    logQ[3] = 0.0;

}



/*
 * Take exp() of a purely imaginary quaternion.
 */
void Lgm_QuatExp( double Q[4], double expQ[4] ) {

    double     Theta, SinTheta, CosTheta;
    Lgm_Vector u;

    u.x = Q[0];
    u.y = Q[1];
    u.z = Q[2];
    Theta = Lgm_NormalizeVector( &u );
    SinTheta = sin( Theta );
    CosTheta = cos( Theta );

    expQ[0] = SinTheta*u.x;
    expQ[1] = SinTheta*u.y;
    expQ[2] = SinTheta*u.z;
    expQ[3] = CosTheta;

}



// Eqn 6.15 of Dam et al.
void Lgm_QuatSquadComputeAuxPoint( double Qim1[4], double Qi[4], double Qip1[4],   double si[4] ) {

    double Qi_inv[4], A[4], B[4], logA[4], logB[4], C[4], D[4];

    Lgm_QuatInverse( Qi, Qi_inv );
    Lgm_QuatMultiply( Qi_inv, Qip1, A ); Lgm_QuatLog( A, logA );
    Lgm_QuatMultiply( Qi_inv, Qim1, B ); Lgm_QuatLog( B, logB );

    Lgm_QuatAdd( logA, logB, C );
    Lgm_QuatScale( C, -0.25 );
    Lgm_QuatExp( C, D );
    Lgm_QuatMultiply( Qi, D, si );

}


long int Lgm_QuatFindTimeIndex( double *T, long int N, double t ){

    long int jl, jh, jm;

    jl = -1;
    jh =  N;
    while ( jh-jl > 1 ) {

        jm = (jh+jl)/2;
        if ( t > T[jm] ) {
            jl = jm;
        } else {
            jh = jm;
        }
    }

    return( jl );

}


/**
 *  \brief
 *      Computes Spherical Spline Quaternion Interpolation (SQUAD) on arrays of quaternions.
 *  \details
 *      Spherical spline interpolation is a method for smoothly
 *      interpolating points on a sphere.
 *
 *       Input arrays to interpolate from:
 *          \param[in]  T: array of input 'time' values.
 *          \param[in]  Q: array of input quaternions (
 *          \param[in]  N: number of quaternions in (T, Q) arrays
 *
 *       Interpolated arrays:
 *          \param[in]  t: array of 'time' values to interpolate to.
 *          \param[out] q: array of resulting quaternions.
 *          \param[in]  n: number of quaternions in (t, q) arrays
 *
 */
int Lgm_QuatSquadInterp( double *T,   double *Q[4], long int N,   double *t, double *q[4], long int n ) {

    long int    i, j, im1, ip1, ip2;
    double      si[4], sip1[4], h;



    for (i=0; i<N; i++){
        printf("\tT[%d] = %g Q[%d] = ", i, T[i], i ); Lgm_PrintQuat( Q[i] );
    }




    /*
     * Check that there are enough points
     */
    if ( N < 2 ) {

        printf( "Lgm_QuatSquadInterp: Not enough points to interpolate: N = %d\n", N );
        return( -1 );

    } else if ( N==2 ) {

        printf( "Lgm_QuatSquadInterp: Not enough points to perform a SQUAD interpolation: N = %d. Do a SLERP instead\n", N );
        // PUT code HERE to PERFORM SLERP
        return( 0 );

    }
     


    /*
     * Loop over all time values contained in the t array
     */ 
    for ( j=0; j<n; j++ ){

        /*
         *  Locate where t[i] falls in the control point time array.
         *
         *      if i<0, t[i] falls below all of the control points and an
         *      extrapolation would be needed on the low side.
         *
         *      if i>=N-1, t[i] falls above all of the control points and an
         *      extrapolation would be needed on the hiugh side.
         *
         */
        i = Lgm_QuatFindTimeIndex( T, N, t[j] );

        if ( i<0 ) {

            printf( "Lgm_QuatSquadInterp: Warning. Extrapolation required on the low side. t[%ld] = %g. i = %ld\n", j, t[j], i );
            return(-1);

        } else if ( i>=N-1 ) {

            printf( "Lgm_QuatSquadInterp: Warning. Extrapolation required on the high side. t[%ld] = %g. i = %ld\n", j, t[j], i );
            return(-1);

        } else if ( i == 0 ) {

            printf( "Lgm_QuatSquadInterp: Warning. i == 0, can only do SLERP on [0,1] interval.\n");
            ip1 = i+1;
            h = (t[j] - T[i])/(T[ip1] - T[i]); // fraction of the way through the interval -- is there abtter way to do this than linear?
if ((h <0.0)||(h>1.0)) {
printf("h not in range [0-1], h = %g   j=%d i=%d t[j] = %g T[i] = %g T[i+1] = %g\n", h, j, i, t[j], T[i], T[i+1]);
//exit(0);
}
            Lgm_QuatSlerp( Q[0], Q[1], h, q[j] );

        } else if ( i == N-2 ) {

            printf( "Lgm_QuatSquadInterp: Warning. i == %ld, can only do SLERP on [%ld,%ld] interval.\n", N-2, N-1 );
            ip1 = i+1;
            h = (t[j] - T[i])/(T[ip1] - T[i]); // fraction of the way through the interval -- is there abtter way to do this than linear?
if ((h <0.0)||(h>1.0)) {
printf("h not in range [0-1], h = %g   j=%d i=%d t[j] = %g T[i] = %g T[i+1] = %g\n", h, j, i, t[j], T[i], T[i+1]);
//exit(0);
}
            Lgm_QuatSlerp( Q[N-2], Q[N-1], h, q[j] );

        } else {

            /*
             * finally we have enough info to do a proper SQUAD.
             */

            /*
             * Compute the auxilliary control points needed for SQUAD
             */
            im1 = i-1;
            ip1 = i+1;
            ip2 = i+2;
            Lgm_QuatSquadComputeAuxPoint( Q[im1], Q[i], Q[ip1],   si ); // s_i
            Lgm_QuatSquadComputeAuxPoint( Q[i], Q[ip1], Q[ip2], sip1 ); // s_i+1

            /*
             * Compute the SQUAD
             */
            h = (t[j] - T[i])/(T[ip1] - T[i]); // fraction of the way through the interval -- is there abtter way to do this than linear?
if ((h <0.0)||(h>1.0)) {
printf("h not in range [0-1], h = %g   j=%d i=%d t[j] = %g T[i] = %g T[i+1] = %g\n", h, j, i, t[j], T[i], T[i+1]);
//exit(0);
}
            Lgm_QuatSquad( Q[i], Q[ip1], si, sip1, h, q[j] );
            

        }
    


    }



    return( 1 );

}







