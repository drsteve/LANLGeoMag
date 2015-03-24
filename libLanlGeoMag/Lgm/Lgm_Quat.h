#include "Lgm_Vec.h"
#ifndef LGM_QUAT_H
#define LGM_QUAT_H

#ifndef DegPerRad
#define DegPerRad       57.29577951308232087680
#endif
#ifndef RadPerDeg
#define RadPerDeg        0.01745329251994329576
#endif

void     Lgm_PrintQuat(double Q[4]);
double   Lgm_QuatNorm( double Q[4] );
void     Lgm_QuatConjugate( double Qin[4], double Qout[4] );
void     Lgm_QuatInverse( double Qin[4], double Qout[4] );
double   Lgm_NormalizeQuat( double Q[4] );
double   Lgm_MatrixTrace( double A[3][3] );
void     Lgm_MatrixToQuat( double A[3][3], double *Q );
void     Lgm_Quat_To_Matrix( double Q[4], double A[3][3] );
void     Lgm_AxisAngleToQuat( Lgm_Vector *u, double Angle, double *Q );
void     Lgm_QuatToAxisAngle( double *Q, double *Angle, Lgm_Vector *u );
void     Lgm_QuatRotateVector( double *Q, Lgm_Vector *v, Lgm_Vector *vp );
double   Lgm_QuatMagnitude( double *Q );
void     Lgm_QuatScale( double *Q, double t );
void     Lgm_QuatAdd( double *Q1, double *Q2, double *Q3 );
double   Lgm_QuatVecLength( double *v );
double   Lgm_QuatVecDot( double *v1, double *v2 );
void     Lgm_QuatVecZero( double *v );
void     Lgm_QuatVecSet( double *v, double x, double y, double z );
void     Lgm_QuatVecAdd( double *a, double *b, double *c );
void     Lgm_QuatVecSub( double *a, double *b, double *c);
void     Lgm_QuatVecCopy( double *v1, double *v2 );
void     Lgm_QuatVecScale( double *v, double f );
void     Lgm_QuatVecNormalize( double *v );
void     Lgm_QuatVecCross( double *a, double *b, double *c );
void     Lgm_QuatCombineQuats( double Q1[4], double Q2[4], double Q[4] );
void     Lgm_QuatMultiply( double Q1[4], double Q2[4], double Q[4] );
void     Lgm_QuatPow( double Qin[4], double t, double Qout[4] );
void     Lgm_QuatSlerp( double p[4], double q[4], double h, double Qout[4] );
void     Lgm_QuatSquad( double Q0[4], double Q1[4], double S0[4], double S1[4], double h, double Qout[4] );
void     Lgm_QuatLog( double Q[4], double logQ[4] );
void     Lgm_QuatExp( double Q[4], double expQ[4] );
void     Lgm_QuatSquadComputeAuxPoint( double Qim1[4], double Qi[4], double Qip1[4],   double si[4] );
long int Lgm_QuatFindTimeIndex( double *T, long int N, double t );
int      Lgm_QuatSquadInterp( double *T,   double *Q[4], long int N,   double *t, double *q[4], long int n );


#endif
