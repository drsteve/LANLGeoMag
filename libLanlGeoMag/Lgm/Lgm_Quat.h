#include "Lgm_Vec.h"
#ifndef LGM_QUAT_H
#define LGM_QUAT_H

#ifndef DegPerRad
#define DegPerRad       57.29577951308232087680
#endif
#ifndef RadPerDeg
#define RadPerDeg        0.01745329251994329576
#endif


double  Lgm_NormalizeQuat( double *Q );
double  Lgm_MatrixTrace( double A[3][3] );
void    Lgm_MatrixToQuat( double A[3][3], double *Q );
void    Lgm_Quat_To_Matrix( double Q[4], double A[3][3] );
void    Lgm_QuatToAxisAngle( double *Q, double *Angle, Lgm_Vector *u );
void    Lgm_AxisAngleToQuat( Lgm_Vector *u, double Angle, double *Q );
void    Lgm_QuatRotateVector( double *Q, Lgm_Vector *v, Lgm_Vector *vp );
double  Lgm_QuatMagnitude( double *Q );
double  Lgm_QuatVecLength( double *v );
double  Lgm_QuatVecDot( double *v1, double *v2 );
void    Lgm_QuatVecZero( double *v );
void    Lgm_QuatVecSet( double *v, double x, double y, double z );
void    Lgm_QuatVecAdd( double *a, double *b, double *c );
void    Lgm_QuatVecSub( double *a, double *b, double *c);
void    Lgm_QuatVecCopy( double *v1, double *v2 );
void    Lgm_QuatVecScale( double *v, double f );
void    Lgm_QuatVecNormalize( double *v );
void    Lgm_QuatVecCross( double *a, double *b, double *result );
void    Lgm_QuatCombineQuats( double Q1[4], double Q2[4], double Q[4] );


#endif
