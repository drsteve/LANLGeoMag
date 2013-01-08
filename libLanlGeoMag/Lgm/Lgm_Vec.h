#ifndef LGM_VEC_H
#define LGM_VEC_H
typedef struct Lgm_Vector {
    double x;
    double y;
    double z;
} Lgm_Vector;

typedef struct LgmPosition {
    double x;
    double y;
    double z;
} LgmPosition;

Lgm_Vector *Lgm_CreateVector( double x, double y, double z );
void        Lgm_FreeVector( Lgm_Vector *v );
void 	    Lgm_CrossProduct(Lgm_Vector *, Lgm_Vector *, Lgm_Vector *);
double 	    Lgm_DotProduct(Lgm_Vector *, Lgm_Vector *);
double 	    Lgm_NormalizeVector(Lgm_Vector *);
void 	    Lgm_ScaleVector(Lgm_Vector *, double);
double 	    Lgm_Magnitude( Lgm_Vector *);
void 	    Lgm_ForceMagnitude(Lgm_Vector *, double);
void 	    Lgm_MatTimesVec(double A[3][3], Lgm_Vector *, Lgm_Vector *);
void        Lgm_MatTimesMat( double A[3][3], double B[3][3], double R[3][3] );
void        Lgm_VecSub(Lgm_Vector *c, Lgm_Vector *a, Lgm_Vector *b );
void        Lgm_VecAdd(Lgm_Vector *c, Lgm_Vector *a, Lgm_Vector *b );
double      Lgm_VecDiffMag(Lgm_Vector *a, Lgm_Vector *b );
void        Lgm_Transpose( double A[3][3], double B[3][3] );
void        Lgm_SphToCartCoords( double Lat, double Lon, double r, Lgm_Vector *c );
void        Lgm_CartToSphCoords( Lgm_Vector *u, double *Lat, double *Lon, double *r);
void        Lgm_SetVecVal( Lgm_Vector *u, double f );
void        Lgm_SetVecElements( Lgm_Vector *u, double x, double y, double z );
void        Lgm_SetArrVal2( double *A, double f );
void        Lgm_SetArrVal3( double *A, double f );
void        Lgm_SetArrVal4( double *A, double f );
void        Lgm_SetArrElements2( double *A, double x, double y );
void        Lgm_SetArrElements3( double *A, double x, double y, double z );
void        Lgm_SetArrElements4( double *A, double a, double b, double c, double d );
void        Lgm_VecToArr( Lgm_Vector *u, double *A );
void        Lgm_ArrToVec( double *A, Lgm_Vector *u);




#endif
