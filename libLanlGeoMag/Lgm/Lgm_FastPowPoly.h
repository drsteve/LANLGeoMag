#ifndef LGM_FASTPOWPOLY_H
#define LGM_FASTPOWPOLY_H

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <xmmintrin.h>

#define EXP_POLY_DEGREE 3
#define POLY0(x, c0) _mm_set1_ps(c0)
#define POLY1(x, c0, c1) _mm_add_ps(_mm_mul_ps(POLY0(x, c1), x), _mm_set1_ps(c0))
#define POLY2(x, c0, c1, c2) _mm_add_ps(_mm_mul_ps(POLY1(x, c1, c2), x), _mm_set1_ps(c0))
#define POLY3(x, c0, c1, c2, c3) _mm_add_ps(_mm_mul_ps(POLY2(x, c1, c2, c3), x), _mm_set1_ps(c0))
#define POLY4(x, c0, c1, c2, c3, c4) _mm_add_ps(_mm_mul_ps(POLY3(x, c1, c2, c3, c4), x), _mm_set1_ps(c0))
#define POLY5(x, c0, c1, c2, c3, c4, c5) _mm_add_ps(_mm_mul_ps(POLY4(x, c1, c2, c3, c4, c5), x), _mm_set1_ps(c0))

#define LOG_POLY_DEGREE 3

typedef union Lgm_myf4 {
    __m128  m;
    float   f[4];
} Lgm_myf4;


typedef union Lgm_f4 {
   int32_t  i[4];
   uint32_t u[4];
   float    f[4];
   __m128   m;
   __m128i  mi;
} Lgm_f4;


__m128  exp2f4(__m128 x);
__m128  log2f4(__m128 x);
__m128  powf4(__m128 x, __m128 y);
float   Lgm_FastPowPoly( double x, double y );
void Lgm_FastPowPoly_v( double *x, double *y, double *result );









typedef struct Lgm_FastPow {

    float         _2p23;
    unsigned int  precision;
    float         ilog2;
    unsigned int  pTable[263000];

} Lgm_FastPow;


Lgm_FastPow *Lgm_InitFastPow( );
void        *Lgm_FreeFastPow( Lgm_FastPow *fp );
void        powFastSetTable( Lgm_FastPow *fp );
float       powFastLookup( float val, Lgm_FastPow *fp );


#endif
