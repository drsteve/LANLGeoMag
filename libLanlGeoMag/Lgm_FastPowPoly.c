#include "Lgm/Lgm_FastPowPoly.h"
#include <stdlib.h>

__m128 exp2f4(__m128 x) {

   __m128i ipart;
   __m128 fpart, expipart, expfpart;

   x = _mm_min_ps(x, _mm_set1_ps( 129.00000f));
   x = _mm_max_ps(x, _mm_set1_ps(-126.99999f));

   /* ipart = int(x - 0.5) */
   ipart = _mm_cvtps_epi32(_mm_sub_ps(x, _mm_set1_ps(0.5f)));

   /* fpart = x - ipart */
   fpart = _mm_sub_ps(x, _mm_cvtepi32_ps(ipart));

   /* expipart = (float) (1 << ipart) */
   expipart = _mm_castsi128_ps(_mm_slli_epi32(_mm_add_epi32(ipart, _mm_set1_epi32(127)), 23));

   /* minimax polynomial fit of 2**x, in range [-0.5, 0.5[ */
#if EXP_POLY_DEGREE == 5
   expfpart = POLY5(fpart, 9.9999994e-1f, 6.9315308e-1f, 2.4015361e-1f, 5.5826318e-2f, 8.9893397e-3f, 1.8775767e-3f);
#elif EXP_POLY_DEGREE == 4
   expfpart = POLY4(fpart, 1.0000026f, 6.9300383e-1f, 2.4144275e-1f, 5.2011464e-2f, 1.3534167e-2f);
#elif EXP_POLY_DEGREE == 3
   expfpart = POLY3(fpart, 9.9992520e-1f, 6.9583356e-1f, 2.2606716e-1f, 7.8024521e-2f);
#elif EXP_POLY_DEGREE == 2
   expfpart = POLY2(fpart, 1.0017247f, 6.5763628e-1f, 3.3718944e-1f);
//#else
//#error
#endif

   return _mm_mul_ps(expipart, expfpart);
}


__m128 log2f4(__m128 x) {
   __m128i exp = _mm_set1_epi32(0x7F800000);
   __m128i mant = _mm_set1_epi32(0x007FFFFF);

   __m128 one = _mm_set1_ps( 1.0f);

   __m128i i = _mm_castps_si128(x);

   __m128 e = _mm_cvtepi32_ps(_mm_sub_epi32(_mm_srli_epi32(_mm_and_si128(i, exp), 23), _mm_set1_epi32(127)));

   __m128 m = _mm_or_ps(_mm_castsi128_ps(_mm_and_si128(i, mant)), one);

   __m128 p;

   /* Minimax polynomial fit of log2(x)/(x - 1), for x in range [1, 2[ */
#if LOG_POLY_DEGREE == 6
   p = POLY5( m, 3.1157899f, -3.3241990f, 2.5988452f, -1.2315303f,  3.1821337e-1f, -3.4436006e-2f);
#elif LOG_POLY_DEGREE == 5
   p = POLY4(m, 2.8882704548164776201f, -2.52074962577807006663f, 1.48116647521213171641f, -0.465725644288844778798f, 0.0596515482674574969533f);
#elif LOG_POLY_DEGREE == 4
   p = POLY3(m, 2.61761038894603480148f, -1.75647175389045657003f, 0.688243882994381274313f, -0.107254423828329604454f);
#elif LOG_POLY_DEGREE == 3
   p = POLY2(m, 2.28330284476918490682f, -1.04913055217340124191f, 0.204446009836232697516f);
//#else
//#error
#endif

   /* This effectively increases the polynomial degree by one, but ensures that log2(1) == 0*/
   p = _mm_mul_ps(p, _mm_sub_ps(m, one));

   return _mm_add_ps(p, e);
}


/*
 * Finally, the fast pow() function
 */
__m128 powf4(__m128 x, __m128 y) {
   return exp2f4(_mm_mul_ps(log2f4(x), y));
}

float Lgm_FastPowPoly( double x, double y ) {

    float       result;
    Lgm_myf4    d;
    __m128      a, b;

    a = (__m128){ (float)x, 0.0f, 0.0f, 0.0f };
    b = (__m128){ (float)y, 0.0f, 0.0f, 0.0f };
    d.m = powf4( a, b );
    result = d.f[0];

    return(result);

}

void Lgm_FastPowPoly_v( double *x, double *y, double *result ) {

    Lgm_myf4    d;
    __m128      a, b;

    a = (__m128){ (float)x[0], (float)x[1], (float)x[2], (float)x[3] };
    b = (__m128){ (float)y[0], (float)y[1], (float)y[2], (float)y[3] };
    d.m = powf4( a, b );
    result[0] = (double)d.f[0];
    result[1] = (double)d.f[1];
    result[2] = (double)d.f[2];
    result[3] = (double)d.f[3];

}





Lgm_FastPow *Lgm_InitFastPow( ){

    Lgm_FastPow *f = (Lgm_FastPow *)calloc( 1, sizeof(Lgm_FastPow) );

    f->_2p23     = 8388608.0f;
    f->precision = 12;                  // >= 0 and <= 18
    f->ilog2     = 3.32192809488736;    // use for 10^x
    powFastSetTable( f );
    

    return( f );

}

void *Lgm_FreeFastPow( Lgm_FastPow *f ){
    free( f );
}



/**
 * Initialize powFast lookup table.
 */
void powFastSetTable( Lgm_FastPow *ff ) {

   /* step along table elements and x-axis positions */
   float zeroToOne = 1.0f / ((float)(1 << ff->precision) * 2.0f);       /* A */
   int   i;                                                             /* B */
   for( i = 0;  i < (1 << ff->precision);  ++i ) {                      /* C */
      /* make y-axis value for table element */
      float f = ((float)pow( 2.0f, zeroToOne ) - 1.0f) * ff->_2p23;
      ff->pTable[i] = (unsigned int)( f < ff->_2p23 ? f : (ff->_2p23 - 1.0f) );
      zeroToOne += 1.0f / (float)(1 << ff->precision);
   }                                                                    /* D */

}

/**
 * Get pow (fast!).
 *
 * @val        power to raise radix to
 */
float powFastLookup( float val, Lgm_FastPow *ff ) {

   /* build float bits */
   int i = (int)( (val * (ff->_2p23 * ff->ilog2)) + (127.0f * ff->_2p23) );

   /* replace mantissa with lookup */
   int it = (i & 0xFF800000) | ff->pTable[(i & 0x7FFFFF) >>        /* E */
      (23 - ff->precision)];                                       /* F */

   /* convert bits to float */
   return *(float*)( &it );
printf("*(float*)( &it ) = %g\n", *(float*)( &it ));

}

