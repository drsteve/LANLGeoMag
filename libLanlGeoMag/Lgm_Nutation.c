#include "Lgm/Lgm_CTrans.h"
#include <math.h>
#include <stdio.h>

static int a[106][5] = {
    { 0,  0,  0,  0,  1}, //   1
    { 0,  0,  2, -2,  2}, //   9
    { 0,  0,  2,  0,  2}, //  31
    { 0,  0,  0,  0,  2}, //   2
    { 0,  1,  0,  0,  0}, //  10
    { 1,  0,  0,  0,  0}, //  32
    { 0,  1,  2, -2,  2}, //  11
    { 0,  0,  2,  0,  1}, //  33
    { 1,  0,  2,  0,  2}, //  34
    { 0, -1,  2, -2,  2}, //  12
    { 1,  0,  0, -2,  0}, //  35
    { 0,  0,  2, -2,  1}, //  13
    {-1,  0,  2,  0,  2}, //  36
    { 1,  0,  0,  0,  1}, //  38
    { 0,  0,  0,  2,  0}, //  37
    {-1,  0,  2,  2,  2}, //  40
    {-1,  0,  0,  0,  1}, //  39
    { 1,  0,  2,  0,  1}, //  41
    { 2,  0,  0, -2,  0}, //  14
    {-2,  0,  2,  0,  1}, //   3
    { 0,  0,  2,  2,  2}, //  42
    { 2,  0,  2,  0,  2}, //  45
    { 2,  0,  0,  0,  0}, //  43
    { 1,  0,  2, -2,  2}, //  44
    { 0,  0,  2,  0,  0}, //  46
    { 0,  0,  2, -2,  0}, //  15
    {-1,  0,  2,  0,  1}, //  47
    { 0,  2,  0,  0,  0}, //  16
    { 0,  2,  2, -2,  2}, //  18
    {-1,  0,  0,  2,  1}, //  48
    { 0,  1,  0,  0,  1}, //  17
    { 1,  0,  0, -2,  1}, //  49
    { 0, -1,  0,  0,  1}, //  19
    { 2,  0, -2,  0,  0}, //   4
    {-1,  0,  2,  2,  1}, //  50
    { 1,  0,  2,  2,  2}, //  54
    { 1,  1,  0, -2,  0}, //  51
    { 0,  1,  2,  0,  2}, //  52
    { 0, -1,  2,  0,  2}, //  53
    { 0,  0,  2,  2,  1}, //  58
    {-2,  0,  0,  2,  1}, //  20
    { 1,  0,  0,  2,  0}, //  55
    { 2,  0,  2, -2,  2}, //  56
    { 0,  0,  0,  2,  1}, //  57
    { 1,  0,  2, -2,  1}, //  59
    { 0, -1,  2, -2,  1}, //  21
    { 0,  0,  0, -2,  1}, //  60
    { 1, -1,  0,  0,  0}, //  61
    { 2,  0,  2,  0,  1}, //  62
    { 2,  0,  0, -2,  1}, //  22
    { 0,  1,  2, -2,  1}, //  23
    { 1,  0,  0, -1,  0}, //  24
    { 0,  1,  0, -2,  0}, //  63
    { 1,  0, -2,  0,  0}, //  64
    { 0,  0,  0,  1,  0}, //  65
    {-2,  0,  2,  0,  2}, //   5
    { 1, -1,  0, -1,  0}, //   6
    { 1,  1,  0,  0,  0}, //  66
    { 1,  0,  2,  0,  0}, //  67
    { 1, -1,  2,  0,  2}, //  68
    {-1, -1,  2,  2,  2}, //  69
    { 3,  0,  2,  0,  2}, //  71
    { 0, -1,  2,  2,  2}, //  72
    { 0, -2,  2, -2,  1}, //   7
    {-2,  0,  0,  0,  1}, //  70
    { 1,  1,  2,  0,  2}, //  73
    {-1,  0,  2, -2,  1}, //  74
    { 2,  0,  0,  0,  1}, //  75
    { 1,  0,  0,  0,  2}, //  76
    { 3,  0,  0,  0,  0}, //  77
    { 0,  0,  2,  1,  2}, //  78
    {-1,  0,  2,  4,  2}, //  82
    { 2,  0, -2,  0,  1}, //   8
    { 2,  1,  0, -2,  0}, //  25
    { 0,  0, -2,  2,  1}, //  26
    { 0,  1, -2,  2,  0}, //  27
    { 0,  1,  0,  0,  2}, //  28
    {-1,  0,  0,  1,  1}, //  29
    { 0,  1,  2, -2,  0}, //  30
    {-1,  0,  0,  0,  2}, //  79
    { 1,  0,  0, -4,  0}, //  80
    {-2,  0,  2,  2,  2}, //  81
    { 2,  0,  0, -4,  0}, //  83
    { 1,  1,  2, -2,  2}, //  84
    { 1,  0,  2,  2,  1}, //  85
    {-2,  0,  2,  4,  2}, //  86
    {-1,  0,  4,  0,  2}, //  87
    { 1, -1,  0, -2,  0}, //  88
    { 2,  0,  2, -2,  1}, //  89
    { 2,  0,  2,  2,  2}, //  90
    { 1,  0,  0,  2,  1}, //  91
    { 0,  0,  4, -2,  2}, //  92
    { 3,  0,  2, -2,  2}, //  93
    { 1,  0,  2, -2,  0}, //  94
    { 0,  1,  2,  0,  1}, //  95
    {-1, -1,  0,  2,  1}, //  96
    { 0,  0, -2,  0,  1}, //  97
    { 0,  0,  2, -1,  2}, //  98
    { 0,  1,  0,  2,  0}, //  99
    { 1,  0, -2, -2,  0}, // 100
    { 0, -1,  2,  0,  1}, // 101
    { 1,  1,  0, -2,  1}, // 102
    { 1,  0, -2,  2,  0}, // 103
    { 2,  0,  0,  2,  0}, // 104
    { 0,  0,  2,  4,  2}, // 105
    { 0,  1,  0,  1,  0}  // 106
};

// 2nd and 4th columns are multiplies by 10 so we can store as long ints.
static long int A[106][4] = {
    { -171996, -1742,  92025,    89 }, //   1
    {  -13187,   -16,   5736,   -31 }, //   9
    {   -2274,    -2,    977,    -5 }, //  31
    {    2062,     2,   -895,     5 }, //   2
    {    1426,   -34,     54,    -1 }, //  10
    {     712,     1,     -7,     0 }, //  32
    {    -517,    12,    224,    -6 }, //  11
    {    -386,    -4,    200,     0 }, //  33
    {    -301,     0,    129,    -1 }, //  34
    {     217,    -5,    -95,     3 }, //  12
    {    -158,     0,     -1,     0 }, //  35
    {     129,     1,    -70,     0 }, //  13
    {     123,     0,    -53,     0 }, //  36
    {      63,     1,    -33,     0 }, //  38
    {      63,     0,     -2,     0 }, //  37
    {     -59,     0,     26,     0 }, //  40
    {     -58,    -1,     32,     0 }, //  39
    {     -51,     0,     27,     0 }, //  41
    {      48,     0,      1,     0 }, //  14
    {      46,     0,    -24,     0 }, //   3
    {     -38,     0,     16,     0 }, //  42
    {     -31,     0,     13,     0 }, //  45
    {      29,     0,     -1,     0 }, //  43
    {      29,     0,    -12,     0 }, //  44
    {      26,     0,     -1,     0 }, //  46
    {     -22,     0,      0,     0 }, //  15
    {      21,     0,    -10,     0 }, //  47
    {      17,    -1,      0,     0 }, //  16
    {     -16,     1,      7,     0 }, //  18
    {      16,     0,     -8,     0 }, //  48
    {     -15,     0,      9,     0 }, //  17
    {     -13,     0,      7,     0 }, //  49
    {     -12,     0,      6,     0 }, //  19
    {      11,     0,      0,     0 }, //   4
    {     -10,     0,      5,     0 }, //  50
    {      -8,     0,      3,     0 }, //  54
    {      -7,     0,      0,     0 }, //  51
    {       7,     0,     -3,     0 }, //  52
    {      -7,     0,      3,     0 }, //  53
    {      -7,     0,      3,     0 }, //  58
    {      -6,     0,      3,     0 }, //  20
    {       6,     0,      0,     0 }, //  55
    {       6,     0,     -3,     0 }, //  56
    {      -6,     0,      3,     0 }, //  57
    {       6,     0,     -3,     0 }, //  59
    {      -5,     0,      3,     0 }, //  21
    {      -5,     0,      3,     0 }, //  60
    {       5,     0,      0,     0 }, //  61
    {      -5,     0,      3,     0 }, //  62
    {       4,     0,     -2,     0 }, //  22
    {       4,     0,     -2,     0 }, //  23
    {      -4,     0,      0,     0 }, //  24
    {      -4,     0,      0,     0 }, //  63
    {       4,     0,      0,     0 }, //  64
    {      -4,     0,      0,     0 }, //  65
    {      -3,     0,      1,     0 }, //   5
    {      -3,     0,      0,     0 }, //   6
    {      -3,     0,      0,     0 }, //  66
    {       3,     0,      0,     0 }, //  67
    {      -3,     0,      1,     0 }, //  68
    {      -3,     0,      1,     0 }, //  69
    {      -3,     0,      1,     0 }, //  71
    {      -3,     0,      1,     0 }, //  72
    {      -2,     0,      1,     0 }, //   7
    {      -2,     0,      1,     0 }, //  70
    {       2,     0,     -1,     0 }, //  73
    {      -2,     0,      1,     0 }, //  74
    {       2,     0,     -1,     0 }, //  75
    {      -2,     0,      1,     0 }, //  76
    {       2,     0,      0,     0 }, //  77
    {       2,     0,     -1,     0 }, //  78
    {      -2,     0,      1,     0 }, //  82
    {       1,     0,      0,     0 }, //   8
    {       1,     0,      0,     0 }, //  25
    {       1,     0,      0,     0 }, //  26
    {      -1,     0,      0,     0 }, //  27
    {       1,     0,      0,     0 }, //  28
    {       1,     0,      0,     0 }, //  29
    {      -1,     0,      0,     0 }, //  30
    {       1,     0,     -1,     0 }, //  79
    {      -1,     0,      0,     0 }, //  80
    {       1,     0,     -1,     0 }, //  81
    {      -1,     0,      0,     0 }, //  83
    {       1,     0,     -1,     0 }, //  84
    {      -1,     0,      1,     0 }, //  85
    {      -1,     0,      1,     0 }, //  86
    {       1,     0,      0,     0 }, //  87
    {       1,     0,      0,     0 }, //  88
    {       1,     0,     -1,     0 }, //  89
    {      -1,     0,      0,     0 }, //  90
    {      -1,     0,      0,     0 }, //  91
    {       1,     0,      0,     0 }, //  92
    {       1,     0,      0,     0 }, //  93
    {      -1,     0,      0,     0 }, //  94
    {       1,     0,      0,     0 }, //  95
    {       1,     0,      0,     0 }, //  96
    {      -1,     0,      0,     0 }, //  97
    {      -1,     0,      0,     0 }, //  98
    {      -1,     0,      0,     0 }, //  99
    {      -1,     0,      0,     0 }, // 100
    {      -1,     0,      0,     0 }, // 101
    {      -1,     0,      0,     0 }, // 102
    {      -1,     0,      0,     0 }, // 103
    {       1,     0,      0,     0 }, // 104
    {      -1,     0,      0,     0 }, // 105
    {       1,     0,      0,     0 }  // 106
};


void Lgm_Nutation( double T_TT, double nTerms, double *dPsi, double *dEps ) {

    int     i;
    double  Mmoon, Msun, uMmoon, Dsun, Nmoon;
    double  ap, PsiSum, EpsSum, T2_TT, T3_TT;

    T2_TT = T_TT*T_TT; T3_TT = T2_TT*T_TT;
    Mmoon  = fmod( 134.96340251 + (1325.0*360.0 + 198.8675605)*T_TT + 0.0088553*T2_TT + 1.4343e-5*T3_TT, 360.0)*RadPerDeg;
    Msun   = fmod( 357.52910918 + (99.0*360.0   + 359.0502911)*T_TT - 0.0001537*T2_TT +    3.8e-8*T3_TT, 360.0)*RadPerDeg;
    uMmoon = fmod(  93.27209062 + (1342.0*360.0 +  82.0174577)*T_TT - 0.0035420*T2_TT -   2.88e-7*T3_TT, 360.0)*RadPerDeg;
    Dsun   = fmod( 297.85019547 + (1236.0*360.0 + 307.1114469)*T_TT - 0.0017696*T2_TT +  1.831e-6*T3_TT, 360.0)*RadPerDeg;
    Nmoon  = fmod( 125.04455501 - (5.0*360.0    + 134.1361851)*T_TT + 0.0020756*T2_TT +  2.139e-6*T3_TT, 360.0)*RadPerDeg;

    if (nTerms > 106) nTerms = 106;
    for (PsiSum = EpsSum = 0.0, i=0; i<nTerms; i++){

        ap = a[i][0]*Mmoon + a[i][1]*Msun + a[i][2]*uMmoon + a[i][3]*Dsun + a[i][4]*Nmoon;
        PsiSum += (A[i][0] + 0.1*A[i][1]*T_TT)*sin( ap );
        EpsSum += (A[i][2] + 0.1*A[i][3]*T_TT)*cos( ap );

    }

    *dPsi = PsiSum*0.0001;
    *dEps = EpsSum*0.0001;

}
