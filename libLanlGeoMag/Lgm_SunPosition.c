#include <stdio.h>
#include <math.h>
#define INCLUDE_PERTURBATIONS 1

double Frac( double x ) {
/*
    x -= floor(x);
    if (x<0.0) x += 1.0;
*/
    x = x-(double)((int)x); if (x<0.0) x += 1.0;
    return( x );
}

void AddThe( double c1, double s1, double c2, double s2, double *c, double *s ){
    *c = c1*c2 - s1*s2;
    *s = s1*c2 + c1*s2;
}

void Term( double T, int i1, int i, int it, double dlc, double dls, double drc, double drs, double dbc, double dbs, 
            double c3[], double s3[], double c[], double s[],  double *dl, double *dr, double *db ) {

    double u, v;

    if (it == 0){
        AddThe( c3[i1], s3[i1], c[i], s[i], &u, &v );
    } else {
        u *= T;
        v *= T;
    }
    *dl += dlc*u + dls*v;
    *dr += drc*u + drs*v;
    *db += dbc*u + dbs*v;
}


void Lgm_SunPosition( double T, double *l, double *r, double *b ) {

    int     i;
    double  p2;
    double  m2, m3, m4, m5, m6;
    double  d, a, uu, c3[9], s3[9], c[9], s[9];
    double  dl, dr, db;


    p2 = 2.0*M_PI;
    dl = dr = db = 0.0;
    m2 = p2*Frac( 0.1387306 + 162.5485917*T );
    m3 = p2*Frac( 0.9931266 + 99.9973604*T );
    m4 = p2*Frac( 0.0543250 + 53.1666028*T );
    m5 = p2*Frac( 0.0551750 + 8.4293972*T );
    m6 = p2*Frac( 0.8816500 + 3.3938722*T );
    d  = p2*Frac( 0.8274 + 1236.8531*T);
    a  = p2*Frac( 0.3749 + 1325.5524*T );
    uu = p2*Frac( 0.2591 + 1342.2278*T );

    c3[1] = 1.0; s3[1] = 0.0;
    c3[2] = cos(m3); s3[2] = sin(m3);
    c3[0] = c3[2]; s3[0] = -s3[2];
    for (i=3; i<9; i++ ) AddThe( c3[i-1], s3[i-1], c3[2], s3[2], &c3[i], &s3[i] );




    if ( INCLUDE_PERTURBATIONS ) {

        // perturbations for Venus
        c[8] = 1.0; s[8] = 0.0;
        c[7] = cos(m2); s[7] = -sin(m2);
        for (i=7; i>=3; i--) AddThe( c[i], s[i], c[7], s[7], &c[i-1], &s[i-1] );
        Term( T, 2, 8, 0, -0.22, 6892.76, -16707.37, -0.54,  0.00,  0.00, c3, s3, c, s, &dl, &dr, &db );
        Term( T, 2, 8, 1, -0.06,  -17.35,     42.04, -0.15,  0.00,  0.00, c3, s3, c, s, &dl, &dr, &db );
        Term( T, 2, 8, 2, -0.01,   -0.05,      0.13, -0.02,  0.00,  0.00, c3, s3, c, s, &dl, &dr, &db );
        Term( T, 3, 8, 0,  0.00,   71.98,   -139.57,  0.00,  0.00,  0.00, c3, s3, c, s, &dl, &dr, &db );
        Term( T, 3, 8, 1,  0.00,   -0.36,      0.70,  0.00,  0.00,  0.00, c3, s3, c, s, &dl, &dr, &db );
        Term( T, 4, 8, 0,  0.00,    1.04,     -1.75,  0.00,  0.00,  0.00, c3, s3, c, s, &dl, &dr, &db );
        Term( T, 1, 7, 0,  0.03,   -0.07,     -0.16, -0.07,  0.02, -0.02, c3, s3, c, s, &dl, &dr, &db );
        Term( T, 2, 7, 0,  2.35,   -4.23,     -4.75, -2.64,  0.00,  0.00, c3, s3, c, s, &dl, &dr, &db );
        Term( T, 2, 6, 0, -0.10,    0.06,      0.12,  0.20,  0.02,  0.00, c3, s3, c, s, &dl, &dr, &db );
        Term( T, 3, 7, 0, -0.06,   -0.03,      0.20, -0.01,  0.01, -0.09, c3, s3, c, s, &dl, &dr, &db );
        Term( T, 3, 6, 0, -4.70,    2.90,      8.28, 13.42,  0.01, -0.01, c3, s3, c, s, &dl, &dr, &db );
        Term( T, 4, 6, 0,  1.80,   -1.74,     -1.44, -1.57,  0.04, -0.06, c3, s3, c, s, &dl, &dr, &db );
        Term( T, 4, 5, 0, -0.67,    0.03,      0.11,  2.43,  0.01,  0.00, c3, s3, c, s, &dl, &dr, &db );
        Term( T, 5, 6, 0,  0.03,   -0.03,      0.10,  0.09,  0.01, -0.01, c3, s3, c, s, &dl, &dr, &db );
        Term( T, 5, 5, 0,  1.51,   -0.40,     -0.88, -3.36,  0.18, -0.10, c3, s3, c, s, &dl, &dr, &db );
        Term( T, 5, 4, 0, -0.19,   -0.09,     -0.38,  0.77,  0.00,  0.00, c3, s3, c, s, &dl, &dr, &db );
        Term( T, 6, 5, 0,  0.76,   -0.68,      0.30,  0.37,  0.01,  0.00, c3, s3, c, s, &dl, &dr, &db );
        Term( T, 6, 4, 0, -0.14,   -0.04,     -0.11,  0.43, -0.03,  0.00, c3, s3, c, s, &dl, &dr, &db );
        Term( T, 6, 3, 0, -0.05,   -0.07,     -0.31,  0.21,  0.00,  0.00, c3, s3, c, s, &dl, &dr, &db );
        Term( T, 7, 4, 0,  0.15,   -0.04,     -0.06, -0.21,  0.01,  0.00, c3, s3, c, s, &dl, &dr, &db );
        Term( T, 7, 3, 0, -0.03,   -0.03,     -0.09,  0.09, -0.01,  0.00, c3, s3, c, s, &dl, &dr, &db );
        Term( T, 7, 2, 0,  0.00,   -0.04,     -0.18,  0.02,  0.00,  0.00, c3, s3, c, s, &dl, &dr, &db );
        Term( T, 8, 3, 0, -0.12,   -0.03,     -0.08,  0.31, -0.02, -0.01, c3, s3, c, s, &dl, &dr, &db );

        // perturbations for Mars
        c[7] = cos(m4); s[7] = -sin(m4);
        for (i=7; i>=1; i--) AddThe( c[i], s[i], c[7], s[7], &c[i-1], &s[i-1] );
        Term( T, 2, 7, 0, -0.22,  0.17, -0.21, -0.27, 0.00, 0.00, c3, s3, c, s, &dl, &dr, &db );
        Term( T, 2, 6, 0, -1.66,  0.62,  0.16,  0.28, 0.00, 0.00, c3, s3, c, s, &dl, &dr, &db );
        Term( T, 3, 6, 0,  1.96,  0.57, -1.32,  4.55, 0.00, 0.01, c3, s3, c, s, &dl, &dr, &db );
        Term( T, 3, 5, 0,  0.40,  0.15, -0.17,  0.46, 0.00, 0.00, c3, s3, c, s, &dl, &dr, &db );
        Term( T, 3, 4, 0,  0.53,  0.26,  0.09, -0.22, 0.00, 0.00, c3, s3, c, s, &dl, &dr, &db );
        Term( T, 4, 5, 0,  0.05,  0.12, -0.35,  0.15, 0.00, 0.00, c3, s3, c, s, &dl, &dr, &db );
        Term( T, 4, 4, 0, -0.13, -0.48,  1.06, -0.29, 0.01, 0.00, c3, s3, c, s, &dl, &dr, &db );
        Term( T, 4, 3, 0, -0.04, -0.20,  0.20, -0.04, 0.00, 0.00, c3, s3, c, s, &dl, &dr, &db );
        Term( T, 5, 4, 0,  0.00, -0.03,  0.10,  0.04, 0.00, 0.00, c3, s3, c, s, &dl, &dr, &db );
        Term( T, 5, 3, 0,  0.05, -0.07,  0.20,  0.14, 0.00, 0.00, c3, s3, c, s, &dl, &dr, &db );
        Term( T, 5, 2, 0, -0.10,  0.11, -0.23, -0.22, 0.00, 0.00, c3, s3, c, s, &dl, &dr, &db );
        Term( T, 6, 1, 0, -0.05,  0.00,  0.01, -0.14, 0.00, 0.00, c3, s3, c, s, &dl, &dr, &db );
        Term( T, 6, 0, 0,  0.05,  0.01, -0.02,  0.10, 0.00, 0.00, c3, s3, c, s, &dl, &dr, &db );

        // perturbations for Jupiter
        c[7] = cos(m5); s[7] = -sin(m5);
        for (i=7; i>=5; i--) AddThe( c[i], s[i], c[7], s[7], &c[i-1], &s[i-1] );
        Term( T, 0, 7, 0,  0.01,  0.07,  0.18,  -0.02,   0.00, -0.02, c3, s3, c, s, &dl, &dr, &db );
        Term( T, 1, 7, 0, -0.31,  2.58,  0.52,   0.34,   0.02,  0.00, c3, s3, c, s, &dl, &dr, &db );
        Term( T, 2, 7, 0, -7.21, -0.06,  0.13, -16.27,   0.00, -0.02, c3, s3, c, s, &dl, &dr, &db );
        Term( T, 2, 6, 0, -0.54, -1.52,  3.09,  -1.12,   0.01, -0.17, c3, s3, c, s, &dl, &dr, &db );
        Term( T, 2, 5, 0, -0.03, -0.21,  0.38,  -0.06,   0.00, -0.02, c3, s3, c, s, &dl, &dr, &db );
        Term( T, 3, 7, 0, -0.16,  0.05, -0.18,  -0.31,   0.01,  0.00, c3, s3, c, s, &dl, &dr, &db );
        Term( T, 3, 6, 0,  0.14, -2.73,  9.23,   0.48,   0.00,  0.00, c3, s3, c, s, &dl, &dr, &db );
        Term( T, 3, 5, 0,  0.07, -0.55,  1.83,   0.25,   0.01,  0.00, c3, s3, c, s, &dl, &dr, &db );
        Term( T, 3, 4, 0,  0.02, -0.08,  0.25,   0.06,   0.00,  0.00, c3, s3, c, s, &dl, &dr, &db );
        Term( T, 4, 6, 0,  0.01, -0.07,  0.16,   0.04,   0.00,  0.00, c3, s3, c, s, &dl, &dr, &db );
        Term( T, 4, 5, 0, -0.16, -0.03,  0.08,  -0.64,   0.00,  0.00, c3, s3, c, s, &dl, &dr, &db );
        Term( T, 4, 4, 0, -0.04, -0.01,  0.03,  -0.17,   0.00,  0.00, c3, s3, c, s, &dl, &dr, &db );

        // perturbations for Saturn
        c[7] = cos(m6); s[7] = -sin(m6);
        AddThe( c[7], s[7], c[7], s[7], &c[6], &s[6] );
        Term( T, 1, 7, 0,  0.00,  0.32,  0.01,   0.00,   0.00,  0.00, c3, s3, c, s, &dl, &dr, &db );
        Term( T, 2, 7, 0, -0.08, -0.41,  0.97,  -0.18,   0.00, -0.01, c3, s3, c, s, &dl, &dr, &db );
        Term( T, 2, 6, 0,  0.04,  0.10, -0.23,   0.10,   0.00,  0.00, c3, s3, c, s, &dl, &dr, &db );
        Term( T, 3, 6, 0,  0.04,  0.10, -0.35,   0.13,   0.00,  0.00, c3, s3, c, s, &dl, &dr, &db );

        // perturbations for Moon
        dl += 6.45*sin(d)  - 0.42*sin(d-a) + 0.18*sin(d+a) + 0.17*sin(d-m3) - 0.06*sin(d+m3);
        dr += 30.76*cos(d) - 3.06*cos(d-a) + 0.85*cos(d+a) - 0.58*cos(d+m3) + 0.57*cos(d-m3);
        db += 0.576*sin(uu);

    }
    

    dl += 6.40*sin(p2*(0.6983 + 0.0561*T)) + 1.87*sin(p2*(0.5764 + 0.4174*T))
            + 0.27*sin(p2*(0.4189 + 0.3306*T)) + 0.20*sin(p2*(0.3581 + 2.4814*T));
    *l = M_PI/180.0 * 360.0*Frac( 0.7859453 + m3/p2 + ((6191.2 + 1.1*T)*T + dl)/1296e3);
    *r = 1.0001398 - 7e-7*T + 1e-6*dr;
    *b = M_PI/180.0 * db/3600.0;


}

