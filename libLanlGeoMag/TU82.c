/* Copyright (c) 1999 Michael G. Henderson <mghenderson@lanl.gov>
 *
 * Permission to use, copy, modify, distribute, and sell this software and its
 * documentation for any purpose is hereby granted without fee, provided that
 * the above copyright notice appear in all copies and that both that
 * copyright notice and this permission notice appear in supporting
 * documentation.  No representations are made about the suitability of this
 * software for any purpose.  It is provided "as is" without express or
 * implied warranty.
 *
 *
 *
 *  The TU82 Magnetic field model. See following reference (note: there are some errors in the
 *  paper);
 *
 *  Tsyganenko, N. A. and A. V. Usmanov, Determination of the Magnetospheric
 *  current system parameters and development of experimental geomagnetic field
 *  models based on data from IMP and HEOS satellites, Planet. Space Sci., 10,
 *  985-998, 1982.
 *
 *
 *
 */
#include "Lgm/Lgm_MagModelInfo.h"

/*
 *
 *  User Selectable Parameters:
 *                              Info->Kp  (integer value between 0-5 selects the Kp model)
 *                              Set this param in the Lgm_MagModelInfo structure before you call routine.
 *
 *
 *  The 26 constant parameters for the 11 TU82 models -- there
 *  is one set of 26 for each Kp level. Specifically;
 *
 *                a[0][]  -> Kp =  0
 *                a[1][]  -> Kp =  0+
 *                a[2][]  -> Kp =  1-
 *                a[3][]  -> Kp =  1
 *                a[4][]  -> Kp =  1+
 *                a[5][]  -> Kp =  2-
 *                a[6][]  -> Kp =  2 
 *                a[7][]  -> Kp =  2+
 *                a[8][]  -> Kp =  3-
 *                a[9][]  -> Kp =  3,3+ 
 *                a[10][] -> Kp >= 3+
 *
 *  Each line contains the following parameters;
 *  {  a1, a2, a3, a4, a5, a6, b1, b2, b3, b4, c1, c2, c3, c4, c5, c6, c7, c8,
 *     BN, DB, B0, Dx, xN, D, Dy, rH };
 *
 */
static double Lgm_TU82_a[11][26] = {
    {   -0.0922, 1.27, 26.5, -14.2, -0.0238, 0.0294, -0.0265, -0.0362, -1.09, 0.00922,
        -29.8, 34.2, 0.157, -0.137, 0.0133, -0.0128, 1.09, 0.68, 46.6, 104.0, -12.5,        /* Kp = 0 */
        20.00, -6.82, 3.13, 10.00, 7.50 },
    {   -0.350, 1.47, 15.2, -7.24, -0.0320, 0.0445, 0.0199, -0.103, -0.441, -1.07,
        -15.4, 26.5, 0.0564, -0.0499, -0.00997, 0.0137, 0.441, 1.44, 43.5, 62.9, -21.0,     /* Kp = 0+ */
        19.44, -6.52, 3.06, 12.0, 7.37 },
    {   -0.467, 1.73, 20.3, -13.4, -0.0160, 0.0425, 0.0281, -0.126, -0.417, -1.39,
        -21.5, 34.4, 0.0928, -0.103, -0.0140, 0.0185, 0.417, 2.08, 51.5, 83.5, -24.9,       /* Kp = 1- */
        19.44, -6.82, 2.76, 12.22, 7.37 },
    {   -0.654, 1.98, 11.6, -7.96, 0.0108, 0.0696, 0.0289, -0.129, -0.116, -1.61,
        -18.2, 33.2, 0.0515, -0.0578, -0.0145, 0.0136, 0.116, 2.02, 56.8, 81.1, -26.6,      /* Kp = 1+ */
        19.44, -6.22, 3.57, 13.11, 7.11 },
    {   -0.235, 1.66, 25.0, -14.7, -0.0257, 0.0621, 0.0334, -0.140, -0.275, -1.61,
        -22.7, 32.6, 0.0673, -0.0643, -0.0167, 0.0275, 0.275, 2.37, 65.0, 121.0, -30.8,     /* Kp = 1  */
        19.44, -6.22, 3.57, 12.89, 6.46 },
    {   -0.666, 2.24, 16.4, -6.94, -0.00335, 0.0748, 0.0552, -0.183, -0.0438, -2.33,
        -18.8, 31.7, 0.0472, -0.0421, -0.0276, 0.0322, 0.0438, 2.69, 63.1, 101.0,           /* Kp = 2- */
        -29.1, 18.89, -6.07, 2.98, 11.79, 7.37 },
    {   -1.01, 2.62, 5.48, -2.37, 0.0375, 0.0867, 0.0355, -0.148, -0.233, -1.73, -19.8,
        36.1, 0.0679, -0.0835, -0.0177, 0.0066, 0.233, 1.85, 61.7, 90.1, -40.0, 19.44,      /* Kp = 2  */
        -5.63, 3.50, 15.11, 5.56 },
    {   -1.07, 2.69, 17.8, -7.32, -0.00876, 0.0247, 0.0260, -0.140, -0.137, -1.68,
        -18.5, 34.4, 0.0416, -0.0478, -0.0130, -0.00107, 0.137, 2.07, 68.8, 98.2,           /* Kp = 2+ */
        -32.3, 19.44, -5.04, 4.17, 17.11, 6.98 },
    {   -0.652, 2.48, 7.25, -1.81, 0.0139, 0.0942, 0.0332, -0.163, -0.630, -1.48,
        -18.4, 32.3, 0.0770, -0.0730, -0.0166, 0.0118, 0.630, 1.58, 61.9, 86.2, -35.2,      /* Kp = 3- */
        17.78, -4.89, 3.20, 16.00, 5.17 },
    {   -1.88, 4.02, -1.45, 1.19, 0.0115, 0.196, 0.00914, -0.156, 0.623, -2.98, -17.8,
        40.1, 0.0595, -0.0792, -0.00457, -0.0347, -0.623, 2.91, 79.1, 105.0, -46.2,         /* Kp = 3,3+ */
        17.78, -5.63, 3.72, 15.33, 4.32 },
    {   -1.17, 3.65, 0.912, 8.01, 0.0758, 0.190, 0.0615, -0.257, 0.168, -3.40, -17.8,
        32.5, 0.0692, -0.0658, -0.0308, 0.00661, -0.168, 2.86, 82.6, 124.0, -48.1,          /* Kp > 3+  */
        15.00, -4.00, 2.24, 12.00, 5.56 }
};









/*
 *   GSM coords are assumed for the input vectors.
 */
int Lgm_Brc_TU82( Lgm_Vector *v, Lgm_Vector *B, Lgm_MagModelInfo *Info ) {

    double      sin_psi, cos_psi, rho, zeta, rho_0;
    double 	    r2, z2, rr;
    Lgm_Vector 	v_sm, Bsm;
    double 	    B_rho, B_zeta;
    double      *p;
    int         indx;

    indx = (int)(Info->fKp*3.0);
    if (indx <  0 ) indx =  0;
    if (indx > 10 ) indx = 10;
    p = Lgm_TU82_a[indx];


    sin_psi = Info->c->sin_psi;
    cos_psi = Info->c->cos_psi;

    /*
     *   Compute the solar-magnetic coordinates (SM)
     */
    v_sm.x = v->x*cos_psi - v->z*sin_psi;
    v_sm.y = v->y;
    v_sm.z = v->x*sin_psi + v->z*cos_psi;



    /*
     *   Compute the cylindrical coords.
     */
    rho_0 = 4.0; 
    rho  = sqrt( v_sm.x*v_sm.x + v_sm.y*v_sm.y ) / rho_0;
    zeta = v_sm.z / rho_0;



    /*
     * Now compute Brc(rho, phi, zeta) i.e. in cylindrical coords.
     */
    r2 = rho*rho;
    z2 = zeta*zeta;
    rr = 4.0 * p[20] * pow( r2 + z2 + 4.0, -2.5 );
    B_rho  = 3.0 * rho * zeta * rr;
    B_zeta = (2.0 * z2 - r2 + 8.0) * rr;



    /*
     *   Now transform Brc(rho, phi, zeta) -> Brc(x_sm, y_sm, z_sm)
     */
    if (rho > 0.0){
	    rr =  B_rho/(rho * rho_0);
        Bsm.x = rr * v_sm.x;
        Bsm.y = rr * v_sm.y;
    }
    else{ 
        Bsm.x = 0.0;
        Bsm.y = 0.0;
    }
    Bsm.z = B_zeta;



    /*
     *    Now transform B1(x_sm, y_sm, z_sm) -> B1(x, y, z)
     */
    B->x =  Bsm.x*cos_psi + Bsm.z*sin_psi;
    B->y =  Bsm.y;
    B->z = -Bsm.x*sin_psi + Bsm.z*cos_psi;


    return(1);

}


int Lgm_Bt_TU82( Lgm_Vector *v, Lgm_Vector *B, Lgm_MagModelInfo *Info ) {

    double      S, z2, D2, z2D2, Sqrt_z2D2, xNmx, xNmxmS, F, xNmx2, xNmxmS2, G;
    double      f, q, z, zs, zr, rH, D, xN, Dy, BN, DB;
    double      *p;
    int         indx;

    indx = (int)(Info->fKp*3.0);
    if (indx <  0 ) indx =  0;
    if (indx > 10 ) indx = 10;
    p = Lgm_TU82_a[indx];


    S  = 50.0;
    BN = p[18];
    DB = p[19];

    xN = p[22];
    D  = p[23];
    Dy = p[24];
    rH = p[25];

    zs = rH * Info->c->sin_psi;
    z  = v->z - zs;
    z2 = z*z; 
    D2 = D*D;
    z2D2 = z2 + D2;
    Sqrt_z2D2 = sqrt( z2D2 );
    xNmx = xN - v->x;
    xNmxmS = xNmx - S;
    F = atan( xNmx / Sqrt_z2D2 ) - atan( xNmxmS/ Sqrt_z2D2 );


    xNmx2 = xNmx*xNmx;
    xNmxmS2 = xNmxmS*xNmxmS;
    G = log( (xNmx2 + z2D2)/(xNmxmS2 + z2D2) );

    f = 1.0/(1.0+v->y*v->y/(Dy*Dy));

    q = BN - xNmx*DB/S;

    B->x = ( z*q*F/Sqrt_z2D2 + 0.5*DB*z*G/S )*f/M_PI;
    B->y = 0.0;
    B->z = ( 0.5*q*G + DB*(1.0-Sqrt_z2D2*F/S) )*f/M_PI;

    return(1);

}


int Lgm_Bmp_TU82( Lgm_Vector *v, Lgm_Vector *B, Lgm_MagModelInfo *Info ) {

    double	sin_psi, cos_psi, e1, e2, x2, y2, z2, zc;
    double  *p;
    int         indx;

    indx = (int)(Info->fKp*3.0);
    if (indx <  0 ) indx =  0;
    if (indx > 10 ) indx = 10;
    p = Lgm_TU82_a[indx];


    sin_psi = Info->c->sin_psi;
    e1      = exp(v->x / p[21]);
    x2      = v->x * v->x;
    y2      = v->y * v->y;
    z2      = v->z * v->z;

    B->x = v->z  *(p[0] + p[1]*e1) + sin_psi*(p[2] + p[3]*e1 + p[4]*y2 + p[5]*z2);
    B->y = v->y*v->z*(p[6] + p[7]*e1) + sin_psi*(p[8] + p[9]*e1);
    B->z = p[10] + p[11]*e1 + y2*(p[12] + p[13]*e1) + z2*(p[14] + p[15]*e1) + v->z*sin_psi*(p[16] + p[17]*e1);


    return(1);
}



int Lgm_B_TU82( Lgm_Vector *v, Lgm_Vector *B, Lgm_MagModelInfo *Info ) {

    Lgm_Vector 	B1, B2, B3, B4;

    Lgm_Brc_TU82( v, &B1,  Info );
    Lgm_Bt_TU82( v, &B2,  Info );
    Lgm_Bmp_TU82( v, &B3,  Info );

    switch ( Info->InternalModel ){

        case LGM_CDIP:
                        Lgm_B_cdip( v, &B4, Info );
                        break;
        case LGM_EDIP:
                        Lgm_B_edip( v, &B4, Info );
                        break;
        case LGM_IGRF:
                        Lgm_B_igrf( v, &B4, Info );
                        break;
        default:
                        fprintf(stderr, "Lgm_B_TU82: Unknown internal model (%d)\n", Info->InternalModel );
                        break;

    }

    B->x = B1.x + B2.x + B3.x + B4.x;
    B->y = B1.y + B2.y + B3.y + B4.y;
    B->z = B1.z + B2.z + B3.z + B4.z;

    ++Info->nFunc;

    return(1);

}
