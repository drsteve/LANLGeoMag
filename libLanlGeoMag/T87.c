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
 *  The T87 long Magnetic field model. See following reference (note: there are some errors in the
 *  paper);
 *
 *  Tsyganenko, N. A., Global quantitative models of the geomagnetic field
 *     in the cislunar magnetosphere for different disturbance levels,
 *     Planet. Space Sci., 35, 1347-1358, 1987.
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
 *  The 36 constant parameters for the 6 T87-long models -- there
 *  is one set of 36 for each Kp level. Specifically;
 *
 *                a[0][] -> Kp =     0, 0+
 *                a[1][] -> Kp = 1-, 1, 1+
 *                a[2][] -> Kp = 2-, 2, 2+
 *                a[3][] -> Kp = 3-, 3, 3+
 *                a[4][] -> Kp = 4-, 4, 4+
 *                a[5][] -> Kp = >= 5-
 *
 *  Each line contains the following parameters;
 *  {  a1,        a2,       a3,        a4,         a5,        a6, 
 *     b1,        b2,       b3,        b4,         b5,        b6, 
 *     c1,        c2,       c3,        c4,         c5,        c6, 
 *     c7,        c8,       c9,        c10,        B0,        B1, 
 *     B2         Brc,      Rrc        xN,         D,         dy, 
 *     RH,        dx1,      dx2,       x1,         x2,        RT         },
 *
 *
 *  Note that it is correct that b6 for Kp = 4-, 4, 4+ is positive (there's
 *  a typo in the original paper).
 *
 */
static double Lgm_T87_a[6][36] = {
    {  -0.09673,  -10.63,   1.210,     34.57,      -0.04502,  -0.06553, 	
       -0.02952,  0.3852,   -0.03665,  -2.084,     0.001795,  0.00638, 
       -23.49,    0.06082,  0.01642,   -0.02137,   32.21,     -0.04373, 	/* Kp =     0, 0+ */
       -0.02311,  -0.2832,  -0.002303, -0.000631,  -6.397,    -967.0, 
       -8650.0,   -20.55,   5.180,     -2.796,     2.715,     13.58,
       8.038,     29.21,    14.605,    4.0,        5.0,       30.0       },
    {  -0.48500,  -12.84,   1.856,     40.06,      -0.0294,   -0.09071, 
       -0.02993,  0.5465,   -0.04928,  -2.453,     0.001587,  0.007402, 
       -29.41,    0.08101,  0.02322,   -0.10910,   40.75,     -0.07995, 	/* Kp = 1-, 1, 1+ */
       -0.03859,  -0.2755,  -0.002759, -0.000408,  -6.189,    -957.8, 
       -7246.0,   -25.51,   5.207,     -4.184,     2.641,     16.56, 
       7.795,     29.36,    14.68,     4.0,        5.0,       30.0       },
    {  -1.13200,  -18.05,   2.625,     48.55,      -0.004868, -0.10870, 
       -0.03824,  0.8514,   -0.05220,  -2.881,     -0.000295, 0.009055, 
       -29.48,    0.06394,  0.03864,   -0.22880,   41.77,     -0.05849, 	/* Kp = 2-, 2, 2+ */
       -0.06443,  -0.4683,  0.001222,  -0.000519,  -3.696,    -991.1, 
       -6955.0,   -31.43,   4.878,     -3.151,     3.277,     19.19, 
       7.248,     28.99,    14.495,    4.0,        5.0,       30.0       },
    {  -1.00300,  -16.98,   3.140,     52.81,      -0.08625,  -0.14780, 
       -0.03501,  0.5500,   -0.07778,  -2.970,     0.002086,  0.01275, 
       -26.79,    0.06328,  0.03622,   0.08345,    39.72,     -0.06009, 	/* Kp = 3-, 3, 3+ */
       -0.07825,  -0.9698,  0.000178,  -0.000573,  -0.9328,   -872.5, 
       -5851.0,   -39.68,   4.902,     -3.848,     2.790,     20.91, 
       6.193,     26.81,    13.405,    4.0,        5.0,       30.0       },
    {  -1.53900,  -14.29,   3.479,     53.36,      -0.004201, -0.20430, 
       -0.03932,  0.6409,   -0.10580,  -3.221,     -0.00114,  0.02166, 
       -30.43,    0.04049,  0.05464,   0.008884,   42.00,     -0.01035, 	/* Kp = 4-, 4, 4+ */
       -0.10530,  -1.630,   0.003802,  -0.001029,  4.204,     -665.6, 
       -1011.0,   -43.49,   4.514,     -2.948,     2.990,     21.59, 
       6.005,     22.00,    11.00,     4.0,        5.0,       30.0       },
    {  -2.581,    -7.726,   5.045,     53.31,      0.02262,   -0.1972, 
       -0.01981,  0.4280,   -0.1055,   -5.075,     0.002762,  0.03277, 
       -27.35,    0.04986,  0.06119,   -0.1211,    47.48,     -0.0502, 		/* Kp >= 5-      */
       -0.1477,   0.838,    -0.01008,  -0.0057,    9.231,     -674.3, 
       -900.0,    -74.43,   4.658,     -3.245,     3.39,      21.80, 
       5.620,     25.17,    12.585,    4.0,        5.0,       30.0       }
};





/*
 *   Routines compute the B1, B2, and B3 terms in the T87 long model.
 *   GSM coords are assumed.
 */
int Lgm_B1_T87( Lgm_Vector *v, Lgm_Vector *B, Lgm_MagModelInfo *Info ) {

    double      sin_psi, cos_psi, rho, zeta;
    double 	    r2, z2, rr;
    Lgm_Vector 	v_sm, Bsm;
    double 	    B_rho, B_zeta;
    double      *p;

    p = Lgm_T87_a[Info->Kp];


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
    rho  = sqrt( v_sm.x*v_sm.x + v_sm.y*v_sm.y ) / p[26];
    zeta = v_sm.z / p[26];



    /*
     * Now compute B1(rho, phi, zeta) i.e. in cylindrical coords.
     */
    r2 = rho*rho;
    z2 = zeta*zeta;
    rr = 4.0 * p[25] * pow( r2 + z2 + 4.0, -2.5 );
    B_rho  = 3.0 * rho * zeta * rr;
    B_zeta = (2.0 * z2 - r2 + 8.0) * rr;



    /*
     *   Now transform B1(rho, phi, zeta) -> B1(x_sm, y_sm, z_sm)
     */
    if (rho > 0.0){
	rr =  B_rho/(rho * p[26]);
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


int Lgm_B2_T87( Lgm_Vector *v, Lgm_Vector *B, Lgm_MagModelInfo *Info ) {

    double 	   xsi_1, xsi_2, xsi_N;
    double 	   z_r, z_p, z_m;
    double 	   beta, beta_p, beta_m;
    double 	   x22mb2, x22mbp2, x22mbm2, xNmx2tgamma_2, xNmx2tgamma_2p, xNmx2tgamma_2m;
    double 	   xsi_12, xsi_22, xsi_N2, D2, beta2, beta_p2, beta_m2;
    double 	   gamma_22, gamma_2p2, gamma_2m2;
    double 	   xN_m_x1, xN_m_x2, xN_m_x1_2, xN_m_x2_2, x2b, x2bp, x2bm;
    double 	   log_xN_m_x1_2, log_xN_m_x2_2, log_x2b, log_x2bp, log_x2bm;
    double 	   gamma_1, gamma_1p, gamma_1m, gamma_2, gamma_2p, gamma_2m;
    double 	   P_1, P_1p, P_1m, P_2, P_2p, P_2m;
    double 	   S_0, S_0p, S_0m, S_1, S_1p, S_1m, S_2, S_2p, S_2m;
    double 	   G_0, G_1, G_1p, G_1m, G_2, G_2p, G_2m;
    double 	   f;
    double    *p;

    p = Lgm_T87_a[Info->Kp];


    xsi_1  = p[33] - v->x;
    xsi_2  = p[34] - v->x;
    xsi_N  = p[27] - v->x;
    xsi_12 = xsi_1*xsi_1;
    xsi_22 = xsi_2*xsi_2;
    xsi_N2 = xsi_N*xsi_N;

    z_r = v->z - p[30] * Info->c->sin_psi;
    z_p = v->z - p[35];
    z_m = v->z + p[35];


    D2      = p[28]*p[28];
    beta2   = z_r*z_r + D2;
    beta    = sqrt( beta2 );
    beta_p2 = z_p*z_p + D2;
    beta_p  = sqrt( beta_p2 );
    beta_m2 = z_m*z_m + D2;
    beta_m  = sqrt( beta_m2 );


    gamma_1   = xsi_12 + beta2;
    gamma_1p  = xsi_12 + beta_p2;
    gamma_1m  = xsi_12 + beta_m2;
    gamma_2   = xsi_22 + beta2;
    gamma_22  = gamma_2*gamma_2;
    gamma_2p  = xsi_22 + beta_p2;
    gamma_2p2 = gamma_2p*gamma_2p;
    gamma_2m  = xsi_22 + beta_m2;
    gamma_2m2 = gamma_2m*gamma_2m;

    xN_m_x1    = p[27] - p[33];
    xN_m_x2    = p[27] - p[34];
    xN_m_x1_2  = xN_m_x1*xN_m_x1;
    xN_m_x2_2  = xN_m_x2*xN_m_x2;
    x2b        = xsi_N2 + beta2;
    x2bp       = xsi_N2 + beta_p2;
    x2bm       = xsi_N2 + beta_m2;

    log_xN_m_x1_2 = log( xN_m_x1_2 );
    log_xN_m_x2_2 = log( xN_m_x2_2 );
    log_x2b       = log( x2b );
    log_x2bp      = log( x2bp );
    log_x2bm      = log( x2bm );

    P_1  = 0.5 / gamma_1   * (log_xN_m_x1_2 - log_x2b  );
    P_1p = 0.5 / gamma_1p  * (log_xN_m_x1_2 - log_x2bp );
    P_1m = 0.5 / gamma_1m  * (log_xN_m_x1_2 - log_x2bm );

    P_2  = 1.0 / gamma_22  * (log_xN_m_x2_2 - log_x2b  );
    P_2p = 1.0 / gamma_2p2 * (log_xN_m_x2_2 - log_x2bp );
    P_2m = 1.0 / gamma_2m2 * (log_xN_m_x2_2 - log_x2bm );

    S_0  = (M_PI_2 + atan(xsi_N / beta  )) / beta;
    S_0p = (M_PI_2 + atan(xsi_N / beta_p)) / beta_p;
    S_0m = (M_PI_2 + atan(xsi_N / beta_m)) / beta_m;

    S_1  = P_1  - (xsi_1 / gamma_1 ) * S_0;
    S_1p = P_1p - (xsi_1 / gamma_1p) * S_0p;
    S_1m = P_1m - (xsi_1 / gamma_1m) * S_0m;

    x22mb2  = xsi_22 - beta2;
    x22mbp2 = xsi_22 - beta_p2;
    x22mbm2 = xsi_22 - beta_m2;
    xNmx2tgamma_2   = xN_m_x2 * gamma_2;
    xNmx2tgamma_2p  = xN_m_x2 * gamma_2p;
    xNmx2tgamma_2m  = xN_m_x2 * gamma_2m;
    S_2  = -xsi_2 * P_2  - 1.0 / xNmx2tgamma_2  + x22mb2  / gamma_22  * S_0;
    S_2p = -xsi_2 * P_2p - 1.0 / xNmx2tgamma_2p + x22mbp2 / gamma_2p2 * S_0p;
    S_2m = -xsi_2 * P_2m - 1.0 / xNmx2tgamma_2m + x22mbm2 / gamma_2m2 * S_0m;

    G_0  = 0.5 * log_x2b - .25 * ( log_x2bp + log_x2bm );

    G_1  = beta2   * S_0  / gamma_1  + xsi_1 * P_1;
    G_1p = beta_p2 * S_0p / gamma_1p + xsi_1 * P_1p;
    G_1m = beta_m2 * S_0m / gamma_1m + xsi_1 * P_1m;

    G_2  = -0.5 * x22mb2  * P_2    - 2.0 * beta2   * xsi_2 * S_0 /gamma_22    - xsi_2/xNmx2tgamma_2;
    G_2p = -0.5 * x22mbp2 * P_2p   - 2.0 * beta_p2 * xsi_2 * S_0p/gamma_2p2   - xsi_2/xNmx2tgamma_2p;
    G_2m = -0.5 * x22mbm2 * P_2m   - 2.0 * beta_m2 * xsi_2 * S_0m/gamma_2m2   - xsi_2/xNmx2tgamma_2m;
	
    f = 1.0 / (M_PI * (1.0 + (v->y/p[29])*(v->y/p[29])));



    /*
     *   Compute B2(x, y, z)  i.e. in GSM coords.
     */
    B->x = f * (p[22] * ( z_r * S_0 - 0.5 * (z_p * S_0p + z_m 
		* S_0m)) + p[23] * (z_r * S_1 - 0.5 * (z_p * S_1p 
		+ z_m * S_1m)) + p[24] * (z_r * S_2 - 0.5 * (z_p 
		* S_2p + z_m * S_2m)));
    B->y = 0.0;
    B->z = f * (p[22] * G_0 + p[23] * (G_1 - 0.5 * (G_1p + G_1m)) 
		+ p[24] * (G_2 - 0.5 * (G_2p + G_2m)));

    return(1);

}


int Lgm_B3_T87( Lgm_Vector *v, Lgm_Vector *B, Lgm_MagModelInfo *Info ) {

    double	sin_psi, cos_psi, e1, e2, x2, y2, z2, zc;
    double  *p;

    p = Lgm_T87_a[Info->Kp];
    

    sin_psi = Info->c->sin_psi;
    cos_psi = Info->c->cos_psi;
    e1      = exp(v->x / p[31]);
    e2      = exp(v->x / p[32]);
    x2      = v->x * v->x;
    y2      = v->y * v->y;
    z2      = v->z * v->z;
    zc      = v->z * cos_psi;

    B->x = e1 * (p[0]*zc + p[1]*sin_psi)
		+ e2 * (p[2]*zc + (p[3] + p[4]*y2 + p[5]*z2) * sin_psi);

    B->y = e1 * (p[6]*v->y*zc + p[7]*v->y*sin_psi) 
		+ e2 * (p[8]*v->y*zc + (p[9]*v->y + p[10]*v->y*y2 + p[11]*v->y*z2) * sin_psi);

    B->z = e1 * ((p[12] + p[13]*y2 + p[14]*z2) * cos_psi + p[15]*v->z*sin_psi) 
		+ e2 * ((p[16] + p[17]*y2 + p[18]*z2) * cos_psi + (p[19]*v->z + p[20]*v->z*y2 + p[21]*v->z*z2) * sin_psi);

    return(1);
}



int Lgm_B_T87( Lgm_Vector *v, Lgm_Vector *B, Lgm_MagModelInfo *Info ) {

    Lgm_Vector 	B1, B2, B3, B4;

    Lgm_B1_T87( v, &B1,  Info );
    Lgm_B2_T87( v, &B2,  Info );
    Lgm_B3_T87( v, &B3,  Info );

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
                        fprintf(stderr, "Lgm_B_T87: Unknown internal model (%d)\n", Info->InternalModel );
                        break;

    }

    B->x = B1.x + B2.x + B3.x + B4.x;
    B->y = B1.y + B2.y + B3.y + B4.y;
    B->z = B1.z + B2.z + B3.z + B4.z;

    ++Info->nFunc;

    return(1);

}


/*
 *   $Id$
 */
