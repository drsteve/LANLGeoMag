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
 *  The T89 Magnetic field model. See following reference (note: there are some
 *  errors in the paper);
 *
 *      Tsyganenko, N. A., A magnetospheric magnetic field model with a warped
 *      tail current sheet, Planetary and Space Science Volume 37, Issue 1,
 *      January 1989, Pages 5-20.
 *
 *
 */

#include "Lgm/Lgm_MagModelInfo.h"

/*
 *  The 39 constant parameters for the 6 T89 Kp models -- there
 *  is one set of 39 for each Kp level. Specifically;
 *
 *                a[0][] -> Kp =     0, 0+
 *                a[1][] -> Kp = 1-, 1, 1+
 *                a[2][] -> Kp = 2-, 2, 2+
 *                a[3][] -> Kp = 3-, 3, 3+
 *                a[4][] -> Kp = 4-, 4, 4+
 *                a[5][] -> Kp = >= 5-
 *
 *  Each line contains the following parameters;
 *  {  c1,        c2,       c3,        c4,         c5,        c6,
 *     c7,        c8,       c9,       c10,        c11,       c12,
 *    c13,       c14,      c15,       c16,        c17,       c18,
 *    c19,       d_x,     a_RC,       D_0,   gamma_RC,       R_c,
 *      G,       a_T,      D_y,       x_0,        L_y,       D_x,
 *   L_RC,       L_T,  gamma_T,     delta,    gamma_1,       R_T,
 *   x_0c,     L2_xc,     D_yc  				   },
 *
 *
 *
 */
static double Lgm_T89_a[6][39] = {
    { -98.72,   -10014.0, 15.03,    76.62,    -10237.0, 1.813, 
      31.10,    -0.07464, -0.07764, 0.003303, -1.129,   0.001663,
      0.000988, 18.21,    -0.03018, -0.03829, -0.1283,  -0.001973,     /* Kp =     0, 0+ */
      0.000717, 24.74,    8.161,    2.08,     -0.8799,  9.084, 
      3.838,    13.55,    26.94,    5.745,    10.0,     13.0, 
      5.0,      6.3,      4.0,      0.01,     1.0,      30.0, 
      4.0,      50.0,     20.0 },
    { -35.64,   -12800.0, 14.37,    124.5,    -13543.0, 2.316, 
      35.64,    -0.0741,  -0.1081,  0.003924, -1.451,   0.00202,
      0.00111,  21.37,    -0.04567, -0.05382, -0.1457,  -0.002742,     /* Kp = 1-, 1, 1+ */
      0.001244, 22.33,    8.119,    1.664,    0.9324,   9.238, 
      2.426,    13.81,    28.83,    6.052,    10.0,     13.0,
      5.0,      6.3,      4.0,      0.01,     1.0,      30.0,
      4.0,      50.0,     20.0 },
    { -77.45,   -14588.0, 64.85,    123.9,    -16229.0, 2.641, 
      42.46,    -0.07611, -0.1579,  0.004078, -1.391,   0.00153,
      0.000727, 21.86,    -0.04199, -0.06523, -0.6412,  -0.000948,     /* Kp = 2-, 2, 2+ */
      0.002276, 20.90,    6.283,    1.541,    4.183,    9.609, 
      6.591,    15.08,    30.57,    7.435,    10.0,     13.0,
      5.0,      6.3,      4.0,      0.01,     1.0,      30.0,
      4.0,      50.0,     20.0 },
    { -70.12,   -16125.0, 90.71,    38.08,    -19630.0, 3.181,    
      47.50,    -0.1327,  -0.1864,  0.01382,  -1.488,   0.002962,
      0.000897, 22.74,    -0.04095, -0.09223, -1.059,   -0.001766,    /* Kp = 3-, 3, 3+ */
      0.003034, 18.64,    6.266,    0.9351,   5.389,    8.573, 
      5.935,    15.63,    31.47,    8.103,    10.0,     13.0,
      5.0,      6.3,      4.0,      0.01,     1.0,      30.0,
      4.0,      50.0,     20.0 },
    { -162.5,   -15806.0, 160.6,    5.888,   -27534.0, 3.607, 
      51.10,    -0.1006,  -0.1927,  0.03353, -1.392,   0.001594,
      0.002439, 22.41,    -0.04925, -0.1153, -1.399,   0.000716,     /* Kp = 4-, 4, 4+ */
      0.002696, 18.31,    6.196,    0.7677,  5.072,    10.06, 
      6.668,    16.11,    30.04,    8.260,    10.0,     13.0,
      5.0,      6.3,      4.0,      0.01,     1.0,      30.0,
      4.0,      50.0,     20.0 },
    { -128.4,   -16184.0, 149.1,    215.5,   -36435.0, 4.090, 
      49.09,    -0.0231,  -0.1359,  0.01989, -2.298,   0.004911,
      0.003421, 21.79,    -0.05447, -0.1149, -0.2214,  -0.01355,     /* Kp >= 5-      */
      0.001185, 19.48,    5.831,    0.3325,  6.472,    10.47, 
      9.081,    15.85,    25.27,    7.976,    10.0,     13.0,
      5.0,      6.3,      4.0,      0.01,     1.0,      30.0,
      4.0,      50.0,     20.0 }
};




int Lgm_BT_T89( Lgm_Vector *v, Lgm_Vector *B, Lgm_MagModelInfo *Info ) {


    double	   sin_psi, cos_psi, tan_psi;
    double 	   rho, zeta;
    double 	   x_sm, y_sm, z_sm;
    double 	   BT_xsm, BT_ysm, BT_zsm;
    double 	   h_1, h_T, D_T, W, z_s, z_r, xi_T, S_T, P, Q_T;
    double 	   W_x, W_y, z_sx, z_sy, D_Tx, D_Ty;
    double	   p28_2, p28_4, x_sm_2, y_sm_2, y_sm_3, y_sm_4, gg, gg2, aa, aa2, hh, bb, bb2, S_T_2;
    double	   p26_2, p29_2, p31_2, cc, cc2, pp;
    double	   zz, ss, ss12, ss32, tt, tt12, tt32, uu, uu12, uu32, nn, oonn, ooS_T, ooP, ooS_T_2, qtzr;
    double    *p;


    p = Lgm_T89_a[Info->Kp];


    sin_psi = Info->c->sin_psi;
    cos_psi = Info->c->cos_psi;
    tan_psi = Info->c->tan_psi;



    /*
     * From the GSM coordinates compute the solar-magnetic coordinates (SM)
     */
    x_sm = v->x*cos_psi - v->z*sin_psi;
    y_sm = v->y;
    z_sm = v->x*sin_psi + v->z*cos_psi;


    p28_2 = p[28]*p[28];
    p28_4 = p28_2*p28_2;
    x_sm_2 = x_sm*x_sm;
    y_sm_2 = y_sm*y_sm;
    y_sm_3 = y_sm_2*y_sm;
    y_sm_4 = y_sm_2*y_sm_2;
    gg = y_sm_4+p28_4;
    gg2 = gg*gg;
    aa = x_sm+16.0;
    aa2 = aa*aa;
    zz = x_sm+p[23];
    hh = sqrt(zz*zz+16.0);
    pp = 3.0/2.0;
    cc = x_sm-p[27];
    cc2 = cc*cc;
    p26_2 = p[26]*p[26];
    p29_2 = p[29]*p[29];
    p31_2 = p[31]*p[31];

    ss = cc2 + p29_2;
    ss12 = sqrt(ss);
    ss32 = ss12*ss;

    tt = aa2+36.0;
    tt12 = sqrt(tt);
    tt32 = tt12*tt;

    uu = x_sm_2+p31_2;
    uu12 = sqrt(uu);
    uu32 = uu12*uu;


    /*
     *   Compute the cylindrical coords.
     */
    rho = sqrt(x_sm_2 + y_sm_2);
    zeta = z_sm;



    /*
     * Now compute BT(x_sm, y_sm, z_sm) i.e. in cartesian SM coords.
     */
    h_1 = 0.5*(1.0-aa/tt12);
    h_T = 0.5*(1.0+x_sm/uu12);
    D_T = p[21] + p[33]*y_sm_2 + p[32]*h_T + p[34]*h_1;
    D_Tx = p[32]*p31_2/(2.0*uu32) -18.0*p[34]/tt32;

    nn = (1.0+y_sm_2/p26_2);
    oonn = 1.0/nn;

    D_Ty = 2.0*p[33]*y_sm;
    W = 0.5*(1.0 - cc/ss12)*oonn;
    W_x = -0.5*p29_2*oonn/ss32;
    W_y = -2.0*y_sm*W/(y_sm_2+p26_2);
    z_s = 0.5*tan_psi*(zz-hh)-p[24]*sin_psi*y_sm_4/(y_sm_4+p28_4);
    z_sx = 0.5*(1.0 - zz/hh)*tan_psi;
    z_sy = -4.0*p[24]*y_sm_3*p28_4*sin_psi/gg2;
    z_r = z_sm - z_s;
    xi_T = sqrt(z_r*z_r + D_T*D_T);
    bb = p[25]+xi_T; bb2 = bb*bb;
    S_T = sqrt(rho*rho + bb2);
    ooS_T = 1.0/S_T;
    P = S_T + p[25] + xi_T; 
    ooP = 1.0/P;
    S_T_2 = S_T*S_T;
    ooS_T_2 = 1.0/S_T_2;
    Q_T = W/(xi_T*S_T) * (p[0]*ooP + p[1]*ooS_T_2);
    qtzr = Q_T*z_r;
    BT_xsm = qtzr * x_sm;
    BT_ysm = qtzr * y_sm;
    BT_zsm = W*ooS_T * (p[0]+p[1]*(p[25]+xi_T)*ooS_T_2) 
		 + (x_sm*W_x+y_sm*W_y)*ooP*(p[0]+p[1]*ooS_T)
		 + BT_xsm*z_sx + BT_ysm*z_sy 
		 - Q_T*D_T*(x_sm*D_Tx+y_sm*D_Ty);
	



    /*
     * Now transform BT(x_sm, y_sm, z_sm) -> BT(x, y, z)
     */
    B->x =  BT_xsm*cos_psi + BT_zsm*sin_psi;
    B->y =  BT_ysm;
    B->z = -BT_xsm*sin_psi + BT_zsm*cos_psi;

    return(1);

}



int Lgm_BRC_T89( Lgm_Vector *v, Lgm_Vector *B, Lgm_MagModelInfo *Info ){


    double 	   sin_psi, cos_psi, tan_psi;
    double 	   rho, zeta;
    double 	   x, y, z, x_sm, y_sm, z_sm;
    double 	   BRC_xsm, BRC_ysm, BRC_zsm;
    double 	   h_1, h_RC, D_RC, z_s, z_r, xi_RC, S_RC, Q_RC;
    double 	   z_sx, z_sy, D_RCx;
    double 	   p28_2, p28_4, y_sm_2, y_sm_3, y_sm_4, gg, gg2, ff, ff2, S_RC2, S_RC3, S_RC_5, hh;
    double	   x_sm_2, ee, ss, ss12, ss32, tt, tt12, tt32, p30_2, qrczr;
    double    *p;


    p = Lgm_T89_a[Info->Kp];


    x = v->x; y = v->y; z = v->z;


    sin_psi = Info->c->sin_psi;
    cos_psi = Info->c->cos_psi;
    tan_psi = Info->c->tan_psi;



    /*
     * From the GSM coordinates compute the solar-magnetic coordinates (SM)
     */
    x_sm = v->x*cos_psi - v->z*sin_psi;
    y_sm = v->y;
    z_sm = v->x*sin_psi + v->z*cos_psi;



    /*
     * From SM coords compute the cylindrical coords
     */
    rho = sqrt(x_sm*x_sm + y_sm*y_sm);
    zeta = z_sm;


    p28_2 = p[28]*p[28];
    p28_4 = p28_2*p28_2;
    x_sm_2 = x_sm*x_sm;
    y_sm_2 = y_sm*y_sm;
    y_sm_3 = y_sm_2*y_sm;
    y_sm_4 = y_sm_2*y_sm_2;
    gg = y_sm_4+p28_4;
    gg2 = gg*gg;
    ee = x_sm+p[23];
    hh = sqrt(ee*ee+16.0);
    
    ff = x_sm+16.0;
    ff2 = ff*ff;

    ss = ff2+36.0;
    ss12 = sqrt(ss);
    ss32 = ss12*ss;

    p30_2 = p[30]*p[30];

    tt = x_sm_2 + p30_2;
    tt12 = sqrt(tt);
    tt32 = tt12*tt;


    /*
     * Now compute BRC(x_sm, y_sm, z_sm) i.e. in cartesian SM coords.
     */
    h_1 = 0.5*(1.0-ff/ss12);
    h_RC = 0.5*(1.0+x_sm/tt12);
    D_RC = p[21] + p[22]*h_RC + p[34]*h_1;
    D_RCx = 0.5*p[22]*p[30]*p[30]/tt32 -18.0*p[34]/ss32;
    z_s = 0.5*tan_psi*(ee-hh) -p[24]*sin_psi*y_sm_4/(y_sm_4+p28_4);
    z_sx = 0.5*(1.0 - ee/hh)*tan_psi;
    z_sy = -4.0*p[24]*y_sm_3*p28_4*sin_psi/gg2;
    z_r = z_sm - z_s;
    xi_RC = sqrt(z_r*z_r + D_RC*D_RC);
    ff = p[20]+xi_RC; ff2 = ff*ff;
    S_RC = sqrt(rho*rho + ff2);
    S_RC2 = S_RC*S_RC; S_RC3 = S_RC2*S_RC; S_RC_5 = S_RC2*S_RC3;
    Q_RC = 3.0*p[4]/(xi_RC*S_RC_5)*(p[20]+xi_RC);
    qrczr = Q_RC*z_r;
    BRC_xsm = qrczr * x_sm;
    BRC_ysm = qrczr * y_sm;
    BRC_zsm = p[4]*(2.0*ff2-rho*rho)/S_RC_5 + BRC_xsm*z_sx + BRC_ysm*z_sy - Q_RC*D_RC*x_sm*D_RCx;





    /*
     * Now transform BRC(x_sm, y_sm, z_sm) -> BRC(x, y, z)
     */
    B->x = BRC_xsm*cos_psi + BRC_zsm*sin_psi;
    B->y = BRC_ysm;
    B->z = -BRC_xsm*sin_psi + BRC_zsm*cos_psi;

    return(1);

}



int Lgm_BM_T89( Lgm_Vector *v, Lgm_Vector *B, Lgm_MagModelInfo *Info ) {

    double 	   sin_psi, cos_psi;
    double	   x, y, y2, z, z2, exod_x;
    double    *p;


    p = Lgm_T89_a[Info->Kp];

    sin_psi = Info->c->sin_psi;
    cos_psi = Info->c->cos_psi;

    x = v->x;
    y = v->y; y2 = y*y;
    z = v->z; z2 = z*z;
    
    exod_x = exp(x/p[19]);

    B->x = exod_x * (p[5] * z * cos_psi + (p[6] + p[7] * y2 + p[8] * z2) * sin_psi);
    B->y = exod_x * (p[9] * y * z * cos_psi + (p[10] * y + p[11] * y*y2 + p[12] * y * z2) *sin_psi);
    B->z = exod_x * ((p[13] + p[14] * y2 + p[15] * z2) * cos_psi + (p[16] * z + p[17] * z * y2 + p[18] * z*z2) * sin_psi);

    return(1);
}


int Lgm_BC_T89( Lgm_Vector *v, Lgm_Vector *B, Lgm_MagModelInfo *Info ) {

    double	   sin_psi, cos_psi, x, y, z;
    double 	   W_c, W_cx, W_cy, S_p, S_m;
    double 	   F_px, F_mx, F_py, F_my, F_pz, F_mz; 
    double	   p38_2, aa, aa2;
    double	   x2, y2, x2py2, zpp35, zmp35, wcx, wcy, cc, dd, ff, gg, ffcc, ggdd, zz;
    double	   psp, ss, ss12, ss32, ee, oop38_2;
    double    *p;


    p = Lgm_T89_a[Info->Kp];


    x = v->x;
    y = v->y;
    z = v->z;

    x2 = x*x;
    y2 = y*y;

    x2py2 = x2 + y2;

    sin_psi = Info->c->sin_psi;
    cos_psi = Info->c->cos_psi;

    p38_2 = p[38]*p[38];
    oop38_2 = 1.0/p38_2;
    ee = 1.0/(1.0+y2*oop38_2);
    aa = x-p[36];
    aa2 = aa*aa;

    zpp35 = z+p[35];
    zmp35 = z-p[35];

    ss = aa2+p[37];
    ss12 = sqrt(ss);
    ss32 = ss*ss12;

    /*
     * Compute BC(x, y, z)  i.e. in GSM coords.
     */
    W_c = 0.5*(1.0-aa/ss12)*ee;
    W_cx = -0.5*p[37]*ee/ss32;
    W_cy = -2.0*y*W_c/(y2+p38_2);
    S_p = sqrt(zpp35*zpp35 + x2py2);
    S_m = sqrt(zmp35*zmp35 + x2py2);

    wcx = W_c*x;
    wcy = W_c*y;

    cc = 1.0/S_p;
    dd = 1.0/S_m;
    ff = 1.0/(S_p+zpp35);
    gg = 1.0/(S_m-zmp35);

    ffcc = ff*cc;
    ggdd = gg*dd;

    zz = x*W_cx+y*W_cy;

    F_px =  wcx*ffcc;
    F_mx = -wcx*ggdd;
    F_py =  wcy*ffcc;
    F_my = -wcy*ggdd;
    F_pz =  W_c*cc + zz*ff;
    F_mz =  W_c*dd + zz*gg;

    psp = p[3]*sin_psi;

    B->x = p[2]*(F_px + F_mx) + psp*(F_px - F_mx);
    B->y = p[2]*(F_py + F_my) + psp*(F_py - F_my);
    B->z = p[2]*(F_pz + F_mz) + psp*(F_pz - F_mz);

    return(1);
}





int Lgm_B_T89( Lgm_Vector *v, Lgm_Vector *B, Lgm_MagModelInfo *Info ) {

    Lgm_Vector      B1, B2, B3, B4, B5;


    Lgm_BT_T89(  v, &B1, Info );
    Lgm_BRC_T89( v, &B2, Info );
    Lgm_BM_T89(  v, &B3, Info );
    Lgm_BC_T89(  v, &B4, Info );
    switch ( Info->InternalModel ){

        case LGM_CDIP:
                        Lgm_B_cdip( v, &B5, Info );
                        break;
        case LGM_EDIP:
                        Lgm_B_edip( v, &B5, Info );
                        break;
        case LGM_IGRF:
                        Lgm_B_igrf( v, &B5, Info );
                        break;
        default:
                        fprintf(stderr, "Lgm_B_T89: Unknown internal model (%d)\n", Info->InternalModel );
                        break;

    }





    B->x = B1.x + B2.x + B3.x + B4.x + B5.x;
    B->y = B1.y + B2.y + B3.y + B4.y + B5.y;
    B->z = B1.z + B2.z + B3.z + B4.z + B5.z;
/*
    B->x = B5.x;
    B->y = B5.y;
    B->z = B5.z;
    */


    ++Info->nFunc;

    return(1);

}


/*
 *   $id$
 */
