/* Copyright (c) 2008 Michael G. Henderson <mghenderson@lanl.gov>
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
 *  The "Simplified Mead" field derived in Shultz and Lanzerotti (see equations 1.44-1.45))
 *
 *
 */
#include "Lgm/Lgm_MagModelInfo.h"

/*
 *
 *  User Selectable Parameters:
 *          SW pressure
 *
 */
int Lgm_SimplifiedMead(Lgm_Vector *v, Lgm_Vector *B, Lgm_MagModelInfo *Info) {

    double      x_sm, y_sm, z_sm;
    double      r, theta, phi;
    Lgm_Vector  Bsm;
    double      B_r, B_theta, B_phi;
    double      B0, B1, B2, Bt, b, b2, b3, b4;
    double      M, r2, r3;
    double      cp, sp, ct, st;

    M = Info->c->M_cd;
    B0 = Info->B0 * M;
    B1 = Info->B1 * M;               	// See page 30 of Schultz and Lanzerotti, [1974]
    B2 = sqrt(3.0) * Info->B2 * M;	// See page 30 of Schultz and Lanzerotti, [1974]

    /*
     *  compute SM coords from GSM coords
     */
    x_sm = v->x*Info->c->cos_psi - v->z*Info->c->sin_psi;
    y_sm = v->y;
    z_sm = v->x*Info->c->sin_psi + v->z*Info->c->cos_psi;


    /*
     *  convert x_sm, y_sm, and, z_sm to spherical coords.
     */
    r     = sqrt(x_sm*x_sm + y_sm*y_sm + z_sm*z_sm);
    phi   = atan2(y_sm, x_sm);
    theta = acos(z_sm / r);

    /*
     *  compute b, the standoff distance
     *
     *  At subsolar point, Ps = Pb 
     *  where Ps = nmv^2
     *  and Pb = B^2/2mu0
     *
     *  For this simple model, B at subsolar (i.e. at r=b) is:
     *      B = Btheta = -(M + B1 + B2)/b^3
     *  Balancing Pb with Psw gives:
     *  
     *          (  (M+B1+B2)^2  ) (1/6)
     *      b = ( --------------)
     *          (  2 mu0 Ps     )
     *
     *  If M+B1+B2 is in nT and Ps is in nPa, then:
     *
     *          (  3.9788736e-4 * (M+B1+B2)^2  ) (1/6)
     *      b = ( ---------------------------- )
     *          (               Ps             )
     *  
     *   Note that eqn 1.43 in Schultz and Lanzerotti, [1974], is for the full Mead field.
     *   we can probably use either one though???
     */
    Bt = B0+B1+B2;
    b = pow( 3.9788736e-4 * Bt*Bt / Info->P, 1.0/6.0 );
    b = 1.068*pow( 3.9788736e-4 * B0*B0/Info->P, 1.0/6.0); //eqn (1.43)
//printf("B0, B1, B2, Psw, Bt, b = %g, %g, %g, %g, %g, %g\n", B0, B1, B2, Info->P, Bt, b);
    b2 = b*b;
    b3 = b*b2;
    b4 = b2*b2;
    


    /*
     *  compute centered dipole field in spherical coords
     *  i.e. (B_r, B_theta, B_phi). In Schultz and Lanzerotti, [1974],
     *  Phi is measured from midnight. Here its from noon.
     */
    r2       = r*r;
    r3       = r*r2;
    cp       = cos( phi );
    sp       = sin( phi );
    ct       = cos( theta );
    st       = sin( theta );
    B_r      = -2.0*B0*ct/r3 + 2.0*B1*ct/b3 + 2.0*B2*r*ct*st*cp/b4;
    B_phi    = -B2*r*ct*sp/b4;
    B_theta  = -1.0*B0*st/r3 - B1*st/b3 - B2*r*(2.0*st*st-1.0)*cp/b4;
//printf("r = %g\n", r);


    /*
     *   Transform (B_r, B_theta, B_phi) -> (B_xsm, B_ysm, B_zsm)  (still SM)
     */
    Bsm.x = B_r*st*cp + B_theta*ct*cp;
    Bsm.y = B_r*st*sp + B_theta*ct*sp;
    Bsm.z = B_r*ct    - B_theta*st;


    /*
     *  Transform (B_xsm, B_ysm, B_zsm) -> (B_x, B_y, B_z) i.e. trans. to GSM
     */
    B->x =  Bsm.x*Info->c->cos_psi + Bsm.z*Info->c->sin_psi;
    B->y =  Bsm.y;
    B->z = -Bsm.x*Info->c->sin_psi + Bsm.z*Info->c->cos_psi;

    return(1);


}



/*
 *   $Id: Lgm_SimplifiedMead.c 45 2010-10-01 20:43:29Z mgh $
 */

