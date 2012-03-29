/* 
 */
#include "Lgm/Lgm_MagModelInfo.h"

/*
 *
 *
 */
int Lgm_B_Dungey(Lgm_Vector *v, Lgm_Vector *B, Lgm_MagModelInfo *Info) {

    double      x_sm, y_sm, z_sm;
    double      r, theta, phi;
    Lgm_Vector  Bsm;
    double      B_r, B_theta, B_phi;
    double      M, r2, r3, dB;
    double      cp, sp, ct, st;

    M  = Info->M_Dungey;  // Chen and Shultyz have 30500.0 nT
    dB = Info->dB_Dungey;  // Constant Delta-B in SM z--direction

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
    B_r      = -2.0*M*ct/r3;
    B_phi    =  0.0;
    B_theta  = -1.0*M*st/r3;


    /*
     *   Transform (B_r, B_theta, B_phi) -> (B_xsm, B_ysm, B_zsm)  (still SM)
     *   (Here we also add the dB value).
     */
    Bsm.x = B_r*st*cp + B_theta*ct*cp;
    Bsm.y = B_r*st*sp + B_theta*ct*sp;
    Bsm.z = B_r*ct    - B_theta*st + dB;
    



    /*
     *  Transform (B_xsm, B_ysm, B_zsm) -> (B_x, B_y, B_z) i.e. trans. to GSM
     */
    B->x =  Bsm.x*Info->c->cos_psi + Bsm.z*Info->c->sin_psi;
    B->y =  Bsm.y;
    B->z = -Bsm.x*Info->c->sin_psi + Bsm.z*Info->c->cos_psi;

    return(1);


}

