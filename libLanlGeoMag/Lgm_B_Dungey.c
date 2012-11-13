/*! \file Lgm_B_Dungey.c
 *
 *  \brief The "Dungey" magnetic field model (i.e. dipole + constant IMF Bz value.)
 *
 *
 *  \author S. K. Morley
 *  \date   2011
 *
 *
 *
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

    Info->M_Dungey = 30500.0;
    Info->dB_Dungey = -14.474;

    M  = Info->M_Dungey;  // Chen and Shultyz have 30500.0 nT
    dB = Info->dB_Dungey;  // Constant Delta-B in SM z--direction

    /*
     *  compute SM coords from GSM coords
     */
    x_sm = v->x*Info->c->cos_psi - v->z*Info->c->sin_psi;
    y_sm = v->y;
    z_sm = v->x*Info->c->sin_psi + v->z*Info->c->cos_psi;

    //Call the CDIP model with M_Dungey
    Info->c->M_cd = M;
    Lgm_B_cdip( v, B, Info);

    //Convert CDIP value to SM
    Lgm_Convert_Coords( B, &Bsm, GSM_TO_SM, Info->c );

    /*
     *   Here we add the dB value.
     */
    Bsm.z = Bsm.z + dB;
    
    /*
     *  Transform (B_xsm, B_ysm, B_zsm) -> (B_x, B_y, B_z) i.e. trans. to GSM
     */
    B->x =  Bsm.x*Info->c->cos_psi + Bsm.z*Info->c->sin_psi;
    B->y =  Bsm.y;
    B->z = -Bsm.x*Info->c->sin_psi + Bsm.z*Info->c->cos_psi;

    return(1);
}

