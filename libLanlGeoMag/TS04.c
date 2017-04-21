#include "Lgm/Lgm_MagModelInfo.h"
int Lgm_B_TS04( Lgm_Vector *v, Lgm_Vector *B, Lgm_MagModelInfo *Info ) {

    Lgm_Vector       B2;
    int		         iopt;
    double	         parmod[11], ps, Bmag, X, Y, Z, Bx, By, Bz;
    double           FEXT, FINT, BBX, BBY, BBZ, OIMFX, OIMFY, OIMFZ;


    parmod[1]  = Info->P; 	    // Pressure in nPa
    parmod[2]  = Info->Dst;     // Dst in nPa
    parmod[3]  = Info->By; 	    // IMF By in nT
    parmod[4]  = Info->Bz; 	    // IMF Bz in nT
    parmod[5]  = Info->W[0];    // W1
    parmod[6]  = Info->W[1];    // W2
    parmod[7]  = Info->W[2];    // W3
    parmod[8]  = Info->W[3];    // W4
    parmod[9]  = Info->W[4];    // W5
    parmod[10] = Info->W[5];    // W6

    iopt = 0;		            // this param supposedly doesnt actually do anything for this model
    ps = Info->c->psi;          // dipole tilt angle
    X = v->x; Y = v->y; Z = v->z;

    /*
     * Notes about Tsyg_TS04():
     *
     *        (1) If the parmod values are not set appropriately, the model can
     *            return bad values. For example, if IMF By and Bz are zero,
     *            the model will return a zero value in the IMF region. This
     *            will screw up FL tracing routines. Also, the model only uses
     *            IMF By and Bz - so if the IMF is purely radial (i.e. all in
     *            Bx), the same problem will occur. I dont know if this ever
     *            occurs in the QinDenton Files -- but its something to watch
     *            out for.
     *
     *        (2) The model separates the final calculations into 3 regions:
     *            (a) INSIDE MAGNETOPAUSE; (b) OUTSIDE MAGNETOPAUSE, BUT INSIDE
     *            a narrow BOUNDARY REGION over which the field transitions to
     *            the IMF field; (c) the IMF region. Across the boundary
     *            region, the model blends the inside field with the IMF field
     *            using the FINT and FEXT factors. At the inner edge of the
     *            boundary, FINT=1, FEXT=0. At the outer edge of the boundary,
     *            FINT=0, FEXT=1. To do the blending properly, it needs to have
     *            the full fields (i.e. B + B_dipole). So in the final calcs,
     *            the model does this sort of thing;
     *
     *                  if ( We are INSIDE ) {
     *
     *                      Bfinal = Bmodel;
     *
     *                  } else if ( We are in BOUNDARY REGION ) {
     *
     *                      Bfinal = (Bmodel + Bdipole) * Fint + Bimf * Fext - Bdipole;
     *
     *                  } else if ( We are in IMF REGION ) {
     *
     *                      Bfinal = Bimf - Bdipole;
     *
     *                  }
     *
     *            Note that the blending is done with the full field (Bmodel +
     *            Bdipole). The dipole field is then subtracted before the
     *            values are returned because Tsyg_TS04() is just supposed to
     *            return the "external field" (i.e. B due to everything but
     *            the internal dipole field).
     *
     *            Why is all of this important? There are two issues: 
     *
     *                  (a) Our choice of internal model is almost certainly
     *                      different from Tsyganenko's and that introduces an
     *                      inconsistency.
     *
     *                  (b) If you somehow run the model with IMF By and Bz =
     *                      0, the differences between Tsyganenko's dipole and ours
     *                      will produce strange residual values in the BOUNDARY
     *                      and IMF regions. You will get very peculiar FLs out there.
     *
     *            To avoid problems, make sure the input params arfe set to
     *            something sane. In addition, I have changed the code so that
     *            the final blending is done here in this code where we have
     *            access to our dipole model.
     *
     */
    Tsyg_TS04( iopt, parmod, ps, Info->c->sin_psi, Info->c->cos_psi, X, Y, Z, &Bx, &By, &Bz, &Info->TS04_Info );

    switch ( Info->InternalModel ){

        case LGM_CDIP:
                        Lgm_B_cdip( v, &B2, Info );
                        break;
        case LGM_EDIP:
                        Lgm_B_edip( v, &B2, Info );
                        break;
        case LGM_IGRF:
                        Lgm_B_igrf( v, &B2, Info );
                        break;
        default:
                        fprintf(stderr, "Lgm_B_TS04: Unknown internal model (%d)\n", Info->InternalModel);
                        break;

    }

    /*
     * Reconstruct the field the way it was done originally, but use our own
     * dipole model instead.
     */
    FINT  = Info->TS04_Info.FINT;
    FEXT  = Info->TS04_Info.FEXT;
    BBX   = Info->TS04_Info.BBX;
    BBY   = Info->TS04_Info.BBY;
    BBZ   = Info->TS04_Info.BBZ;
    OIMFX = Info->TS04_Info.OIMFX;
    OIMFY = Info->TS04_Info.OIMFY;
    OIMFZ = Info->TS04_Info.OIMFZ;

    if ( Info->TS04_Info.Region == LGM_TS04_MAGNETOSPHERE ) {

        /*
         * Same formula as in the EXTERN model. Except we add dipole.
         */
        B->x = BBX  + B2.x;
        B->y = BBY  + B2.y;
        B->z = BBZ  + B2.z;

    } else if ( Info->TS04_Info.Region == LGM_TS04_BOUNDARY ) {

        /*
         * Same formula as in the EXTERN model. Except we dont subtract dipole.
         */
        B->x = (BBX + B2.x)*FINT + OIMFX*FEXT;
        B->y = (BBY + B2.y)*FINT + OIMFY*FEXT;
        B->z = (BBZ + B2.z)*FINT + OIMFZ*FEXT;

    } else {

        /*
         * Same formula as in the EXTERN model. Except we dont subtract dipole.
         */
        B->x = OIMFX;
        B->y = OIMFY;
        B->z = OIMFZ;

    }


    ++Info->nFunc;

    return(1);

}

