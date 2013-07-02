#include "Lgm/Lgm_MagModelInfo.h"
int Lgm_B_T02( Lgm_Vector *v, Lgm_Vector *B, Lgm_MagModelInfo *Info ) {

    Lgm_Vector  B2;
    int		    iopt;
    double	    parmod[11], ps, X, Y, Z, Bx, By, Bz;


    parmod[1]  = Info->P; 	// Pressure in nPa
    parmod[2]  = Info->Dst; // Dst in nPa
    parmod[3]  = Info->By; 	// IMF By in nT
    parmod[4]  = Info->Bz; 	// IMF Bz in nT
    parmod[5]  = Info->G1;  // G1
    parmod[6]  = Info->G2;  // G2

    iopt = 0;		// this param supposedly doesnt actually do anything for this model

    ps = Info->c->psi;  // dipole tilt angle
    X = v->x; Y = v->y; Z = v->z;

    Tsyg_T02( iopt, parmod, ps, Info->c->sin_psi, Info->c->cos_psi, X, Y, Z, &Bx, &By, &Bz, &Info->T01_Info );
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
                        fprintf(stderr, "Lgm_B_T02: Unknown internal model (%d)\n", Info->InternalModel);
                        break;

    }

    B->x = Bx + B2.x;
    B->y = By + B2.y;
    B->z = Bz + B2.z;
    ++Info->nFunc;

    return(1);

}


/*
 *   $id$
 */
