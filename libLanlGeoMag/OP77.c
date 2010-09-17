#include "Lgm/Lgm_MagModelInfo.h"
int Lgm_B_OP77( Lgm_Vector *v, Lgm_Vector *B, Lgm_MagModelInfo *m ) {

    Lgm_Vector  B1, B2, B3, u;
    double	    ps, XX[4], BF[4];

    ps = m->c->psi;  // dipole tilt angle


    /*
     *  v is in GSM, but OP77 needs it in SM
     */
    Lgm_Convert_Coords( v, &u, GSM_TO_SM, m->c );
    XX[1] = u.x; XX[2] = u.y; XX[3] = u.z;

    OlsenPfitzerStatic( XX, BF, ps*DegPerRad, m );

    /*
     *  BFv is in SM, but we need it in GSM
     */
    B1.x = BF[1]; B1.y = BF[2]; B1.z = BF[3];
    Lgm_Convert_Coords( &B1, &B2, SM_TO_GSM, m->c );

    // printf("Bop77 =  (%g, %g, %g)\n", B2.x, B2.y, B2.z);
    switch ( m->InternalModel ){

        case LGM_CDIP:
                        Lgm_B_cdip( v, &B3, m );
                        break;
        case LGM_EDIP:
                        Lgm_B_edip( v, &B3, m );
                        break;
        case LGM_IGRF:
                        Lgm_B_igrf( v, &B3, m );
                        break;
        default:
                        fprintf(stderr, "Lgm_B_OP77: Unknown internal model (%d)\n", m->InternalModel);
                        break;

    }

    B->x = B2.x + B3.x;
    B->y = B2.y + B3.y;
    B->z = B2.z + B3.z;

    ++m->nFunc;

    return(1);

}
/*
 *   $id$
 */
