#include "Lgm/Lgm_MagModelInfo.h"
int Lgm_B_TS07( Lgm_Vector *v, Lgm_Vector *B, Lgm_MagModelInfo *Info ) {

    Lgm_Vector       B2;
    int		         iopt;
    double	         parmod[11], ps, X, Y, Z, Bx, By, Bz;


    parmod[1]  = Info->TS07_Info.Pdyn; 	// Pressure in nPa
//printf("parmod[1] = %g\n", parmod[1]);

    iopt = 0;		// this param supposedly doesnt actually do anything for this model

    ps = Info->c->psi;  // dipole tilt angle
    X = v->x; Y = v->y; Z = v->z;

    Tsyg_TS07( iopt, parmod, ps, Info->c->sin_psi, Info->c->cos_psi, X, Y, Z, &Bx, &By, &Bz, &Info->TS07_Info );
    /*
    printf("Bts07 =  (%g, %g, %g)\n", Bx, By, Bz);
    Lgm_B_cdip(  v, &B2, Info );
    printf("Bcdip =  (%f, %f, %f)\n", B2.x, B2.y, B2.z);
    */
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
                        fprintf(stderr, "Lgm_B_TS07: Unknown internal model (%d)\n", Info->InternalModel);
                        break;

    }

    B->x = Bx + B2.x;
    B->y = By + B2.y;
    B->z = Bz + B2.z;
/*
    B->x = B2.x;
    B->y = B2.y;
    B->z = B2.z;
    printf("Bigrf =  (%f, %f, %f)\n", B2.x, B2.y, B2.z);
*/


    ++Info->nFunc;

    return(1);

}

