#include "Lgm/Lgm_MagModelInfo.h"
int Lgm_B_TS04( Lgm_Vector *v, Lgm_Vector *B, Lgm_MagModelInfo *Info ) {

    Lgm_Vector  B2;
    int		    iopt;
    double	    parmod[11], ps, X, Y, Z, Bx, By, Bz;


    parmod[1]  = Info->P; 	// Pressure in nPa
    parmod[2]  = Info->Dst;     // Dst in nPa
    parmod[3]  = Info->By; 	// IMF By in nT
    parmod[4]  = Info->Bz; 	// IMF Bz in nT
    parmod[5]  = Info->W[0];   // W1
    parmod[6]  = Info->W[1];   // W2
    parmod[7]  = Info->W[2];   // W3
    parmod[8]  = Info->W[3];   // W4
    parmod[9]  = Info->W[4];   // W5
    parmod[10] = Info->W[5];   // W6

    iopt = 0;		// this param supposedly doesnt actually do anything for this model

    ps = Info->c->psi;  // dipole tilt angle
    X = v->x; Y = v->y; Z = v->z;

    Tsyg_TS04( iopt, parmod, ps, Info->c->sin_psi, Info->c->cos_psi, X, Y, Z, &Bx, &By, &Bz );
//    t04s_( &iopt, parmod+1, &ps, &X, &Y, &Z, &Bx, &By, &Bz );
    /*
    printf("Bts04 =  (%g, %g, %g)\n", Bx, By, Bz);
    B_cdip(  v, &B2, Info );
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
                        fprintf(stderr, "Lgm_B_TS04: Unknown internal model (%d)\n", Info->InternalModel);
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

