#include "Lgm/Lgm_MagModelInfo.h"
int Lgm_B_OP88( Lgm_Vector *v, Lgm_Vector *B, Lgm_MagModelInfo *Info ) {

    Lgm_Vector       B2;
    double	         DEN, VEL, DST, Bmag, X, Y, Z, Bx, By, Bz;


    X = v->x; Y = v->y; Z = v->z;
    DEN = Info->Den;
    VEL = Info->V;
    DST = Info->Dst;
    Lgm_OP88_BDYN( DEN, VEL, DST, X, Y, Z, &Bx, &By, &Bz );

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
    Bmag = Lgm_Magnitude( B );
/*
    B->x = B2.x;
    B->y = B2.y;
    B->z = B2.z;
    printf("Bigrf =  (%f, %f, %f)\n", B2.x, B2.y, B2.z);
*/





    ++Info->nFunc;

    return(1);

}

