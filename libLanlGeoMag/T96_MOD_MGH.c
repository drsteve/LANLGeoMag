#include "Lgm/Lgm_MagModelInfo.h"
int Lgm_B_T96MOD_MGH( Lgm_Vector *v, Lgm_Vector *B, Lgm_MagModelInfo *MagInfo ) {

    Lgm_Vector  B2;
    int		    IYEAR, IDAY, IH, IM;
    double	    UTC, X, Y, Z, Bx, By, Bz, SEC, parmod[11], amdf[11];


    X = v->x; Y = v->y; Z = v->z;
    IYEAR = MagInfo->c->UTC.Year;
    IDAY  = MagInfo->c->UTC.Doy;
    UTC    = MagInfo->c->UTC.Time;
    IH    = (int)UTC;
    IM    = (int)((UTC - (double)IH)*60.0);
    SEC   = (UTC - (double)IH - (double)IM/60.0)/3600.0;


    /* 
     *  These are the extra variable params for the T96MOD
     *  We do it this way, so we can fit them more easily.
     */
    amdf[1] = MagInfo->T96MOD_V[0]; // multiplier for B_RC
    amdf[2] = MagInfo->T96MOD_V[1]; // multiplier for B_Tail2
    amdf[3] = MagInfo->T96MOD_V[2]; // multiplier for B_Tail3
    amdf[4] = MagInfo->T96MOD_V[3]; // multiplier for B_Tail3_Additional
    amdf[5] = MagInfo->T96MOD_V[4]; // DMODIF in THIN CURRENT SHEET THICKNESS (set to zero if you dont want the extra term)

    parmod[1] = MagInfo->P;
    parmod[2] = MagInfo->Dst;
    parmod[3] = MagInfo->By;
    parmod[4] = MagInfo->Bz;
    parmod[5] = 0.0;
    parmod[6] = 0.0;


    lgm_field_t96mod_mgh_( &parmod[1], &amdf[1], &IYEAR, &IDAY, &IH, &IM, &SEC, &X, &Y, &Z, &Bx, &By, &Bz );

    Lgm_B_igrf(  v, &B2, MagInfo );
    B->x = Bx + B2.x;
    B->y = By + B2.y;
    B->z = Bz + B2.z;

    ++MagInfo->nFunc;

    return(1);

}
