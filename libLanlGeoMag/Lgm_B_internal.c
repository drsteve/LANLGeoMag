#include "Lgm/Lgm_MagModelInfo.h"

int Lgm_B_igrf(Lgm_Vector *v, Lgm_Vector *B, Lgm_MagModelInfo *MagInfo) {
    Lgm_B_igrf_ctrans( v, B, MagInfo->c );
    return(1);
}

int Lgm_B_cdip(Lgm_Vector *v, Lgm_Vector *B, Lgm_MagModelInfo *MagInfo) {
    Lgm_B_cdip_ctrans( v, B, MagInfo->c );
    return(1);
}

int Lgm_B_edip(Lgm_Vector *v, Lgm_Vector *B, Lgm_MagModelInfo *MagInfo) {
    Lgm_B_edip_ctrans( v, B, MagInfo->c );
    return(1);
}

int Lgm_B_JensenCain1960(Lgm_Vector *v, Lgm_Vector *B, Lgm_MagModelInfo *MagInfo) {
    Lgm_B_JensenCain1960_ctrans( v, B, MagInfo->c );
    return(1);
}
