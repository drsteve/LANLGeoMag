#ifndef LGM_TA2016_H
#define LGM_TA2016_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/*
 * Define a structure to hold all of the info needed in TA16
 */
typedef struct LgmTA16_Info {

    int     SetupDone;
    int     ArraysAlloced;
    int     lenA;
    int     lenL;
    double  *A;
    double  *XX;
    double  *YY;
    double  *ZZ;
    double  *ST;
    double  *RHO;
    double  *ZSP;
    double  *ZCP;
    double  *RHBR;

    double Pdyn;
    double SymHc_avg;
    double Xind_avg;
    double By_avg;

    Lgm_CTrans *c;


} LgmTA16_Info;







/*
 *  Function declarations
 */
void Lgm_Init_TA16( LgmTA16_Info *ta );

void Lgm_DeAllocate_TA16( LgmTA16_Info *t );

int Lgm_Copy_TA16_Info( LgmTA16_Info *targ, LgmTA16_Info *src );

int Lgm_GetData_TA16( LgmTA16_Info *ta );

int Lgm_SetCoeffs_TA16( long int Date, double UTC, LgmTA16_Info *ta );

int TA2016_SetGrid( LgmTA16_Info *t );

int TA2016( Lgm_Vector *posGSM, double *PARMOD, Lgm_CTrans *ctrans, Lgm_Vector *BvecGSM, LgmTA16_Info *Info );

#endif

