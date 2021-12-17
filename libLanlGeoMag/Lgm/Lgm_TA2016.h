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
    int     AAlloced;
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

    Lgm_CTrans *c;


} LgmTA16_Info;







/*
 *  Function declarations
 */
void Lgm_Init_TA16( LgmTA16_Info *t );

int Lgm_Copy_TA16_Info( LgmTsyg2007_Info *t, LgmTsyg2007_Info *s );

int Lgm_SetCoeffs_TA16( long int Date, double UTC, LgmTA16_Info *t );


#endif

