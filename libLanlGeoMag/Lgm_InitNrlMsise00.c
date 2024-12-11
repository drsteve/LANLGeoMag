#include "Lgm/Lgm_DynamicMemory.h"
#include "Lgm/Lgm_NrlMsise00.h"
#include "Lgm/Lgm_NrlMsise00_Data.h"


Lgm_Msis00Info *InitMsis00( ) {

    int             i, j;
    Lgm_Msis00Info *p;

    p = (Lgm_Msis00Info *)calloc( 1, sizeof(*p) );
    

    /*
     * Set Name, DATE andm TIME identifiers...
     */
//    strcpy( p->NAME,   NAME );
//    strcpy( p->ISDATE, ISDATE );
//    strcpy( p->ISTIME, ISTIME );

    p->IMR = 0;


    /*
     * Copy over PT array
     */
    //LGM_ARRAY_1D( p->PT, 150+2, double );
    LGM_ARRAY_1D( p->PT, 152, double );
    for ( i=1; i<=150; i++ ) p->PT[i] = PT[i];


    /*
     * Set up the PD array in correct order.
     */
    //LGM_ARRAY_2D( p->PD, 150+2, 9+2, double );
    LGM_ARRAY_2D( p->PD, 152, 11, double );
    for ( i=1; i<=150; i++ ) {
        p->PD[i][1] = PA[i];
//printf("p->PD[%d][1] = %g  PA[%d] = %g\n", i, p->PD[i][1], i, PA[i]);
    }

    for ( i=1; i<=150; i++ ) p->PD[i][1] = PA[i];
    for ( i=1; i<=150; i++ ) p->PD[i][2] = PB[i];
    for ( i=1; i<=150; i++ ) p->PD[i][3] = PC[i];
    for ( i=1; i<=150; i++ ) p->PD[i][4] = PD[i];
    for ( i=1; i<=150; i++ ) p->PD[i][5] = PE[i];
    for ( i=1; i<=150; i++ ) p->PD[i][6] = PF[i];
    for ( i=1; i<=150; i++ ) p->PD[i][7] = PG[i];
    for ( i=1; i<=150; i++ ) p->PD[i][8] = PH[i];
    for ( i=1; i<=150; i++ ) p->PD[i][9] = PI[i];
    // Also set up a C-style row-major version of this array.
    //LGM_ARRAY_2D( p->PMA, 9+2, 150+2, double );
    LGM_ARRAY_2D( p->PD_rc, 9+2, 150+2, double );
    for ( i=1; i<=150; i++ ) {
        for ( j=1; j<=9; j++ ) {
            p->PD_rc[j][i] = p->PD[i][j];
        }
    }


    /*
     * Copy over PS array
     */
    LGM_ARRAY_1D( p->PS, 150+2, double );
    for ( i=1; i<=150; i++ ) p->PS[i] = PJ[i];



    /*
     * Set up the PDL array in correct order.
     */
    LGM_ARRAY_2D( p->PDL, 25+1, 2+1, double );
    for ( i=1; i<=25; i++ ) p->PDL[i][1] = PK1[i];
    for ( i=1; i<=25; i++ ) p->PDL[i][2] = PK1[i+25];


    /*
     * Set up the PD array in correct order.
     */
    LGM_ARRAY_2D( p->PTL, 150+1, 9+1, double );
    for ( i=1; i<=150; i++ ) p->PTL[i][1] = PL[i];
    for ( i=1; i<=150; i++ ) p->PTL[i][2] = PM[i];
    for ( i=1; i<=150; i++ ) p->PTL[i][3] = PN[i];
    for ( i=1; i<=150; i++ ) p->PTL[i][4] = PO[i];
    // Also set up a C-style row-major version of this array.
    LGM_ARRAY_2D( p->PTL_rc, 9+1, 150+1, double );
    for ( i=1; i<=150; i++ ) {
        for ( j=1; j<=9; j++ ) {
            p->PTL_rc[j][i] = p->PTL[i][j];
        }
    }


    /*
     * Set up the PD array in correct order.
     */
    LGM_ARRAY_2D( p->PMA, 100+1, 10+1, double );
    for ( i=1; i<=100; i++ ) p->PMA[i][1]  = PP[i];
    for ( i=1; i<=100; i++ ) p->PMA[i][2]  = PQ[i];
    for ( i=1; i<=100; i++ ) p->PMA[i][3]  = PR[i];
    for ( i=1; i<=100; i++ ) p->PMA[i][4]  = PS[i];
    for ( i=1; i<=100; i++ ) p->PMA[i][5]  = PU[i];
    for ( i=1; i<=100; i++ ) p->PMA[i][6]  = PV[i];
    for ( i=1; i<=100; i++ ) p->PMA[i][7]  = PW[i];
    for ( i=1; i<=100; i++ ) p->PMA[i][8]  = PX[i];
    for ( i=1; i<=100; i++ ) p->PMA[i][9]  = PY[i];
    for ( i=1; i<=100; i++ ) p->PMA[i][10] = PZ[i];
    // Also set up a C-style row-major version of this array.
    LGM_ARRAY_2D( p->PMA_rc, 10+1, 100+1, double );
    for ( i=1; i<=100; i++ ) {
        for ( j=1; j<=10; j++ ) {
            p->PMA_rc[j][i] = p->PMA[i][j];
        }
    }
    
    


    /*
     * Copy over SAM array
     */
    LGM_ARRAY_1D( p->SAM, 100+1, double );
    for ( i=1; i<=100; i++ ) p->SAM[i] = PAA[i];


    /*
     * Copy over PAVGM array
     */
    LGM_ARRAY_1D( p->PAVGM, 10+1, double );
    for ( i=1; i<=10; i++ ) p->PAVGM[i] = PAVGM[i];


    /*
     * Copy over PTM array
     */
    LGM_ARRAY_1D( p->PTM, 100+1, double );
    for ( i=1; i<=100; i++ ) p->PTM[i] = PTM[i];


    /*
     * Set up the PDM array in correct order.
     */
    LGM_ARRAY_2D( p->PDM, 10+1, 8+1, double );
    for ( i=1; i<=10; i++ ) p->PDM[i][1] = PDM[1][i];
    for ( i=1; i<=10; i++ ) p->PDM[i][2] = PDM[2][i];
    for ( i=1; i<=10; i++ ) p->PDM[i][3] = PDM[3][i];
    for ( i=1; i<=10; i++ ) p->PDM[i][4] = PDM[4][i];
    for ( i=1; i<=10; i++ ) p->PDM[i][5] = PDM[5][i];
    for ( i=1; i<=10; i++ ) p->PDM[i][6] = PDM[6][i];
    for ( i=1; i<=10; i++ ) p->PDM[i][7] = PDM[7][i];
    for ( i=1; i<=10; i++ ) p->PDM[i][8] = PDM[8][i];


    // SWITCH Arrays and vars
    p->ISW = 0;
    LGM_ARRAY_1D( p->SW,  25+1, double );
    LGM_ARRAY_1D( p->SWC, 25+1, double );
    for ( i=1; i<=25; i++ ) {
        p->SW[i]  = 1.0;
        p->SWC[i] = 1.0;
    }


    p->DAYL = -1.0;
    p->XL   = 1000.0;
    p->TLL  = 1000.0;
    for ( i=1; i<=25; i++ ) {
        p->SV[i] = 1.0;
    }


    // Init some vars in VTST7
    for (i=1; i<=2; i++ ) {
        p->IYDL[i]  = -999;
        p->SECL[i]  = -999.0;
        p->GLATL[i] = -999.0;
        p->GLL[i]   = -999.0;
        p->STLL[i]  = -999.0;
        p->FAL[i]   = -999.0;
        p->FL[i]    = -999.0;
        for (j=1; j<=7; j++ ) {
            p->APL[j][i]    = -999.0;
        }
        for (j=1; j<=25; j++ ) {
            p->SWL[j][i]   = -999.0;
            p->SWCL[j][i]  = -999.0;
        }
    }





    /*
     * Some variables that need to be saved/initialized.
     */
    p->GTD7_ALAST = 99999.0;
    p->GTD7_MSSL  = -999;

    p->GTS7_ALAST = -999.0;

    p->GLOBE7_P14 = -1000.0;
    p->GLOBE7_P18 = -1000.0;
    p->GLOBE7_P32 = -1000.0;
    p->GLOBE7_P39 = -1000.0;

    p->GLOBE7_CD14 = -9e99;
    p->GLOBE7_CD18 = -9e99;
    p->GLOBE7_CD32 = -9e99;
    p->GLOBE7_CD39 = -9e99;


    p->GLOB7S_P14 = -1000.0;
    p->GLOB7S_P18 = -1000.0;
    p->GLOB7S_P32 = -1000.0;
    p->GLOB7S_P39 = -1000.0;

    p->GLOB7S_CD14 = -9e99;
    p->GLOB7S_CD18 = -9e99;
    p->GLOB7S_CD32 = -9e99;
    p->GLOB7S_CD39 = -9e99;











    return( p );

}


void Lgm_FreeMsis00( Lgm_Msis00Info *p ) {

    LGM_ARRAY_1D_FREE( p->PT );
    LGM_ARRAY_2D_FREE( p->PD );
    LGM_ARRAY_2D_FREE( p->PD_rc );
    LGM_ARRAY_1D_FREE( p->PS );
    LGM_ARRAY_2D_FREE( p->PDL );
    LGM_ARRAY_2D_FREE( p->PTL );
    LGM_ARRAY_2D_FREE( p->PTL_rc );
    LGM_ARRAY_2D_FREE( p->PMA );
    LGM_ARRAY_2D_FREE( p->PMA_rc );
    LGM_ARRAY_1D_FREE( p->SAM );
    LGM_ARRAY_1D_FREE( p->PAVGM );
    LGM_ARRAY_1D_FREE( p->PTM );
    LGM_ARRAY_2D_FREE( p->PDM );
    LGM_ARRAY_1D_FREE( p->SW );
    LGM_ARRAY_1D_FREE( p->SWC );

    free( p );

    return;
}










