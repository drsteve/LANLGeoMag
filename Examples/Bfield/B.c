#include <Lgm_MagModelInfo.h>

/* BAL 04-Jan-2011 modified for more cases */

int main(){


    long int            Date;
    double              UTC;
    Lgm_Vector          u, B, ugsm;
    Lgm_MagModelInfo    *mInfo;
    int                 i, j;


    Date = 20050831;
    UTC  = 9.0;

    mInfo = Lgm_InitMagInfo( );
    Lgm_Set_Coord_Transforms( Date, UTC, mInfo->c );


    /*
     * Can set model manually ...
     */
    //mInfo->Bfield = Lgm_B_T89;
    //mInfo->Bfield = Lgm_B_TS04;
    //mInfo->Bfield = Lgm_B_TS04;
    //mInfo->Bfield = Lgm_B_T87;
    //mInfo->Bfield = Lgm_B_T96;
    mInfo->Bfield = Lgm_B_TS07;



    /*
     * Or there is a setter routine ...
     */
    //Lgm_MagModelInfo_Set_MagModel( LGM_EDIP, LGM_EXTMODEL_T89, mInfo );



    /*
     * For TS07, the coeffs need to be initialized for each new time...
     */
    Lgm_SetCoeffs_TS07( Date, UTC, &mInfo->TS07_Info );




    /*
     *  Can set/over-ride Qin-Denton parameters manually ....
     */
    mInfo->P      = 4.1011111111111118;
    mInfo->Dst    = 7.7777777777777777;
    mInfo->By     = 3.7244444444444444;
    mInfo->Bz     = -0.12666666666666665;
    mInfo->W[0]   = 0.12244444444444445;
    mInfo->W[1]   = 0.2514;
    mInfo->W[2]   = 0.089266666666666661;
    mInfo->W[3]   = 0.047866666666666668;
    mInfo->W[4]   = 0.22586666666666666;
    mInfo->W[5]   = 1.0461333333333334;


    /*
     *  Or Qin-Denton parameters can be obtained automatically by date/time ....
     */
//    JD = Lgm_Date_to_JD( Date, UTC, mInfo->c );     // Compute JD.
//    Lgm_get_QinDenton_at_JD( JD, &p, 1 );           // Get (interpolate) the QinDenton vals 
//                                                    // from the values in the file at the 
//                                                    // given Julian Date.
//    Lgm_set_QinDenton( &p, mInfo );                 // Set params in mInfo structure.


    printf("%13s", "Kp");
    printf("%13s", "Ugsmx (Re)");
    printf("%13s", "Ugsmy (Re)");
    printf("%13s", "Ugsmz (Re)");
    printf("%13s", "Bgsmx (nT)");
    printf("%13s", "Bgsmy (nT)");
    printf("%13s", "Bgsmz (nT)");
    printf("%13s", "Bmag (nT)\n");


    for (i=0; i<=5; i++) {
        mInfo->Kp = i;
        u.x = -6.6; u.y =  0.0;  u.z =  0.0;
        //Lgm_Convert_Coords( &u, &ugsm, GEO_TO_GSM, mInfo->c );
        Lgm_Convert_Coords( &u, &ugsm, SM_TO_GSM, mInfo->c );
        mInfo->Bfield( &ugsm, &B, mInfo );
        printf( "%13i", mInfo->Kp);
        printf( "%13g%13g%13g", ugsm.x, ugsm.y, ugsm.z );
        printf( "%13g%13g%13g", B.x, B.y, B.z );
        printf( "%13g\n", Lgm_Magnitude( &B ) );
    }





    for (j=0; j<100; j++){
        mInfo->Kp = 3;
        for (i=0; i<13; i++) {
            u.x = -1.0 - (double)i * 0.5;
            u.y =  0.0;  u.z =  0.0;
            mInfo->Bfield( &u, &B, mInfo );
            printf( "%13i", mInfo->Kp);
            printf( "%13g%13g%13g", u.x, u.y, u.z );
            printf( "%13g%13g%13g", B.x, B.y, B.z );
            printf( "%13g\n", Lgm_Magnitude( &B ) );
        }
    }






    Lgm_FreeMagInfo( mInfo );


    exit(0);
}
