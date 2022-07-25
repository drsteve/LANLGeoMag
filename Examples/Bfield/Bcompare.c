#include <stdio.h>
#include <Lgm_MagModelInfo.h>


int main(){

    long int            Date;
    double              UTC;
    Lgm_Vector          B, ugsm;
    Lgm_MagModelInfo    *mInfo;

    Date = 20080102;
    UTC  = 0.04; //22.0;

    mInfo = Lgm_InitMagInfo( );

    /*
     * Set TA2016 using setter routine ...
     */
    Lgm_MagModelInfo_Set_MagModel( LGM_IGRF, LGM_EXTMODEL_TA16, mInfo );

    Lgm_Init_TA16( &mInfo->TA16_Info, 1 );  // use verbose output

    for (int offset=0; offset<6; offset++) {
        Date += offset;
        UTC += 0.01*(double)offset;
        Lgm_Set_Coord_Transforms( Date, UTC, mInfo->c );

        /*
         * For TA16, the coeffs need to be initialized for each new time...
         */
        Lgm_SetCoeffs_TA16( Date, UTC, &mInfo->TA16_Info );

        // Like other Tsyganenko models, input is in GSM
        ugsm.x = -6.5;
        ugsm.y = 1.5;
        ugsm.z = 0.5;

        mInfo->Bfield( &ugsm, &B, mInfo );

        // Set header line for console output
        printf("%13s", "SymhHc (nT)");
        printf("%13s", "Ugsmx (Re)");
        printf("%13s", "Ugsmy (Re)");
        printf("%13s", "Ugsmz (Re)");
        printf("%13s", "Bgsmx (nT)");
        printf("%13s", "Bgsmy (nT)");
        printf("%13s", "Bgsmz (nT)");
        printf("%13s", "Bmag (nT)\n");
        printf( "%13g", mInfo->TA16_Info.SymHc_avg);
        printf( "%13g%13g%13g", ugsm.x, ugsm.y, ugsm.z );
        printf( "%13g%13g%13g", B.x, B.y, B.z );
        printf( " %13.2lf\n", Lgm_Magnitude( &B ) );
    }
Lgm_FreeMagInfo( mInfo );
exit(0);
}
