#include <Lgm_MagModelInfo.h>
#include <Lgm_QinDenton.h>

int main(){


    long int            Date;
    double              Time;
    Lgm_Vector          u, v1, v2, v3;
    Lgm_MagModelInfo    *MagInfo;
    Lgm_DateTime        UTC;
    Lgm_QinDentonOne    p;
    
    Date = 20100203;
    Time = 12.34567;
    
    Date = 20120903;
    Time = 14.0 + 57.0/60.0 + 30.0/3600.0;

    MagInfo = Lgm_InitMagInfo( );
    Lgm_Set_Coord_Transforms( Date, Time, MagInfo->c );

    //Lgm_MagModelInfo_Set_MagModel( LGM_IGRF, LGM_EXTMODEL_TS07, MagInfo );
    Lgm_MagModelInfo_Set_MagModel( LGM_IGRF, LGM_EXTMODEL_TS04, MagInfo );

    /*
     * Fill UTC DateTime structure
     */
    Lgm_Make_UTC( Date, Time, &UTC, MagInfo->c );


    /*
     * Get QinDenton Parameters for the given Date and Time
     */
    Lgm_get_QinDenton_at_JD( UTC.JD, &p, 1, 0 );
    Lgm_set_QinDenton( &p, MagInfo );

//    Lgm_SetCoeffs_TS07( Date, Time, &(MagInfo->TS07_Info) );




    u.x = 2.74384; u.y =  5.96585;  u.z =  -0.815308; // Re

    Lgm_Trace( &u, &v1, &v2, &v3, 120.0, 1e-7, 1e-7, MagInfo );
    printf( "u = %g %g %g Re\n", u.x, u.y, u.z );
    printf( "v1 = %g %g %g Re\n", v1.x, v1.y, v1.z );
    printf( "v2 = %g %g %g Re\n", v2.x, v2.y, v2.z );
    printf( "v3 = %g %g %g Re\n", v3.x, v3.y, v3.z );



    Lgm_FreeMagInfo( MagInfo );


    exit(0);
}
