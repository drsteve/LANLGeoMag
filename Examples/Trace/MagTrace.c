#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <Lgm_CTrans.h>
#include <Lgm_QinDenton.h>
#include <Lgm_ElapsedTime.h>
#include <Lgm_FluxToPsd.h>
#include <omp.h>

int ClassifyFL_Enhanced( int Flag, Lgm_Vector *w ) {

    int NewFlag = 0;

    if ( Flag < 0 ) return( Flag );

    switch (Flag) {

            case LGM_OPEN_IMF:
                                break;
            case LGM_CLOSED:
                                if ( w->x > 0.0 ) {
                                    NewFlag = 1;
                                } else {
                                    NewFlag = ( fabs(w->y) > 10.0 ) ? 2 : 3;
                                    //NewFlag = 2+fabs( w->x )*10.0;
                                }
                                break;
            case LGM_OPEN_N_LOBE:
                                NewFlag = ( (fabs(w->z) < 20.0) && (fabs(w->y) > 10.0 ) ) ? 5 : 4;
                                NewFlag = ( fabs(w->z) < 10.0 ) ? 5 : 4;
//printf("w = %g %g %g NewFlag = %d\n", w->x, w->y, w->z, NewFlag);
                                break;
            case LGM_OPEN_S_LOBE:
                                NewFlag = ( fabs(w->z) < 10.0 ) ? 6 : 4;
                                break;


    }

    return( NewFlag );

}

int ClassifyFL_Enhanced2( int Flag, Lgm_Vector *v1, Lgm_Vector *v2, Lgm_Vector *v3, Lgm_MagModelInfo *mInfo ) {

    Lgm_Vector w;
    int NewFlag = 0;

    if ( Flag < 0 ) return( Flag );

    Lgm_Convert_Coords( v3, &w, GSM_TO_SM, mInfo->c );

    switch (Flag) {

            case LGM_OPEN_IMF:
                                break;
            case LGM_CLOSED:
                                if ( w.x > 0.0 ) {
                                    NewFlag = 1;
                                } else {
                                    NewFlag = ( fabs(w.y) > 10.0 ) ? 2 : 3;
                                    //NewFlag = 2+fabs( w.x )*10.0;
                                }
                                break;
            case LGM_OPEN_N_LOBE:
                                //NewFlag = ( (fabs(w.z) < 20.0) && (fabs(w.y) > 10.0 ) ) ? 5 : 4;
                                //NewFlag = ( fabs(w.z) < 10.0 ) ? 5 : 4;
                                NewFlag = (( fabs(v3->z) < 10.0 ) && ( mInfo->v1_final.x > -50.0 )) ? 5 : 4;
if (NewFlag == 5){
    printf("mInfo->v1_final = %g %g %g\n", mInfo->v1_final.x, mInfo->v1_final.y, mInfo->v1_final.z );
}
                                break;
            case LGM_OPEN_S_LOBE:
                                NewFlag = ( fabs(w.z) < 10.0 ) ? 6 : 4;
                                break;


    }

    return( NewFlag );

}

int main(){

    int                 NX, NY, i, j;
    double              LX_MIN, LX_MAX;
    double              LY_MIN, LY_MAX;
    double              x, y, GeodLat, GeodLong, GeodHeight ;
    double              UTC, JD;
    long int            Date;
    int                 Flag;
    double              **Image;
    Lgm_QinDentonOne    p;
    Lgm_Vector          u, v, v1, v2, v3, w;
    Lgm_MagModelInfo    *mInfo2;
    Lgm_MagModelInfo    *mInfo = Lgm_InitMagInfo();
    

    NX     = 1000;
    LX_MIN = -30.0;
    LX_MAX =  30.0;

    NY     = 1000;
    LY_MIN = -30.0;
    LY_MAX =  30.0;

    GeodHeight = 100.0;  //km

    Date = 20121210;
    UTC  = 16.0 + 27.0/60.0 + 29.99/3600.0;
    JD = Lgm_Date_to_JD( Date, UTC, mInfo->c );
    Lgm_Set_Coord_Transforms( Date, UTC, mInfo->c );
    printf("Tilt = %g\n", mInfo->c->psi*DegPerRad );
    Lgm_get_QinDenton_at_JD( JD, &p, 1 );
    Lgm_set_QinDenton( &p, mInfo );

    Lgm_MagModelInfo_Set_MagModel( LGM_IGRF, LGM_EXTMODEL_TS04, mInfo );
    Lgm_MagModelInfo_Set_MagModel( LGM_IGRF, LGM_EXTMODEL_T89, mInfo );
    Lgm_MagModelInfo_Set_MagModel( LGM_IGRF, LGM_EXTMODEL_T96, mInfo );
    
    Lgm_Set_Open_Limits( mInfo, -100.0, 30.0, -40.0, 40.0, -40.0, 40.0 );


    LGM_ARRAY_2D( Image, NX, NY, double );

    { // BEGIN PARALLEL EXECUTION
        #pragma omp parallel private(x,y,j,GeodLat,GeodLong,u,v,v1,v2,v3,Flag,mInfo2)
        #pragma omp for schedule(dynamic, 1)
        for ( i=0; i<NX; i++ ) {
        //for ( i=25; i<=25; i++ ) {
            x = (LX_MAX-LX_MIN) * i / ((double)(NX-1)) + LX_MIN;

            mInfo2 = Lgm_CopyMagInfo( mInfo );

            printf("i=%d\n", i);
            
            for ( j=0; j<NY; j++ ) {
            //for ( j=NY/2; j<=NY/2; j++ ) {
                y = (LY_MAX-LY_MIN) * j / ((double)(NY-1)) + LY_MIN;

                GeodLat  = 90.0 - sqrt( x*x + y*y );
                GeodLong = atan2( y, x )*DegPerRad;
                Lgm_GEOD_to_WGS84( GeodLat, GeodLong, GeodHeight, &v );
                Lgm_Convert_Coords( &v, &u, WGS84_TO_GSM, mInfo->c );



//mInfo2->VerbosityLevel = 5;
                Flag = Lgm_Trace( &u, &v1, &v2, &v3, GeodHeight, 1e-7, 1e-7, mInfo2 );
                Lgm_Convert_Coords( &v3, &w, GSM_TO_SM, mInfo->c );
                //Image[i][j] = (double)ClassifyFL_Enhanced(Flag, &w)+100;
                Image[i][j] = (double)ClassifyFL_Enhanced2( Flag, &v1, &v2, &v3, mInfo2 )+100;

                //printf("x, y, = %g %g GeodLat, GeodLong, GeodHeight = %g %g %g     u = %g %g %g   Flag = %d\n", x, y, GeodLat, GeodLong, GeodHeight, u.x, u.y, u.z, Flag);

            }

            Lgm_FreeMagInfo( mInfo2 );

        }
    } // END PARALLEL EXECUTION



    DumpGif( "Image_T96_3", NX, NY, Image );
    LGM_ARRAY_2D_FREE( Image );

    Lgm_FreeMagInfo( mInfo );
    return(0);


}
