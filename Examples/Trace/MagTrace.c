#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <Lgm_CTrans.h>
#include <Lgm_QinDenton.h>
#include <Lgm_ElapsedTime.h>
#include <Lgm_FluxToPsd.h>
#include <omp.h>


void DumpImage( char *FilenameBase, int W, int H, double **Image ){

    double           Val, Val2, Min, Max, dVal;
    int              w, h;
    unsigned char   *uImage, *Red, *Grn, *Blu, uVal;
    char            Filename[1024];
    FILE            *fp_gif, *fp_info;

    int             LogScale;

    uImage = (unsigned char *)calloc( W*H, sizeof(unsigned char) );
    Red = (unsigned char *)calloc( 256, sizeof(unsigned char) );
    Grn = (unsigned char *)calloc( 256, sizeof(unsigned char) );
    Blu = (unsigned char *)calloc( 256, sizeof(unsigned char) );
    Red[1] =   0; Grn[1] =  105;  Blu[1] =   27;
    Red[2] = 152/2; Grn[2] =  152/2;  Blu[2] =  122/2;
    Red[3] =  98/2; Grn[3] =   98/2;  Blu[3] =   77/2;
    Red[4] =  19; Grn[4] =   83;  Blu[4] =  194;
    Red[5] = 150; Grn[5] =   20;  Blu[5] =   20;
    Red[6] = 230; Grn[6] =  230;  Blu[6] =   19;
    Red[7] =  88; Grn[7] =    0;  Blu[7] =   88;
    Red[8] =   0; Grn[8] =    0;  Blu[8] =    0;
    Red[80] =   90; Grn[80] =    255;  Blu[80] =    255;

    for ( w=0; w<W; w++ ){
        for ( h=0; h<H; h++ ) {

            *(uImage + W*(H-1-h) + w) = (unsigned char)( (int)(Image[h][w]+.1));

        }
    }

    sprintf( Filename, "%s.gif", FilenameBase);
    fp_gif = fopen(Filename, "w");
    WriteGIF(fp_gif, (byte *)uImage, 0, W, H, Red, Grn, Blu, 256, 0, "");
    fclose(fp_gif);

    free( uImage );
    free( Red );
    free( Grn );
    free( Blu );



}
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
            default:
                                break;


    }

    return( NewFlag );

}

int ClassifyFL_Enhanced2( int Flag, Lgm_Vector *v1, Lgm_Vector *v2, Lgm_Vector *v3, Lgm_MagModelInfo *mInfo ) {

    Lgm_Vector w;
    int NewFlag = 0;


    if ( Flag < 0 ) return( 8 );

    Lgm_Convert_Coords( v3, &w, GSM_TO_SM, mInfo->c );

    switch (Flag) {

            case LGM_OPEN_IMF:
                                NewFlag = 0;
                                break;
            case LGM_CLOSED:
                                if ( w.x > 0.0 ) {
                                //if ( v3->x > 0.0 ) {
                                    //printf("v3->x = %g\n", v3->x );
                                    NewFlag = 1;
                                } else {
                                    NewFlag = ( fabs(w.y) > 20.0 ) ? 2 : 3;
                                    //NewFlag = 2+fabs( w.x )*10.0;
                                }
                                break;
            case LGM_OPEN_N_LOBE:
                                //NewFlag = ( fabs(w.z) < 10.0 ) ? 5 : 4;
                                //NewFlag = (( fabs(v3->z) < 10.0 ) && ( mInfo->v1_final.x > 0.0 )) ? 5 : 4;
                                NewFlag = ( mInfo->v1_final.x > 0.0 ) ? 5 : 4;
                                //NewFlag = (( mInfo->v1_final.x > -50.0 )) ? 5 : 4;
                                break;
            case LGM_OPEN_S_LOBE:
                                //NewFlag = ( (fabs(w.z) < 20.0) && (fabs(w.y) > 10.0 ) ) ? 5 : 4;
                                //NewFlag = ( fabs(w.z) < 10.0 ) ? 6 : 4;
                                //NewFlag = (( fabs(v3->z) < 10.0 ) && ( mInfo->v1_final.x > 0.0 )) ? 6 : 4;
                                NewFlag = ( mInfo->v2_final.x > 0.0 ) ? 6 : 4;
                                //NewFlag = (( mInfo->v1_final.x > -50.0 )) ? 6 : 4;
                                break;
            default:
                                NewFlag = 8;
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
    double              **Image, **ImageSouth, **ImageEq, **ImageYZ15, **ImageYZ30, **ImageYZ45, R, MLAT, MLT;

    int                 EQ_NX, EQ_NY, ii, jj;
    double              EQ_XMIN, EQ_XMAX, EQ_YMIN, EQ_YMAX;

    int                 YZ_NY, YZ_NZ;
    double              YZ_YMIN, YZ_YMAX, YZ_ZMIN, YZ_ZMAX;

    double              Hdid, Hnext, s;
    int                 reset=0;
    Lgm_Vector          u_scale;


    Lgm_QinDentonOne    p;
    Lgm_Vector          u, v, v1, v2, v3, w, ww, ww2;
    Lgm_MagModelInfo    *mInfo2;
    Lgm_MagModelInfo    *mInfo = Lgm_InitMagInfo();
    

    NX     = 50; LX_MIN = -30.0; LX_MAX =  30.0;
    NY     = 50; LY_MIN = -30.0; LY_MAX =  30.0;



    EQ_NX   = 200; EQ_XMIN = -60.0; EQ_XMAX = 20.0;
    EQ_NY   = 200; EQ_YMIN = -40.0; EQ_YMAX = 40.0;


    YZ_NY   = 200; YZ_YMIN = -40.0; YZ_YMAX = 40.0;
    YZ_NZ   = 200; YZ_ZMIN = -40.0; YZ_ZMAX = 40.0;



    u_scale.x =  100.0;  u_scale.y = 100.0; u_scale.z = 100.0;

    GeodHeight = 100.0;  //km

    Date = 20020417;
    UTC  = 10.0 + 50.0/60.0 + 0.0/3600.0;
Date = 20080310;
UTC = 0.0;
    JD = Lgm_Date_to_JD( Date, UTC, mInfo->c );
    Lgm_Set_Coord_Transforms( Date, UTC, mInfo->c );
    printf("Tilt = %g\n", mInfo->c->psi*DegPerRad );
//exit(0);
    Lgm_get_QinDenton_at_JD( JD, &p, 1 );
    Lgm_set_QinDenton( &p, mInfo );

    Lgm_MagModelInfo_Set_MagModel( LGM_IGRF, LGM_EXTMODEL_T01S, mInfo );
    Lgm_MagModelInfo_Set_MagModel( LGM_IGRF, LGM_EXTMODEL_T89, mInfo );
    Lgm_MagModelInfo_Set_MagModel( LGM_IGRF, LGM_EXTMODEL_T87, mInfo );
    Lgm_MagModelInfo_Set_MagModel( LGM_IGRF, LGM_EXTMODEL_OP77, mInfo );
    Lgm_MagModelInfo_Set_MagModel( LGM_IGRF, LGM_EXTMODEL_T96, mInfo );
    Lgm_MagModelInfo_Set_MagModel( LGM_IGRF, LGM_EXTMODEL_TS04, mInfo );
    Lgm_MagModelInfo_Set_MagModel( LGM_IGRF, LGM_EXTMODEL_TU82, mInfo );
    Lgm_MagModelInfo_Set_MagModel( LGM_IGRF, LGM_EXTMODEL_OP88, mInfo );
Lgm_MagModelInfo_Set_MagModel( LGM_IGRF, LGM_EXTMODEL_TS07, mInfo );
Lgm_SetCoeffs_TS07( Date, UTC, &mInfo->TS07_Info );
    Lgm_MagModelInfo_Set_MagModel( LGM_IGRF, LGM_EXTMODEL_T02, mInfo );
    
    Lgm_Set_Open_Limits( mInfo, -60.0, 30.0, -40.0, 40.0, -40.0, 40.0 );


    LGM_ARRAY_2D( Image, NX, NY, double );
    LGM_ARRAY_2D( ImageSouth, NX, NY, double );
    LGM_ARRAY_2D( ImageEq, EQ_NX, EQ_NY, double );
    LGM_ARRAY_2D( ImageYZ15, YZ_NY, YZ_NZ, double );
    LGM_ARRAY_2D( ImageYZ30, YZ_NY, YZ_NZ, double );
    LGM_ARRAY_2D( ImageYZ45, YZ_NY, YZ_NZ, double );

    { // BEGIN PARALLEL EXECUTION
        #pragma omp parallel private(ii,jj,x,y,j,GeodLat,GeodLong,u,v,v1,v2,v3,ww,ww2,Flag,mInfo2)
        #pragma omp for schedule(dynamic, 1)
        for ( i=0; i<NX; i++ ) {
            x = (LX_MAX-LX_MIN) * i / ((double)(NX-1)) + LX_MIN;

            mInfo2 = Lgm_CopyMagInfo( mInfo );

            printf("i=%d\n", i);
            
            for ( j=0; j<NY; j++ ) {
                y = (LY_MAX-LY_MIN) * j / ((double)(NY-1)) + LY_MIN;




                /*
                 * Trace From the North
                 */
                R    = 120.0/Re + 1.0;
                MLAT = 90.0 - sqrt( x*x + y*y );
                MLT  = atan2( y, x )*DegPerRad/15.0+12.0;
                Lgm_R_MLAT_MLT_to_CDMAG( R, MLAT, MLT, &v, mInfo->c );
                Lgm_Convert_Coords( &v, &u, CDMAG_TO_GSM, mInfo->c );

                v3.x = v3.y = v3.z = -1e31;
                Flag = Lgm_Trace( &u, &v1, &v2, &v3, GeodHeight, 1e-7, 1e-7, mInfo2 );

                Lgm_Convert_Coords( &v2, &w, GSM_TO_SM, mInfo->c );
                if ( w.x > 0.0 ) Lgm_TraceToSMEquat( &v2, &v3, 1e-6, mInfo2 );

                Lgm_Convert_Coords( &v3, &w, GSM_TO_SM, mInfo->c );
                Image[i][j] = (double)ClassifyFL_Enhanced2( Flag, &v1, &v2, &v3, mInfo2 );
                if (Image[i][j] < 0.0) Image[i][j] = 8.0;

                /*
                 *  Add point to the Equatorial plane image.
                 */
                if ( ( Flag == LGM_CLOSED ) || ( Image[i][j] >= 5 ) ) {
                    ii = (v3.y - EQ_YMIN)/(EQ_YMAX-EQ_YMIN) * (EQ_NY-1);
                    jj = (v3.x - EQ_XMIN)/(EQ_XMAX-EQ_XMIN) * (EQ_NX-1);
                    //printf("ii, jj = %d %d\n", ii, jj );
                    if ( (ii>=0)&&(ii<EQ_NY)&&(jj>=0)&&(jj<EQ_NX) ) {
                        ImageEq[ii][jj] = Image[i][j];
                    }
                }

                /*
                 *  Trace to X = -15Re
                 */
                if ( Lgm_TraceToYZPlane( &u, &ww, -15.0, -1.0, 1e-7, mInfo2 ) > 0 ){
                    //printf( "ww = %g %g %g\n", ww.x, ww.y, ww.z);
                    ii = (ww.z - YZ_ZMIN)/(YZ_ZMAX-YZ_ZMIN) * (YZ_NZ-1);
                    jj = (ww.y - YZ_YMIN)/(YZ_YMAX-YZ_YMIN) * (YZ_NY-1);
                    if ( (ii>=0)&&(ii<YZ_NZ)&&(jj>=0)&&(jj<YZ_NY) ) {
                        ImageYZ15[ii][jj] = Image[i][j];
                    }

                    /*
                     *  Trace to X = -15Re for second hit
                     */
                    Lgm_MagStep( &ww, &u_scale, 1e-2, &Hdid, &Hnext, -1.0, &s, &reset, mInfo2->Bfield, mInfo2 );
                    if ( Lgm_TraceToYZPlane( &ww, &ww2, -15.0, -1.0, 1e-4, mInfo2 ) > 0 ){
                        //printf( "ww = %g %g %g\n", ww.x, ww.y, ww.z);
                        ii = (ww2.z - YZ_ZMIN)/(YZ_ZMAX-YZ_ZMIN) * (YZ_NZ-1);
                        jj = (ww2.y - YZ_YMIN)/(YZ_YMAX-YZ_YMIN) * (YZ_NY-1);
                        if ( (ii>=0)&&(ii<YZ_NZ)&&(jj>=0)&&(jj<YZ_NY) ) {
                            ImageYZ15[ii][jj] = Image[i][j];
                        }
                    } 

                } 

                /*
                 *  Trace to X = -30Re
                 */
                if ( Lgm_TraceToYZPlane( &u, &ww, -30.0, -1.0, 1e-7, mInfo2 ) > 0 ){
                    //printf( "ww = %g %g %g\n", ww.x, ww.y, ww.z);
                    ii = (ww.z - YZ_ZMIN)/(YZ_ZMAX-YZ_ZMIN) * (YZ_NZ-1);
                    jj = (ww.y - YZ_YMIN)/(YZ_YMAX-YZ_YMIN) * (YZ_NY-1);
                    if ( (ii>=0)&&(ii<YZ_NZ)&&(jj>=0)&&(jj<YZ_NY) ) {
                        ImageYZ30[ii][jj] = Image[i][j];
                    }

                    /*
                     *  Trace to X = -30Re for second hit
                     */
                    Lgm_MagStep( &ww, &u_scale, 1e-2, &Hdid, &Hnext, -1.0, &s, &reset, mInfo2->Bfield, mInfo2 );
                    if ( Lgm_TraceToYZPlane( &ww, &ww2, -30.0, -1.0, 1e-4, mInfo2 ) > 0 ){
                        //printf( "ww = %g %g %g\n", ww.x, ww.y, ww.z);
                        ii = (ww2.z - YZ_ZMIN)/(YZ_ZMAX-YZ_ZMIN) * (YZ_NZ-1);
                        jj = (ww2.y - YZ_YMIN)/(YZ_YMAX-YZ_YMIN) * (YZ_NY-1);
                        if ( (ii>=0)&&(ii<YZ_NZ)&&(jj>=0)&&(jj<YZ_NY) ) {
                            ImageYZ30[ii][jj] = Image[i][j];
                        }
                    } 

                } 

                /*
                 *  Trace to X = -45Re
                 */
                if ( Lgm_TraceToYZPlane( &u, &ww, -45.0, -1.0, 1e-7, mInfo2 ) > 0 ){
                    //printf( "ww = %g %g %g\n", ww.x, ww.y, ww.z);
                    ii = (ww.z - YZ_ZMIN)/(YZ_ZMAX-YZ_ZMIN) * (YZ_NZ-1);
                    jj = (ww.y - YZ_YMIN)/(YZ_YMAX-YZ_YMIN) * (YZ_NY-1);
                    if ( (ii>=0)&&(ii<YZ_NZ)&&(jj>=0)&&(jj<YZ_NY) ) {
                        ImageYZ45[ii][jj] = Image[i][j];
                    }

                    /*
                     *  Trace to X = -45Re for second hit
                     */
                    Lgm_MagStep( &ww, &u_scale, 1e-2, &Hdid, &Hnext, -1.0, &s, &reset, mInfo2->Bfield, mInfo2 );
                    if ( Lgm_TraceToYZPlane( &ww, &ww2, -45.0, -1.0, 1e-4, mInfo2 ) > 0 ){
                        //printf( "ww = %g %g %g\n", ww.x, ww.y, ww.z);
                        ii = (ww2.z - YZ_ZMIN)/(YZ_ZMAX-YZ_ZMIN) * (YZ_NZ-1);
                        jj = (ww2.y - YZ_YMIN)/(YZ_YMAX-YZ_YMIN) * (YZ_NY-1);
                        if ( (ii>=0)&&(ii<YZ_NZ)&&(jj>=0)&&(jj<YZ_NY) ) {
                            ImageYZ45[ii][jj] = Image[i][j];
                        }
                    } 

                } 



                /*
                 * Trace From the South
                 */
                R    = 120.0/Re + 1.0;
                MLAT = 90.0 - sqrt( x*x + y*y );
                MLAT *= -1.0;
                MLT  = atan2( y, x )*DegPerRad/15.0+12.0;
                Lgm_R_MLAT_MLT_to_CDMAG( R, MLAT, MLT, &v, mInfo->c );
                Lgm_Convert_Coords( &v, &u, CDMAG_TO_GSM, mInfo->c );

                v3.x = v3.y = v3.z = -1e31;
                Flag = Lgm_Trace( &u, &v1, &v2, &v3, GeodHeight, 1e-7, 1e-7, mInfo2 );

                Lgm_Convert_Coords( &v2, &w, GSM_TO_SM, mInfo->c );
                if ( w.x > 0.0 ) Lgm_TraceToSMEquat( &v2, &v3, 1e-6, mInfo2 );

                Lgm_Convert_Coords( &v3, &w, GSM_TO_SM, mInfo->c );
                ImageSouth[i][j] = (double)ClassifyFL_Enhanced2( Flag, &v1, &v2, &v3, mInfo2 );
                if (ImageSouth[i][j] < 0.0) ImageSouth[i][j] = 8.0;

                /*
                 *  Add point to the Equatorial plane image.
                 */
                if ( ( Flag == LGM_CLOSED ) || ( ImageSouth[i][j] >= 5 ) ) {
                    ii = (v3.y - EQ_YMIN)/(EQ_YMAX-EQ_YMIN) * (EQ_NY-1);
                    jj = (v3.x - EQ_XMIN)/(EQ_XMAX-EQ_XMIN) * (EQ_NX-1);
                    //printf("ii, jj = %d %d\n", ii, jj );
                    if ( (ii>=0)&&(ii<EQ_NY)&&(jj>=0)&&(jj<EQ_NX) ) {
                        ImageEq[ii][jj] = ImageSouth[i][j];
                    }
                }

                /*
                 *  Trace to X = -15Re
                 */
                if ( Lgm_TraceToYZPlane( &u, &ww, -15.0, 1.0, 1e-7, mInfo2 ) > 0 ){
                    //printf( "ww = %g %g %g\n", ww.x, ww.y, ww.z);
                    ii = (ww.z - YZ_ZMIN)/(YZ_ZMAX-YZ_ZMIN) * (YZ_NZ-1);
                    jj = (ww.y - YZ_YMIN)/(YZ_YMAX-YZ_YMIN) * (YZ_NY-1);
                    if ( (ii>=0)&&(ii<YZ_NZ)&&(jj>=0)&&(jj<YZ_NY) ) {
                        ImageYZ15[ii][jj] = ImageSouth[i][j];
                    }

                    /*
                     *  Trace to X = -15Re for second hit
                     */
                    Lgm_MagStep( &ww, &u_scale, 1e-2, &Hdid, &Hnext, 1.0, &s, &reset, mInfo2->Bfield, mInfo2 );
                    if ( Lgm_TraceToYZPlane( &ww, &ww2, -15.0, 1.0, 1e-4, mInfo2 ) > 0 ){
                        //printf( "ww = %g %g %g\n", ww.x, ww.y, ww.z);
                        ii = (ww2.z - YZ_ZMIN)/(YZ_ZMAX-YZ_ZMIN) * (YZ_NZ-1);
                        jj = (ww2.y - YZ_YMIN)/(YZ_YMAX-YZ_YMIN) * (YZ_NY-1);
                        if ( (ii>=0)&&(ii<YZ_NZ)&&(jj>=0)&&(jj<YZ_NY) ) {
                            ImageYZ15[ii][jj] = ImageSouth[i][j];
                        }
                    } 

                } 

                /*
                 *  Trace to X = -30Re
                 */
                if ( Lgm_TraceToYZPlane( &u, &ww, -30.0, 1.0, 1e-7, mInfo2 ) > 0 ){
                    //printf( "ww = %g %g %g\n", ww.x, ww.y, ww.z);
                    ii = (ww.z - YZ_ZMIN)/(YZ_ZMAX-YZ_ZMIN) * (YZ_NZ-1);
                    jj = (ww.y - YZ_YMIN)/(YZ_YMAX-YZ_YMIN) * (YZ_NY-1);
                    if ( (ii>=0)&&(ii<YZ_NZ)&&(jj>=0)&&(jj<YZ_NY) ) {
                        ImageYZ30[ii][jj] = ImageSouth[i][j];
                    }

                    /*
                     *  Trace to X = -30Re for second hit
                     */
                    Lgm_MagStep( &ww, &u_scale, 1e-2, &Hdid, &Hnext, 1.0, &s, &reset, mInfo2->Bfield, mInfo2 );
                    if ( Lgm_TraceToYZPlane( &ww, &ww2, -30.0, 1.0, 1e-4, mInfo2 ) > 0 ){
                        //printf( "ww = %g %g %g\n", ww.x, ww.y, ww.z);
                        ii = (ww2.z - YZ_ZMIN)/(YZ_ZMAX-YZ_ZMIN) * (YZ_NZ-1);
                        jj = (ww2.y - YZ_YMIN)/(YZ_YMAX-YZ_YMIN) * (YZ_NY-1);
                        if ( (ii>=0)&&(ii<YZ_NZ)&&(jj>=0)&&(jj<YZ_NY) ) {
                            ImageYZ30[ii][jj] = ImageSouth[i][j];
                        }
                    } 

                } 


                /*
                 *  Trace to X = -45Re
                 */
                if ( Lgm_TraceToYZPlane( &u, &ww, -45.0, 1.0, 1e-7, mInfo2 ) > 0 ){
                    //printf( "ww = %g %g %g\n", ww.x, ww.y, ww.z);
                    ii = (ww.z - YZ_ZMIN)/(YZ_ZMAX-YZ_ZMIN) * (YZ_NZ-1);
                    jj = (ww.y - YZ_YMIN)/(YZ_YMAX-YZ_YMIN) * (YZ_NY-1);
                    if ( (ii>=0)&&(ii<YZ_NZ)&&(jj>=0)&&(jj<YZ_NY) ) {
                        ImageYZ45[ii][jj] = ImageSouth[i][j];
                    }

                    /*
                     *  Trace to X = -45Re for second hit
                     */
                    Lgm_MagStep( &ww, &u_scale, 1e-2, &Hdid, &Hnext, 1.0, &s, &reset, mInfo2->Bfield, mInfo2 );
                    if ( Lgm_TraceToYZPlane( &ww, &ww2, -45.0, 1.0, 1e-4, mInfo2 ) > 0 ){
                        //printf( "ww = %g %g %g\n", ww.x, ww.y, ww.z);
                        ii = (ww2.z - YZ_ZMIN)/(YZ_ZMAX-YZ_ZMIN) * (YZ_NZ-1);
                        jj = (ww2.y - YZ_YMIN)/(YZ_YMAX-YZ_YMIN) * (YZ_NY-1);
                        if ( (ii>=0)&&(ii<YZ_NZ)&&(jj>=0)&&(jj<YZ_NY) ) {
                            ImageYZ45[ii][jj] = ImageSouth[i][j];
                        }
                    } 

                } 















            }

            Lgm_FreeMagInfo( mInfo2 );

        }
    } // END PARALLEL EXECUTION



    DumpImage( "ImageNorth_TS07", NX, NY, Image );
    DumpImage( "ImageSouth_TS07", NX, NY, ImageSouth );
    DumpImage( "ImageEq_TS07", EQ_NX, EQ_NY, ImageEq );
    DumpImage( "ImageYZ15_TS07", YZ_NY, YZ_NZ, ImageYZ15 );
    DumpImage( "ImageYZ30_TS07", YZ_NY, YZ_NZ, ImageYZ30 );
    DumpImage( "ImageYZ45_TS07", YZ_NY, YZ_NZ, ImageYZ45 );
    LGM_ARRAY_2D_FREE( Image );
    LGM_ARRAY_2D_FREE( ImageSouth );
    LGM_ARRAY_2D_FREE( ImageEq );
    LGM_ARRAY_2D_FREE( ImageYZ15 );
    LGM_ARRAY_2D_FREE( ImageYZ30 );
    LGM_ARRAY_2D_FREE( ImageYZ45 );

    Lgm_FreeMagInfo( mInfo );
    return(0);


}
