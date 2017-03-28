#include <math.h>
#include <Lgm_FluxToPsd.h>
#include <Lgm_QuadPack.h>

typedef struct Fs_Info {
    double  n;
    double  Beta;
} Fs_Info;
double  Fs_Integrand( double Theta, _qpInfo *qpInfo ) {
    double val, CosTmB;
    Fs_Info *fi;
    fi = (Fs_Info *)qpInfo;
    CosTmB = cos( Theta - fi->Beta );
    val = 0.5*pow( CosTmB, fi->n )*cos( Theta );
    return( val );
}
unsigned char Ctab_Red[] = { 0, 76, 78, 79, 80, 81, 83, 84, 85, 86, 88, 89, 90, 92, 93, 94, 95, 97, 98, 99, 100, 102, 103, 104, 106, 107, 108, 109, 111, 112, 113, 114, 116, 117, 118, 119, 121, 122, 123, 125, 126, 127, 128, 130, 131, 132, 133, 135, 136, 137, 139, 140, 141, 142, 144, 145, 146, 147, 149, 150, 151, 152, 154, 155, 156, 158, 159, 160, 161, 163, 164, 165, 166, 168, 169, 170, 172, 173, 174, 175, 177, 178, 179, 180, 182, 183, 184, 186, 187, 188, 189, 191, 192, 193, 194, 196, 197, 198, 199, 201, 202, 203, 205, 206, 207, 208, 210, 211, 212, 213, 215, 216, 217, 219, 220, 221, 222, 224, 225, 226, 227, 229, 230, 231, 232, 234, 235, 236, 238, 238, 238, 238, 238, 238, 238, 239, 239, 239, 239, 239, 239, 239, 240, 240, 240, 240, 240, 240, 240, 241, 241, 241, 241, 241, 241, 241, 242, 242, 242, 242, 242, 242, 242, 243, 243, 243, 243, 243, 243, 243, 243, 244, 244, 244, 244, 244, 244, 244, 245, 245, 245, 245, 245, 245, 245, 246, 246, 246, 246, 246, 246, 246, 247, 247, 247, 247, 247, 247, 247, 248, 248, 248, 248, 248, 248, 248, 249, 249, 249, 249, 249, 249, 249, 250, 250, 250, 250, 250, 250, 250, 251, 251, 251, 251, 251, 251, 251, 252, 252, 252, 252, 252, 252, 252, 253, 253, 253, 252, 251, 250, 250, 249, 248, 247, 246, 245, 245, 244, 243, 242, 241, 240, 240, 239, 238, 237 };
unsigned char Ctab_Grn[] = { 0, 23, 23, 24, 24, 25, 25, 26, 26, 27, 27, 28, 29, 29, 30, 30, 31, 32, 32, 33, 33, 34, 35, 35, 36, 37, 37, 38, 38, 39, 40, 40, 41, 42, 43, 43, 44, 45, 45, 46, 47, 47, 48, 49, 50, 50, 51, 52, 53, 53, 54, 55, 56, 56, 57, 58, 59, 60, 60, 61, 62, 63, 64, 64, 65, 66, 67, 68, 69, 69, 70, 71, 72, 73, 74, 75, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 124, 125, 126, 127, 128, 129, 130, 132, 133, 134, 135, 136, 137, 139, 140, 141, 142, 143, 144, 146, 147, 148, 149, 150, 151, 153, 154, 155, 156, 157, 158, 160, 161, 162, 163, 164, 165, 167, 168, 169, 170, 171, 172, 174, 175, 176, 177, 178, 179, 181, 182, 183, 184, 185, 186, 188, 189, 190, 191, 192, 193, 194, 196, 197, 198, 199, 200, 201, 203, 204, 205, 206, 207, 208, 209, 211, 212, 213, 214, 215, 216, 217, 219, 220, 221, 222, 223, 224, 225, 227, 228, 229, 230, 231, 232, 233, 235, 236, 237, 238, 239, 240, 241, 242, 244, 245, 246, 247, 248, 249, 250, 251, 253, 253, 253, 253, 253, 254, 254, 254, 254, 254, 254, 254, 255, 255, 255, 255, 255, 255, 255, 255 };
unsigned char Ctab_Blu[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 10, 10, 10, 10, 10, 11, 11, 11, 11, 11, 12, 12, 12, 12, 12, 13, 13, 13, 13, 13, 14, 14, 14, 14, 14, 15, 15, 15, 15, 15, 16, 16, 16, 16, 16, 17, 17, 17, 17, 17, 18, 18, 18, 18, 19, 19, 19, 19, 19, 20, 20, 20, 20, 20, 21, 21, 21, 21, 21, 22, 22, 22, 22, 23, 23, 23, 23, 23, 24, 24, 24, 24, 24, 25, 25, 25 };

/*
 *   Routine to write out a GIF image
 */
void DumpGif3( char *FilenameBase, int W, int H, double **Image ){
double           Val, Val2, Min, Max, dVal;
    int              w, h;
    unsigned char   *uImage, uVal, Filename[1024];
    FILE            *fp_gif, *fp_info;

    int             LogScale;

    LogScale = TRUE;
    LogScale = FALSE;


    // Determine Min/Max values...
    Min =  9e99;
    Max = -9e99;
    for ( w=0; w<W; w++ ){
        for ( h=0; h<H; h++ ) {

            Val = Image[h][w];
            if ( LogScale ) {
                Val2 = (Val > 0.0) ? log10( Val ) : -9e99;
            } else {
                Val2 = Image[h][w];
            }
            if (Val2 > Max) Max = Val2;
            if ((Val2 < Min)&&(Val > 0.0)) Min = Val2;

        }
    }

    printf("Min, Max = %g %g\n", Min, Max);
//Min = 1.9;
//Max = 3.5;

//Min = 1.1;
//Max = 2.0;

    sprintf( Filename, "%s.info", FilenameBase);
    fp_info = fopen( Filename, "w" );
    fprintf( fp_info, "Min: %g\n", Min );
    fprintf( fp_info, "Max: %g\n", Max );
    fclose( fp_info );



    uImage = (unsigned char *)calloc( W*H, sizeof(unsigned char) );

    for ( w=0; w<W; w++ ){
        for ( h=0; h<H; h++ ) {

            if ( LogScale ) {
                Val = Image[h][w] > 0.0 ? log10( Image[h][w] ) : -1e31;
            } else {
                Val = Image[h][w];
            }
            if ( Val < Min ) {
                uVal = 0;
            } else {
                dVal = (Val - Min)/(Max-Min)*255.0;
                uVal = (dVal > 255.0) ? 255 : (unsigned char)dVal;
            }

            *(uImage + W*(H-1-h) + w) = uVal;

        }
    }

    sprintf( Filename, "%s.gif", FilenameBase);
    fp_gif = fopen(Filename, "w");
    WriteGIF(fp_gif, (byte *)uImage, 0, W, H, Ctab_Red, Ctab_Grn, Ctab_Blu, 256, 0, "");
    fclose(fp_gif);

    free( uImage );


    // dump a colorbar image
    W = 10; H = 256;
    uImage = (unsigned char *)calloc( W*H, sizeof(unsigned char) );
    for ( w=0; w<W; w++ ){
        for ( h=0; h<H; h++ ) {
            *(uImage + W*(H-1-h) + w) = h;
        }
    }
    sprintf( Filename, "%s_Bar.gif", FilenameBase);
    fp_gif = fopen(Filename, "w");
    WriteGIF(fp_gif, (byte *)uImage, 0, W, H, Ctab_Red, Ctab_Grn, Ctab_Blu, 256, 0, "");
    fclose(fp_gif);
    free( uImage );

}

int GenIllumTextures( int N, float **FdImage, float **FsImage ) {

    int     i, j;
    double  CosAlpha_Min, CosAlpha_Max, s1;
    double  CosBeta_Min, CosBeta_Max;
    double  LT_Min, LT_Max, s2;
    double  CosAlpha, SinAlpha, Alpha, LT, n, Fd;
    double  CosBeta, Beta, Fs;
    double  **Image;
    double  T, a, b;
    double  epsabs, epsrel, abserr, work[2002];
    int     limit=500, lenw=4*limit, iwork[502], last, ier, neval;
    int     VerbosityLevel = 0;
    Fs_Info fi;

    float   *FsData;
    float   *FdData;

    FsData = (float *)calloc( N*N, sizeof( float) );
    FdData = (float *)calloc( N*N, sizeof( float) );


    n = 30.0;

    CosAlpha_Min = -1.0;
    CosAlpha_Max =  1.0;
    s1 = (CosAlpha_Max - CosAlpha_Min)/(double)(N-1);

    LT_Min = -1.0;
    LT_Max =  1.0;
    s2 = (LT_Max - LT_Min)/(double)(N-1);

    LGM_ARRAY_2D( Image, N, N, double );

    for ( j=0; j<N; ++j ){

        CosAlpha = s1*j + CosAlpha_Min;
        SinAlpha = sqrt( 1.0 - CosAlpha*CosAlpha );
        Alpha = acos( CosAlpha );

        for ( i=0; i<N; ++i ){

            LT = s2*i + LT_Min;

            Fd = 0.25*sqrt( 1.0-LT*LT ) * ( SinAlpha + (M_PI-Alpha)*CosAlpha );
            Image[i][j] = Fd;


        }
    }
//DumpGif3( "Fd", N, N, Image );

    for ( j=0; j<N; ++j ){
        for ( i=0; i<N; ++i ){
            *(FdData + N*i + j) = (float)Image[i][j];
        }
    }
    



    CosBeta_Min = -1.0;
    CosBeta_Max =  1.0;
    s2 = (CosBeta_Max - CosBeta_Min)/(double)(N-1);

    for ( j=0; j<N; ++j ){

        CosAlpha = s1*j + CosAlpha_Min;
        SinAlpha = sqrt( 1.0 - CosAlpha*CosAlpha );
        Alpha = acos( CosAlpha );

        for ( i=0; i<N; ++i ){

            CosBeta = s2*i + CosBeta_Min;
            Beta = acos( CosBeta );


            a = Alpha - M_PI/2.0;
            b = M_PI/2.0;
            epsabs = 0.0;
            epsrel = 1e-5;

            dqags( Fs_Integrand, (_qpInfo *)&fi, a, b, epsabs, epsrel, &Fs, &abserr, &neval, &ier, limit, lenw, &last, iwork, work, VerbosityLevel );


            Image[i][j] = Fs;


        }
    }

    for ( j=0; j<N; ++j ){
        for ( i=0; i<N; ++i ){
            *(FsData + N*i + j) = (float)Image[i][j];
        }
    }






    *FdImage = FdData;
    *FsImage = FsData;
printf("Here\n");
    LGM_ARRAY_2D_FREE( Image );

    return( 1 );


}
