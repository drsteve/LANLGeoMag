#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <Lgm_CTrans.h>
#include <Lgm_Sgp.h>
#include <Lgm_MagEphemInfo.h>
#include <Lgm_DynamicMemory.h>
#include <Lgm_ElapsedTime.h>
#include <Lgm_Eop.h>
#include <ctype.h>


#define KP_DEFAULT 0



/*
 *  Compute position of S/C from TLE (Two Line Elements).
 *  Then use this to compute Magnetic Ephemerides...
 */

int main( int argc, char *argv[] ){

    long int        Date;
    double          UTC;
    double          r, lat, lon;
    int             i;
    Lgm_CTrans      *c = Lgm_init_ctrans( 0 );
    Lgm_Vector      Ugsm, Ugeo;
    char            IntModel[20], ExtModel[20];
    int             Colorize;
    FILE            *fp_MagEphem;


    double           Alpha[1000], Kp, Dst;
    int              nAlpha;
    Lgm_MagEphemInfo *MagEphemInfo;
    Lgm_ElapsedTimeInfo  tInfo;

    Lgm_ElapsedTimeInit( &tInfo, 255, 95, 0 );

    nAlpha   = 36;
    for (i=0; i<nAlpha; i++) {
        Alpha[i] = 90.0 -  i*2.5;
        printf("Alpha = %lf\n", Alpha[i]);
    }

    if ( nAlpha > 0 ){
        MagEphemInfo = Lgm_InitMagEphemInfo(0, nAlpha);
    } else {
        // doesnt seem to like allocating zero size...
        MagEphemInfo = Lgm_InitMagEphemInfo(0, 1);
    }

    // Settings for Lstar calcs
    MagEphemInfo->ComputeVgc = TRUE;
    MagEphemInfo->ComputeVgc = FALSE;
    MagEphemInfo->LstarQuality = 3;
    MagEphemInfo->SaveShellLines = TRUE;
    MagEphemInfo->LstarInfo->LSimpleMax = 20.0;
    MagEphemInfo->LstarInfo->VerbosityLevel = 4;
    MagEphemInfo->LstarInfo->mInfo->VerbosityLevel = 2;
    MagEphemInfo->LstarInfo->mInfo->Lgm_LossConeHeight = 100.0;


//    MagEphemInfo->LstarInfo->mInfo->Bfield = Lgm_B_TS04;
    MagEphemInfo->LstarInfo->mInfo->Bfield = Lgm_B_TS04_opt;
//MagEphemInfo->LstarInfo->mInfo->Bfield = Lgm_B_OP77;
//MagEphemInfo->LstarInfo->mInfo->Bfield = Lgm_B_cdip;
//MagEphemInfo->LstarInfo->mInfo->Bfield = Lgm_B_T89;
    MagEphemInfo->LstarInfo->mInfo->P      = 4.1011111111111118;
    MagEphemInfo->LstarInfo->mInfo->Dst    = 7.7777777777777777;
    MagEphemInfo->LstarInfo->mInfo->By     = 3.7244444444444444;
    MagEphemInfo->LstarInfo->mInfo->Bz     = -0.12666666666666665;
    MagEphemInfo->LstarInfo->mInfo->W[0]   = 0.12244444444444445;
    MagEphemInfo->LstarInfo->mInfo->W[1]   = 0.2514;
    MagEphemInfo->LstarInfo->mInfo->W[2]   = 0.089266666666666661;
    MagEphemInfo->LstarInfo->mInfo->W[3]   = 0.047866666666666668;
    MagEphemInfo->LstarInfo->mInfo->W[4]   = 0.22586666666666666;
    MagEphemInfo->LstarInfo->mInfo->W[5]   = 1.0461333333333334;

    //vals = magcoords.Lvalue([-4, 0, 1], datetime.datetime(2009,1,1))
    Date = 19960106;
    //UTC = 1.0 + 40.0/60.0 + 59.145/3600.0;
    UTC = 1.2444444444444445;




    Colorize = TRUE;
//    MagEphemInfo->LstarInfo->mInfo->Kp = Kp;
    r = 4.83415065;
    lon = -40.26632902*RadPerDeg;
    lat = 36.44696388*RadPerDeg;
    Ugsm.x = r*cos(lat)*cos(lon);
    Ugsm.y = r*cos(lat)*sin(lon);
    Ugsm.z = r*sin(lat);
    printf("Ugsm = %.15lf %.15lf %.15lf\n", Ugsm.x, Ugsm.y, Ugsm.z);



    MagEphemInfo->nAlpha = nAlpha;
    for (i=0; i<nAlpha; i++){
        MagEphemInfo->Alpha[i] = Alpha[i];
        printf("Alpha[%d] = %g\n", i, Alpha[i]);
    }

    /*
     * Open Mag Ephem file for writing
     */
    fp_MagEphem = fopen( "OutputFile.txt", "wb" );
    Lgm_WriteMagEphemHeader( fp_MagEphem, "", 0, "", IntModel, ExtModel, MagEphemInfo );



    /*
     * Compute L*s, Is, Bms, Footprints, etc...
     * These quantities are stored in the MagEphemInfo Structure
     */
    printf("\n\n\nDate, ut = %ld %g   Ugsm = %g %g %g \n", Date, UTC, Ugsm.x, Ugsm.y, Ugsm.z );
    Lgm_ComputeLstarVersusPA( Date, UTC, &Ugsm, nAlpha, Alpha, MagEphemInfo->LstarQuality, Colorize, MagEphemInfo );
    Lgm_WriteMagEphemData( fp_MagEphem, IntModel, ExtModel, Kp, Dst, MagEphemInfo );
    Lgm_PrintElapsedTime( &tInfo );


    fclose(fp_MagEphem);
    Lgm_free_ctrans( c );
    Lgm_FreeMagEphemInfo( MagEphemInfo );

    char Filename[256];
    FILE *fp;
    Lgm_SetCurrentTimeStr( &tInfo );
    sprintf( Filename, "MagEphemFromPosition_Timing_%s_%s.txt", getenv("HOSTNAME"), tInfo.CurrentTimeStr2 );
    fp = fopen( Filename, "w");
    fprintf( fp, "*****  Elapsed Time (DD:HH:MM:SS): %s  *****", tInfo.ElapsedTimeStr );
    fclose(fp);



    return(0);
}

