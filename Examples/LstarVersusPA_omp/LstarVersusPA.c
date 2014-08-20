#define MAIN
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <Lgm_CTrans.h>
#include <Lgm_MagEphemInfo.h>

#define KP_DEFAULT 0

void LstarVersusPA( long int Date, double UTC, Lgm_Vector *u, int nAlpha, double *Alpha, int Quality, Lgm_MagEphemInfo *MagEphemInfo );

int main( int argc, char *argv[] ){

    double           UTC, Alpha[1000], a;
    long int         Date;
    int              nAlpha, Kp, i;
    char             Filename2[1024];
    Lgm_Vector       Psm, P;
    Lgm_CTrans       *c = Lgm_init_ctrans(0);
    Lgm_MagModelInfo *mmi = Lgm_InitMagInfo();
    Lgm_MagEphemInfo *MagEphemInfo = Lgm_InitMagEphemInfo(0, 20);
    FILE             *fp;

    mmi->Lgm_MagStep_Integrator = LGM_MAGSTEP_ODE_BS;

    // Date and UTC
    Date       = 20000701;
    Date       = 20020901;
    UTC        = 4.0 + 0.0/60.0 + 0.0/3600.0;
    Lgm_Set_Coord_Transforms( Date, UTC, c );

    // Position in SM
//    Psm.x = -6.6; Psm.y = 0.0; Psm.z = 0.0;
//    Lgm_Convert_Coords( &Psm, &P, SM_TO_GSM, c );

MagEphemInfo->LstarInfo->LSimpleMax = 20.0;
MagEphemInfo->LstarInfo->mInfo->Lgm_MagStep_Integrator = LGM_MAGSTEP_ODE_BS;

P.x = -18.53;
P.x = 10.815;
P.y = 0.0;
P.z = 0.0;
double xxx;
for (xxx=0.0; xxx<3.0; xxx += 0.001 ){





    // Create array of Pitch Angles to compute
    for (nAlpha=0,a=50.0; a<=90.0; a+=10.0,++nAlpha) {
        Alpha[nAlpha] = a ;
        //printf("Alpha[%d] = %g\n", nAlpha, Alpha[nAlpha]);
    }
nAlpha = 1;
Alpha[0] = 90.0;
Alpha[1] = 68.0;
Alpha[2] = 15.0;
Alpha[3] = 2.0;


//P.x =  2.55278;
//P.y = -4.77644;
//P.z =  3.4731;
//Alpha[0] =  81.7844;




//USER INPUT STUFF
    MagEphemInfo->LstarQuality   = 0;
    MagEphemInfo->SaveShellLines = TRUE;
    MagEphemInfo->LstarInfo->VerbosityLevel = 1;
    MagEphemInfo->LstarInfo->mInfo->VerbosityLevel = 0;

    Kp = 2;
    MagEphemInfo->LstarInfo->mInfo->Bfield        = Lgm_B_cdip;
    MagEphemInfo->LstarInfo->mInfo->Bfield        = Lgm_B_T89;
    MagEphemInfo->LstarInfo->mInfo->InternalModel = LGM_CDIP;
    MagEphemInfo->LstarInfo->mInfo->InternalModel = LGM_IGRF;
    MagEphemInfo->LstarInfo->mInfo->Kp = ( Kp >= 0 ) ? Kp : KP_DEFAULT;
    if ( MagEphemInfo->LstarInfo->mInfo->Kp > 5 ) MagEphemInfo->LstarInfo->mInfo->Kp = 5;

MagEphemInfo->LstarInfo->mInfo->Bfield        = Lgm_B_OP77;
MagEphemInfo->LstarInfo->mInfo->Bfield        = Lgm_B_cdip;
MagEphemInfo->LstarInfo->mInfo->Bfield        = Lgm_B_T89c;

Lgm_Vector v1, v2, v3;
Lgm_Set_Coord_Transforms( Date, UTC, MagEphemInfo->LstarInfo->mInfo->c );
Lgm_Trace( &P, &v1, &v2, &v3, 120.0, 0.01, 1e-8, MagEphemInfo->LstarInfo->mInfo);
MagEphemInfo->LstarInfo->mInfo->Bm = MagEphemInfo->LstarInfo->mInfo->Bmin;
printf("ZZZZZZZZZ: MagEphemInfo->LstarInfo->mInfo->Bm = %g\n", MagEphemInfo->LstarInfo->mInfo->Bm);
//v3 = P;

    /*
     * Compute L*s, Is, Bms, Footprints, etc...
     * These quantities are stored in the MagEphemInfo Structure
     */
//P.x = -6.6; P.y = P.z = 0.0;
printf("P = %g %g %g\n", P.x, P.y, P.z);
//    ComputeLstarVersusPA( Date, UTC, &P, nAlpha, Alpha, MagEphemInfo->LstarQuality, MagEphemInfo );
Lgm_Convert_Coords( &P, &v3, SM_TO_GSM, MagEphemInfo->LstarInfo->mInfo->c );




    double Kin = 1.8732;
    Lgm_DateTime DT_UTC;
    Lgm_Make_UTC( Date, UTC, &DT_UTC, c );
    if ( Lgm_Setup_AlphaOfK( &DT_UTC, &v3, MagEphemInfo->LstarInfo->mInfo ) > 0 ) {
        Alpha[0] = Lgm_AlphaOfK( Kin, MagEphemInfo->LstarInfo->mInfo );
        Lgm_TearDown_AlphaOfK( MagEphemInfo->LstarInfo->mInfo );
    } else {
        Alpha[0] = LGM_FILL_VALUE;
    }










    Lgm_ComputeLstarVersusPA( Date, UTC, &v3, nAlpha, Alpha, TRUE, MagEphemInfo );                           
//for (i=0; i<24; i++)printf("%g ", MagEphemInfo->I[i] - MagEphemInfo->I[0]);
//printf("\n");





P.x += 0.0001;
}
WriteMagEphemInfoStruct( "test.dat", nAlpha, MagEphemInfo );

    fp = fopen("Lstar.dat", "w");
    for ( i=0; i<nAlpha; ++i ) {
        fprintf( fp, "%g %.10lf\n", Alpha[i], MagEphemInfo->Lstar[i] );
    }
    fclose(fp);

    Lgm_free_ctrans( c );
    Lgm_FreeMagEphemInfo( MagEphemInfo );

    return(0);

}





