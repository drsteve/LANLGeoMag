#define MAIN
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <Lgm_CTrans.h>
#include <Lgm_MagEphemInfo.h>

#define KP_DEFAULT 0

void LstarVersusPA( long int Date, double UTC, Lgm_Vector *u, int nAlpha, double *Alpha, int Quality, Lgm_MagEphemInfo *MagEphemInfo );


double LS( long int Date, double Time, double Kin, Lgm_Vector *Pgsm, Lgm_MagEphemInfo *MagEphemInfo  ) {

    int             nAlpha = 1, Flag;
    double          Alpha[1];
    Lgm_Vector      v1, v2, v3;
    Lgm_DateTime    DT_UTC;
printf("1. Pgsm = %g %g %g \n", Pgsm->x, Pgsm->y, Pgsm->z);

    Lgm_Set_Coord_Transforms( Date, Time, MagEphemInfo->LstarInfo->mInfo->c );
    Lgm_Make_UTC( Date, Time, &DT_UTC, MagEphemInfo->LstarInfo->mInfo->c );

    Flag = Lgm_Trace( Pgsm, &v1, &v2, &v3, 120.0, 1e-4, 1e-8, MagEphemInfo->LstarInfo->mInfo );
    MagEphemInfo->LstarInfo->mInfo->Bm = MagEphemInfo->LstarInfo->mInfo->Bmin;
printf("2. Pgsm = %g %g %g Flag, Alpha[0] = %d\n", Pgsm->x, Pgsm->y, Pgsm->z, Flag );
printf("3. v1 = %g %g %g\n", v1.x, v1.y, v1.z );
printf("4. v2 = %g %g %g\n", v2.x, v2.y, v2.z );
printf("5. v3 = %g %g %g\n", v3.x, v3.y, v3.z );
    if ( Flag != LGM_CLOSED ) return( LGM_FILL_VALUE );

    if ( Lgm_Setup_AlphaOfK( &DT_UTC, &v3, MagEphemInfo->LstarInfo->mInfo ) > 0 ) {
        Alpha[0] = Lgm_AlphaOfK( Kin, MagEphemInfo->LstarInfo->mInfo );
        Lgm_TearDown_AlphaOfK( MagEphemInfo->LstarInfo->mInfo );
    } else {
        Alpha[0] = LGM_FILL_VALUE;
    }

printf("6. Pgsm = %g %g %g Flag, Alpha[0] = %d %g\n", Pgsm->x, Pgsm->y, Pgsm->z, Flag, Alpha[0]);
    Lgm_ComputeLstarVersusPA( Date, Time, &v3, nAlpha, Alpha, TRUE, MagEphemInfo );                           

printf("\n\n");
    return( MagEphemInfo->Lstar[0] );

}

int main( int argc, char *argv[] ){

    double           UTC, Alpha[1000], a;
    long int         Date;
    int              nAlpha, Kp, i;
    Lgm_Vector       P, Pgsm;
    Lgm_CTrans       *c = Lgm_init_ctrans(0);
    Lgm_MagModelInfo *mmi = Lgm_InitMagInfo();
    Lgm_MagEphemInfo *MagEphemInfo = Lgm_InitMagEphemInfo(0, 20);
    FILE             *fp;

    mmi->Lgm_MagStep_Integrator = LGM_MAGSTEP_ODE_BS;

    // Date and UTC
    Date       = 20130327;
    Date       = 20000701;
    Date       = 20130101;
    UTC        = 0.0 + 0.0/60.0 + 0.0/3600.0;
    Lgm_Set_Coord_Transforms( Date, UTC, c );
    MagEphemInfo->LstarInfo->LSimpleMax = 35.0;
    MagEphemInfo->LstarInfo->mInfo->Lgm_MagStep_Integrator = LGM_MAGSTEP_ODE_BS;
    MagEphemInfo->LstarInfo->nFLsInDriftShell = 96;
    MagEphemInfo->nFLsInDriftShell = 96;

    //USER INPUT STUFF
    MagEphemInfo->LstarQuality   = 3;
Lgm_SetLstarTolerances( 3, 96, MagEphemInfo->LstarInfo );
    MagEphemInfo->SaveShellLines = TRUE;
    MagEphemInfo->LstarInfo->VerbosityLevel = 1;
    MagEphemInfo->LstarInfo->mInfo->VerbosityLevel = 0;

    Kp = 2;
    MagEphemInfo->LstarInfo->mInfo->Bfield        = Lgm_B_cdip;
    MagEphemInfo->LstarInfo->mInfo->Bfield        = Lgm_B_T89;
    MagEphemInfo->LstarInfo->mInfo->Kp = ( Kp >= 0 ) ? Kp : KP_DEFAULT;
    if ( MagEphemInfo->LstarInfo->mInfo->Kp > 5 ) MagEphemInfo->LstarInfo->mInfo->Kp = 5;

    MagEphemInfo->LstarInfo->mInfo->Bfield        = Lgm_B_cdip;
    MagEphemInfo->LstarInfo->mInfo->Bfield        = Lgm_B_Dungey;
    MagEphemInfo->LstarInfo->mInfo->Bfield        = Lgm_B_OP77;
    MagEphemInfo->LstarInfo->mInfo->Bfield        = Lgm_B_T89c;
    MagEphemInfo->LstarInfo->mInfo->InternalModel = LGM_IGRF;
    MagEphemInfo->LstarInfo->mInfo->InternalModel = LGM_CDIP;


    FILE         *fpjunk;
    int          done;
    char         fname[1024];
    double       Ra, Rc, Fa, Fc, R, F, Phi, Lat, MLT, Fmax, Rmax;
    double       Kin = 0.5;

    sprintf(fname, "junk_K_v5_24.txt" );
    fpjunk = fopen(fname, "w");





    for (MLT=0.0; MLT<=24.0; MLT+=0.05){


        done = FALSE;

        Ra  = 5.0;   
        Phi = (MLT*15.0 - 180.0)*RadPerDeg; Lat = 0.0*RadPerDeg;
        P.x = Ra*cos( Phi )*cos(Lat); P.y = Ra*sin( Phi )*cos(Lat); P.z = Ra*sin(Lat);
        Lgm_Convert_Coords( &P, &Pgsm, SM_TO_GSM, c );
        Fa  =  LS(  Date, UTC, Kin, &Pgsm, MagEphemInfo  );
//    WriteMagEphemInfoStruct( "test.dat", 1, MagEphemInfo );
//exit(0);

        Rc  = 30.0;   
        Phi = (MLT*15.0 - 180.0)*RadPerDeg; Lat = 0.0*RadPerDeg;
        P.x = Rc*cos( Phi )*cos(Lat); P.y = Rc*sin( Phi )*cos(Lat); P.z = Rc*sin(Lat);
        Lgm_Convert_Coords( &P, &Pgsm, SM_TO_GSM, c );
        Fc  =  LS(  Date, UTC, Kin, &Pgsm, MagEphemInfo  );

        Fmax = -1.0;

        if ( (Fc < 0.0) && (Fa>0.0) ) {

            while (!done ) {

                R   = (Ra+Rc)/2.0;
printf("R = %g\n", R);
                Phi = (MLT*15.0 - 180.0)*RadPerDeg; Lat = 0.0*RadPerDeg;
                P.x = R*cos( Phi )*cos(Lat); P.y = R*sin( Phi )*cos(Lat); P.z = R*sin(Lat);
                Lgm_Convert_Coords( &P, &Pgsm, SM_TO_GSM, c );
                F  =  LS(  Date, UTC, Kin, &Pgsm, MagEphemInfo  );

                if ( F < 0.0 ) {
                    Rc = R;
                    Fc = F;
                } else {
                    if ( F > Fmax ) { Fmax = F; Rmax = R; };
                    Ra = R;
                    Fa = F;
                }

                if ( fabs(Rc-Ra) < 1e-3 ) done = TRUE;

            }

        } else {
            printf("No bracket!\n");
        }


if (Fmax > Fa) printf("************************* Fmax, Fa = %g %g\n", Fmax, Fa);
        fprintf( fpjunk, "%g %g\n", MLT, Fa );
        fflush(fpjunk);

    }

    fclose(fpjunk);



    MLT = 3.00;
    Phi = (MLT*15.0 - 180.0)*RadPerDeg; Lat = 0.0*RadPerDeg;
    P.x = Ra*cos( Phi )*cos(Lat); P.y = Ra*sin( Phi )*cos(Lat); P.z = Ra*sin(Lat);
    Lgm_Convert_Coords( &P, &Pgsm, SM_TO_GSM, c );
    F  =  LS(  Date, UTC, Kin, &Pgsm, MagEphemInfo  );
    WriteMagEphemInfoStruct( "test.dat", 1, MagEphemInfo );

    fp = fopen("Lstar.dat", "w");
    for ( i=0; i<nAlpha; ++i ) {
        fprintf( fp, "%g %.10lf\n", Alpha[i], MagEphemInfo->Lstar[i] );
    }
    fclose(fp);

    Lgm_free_ctrans( c );
    Lgm_FreeMagEphemInfo( MagEphemInfo );

    return(0);

}





