#define MAIN
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <Lgm_CTrans.h>
#include <Lgm_MagEphemInfo.h>

#define KP_DEFAULT 0

void LstarVersusPA( long int Date, double UTC, Lgm_Vector *u, int nAlpha, double *Alpha, int Quality, Lgm_MagEphemInfo *MagEphemInfo );


double LS( long int Date, double Time, double Kin, Lgm_Vector *P, Lgm_MagEphemInfo *MagEphemInfo  ) {

    int             nAlpha = 1;
    double          Alpha[1];
    Lgm_Vector      v1, v2, v3;
    Lgm_DateTime    DT_UTC;

    Lgm_Set_Coord_Transforms( Date, Time, MagEphemInfo->LstarInfo->mInfo->c );
    Lgm_Make_UTC( Date, Time, &DT_UTC, MagEphemInfo->LstarInfo->mInfo->c );

    Lgm_Trace( P, &v1, &v2, &v3, 120.0, 0.001, 1e-8, MagEphemInfo->LstarInfo->mInfo );
    MagEphemInfo->LstarInfo->mInfo->Bm = MagEphemInfo->LstarInfo->mInfo->Bmin;
    Lgm_Convert_Coords( P, &v3, SM_TO_GSM, MagEphemInfo->LstarInfo->mInfo->c );

    if ( Lgm_Setup_AlphaOfK( &DT_UTC, &v3, MagEphemInfo->LstarInfo->mInfo ) > 0 ) {
        Alpha[0] = Lgm_AlphaOfK( Kin, MagEphemInfo->LstarInfo->mInfo );
        Lgm_TearDown_AlphaOfK( MagEphemInfo->LstarInfo->mInfo );
    } else {
        Alpha[0] = LGM_FILL_VALUE;
    }

    Lgm_ComputeLstarVersusPA( Date, Time, &v3, nAlpha, Alpha, TRUE, MagEphemInfo );                           

    return( MagEphemInfo->Lstar[0] );

}

int main( int argc, char *argv[] ){

    double           UTC, Alpha[1000], a;
    long int         Date;
    int              nAlpha, Kp, i;
    Lgm_Vector       P;
    Lgm_CTrans       *c = Lgm_init_ctrans(0);
    Lgm_MagModelInfo *mmi = Lgm_InitMagInfo();
    Lgm_MagEphemInfo *MagEphemInfo = Lgm_InitMagEphemInfo(0, 20);
    FILE             *fp;

    mmi->Lgm_MagStep_Integrator = LGM_MAGSTEP_ODE_BS;

    // Date and UTC
    Date       = 20000701;
    Date       = 20130327;
    Date       = 20130727;
    UTC        = 0.0 + 0.0/60.0 + 0.0/3600.0;
    Lgm_Set_Coord_Transforms( Date, UTC, c );
    MagEphemInfo->LstarInfo->LSimpleMax = 25.0;
    MagEphemInfo->LstarInfo->mInfo->Lgm_MagStep_Integrator = LGM_MAGSTEP_ODE_BS;
    MagEphemInfo->LstarInfo->nFLsInDriftShell = 48;

    //USER INPUT STUFF
    MagEphemInfo->LstarQuality   = 3;
    MagEphemInfo->SaveShellLines = TRUE;
    MagEphemInfo->LstarInfo->VerbosityLevel = 3;
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
    MagEphemInfo->LstarInfo->mInfo->InternalModel = LGM_CDIP;
    MagEphemInfo->LstarInfo->mInfo->InternalModel = LGM_IGRF;


    FILE         *fpjunk;
    int          done;
    char         fname[1024];
    double       Ra, Rc, Fa, Fc, R, F, Phi, Lat, MLT, Fmax, Rmax;
    double       Kin = 2.0;

    sprintf(fname, "junk_K_v5_24.txt" );
    fpjunk = fopen(fname, "w");





    for (MLT=0.00; MLT<=-1.0; MLT+=0.05){


        done = FALSE;

        Ra  = 6.0;   
        Phi = (MLT*15.0 - 180.0)*RadPerDeg; Lat = 0.0*RadPerDeg;
        P.x = Ra*cos( Phi )*cos(Lat); P.y = Ra*sin( Phi )*cos(Lat); P.z = Ra*sin(Lat);
        Fa  =  LS(  Date, UTC, Kin, &P, MagEphemInfo  );

        Rc  = 20.0;   
        Phi = (MLT*15.0 - 180.0)*RadPerDeg; Lat = 0.0*RadPerDeg;
        P.x = Rc*cos( Phi )*cos(Lat); P.y = Rc*sin( Phi )*cos(Lat); P.z = Rc*sin(Lat);
        Fc  =  LS(  Date, UTC, Kin, &P, MagEphemInfo  );

        Fmax = -1.0;

        if ( (Fc < 0.0) && (Fa>0.0) ) {

            while (!done ) {

                R   = (Ra+Rc)/2.0;
                Phi = (MLT*15.0 - 180.0)*RadPerDeg; Lat = 0.0*RadPerDeg;
                P.x = R*cos( Phi )*cos(Lat); P.y = R*sin( Phi )*cos(Lat); P.z = R*sin(Lat);
                F  =  LS(  Date, UTC, Kin, &P, MagEphemInfo  );

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



MLT = 0.0;
Rmax= 4.1;
Kin = 0.1;
    Phi = (MLT*15.0 - 180.0)*RadPerDeg; Lat = 0.0*RadPerDeg;
    P.x = Rmax*cos( Phi )*cos(Lat); P.y = Rmax*sin( Phi )*cos(Lat); P.z = Rmax*sin(Lat);
    printf( "Final = %g\n", LS(  Date, UTC, Kin, &P, MagEphemInfo  ));
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





