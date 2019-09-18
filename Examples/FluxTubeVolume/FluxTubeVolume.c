#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <Lgm_CTrans.h>
#include <Lgm_QinDenton.h>
#include <Lgm_ElapsedTime.h>
#include <Lgm_FluxToPsd.h>
#include <omp.h>


int main(){

    double              x, y, GeodLat, GeodLong, GeodHeight ;
    double              UTC, JD;
    long int            Date;
    int                 Flag;
    double              V, r, mlat, MLT, cl, sl, Phi;


    double              Hdid, Hnext, s, Inc;
    int                 reset=0, Ni;
    Lgm_Vector          u_scale;
    int                 EnhancedFlag;
    Lgm_QinDentonOne    p;
    Lgm_Vector          u, u_sm, v1, v2, v3, Bvec;
    Lgm_MagModelInfo    *mInfo2, *mInfo = Lgm_InitMagInfo();
    char    Filename[1024];
    int     i, ii, nDivs, Status;
    FILE    *fpout;
    FILE    *fpout2 = fopen( "FL.txt", "w" );
    

    Date       = 20130324;
    UTC        = 23.462;
    UTC        = 23.662;
    JD = Lgm_Date_to_JD( Date, UTC, mInfo->c );
    Lgm_Set_Coord_Transforms( Date, UTC, mInfo->c );
    printf("Tilt = %g\n", mInfo->c->psi*DegPerRad );
    //exit(0);

    Lgm_get_QinDenton_at_JD( JD, &p, 1, 1 );
    Lgm_set_QinDenton( &p, mInfo );

    Lgm_MagModelInfo_Set_MagModel( LGM_IGRF, LGM_EXTMODEL_TS04, mInfo );
    Lgm_MagModelInfo_Set_MagModel( LGM_CDIP, LGM_EXTMODEL_T87, mInfo );
    Lgm_MagModelInfo_Set_MagModel( LGM_IGRF, LGM_EXTMODEL_T96, mInfo );
    Lgm_MagModelInfo_Set_MagModel( LGM_IGRF, LGM_EXTMODEL_T89, mInfo );
    Lgm_MagModelInfo_Set_MagModel( LGM_IGRF, LGM_EXTMODEL_TS04, mInfo );
    Lgm_MagModelInfo_Set_MagModel( LGM_IGRF, LGM_EXTMODEL_OP77, mInfo );
    Lgm_MagModelInfo_Set_MagModel( LGM_IGRF, LGM_EXTMODEL_T01S, mInfo );

    char *s1, *s2, *s3, *s4;
    Lgm_Get_ExtMagModelStrings( &s1, &s2, &s3, &s4, mInfo );
    printf( "s1 = %s\n", s1 );
    printf( "s2 = %s\n", s2 );
    printf( "s3 = %s\n", s3 );
    printf( "s4 = %s\n", s4 );
    sprintf( Filename, "V_%s.txt", mInfo->ExtMagModelStr1 );
    fpout = fopen( Filename, "w" );
    

    //mInfo->fp = fpout2;
    mInfo->SavePoints = FALSE;
    MLT = 0.0;
    r = 1.0 + 200.0/Re;
    GeodHeight = 200.0;
    Inc = 0.001;
    Ni = (80.0-60.0)/Inc;
mInfo->Lgm_MagStep_BS_atol = 1e-6;
mInfo->Lgm_MagStep_BS_rtol = 0.0;
mInfo->VerbosityLevel = 0;

    #pragma omp parallel private(Status,mlat,Phi,cl,sl,u_sm,u,v1,v2,v3,mInfo2,Flag,nDivs,V)
    #pragma omp for schedule(dynamic, 20)
    for ( i=0; i<=Ni; i++ ) {

        mlat = 60.0 + Inc*i;

        mInfo2 = Lgm_CopyMagInfo( mInfo );

        /*
         * Trace from given SM position
         */
        //printf("mlat, MLT = %g %g     ", mlat, MLT);
        Phi = 15.0*(MLT-12.0)*RadPerDeg;
        cl = cos( mlat * RadPerDeg ); sl = sin( mlat * RadPerDeg );
        u_sm.x = r*cl*cos(Phi); u_sm.y = r*cl*sin(Phi); u_sm.z = r*sl;
        Lgm_Convert_Coords( &u_sm, &u, SM_TO_GSM, mInfo2->c );
        //printf("u_sm  = %g %g %g\n", u_sm.x, u_sm.y, u_sm.z );
        //printf("u     = %g %g %g\n", u.x, u.y, u.z );

        Flag = Lgm_Trace( &u, &v1, &v2, &v3, GeodHeight, 1e-7, 1e-7, mInfo2 );


        if ( Flag == LGM_CLOSED ) {
            nDivs = mInfo2->Stotal/0.1;
            if ( nDivs < 200 ) nDivs = 200;
            mInfo2->Hmax = mInfo2->Stotal/((double)(nDivs));

            Status = Lgm_TraceLine3( &v1, mInfo2->Stotal, nDivs, 1.0, 1e-7, FALSE, mInfo2 );

            if ( !InitSpline( mInfo2 ) ) {
                printf("Failed to init spline\n");
                exit(1);
            } else {

                V = FluxTubeVolume( mInfo2 );
                printf( "mlat, MLT, FluxTubeVolume = %g %g %g\n", mlat, MLT, V );
                fprintf( fpout, "%g %g\n", mlat, V );
                fflush( fpout );
                FreeSpline( mInfo2 );

            }
        } else {
            printf("\n");
        }
        
        Lgm_FreeMagInfo( mInfo2 );


    }

    fclose( fpout );
    fclose( fpout2 );
    Lgm_FreeMagInfo( mInfo );




}
