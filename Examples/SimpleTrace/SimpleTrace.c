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
    int                 Flag, iMLT;
    double              r, mlat, MLT, cl, sl, Phi;


    double              Hdid, Hnext, s;
    int                 reset=0;
    Lgm_Vector          u_scale;
    int                 EnhancedFlag;
    Lgm_QinDentonOne    p;
    Lgm_Vector          u, u_sm, v1, v2, v3, Bvec;
    Lgm_MagModelInfo    *mInfo = Lgm_InitMagInfo();
    int     i;
    FILE    *fpout = fopen( "FieldLine.txt", "w" );
    FILE    *fpout2 = fopen( "Bmin_Points.txt", "w" );
    

    Date       = 20130410;
    UTC        = 7.0;
//    UTC        = 23.562;
    JD = Lgm_Date_to_JD( Date, UTC, mInfo->c );
    Lgm_Set_Coord_Transforms( Date, UTC, mInfo->c );
    printf("Tilt = %g\n", mInfo->c->psi*DegPerRad );
    //exit(0);

    Lgm_get_QinDenton_at_JD( JD, &p, 1, 1 );
    Lgm_set_QinDenton( &p, mInfo );

    Lgm_MagModelInfo_Set_MagModel( LGM_IGRF, LGM_EXTMODEL_T01S, mInfo );
    Lgm_MagModelInfo_Set_MagModel( LGM_CDIP, LGM_EXTMODEL_TS04, mInfo );
    Lgm_MagModelInfo_Set_MagModel( LGM_CDIP, LGM_EXTMODEL_T89, mInfo );

    char *s1, *s2, *s3, *s4;
    Lgm_Get_ExtMagModelStrings( &s1, &s2, &s3, &s4, mInfo );
    printf( "s1 = %s\n", s1 );
    printf( "s2 = %s\n", s2 );
    printf( "s3 = %s\n", s3 );
    printf( "s4 = %s\n", s4 );
    
    mInfo->fp = fopen( "FL.txt", "w");
    mInfo->SavePoints = TRUE;

    for ( iMLT = 12; iMLT <= 12; iMLT += 12 ) {
        MLT = (double)iMLT;

        for ( mlat = 0.0; mlat <= 90.0; mlat += 0.001 ) {

            /*
             * Trace from given SM position
             */
            printf("mlat, MLT = %g %g\n", mlat, MLT);
            r = 1.0 + 120.0/Re;
            Phi = 15.0*(MLT-12.0)*RadPerDeg;
            cl = cos( mlat * RadPerDeg ); sl = sin( mlat * RadPerDeg );
            u_sm.x = r*cl*cos(Phi); u_sm.y = r*cl*sin(Phi); u_sm.z = r*sl;
            Lgm_Convert_Coords( &u_sm, &u, SM_TO_GSM, mInfo->c );
            printf("u_sm  = %g %g %g\n", u_sm.x, u_sm.y, u_sm.z );
            printf("u     = %g %g %g\n", u.x, u.y, u.z );

            GeodHeight = 120.0;
            //mInfo->VerbosityLevel = 6;
// this just about nails it!
mInfo->Lgm_MagStep_BS_atol = 1e-9;
mInfo->Lgm_MagStep_BS_rtol = 0.0;
mInfo->nFunc = 0;
            Flag = Lgm_Trace( &u, &v1, &v2, &v3, GeodHeight, 1e-7, 1e-7, mInfo );
            printf("Flag = %d, nFunc = %d\n", Flag, mInfo->nFunc);

            mInfo->Bfield( &u, &Bvec, mInfo );
            printf("u  = %g %g %g    |B| = %g\n", u.x, u.y, u.z, Lgm_Magnitude( &Bvec ) );

            mInfo->Bfield( &v1, &Bvec, mInfo );
            printf("v1 = %g %g %g    |B| = %g\n", v1.x, v1.y, v1.z, Lgm_Magnitude( &Bvec ) );

            mInfo->Bfield( &v2, &Bvec, mInfo );
            printf("v2 = %g %g %g    |B| = %g\n", v2.x, v2.y, v2.z, Lgm_Magnitude( &Bvec ) );

            mInfo->Bfield( &v3, &Bvec, mInfo );
            printf("v3 = %g %g %g    |B| = %g\n\n", v3.x, v3.y, v3.z, Lgm_Magnitude( &Bvec ) );

            fprintf( fpout2, "%.10lf %.10lf %.10lf %.10lf %.10lf\n", mInfo->Smin, mInfo->Pmin.x, mInfo->Pmin.y, mInfo->Pmin.z, mInfo->Bmin);
            fprintf( fpout2, "%.10lf %.10lf %.10lf %.10lf %.10lf\n", mInfo->Smin, v3.x, v3.y, v3.z, mInfo->Bmin);

            if ( (Flag == LGM_CLOSED) || (Flag == LGM_OPEN_S_LOBE) ) {

//                Lgm_TraceLine3( &v1, mInfo->Stotal, 4000, 1.0, 1e-7, FALSE, mInfo );
//                for (i=0; i<mInfo->nPnts; i++){ fprintf( fpout, "%.10lf %.10lf %.10lf %.10lf %.10lf\n", mInfo->s[i], mInfo->Px[i], mInfo->Py[i], mInfo->Pz[i], mInfo->Bmag[i]); }
/*
                Lgm_TraceLine3( &v2, mInfo->Stotal, 8000, -1.0, 1e-1, FALSE, mInfo );
                for (i=0; i<mInfo->nPnts; i++){ fprintf( fpout, "%.10lf %.10lf %.10lf %.10lf %.10lf\n", mInfo->s[i], mInfo->Px[i], mInfo->Py[i], mInfo->Pz[i], mInfo->Bmag[i]); }
                Lgm_TraceLine3( &v3, mInfo->Stotal, 8000, -1.0, 1e-1, FALSE, mInfo );
                for (i=0; i<mInfo->nPnts; i++){ fprintf( fpout, "%.10lf %.10lf %.10lf %.10lf %.10lf\n", mInfo->s[i], mInfo->Px[i], mInfo->Py[i], mInfo->Pz[i], mInfo->Bmag[i]); }
*/

            } else if ( Flag == LGM_OPEN_N_LOBE ) {

//                Lgm_TraceLine3( &v2, mInfo->Stotal, 1000, -1.0, 1e-1, TRUE, mInfo );
//                for (i=0; i<mInfo->nPnts; i++){ fprintf( fpout, "%.10lf %.10lf %.10lf %.10lf %.10lf\n", mInfo->s[i], mInfo->Px[i], mInfo->Py[i], mInfo->Pz[i], mInfo->Bmag[i]); }

            } 

        }
        }

        fclose( fpout );
    fclose( fpout2 );




}
