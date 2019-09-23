#include <check.h>
#include "../libLanlGeoMag/Lgm/Lgm_CTrans.h"
#include "../libLanlGeoMag/Lgm/Lgm_QinDenton.h"
#include "../libLanlGeoMag/Lgm/Lgm_MagModelInfo.h"
#include <stdio.h>
#include <stdlib.h>
#define TRUE    1
#define FALSE   0

/*
 *  Unit and Regression tests for FluxTubeVolume solvers
 */


START_TEST(test_FluxTubeVolume_01) {

    double              UTC, JD;
    long int            Date;
    int                 Passed, Flag, nDivs, Status;
    double              V, V_expected, GeodHeight, r, mlat, MLT, cl, sl, Phi;
    Lgm_QinDentonOne    p;
    Lgm_Vector          u, u_sm, v1, v2, v3;
    Lgm_MagModelInfo    *mInfo = Lgm_InitMagInfo();
    FILE                *fp_expected;
    FILE                *fp_got;

    if ( (fp_expected = fopen( "check_FluxTubeVolume_01.expected", "r" )) != NULL ) {

        // real and imaginary parts of solution
        fscanf( fp_expected, "%lf %lf %lf", &mlat, &MLT, &V_expected );
        fclose( fp_expected );
    } else {
        printf("Cant open file: check_FluxTubeVolume_01.expected\n" );
    }

    Date       = 20130324;
    UTC        = 23.662;
    JD = Lgm_Date_to_JD( Date, UTC, mInfo->c );
    Lgm_Set_Coord_Transforms( Date, UTC, mInfo->c );
    printf("Tilt = %g\n", mInfo->c->psi*DegPerRad );
    //exit(0);

    //Lgm_get_QinDenton_at_JD( JD, &p, 1, 1 );
    //Lgm_set_QinDenton( &p, mInfo );
    mInfo->Kp = 4;
    Lgm_MagModelInfo_Set_MagModel( LGM_IGRF, LGM_EXTMODEL_T89, mInfo );
    

    GeodHeight = 120.0;
    r    = 1.0 + GeodHeight/Re;
    mlat = 70.0;
    MLT  = 3.2;
    Phi  = 15.0*(MLT-12.0)*RadPerDeg;
    cl = cos( mlat * RadPerDeg ); sl = sin( mlat * RadPerDeg );
    u_sm.x = r*cl*cos(Phi); u_sm.y = r*cl*sin(Phi); u_sm.z = r*sl;
    Lgm_Convert_Coords( &u_sm, &u, SM_TO_GSM, mInfo->c );
    Flag = Lgm_Trace( &u, &v1, &v2, &v3, GeodHeight, 1e-7, 1e-7, mInfo );
    
    if ( Flag == LGM_CLOSED ) {
        nDivs = mInfo->Stotal/0.1;
        if ( nDivs < 200 ) nDivs = 200;
        mInfo->Hmax = mInfo->Stotal/((double)(nDivs));

        Status = Lgm_TraceLine3( &v1, mInfo->Stotal, nDivs, 1.0, 1e-7, FALSE, mInfo );

        if ( !InitSpline( mInfo ) ) {

            printf("Failed to init spline\n");
            exit(1);

        } else {

            V = Lgm_FluxTubeVolume( mInfo );
            printf( "mlat, MLT, FluxTubeVolume = %g %g %g\n", mlat, MLT, V );
            FreeSpline( mInfo );

        }
    } else {
        printf("Field Line not closed?\n");
    }

    Lgm_FreeMagInfo( mInfo );

    Passed = FALSE;
    if (  (fabs( V_expected - V ) < 1e-8 )  ) Passed = TRUE;

    if ( !Passed ){
        printf("\nTest 01,  Lgm_FluxTubeVolume(): mlat, MLT = %g %g\n", mlat, MLT );
        printf("       Expected V: %.15lf\n", V_expected );
        printf("            Got V: %.15lf\n", V );
        printf("             Diff: %.15g\n\n\n", fabs( V_expected - V ) );
    }

    if ( (fp_got = fopen( "check_FluxTubeVolume_01.got", "w" )) != NULL ) {
        fprintf( fp_got, "%g %g %.15lf\n", mlat, MLT, V );
        fclose( fp_got );
    } else {
        printf("Could not open file: check_FluxTubeVolume_01.got\\n");
    }

    fail_unless( Passed, "Lgm_FluxTubeVolume(): Regression test failed. Compare 'got' and 'expected' files: check_FluxTubeVolume_01.got check_FluxTubeVolume_01.expected\n" );


    return;
}
END_TEST



Suite *FluxTubeVolume_suite(void) {

  Suite *s = suite_create("FLUXTUBEVOLUME_TESTS");

  TCase *tc_FluxTubeVolume = tcase_create("FluxTubeVolume");

  tcase_add_test(tc_FluxTubeVolume, test_FluxTubeVolume_01);

  suite_add_tcase(s, tc_FluxTubeVolume);

  return s;

}

int main(void) {

    int      number_failed;
    Suite   *s  = FluxTubeVolume_suite();
    SRunner *sr = srunner_create(s);

    printf("\n\n");
    srunner_run_all(sr, CK_NORMAL);
    number_failed = srunner_ntests_failed(sr);
    srunner_free(sr);

    return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;

}
