#include <check.h>
#include "../libLanlGeoMag/Lgm/Lgm_MagModelInfo.h"

/*
 *  Regression tests for McIlwain L-shell calculations
 */


Lgm_MagModelInfo    *mInfo;
Lgm_CTrans          *c;

void McIlwain_L_Setup(void) {
    mInfo = Lgm_InitMagInfo();
    return;
}

void McIlwain_L_TearDown(void) {
    Lgm_FreeMagInfo( mInfo );
    return;
}


START_TEST(test_MCILWAIN_01) {

    int                 Passed = FALSE;
    long int            Date;
    double              L, I, Bm, M, UTC;
    double              L_expected=-9e99, I_expected=-9e99, Bm_expected=-9e99, M_expected=-9e99;
    Lgm_Vector          u, v;
    FILE                *fp_expected;
    FILE                *fp_got;

    if ( (fp_expected = fopen( "check_McIlwain_L_01.expected", "r" )) != NULL ) {
        fscanf( fp_expected, "%lf", &L_expected );
        fscanf( fp_expected, "%lf", &I_expected );
        fscanf( fp_expected, "%lf", &Bm_expected );
        fscanf( fp_expected, "%lf", &M_expected );
        fclose( fp_expected );
    } else {
        printf("Lgm_McIlwain_L(): Cant open file: check_McIlwain_L_01.expected\n" );
    }


    Date = 20090101; UTC  = 0.0; mInfo->Kp = 2; Lgm_Set_Coord_Transforms( Date, UTC, mInfo->c );
    u.x = -4.0; u.y = 0.0; u.z = 0.0;
    Lgm_Convert_Coords( &u, &v, SM_TO_GSM, mInfo->c );
    Lgm_MagModelInfo_Set_MagModel( LGM_CDIP, LGM_EXTMODEL_NULL, mInfo );
    L = Lgm_McIlwain_L( Date, UTC, &v, 90.0, 0, &I, &Bm, &M, mInfo );

    if (    (fabs( L-L_expected ) < 1e-7) && (fabs( I-I_expected ) < 1e-7) && (fabs( Bm-Bm_expected ) < 1e-7) && (fabs( M-M_expected ) < 1e-7) ) Passed = TRUE;
    if ( !Passed ){
        printf("\nTest 01, Lgm_McIlwain_L(): %15s    %15s    %15s    %15s\n", "       L       ", "       I       ", "       Bm       ", "       M       ");
        printf("                   Expected: %.15g   %.15g   %.15g   %.15g\n", L_expected,  I_expected, Bm_expected, M_expected );
        printf("                        Got: %.15g   %.15g   %.15g   %.15g\n\n\n", L,  I, Bm, M);
    }

    if ( (fp_got = fopen( "check_McIlwain_L_01.got", "w" )) != NULL ) {
        fprintf( fp_got, "%.15lf\n", L );
        fprintf( fp_got, "%.15lf\n", I );
        fprintf( fp_got, "%.15lf\n", Bm );
        fprintf( fp_got, "%.15lf\n", M );
        fclose( fp_got );
    } else {
        printf("Could not open file: check_McIlwain_L_01.got\\n");
    }

    fail_unless( Passed, "Lgm_McIlwain_L(): Regression test failed. Compare 'expected' and 'got' files: check_McIlwain_L_01.expected check_McIlwain_L_01.got\n" );


    return;
}
END_TEST

START_TEST(test_MCILWAIN_02) {

    int                 Passed = FALSE;
    long int            Date;
    double              L, I, Bm, M, UTC;
    double              L_expected=-9e99, I_expected=-9e99, Bm_expected=-9e99, M_expected=-9e99;
    Lgm_Vector          u, v;
    FILE                *fp_expected;
    FILE                *fp_got;

    if ( (fp_expected = fopen( "check_McIlwain_L_02.expected", "r" )) != NULL ) {
        fscanf( fp_expected, "%lf", &L_expected );
        fscanf( fp_expected, "%lf", &I_expected );
        fscanf( fp_expected, "%lf", &Bm_expected );
        fscanf( fp_expected, "%lf", &M_expected );
        fclose( fp_expected );
    } else {
        printf("Lgm_McIlwain_L(): Cant open file: check_McIlwain_L_02.expected\n" );
    }

    Date = 20090101; UTC  = 0.0; mInfo->Kp = 2; Lgm_Set_Coord_Transforms( Date, UTC, mInfo->c );
    u.x = -6.6; u.y = 2.3; u.z = 0.4;
    Lgm_MagModelInfo_Set_MagModel( LGM_IGRF, LGM_EXTMODEL_T89, mInfo );
    L = Lgm_McIlwain_L( Date, UTC, &u, 90.0, 1, &I, &Bm, &M, mInfo );

    if (    (fabs( L-L_expected ) < 1e-7) && (fabs( I-I_expected ) < 1e-7) && (fabs( Bm-Bm_expected ) < 1e-7) && (fabs( M-M_expected ) < 1e-7) ) Passed = TRUE;
    if ( !Passed ){
        printf("\nTest 02, Lgm_McIlwain_L(): %15s    %15s    %15s    %15s\n", "       L       ", "       I       ", "       Bm       ", "       M       ");
        printf("                   Expected: %.15g   %.15g   %.15g   %.15g\n", L_expected,  I_expected, Bm_expected, M_expected );
        printf("                        Got: %.15g   %.15g   %.15g   %.15g\n\n\n", L,  I, Bm, M);
    }

    if ( (fp_got = fopen( "check_McIlwain_L_02.got", "w" )) != NULL ) {
        fprintf( fp_got, "%.15lf\n", L );
        fprintf( fp_got, "%.15lf\n", I );
        fprintf( fp_got, "%.15lf\n", Bm );
        fprintf( fp_got, "%.15lf\n", M );
        fclose( fp_got );
    } else {
        printf("Could not open file: check_McIlwain_L_02.got\\n");
    }

    fail_unless( Passed, "Lgm_McIlwain_L(): Regression test failed. Compare 'expected' and 'got' files: check_McIlwain_L_02.expected check_McIlwain_L_02.got\n" );


    return;
}
END_TEST


START_TEST(test_MCILWAIN_03) {

    int                 Passed = FALSE;
    long int            Date;
    double              L, I, Bm, M, UTC;
    double              L_expected=-9e99, I_expected=-9e99, Bm_expected=-9e99, M_expected=-9e99;
    Lgm_Vector          u, v;
    FILE                *fp_expected;
    FILE                *fp_got;

    if ( (fp_expected = fopen( "check_McIlwain_L_03.expected", "r" )) != NULL ) {

        fscanf( fp_expected, "%lf", &L_expected );
        fscanf( fp_expected, "%lf", &I_expected );
        fscanf( fp_expected, "%lf", &Bm_expected );
        fscanf( fp_expected, "%lf", &M_expected );
        fclose( fp_expected );
    } else {
        printf("Lgm_McIlwain_L(): Cant open file: check_McIlwain_L_03.expected\n" );
    }

    Date = 20060823; UTC  = 13.213; mInfo->Kp = 5; Lgm_Set_Coord_Transforms( Date, UTC, mInfo->c );
    u.x = 2.4; u.y = 1.2; u.z = 0.4;
    Lgm_MagModelInfo_Set_MagModel( LGM_IGRF, LGM_EXTMODEL_T89c, mInfo );
    L = Lgm_McIlwain_L( Date, UTC, &u, 90.0, 1, &I, &Bm, &M, mInfo );

    if (    (fabs( L-L_expected ) < 1e-7) && (fabs( I-I_expected ) < 1e-7) && (fabs( Bm-Bm_expected ) < 1e-7) && (fabs( M-M_expected ) < 1e-7) ) Passed = TRUE;
    if ( !Passed ){
        printf("\nTest 03, Lgm_McIlwain_L(): %15s    %15s    %15s    %15s\n", "       L       ", "       I       ", "       Bm       ", "       M       ");
        printf("                   Expected: %.15g   %.15g   %.15g   %.15g\n", L_expected,  I_expected, Bm_expected, M_expected );
        printf("                        Got: %.15g   %.15g   %.15g   %.15g\n\n\n", L,  I, Bm, M);
    }

    if ( (fp_got = fopen( "check_McIlwain_L_03.got", "w" )) != NULL ) {
        fprintf( fp_got, "%.15lf\n", L );
        fprintf( fp_got, "%.15lf\n", I );
        fprintf( fp_got, "%.15lf\n", Bm );
        fprintf( fp_got, "%.15lf\n", M );
        fclose( fp_got );
    } else {
        printf("Could not open file: check_McIlwain_L_03.got\\n");
    }

    fail_unless( Passed, "Lgm_McIlwain_L(): Regression test failed. Compare 'expected' and 'got' files: check_McIlwain_L_03.expected check_McIlwain_L_03.got\n" );


    return;
}
END_TEST

START_TEST(test_MCILWAIN_04) {

    int                 Passed = FALSE;
    long int            Date;
    double              L, I, Bm, M, UTC;
    double              L_expected=-9e99, I_expected=-9e99, Bm_expected=-9e99, M_expected=-9e99;
    Lgm_Vector          u, v;
    FILE                *fp_expected;
    FILE                *fp_got;

    if ( (fp_expected = fopen( "check_McIlwain_L_04.expected", "r" )) != NULL ) {

        fscanf( fp_expected, "%lf", &L_expected );
        fscanf( fp_expected, "%lf", &I_expected );
        fscanf( fp_expected, "%lf", &Bm_expected );
        fscanf( fp_expected, "%lf", &M_expected );
        fclose( fp_expected );
    } else {
        printf("Lgm_McIlwain_L(): Cant open file: check_McIlwain_L_04.expected\n" );
    }


    Date = 20090101; UTC  = 0.0; mInfo->Kp = 2; Lgm_Set_Coord_Transforms( Date, UTC, mInfo->c );
    u.x = -7.0; u.y = 1.0; u.z = 2.0;
    //Lgm_Convert_Coords( &u, &v, SM_TO_GSM, mInfo->c );
    Lgm_MagModelInfo_Set_MagModel( LGM_IGRF, LGM_EXTMODEL_T89, mInfo );
    L = Lgm_McIlwain_L( Date, UTC, &u, 90.0, 1, &I, &Bm, &M, mInfo );

    if (    (fabs( L-L_expected ) < 1e-7) && (fabs( I-I_expected ) < 1e-7) && (fabs( Bm-Bm_expected ) < 1e-7) && (fabs( M-M_expected ) < 1e-7) ) Passed = TRUE;
    if ( !Passed ){

        printf("\nTest 04, Lgm_McIlwain_L(): %15s    %15s    %15s    %15s\n", "       L       ", "       I       ", "       Bm       ", "       M       ");
        printf("                   Expected: %.15g   %.15g   %.15g   %.15g\n", L_expected,  I_expected, Bm_expected, M_expected );
        printf("                        Got: %.15g   %.15g   %.15g   %.15g\n\n\n", L,  I, Bm, M);
    }

    if ( (fp_got = fopen( "check_McIlwain_L_04.got", "w" )) != NULL ) {
        fprintf( fp_got, "%.15lf\n", L );
        fprintf( fp_got, "%.15lf\n", I );
        fprintf( fp_got, "%.15lf\n", Bm );
        fprintf( fp_got, "%.15lf\n", M );
        fclose( fp_got );
    } else {
        printf("Could not open file: check_McIlwain_L_04.got\\n");
    }

    fail_unless( Passed, "Lgm_McIlwain_L(): Regression test failed. Compare 'expected' and 'got' files: check_McIlwain_L_04.expected check_McIlwain_L_04.got\n" );


    return;
}
END_TEST

START_TEST(test_MCILWAIN_05) {

    int                 Passed = FALSE;
    long int            Date;
    double              L, I, Bm, M, UTC;
    double              L_expected=-9e99, I_expected=-9e99, Bm_expected=-9e99, M_expected=-9e99;
    Lgm_Vector          u, v;
    FILE                *fp_expected;
    FILE                *fp_got;

    if ( (fp_expected = fopen( "check_McIlwain_L_05.expected", "r" )) != NULL ) {

        fscanf( fp_expected, "%lf", &L_expected );
        fscanf( fp_expected, "%lf", &I_expected );
        fscanf( fp_expected, "%lf", &Bm_expected );
        fscanf( fp_expected, "%lf", &M_expected );
        fclose( fp_expected );
    } else {
        printf("Lgm_McIlwain_L(): Cant open file: check_McIlwain_L_05.expected\n" );
    }


    Date = 20090101; UTC  = 0.0;
    Lgm_Set_Coord_Transforms( Date, UTC, mInfo->c );
    u.x = -7.0; u.y = 0.0; u.z = 0.0;
    Lgm_Convert_Coords( &u, &v, SM_TO_GSM, mInfo->c );
    Lgm_MagModelInfo_Set_MagModel( LGM_CDIP, LGM_EXTMODEL_NULL, mInfo );
    L = Lgm_McIlwain_L( Date, UTC, &v, 90.0, 1, &I, &Bm, &M, mInfo );

    if (    (fabs( L-L_expected ) < 1e-7) && (fabs( I-I_expected ) < 1e-7) && (fabs( Bm-Bm_expected ) < 1e-7) && (fabs( M-M_expected ) < 1e-7) ) Passed = TRUE;
    if ( !Passed ){

        printf("\nTest 05, Lgm_McIlwain_L(): %15s    %15s    %15s    %15s\n", "       L       ", "       I       ", "       Bm       ", "       M       ");
        printf("                   Expected: %.15g   %.15g   %.15g   %.15g\n", L_expected,  I_expected, Bm_expected, M_expected );
        printf("                        Got: %.15g   %.15g   %.15g   %.15g\n\n\n", L,  I, Bm, M);
    }

    if ( (fp_got = fopen( "check_McIlwain_L_05.got", "w" )) != NULL ) {
        fprintf( fp_got, "%.15lf\n", L );
        fprintf( fp_got, "%.15lf\n", I );
        fprintf( fp_got, "%.15lf\n", Bm );
        fprintf( fp_got, "%.15lf\n", M );
        fclose( fp_got );
    } else {
        printf("Could not open file: check_McIlwain_L_05.got\\n");
    }

    fail_unless( Passed, "Lgm_McIlwain_L(): Regression test failed. Compare 'expected' and 'got' files: check_McIlwain_L_05.expected check_McIlwain_L_05.got\n" );


    return;
}
END_TEST

START_TEST(test_MCILWAIN_06) {

    int                 Passed = FALSE;
    long int            Date;
    double              L, I, Bm, M, UTC, Bmin, Blocal;
    double              L_expected=-9e99, I_expected=-9e99, Bm_expected=-9e99, M_expected=-9e99, Bmin_expected=-9e99, Blocal_expected=-9e99;
    Lgm_Vector          u, v;
    FILE                *fp_expected;
    FILE                *fp_got;

    if ( (fp_expected = fopen( "check_McIlwain_L_06.expected", "r" )) != NULL ) {

        fscanf( fp_expected, "%lf", &L_expected );
        fscanf( fp_expected, "%lf", &I_expected );
        fscanf( fp_expected, "%lf", &Bm_expected );
        fscanf( fp_expected, "%lf", &M_expected );
        fscanf( fp_expected, "%lf", &Blocal_expected );
        fscanf( fp_expected, "%lf", &Bmin_expected );
        fclose( fp_expected );
    } else {
        printf("Lgm_McIlwain_L(): Cant open file: check_McIlwain_L_06.expected\n" );
    }


    Date = 20090101; UTC  = 0.0;
    Lgm_Set_Coord_Transforms( Date, UTC, mInfo->c );
    u.x = -4.0; u.y = 0.0; u.z = 1.0;
    Lgm_Convert_Coords( &u, &v, GSM_TO_GSM, mInfo->c );
    Lgm_MagModelInfo_Set_MagModel( LGM_IGRF, LGM_EXTMODEL_T89, mInfo );
    L = Lgm_McIlwain_L( Date, UTC, &v, 90.0, 1, &I, &Bm, &M, mInfo );
    Blocal = mInfo->Blocal;
    Bmin = mInfo-> Bmin;

    if (    (fabs( L-L_expected ) < 1e-7) && (fabs( I-I_expected ) < 1e-7) && (fabs( Bm-Bm_expected ) < 1e-7) && (fabs( M-M_expected ) < 1e-7) && (fabs( Blocal-Blocal_expected ) < 1e-7) && (fabs( Bmin-Bmin_expected ) < 1e-7)) Passed = TRUE;
    if ( !Passed ){

        printf("\nTest 06, Lgm_McIlwain_L(): %15s    %15s    %15s    %15s\n", "       L       ", "       I       ", "       Bm       ", "       M       ");
        printf("                   Expected: %.15g   %.15g   %.15g   %.15g   %.15g   %.15g\n", L_expected,  I_expected, Bm_expected, M_expected, Blocal_expected, Bmin_expected);
        printf("                        Got: %.15g   %.15g   %.15g   %.15g   %.15g   %.15g\n\n\n", L,  I, Bm, M, Blocal, Bmin);
    }

    if ( (fp_got = fopen( "check_McIlwain_L_06.got", "w" )) != NULL ) {
        fprintf( fp_got, "%.15lf\n", L );
        fprintf( fp_got, "%.15lf\n", I );
        fprintf( fp_got, "%.15lf\n", Bm );
        fprintf( fp_got, "%.15lf\n", M );
        fprintf( fp_got, "%.15lf\n", Blocal );
        fprintf( fp_got, "%.15lf\n", Bmin );
        fclose( fp_got );
    } else {
        printf("Could not open file: check_McIlwain_L_06.got\\n");
    }

    fflush(stdout);
    fail_unless( Passed, "Lgm_McIlwain_L(): Regression test failed. Compare 'expected' and 'got' files: check_McIlwain_L_06.expected check_McIlwain_L_06.got\n" );


    return;
}
END_TEST

START_TEST(test_MCILWAIN_07) {

    int                 Passed = FALSE;
    long int            Date;
    double              L, I, Bm, M, UTC;
    double              L_expected=-9e99, I_expected=-9e99, Bm_expected=-9e99, M_expected=-9e99;
    Lgm_Vector          u, v;
    FILE                *fp_expected;
    FILE                *fp_got;

    if ( (fp_expected = fopen( "check_McIlwain_L_07.expected", "r" )) != NULL ) {

        fscanf( fp_expected, "%lf", &L_expected );
        fscanf( fp_expected, "%lf", &I_expected );
        fscanf( fp_expected, "%lf", &Bm_expected );
        fscanf( fp_expected, "%lf", &M_expected );
        fclose( fp_expected );
    } else {
        printf("Lgm_McIlwain_L(): Cant open file: check_McIlwain_L_07.expected\n" );
    }


    Date = 20090101; UTC  = 0.0;
    Lgm_Set_Coord_Transforms( Date, UTC, mInfo->c );
    u.x = -4.0; u.y = 0.0; u.z = 1.0;
    Lgm_Convert_Coords( &u, &v, SM_TO_GSM, mInfo->c );
    Lgm_MagModelInfo_Set_MagModel( LGM_IGRF, LGM_EXTMODEL_T89, mInfo );
    L = Lgm_McIlwain_L( Date, UTC, &v, 90.0, 1, &I, &Bm, &M, mInfo );

    if (    (fabs( L-L_expected ) < 1e-7) && (fabs( I-I_expected ) < 1e-7) && (fabs( Bm-Bm_expected ) < 1e-7) && (fabs( M-M_expected ) < 1e-7) ) Passed = TRUE;
    if ( !Passed ){

        printf("\nTest 07, Lgm_McIlwain_L(): %15s    %15s    %15s    %15s\n", "       L       ", "       I       ", "       Bm       ", "       M       ");
        printf("                   Expected: %.15g   %.15g   %.15g   %.15g\n", L_expected,  I_expected, Bm_expected, M_expected );
        printf("                        Got: %.15g   %.15g   %.15g   %.15g\n\n\n", L,  I, Bm, M);
    }

    if ( (fp_got = fopen( "check_McIlwain_L_07.got", "w" )) != NULL ) {
        fprintf( fp_got, "%.15lf\n", L );
        fprintf( fp_got, "%.15lf\n", I );
        fprintf( fp_got, "%.15lf\n", Bm );
        fprintf( fp_got, "%.15lf\n", M );
        fclose( fp_got );
    } else {
        printf("Could not open file: check_McIlwain_L_07.got\\n");
    }

    fflush(stdout);
    fail_unless( Passed, "Lgm_McIlwain_L(): Regression test failed. Compare 'expected' and 'got' files: check_McIlwain_L_07.expected check_McIlwain_L_07.got\n" );


    return;
}
END_TEST

Suite *McIlwain_L_suite(void) {

  Suite *s = suite_create("MCILWAIN_L_TESTS");

  TCase *tc_McIlwain_L = tcase_create("McIlwain L-Shell Values");
  tcase_add_checked_fixture(tc_McIlwain_L, McIlwain_L_Setup, McIlwain_L_TearDown);

  tcase_add_test(tc_McIlwain_L, test_MCILWAIN_01);
  tcase_add_test(tc_McIlwain_L, test_MCILWAIN_02);
  tcase_add_test(tc_McIlwain_L, test_MCILWAIN_03);
  tcase_add_test(tc_McIlwain_L, test_MCILWAIN_04);
  tcase_add_test(tc_McIlwain_L, test_MCILWAIN_05);
  tcase_add_test(tc_McIlwain_L, test_MCILWAIN_06);
  tcase_add_test(tc_McIlwain_L, test_MCILWAIN_07);

  suite_add_tcase(s, tc_McIlwain_L);

  return s;

}

int main(void) {

    int      number_failed;
    Suite   *s  = McIlwain_L_suite();
    SRunner *sr = srunner_create(s);

    printf("\n\n");
    srunner_run_all(sr, CK_NORMAL);
    number_failed = srunner_ntests_failed(sr);
    srunner_free(sr);

    return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;

}
