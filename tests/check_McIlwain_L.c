#include <check.h>
#include "../libLanlGeoMag/Lgm/Lgm_MagModelInfo.h"

/*
 *  Regression tests for McIlwain L-shell calculations
 */


Lgm_MagModelInfo    *mInfo;
Lgm_CTrans          *c;

void McIlwain_L_Setup(void) {
    c = Lgm_init_ctrans( 0 );
    mInfo = Lgm_InitMagInfo();
    return;
}

void McIlwain_L_TearDown(void) {
    Lgm_FreeMagInfo( mInfo );
    Lgm_free_ctrans( c ) ;
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
    } else {
        printf("Lgm_McIlwain_L(): Cant open file: check_McIlwain_L_01.expected\n" );
    }


    Date = 20090101; UTC  = 0.0; mInfo->Kp = 2; Lgm_Set_Coord_Transforms( Date, UTC, c );
    u.x = -4.0; u.y = 0.0; u.z = 1.0;
    L = Lgm_McIlwain_L( Date, UTC, &u, 90.0, 1, &I, &Bm, &M, mInfo );

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

    fclose( fp_expected );
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
    } else {
        printf("Lgm_McIlwain_L(): Cant open file: check_McIlwain_L_02.expected\n" );
    }

    Date = 20090101; UTC  = 0.0; mInfo->Kp = 2; Lgm_Set_Coord_Transforms( Date, UTC, c );
    u.x = -6.6; u.y = 2.3; u.z = 0.4;
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

    fclose( fp_expected );
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
    } else {
        printf("Lgm_McIlwain_L(): Cant open file: check_McIlwain_L_03.expected\n" );
    }

    Date = 20060823; UTC  = 13.213; mInfo->Kp = 5; Lgm_Set_Coord_Transforms( Date, UTC, c );
    u.x = 2.4; u.y = 1.2; u.z = 0.4;
    L = Lgm_McIlwain_L( Date, UTC, &u, 90.0, 1, &I, &Bm, &M, mInfo );

    if (    (fabs( L-L_expected ) < 1e-7) && (fabs( I-I_expected ) < 1e-7) && (fabs( Bm-Bm_expected ) < 1e-7) && (fabs( M-M_expected ) < 1e-7) ) Passed = TRUE;
    if ( !Passed ){
        printf("\nTest 02, Lgm_McIlwain_L(): %15s    %15s    %15s    %15s\n", "       L       ", "       I       ", "       Bm       ", "       M       ");
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

    fclose( fp_expected );
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
    } else {
        printf("Lgm_McIlwain_L(): Cant open file: check_McIlwain_L_04.expected\n" );
    }


    Date = 20090101; UTC  = 0.0; mInfo->Kp = 2; Lgm_Set_Coord_Transforms( Date, UTC, c );
    u.x = 7.0; u.y = 1.0; u.z = 3.0;
    L = Lgm_McIlwain_L( Date, UTC, &u, 90.0, 1, &I, &Bm, &M, mInfo );

    if (    (fabs( L-L_expected ) < 1e-7) && (fabs( I-I_expected ) < 1e-7) && (fabs( Bm-Bm_expected ) < 1e-7) && (fabs( M-M_expected ) < 1e-7) ) Passed = TRUE;
    if ( !Passed ){
        printf("\nTest 01, Lgm_McIlwain_L(): %15s    %15s    %15s    %15s\n", "       L       ", "       I       ", "       Bm       ", "       M       ");
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

    fclose( fp_expected );
    fail_unless( Passed, "Lgm_McIlwain_L(): Regression test failed. Compare 'expected' and 'got' files: check_McIlwain_L_04.expected check_McIlwain_L_04.got\n" );


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
