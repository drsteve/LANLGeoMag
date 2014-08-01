#include <check.h>
#include "../libLanlGeoMag/Lgm/Lgm_MagModelInfo.h"

/*
 *  Validation and Regression tests for Lstar
 */


Lgm_MagModelInfo    *mInfo;
Lgm_CTrans          *c;

void Lstar_Setup(void) {
    c = Lgm_init_ctrans( 0 );
    mInfo = Lgm_InitMagInfo();
    return;
}

void Lstar_TearDown(void) {
    Lgm_FreeMagInfo( mInfo );
    Lgm_free_ctrans( c ) ;
    return;
}


START_TEST(test_LSTAR_01) {

    int                 Passed = FALSE;
    long int            Date;
    double              L, I, Bm, M, UTC;
    double              L_expected=-9e99, I_expected=-9e99, Bm_expected=-9e99, M_expected=-9e99;
    Lgm_Vector          u, v;
    FILE                *fp_expected;
    FILE                *fp_got;

    if ( (fp_expected = fopen( "check_Lstar_01.expected", "r" )) != NULL ) {

        fscanf( fp_expected, "%lf", &L_expected );
        fscanf( fp_expected, "%lf", &I_expected );
        fscanf( fp_expected, "%lf", &Bm_expected );
        fscanf( fp_expected, "%lf", &M_expected );
    } else {
        printf("Lgm_Lstar(): Cant open file: check_Lstar_01.expected\n" );
    }


    Date = 20090101; UTC  = 0.0; mInfo->Kp = 2; Lgm_Set_Coord_Transforms( Date, UTC, c );
    u.x = -4.0; u.y = 0.0; u.z = 1.0;
    L = Lgm_Lstar( Date, UTC, &u, 90.0, 1, &I, &Bm, &M, mInfo );

    if (    (fabs( L-L_expected ) < 1e-7) && (fabs( I-I_expected ) < 1e-7) && (fabs( Bm-Bm_expected ) < 1e-7) && (fabs( M-M_expected ) < 1e-7) ) Passed = TRUE;
    if ( !Passed ){
        printf("\nTest 01, Lgm_Lstar(): %15s    %15s    %15s    %15s\n", "       L       ", "       L       ", "       Bm       ", "       M       ");
        printf("                   Expected: %.15g   %.15g   %.15g   %.15g\n", L_expected,  I_expected, Bm_expected, M_expected );
        printf("                        Got: %.15g   %.15g   %.15g   %.15g\n\n\n", L,  I, Bm, M);
    }

    if ( (fp_got = fopen( "check_Lstar_01.got", "w" )) != NULL ) {
        fprintf( fp_got, "%.15lf\n", L );
        fprintf( fp_got, "%.15lf\n", I );
        fprintf( fp_got, "%.15lf\n", Bm );
        fprintf( fp_got, "%.15lf\n", M );
        fclose( fp_got );
    } else {
        printf("Could not open file: check_Lstar_01.got\\n");
    }

    fclose( fp_expected );
    fail_unless( Passed, "Lgm_Lstar(): Regression test failed. Compare 'expected' and 'got' files: check_Lstar_01.expected check_Lstar_01.got\n" );


    return;
}
END_TEST

START_TEST(test_LSTAR_02) {

    int                 Passed = FALSE;
    long int            Date;
    double              L, I, Bm, M, UTC;
    double              L_expected=-9e99, I_expected=-9e99, Bm_expected=-9e99, M_expected=-9e99;
    Lgm_Vector          u, v;
    FILE                *fp_expected;
    FILE                *fp_got;

    if ( (fp_expected = fopen( "check_Lstar_02.expected", "r" )) != NULL ) {

        fscanf( fp_expected, "%lf", &L_expected );
        fscanf( fp_expected, "%lf", &I_expected );
        fscanf( fp_expected, "%lf", &Bm_expected );
        fscanf( fp_expected, "%lf", &M_expected );
    } else {
        printf("Lgm_Lstar(): Cant open file: check_Lstar_02.expected\n" );
    }

    Date = 20090101; UTC  = 0.0; mInfo->Kp = 2; Lgm_Set_Coord_Transforms( Date, UTC, c );
    u.x = -6.6; u.y = 2.3; u.z = 0.4; 
    L = Lgm_Lstar( Date, UTC, &u, 90.0, 1, &I, &Bm, &M, mInfo );

    if (    (fabs( L-L_expected ) < 1e-7) && (fabs( I-I_expected ) < 1e-7) && (fabs( Bm-Bm_expected ) < 1e-7) && (fabs( M-M_expected ) < 1e-7) ) Passed = TRUE;
    if ( !Passed ){
        printf("\nTest 02, Lgm_Lstar(): %15s    %15s    %15s    %15s\n", "       L       ", "       L       ", "       Bm       ", "       M       ");
        printf("                   Expected: %.15g   %.15g   %.15g   %.15g\n", L_expected,  I_expected, Bm_expected, M_expected );
        printf("                        Got: %.15g   %.15g   %.15g   %.15g\n\n\n", L,  I, Bm, M);
    }

    if ( (fp_got = fopen( "check_Lstar_02.got", "w" )) != NULL ) {
        fprintf( fp_got, "%.15lf\n", L );
        fprintf( fp_got, "%.15lf\n", I );
        fprintf( fp_got, "%.15lf\n", Bm );
        fprintf( fp_got, "%.15lf\n", M );
        fclose( fp_got );
    } else {
        printf("Could not open file: check_Lstar_02.got\\n");
    }

    fclose( fp_expected );
    fail_unless( Passed, "Lgm_Lstar(): Regression test failed. Compare 'expected' and 'got' files: check_Lstar_02.expected check_Lstar_02.got\n" );


    return;
}
END_TEST


START_TEST(test_LSTAR_03) {

    int                 Passed = FALSE;
    long int            Date;
    double              L, I, Bm, M, UTC;
    double              L_expected=-9e99, I_expected=-9e99, Bm_expected=-9e99, M_expected=-9e99;
    Lgm_Vector          u, v;
    FILE                *fp_expected;
    FILE                *fp_got;

    if ( (fp_expected = fopen( "check_Lstar_03.expected", "r" )) != NULL ) {

        fscanf( fp_expected, "%lf", &L_expected );
        fscanf( fp_expected, "%lf", &I_expected );
        fscanf( fp_expected, "%lf", &Bm_expected );
        fscanf( fp_expected, "%lf", &M_expected );
    } else {
        printf("Lgm_Lstar(): Cant open file: check_Lstar_03.expected\n" );
    }

    Date = 20060823; UTC  = 13.213; mInfo->Kp = 5; Lgm_Set_Coord_Transforms( Date, UTC, c );
    u.x = 2.4; u.y = 1.2; u.z = 0.4; 
    L = Lgm_Lstar( Date, UTC, &u, 90.0, 1, &I, &Bm, &M, mInfo );

    if (    (fabs( L-L_expected ) < 1e-7) && (fabs( I-I_expected ) < 1e-7) && (fabs( Bm-Bm_expected ) < 1e-7) && (fabs( M-M_expected ) < 1e-7) ) Passed = TRUE;
    if ( !Passed ){
        printf("\nTest 02, Lgm_Lstar(): %15s    %15s    %15s    %15s\n", "       L       ", "       L       ", "       Bm       ", "       M       ");
        printf("                   Expected: %.15g   %.15g   %.15g   %.15g\n", L_expected,  I_expected, Bm_expected, M_expected );
        printf("                        Got: %.15g   %.15g   %.15g   %.15g\n\n\n", L,  I, Bm, M);
    }

    if ( (fp_got = fopen( "check_Lstar_03.got", "w" )) != NULL ) {
        fprintf( fp_got, "%.15lf\n", L );
        fprintf( fp_got, "%.15lf\n", I );
        fprintf( fp_got, "%.15lf\n", Bm );
        fprintf( fp_got, "%.15lf\n", M );
        fclose( fp_got );
    } else {
        printf("Could not open file: check_Lstar_03.got\\n");
    }

    fclose( fp_expected );
    fail_unless( Passed, "Lgm_Lstar(): Regression test failed. Compare 'expected' and 'got' files: check_Lstar_03.expected check_Lstar_03.got\n" );


    return;
}
END_TEST



Suite *Lstar_suite(void) {

  Suite *s = suite_create("LSTAR_L_TESTS");

  TCase *tc_Lstar = tcase_create("McIlwain L-Shell Values");
  tcase_add_checked_fixture(tc_Lstar, Lstar_Setup, Lstar_TearDown);

  tcase_add_test(tc_Lstar, test_LSTAR_01);
  tcase_add_test(tc_Lstar, test_LSTAR_02);
  tcase_add_test(tc_Lstar, test_LSTAR_03);

  suite_add_tcase(s, tc_Lstar);

  return s;

}

int main(void) {

    int      number_failed;
    Suite   *s  = Lstar_suite();
    SRunner *sr = srunner_create(s);

    printf("\n\n");
    srunner_run_all(sr, CK_NORMAL);
    number_failed = srunner_ntests_failed(sr);
    srunner_free(sr);

    return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;

}
