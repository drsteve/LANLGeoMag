#include <check.h>
#include "../libLanlGeoMag/Lgm/Lgm_PolyRoots.h"
#include <stdio.h>
#include <stdlib.h>
#define TRUE    1
#define FALSE   0

/*
 *  Unit and Regression tests for PolyRoots solvers
 */


START_TEST(test_PolyRoots_01) {

    int             nReal, nReal_expected, Passed=FALSE;
    double          a, b, c;
    double          z1_r_expected, z1_i_expected;
    double          z2_r_expected, z2_i_expected;
    double complex  z1, z2, Residual1, Residual2;
    FILE            *fp_expected;
    FILE            *fp_got;

    if ( (fp_expected = fopen( "check_PolyRoots_01.expected", "r" )) != NULL ) {

        // real and imaginary parts of solution
        fscanf( fp_expected, "%d", &nReal_expected );
        fscanf( fp_expected, "%lf %lf", &z1_r_expected, &z1_i_expected );
        fscanf( fp_expected, "%lf %lf", &z2_r_expected, &z2_i_expected );
        fclose( fp_expected );
    } else {
        printf("Cant open file: check_PolyRoots_01.expected\n" );
    }

    a = 3.0; b = 6.0; c = 4.0;
    nReal = Lgm_QuadraticRoots( a, b, c, &z1, &z2 );
    Residual1 = a*z1*z1 + b*z1 + c;
    Residual2 = a*z2*z2 + b*z2 + c;

    if (   (fabs( creal(z1)-z1_r_expected) < 1e-10) && (fabs( cimag(z1)-z1_i_expected) < 1e-10)
        && (fabs( creal(z2)-z2_r_expected) < 1e-10) && (fabs( cimag(z2)-z2_i_expected) < 1e-10)
        && (cabs(Residual1)< 1e-10) && (cabs(Residual2) < 1e-10)
        && (nReal == nReal_expected) ) Passed = TRUE;

    if ( !Passed ){
        printf("\nTest 01,  Lgm_QuadraticRoots(): Expected %d real roots, got %d\n\n", nReal_expected, nReal );
        printf("       Root 1, Expected: %.15lf + %.15lf I\n",   z1_r_expected, z1_i_expected );
        printf("                    Got: %.15lf + %.15lf I\n", creal(z1),     cimag(z1) );
        printf("               Residual: %.15g + %.15g    |Residual| = %.15g\n\n", creal(Residual1), cimag(Residual1), cabs(Residual1) );
        printf("       Root 2, Expected: %.15lf + %.15lf I\n",   z2_r_expected, z2_i_expected );
        printf("                    Got: %.15lf + %.15lf I\n\n", creal(z2),     cimag(z2) );
        printf("               Residual: %.15g + %.15g    |Residual| = %.15g\n\n\n", creal(Residual2), cimag(Residual2), cabs(Residual2) );
    }

    if ( (fp_got = fopen( "check_PolyRoots_01.got", "w" )) != NULL ) {
        fprintf( fp_got, "%d\n", nReal );
        fprintf( fp_got, "%.15lf %.15lf\n", creal(z1), cimag(z1) );
        fprintf( fp_got, "%.15lf %.15lf\n", creal(z2), cimag(z2) );
        fclose( fp_got );
    } else {
        printf("Could not open file: check_PolyRoots_01.got\\n");
    }

    fail_unless( Passed, "Lgm_QuadraticRoots(): Regression test failed. Compare 'got' and 'expected' files: check_PolyRoots_01.got check_PolyRoots_01.expected\n" );


    return;
}
END_TEST

START_TEST(test_PolyRoots_02) {

    int             nReal, nReal_expected, Passed=FALSE;
    double          a, b, c;
    double          x1_expected;
    double          z2_r_expected, z2_i_expected;
    double          z3_r_expected, z3_i_expected;
    double          x1;
    double complex  z2, z3;
    double          Residual1;
    double complex  Residual2, Residual3;
    FILE            *fp_expected;
    FILE            *fp_got;

    if ( (fp_expected = fopen( "check_PolyRoots_02.expected", "r" )) != NULL ) {

        // real and imaginary parts of solution
        fscanf( fp_expected, "%d", &nReal_expected );
        fscanf( fp_expected, "%lf", &x1_expected );
        fscanf( fp_expected, "%lf %lf", &z2_r_expected, &z2_i_expected );
        fscanf( fp_expected, "%lf %lf", &z3_r_expected, &z3_i_expected );
        fclose( fp_expected );
    } else {
        printf("Cant open file: check_PolyRoots_02.expected\n" );
    }

    a = 7.0; b = -1.0; c = 36.0;
    nReal = Lgm_CubicRoots( a, b, c, &x1, &z2, &z3 );
    Residual1 = x1*x1*x1 + a*x1*x1 + b*x1 + c;
    Residual2 = z2*z2*z2 + a*z2*z2 + b*z2 + c;
    Residual3 = z3*z3*z3 + a*z3*z3 + b*z3 + c;

    if (   (fabs( creal(z2)-z2_r_expected) < 1e-10) && (fabs( cimag(z2)-z2_i_expected) < 1e-10)
        && (fabs( creal(z3)-z3_r_expected) < 1e-10) && (fabs( cimag(z3)-z3_i_expected) < 1e-10)
        && (fabs( x1-x1_expected) < 1e-10)
        && (Residual1 < 1e-10) && (cabs(Residual2)< 1e-10) && (cabs(Residual3) < 1e-10)
        && (nReal == nReal_expected) ) Passed = TRUE;

    if ( !Passed ){
        printf("\nTest 02,  Lgm_CubicRoots(): Expected %d real roots, got %d\n\n", nReal_expected, nReal );
        printf("       Root 1, Expected: %.15lf\n",   x1_expected );
        printf("                    Got: %.15lf\n", x1 );
        printf("               Residual: %.15g\n\n", Residual1 );
        printf("       Root 2, Expected: %.15lf + %.15lf I\n",   z2_r_expected, z2_i_expected );
        printf("                    Got: %.15lf + %.15lf I\n", creal(z2),     cimag(z2) );
        printf("               Residual: %.15g + %.15g    |Residual| = %.15g\n\n", creal(Residual2), cimag(Residual2), cabs(Residual2) );
        printf("       Root 3, Expected: %.15lf + %.15lf I\n",   z3_r_expected, z3_i_expected );
        printf("                    Got: %.15lf + %.15lf I\n", creal(z3),     cimag(z3) );
        printf("               Residual: %.15g + %.15g    |Residual| = %.15g\n\n", creal(Residual3), cimag(Residual3), cabs(Residual3) );
    }

    if ( (fp_got = fopen( "check_PolyRoots_02.got", "w" )) != NULL ) {
        fprintf( fp_got, "%d\n", nReal );
        fprintf( fp_got, "%.15lf\n", x1 );
        fprintf( fp_got, "%.15lf %.15lf\n", creal(z2), cimag(z2) );
        fprintf( fp_got, "%.15lf %.15lf\n", creal(z3), cimag(z3) );
        fclose( fp_got );
    } else {
        printf("Could not open file: check_PolyRoots_02.got\\n");
    }

    fail_unless( Passed, "Lgm_CubicRoots(): Regression test failed. Compare 'got' and 'expected' files: check_PolyRoots_02.got check_PolyRoots_02.expected\n" );


    return;
}
END_TEST

START_TEST(test_PolyRoots_03) {

    int             nReal, nReal_expected, Passed=FALSE;
    double          b, c, d, e;
    double          z1_r_expected, z1_i_expected;
    double          z2_r_expected, z2_i_expected;
    double          z3_r_expected, z3_i_expected;
    double          z4_r_expected, z4_i_expected;
    double complex  z1, z2, z3, z4, ztmp;
    double complex  Residual1, Residual2, Residual3, Residual4;
    FILE            *fp_expected;
    FILE            *fp_got;

    if ( (fp_expected = fopen( "check_PolyRoots_03.expected", "r" )) != NULL ) {

        // real and imaginary parts of solution
        fscanf( fp_expected, "%d", &nReal_expected );
        fscanf( fp_expected, "%lf %lf", &z1_r_expected, &z1_i_expected );
        fscanf( fp_expected, "%lf %lf", &z2_r_expected, &z2_i_expected );
        fscanf( fp_expected, "%lf %lf", &z3_r_expected, &z3_i_expected );
        fscanf( fp_expected, "%lf %lf", &z4_r_expected, &z4_i_expected );
        fclose( fp_expected );
    } else {
        printf("Cant open file: check_PolyRoots_03.expected\n" );
    }

    b = -37.0; c = -1.0; d = 7.0; e = 60.0;
    nReal = Lgm_QuarticRoots( b, c, d, e, &z1, &z2, &z3, &z4 );
    /*
     * The roots dshould be:
     * 37.020721858115365 0.000000000000000
     * 1.233413878980077 0.000000000000000
     * -0.627067868547723 -0.959579313119009
     * -0.627067868547723 0.959579313119009
     * On some machines the last two get swapped.
     * Do a swap on them if cimag(z4) < cimag(z3)
     */
    if ( cimag(z4) < cimag(z3) ) {
        ztmp = z4;
        z4 = z3;
        z3 = ztmp;
    }

    Residual1 = z1*z1*z1*z1 + b*z1*z1*z1 + c*z1*z1 + d*z1 + e;
    Residual2 = z2*z2*z2*z2 + b*z2*z2*z2 + c*z2*z2 + d*z2 + e;
    Residual3 = z3*z3*z3*z3 + b*z3*z3*z3 + c*z3*z3 + d*z3 + e;
    Residual4 = z4*z4*z4*z4 + b*z4*z4*z4 + c*z4*z4 + d*z4 + e;

    if (   (fabs( creal(z1)-z1_r_expected) < 1e-10) && (fabs( cimag(z1)-z1_i_expected) < 1e-10)
        && (fabs( creal(z2)-z2_r_expected) < 1e-10) && (fabs( cimag(z2)-z2_i_expected) < 1e-10)
        && (fabs( creal(z3)-z3_r_expected) < 1e-10) && (fabs( cimag(z3)-z3_i_expected) < 1e-10)
        && (fabs( creal(z4)-z4_r_expected) < 1e-10) && (fabs( cimag(z4)-z4_i_expected) < 1e-10)
        && (cabs(Residual1) < 1e-8) && (cabs(Residual2)< 1e-8) && (cabs(Residual3) < 1e-8) && (cabs(Residual4) < 1e-8)
        && (nReal == nReal_expected) ) Passed = TRUE;

    if ( !Passed ){
        printf("\nTest 03,  Lgm_QuarticRoots(): Expected %d real roots, got %d\n\n", nReal_expected, nReal );
        printf("       Root 1, Expected: %.15lf + %.15lf I\n",   z1_r_expected, z1_i_expected );
        printf("                    Got: %.15lf + %.15lf I\n", creal(z1),     cimag(z1) );
        printf("               Residual: %.15g + %.15g    |Residual| = %.15g\n\n", creal(Residual1), cimag(Residual1), cabs(Residual1) );
        printf("       Root 2, Expected: %.15lf + %.15lf I\n",   z2_r_expected, z2_i_expected );
        printf("                    Got: %.15lf + %.15lf I\n", creal(z2),     cimag(z2) );
        printf("               Residual: %.15g + %.15g    |Residual| = %.15g\n\n", creal(Residual2), cimag(Residual2), cabs(Residual2) );
        printf("       Root 3, Expected: %.15lf + %.15lf I\n",   z3_r_expected, z3_i_expected );
        printf("                    Got: %.15lf + %.15lf I\n", creal(z3),     cimag(z3) );
        printf("               Residual: %.15g + %.15g    |Residual| = %.15g\n\n", creal(Residual3), cimag(Residual3), cabs(Residual3) );
        printf("       Root 4, Expected: %.15lf + %.15lf I\n",   z4_r_expected, z4_i_expected );
        printf("                    Got: %.15lf + %.15lf I\n", creal(z4),     cimag(z4) );
        printf("               Residual: %.15g + %.15g    |Residual| = %.15g\n\n", creal(Residual4), cimag(Residual4), cabs(Residual4) );
    }

    if ( (fp_got = fopen( "check_PolyRoots_03.got", "w" )) != NULL ) {
        fprintf( fp_got, "%d\n", nReal );
        fprintf( fp_got, "%.15lf %.15lf\n", creal(z1), cimag(z1) );
        fprintf( fp_got, "%.15lf %.15lf\n", creal(z2), cimag(z2) );
        fprintf( fp_got, "%.15lf %.15lf\n", creal(z3), cimag(z3) );
        fprintf( fp_got, "%.15lf %.15lf\n", creal(z4), cimag(z4) );
        fclose( fp_got );
    } else {
        printf("Could not open file: check_PolyRoots_03.got\\n");
    }

    fail_unless( Passed, "Lgm_QuarticRoots(): Regression test failed. Compare 'got' and 'expected' files: check_PolyRoots_03.got check_PolyRoots_03.expected\n" );


    return;
}
END_TEST

START_TEST(test_PolyRoots_04) {

    int             nReal, nReal_expected, Passed=FALSE;
    double          a, b, c;
    double          x1_expected;
    double          x1;
    double          Residual1;
    FILE            *fp_expected;
    FILE            *fp_got;

    if ( (fp_expected = fopen( "check_PolyRoots_04.expected", "r" )) != NULL ) {

        // real and imaginary parts of solution
        fscanf( fp_expected, "%d", &nReal_expected );
        fscanf( fp_expected, "%lf", &x1_expected );
        fclose( fp_expected );
    } else {
        printf("Cant open file: check_PolyRoots_04.expected\n" );
    }

    a = 7.0; b = -1.0; c = 36.0;
    x1 = Lgm_CubicRealRoot( a, b, c );
    Residual1 = x1*x1*x1 + a*x1*x1 + b*x1 + c;

    if (   (fabs( x1-x1_expected) < 1e-10) && (Residual1 < 1e-10) ) Passed = TRUE;

    if ( !Passed ){
        printf("\nTest 04,  Lgm_CubicRealRoot(): Expected 1 real root\n\n" );
        printf("       Root 1, Expected: %.15lf\n",   x1_expected );
        printf("                    Got: %.15lf\n", x1 );
        printf("               Residual: %.15g\n\n\n", Residual1 );
    }

    if ( (fp_got = fopen( "check_PolyRoots_04.got", "w" )) != NULL ) {
        fprintf( fp_got, "%d\n", nReal );
        fprintf( fp_got, "%.15lf\n", x1 );
        fclose( fp_got );
    } else {
        printf("Could not open file: check_PolyRoots_04.got\\n");
    }

    fail_unless( Passed, "Lgm_CubicRealRoot(): Regression test failed. Compare 'got' and 'expected' files: check_PolyRoots_04.got check_PolyRoots_04.expected\n" );


    return;
}
END_TEST

START_TEST(test_PolyRoots_05) {

    int             nErrors, nReal, nPolys=0, Passed=FALSE;
    double          a, b, c;
    double complex  z1, z2, f1, f2;

    nErrors = 0;
    for (a = -10.0; a< 10.0; a += 1.0 ){
        for (b = -10.0; b< 10.0; b += 1.0 ){
            for (c = -10.0; c< 10.0; c += 1.0 ){
                nReal = Lgm_QuadraticRoots( a, b, c, &z1, &z2 );
                f1 = a*z1*z1 + b*z1 + c;
                f2 = a*z2*z2 + b*z2 + c;
                if ( ( cabs(f1) > 1e-10 ) || ( cabs(f2) > 1e-10 ) ) {
                    ++nErrors;
                    printf("a, b, c = %g %g %g   z1 = %lf + %lf I   f1 = %g + %g I    |f1| = %g\n",   a, b, c, creal(z1), cimag(z1), creal(f1), cimag(f1), cabs(f1) );
                    printf("a, b, c = %g %g %g   z2 = %lf + %lf I   f1 = %g + %g I    |f1| = %g\n\n", a, b, c, creal(z2), cimag(z2), creal(f2), cimag(f2), cabs(f2) );
                }
                ++nPolys;
            }
        }
    }

    if ( nErrors == 0 ) Passed = TRUE;

    if ( !Passed ){
        printf("\nTest 05,  Lgm_QuadraticRoots(): %d Polynomials solved, %d errors\n\n\n", nPolys, nErrors );
    }

    fail_unless( Passed, "Lgm_QuadraticRoots(): Unit test failed. Residuals out of acceptable limits.\n" );

    return;
}
END_TEST

START_TEST(test_PolyRoots_06) {

    int             nErrors, nReal, nPolys=0, Passed=FALSE;
    double          b, c, d, x1, f1;
    double complex  z2, z3, f2, f3;

    nErrors = 0;
    for (b = -10.0; b< 10.0; b += 1.0 ){
        for (c = -10.0; c< 10.0; c += 1.0 ){
            for (d = -10.0; d< 10.0; d += 1.0 ){
                nReal = Lgm_CubicRoots( b, c, d, &x1, &z2, &z3 );
                f1 = x1*x1*x1 + b*x1*x1 + c*x1 + d;
                f2 = z2*z2*z2 + b*z2*z2 + c*z2 + d;
                f3 = z3*z3*z3 + b*z3*z3 + c*z3 + d;
                if (( fabs(f1) > 1e-10 ) ||  ( cabs(f1) > 1e-10 ) || ( cabs(f2) > 1e-10 ) ) {
                    ++nErrors;
                    printf("a, b, c = %g %g %g   x1 = %lf    |f1| = %g\n",   b, c, d, x1, fabs(f1) );
                    printf("a, b, c = %g %g %g   z2 = %lf + %lf I   f1 = %g + %g I    |f1| = %g\n\n", b, c, d, creal(z2), cimag(z2), creal(f2), cimag(f2), cabs(f2) );
                    printf("a, b, c = %g %g %g   z3 = %lf + %lf I   f3 = %g + %g I    |f1| = %g\n",   b, c, d, creal(z3), cimag(z3), creal(f3), cimag(f3), cabs(f3) );
                }
                ++nPolys;
            }
        }
    }

    if ( nErrors == 0 ) Passed = TRUE;

    if ( !Passed ){
        printf("\nTest 06,  Lgm_CubicRoots(): %d Polynomials solved, %d errors\n\n\n", nPolys, nErrors );
    }

    fail_unless( Passed, "Lgm_CubicRoots(): Unit test failed. Residuals out of acceptable limits.\n" );

    return;
}
END_TEST

START_TEST(test_PolyRoots_07) {

    int             nErrors, nReal, nPolys=0, Passed=FALSE;
    double          b, c, d, e;
    double complex  z1, z2, z3, z4, f1, f2, f3, f4;

    nErrors = 0;
    for (b = -10.0; b< 10.0; b += 1.0 ){
        for (c = -10.0; c< 10.0; c += 1.0 ){
            for (d = -10.0; d< 10.0; d += 1.0 ){
                for (e = -10.0; e< 10.0; e += 1.0 ){
                    nReal = Lgm_QuarticRoots( b, c, d, e, &z1, &z2, &z3, &z4 );
                    f1 = z1*z1*z1*z1 + b*z1*z1*z1 + c*z1*z1 + d*z1 + e;
                    f2 = z2*z2*z2*z2 + b*z2*z2*z2 + c*z2*z2 + d*z2 + e;
                    f3 = z3*z3*z3*z3 + b*z3*z3*z3 + c*z3*z3 + d*z3 + e;
                    f4 = z4*z4*z4*z4 + b*z4*z4*z4 + c*z4*z4 + d*z4 + e;
                    if (( cabs(f1) > 1e-6 ) ||  ( cabs(f2) > 1e-6 ) || ( cabs(f3) > 1e-6 ) || ( cabs(f4) > 1e-6 ) ) {
                        ++nErrors;
                        printf("b, c, d, e = %g %g %g %g   z1 = %lf + %lf I   f1 = %g + %g I    |f1| = %g\n\n", b, c, d, e, creal(z1), cimag(z1), creal(f1), cimag(f1), cabs(f1) );
                        printf("b, c, d, e = %g %g %g %g   z2 = %lf + %lf I   f2 = %g + %g I    |f2| = %g\n",   b, c, d, e, creal(z2), cimag(z2), creal(f2), cimag(f2), cabs(f2) );
                        printf("b, c, d, e = %g %g %g %g   z3 = %lf + %lf I   f3 = %g + %g I    |f3| = %g\n\n", b, c, d, e, creal(z3), cimag(z3), creal(f3), cimag(f3), cabs(f3) );
                        printf("b, c, d, e = %g %g %g %g   z4 = %lf + %lf I   f4 = %g + %g I    |f4| = %g\n",   b, c, d, e, creal(z4), cimag(z4), creal(f4), cimag(f4), cabs(f4) );
                    }
                    ++nPolys;
                }
            }
        }
    }

    if ( nErrors == 0 ) Passed = TRUE;

    if ( !Passed ){
        printf("\nTest 07,  Lgm_QuarticRoots(): %d Polynomials solved, %d errors\n\n\n", nPolys, nErrors );
    }

    fail_unless( Passed, "Lgm_QuarticRoots(): Unit test failed. Residuals out of acceptable limits.\n" );

    return;
}
END_TEST


Suite *PolyRoots_suite(void) {

  Suite *s = suite_create("POLY_ROOTS_TESTS");

  TCase *tc_PolyRoots = tcase_create("PolyRoots");
  //tcase_add_checked_fixture(tc_PolyRoots, PolyRoots_Setup, PolyRoots_TearDown);

  tcase_add_test(tc_PolyRoots, test_PolyRoots_01);
  tcase_add_test(tc_PolyRoots, test_PolyRoots_02);
  tcase_add_test(tc_PolyRoots, test_PolyRoots_03);
  tcase_add_test(tc_PolyRoots, test_PolyRoots_04);
  tcase_add_test(tc_PolyRoots, test_PolyRoots_05);
  tcase_add_test(tc_PolyRoots, test_PolyRoots_06);
  tcase_add_test(tc_PolyRoots, test_PolyRoots_07);

  suite_add_tcase(s, tc_PolyRoots);

  return s;

}

int main(void) {

    int      number_failed;
    Suite   *s  = PolyRoots_suite();
    SRunner *sr = srunner_create(s);

    printf("\n\n");
    srunner_run_all(sr, CK_NORMAL);
    number_failed = srunner_ntests_failed(sr);
    srunner_free(sr);

    return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;

}
