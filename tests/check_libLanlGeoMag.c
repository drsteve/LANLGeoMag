#include <check.h>
#include "../libLanlGeoMag/Lgm/Lgm_CTrans.h"

#define TRUE    1
#define FALSE   0

/*BEGIN Leap seconds test case*/
Lgm_CTrans *c;

void leapsecond_setup(void) {
    c = Lgm_init_ctrans( 0 );
    Lgm_LoadLeapSeconds( c );  
    return;
}

void leapsecond_teardown(void) {
    Lgm_free_ctrans( c ) ;
    return;
}


START_TEST(test_IsLeapSecondDay_01) {

    int    Result;
    double sec_in_day;

    printf("Testing to see if 19920630 is a leap second date (it should be)\n");
    Result = Lgm_IsLeapSecondDay( 19920630, &sec_in_day, c );
    ck_assert_msg( ((Result == 1) && (sec_in_day == 86401)), "1992/6/30 should be a leap second day, with sec_in_day = 86401 (got %g)", sec_in_day);

    return;
}
END_TEST

START_TEST(test_IsLeapSecondDay_02) {

    int    Result;
    double sec_in_day;

    printf("Testing to see if 19920629 is a leap second date (it should not be)\n");
    Result = Lgm_IsLeapSecondDay( 19920629, &sec_in_day, c );
    ck_assert_msg( ((Result == 0) && (sec_in_day == 86400)), "1992/6/29 should be a leap second day, with sec_in_day = 86400 (got %g)", sec_in_day);

    return;
}
END_TEST



START_TEST(test_GetLeapSeconds) {
    double n_leap = Lgm_GetLeapSeconds(2454984.0, c);
    printf("Check to see that number of leap seconds by 2009/6/1 is 34\n");
    ck_assert_msg( n_leap == 34.0, "Should be 34 leap seconds by 2009/6/1, not %f", n_leap);
  return;
}
END_TEST // END Leap seconds test case



START_TEST(test_VectorMagn) {
    Lgm_Vector  Utest;
    int         Passed;
    double      ans=3.7416573867739413, res;

    Utest.x = 1; Utest.y = 2; Utest.z = 3;
    printf("Check that Lgm_Magnitude gives known answer: ");

    res = Lgm_Magnitude(&Utest);
    Passed = (fabs(res-ans) <= 1e-12) ? TRUE : FALSE;
    if (Passed) {
        printf("Passed.\n");
    } else {
        printf("Failed.\n");
    }
    ck_assert(Passed);
  return;
}
END_TEST


START_TEST(test_VectorNorm) {
    Lgm_Vector  Utest;
    int         Passed;
    double      ans=1.0, res;

    Utest.x = 1; Utest.y = 2; Utest.z = 3;
    printf("Check that Lgm_NormalizeVector result has magnitude of 1: ");

    Lgm_NormalizeVector(&Utest);
    res = Lgm_Magnitude(&Utest);
    Passed = (fabs(res-ans) <= 1e-12) ? TRUE : FALSE;
    if (Passed) {
        printf("Passed.\n");
    } else {
        printf("Failed.\n");
    }
    ck_assert(Passed);
  return;
}
END_TEST




Suite *make_master_suite (void) {
    Suite *s = suite_create("LGM Master Suite");

    return s;
}

Suite *lgm_ls_suite(void) {

    Suite *s = suite_create("LEAP_SECOND_TESTS");
  
    TCase *tc_leapseconds = tcase_create("Leap seconds");
    tcase_add_checked_fixture(tc_leapseconds, leapsecond_setup, leapsecond_teardown);
    tcase_add_test(tc_leapseconds, test_IsLeapSecondDay_01);
    tcase_add_test(tc_leapseconds, test_IsLeapSecondDay_02);
    tcase_add_test(tc_leapseconds, test_GetLeapSeconds);
    suite_add_tcase(s, tc_leapseconds);
  
    return s;

}

Suite *lgm_vec_suite(void) {

    Suite *s = suite_create("VECTOR_OPERATION_TESTS");

    TCase *tc_vector = tcase_create("Test Lgm_Vector operations");
    tcase_add_test(tc_vector, test_VectorMagn);
    tcase_add_test(tc_vector, test_VectorNorm);
    suite_add_tcase(s, tc_vector);

    return s;

}


int main(void) {

    int      number_failed;
    Suite   *s = make_master_suite();
    SRunner *sr = srunner_create(s);
    srunner_add_suite (sr, lgm_ls_suite());
    srunner_add_suite (sr, lgm_vec_suite());

    printf("\n\n\n======================================================\n");
    printf("\n               UNIT/REGRESSION TESTS                   \n");
    printf("    Note: Some reported errors may be normal. For\n");
    printf("    example, we may be testing to see if the code \n");
    printf("    properly catches errors (e.g. like invalid \n");
    printf("    dates etc.)\n\n");
    printf("\n======================================================\n\n");
    srunner_run_all(sr, CK_ENV);
    number_failed = srunner_ntests_failed(sr);
    srunner_free(sr);

    return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;

}
