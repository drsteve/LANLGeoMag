// REALLY NEED TO CHECK ALL THE PARAMS
#include <check.h>
#include "../libLanlGeoMag/Lgm/Lgm_Misc.h"

/*BEGIN String case test case*/

START_TEST(test_StrToLower) {

    char*        Result;
    char         TestString[128];

    printf("Starting test_StrToUpper\n");
    sprintf(TestString, "AN_UPPER_CASE_STRING");
    printf("Converting string to  lower: %s\n", TestString);
    Result = Lgm_StrToLower(TestString, 128);
    printf("Result: %s\n\n", TestString);
    fflush(stdout);
    fail_unless(!strcmp(TestString, "an_upper_case_string"));
    return;
}
END_TEST

START_TEST(test_StrToUpper) {

    char*        Result;
    char         TestString[128];

    printf("Starting test_StrToUpper\n");
    sprintf(TestString, "a_lower_case_string");
    printf("Converting string to upper: %s\n", TestString);
    Result = Lgm_StrToUpper(TestString, 128);
    printf("Result: %s\n\n", TestString);
    fflush(stdout);
    fail_unless(!strcmp(TestString, "A_LOWER_CASE_STRING"));
    return;
}
END_TEST

START_TEST(test_ReplaceSubString) {

    int          Result;
    char         TestString[128], NewString[128], Substr[5];

    printf("Starting test_ReplaceSubString\n");
    sprintf(Substr, "%02d", 4);
    sprintf(TestString, "2024-%MM-01T12:00:00Z");
    printf("Replacing wildcard substring: %s\n", TestString);
    Lgm_ReplaceSubString(NewString, TestString, "%MM", Substr);
    printf("Result: %s\n\n", NewString);
    fflush(stdout);
    fail_unless(!strcmp(NewString, "2024-04-01T12:00:00Z"));
    return;
}
END_TEST

START_TEST(test_ReplaceSubString2) {

    char         TestString[128], Substr[5];
    char*        NewString;

    printf("Starting test_ReplaceSubString2\n");
    sprintf(Substr, "%02d", 4);
    sprintf(TestString, "2024-%MM-01T12:00:00Z");
    printf("Replacing wildcard substring: %s\n", TestString);
    Lgm_ReplaceSubString2(&NewString, TestString, "%MM", Substr);
    printf("Result: %s\n\n", NewString);
    fflush(stdout);
    fail_unless(!strcmp(NewString, "2024-04-01T12:00:00Z"));
    free(NewString);
    return;
}
END_TEST


Suite *ManipStr_suite(void) {

  Suite *s = suite_create("STRING_FUNCTION_TESTS");

  TCase *tc_manipStr = tcase_create("ISO Time Parser");
  tcase_add_test(tc_manipStr, test_StrToLower);
  tcase_add_test(tc_manipStr, test_StrToUpper);
  tcase_add_test(tc_manipStr, test_ReplaceSubString);
  tcase_add_test(tc_manipStr, test_ReplaceSubString2);
  suite_add_tcase(s, tc_manipStr);

  return s;

}

int main(void) {

    int      number_failed;
    Suite   *s  = ManipStr_suite();
    SRunner *sr = srunner_create(s);

    printf("\n\n");
    srunner_run_all(sr, CK_ENV);
    number_failed = srunner_ntests_failed(sr);
    srunner_free(sr);

    return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;

}
