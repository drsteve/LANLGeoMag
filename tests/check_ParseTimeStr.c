#include <check.h>
#include "../libLanlGeoMag/Lgm/Lgm_CTrans.h"

/*BEGIN Leap seconds test case*/
Lgm_CTrans *c;

void ParseTimeStr_setup(void) {
    c = Lgm_init_ctrans( 0 );
    Lgm_LoadLeapSeconds( c );  
    return;
}

void ParseTimeStr_teardown(void) {
    Lgm_free_ctrans( c ) ;
    return;
}


START_TEST(test_ISO_01) {

    int          Result;
    char         TimeString[128];
    Lgm_DateTime d;

    sprintf(TimeString, "2009-06-21T05:12:34.123456789Z");
    Result = ParseTimeString( TimeString, &d, c );

    fail_unless( (Result   == 1),    "Error returned by ParseTimeString()");
    fail_unless( (d.Year   == 2009), "For string %s Year should be 2009, got %d", TimeString, d.Year );
    fail_unless( (d.Month  == 6),    "For string %s Month should be 6, got %d", TimeString, d.Month );
    fail_unless( (d.Day    == 21),   "For string %s Day should be 21, got %d", TimeString, d.Day );
    fail_unless( (d.Hour   == 5),    "For string %s Hour should be 5, got %d", TimeString, d.Hour );
    fail_unless( (d.Minute == 12),   "For string %s Minute should be 12, got %d", TimeString, d.Minute );
    fail_unless( (d.Second == 34.123456789),   "For string %s Second should be 34.123456789, got %lf", TimeString, d.Second );

    return;
}
END_TEST

START_TEST(test_ISO_02) {

    int          Result;
    char         TimeString[128];
    Lgm_DateTime d;

    sprintf(TimeString, "2010-W52-6T00:00:00Z");
    Result = ParseTimeString( TimeString, &d, c );

    fail_unless( (Result   == 1),    "Error returned by ParseTimeString()");
    fail_unless( (d.Year   == 2011), "For string %s Year should be 2011, got %d", TimeString, d.Year );
    fail_unless( (d.Month  == 1),    "For string %s Month should be 1, got %d", TimeString, d.Month );
    fail_unless( (d.Day    == 1),    "For string %s Day should be 1, got %d", TimeString, d.Day );
    fail_unless( (d.Hour   == 0),    "For string %s Hour should be 0, got %d", TimeString, d.Hour );
    fail_unless( (d.Minute == 0),    "For string %s Minute should be 0, got %d", TimeString, d.Minute );
    fail_unless( (d.Second == 0.0),  "For string %s Second should be 0, got %d", TimeString, d.Second );

    return;
}
END_TEST


Suite *ParseTimeStr_suite(void) {

  Suite *s = suite_create("ParseTimeStr");

  TCase *tc_ParseTimeStr = tcase_create("ISO Time Parser");
  tcase_add_checked_fixture(tc_ParseTimeStr, ParseTimeStr_setup, ParseTimeStr_teardown);
  tcase_add_test(tc_ParseTimeStr, test_ISO_01);
  tcase_add_test(tc_ParseTimeStr, test_ISO_02);
  suite_add_tcase(s, tc_ParseTimeStr);

  return s;

}

int main(void) {

  int      number_failed;
  Suite   *s  = ParseTimeStr_suite();
  SRunner *sr = srunner_create(s);

  srunner_run_all(sr, CK_NORMAL);
  number_failed = srunner_ntests_failed(sr);
  srunner_free(sr);

  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;

}
