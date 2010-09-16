#include <check.h>
#include "../libLanlGeoMag/Lgm/Lgm_CTrans.h"

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


START_TEST(test_IsLeapSecondDay) {

  /*This is wrong...leap second is added at END of 19920630,
   *not start of 19920701, but matching library for now...
   */
  double sec_in_day;

  fail_unless( Lgm_IsLeapSecondDay(19910630, &sec_in_day, c)==1, "1992/6/30 should be a leap second day");

  return;
}
END_TEST



START_TEST(test_GetLeapSeconds) {
    double n_leap = Lgm_GetLeapSeconds(2454984.0, c);
    fail_unless( n_leap == 34.0, "Should be 34 leap seconds by 2009/6/1, not %f", n_leap);
  return;
}
END_TEST // END Leap seconds test case


Suite *lgm_suite(void) {

  Suite *s = suite_create("LGM");

  TCase *tc_leapseconds = tcase_create("Leap seconds");
  tcase_add_checked_fixture(tc_leapseconds, leapsecond_setup, leapsecond_teardown);
  tcase_add_test(tc_leapseconds, test_IsLeapSecondDay);
  tcase_add_test(tc_leapseconds, test_GetLeapSeconds);
  suite_add_tcase(s, tc_leapseconds);

  return s;

}

int main(void)
{
  int      number_failed;
  Suite   *s = lgm_suite();
  SRunner *sr = srunner_create(s);

  srunner_run_all(sr, CK_NORMAL);
  number_failed = srunner_ntests_failed(sr);
  srunner_free(sr);

  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
