// REALLY NEED TO CHECK ALL THE PARAMS
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
    double       t, fYear, JD, T;
    Lgm_DateTime d;

    sprintf(TimeString, "2009-06-21T05:12:34.123456789Z");
    printf("Converting ISO string: %s\n", TimeString);
    Result = IsoTimeStringToDateTime( TimeString, &d, c );

    fail_unless( (Result   == 1),        "Error returned by IsoTimeStringToDateTime()");
    fail_unless( (d.Date   == 20090621), "For string %s Year should be 20090621, got %ld", TimeString, d.Date );
    fail_unless( (d.Year   == 2009),     "For string %s Year should be 2009, got %d", TimeString, d.Year );
    fail_unless( (d.Month  == 6),        "For string %s Month should be 6, got %d", TimeString, d.Month );
    fail_unless( (d.Day    == 21),       "For string %s Day should be 21, got %d", TimeString, d.Day );
    fail_unless( (d.Doy    == 172),      "For string %s Doy should be 172, got %d", TimeString, d.Doy );
    t = 5.0 + 12.0/60.0 + 34.123456789/3600.0;
    fail_unless( (fabs(d.Time-t)<1e-12), "For string %s Time should be %.10lf, got %.10lf", TimeString, t, d.Time );
    fail_unless( (d.Hour   == 5),        "For string %s Hour should be 5, got %d", TimeString, d.Hour );
    fail_unless( (d.Minute == 12),       "For string %s Minute should be 12, got %d", TimeString, d.Minute );
    fail_unless( (fabs(d.Second-34.123456789)<1e-12),   "For string %s Second should be 34.123456789, got %.10lf", TimeString, d.Second );
    fail_unless( (d.Week   == 25),       "For string %s Week should be 25, got %d", TimeString, d.Week );
    fail_unless( (d.wYear  == 2009),     "For string %s Iso Week Year should be 2009, got %d", TimeString, d.wYear );
    fail_unless( (d.Dow    == 7),        "For string %s Day of Week should be 7, got %d", TimeString, d.Dow );
    fail_unless( !strcmp(d.DowStr,"Sun"),"For string %s DowStr should be \"Sun\", got \"%s\"", TimeString, d.DowStr );
    fYear = 2009.0 + (172.0 + t/24.0)/365.0;
    fail_unless( (fabs(d.fYear-fYear)<1e-12), "For string %s fYear should be %.10lf, got %.10lf", TimeString, fYear, d.fYear );
    JD = 2455003.5+t/24.0;
    fail_unless( (fabs(d.JD-JD)<1e-12),       "For string %s JD should be %.10lf, got %.10lf", TimeString, JD, d.JD );
    T = (JD-2451545.0)/36525.0;
    fail_unless( (fabs(d.T-T)<1e-12),         "For string %s T should be %.10lf, got %.10lf", TimeString, T, d.T );
    fail_unless( (fabs(d.DaySeconds-86400.0)<1e-10), "For string %s DaySeconds should be 86400, got %g", TimeString, d.DaySeconds );
    fail_unless( (d.TZD_sgn == 1),       "For string %s TZD_sgn should be 1, got %d", TimeString, d.TZD_sgn );
    fail_unless( (d.TZD_hh  == 0),       "For string %s TZD_hh should be 0, got %d", TimeString, d.TZD_hh );
    fail_unless( (d.TZD_mm  == 0),       "For string %s TZD_mm should be 0, got %d", TimeString, d.TZD_mm );
    fail_unless( (d.TimeSystem == LGM_TIME_SYS_UTC), "For string %s TimeSystem should be %d, got %d", TimeString, LGM_TIME_SYS_UTC, d.TimeSystem );

    return;
}
END_TEST


START_TEST(test_ISO_02) {

    int          Result;
    char         TimeString[128];
    double       t, fYear, JD, T;
    Lgm_DateTime d;

    sprintf(TimeString, "2009-06-21T05:12:34.123456789-03:30");
    printf("Converting ISO string: %s\n", TimeString);
    Result = IsoTimeStringToDateTime( TimeString, &d, c );

    fail_unless( (Result   == 1),        "Error returned by IsoTimeStringToDateTime()");
    fail_unless( (d.Date   == 20090621), "For string %s Year should be 20090621, got %ld", TimeString, d.Date );
    fail_unless( (d.Year   == 2009),     "For string %s Year should be 2009, got %d", TimeString, d.Year );
    fail_unless( (d.Month  == 6),        "For string %s Month should be 6, got %d", TimeString, d.Month );
    fail_unless( (d.Day    == 21),       "For string %s Day should be 21, got %d", TimeString, d.Day );
    fail_unless( (d.Doy    == 172),      "For string %s Doy should be 172, got %d", TimeString, d.Doy );
    t = 8.0 + 42.0/60.0 + 34.123456789/3600.0;
    fail_unless( (fabs(d.Time-t)<1e-12), "For string %s Time should be %.10lf, got %.10lf", TimeString, t, d.Time );
    fail_unless( (d.Hour   == 8),        "For string %s Hour should be 8, got %d", TimeString, d.Hour );
    fail_unless( (d.Minute == 42),       "For string %s Minute should be 42, got %d", TimeString, d.Minute );
    fail_unless( (fabs(d.Second-34.123456789)<1e-12),   "For string %s Second should be 34.123456789, got %lf", TimeString, d.Second );
    fail_unless( (d.Week   == 25),       "For string %s Week should be 25, got %d", TimeString, d.Week );
    fail_unless( (d.wYear  == 2009),     "For string %s Iso Week Year should be 2009, got %d", TimeString, d.wYear );
    fail_unless( (d.Dow    == 7),        "For string %s Day of Week should be 7, got %d", TimeString, d.Dow );
    fail_unless( !strcmp(d.DowStr,"Sun"),"For string %s DowStr should be \"Sun\", got \"%s\"", TimeString, d.DowStr );
    fYear = 2009.0 + (172.0 + t/24.0)/365.0;
    fail_unless( (fabs(d.fYear-fYear)<1e-12), "For string %s fYear should be %.10lf, got %.10lf", TimeString, fYear, d.fYear );
    JD = 2455003.5+t/24.0;
    fail_unless( (fabs(d.JD-JD)<1e-12),       "For string %s JD should be %.10lf, got %.10lf", TimeString, JD, d.JD );
    T = (JD-2451545.0)/36525.0;
    fail_unless( (fabs(d.T-T)<1e-12),         "For string %s T should be %.10lf, got %.10lf", TimeString, T, d.T );
    fail_unless( (fabs(d.DaySeconds-86400.0)<1e-10), "For string %s DaySeconds should be 86400, got %g", TimeString, d.DaySeconds );
    fail_unless( (d.TZD_sgn == -1),       "For string %s TZD_sgn should be -1, got %d", TimeString, d.TZD_sgn );
    fail_unless( (d.TZD_hh  == 3),       "For string %s TZD_hh should be 3, got %d", TimeString, d.TZD_hh );
    fail_unless( (d.TZD_mm  == 30),       "For string %s TZD_mm should be 30, got %d", TimeString, d.TZD_mm );
    fail_unless( (d.TimeSystem == LGM_TIME_SYS_UTC), "For string %s TimeSystem should be %d, got %d", TimeString, LGM_TIME_SYS_UTC, d.TimeSystem );

    return;
}
END_TEST


START_TEST(test_ISO_03) {

    int          Result;
    char         TimeString[128];
    double       t, fYear, JD, T;
    Lgm_DateTime d;

    // test TZ offset.
    sprintf(TimeString, "2009-06-21T05:12:34.123456789+03:30");
    printf("Converting ISO string: %s\n", TimeString);
    Result = IsoTimeStringToDateTime( TimeString, &d, c );

    fail_unless( (Result   == 1),        "Error returned by IsoTimeStringToDateTime()");
    fail_unless( (d.Date   == 20090621), "For string %s Year should be 20090621, got %ld", TimeString, d.Date );
    fail_unless( (d.Year   == 2009),     "For string %s Year should be 2009, got %d", TimeString, d.Year );
    fail_unless( (d.Month  == 6),        "For string %s Month should be 6, got %d", TimeString, d.Month );
    fail_unless( (d.Day    == 21),       "For string %s Day should be 21, got %d", TimeString, d.Day );
    fail_unless( (d.Doy    == 172),      "For string %s Doy should be 172, got %d", TimeString, d.Doy );
    t = 1.0 + 42.0/60.0 + 34.123456789/3600.0;
    fail_unless( (fabs(d.Time-t)<1e-12), "For string %s Time should be %.10lf, got %.10lf", TimeString, t, d.Time );
    fail_unless( (d.Hour   == 1),        "For string %s Hour should be 1, got %d", TimeString, d.Hour );
    fail_unless( (d.Minute == 42),       "For string %s Minute should be 42, got %d", TimeString, d.Minute );
    fail_unless( (fabs(d.Second-34.123456789)<1e-12),   "For string %s Second should be 34.123456789, got %lf", TimeString, d.Second );
    fail_unless( (d.Week   == 25),       "For string %s Week should be 25, got %d", TimeString, d.Week );
    fail_unless( (d.wYear  == 2009),     "For string %s Iso Week Year should be 2009, got %d", TimeString, d.wYear );
    fail_unless( (d.Dow    == 7),        "For string %s Day of Week should be 7, got %d", TimeString, d.Dow );
    fail_unless( !strcmp(d.DowStr,"Sun"),"For string %s DowStr should be \"Sun\", got \"%s\"", TimeString, d.DowStr );
    fYear = 2009.0 + (172.0 + t/24.0)/365.0;
    fail_unless( (fabs(d.fYear-fYear)<1e-12), "For string %s fYear should be %.10lf, got %.10lf", TimeString, fYear, d.fYear );
    JD = 2455003.5+t/24.0;
    fail_unless( (fabs(d.JD-JD)<1e-12),       "For string %s JD should be %.10lf, got %.10lf", TimeString, JD, d.JD );
    T = (JD-2451545.0)/36525.0;
    fail_unless( (fabs(d.T-T)<1e-12),         "For string %s T should be %.10lf, got %.10lf", TimeString, T, d.T );
    fail_unless( (fabs(d.DaySeconds-86400.0)<1e-10), "For string %s DaySeconds should be 86400, got %g", TimeString, d.DaySeconds );
    fail_unless( (d.TZD_sgn == 1),       "For string %s TZD_sgn should be 1, got %d", TimeString, d.TZD_sgn );
    fail_unless( (d.TZD_hh  == 3),       "For string %s TZD_hh should be 3, got %d", TimeString, d.TZD_hh );
    fail_unless( (d.TZD_mm  == 30),       "For string %s TZD_mm should be 30, got %d", TimeString, d.TZD_mm );
    fail_unless( (d.TimeSystem == LGM_TIME_SYS_UTC), "For string %s TimeSystem should be %d, got %d", TimeString, LGM_TIME_SYS_UTC, d.TimeSystem );

    return;
}
END_TEST


START_TEST(test_ISO_04) {

    int          Result;
    char         TimeString[128];
    double       t, fYear, JD, T;
    Lgm_DateTime d;

    // test roll-over to previous day due to TZ info.
    sprintf(TimeString, "2009-06-21T01:12:34.123456789+03:30");
    printf("Converting ISO string: %s\n", TimeString);
    Result = IsoTimeStringToDateTime( TimeString, &d, c );

    fail_unless( (Result   == 1),        "Error returned by IsoTimeStringToDateTime()");
    fail_unless( (d.Date   == 20090620), "For string %s Year should be 20090620, got %ld", TimeString, d.Date );
    fail_unless( (d.Year   == 2009),     "For string %s Year should be 2009, got %d", TimeString, d.Year );
    fail_unless( (d.Month  == 6),        "For string %s Month should be 6, got %d", TimeString, d.Month );
    fail_unless( (d.Day    == 20),       "For string %s Day should be 20, got %d", TimeString, d.Day );
    fail_unless( (d.Doy    == 171),      "For string %s Doy should be 171, got %d", TimeString, d.Doy );
    t = 21.0 + 42.0/60.0 + 34.123456789/3600.0;
    fail_unless( (fabs(d.Time-t)<1e-12), "For string %s Time should be %.10lf, got %.10lf", TimeString, t, d.Time );
    fail_unless( (d.Hour   == 21),        "For string %s Hour should be 21, got %d", TimeString, d.Hour );
    fail_unless( (d.Minute == 42),       "For string %s Minute should be 42, got %d", TimeString, d.Minute );
    fail_unless( (fabs(d.Second-34.123456789)<1e-12),   "For string %s Second should be 34.123456789, got %lf", TimeString, d.Second );
    fail_unless( (d.Week   == 25),       "For string %s Week should be 25, got %d", TimeString, d.Week );
    fail_unless( (d.wYear  == 2009),     "For string %s Iso Week Year should be 2009, got %d", TimeString, d.wYear );
    fail_unless( (d.Dow    == 6),        "For string %s Day of Week should be 6, got %d", TimeString, d.Dow );
    fail_unless( !strcmp(d.DowStr,"Sat"),"For string %s DowStr should be \"Sat\", got \"%s\"", TimeString, d.DowStr );
    fYear = 2009.0 + (171.0 + t/24.0)/365.0;
    fail_unless( (fabs(d.fYear-fYear)<1e-12), "For string %s fYear should be %.10lf, got %.10lf", TimeString, fYear, d.fYear );
    JD = 2455002.5+t/24.0;
    fail_unless( (fabs(d.JD-JD)<1e-12),       "For string %s JD should be %.10lf, got %.10lf", TimeString, JD, d.JD );
    T = (JD-2451545.0)/36525.0;
    fail_unless( (fabs(d.T-T)<1e-12),         "For string %s T should be %.10lf, got %.10lf", TimeString, T, d.T );
    fail_unless( (fabs(d.DaySeconds-86400.0)<1e-10), "For string %s DaySeconds should be 86400, got %g", TimeString, d.DaySeconds );
    fail_unless( (d.TZD_sgn == 1),       "For string %s TZD_sgn should be 1, got %d", TimeString, d.TZD_sgn );
    fail_unless( (d.TZD_hh  == 3),       "For string %s TZD_hh should be 3, got %d", TimeString, d.TZD_hh );
    fail_unless( (d.TZD_mm  == 30),       "For string %s TZD_mm should be 30, got %d", TimeString, d.TZD_mm );
    fail_unless( (d.TimeSystem == LGM_TIME_SYS_UTC), "For string %s TimeSystem should be %d, got %d", TimeString, LGM_TIME_SYS_UTC, d.TimeSystem );

    return;
}
END_TEST


START_TEST(test_ISO_05) {

    int          Result;
    char         TimeString[128];
    double       t, fYear, JD, T;
    Lgm_DateTime d;

    /*
     * More complicatehd test. The date rolls back over a year boundary.  I.e.
     * from Jan 1, 2009 local time to Dec 31, 2008 UTC. And that was a leap
     * seconds day so DaySeconds should be 86401.0. The Week Year doesnt roll
     * back -- it should still be 2009. Plus, 2008 was a leap year, so we test
     * lots of things here....
     */
    sprintf(TimeString, "2009-01-01T01:52:23+07");
    printf("Converting ISO string: %s\n", TimeString);
    Result = IsoTimeStringToDateTime( TimeString, &d, c );

    fail_unless( (Result   == 1),        "Error returned by IsoTimeStringToDateTime()");
    fail_unless( (d.Date   == 20081231), "For string %s Year should be 20081231, got %ld", TimeString, d.Date );
    fail_unless( (d.Year   == 2008),     "For string %s Year should be 2008, got %d", TimeString, d.Year );
    fail_unless( (d.Month  == 12),        "For string %s Month should be 12, got %d", TimeString, d.Month );
    fail_unless( (d.Day    == 31),       "For string %s Day should be 31, got %d", TimeString, d.Day );
    fail_unless( (d.Doy    == 366),      "For string %s Doy should be 366, got %d", TimeString, d.Doy );
    t = 18.0 + 52.0/60.0 + 23.0/3600.0;
    fail_unless( (fabs(d.Time-t)<1e-12), "For string %s Time should be %.10lf, got %.10lf", TimeString, t, d.Time );
    fail_unless( (d.Hour   == 18),        "For string %s Hour should be 18, got %d", TimeString, d.Hour );
    fail_unless( (d.Minute == 52),       "For string %s Minute should be 52, got %d", TimeString, d.Minute );
    fail_unless( (fabs(d.Second-23.0)<1e-12),   "For string %s Second should be 23, got %lf", TimeString, d.Second );
    fail_unless( (d.Week   ==  1),       "For string %s Week should be 1, got %d", TimeString, d.Week );
    fail_unless( (d.wYear  == 2009),     "For string %s Iso Week Year should be 2009, got %d", TimeString, d.wYear );
    fail_unless( (d.Dow    == 3),        "For string %s Day of Week should be 3, got %d", TimeString, d.Dow );
    fail_unless( !strcmp(d.DowStr,"Wed"),"For string %s DowStr should be \"Wed\", got \"%s\"", TimeString, d.DowStr );
    fYear = 2008.0 + (366.0 + t/24.0)/366.0;
    fail_unless( (fabs(d.fYear-fYear)<1e-12), "For string %s fYear should be %.10lf, got %.10lf", TimeString, fYear, d.fYear );
    JD = 2454831.5+t/24.0;
    fail_unless( (fabs(d.JD-JD)<1e-12),       "For string %s JD should be %.10lf, got %.10lf", TimeString, JD, d.JD );
    T = (JD-2451545.0)/36525.0;
    fail_unless( (fabs(d.T-T)<1e-12),         "For string %s T should be %.10lf, got %.10lf", TimeString, T, d.T );
    fail_unless( (fabs(d.DaySeconds-86401.0)<1e-10), "For string %s DaySeconds should be 86401, got %g", TimeString, d.DaySeconds );
    fail_unless( (d.TZD_sgn == 1),       "For string %s TZD_sgn should be 1, got %d", TimeString, d.TZD_sgn );
    fail_unless( (d.TZD_hh  == 7),       "For string %s TZD_hh should be 3, got %d", TimeString, d.TZD_hh );
    fail_unless( (d.TZD_mm  == 0),       "For string %s TZD_mm should be 30, got %d", TimeString, d.TZD_mm );
    fail_unless( (d.TimeSystem == LGM_TIME_SYS_UTC), "For string %s TimeSystem should be %d, got %d", TimeString, LGM_TIME_SYS_UTC, d.TimeSystem );

    return;
}
END_TEST


START_TEST(test_ISO_06) {

    int          Result;
    char         TimeString[128];
    double       t, fYear, JD, T;
    Lgm_DateTime d;

    // Like a previous test but truncated form
    sprintf(TimeString, "20090621T051234.123456789-0330");
    printf("Converting ISO string: %s\n", TimeString);
    Result = IsoTimeStringToDateTime( TimeString, &d, c );

    fail_unless( (Result   == 1),        "Error returned by IsoTimeStringToDateTime()");
    fail_unless( (d.Date   == 20090621), "For string %s Year should be 20090621, got %ld", TimeString, d.Date );
    fail_unless( (d.Year   == 2009),     "For string %s Year should be 2009, got %d", TimeString, d.Year );
    fail_unless( (d.Month  == 6),        "For string %s Month should be 6, got %d", TimeString, d.Month );
    fail_unless( (d.Day    == 21),       "For string %s Day should be 21, got %d", TimeString, d.Day );
    fail_unless( (d.Doy    == 172),      "For string %s Doy should be 172, got %d", TimeString, d.Doy );
    t = 8.0 + 42.0/60.0 + 34.123456789/3600.0;
    fail_unless( (fabs(d.Time-t)<1e-12), "For string %s Time should be %.10lf, got %.10lf", TimeString, t, d.Time );
    fail_unless( (d.Hour   == 8),        "For string %s Hour should be 8, got %d", TimeString, d.Hour );
    fail_unless( (d.Minute == 42),       "For string %s Minute should be 42, got %d", TimeString, d.Minute );
    fail_unless( (fabs(d.Second-34.123456789)<1e-12),   "For string %s Second should be 34.123456789, got %lf", TimeString, d.Second );
    fail_unless( (d.Week   == 25),       "For string %s Week should be 25, got %d", TimeString, d.Week );
    fail_unless( (d.wYear  == 2009),     "For string %s Iso Week Year should be 2009, got %d", TimeString, d.wYear );
    fail_unless( (d.Dow    == 7),        "For string %s Day of Week should be 7, got %d", TimeString, d.Dow );
    fail_unless( !strcmp(d.DowStr,"Sun"),"For string %s DowStr should be \"Sun\", got \"%s\"", TimeString, d.DowStr );
    fYear = 2009.0 + (172.0 + t/24.0)/365.0;
    fail_unless( (fabs(d.fYear-fYear)<1e-12), "For string %s fYear should be %.10lf, got %.10lf", TimeString, fYear, d.fYear );
    JD = 2455003.5+t/24.0;
    fail_unless( (fabs(d.JD-JD)<1e-12),       "For string %s JD should be %.10lf, got %.10lf", TimeString, JD, d.JD );
    T = (JD-2451545.0)/36525.0;
    fail_unless( (fabs(d.T-T)<1e-12),         "For string %s T should be %.10lf, got %.10lf", TimeString, T, d.T );
    fail_unless( (fabs(d.DaySeconds-86400.0)<1e-10), "For string %s DaySeconds should be 86400, got %g", TimeString, d.DaySeconds );
    fail_unless( (d.TZD_sgn == -1),       "For string %s TZD_sgn should be -1, got %d", TimeString, d.TZD_sgn );
    fail_unless( (d.TZD_hh  == 3),       "For string %s TZD_hh should be 3, got %d", TimeString, d.TZD_hh );
    fail_unless( (d.TZD_mm  == 30),       "For string %s TZD_mm should be 30, got %d", TimeString, d.TZD_mm );
    fail_unless( (d.TimeSystem == LGM_TIME_SYS_UTC), "For string %s TimeSystem should be %d, got %d", TimeString, LGM_TIME_SYS_UTC, d.TimeSystem );

    return;
}
END_TEST


START_TEST(test_ISO_07) {

    int          Result;
    char         TimeString[128];
    double       t, fYear, JD, T;
    Lgm_DateTime d;

    // Like a previous test but truncated form and no T
    sprintf(TimeString, "20090621 051234.123456789-0330");
    printf("Converting ISO string: %s\n", TimeString);
    Result = IsoTimeStringToDateTime( TimeString, &d, c );

    fail_unless( (Result   == 1),        "Error returned by IsoTimeStringToDateTime()");
    fail_unless( (d.Date   == 20090621), "For string %s Year should be 20090621, got %ld", TimeString, d.Date );
    fail_unless( (d.Year   == 2009),     "For string %s Year should be 2009, got %d", TimeString, d.Year );
    fail_unless( (d.Month  == 6),        "For string %s Month should be 6, got %d", TimeString, d.Month );
    fail_unless( (d.Day    == 21),       "For string %s Day should be 21, got %d", TimeString, d.Day );
    fail_unless( (d.Doy    == 172),      "For string %s Doy should be 172, got %d", TimeString, d.Doy );
    t = 8.0 + 42.0/60.0 + 34.123456789/3600.0;
    fail_unless( (fabs(d.Time-t)<1e-12), "For string %s Time should be %.10lf, got %.10lf", TimeString, t, d.Time );
    fail_unless( (d.Hour   == 8),        "For string %s Hour should be 8, got %d", TimeString, d.Hour );
    fail_unless( (d.Minute == 42),       "For string %s Minute should be 42, got %d", TimeString, d.Minute );
    fail_unless( (fabs(d.Second-34.123456789)<1e-12),   "For string %s Second should be 34.123456789, got %lf", TimeString, d.Second );
    fail_unless( (d.Week   == 25),       "For string %s Week should be 25, got %d", TimeString, d.Week );
    fail_unless( (d.wYear  == 2009),     "For string %s Iso Week Year should be 2009, got %d", TimeString, d.wYear );
    fail_unless( (d.Dow    == 7),        "For string %s Day of Week should be 7, got %d", TimeString, d.Dow );
    fail_unless( !strcmp(d.DowStr,"Sun"),"For string %s DowStr should be \"Sun\", got \"%s\"", TimeString, d.DowStr );
    fYear = 2009.0 + (172.0 + t/24.0)/365.0;
    fail_unless( (fabs(d.fYear-fYear)<1e-12), "For string %s fYear should be %.10lf, got %.10lf", TimeString, fYear, d.fYear );
    JD = 2455003.5+t/24.0;
    fail_unless( (fabs(d.JD-JD)<1e-12),       "For string %s JD should be %.10lf, got %.10lf", TimeString, JD, d.JD );
    T = (JD-2451545.0)/36525.0;
    fail_unless( (fabs(d.T-T)<1e-12),         "For string %s T should be %.10lf, got %.10lf", TimeString, T, d.T );
    fail_unless( (fabs(d.DaySeconds-86400.0)<1e-10), "For string %s DaySeconds should be 86400, got %g", TimeString, d.DaySeconds );
    fail_unless( (d.TZD_sgn == -1),       "For string %s TZD_sgn should be -1, got %d", TimeString, d.TZD_sgn );
    fail_unless( (d.TZD_hh  == 3),       "For string %s TZD_hh should be 3, got %d", TimeString, d.TZD_hh );
    fail_unless( (d.TZD_mm  == 30),       "For string %s TZD_mm should be 30, got %d", TimeString, d.TZD_mm );
    fail_unless( (d.TimeSystem == LGM_TIME_SYS_UTC), "For string %s TimeSystem should be %d, got %d", TimeString, LGM_TIME_SYS_UTC, d.TimeSystem );

    return;
}
END_TEST


START_TEST(test_ISO_022) {

    int          Result;
    char         TimeString[128];
    Lgm_DateTime d;

    sprintf(TimeString, "2010-W52-6T00:00:00Z");
    Result = IsoTimeStringToDateTime( TimeString, &d, c );

    fail_unless( (Result   == 1),    "Error returned by IsoTimeStringToDateTime()");
    fail_unless( (d.Year   == 2011), "For string %s Year should be 2011, got %d", TimeString, d.Year );
    fail_unless( (d.Month  == 1),    "For string %s Month should be 1, got %d", TimeString, d.Month );
    fail_unless( (d.Day    == 1),    "For string %s Day should be 1, got %d", TimeString, d.Day );
    fail_unless( (d.Hour   == 0),    "For string %s Hour should be 0, got %d", TimeString, d.Hour );
    fail_unless( (d.Minute == 0),    "For string %s Minute should be 0, got %d", TimeString, d.Minute );
    fail_unless( (d.Second == 0.0),  "For string %s Second should be 0, got %d", TimeString, d.Second );

    return;
}
END_TEST


START_TEST(test_ISO_033) {

    int          Result;
    char         TimeString[128];
    Lgm_DateTime d;

    sprintf(TimeString, "2009-W53-6");
    Result = IsoTimeStringToDateTime( TimeString, &d, c );

    fail_unless( (Result   == 1),    "Error returned by IsoTimeStringToDateTime()");
    fail_unless( (d.Year   == 2010), "For string %s Year should be 2010, got %d", TimeString, d.Year );
    fail_unless( (d.Month  == 1),    "For string %s Month should be 1, got %d", TimeString, d.Month );
    fail_unless( (d.Day    == 2),    "For string %s Day should be 2, got %d", TimeString, d.Day );
    fail_unless( (d.Hour   == 0),    "For string %s Hour should be 0, got %d", TimeString, d.Hour );
    fail_unless( (d.Minute == 0),    "For string %s Minute should be 0, got %d", TimeString, d.Minute );
    fail_unless( (d.Second == 0.0),  "For string %s Second should be 0, got %d", TimeString, d.Second );

    return;
}
END_TEST


START_TEST(test_ISO_040) {

    int          Result;
    char         TimeString[128];
    Lgm_DateTime d;

    // Invalid date test
    sprintf(TimeString, "2009-02-29");
    Result = IsoTimeStringToDateTime( TimeString, &d, c );
    fail_unless( (Result   == -1),    "%s is an invalid date and IsoTimeStringToDateTime() should have returned an error", TimeString );

    return;
}
END_TEST

START_TEST(test_ISO_055) {

    int          Result;
    char         TimeString[128];
    Lgm_DateTime d;

    // Invalid date test
    sprintf(TimeString, "2010-W53-2");
    Result = IsoTimeStringToDateTime( TimeString, &d, c );
    fail_unless( (Result   == -1),    "%s is an invalid date and IsoTimeStringToDateTime() should have returned an error", TimeString );

    return;
}
END_TEST

START_TEST(test_ISO_066) {

    int          Result;
    char         TimeString[128];
    Lgm_DateTime d;

    sprintf(TimeString, "2009W536");
    Result = IsoTimeStringToDateTime( TimeString, &d, c );

    fail_unless( (Result   == 1),    "Error returned by IsoTimeStringToDateTime()");
    fail_unless( (d.Year   == 2010), "For string %s Year should be 2010, got %d", TimeString, d.Year );
    fail_unless( (d.Month  == 1),    "For string %s Month should be 1, got %d", TimeString, d.Month );
    fail_unless( (d.Day    == 2),    "For string %s Day should be 2, got %d", TimeString, d.Day );
    fail_unless( (d.Hour   == 0),    "For string %s Hour should be 0, got %d", TimeString, d.Hour );
    fail_unless( (d.Minute == 0),    "For string %s Minute should be 0, got %d", TimeString, d.Minute );
    fail_unless( (d.Second == 0.0),  "For string %s Second should be 0, got %d", TimeString, d.Second );

    return;
}
END_TEST


Suite *ParseTimeStr_suite(void) {

  Suite *s = suite_create("ISO_TIME_STRING_PARSE_TESTS");

  TCase *tc_ParseTimeStr = tcase_create("ISO Time Parser");
  tcase_add_checked_fixture(tc_ParseTimeStr, ParseTimeStr_setup, ParseTimeStr_teardown);
  tcase_add_test(tc_ParseTimeStr, test_ISO_01);
  tcase_add_test(tc_ParseTimeStr, test_ISO_02);
  tcase_add_test(tc_ParseTimeStr, test_ISO_03);
  tcase_add_test(tc_ParseTimeStr, test_ISO_04);
  tcase_add_test(tc_ParseTimeStr, test_ISO_05);
  tcase_add_test(tc_ParseTimeStr, test_ISO_06);
  tcase_add_test(tc_ParseTimeStr, test_ISO_07);
  suite_add_tcase(s, tc_ParseTimeStr);

  return s;

}

int main(void) {

    int      number_failed;
    Suite   *s  = ParseTimeStr_suite();
    SRunner *sr = srunner_create(s);

    printf("\n\n");
    srunner_run_all(sr, CK_NORMAL);
    number_failed = srunner_ntests_failed(sr);
    srunner_free(sr);

    return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;

}
