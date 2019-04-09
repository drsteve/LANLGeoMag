#include <check.h>
#include "../libLanlGeoMag/Lgm/Lgm_MagModelInfo.h"

/*
 *  Regression tests for Tsyganenko model calculations
 */


Lgm_MagModelInfo    *mInfo;

void Tsyganenko_Setup(void) {
    mInfo = Lgm_InitMagInfo();
    return;
}

void Tsyganenko_TearDown(void) {
    Lgm_FreeMagInfo( mInfo );
    return;
}


START_TEST(test_Tsyganenko_01) {
    Lgm_Vector        Bexpect, Pos, Btest, Udiff;
    int               nTests, nPass, nFail, transflag, Passed=FALSE;
    int               retVal;
    double            del, Kp;
    char              buff[262], extModel[10], intModel[10];
    char              IsoDate[80];
    long long         TT2000;
    FILE              *testfile, *outfile;
    Lgm_DateTime      d;

    int makeNew = 1;

    /* read test file */
    testfile = fopen("check_Tsyganenko.expected","r");
    if (makeNew) outfile = fopen("check_Tsyganenko.got", "w");

    /* step through test cases one line at a time */
    nTests = 0;
    nPass = 0;
    nFail = 0;
    Passed = TRUE;
    while( fgets(buff,260,testfile) != NULL) {
        if (buff[0]!='#') {
            // read line
            sscanf(buff,"%s %s %s %lf %lf %lf %lf %lf %lf %lf",
		   &IsoDate[0], &extModel[0], &intModel[0], &Kp, 
                   &Pos.x, &Pos.y, &Pos.z, &Bexpect.x, &Bexpect.y, &Bexpect.z);

	    Lgm_InitMagInfoDefaults( mInfo );
            IsoTimeStringToDateTime( IsoDate, &d, mInfo->c );
            Lgm_Set_Coord_Transforms( d.Date, d.Time, mInfo->c );
	    if (!strncmp(intModel, "IGRF", 10)) {
	        Lgm_Set_Lgm_B_IGRF_InternalModel( mInfo );
	    }
	    else if (!strncmp(intModel, "CDIP", 10)) {
	        Lgm_Set_Lgm_B_cdip_InternalModel( mInfo );
	    }
	    else if (!strncmp(intModel, "EDIP", 10)) {
	        Lgm_Set_Lgm_B_edip_InternalModel( mInfo );
	    }
	    else {
	        nFail++;
	        printf("Test %d bad internal model %s\n", nTests, intModel);
		continue;
	    }
	    if (!strncmp(extModel, "T89", 10)) {
	        Lgm_Set_Lgm_B_T89(mInfo);
	    }
	    else {
	        nFail++;
	        printf("Test %d bad external model %s\n", nTests, extModel);
		continue;
	    }
	    mInfo->Kp = Kp;
	    retVal = mInfo->Bfield(&Pos, &Btest, mInfo);
	    if (retVal != 1) {
	        nFail++;
	        printf("Test %d odd return from Lgm_B_T89\n", nTests);
	        continue;
		}
            //printf("IsoDate = %s; d.Date, d.time = %ld, %lf \n", IsoDate, d.Date, d.Time);
            Udiff.x = Btest.x - Bexpect.x;
            Udiff.y = Btest.y - Bexpect.y;
            Udiff.z = Btest.z - Bexpect.z;
            del = Lgm_Magnitude(&Udiff);
            nTests++;
            if (fabs(del) <= 1.0e-5) {
                nPass++;
                printf("Test %d passed\n", nTests);
                }
            else {
                nFail++;
                printf("*****  warning : difference >= 1.0e-5 nT  *****\n");
                printf("Test %d failed (diff: %g %g %g   %g)\n", nTests, Udiff.x, Udiff.y, Udiff.z, fabs(del));
                }
            if (makeNew) fprintf(outfile, "%s %s %s %lf %lf %lf %lf %lf %lf %lf\n", IsoDate, extModel, intModel, Kp, Pos.x, Pos.y, Pos.z, Btest.x, Btest.y, Btest.z);

            }
        else {
            if (makeNew) fprintf(outfile, "%s", buff);
            }
        }
    if (nFail>0) Passed = FALSE;
    fclose(testfile);
    if (makeNew) fclose(outfile);
    printf("Result: %d tests pass; %d tests fail (Precision=1.0e-5 nT)\n", nPass, nFail);
    fflush(stdout);

    ck_assert_msg( Passed, "Tsyganenko tests failed.\n" );

}
END_TEST


Suite *Tsyganenko_suite(void) {

  Suite *s = suite_create("TSYGANENKO_TESTS");

  TCase *tc_Tsyganenko = tcase_create("Tsyganenko Models");
  tcase_add_checked_fixture(tc_Tsyganenko, Tsyganenko_Setup, Tsyganenko_TearDown);

  tcase_add_test(tc_Tsyganenko, test_Tsyganenko_01);

  suite_add_tcase(s, tc_Tsyganenko);

  return s;

}

int main(void) {

    int      number_failed;
    Suite   *s  = Tsyganenko_suite();
    SRunner *sr = srunner_create(s);

    printf("\n\n");
    srunner_run_all(sr, CK_NORMAL);
    number_failed = srunner_ntests_failed(sr);
    srunner_free(sr);

    return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;

}
