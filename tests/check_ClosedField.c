#include <check.h>
#include "../libLanlGeoMag/Lgm/Lgm_MagModelInfo.h"

/*
 *  Regression tests for closed field line calculations
 */


Lgm_MagModelInfo    *mInfo;

void ClosedField_Setup(void) {
    mInfo = Lgm_InitMagInfo();
    return;
}

void ClosedField_TearDown(void) {
    Lgm_FreeMagInfo( mInfo );
    return;
}


START_TEST(test_ClosedField_01) {
    Lgm_Vector        MinBexpect, Pos, MinBtest, Udiff, NFPexpect, SFPexpect, \
        NFPtest, SFPtest;
    int               nTests, nPass, nFail, transflag, SubtestPassed, \
        Passed=FALSE;
    int               retVal;
    double            del, Kp;
    char              buff[262], extModel[10], intModel[10], Resultexpect[20];
    char              *Resulttest;
    char              IsoDate[80];
    long long         TT2000;
    FILE              *testfile, *outfile;
    Lgm_DateTime      d;

    int makeNew = 1;

    /* read test file */
    testfile = fopen("check_ClosedField_01.expected","r");
    if (makeNew) outfile = fopen("check_ClosedField_01.got", "w");

    /* step through test cases one line at a time */
    nTests = 0;
    nPass = 0;
    nFail = 0;
    Passed = TRUE;
    while( fgets(buff,260,testfile) != NULL) {
        if (buff[0]!='#') {
	    SubtestPassed = TRUE;
            // read line
            sscanf(buff,
		   "%s %s %s %lf "
		   "%lf %lf %lf %s "
		   "%lf %lf %lf "
		   "%lf %lf %lf "
		   "%lf %lf %lf",
		   IsoDate, extModel, intModel, &Kp, 
                   &Pos.x, &Pos.y, &Pos.z, Resultexpect,
		   &NFPexpect.x, &NFPexpect.y, &NFPexpect.z,
		   &SFPexpect.x, &SFPexpect.y, &SFPexpect.z,
		   &MinBexpect.x, &MinBexpect.y, &MinBexpect.z);
	    Lgm_InitMagInfoDefaults( mInfo );
            IsoTimeStringToDateTime( IsoDate, &d, mInfo->c );
            Lgm_Set_Coord_Transforms( d.Date, d.Time, mInfo->c );
            nTests++;
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
	    else if (!strncmp(extModel, "OP77", 10)) {
	        Lgm_Set_Lgm_B_OP77(mInfo);
	    }
	    else {
	        nFail++;
	        printf("Test %d bad external model %s\n", nTests, extModel);
		continue;
	    }
	    mInfo->Kp = Kp;
	    retVal = Lgm_Trace(&Pos, &SFPtest, &NFPtest, &MinBtest,
			       100., 0.01, 1.e-7, mInfo);
            //printf("IsoDate = %s; d.Date, d.time = %ld, %lf \n", IsoDate, d.Date, d.Time);
	    switch (retVal) {
	        case LGM_CLOSED:
		    Resulttest = "LGM_CLOSED";
		    break;
	        default:
		    Resulttest = "UNKNOWN";
	    }	    
	    if (strncmp(Resultexpect, Resulttest, 19)) {
	        SubtestPassed = FALSE;
		printf("Test %d result %s, expected %s\n",
		       nTests, Resulttest, Resultexpect);
	    }
            Udiff.x = NFPtest.x - NFPexpect.x;
            Udiff.y = NFPtest.y - NFPexpect.y;
            Udiff.z = NFPtest.z - NFPexpect.z;
            del = Lgm_Magnitude(&Udiff);
            if (fabs(del) > 1.0e-5) {
	        SubtestPassed = FALSE;
                printf("*****  warning : NFP difference >= 1.0e-5 km  *****\n");
                printf("Test %d failed (diff: %g %g %g   %g)\n",
		       nTests, Udiff.x, Udiff.y, Udiff.z, fabs(del));
            }
            Udiff.x = SFPtest.x - SFPexpect.x;
            Udiff.y = SFPtest.y - SFPexpect.y;
            Udiff.z = SFPtest.z - SFPexpect.z;
            del = Lgm_Magnitude(&Udiff);
            if (fabs(del) > 1.0e-5) {
	        SubtestPassed = FALSE;
                printf("*****  warning : SFP difference >= 1.0e-5 km  *****\n");
                printf("Test %d failed (diff: %g %g %g   %g)\n",
		       nTests, Udiff.x, Udiff.y, Udiff.z, fabs(del));
            }
            Udiff.x = MinBtest.x - MinBexpect.x;
            Udiff.y = MinBtest.y - MinBexpect.y;
            Udiff.z = MinBtest.z - MinBexpect.z;
            del = Lgm_Magnitude(&Udiff);
            if (fabs(del) > 1.0e-5) {
	        SubtestPassed = FALSE;
                printf("*****  warning : MinB difference >= 1.0e-5 nT  *****\n");
                printf("Test %d failed (diff: %g %g %g   %g)\n",
		       nTests, Udiff.x, Udiff.y, Udiff.z, fabs(del));
            }
	    if (SubtestPassed) {
	        nPass++;
                printf("Test %d passed\n", nTests);
            }
	    else {
	        nFail++;
	    }
            if (makeNew) fprintf(
		   outfile,
		   "%s %s %s %lf "
		   "%lf %lf %lf %s "
		   "%lf %lf %lf "
		   "%lf %lf %lf "
		   "%lf %lf %lf\n",
		   IsoDate, extModel, intModel, Kp,
		   Pos.x, Pos.y, Pos.z, Resulttest,
		   NFPtest.x, NFPtest.y, NFPtest.z,
		   SFPtest.x, SFPtest.y, SFPtest.z,
		   MinBtest.x, MinBtest.y, MinBtest.z);
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

    ck_assert_msg( Passed, "ClosedField tests failed.\n" );

}
END_TEST


Suite *ClosedField_suite(void) {

  Suite *s = suite_create("CLOSEDFIELD_TESTS");

  TCase *tc_ClosedField = tcase_create("ClosedField Models");
  tcase_add_checked_fixture(tc_ClosedField, ClosedField_Setup,
			    ClosedField_TearDown);

  tcase_add_test(tc_ClosedField, test_ClosedField_01);

  suite_add_tcase(s, tc_ClosedField);

  return s;

}

int main(void) {

    int      number_failed;
    Suite   *s  = ClosedField_suite();
    SRunner *sr = srunner_create(s);

    printf("\n\n");
    srunner_run_all(sr, CK_ENV);
    number_failed = srunner_ntests_failed(sr);
    srunner_free(sr);

    return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;

}
