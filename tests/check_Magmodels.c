#include <check.h>
#include "../libLanlGeoMag/Lgm/Lgm_MagModelInfo.h"

/*
 *  Regression tests for magnetic field model calculations
 */


Lgm_MagModelInfo    *mInfo;

void Magmodels_Setup(void) {
    mInfo = Lgm_InitMagInfo();
    return;
}

void Magmodels_TearDown(void) {
    Lgm_FreeMagInfo( mInfo );
    return;
}


START_TEST(test_Magmodels_01) {
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
    testfile = fopen("check_Magmodels_01.expected","r");
    if (makeNew) outfile = fopen("check_Magmodels_01.got", "w");

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

	    nTests++;
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
	    else if (!strncmp(extModel, "OP77", 10)) {
	        Lgm_Set_Lgm_B_OP77(mInfo);
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
	        printf("Test %d odd return from Bfield\n", nTests);
	        continue;
		}
            //printf("IsoDate = %s; d.Date, d.time = %ld, %lf \n", IsoDate, d.Date, d.Time);
            Udiff.x = Btest.x - Bexpect.x;
            Udiff.y = Btest.y - Bexpect.y;
            Udiff.z = Btest.z - Bexpect.z;
            del = Lgm_Magnitude(&Udiff);
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

    ck_assert_msg( Passed, "Magmodels tests failed.\n" );

}
END_TEST


START_TEST(test_Magmodels_02) {
    Lgm_Vector        Bexpect, Pos, Btest, Udiff;
    int               nTests, nPass, nFail, transflag, Passed=FALSE;
    int               retVal;
    double            del, Pdyn, SymHc, Nidx, Byavg;
    char              buff[262], extModel[10], intModel[10];
    char              IsoDate[80];
    long long         TT2000;
    FILE              *testfile, *outfile;
    Lgm_DateTime      d;

    int makeNew = 1;

    /* read test file for TA2016*/
    testfile = fopen("check_Magmodels_02.expected","r");
    if (makeNew) outfile = fopen("check_Magmodels_02.got", "w");

    /* step through test cases one line at a time */
    nTests = 0;
    nPass = 0;
    nFail = 0;
    Passed = TRUE;
	Lgm_InitMagInfoDefaults( mInfo );
    Lgm_Init_TA16( &mInfo->TA16_Info, 1 );  // use verbose output
    while( fgets(buff,260,testfile) != NULL) {
        if (buff[0]!='#') {
            // read line
            sscanf(buff, "%s %s %s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
            //sscanf(buff,"%s %s %s %lf %lf %lf %lf %lf %lf %lf",
		       &IsoDate[0], &extModel[0], &intModel[0],
               &Pdyn, &SymHc, &Nidx, &Byavg,
               &Pos.x, &Pos.y, &Pos.z, &Bexpect.x, &Bexpect.y, &Bexpect.z);

	        nTests++;
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
	        if (!strncmp(extModel, "TA16", 10)) {
	            Lgm_Set_Lgm_B_TA16(mInfo);
                // Lgm_SetCoeffs_TA16( d.Date, d.Time, &mInfo->TA16_Info );

                mInfo->TA16_Info.Pdyn = Pdyn;
                mInfo->TA16_Info.SymHc_avg = SymHc;
                mInfo->TA16_Info.Xind_avg = Nidx;
                mInfo->TA16_Info.By_avg = Byavg;
	        }
	        else {
	            nFail++;
	            printf("Test %d bad external model %s\n", nTests, extModel);
		    continue;
	        }
	        retVal = mInfo->Bfield(&Pos, &Btest, mInfo);
            Udiff.x = Btest.x - Bexpect.x;
            Udiff.y = Btest.y - Bexpect.y;
            Udiff.z = Btest.z - Bexpect.z;
            del = Lgm_Magnitude(&Udiff);
            if (fabs(del) <= 1.0e-5) {
                nPass++;
                printf("Test %d passed\n", nTests);
                }
            else {
                nFail++;
                printf("*****  warning : difference >= 1.0e-5 nT  *****\n");
                printf("Test %d failed (diff: %g %g %g   %g)\n", nTests, Udiff.x, Udiff.y, Udiff.z, fabs(del));
            }
            if (makeNew) {
                fprintf(outfile, "%s %s %s %lf %lf %lf %lf %lf %lf %lf %11.8g %11.8g %11.8g\n",
                        IsoDate, extModel, intModel,
                        Pdyn, SymHc, Nidx, Byavg,
                        Pos.x, Pos.y, Pos.z, Btest.x, Btest.y, Btest.z);
            }
        }
        else {
            if (makeNew) {
                fprintf(outfile, "%s", buff);
            }
        }
    }
    if (nFail>0) Passed = FALSE;
    fclose(testfile);
    if (makeNew) fclose(outfile);
    printf("Result: %d tests pass; %d tests fail (Precision=1.0e-5 nT)\n", nPass, nFail);
    fflush(stdout);

    ck_assert_msg( Passed, "Magmodels tests failed.\n" );

}
END_TEST



Suite *Magmodels_suite(void) {

  Suite *s = suite_create("MAGMODELS_TESTS");

  TCase *tc_Magmodels = tcase_create("Magmodels Models");
  tcase_add_checked_fixture(tc_Magmodels, Magmodels_Setup, Magmodels_TearDown);

  tcase_add_test(tc_Magmodels, test_Magmodels_01);
  tcase_add_test(tc_Magmodels, test_Magmodels_02);

  suite_add_tcase(s, tc_Magmodels);

  return s;

}

int main(void) {

    int      number_failed;
    Suite   *s  = Magmodels_suite();
    SRunner *sr = srunner_create(s);

    printf("\n\n");
    srunner_run_all(sr, CK_ENV);
    number_failed = srunner_ntests_failed(sr);
    srunner_free(sr);

    return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;

}
