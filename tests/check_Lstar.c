#include <check.h>
#include "../libLanlGeoMag/Lgm/Lgm_LstarInfo.h"
#include "../libLanlGeoMag/Lgm/Lgm_MagEphemInfo.h"


Lgm_LstarInfo *LstarInfo;


void Lstar_Setup(void) {
    LstarInfo = InitLstarInfo(1);
    return;
}

void Lstar_TearDown(void) {
    FreeLstarInfo(LstarInfo);
    return;
}

START_TEST(test_Lstar_CDIP){
    /*
     *  Lstar in centered dipole should give known result
     */

    double           UTC, Blocal, sa, sa2, LstarDiff, tol;
    long int         Date;
    int              Kp, quality;
    Lgm_Vector       Psm, P, Bvec, v1, v2, v3;

    // Date and UTC
    Date       = 19891020;
    UTC        = 21.0;
    Lgm_Set_Coord_Transforms( Date, UTC, LstarInfo->mInfo->c );

    // Position in SM
    Psm.x = 0.0; Psm.z = 0.0;
    Psm.y = 7.0;
    Lgm_Convert_Coords( &Psm, &P, SM_TO_GSM, LstarInfo->mInfo->c );

    // pitch angles to compute
    LstarInfo->PitchAngle = 90.0;
    LstarInfo->LstarMoment = LGM_LSTAR_MOMENT_CDIP;
    Lgm_MagModelInfo_Set_MagModel( LGM_CDIP, LGM_EXTMODEL_NULL, LstarInfo->mInfo );
    quality = 3;
    Lgm_SetLstarTolerances( quality, 48, LstarInfo );

    Lstar( &P, LstarInfo );

    LstarDiff = fabs(LstarInfo->LS - Psm.y);

    tol = pow(10.0, (double) -quality );
    ck_assert((LstarDiff<tol));

    return;

}END_TEST


START_TEST(test_Lstar_CDIPapprox){
    /*
     *  Lstar in centered dipole should give same result with full field and dipole approx
     */

    double           UTC, Blocal, sa, sa2, LstarDiff, tol;
    long int         Date;
    int              Kp, quality;
    Lgm_Vector       Psm, P, Bvec, v1, v2, v3;

    // Date and UTC
    Date       = 20110123;
    UTC        = 3.14;
    Lgm_Set_Coord_Transforms( Date, UTC, LstarInfo->mInfo->c );

    // Position in SM
    Psm.x = 0.0; Psm.z = 0.0;
    Psm.y = 6.0;
    Lgm_Convert_Coords( &Psm, &P, SM_TO_GSM, LstarInfo->mInfo->c );

    // pitch angles to compute
    LstarInfo->PitchAngle = 75.0;
    LstarInfo->LstarMoment = LGM_LSTAR_MOMENT_CDIP;
    Lgm_MagModelInfo_Set_MagModel( LGM_CDIP, LGM_EXTMODEL_NULL, LstarInfo->mInfo );
    quality = 3;
    Lgm_SetLstarTolerances( quality, 36, LstarInfo );

    Lstar( &P, LstarInfo );

    LstarDiff = fabs(LstarInfo->LS - LstarInfo->LS_dip_approx);

    tol = pow(20.0, (double) -quality );
    ck_assert_msg((LstarDiff<tol), "LS = %g; LS_dip_approx = %g; LstarDiff = %g; Mused = %g; Mcd = %g; Mcd2010 = %g\n", LstarInfo->LS, LstarInfo->LS_dip_approx, LstarDiff, LstarInfo->Mused, LstarInfo->mInfo->c->M_cd, LstarInfo->mInfo->c->M_cd_2010);

    return;

}END_TEST


START_TEST(test_Lstar_CDIPalpha){
    /*
     *  Same pitch angle in centered dipole should give same result
     */

    double           UTC, Blocal, sa, sa2, LstarDiff, ans1, ans2, tol;
    long int         Date;
    int              Kp, quality;
    Lgm_Vector       Psm, P, Bvec, v1, v2, v3;

    // Date and UTC
    Date       = 19991122;
    UTC        = 19.0;
    Lgm_Set_Coord_Transforms( Date, UTC, LstarInfo->mInfo->c );

    // Position in SM
    Psm.x = 0.0; Psm.y = 6.6; Psm.z = 0.0;
    Lgm_Convert_Coords( &Psm, &P, SM_TO_GSM, LstarInfo->mInfo->c );

    // set up Lstar calcs
    LstarInfo->LstarMoment = LGM_LSTAR_MOMENT_CDIP_2010;
    Lgm_MagModelInfo_Set_MagModel( LGM_CDIP, LGM_EXTMODEL_NULL, LstarInfo->mInfo );
    Lgm_SetLstarTolerances( 5, 24, LstarInfo );

    // run 1
    LstarInfo->PitchAngle = 87.5;
    Lstar( &P, LstarInfo );
    ans1 = LstarInfo->LS;

    // run 2
    LstarInfo->PitchAngle = 37.5;
    Lstar( &P, LstarInfo );
    ans2 = LstarInfo->LS;

    // compare
    LstarDiff = fabs(ans2-ans1);
    tol = pow(20.0, (double) -quality );

    ck_assert_msg((LstarDiff<tol), "Roederer L differs in CDIP for different pitch angles. %g - %g = %g\n", ans2, ans1, LstarDiff);


    return;

}END_TEST


START_TEST(test_Lstar_CDIPalpha2){
    /*
     *  Should get same result from Lstar and ComputeLstarVersusPA
     */

    double           UTC, Blocal, sa, sa2, LstarDiff, ans2, tol;
    double           PAs[2] = {87.5, 37.5};
    long int         Date;
    int              Kp, quality=3;
    Lgm_Vector       Psm, P, Bvec;
    Lgm_MagEphemInfo *MagEphemInfo = Lgm_InitMagEphemInfo(1, 2);

    // Date, UTC, position
    Date       = 19991122;
    UTC        = 19.0;
    Psm.x = 0.0; Psm.y = 6.6; Psm.z = 0.0;

    /*
     * Calculate using Lstar directly
     */
    Lgm_Set_Coord_Transforms( Date, UTC, LstarInfo->mInfo->c );

    // Position in SM
    Lgm_Convert_Coords( &Psm, &P, SM_TO_GSM, LstarInfo->mInfo->c );

    // set up Lstar calcs
    LstarInfo->LstarMoment = LGM_LSTAR_MOMENT_CDIP_2010;
    Lgm_MagModelInfo_Set_MagModel( LGM_EDIP, LGM_EXTMODEL_NULL, LstarInfo->mInfo );
    Lgm_SetLstarTolerances( quality, 72, LstarInfo );

    // Get Lstar
    LstarInfo->PitchAngle = PAs[1];
    Lstar( &P, LstarInfo );
    ans2 = LstarInfo->LS;

    /*
     * Calculate using Lgm_ComputeLstarVersusPA
     */

    Lgm_SetMagEphemLstarQuality( quality, 72, MagEphemInfo );
    MagEphemInfo->LstarInfo->LstarMoment = LGM_LSTAR_MOMENT_CDIP_2010;
    Lgm_MagModelInfo_Set_MagModel( LGM_EDIP, LGM_EXTMODEL_NULL, MagEphemInfo->LstarInfo->mInfo );

    Lgm_ComputeLstarVersusPA( Date, UTC, &P, 2, PAs, FALSE, MagEphemInfo );

    printf("Lstar (from Lstar) = %g\nLstar (from ComputeLstarVersusPA) = %g\n", ans2, MagEphemInfo->Lstar[1]);

    // compare
    LstarDiff = fabs(ans2 - MagEphemInfo->Lstar[1]);
    tol = pow(20.0, (double) -quality );

    ck_assert_msg((LstarDiff<tol), "Roederer L differs between Lstar and ComputeLstarVersusPA. %g - %g = %g\n", ans2, MagEphemInfo->Lstar[1], LstarDiff);


    return;

}END_TEST


START_TEST(test_Lstar_Regressions) {
    /* Regression tests against previous L* results */
    Lgm_Vector        Pos, PosGSM;
    int               nTests, nPass, nFail, transflag, SubtestPassed, \
                      Passed=FALSE;
    double            del, Kp, PA, HiltonExpect, HiltonTest, McIlwainExpect, \
                      McIlwainTest, RoedererExpect, RoedererTest, I, Bm, M;
    char              buff[262], extModel[10], intModel[10];
    char              IsoDate[80];
    FILE              *testfile, *outfile;
    Lgm_DateTime      d;

    int makeNew = 1;

    /* read test file */
    testfile = fopen("check_Lstar.expected","r");
    if (makeNew) outfile = fopen("check_Lstar.got", "w");

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
		   "%lf %lf %lf "
		   "%lf %lf %lf %lf",
		   IsoDate, extModel, intModel, &Kp, 
                   &Pos.x, &Pos.y, &Pos.z,
		   &PA, &HiltonExpect, &McIlwainExpect, &RoedererExpect);
            IsoTimeStringToDateTime( IsoDate, &d, LstarInfo->mInfo->c );
            //printf("IsoDate = %s; d.Date, d.Time = %ld, %lf \n", IsoDate, d.Date, d.Time);
            Lgm_Set_Coord_Transforms( d.Date, d.Time, LstarInfo->mInfo->c );
	    //Assume input is SM (since it is in Python tests)
	    Lgm_Convert_Coords(&Pos, &PosGSM, SM_TO_GSM, LstarInfo->mInfo->c );
            nTests++;
	    if (!strncmp(intModel, "IGRF", 10)) {
	        Lgm_Set_Lgm_B_IGRF_InternalModel( LstarInfo->mInfo );
	    }
	    else if (!strncmp(intModel, "CDIP", 10)) {
	        Lgm_Set_Lgm_B_cdip_InternalModel( LstarInfo->mInfo );
	    }
	    else if (!strncmp(intModel, "EDIP", 10)) {
	        Lgm_Set_Lgm_B_edip_InternalModel( LstarInfo->mInfo );
	    }
	    else {
	        nFail++;
	        printf("Test %d bad internal model %s\n", nTests, intModel);
		continue;
	    }
	    if (!strncmp(extModel, "T89", 10)) {
	        Lgm_Set_Lgm_B_T89(LstarInfo->mInfo);
	    }
	    else if (!strncmp(extModel, "OP77", 10)) {
	        Lgm_Set_Lgm_B_OP77(LstarInfo->mInfo);
	    }
	    else if (!strncmp(extModel, "NULL", 10)) {
	        LstarInfo->mInfo->ExternalModel = LGM_EXTMODEL_NULL;
	    }
	    else {
	        nFail++;
	        printf("Test %d bad external model %s\n", nTests, extModel);
		continue;
	    }
	    LstarInfo->mInfo->Kp = Kp;
	    LstarInfo->PitchAngle = PA;
	    Lgm_SetLstarTolerances( 1, 24, LstarInfo );
	    //	    LstarInfo->ShabanskyHandling = LGM_SHABANSKY_IGNORE;
	    McIlwainTest = Lgm_McIlwain_L(d.Date, d.Time, &PosGSM,
					  PA, 0, &I, &Bm, &M, LstarInfo->mInfo);
	    HiltonTest = Lgm_McIlwain_L(d.Date, d.Time, &PosGSM,
					PA, 1, &I, &Bm, &M, LstarInfo->mInfo);
	    Lstar(&PosGSM, LstarInfo);
	    RoedererTest = LstarInfo->LS;
	    del = HiltonTest - HiltonExpect;
	    if (fabs(del) > 1.0e-5) {
	        SubtestPassed = FALSE;
                printf("*****  warning : Hilton difference >= 1.0e-5 *****\n");
                printf("Test %d failed (diff: %g)\n", nTests, del);
            }
	    del = McIlwainTest - McIlwainExpect;
	    if (fabs(del) > 1.0e-5) {
	        SubtestPassed = FALSE;
                printf("*****  warning : McIlwain difference >= 1.0e-5 *****\n");
                printf("Test %d failed (diff: %g)\n", nTests, del);
            }
	    del = RoedererTest - RoedererExpect;
	    if (fabs(del) > 1.0e-5) {
	        SubtestPassed = FALSE;
                printf("*****  warning : Roederer difference >= 1.0e-5 *****\n");
                printf("Test %d failed (diff: %g)\n", nTests, del);
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
		   "%lf %lf %lf "
		   "%lf %lf %lf %lf\n",
		   IsoDate, extModel, intModel, Kp,
		   Pos.x, Pos.y, Pos.z,
		   PA, HiltonTest, McIlwainTest, RoedererTest);
            }
        else {
            if (makeNew) fprintf(outfile, "%s", buff);
            }
        }
    if (nFail>0) Passed = FALSE;
    fclose(testfile);
    if (makeNew) fclose(outfile);
    printf("Result: %d tests pass; %d tests fail (Precision=1.0e-5)\n", nPass, nFail);
    fflush(stdout);

    ck_assert_msg( Passed, "LstarRegressions tests failed.\n" );

}
END_TEST


Suite *Lstar_suite(void) {

  Suite *s = suite_create("ROEDERER_L_TESTS");

  TCase *tc_Lstar = tcase_create("Roederer L tests");
  tcase_add_checked_fixture(tc_Lstar, Lstar_Setup, Lstar_TearDown);

  tcase_add_test(tc_Lstar, test_Lstar_CDIP);
  tcase_add_test(tc_Lstar, test_Lstar_CDIPapprox);
  tcase_add_test(tc_Lstar, test_Lstar_CDIPalpha);
  tcase_add_test(tc_Lstar, test_Lstar_CDIPalpha2);
  tcase_add_test(tc_Lstar, test_Lstar_Regressions);

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
