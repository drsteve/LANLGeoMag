#include <check.h>
#include "../libLanlGeoMag/Lgm/Lgm_LstarInfo.h"


Lgm_LstarInfo *LstarInfo;


void Lstar_Setup(void) {
    LstarInfo = InitLstarInfo(0);
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
    LstarInfo->VerbosityLevel = 1;
    LstarInfo->mInfo->VerbosityLevel = 0;

    /*
     *  Blocal at given location.
     */
    LstarInfo->mInfo->Bfield( &P, &Bvec, LstarInfo->mInfo );
    Blocal = Lgm_Magnitude( &Bvec );
    sa = sin( LstarInfo->PitchAngle*RadPerDeg ); sa2 = sa*sa;
    LstarInfo->mInfo->Bm = Blocal/sa2;
    printf("Bm = %g\n", LstarInfo->mInfo->Bm);


    Lstar( &P, LstarInfo );
    printf("L* = %g\n", LstarInfo->LS);

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
    LstarInfo->VerbosityLevel = 1;
    LstarInfo->mInfo->VerbosityLevel = 0;

    /*
     *  Blocal at given location.
     */
    LstarInfo->mInfo->Bfield( &P, &Bvec, LstarInfo->mInfo );
    Blocal = Lgm_Magnitude( &Bvec );
    sa = sin( LstarInfo->PitchAngle*RadPerDeg ); sa2 = sa*sa;
    LstarInfo->mInfo->Bm = Blocal/sa2;

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
    LstarInfo->VerbosityLevel = 1;
    LstarInfo->mInfo->VerbosityLevel = 0;

    // run 1
    LstarInfo->PitchAngle = 87.5;
    LstarInfo->mInfo->Bfield( &P, &Bvec, LstarInfo->mInfo );
    Blocal = Lgm_Magnitude( &Bvec );
    sa = sin( LstarInfo->PitchAngle*RadPerDeg ); sa2 = sa*sa;
    LstarInfo->mInfo->Bm = Blocal/sa2;
    Lstar( &P, LstarInfo );
    ans1 = LstarInfo->LS;

    // run 2
    LstarInfo->PitchAngle = 37.5;
    LstarInfo->mInfo->Bfield( &P, &Bvec, LstarInfo->mInfo );
    Blocal = Lgm_Magnitude( &Bvec );
    sa = sin( LstarInfo->PitchAngle*RadPerDeg ); sa2 = sa*sa;
    LstarInfo->mInfo->Bm = Blocal/sa2;
    Lstar( &P, LstarInfo );
    ans2 = LstarInfo->LS;

    // compare
    LstarDiff = fabs(ans2-ans1);
    tol = pow(20.0, (double) -quality );

    ck_assert_msg((LstarDiff<tol), "Roederer L differs in CDIP for different pitch angles. %g - %g = %g\n", ans2, ans1, LstarDiff);


    return;

}END_TEST


Suite *Lstar_suite(void) {

  Suite *s = suite_create("ROEDERER_L_TESTS");

  TCase *tc_Lstar = tcase_create("Roederer L tests");
  tcase_add_checked_fixture(tc_Lstar, Lstar_Setup, Lstar_TearDown);

  tcase_add_test(tc_Lstar, test_Lstar_CDIP);
  tcase_add_test(tc_Lstar, test_Lstar_CDIPapprox);
  tcase_add_test(tc_Lstar, test_Lstar_CDIPalpha);

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
