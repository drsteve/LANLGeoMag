#include <stdio.h>
#include <string.h>
#include <check.h>
#include "../libLanlGeoMag/Lgm/Lgm_CTrans.h"
#include "../libLanlGeoMag/Lgm/Lgm_Vec.h"

#define TRUE    1
#define FALSE   0


int getSys( char *sys );


START_TEST(test_CoordTrans) {

    Lgm_CTrans        *c = Lgm_init_ctrans( 0 ); 
    Lgm_Vector        Utarg, Ucent, Utest, Udiff;
    int               nTests, line, nPass, nFail, transflag, Passed=FALSE;
    double            del;
    char              buff[262], sysIn[10], sysOut[10];
    char              IsoDate[80];
    long long         TT2000;
    FILE              *testfile, *outfile;
    Lgm_DateTime      d;

    int makeNew = 1;
    Lgm_Set_CTrans_Options(LGM_EPH_DE, LGM_PN_IAU76, c);

    /* read test file */
    testfile = fopen("check_CoordTrans.expected","r");
    if (makeNew) outfile = fopen("check_CoordTrans.got", "w");

    /* step through test cases one line at a time */
    line = 0;
    nTests = 0;
    nPass = 0;
    nFail = 0;
    Passed = TRUE;
    while( fgets(buff,260,testfile) != NULL) {
        //if (line>=15) exit(0);
        if (buff[0]!='#') {
            // read line
            sscanf(buff,"%s %lld %s %s %lf %lf %lf %lf %lf %lf", &IsoDate[0], &TT2000, &sysIn[0], &sysOut[0], 
                                                                 &Ucent.x, &Ucent.y, &Ucent.z, &Utarg.x, &Utarg.y, &Utarg.z);
            line++;

            // Set up all the necessary variables to do transformations for this Date and UTC
            IsoTimeStringToDateTime( IsoDate, &d, c );
            Lgm_Set_Coord_Transforms( d.Date, d.Time, c );
            //printf("IsoDate = %s; d.Date, d.time = %ld, %lf \n", IsoDate, d.Date, d.Time);
            //get transformed coordinate
            transflag = getSys(sysIn)*100 + getSys(sysOut);
            Lgm_Convert_Coords( &Ucent, &Utest, transflag, c );

            //then test difference
            Udiff.x = Utest.x - Utarg.x;
            Udiff.y = Utest.y - Utarg.y;
            Udiff.z = Utest.z - Utarg.z;
            del = Lgm_Magnitude(&Udiff);
            nTests++;
            if (fabs(del) <= 1.0e-5) {
                nPass++;
                printf("Test %d passed\n", nTests);
                }
            else {
                nFail++;
                printf("*****  warning : difference >= 1.0e-5 km (1 cm)  *****\n");
                printf("Test %d failed (diff: %g %g %g   %g)\n", nTests, Udiff.x, Udiff.y, Udiff.z, fabs(del));
                }
            if (makeNew) fprintf(outfile, "%s %lld %s %s %lf %lf %lf %lf %lf %lf\n", IsoDate, TT2000, sysIn, sysOut, Ucent.x, Ucent.y, Ucent.z, Utest.x, Utest.y, Utest.z);

            }
        else {
            if (makeNew) fprintf(outfile, "%s", buff);
            }
        }
    if (nFail>0) Passed = FALSE;
    fclose(testfile);
    if (makeNew) fclose(outfile);
    printf("Result: %d tests pass; %d tests fail (Precision=1.0e-5 km [=1cm])\n", nPass, nFail);
    Lgm_free_ctrans( c ); // free the structure

    fail_unless( Passed, "CoordTrans test failed. Have the leap-seconds changed?\n" );

    return;
    } END_TEST


START_TEST(test_CoordRoundtrip) {

    Lgm_CTrans        *c = Lgm_init_ctrans( 0 ); 
    Lgm_Vector        Utarg, Umid, Umid2, Utest, Udiff;
    int               nTests, nPass, nFail, transflag, ntest, Passed=FALSE;
    int               n=0, sysIn[5]={1,3,4,5,7}, sysMid[5]={2,6,6,4,1};
    long int          Date;
    double            UTC, del;

    Lgm_Set_CTrans_Options(LGM_EPH_DE, LGM_PN_IAU76, c);

    /* step through test cases one line at a time */
    ntest = 0;
    nPass = 0;
    nFail = 0;
    Passed = TRUE;

    Date = 20000101;    // Jan 1, 2000
    UTC  = 0.0+0.0/60.0;         // Universal Time Coordinated (in decimal hours)

    Utarg.x = -25.275; Utarg.y = 9.735; Utarg.z = 0.000; // Set a vector in GSM coordinates
    Utarg.x =   5.271; Utarg.y = 3.535; Utarg.z = 1.362; // Set a vector in GSM coordinates
    Utarg.x =   4.427; Utarg.y = 0.773; Utarg.z = 3.608; // Set a vector in GSM coordinates

    printf("\nCoordinate roundtrip tests:\n");
    printf("Test %d. %d -> %d -> %d\n", ntest, sysIn[ntest], sysMid[ntest], sysIn[ntest]);
    // Set up all the necessary variables to do transformations for this Date and UTC
    Lgm_Set_CTrans_Options(LGM_EPH_HIGH_ACCURACY, LGM_PN_IAU76, c); /* Uses LGM high accuracy analytic solution */
    Lgm_Set_Coord_Transforms( Date, UTC, c );

    // Do the transformation to intermediate system and back out
    transflag = sysIn[ntest]*100 + sysMid[ntest];
    Lgm_Convert_Coords( &Utarg, &Umid, transflag, c );
    transflag = sysMid[ntest]*100 + sysIn[ntest];
    Lgm_Convert_Coords( &Umid, &Utest, transflag, c );
    printf("Roundtrip from system %d to %d to %d\n", sysIn[ntest], sysMid[ntest], sysIn[ntest]);
    nFail = nFail + (testDiff(Utest, Utarg, Umid));
    ntest++;


    // **** Test 2 **** //
    printf("Test %d. %d -> %d -> %d\n", ntest, sysIn[ntest], sysMid[ntest], sysIn[ntest]);
    // Set up all the necessary variables to do transformations for this Date and UTC
    Lgm_Set_CTrans_Options(LGM_EPH_DE, LGM_PN_IAU76, c); /* Uses JPL Development Ephemeris */
    Lgm_Set_Coord_Transforms( Date, UTC, c );

    // Do the transformation to intermediate system and back out
    transflag = sysIn[ntest]*100 + sysMid[ntest];
    Lgm_Convert_Coords( &Utarg, &Umid, transflag, c );
    transflag = sysMid[ntest]*100 + sysIn[ntest];
    Lgm_Convert_Coords( &Umid, &Utest, transflag, c );
    printf("Roundtrip from system %d to %d to %d\n", sysIn[ntest], sysMid[ntest], sysIn[ntest]);
    nFail = nFail + (testDiff(Utest, Utarg, Umid));
    ntest++;


    // **** Test 3 (multiple roundtrips) **** //
    printf("Test %d. %d -> %d -> %d -> %d -> %d\n", ntest, sysIn[ntest], sysMid[ntest], sysIn[ntest], sysMid[ntest], sysIn[ntest]);
    // Set up all the necessary variables to do transformations for this Date and UTC
    Lgm_Set_CTrans_Options(LGM_EPH_HIGH_ACCURACY, LGM_PN_IAU76, c); /* Uses LGM high accuracy analytic solution */
    Lgm_Set_Coord_Transforms( Date, UTC, c );

    // Do the transformation to intermediate system and back out
    transflag = sysIn[ntest]*100 + sysMid[ntest];
    Lgm_Convert_Coords( &Utarg, &Umid, transflag, c );
    transflag = sysMid[ntest]*100 + sysIn[ntest];
    Lgm_Convert_Coords( &Umid, &Utest, transflag, c );
    transflag = sysIn[ntest]*100 + sysMid[ntest];
    Lgm_Convert_Coords( &Utest, &Umid, transflag, c );
    transflag = sysMid[ntest]*100 + sysIn[ntest];
    Lgm_Convert_Coords( &Umid, &Utest, transflag, c );
    printf("Roundtrip from system %d to %d to %d to %d to %d\n", sysIn[ntest], sysMid[ntest], sysIn[ntest], sysMid[ntest], sysIn[ntest]);
    nFail = nFail + (testDiff(Utest, Utarg, Umid));
    ntest++;

    // **** Test 4 (via multiple) **** //
    printf("Test %d. %d -> %d -> 1 -> %d\n", ntest, sysIn[ntest], sysMid[ntest], sysIn[ntest]);
    // Set up all the necessary variables to do transformations for this Date and UTC
    Lgm_Set_CTrans_Options(LGM_EPH_HIGH_ACCURACY, LGM_PN_IAU76, c); /* Uses LGM high accuracy analytic solution */
    Lgm_Set_Coord_Transforms( Date, UTC, c );

    // Do the transformation to intermediate system and back out
    transflag = sysIn[ntest]*100 + sysMid[ntest];
    Lgm_Convert_Coords( &Utarg, &Umid, transflag, c );
    transflag = sysMid[ntest]*100 + 1;
    Lgm_Convert_Coords( &Umid, &Umid2, transflag, c );
    transflag = 100 + sysIn[ntest];
    Lgm_Convert_Coords( &Umid2, &Utest, transflag, c );
    printf("Roundtrip from system %d to %d to 1 to %d\n", sysIn[ntest], sysMid[ntest], sysIn[ntest]);
    nFail = nFail + (testDiff(Utest, Utarg, Umid));
    ntest++;

    // **** Test 5 **** //
    printf("Test %d. %d -> %d -> %d\n", ntest, sysIn[ntest], sysMid[ntest], sysIn[ntest]);
    // Set up all the necessary variables to do transformations for this Date and UTC
    Lgm_Set_CTrans_Options(LGM_EPH_HIGH_ACCURACY, LGM_PN_IAU76, c); /* Uses LGM high accuracy analytic solution */
    Lgm_Set_Coord_Transforms( Date, UTC, c );

    // Do the transformation to intermediate system and back out
    transflag = sysIn[ntest]*100 + sysMid[ntest];
    Lgm_Convert_Coords( &Utarg, &Umid, transflag, c );
    transflag = sysMid[ntest]*100 + sysIn[ntest];
    Lgm_Convert_Coords( &Umid, &Utest, transflag, c );
    printf("Roundtrip from system %d to %d to %d\n", sysIn[ntest], sysMid[ntest], sysIn[ntest]);
    nFail = nFail + (testDiff(Utest, Utarg, Umid));
    ntest++;


    if (nFail>0) Passed = FALSE;
    printf("Result: %d tests pass; %d tests fail \n", nPass, nFail);
    Lgm_free_ctrans( c ); // free the structure

    fail_unless( Passed, "CoordRoundtrip test failed.\n" );

    return;
    } END_TEST


int testDiff(Lgm_Vector Utest, Lgm_Vector Utarg, Lgm_Vector Umid) {
    Lgm_Vector  Udiff;
    double      del;
    int         Fail = FALSE;

    //then test difference
    Udiff.x = Utest.x - Utarg.x;
    Udiff.y = Utest.y - Utarg.y;
    Udiff.z = Utest.z - Utarg.z;
    del = Lgm_Magnitude(&Udiff);
    printf("X: In (%g), Out(%g), Diff (%g)\n", Utarg.x, Utest.x, Udiff.x);
    printf("Y: In (%g), Out(%g), Diff (%g)\n", Utarg.y, Utest.y, Udiff.y);
    printf("Z: In (%g), Out(%g), Diff (%g)\n", Utarg.z, Utest.z, Udiff.z);

    if (fabs(del) <= 1.0e-9) {
        printf("Test passed\n");
        }
    else {
        Fail = TRUE;
        printf("Test failed (diff: %g %g %g   %g)\n", Udiff.x, Udiff.y, Udiff.z, fabs(del));
    }

    return Fail;
    }


int getSys( char *sys ) {
    int out;
    if (!strcmp(sys, "GEI2000")){
        out = GEI2000_COORDS;
        } 
    else if (!strcmp(sys, "GEO")) {
        out = GEO_COORDS;
        }
    else if (!strcmp(sys, "GSE")) {
        out = GSE_COORDS;
        }
    else if (!strcmp(sys, "GSM")) {
        out = GSM_COORDS;
        }
    else if (!strcmp(sys, "SM")) {
        out = SM_COORDS;
        }
    else if (!strcmp(sys, "GSE2000")) {
        out = GSE2000_COORDS;
        }
    else if (!strcmp(sys, "CDMAG")) {
        out = CDMAG_COORDS;
        }
    else if (!strcmp(sys, "EDMAG")) {
        out = EDMAG_COORDS;
        }
    return(out);
    }


Suite *CT_suite(void) {

  Suite *s = suite_create("CoordTrans_TESTS");

  TCase *tc_CoordTrans = tcase_create("Coordinate Transformations");

  tcase_add_test(tc_CoordTrans, test_CoordTrans);
  tcase_add_test(tc_CoordTrans, test_CoordRoundtrip);

  suite_add_tcase(s, tc_CoordTrans);

  return s;

}

int main(void) {

    int      number_failed;
    Suite   *s  = CT_suite();
    SRunner *sr = srunner_create(s);

    printf("\n\n");
    srunner_run_all(sr, CK_NORMAL);
    number_failed = srunner_ntests_failed(sr);
    srunner_free(sr);

    return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;

}
