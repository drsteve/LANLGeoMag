#include <stdio.h>
#include <string.h>
#include <check.h>
#include "../libLanlGeoMag/Lgm/Lgm_CTrans.h"
#include "../libLanlGeoMag/Lgm/Lgm_Vec.h"

#define TRUE    1
#define FALSE   0


int getSys( char *sys );
int testDiff( Lgm_Vector Utest, Lgm_Vector Utarg, double tol );


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

    ck_assert_msg( Passed, "CoordTrans test failed. Have the leap-seconds changed?\n" );

    return;
    } END_TEST


START_TEST(test_CoordTransNoEph) {
    /* Coordinate transformations without using the JPL DE421*/
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

    /* read test file */
    testfile = fopen("check_CoordTransNoEph.expected","r");
    if (makeNew) outfile = fopen("check_CoordTransNoEph.got", "w");

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
    fflush(stdout); // get all that status info out before ck_assert
    Lgm_free_ctrans( c ); // free the structure

    ck_assert_msg( Passed, "CoordTransNoEph test failed. Have the leap-seconds changed?\n" );

    return;
    } END_TEST


START_TEST(test_CoordRoundtrip) {

    Lgm_CTrans        *c = Lgm_init_ctrans( 0 ); 
    Lgm_Vector        Utarg, Umid, Umid2, Utest, Udiff;
    int               nTests, nPass, nFail, transflag, ntest, Passed=FALSE;
    int               n=0, sysIn[6]={1,3,4,5,7,11}, sysMid[6]={2,6,6,4,1,8};
    long int          Date;
    double            UTC, del, tol=1e-10;

    Lgm_Set_CTrans_Options(LGM_EPH_DE, LGM_PN_IAU76, c);

    /* step through test cases one line at a time */
    ntest = 0;
    nPass = 0;
    nFail = 0;
    Passed = TRUE;

    Date = 20120501;    // May 1, 2021
    UTC  = 0.0+0.0/60.0;         // Universal Time Coordinated (in decimal hours)

    Utarg.x =   4.427; Utarg.y = 0.773; Utarg.z = 3.608; // Set a vector in input coordinates


    // **** Test 1 **** //
    printf("\nCoordinate roundtrip tests:\n");
    printf("Test %d. %d -> %d -> %d\n", ntest+1, sysIn[ntest], sysMid[ntest], sysIn[ntest]);
    // Set up all the necessary variables to do transformations for this Date and UTC
    Lgm_Set_CTrans_Options(LGM_EPH_HIGH_ACCURACY, LGM_PN_IAU76, c); /* Uses LGM high accuracy analytic solution */
    Lgm_Set_Coord_Transforms( Date, UTC, c );

    // Do the transformation to intermediate system and back out
    transflag = sysIn[ntest]*100 + sysMid[ntest];
    Lgm_Convert_Coords( &Utarg, &Umid, transflag, c );
    transflag = sysMid[ntest]*100 + sysIn[ntest];
    Lgm_Convert_Coords( &Umid, &Utest, transflag, c );
    printf("Roundtrip from system %d to %d to %d\n", sysIn[ntest], sysMid[ntest], sysIn[ntest]);
    nFail = nFail + (testDiff(Utest, Utarg, tol));
    ntest++;


    // **** Test 2 **** //
    printf("Test %d. %d -> %d -> %d\n", ntest+1, sysIn[ntest], sysMid[ntest], sysIn[ntest]);
    // Set up all the necessary variables to do transformations for this Date and UTC
    Lgm_Set_CTrans_Options(LGM_EPH_DE, LGM_PN_IAU76, c); /* Uses JPL Development Ephemeris */
    Lgm_Set_Coord_Transforms( Date, UTC, c );

    // Do the transformation to intermediate system and back out
    transflag = sysIn[ntest]*100 + sysMid[ntest];
    Lgm_Convert_Coords( &Utarg, &Umid, transflag, c );
    transflag = sysMid[ntest]*100 + sysIn[ntest];
    Lgm_Convert_Coords( &Umid, &Utest, transflag, c );
    printf("Roundtrip from system %d to %d to %d\n", sysIn[ntest], sysMid[ntest], sysIn[ntest]);
    nFail = nFail + (testDiff(Utest, Utarg, tol));
    ntest++;


    // **** Test 3 (multiple roundtrips) **** //
    printf("Test %d. %d -> %d -> %d -> %d -> %d\n", ntest+1, sysIn[ntest], sysMid[ntest], sysIn[ntest], sysMid[ntest], sysIn[ntest]);
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
    nFail = nFail + (testDiff(Utest, Utarg, tol));
    ntest++;

    // **** Test 4 (via multiple) **** //
    printf("Test %d. %d -> %d -> 1 -> %d\n", ntest+1, sysIn[ntest], sysMid[ntest], sysIn[ntest]);
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
    nFail = nFail + (testDiff(Utest, Utarg, tol));
    ntest++;

    // **** Test 5 **** //
    printf("Test %d. %d -> %d -> %d\n", ntest+1, sysIn[ntest], sysMid[ntest], sysIn[ntest]);
    // Set up all the necessary variables to do transformations for this Date and UTC
    Lgm_Set_CTrans_Options(LGM_EPH_HIGH_ACCURACY, LGM_PN_IAU76, c); /* Uses LGM high accuracy analytic solution */
    Lgm_Set_Coord_Transforms( Date, UTC, c );

    // Do the transformation to intermediate system and back out
    transflag = sysIn[ntest]*100 + sysMid[ntest];
    Lgm_Convert_Coords( &Utarg, &Umid, transflag, c );
    transflag = sysMid[ntest]*100 + sysIn[ntest];
    Lgm_Convert_Coords( &Umid, &Utest, transflag, c );
    printf("Roundtrip from system %d to %d to %d\n", sysIn[ntest], sysMid[ntest], sysIn[ntest]);
    nFail = nFail + (testDiff(Utest, Utarg, tol));
    ntest++;

    // **** Test 6 **** //
    printf("Test %d. %d -> %d -> 2 -> %d\n", ntest+1, sysIn[ntest], sysMid[ntest], sysIn[ntest]);
    // Set up all the necessary variables to do transformations for this Date and UTC
    Lgm_Set_CTrans_Options(LGM_EPH_LOW_ACCURACY, LGM_PN_IAU76, c); /* Uses LGM low accuracy analytic solution */
    Lgm_Set_Coord_Transforms( Date, UTC, c );

    // Do the transformation to intermediate system and back out
    transflag = sysIn[ntest]*100 + sysMid[ntest];
    Lgm_Convert_Coords( &Utarg, &Umid, transflag, c );
    transflag = sysMid[ntest]*100 + 2;
    Lgm_Convert_Coords( &Umid, &Umid2, transflag, c );
    transflag = 200 + sysIn[ntest];
    Lgm_Convert_Coords( &Umid2, &Utest, transflag, c );
    printf("Roundtrip from system %d to %d to 2 to %d\n", sysIn[ntest], sysMid[ntest], sysIn[ntest]);
    nFail = nFail + (testDiff(Utest, Utarg, tol));
    ntest++;
    

    // **** Tests complete *** //
    if (nFail>0) Passed = FALSE;
    printf("Result: %d tests pass; %d tests fail \n", ntest-nFail, nFail);
    Lgm_free_ctrans( c ); // free the structure

    ck_assert_msg( Passed, "CoordRoundtrip test failed.\n" );

    return;
    } END_TEST


START_TEST(test_CoordGSE_equiv) {
    Lgm_CTrans        *c = Lgm_init_ctrans( 0 ); 
    Lgm_Vector        Ufrom, Utest1, Utest2, Udiff;
    int               nTests, nPass, nFail, ntest, Passed=FALSE;
    int               n=0, sysFrom[3]={2,6,9};
    long int          Date;
    double            UTC, del, tol=1e-9;

    Lgm_Set_CTrans_Options(LGM_EPH_DE, LGM_PN_IAU76, c);

    /* step through test cases one line at a time */
    ntest = 0;
    nPass = 0;
    nFail = 0;
    Passed = TRUE;

    Ufrom.x =   5.271324; Ufrom.y = 3.535221; Ufrom.z = 1.362777; // Set a vector in input coordinates

    printf("\nGSE/GSE2000 tests (equivalence at J2000):\n");

    // **** Test 1 **** //
    printf("Test %d. Equivalence of GSE and GSE2000 at J2000 (DE)\n", ntest+1);
    // Set up all the necessary variables to do transformations for this Date and UTC
    Lgm_Set_CTrans_Options(LGM_EPH_DE, LGM_PN_IAU76, c); /* Uses LGM high accuracy analytic solution */
    Lgm_Set_Coord_Transforms( 20000101, 12.0, c );

    // Do the transformations from input system to both target systems
    Lgm_Convert_Coords( &Ufrom, &Utest1, MOD_TO_GSE, c );
    Lgm_Convert_Coords( &Ufrom, &Utest2, MOD_TO_GSE2000, c );
    nFail = nFail + (testDiff(Utest1, Utest2, tol));
    ntest++;

    // **** Test 2 **** //
    printf("Test %d. Equivalence of GSE and GSE2000 at J2000 (HA)\n", ntest+1);
    // Set up all the necessary variables to do transformations for this Date and UTC
    Lgm_Set_CTrans_Options(LGM_EPH_HIGH_ACCURACY, LGM_PN_IAU76, c); /* Uses LGM high accuracy analytic solution */
    Lgm_Set_Coord_Transforms( 20000101, 12.0, c );

    // Do the transformation to intermediate system and back out
    Lgm_Convert_Coords( &Ufrom, &Utest1, MOD_TO_GSE, c );
    Lgm_Convert_Coords( &Ufrom, &Utest2, MOD_TO_GSE2000, c );
    nFail = nFail + (testDiff(Utest1, Utest2, tol));
    ntest++;

    // **** Test 3 **** //
    printf("Test %d. Equivalence of GSE and GSE2000 at J2000 (LA)\n", ntest+1);
    // Set up all the necessary variables to do transformations for this Date and UTC
    Lgm_Set_CTrans_Options(LGM_EPH_LOW_ACCURACY, LGM_PN_IAU76, c); /* Uses LGM low accuracy analytic solution */
    Lgm_Set_Coord_Transforms( 20000101, 12.0, c );

    // Do the transformation to intermediate system and back out
    Lgm_Convert_Coords( &Ufrom, &Utest1, MOD_TO_GSE, c );
    Lgm_Convert_Coords( &Ufrom, &Utest2, MOD_TO_GSE2000, c );
    nFail = nFail + (testDiff(Utest1, Utest2, tol));
    ntest++;


    // **** Tests complete *** //
    if (nFail>0) Passed = FALSE;
    printf("Result: %d tests pass; %d tests fail \n", ntest-nFail, nFail);
    Lgm_free_ctrans( c ); // free the structure

    ck_assert_msg( Passed, "GSE test failed. GSE and GSE2000 do not appear equivalent at J2000.\n" );

    return;
    } END_TEST


START_TEST(test_CoordGSE_fail) {
    Lgm_CTrans        *c = Lgm_init_ctrans( 0 ); 
    Lgm_Vector        Ufrom, Utest1, Utest2, Udiff;
    int               nTests, nPass, nFail, ntest, Passed=FALSE;
    int               n=0, sysFrom[3]={2,6,9};
    long int          Date;
    double            UTC, del, tol=1e-9;

    Lgm_Set_CTrans_Options(LGM_EPH_DE, LGM_PN_IAU76, c);

    /* step through test cases one line at a time */
    ntest = 0;
    nPass = 0;
    nFail = 0;
    Passed = TRUE;

    Ufrom.x =   5.271324; Ufrom.y = 3.535221; Ufrom.z = 1.362777; // Set a vector in input coordinates

    printf("\nGSE/GSE2000 tests (non-equiv. away from J2000):\n");

    // **** Test **** //
    printf("Test %d. Non-equivalence of GSE and GSE2000 at arbitrary time (DE)\n", ntest+1);
    // Set up all the necessary variables to do transformations for this Date and UTC
    Lgm_Set_CTrans_Options(LGM_EPH_DE, LGM_PN_IAU76, c); /* Uses LGM high accuracy analytic solution */
    Lgm_Set_Coord_Transforms( 19940611, 1.0, c );

    // Do the transformations from input system to both target systems
    Lgm_Convert_Coords( &Ufrom, &Utest1, MOD_TO_GSE, c );
    Lgm_Convert_Coords( &Ufrom, &Utest2, MOD_TO_GSE2000, c );
    nFail = nFail + (testDiff(Utest1, Utest2, tol));
    ntest++;

    // **** Test 2 **** //
    printf("Test %d. Non-equivalence of GSE and GSE2000 at J2000 (HA)\n", ntest+1);
    // Set up all the necessary variables to do transformations for this Date and UTC
    Lgm_Set_CTrans_Options(LGM_EPH_HIGH_ACCURACY, LGM_PN_IAU76, c); /* Uses LGM high accuracy analytic solution */
    Lgm_Set_Coord_Transforms( 20061120, 15.0, c );

    // Do the transformation to intermediate system and back out
    Lgm_Convert_Coords( &Ufrom, &Utest1, MOD_TO_GSE, c );
    Lgm_Convert_Coords( &Ufrom, &Utest2, MOD_TO_GSE2000, c );
    nFail = nFail + (testDiff(Utest1, Utest2, tol));
    ntest++;


    // **** Tests complete *** //
    printf("Result: %d tests pass; %d tests fail \n", nFail, ntest-nFail);
    Lgm_free_ctrans( c ); // free the structure

    ck_assert( nFail>0 );

    return;
    } END_TEST


START_TEST(test_CoordDipoleTilt) {
    /* Dipole tilt*/
    Lgm_CTrans        *c = Lgm_init_ctrans( 0 );
    int               nTests, line, nPass, nFail, Passed=FALSE;
    double            TiltTest, TiltExpected, del;
    char              buff[262];
    char              IsoDate[80];
    FILE              *testfile, *outfile;
    Lgm_DateTime      d;

    int makeNew = 1;

    /* read test file */
    testfile = fopen("check_CoordDipoleTilt.expected","r");
    if (makeNew) outfile = fopen("check_CoordDipoleTilt.got", "w");

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
            sscanf(buff,"%s %lf", IsoDate, &TiltExpected);
            line++;

            // Set up all the necessary variables to do transformations for this Date and UTC
            IsoTimeStringToDateTime( IsoDate, &d, c );
            //printf("IsoDate = %s; d.Date, d.time = %ld, %lf \n", IsoDate, d.Date, d.Time);
	    TiltTest = Lgm_Dipole_Tilt(d.Date, d.Time);
            del = TiltTest - TiltExpected;
            nTests++;
            if (fabs(del) <= 1.0e-5) {
                nPass++;
                printf("Test %d passed\n", nTests);
                }
            else {
                nFail++;
                printf("*****  warning : difference >= 1.0e-5 r  *****\n");
                printf("Test %d failed (diff: %g)\n", nTests, del);
                }
            if (makeNew) fprintf(outfile, "%s %lf\n", IsoDate, TiltTest);
            }
        else {
            if (makeNew) fprintf(outfile, "%s", buff);
            }
        }
    if (nFail>0) Passed = FALSE;
    fclose(testfile);
    if (makeNew) fclose(outfile);
    printf("Result: %d tests pass; %d tests fail (Precision=1.0e-5)\n", nPass, nFail);
    fflush(stdout); // get all that status info out before ck_assert
    Lgm_free_ctrans( c ); // free the structure

    ck_assert_msg( Passed, "CoordDipoleTilt test failed.\n" );

    return;
    } END_TEST


int testDiff(Lgm_Vector Utest, Lgm_Vector Utarg, double tol) {
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

    if ((fabs(del) <= tol) && (fabs(Udiff.x) <= tol)  && (fabs(Udiff.y) <= tol)  && (fabs(Udiff.z) <= tol)){
        printf("Test passed\n\n");
        }
    else {
        Fail = TRUE;
        printf("Test failed (diff: %g %g %g   %g)\n\n", Udiff.x, Udiff.y, Udiff.z, fabs(del));
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
    else if (!strcmp(sys, "WGS84")) {
        out = WGS84_COORDS;
        }
    return(out);
    }


Suite *CT_suite(void) {

  Suite *s = suite_create("CoordTrans_TESTS");

  TCase *tc_CoordTrans = tcase_create("Coordinate Transformations");

  tcase_add_test(tc_CoordTrans, test_CoordTrans);
  tcase_add_test(tc_CoordTrans, test_CoordTransNoEph);
  tcase_add_test(tc_CoordTrans, test_CoordRoundtrip);
  tcase_add_test(tc_CoordTrans, test_CoordGSE_equiv);
  tcase_add_test(tc_CoordTrans, test_CoordGSE_fail);
  tcase_add_test(tc_CoordTrans, test_CoordDipoleTilt);

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
