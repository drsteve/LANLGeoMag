#include <check.h>
#include "../libLanlGeoMag/Lgm/Lgm_CTrans.h"
#include "../libLanlGeoMag/Lgm/Lgm_Vec.h"

#define TRUE    1
#define FALSE   0


int target( int targ );


START_TEST(test_DE421) {

    Lgm_CTrans        *c = Lgm_init_ctrans( 0 ); 
    Lgm_JPLephemInfo  *jpl = Lgm_InitJPLephemInfo( 421, LGM_DE_SUN|LGM_DE_EARTHMOON|LGM_DE_INNERPLANETS|LGM_DE_OUTERPLANETS, 1);
    Lgm_Vector        Utarg, Ucent, Utest;
    long int          Date;
    int               targObj, ntarg, nctr, ncoord, nTests, line, Year, Month, Day, nPass, nFail, Passed=FALSE;
    double            xi, UTC, JD, del, rval;
    char              buff[102];
    FILE              *testfile;

    Lgm_ReadJPLephem( jpl );

    /* read test file */
    testfile = fopen("testpo.421","r");

    while(1) {
     fgets(buff,100,testfile);
     buff[3]='\0';
     if(strcmp(buff,"EOT")==0) break;
    }

    /* step through test cases one line at a time */
    printf("DE %d\n",jpl->DEnum);
    line = 0;
    nTests = 0;
    nPass = 0;
    nFail = 0;
    Passed = TRUE;
    while( fgets(buff,100,testfile) != NULL) {
        // read line
        sscanf(buff+15," %lf %d %d %d %lf",&JD,&ntarg,&nctr,&ncoord,&xi);
        line++;

        if (ntarg<=13) { //exclude nutations/librations for now
            // Set up all the necessary variables to do transformations for this Date and UTC
            Date = Lgm_JD_to_Date( JD, &Year, &Month, &Day, &UTC );
            Lgm_Set_Coord_Transforms( Date, UTC, c );
            targObj = target(ntarg);
            if (targObj > 0) {
                if (ncoord<=3) {
                    Lgm_JPLephem_position( JD, targObj, jpl, &Utarg);
                    }
                else {
                    Lgm_JPLephem_velocity( JD, targObj, jpl, &Utarg);
                    }
                }
            else {
                Utarg.x = 0.0;
                Utarg.y = 0.0;
                Utarg.z = 0.0;
                }
        } else {
            continue;
            }
        if (nctr<=13) {
            targObj = target(nctr);
            if (targObj > 0) {
                if (ncoord<=3) {
                    Lgm_JPLephem_position( JD, targObj, jpl, &Ucent);
                    }
                else {
                    Lgm_JPLephem_velocity( JD, targObj, jpl, &Ucent);
                    }
                }
            else {
                Ucent.x = 0.0;
                Ucent.y = 0.0;
                Ucent.z = 0.0;
                }
        } else {
            continue;
            }
        //get target relative to center
        Utest.x = Utarg.x - Ucent.x;
        Utest.y = Utarg.y - Ucent.y;
        Utest.z = Utarg.z - Ucent.z;
        //then test difference
        switch (ncoord) {
            case 1:
                rval = Utest.x;
                break;
            case 2:
                rval = Utest.y;
                break;
            case 3:
                rval = Utest.z;
                break;
            case 4:
                rval = Utest.x;
                break;
            case 5:
                rval = Utest.y;
                break;
            case 6:
                rval = Utest.z;
                break;
            }

        xi  = xi*LGM_DE421_AU; //pos are AU, vel are AU/day (??)
        del = rval - xi;
        nTests++;
        if (fabs(del) <= 2.0e-5) {
            nPass++;
            printf("%7d %7d %10.1f %2d %2d %2d %25.13f %25.13f %13.5e\n",
                   nTests, line, JD, ntarg, nctr, ncoord, xi, rval, del);
            }
        else {
            nFail++;
            printf("*****  warning : difference >= 1.0e-5 km (1 cm)  *****\n");
            printf("%7d %7d %10.1f %2d %2d %2d %25.13f %25.13f %13.5e\n",
                   nTests, line, JD, ntarg, nctr, ncoord, xi, rval, del);
            }

        }
    if (nFail>0) Passed = FALSE;
    fclose(testfile);
    printf("Result: %d tests pass; %d tests fail (Precision=1.0e-5 km [=1cm])\n", nPass, nFail);
    printf("Caveat: No tests on libration or nutation\n\n");
    Lgm_free_ctrans( c ); // free the structure
    Lgm_FreeJPLephemInfo( jpl ); // free the structure

    fail_unless( Passed, "DE421 test failed. LANLGeoMag does not match the test file from JPL.\n" );
   
    return;
    } END_TEST

int target( int targ ) {
    int out;
    switch (targ) {
        case 1:
            out = LGM_DE_MERCURY;
            break;
        case 2:
            out = LGM_DE_VENUS;
            break;
        case 3:
            out = LGM_DE_EARTH;
            break;
        case 4:
            out = LGM_DE_MARS;
            break;
        case 5:
            out = LGM_DE_JUPITER;
            break;
        case 6:
            out = LGM_DE_SATURN;
            break;
        case 7:
            out = LGM_DE_URANUS;
            break;
        case 8:
            out = LGM_DE_NEPTUNE;
            break;
        case 9:
            out = LGM_DE_PLUTO;
            break;
        case 10:
            out = LGM_DE_MOON_ICRF;
            break;
        case 11:
            out = LGM_DE_SUN;
            break;
        case 12:
            out = -1;
            break;
        case 13:
            out = LGM_DE_EARTHMOON;
            break;
        }
    return(out);
    }



Suite *DE421_suite(void) {

  Suite *s = suite_create("DE421_TESTS");

  TCase *tc_DE421 = tcase_create("DE421");

  tcase_add_test(tc_DE421, test_DE421);

  suite_add_tcase(s, tc_DE421);

  return s;

}

int main(void) {

    int      number_failed;
    Suite   *s  = DE421_suite();
    SRunner *sr = srunner_create(s);

    printf("\n\n");
    srunner_run_all(sr, CK_NORMAL);
    number_failed = srunner_ntests_failed(sr);
    srunner_free(sr);

    return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;

}
