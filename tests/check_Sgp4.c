#include <check.h>
#include "../libLanlGeoMag/Lgm/Lgm_Sgp.h"
#include <stdio.h>
#include <stdlib.h>
#define TRUE    1
#define FALSE   0

/*
 *  Unit and Regression tests for Sgp4 solvers
 */


START_TEST(test_Sgp4_01) {

    double          tsince, Xe, Ye, Ze, VXe, VYe, VZe, X, Y, Z, VX, VY, VZ;
    char            Line[5000], Line0[120], Line1[120], Line2[120];
    int             nTLEs, SatNum, Passed=FALSE;
    FILE            *fp_expected;
    FILE            *fp_got;
    _SgpTLE         *TLEs;  // pointer to a struct
    _SgpInfo        *s;

    // create a place to hold an array of TLEs (1 element)
    TLEs  =  (_SgpTLE *)calloc( 100, sizeof(_SgpTLE) );
    s     = (_SgpInfo *)calloc( 1, sizeof(_SgpInfo) );
    nTLEs = 0;

    Passed = TRUE;
    if ( (fp_expected = fopen( "check_Sgp4_01.expected", "r" )) != NULL ) {

        while ( fgets( Line, 4096, fp_expected ) != NULL ) {

            if ( Line[0] == '#' ) {
                // a comment line - do nothing.
            } else if ( strstr( Line, "xx" ) != NULL ) {
                sscanf( Line, "%d xx", &SatNum );
                printf("New Sat: %d\n", SatNum);


                // read in the TLE
                fgets( Line1, 100, fp_expected );
                fgets( Line2, 100, fp_expected );
                sprintf( Line0, "%d xx", SatNum );


                // Add the three lines of the TLE to the list of TLEs
                LgmSgp_ReadTlesFromStrings( Line0, Line1, Line2, &nTLEs, TLEs, 0 );

                // initialize the SGP4 propagator
                LgmSgp_SGP4_Init( s, &TLEs[nTLEs-1] );
                
            } else {
                // Its a line of expected results ( Time since Epoch in minutes, X, Y, Z, VX, VY, VZ )
                sscanf( Line, "%lf %lf %lf %lf %lf %lf %lf", &tsince, &Xe, &Ye, &Ze, &VXe, &VYe, &VZe );

                // try to duplicate -- do the propagation
                LgmSgp_SGP4( tsince, s ); X = s->X; Y = s->Y; Z = s->Z; VX = s->VX; VY = s->VY; VZ = s->VZ;

                if ( (SatNum == 33334) ) {
                    if (s->error != 3)  Passed = FALSE;
                } else {
                    if ( ( fabs(Xe-X) > 1e-7 )   || ( fabs(Ye-Y) > 1e-7 )  || ( fabs(Ze-Z) > 1e-7 )  
                            || ( fabs(VXe-VX) > 1e-7 ) || ( fabs(VYe-VY) > 1e-7 ) || ( fabs(VZe-VZ) > 1e-7 ) ) {

                        Passed = FALSE;
                        printf("\tTime: %-16.8lf\n", tsince );
                        printf("\t\t   Expected: %16.8lf %16.8lf %16.8lf     %12.9lf %12.9lf %12.9lf\n", Xe, Ye, Ze, VXe, VYe, VZe );  
                        printf("\t\t        Got: %16.8lf %16.8lf %16.8lf     %12.9lf %12.9lf %12.9lf\n\n", X, Y, Z, VX, VY, VZ );  
                    }
                }
                
            }
            
        }

        fclose( fp_expected );

    } else {
        printf("Cant open file: check_Sgp4_01.expected\n" );
    }


//    if (   (fabs( creal(z1)-z1_r_expected) < 1e-10) && (fabs( cimag(z1)-z1_i_expected) < 1e-10)
//        && (fabs( creal(z2)-z2_r_expected) < 1e-10) && (fabs( cimag(z2)-z2_i_expected) < 1e-10)
//        && (cabs(Residual1)< 1e-10) && (cabs(Residual2) < 1e-10)
//        && (nReal == nReal_expected) ) Passed = TRUE;
//
//    if ( !Passed ){
//        printf("\nTest 01,  Lgm_QuadraticRoots(): Expected %d real roots, got %d\n\n", nReal_expected, nReal );
//        printf("       Root 1, Expected: %.15lf + %.15lf I\n",   z1_r_expected, z1_i_expected );
//        printf("                    Got: %.15lf + %.15lf I\n", creal(z1),     cimag(z1) );
//        printf("               Residual: %.15g + %.15g    |Residual| = %.15g\n\n", creal(Residual1), cimag(Residual1), cabs(Residual1) );
//        printf("       Root 2, Expected: %.15lf + %.15lf I\n",   z2_r_expected, z2_i_expected );
//        printf("                    Got: %.15lf + %.15lf I\n\n", creal(z2),     cimag(z2) );
//        printf("               Residual: %.15g + %.15g    |Residual| = %.15g\n\n\n", creal(Residual2), cimag(Residual2), cabs(Residual2) );
//    }

//    if ( (fp_got = fopen( "check_Sgp4_01.got", "w" )) != NULL ) {
//        fprintf( fp_got, "%d\n", nReal );
//        fprintf( fp_got, "%.15lf %.15lf\n", creal(z1), cimag(z1) );
//        fprintf( fp_got, "%.15lf %.15lf\n", creal(z2), cimag(z2) );
//        fclose( fp_got );
//    } else {
//        printf("Could not open file: check_Sgp4_01.got\\n");
//    }



    fail_unless( Passed, "Lgm_Sgp(): Regression test failed. Compare 'got' and 'expected' files: check_Sgp4_01.got check_Sgp4_01.expected\n" );


    free( s );
    free( TLEs );
    return;
}
END_TEST







Suite *Sgp4_suite(void) {

  Suite *s = suite_create("SGP4_TESTS");

  TCase *tc_Sgp4 = tcase_create("Sgp4");
  //tcase_add_checked_fixture(tc_Sgp4, Sgp4_Setup, Sgp4_TearDown);

  tcase_add_test(tc_Sgp4, test_Sgp4_01);

  suite_add_tcase(s, tc_Sgp4);

  return s;

}

int main(void) {

    int      number_failed;
    Suite   *s  = Sgp4_suite();
    SRunner *sr = srunner_create(s);

    printf("\n\n");
    srunner_run_all(sr, CK_NORMAL);
    number_failed = srunner_ntests_failed(sr);
    srunner_free(sr);

    return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;

}
