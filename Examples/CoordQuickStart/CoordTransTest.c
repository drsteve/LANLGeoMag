#include <stdio.h>
#include <string.h>
#include <Lgm_CTrans.h>
#include <Lgm_Vec.h>

int getSys( char *sys );


int main( ) {

    Lgm_CTrans        *c = Lgm_init_ctrans( 0 ); 
    Lgm_Vector        Utarg, Ucent, Utest, Udiff;
    int               nTests, line, nPass, nFail, transflag;
    double            del;
    char              buff[262], sysIn[10], sysOut[10];
    char              IsoDate[80];
    long long         TT2000;
    FILE              *testfile, *outfile;
    Lgm_DateTime      d;

    int makeNew = 0;
    Lgm_Set_CTrans_Options(LGM_EPH_DE, LGM_PN_IAU76, c);

    /* read test file */
    testfile = fopen("LgmTestTransforms_MMS.txt","r");
    if (makeNew) outfile = fopen("LgmTestTransforms_MMS_new.txt", "w");

    /* step through test cases one line at a time */
    line = 0;
    nTests = 0;
    nPass = 0;
    nFail = 0;
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
    fclose(testfile);
    if (makeNew) fclose(outfile);
    printf("Result: %d tests pass; %d tests fail (Precision=1.0e-5 km [=1cm])\n", nPass, nFail);
    Lgm_free_ctrans( c ); // free the structure

    return(0);
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
