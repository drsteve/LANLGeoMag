#include <Lgm_CTrans.h>
#include <Lgm_Vec.h>
#include <Lgm_JPLeph.h>

int main( ) {

    Lgm_CTrans        *c = Lgm_init_ctrans( 0 ); 
    Lgm_JPLephemInfo  *jpl = Lgm_InitJPLephemInfo( 421, LGM_DE_SUN|LGM_DE_EARTHMOON|LGM_DE_OUTERPLANETS, 1);
    Lgm_Vector        Uicrf, pyResult1, pyResult2, pyResult3;
    long int          Date;
    double            UTC, JD;
    
    Date = 20170101;              // Jan 1, 2017
    UTC  = 12.0+0.0/60.0;         // Universal Time Coordinated (in decimal hours)
    pyResult1.x =  506274.15309483; //Sun (wrt solar sys barycenter)
    pyResult1.y =  550063.81126863; 
    pyResult1.z =  213030.58611033;
    pyResult2.x =  -27645517.94343229; //Earth-Moon barycenter (wrt solar sys barycenter)
    pyResult2.y =  133019570.59230532;
    pyResult2.z =  57639804.30457275;
    pyResult3.x =  1445862314.051871; //Pluto
    pyResult3.y =  -4400572504.816436;
    pyResult3.z =  -1808918123.981852;

    // Set up all the necessary variables to do transformations for this Date and UTC
    Lgm_Set_Coord_Transforms( Date, UTC, c );
    JD = Lgm_Date_to_JD( Date, UTC, c);
    Lgm_ReadJPLephem( jpl );
    Lgm_JPLephem_position( JD, LGM_DE_SUN, jpl, &Uicrf);

    // Test a few planets for this time
    printf("DE %d\n",jpl->DEnum);
    printf("Date, UTC = %8ld %lf\n", Date, UTC);
    printf("JD  = %lf\n", JD);
    printf("Sun Position (ICRF): ");
    Lgm_PrintVector(&Uicrf);
    printf("Difference from expected result = %g %g %g\n\n", Uicrf.x-pyResult1.x, Uicrf.y-pyResult1.y, Uicrf.z-pyResult1.z);

    Lgm_JPLephem_position( JD, LGM_DE_EARTHMOON, jpl, &Uicrf);
    printf("Earth-Moon Barycenter Position (ICRF): ");
    Lgm_PrintVector(&Uicrf);
    printf("Difference from expected result = %g %g %g\n\n", Uicrf.x-pyResult2.x, Uicrf.y-pyResult2.y, Uicrf.z-pyResult2.z);

    Lgm_JPLephem_position( JD, LGM_DE_PLUTO, jpl, &Uicrf);
    printf("Pluto Position (ICRF): ");
    Lgm_PrintVector(&Uicrf);
    printf("Difference from expected result = %g %g %g\n\n", Uicrf.x-pyResult3.x, Uicrf.y-pyResult3.y, Uicrf.z-pyResult3.z);

    //Test Sun position for different time
    Date = 19821111;              // Jan 1, 2017
    UTC  = 11.0+0.0/60.0;         // Universal Time Coordinated (in decimal hours)
    pyResult1.x =  922949.0261094079; //Sun (wrt solar sys barycenter)
    pyResult1.y = 1047536.3872436787;
    pyResult1.z =  414153.2146306878;

    Lgm_Set_Coord_Transforms( Date, UTC, c );
    JD = Lgm_Date_to_JD( Date, UTC, c);
    Lgm_JPLephem_position( JD, LGM_DE_SUN, jpl, &Uicrf);

    printf("Date, UTC = %8ld %lf\n", Date, UTC);
    printf("JD  = %lf\n", JD);
    printf("Sun Position (ICRF): ");
    Lgm_PrintVector(&Uicrf);
    printf("Difference from expected result = %g %g %g\n\n", Uicrf.x-pyResult1.x, Uicrf.y-pyResult1.y, Uicrf.z-pyResult1.z);

    Lgm_free_ctrans( c ); // free the structure
    Lgm_FreeJPLephemInfo( jpl ); // free the structure

    return(0);
    }


