#include <Lgm_CTrans.h>
#include <Lgm_Vec.h>
#include <Lgm_JPLeph.h>

int main( ) {

    Lgm_CTrans        *c = Lgm_init_ctrans( 0 ); 
    Lgm_JPLephemInfo  *jpl = Lgm_InitJPLephemInfo( 421, LGM_DE_SUN|LGM_DE_EARTHMOON, 1);
    Lgm_Vector        Uicrf, pyResult1, pyResult2;
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

    // Set up all the necessary variables to do transformations for this Date and UTC
    Lgm_Set_Coord_Transforms( Date, UTC, c );
    JD = Lgm_Date_to_JD( Date, UTC, c);
    Lgm_ReadJPLephem( jpl );
    Lgm_JPLephem_position( JD, LGM_DE_SUN, jpl, &Uicrf);

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

    Lgm_free_ctrans( c ); // free the structure
    Lgm_FreeJPLephemInfo( jpl ); // free the structure

    return(0);
    }


