#include <Lgm_CTrans.h>
#include <Lgm_Vec.h>
#include <Lgm_JPLeph.h>
#include <Lgm_WGS84.h>

int main( ) {

    Lgm_CTrans        *c = Lgm_init_ctrans( 0 ); 
    Lgm_JPLephemInfo  *jpl = Lgm_InitJPLephemInfo( 421, LGM_DE_SUN|LGM_DE_EARTHMOON, 1);
    Lgm_Vector        Uicrf, Umod, Ulgm, diff;
    long int          Date;
    double            UTC, JD, Dec, RA, r;
    
    Date = 20170101;              // Jan 1, 2017
    UTC  = 12.0+0.0/60.0;         // Universal Time Coordinated (in decimal hours)

    // Set up all the necessary variables to do transformations for this Date and UTC
    Lgm_Set_Coord_Transforms( Date, UTC, c );
    JD = Lgm_Date_to_JD( Date, UTC, c);
    Lgm_ReadJPLephem( jpl );
    Lgm_JPL_getSunVector( JD, jpl, &Uicrf);
    printf("DE421 Earth-Sun Vector: ");
    Lgm_PrintVector(&Uicrf);
    //Lgm_NormalizeVector(&Uicrf);
    //printf("DE421 E-S Vec (Normed): ");
    //Lgm_PrintVector(&Uicrf);
    
    Ulgm.x = c->earth_sun_dist*WGS84_A * c->Sun.x ;
    Ulgm.y = c->earth_sun_dist*WGS84_A * c->Sun.y ;
    Ulgm.z = c->earth_sun_dist*WGS84_A * c->Sun.z ;
    printf("Lgm GEI2000 Sun Vector: ");
    Lgm_PrintVector(&Ulgm);

    diff.x = Uicrf.x - Ulgm.x;
    diff.y = Uicrf.y - Ulgm.y;
    diff.z = Uicrf.z - Ulgm.z;
    printf("Difference DE421-Lgm  : ");
    Lgm_PrintVector(&diff);

    printf("\nRA-Dec differences in MOD\n");
    Lgm_Convert_Coords( &Uicrf, &Umod, GEI2000_TO_MOD, c );
    Lgm_CartToSphCoords( &Umod, &Dec, &RA, &r);
    RA += 360.0;
    printf("DE421 Sun Vec (RA, Dec): %g %g\n", RA, Dec);
    printf("Lgm   Sun Vec (RA, Dec): %g %g\n", c->RA_sun, c->DEC_sun);

    Lgm_free_ctrans( c ); // free the structure
    Lgm_FreeJPLephemInfo( jpl ); // free the structure

    return(0);
    }


