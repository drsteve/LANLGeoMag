#include <Lgm_CTrans.h>
#include <Lgm_Vec.h>
#include <Lgm_DynamicMemory.h>
#include <Lgm_WGS84.h>

/* Example program to show how to set up coordinate transforms
 * and to compare the use of different methods of finding the
 * Sun direction.
 */


int main( ) {

    Lgm_CTrans        *c = Lgm_init_ctrans( 0 ); 
    Lgm_Vector        EcPole_J2000;
    Lgm_DateTime      DateTime;
    long int          Date;
    double            UTC, ang_deg;
    
    /* Set date and time for LGM */
    Date = 20000101;              // Jan 1, 2017
    UTC  = 12.0+0.0/60.0;         // Universal Time Coordinated (in decimal hours)
    printf("Date and Time for example calculations:\n");
    Lgm_Make_UTC( 20040514, 16.0 + 43.0/60.0, &DateTime, c );
    Lgm_Print_DateTime( &DateTime, 4, 8 ); printf("\n\n\n");

    /* Set up all the necessary variables to do transformations for this Date and UTC
     * The CTrans options let you select the use of JPL Development Ephmerides, a
     * high-accuracy algortihm, or a lower accuracy algorithm, for finding the Sun
     * direction.
     */
    Lgm_Set_CTrans_Options( LGM_EPH_DE, LGM_PN_IAU76, c );
    Lgm_Set_Coord_Transforms( Date, UTC, c );
    printf("DE421 Earth-Sun Vector in MOD: ");
    Lgm_PrintVector(&c->Sun);
    
    //z_vec of GSE in MOD is c->EcPole;
    ang_deg = Lgm_VectorAngle(&c->EcPole, &c->Sun);

    printf("Angle between DE421 Sun vector and Ecliptic Pole = %g\n", ang_deg);
    printf("Difference from 90deg = %g\n", ang_deg-90.0);

    printf("\n");
    printf("Transformation Matrix:\n");
    printf("                        [ %15.8lf  %15.8lf  %15.8lf ]\n",   c->Amod_to_gse[0][0], c->Amod_to_gse[1][0], c->Amod_to_gse[2][0]);
    printf("  Amod_to_gse (DE) =    [ %15.8lf  %15.8lf  %15.8lf ]\n",   c->Amod_to_gse[0][1], c->Amod_to_gse[1][1], c->Amod_to_gse[2][1]);
    printf("                        [ %15.8lf  %15.8lf  %15.8lf ]\n\n", c->Amod_to_gse[0][2], c->Amod_to_gse[1][2], c->Amod_to_gse[2][2]);


    /* Redo the setup for the high-accuracy Sun direction alogrithm */
    Lgm_Set_CTrans_Options( LGM_EPH_HIGH_ACCURACY, LGM_PN_IAU76, c );
    Lgm_Set_Coord_Transforms( Date, UTC, c );
    ang_deg = Lgm_VectorAngle(&c->EcPole, &c->Sun);
    printf("\nHigh-accuracy Earth-Sun Vector in MOD: ");
    Lgm_PrintVector(&c->Sun);

    printf("Angle between LGM high-Q Sun vector and Ecliptic Pole = %g\n", ang_deg);
    printf("Difference from 90deg = %g\n", ang_deg-90.0);
    
    printf("                        [ %15.8lf  %15.8lf  %15.8lf ]\n",   c->Amod_to_gse[0][0], c->Amod_to_gse[1][0], c->Amod_to_gse[2][0]);
    printf("  Amod_to_gse (HA)   =  [ %15.8lf  %15.8lf  %15.8lf ]\n",   c->Amod_to_gse[0][1], c->Amod_to_gse[1][1], c->Amod_to_gse[2][1]);
    printf("                        [ %15.8lf  %15.8lf  %15.8lf ]\n\n", c->Amod_to_gse[0][2], c->Amod_to_gse[1][2], c->Amod_to_gse[2][2]);
 

    /* Redo the setup for the lower-accuracy Sun direction alogrithm */
    Lgm_Set_CTrans_Options( LGM_EPH_LOW_ACCURACY, LGM_PN_IAU76, c );
    Lgm_Set_Coord_Transforms( Date, UTC, c );
    ang_deg = Lgm_VectorAngle(&c->EcPole, &c->Sun);
    printf("\nLower-accuracy Earth-Sun Vector in MOD: ");
    Lgm_PrintVector(&c->Sun);

    printf("Angle between LGM low-Q Sun vector and Z-axis = %g\n", ang_deg);
    printf("Difference from 90deg = %g\n", ang_deg-90.0);

    printf("                        [ %15.8lf  %15.8lf  %15.8lf ]\n",   c->Amod_to_gse[0][0], c->Amod_to_gse[1][0], c->Amod_to_gse[2][0]);
    printf("  Amod_to_gse (LA)   =  [ %15.8lf  %15.8lf  %15.8lf ]\n",   c->Amod_to_gse[0][1], c->Amod_to_gse[1][1], c->Amod_to_gse[2][1]);
    printf("                        [ %15.8lf  %15.8lf  %15.8lf ]\n\n", c->Amod_to_gse[0][2], c->Amod_to_gse[1][2], c->Amod_to_gse[2][2]);



    Lgm_Convert_Coords( &c->EcPole, &EcPole_J2000, MOD_TO_GEI2000, c );
    Lgm_NormalizeVector(&EcPole_J2000);
    printf("Ecliptic Pole Vector at J2000: %15.8lf, %15.8lf, %15.8lf \n", EcPole_J2000.x, EcPole_J2000.y, EcPole_J2000.z);

    Lgm_free_ctrans( c ); // free the structure

    return(0);
    }


