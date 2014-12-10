#include <Lgm_CTrans.h>
#include <Lgm_Vec.h>
#include <Lgm_JPLeph.h>

int main( ) {

    Lgm_CTrans  *c = Lgm_init_ctrans( 0 ); 
    Lgm_JPLeph  *jpl = Lgm_InitJPLeph( 421, LGM_DE_SUN, 1);
    Lgm_Vector  Uicrf;
    long int    Date;
    double      UTC, JD;
    
    Date = 20170101;    // Jan 1, 2000
    UTC  = 12.0+0.0/60.0;         // Universal Time Coordinated (in decimal hours)


    // Set up all the necessary variables to do transformations for this Date and UTC
    Lgm_Set_Coord_Transforms( Date, UTC, c );
    JD = Lgm_Date_to_JD( Date, UTC, c);
    Lgm_ReadJPLephem( jpl );
    Lgm_JPLephem_position( JD, LGM_DE_SUN, jpl, &Uicrf);

    printf("DE %d\n",jpl->DEnum);
    printf("Date = %8ld\n", Date);
    printf("UTC  = %lf\n", UTC);
    printf("JD  = %lf\n", JD);
    Lgm_PrintVector(&Uicrf);
    printf("\n");
    printf("(Should read [ 506274.15309483, 550063.81126863, 213030.58611033])\n");

    Lgm_free_ctrans( c ); // free the structure
    Lgm_FreeJPLeph( jpl ); // free the structure

    return(0);
    }


