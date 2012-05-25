#include <Lgm_CTrans.h>
int main( ) {
    Lgm_CTrans      *c = Lgm_init_ctrans( 1 ); // more compact declaration
    Lgm_Vector      Ugsm, Usm;
    long int        Date;
    double      UTC;
    
    Date = 20040812;    // August 12, 2004

    UTC  = 12.34567;  // Universal Time Coordinated (in decimal hours)
    Ugsm.x = -6.6; Ugsm.y = 3.4; Ugsm.z = -2.3; // Set a vector in GSM coordinates

    // Set up all the necessary variables to do transformations for this
    // Date and UTC
    Lgm_Set_Coord_Transforms( Date, UTC, c );

    // Do the transformation from GSM->SM
    Lgm_Convert_Coords( &Ugsm, &Usm, GSM_TO_SM, c );

    // Print out the final results
    printf("Date = %8ld\n", Date);
    printf("UTC  = %lf\n", UTC);
    printf("Ugsm = %.8lf %.8lf %.8lf Re\n", Ugsm.x, Ugsm.y, Ugsm.z);
    printf("Usm  = %.8lf %.8lf %.8lf Re\n", Usm.x, Usm.y, Usm.z);



    Lgm_free_ctrans( c ); // free the structure

    return(0);
}

