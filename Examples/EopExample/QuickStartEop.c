#include <Lgm_CTrans.h>
#include <Lgm_Eop.h>
int main( ) {
    Lgm_CTrans      *c = Lgm_init_ctrans( 1 ); // more compact declaration
    Lgm_Eop         *e = Lgm_init_eop( 1 );
    Lgm_EopOne      eop;
    Lgm_Vector      Ugsm, Usm;
    long int        Date;
    double          UTC, JD;
    
    Date = 20040812;                  // August 12, 2004
    UTC  = 12.34567;                  // Universal Time Coordinated (in decimal hours)
Date = 20140829;    // Jan 1, 2000
UTC  = 15.0+32.0/60.0+13.811285/3600.0;

    JD = Lgm_Date_to_JD( Date, UTC, c ); // Compute JD
    Ugsm.x = -6.6; Ugsm.y = 3.4; Ugsm.z = -2.3; // Set a vector in GSM coordinates

    // Read in the EOP vals
    Lgm_read_eop( e );

    // Get (interpolate) the EOP vals from the values in the file at the given Julian Date
    Lgm_get_eop_at_JD( JD, &eop, e );

    // Set the EOP vals in the CTrans structure.
    Lgm_set_eop( &eop, c );


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
    Lgm_destroy_eop( e );
    return(0);

}

