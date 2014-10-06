#include <Lgm_CTrans.h>

#include <Lgm_Vec.h>

int main( ) {

    Lgm_CTrans  *c = Lgm_init_ctrans( 1 ); 
    Lgm_Vector  Ugsm, Usm, Ugse1, Ugse2, tmp, Uj2000;
    long int    Date;
    double      UTC, Lat, Lon, r;
    
    Date = 20170101;    // Jan 1, 2000
    UTC  = 12.0+0.0/60.0;         // Universal Time Coordinated (in decimal hours)
    Ugsm.x = -6.6; Ugsm.y = 3.4; Ugsm.z = -2.3; // Set a vector in GSM coordinates

    // Set up all the necessary variables to do transformations for this Date and UTC
    Lgm_Set_Coord_Transforms( Date, UTC, c );

    // Do the transformation from GSM->SM
    Lgm_Convert_Coords( &Ugsm, &Usm, GSM_TO_SM, c );
    // Do the transformation from GSM->GSE
    Lgm_Convert_Coords( &Ugsm, &Ugse1, GSM_TO_GSE, c );
    // Do the transformation from SM->GSE
    Lgm_Convert_Coords( &Usm, &Ugse2, SM_TO_GSE, c );


    // Print out the final results
    printf("Date = %8ld\n", Date);
    printf("UTC  = %lf\n", UTC);
    printf("Ugsm = %.8lf %.8lf %.8lf Re\n", Ugsm.x, Ugsm.y, Ugsm.z);
    printf("Usm  = %.8lf %.8lf %.8lf Re\n", Usm.x, Usm.y, Usm.z);
    printf("\nGoing to GSE from SM and GSM\n");
    printf("Ugse  = %.8lf %.8lf %.8lf Re\n", Ugse1.x, Ugse1.y, Ugse1.z);
    printf("Ugse  = %.8lf %.8lf %.8lf Re\n", Ugse2.x, Ugse2.y, Ugse2.z);
    printf("They are different by\n");
    Lgm_VecSub(&tmp, &Ugse1, &Ugse2 );

    Lgm_PrintVector(&tmp);
    printf("\n");


    // compute the ground track
    // Do the transformation from SM->J2000
    Lgm_Convert_Coords( &Usm, &Uj2000, SM_TO_GEI2000, c );
    printf("The ground track point:\n");
    Lgm_CartToSphCoords(&Uj2000, &Lat, &Lon, &r);
    printf("Lat:%lf Lon:%lf\n", Lat, Lon);


    Lgm_free_ctrans( c ); // free the structure

    return(0);
}

