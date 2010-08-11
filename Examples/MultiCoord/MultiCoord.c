#include <Lgm_CTrans.h>
int main( ) {
    Lgm_CTrans   *c = Lgm_init_ctrans( 1 ); // more compact declaration
    long int    Date;
    double      UTC, R, MLAT, MLT, MLON;
    Lgm_Vector  Ugsm, v, w;
    
    Date = 20040812;    // August 12, 2004
    UTC  = 12.0+20.0/60.0+44.0/3600.0;  // Universal Time Coordinated (in decimal hours)
    Ugsm.x = -6.6; Ugsm.y =  6.6; Ugsm.z = 6.6; // Set a vector in GSM coordinates:wq


    // Set up all the necessary variables to do transformations for this
    // Date and UTC
    Lgm_Set_Coord_Transforms( Date, UTC, c);


    // Print out the final results
    printf("Date = %8ld\n", Date);
    printf("UTC  = %lf\n", UTC);
    printf("Ugsm = %.8lf %.8lf %.8lf\n", Ugsm.x, Ugsm.y, Ugsm.z);

    Lgm_Convert_Coords( &Ugsm, &v, GSM_TO_SM, c );
    printf("Usm  = %.8lf %.8lf %.8lf\n", v.x, v.y, v.z);

    Lgm_Convert_Coords( &Ugsm, &v, GSM_TO_GSE, c );
    printf("Ugse  = %.8lf %.8lf %.8lf\n", v.x, v.y, v.z);
    
    Lgm_Convert_Coords( &Ugsm, &v, GSM_TO_WGS84, c );
    printf("Uwgs84  = %.8lf %.8lf %.8lf\n", v.x, v.y, v.z);
    
    Lgm_Convert_Coords( &Ugsm, &v, GSM_TO_GEI2000, c );
    printf("Ugei2000  = %.8lf %.8lf %.8lf\n", v.x, v.y, v.z);

    Lgm_Convert_Coords( &Ugsm, &v, GSM_TO_MOD, c );
    printf("Umod  = %.8lf %.8lf %.8lf\n", v.x, v.y, v.z);

    Lgm_Convert_Coords( &Ugsm, &v, GSM_TO_MOD, c );
    printf("Umod  = %.8lf %.8lf %.8lf\n", v.x, v.y, v.z);
    
    Lgm_Convert_Coords( &Ugsm, &v, GSM_TO_CDMAG, c );
    printf("Ucdmag  = %.8lf %.8lf %.8lf\n", v.x, v.y, v.z);

    Lgm_Convert_Coords( &Ugsm, &v, GSM_TO_EDMAG, c );
    printf("Uedmag  = %.8lf %.8lf %.8lf\n", v.x, v.y, v.z);

    Lgm_Convert_Coords( &v, &w, EDMAG_TO_GSM, c );
    printf("Ugsm  = %.8lf %.8lf %.8lf (testing ED back to GSM)\n", w.x, w.y, w.z);





    Lgm_CDMAG_to_R_MLAT_MLON_MLT( &v, &R, &MLAT, &MLON, &MLT, c );
    printf("Cdmag R, MLAT, MLON, MLT  = %.8lf %.8lf %.8lf %.8lf\n", R, MLAT, MLON, MLT);

    double GLAT = 35.00, GLON = -106.0;
    printf("\nGiven: GLAT, GLON = %g %g\n", GLAT, GLON);
    Lgm_GLATLON_TO_CDMLATLONMLT( GLAT, GLON, &MLAT, &MLON, &MLT, c );
    printf("Cdmag MLAT, MLON, MLT  = %.8lf %.8lf %.8lf\n", MLAT, MLON, MLT);

    Lgm_GLATLON_TO_EDMLATLONMLT( GLAT, GLON, &MLAT, &MLON, &MLT, c );
    printf("Edmag MLAT, MLON, MLT  = %.8lf %.8lf %.8lf\n", MLAT, MLON, MLT);
    

    
    
    
    
    Lgm_free_ctrans( c ); // free the structure
    return(0);
}

