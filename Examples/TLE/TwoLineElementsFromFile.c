#include <stdlib.h>
#include <stdio.h>
#include <Lgm_CTrans.h>
#include <Lgm_Sgp.h>

/*
 * This example shows how to compute positions of objects from TLEs (Two Line
 * Elements). In this example, we read multiple TLEs from a file
 */
int main( int argc, char *argv[] ){

    int         Year, Month, Day;
    int         StartYear, StartMonth, StartDay, StartDoy;
    int         EndYear, EndMonth, EndDay, EndDoy;
    long int    Date, StartDate, EndDate;
    double      UT, tsince, JD, StartUT, EndUT, StartJD, EndJD, JDinc;
    Lgm_CTrans  *c = Lgm_init_ctrans( 0 );
    Lgm_Vector  Ugse, Uteme;
    int         nTLEs, i;
    char        TleFile[512];
    char        *InputFile  = "input.txt";
    char        *OutputFile = "output.txt";
    FILE        *fp;


    /* 
     * Open input file and extract:
     *   1) name of a file containing TLEs
     *   2) Start Date and Time
     *   3) End Date and Time
     */
    if ( (fp = fopen( InputFile, "r" )) != NULL ) {
        fscanf( fp, "%s", TleFile );
        fscanf( fp, "%ld %lf", &StartDate, &StartUT );
        fscanf( fp, "%ld %lf", &EndDate, &EndUT );
    } else {
        printf( "Couldnt open file %s for reading\n", InputFile );
        exit( 1 );
    }
    fclose( fp );


    /*
     * Alloc some memory for the SgpInfo structure and the TLEs array (here we
     * alloc 20000 elements -- in more robust applications you'd want to
     * dynamically alloc and realloc as you go so you dont waste memory. )
     */
    _SgpInfo *s   = (_SgpInfo *)calloc( 1, sizeof(_SgpInfo) );
    _SgpTLE *TLEs = (_SgpTLE *)calloc( 20000, sizeof(_SgpTLE) );

    
    /*
     * Read in TLEs.
     */
    nTLEs = 0;
    LgmSgp_ReadTlesFromFile( TleFile, &nTLEs, TLEs, 1 );
    printf("Found %d Valid TLEs in file %s\n", nTLEs, TleFile );


    /*
     * All the TLEs have their own epoch times in them. And the propagator
     * (sgp4) uses the "time since (in minutes)". So for a given time of
     * interest, we need to compute the tsince needed. Convert Start/End Dates
     * to Julian dates -- they are easier to loop over contiguously.
     */
    Lgm_Doy( StartDate, &StartYear, &StartMonth, &StartDay, &StartDoy );
    StartJD = Lgm_JD( StartYear, StartMonth, StartDay, StartUT, LGM_TIME_SYS_UTC, c );

    Lgm_Doy( EndDate, &EndYear, &EndMonth, &EndDay, &EndDoy );
    EndJD = Lgm_JD( EndYear, EndMonth, EndDay, EndUT, LGM_TIME_SYS_UTC, c );

    JDinc = 5.0/1440.0; // 5-min increment

    if ( (fp = fopen( OutputFile, "w" )) == NULL ) {
        printf( "Couldnt open file %s for writing\n", OutputFile );
        exit( 1 );
    }

    // loop through each TLE found
    for (i=0; i<nTLEs; i++){

        // init SGP4
        LgmSgp_SGP4_Init( s, &TLEs[i] );
        fprintf(fp, "%%%s\n", TLEs[i].Line0 );
        fprintf(fp, "%%%s\n", TLEs[i].Line1 );
        fprintf(fp, "%%%s\n", TLEs[i].Line2 );

        // loop over specified time range
        for ( JD = StartJD; JD <= EndJD; JD += JDinc ) {

            // Convert the current JD back to Date/UT etc..
            Lgm_jd_to_ymdh ( JD, &Date, &Year, &Month, &Day, &UT );

            // Set up the trans matrices
            Lgm_Set_Coord_Transforms( Date, UT, c );
        
            // "time since" in minutes (thats what SGP4 wants)
            tsince = (JD - TLEs[i].JD)*1440.0; 

            // Call SGP4. Coords are in TEME. 
            LgmSgp_SGP4( tsince, s );
            Uteme.x = s->X/WGS84_A; Uteme.y = s->Y/WGS84_A; Uteme.z = s->Z/WGS84_A;

            // Example of converting TEME->GSM coords.
            Lgm_Convert_Coords( &Uteme, &Ugse, TEME_TO_GSM, c );

            // Dump result to file
            fprintf(fp, "%8ld %.10lf %10.6lf %10.6lf %10.6lf\n", Date, UT, Ugse.x, Ugse.y, Ugse.z  );

        }
    }

    fclose(fp);


    Lgm_free_ctrans( c );
    free( s );
    free( TLEs );


    return(0);
}
