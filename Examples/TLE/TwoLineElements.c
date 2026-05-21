#include <stdio.h>
#include <stdlib.h>
#include <Lgm_CTrans.h>
#include <Lgm_Eop.h>
#include <Lgm_Sgp.h>

/*
 * This example shows how to compute positions of objects from TLEs
 * (Two Line Elements). 
 */

int main( int argc, char *argv[] ){

    int         Year, Month, Day;
    int         StartYear, StartMonth, StartDay, StartDoy;
    int         EndYear, EndMonth, EndDay, EndDoy;
    long int    Date, StartDate, EndDate;
    double      UT, tsince, JD, StartUT, EndUT, StartJD, EndJD, JDinc;
    Lgm_CTrans  *c = Lgm_init_ctrans( 0 );
    Lgm_Vector  Ugse, Uteme;
    Lgm_Eop     *e = Lgm_init_eop( 0 );
    Lgm_EopOne  eop;
    int         nTLEs;
    char        Line0[100], Line1[100], Line2[100], *ptr;
    char        *InputFile  = "SSA_Forte.txt";
    char        *OutputFile = "Output.txt";
    FILE        *fp;


    /* 
     * Open input file and extract:
     *   1) the 3 TLE lines
     *   2) Start Date and Time
     *   3) End Date and Time
     */
    if ( (fp = fopen( InputFile, "r" )) != NULL ) {
        fgets( Line0, 99, fp );
        fgets( Line1, 99, fp );
        fgets( Line2, 99, fp );
        fscanf( fp, "%ld %lf", &StartDate, &StartUT );
        fscanf( fp, "%ld %lf", &EndDate, &EndUT );
    } else {
        printf( "Couldnt open file %s for reading\n", InputFile );
        exit( 1 );
    }
    fclose( fp );

    /*
     * Remove any extraneous newline and/or linefeeds at the end of the strings.
     */
    if ( (ptr = strstr(Line0, "\n")) != NULL ) *ptr = '\0'; 
    if ( (ptr = strstr(Line0, "\r")) != NULL ) *ptr = '\0';
    if ( (ptr = strstr(Line1, "\n")) != NULL ) *ptr = '\0'; 
    if ( (ptr = strstr(Line1, "\r")) != NULL ) *ptr = '\0';
    if ( (ptr = strstr(Line2, "\n")) != NULL ) *ptr = '\0'; 
    if ( (ptr = strstr(Line2, "\r")) != NULL ) *ptr = '\0';


    /*
     * Alloc some memory for the SgpInfo structure and the TLEs array (here we
     * only have a single element in the array)
     */
    _SgpInfo *s   = (_SgpInfo *)calloc( 1, sizeof(_SgpInfo) );
    _SgpTLE *TLEs = (_SgpTLE *)calloc( 1, sizeof(_SgpTLE) );

    
    /*
     * Read in TLEs from the Line0, Line1 and Line2 strings. nTLEs must be
     * initialized to zero by the user.
     */
    nTLEs = 0;
    LgmSgp_ReadTlesFromStrings( Line0, Line1, Line2, &nTLEs, TLEs, 4 );


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

    JDinc = 60.0/1440.0; // 60-min increment

    if ( (fp = fopen( OutputFile, "w" )) == NULL ) {
        printf( "Couldnt open file %s for writing\n", OutputFile );
        exit( 1 );
    }

    // init SGP4
    LgmSgp_SGP4_Init( s, &TLEs[0] );
    printf("%%%s\n", TLEs[0].Line0 );
    printf("%%%s\n", TLEs[0].Line1 );
    printf("%%%s\n", TLEs[0].Line2 );
    fprintf(fp, "%%%s\n", TLEs[0].Line0 );
    fprintf(fp, "%%%s\n", TLEs[0].Line1 );
    fprintf(fp, "%%%s\n", TLEs[0].Line2 );

    // Read in the EOP vals
    Lgm_read_eop( e );

    // loop over specified time range
    for ( JD = StartJD; JD <= EndJD; JD += JDinc ) {

        // Convert the current JD back to Date/UT etc..
        Lgm_jd_to_ymdh ( JD, &Date, &Year, &Month, &Day, &UT );

        // Get (interpolate) the EOP vals from the values in the file at the given Julian Date
        Lgm_get_eop_at_JD( JD, &eop, e );

        // Set the EOP vals in the CTrans structure.
        Lgm_set_eop( &eop, c );

        // Set up the trans matrices
        Lgm_Set_Coord_Transforms( Date, UT, c );
    
        // "time since" in minutes (thats what SGP4 wants)
        tsince = (JD - TLEs[0].JD)*1440.0; 

        // Call SGP4. Coords are in TEME. 
        LgmSgp_SGP4( tsince, s );
        Uteme.x = s->X/WGS84_A; Uteme.y = s->Y/WGS84_A; Uteme.z = s->Z/WGS84_A;
printf("*%8ld %.10lf %15.9lf %15.9lf %15.9lf\n", Date, UT, Uteme.x, Uteme.y, Uteme.z  );

        // Example of converting TEME->GSE coords.
        Lgm_Convert_Coords( &Uteme, &Ugse, TEME_TO_WGS84, c );
        //Lgm_Convert_Coords( &Uteme, &Ugse, TEME_TO_GEO, c );

        // Dump result to file
        printf("%8ld %.10lf %15.9lf %15.9lf %15.9lf\n", Date, UT, Ugse.x, Ugse.y, Ugse.z  );
        fprintf(fp, "%8ld %.10lf %10.6lf %10.6lf %10.6lf\n", Date, UT, Ugse.x, Ugse.y, Ugse.z  );

    }
    fclose(fp);

    Lgm_free_ctrans( c );
    free( s );
    free( TLEs );
    Lgm_destroy_eop( e );

    return(0);
}
