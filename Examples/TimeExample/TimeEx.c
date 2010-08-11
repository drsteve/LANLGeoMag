#include <Lgm_CTrans.h>
#include <Lgm_Eop.h>
int main( ) {
    Lgm_CTrans      *c = Lgm_init_ctrans( 0 );
    Lgm_Eop         *e = Lgm_init_eop( 0 );
    Lgm_EopOne      eop;
    int             ny, nm, nd;
    long int        Date;
    double          JD, JDs, JDe, UTC, TAI_minus_TAI, UT1_minus_TAI, UTC_minus_TAI, TT_minus_TAI;
    FILE            *fp;


    // Read in the EOP vals
    Lgm_read_eop( e );

    JDs = Lgm_Date_to_JD( 19620101, 0.0, c ); 
    JDe = Lgm_Date_to_JD( 20100101, 0.0, c ); 

    fp = fopen("Time.dat", "w");
    for ( JD = JDs; JD < JDe; JD += 1.0 ) {

        // Get (interpolate) the EOP vals from the values in the file at the given Julian Date
        Lgm_get_eop_at_JD( JD, &eop, e );

        // Set the EOP vals in the CTrans structure.
        Lgm_set_eop( &eop, c );

        printf("DUT1, DAT = %lf %lf\n", c->DUT1, c->DAT );
        TAI_minus_TAI = 0.0;
        UT1_minus_TAI = c->DUT1 - c->DAT;
        UTC_minus_TAI =  -c->DAT;
        TT_minus_TAI  =  32.184;
        Date = Lgm_JD_to_Date( JD, &ny, &nm, &nd, &UTC );
        fprintf(fp, "%8ld, %lf,   %lf, %lf,  %lf, %lf, %lf, %lf, %lf, %lf\n", Date, UTC, JD, 1962.0+(JD-JDs)/365.0, c->DUT1, c->DAT, TAI_minus_TAI, UT1_minus_TAI, UTC_minus_TAI, TT_minus_TAI );

    }
    fclose(fp);
    
    // free structures
    Lgm_free_ctrans( c ); 
    Lgm_destroy_eop( e );

    return(0);

}

