#include <Lgm_CTrans.h>
#include <Lgm_Eop.h>

int main(){

    Lgm_CTrans      *c = Lgm_init_ctrans( 0 ); // more compact declaration
    Lgm_Eop         *e = Lgm_init_eop( 0 );
    Lgm_EopOne      eop;
    Lgm_Vector      u_mod, u_gei;
    long int        Date;
    double          UTC, UT, JD, JD0;
    int             HH, ny, nm, nd;
    double          RA_JPL, DEC_JPL, RA_LOW, DEC_LOW, RA_HIGH, DEC_HIGH;
    double          RA_hh, RA_mm, RA_ss, DEC_dd, DEC_mm, DEC_ss, sgn;
    char            Line[300];
    FILE            *fp, *fpout;


    Date = 20070101;
    UTC  = 0.0;
    JD   = Lgm_Date_to_JD( Date, UTC, c ); // Compute JD
    JD0  = JD;


    // Read in the EOP vals
    Lgm_read_eop( e );



    fp    = fopen("jpl.dat", "r");
    fpout = fopen("SUN.dat", "w");
    fprintf(fpout, "%15s  %15s  %15s  %15s  %15s  %15s  %15s  %15s  %15s  %15s  %15s  %15s\n", 
                "Julian Date", "UTC", "RA_JPL", "DEC_JPL", "RA_LOW", "DEC_LOW", "RA_HIGH", "DEC_HIGH", 
                "RA_JPL-RA_LOW (arcsec)", "DEC_JPL-DEC_LOW (arcsec)", "RA_JPL-RA_HIGH (arcsec)", "DEC_JPL-DEC_HIGH (arcsec)" );
    while ( fgets( Line, 256, fp ) != NULL ) {

        sscanf( Line, "%lf %lf %lf %lf %lf %lf", &RA_hh, &RA_mm, &RA_ss, &DEC_dd, &DEC_mm, &DEC_ss );

        /*
         * Position of Sun according to JPL Ephemerides
         */
        RA_JPL  = RA_hh + RA_mm/60.0 + RA_ss/3600.0;
        // the sign is carried by the degrees field. But if its -00 this gets read as
        // 0. So we actually have to test for the - sign. PITA.
        sgn = (Line[12] == '-') ? -1.0 : 1.0; 
        DEC_JPL = sgn*( fabs(DEC_dd) + DEC_mm/60.0 + DEC_ss/3600.0);


        /*
         * Compute low- and high-accuracy positions of Sun with LGM
         */
        Date = Lgm_JD_to_Date( JD, &ny, &nm, &nd, &UTC );
        Lgm_get_eop_at_JD( JD, &eop, e ); // Get (interpolate) the EOP vals from the values in the file at the given Julian Date
        Lgm_set_eop( &eop, c );           // Set the EOP vals in the CTrans structure.
        Lgm_Set_Coord_Transforms( Date, UTC, c );

        Lgm_Radec_to_Cart( c->RA_sun, c->DEC_sun, &u_mod );
        Lgm_Convert_Coords( &u_mod, &u_gei, MOD_TO_GEI2000, c );
        RA_LOW = atan2( u_gei.y, u_gei.x)*DegPerRad; if (RA_LOW<0.0) RA_LOW += 360.0;
        RA_LOW /= 15.0;
        DEC_LOW = asin( u_gei.z )*DegPerRad;

        Lgm_Radec_to_Cart( c->RA_sun_ha, c->DEC_sun_ha, &u_mod );
        Lgm_Convert_Coords( &u_mod, &u_gei, MOD_TO_GEI2000, c );
        RA_HIGH = atan2( u_gei.y, u_gei.x)*DegPerRad; if (RA_HIGH<0.0) RA_HIGH += 360.0;
        RA_HIGH /= 15.0;
        DEC_HIGH = asin( u_gei.z )*DegPerRad;


        fprintf(fpout, "%15.9lf, %15.9lf, %15.9lf, %15.9lf, %15.9lf, %15.9lf, %15.9lf, %15.9lf, %15.9lf, %15.9lf, %15.9lf, %15.9lf\n", 
                        JD, UTC, RA_JPL, DEC_JPL, RA_LOW, DEC_LOW, RA_HIGH, DEC_HIGH, 
                        (RA_JPL-RA_LOW)*3600.0*15.0, (DEC_JPL-DEC_LOW)*3600.0, (RA_JPL-RA_HIGH)*3600.0*15.0, (DEC_JPL-DEC_HIGH)*3600.0 );

        JD  += 1.0/24.0;

    }
    fclose(fp);
    fclose(fpout);

    Lgm_free_ctrans( c );
    Lgm_destroy_eop( e );



    exit(0);

}
