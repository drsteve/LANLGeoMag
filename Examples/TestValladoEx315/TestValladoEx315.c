#include <Lgm_CTrans.h>

/*
 *  Verify that we get the same results as those obtained in Example 3-15 of
 *  Vallado's "Fundamentals of Astrodynamics and Applications", 3rd Ed., Space
 *  Technology Library, Microcosm Press, 2007. Note that there is an
 *  inconsistency in the last digit of the ITRF vector between examples 3-14
 *  and 3-15 even though they are intended to be to same. I don't know to what
 *  extent this means there are bugs in Vallado's results...
 */

int main( int argc, char *argv[] ){

    long int    Date;
    double      UTC;
    Lgm_Vector  u, v;
    Lgm_Vector  Vpef, Vtod, Vmod, Vgcrf;
    Lgm_CTrans  *c = Lgm_init_ctrans( 1 );


    printf("\n\n*******************************************************************************************\n");
    printf("\"Test of Vallado's Example 3-15. Performing IAU-76/FK5 Reduction.\"\n");
    printf("  (See p. 235, Vallado, \"Fundamentals of Astrodynamics and Applications\", 3rd Ed., 2007.)\n");
    printf("\n\nGiven:\n");
    printf("        r_itrf = -1033.4793830 Ihat + u.y =  7901.2952754 Jhat + 6380.3565958 Khat\n");
    printf("Find:\n");
    printf("        r_gcrf  on April 6, 2004, 07:51:28.386009UTC\n");
    printf("\n\n*******************************************************************************************\n");


    printf("\n\n");
    printf("    Setting Date and Time.\n");
    printf("   ------------------------\n");
    Date = 20040406;
    UTC  = 7.0 + 51.0/60.0 + 28.386009/3600.0;
    printf("    Date: %8ld\n", Date);
    printf("     UTC: %15.12lf   ( ", UTC); Lgm_Print_HMSd( UTC ); printf(" )\n");


    printf("\n\n");
    printf("    Setting Earth Orientation Parameters.\n");
    printf("   --------------------------------------\n");
    c->DUT1  = -0.4399619;
    c->xp    = -0.140682; // arcsec
    c->yp    =  0.333309; // arcsec
    if ( ( argc > 1 ) && ( !strcmp(argv[1], "-noeop") ) ){
        c->ddPsi =  0.0; // arcsec
        c->ddEps =  0.0; // arcsec
        Vpef.x = -1033.4750313; Vpef.y = 7901.3055856; Vpef.z = 6380.3445328;
        Vtod.x =  5094.5147804; Vtod.y = 6127.3664612; Vtod.z = 6380.3445328;
        Vmod.x =  5094.0290167; Vmod.y = 6127.8709363; Vmod.z = 6380.2478885;
        Vgcrf.x =  5102.5096; Vgcrf.y = 6123.01152; Vgcrf.z = 6378.1363;
    } else {
        c->ddPsi = -0.052195; // arcsec
        c->ddEps = -0.003875; // arcsec
        Vpef.x = -1033.4750313; Vpef.y = 7901.3055856; Vpef.z = 6380.3445328;
        Vtod.x =  5094.5147804; Vtod.y = 6127.3664612; Vtod.z = 6380.3445328;
        Vmod.x =  5094.0283745; Vmod.y = 6127.8708164; Vmod.z = 6380.2485164;
        Vgcrf.x =  5102.508953; Vgcrf.y = 6123.011396; Vgcrf.z = 6378.136937;
    }
    printf("    DUT1:  %15.12lf   ( ", c->DUT1);         Lgm_Print_HMSd( c->DUT1/3600.0 );  printf(" )\n");
    printf("      xp:  %15.12lf   ( ", c->xp/3600.0);    Lgm_Print_DMSd( c->xp/3600.0 );    printf(" )\n");
    printf("      yp:  %15.12lf   ( ", c->yp/3600.0);    Lgm_Print_DMSd( c->yp/3600.0 );    printf(" )\n");
    printf("   ddPsi:  %15.12lf   ( ", c->ddPsi/3600.0); Lgm_Print_DMSd( c->ddPsi/3600.0 ); printf(" )\n");
    printf("   ddEps:  %15.12lf   ( ", c->ddEps/3600.0); Lgm_Print_DMSd( c->ddEps/3600.0 ); printf(" )\n");

    printf("\n\n");
    printf("    Computing Coordinate Transformations.\n");
    printf("   --------------------------------------\n");
    Lgm_Set_Coord_Transforms( Date, UTC, c);


    printf("\n\n");
    printf("    Setting ITRF Coordinates (km).\n");
    printf("   --------------------------------\n");
    u.x = -1033.4793830;
    u.y =  7901.2952754;
    u.z =  6380.3565958;
    printf("    u_itrf: %.9lf  %.9lf  %.9lf\n", u.x, u.y, u.z);

    printf("\n\n");
    printf("    Transforming to PEF Coordinates (km).\n");
    printf("   ---------------------------------------\n");
    Lgm_Convert_Coords( &u, &v, WGS84_TO_PEF, c);
    printf("     u_pef: %15.9lf  %15.9lf  %15.9lf\n", v.x, v.y, v.z);
    printf("     u_pef: %15.9lf  %15.9lf  %15.9lf  (Vallado's result)\n", Vpef.x, Vpef.y, Vpef.z);
    printf("      DIFF: %15.9lf  %15.9lf  %15.9lf  (LGM - Vallado's result)\n", v.x - Vpef.x, v.y - Vpef.y, v.z - Vpef.z);

    printf("\n\n");
    printf("    Transforming to TOD Coordinates (km).\n");
    printf("   ---------------------------------------\n");
    Lgm_Convert_Coords( &u, &v, WGS84_TO_TOD, c);
    printf("     u_tod: %15.9lf  %15.9lf  %15.9lf\n", v.x, v.y, v.z);
    printf("     u_tod: %15.9lf  %15.9lf  %15.9lf  (Vallado's result)\n", Vtod.x, Vtod.y, Vtod.z);
    printf("      DIFF: %15.9lf  %15.9lf  %15.9lf  (LGM - Vallado's result)\n", v.x - Vtod.x, v.y - Vtod.y, v.z - Vtod.z);

    printf("\n\n");
    printf("    Transforming to MOD Coordinates (km).\n");
    printf("   ---------------------------------------\n");
    Lgm_Convert_Coords( &u, &v, WGS84_TO_MOD, c);
    printf("     u_mod: %.8lf  %.8lf  %.8lf\n", v.x, v.y, v.z);
    printf("     u_mod: %15.9lf  %15.9lf  %15.9lf  (Vallado's result)\n", Vmod.x, Vmod.y, Vmod.z);
    printf("      DIFF: %15.9lf  %15.9lf  %15.9lf  (LGM - Vallado's result)\n", v.x - Vmod.x, v.y - Vmod.y, v.z - Vmod.z);


    printf("\n\n");
    printf("    Transforming to GCRF Coordinates (km).\n");
    printf("   ---------------------------------------\n");
    Lgm_Convert_Coords( &u, &v, WGS84_TO_EME2000, c);
    printf("     u_gcrf: %15.9lf  %15.9lf  %15.9lf\n", v.x, v.y, v.z);
    printf("     u_gcrf: %15.9lf  %15.9lf  %15.9lf  (Vallado's result)\n", Vgcrf.x, Vgcrf.y, Vgcrf.z);
    printf("      DIFF: %15.9lf  %15.9lf  %15.9lf  (LGM - Vallado's result)\n", v.x - Vgcrf.x, v.y - Vgcrf.y, v.z - Vgcrf.z);


if (0==1){
    u.x = -9060.47373569/WGS84_A;
    u.y =  4658.70952502/WGS84_A;
    u.z =  813.68673153/WGS84_A;
    Date = 2000182;
    UTC  = 0.78495062*24.0;
c->DUT1 = 0.0;
c->xp = 0.0; // radians
c->yp = 0.0; // radians
c->ddPsi = 0.0;
c->ddEps = 0.0;
    Lgm_Set_Coord_Transforms( Date, UTC, c);
    Lgm_Convert_Coords( &u, &v, TEME_TO_EME2000, c);
    printf("u = %.8lf  %.8lf  %.8lf\n", u.x*WGS84_A, u.y*WGS84_A, u.z*WGS84_A);
    printf("v = %.8lf  %.8lf  %.8lf\n", v.x*WGS84_A, v.y*WGS84_A, v.z*WGS84_A);
}

    Lgm_free_ctrans( c );

    return(0);

}
