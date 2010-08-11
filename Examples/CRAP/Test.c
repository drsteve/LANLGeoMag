#include <Lgm_CTrans.h>
const char *sMonth[] = { "", "January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December" };
main(){

WritePsdHeader( "puke.txt", "mgh", "mithril" );

}


void WritePsdHeader( char *Filename, char *UserName, char *Machine ){

    int         i, j, Year, Month, Day, HH, MM, SS;
    char        Str[80];
    FILE        *fp;
    long int    CreateDate;
    double      JD, UTC;
    Lgm_CTrans  *c = Lgm_init_ctrans(0);

    int         nAlpha;
    double      a, Alpha[100];

    for (nAlpha=0,a=5.0; a<=90.0; a+=5,++nAlpha) Alpha[nAlpha] = a;


    JD = Lgm_GetCurrentJD(c);
    CreateDate = Lgm_JD_to_Date( JD, &Year, &Month, &Day, &UTC );
    Lgm_UT_to_HMS( UTC, &HH, &MM, &SS );

    /*
     * Write Header
     */
    fp = fopen(Filename, "wb");
    fprintf( fp, "%% Spacecraft:  %s\n", "RBSPA" );
    fprintf( fp, "%% Field Model:  %s\n", "T89" );
    fprintf( fp, "%% nAlpha:  %d; ", nAlpha ); for (i=0; i<nAlpha; i++) fprintf(fp, " %g", Alpha[i]); fprintf( fp, ";   Units: Degrees\n");
    fprintf( fp, "%%\n");
    fprintf( fp, "%% File Contents    :  L* values and Phase Space Densities at Constant Mu and K\n");
    fprintf( fp, "%% File Created at  :  %02d:%02d:%02d UTC  %s %02d %4d\n", HH, MM, SS, sMonth[Month], Day, Year );
    fprintf( fp, "%% File Created by  :  %s\n", UserName );
    fprintf( fp, "%% File Created on  :  %s\n", Machine );
    fprintf( fp, "%% File Fmt Version :  %s\n", "" );
    fprintf( fp, "%%\n");
    fprintf( fp, "%% Description of Variables:\n");
    fprintf( fp, "%%     nK:             Number, of K values used.\n");
    fprintf( fp, "%%\n");
    fprintf( fp, "%%     nMu:            Number, of Mu values used.\n");
    fprintf( fp, "%%\n");
    fprintf( fp, "%%     Date:           The date. In YYYMMDD format.\n");
    fprintf( fp, "%%\n");
    fprintf( fp, "%%     UTC:            Universal Time (Coordinated). In decimal hours.\n");
    fprintf( fp, "%%\n");
    fprintf( fp, "%%     Lat:            Geographic latitude of S/C. In units of Deg.\n");
    fprintf( fp, "%%\n");
    fprintf( fp, "%%     Lon:            Geographic Longitude of S/C. In units of Deg.\n");
    fprintf( fp, "%%\n");
    fprintf( fp, "%%     Rad:            Geocentric Radial Distance of S/C. In units of km.\n");
    fprintf( fp, "%%\n");
    fprintf( fp, "%%     Xgsm:           X-compoentent of GSM position vector. In units of Re.\n");
    fprintf( fp, "%%\n");
    fprintf( fp, "%%     Ygsm:           Y-compoentent of GSM position vector. In units of Re.\n");
    fprintf( fp, "%%\n");
    fprintf( fp, "%%     Zgsm:           Z-compoentent of GSM position vector. In units of Re.\n");
    fprintf( fp, "%%\n");
    fprintf( fp, "%%     B_Used:         The Magnetic Field strength that was used for\n");
    fprintf( fp, "%%                     computations (it will be either Bmod or Bobs). In units\n");
    fprintf( fp, "%%                     of nT.\n");
    fprintf( fp, "%%\n");
    fprintf( fp, "%%     B_Mod:          Magnetic field strength given by the field model (in nT)\n");
    fprintf( fp, "%%                     In units of nT.\n");
    fprintf( fp, "%%\n");
    fprintf( fp, "%%     B_Obs:          Magnetic field strength observed at the S/C (No Obs if <\n");
    fprintf( fp, "%%                     0). In units of nT.\n");
    fprintf( fp, "%%\n");
    fprintf( fp, "%%     M_Used:         The Magnetic Dipole Moment that was used to convert\n");
    fprintf( fp, "%%                     magnetic flux to L*. In units of nT.\n");
    fprintf( fp, "%%\n");
    fprintf( fp, "%%     M_Ref:          A fixed reference magnetic dipole moment for converting\n");
    fprintf( fp, "%%                     magnetic flux to L*. In units of nT.\n");
    fprintf( fp, "%%\n");
    fprintf( fp, "%%     M_IGRF:         Time-dependant magnetic dipole moment (probably shouldn't\n");
    fprintf( fp, "%%                     be used for converting Magnetic Flux to L*, but it\n");
    fprintf( fp, "%%                     sometimes is). In units of nT.\n");
    fprintf( fp, "%%\n");
    fprintf( fp, "%%     L*0-(nK-1):     The nK L* values obtained. (One associated with each of\n");
    fprintf( fp, "%%                     the nK K values). L* values are dimensionless quantities.\n");
    fprintf( fp, "%%\n");
    fprintf( fp, "%%     Bm0-(nK-1):     The nK Bmirror values computed. (One associated with each\n");
    fprintf( fp, "%%                     of the nK K values). In units of nT.\n");
    fprintf( fp, "%%\n");
    fprintf( fp, "%%     I0-(nK-1):      The nK Integral Invariant values computed. (One\n");
    fprintf( fp, "%%                     associated with each of the nK K values). In units of Re.\n");
    fprintf( fp, "%%\n");
    fprintf( fp, "%%     Alpha0-(nK-1):  The nK local pitch angles computed. (One associated with\n");
    fprintf( fp, "%%                     each of the nK K values). In degrees.\n");
    fprintf( fp, "%% \n");
    fprintf( fp, "%%     fp(Ki,Muj):     The phase space density computed at the ith K and jth Mu\n");
    fprintf( fp, "%%                     value. In units of (c/cm/MeV)^3. \n");

    // Column Header
    fprintf( fp, "%% %-10s", "Date" );
    fprintf( fp, " %13s", "UT" );
    fprintf( fp, " %12s", "Lat" );
    fprintf( fp, " %12s", "Lon" );
    fprintf( fp, " %12s", "Rad" );
    fprintf( fp, " %12s", "Xgsm" );
    fprintf( fp, " %12s", "Ygsm" );
    fprintf( fp, " %12s", "Zgsm" );
    fprintf( fp, " %12s", "B_Used" );
    fprintf( fp, " %12s", "B_Mod" );
    fprintf( fp, " %12s", "B_Obs" );
    fprintf( fp, " %12s", "M_Used" );
    fprintf( fp, " %12s", "M_Ref" );
    fprintf( fp, " %12s", "M_IGRF" );
    fprintf(fp, "    ");
    for (i=0; i<nAlpha; i++) { sprintf( Str, "L*%d", i ); fprintf(fp, " %12s", Str ); }
    fprintf(fp, "    ");
    for (i=0; i<nAlpha; i++) { sprintf( Str, "Bm%d", i ); fprintf(fp, " %12s", Str ); }
    fprintf(fp, "    ");
    for (i=0; i<nAlpha; i++) { sprintf( Str, "I%d", i ); fprintf(fp, " %12s", Str ); }
/*
    fprintf(fp, "    ");
    for (i=0; i<nK; i++) {
        for (j=0; j<PsdAssim->nMu; j++) {
            sprintf( Str, "fp(K%d,Mu%d)", i, j );
            fprintf(fp, " %14s", Str);
        }
    }
*/
    fprintf(fp, "\n");
    // Units/Format
    fprintf( fp, "%% %-10s", "YYYYMMDD" );
    fprintf( fp, " %13s", "Hours" );
    fprintf( fp, " %12s", "Deg." );
    fprintf( fp, " %12s", "Deg." );
    fprintf( fp, " %12s", "km" );
    fprintf( fp, " %12s", "Re" );
    fprintf( fp, " %12s", "Re" );
    fprintf( fp, " %12s", "Re" );
    fprintf( fp, " %12s", "nT" );
    fprintf( fp, " %12s", "nT" );
    fprintf( fp, " %12s", "nT" );
    fprintf( fp, " %12s", "nT" );
    fprintf( fp, " %12s", "nT" );
    fprintf( fp, " %12s", "nT" );
    fprintf(fp, "    ");
    for (i=0; i<nAlpha; i++) { sprintf( Str, "dimless" ); fprintf(fp, " %12s", Str ); }
    fprintf(fp, "    ");
    for (i=0; i<nAlpha; i++) { sprintf( Str, "nT" ); fprintf(fp, " %12s", Str ); }
    fprintf(fp, "    ");
    for (i=0; i<nAlpha; i++) { sprintf( Str, "Re" ); fprintf(fp, " %12s", Str ); }
/*
    fprintf(fp, "    ");
    for (i=0; i<nAlpha; i++) {
        for (j=0; j<PsdAssim->nMu; j++) {
            sprintf( Str, "(c/cm/MeV)^3" );
            fprintf(fp, " %14s", Str);
        }
    }
*/
    fprintf(fp, "\n");
    fclose( fp );

    free(c);

}

