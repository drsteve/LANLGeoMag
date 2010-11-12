#include <Lgm_CTrans.h>
#include <Lgm_MagEphemInfo.h>
const char *sMonth[] = { "", "January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December" };

/*
main(){

    FILE        *fp;

    fp = fopen( "puke.txt", "wb");
    WriteMagEphemHeader( fp, "mgh", "mithril" );
    WriteMagEphemData( fp );

}
*/


void WriteMagEphemHeader( FILE *fp, char *Spacecraft, char *IntModel, char *ExtModel, double Kp, double Dst, Lgm_MagEphemInfo *m ){

    int         i, Year, Month, Day, HH, MM, SS;
    char        Str[80];
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
    fprintf( fp, "# Spacecraft:  %s\n", Spacecraft );
//    fprintf( fp, "# Field Model:  %s\n", MagModel );
    fprintf( fp, "# nAlpha:  %d; Alphas: ", m->nAlpha ); for (i=0; i<m->nAlpha; i++) fprintf(fp, " %g", m->Alpha[i]); fprintf( fp, ";   Units: Degrees\n");
    fprintf( fp, "#\n");
    fprintf( fp, "# File Contents    :  Magnetic Empherii for spacecraft trajectory.\n");
    fprintf( fp, "# File Created at  :  %02d:%02d:%02d UTC  %s %02d %4d\n", HH, MM, SS, sMonth[Month], Day, Year );
//    fprintf( fp, "# File Created by  :  %s\n", UserName );
//    fprintf( fp, "# File Created on  :  %s\n", Machine );
//    fprintf( fp, "# File Fmt Version :  %s\n", "" );
    fprintf( fp, "#\n");
    fprintf( fp, "# Description of Variables:\n");
    fprintf( fp, "#          DateTime:  The date and time in ISO 8601 compliant format.\n");
    fprintf( fp, "#\n");
    fprintf( fp, "#              Date:  The date. In YYYMMDD format.\n");
    fprintf( fp, "#\n");
    fprintf( fp, "#               DOY:  Day of Yeay. DDD\n");
    fprintf( fp, "#\n");
    fprintf( fp, "#               UTC:  Universal Time (Coordinated). In decimal hours.\n");
    fprintf( fp, "#\n");
    fprintf( fp, "#       Julian Date:  Julian Date. In decimal days.\n");
    fprintf( fp, "#\n");
    fprintf( fp, "#          GPS Time:  Number of SI seconds since 0h Jan 6, 1980.\n");
    fprintf( fp, "#\n");
    fprintf( fp, "#              Xgeo:  X-component of Geographic position of S/C. In units of Re.\n");
    fprintf( fp, "#              Ygeo:  Y-component of Geographic position of S/C. In units of Re.\n");
    fprintf( fp, "#              Zgeo:  Z-component of Geographic position of S/C. In units of Re.\n");
    fprintf( fp, "#\n");
    fprintf( fp, "#              Xgsm:  X-compontent of GSM position vector. In units of Re.\n");
    fprintf( fp, "#              Ygsm:  Y-compontent of GSM position vector. In units of Re.\n");
    fprintf( fp, "#              Zgsm:  Z-compontent of GSM position vector. In units of Re.\n");
    fprintf( fp, "#\n");
    fprintf( fp, "#               Xsm:  X-compontent of SM position vector. In units of Re.\n");
    fprintf( fp, "#               Ysm:  Y-compontent of SM position vector. In units of Re.\n");
    fprintf( fp, "#               Zsm:  Z-compontent of SM position vector. In units of Re.\n");
    fprintf( fp, "#\n");
    fprintf( fp, "#              Xgei:  X-compontent of GEI position vector. In units of Re.\n");
    fprintf( fp, "#              Ygei:  Y-compontent of GEI position vector. In units of Re.\n");
    fprintf( fp, "#              Zgei:  Z-compontent of GEI position vector. In units of Re.\n");
    fprintf( fp, "#\n");
    fprintf( fp, "#              Xgse:  X-compontent of GSE position vector. In units of Re.\n");
    fprintf( fp, "#              Ygse:  Y-compontent of GSE position vector. In units of Re.\n");
    fprintf( fp, "#              Zgse:  Z-compontent of GSE position vector. In units of Re.\n");
    fprintf( fp, "#\n");
    fprintf( fp, "#          IntModel:  Internal magnetic field model used.\n");
    fprintf( fp, "#          ExtModel:  External magnetic field model used.\n");
    fprintf( fp, "#                Kp:  Kp index value.\n");
    fprintf( fp, "#               Dst:  Dst index value. In units of nT.\n");
    fprintf( fp, "#\n");
    fprintf( fp, "#             Bsc_x:  X-component of magnetic field vector at S/C (in GSM coords). In units of nT.\n");
    fprintf( fp, "#             Bsc_y:  y-component of magnetic field vector at S/C (in GSM coords). In units of nT.\n");
    fprintf( fp, "#             Bsc_z:  z-component of magnetic field vector at S/C (in GSM coords). In units of nT.\n");
    fprintf( fp, "#               Bsc:  Magnitude of magnetic field vector at S/C. In units of nT.\n");
    fprintf( fp, "#\n");
    fprintf( fp, "#     FieldLineType:  Description of the type of field line the S/C is on.\n");
    fprintf( fp, "#                     Can be one of 4 types:\n");
    fprintf( fp, "#                         LGM_CLOSED      - FL hits Earth at both ends.\n");
    fprintf( fp, "#                         LGM_OPEN_N_LOBE - FL is an OPEN field line rooted in the Northern polar cap.\n");
    fprintf( fp, "#                         LGM_OPEN_S_LOBE - FL is an OPEN field line rooted in the Southern polar cap.\n");
    fprintf( fp, "#                         LGM_OPEN_IMF    - FL does not hit Earth at eitrher end.\n");
    fprintf( fp, "#\n");
    fprintf( fp, "#       Foot_N_Xgsm:  X-compontent of GSM position vector of Northern Magnetic Footpoint. (Re).\n");
    fprintf( fp, "#       Foot_N_Ygsm:  Y-compontent of GSM position vector of Northern Magnetic Footpoint. (Re).\n");
    fprintf( fp, "#       Foot_N_Zgsm:  Z-compontent of GSM position vector of Northern Magnetic Footpoint. (Re).\n");
    fprintf( fp, "#\n");
    fprintf( fp, "#       Foot_N_Xgeo:  X-compontent of GEO position vector of Northern Magnetic Footpoint. (Re).\n");
    fprintf( fp, "#       Foot_N_Ygeo:  Y-compontent of GEO position vector of Northern Magnetic Footpoint. (Re).\n");
    fprintf( fp, "#       Foot_N_Zgeo:  Z-compontent of GEO position vector of Northern Magnetic Footpoint. (Re).\n");
    fprintf( fp, "#\n");
    fprintf( fp, "#    Foot_N_GeodLat:  Geodetic Latitude of Northern Magnetic Footpoint. (Degrees).\n");
    fprintf( fp, "#    Foot_N_GeodLon:  Geodetic Longitude of Northern Magnetic Footpoint. (Degrees).\n");
    fprintf( fp, "# Foot_N_GeodHeight:  Geodetic Height of Northern Magnetic Footpoint. (km).\n");
    fprintf( fp, "#\n");
    fprintf( fp, "#         Bfn_geo_x:  X-component of magnetic field vector at Northern Footpoint (in GEO coords). In units of nT.\n");
    fprintf( fp, "#         Bfn_geo_y:  y-component of magnetic field vector at Northern Footpoint (in GEO coords). In units of nT.\n");
    fprintf( fp, "#         Bfn_geo_z:  z-component of magnetic field vector at Northern Footpoint (in GEO coords). In units of nT.\n");
    fprintf( fp, "#           Bfn_geo:  Magnitude of magnetic field vector at Northern Footpoint. In units of nT.\n");
    fprintf( fp, "#\n");
    fprintf( fp, "#             Bfn_x:  X-component of magnetic field vector at Northern Footpoint (in GSM coords). In units of nT.\n");
    fprintf( fp, "#             Bfn_y:  y-component of magnetic field vector at Northern Footpoint (in GSM coords). In units of nT.\n");
    fprintf( fp, "#             Bfn_z:  z-component of magnetic field vector at Northern Footpoint (in GSM coords). In units of nT.\n");
    fprintf( fp, "#               Bfn:  Magnitude of magnetic field vector at Northern Footpoint. In units of nT.\n");
    fprintf( fp, "#\n");
    fprintf( fp, "# Alpha_Loss_Cone_n:  Value of Northern Loss Cone angle. asin( sqrt(Bsc/Bfn) ). In units of Degrees.\n");
    fprintf( fp, "#\n");
    fprintf( fp, "#       Foot_S_Xgsm:  X-compontent of GSM position vector of Southern Magnetic Footpoint. (Re).\n");
    fprintf( fp, "#       Foot_S_Ygsm:  Y-compontent of GSM position vector of Southern Magnetic Footpoint. (Re).\n");
    fprintf( fp, "#       Foot_S_Zgsm:  Z-compontent of GSM position vector of Southern Magnetic Footpoint. (Re).\n");
    fprintf( fp, "#\n");
    fprintf( fp, "#       Foot_S_Xgeo:  X-compontent of GEO position vector of Southern Magnetic Footpoint. (Re).\n");
    fprintf( fp, "#       Foot_S_Ygeo:  Y-compontent of GEO position vector of Southern Magnetic Footpoint. (Re).\n");
    fprintf( fp, "#       Foot_S_Zgeo:  Z-compontent of GEO position vector of Southern Magnetic Footpoint. (Re).\n");
    fprintf( fp, "#\n");
    fprintf( fp, "#    Foot_S_GeodLat:  Geodetic Latitude of Southern Magnetic Footpoint. (Degrees).\n");
    fprintf( fp, "#    Foot_S_GeodLon:  Geodetic Longitude of Southern Magnetic Footpoint. (Degrees).\n");
    fprintf( fp, "# Foot_S_GeodHeight:  Geodetic Height of Southern Magnetic Footpoint. (km).\n");
    fprintf( fp, "#\n");
    fprintf( fp, "#         Bfs_geo_x:  X-component of magnetic field vector at Southern Footpoint (in GEO coords). In units of nT.\n");
    fprintf( fp, "#         Bfs_geo_y:  y-component of magnetic field vector at Southern Footpoint (in GEO coords). In units of nT.\n");
    fprintf( fp, "#         Bfs_geo_z:  z-component of magnetic field vector at Southern Footpoint (in GEO coords). In units of nT.\n");
    fprintf( fp, "#           Bfn_geo:  Magnitude of magnetic field vector at Southern Footpoint. In units of nT.\n");
    fprintf( fp, "#\n");
    fprintf( fp, "#             Bfs_x:  X-component of magnetic field vector at Southern Footpoint (in GSM coords). In units of nT.\n");
    fprintf( fp, "#             Bfs_y:  y-component of magnetic field vector at Southern Footpoint (in GSM coords). In units of nT.\n");
    fprintf( fp, "#             Bfs_z:  z-component of magnetic field vector at Southern Footpoint (in GSM coords). In units of nT.\n");
    fprintf( fp, "#               Bfs:  Magnitude of magnetic field vector at Southern Footpoint. In units of nT.\n");
    fprintf( fp, "#\n");
    fprintf( fp, "# Alpha_Loss_Cone_s:  Value of Southern Loss Cone angle. asin( sqrt(Bsc/Bfn) ). In units of Degrees.\n");
    fprintf( fp, "#\n");
    fprintf( fp, "#            Pmin_x:  X-component of GSM position vector of minimum-|B| point. In units of Re.\n");
    fprintf( fp, "#            Pmin_y:  Y-component of GSM position vector of minimum-|B| point. In units of Re.\n");
    fprintf( fp, "#            Pmin_z:  Z-component of GSM position vector of minimum-|B| point. In units of Re.\n");
    fprintf( fp, "#\n");
    fprintf( fp, "#            Bmin_x:  X-component of magnetic field vector at Pmin (in GSM coords). In units of nT.\n");
    fprintf( fp, "#            Bmin_y:  y-component of magnetic field vector at Pmin (in GSM coords). In units of nT.\n");
    fprintf( fp, "#            Bmin_z:  z-component of magnetic field vector at Pmin (in GSM coords). In units of nT.\n");
    fprintf( fp, "#              Bmin:  Magnitude of magnetic field vector at Pmin. In units of nT.\n");
    fprintf( fp, "#\n");
    fprintf( fp, "#            M_used:  The magnetic dipole moment that was used to convert\n");
    fprintf( fp, "#                     magnetic flux to L*. In units of nT.\n");
    fprintf( fp, "#\n");
    fprintf( fp, "#             M_ref:  T fixed reference magnetic dipole moment for converting\n");
    fprintf( fp, "#                     magnetic flux to L*. In units of nT.\n");
    fprintf( fp, "#\n");
    fprintf( fp, "#            M_igrf:  Time-dependant magnetic dipole moment (probably shouldn't\n");
    fprintf( fp, "#                     be used for converting magnetic flux to L*, but it\n");
    fprintf( fp, "#                     sometimes is). In units of nT.\n");
    fprintf( fp, "#\n");
    fprintf( fp, "#    L*0-(nAlpha-1):  The nAlpha L* values obtained. L* values are dimensionless quantities.\n");
    fprintf( fp, "#\n");
    fprintf( fp, "#     L0-(nAlpha-1):  The nAlpha mcilwain L values obtained. L values are dimensionless quantities.\n");
    fprintf( fp, "#\n");
    fprintf( fp, "#    Bm0-(nAlpha-1):  The nAlpha bmirror values computed. In units of nT.\n");
    fprintf( fp, "#\n");
    fprintf( fp, "#     I0-(nAlpha-1):  The nAlpha integral invariant values computed.  In units of Re.\n");
    fprintf( fp, "#\n");
    fprintf( fp, "# Fill Value: %g\n", LGM_FILL_VALUE );
    fprintf( fp, "#\n");
    fprintf( fp, "#\n");
    fprintf( fp, "#\n");

    // column header
    fprintf( fp, "%91s",  "#  +------------------------------------ Date and Time -----------------------------------+" );
    fprintf( fp, " %38s",  " +----- Geographic Coordinates ------+" );
    fprintf( fp, " %38s",  " +--------- GSM Coordinates ---------+" );
    fprintf( fp, " %38s",  " +---------- SM Coordinates ---------+" );
    fprintf( fp, " %38s",  " +------- GEI 2000 Coordinates ------+" );
    fprintf( fp, " %38s",  " +---------- GSE Coordinates --------+" );
    fprintf( fp, " %13s",  " +-Int Model-+" );
    fprintf( fp, " %13s",  " +-Ext Model-+" );
    fprintf( fp, " %6s",   " +-Kp-+" );
    fprintf( fp, " %7s",   " +-Dst-+" );
    fprintf( fp, " %51s",  " +--------- Magnetic Field at SpaceCraft ---------+" );
    fprintf( fp, " %29s",  " +----- Field Line Type ----+" );
    fprintf( fp, " %38s",  " +---- North Mag. Footpoint GSM -----+" );
    fprintf( fp, " %38s",  " +- North Mag. Footpoint Geographic -+" );
    fprintf( fp, " %38s",  " +-- North Mag. Footpoint Geodetic --+" );
    fprintf( fp, " %51s",  " +---- Mag. Field at North Mag. Footpoint GEO ----+" );
    fprintf( fp, " %51s",  " +---- Mag. Field at North Mag. Footpoint GSM ----+" );
    fprintf( fp, " %12s",  " +-N.L.Cone-+" );
    fprintf( fp, " %38s",  " +---- South Mag. Footpoint GSM -----+" );
    fprintf( fp, " %38s",  " +- South Mag. Footpoint Geographic -+" );
    fprintf( fp, " %38s",  " +-- South Mag. Footpoint Geodetic --+" );
    fprintf( fp, " %51s",  " +---- Mag. Field at South Mag. Footpoint GEO ----+" );
    fprintf( fp, " %51s",  " +---- Mag. Field at South Mag. Footpoint GSM ----+" );
    fprintf( fp, " %12s",  " +-S.L.Cone-+" );
    fprintf( fp, " %38s",  " +----- Minimum |B| Point GSM -------+" );
    fprintf( fp, " %51s",  " +---- Magnetic Field at Minimum |B| Pointint ----+" );

    fprintf(fp, "\n");


    // column header
    fprintf( fp, "# %25s", "DateTime" );
    fprintf( fp, " %10s", "Date" );
    fprintf( fp, " %5s",  "DOY" );
    fprintf( fp, " %13s", "UTC" );
    fprintf( fp, " %16s", "Julian Date" );
    fprintf( fp, " %15s", "GPS Time" );

    fprintf( fp, " %12s", "Xgeo" );
    fprintf( fp, " %12s", "Ygeo" );
    fprintf( fp, " %12s", "Zgeo" );

    fprintf( fp, " %12s", "Xgsm" );
    fprintf( fp, " %12s", "Ygsm" );
    fprintf( fp, " %12s", "Zgsm" );

    fprintf( fp, " %12s", "Xsm" );
    fprintf( fp, " %12s", "Ysm" );
    fprintf( fp, " %12s", "Zsm" );

    fprintf( fp, " %12s", "Xgei" );
    fprintf( fp, " %12s", "Ygei" );
    fprintf( fp, " %12s", "Zgei" );

    fprintf( fp, " %12s", "Xgse" );
    fprintf( fp, " %12s", "Ygse" );
    fprintf( fp, " %12s", "Zgse" );

    fprintf( fp, " %14s",  "              " );
    fprintf( fp, " %14s",  "              " );
    fprintf( fp, " %7s",   "       " );
    fprintf( fp, " %8s",   "        " );

    fprintf( fp, " %12s", "Bsc_x" ); // Bsc  gsm
    fprintf( fp, " %12s", "Bsc_y" );
    fprintf( fp, " %12s", "Bsc_z" );
    fprintf( fp, " %12s", "Bsc" );

    fprintf( fp, " %29s",  "" );

    fprintf( fp, " %12s", "Xgsm" ); // n. foot
    fprintf( fp, " %12s", "Ygsm" );
    fprintf( fp, " %12s", "Zgsm" );

    fprintf( fp, " %12s", "Xgeo" );
    fprintf( fp, " %12s", "Ygeo" );
    fprintf( fp, " %12s", "Zgeo" );

    fprintf( fp, " %12s", "GeodLat" );
    fprintf( fp, " %12s", "GeodLon" );
    fprintf( fp, " %12s", "GeodHeight" );

    fprintf( fp, " %12s", "Bfn_geo_x" ); // Bfn  geo
    fprintf( fp, " %12s", "Bfn_geo_y" );
    fprintf( fp, " %12s", "Bfn_geo_z" );
    fprintf( fp, " %12s", "Bfn_geo" );

    fprintf( fp, " %12s", "Bfn_x" ); // Bfn  gsm
    fprintf( fp, " %12s", "Bfn_y" );
    fprintf( fp, " %12s", "Bfn_z" );
    fprintf( fp, " %12s", "Bfn" );

    fprintf( fp, " %12s", "Alpha_LC_N" ); // Northern Loss Cone

    fprintf( fp, " %12s", "Xgsm" ); // s. foot
    fprintf( fp, " %12s", "Ygsm" );
    fprintf( fp, " %12s", "Zgsm" );

    fprintf( fp, " %12s", "Xgeo" );
    fprintf( fp, " %12s", "Ygeo" );
    fprintf( fp, " %12s", "Zgeo" );

    fprintf( fp, " %12s", "GeodLat" );
    fprintf( fp, " %12s", "GeodLon" );
    fprintf( fp, " %12s", "GeodHeight" );

    fprintf( fp, " %12s", "Bfs_geo_x" ); // Bfs  geo
    fprintf( fp, " %12s", "Bfs_geo_y" );
    fprintf( fp, " %12s", "Bfs_geo_z" );
    fprintf( fp, " %12s", "Bfs_geo" );

    fprintf( fp, " %12s", "Bfs_x" ); // Bfs  gsm
    fprintf( fp, " %12s", "Bfs_y" );
    fprintf( fp, " %12s", "Bfs_z" );
    fprintf( fp, " %12s", "Bfs" );

    fprintf( fp, " %12s", "Alpha_LC_N" ); // Southern Loss Cone

    fprintf( fp, " %12s", "Xgsm" ); // Pmin gsm
    fprintf( fp, " %12s", "Ygsm" );
    fprintf( fp, " %12s", "Zgsm" );

    fprintf( fp, " %12s", "Bmin_x" ); // |B|min gsm
    fprintf( fp, " %12s", "Bmin_y" );
    fprintf( fp, " %12s", "Bmin_z" );
    fprintf( fp, " %12s", "Bmin" );

    fprintf( fp, " %12s", "M_used" );
    fprintf( fp, " %12s", "M_ref" );
    fprintf( fp, " %12s", "M_igrf" );
    fprintf(fp, "    ");
    for (i=0; i<m->nAlpha; i++) { sprintf( Str, "L*%d", i ); fprintf(fp, " %12s", Str ); }
    fprintf(fp, "    ");
    for (i=0; i<m->nAlpha; i++) { sprintf( Str, "L%d", i ); fprintf(fp, " %12s", Str ); }
    fprintf(fp, "    ");
    for (i=0; i<m->nAlpha; i++) { sprintf( Str, "Bm%d", i ); fprintf(fp, " %12s", Str ); }
    fprintf(fp, "    ");
    for (i=0; i<m->nAlpha; i++) { sprintf( Str, "I%d", i ); fprintf(fp, " %12s", Str ); }
    fprintf(fp, "\n");
    // units/format
    fprintf( fp, "# %25s", "YYYY-MM-DDTHH:MM:SS.SSSSZ" );
    fprintf( fp, " %10s", "YYYYMMDD" );
    fprintf( fp, " %5s",  "DDD" );
    fprintf( fp, " %13s", "Hours" );
    fprintf( fp, " %16s", "Days" );
    fprintf( fp, " %15s", "Seconds" );

    fprintf( fp, " %12s", "Re" );
    fprintf( fp, " %12s", "Re" );
    fprintf( fp, " %12s", "Re" );

    fprintf( fp, " %12s", "Re" );
    fprintf( fp, " %12s", "Re" );
    fprintf( fp, " %12s", "Re" );

    fprintf( fp, " %12s", "Re" );
    fprintf( fp, " %12s", "Re" );
    fprintf( fp, " %12s", "Re" );

    fprintf( fp, " %12s", "Re" );
    fprintf( fp, " %12s", "Re" );
    fprintf( fp, " %12s", "Re" );

    fprintf( fp, " %12s", "Re" );
    fprintf( fp, " %12s", "Re" );
    fprintf( fp, " %12s", "Re" );

    fprintf( fp, " %14s",  " " );  // Int Model
    fprintf( fp, " %14s",  " " );  // Ext Model
    fprintf( fp, " %7s",   " " );  // Kp
    fprintf( fp, " %8s",   "nT" ); // Dst

    fprintf( fp, " %12s", "nT" );   // Bsc
    fprintf( fp, " %12s", "nT" );
    fprintf( fp, " %12s", "nT" );
    fprintf( fp, " %12s", "nT" );

    fprintf( fp, " %29s",  "" );    // Field Line Type

    fprintf( fp, " %12s", "Re" );   // GSM
    fprintf( fp, " %12s", "Re" );
    fprintf( fp, " %12s", "Re" );

    fprintf( fp, " %12s", "Re" );   // GEO
    fprintf( fp, " %12s", "Re" );
    fprintf( fp, " %12s", "Re" );

    fprintf( fp, " %12s", "Deg." ); // Geodetic
    fprintf( fp, " %12s", "Deg." );
    fprintf( fp, " %12s", "km" );

    fprintf( fp, " %12s", "nT" );   // Bfn geo
    fprintf( fp, " %12s", "nT" );
    fprintf( fp, " %12s", "nT" );
    fprintf( fp, " %12s", "nT" );

    fprintf( fp, " %12s", "nT" );   // Bfn gsm
    fprintf( fp, " %12s", "nT" );
    fprintf( fp, " %12s", "nT" );
    fprintf( fp, " %12s", "nT" );

    fprintf( fp, " %12s", "Deg." ); // Loss Cone North

    fprintf( fp, " %12s", "Re" );   // GSM
    fprintf( fp, " %12s", "Re" );
    fprintf( fp, " %12s", "Re" );

    fprintf( fp, " %12s", "Re" );   // GEO
    fprintf( fp, " %12s", "Re" );
    fprintf( fp, " %12s", "Re" );

    fprintf( fp, " %12s", "Deg." ); // Geodetic
    fprintf( fp, " %12s", "Deg." );
    fprintf( fp, " %12s", "km" );

    fprintf( fp, " %12s", "nT" );   // Bfs geo
    fprintf( fp, " %12s", "nT" );
    fprintf( fp, " %12s", "nT" );
    fprintf( fp, " %12s", "nT" );

    fprintf( fp, " %12s", "nT" );   // Bfs gsm
    fprintf( fp, " %12s", "nT" );
    fprintf( fp, " %12s", "nT" );
    fprintf( fp, " %12s", "nT" );

    fprintf( fp, " %12s", "Deg." ); // Loss Cone South

    fprintf( fp, " %12s", "Re" );   // Pmin
    fprintf( fp, " %12s", "Re" );
    fprintf( fp, " %12s", "Re" );

    fprintf( fp, " %12s", "nT" );   // Bmin
    fprintf( fp, " %12s", "nT" );
    fprintf( fp, " %12s", "nT" );
    fprintf( fp, " %12s", "nT" );

    fprintf( fp, " %12s", "nT" );   // M_Used
    fprintf( fp, " %12s", "nT" );   // M_Ref
    fprintf( fp, " %12s", "nT" );   // M_IGRF

    fprintf(fp, "    ");
    for (i=0; i<m->nAlpha; i++) { sprintf( Str, "Dimless" ); fprintf(fp, " %12s", Str ); }
    fprintf(fp, "    ");
    for (i=0; i<m->nAlpha; i++) { sprintf( Str, "Dimless" ); fprintf(fp, " %12s", Str ); }
    fprintf(fp, "    ");
    for (i=0; i<m->nAlpha; i++) { sprintf( Str, "nT" ); fprintf(fp, " %12s", Str ); }
    fprintf(fp, "    ");
    for (i=0; i<m->nAlpha; i++) { sprintf( Str, "Re" ); fprintf(fp, " %12s", Str ); }
    fprintf(fp, "\n");

    Lgm_free_ctrans(c);

}


void WriteMagEphemData( FILE *fp, char *IntModel, char *ExtModel, double Kp, double Dst, Lgm_MagEphemInfo *m ){

    int             i;
    char            Str[128];
    double          GeodLat, GeodLong, GeodHeight, L;
    double          Bsc_mag, Bfn_mag, Bfs_mag, Bmin_mag, Alpha_Loss_Cone_n, Alpha_Loss_Cone_s;
    Lgm_DateTime    DT_UTC;
    Lgm_CTrans      *c = Lgm_init_ctrans(0);
    Lgm_Vector      v, Bsc, Bfn, Bfn_geo, Bfs, Bfs_geo, Bmin;

    Lgm_Set_Coord_Transforms( m->Date, m->UTC, c );
    Lgm_Make_UTC( m->Date, m->UTC, &DT_UTC, c );
    Lgm_DateTimeToString( Str, DT_UTC, 0, 4 );
    

    fprintf( fp, "%25s",     Str );          // Date+Time in ISO 8601 format
    fprintf( fp, "   %10ld", c->UTC.Date );  // Date
    fprintf( fp, " %5d",     c->UTC.Doy );   // DOY
    fprintf( fp, " %13.8lf", c->UTC.Time );  // UTC
    fprintf( fp, " %16.8lf", c->UTC.JD );    // Julian Date
    fprintf( fp, " %15.3lf", Lgm_UTC_to_GpsSeconds( &c->UTC, c ) ); // GpsTime


    Lgm_Convert_Coords( &m->P, &v, GSM_TO_GEO, c );
    fprintf( fp, " %12g", v.x );     // Xgeo
    fprintf( fp, " %12g", v.y );     // Ygeo
    fprintf( fp, " %12g", v.z );     // Zgeo

    fprintf( fp, " %12g", m->P.x );  // Xgsm
    fprintf( fp, " %12g", m->P.y );  // Ygsm
    fprintf( fp, " %12g", m->P.z );  // Zgsm

    Lgm_Convert_Coords( &m->P, &v, GSM_TO_SM, c );
    fprintf( fp, " %12g", v.x );        // Xsm
    fprintf( fp, " %12g", v.y );        // Ysm
    fprintf( fp, " %12g", v.z );        // Zsm

    Lgm_Convert_Coords( &m->P, &v, GSM_TO_GEI2000, c );
    fprintf( fp, " %12g", v.x );        // Xgei
    fprintf( fp, " %12g", v.y );        // Ygei
    fprintf( fp, " %12g", v.z );        // Zgei

    Lgm_Convert_Coords( &m->P, &v, GSM_TO_GSE, c );
    fprintf( fp, " %12g", v.x );        // Xgse
    fprintf( fp, " %12g", v.y );        // Ygse
    fprintf( fp, " %12g", v.z );        // Zgse

    if ( !strcmp( ExtModel, "IGRF" ) || !strcmp( ExtModel, "CDIP" ) || !strcmp( ExtModel, "EDIP" ) ) {
        // If our "external model is just a dipole or igrf, then "internal doesnt really mean anything...)
        fprintf( fp, " %14s", "N/A" );       // Int model is Not applicable
    } else {
        fprintf( fp, " %14s", IntModel );    // Int Model
    }
    fprintf( fp, " %14s", ExtModel );        // Ext Model
    fprintf( fp, " %7.1f",  Kp );            // Kp
    fprintf( fp, " %8d",  (int)Dst );             // Dst
    
    m->LstarInfo->mInfo->Bfield( &m->P, &Bsc, m->LstarInfo->mInfo );
    fprintf( fp, " %12g", Bsc.x );  // Bsc_x_gsm
    fprintf( fp, " %12g", Bsc.y );  // Bsc_y_gsm
    fprintf( fp, " %12g", Bsc.z );  // Bsc_z_gsm
    fprintf( fp, " %12g", (Bsc_mag = Lgm_Magnitude( &Bsc )) );    // |B|

    switch ( m->FieldLineType ) {
        case LGM_OPEN_IMF:
                            fprintf( fp, " %29s",  "LGM_OPEN_IMF" ); // FL Type
                            break;
        case LGM_CLOSED:
                            fprintf( fp, " %29s",  "LGM_CLOSED" ); // FL Type
                            break;
        case LGM_OPEN_N_LOBE:
                            fprintf( fp, " %29s",  "LGM_OPEN_N_LOBE" ); // FL Type
                            break;
        case LGM_OPEN_S_LOBE:
                            fprintf( fp, " %29s",  "LGM_OPEN_S_LOBE" ); // FL Type
                            break;
        case LGM_INSIDE_EARTH:
                            fprintf( fp, " %29s",  "LGM_INSIDE_EARTH" ); // FL Type
                            break;
        case LGM_TARGET_HEIGHT_UNREACHABLE:
                            fprintf( fp, " %29s",  "LGM_TARGET_HEIGHT_UNREACHABLE" ); // FL Type
                            break;
        default:
                            fprintf( fp, " %29s",  "UNKNOWN FIELD TYPE" ); // FL Type
                            break;
    }

    if ( (m->FieldLineType == LGM_CLOSED) || (m->FieldLineType == LGM_OPEN_N_LOBE) ) {

        fprintf( fp, " %12g", m->Ellipsoid_Footprint_Pn.x );     // Xgsm   North Foot
        fprintf( fp, " %12g", m->Ellipsoid_Footprint_Pn.y );     // Ygsm
        fprintf( fp, " %12g", m->Ellipsoid_Footprint_Pn.z );     // Zgsm

        Lgm_Convert_Coords( &m->Ellipsoid_Footprint_Pn, &v, GSM_TO_GEO, c );
        fprintf( fp, " %12g", v.x );                    // Xgeo   North Foot
        fprintf( fp, " %12g", v.y );                    // Ygeo
        fprintf( fp, " %12g", v.z );                    // Zgeo

        Lgm_WGS84_to_GEOD( &v, &GeodLat, &GeodLong, &GeodHeight );
        fprintf( fp, " %12g", GeodLat );                // Geod Lat   North Foot
        fprintf( fp, " %12g", GeodLong );               // Geod Long
        fprintf( fp, " %12g", GeodHeight );             // Geod Height

        m->LstarInfo->mInfo->Bfield( &m->Ellipsoid_Footprint_Pn, &Bfn, m->LstarInfo->mInfo );
        Lgm_Convert_Coords( &Bfn, &Bfn_geo, GSM_TO_WGS84, c );
        fprintf( fp, " %12g", Bfn_geo.x );                                  // Bfn_x_geo
        fprintf( fp, " %12g", Bfn_geo.y );                                  // Bfn_y_geo
        fprintf( fp, " %12g", Bfn_geo.z );                                  // Bfn_z_geo
        fprintf( fp, " %12g", (Bfn_mag = Lgm_Magnitude( &Bfn_geo )) );      // |B|

        fprintf( fp, " %12g", Bfn.x );                                  // Bfn_x_gsm
        fprintf( fp, " %12g", Bfn.y );                                  // Bfn_y_gsm
        fprintf( fp, " %12g", Bfn.z );                                  // Bfn_z_gsm
        fprintf( fp, " %12g", (Bfn_mag = Lgm_Magnitude( &Bfn )) );      // |B|

        Alpha_Loss_Cone_n = asin( sqrt( Bsc_mag/Bfn_mag ) )*DegPerRad;
        fprintf( fp, " %12g", Alpha_Loss_Cone_n );                      // Northern Loss Cone Angle


    } else {

        fprintf( fp, " %12g", LGM_FILL_VALUE ); 
        fprintf( fp, " %12g", LGM_FILL_VALUE ); 
        fprintf( fp, " %12g", LGM_FILL_VALUE );

        fprintf( fp, " %12g", LGM_FILL_VALUE );
        fprintf( fp, " %12g", LGM_FILL_VALUE );
        fprintf( fp, " %12g", LGM_FILL_VALUE );

        fprintf( fp, " %12g", LGM_FILL_VALUE );
        fprintf( fp, " %12g", LGM_FILL_VALUE );
        fprintf( fp, " %12g", LGM_FILL_VALUE );

        fprintf( fp, " %12g", LGM_FILL_VALUE );
        fprintf( fp, " %12g", LGM_FILL_VALUE );
        fprintf( fp, " %12g", LGM_FILL_VALUE );
        fprintf( fp, " %12g", LGM_FILL_VALUE );

        fprintf( fp, " %12g", LGM_FILL_VALUE );

    }


    if ( (m->FieldLineType == LGM_CLOSED) || (m->FieldLineType == LGM_OPEN_S_LOBE) ) {

        fprintf( fp, " %12g", m->Ellipsoid_Footprint_Ps.x );     // Xgsm   South Foot
        fprintf( fp, " %12g", m->Ellipsoid_Footprint_Ps.y );     // Ygsm
        fprintf( fp, " %12g", m->Ellipsoid_Footprint_Ps.z );     // Zgsm

        Lgm_Convert_Coords( &m->Ellipsoid_Footprint_Ps, &v, GSM_TO_GEO, c );
        fprintf( fp, " %12g", v.x );                    // Xgeo   South Foot
        fprintf( fp, " %12g", v.y );                    // Ygeo
        fprintf( fp, " %12g", v.z );                    // Zgeo


        Lgm_WGS84_to_GEOD( &v, &GeodLat, &GeodLong, &GeodHeight );
        fprintf( fp, " %12g", GeodLat );                // Geod Lat   South Foot
        fprintf( fp, " %12g", GeodLong );               // Geod Long
        fprintf( fp, " %12g", GeodHeight );             // Geod Height

        m->LstarInfo->mInfo->Bfield( &m->Ellipsoid_Footprint_Ps, &Bfs, m->LstarInfo->mInfo );
        Lgm_Convert_Coords( &Bfs, &Bfs_geo, GSM_TO_WGS84, c );
        fprintf( fp, " %12g", Bfs_geo.x );                                  // Bfs_x_geo
        fprintf( fp, " %12g", Bfs_geo.y );                                  // Bfs_y_geo
        fprintf( fp, " %12g", Bfs_geo.z );                                  // Bfs_z_geo
        fprintf( fp, " %12g", (Bfs_mag = Lgm_Magnitude( &Bfs_geo )) );      // |B|

        fprintf( fp, " %12g", Bfs.x );                                  // Bfs_x_gsm
        fprintf( fp, " %12g", Bfs.y );                                  // Bfs_y_gsm
        fprintf( fp, " %12g", Bfs.z );                                  // Bfs_z_gsm
        fprintf( fp, " %12g", (Bfs_mag = Lgm_Magnitude( &Bfs )) );      // |B|

        Alpha_Loss_Cone_s = asin( sqrt( Bsc_mag/Bfs_mag ) )*DegPerRad;
        fprintf( fp, " %12g", Alpha_Loss_Cone_s );                      // Southern Loss Cone Angle


    } else {
        fprintf( fp, " %12g", LGM_FILL_VALUE ); 
        fprintf( fp, " %12g", LGM_FILL_VALUE ); 
        fprintf( fp, " %12g", LGM_FILL_VALUE );

        fprintf( fp, " %12g", LGM_FILL_VALUE );
        fprintf( fp, " %12g", LGM_FILL_VALUE );
        fprintf( fp, " %12g", LGM_FILL_VALUE );

        fprintf( fp, " %12g", LGM_FILL_VALUE );
        fprintf( fp, " %12g", LGM_FILL_VALUE );
        fprintf( fp, " %12g", LGM_FILL_VALUE );

        fprintf( fp, " %12g", LGM_FILL_VALUE );
        fprintf( fp, " %12g", LGM_FILL_VALUE );
        fprintf( fp, " %12g", LGM_FILL_VALUE );
        fprintf( fp, " %12g", LGM_FILL_VALUE );

        fprintf( fp, " %12g", LGM_FILL_VALUE );

    }

    

    if ( m->FieldLineType == LGM_CLOSED ) {
        fprintf( fp, " %12g", m->Pmin.x );          // Xgsm  Pmin
        fprintf( fp, " %12g", m->Pmin.y );          // Ygsm
        fprintf( fp, " %12g", m->Pmin.z );          // Zgsm

        m->LstarInfo->mInfo->Bfield( &m->Pmin, &Bmin, m->LstarInfo->mInfo );
        fprintf( fp, " %12g", Bmin.x );  // Bmin_x_gsm
        fprintf( fp, " %12g", Bmin.y );  // Bmin_y_gsm
        fprintf( fp, " %12g", Bmin.z );  // Bmin_z_gsm
        fprintf( fp, " %12g", (Bmin_mag = Lgm_Magnitude( &Bmin )) );    // |B|
    } else {
        fprintf( fp, " %12g", LGM_FILL_VALUE );
        fprintf( fp, " %12g", LGM_FILL_VALUE );
        fprintf( fp, " %12g", LGM_FILL_VALUE );

        fprintf( fp, " %12g", LGM_FILL_VALUE );
        fprintf( fp, " %12g", LGM_FILL_VALUE );
        fprintf( fp, " %12g", LGM_FILL_VALUE );
        fprintf( fp, " %12g", LGM_FILL_VALUE );
    }

    
    fprintf( fp, " %12g", m->Mused );   // M_Used
    fprintf( fp, " %12g", m->Mref );    // M_Ref
    fprintf( fp, " %12g", m->Mcurr );   // M_IGRF

    // L*'s
    fprintf(fp, "    ");
    for (i=0; i<m->nAlpha; i++) { fprintf(fp, " %12g", m->Lstar[i] ); }

    // McIlwain L (computed from I, Bm, M)
    fprintf(fp, "    ");
    for (i=0; i<m->nAlpha; i++) { 
        L = ( m->I[i] > 0.0 ) ? LFromIBmM_McIlwain(m->I[i], m->Bm[i], m->Mused ) : 0.0;
        fprintf(fp, " %12g", L);
    }

    // Bms's
    fprintf(fp, "    ");
    for (i=0; i<m->nAlpha; i++) { fprintf(fp, " %12g", m->Bm[i] ); }

    // I's
    fprintf(fp, "    ");
    for (i=0; i<m->nAlpha; i++) { fprintf(fp, " %12g", m->I[i] ); }


    fprintf(fp, "\n");

    
    Lgm_free_ctrans( c );
}
