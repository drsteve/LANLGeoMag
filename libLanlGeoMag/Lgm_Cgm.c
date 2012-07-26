#include <stdio.h>
#include <math.h>
#include "Lgm/Lgm_QuadPack.h"
#include "Lgm/Lgm_MagModelInfo.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <string.h>

#define TRACE_TOL   1e-8


/*
 *
 * Input:
 *
 *
 * Output:
 *
 *
 */
int  Lgm_GSM_TO_CBM( Lgm_Vector *u, double *CgmLat, double *CgmLon, double *CgmRad, Lgm_MagModelInfo *m ) {

    int         TraceFlag;
    double      r, R, L, Lat, SinLat, SinLat2, CosLat2;
    Lgm_Vector  v1, v2, v3, v3_cdmag, Bvec;


    /*
     * Trace the field line for the given position. The Min-B point is at v3.
     */
    r = Lgm_Magnitude( u );
    TraceFlag = Lgm_Trace( u, &v1, &v2, &v3, 100.0, TRACE_TOL, TRACE_TOL, m );


    /*
     * Get point to trace back from.
     */
    if ( TraceFlag == LGM_CLOSED ) {

        Lgm_Convert_Coords( &v3, &v3_cdmag, GSM_TO_CDMAG, m->c );

    } else if ( TraceFlag == LGM_OPEN_N_LOBE ) {

        // should really trace to a set R
        Lgm_Convert_Coords( &v2, &v3_cdmag, GSM_TO_CDMAG, m->c );

    } else if ( TraceFlag == LGM_OPEN_S_LOBE ) {

        // should really trace to a set R
        Lgm_Convert_Coords( &v1, &v3_cdmag, GSM_TO_CDMAG, m->c );

    }


    /*
     * The Corrected Geomagnetic Coordinates of the original point are obtained
     * by now tracing back -- with a dipole field -- from v3 to the same radius
     * as we started from.  However, v3_cdmag is in the cdmag coord system,
     * where the equation of a dipole FL there is trivial: R = L cos^2(Lat).
     * 
     * First, we can compute L bsed on the known point v3_cdmag.
     */
    R   = Lgm_Magnitude( &v3_cdmag );
    SinLat = v3_cdmag.z/R;
    SinLat2 = SinLat*SinLat;
    CosLat2 = 1.0 - SinLat2;
    L = R/CosLat2;

    /*
     * Now using the L value and the original r, get CgmLat.
     */
    *CgmRad = r;
    *CgmLat = DegPerRad*acos( sqrt(r/L) );


    /*
     * CgmLon should be the same as the longitude of v3_cdmag.
     */
    *CgmLon = DegPerRad*atan2( v3_cdmag.y, v3_cdmag.x );
    if ( *CgmLon < 0.0 ) *CgmLon += 360.0;
    
    
    return( TraceFlag );

}


int  Lgm_GEOD_TO_CGM( double geoLat, double geoLon, double geoAlt, double *CgmLat, 
                     double *CgmLon, double *CgmRad, Lgm_MagModelInfo *m ) {

    int         res, TraceFlag;
    double      r, R, L, Lat, SinLat, SinLat2, CosLat2;
    Lgm_Vector  wgsVec, gsmVec, vSME, v3_cdmag, Bvec;


    /*
     * Trace the field line for the given position. The Min-B point is at v3.
     */
    // First convert geo triad to GSM (x,y,z)
//    Lgm_GEOD_to_WGS84( geoLat, geoLon, geoAlt, &wgsVec);
Lgm_SphToCartCoords( geoLat , geoLon, geoAlt, &wgsVec );
    Lgm_Convert_Coords( &wgsVec, &gsmVec, WGS84_TO_GSM, m->c );
    //Trace to SM equator
    if ( Lgm_TraceToSMEquat( &gsmVec, &vSME, TRACE_TOL, m)  < 1) return( -1 );
    
    //Now vSME has position of SM equator (in GSM)
    //So convert to CDMAG
    Lgm_Convert_Coords( &vSME, &v3_cdmag, GSM_TO_CDMAG, m->c );
    //Lgm_Convert_Coords( &vSME, &v3_cdmag, GSM_TO_SM, m->c );


    /*
    * The Corrected Geomagnetic Coordinates of the original point are obtained
    * by now tracing back -- with a dipole field -- from v3 to the same radius
    * as we started from.  However, v3_cdmag is in the cdmag coord system,
    * where the equation of a dipole FL there is trivial: R = L cos^2(Lat).
    * 
    * First, we can compute L bsed on the known point v3_cdmag.
    */

printf("v3_cdmag-> = %.15lf %.15lf %.15lf %.15lf\n", v3_cdmag.x, v3_cdmag.y, v3_cdmag.z, Lgm_Magnitude( &v3_cdmag)  );
    R   = Lgm_Magnitude( &v3_cdmag );

    SinLat = v3_cdmag.z/R;
    SinLat2 = SinLat*SinLat;
    CosLat2 = 1.0 - SinLat2;
    L = R/CosLat2;
printf("L-> = %.15lf\n", L);

    /*
     * CgmLon should be the same as the longitude of v3_cdmag.
     */
    *CgmLon = DegPerRad*atan2( v3_cdmag.y, v3_cdmag.x );
    if ( *CgmLon < 0.0 ) *CgmLon += 360.0;

    *CgmRad = geoAlt; 
    *CgmLat = DegPerRad*acos( sqrt(geoAlt/L) );

    return( 1 );

}

int  Lgm_CGM_TO_GEOD( double CgmLat, double CgmLon, double CgmRad, double *geoLat, 
                     double *geoLon, double *geoAlt, Lgm_MagModelInfo *m ) {

    int         res, TraceFlag;
    double      r, R, L, Lat, SinLat, SinLat2, CosLat2, theta, phi;
    Lgm_Vector  wgsVec, gsmVec, vSME, v3_cdmag, Bvec;
//printf("CGM  lat %g, lon %g, alt %g\n", CgmLat, CgmLon, CgmAlt );

    //First get R, L, etc.
//    R = 1.0+CgmAlt/WGS84_A;
R = CgmRad;
    SinLat = sin(RadPerDeg*CgmLat);
    SinLat2 = SinLat*SinLat;
    CosLat2 = 1.0 - SinLat2;
    L = R/CosLat2;
printf("L<- = %.15lf\n", L);
//printf("R: %g, L: %g\n", R, L );
    //construct location of Pmin (for dipole FL)
//    Lgm_SphToCartCoords( 0.0 , CgmLon, L, &v3_cdmag); //Is this right??
    Lgm_SphToCartCoords( 0.0 , CgmLon, L, &v3_cdmag); //Is this right??
    Lgm_Convert_Coords( &v3_cdmag, &vSME, CDMAG_TO_GSM, m->c );

printf("v3_cdmag<- = %.15lf %.15lf %.15lf     %.15lf\n", v3_cdmag.x, v3_cdmag.y, v3_cdmag.z, Lgm_Magnitude( &v3_cdmag)  );

    
    //Now trace from vSME to requested altitude
//Lgm_TraceToEarth( &vSME, &gsmVec, CgmAlt, 1, TRACE_TOL, m );
    Lgm_TraceToSphericalEarth( &vSME, &gsmVec, (CgmRad-1.0)*WGS84_A, 1, TRACE_TOL, m );
Lgm_Vector vSME2;
    if ( Lgm_TraceToSMEquat( &gsmVec, &vSME2, TRACE_TOL, m)  < 1) return( -1 );
printf("vSME  = %g %g %g\n", vSME.x, vSME.y, vSME.z);
printf("vSME2 = %g %g %g\n", vSME2.x, vSME2.y, vSME2.z);
    //Lgm_TraceToSphericalEarth( &vSME, &gsmVec, 1.0+CgmRad/WGS84_A, 1, TRACE_TOL, m );
    //Lgm_TraceToSphericalEarth( &vSME, &gsmVec, 100.0, 1, TRACE_TOL, m );
//printf("GSM  x %g, y %g, z %g\n", gsmVec.x, gsmVec.y, gsmVec.z );
printf("HHHHHH = %g   (CgmRad-1.0)*WGS84_A = %g\n", (Lgm_Magnitude( &gsmVec ) -1.0)*WGS84_A, (CgmRad-1.0)*WGS84_A);
    
    //Get coords of trace point in WGS84
    Lgm_Convert_Coords( &gsmVec, &wgsVec, GSM_TO_WGS84, m->c );
*geoAlt = Lgm_Magnitude( &wgsVec );
    //Lgm_WGS84_to_GEOD( &wgsVec, geoLat, geoLon, geoAlt );
    Lgm_CartToSphCoords( &wgsVec, geoLat, geoLon, geoAlt );
    *geoAlt = (*geoAlt-1.0)*WGS84_A;
//*geoAlt = (*geoAlt-1.0)*WGS84_A;
//printf("Radius WGS: %g\n", Lgm_Magnitude( &wgsVec ) );
//printf("WGS x %g, y %g, z %g\n", wgsVec.x, wgsVec.y, wgsVec.z );
//printf("Lat %g, Lon %g, Alt %g\n", *geoLat, *geoLon, *geoAlt );

    return( 1 );
}

