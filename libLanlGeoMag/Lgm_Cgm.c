#include <stdio.h>
#include <math.h>
#include "Lgm/Lgm_QuadPack.h"
#include "Lgm/Lgm_MagModelInfo.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <string.h>

#define TRACE_TOL   1e-7


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


int  Lgm_GEO_TO_CGM( double geoLat, double geoLon, double geoAlt, double *CgmLat, 
                     double *CgmLon, double *CgmRad, Lgm_MagModelInfo *m ) {

    int         res, TraceFlag;
    double      r, R, L, Lat, SinLat, SinLat2, CosLat2;
    Lgm_Vector  wgsVec, gsmVec, vSME, v3_cdmag, Bvec;


    /*
     * Trace the field line for the given position. The Min-B point is at v3.
     */
    // First convert geo triad to GSM (x,y,z)
    Lgm_GEOD_to_WGS84( geoLat, geoLon, geoAlt, &wgsVec);
    Lgm_Convert_Coords( &wgsVec, &gsmVec, WGS84_TO_GSM, m->c );
    //Trace to SM equator
    if ( Lgm_TraceToSMEquat( &gsmVec, &vSME, TRACE_TOL, m)  < 1) return( -1 );
    
    //Now vSME has position of SM equator (in GSM)
    //So convert to CDMAG
    Lgm_Convert_Coords( &vSME, &v3_cdmag, GSM_TO_CDMAG, m->c );

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
     * CgmLon should be the same as the longitude of v3_cdmag.
     */
    *CgmLon = DegPerRad*atan2( v3_cdmag.y, v3_cdmag.x );
    if ( *CgmLon < 0.0 ) *CgmLon += 360.0;

    *CgmRad = geoAlt; 
    *CgmLat = DegPerRad*acos( sqrt(geoAlt/L) );

    return( 1 );

}

int  Lgm_CGM_TO_GEO( double CgmLat, double CgmLon, double CgmAlt, double *geoLat, 
                     double *geoLon, double *geoAlt, Lgm_MagModelInfo *m ) {

    int         res, TraceFlag;
    double      r, R, L, Lat, SinLat, SinLat2, CosLat2;
    Lgm_Vector  wgsVec, gsmVec, vSME, v3_cdmag, Bvec;

    //First get R, L, etc.
    R = CgmAlt+WGS84_A;
    SinLat = v3_cdmag.z/R;
    SinLat2 = SinLat*SinLat;
    CosLat2 = 1.0 - SinLat2;
    L = R/CosLat2;

    //construct location of Pmin (for dipole FL)
    Lgm_SphToCartCoords( 0.0 , CgmLon, L, &v3_cdmag ); //Is this right??
    Lgm_Convert_Coords( &v3_cdmag, &vSME, CDMAG_TO_GSM, m->c );
    
    //Now trace from vSME to requested altitude
    Lgm_TraceToEarth( &vSME, &gsmVec, CgmAlt, 1, TRACE_TOL, m );
    
    //Get coords of trace point in WGS84
    Lgm_Convert_Coords( &gsmVec, &wgsVec, GSM_TO_WGS84, m->c );


}

