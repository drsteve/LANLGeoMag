#include "Lgm/Lgm_MagModelInfo.h"

/*
 *  Compute the radial component of B. MLT/mlat assumed to be in SM coords.
 */
double Lgm_Bradial( double MLT, double mlat, Lgm_MagModelInfo *m ) {

    double      Phi, Lat, r, CosLat, SinLat, CosPhi, SinPhi, Br;
    Lgm_Vector  v1, u, B;

    // Get radius -- assume 100km above a spherical Earth
    r = 1.0 + 100.0/Re;

    // get azimuth angle in SM coords.
    Phi = (MLT*15.0 - 180.0)*RadPerDeg;

    // get latitude
    Lat = mlat*RadPerDeg;


    // Vector to point in question
    CosLat = cos( Lat ); SinLat = sin( Lat );
    CosPhi = cos( Phi ); SinPhi = sin( Phi );
    v1.x = r*CosLat*CosPhi;
    v1.y = r*CosLat*SinPhi;
    v1.z = r*SinLat;

    // Convert v1 from SM to GSM.
    Lgm_Convert_Coords( &v1, &u, SM_TO_GSM, m->c );

    // Get B-field vector at the point
    m->Bfield( &u, &B, m );


    // compute radial component of field
    Br = CosLat*CosPhi*B.x + CosLat*SinPhi*B.y + SinLat*B.z;



    //return( (Br<0.0) ? -Ang : Ang );
    return( Br );


}


/*
 *    Find magnetic equator and last closed FL in SM MLT/mlat coords.
 */
int Lgm_FindDipEquator( double MLT, double *mlat, Lgm_MagModelInfo *m ) {

    int     done;
    double  a, b, c, Da, Db, Dc;


    /*
     * Use bissection to find the DIP equator. I.e. where Bradial is zero.
     */
    a = -30.0; Da = Lgm_Bradial( MLT, a, m );
    c =  30.0; Dc = Lgm_Bradial( MLT, c, m );
    if ( Da*Db > 0.0 ) return(-1); // no bracket...

    done = FALSE;
    while ( !done ) {

        if ( fabs(c-a) < 1e-8 ) {

            // Our bracket is pretty small and we still havent converged. Return what we have...
            *mlat = b;
            return(-1);

        } else {

            b  = 0.5*(a+c); // take midpoint
            Db = Lgm_Bradial( MLT, b, m );

            if ( fabs(Db) < 1e-6 ) {
                done = TRUE;
                *mlat = b;
            } else if ( Db < 0.0 ) {
                c  = b;
                Dc = Db;
            } else {
                a  = b;
                Da = Db;
            }
        }

    }

    return(1);

}
