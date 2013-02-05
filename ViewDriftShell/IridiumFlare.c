#include <Lgm/Lgm_CTrans.h>
#include "Objects.h"
#include <Lgm/Lgm_Vec.h>
#include <Lgm/Lgm_Quat.h>
#include <Lgm/Lgm_Sgp.h>
#include <stdio.h>

/*
 * Compute the ground intersection points where the reflected sunlight will hit
 * Flag1, Flag2, and Flag3 indicate if the reflected ray actually hits the Earth
 */
void IridiumFlare( double JD, _SgpTLE *TLE, Lgm_Vector *EarthToSun, int *Flag1, int *Flag2, int *Flag3, Lgm_Vector *g1, Lgm_Vector *g2, Lgm_Vector *g3, Lgm_Vector *P ) {

    double          RotAngle, Q[4], f, tsince, tMin, tMax, tMinFwd, UTC;
    int             Year, Month, Day;
    long int        Date;
    Lgm_Vector      v, zenith, RotAxis, n1, n2, n3;
    Lgm_Vector      sun, vv, r1, r2, r3, Pteme, Vteme;
    RayType         SatToSunRay, Ray;
    EllipsoidType   Earth;
    _SgpInfo        sgp;
    Lgm_CTrans      *c = Lgm_init_ctrans( 0 );

    Lgm_JD_to_Date(JD, &Year, &Month, &Day, &UTC );
    Date = Year*10000+Month*100+Day;

    /*
     * Define Earth
     */
    Earth.Origin.x = Earth.Origin.y = Earth.Origin.z = 0.0;
    Earth.Radius_a  = 1.0;
    Earth.Radius2_a = Earth.Radius_a*Earth.Radius_a;
    Earth.Radius_b  = WGS84_B/WGS84_A;
    Earth.Radius2_b = Earth.Radius_b*Earth.Radius_b;

    /*
     * Get sat position in MOD coords (Sun will be computed in MOD)
     */
    tsince = (JD - TLE->JD)*1440.0;
    LgmSgp_SGP4_Init( &sgp, TLE );
    LgmSgp_SGP4( tsince, &sgp ); Pteme.x = sgp.X/SGP_XKMPER; Pteme.y = sgp.Y/SGP_XKMPER; Pteme.z = sgp.Z/SGP_XKMPER;
    Lgm_Set_Coord_Transforms( Date, UTC, c ); Lgm_Convert_Coords( &Pteme, P, TEME_TO_TOD, c );


    /*
     * construct ray from sat to sun
     */
    sun.x = EarthToSun->x - P->x; sun.y = EarthToSun->y - P->y; sun.z = EarthToSun->z - P->z;
    Lgm_NormalizeVector( &sun );
    SatToSunRay.Origin.x = P->x; SatToSunRay.Origin.y = P->y; SatToSunRay.Origin.z = P->z;
    SatToSunRay.Direction.x = sun.x; SatToSunRay.Direction.y = sun.y; SatToSunRay.Direction.z = sun.z;

    /*
     * Test for intersection with Earth
     */
    *Flag1 = *Flag2 = *Flag3 = FALSE;
    if ( !EllipsoidIntersect( &Earth, &SatToSunRay, &tMin, &tMax, &tMinFwd ) ) {

        // There is a clear path from the sat to the sun


        /*
         * Set zenith vector (vector that points from center of earth to sat)
         */
        zenith = *P;
        Lgm_NormalizeVector( &zenith );


        /*
         * Compute trajectory of the sat at time JD
         */
//tsince = (JD - TLE->JD)*1440.0+1.0/60.0;
//LgmSgp_SGP4( tsince, &sgp ); e.x = sgp.X/SGP_XKMPER; e.y = sgp.Y/SGP_XKMPER; e.z = sgp.Z/SGP_XKMPER;
//tsince = (JD - TLE->JD)*1440.0-1.0/60.0;
//LgmSgp_SGP4( tsince, &sgp ); s.x = sgp.X/SGP_XKMPER; s.y = sgp.Y/SGP_XKMPER; s.z = sgp.Z/SGP_XKMPER;

//Vteme.x = e.x-s.x; Vteme.y = e.y-s.y; Vteme.z = e.z-s.z;

        Vteme.x = sgp.VX/SGP_XKMPER; Vteme.y = sgp.VY/SGP_XKMPER; Vteme.z = sgp.VZ/SGP_XKMPER;
        Lgm_Convert_Coords( &Vteme, &v, TEME_TO_TOD, c );
        Lgm_NormalizeVector( &v );


        /*
         * Get vector thats perp to both of these vectors
         */
        Lgm_CrossProduct( &zenith, &v, &RotAxis );
        Lgm_NormalizeVector( &RotAxis );

        /*
         * Find the normal to the leading mirror (Main Mirror Array or MMA)
         * Start the vector pointing to zenith. Then rotate by 40 deg. about RotAxis1
         */
        RotAngle = 130.0;
        Lgm_AxisAngleToQuat( &RotAxis, RotAngle, Q );
        Lgm_QuatRotateVector( Q, &zenith, &n1 );


        /*
         * Now find the other two mirror normals. Each are found by rotating n1 by +/-120deg 
         * around zenith.
         */
        RotAngle = 120.0;
        Lgm_AxisAngleToQuat( &zenith, RotAngle, Q );
        Lgm_QuatRotateVector( Q, &n1, &n2 );

        Lgm_AxisAngleToQuat( &zenith, -RotAngle, Q );
        Lgm_QuatRotateVector( Q, &n1, &n3 );



        /*
         * compute reflected vectors -- one from each mirror
         * and compute intersection of each reflected ray with the ground
         */
        Ray.Origin.x = P->x; Ray.Origin.y = P->y; Ray.Origin.z = P->z;
        g1->x = g1->y = g1->z = 0.0;
        g2->x = g2->y = g2->z = 0.0;
        g3->x = g3->y = g3->z = 0.0;


        // n1
        f = 2.0*Lgm_DotProduct( &sun, &n1 );
        if ( f > 0.0 ) {
            vv = n1; Lgm_ScaleVector( &vv, f );
            r1.x = vv.x - sun.x; r1.y = vv.y - sun.y; r1.z = vv.z - sun.z;
            Lgm_NormalizeVector( &r1 );

            Ray.Direction.x = r1.x; Ray.Direction.y = r1.y; Ray.Direction.z = r1.z;
            if ( EllipsoidIntersect( &Earth, &Ray, &tMin, &tMax, &tMinFwd ) ){
                g1->x = Ray.Origin.x + tMinFwd*Ray.Direction.x;
                g1->y = Ray.Origin.y + tMinFwd*Ray.Direction.y;
                g1->z = Ray.Origin.z + tMinFwd*Ray.Direction.z;
                *Flag1 = TRUE;
            }
        }




        // n2
        f = 2.0*Lgm_DotProduct( &sun, &n2 );
        if ( f > 0.0 ) {
            vv = n2; Lgm_ScaleVector( &vv, f );
            r2.x = vv.x - sun.x; r2.y = vv.y - sun.y; r2.z = vv.z - sun.z;
            Lgm_NormalizeVector( &r2 );

            Ray.Direction.x = r2.x; Ray.Direction.y = r2.y; Ray.Direction.z = r2.z;
            if ( EllipsoidIntersect( &Earth, &Ray, &tMin, &tMax, &tMinFwd ) ){
                g2->x = Ray.Origin.x + tMinFwd*Ray.Direction.x;
                g2->y = Ray.Origin.y + tMinFwd*Ray.Direction.y;
                g2->z = Ray.Origin.z + tMinFwd*Ray.Direction.z;
                *Flag2 = TRUE;
            }
        }



        f = 2.0*Lgm_DotProduct( &sun, &n3 );
        if ( f > 0.0 ) {
            vv = n3; Lgm_ScaleVector( &vv, f );
            r3.x = vv.x - sun.x; r3.y = vv.y - sun.y; r3.z = vv.z - sun.z;
            Lgm_NormalizeVector( &r3 );

            Ray.Direction.x = r3.x; Ray.Direction.y = r3.y; Ray.Direction.z = r3.z;
            if ( EllipsoidIntersect( &Earth, &Ray, &tMin, &tMax, &tMinFwd ) ){
                g3->x = Ray.Origin.x + tMinFwd*Ray.Direction.x;
                g3->y = Ray.Origin.y + tMinFwd*Ray.Direction.y;
                g3->z = Ray.Origin.z + tMinFwd*Ray.Direction.z;
                *Flag3 = TRUE;
            }
        }



    }

    free( c );

}



