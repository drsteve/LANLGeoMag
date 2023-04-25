// Geodetic.c
#include <Lgm_CTrans.h>
int main(){

    double      GeodLat, GeodLong, GeodHeight, r;
    double      GeoLat, GeoLong, GeoHeight;
    Lgm_Vector  v;

    GeodLat    = 35.0 + 53.0/60.0 + 17.0/3600.0;   // degrees
    GeodLong   = -106.0 - 18.0/60.0 - 23.0/3600.0; // degrees
    GeodHeight = 20.0; // km
GeoLat    = 56.02;   // degrees
GeoLong   = -141.42; // degrees
GeoHeight = 850.0; // km
v.x = (1.0+GeoHeight/WGS84_A)*cos( GeoLong * RadPerDeg )*cos( GeoLat*RadPerDeg );
v.y = (1.0+GeoHeight/WGS84_A)*sin( GeoLong * RadPerDeg )*cos( GeoLat*RadPerDeg );
v.z = (1.0+GeoHeight/WGS84_A)*sin( GeoLat*RadPerDeg );
    Lgm_WGS84_to_GEOD( &v, &GeodLat, &GeodLong, &GeodHeight );

//    Lgm_GEOD_to_WGS84( GeodLat, GeodLong, GeodHeight, &v );

    printf( "            Geodetic Position\n");
    printf( "   ===================================\n");
    printf( "      GeodLat    : %lf degrees\n", GeodLat );
    printf( "      GeodLong   : %lf degrees\n", GeodLong );
    printf( "      GeodHeight : %lf km\n\n", GeodHeight );

    printf( "            Geocentric Position\n");
    printf( "   ===================================\n");
    printf( "      X         : %-15lf meters\n", v.x*WGS84_A*1000.0);   // meters
    printf( "      Y         : %-15lf meters\n", v.y*WGS84_A*1000.0);   // meters
    printf( "      Z         : %-15lf meters\n", v.z*WGS84_A*1000.0); // meters
    r = sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
    printf( "      Radius    : %lf Re (%lf km)\n", r, r*WGS84_A );
    printf( "      Latitude  : %lf degrees\n", DegPerRad*asin( v.z/r ) );
    printf( "      Longitude : %lf degrees\n\n", DegPerRad*atan2( v.y, v.x ) );

    exit(0);
}
