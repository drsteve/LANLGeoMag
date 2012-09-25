#include "Lgm/Lgm_CTrans.h"
#include "Lgm/Lgm_MagEphemInfo.h"
const char *sMonth[] = { "", "January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December" };


void Lgm_WriteMagEphemHeader( FILE *fp, int SpiceBody,  char *Spacecraft, int IdNumber, char *IntDesig, char *CmdLine, int nAscend, Lgm_DateTime *Ascend_UTC, Lgm_Vector *Ascend_U, int nPerigee, Lgm_DateTime *Perigee_UTC, Lgm_Vector *Perigee_U, int nApogee, Lgm_DateTime *Apogee_UTC, Lgm_Vector *Apogee_U, Lgm_MagEphemInfo *m ){

    int         i, Year, Month, Day, HH, MM, SS, n, tsl, n2;
    char        Str[80], *Str2;
    char        *p, TextStr[80], IsoTimeString[80];
    long int    CreateDate;
    double      JD, UTC, R;
    double      GeodLat, GeodLong, GeodHeight;
    Lgm_CTrans  *c = Lgm_init_ctrans(0);
    Lgm_Vector  w;


    JD = Lgm_GetCurrentJD(c);
    CreateDate = Lgm_JD_to_Date( JD, &Year, &Month, &Day, &UTC );
    Lgm_UT_to_HMS( UTC, &HH, &MM, &SS );

    /*
     * Write Header
     */
    int nCol = 0;

    fprintf( fp, "# {\n");
    if ( m->nAlpha > 0 ) {
        fprintf( fp, "#  \"Alpha\":            { \"DESCRIPTION\": \"Pitch Angles.\",\n");
        fprintf( fp, "#                               \"NAME\": \"Pitch Angle\",\n");
        fprintf( fp, "#                              \"TITLE\": \"Pitch Angle\",\n");
        fprintf( fp, "#                              \"LABEL\": \"Pitch Angle\",\n");
        fprintf( fp, "#                          \"DIMENSION\": [ %d ],\n", m->nAlpha );
        fprintf( fp, "#                             \"VALUES\": [ ");
        for (i=0; i<m->nAlpha-1; i++) fprintf(fp, "%g, ", m->Alpha[i] );
        fprintf(fp, "%g ],\n", m->Alpha[i] ); 

        fprintf( fp, "#                      \"ELEMENT_NAMES\": [ ");
        for (i=0; i<m->nAlpha-1; i++) fprintf(fp, "\"PA%d\", ", i );
        fprintf(fp, "\"PA%d\" ],\n", i ); 

        fprintf( fp, "#                     \"ELEMENT_LABELS\": [ ");
        for (i=0; i<m->nAlpha-1; i++) fprintf(fp, "\"%g (Deg.)\", ", m->Alpha[i] );
        fprintf(fp, "\"%g Deg.\" ],\n", m->Alpha[i] ); 
        fprintf( fp, "#                              \"UNITS\": \"Degrees\",\n");
        fprintf( fp, "#                          \"VALID_MIN\":  0.0,\n");
        fprintf( fp, "#                          \"VALID_MAX\": 90.0,\n");
        fprintf( fp, "#                         \"FILL_VALUE\": -1e31 },\n");
        fprintf( fp, "#\n");
    }


    if ( nAscend > 0 ) {
        fprintf( fp, "#  \"AscendingNodeTimes\":{ \"DESCRIPTION\": \"Ascending Node Crossing Times on this day.\",\n");
        fprintf( fp, "#                                \"NAME\": \"Ascending Node Crossing Times\",\n");
        fprintf( fp, "#                               \"TITLE\": \"Ascending Node Crossing Times\",\n");
        fprintf( fp, "#                               \"LABEL\": \"Ascending Node Crossing Times\",\n");
        fprintf( fp, "#                           \"DIMENSION\": [ %d ],\n", nAscend );

        fprintf( fp, "#                              \"VALUES\": [ ");
        for (i=0; i<nAscend-1; i++) {
            Lgm_DateTimeToString( IsoTimeString, &Ascend_UTC[i], 0, 3 );
            fprintf(fp, "\"%s\", ", IsoTimeString );
        }
        Lgm_DateTimeToString( IsoTimeString, &Ascend_UTC[i], 0, 3 );
        fprintf(fp, "\"%s\" ],\n", IsoTimeString ); 


        fprintf( fp, "#                       \"ELEMENT_NAMES\": [ ");
        for (i=0; i<nAscend-1; i++) fprintf(fp, "\"AscendTime%d\", ", i );
        fprintf(fp, "\"AscendTime%d\" ],\n", i ); 

        fprintf( fp, "#                              \"UNITS\": \"UTC\" },\n");
        fprintf( fp, "#\n");



        fprintf( fp, "#  \"AscendingPosGeod\": { \"DESCRIPTION\": \"Ascending Node Crossing Positions on this day in Geodetic Coords (lat/lon/alt).\",\n");
        fprintf( fp, "#                               \"NAME\": \"Ascending Node Crossing Pos Geodetic\",\n");
        fprintf( fp, "#                              \"TITLE\": \"Ascending Node Crossing Pos Geodetic\",\n");
        fprintf( fp, "#                              \"LABEL\": \"Ascending Node Crossing Pos Geodetic\",\n");
        fprintf( fp, "#                          \"DIMENSION\": [ %d ],\n", nAscend );

        fprintf( fp, "#                             \"VALUES\": [ ");
        for (i=0; i<nAscend-1; i++) {
            Lgm_Set_Coord_Transforms( Ascend_UTC[i].Date, Ascend_UTC[i].Time, c );
            Lgm_Convert_Coords( &Ascend_U[i], &w, GEI2000_TO_WGS84, c );
            Lgm_WGS84_to_GEOD( &w, &GeodLat, &GeodLong, &GeodHeight );
            fprintf(fp, "\" %.6lf, %.6lf, %.6lf\", ", GeodLat, GeodLong, GeodHeight );
        }
        Lgm_Set_Coord_Transforms( Ascend_UTC[i].Date, Ascend_UTC[i].Time, c );
        Lgm_Convert_Coords( &Ascend_U[i], &w, GEI2000_TO_WGS84, c );
        Lgm_WGS84_to_GEOD( &w, &GeodLat, &GeodLong, &GeodHeight );
        fprintf(fp, "\" %.6lf, %.6lf, %.6lf\" ],\n", GeodLat, GeodLong, GeodHeight );


        fprintf( fp, "#                      \"ELEMENT_NAMES\": [ ");
        for (i=0; i<nAscend-1; i++) fprintf(fp, "\"AscendPosGeod%d\", ", i );
        fprintf(fp, "\"AscendPosGeod%d\" ],\n", i ); 

        fprintf( fp, "#                              \"UNITS\": \"Deg., Deg., km\" },\n");
        fprintf( fp, "#\n");
    }
    if ( nPerigee > 0 ) {
        fprintf( fp, "#  \"PerigeeTimes\":     { \"DESCRIPTION\": \"Perigee Times on this day.\",\n");
        fprintf( fp, "#                               \"NAME\": \"Perigee Times\",\n");
        fprintf( fp, "#                              \"TITLE\": \"Perigee Times\",\n");
        fprintf( fp, "#                              \"LABEL\": \"Perigee Times\",\n");
        fprintf( fp, "#                          \"DIMENSION\": [ %d ],\n", nPerigee );

        fprintf( fp, "#                             \"VALUES\": [ ");
        for (i=0; i<nPerigee-1; i++) {
            Lgm_DateTimeToString( IsoTimeString, &Perigee_UTC[i], 0, 3 );
            fprintf(fp, "\"%s\", ", IsoTimeString );
        }
        Lgm_DateTimeToString( IsoTimeString, &Perigee_UTC[i], 0, 3 );
        fprintf(fp, "\"%s\" ],\n", IsoTimeString ); 


        fprintf( fp, "#                      \"ELEMENT_NAMES\": [ ");
        for (i=0; i<nPerigee-1; i++) fprintf(fp, "\"PerigeeTime%d\", ", i );
        fprintf(fp, "\"PerigeeTime%d\" ],\n", i ); 

        fprintf( fp, "#                              \"UNITS\": \"UTC\" },\n");
        fprintf( fp, "#\n");



        fprintf( fp, "#  \"PerigeePosGeod\":   { \"DESCRIPTION\": \"Perigee Positions on this day in Geodetic Coords (lat/lon/alt).\",\n");
        fprintf( fp, "#                               \"NAME\": \"Perigee Pos Geodetic\",\n");
        fprintf( fp, "#                              \"TITLE\": \"Perigee Pos Geodetic\",\n");
        fprintf( fp, "#                              \"LABEL\": \"Perigee Pos Geodetic\",\n");
        fprintf( fp, "#                          \"DIMENSION\": [ %d ],\n", nPerigee );

        fprintf( fp, "#                             \"VALUES\": [ ");
        for (i=0; i<nPerigee-1; i++) {
            Lgm_Set_Coord_Transforms( Perigee_UTC[i].Date, Perigee_UTC[i].Time, c );
            Lgm_Convert_Coords( &Perigee_U[i], &w, GEI2000_TO_WGS84, c );
            Lgm_WGS84_to_GEOD( &w, &GeodLat, &GeodLong, &GeodHeight );
            fprintf(fp, "\" %.6lf, %.6lf, %.6lf\", ", GeodLat, GeodLong, GeodHeight );
        }
        Lgm_Set_Coord_Transforms( Perigee_UTC[i].Date, Perigee_UTC[i].Time, c );
        Lgm_Convert_Coords( &Perigee_U[i], &w, GEI2000_TO_WGS84, c );
        Lgm_WGS84_to_GEOD( &w, &GeodLat, &GeodLong, &GeodHeight );
        fprintf(fp, "\" %.6lf, %.6lf, %.6lf\" ],\n", GeodLat, GeodLong, GeodHeight );


        fprintf( fp, "#                      \"ELEMENT_NAMES\": [ ");
        for (i=0; i<nPerigee-1; i++) fprintf(fp, "\"PerigeePosGeod%d\", ", i );
        fprintf(fp, "\"PerigeePosGeod%d\" ],\n", i ); 

        fprintf( fp, "#                              \"UNITS\": \"Deg., Deg., km\" },\n");
        fprintf( fp, "#\n");
    }

    if ( nApogee > 0 ) {
        fprintf( fp, "#  \"ApogeeTimes\":      { \"DESCRIPTION\": \"Apogee Times on this day.\",\n");
        fprintf( fp, "#                               \"NAME\": \"Apogee Times\",\n");
        fprintf( fp, "#                              \"TITLE\": \"Apogee Times\",\n");
        fprintf( fp, "#                              \"LABEL\": \"Apogee Times\",\n");
        fprintf( fp, "#                          \"DIMENSION\": [ %d ],\n", nApogee );

        fprintf( fp, "#                             \"VALUES\": [ ");
        for (i=0; i<nApogee-1; i++) {
            Lgm_DateTimeToString( IsoTimeString, &Apogee_UTC[i], 0, 3 );
            fprintf(fp, "\"%s\", ", IsoTimeString );
        }
        Lgm_DateTimeToString( IsoTimeString, &Apogee_UTC[i], 0, 3 );
        fprintf(fp, "\"%s\" ],\n", IsoTimeString ); 


        fprintf( fp, "#                      \"ELEMENT_NAMES\": [ ");
        for (i=0; i<nApogee-1; i++) fprintf(fp, "\"ApogeeTime%d\", ", i );
        fprintf(fp, "\"ApogeeTime%d\" ],\n", i ); 

        fprintf( fp, "#                              \"UNITS\": \"UTC\" },\n");
        fprintf( fp, "#\n");

        fprintf( fp, "#  \"ApogeePosGeod\":   { \"DESCRIPTION\": \"Apogee Positions on this day in Geodetic Coords (lat/lon/alt).\",\n");
        fprintf( fp, "#                               \"NAME\": \"Apogee Pos Geodetic\",\n");
        fprintf( fp, "#                              \"TITLE\": \"Apogee Pos Geodetic\",\n");
        fprintf( fp, "#                              \"LABEL\": \"Apogee Pos Geodetic\",\n");
        fprintf( fp, "#                          \"DIMENSION\": [ %d ],\n", nApogee );

        fprintf( fp, "#                             \"VALUES\": [ ");
        for (i=0; i<nApogee-1; i++) {
            Lgm_Set_Coord_Transforms( Apogee_UTC[i].Date, Apogee_UTC[i].Time, c );
            Lgm_Convert_Coords( &Apogee_U[i], &w, GEI2000_TO_WGS84, c );
            Lgm_WGS84_to_GEOD( &w, &GeodLat, &GeodLong, &GeodHeight );
            fprintf(fp, "\" %.6lf, %.6lf, %.6lf\", ", GeodLat, GeodLong, GeodHeight );
        }
        Lgm_Set_Coord_Transforms( Apogee_UTC[i].Date, Apogee_UTC[i].Time, c );
        Lgm_Convert_Coords( &Apogee_U[i], &w, GEI2000_TO_WGS84, c );
        Lgm_WGS84_to_GEOD( &w, &GeodLat, &GeodLong, &GeodHeight );
        fprintf(fp, "\" %.6lf, %.6lf, %.6lf\" ],\n", GeodLat, GeodLong, GeodHeight );


        fprintf( fp, "#                      \"ELEMENT_NAMES\": [ ");
        for (i=0; i<nApogee-1; i++) fprintf(fp, "\"ApogeePosGeod%d\", ", i );
        fprintf(fp, "\"ApogeePosGeod%d\" ],\n", i ); 

        fprintf( fp, "#                              \"UNITS\": \"Deg., Deg., km\" },\n");
        fprintf( fp, "#\n");
    }


    fprintf( fp, "#  \"DateTime\":         { \"DESCRIPTION\": \"The date and time in ISO 8601 compliant format.\",\n");
    fprintf( fp, "#                               \"NAME\": \"IsoDateTime\",\n");
    fprintf( fp, "#                              \"TITLE\": \"ISO DateTime\",\n");
    fprintf( fp, "#                              \"LABEL\": \"Time\",\n");
    fprintf( fp, "#                              \"UNITS\": \"UTC\",\n");
    fprintf( fp, "#                       \"START_COLUMN\": %d },\n", nCol++);
    fprintf( fp, "#\n");

    fprintf( fp, "#  \"Date\":             { \"DESCRIPTION\": \"The date. In YYYMMDD format.\",\n");
    fprintf( fp, "#                               \"NAME\": \"Date\",\n");
    fprintf( fp, "#                              \"TITLE\": \"Date\",\n");
    fprintf( fp, "#                              \"LABEL\": \"Date\",\n");
    fprintf( fp, "#                       \"START_COLUMN\": %d },\n", nCol++);
    fprintf( fp, "#\n");

    fprintf( fp, "#  \"DOY\":              { \"DESCRIPTION\": \"Ordinal Day of Year.\",\n");
    fprintf( fp, "#                               \"NAME\": \"DOY\",\n");
    fprintf( fp, "#                              \"TITLE\": \"Day Of Year\",\n");
    fprintf( fp, "#                              \"LABEL\": \"DOY\",\n");
    fprintf( fp, "#                       \"START_COLUMN\": %d,\n", nCol++);
    fprintf( fp, "#                              \"UNITS\": \"Days\",\n");
    fprintf( fp, "#                          \"VALID_MIN\": 0,\n");
    fprintf( fp, "#                          \"VALID_MAX\": 366 },\n");
    fprintf( fp, "#\n");

    fprintf( fp, "#  \"UTC\":              { \"DESCRIPTION\": \"Universal Time (Coordinated). In decimal hours.\",\n");
    fprintf( fp, "#                               \"NAME\": \"UTC\",\n");
    fprintf( fp, "#                              \"TITLE\": \"UTC\",\n");
    fprintf( fp, "#                              \"LABEL\": \"UTC\",\n");
    fprintf( fp, "#                       \"START_COLUMN\": %d,\n", nCol++);
    fprintf( fp, "#                              \"UNITS\": \"Hours\",\n");
    fprintf( fp, "#                          \"VALID_MIN\": 0.0,\n");
    fprintf( fp, "#                          \"VALID_MAX\": 24.0 },\n");
    fprintf( fp, "#\n");





    fprintf( fp, "#  \"JulianDate\":       { \"DESCRIPTION\": \"Julian Date. In decimal days.\",\n");
    fprintf( fp, "#                               \"NAME\": \"JulianDate\",\n");
    fprintf( fp, "#                              \"TITLE\": \"Julian Date\",\n");
    fprintf( fp, "#                              \"LABEL\": \"Julian Date\",\n");
    fprintf( fp, "#                       \"START_COLUMN\": %d,\n", nCol++);
    fprintf( fp, "#                              \"UNITS\": \"Days\" },\n");
    fprintf( fp, "#\n");



    fprintf( fp, "#  \"GpsTime\":          { \"DESCRIPTION\": \"Number of SI seconds since 0h Jan 6, 1980 UTC.\",\n");
    fprintf( fp, "#                               \"NAME\": \"GpsTime\",\n");
    fprintf( fp, "#                              \"TITLE\": \"Gps Time\",\n");
    fprintf( fp, "#                              \"LABEL\": \"Gps Time (Seconds)\",\n");
    fprintf( fp, "#                       \"START_COLUMN\": %d,\n", nCol++);
    fprintf( fp, "#                              \"UNITS\": \"Seconds\" },\n");
    fprintf( fp, "#\n");



    fprintf( fp, "#  \"DipoleTiltAngle\":  { \"DESCRIPTION\": \"Angle between Zgsm and Zsm (i.e. between Zgsm and dipole axis direction).\",\n");
    fprintf( fp, "#                               \"NAME\": \"DipoleTiltAngle\",\n");
    fprintf( fp, "#                              \"TITLE\": \"Dipole Tilt Angle\",\n");
    fprintf( fp, "#                              \"LABEL\": \"Dipole Tilt Angle (Degrees)\",\n");
    fprintf( fp, "#                       \"START_COLUMN\": %d,\n", nCol++);
    fprintf( fp, "#                              \"UNITS\": \"Degrees\" },\n");
    fprintf( fp, "#\n");

    fprintf( fp, "#  \"InOut\":            { \"DESCRIPTION\": \"Flag indicating whether we are inbound (-1) or outbound (+1).\",\n");
    fprintf( fp, "#                               \"NAME\": \"InOut\",\n");
    fprintf( fp, "#                              \"TITLE\": \"InOut\",\n");
    fprintf( fp, "#                              \"LABEL\": \"InOut\",\n");
    fprintf( fp, "#                       \"START_COLUMN\": %d,\n", nCol++);
    fprintf( fp, "#                              \"UNITS\": \"dimless\" },\n");
    fprintf( fp, "#\n");

    fprintf( fp, "#  \"OrbitNumber\":      { \"DESCRIPTION\": \"Orbit Number.\",\n");
    fprintf( fp, "#                               \"NAME\": \"OrbitNumber\",\n");
    fprintf( fp, "#                              \"TITLE\": \"Orbit Number\",\n");
    fprintf( fp, "#                              \"LABEL\": \"Orbit Number\",\n");
    fprintf( fp, "#                       \"START_COLUMN\": %d,\n", nCol++);
    fprintf( fp, "#                              \"UNITS\": \"dimless\" },\n");
    fprintf( fp, "#\n");


    fprintf( fp, "#  \"Rgeo\":             { \"DESCRIPTION\": \"Geocentric Geographic position vector of S/C.\",\n");
    fprintf( fp, "#                               \"NAME\": \"Rgeo\",\n");
    fprintf( fp, "#                              \"TITLE\": \"Geocentric Geographic Position\",\n");
    fprintf( fp, "#                              \"LABEL\": \"Distance (R!BE)\",\n");
    fprintf( fp, "#                              \"UNITS\": \"R!BE!N\",\n");
    fprintf( fp, "#                          \"DIMENSION\": [ 3 ],\n");
    fprintf( fp, "#                       \"START_COLUMN\": %d,\n", nCol); nCol += 3;
    fprintf( fp, "#                      \"ELEMENT_NAMES\": [\"Rgeo_x\", \"Rgeo_y\", \"Rgeo_z\" ],\n");
    fprintf( fp, "#                     \"ELEMENT_LABELS\": [\"Rgeo!Bx!N ,  R!BE\", \"Rgeo!By!N ,  R!BE\", \"Rgeo!Bz!N ,  R!BE\" ],\n");
    fprintf( fp, "#                          \"VALID_MIN\": -1000.0,\n");
    fprintf( fp, "#                          \"VALID_MAX\":  1000.0,\n");
    fprintf( fp, "#                         \"FILL_VALUE\": -1e31 },\n");
    fprintf( fp, "#\n");


    fprintf( fp, "#  \"Rgeod_LatLon\":     { \"DESCRIPTION\": \"Geodetic Geographic Latitude and Longitude of S/C.\",\n");
    fprintf( fp, "#                               \"NAME\": \"Rgeod_LatLon\",\n");
    fprintf( fp, "#                              \"TITLE\": \"Geodetic Latitude and Longitude\",\n");
    fprintf( fp, "#                              \"LABEL\": \"R!Bgeod!N Latitude and Longitude (Degrees)\",\n");
    fprintf( fp, "#                              \"UNITS\": \"Degrees\",\n");
    fprintf( fp, "#                          \"DIMENSION\": [ 2 ],\n");
    fprintf( fp, "#                       \"START_COLUMN\": %d,\n", nCol); nCol += 2;
    fprintf( fp, "#                      \"ELEMENT_NAMES\": [ \"Rgeod_Lat\", \"Rgeod_Lon\" ],\n");
    fprintf( fp, "#                      \"ELEMENT_LABELS\": [ \"Geodetic Latitude\", \"Geodetic Longitude\" ],\n");
    fprintf( fp, "#                          \"VALID_MIN\": -360.0,\n");
    fprintf( fp, "#                          \"VALID_MAX\":  360.0,\n");
    fprintf( fp, "#                         \"FILL_VALUE\": -1e31 },\n");
    fprintf( fp, "#\n");

    fprintf( fp, "#  \"Rgeod_Height\":     { \"DESCRIPTION\": \"Geodetic Geographic Height (Above WGS84 Ellipsoid) of S/C.\",\n");
    fprintf( fp, "#                               \"NAME\": \"Rgeod_Height\",\n");
    fprintf( fp, "#                              \"TITLE\": \"Geodetic Height\",\n");
    fprintf( fp, "#                              \"LABEL\": \"R!Bgeod!N Height (km)\",\n");
    fprintf( fp, "#                              \"UNITS\": \"km\",\n");
    fprintf( fp, "#                       \"START_COLUMN\": %d,\n", nCol++);
    fprintf( fp, "#                          \"VALID_MIN\": -100.0,\n");
    fprintf( fp, "#                          \"VALID_MAX\":  1e7,\n");
    fprintf( fp, "#                         \"FILL_VALUE\": -1e31 },\n");
    fprintf( fp, "#\n");

    fprintf( fp, "#  \"Rgsm\":             { \"DESCRIPTION\": \"Geocentric Solar Magnetospheric position vector of S/C.\",\n");
    fprintf( fp, "#                               \"NAME\": \"Rgsm\",\n");
    fprintf( fp, "#                              \"TITLE\": \"Geocentric Solar Magnetospheric Position\",\n");
    fprintf( fp, "#                              \"LABEL\": \"Distance (R!BE)\",\n");
    fprintf( fp, "#                              \"UNITS\": \"R!BE!N\",\n");
    fprintf( fp, "#                          \"DIMENSION\": [ 3 ],\n");
    fprintf( fp, "#                       \"START_COLUMN\": %d,\n", nCol); nCol += 3;
    fprintf( fp, "#                      \"ELEMENT_NAMES\": [\"Rgsm_x\", \"Rgsm_y\", \"Rgsm_z\" ],\n");
    fprintf( fp, "#                     \"ELEMENT_LABELS\": [\"Rgsm!Bx!N ,  R!BE\", \"Rgsm!By!N ,  R!BE\", \"Rgsm!Bz!N ,  R!BE\" ],\n");
    fprintf( fp, "#                          \"VALID_MIN\": -1000.0,\n");
    fprintf( fp, "#                          \"VALID_MAX\":  1000.0,\n");
    fprintf( fp, "#                         \"FILL_VALUE\": -1e31 },\n");
    fprintf( fp, "#\n");

    fprintf( fp, "#  \"Rsm\":              { \"DESCRIPTION\": \"Geocentric Solar Magnetic position vector of S/C.\",\n");
    fprintf( fp, "#                               \"NAME\": \"Rsm\",\n");
    fprintf( fp, "#                              \"TITLE\": \"Geocentric Solar Magnetic Position\",\n");
    fprintf( fp, "#                              \"LABEL\": \"Distance (R!BE)\",\n");
    fprintf( fp, "#                              \"UNITS\": \"R!BE!N\",\n");
    fprintf( fp, "#                          \"DIMENSION\": [ 3 ],\n");
    fprintf( fp, "#                       \"START_COLUMN\": %d,\n", nCol); nCol += 3;
    fprintf( fp, "#                      \"ELEMENT_NAMES\": [\"Rsm_x\", \"Rsm_y\", \"Rsm_z\" ],\n");
    fprintf( fp, "#                     \"ELEMENT_LABELS\": [\"Rsm!Bx!N ,  R!BE\", \"Rsm!By!N ,  R!BE\", \"Rsm!Bz!N ,  R!BE\" ],\n");
    fprintf( fp, "#                          \"VALID_MIN\": -1000.0,\n");
    fprintf( fp, "#                          \"VALID_MAX\":  1000.0,\n");
    fprintf( fp, "#                         \"FILL_VALUE\": -1e31 },\n");
    fprintf( fp, "#\n");

    fprintf( fp, "#  \"Rgei\":             { \"DESCRIPTION\": \"Geocentric Equatorial Inertial position vector of S/C.\",\n");
    fprintf( fp, "#                               \"NAME\": \"Rgei\",\n");
    fprintf( fp, "#                              \"TITLE\": \"Geocentric Equatorial Inertial Position\",\n");
    fprintf( fp, "#                              \"LABEL\": \"Distance (R!BE)\",\n");
    fprintf( fp, "#                              \"UNITS\": \"R!BE!N\",\n");
    fprintf( fp, "#                          \"DIMENSION\": [ 3 ],\n");
    fprintf( fp, "#                       \"START_COLUMN\": %d,\n", nCol); nCol += 3;
    fprintf( fp, "#                      \"ELEMENT_NAMES\": [\"Rgei_x\", \"Rgei_y\", \"Rgei_z\" ],\n");
    fprintf( fp, "#                     \"ELEMENT_LABELS\": [\"Rgei!Bx!N ,  R!BE\", \"Rgei!By!N ,  R!BE\", \"Rgei!Bz!N ,  R!BE\" ],\n");
    fprintf( fp, "#                          \"VALID_MIN\": -1000.0,\n");
    fprintf( fp, "#                          \"VALID_MAX\":  1000.0,\n");
    fprintf( fp, "#                         \"FILL_VALUE\": -1e31 },\n");
    fprintf( fp, "#\n");

    fprintf( fp, "#  \"Rgse\":             { \"DESCRIPTION\": \"Geocentric Solar Ecliptic position vector of S/C.\",\n");
    fprintf( fp, "#                               \"NAME\": \"Rgse\",\n");
    fprintf( fp, "#                              \"TITLE\": \"Geocentric Solar Ecliptic Position\",\n");
    fprintf( fp, "#                              \"LABEL\": \"Distance (R!BE)\",\n");
    fprintf( fp, "#                              \"UNITS\": \"R!BE!N\",\n");
    fprintf( fp, "#                          \"DIMENSION\": [ 3 ],\n");
    fprintf( fp, "#                       \"START_COLUMN\": %d,\n", nCol); nCol += 3;
    fprintf( fp, "#                      \"ELEMENT_NAMES\": [\"Rgse_x\", \"Rgse_y\", \"Rgse_z\" ],\n");
    fprintf( fp, "#                     \"ELEMENT_LABELS\": [\"Rgse!Bx!N ,  R!BE\", \"Rgse!By!N ,  R!BE\", \"Rgse!Bz!N ,  R!BE\" ],\n");
    fprintf( fp, "#                          \"VALID_MIN\": -1000.0,\n");
    fprintf( fp, "#                          \"VALID_MAX\":  1000.0,\n");
    fprintf( fp, "#                         \"FILL_VALUE\": -1e31 },\n");
    fprintf( fp, "#\n");

    fprintf( fp, "#  \"CDMAG_MLAT\":       { \"DESCRIPTION\": \"Magnetic Latitude of S/C in Centerted Dipole Coordinates.\",\n");
    fprintf( fp, "#                               \"NAME\": \"CDMAG_MLAT\",\n");
    fprintf( fp, "#                              \"TITLE\": \"S/C Centered Dipole Magnetic Latitude\",\n");
    fprintf( fp, "#                              \"LABEL\": \"CD MLAT (Degrees)\",\n");
    fprintf( fp, "#                              \"UNITS\": \"Degrees\",\n");
    fprintf( fp, "#                       \"START_COLUMN\": %d,\n", nCol++);
    fprintf( fp, "#                          \"VALID_MIN\": -90.0,\n");
    fprintf( fp, "#                          \"VALID_MAX\":  90.0,\n");
    fprintf( fp, "#                         \"FILL_VALUE\": -1e31 },\n");
    fprintf( fp, "#\n");

    fprintf( fp, "#  \"CDMAG_MLON\":       { \"DESCRIPTION\": \"Magnetic Longitude of S/C Centerted Dipole Coordinates.\",\n");
    fprintf( fp, "#                               \"NAME\": \"CDMAG_MLON\",\n");
    fprintf( fp, "#                              \"TITLE\": \"S/C Centered Dipole Magnetic Longitude\",\n");
    fprintf( fp, "#                              \"LABEL\": \"CD MLON (Degrees)\",\n");
    fprintf( fp, "#                              \"UNITS\": \"Degrees\",\n");
    fprintf( fp, "#                       \"START_COLUMN\": %d,\n", nCol++);
    fprintf( fp, "#                          \"VALID_MIN\":   0.0,\n");
    fprintf( fp, "#                          \"VALID_MAX\": 360.0,\n");
    fprintf( fp, "#                         \"FILL_VALUE\": -1e31 },\n");
    fprintf( fp, "#\n");

    fprintf( fp, "#  \"CDMAG_MLT\":        { \"DESCRIPTION\": \"Magnetic Local Time of S/C in Centerted Dipole Coordinates.\",\n");
    fprintf( fp, "#                               \"NAME\": \"CDMAG_MLT\",\n");
    fprintf( fp, "#                              \"TITLE\": \"S/C Centered Dipole Magnetic Local Time (MLT)\",\n");
    fprintf( fp, "#                              \"LABEL\": \"CD MLT (Hours)\",\n");
    fprintf( fp, "#                              \"UNITS\": \"Hours\",\n");
    fprintf( fp, "#                       \"START_COLUMN\": %d,\n", nCol++);
    fprintf( fp, "#                          \"VALID_MIN\": 0.0,\n");
    fprintf( fp, "#                          \"VALID_MAX\": 24.0,\n");
    fprintf( fp, "#                         \"FILL_VALUE\": -1e31 },\n");
    fprintf( fp, "#\n");

    fprintf( fp, "#  \"CDMAG_R\":          { \"DESCRIPTION\": \"Radial distance of S/C from center of CDMAG coordinate system.\",\n");
    fprintf( fp, "#                               \"NAME\": \"CDMAG_R\",\n");
    fprintf( fp, "#                              \"TITLE\": \"S/C Centered Dipole Radial Distance\",\n");
    fprintf( fp, "#                              \"LABEL\": \"CD R (Re)\",\n");
    fprintf( fp, "#                              \"UNITS\": \"Re\",\n");
    fprintf( fp, "#                       \"START_COLUMN\": %d,\n", nCol++);
    fprintf( fp, "#                          \"VALID_MIN\": 0.0,\n");
    fprintf( fp, "#                          \"VALID_MAX\": 1000.0,\n");
    fprintf( fp, "#                         \"FILL_VALUE\": -1e31 },\n");
    fprintf( fp, "#\n");




    fprintf( fp, "#  \"EDMAG_MLAT\":       { \"DESCRIPTION\": \"Magnetic Latitude of S/C in Eccentric Dipole Coordinates.\",\n");
    fprintf( fp, "#                               \"NAME\": \"EDMAG_MLAT\",\n");
    fprintf( fp, "#                              \"TITLE\": \"S/C Eccentric Dipole Magnetic Latitude\",\n");
    fprintf( fp, "#                              \"LABEL\": \"ED MLAT (Degrees)\",\n");
    fprintf( fp, "#                              \"UNITS\": \"Degrees\",\n");
    fprintf( fp, "#                       \"START_COLUMN\": %d,\n", nCol++);
    fprintf( fp, "#                          \"VALID_MIN\": -90.0,\n");
    fprintf( fp, "#                          \"VALID_MAX\":  90.0,\n");
    fprintf( fp, "#                         \"FILL_VALUE\": -1e31 },\n");
    fprintf( fp, "#\n");

    fprintf( fp, "#  \"EDMAG_MLON\":       { \"DESCRIPTION\": \"Magnetic Longitude of S/C Eccentric Dipole Coordinates.\",\n");
    fprintf( fp, "#                               \"NAME\": \"EDMAG_MLON\",\n");
    fprintf( fp, "#                              \"TITLE\": \"S/C Eccentric Dipole Magnetic Longitude\",\n");
    fprintf( fp, "#                              \"LABEL\": \"ED MLON (Degrees)\",\n");
    fprintf( fp, "#                              \"UNITS\": \"Degrees\",\n");
    fprintf( fp, "#                       \"START_COLUMN\": %d,\n", nCol++);
    fprintf( fp, "#                          \"VALID_MIN\":   0.0,\n");
    fprintf( fp, "#                          \"VALID_MAX\": 360.0,\n");
    fprintf( fp, "#                         \"FILL_VALUE\": -1e31 },\n");
    fprintf( fp, "#\n");

    fprintf( fp, "#  \"EDMAG_MLT\":        { \"DESCRIPTION\": \"Magnetic Local Time of S/C in Eccentric Dipole Coordinates.\",\n");
    fprintf( fp, "#                               \"NAME\": \"EDMAG_MLT\",\n");
    fprintf( fp, "#                              \"TITLE\": \"S/C Eccentric Dipole Magnetic Local Time (MLT)\",\n");
    fprintf( fp, "#                              \"LABEL\": \"ED MLT (Hours)\",\n");
    fprintf( fp, "#                              \"UNITS\": \"Hours\",\n");
    fprintf( fp, "#                       \"START_COLUMN\": %d,\n", nCol++);
    fprintf( fp, "#                          \"VALID_MIN\": 0.0,\n");
    fprintf( fp, "#                          \"VALID_MAX\": 24.0,\n");
    fprintf( fp, "#                         \"FILL_VALUE\": -1e31 },\n");
    fprintf( fp, "#\n");

    fprintf( fp, "#  \"EDMAG_R\":          { \"DESCRIPTION\": \"Radial distance of S/C from center of EDMAG coordinate system.\",\n");
    fprintf( fp, "#                               \"NAME\": \"EDMAG_R\",\n");
    fprintf( fp, "#                              \"TITLE\": \"S/C Eccentric Dipole Radial Distance\",\n");
    fprintf( fp, "#                              \"LABEL\": \"ED R (Re)\",\n");
    fprintf( fp, "#                              \"UNITS\": \"Re\",\n");
    fprintf( fp, "#                       \"START_COLUMN\": %d,\n", nCol++);
    fprintf( fp, "#                          \"VALID_MIN\": 0.0,\n");
    fprintf( fp, "#                          \"VALID_MAX\": 1000.0,\n");
    fprintf( fp, "#                         \"FILL_VALUE\": -1e31 },\n");
    fprintf( fp, "#\n");




    fprintf( fp, "#  \"IntModel\":         { \"DESCRIPTION\": \"Internal magnetic field model.\",\n");
    fprintf( fp, "#                               \"NAME\": \"IntModel\",\n");
    fprintf( fp, "#                              \"TITLE\": \"Internal Magnetic Field Model\",\n");
    fprintf( fp, "#                              \"LABEL\": \"Internal Field Model\",\n");
    fprintf( fp, "#                       \"START_COLUMN\": %d,\n", nCol++);
    fprintf( fp, "#                               \"ENUM\": [ \"CDIP\", \"EDIP\", \"IGRF\" ]  },\n");
    fprintf( fp, "#\n");

    fprintf( fp, "#  \"ExtModel\":         { \"DESCRIPTION\": \"External magnetic field model.\",\n");
    fprintf( fp, "#                               \"NAME\": \"ExtModel\",\n");
    fprintf( fp, "#                              \"TITLE\": \"External Magnetic Field Model\",\n");
    fprintf( fp, "#                              \"LABEL\": \"External Field Model\",\n");
    fprintf( fp, "#                       \"START_COLUMN\": %d,\n", nCol++);
    fprintf( fp, "#                               \"ENUM\": [ \"OP77\", \"T87\", \"T89\" ]  },\n");
    fprintf( fp, "#\n");

    fprintf( fp, "#  \"Kp\":               { \"DESCRIPTION\": \"Kp index value.\",\n");
    fprintf( fp, "#                               \"NAME\": \"Kp\",\n");
    fprintf( fp, "#                              \"TITLE\": \"Kp Index\",\n");
    fprintf( fp, "#                              \"LABEL\": \"Kp\",\n");
    fprintf( fp, "#                       \"START_COLUMN\": %d,\n", nCol++);
    fprintf( fp, "#                          \"VALID_MIN\": 0.0,\n");
    fprintf( fp, "#                          \"VALID_MAX\": 9.0,\n");
    fprintf( fp, "#                         \"FILL_VALUE\": -1e31 },\n");
    fprintf( fp, "#\n");

    fprintf( fp, "#  \"Dst\":              { \"DESCRIPTION\": \"Dst index value.\",\n");
    fprintf( fp, "#                               \"NAME\": \"Dst\",\n");
    fprintf( fp, "#                              \"TITLE\": \"Dst Index\",\n");
    fprintf( fp, "#                              \"LABEL\": \"Dst (nT)\",\n");
    fprintf( fp, "#                              \"UNITS\": \"nT\",\n");
    fprintf( fp, "#                       \"START_COLUMN\": %d,\n", nCol++);
    fprintf( fp, "#                          \"VALID_MIN\": -10000.0,\n");
    fprintf( fp, "#                          \"VALID_MAX\":  10000.0,\n");
    fprintf( fp, "#                         \"FILL_VALUE\": -1e31 },\n");
    fprintf( fp, "#\n");


    fprintf( fp, "#  \"Bsc_gsm\":          { \"DESCRIPTION\": \"Magnetic field vector at S/C (in GSM coords).\",\n");
    fprintf( fp, "#                               \"NAME\": \"Bsc_gsm\",\n");
    fprintf( fp, "#                              \"TITLE\": \"Magnetic field vector at S/C\",\n");
    fprintf( fp, "#                              \"LABEL\": \"Bsc_gsm (nT)\",\n");
    fprintf( fp, "#                              \"UNITS\": \"nT\",\n");
    fprintf( fp, "#                          \"DIMENSION\": [ 4 ],\n");
    fprintf( fp, "#                       \"START_COLUMN\": %d,\n", nCol); nCol += 4;
    fprintf( fp, "#                      \"ELEMENT_NAMES\": [ \"Bsc_gsm_x\", \"Bsc_gsm_y\", \"Bsc_gsm_z\", \"Bsc_gsm_mag\" ],\n");
    fprintf( fp, "#                     \"ELEMENT_LABELS\": [ \"Bsc!Bgsm_x!N , nT\", \"Bsc!Bgsm_y!N , nT\", \"Bsc!Bgsm_z!N , nT\", \"Bsc_gsm_mag , nT\" ],\n");
    fprintf( fp, "#                          \"VALID_MIN\": -70000.0,\n");
    fprintf( fp, "#                          \"VALID_MAX\":  70000.0,\n");
    fprintf( fp, "#                         \"FILL_VALUE\": -1e31 },\n");
    fprintf( fp, "#\n");


    fprintf( fp, "#  \"FieldLineType\":    { \"DESCRIPTION\": \"Description of the type of field line the S/C is on., Can be one of 4 types: LGM_CLOSED      - FL hits Earth at both ends. LGM_OPEN_N_LOBE - FL is an OPEN field line rooted in the Northern polar cap. LGM_OPEN_S_LOBE - FL is an OPEN field line rooted in the Southern polar cap. LGM_OPEN_IMF    - FL does not hit Earth at eitrher end.\",\n");
    fprintf( fp, "#                               \"NAME\": \"FieldLineType\",\n");
    fprintf( fp, "#                              \"TITLE\": \"Field Line Type\",\n");
    fprintf( fp, "#                              \"LABEL\": \"Field Line Type\",\n");
    fprintf( fp, "#                       \"START_COLUMN\": %d,\n", nCol++);
    fprintf( fp, "#                               \"ENUM\": [ \"LGM_CLOSED\", \"LGM_OPEN_N_LOBE\", \"LGM_OPEN_S_LOBE\", \"LGM_OPEN_IMF\" ]  },\n");
    fprintf( fp, "#\n");


    fprintf( fp, "# \"S_sc_to_pfn\":       { \"DESCRIPTION\": \"Distance between S/C and Northern Footpoint along field line.\",\n");
    fprintf( fp, "#                               \"NAME\": \"S_sc_to_pfn\",\n");
    fprintf( fp, "#                              \"TITLE\": \"Distance to N. Foot. along FL\",\n");
    fprintf( fp, "#                              \"LABEL\": \"Distance to N. Foot. along FL (Re)\",\n");
    fprintf( fp, "#                              \"UNITS\": \"Re\",\n");
    fprintf( fp, "#                       \"START_COLUMN\": %d,\n", nCol++);
    //fprintf( fp, "#                          \"VALID_MIN\":  0.0,\n");
    //fprintf( fp, "#                          \"VALID_MAX\": 90.0,\n");
    fprintf( fp, "#                         \"FILL_VALUE\": -1e31 },\n");
    fprintf( fp, "#\n");

    fprintf( fp, "# \"S_sc_to_pfs\":       { \"DESCRIPTION\": \"Distance between S/C and Southern Footpoint along field line.\",\n");
    fprintf( fp, "#                               \"NAME\": \"S_sc_to_pfs\",\n");
    fprintf( fp, "#                              \"TITLE\": \"Distance to S. Foot. along FL\",\n");
    fprintf( fp, "#                              \"LABEL\": \"Distance to S. Foot. along FL (Re)\",\n");
    fprintf( fp, "#                              \"UNITS\": \"Re\",\n");
    fprintf( fp, "#                       \"START_COLUMN\": %d,\n", nCol++);
    //fprintf( fp, "#                          \"VALID_MIN\":  0.0,\n");
    //fprintf( fp, "#                          \"VALID_MAX\": 90.0,\n");
    fprintf( fp, "#                         \"FILL_VALUE\": -1e31 },\n");
    fprintf( fp, "#\n");

    fprintf( fp, "# \"S_pfs_to_Bmin\":     { \"DESCRIPTION\": \"Distance between Southern Footpoint and Bmin point along field line.\",\n");
    fprintf( fp, "#                               \"NAME\": \"S_pfs_to_Bmin\",\n");
    fprintf( fp, "#                              \"TITLE\": \"Distance from S. Foot. to Pmin along FL\",\n");
    fprintf( fp, "#                              \"LABEL\": \"Distance from S. Foot. to Pmin along FL (Re)\",\n");
    fprintf( fp, "#                              \"UNITS\": \"Re\",\n");
    fprintf( fp, "#                       \"START_COLUMN\": %d,\n", nCol++);
    //fprintf( fp, "#                          \"VALID_MIN\":  0.0,\n");
    //fprintf( fp, "#                          \"VALID_MAX\": 90.0,\n");
    fprintf( fp, "#                         \"FILL_VALUE\": -1e31 },\n");
    fprintf( fp, "#\n");

    fprintf( fp, "# \"S_Bmin_to_sc\":     { \"DESCRIPTION\": \"Distance between Bmin point and S/C along field line (positive if north of Bmin).\",\n");
    fprintf( fp, "#                               \"NAME\": \"S_Bmin_to_sc\",\n");
    fprintf( fp, "#                              \"TITLE\": \"Distance from Pmin to S/C along FL\",\n");
    fprintf( fp, "#                              \"LABEL\": \"Distance from Pmin to S/C along FL (Re)\",\n");
    fprintf( fp, "#                              \"UNITS\": \"Re\",\n");
    fprintf( fp, "#                       \"START_COLUMN\": %d,\n", nCol++);
    //fprintf( fp, "#                          \"VALID_MIN\":  0.0,\n");
    //fprintf( fp, "#                          \"VALID_MAX\": 90.0,\n");
    fprintf( fp, "#                         \"FILL_VALUE\": -1e31 },\n");
    fprintf( fp, "#\n");

    fprintf( fp, "# \"S_total\":           { \"DESCRIPTION\": \"Total Field Line length (along field line).\",\n");
    fprintf( fp, "#                               \"NAME\": \"S_total\",\n");
    fprintf( fp, "#                              \"TITLE\": \"Field Line length\",\n");
    fprintf( fp, "#                              \"LABEL\": \"Field Line length (Re)\",\n");
    fprintf( fp, "#                              \"UNITS\": \"Re\",\n");
    fprintf( fp, "#                       \"START_COLUMN\": %d,\n", nCol++);
    //fprintf( fp, "#                          \"VALID_MIN\":  0.0,\n");
    //fprintf( fp, "#                          \"VALID_MAX\": 90.0,\n");
    fprintf( fp, "#                         \"FILL_VALUE\": -1e31 },\n");
    fprintf( fp, "#\n");

    fprintf( fp, "# \"d2B_ds2\":           { \"DESCRIPTION\": \"Second derivative of |B| with respect to s (dist along FL) at minimum |B| point.\",\n");
    fprintf( fp, "#                               \"NAME\": \"d2B_ds2\",\n");
    fprintf( fp, "#                              \"TITLE\": \"d!A2!NB/ds!A2!N\",\n");
    fprintf( fp, "#                              \"LABEL\": \"d!A2!NB/ds!A2!N (nT!A2!N/Re!A2!N)\",\n");
    fprintf( fp, "#                              \"UNITS\": \"nT!A2!N/Re!A2!N\",\n");
    fprintf( fp, "#                       \"START_COLUMN\": %d,\n", nCol++);
    //fprintf( fp, "#                          \"VALID_MIN\":  0.0,\n");
    //fprintf( fp, "#                          \"VALID_MAX\": 90.0,\n");
    fprintf( fp, "#                         \"FILL_VALUE\": -1e31 },\n");
    fprintf( fp, "#\n");

    fprintf( fp, "# \"Sb0\":               { \"DESCRIPTION\": \"Value of the 'Sb Integral' for equatorially mirroring particles (not generally zero).\",\n");
    fprintf( fp, "#                               \"NAME\": \"Sb0\",\n");
    fprintf( fp, "#                              \"TITLE\": \"Sb0\",\n");
    fprintf( fp, "#                              \"LABEL\": \"Sb0 (Re)\",\n");
    fprintf( fp, "#                              \"UNITS\": \"Re\",\n");
    fprintf( fp, "#                       \"START_COLUMN\": %d,\n", nCol++);
    //fprintf( fp, "#                          \"VALID_MIN\":  0.0,\n");
    //fprintf( fp, "#                          \"VALID_MAX\": 90.0,\n");
    fprintf( fp, "#                         \"FILL_VALUE\": -1e31 },\n");
    fprintf( fp, "#\n");

    fprintf( fp, "# \"RadiusOfCurv\":      { \"DESCRIPTION\": \"Field line radius of curvature at minimum |B| point.\",\n");
    fprintf( fp, "#                               \"NAME\": \"RadiusOfCurv\",\n");
    fprintf( fp, "#                              \"TITLE\": \"FL Radius of Curvature at Bmin\",\n");
    fprintf( fp, "#                              \"LABEL\": \"Radius of Curvature (R!BE!N)\",\n");
    fprintf( fp, "#                              \"UNITS\": \"R!BE!N\",\n");
    fprintf( fp, "#                       \"START_COLUMN\": %d,\n", nCol++);
    //fprintf( fp, "#                          \"VALID_MIN\":  0.0,\n");
    //fprintf( fp, "#                          \"VALID_MAX\": 90.0,\n");
    fprintf( fp, "#                         \"FILL_VALUE\": -1e31 },\n");
    fprintf( fp, "#\n");






    fprintf( fp, "#  \"Pfn_geo\":          { \"DESCRIPTION\": \"Location of Northern Footpoint (in GEO coords).\",\n");
    fprintf( fp, "#                               \"NAME\": \"Pfn_geo\",\n");
    fprintf( fp, "#                              \"TITLE\": \"Geographic Northern Footpoint\",\n");
    fprintf( fp, "#                              \"LABEL\": \"Pfn_geo (R!BE)\",\n");
    fprintf( fp, "#                              \"UNITS\": \"R!BE!N\",\n");
    fprintf( fp, "#                          \"DIMENSION\": [ 3 ],\n");
    fprintf( fp, "#                       \"START_COLUMN\": %d,\n", nCol); nCol += 3;
    fprintf( fp, "#                      \"ELEMENT_NAMES\": [\"Pfn_geo_x\", \"Pfn_geo_y\", \"Pfn_geo_z\"],\n");
    fprintf( fp, "#                     \"ELEMENT_LABELS\": [\"Pfn_geo!Bx!N , R!BE\", \"Pfn_geo!By!N , R!BE\", \"Pfn_geo!Bz!N , R!BE\"],\n");
    fprintf( fp, "#                          \"VALID_MIN\": -1000.0,\n");
    fprintf( fp, "#                          \"VALID_MAX\":  1000.0,\n");
    fprintf( fp, "#                         \"FILL_VALUE\": -1e31 },\n");
    fprintf( fp, "#\n");

    fprintf( fp, "#  \"Pfn_gsm\":          { \"DESCRIPTION\": \"Location of Northern Footpoint (in GSM coords).\",\n");
    fprintf( fp, "#                               \"NAME\": \"Pfn_gsm\",\n");
    fprintf( fp, "#                              \"TITLE\": \"GSM Position of Northern Footpoint\",\n");
    fprintf( fp, "#                              \"LABEL\": \"Pfn_gsm (R!BE)\",\n");
    fprintf( fp, "#                              \"UNITS\": \"R!BE!N\",\n");
    fprintf( fp, "#                          \"DIMENSION\": [ 3 ],\n");
    fprintf( fp, "#                       \"START_COLUMN\": %d,\n", nCol); nCol += 3;
    fprintf( fp, "#                      \"ELEMENT_NAMES\": [ \"Pfn_gsm_x\", \"Pfn_gsm_y\", \"Pfn_gsm_z\" ],\n");
    fprintf( fp, "#                     \"ELEMENT_LABELS\": [\"Pfn_gsm!Bx!N , R!BE\", \"Pfn_gsm!By!N , R!BE\", \"Pfn_gsm!Bz!N , R!BE\"],\n");
    fprintf( fp, "#                          \"VALID_MIN\": -1000.0,\n");
    fprintf( fp, "#                          \"VALID_MAX\":  1000.0,\n");
    fprintf( fp, "#                         \"FILL_VALUE\": -1e31 },\n");
    fprintf( fp, "#\n");

    fprintf( fp, "#  \"Pfn_geod_LatLon\":  { \"DESCRIPTION\": \"Geodetic Latitude and Longitude of Northern Footpoint.\",\n");
    fprintf( fp, "#                               \"NAME\": \"Pfn_geod_LatLon\",\n");
    fprintf( fp, "#                              \"TITLE\": \"Geodetic Latitude and Longitude of Northern Footpoint\",\n");
    fprintf( fp, "#                              \"LABEL\": \"Pfn_geod_LatLon (Degrees)\",\n");
    fprintf( fp, "#                              \"UNITS\": \"Degrees\",\n");
    fprintf( fp, "#                          \"DIMENSION\": [ 2 ],\n");
    fprintf( fp, "#                       \"START_COLUMN\": %d,\n", nCol); nCol += 2;
    fprintf( fp, "#                      \"ELEMENT_NAMES\": [ \"Pfn_geod_Lat\", \"Pfn_geod_Lon\" ],\n");
    fprintf( fp, "#                     \"ELEMENT_LABELS\": [ \"Pfn_geod Latitude, Degrees\", \"Pfn_geod Longitude, Degrees\" ],\n");
    fprintf( fp, "#                          \"VALID_MIN\": -360.0,\n");
    fprintf( fp, "#                          \"VALID_MAX\":  360.0,\n");
    fprintf( fp, "#                         \"FILL_VALUE\": -1e31 },\n");
    fprintf( fp, "#\n");

    fprintf( fp, "#  \"Pfn_geod_Height\":  { \"DESCRIPTION\": \"Geodetic Height of Northern Footpoint.\",\n");
    fprintf( fp, "#                               \"NAME\": \"Pfn_geod_Height\",\n");
    fprintf( fp, "#                              \"TITLE\": \"Geodetic Height of Northern Footpoint\",\n");
    fprintf( fp, "#                              \"LABEL\": \"Pfn_geod_Height (km)\",\n");
    fprintf( fp, "#                              \"UNITS\": \"km\",\n");
    fprintf( fp, "#                       \"START_COLUMN\": %d,\n", nCol++);
    fprintf( fp, "#                          \"VALID_MIN\": -100.0,\n");
    fprintf( fp, "#                          \"VALID_MAX\": 1e7,\n");
    fprintf( fp, "#                         \"FILL_VALUE\": -1e31 },\n");
    fprintf( fp, "#\n");

    fprintf( fp, "#  \"Pfn_CD_MLAT\":      { \"DESCRIPTION\": \"Magnetic Latitude of Northern Footpoint in Centerted Dipole Coordinates.\",\n");
    fprintf( fp, "#                               \"NAME\": \"Pfn_CD_MLAT\",\n");
    fprintf( fp, "#                              \"TITLE\": \"CD MLAT of Northern Footpoint\",\n");
    fprintf( fp, "#                              \"LABEL\": \"CD MLAT (Degrees)\",\n");
    fprintf( fp, "#                              \"UNITS\": \"Degrees\",\n");
    fprintf( fp, "#                       \"START_COLUMN\": %d,\n", nCol++);
    fprintf( fp, "#                          \"VALID_MIN\": -90.0,\n");
    fprintf( fp, "#                          \"VALID_MAX\":  90.0,\n");
    fprintf( fp, "#                         \"FILL_VALUE\": -1e31 },\n");
    fprintf( fp, "#\n");

    fprintf( fp, "#  \"Pfn_CD_MLON\":      { \"DESCRIPTION\":\"Magnetic Longitude of Northern Footpoint Centerted Dipole Coordinates.\",\n");
    fprintf( fp, "#                               \"NAME\": \"Pfn_CD_MLON\",\n");
    fprintf( fp, "#                              \"TITLE\": \"CD MLON of Northern Footpoint\",\n");
    fprintf( fp, "#                              \"LABEL\": \"CD MLON (Degrees)\",\n");
    fprintf( fp, "#                              \"UNITS\": \"Degrees\",\n");
    fprintf( fp, "#                       \"START_COLUMN\": %d,\n", nCol++);
    fprintf( fp, "#                          \"VALID_MIN\":   0.0,\n");
    fprintf( fp, "#                          \"VALID_MAX\": 360.0,\n");
    fprintf( fp, "#                         \"FILL_VALUE\": -1e31 },\n");
    fprintf( fp, "#\n");

    fprintf( fp, "#  \"Pfn_CD_MLT\":       { \"DESCRIPTION\": \"Magnetic Local Time of Northern Footpoint in Centerted Dipole Coordinates.\",\n");
    fprintf( fp, "#                               \"NAME\": \"Pfn_CD_MLT\",\n");
    fprintf( fp, "#                              \"TITLE\": \"CD MLT of Northern Footpoint\",\n");
    fprintf( fp, "#                              \"LABEL\": \"CD MLT (Hours)\",\n");
    fprintf( fp, "#                              \"UNITS\": \"Hours\",\n");
    fprintf( fp, "#                       \"START_COLUMN\": %d,\n", nCol++);
    fprintf( fp, "#                          \"VALID_MIN\": 0.0,\n");
    fprintf( fp, "#                          \"VALID_MAX\": 24.0,\n");
    fprintf( fp, "#                         \"FILL_VALUE\": -1e31 },\n");
    fprintf( fp, "#\n");

    fprintf( fp, "#  \"Pfn_ED_MLAT\":      { \"DESCRIPTION\": \"Magnetic Latitude of Northern Footpoint in Eccentric Dipole Coordinates.\",\n");
    fprintf( fp, "#                               \"NAME\": \"Pfn_ED_MLAT\",\n");
    fprintf( fp, "#                              \"TITLE\": \"ED MLAT of Northern Footpoint\",\n");
    fprintf( fp, "#                              \"LABEL\": \"ED MLAT (Degrees)\",\n");
    fprintf( fp, "#                              \"UNITS\": \"Degrees\",\n");
    fprintf( fp, "#                       \"START_COLUMN\": %d,\n", nCol++);
    fprintf( fp, "#                          \"VALID_MIN\": -90.0,\n");
    fprintf( fp, "#                          \"VALID_MAX\":  90.0,\n");
    fprintf( fp, "#                         \"FILL_VALUE\": -1e31 },\n");
    fprintf( fp, "#\n");

    fprintf( fp, "#  \"Pfn_ED_MLON\":      { \"DESCRIPTION\":  \"Magnetic Longitude of Northern Footpoint Eccentric Dipole Coordinates.\",\n");
    fprintf( fp, "#                               \"NAME\": \"Pfn_ED_MLON\",\n");
    fprintf( fp, "#                              \"TITLE\": \"ED MLON of Northern Footpoint\",\n");
    fprintf( fp, "#                              \"LABEL\": \"ED MLON (Degrees)\",\n");
    fprintf( fp, "#                              \"UNITS\": \"Degrees\",\n");
    fprintf( fp, "#                       \"START_COLUMN\": %d,\n", nCol++);
    fprintf( fp, "#                          \"VALID_MIN\":   0.0,\n");
    fprintf( fp, "#                          \"VALID_MAX\": 360.0,\n");
    fprintf( fp, "#                         \"FILL_VALUE\": -1e31 },\n");
    fprintf( fp, "#\n");

    fprintf( fp, "#  \"Pfn_ED_MLT\":       { \"DESCRIPTION\":\"Magnetic Local Time of Northern Footpoint in Eccentric Dipole Coordinates.\",\n");
    fprintf( fp, "#                               \"NAME\": \"Pfn_ED_MLT\",\n");
    fprintf( fp, "#                              \"TITLE\": \"ED MLT of Northern Footpoint\",\n");
    fprintf( fp, "#                              \"LABEL\": \"ED MLT (Hours)\",\n");
    fprintf( fp, "#                              \"UNITS\": \"Hours\",\n");
    fprintf( fp, "#                       \"START_COLUMN\": %d,\n", nCol++);
    fprintf( fp, "#                          \"VALID_MIN\": 0.0,\n");
    fprintf( fp, "#                          \"VALID_MAX\": 24.0,\n");
    fprintf( fp, "#                         \"FILL_VALUE\": -1e31 },\n");
    fprintf( fp, "#\n");


    fprintf( fp, "#  \"Bfn_geo\":          { \"DESCRIPTION\": \"Magnetic field vector at Northern Footpoint (in GEO coords).\",\n");
    fprintf( fp, "#                               \"NAME\": \"Bfn_geo\",\n");
    fprintf( fp, "#                              \"TITLE\": \"Magnetic field vector at Northern Footpoint (in GEO coords)\",\n");
    fprintf( fp, "#                              \"LABEL\": \"Bfn_geo (nT)\",\n");
    fprintf( fp, "#                              \"UNITS\": \"nT\",\n");
    fprintf( fp, "#                          \"DIMENSION\": [ 4 ],\n");
    fprintf( fp, "#                       \"START_COLUMN\": %d,\n", nCol); nCol += 4;
    fprintf( fp, "#                      \"ELEMENT_NAMES\": [ \"Bfn_geo_x\", \"Bfn_geo_y\", \"Bfn_geo_z\", \"Bfn_geo_mag\" ],\n");
    fprintf( fp, "#                     \"ELEMENT_LABELS\": [ \"Bfn_geo!Bx!N , nT\", \"Bfn_geo!By!N , nT\", \"Bfn_geo!Bz!N , nT\", \"Bfn_geo_mag , nT\" ],\n");
    fprintf( fp, "#                          \"VALID_MIN\": -70000.0,\n");
    fprintf( fp, "#                          \"VALID_MAX\":  70000.0,\n");
    fprintf( fp, "#                         \"FILL_VALUE\": -1e31 },\n");
    fprintf( fp, "#\n");

    fprintf( fp, "#  \"Bfn_gsm\":          { \"DESCRIPTION\": \"Magnetic field vector at Northern Footpoint (in GSM coords).\",\n");
    fprintf( fp, "#                               \"NAME\": \"Bfn_gsm\",\n");
    fprintf( fp, "#                              \"TITLE\": \"Magnetic field vector at Northern Footpoint (in GSM coords)\",\n");
    fprintf( fp, "#                              \"LABEL\": \"Bfn_gsm (nT)\",\n");
    fprintf( fp, "#                              \"UNITS\": \"nT\",\n");
    fprintf( fp, "#                          \"DIMENSION\": [ 4 ],\n");
    fprintf( fp, "#                       \"START_COLUMN\": %d,\n", nCol); nCol += 4;
    fprintf( fp, "#                      \"ELEMENT_NAMES\": [ \"Bfn_gsm_x\", \"Bfn_gsm_y\", \"Bfn_gsm_z\", \"Bfn_gsm_mag\" ],\n");
    fprintf( fp, "#                     \"ELEMENT_LABELS\": [ \"Bfn_gsm!Bx!N , nT\", \"Bfn_gsm!By!N , nT\", \"Bfn_gsm!Bz!N , nT\", \"Bfn_gsm_mag , nT\" ],\n");
    fprintf( fp, "#                          \"VALID_MIN\": -70000.0,\n");
    fprintf( fp, "#                          \"VALID_MAX\":  70000.0,\n");
    fprintf( fp, "#                         \"FILL_VALUE\": -1e31 },\n");
    fprintf( fp, "#\n");

    fprintf( fp, "# \"Loss_Cone_Alpha_n\": { \"DESCRIPTION\": \"Value of Northern Loss Cone angle. asin( sqrt(Bsc/Bfn) ). Bsc is the magntiude of B at spacecraft,"
                                                               " Bfn is the magnitude of B at the northern footpoint. Note, loss-cone height is taken to be equal"
                                                               " to the footpoint height which is %g km above the WGS84 geoid.\",\n", m->LstarInfo->mInfo->Lgm_LossConeHeight );
    fprintf( fp, "#                               \"NAME\": \"Loss_Cone_Alpha_n\",\n");
    fprintf( fp, "#                              \"TITLE\": \"Northern Loss Cone Angle\",\n");
    fprintf( fp, "#                              \"LABEL\": \"Northern Loss Cone Angle (Degrees)\",\n");
    fprintf( fp, "#                              \"UNITS\": \"Degrees\",\n");
    fprintf( fp, "#                       \"START_COLUMN\": %d,\n", nCol++);
    //fprintf( fp, "#                          \"VALID_MIN\":  0.0,\n");
    //fprintf( fp, "#                          \"VALID_MAX\": 90.0,\n");
    fprintf( fp, "#                         \"FILL_VALUE\": -1e31 },\n");
    fprintf( fp, "#\n");











    fprintf( fp, "#  \"Pfs_geo\":          { \"DESCRIPTION\": \"Location of Southern Footpoint (in GEO coords).\",\n");
    fprintf( fp, "#                               \"NAME\": \"Pfs_geo\",\n");
    fprintf( fp, "#                              \"TITLE\": \"Geographic Southern Footpoint\",\n");
    fprintf( fp, "#                              \"LABEL\": \"Pfs_geo (R!BE)\",\n");
    fprintf( fp, "#                              \"UNITS\": \"R!BE!N\",\n");
    fprintf( fp, "#                          \"DIMENSION\": [ 3 ],\n");
    fprintf( fp, "#                       \"START_COLUMN\": %d,\n", nCol); nCol += 3;
    fprintf( fp, "#                      \"ELEMENT_NAMES\": [ \"Pfs_geo_x\", \"Pfs_geo_y\", \"Pfs_geo_z\" ],\n");
    fprintf( fp, "#                     \"ELEMENT_LABELS\": [ \"Pfs_geo!Bx!N , R!BE\", \"Pfs_geo!By!N , R!BE\", \"Pfs_geo!Bz!N , R!BE\" ],\n");
    fprintf( fp, "#                          \"VALID_MIN\": -1000.0,\n");
    fprintf( fp, "#                          \"VALID_MAX\":  1000.0,\n");
    fprintf( fp, "#                         \"FILL_VALUE\": -1e31 },\n");
    fprintf( fp, "#\n");

    fprintf( fp, "#  \"Pfs_gsm\":          { \"DESCRIPTION\": \"Location of Southern Footpoint (in GSM coords).\",\n");
    fprintf( fp, "#                               \"NAME\": \"Pfs_gsm\",\n");
    fprintf( fp, "#                              \"TITLE\": \"GSM Southern Footpoint\",\n");
    fprintf( fp, "#                              \"LABEL\": \"Pfs_gsm (R!BE)\",\n");
    fprintf( fp, "#                              \"UNITS\": \"R!BE!N\",\n");
    fprintf( fp, "#                          \"DIMENSION\": [ 3 ],\n");
    fprintf( fp, "#                       \"START_COLUMN\": %d,\n", nCol); nCol += 3;
    fprintf( fp, "#                      \"ELEMENT_NAMES\": [ \"Pfs_gsm_x\", \"Pfs_gsm_y\", \"Pfs_gsm_z\" ],\n");
    fprintf( fp, "#                     \"ELEMENT_LABELS\": [ \"Pfs_gsm!Bx!N , R!BE\", \"Pfs_gsm!By!N , R!BE\", \"Pfs_gsm!Bz!N , R!BE\" ],\n");
    fprintf( fp, "#                          \"VALID_MIN\": -1000.0,\n");
    fprintf( fp, "#                          \"VALID_MAX\":  1000.0,\n");
    fprintf( fp, "#                         \"FILL_VALUE\": -1e31 },\n");
    fprintf( fp, "#\n");

    fprintf( fp, "#  \"Pfs_geod_LatLon\":  { \"DESCRIPTION\": \"Geodetic Latitude and Longitude of Southern Footpoint.\",\n");
    fprintf( fp, "#                               \"NAME\": \"Pfs_geod_LatLon\",\n");
    fprintf( fp, "#                              \"TITLE\": \"Geodetic Latitude and Longitude of Southern Footpoint\",\n");
    fprintf( fp, "#                              \"LABEL\": \"Pfs_geod_LatLon (Degrees)\",\n");
    fprintf( fp, "#                              \"UNITS\": \"Degrees\",\n");
    fprintf( fp, "#                          \"DIMENSION\": [ 2 ],\n");
    fprintf( fp, "#                       \"START_COLUMN\": %d,\n", nCol); nCol += 2;
    fprintf( fp, "#                      \"ELEMENT_NAMES\": [ \"Pfs_geod_Lat\", \"Pfs_geod_Lon\" ],\n");
    fprintf( fp, "#                     \"ELEMENT_LABELS\": [ \"Pfs_geod Latitude, Degrees\", \"Pfs_geod Longitude, Degrees\" ],\n");
    fprintf( fp, "#                          \"VALID_MIN\": -360.0,\n");
    fprintf( fp, "#                          \"VALID_MAX\":  360.0,\n");
    fprintf( fp, "#                         \"FILL_VALUE\": -1e31 },\n");
    fprintf( fp, "#\n");

    fprintf( fp, "#  \"Pfs_geod_Height\":  { \"DESCRIPTION\": \"Geodetic Height of Southern Footpoint.\",\n");
    fprintf( fp, "#                               \"NAME\": \"Pfs_geod_Height\",\n");
    fprintf( fp, "#                              \"TITLE\": \"Geodetic Height of Southern Footpoint\",\n");
    fprintf( fp, "#                              \"LABEL\": \"Pfs_geod_Height (km)\",\n");
    fprintf( fp, "#                              \"UNITS\": \"km\",\n");
    fprintf( fp, "#                       \"START_COLUMN\": %d,\n", nCol++);
    fprintf( fp, "#                          \"VALID_MIN\": -100.0,\n");
    fprintf( fp, "#                          \"VALID_MAX\": 1e7,\n");
    fprintf( fp, "#                         \"FILL_VALUE\": -1e31 },\n");
    fprintf( fp, "#\n");

    fprintf( fp, "#  \"Pfs_CD_MLAT\":      { \"DESCRIPTION\": \"Magnetic Latitude of Southern Footpoint in Centerted Dipole Coordinates.\",\n");
    fprintf( fp, "#                               \"NAME\": \"Pfs_CD_MLAT\",\n");
    fprintf( fp, "#                              \"TITLE\": \"CD MLAT of Southern Footpoint\",\n");
    fprintf( fp, "#                              \"LABEL\": \"CD MLAT (Degrees)\",\n");
    fprintf( fp, "#                              \"UNITS\": \"Degrees\",\n");
    fprintf( fp, "#                       \"START_COLUMN\": %d,\n", nCol++);
    fprintf( fp, "#                          \"VALID_MIN\": -90.0,\n");
    fprintf( fp, "#                          \"VALID_MAX\":  90.0,\n");
    fprintf( fp, "#                         \"FILL_VALUE\": -1e31 },\n");
    fprintf( fp, "#\n");

    fprintf( fp, "#  \"Pfs_CD_MLON\":      { \"DESCRIPTION\": \"Magnetic Longitude of Southern Footpoint Centerted Dipole Coordinates.\",\n");
    fprintf( fp, "#                               \"NAME\": \"Pfs_CD_MLON\",\n");
    fprintf( fp, "#                              \"TITLE\": \"CD MLON of Southern Footpoint\",\n");
    fprintf( fp, "#                              \"LABEL\": \"CD MLON (Degrees)\",\n");
    fprintf( fp, "#                              \"UNITS\": \"Degrees\",\n");
    fprintf( fp, "#                       \"START_COLUMN\": %d,\n", nCol++);
    fprintf( fp, "#                          \"VALID_MIN\":   0.0,\n");
    fprintf( fp, "#                          \"VALID_MAX\": 360.0,\n");
    fprintf( fp, "#                         \"FILL_VALUE\": -1e31 },\n");
    fprintf( fp, "#\n");

    fprintf( fp, "#  \"Pfs_CD_MLT\":       { \"DESCRIPTION\": \"Magnetic Local Time of Southern Footpoint in Centerted Dipole Coordinates.\",\n");
    fprintf( fp, "#                               \"NAME\": \"Pfs_CD_MLT\",\n");
    fprintf( fp, "#                              \"TITLE\": \"CD MLT of Southern Footpoint\",\n");
    fprintf( fp, "#                              \"LABEL\": \"CD MLT (Hours)\",\n");
    fprintf( fp, "#                              \"UNITS\": \"Hours\",\n");
    fprintf( fp, "#                       \"START_COLUMN\": %d,\n", nCol++);
    fprintf( fp, "#                          \"VALID_MIN\": 0.0,\n");
    fprintf( fp, "#                          \"VALID_MAX\": 24.0,\n");
    fprintf( fp, "#                         \"FILL_VALUE\": -1e31 },\n");
    fprintf( fp, "#\n");

    fprintf( fp, "#  \"Pfs_ED_MLAT\":      { \"DESCRIPTION\":  \"Magnetic Latitude of Southern Footpoint in Eccentric Dipole Coordinates.\",\n");
    fprintf( fp, "#                               \"NAME\": \"Pfs_ED_MLAT\",\n");
    fprintf( fp, "#                              \"TITLE\": \"ED MLAT of Southern Footpoint\",\n");
    fprintf( fp, "#                              \"LABEL\": \"ED MLAT (Degrees)\",\n");
    fprintf( fp, "#                              \"UNITS\": \"Degrees\",\n");
    fprintf( fp, "#                       \"START_COLUMN\": %d,\n", nCol++);
    fprintf( fp, "#                          \"VALID_MIN\": -90.0,\n");
    fprintf( fp, "#                          \"VALID_MAX\":  90.0,\n");
    fprintf( fp, "#                         \"FILL_VALUE\": -1e31 },\n");
    fprintf( fp, "#\n");

    fprintf( fp, "#  \"Pfs_ED_MLON\":      { \"DESCRIPTION\":  \"Magnetic Longitude of Southern Footpoint Eccentric Dipole Coordinates.\",\n");
    fprintf( fp, "#                               \"NAME\": \"Pfs_ED_MLON\",\n");
    fprintf( fp, "#                              \"TITLE\": \"ED MLON of Southern Footpoint\",\n");
    fprintf( fp, "#                              \"LABEL\": \"ED MLON (Degrees)\",\n");
    fprintf( fp, "#                              \"UNITS\": \"Degrees\",\n");
    fprintf( fp, "#                       \"START_COLUMN\": %d,\n", nCol++);
    fprintf( fp, "#                          \"VALID_MIN\":   0.0,\n");
    fprintf( fp, "#                          \"VALID_MAX\": 360.0,\n");
    fprintf( fp, "#                         \"FILL_VALUE\": -1e31 },\n");
    fprintf( fp, "#\n");

    fprintf( fp, "#  \"Pfs_ED_MLT\":       { \"DESCRIPTION\": \"Magnetic Local Time of Southern Footpoint in Eccentric Dipole Coordinates.\",\n");
    fprintf( fp, "#                               \"NAME\": \"Pfs_ED_MLT\",\n");
    fprintf( fp, "#                              \"TITLE\": \"ED MLT of Southern Footpoint\",\n");
    fprintf( fp, "#                              \"LABEL\": \"ED MLT (Hours)\",\n");
    fprintf( fp, "#                              \"UNITS\": \"Hours\",\n");
    fprintf( fp, "#                       \"START_COLUMN\": %d,\n", nCol++);
    fprintf( fp, "#                          \"VALID_MIN\": 0.0,\n");
    fprintf( fp, "#                          \"VALID_MAX\": 24.0,\n");
    fprintf( fp, "#                         \"FILL_VALUE\": -1e31 },\n");
    fprintf( fp, "#\n");

    fprintf( fp, "#  \"Bfs_geo\":          { \"DESCRIPTION\": \"Magnetic field vector at Southern Footpoint (in GEO coords).\",\n");
    fprintf( fp, "#                               \"NAME\": \"Bfs_geo\",\n");
    fprintf( fp, "#                              \"TITLE\": \"Magnetic field vector at Southern Footpoint (in GEO coords)\",\n");
    fprintf( fp, "#                              \"LABEL\": \"Bfs_geo (nT)\",\n");
    fprintf( fp, "#                              \"UNITS\": \"nT\",\n");
    fprintf( fp, "#                          \"DIMENSION\": [ 4 ],\n");
    fprintf( fp, "#                       \"START_COLUMN\": %d,\n", nCol); nCol += 4;
    fprintf( fp, "#                      \"ELEMENT_NAMES\": [ \"Bfs_geo_x\", \"Bfs_geo_y\", \"Bfs_geo_z\", \"Bfs_geo_mag\" ],\n");
    fprintf( fp, "#                     \"ELEMENT_LABELS\": [ \"Bfs_geo!Bx!N , nT\", \"Bfs_geo!By!N , nT\", \"Bfs_geo!Bz!N , nT\", \"Bfs_geo_mag , nT\" ],\n");
    //fprintf( fp, "#                          \"VALID_MIN\": -70000.0,\n");
    //fprintf( fp, "#                          \"VALID_MAX\":  70000.0,\n");
    fprintf( fp, "#                         \"FILL_VALUE\": -1e31 },\n");
    fprintf( fp, "#\n");

    fprintf( fp, "#  \"Bfs_gsm\":          { \"DESCRIPTION\": \"Magnetic field vector at Southern Footpoint (in GSM coords).\",\n");
    fprintf( fp, "#                               \"NAME\": \"Bfs_gsm\",\n");
    fprintf( fp, "#                              \"TITLE\": \"Magnetic field vector at Southern Footpoint (in GSM coords)\",\n");
    fprintf( fp, "#                              \"LABEL\": \"Bfs_gsm (nT)\",\n");
    fprintf( fp, "#                              \"UNITS\": \"nT\",\n");
    fprintf( fp, "#                          \"DIMENSION\": [ 4 ],\n");
    fprintf( fp, "#                       \"START_COLUMN\": %d,\n", nCol); nCol += 4;
    fprintf( fp, "#                      \"ELEMENT_NAMES\": [ \"Bfs_gsm_x\", \"Bfs_gsm_y\", \"Bfs_gsm_z\", \"Bfs_gsm_mag\" ],\n");
    fprintf( fp, "#                     \"ELEMENT_LABELS\": [ \"Bfs_gsm!Bx!N , nT\", \"Bfs_gsm!By!N , nT\", \"Bfs_gsm!Bz!N , nT\", \"Bfs_gsm_mag , nT\" ],\n");
    //fprintf( fp, "#                          \"VALID_MIN\": -70000.0,\n");
    //fprintf( fp, "#                          \"VALID_MAX\":  70000.0,\n");
    fprintf( fp, "#                         \"FILL_VALUE\": -1e31 },\n");
    fprintf( fp, "#\n");

    fprintf( fp, "# \"Loss_Cone_Alpha_s\": { \"DESCRIPTION\": \"Value of Southern Loss Cone angle. asin( sqrt(Bsc/Bfs) ). Bsc is the magntiude of B at spacecraft,"
                                                               " Bfs is the magnitude of B at the southern footpoint. Note, loss-cone height is taken to be equal"
                                                               " to the footpoint height which is %g km above the WGS84 geoid.\",\n", m->LstarInfo->mInfo->Lgm_LossConeHeight );
    fprintf( fp, "#                               \"NAME\": \"Loss_Cone_Alpha_s\",\n");
    fprintf( fp, "#                              \"TITLE\": \"Southern Loss Cone Angle\",\n");
    fprintf( fp, "#                              \"LABEL\": \"Southern Loss Cone Angle (Degrees)\",\n");
    fprintf( fp, "#                              \"UNITS\": \"Degrees\",\n");
    fprintf( fp, "#                       \"START_COLUMN\": %d,\n", nCol++);
    //fprintf( fp, "#                          \"VALID_MIN\":  0.0,\n");
    //fprintf( fp, "#                          \"VALID_MAX\": 90.0,\n");
    fprintf( fp, "#                         \"FILL_VALUE\": -1e31 },\n");
    fprintf( fp, "#\n");












    fprintf( fp, "#  \"Pmin_gsm\":         { \"DESCRIPTION\": \"Location of minimum-|B| point (in GSM coords).\",\n");
    fprintf( fp, "#                               \"NAME\": \"Pmin_gsm\",\n");
    fprintf( fp, "#                              \"TITLE\": \"Minimum-|B| point (in GSM Coordinates)\",\n");
    fprintf( fp, "#                              \"LABEL\": \"Pmin_gsm (R!BE)\",\n");
    fprintf( fp, "#                              \"UNITS\": \"R!BE!N\",\n");
    fprintf( fp, "#                          \"DIMENSION\": [ 3 ],\n");
    fprintf( fp, "#                       \"START_COLUMN\": %d,\n", nCol); nCol += 3;
    fprintf( fp, "#                      \"ELEMENT_NAMES\": [ \"Pmin_gsm_x\", \"Pmin_gsm_y\", \"Pmin_gsm_z\" ],\n");
    fprintf( fp, "#                     \"ELEMENT_LABELS\": [ \"Pmin_gsm!Bx!N , R!BE\", \"Pmin_gsm!By!N , R!BE\", \"Pmin_gsm!Bz!N , R!BE\" ],\n");
    //fprintf( fp, "#                          \"VALID_MIN\": -1000.0,\n");
    //fprintf( fp, "#                          \"VALID_MAX\":  1000.0,\n");
    fprintf( fp, "#                         \"FILL_VALUE\": -1e31 },\n");
    fprintf( fp, "#\n");


    fprintf( fp, "#  \"Bmin_gsm\":         { \"DESCRIPTION\": \"B-field at minimum-|B| point (in GSM coords).\",\n");
    fprintf( fp, "#                               \"NAME\": \"Bmin_gsm\",\n");
    fprintf( fp, "#                              \"TITLE\": \"Minimum-|B| (in GSM Coordinates)\",\n");
    fprintf( fp, "#                              \"LABEL\": \"Bmin_gsm (nT)\",\n");
    fprintf( fp, "#                              \"UNITS\": \"nT\",\n");
    fprintf( fp, "#                          \"DIMENSION\": [ 4 ],\n");
    fprintf( fp, "#                       \"START_COLUMN\": %d,\n", nCol); nCol += 4;
    fprintf( fp, "#                      \"ELEMENT_NAMES\": [ \"Bmin_gsm_x\", \"Bmin_gsm_y\", \"Bmin_gsm_z\", \"Bmin_gsm_mag\" ],\n");
    fprintf( fp, "#                     \"ELEMENT_LABELS\": [ \"Bmin_gsm!Bx!N , nT\", \"Bmin_gsm!By!N , nT\", \"Bmin_gsm!Bz!N , nT\", \"Bmin_gsm_mag , nT\" ],\n");
    //fprintf( fp, "#                          \"VALID_MIN\": -70000.0,\n");
    //fprintf( fp, "#                          \"VALID_MAX\":  70000.0,\n");
    fprintf( fp, "#                         \"FILL_VALUE\": -1e31 },\n");
    fprintf( fp, "#\n");

    fprintf( fp, "#  \"Lsimple\":          { \"DESCRIPTION\": \"Geocentric distance to Bmin point for FL threading vehicle (i.e. |Pmin|).\",\n");
    fprintf( fp, "#                               \"NAME\": \"Lsimple\",\n");
    fprintf( fp, "#                              \"TITLE\": \"Lsimple\",\n");
    fprintf( fp, "#                              \"LABEL\": \"Lsimple\",\n");
    fprintf( fp, "#                              \"UNITS\": \"Dimless\",\n");
    fprintf( fp, "#                       \"START_COLUMN\": %d,\n", nCol++);
    //fprintf( fp, "#                          \"VALID_MIN\":  0.0,\n");
    //fprintf( fp, "#                          \"VALID_MAX\": 90.0,\n");
    fprintf( fp, "#                         \"FILL_VALUE\": -1e31 },\n");
    fprintf( fp, "#\n");

    fprintf( fp, "#  \"InvLat\":           { \"DESCRIPTION\": \"Invariant latitude of vehicle computed from Lambda=acos(sqrt(1/Lsimple)).\",\n");
    fprintf( fp, "#                               \"NAME\": \"InvLat\",\n");
    fprintf( fp, "#                              \"TITLE\": \"Invariant Latitude\",\n");
    fprintf( fp, "#                              \"LABEL\": \"Invariant Latitude (Degrees)\",\n");
    fprintf( fp, "#                              \"UNITS\": \"Degrees\",\n");
    fprintf( fp, "#                       \"START_COLUMN\": %d,\n", nCol++);
    //fprintf( fp, "#                          \"VALID_MIN\":  0.0,\n");
    //fprintf( fp, "#                          \"VALID_MAX\": 90.0,\n");
    fprintf( fp, "#                         \"FILL_VALUE\": -1e31 },\n");
    fprintf( fp, "#\n");

    fprintf( fp, "#  \"Lm_eq\":            { \"DESCRIPTION\": \"McIlwain L of an eq. mirroring particle on same FL as vehicle (computed from L=Lm_eq, I=0, and Bm=|Bmin_gsm|, M=M_igrf).\",\n");
    fprintf( fp, "#                               \"NAME\": \"Lm_eq\",\n");
    fprintf( fp, "#                              \"TITLE\": \"McIlwain L for Equatorially Mirroring Particles\",\n");
    fprintf( fp, "#                              \"LABEL\": \"Lm for Eq. Mirr. Particles\",\n");
    fprintf( fp, "#                              \"UNITS\": \"Dimless\",\n");
    fprintf( fp, "#                       \"START_COLUMN\": %d,\n", nCol++);
    //fprintf( fp, "#                          \"VALID_MIN\":  0.0,\n");
    //fprintf( fp, "#                          \"VALID_MAX\": 90.0,\n");
    fprintf( fp, "#                         \"FILL_VALUE\": -1e31 },\n");
    fprintf( fp, "#\n");

    fprintf( fp, "#  \"InvLat_eq\":        { \"DESCRIPTION\": \"Invariant latitude of vehicle computed from Lambda=acos(sqrt(1.0/Lm_eq)).\",\n");
    fprintf( fp, "#                               \"NAME\": \"InvLat_eq\",\n");
    fprintf( fp, "#                              \"TITLE\": \"Lm_eq Invariant Latitude\",\n");
    fprintf( fp, "#                              \"LABEL\": \"Lm_eq Invariant Latitude (Degrees)\",\n");
    fprintf( fp, "#                              \"UNITS\": \"Degrees\",\n");
    fprintf( fp, "#                       \"START_COLUMN\": %d,\n", nCol++);
    //fprintf( fp, "#                          \"VALID_MIN\":  0.0,\n");
    //fprintf( fp, "#                          \"VALID_MAX\": 90.0,\n");
    fprintf( fp, "#                         \"FILL_VALUE\": -1e31 },\n");
    fprintf( fp, "#\n");


    fprintf( fp, "#  \"BoverBeq\":         { \"DESCRIPTION\": \"Magntiude of Bsc over magnitude of Bmin.\",\n");
    fprintf( fp, "#                               \"NAME\": \"BoverBeq\",\n");
    fprintf( fp, "#                              \"TITLE\": \"BoverBeq\",\n");
    fprintf( fp, "#                              \"LABEL\": \"BoverBeq\",\n");
    fprintf( fp, "#                              \"UNITS\": \"Dimless\",\n");
    fprintf( fp, "#                       \"START_COLUMN\": %d,\n", nCol++);
    //fprintf( fp, "#                          \"VALID_MIN\":  0.0,\n");
    //fprintf( fp, "#                          \"VALID_MAX\": 90.0,\n");
    fprintf( fp, "#                         \"FILL_VALUE\": -1e31 },\n");
    fprintf( fp, "#\n");

    fprintf( fp, "#  \"MlatFromBoverBeq\": { \"DESCRIPTION\": \"Dipole latitude where (B/Beq)_dipole == BoverBeq.\",\n");
    fprintf( fp, "#                               \"NAME\": \"MlatFromBoverBeq\",\n");
    fprintf( fp, "#                              \"TITLE\": \"MlatFromBoverBeq\",\n");
    fprintf( fp, "#                              \"LABEL\": \"MlatFromBoverBeq (Degrees)\",\n");
    fprintf( fp, "#                              \"UNITS\": \"Degrees\",\n");
    fprintf( fp, "#                       \"START_COLUMN\": %d,\n", nCol++);
    //fprintf( fp, "#                          \"VALID_MIN\":  0.0,\n");
    //fprintf( fp, "#                          \"VALID_MAX\": 90.0,\n");
    fprintf( fp, "#                         \"FILL_VALUE\": -1e31 },\n");
    fprintf( fp, "#\n");



    fprintf( fp, "#  \"M_used\":           { \"DESCRIPTION\": \"The magnetic dipole moment that was used to convert magnetic flux to L*. In units of nT.\",\n");
    fprintf( fp, "#                               \"NAME\": \"M_used\",\n");
    fprintf( fp, "#                              \"TITLE\": \"M_used\",\n");
    fprintf( fp, "#                              \"LABEL\": \"M_used (nT)\",\n");
    fprintf( fp, "#                              \"UNITS\": \"nT\",\n");
    fprintf( fp, "#                       \"START_COLUMN\": %d,\n", nCol++);
    //fprintf( fp, "#                          \"VALID_MIN\":     0.0,\n");
    //fprintf( fp, "#                          \"VALID_MAX\": 70000.0,\n");
    fprintf( fp, "#                         \"FILL_VALUE\": -1e31 },\n");
    fprintf( fp, "#\n");

    fprintf( fp, "#  \"M_ref\":            { \"DESCRIPTION\": \"The fixed reference magnetic dipole moment for converting magnetic flux to L*. In units of nT.\",\n");
    fprintf( fp, "#                               \"NAME\": \"M_ref\",\n");
    fprintf( fp, "#                              \"TITLE\": \"M_ref\",\n");
    fprintf( fp, "#                              \"LABEL\": \"M_ref (nT)\",\n");
    fprintf( fp, "#                              \"UNITS\": \"nT\",\n");
    fprintf( fp, "#                       \"START_COLUMN\": %d,\n", nCol++);
    //fprintf( fp, "#                          \"VALID_MIN\":     0.0,\n");
    //fprintf( fp, "#                          \"VALID_MAX\": 70000.0,\n");
    fprintf( fp, "#                         \"FILL_VALUE\": -1e31 },\n");
    fprintf( fp, "#\n");

    fprintf( fp, "#  \"M_igrf\":           { \"DESCRIPTION\": \"Time-dependant magnetic dipole moment (probably shouldn't be used for converting magnetic flux to L*, but it sometimes is). In units of nT.\",\n");
    fprintf( fp, "#                               \"NAME\": \"M_igrf\",\n");
    fprintf( fp, "#                              \"TITLE\": \"M_igrf\",\n");
    fprintf( fp, "#                              \"LABEL\": \"M_igrf (nT)\",\n");
    fprintf( fp, "#                              \"UNITS\": \"nT\",\n");
    fprintf( fp, "#                       \"START_COLUMN\": %d,\n", nCol++);
    //fprintf( fp, "#                          \"VALID_MIN\":     0.0,\n");
    //fprintf( fp, "#                          \"VALID_MAX\": 70000.0,\n");
    if ( m->nAlpha > 0 ){
        fprintf( fp, "#                         \"FILL_VALUE\": -1e31 },\n");
    } else {
        fprintf( fp, "#                         \"FILL_VALUE\": -1e31 }\n");
    }
    fprintf( fp, "#\n");






    if ( m->nAlpha > 0 ) {
        fprintf( fp, "#  \"Lstar\":            { \"DESCRIPTION\": \"Generalized Roederer L-shell value (also known as L*).\",\n");
        fprintf( fp, "#                               \"NAME\": \"Lstar\",\n");
        fprintf( fp, "#                              \"TITLE\": \"L!A*!N\",\n");
        fprintf( fp, "#                              \"LABEL\": \"L!A*!N\",\n");
        fprintf( fp, "#                              \"UNITS\": \"Dimensionless\",\n");
        fprintf( fp, "#                          \"DIMENSION\": [ %d ],\n", m->nAlpha );
        fprintf( fp, "#                       \"START_COLUMN\": %d,\n", nCol); nCol += m->nAlpha;
        fprintf( fp, "#                      \"ELEMENT_NAMES\": [ ");
        for (i=0; i<m->nAlpha-1; i++) fprintf(fp, "\"Lstar_%g\", ", m->Alpha[i] );
        fprintf(fp, "\"Lstar_%g\" ],\n", m->Alpha[i] ); 
        fprintf( fp, "#                     \"ELEMENT_LABELS\": [ ");
        for (i=0; i<m->nAlpha-1; i++) fprintf(fp, "\"L!A*!N %g!Ao!N\", ", m->Alpha[i] );
        fprintf(fp, "\"L!A*!N %g!Ao!N\" ],\n", m->Alpha[i] ); 
        fprintf( fp, "#                           \"DEPEND_1\": \"Alpha\",\n");
        //fprintf( fp, "#                          \"VALID_MIN\": 0.0,\n");
        //fprintf( fp, "#                          \"VALID_MAX\": 1000.0,\n");
        fprintf( fp, "#                         \"FILL_VALUE\": -1e31 },\n");
        fprintf( fp, "#\n");

        fprintf( fp, "#  \"DriftShellType\":   { \"DESCRIPTION\": \"Type of Drift Shell (e.g. %d=CLOSED, %d=CLOSED_SHABANSKY, %d=OPEN, %d=OPEN_SHABANSKY)\",\n", LGM_DRIFT_ORBIT_CLOSED, LGM_DRIFT_ORBIT_CLOSED_SHABANSKY, LGM_DRIFT_ORBIT_OPEN, LGM_DRIFT_ORBIT_OPEN_SHABANSKY);
        fprintf( fp, "#                               \"NAME\": \"DriftShellType\",\n");
        fprintf( fp, "#                              \"TITLE\": \"Drift Shell Type\",\n");
        fprintf( fp, "#                              \"LABEL\": \"DS Type\",\n");
        fprintf( fp, "#                              \"UNITS\": \"Dimensionless\",\n");
        fprintf( fp, "#                          \"DIMENSION\": [ %d ],\n", m->nAlpha );
        fprintf( fp, "#                       \"START_COLUMN\": %d,\n", nCol); nCol += m->nAlpha;
        fprintf( fp, "#                      \"ELEMENT_NAMES\": [ ");
        for (i=0; i<m->nAlpha-1; i++) fprintf(fp, "\"DSType_%g\", ", m->Alpha[i] );
        fprintf(fp, "\"DSType_%g\" ],\n", m->Alpha[i] ); 
        fprintf( fp, "#                     \"ELEMENT_LABELS\": [ ");
        for (i=0; i<m->nAlpha-1; i++) fprintf(fp, "\"DSType %g!Ao!N\", ", m->Alpha[i] );
        fprintf(fp, "\"DSType %g!Ao!N\" ],\n", m->Alpha[i] ); 
        fprintf( fp, "#                           \"DEPEND_1\": \"Alpha\",\n");
        //fprintf( fp, "#                          \"VALID_MIN\": 0.0,\n");
        //fprintf( fp, "#                          \"VALID_MAX\": 1000.0,\n");
        fprintf( fp, "#                         \"FILL_VALUE\": -1e31 },\n");
        fprintf( fp, "#\n");


        fprintf( fp, "#  \"L\":                { \"DESCRIPTION\": \"McIlwain L-shell value.\",\n");
        fprintf( fp, "#                               \"NAME\": \"L\",\n");
        fprintf( fp, "#                              \"TITLE\": \"L\",\n");
        fprintf( fp, "#                              \"LABEL\": \"L\",\n");
        fprintf( fp, "#                              \"UNITS\": \"Dimensionless\",\n");
        fprintf( fp, "#                          \"DIMENSION\": [ %d ],\n", m->nAlpha );
        fprintf( fp, "#                       \"START_COLUMN\": %d,\n", nCol); nCol += m->nAlpha;
        fprintf( fp, "#                      \"ELEMENT_NAMES\": [ ");
        for (i=0; i<m->nAlpha-1; i++) fprintf(fp, "\"L_%g\", ", m->Alpha[i] );
        fprintf(fp, "\"L_%g\" ],\n", m->Alpha[i] ); 
        fprintf( fp, "#                     \"ELEMENT_LABELS\": [ ");
        for (i=0; i<m->nAlpha-1; i++) fprintf(fp, "\"L %g!Ao!N\", ", m->Alpha[i] );
        fprintf(fp, "\"L %g!Ao!N\" ],\n", m->Alpha[i] ); 
        fprintf( fp, "#                           \"DEPEND_1\": \"Alpha\",\n");
        //fprintf( fp, "#                          \"VALID_MIN\": 0.0,\n");
        //fprintf( fp, "#                          \"VALID_MAX\": 1000.0,\n");
        fprintf( fp, "#                         \"FILL_VALUE\": -1e31 },\n");
        fprintf( fp, "#\n");


        fprintf( fp, "#  \"Bm\":               { \"DESCRIPTION\": \"Magnetic field strength at mirror points for each pitch angle.\",\n");
        fprintf( fp, "#                               \"NAME\": \"Bm\",\n");
        fprintf( fp, "#                              \"TITLE\": \"Mirror Magnetic Field Strength\",\n");
        fprintf( fp, "#                              \"LABEL\": \"Bm (nT)\",\n");
        fprintf( fp, "#                              \"UNITS\": \"nT\",\n");
        fprintf( fp, "#                          \"DIMENSION\": [ %d ],\n", m->nAlpha );
        fprintf( fp, "#                       \"START_COLUMN\": %d,\n", nCol); nCol += m->nAlpha;
        fprintf( fp, "#                      \"ELEMENT_NAMES\": [ ");
        for (i=0; i<m->nAlpha-1; i++) fprintf(fp, "\"Bm_%g\", ", m->Alpha[i] );
        fprintf(fp, "\"Bm_%g\" ],\n", m->Alpha[i] ); 
        fprintf( fp, "#                     \"ELEMENT_LABELS\": [ ");
        for (i=0; i<m->nAlpha-1; i++) fprintf(fp, "\"B!Dm!N %g!Ao!N\", ", m->Alpha[i] );
        fprintf(fp, "\"B!Dm!N %g!Ao!N\" ],\n", m->Alpha[i] ); 
        fprintf( fp, "#                           \"DEPEND_1\": \"Alpha\",\n");
        //fprintf( fp, "#                          \"VALID_MIN\": 0.0,\n");
        //fprintf( fp, "#                          \"VALID_MAX\": 70000.0,\n");
        fprintf( fp, "#                         \"FILL_VALUE\": -1e31 },\n");
        fprintf( fp, "#\n");


        fprintf( fp, "#  \"I\":                { \"DESCRIPTION\": \"Integral invariant for each pitch angle.\",\n");
        fprintf( fp, "#                               \"NAME\": \"I\",\n");
        fprintf( fp, "#                              \"TITLE\": \"Integral Invariant\",\n");
        fprintf( fp, "#                              \"LABEL\": \"I (R!BE)\",\n");
        fprintf( fp, "#                              \"UNITS\": \"R!BE!N\",\n");
        fprintf( fp, "#                          \"DIMENSION\": [ %d ],\n", m->nAlpha );
        fprintf( fp, "#                       \"START_COLUMN\": %d,\n", nCol); nCol += m->nAlpha;
        fprintf( fp, "#                      \"ELEMENT_NAMES\": [ ");
        for (i=0; i<m->nAlpha-1; i++) fprintf(fp, "\"I_%g\", ", m->Alpha[i] );
        fprintf(fp, "\"I_%g\" ],\n", m->Alpha[i] ); 
        fprintf( fp, "#                     \"ELEMENT_LABELS\": [ ");
        for (i=0; i<m->nAlpha-1; i++) fprintf(fp, "\"I %g!Ao!N\", ", m->Alpha[i] );
        fprintf(fp, "\"I %g!Ao!N\" ],\n", m->Alpha[i] ); 
        fprintf( fp, "#                           \"DEPEND_1\": \"Alpha\",\n");
        //fprintf( fp, "#                          \"VALID_MIN\": 0.0,\n");
        //fprintf( fp, "#                          \"VALID_MAX\": 1000.0,\n");
        fprintf( fp, "#                         \"FILL_VALUE\": -1e31 },\n");
        fprintf( fp, "#\n");

        fprintf( fp, "#  \"K\":                { \"DESCRIPTION\": \"Second Invariant ( I*sqrt(Bm) ) for each pitch angle.\",\n");
        fprintf( fp, "#                               \"NAME\": \"K\",\n");
        fprintf( fp, "#                              \"TITLE\": \"Second Invariant\",\n");
        fprintf( fp, "#                              \"LABEL\": \"K (R!BE!N G!A0.5!N)\",\n");
        fprintf( fp, "#                              \"UNITS\": \"R!BE!N G!A0.5!N\",\n");
        fprintf( fp, "#                          \"DIMENSION\": [ %d ],\n", m->nAlpha );
        fprintf( fp, "#                       \"START_COLUMN\": %d,\n", nCol); nCol += m->nAlpha;
        fprintf( fp, "#                      \"ELEMENT_NAMES\": [ ");
        for (i=0; i<m->nAlpha-1; i++) fprintf(fp, "\"K_%g\", ", m->Alpha[i] );
        fprintf(fp, "\"K_%g\" ],\n", m->Alpha[i] ); 
        fprintf( fp, "#                     \"ELEMENT_LABELS\": [ ");
        for (i=0; i<m->nAlpha-1; i++) fprintf(fp, "\"K %g!Ao!N\", ", m->Alpha[i] );
        fprintf(fp, "\"K %g!Ao!N\" ],\n", m->Alpha[i] ); 
        fprintf( fp, "#                           \"DEPEND_1\": \"Alpha\",\n");
        //fprintf( fp, "#                          \"VALID_MIN\": 0.0,\n");
        //fprintf( fp, "#                          \"VALID_MAX\": 100000.0,\n");
        fprintf( fp, "#                         \"FILL_VALUE\": -1e31 },\n");
        fprintf( fp, "#\n");

        fprintf( fp, "#  \"Sb\":               { \"DESCRIPTION\": \"Sb Integral for each pitch angle.\",\n");
        fprintf( fp, "#                               \"NAME\": \"Sb\",\n");
        fprintf( fp, "#                              \"TITLE\": \"Sb Integral\",\n");
        fprintf( fp, "#                              \"LABEL\": \"Sb (R!BE!N)\",\n");
        fprintf( fp, "#                              \"UNITS\": \"R!BE!N\",\n");
        fprintf( fp, "#                          \"DIMENSION\": [ %d ],\n", m->nAlpha );
        fprintf( fp, "#                       \"START_COLUMN\": %d,\n", nCol); nCol += m->nAlpha;
        fprintf( fp, "#                      \"ELEMENT_NAMES\": [ ");
        for (i=0; i<m->nAlpha-1; i++) fprintf(fp, "\"Sb_%g\", ", m->Alpha[i] );
        fprintf(fp, "\"Sb_%g\" ],\n", m->Alpha[i] ); 
        fprintf( fp, "#                     \"ELEMENT_LABELS\": [ ");
        for (i=0; i<m->nAlpha-1; i++) fprintf(fp, "\"Sb %g!Ao!N\", ", m->Alpha[i] );
        fprintf(fp, "\"Sb %g!Ao!N\" ],\n", m->Alpha[i] ); 
        fprintf( fp, "#                           \"DEPEND_1\": \"Alpha\",\n");
        //fprintf( fp, "#                          \"VALID_MIN\": 0.0,\n");
        //fprintf( fp, "#                          \"VALID_MAX\": 100000.0,\n");
        fprintf( fp, "#                         \"FILL_VALUE\": -1e31 },\n");
        fprintf( fp, "#\n");

        fprintf( fp, "#  \"Tb\":               { \"DESCRIPTION\": \"Bounce period for 1 MeV electrons.\",\n");
        fprintf( fp, "#                               \"NAME\": \"Tb\",\n");
        fprintf( fp, "#                              \"TITLE\": \"Bounce period for 1 MeV electrons\",\n");
        fprintf( fp, "#                              \"LABEL\": \"Tb (s)\",\n");
        fprintf( fp, "#                              \"UNITS\": \"s\",\n");
        fprintf( fp, "#                          \"DIMENSION\": [ %d ],\n", m->nAlpha );
        fprintf( fp, "#                       \"START_COLUMN\": %d,\n", nCol); nCol += m->nAlpha;
        fprintf( fp, "#                      \"ELEMENT_NAMES\": [ ");
        for (i=0; i<m->nAlpha-1; i++) fprintf(fp, "\"Tb_%g\", ", m->Alpha[i] );
        fprintf(fp, "\"Tb_%g\" ],\n", m->Alpha[i] ); 
        fprintf( fp, "#                     \"ELEMENT_LABELS\": [ ");
        for (i=0; i<m->nAlpha-1; i++) fprintf(fp, "\"Tb %g!Ao!N\", ", m->Alpha[i] );
        fprintf(fp, "\"Tb %g!Ao!N\" ],\n", m->Alpha[i] ); 
        fprintf( fp, "#                           \"DEPEND_1\": \"Alpha\",\n");
        //fprintf( fp, "#                          \"VALID_MIN\": 0.0,\n");
        //fprintf( fp, "#                          \"VALID_MAX\": 100000.0,\n");
        fprintf( fp, "#                         \"FILL_VALUE\": -1e31 },\n");
        fprintf( fp, "#\n");

        fprintf( fp, "#  \"Kappa\":            { \"DESCRIPTION\": \"Kappa parameter for 1MeV electrons -- sqrt( (Minimum Radius of Curvature)/Maximum gyroradius)) (see Büchner, J., and L. M. Zelenyi (1989), Regular and Chaotic Charged Particle Motion in Magnetotaillike Field Reversals, 1. Basic Theory of Trapped Motion, J. Geophys. Res., 94(A9), 11,821–11,842, doi:10.1029/JA094iA09p11821.\",\n");
        fprintf( fp, "#                               \"NAME\": \"Kappa\",\n");
        fprintf( fp, "#                              \"TITLE\": \"Kappa parameter for 1MeV electrons\",\n");
        fprintf( fp, "#                              \"LABEL\": \"Kappa\",\n");
        fprintf( fp, "#                              \"UNITS\": \"dimensionless\",\n");
        fprintf( fp, "#                          \"DIMENSION\": [ %d ],\n", m->nAlpha );
        fprintf( fp, "#                       \"START_COLUMN\": %d,\n", nCol); nCol += m->nAlpha;
        fprintf( fp, "#                      \"ELEMENT_NAMES\": [ ");
        for (i=0; i<m->nAlpha-1; i++) fprintf(fp, "\"Kappa_%g\", ", m->Alpha[i] );
        fprintf(fp, "\"Kappa_%g\" ],\n", m->Alpha[i] ); 
        fprintf( fp, "#                     \"ELEMENT_LABELS\": [ ");
        for (i=0; i<m->nAlpha-1; i++) fprintf(fp, "\"Kappa %g!Ao!N\", ", m->Alpha[i] );
        fprintf(fp, "\"Kappa %g!Ao!N\" ],\n", m->Alpha[i] ); 
        fprintf( fp, "#                           \"DEPEND_1\": \"Alpha\",\n");
        //fprintf( fp, "#                          \"VALID_MIN\": 0.0,\n");
        //fprintf( fp, "#                          \"VALID_MAX\": 100000.0,\n");
        fprintf( fp, "#                         \"FILL_VALUE\": -1e31 }\n");
//        fprintf( fp, "#\n");

    }


    fprintf( fp, "# } end JSON\n");
    fprintf( fp, "#\n");



    fprintf( fp, "#\n");
    fprintf( fp, "#  \"Spacecraft\":       { \"DESCRIPTION\": \"Spacecraft identifier.\",\n");
    fprintf( fp, "#                        \"COMMON_NAME\": \"%s\",\n", Spacecraft);
    fprintf( fp, "#                          \"ID_NUMBER\": \"%d\",\n", IdNumber);
    fprintf( fp, "#                          \"INT_DESIG\": \"%s\",\n", IntDesig);
    fprintf( fp, "#                      \"SPICE_BODY_ID\": \"%d\"\n",  SpiceBody);
    fprintf( fp, "#               \"SPICE_KERNELS_LOADED\": \"%s\"\n", "SPICE_KERNEL_FILES_LOADED" );
    fprintf( fp, "#  },\n");

    fprintf( fp, "#  \"File\":              { \"DESCRIPTION\": \"Description of file contents. The format for ElapsedTime is DDD:HH:MM:SS.\",\n");
    fprintf( fp, "#                       \"CreationTime\": \"%02d:%02d:%02d UTC  %s %02d %4d\",\n", HH, MM, SS, sMonth[Month], Day, Year);
    fprintf( fp, "#                        \"ElapsedTime\": \"%s\",\n", "ELAPSED_TIME" );
    fprintf( fp, "#                          \"CreatedBy\": \"%s\",\n", getenv( "USER" ) );
    fprintf( fp, "#                          \"CreatedOn\": \"%s\",\n", getenv("HOSTNAME") );
    fprintf( fp, "#                        \"CommandLine\": \"%s\"\n", CmdLine );
    fprintf( fp, "#  },\n");

    fprintf( fp, "#\n");
    fprintf( fp, "#\n");



    // column header
    fprintf( fp, "%91s",  "#  +------------------------------------ Date and Time -----------------------------------+" );
    fprintf( fp, " %13s",  " +- TiltAng +" );
    fprintf( fp, " %11s",  " +- InOut +" );
    fprintf( fp, " %11s",  " +-OrbNum-+" );
    fprintf( fp, " %41s",  " +--- Geocentric Geographic Coords --+" );
    fprintf( fp, " %41s",  " +---- Geodetic Geographic Coords ---+" );
    fprintf( fp, " %41s",  " +--------- GSM Coordinates ---------+" );
    fprintf( fp, " %41s",  " +---------- SM Coordinates ---------+" );
    fprintf( fp, " %41s",  " +------- GEI 2000 Coordinates ------+" );
    fprintf( fp, " %41s",  " +---------- GSE Coordinates --------+" );
    fprintf( fp, " %41s",  " +----------------- CD MLAT/MLON/MLT/R ---------------+" );
    fprintf( fp, " %41s",  " +----------------- ED MLAT/MLON/MLT/R ---------------+" );
    fprintf( fp, " %13s",  " +-Int Model-+" );
    fprintf( fp, " %13s",  " +-Ext Model-+" );
    fprintf( fp, " %6s",   " +-Kp-+" );
    fprintf( fp, " %7s",   " +-Dst-+" );
    fprintf( fp, " %51s",  " +--------- Magnetic Field at SpaceCraft ---------+" );
    fprintf( fp, " %29s",  " +----- Field Line Type ----+" );
    fprintf( fp, " %16s",  " +-S_sc_to_pfn-+" );
    fprintf( fp, " %16s",  " +-S_sc_to_pfs-+" );
    fprintf( fp, " %16s",  " +S_pfs_to_Bmin+" );
    fprintf( fp, " %16s",  " +-S_sc_to_Bmin+" );
    fprintf( fp, " %16s",  " +-- S_total --+" );
    fprintf( fp, " %16s",  " +-- d2B_ds2 --+" );
    fprintf( fp, " %16s",  " +---- Sb0 ----+" );
    fprintf( fp, " %16s",  " +- RadOfCurv -+" );
    fprintf( fp, " %38s",  " +---- North Mag. Footpoint GSM -----+" );
    fprintf( fp, " %38s",  " +- North Mag. Footpoint Geographic -+" );
    fprintf( fp, " %38s",  " +-- North Mag. Footpoint Geodetic --+" );
    fprintf( fp, " %38s",  " +--- North Mag. Footpoint CD_MAG  --+" );
    fprintf( fp, " %38s",  " +--- North Mag. Footpoint ED_MAG  --+" );
    fprintf( fp, " %51s",  " +---- Mag. Field at North Mag. Footpoint GEO ----+" );
    fprintf( fp, " %51s",  " +---- Mag. Field at North Mag. Footpoint GSM ----+" );
    fprintf( fp, " %12s",  " +-N.L.Cone+" );
    fprintf( fp, " %38s",  " +---- South Mag. Footpoint GSM -----+" );
    fprintf( fp, " %38s",  " +- South Mag. Footpoint Geographic -+" );
    fprintf( fp, " %38s",  " +-- South Mag. Footpoint Geodetic --+" );
    fprintf( fp, " %38s",  " +--- South Mag. Footpoint CD_MAG  --+" );
    fprintf( fp, " %38s",  " +--- South Mag. Footpoint ED_MAG  --+" );
    fprintf( fp, " %51s",  " +---- Mag. Field at South Mag. Footpoint GEO ----+" );
    fprintf( fp, " %51s",  " +---- Mag. Field at South Mag. Footpoint GSM ----+" );
    fprintf( fp, " %12s",  " +-S.L.Cone+" );
    fprintf( fp, " %38s",  " +----- Minimum |B| Point GSM -------+" );
    fprintf( fp, " %51s",  " +---- Magnetic Field at Minimum |B| Point    ----+" );
    fprintf( fp, " %12s",  " +-Lsimple-+" );
    fprintf( fp, " %12s",  " +-InvLat -+" );
    fprintf( fp, " %12s",  " +- Lm_eq -+" );
    fprintf( fp, " %12s",  " + InvLatEq+" );
    fprintf( fp, " %12s",  " +-BoverBeq+" );
    fprintf( fp, " %12s",  " +-MlatBBeq+" );
    fprintf( fp, " %12s",  " +- M_used +" );
    fprintf( fp, " %12s",  " +- M_ref -+" );
    fprintf( fp, " %12s",  " +- M_igrf +" );
    if (m->nAlpha>0){
        n = 13*m->nAlpha;
        LGM_ARRAY_1D( Str2, n+10, char );

        // Lstar header
        p = Str2; *p++ = ' '; *p++ = '+'; for (i=0; i<n; i++) *p++ = '-'; *p++ = '+'; *p++ = '\0';
        // add text in center
        if ( m->nAlpha > 1 ) {
            sprintf( TextStr, " Lstar (%d Pitch Angles) ", m->nAlpha ); tsl = strlen( TextStr );
        } else {
            sprintf( TextStr, " Lstar " ); tsl = strlen( TextStr );
        }
        n2 = (n+3 - tsl)/2; if (n2<0) n2 = 0;
        for (p = Str2+n2, i=0; i<tsl; i++) *p++ = TextStr[i];
        fprintf( fp, " %*s", n+3,  Str2 );

        // DS Type header
        p = Str2; *p++ = ' '; *p++ = '+'; for (i=0; i<n; i++) *p++ = '-'; *p++ = '+'; *p++ = '\0';
        // add text in center
        if ( m->nAlpha > 1 ) {
            sprintf( TextStr, " Drift Shell Type (%d Pitch Angles) ", m->nAlpha ); tsl = strlen( TextStr );
        } else {
            sprintf( TextStr, " DSType " ); tsl = strlen( TextStr );
        }
        n2 = (n+3 - tsl)/2; if (n2<0) n2 = 0;
        for (p = Str2+n2, i=0; i<tsl; i++) *p++ = TextStr[i];
        fprintf( fp, " %*s", n+3,  Str2 );

        // L header
        p = Str2; *p++ = ' '; *p++ = '+'; for (i=0; i<n; i++) *p++ = '-'; *p++ = '+'; *p++ = '\0';
        // add text in center
        if ( m->nAlpha > 1 ) {
            sprintf( TextStr, " L (%d Pitch Angles) ", m->nAlpha ); tsl = strlen( TextStr );
        } else {
            sprintf( TextStr, " L " ); tsl = strlen( TextStr );
        }
        n2 = (n+3 - tsl)/2; if (n2<0) n2 = 0;
        for (p = Str2+n2, i=0; i<tsl; i++) *p++ = TextStr[i];
        fprintf( fp, " %*s", n+3,  Str2 );

        // Bm header
        p = Str2; *p++ = ' '; *p++ = '+'; for (i=0; i<n; i++) *p++ = '-'; *p++ = '+'; *p++ = '\0';
        // add text in center
        if ( m->nAlpha > 1 ) {
            sprintf( TextStr, " Bm (%d Pitch Angles) ", m->nAlpha ); tsl = strlen( TextStr );
        } else {
            sprintf( TextStr, " Bm " ); tsl = strlen( TextStr );
        }
        n2 = (n+3 - tsl)/2; if (n2<0) n2 = 0;
        for (p = Str2+n2, i=0; i<tsl; i++) *p++ = TextStr[i];
        fprintf( fp, " %*s", n+3,  Str2 );

        // I header
        p = Str2; *p++ = ' '; *p++ = '+'; for (i=0; i<n; i++) *p++ = '-'; *p++ = '+'; *p++ = '\0';
        // add text in center
        if ( m->nAlpha > 1 ) {
            sprintf( TextStr, " I (%d Pitch Angles) ", m->nAlpha ); tsl = strlen( TextStr );
        } else {
            sprintf( TextStr, " I " ); tsl = strlen( TextStr );
        }
        n2 = (n+3 - tsl)/2; if (n2<0) n2 = 0;
        for (p = Str2+n2, i=0; i<tsl; i++) *p++ = TextStr[i];
        fprintf( fp, " %*s", n+3,  Str2 );

        // K header
        p = Str2; *p++ = ' '; *p++ = '+'; for (i=0; i<n; i++) *p++ = '-'; *p++ = '+'; *p++ = '\0';
        // add text in center
        if ( m->nAlpha > 1 ) {
            sprintf( TextStr, " K (%d Pitch Angles) ", m->nAlpha ); tsl = strlen( TextStr );
        } else {
            sprintf( TextStr, " K " ); tsl = strlen( TextStr );
        }
        n2 = (n+3 - tsl)/2; if (n2<0) n2 = 0;
        for (p = Str2+n2, i=0; i<tsl; i++) *p++ = TextStr[i];
        fprintf( fp, " %*s", n+3,  Str2 );

        // Sb header
        p = Str2; *p++ = ' '; *p++ = '+'; for (i=0; i<n; i++) *p++ = '-'; *p++ = '+'; *p++ = '\0';
        // add text in center
        if ( m->nAlpha > 1 ) {
            sprintf( TextStr, " Sb (%d Pitch Angles) ", m->nAlpha ); tsl = strlen( TextStr );
        } else {
            sprintf( TextStr, " Sb " ); tsl = strlen( TextStr );
        }
        n2 = (n+3 - tsl)/2; if (n2<0) n2 = 0;
        for (p = Str2+n2, i=0; i<tsl; i++) *p++ = TextStr[i];
        fprintf( fp, " %*s", n+3,  Str2 );

        // Kappa header
        p = Str2; *p++ = ' '; *p++ = '+'; for (i=0; i<n; i++) *p++ = '-'; *p++ = '+'; *p++ = '\0';
        // add text in center
        if ( m->nAlpha > 1 ) {
            sprintf( TextStr, " Kappa (%d Pitch Angles) ", m->nAlpha ); tsl = strlen( TextStr );
        } else {
            sprintf( TextStr, " Kappa " ); tsl = strlen( TextStr );
        }
        n2 = (n+3 - tsl)/2; if (n2<0) n2 = 0;
        for (p = Str2+n2, i=0; i<tsl; i++) *p++ = TextStr[i];
        fprintf( fp, " %*s", n+3,  Str2 );

        // Tb header
        p = Str2; *p++ = ' '; *p++ = '+'; for (i=0; i<n; i++) *p++ = '-'; *p++ = '+'; *p++ = '\0';
        // add text in center
        if ( m->nAlpha > 1 ) {
            sprintf( TextStr, " Tb (%d Pitch Angles) ", m->nAlpha ); tsl = strlen( TextStr );
        } else {
            sprintf( TextStr, " Tb " ); tsl = strlen( TextStr );
        }
        n2 = (n+3 - tsl)/2; if (n2<0) n2 = 0;
        for (p = Str2+n2, i=0; i<tsl; i++) *p++ = TextStr[i];
        fprintf( fp, " %*s", n+3,  Str2 );

        LGM_ARRAY_1D_FREE( Str2 );
    }

    fprintf(fp, "\n");


    // column header
    fprintf( fp, "# %25s", "Time" );
    fprintf( fp, " %10s", "Date" );
    fprintf( fp, " %5s",  "DOY" );
    fprintf( fp, " %13s", "UTC" );
    fprintf( fp, " %16s", "Julian Date" );
    fprintf( fp, " %15s", "GPS Time" );
    fprintf( fp, " %13s", "    " );
    fprintf( fp, " %11s", "    " );
    fprintf( fp, " %11s", "    " );

    fprintf( fp, " %13s", "Xgeo" );
    fprintf( fp, " %13s", "Ygeo" );
    fprintf( fp, " %13s", "Zgeo" );

    fprintf( fp, " %13s", "GeodLat" );
    fprintf( fp, " %13s", "GeodLon" );
    fprintf( fp, " %13s", "GeodHeight" );

    fprintf( fp, " %13s", "Xgsm" );
    fprintf( fp, " %13s", "Ygsm" );
    fprintf( fp, " %13s", "Zgsm" );

    fprintf( fp, " %13s", "Xsm" );
    fprintf( fp, " %13s", "Ysm" );
    fprintf( fp, " %13s", "Zsm" );

    fprintf( fp, " %13s", "Xgei" );
    fprintf( fp, " %13s", "Ygei" );
    fprintf( fp, " %13s", "Zgei" );

    fprintf( fp, " %13s", "Xgse" );
    fprintf( fp, " %13s", "Ygse" );
    fprintf( fp, " %13s", "Zgse" );

    fprintf( fp, " %13s", "CD_MLAT" );
    fprintf( fp, " %13s", "CD_MLON" );
    fprintf( fp, " %13s", "CD_MLT" );
    fprintf( fp, " %13s", "CD_R" );

    fprintf( fp, " %13s", "ED_MLAT" );
    fprintf( fp, " %13s", "ED_MLON" );
    fprintf( fp, " %13s", "ED_MLT" );
    fprintf( fp, " %13s", "ED_R" );

    fprintf( fp, " %14s",  "Int Model" );
    fprintf( fp, " %14s",  "Ext Model" );
    fprintf( fp, " %7s",   "Kp" );
    fprintf( fp, " %8s",   "Dst" );

    fprintf( fp, " %12s", "Bsc_x" ); // Bsc  gsm
    fprintf( fp, " %12s", "Bsc_y" );
    fprintf( fp, " %12s", "Bsc_z" );
    fprintf( fp, " %12s", "Bsc" );

    fprintf( fp, " %29s",  "Field Line Type" );

    fprintf( fp, " %16s",  "S_sc_to_pfn" );
    fprintf( fp, " %16s",  "S_sc_to_pfs" );
    fprintf( fp, " %16s",  "S_pfs_to_Bmin" );
    fprintf( fp, " %16s",  "S_sc_to_Bmin" );
    fprintf( fp, " %16s",  "S_total" );

    fprintf( fp, " %16s",  "d2B_ds2" );
    fprintf( fp, " %16s",  "Sb0" );
    fprintf( fp, " %16s",  "RadOfCurv");

    fprintf( fp, " %12s", "Xgsm" ); // n. foot
    fprintf( fp, " %12s", "Ygsm" );
    fprintf( fp, " %12s", "Zgsm" );

    fprintf( fp, " %12s", "Xgeo" );
    fprintf( fp, " %12s", "Ygeo" );
    fprintf( fp, " %12s", "Zgeo" );

    fprintf( fp, " %12s", "GeodLat" );
    fprintf( fp, " %12s", "GeodLon" );
    fprintf( fp, " %12s", "GeodHeight" );

    fprintf( fp, " %12s", "Pfn_CD_MLAT" );
    fprintf( fp, " %12s", "Pfn_CD_MLON" );
    fprintf( fp, " %12s", "Pfn_CD_MLT" );

    fprintf( fp, " %12s", "Pfn_ED_MLAT" );
    fprintf( fp, " %12s", "Pfn_ED_MLON" );
    fprintf( fp, " %12s", "Pfn_ED_MLT" );

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

    fprintf( fp, " %12s", "Pfs_CD_MLAT" );
    fprintf( fp, " %12s", "Pfs_CD_MLON" );
    fprintf( fp, " %12s", "Pfs_CD_MLT" );

    fprintf( fp, " %12s", "Pfs_ED_MLAT" );
    fprintf( fp, " %12s", "Pfs_ED_MLON" );
    fprintf( fp, " %12s", "Pfs_ED_MLT" );

    fprintf( fp, " %12s", "Bfs_geo_x" ); // Bfs  geo
    fprintf( fp, " %12s", "Bfs_geo_y" );
    fprintf( fp, " %12s", "Bfs_geo_z" );
    fprintf( fp, " %12s", "Bfs_geo" );

    fprintf( fp, " %12s", "Bfs_x" ); // Bfs  gsm
    fprintf( fp, " %12s", "Bfs_y" );
    fprintf( fp, " %12s", "Bfs_z" );
    fprintf( fp, " %12s", "Bfs" );

    fprintf( fp, " %12s", "Alpha_LC_S" ); // Southern Loss Cone

    fprintf( fp, " %12s", "Xgsm" ); // Pmin gsm
    fprintf( fp, " %12s", "Ygsm" );
    fprintf( fp, " %12s", "Zgsm" );

    fprintf( fp, " %12s", "Bmin_x" ); // |B|min gsm
    fprintf( fp, " %12s", "Bmin_y" );
    fprintf( fp, " %12s", "Bmin_z" );
    fprintf( fp, " %12s", "Bmin" );

    fprintf( fp, " %12s", "Lsimple" );
    fprintf( fp, " %12s", "InvLat" );
    fprintf( fp, " %12s", "Lm_eq" );
    fprintf( fp, " %12s", "InvLatEq" );
    fprintf( fp, " %12s", "BoverBeq" );
    fprintf( fp, " %12s", "MlatBBeq" );

    fprintf( fp, " %12s", "M_used" );
    fprintf( fp, " %12s", "M_ref" );
    fprintf( fp, " %12s", "M_igrf" );
    fprintf(fp, "    ");
    for (i=0; i<m->nAlpha; i++) { sprintf( Str, "L*%d", i ); fprintf(fp, " %12s", Str ); }
    fprintf(fp, "    ");
    for (i=0; i<m->nAlpha; i++) { sprintf( Str, "DSType%d", i ); fprintf(fp, " %12s", Str ); }
    fprintf(fp, "    ");
    for (i=0; i<m->nAlpha; i++) { sprintf( Str, "L%d", i ); fprintf(fp, " %12s", Str ); }
    fprintf(fp, "    ");
    for (i=0; i<m->nAlpha; i++) { sprintf( Str, "Bm%d", i ); fprintf(fp, " %12s", Str ); }
    fprintf(fp, "    ");
    for (i=0; i<m->nAlpha; i++) { sprintf( Str, "I%d", i ); fprintf(fp, " %12s", Str ); }
    fprintf(fp, "    ");
    for (i=0; i<m->nAlpha; i++) { sprintf( Str, "K%d", i ); fprintf(fp, " %12s", Str ); }
    fprintf(fp, "    ");
    for (i=0; i<m->nAlpha; i++) { sprintf( Str, "Sb%d", i ); fprintf(fp, " %12s", Str ); }
    fprintf(fp, "    ");
    for (i=0; i<m->nAlpha; i++) { sprintf( Str, "Tb%d", i ); fprintf(fp, " %12s", Str ); }
    fprintf(fp, "    ");
    for (i=0; i<m->nAlpha; i++) { sprintf( Str, "Kappa%d", i ); fprintf(fp, " %12s", Str ); }
    fprintf(fp, "\n");
    // units/format
    fprintf( fp, "# %25s", "YYYY-MM-DDTHH:MM:SS.SSSSZ" );
    fprintf( fp, " %10s", "YYYYMMDD" );
    fprintf( fp, " %5s",  "DDD" );
    fprintf( fp, " %13s", "Hours" );
    fprintf( fp, " %16s", "Days" );
    fprintf( fp, " %15s", "Seconds" );
    fprintf( fp, " %13s", "Degrees" );
    fprintf( fp, " %11s", "       " );
    fprintf( fp, " %11s", "       " );

    fprintf( fp, " %13s", "Re" ); // Geocentric GEO
    fprintf( fp, " %13s", "Re" );
    fprintf( fp, " %13s", "Re" );

    fprintf( fp, " %13s", "Deg." ); // Geodetic GEO
    fprintf( fp, " %13s", "Deg." );
    fprintf( fp, " %13s", "km" );

    fprintf( fp, " %13s", "Re" );
    fprintf( fp, " %13s", "Re" );
    fprintf( fp, " %13s", "Re" );

    fprintf( fp, " %13s", "Re" );
    fprintf( fp, " %13s", "Re" );
    fprintf( fp, " %13s", "Re" );

    fprintf( fp, " %13s", "Re" );
    fprintf( fp, " %13s", "Re" );
    fprintf( fp, " %13s", "Re" );

    fprintf( fp, " %13s", "Re" );
    fprintf( fp, " %13s", "Re" );
    fprintf( fp, " %13s", "Re" );

    fprintf( fp, " %13s", "Deg." );
    fprintf( fp, " %13s", "Deg." );
    fprintf( fp, " %13s", "Hours" );
    fprintf( fp, " %13s", "Re" );

    fprintf( fp, " %13s", "Deg." );
    fprintf( fp, " %13s", "Deg." );
    fprintf( fp, " %13s", "Hours" );
    fprintf( fp, " %13s", "Re" );

    fprintf( fp, " %14s",  " " );  // Int Model
    fprintf( fp, " %14s",  " " );  // Ext Model
    fprintf( fp, " %7s",   " " );  // Kp
    fprintf( fp, " %8s",   "nT" ); // Dst

    fprintf( fp, " %12s", "nT" );   // Bsc
    fprintf( fp, " %12s", "nT" );
    fprintf( fp, " %12s", "nT" );
    fprintf( fp, " %12s", "nT" );

    fprintf( fp, " %29s",  "" );    // Field Line Type

    fprintf( fp, " %16s",  "Re" ); // S
    fprintf( fp, " %16s",  "Re" ); // S
    fprintf( fp, " %16s",  "Re" ); // S
    fprintf( fp, " %16s",  "Re" ); // S
    fprintf( fp, " %16s",  "Re" ); // S_total

    fprintf( fp, " %16s",  "nT^2/Re^2" ); // d2B_ds2
    fprintf( fp, " %16s",  "Re" ); // Sb0
    fprintf( fp, " %16s",  "Re" ); // Radius of Curvature

    fprintf( fp, " %12s", "Re" );   // GSM
    fprintf( fp, " %12s", "Re" );
    fprintf( fp, " %12s", "Re" );

    fprintf( fp, " %12s", "Re" );   // GEO
    fprintf( fp, " %12s", "Re" );
    fprintf( fp, " %12s", "Re" );

    fprintf( fp, " %12s", "Deg." ); // Geodetic
    fprintf( fp, " %12s", "Deg." );
    fprintf( fp, " %12s", "km" );

    fprintf( fp, " %12s", "Deg." ); // Pfn_CD_MLAT
    fprintf( fp, " %12s", "Deg." ); // Pfn_CD_MLON
    fprintf( fp, " %12s", "Hours" );// Pfn_CD_MLT

    fprintf( fp, " %12s", "Deg." ); // Pfn_ED_MLAT
    fprintf( fp, " %12s", "Deg." ); // Pfn_ED_MLON
    fprintf( fp, " %12s", "Hours" );// Pfn_ED_MLT

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

    fprintf( fp, " %12s", "Deg." ); // Pfs_CD_MLAT
    fprintf( fp, " %12s", "Deg." ); // Pfs_CD_MLON
    fprintf( fp, " %12s", "Hours" );// Pfs_CD_MLT

    fprintf( fp, " %12s", "Deg." ); // Pfs_ED_MLAT
    fprintf( fp, " %12s", "Deg." ); // Pfs_ED_MLON
    fprintf( fp, " %12s", "Hours" );// Pfs_ED_MLT

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

    fprintf( fp, " %12s", "Dimless" );   // Lsimple units
    fprintf( fp, " %12s", "Degrees" );   // InvLat units

    fprintf( fp, " %12s", "Dimless" );   // Lm_eq units
    fprintf( fp, " %12s", "Degrees" );   // InvLat_eq units

    fprintf( fp, " %12s", "Dimless" );   // BoverBeq
    fprintf( fp, " %12s", "Degrees" );   // MlatFromBoverBeq

    fprintf( fp, " %12s", "nT" );   // M_Used
    fprintf( fp, " %12s", "nT" );   // M_Ref
    fprintf( fp, " %12s", "nT" );   // M_IGRF

    fprintf(fp, "    ");
    for (i=0; i<m->nAlpha; i++) { sprintf( Str, "Dimless" ); fprintf(fp, " %12s", Str ); }
    fprintf(fp, "    ");
    for (i=0; i<m->nAlpha; i++) { sprintf( Str, "Dimless" ); fprintf(fp, " %12s", Str ); }
    fprintf(fp, "    ");
    for (i=0; i<m->nAlpha; i++) { sprintf( Str, "Dimless" ); fprintf(fp, " %12s", Str ); }
    fprintf(fp, "    ");
    for (i=0; i<m->nAlpha; i++) { sprintf( Str, "nT" ); fprintf(fp, " %12s", Str ); }
    fprintf(fp, "    ");
    for (i=0; i<m->nAlpha; i++) { sprintf( Str, "Re" ); fprintf(fp, " %12s", Str ); }
    fprintf(fp, "    ");
    for (i=0; i<m->nAlpha; i++) { sprintf( Str, "Re G^0.5" ); fprintf(fp, " %12s", Str ); }
    fprintf(fp, "    ");
    for (i=0; i<m->nAlpha; i++) { sprintf( Str, "Re" ); fprintf(fp, " %12s", Str ); }
    fprintf(fp, "    ");
    for (i=0; i<m->nAlpha; i++) { sprintf( Str, "s" ); fprintf(fp, " %12s", Str ); }
    fprintf(fp, "    ");
    for (i=0; i<m->nAlpha; i++) { sprintf( Str, "dimless" ); fprintf(fp, " %12s", Str ); }
    fprintf(fp, "\n");

    Lgm_free_ctrans(c);

}


void Lgm_WriteMagEphemData( FILE *fp, char *IntModel, char *ExtModel, double Kp, double Dst, Lgm_MagEphemInfo *m ){

    int             i;
    char            Str[128];
    double          GeodLat, GeodLong, GeodHeight, L;
    double          Bsc_mag, Bfn_mag, Bfs_mag, Bmin_mag, Alpha_Loss_Cone_n, Alpha_Loss_Cone_s;
    double          R, MLAT, MLON, MLT;
    double          Lm_eq, InvLat_eq, s, cl, BoverBeq, MagLatFromBoverBeq, Lsimple, InvLat, S_Bmin_to_sc;
    double          Ek, E, p2c2, Beta2, Beta, T, vel, p, rg, Kappa;
    Lgm_DateTime    DT_UTC;
    Lgm_CTrans      *c = Lgm_init_ctrans(0);
    Lgm_Vector      v, vv, Bsc, Bfn, Bfn_geo, Bfs, Bfs_geo, Bmin;

    Lgm_Set_Coord_Transforms( m->Date, m->UTC, c );
    Lgm_Make_UTC( m->Date, m->UTC, &DT_UTC, c );
    Lgm_DateTimeToString( Str, &DT_UTC, 0, 4 );
    

    fprintf( fp, "%25s",     Str );          // Date+Time in ISO 8601 format
    fprintf( fp, "   %10ld", c->UTC.Date );  // Date
    fprintf( fp, " %5d",     c->UTC.Doy );   // DOY
    fprintf( fp, " %13.8lf", c->UTC.Time );  // UTC
    fprintf( fp, " %16.8lf", c->UTC.JD );    // Julian Date
    fprintf( fp, " %15.3lf", Lgm_UTC_to_GpsSeconds( &c->UTC, c ) ); // GpsTime
    fprintf( fp, " %13.8lf", c->psi*DegPerRad ); // DipoleTiltAngle
    fprintf( fp, " %11d", m->InOut );        // InOut
    fprintf( fp, " %11d", m->OrbitNumber );        // InOut


    Lgm_Convert_Coords( &m->P, &v, GSM_TO_GEO, c );
    fprintf( fp, " %13.6lf", v.x );     // Xgeo
    fprintf( fp, " %13.6lf", v.y );     // Ygeo
    fprintf( fp, " %13.6lf", v.z );     // Zgeo

    Lgm_WGS84_to_GEOD( &v, &GeodLat, &GeodLong, &GeodHeight );
    fprintf( fp, " %13.6lf", GeodLat );                // Geod Lat   of SC
    fprintf( fp, " %13.6lf", GeodLong );               // Geod Long
    fprintf( fp, " %13.4lf", GeodHeight );             // Geod Height

    fprintf( fp, " %13.6lf", m->P.x );  // Xgsm
    fprintf( fp, " %13.6lf", m->P.y );  // Ygsm
    fprintf( fp, " %13.6lf", m->P.z );  // Zgsm

    Lgm_Convert_Coords( &m->P, &v, GSM_TO_SM, c );
    fprintf( fp, " %13.6lf", v.x );        // Xsm
    fprintf( fp, " %13.6lf", v.y );        // Ysm
    fprintf( fp, " %13.6lf", v.z );        // Zsm

    Lgm_Convert_Coords( &m->P, &v, GSM_TO_GEI2000, c );
    fprintf( fp, " %13.6lf", v.x );        // Xgei
    fprintf( fp, " %13.6lf", v.y );        // Ygei
    fprintf( fp, " %13.6lf", v.z );        // Zgei

    Lgm_Convert_Coords( &m->P, &v, GSM_TO_GSE, c );
    fprintf( fp, " %13.6lf", v.x );        // Xgse
    fprintf( fp, " %13.6lf", v.y );        // Ygse
    fprintf( fp, " %13.6lf", v.z );        // Zgse


    Lgm_Convert_Coords( &m->P, &v, GSM_TO_CDMAG, c );   // Convert to cartesian CDMAG coords.
    Lgm_CDMAG_to_R_MLAT_MLON_MLT( &v, &R, &MLAT, &MLON, &MLT, c );
    fprintf( fp, " %13.6lf", MLAT );                    // CD MLAT
    fprintf( fp, " %13.6lf", MLON );                    // CD MLON
    fprintf( fp, " %13.6lf", MLT );                     // CD MLT
    fprintf( fp, " %13.6lf", R );                       // CD R (since CDMAG is geocentric, this should be same as |v| )

    Lgm_Convert_Coords( &m->P, &v, GSM_TO_EDMAG, c );   // Convert to cartesian EDMAG coords.
    Lgm_EDMAG_to_R_MLAT_MLON_MLT( &v, &R, &MLAT, &MLON, &MLT, c );
    fprintf( fp, " %13.6lf", MLAT );                    // ED MLAT
    fprintf( fp, " %13.6lf", MLON );                    // ED MLON
    fprintf( fp, " %13.6lf", MLT );                     // ED MLT
    fprintf( fp, " %13.6lf", R );                       // ED R (since EDMAG is NOT geocentric, this should NOT be same as |v| (in general))



    if ( !strcmp( ExtModel, "IGRF" ) || !strcmp( ExtModel, "CDIP" ) || !strcmp( ExtModel, "EDIP" ) ) {
        // If our "external model is just a dipole or igrf, then "internal doesnt really mean anything...)
        fprintf( fp, " %14s", "N/A" );       // Int model is Not applicable
    } else {
        fprintf( fp, " %14s", IntModel );    // Int Model
    }
    fprintf( fp, " %14s", ExtModel );        // Ext Model
    fprintf( fp, " %7.1f",  Kp );            // Kp
    fprintf( fp, " %8.3lf",  Dst );             // Dst
    
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


    fprintf( fp, " %16g", (m->Snorth > 0.0) ? m->Snorth : LGM_FILL_VALUE ); // S_sc_to_pfn
    fprintf( fp, " %16g", (m->Ssouth > 0.0) ? m->Ssouth : LGM_FILL_VALUE ); // S_sc_to_pfs
    fprintf( fp, " %16g", (m->Smin > 0.0) ? m->Smin : LGM_FILL_VALUE ); // S_pfs_to_Bmin
    S_Bmin_to_sc = ((m->Ssouth>0.0)&&(m->Smin > 0.0)) ? m->Ssouth-m->Smin : LGM_FILL_VALUE;
    fprintf( fp, " %16g", S_Bmin_to_sc ); // S_Bmin_to_sc
    fprintf( fp, " %16g", ((m->Snorth > 0.0)&&(m->Ssouth > 0.0)) ? m->Snorth + m->Ssouth : LGM_FILL_VALUE ); // S_total


    fprintf( fp, " %16g", m->d2B_ds2 ); // d2B_ds2
    fprintf( fp, " %16g", m->Sb0 );     // Sb0
    fprintf( fp, " %16g", m->RofC );    // Radius of curvature at Bmin point


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


        Lgm_Convert_Coords( &m->Ellipsoid_Footprint_Pn, &vv, GSM_TO_CDMAG, c );   // Convert to cartesian CDMAG coords.
        Lgm_CDMAG_to_R_MLAT_MLON_MLT( &vv, &R, &MLAT, &MLON, &MLT, c );
        fprintf( fp, " %12g", MLAT );                    // CD MLAT
        fprintf( fp, " %12g", MLON );                    // CD MLON
        fprintf( fp, " %12g", MLT );                     // CD MLT

        Lgm_Convert_Coords( &m->Ellipsoid_Footprint_Pn, &vv, GSM_TO_EDMAG, c );   // Convert to cartesian EDMAG coords.
        Lgm_EDMAG_to_R_MLAT_MLON_MLT( &vv, &R, &MLAT, &MLON, &MLT, c );
        fprintf( fp, " %12g", MLAT );                    // ED MLAT
        fprintf( fp, " %12g", MLON );                    // ED MLON
        fprintf( fp, " %12g", MLT );                     // ED MLT


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

        Lgm_Convert_Coords( &m->Ellipsoid_Footprint_Ps, &vv, GSM_TO_CDMAG, c );   // Convert to cartesian CDMAG coords.
        Lgm_CDMAG_to_R_MLAT_MLON_MLT( &vv, &R, &MLAT, &MLON, &MLT, c );
        fprintf( fp, " %12g", MLAT );                    // CD MLAT
        fprintf( fp, " %12g", MLON );                    // CD MLON
        fprintf( fp, " %12g", MLT );                     // CD MLT

        Lgm_Convert_Coords( &m->Ellipsoid_Footprint_Ps, &vv, GSM_TO_EDMAG, c );   // Convert to cartesian EDMAG coords.
        Lgm_EDMAG_to_R_MLAT_MLON_MLT( &vv, &R, &MLAT, &MLON, &MLT, c );
        fprintf( fp, " %12g", MLAT );                    // ED MLAT
        fprintf( fp, " %12g", MLON );                    // ED MLON
        fprintf( fp, " %12g", MLT );                     // ED MLT

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

    Lsimple = (Bmin_mag > 0.0) ? Lgm_Magnitude( &m->Pmin ) : LGM_FILL_VALUE;
    fprintf( fp, " %12g", Lsimple );      // Lsimple (equatorial distance to Bmin point)
    if (Lsimple > 0.0) {
        InvLat = DegPerRad*acos(sqrt(1.0/Lsimple));
        //if (S_Bmin_to_sc<0.0) InvLat *= -1.0;
    } else {
        InvLat = LGM_FILL_VALUE;
    }
    fprintf( fp, " %12g", InvLat );   // InvLat


    
    Lm_eq = (Bmin_mag > 0.0) ? LFromIBmM_McIlwain( 0.0, Bmin_mag, m->Mcurr ) : LGM_FILL_VALUE;
    fprintf( fp, " %12g", Lm_eq );      // Lm_eq
    if (Lm_eq > 0.0) {
        InvLat_eq = DegPerRad*acos(sqrt(1.0/Lm_eq));
        //if (S_Bmin_to_sc<0.0) InvLat_eq *= -1.0;
    } else {
        InvLat_eq = LGM_FILL_VALUE;
    }
    fprintf( fp, " %12g", InvLat_eq );   // InvLat_eq

    
    BoverBeq = ( Bmin_mag > 0.0) ? Bsc_mag / Bmin_mag : LGM_FILL_VALUE;
    fprintf( fp, " %12g", BoverBeq );   // BoverBeq

    
    
    /* 
     * In a dipole BoverBeq = sqrt(4-3cos^2(lambda))/cos^6(lambda) Define
     * MagLatFromBoverBeq as the dipole latitude that has the same BoverBeq as
     * the realistic model field.  You have to solve the same dipole mirror lat
     * equation that is done in Lgm_CdipMirrorLat where s = sqrt( 1.0/BoverBeq )
     */
    if ( BoverBeq > 0.0 ) {
        s = sqrt( 1.0/BoverBeq ); 
        cl = Lgm_CdipMirrorLat( s );
        if ( fabs(cl) <= 1.0 ){
            MagLatFromBoverBeq = DegPerRad*acos( cl );
            if (S_Bmin_to_sc<0.0) MagLatFromBoverBeq *= -1.0;
        } else {
            MagLatFromBoverBeq = LGM_FILL_VALUE;
        }
    } else {
        MagLatFromBoverBeq = LGM_FILL_VALUE;
    }
    fprintf( fp, " %12g", MagLatFromBoverBeq );   // MagLatFromBoverBeq


    fprintf( fp, " %12g", m->Mused );   // M_Used
    fprintf( fp, " %12g", m->Mref );    // M_Ref
    fprintf( fp, " %12g", m->Mcurr );   // M_IGRF

    // L*'s
    fprintf(fp, "    ");
    for (i=0; i<m->nAlpha; i++) { fprintf(fp, " %12g", m->Lstar[i] ); }

    // Drift Shell Types
    fprintf(fp, "    ");
    for (i=0; i<m->nAlpha; i++) { fprintf(fp, " %12d", m->DriftOrbitType[i] ); }

    // McIlwain L (computed from I, Bm, M)
    fprintf(fp, "    ");
    for (i=0; i<m->nAlpha; i++) { 
        L = ( m->I[i] >= 0.0 ) ? LFromIBmM_McIlwain(m->I[i], m->Bm[i], m->Mused ) : LGM_FILL_VALUE;
        fprintf(fp, " %12g", L);
    }

    // Bms's
    fprintf(fp, "    ");
    for (i=0; i<m->nAlpha; i++) { fprintf(fp, " %12g", m->Bm[i] ); }

    // I's
    fprintf(fp, "    ");
    for (i=0; i<m->nAlpha; i++) { fprintf(fp, " %12g", m->I[i] ); }

    // K's
    fprintf(fp, "    ");
    for (i=0; i<m->nAlpha; i++) { 
        if ( (m->I[i] > 0.0) && (m->Bm[i] > 0.0) ) {
            fprintf(fp, " %12g", m->I[i]*sqrt(m->Bm[i]*1e-5) ); 
        } else {
            fprintf(fp, " %12g", LGM_FILL_VALUE ); 
        }
    }

    // Sb's
    fprintf(fp, "    ");
    for (i=0; i<m->nAlpha; i++) { fprintf(fp, " %12g", m->Sb[i] ); }

    // Tb's
    fprintf(fp, "    ");
    for (i=0; i<m->nAlpha; i++) { 


        Ek    = 1.0; // MeV
        E     = Ek + LGM_Ee0; // total energy, MeV
        p2c2  = Ek*(Ek+2.0*LGM_Ee0); // p^2c^2,  MeV^2
        Beta2 = p2c2/(E*E); // beta^2 = v^2/c^2   (dimensionless)
        Beta  = sqrt( Beta2 );
        vel   = Beta*LGM_c;  // velocity in m/s
        vel  /= (Re*1000.0); // Re/s
        T     = ( m->Sb[i] > 0.0 ) ? 2.0*m->Sb[i]/vel : LGM_FILL_VALUE;

        fprintf(fp, " %12g", T ); 
    }

    // Kappa's
    fprintf(fp, "    ");
    for (i=0; i<m->nAlpha; i++) { 


        Ek    = 1.0; // MeV
        E     = Ek + LGM_Ee0; // total energy, MeV
        p2c2  = Ek*(Ek+2.0*LGM_Ee0); // p^2c^2,  MeV^2
        p     = sqrt(p2c2)*1.60217646e-13/LGM_c;  // mks
        rg    = sin(m->Alpha[i]*RadPerDeg)*p/(LGM_e*Bmin_mag*1e-9); // m. Bmin_mag calced above

        Kappa = sqrt( m->RofC*Re*1e3/rg );

        fprintf(fp, " %12g", Kappa ); 
    }


    fprintf(fp, "\n");

    
    Lgm_free_ctrans( c );
}
