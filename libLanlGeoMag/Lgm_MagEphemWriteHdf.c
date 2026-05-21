#include "Lgm/Lgm_CTrans.h"
#include "Lgm/Lgm_HDF5.h"
#include "Lgm/Lgm_MagEphemInfo.h"
//const char *sMonth[] = { "", "January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December" };


void Lgm_WriteMagEphemHeaderHdf( hid_t file, char *CodeVersion, char *ExtModel, int SpiceBody,  char *Spacecraft, int IdNumber, char *IntDesig, char *CmdLine, int nAscend, Lgm_DateTime *Ascend_UTC, Lgm_Vector *Ascend_U, int nPerigee, Lgm_DateTime *Perigee_UTC, Lgm_Vector *Perigee_U, int nApogee, Lgm_DateTime *Apogee_UTC, Lgm_Vector *Apogee_U, Lgm_MagEphemInfo *m, Lgm_MagEphemData *med ){

    int             Rank;
    hsize_t         Dims[4], MaxDims[4], ChunkDims[4];
    hid_t           cparms;
    herr_t          status;
    hid_t           space;
    hid_t           atype;
    hid_t           DataSet;
    char            TmpStr[1024];



    // Create and Write Alpha Dataset
    DataSet = CreateSimpleRank1DataSet( file, "Alpha", m->nAlpha, H5T_NATIVE_DOUBLE, &space );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "Degrees" );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );
    LGM_HDF5_WRITE_SIMPLE_RANK1_DATASET( file, "Alpha", med->H5_nAlpha, H5T_NATIVE_DOUBLE, &med->H5_Alpha[0] );

    // Create and Write AscendTimes Dataset
    atype = CreateStrType( 32 );
    DataSet = CreateSimpleRank1DataSet( file, "AscendTimes", nAscend, atype, &space );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "UTC" );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "The times of Ascending Node Crossings  in ISO 8601 compliant format." );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );
    LGM_HDF5_WRITE_SIMPLE_RANK1_STRING_DATASET( file, "AscendTimes", med->H5_nAscend, atype, med->H5_Ascend_IsoTimes );
    status  = H5Tclose( atype );


    // Create and Write AscendPosGeod Dataset
    DataSet = CreateSimpleRank2DataSet( file, "AscendPosGeod", med->H5_nAscend, 3, H5T_NATIVE_DOUBLE, &space );
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "AscendTimes" );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Geodetic position of ascending node croissing (lat/lon/rad)." );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "Deg./Deg./km" );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );
    //printf("med->H5_Ascend_Geod[0][0] = %g\n", med->H5_Ascend_Geod[0][0]);
    LGM_HDF5_WRITE_SIMPLE_RANK2_DATASET( file, "AscendPosGeod", med->H5_nAscend, 3, H5T_NATIVE_DOUBLE, &med->H5_Ascend_Geod[0][0] );

    // Create and Write PerigeeTimes Dataset
    atype = CreateStrType( 32 );
    DataSet = CreateSimpleRank1DataSet( file, "PerigeeTimes", nPerigee, atype, &space );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "UTC" );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "The times of Perigee in ISO 8601 compliant format. Defined as smallest geocentric distance." );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );
    LGM_HDF5_WRITE_SIMPLE_RANK1_STRING_DATASET( file, "PerigeeTimes", med->H5_nPerigee, atype, med->H5_Perigee_IsoTimes );
    status  = H5Tclose( atype );

    // Create and Write PerigeePosGeod Dataset
    DataSet = CreateSimpleRank2DataSet( file, "PerigeePosGeod", med->H5_nPerigee, 3, H5T_NATIVE_DOUBLE, &space );
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "PerigeeTimes" );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Geodetic position of perigee (lat/lon/rad)." );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "Deg./Deg./km" );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );
    LGM_HDF5_WRITE_SIMPLE_RANK2_DATASET( file, "PerigeePosGeod", med->H5_nPerigee, 3, H5T_NATIVE_DOUBLE, &med->H5_Perigee_Geod[0][0] );


    // Create and Write ApogeeTimes Dataset
    atype = CreateStrType( 32 );
    DataSet = CreateSimpleRank1DataSet( file, "ApogeeTimes", nApogee, atype, &space );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "The times of Apogee in ISO 8601 compliant format. Defined as largest geocentric distance." );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "UTC" );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );
    LGM_HDF5_WRITE_SIMPLE_RANK1_STRING_DATASET( file, "ApogeeTimes", med->H5_nApogee, atype, med->H5_Apogee_IsoTimes );
    status  = H5Tclose( atype );



    // Create and Write ApogeePosGeod Dataset
    DataSet = CreateSimpleRank2DataSet( file, "ApogeePosGeod", med->H5_nApogee, 3, H5T_NATIVE_DOUBLE, &space );
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "ApogeeTimes" );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Geodetic position of apogee (lat/lon/rad)." );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "Deg./Deg./km" );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );
    LGM_HDF5_WRITE_SIMPLE_RANK2_DATASET( file, "ApogeePosGeod", med->H5_nApogee, 3, H5T_NATIVE_DOUBLE, &med->H5_Apogee_Geod[0][0] );


    // Create IsoTime Dataset
    atype = CreateStrType( 32 );
    DataSet = CreateExtendableRank1DataSet( file, "IsoTime", atype, &space );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "The date and time in ISO 8601 compliant format." );
    Lgm_WriteStringAttr( DataSet, "UNITS",       "UTC" );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",    "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",     "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",    "data" );
    status  = H5Tclose( atype );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );

    // Create Date Dataset
    DataSet = CreateExtendableRank1DataSet( file, "Date", H5T_NATIVE_LONG, &space );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "The date. In YYYMMDD format." );
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "YYYYMMDD" );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );


    // Create Doy Dataset
    DataSet = CreateExtendableRank1DataSet( file, "Doy", H5T_NATIVE_INT, &space );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Ordinal Day of Year." );
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "DDD" );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );


    // Create UTC (hours) Dataset
    DataSet = CreateExtendableRank1DataSet( file, "UTC", H5T_NATIVE_DOUBLE, &space );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Universal Time (Coordinated). In decimal hours." );
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "Hours" );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );


    // Create JD Dataset
    DataSet = CreateExtendableRank1DataSet( file, "JulianDate", H5T_NATIVE_DOUBLE, &space );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Julian Date. In decimal days." );
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "Days" );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );


    // Create GpsTime Dataset
    DataSet = CreateExtendableRank1DataSet( file, "GpsTime", H5T_NATIVE_DOUBLE, &space );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Number of SI seconds since 0h Jan 6, 1980 UTC." );
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "Seconds" );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );


    // Create DipoleTiltAngle Dataset
    DataSet = CreateExtendableRank1DataSet( file, "DipoleTiltAngle", H5T_NATIVE_DOUBLE, &space );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Angle between Zgsm and Zsm (i.e. between Zgsm and dipole axis direction)." );
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "Degrees" );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );


    // Create InOut Dataset
    DataSet = CreateExtendableRank1DataSet( file, "InOut", H5T_NATIVE_INT, &space );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Flag indicating whether we are inbound (-1) or outbound (+1)" );
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "dimless" );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );

    // Create OrbitNumber Dataset
    DataSet = CreateExtendableRank1DataSet( file, "OrbitNumber", H5T_NATIVE_INT, &space );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Orbit Number" );
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "dimless" );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );


    // Create SunAngle Dataset
    DataSet = CreateExtendableRank1DataSet( file, "SunAngle", H5T_NATIVE_DOUBLE, &space );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "SunAngle" );
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "degrees" );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );


    // Create Eclipse Flag Dataset
    DataSet = CreateExtendableRank1DataSet( file, "Eclipse", H5T_NATIVE_INT, &space );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Eclipse Flag" );
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "dimless" );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );


    // Create Rgeo Dataset
    DataSet = CreateExtendableRank2DataSet( file, "Rgeo", 3, H5T_NATIVE_DOUBLE, &space );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Geocentric Geographic position vector of S/C." );
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "Re" );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );

    // Create Rgeod_LatLon Dataset
    DataSet = CreateExtendableRank2DataSet( file, "Rgeod_LatLon", 2, H5T_NATIVE_DOUBLE, &space );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Geodetic Geographic Latitude and Longitude of S/C." );
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "Deg." );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );

    // Create Rgeod_Height Dataset
    DataSet = CreateExtendableRank1DataSet( file, "Rgeod_Height", H5T_NATIVE_DOUBLE, &space );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Geodetic Geographic Height (Above WGS84 Ellipsoid) of S/C." );
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "km" );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );


    // Create Rgsm Dataset
    DataSet = CreateExtendableRank2DataSet( file, "Rgsm", 3, H5T_NATIVE_DOUBLE, &space );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Geocentric Solar Magnetospheric position vector of S/C." );
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "Re" );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );

    // Create Rsm Dataset
    DataSet = CreateExtendableRank2DataSet( file, "Rsm", 3, H5T_NATIVE_DOUBLE, &space );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Geocentric Solar Magnetic position vector of S/C." );
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "Re" );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );

    // Create Rgei Dataset
    DataSet = CreateExtendableRank2DataSet( file, "Rgei", 3, H5T_NATIVE_DOUBLE, &space );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Geocentric Equatorial Inertial position vector of S/C." );
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "Re" );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );

    // Create Rgse Dataset
    DataSet = CreateExtendableRank2DataSet( file, "Rgse", 3, H5T_NATIVE_DOUBLE, &space );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Geocentric Solar Ecliptic position vector of S/C." );
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "Re" );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );



    // Create CDMAG_MLAT Dataset
    DataSet = CreateExtendableRank1DataSet( file, "CDMAG_MLAT", H5T_NATIVE_DOUBLE, &space );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Magnetic Latitude of S/C in Centerted Dipole Coordinates." );
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "Deg." );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );

    // Create CDMAG_MLON Dataset
    DataSet = CreateExtendableRank1DataSet( file, "CDMAG_MLON", H5T_NATIVE_DOUBLE, &space );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Magnetic Longitude of S/C in Centerted Dipole Coordinates." );
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "Deg." );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );

    // Create CDMAG_MLT Dataset
    DataSet = CreateExtendableRank1DataSet( file, "CDMAG_MLT", H5T_NATIVE_DOUBLE, &space );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Magnetic Local Time of S/C in Centerted Dipole Coordinates." );
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "Deg." );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );

    // Create CDMAG_R Dataset
    DataSet = CreateExtendableRank1DataSet( file, "CDMAG_R", H5T_NATIVE_DOUBLE, &space );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Radial distance of S/C from center of CDMAG coordinate system." );
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "Deg." );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );


    // Create EDMAG_MLAT Dataset
    DataSet = CreateExtendableRank1DataSet( file, "EDMAG_MLAT", H5T_NATIVE_DOUBLE, &space );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Magnetic Latitude of S/C in Eccentric Dipole Coordinates." );
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "Deg." );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );

    // Create EDMAG_MLON Dataset
    DataSet = CreateExtendableRank1DataSet( file, "EDMAG_MLON", H5T_NATIVE_DOUBLE, &space );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Magnetic Longitude of S/C in Eccentric Dipole Coordinates." );
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "Deg." );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );

    // Create EDMAG_MLT Dataset
    DataSet = CreateExtendableRank1DataSet( file, "EDMAG_MLT", H5T_NATIVE_DOUBLE, &space );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Magnetic Local Time of S/C in Eccentric Dipole Coordinates." );
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "Deg." );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );

    // Create EDMAG_R Dataset
    DataSet = CreateExtendableRank1DataSet( file, "EDMAG_R", H5T_NATIVE_DOUBLE, &space );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Radial distance of S/C from center of EDMAG coordinate system." );
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "Deg." );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );




    // Create IntModel Dataset
    atype = CreateStrType( 32 );
    DataSet = CreateExtendableRank1DataSet( file, "IntModel", atype, &space );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Internal magnetic field model." );
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );
    status  = H5Tclose( atype );


    // Create ExtModel Dataset
    atype = CreateStrType( 32 );
    DataSet = CreateExtendableRank1DataSet( file, "ExtModel", atype, &space );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "External magnetic field model." );
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );
    status  = H5Tclose( atype );



    // Create Kp Dataset
    DataSet = CreateExtendableRank1DataSet( file, "Kp", H5T_NATIVE_DOUBLE, &space );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Kp index value." );
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "Dimensionless" );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );

    // Create Dst Dataset
    DataSet = CreateExtendableRank1DataSet( file, "Dst", H5T_NATIVE_DOUBLE, &space );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Dst index value." );
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "nT" );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );


    // Create Bsc_gsm Dataset
    DataSet = CreateExtendableRank2DataSet( file, "Bsc_gsm", 4, H5T_NATIVE_DOUBLE, &space );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Magnetic field vector at S/C (in GSM coords)." );
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "nT" );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );




    // Create FieldLineType Dataset
    atype = CreateStrType( 32 );
    DataSet = CreateExtendableRank1DataSet( file, "FieldLineType", atype, &space );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Description of the type of field line the S/C is on., Can be one of 4 types: LGM_CLOSED      - FL hits Earth at both ends. LGM_OPEN_N_LOBE - FL is an OPEN field line rooted in the Northern polar cap. LGM_OPEN_S_LOBE - FL is an OPEN field line rooted in the Southern polar cap. LGM_OPEN_IMF    - FL does not hit Earth at eitrher end." );
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );
    status  = H5Tclose( atype );


    // Create S_sc_to_pfn Dataset
    DataSet = CreateExtendableRank1DataSet( file, "S_sc_to_pfn", H5T_NATIVE_DOUBLE, &space );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Distance between S/C and Northern Footpoint along field line." );
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "Re" );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );

    // Create S_sc_to_pfs Dataset
    DataSet = CreateExtendableRank1DataSet( file, "S_sc_to_pfs", H5T_NATIVE_DOUBLE, &space );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Distance between S/C and Southern Footpoint along field line." );
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "Re" );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );

    // Create S_pfs_to_Bmin Dataset
    DataSet = CreateExtendableRank1DataSet( file, "S_pfs_to_Bmin", H5T_NATIVE_DOUBLE, &space );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Distance between Southern Footpoint and Bmin point along field line.");
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "Re" );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );

    // Create S_Bmin_to_sc Dataset
    DataSet = CreateExtendableRank1DataSet( file, "S_Bmin_to_sc", H5T_NATIVE_DOUBLE, &space );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Distance between Bmin point and S/C along field line (positive if north of Bmin).");
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "Re" );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );

    // Create S_total Dataset
    DataSet = CreateExtendableRank1DataSet( file, "S_total", H5T_NATIVE_DOUBLE, &space );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Total Field Line length (along field line).");
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "Re" );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );

    // Create d2B_ds2 Dataset
    DataSet = CreateExtendableRank1DataSet( file, "d2B_ds2", H5T_NATIVE_DOUBLE, &space );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Second derivative of |B| with respect to s (dist along FL) at minimum |B| point.");
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "nT^2/Re^2" );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );

    // Create Sb0 Dataset
    DataSet = CreateExtendableRank1DataSet( file, "Sb0", H5T_NATIVE_DOUBLE, &space );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Value of the 'Sb Integral' for equatorially mirroring particles (not generally zero).");
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "Re" );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );

    // Create RadiusOfCurv Dataset
    DataSet = CreateExtendableRank1DataSet( file, "RadiusOfCurv", H5T_NATIVE_DOUBLE, &space );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Field line radius of curvature at minimum |B| point.");
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "Re" );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );





    // Create Pfn_geo Dataset
    DataSet = CreateExtendableRank2DataSet( file, "Pfn_geo", 3, H5T_NATIVE_DOUBLE, &space );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Location of Northern Footpoint (in GEO coords).");
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "Re" );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );

    // Create Pfn_gsm Dataset
    DataSet = CreateExtendableRank2DataSet( file, "Pfn_gsm", 3, H5T_NATIVE_DOUBLE, &space );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Location of Northern Footpoint (in GSM coords).");
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "Re" );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );



    // Create Pfn_geod_LatLon Dataset
    DataSet = CreateExtendableRank2DataSet( file, "Pfn_geod_LatLon", 2, H5T_NATIVE_DOUBLE, &space );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Geodetic Latitude and Longitude of Northern Footpoint.");
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "Deg." );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );

    // Create Pfn_geod_Height Dataset
    DataSet = CreateExtendableRank1DataSet( file, "Pfn_geod_Height", H5T_NATIVE_DOUBLE, &space );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Geodetic Height of Northern Footpoint.");
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "km" );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );


    // Create Pfn_CD_MLAT Dataset
    DataSet = CreateExtendableRank1DataSet( file, "Pfn_CD_MLAT", H5T_NATIVE_DOUBLE, &space );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Magnetic Latitude of Northern Footpoint in Centerted Dipole Coordinates.");
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "Deg." );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );

    // Create Pfn_CD_MLON Dataset
    DataSet = CreateExtendableRank1DataSet( file, "Pfn_CD_MLON", H5T_NATIVE_DOUBLE, &space );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Magnetic Longitude of Northern Footpoint in Centerted Dipole Coordinates.");
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "Deg." );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );

    // Create Pfn_CD_MLT Dataset
    DataSet = CreateExtendableRank1DataSet( file, "Pfn_CD_MLT", H5T_NATIVE_DOUBLE, &space );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Magnetic Local Time of Northern Footpoint in Centerted Dipole Coordinates.");
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "Hours" );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );

    // Create Pfn_ED_MLAT Dataset
    DataSet = CreateExtendableRank1DataSet( file, "Pfn_ED_MLAT", H5T_NATIVE_DOUBLE, &space );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Magnetic Latitude of Northern Footpoint in Eccentric Dipole Coordinates.");
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "Deg." );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );

    // Create Pfn_ED_MLON Dataset
    DataSet = CreateExtendableRank1DataSet( file, "Pfn_ED_MLON", H5T_NATIVE_DOUBLE, &space );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Magnetic Longitude of Northern Footpoint in Eccentric Dipole Coordinates.");
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "Deg." );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );

    // Create Pfn_ED_MLT Dataset
    DataSet = CreateExtendableRank1DataSet( file, "Pfn_ED_MLT", H5T_NATIVE_DOUBLE, &space );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Magnetic Local Time of Northern Footpoint in Eccentric Dipole Coordinates.");
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "Hours" );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );

    // Create Bfn_geo Dataset
    DataSet = CreateExtendableRank2DataSet( file, "Bfn_geo", 4, H5T_NATIVE_DOUBLE, &space );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Magnetic field vector at Northern Footpoint (in GEO coords).");
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "nT" );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );

    // Create Bfn_gsm Dataset
    DataSet = CreateExtendableRank2DataSet( file, "Bfn_gsm", 4, H5T_NATIVE_DOUBLE, &space );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Magnetic field vector at Northern Footpoint (in GSM coords).");
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "nT" );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );

    // Create Loss_Cone_Alpha_n Dataset
    DataSet = CreateExtendableRank1DataSet( file, "Loss_Cone_Alpha_n", H5T_NATIVE_DOUBLE, &space );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Value of Northern Loss Cone angle. asin( sqrt(Bsc/Bfn) ).");
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "Deg." );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );







    // Create Pfs_geo Dataset
    DataSet = CreateExtendableRank2DataSet( file, "Pfs_geo", 3, H5T_NATIVE_DOUBLE, &space );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Location of Southern Footpoint (in GEO coords).");
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "Re" );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );

    // Create Pfs_gsm Dataset
    DataSet = CreateExtendableRank2DataSet( file, "Pfs_gsm", 3, H5T_NATIVE_DOUBLE, &space );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Location of Southern Footpoint (in GSM coords).");
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "Re" );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );



    // Create Pfs_geod_LatLon Dataset
    DataSet = CreateExtendableRank2DataSet( file, "Pfs_geod_LatLon", 2, H5T_NATIVE_DOUBLE, &space );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Geodetic Latitude and Longitude of Southern Footpoint.");
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "Deg." );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );

    // Create Pfs_geod_Height Dataset
    DataSet = CreateExtendableRank1DataSet( file, "Pfs_geod_Height", H5T_NATIVE_DOUBLE, &space );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Geodetic Height of Southern Footpoint.");
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "km" );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );


    // Create Pfs_CD_MLAT Dataset
    DataSet = CreateExtendableRank1DataSet( file, "Pfs_CD_MLAT", H5T_NATIVE_DOUBLE, &space );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Magnetic Latitude of Southern Footpoint in Centerted Dipole Coordinates.");
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "Deg." );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );

    // Create Pfs_CD_MLON Dataset
    DataSet = CreateExtendableRank1DataSet( file, "Pfs_CD_MLON", H5T_NATIVE_DOUBLE, &space );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Magnetic Longitude of Southern Footpoint in Centerted Dipole Coordinates.");
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "Deg." );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );

    // Create Pfs_CD_MLT Dataset
    DataSet = CreateExtendableRank1DataSet( file, "Pfs_CD_MLT", H5T_NATIVE_DOUBLE, &space );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Magnetic Local Time of Southern Footpoint in Centerted Dipole Coordinates.");
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "Hours" );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );

    // Create Pfs_ED_MLAT Dataset
    DataSet = CreateExtendableRank1DataSet( file, "Pfs_ED_MLAT", H5T_NATIVE_DOUBLE, &space );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Magnetic Latitude of Southern Footpoint in Eccentric Dipole Coordinates.");
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "Deg." );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );

    // Create Pfs_ED_MLON Dataset
    DataSet = CreateExtendableRank1DataSet( file, "Pfs_ED_MLON", H5T_NATIVE_DOUBLE, &space );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Magnetic Longitude of Southern Footpoint in Eccentric Dipole Coordinates.");
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "Deg." );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );

    // Create Pfs_ED_MLT Dataset
    DataSet = CreateExtendableRank1DataSet( file, "Pfs_ED_MLT", H5T_NATIVE_DOUBLE, &space );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Magnetic Local Time of Southern Footpoint in Eccentric Dipole Coordinates.");
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "Hours" );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );

    // Create Bfs_geo Dataset
    DataSet = CreateExtendableRank2DataSet( file, "Bfs_geo", 4, H5T_NATIVE_DOUBLE, &space );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Magnetic field vector at Southern Footpoint (in GEO coords).");
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "nT" );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );

    // Create Bfs_gsm Dataset
    DataSet = CreateExtendableRank2DataSet( file, "Bfs_gsm", 4, H5T_NATIVE_DOUBLE, &space );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Magnetic field vector at Southern Footpoint (in GSM coords).");
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "nT" );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );

    // Create Loss_Cone_Alpha_s Dataset
    DataSet = CreateExtendableRank1DataSet( file, "Loss_Cone_Alpha_s", H5T_NATIVE_DOUBLE, &space );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Value of Southern Loss Cone angle. asin( sqrt(Bsc/Bfs) ).");
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "Deg." );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );


    // Create Pmin_gsm Dataset
    DataSet = CreateExtendableRank2DataSet( file, "Pmin_gsm", 3, H5T_NATIVE_DOUBLE, &space );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Location of minimum-|B| point (in GSM coords).");
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "Re" );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );

    // Create Bmin_gsm Dataset
    DataSet = CreateExtendableRank2DataSet( file, "Bmin_gsm", 4, H5T_NATIVE_DOUBLE, &space );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "B-field at minimum-|B| point (in GSM coords).");
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "Re" );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );



    // Create Lsimple Dataset
    DataSet = CreateExtendableRank1DataSet( file, "Lsimple", H5T_NATIVE_DOUBLE, &space );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Geocentric distance to Bmin point for FL threading vehicle (i.e. |Pmin|).");
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "Dimless" );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );

    // Create InvLat Dataset
    DataSet = CreateExtendableRank1DataSet( file, "InvLat", H5T_NATIVE_DOUBLE, &space );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Invariant latitude of vehicle computed from Lambda=acos(sqrt(1/Lsimple)).");
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "Deg." );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );


    // Create Lm_eq Dataset
    DataSet = CreateExtendableRank1DataSet( file, "Lm_eq", H5T_NATIVE_DOUBLE, &space );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "McIlwain L of an eq. mirroring particle on same FL as vehicle (computed from L=Lm_eq, I=0, and Bm=|Bmin_gsm|, M=M_igrf).");
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "Dimless" );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );

    // Create InvLat_eq Dataset
    DataSet = CreateExtendableRank1DataSet( file, "InvLat_eq", H5T_NATIVE_DOUBLE, &space );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Invariant latitude of vehicle computed from Lambda=acos(sqrt(1.0/Lm_eq)).");
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "Deg." );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );


    // Create BoverBeq Dataset
    DataSet = CreateExtendableRank1DataSet( file, "BoverBeq", H5T_NATIVE_DOUBLE, &space );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Magntiude of Bsc over magnitude of Bmin.");
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "Dimless" );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );

    // Create MlatFromBoverBeq Dataset
    DataSet = CreateExtendableRank1DataSet( file, "MlatFromBoverBeq", H5T_NATIVE_DOUBLE, &space );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Dipole latitude where (B/Beq)_dipole == BoverBeq.");
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "Deg." );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );



    // Create M_used Dataset
    DataSet = CreateExtendableRank1DataSet( file, "M_used", H5T_NATIVE_DOUBLE, &space );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "The magnetic dipole moment that was used to convert magnetic flux to L*. In units of nT.");
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "nT" );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );

    // Create M_ref Dataset
    DataSet = CreateExtendableRank1DataSet( file, "M_ref", H5T_NATIVE_DOUBLE, &space );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "The fixed reference magnetic dipole moment for converting magnetic flux to L*. In units of nT.");
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "nT" );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );

    // Create M_igrf Dataset
    DataSet = CreateExtendableRank1DataSet( file, "M_igrf", H5T_NATIVE_DOUBLE, &space );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Time-dependant magnetic dipole moment (probably shouldn't be used for converting magnetic flux to L*, but it sometimes is). In units of nT.");
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "nT" );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );


    // Create Lstar Dataset
    DataSet = CreateExtendableRank2DataSet( file, "Lstar", m->nAlpha, H5T_NATIVE_DOUBLE, &space );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Generalized Roederer L-shell value (also known as L*).");
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",       "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "DEPEND_1",       "Alpha" );
    sprintf( TmpStr, "%s Lstar, Dimensionless",     ExtModel );
    Lgm_WriteStringAttr( DataSet, "UNITS",          TmpStr );
    Lgm_WriteStringAttr( DataSet, "LABEL",          TmpStr );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",       "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",        "-1E31" );
    Lgm_WriteStringAttr( DataSet, "TYPICAL_MIN",    "1.0" );
    Lgm_WriteStringAttr( DataSet, "TYPICAL_MAX",    "10.0" );
    Lgm_WriteStringAttr( DataSet, "RENDER_TYPE",    "time_series" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",       "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );

    // Create L Dataset
    DataSet = CreateExtendableRank2DataSet( file, "L", m->nAlpha, H5T_NATIVE_DOUBLE, &space );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "McIlwain L-shell value.");
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "DEPEND_1",   "Alpha" );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "Dimensionless" );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );

    // Create Bm Dataset
    DataSet = CreateExtendableRank2DataSet( file, "Bm", m->nAlpha, H5T_NATIVE_DOUBLE, &space );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Magnetic field strength at mirror points for each pitch angle.");
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "DEPEND_1",   "Alpha" );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "nT" );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "log" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );

    // Create I Dataset
    DataSet = CreateExtendableRank2DataSet( file, "I", m->nAlpha, H5T_NATIVE_DOUBLE, &space );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Integral invariant for each pitch angle.");
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "DEPEND_1",   "Alpha" );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "Re" );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );

    // Create K Dataset
    DataSet = CreateExtendableRank2DataSet( file, "K", m->nAlpha, H5T_NATIVE_DOUBLE, &space );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Second Invariant ( I*sqrt(Bm) ) for each pitch angle.");
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "DEPEND_1",   "Alpha" );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "Re G^.5" );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "log" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    Lgm_WriteStringAttr( DataSet, "VALIDMIN",   "0.0" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );

    // Create Sb Dataset
    DataSet = CreateExtendableRank2DataSet( file, "Sb", m->nAlpha, H5T_NATIVE_DOUBLE, &space );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Value of the 'Sb Integral' for equatorially mirroring particles (not generally zero).");
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "DEPEND_1",   "Alpha" );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "Re" );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );

    // Create Tb Dataset
    DataSet = CreateExtendableRank2DataSet( file, "Tb", m->nAlpha, H5T_NATIVE_DOUBLE, &space );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Bounce period for 1 MeV electrons.");
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "DEPEND_1",   "Alpha" );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "s" );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );


    // Create Kappa Dataset
    DataSet = CreateExtendableRank2DataSet( file, "Kappa", m->nAlpha, H5T_NATIVE_DOUBLE, &space );
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", "Kappa parameter for 1MeV electrons -- sqrt( (Minimum Radius of Curvature)/Maximum gyroradius)) (see Büchner, J., and L. M. Zelenyi (1989), Regular and Chaotic Charged Particle Motion in Magnetotaillike Field Reversals, 1. Basic Theory of Trapped Motion, J. Geophys. Res., 94(A9), 11,821-11,842, doi:10.1029/JA094iA09p11821.");
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "DEPEND_1",   "Alpha" );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "dimlesss" );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );

    // Create DriftShellType Dataset
    DataSet = CreateExtendableRank2DataSet( file, "DriftShellType", m->nAlpha, H5T_NATIVE_INT, &space );
    sprintf( TmpStr, "Type of Drift Shell (e.g. %d=CLOSED, %d=CLOSED_SHABANSKY, %d=OPEN, %d=OPEN_SHABANSKY)", LGM_DRIFT_ORBIT_CLOSED, LGM_DRIFT_ORBIT_CLOSED_SHABANSKY, LGM_DRIFT_ORBIT_OPEN, LGM_DRIFT_ORBIT_OPEN_SHABANSKY);
    Lgm_WriteStringAttr( DataSet, "DESCRIPTION", TmpStr );
    Lgm_WriteStringAttr( DataSet, "DEPEND_0",   "IsoTime" );
    Lgm_WriteStringAttr( DataSet, "DEPEND_1",   "Alpha" );
    Lgm_WriteStringAttr( DataSet, "UNITS",      "dimlesss" );
    Lgm_WriteStringAttr( DataSet, "SCALETYP",   "linear" );
    Lgm_WriteStringAttr( DataSet, "FILLVAL",    "-1E31" );
    Lgm_WriteStringAttr( DataSet, "VAR_TYPE",   "data" );
    status  = H5Sclose( space );
    status  = H5Dclose( DataSet );


}


void Lgm_WriteMagEphemDataHdf( hid_t file, int iRow, int i, Lgm_MagEphemData *med ) {

    hid_t   atype;
    herr_t  status;





    // Write String variables
    atype = CreateStrType( 32 );
    LGM_HDF5_EXTEND_RANK1_DATASET( file, "IsoTime",           iRow,                  atype,             &med->H5_IsoTimes[i][0] );         // Write IsoTime
    LGM_HDF5_EXTEND_RANK1_DATASET( file, "FieldLineType",     iRow,                  atype,             &med->H5_FieldLineType[i][0] );    // Write H5_FieldLineType
    LGM_HDF5_EXTEND_RANK1_DATASET( file, "IntModel",          iRow,                  atype,             &med->H5_IntModel[i][0] );         // Write IntModel
    LGM_HDF5_EXTEND_RANK1_DATASET( file, "ExtModel",          iRow,                  atype,             &med->H5_ExtModel[i][0] );         // Write ExtModel
    status  = H5Tclose( atype );

    // Write Non-String variables
    LGM_HDF5_EXTEND_RANK1_DATASET( file, "Date",              iRow,                  H5T_NATIVE_LONG,   &med->H5_Date[i] );                // Write Date
    LGM_HDF5_EXTEND_RANK1_DATASET( file, "Doy",               iRow,                  H5T_NATIVE_INT,    &med->H5_Doy[i] );                 // Write Doy
    LGM_HDF5_EXTEND_RANK1_DATASET( file, "UTC",               iRow,                  H5T_NATIVE_DOUBLE, &med->H5_UTC[i] );                 // Write UTC (hours)
    LGM_HDF5_EXTEND_RANK1_DATASET( file, "JulianDate",        iRow,                  H5T_NATIVE_DOUBLE, &med->H5_JD[i] );                  // Write JD
    LGM_HDF5_EXTEND_RANK1_DATASET( file, "GpsTime",           iRow,                  H5T_NATIVE_DOUBLE, &med->H5_GpsTime[i] );             // Write GpsTime
    LGM_HDF5_EXTEND_RANK1_DATASET( file, "DipoleTiltAngle",   iRow,                  H5T_NATIVE_DOUBLE, &med->H5_TiltAngle[i] );           // Write DipoleTiltAngle
    LGM_HDF5_EXTEND_RANK1_DATASET( file, "InOut",             iRow,                  H5T_NATIVE_INT,    &med->H5_InOut[i] );               // Write InOut
    LGM_HDF5_EXTEND_RANK1_DATASET( file, "OrbitNumber",       iRow,                  H5T_NATIVE_INT,    &med->H5_OrbitNumber[i] );         // Write OrbitNumber
    LGM_HDF5_EXTEND_RANK1_DATASET( file, "SunAngle",          iRow,                  H5T_NATIVE_DOUBLE, &med->H5_SunAngle[i] );            // Write SunAngle

    /*
     * Kludge....
     * Add in Eclipse Flag
     */
    Lgm_Vector u_gei, u_mod;
    Lgm_CTrans *c = Lgm_init_ctrans(0);
// c->ephModel = LGM_EPH_DE;
// TODO:  make the PN model selectable
c->ephModel = LGM_EPH_LOW_ACCURACY;
    int Eclipse_Flag;
    Lgm_Set_Coord_Transforms( med->H5_Date[i], med->H5_UTC[i], c );
    u_gei.x = med->H5_Rgei[i][0]; u_gei.y = med->H5_Rgei[i][1]; u_gei.z = med->H5_Rgei[i][2];
    Lgm_Convert_Coords( &u_gei, &u_mod, GEI2000_TO_MOD, c );
    Eclipse_Flag = Lgm_EarthEclipse( &u_mod, c );
// printf("med->H5_Date[i], med->H5_UTC[i] = %8ld %g\n", med->H5_Date[i], med->H5_UTC[i] );
// printf("u_gei = %g %g %g\n", u_gei.x, u_gei.y, u_gei.z);
// printf("u_mod = %g %g %g\n", u_mod.x, u_mod.y, u_mod.z);
// printf("Eclipse_Flag = %d\n", Eclipse_Flag );
// //double radius = Lgm_Magnitude( &u_gei );
// //if ( (radius > 2.0)&&(Eclipse_Flag>0) ) exit(0);
    Lgm_free_ctrans( c );
    LGM_HDF5_EXTEND_RANK1_DATASET( file, "Eclipse",           iRow, H5T_NATIVE_INT, &Eclipse_Flag );                // Write Eclipse_Flag

    LGM_HDF5_EXTEND_RANK2_DATASET( file, "Rgsm",              iRow, 3,               H5T_NATIVE_DOUBLE, &med->H5_Rgsm[i][0] );             // Write Rgsm
    LGM_HDF5_EXTEND_RANK2_DATASET( file, "Rgeo",              iRow, 3,               H5T_NATIVE_DOUBLE, &med->H5_Rgeo[i][0] );             // Write Rgeo
    LGM_HDF5_EXTEND_RANK2_DATASET( file, "Rsm",               iRow, 3,               H5T_NATIVE_DOUBLE, &med->H5_Rsm[i][0] );              // Write Rsm
    LGM_HDF5_EXTEND_RANK2_DATASET( file, "Rgei",              iRow, 3,               H5T_NATIVE_DOUBLE, &med->H5_Rgei[i][0] );             // Write Rgei
    LGM_HDF5_EXTEND_RANK2_DATASET( file, "Rgse",              iRow, 3,               H5T_NATIVE_DOUBLE, &med->H5_Rgse[i][0] );             // Write Rgse
    LGM_HDF5_EXTEND_RANK2_DATASET( file, "Rgeod_LatLon",      iRow, 2,               H5T_NATIVE_DOUBLE, &med->H5_Rgeod_LatLon[i][0] );     // Write Rgeod_LatLon
    LGM_HDF5_EXTEND_RANK1_DATASET( file, "Rgeod_Height",      iRow,                  H5T_NATIVE_DOUBLE, &med->H5_Rgeod_Height[i] );        // Write Rgeod_Height
    LGM_HDF5_EXTEND_RANK1_DATASET( file, "CDMAG_MLAT",        iRow,                  H5T_NATIVE_DOUBLE, &med->H5_CDMAG_MLAT[i] );          // Write CDMAG_MLAT
    LGM_HDF5_EXTEND_RANK1_DATASET( file, "CDMAG_MLON",        iRow,                  H5T_NATIVE_DOUBLE, &med->H5_CDMAG_MLON[i] );          // Write CDMAG_MLAT
    LGM_HDF5_EXTEND_RANK1_DATASET( file, "CDMAG_MLT",         iRow,                  H5T_NATIVE_DOUBLE, &med->H5_CDMAG_MLT[i] );           // Write CDMAG_MLAT
    LGM_HDF5_EXTEND_RANK1_DATASET( file, "CDMAG_R",           iRow,                  H5T_NATIVE_DOUBLE, &med->H5_CDMAG_R[i] );             // Write CDMAG_MLAT
    LGM_HDF5_EXTEND_RANK1_DATASET( file, "EDMAG_MLAT",        iRow,                  H5T_NATIVE_DOUBLE, &med->H5_EDMAG_MLAT[i] );          // Write CDMAG_MLAT
    LGM_HDF5_EXTEND_RANK1_DATASET( file, "EDMAG_MLON",        iRow,                  H5T_NATIVE_DOUBLE, &med->H5_EDMAG_MLON[i] );          // Write CDMAG_MLAT
    LGM_HDF5_EXTEND_RANK1_DATASET( file, "EDMAG_MLT",         iRow,                  H5T_NATIVE_DOUBLE, &med->H5_EDMAG_MLT[i] );           // Write CDMAG_MLAT
    LGM_HDF5_EXTEND_RANK1_DATASET( file, "EDMAG_R",           iRow,                  H5T_NATIVE_DOUBLE, &med->H5_EDMAG_R[i] );             // Write CDMAG_MLAT
    LGM_HDF5_EXTEND_RANK1_DATASET( file, "Kp",                iRow,                  H5T_NATIVE_DOUBLE, &med->H5_Kp[i] );                  // Write Kp
    LGM_HDF5_EXTEND_RANK1_DATASET( file, "Dst",               iRow,                  H5T_NATIVE_DOUBLE, &med->H5_Dst[i] );                 // Write Dst
    LGM_HDF5_EXTEND_RANK2_DATASET( file, "Bsc_gsm",           iRow, 4,               H5T_NATIVE_DOUBLE, &med->H5_Bsc_gsm[i][0] );          // Write Bsc_gsm
    LGM_HDF5_EXTEND_RANK1_DATASET( file, "S_sc_to_pfn",       iRow,                  H5T_NATIVE_DOUBLE, &med->H5_S_sc_to_pfn[i] );         // Write S_sc_to_pfn
    LGM_HDF5_EXTEND_RANK1_DATASET( file, "S_sc_to_pfs",       iRow,                  H5T_NATIVE_DOUBLE, &med->H5_S_sc_to_pfs[i] );         // Write S_sc_to_pfs
    LGM_HDF5_EXTEND_RANK1_DATASET( file, "S_pfs_to_Bmin",     iRow,                  H5T_NATIVE_DOUBLE, &med->H5_S_pfs_to_Bmin[i] );       // Write S_pfs_to_Bmin
    LGM_HDF5_EXTEND_RANK1_DATASET( file, "S_Bmin_to_sc",      iRow,                  H5T_NATIVE_DOUBLE, &med->H5_S_Bmin_to_sc[i] );        // Write S_Bmin_to_sc
    LGM_HDF5_EXTEND_RANK1_DATASET( file, "S_total",           iRow,                  H5T_NATIVE_DOUBLE, &med->H5_S_total[i] );             // Write S_total
    LGM_HDF5_EXTEND_RANK1_DATASET( file, "d2B_ds2",           iRow,                  H5T_NATIVE_DOUBLE, &med->H5_d2B_ds2[i] );             // Write d2B_ds2
    LGM_HDF5_EXTEND_RANK1_DATASET( file, "Sb0",               iRow,                  H5T_NATIVE_DOUBLE, &med->H5_Sb0[i] );                 // Write Sb0
    LGM_HDF5_EXTEND_RANK1_DATASET( file, "RadiusOfCurv",      iRow,                  H5T_NATIVE_DOUBLE, &med->H5_RadiusOfCurv[i] );        // Write RadiusOfCurv
    LGM_HDF5_EXTEND_RANK2_DATASET( file, "Pfn_geo",           iRow, 3,               H5T_NATIVE_DOUBLE, &med->H5_Pfn_geo[i][0] );          // Write Pfn_geo
    LGM_HDF5_EXTEND_RANK2_DATASET( file, "Pfn_gsm",           iRow, 3,               H5T_NATIVE_DOUBLE, &med->H5_Pfn_gsm[i][0] );          // Write Pfn_gsm
    LGM_HDF5_EXTEND_RANK2_DATASET( file, "Pfn_geod_LatLon",   iRow, 2,               H5T_NATIVE_DOUBLE, &med->H5_Pfn_geod_LatLon[i][0] );  // Write Pfn_geod_LatLon
    LGM_HDF5_EXTEND_RANK1_DATASET( file, "Pfn_geod_Height",   iRow,                  H5T_NATIVE_DOUBLE, &med->H5_Pfn_geod_Height[i] );     // Write Pfn_geod_Height
    LGM_HDF5_EXTEND_RANK1_DATASET( file, "Pfn_CD_MLAT",       iRow,                  H5T_NATIVE_DOUBLE, &med->H5_Pfn_CD_MLAT[i] );         // Write Pfn_CD_MLAT
    LGM_HDF5_EXTEND_RANK1_DATASET( file, "Pfn_CD_MLON",       iRow,                  H5T_NATIVE_DOUBLE, &med->H5_Pfn_CD_MLON[i] );         // Write Pfn_CD_MLON
    LGM_HDF5_EXTEND_RANK1_DATASET( file, "Pfn_CD_MLT",        iRow,                  H5T_NATIVE_DOUBLE, &med->H5_Pfn_CD_MLT[i] );          // Write Pfn_CD_MLT
    LGM_HDF5_EXTEND_RANK1_DATASET( file, "Pfn_ED_MLAT",       iRow,                  H5T_NATIVE_DOUBLE, &med->H5_Pfn_ED_MLAT[i] );         // Write Pfn_ED_MLAT
    LGM_HDF5_EXTEND_RANK1_DATASET( file, "Pfn_ED_MLON",       iRow,                  H5T_NATIVE_DOUBLE, &med->H5_Pfn_ED_MLON[i] );         // Write Pfn_ED_MLON
    LGM_HDF5_EXTEND_RANK1_DATASET( file, "Pfn_ED_MLT",        iRow,                  H5T_NATIVE_DOUBLE, &med->H5_Pfn_ED_MLT[i] );          // Write Pfn_ED_MLT
    LGM_HDF5_EXTEND_RANK2_DATASET( file, "Bfn_geo",           iRow, 4,               H5T_NATIVE_DOUBLE, &med->H5_Bfn_geo[i][0] );          // Write Bfn_geo
    LGM_HDF5_EXTEND_RANK2_DATASET( file, "Bfn_gsm",           iRow, 4,               H5T_NATIVE_DOUBLE, &med->H5_Bfn_gsm[i][0] );          // Write Bfn_gsm
    LGM_HDF5_EXTEND_RANK1_DATASET( file, "Loss_Cone_Alpha_n", iRow,                  H5T_NATIVE_DOUBLE, &med->H5_LossConeAngleN[i] );      // Write Loss_Cone_Alpha_n
    LGM_HDF5_EXTEND_RANK2_DATASET( file, "Pfs_geo",           iRow, 3,               H5T_NATIVE_DOUBLE, &med->H5_Pfs_geo[i][0] );          // Write Pfs_geo
    LGM_HDF5_EXTEND_RANK2_DATASET( file, "Pfs_gsm",           iRow, 3,               H5T_NATIVE_DOUBLE, &med->H5_Pfs_gsm[i][0] );          // Write Pfs_gsm
    LGM_HDF5_EXTEND_RANK2_DATASET( file, "Pfs_geod_LatLon",   iRow, 2,               H5T_NATIVE_DOUBLE, &med->H5_Pfs_geod_LatLon[i][0] );  // Write Pfs_geod_LatLon
    LGM_HDF5_EXTEND_RANK1_DATASET( file, "Pfs_geod_Height",   iRow,                  H5T_NATIVE_DOUBLE, &med->H5_Pfs_geod_Height[i] );     // Write Pfs_geod_Height
    LGM_HDF5_EXTEND_RANK1_DATASET( file, "Pfs_CD_MLAT",       iRow,                  H5T_NATIVE_DOUBLE, &med->H5_Pfs_CD_MLAT[i] );         // Write Pfs_CD_MLAT
    LGM_HDF5_EXTEND_RANK1_DATASET( file, "Pfs_CD_MLON",       iRow,                  H5T_NATIVE_DOUBLE, &med->H5_Pfs_CD_MLON[i] );         // Write Pfs_CD_MLON
    LGM_HDF5_EXTEND_RANK1_DATASET( file, "Pfs_CD_MLT",        iRow,                  H5T_NATIVE_DOUBLE, &med->H5_Pfs_CD_MLT[i] );          // Write Pfs_CD_MLT
    LGM_HDF5_EXTEND_RANK1_DATASET( file, "Pfs_ED_MLAT",       iRow,                  H5T_NATIVE_DOUBLE, &med->H5_Pfs_ED_MLAT[i] );         // Write Pfs_ED_MLAT
    LGM_HDF5_EXTEND_RANK1_DATASET( file, "Pfs_ED_MLON",       iRow,                  H5T_NATIVE_DOUBLE, &med->H5_Pfs_ED_MLON[i] );         // Write Pfs_ED_MLON
    LGM_HDF5_EXTEND_RANK1_DATASET( file, "Pfs_ED_MLT",        iRow,                  H5T_NATIVE_DOUBLE, &med->H5_Pfs_ED_MLT[i] );          // Write Pfs_ED_MLT
    LGM_HDF5_EXTEND_RANK2_DATASET( file, "Bfs_geo",           iRow, 4,               H5T_NATIVE_DOUBLE, &med->H5_Bfs_geo[i][0] );          // Write Bfs_geo
    LGM_HDF5_EXTEND_RANK2_DATASET( file, "Bfs_gsm",           iRow, 4,               H5T_NATIVE_DOUBLE, &med->H5_Bfs_gsm[i][0] );          // Write Bfs_gsm
    LGM_HDF5_EXTEND_RANK1_DATASET( file, "Loss_Cone_Alpha_s", iRow,                  H5T_NATIVE_DOUBLE, &med->H5_LossConeAngleS[i] );      // Write Loss_Cone_Alpha_s
    LGM_HDF5_EXTEND_RANK2_DATASET( file, "Pmin_gsm",          iRow, 3,               H5T_NATIVE_DOUBLE, &med->H5_Pmin_gsm[i][0] );         // Write Pmin_gsm
    LGM_HDF5_EXTEND_RANK2_DATASET( file, "Bmin_gsm",          iRow, 4,               H5T_NATIVE_DOUBLE, &med->H5_Bmin_gsm[i][0] );         // Write Bmin_gsm
    LGM_HDF5_EXTEND_RANK1_DATASET( file, "Lsimple",           iRow,                  H5T_NATIVE_DOUBLE, &med->H5_Lsimple[i] );             // Write Lsimple
    LGM_HDF5_EXTEND_RANK1_DATASET( file, "InvLat",            iRow,                  H5T_NATIVE_DOUBLE, &med->H5_InvLat[i] );              // Write InvLat
    LGM_HDF5_EXTEND_RANK1_DATASET( file, "Lm_eq",             iRow,                  H5T_NATIVE_DOUBLE, &med->H5_Lm_eq[i] );               // Write Lm_eq
    LGM_HDF5_EXTEND_RANK1_DATASET( file, "InvLat_eq",         iRow,                  H5T_NATIVE_DOUBLE, &med->H5_InvLat_eq[i] );           // Write InvLat_eq
    LGM_HDF5_EXTEND_RANK1_DATASET( file, "BoverBeq",          iRow,                  H5T_NATIVE_DOUBLE, &med->H5_BoverBeq[i] );            // Write BoverBeq
    LGM_HDF5_EXTEND_RANK1_DATASET( file, "MlatFromBoverBeq",  iRow,                  H5T_NATIVE_DOUBLE, &med->H5_MlatFromBoverBeq[i] );    // Write MlatFromBoverBeq
    LGM_HDF5_EXTEND_RANK1_DATASET( file, "M_used",            iRow,                  H5T_NATIVE_DOUBLE, &med->H5_M_used[i] );              // Write M_used
    LGM_HDF5_EXTEND_RANK1_DATASET( file, "M_ref",             iRow,                  H5T_NATIVE_DOUBLE, &med->H5_M_ref[i] );               // Write M_ref
    LGM_HDF5_EXTEND_RANK1_DATASET( file, "M_igrf",            iRow,                  H5T_NATIVE_DOUBLE, &med->H5_M_igrf[i] );              // Write M_igrf
    LGM_HDF5_EXTEND_RANK2_DATASET( file, "Lstar",             iRow, med->H5_nAlpha,  H5T_NATIVE_DOUBLE, &med->H5_Lstar[i][0] );            // Write Lstar
    LGM_HDF5_EXTEND_RANK2_DATASET( file, "Sb",                iRow, med->H5_nAlpha,  H5T_NATIVE_DOUBLE, &med->H5_Sb[i][0] );               // Write Sb
    LGM_HDF5_EXTEND_RANK2_DATASET( file, "Tb",                iRow, med->H5_nAlpha,  H5T_NATIVE_DOUBLE, &med->H5_Tb[i][0] );               // Write Tb
    LGM_HDF5_EXTEND_RANK2_DATASET( file, "Kappa",             iRow, med->H5_nAlpha,  H5T_NATIVE_DOUBLE, &med->H5_Kappa[i][0] );            // Write Kappa
    LGM_HDF5_EXTEND_RANK2_DATASET( file, "DriftShellType",    iRow, med->H5_nAlpha,  H5T_NATIVE_INT,    &med->H5_DriftShellType[i][0] );   // Write DriftShellType
    LGM_HDF5_EXTEND_RANK2_DATASET( file, "L",                 iRow, med->H5_nAlpha,  H5T_NATIVE_DOUBLE, &med->H5_L[i][0] );                // Write L
    LGM_HDF5_EXTEND_RANK2_DATASET( file, "Bm",                iRow, med->H5_nAlpha,  H5T_NATIVE_DOUBLE, &med->H5_Bm[i][0] );               // Write Bm
    LGM_HDF5_EXTEND_RANK2_DATASET( file, "I",                 iRow, med->H5_nAlpha,  H5T_NATIVE_DOUBLE, &med->H5_I[i][0] );                // Write I
    LGM_HDF5_EXTEND_RANK2_DATASET( file, "K",                 iRow, med->H5_nAlpha,  H5T_NATIVE_DOUBLE, &med->H5_K[i][0] );                // Write K

}
