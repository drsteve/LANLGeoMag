#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <string.h>
#include "Lgm/Lgm_HDF5.h"

/**
 *  \file
 *      This file contains simple convenience routines to read/write
 *      datasets/attributes to/from HDF5 files.
 *
 *
 */
void Lgm_WriteStringAttr( hid_t DataSet, char *AttrNAme, char *Str ) {

    hid_t   space, type, Attr;
    hsize_t Dims[1] = {1};
    herr_t  status;

    space = H5Screate_simple(1, Dims, NULL);
    type    = H5Tcopy( H5T_C_S1 );
    status  = H5Tset_size( type, strlen(Str) );
    status  = H5Tset_strpad( type, H5T_STR_NULLPAD );
    status  = H5Tset_cset( type, H5T_CSET_ASCII );
    Attr    = H5Acreate( DataSet, AttrNAme, type, space, H5P_DEFAULT, H5P_DEFAULT);
    status  = H5Awrite( Attr, type, &Str[0] );
    status  = H5Aclose( Attr );
    status  = H5Tclose( type );
    status  = H5Sclose( space );

    return;
}

/** here is a C++ version just so we have it **/
/* #include "H5Cpp.h" */
/* void h5_addStringAttribute(DataSet ds, std::string attrName, std::string attrValue) { */
/*   StrType str_type(0, attrValue.length()); */
/*   DataSpace att_space(H5S_SCALAR); */
/*   Attribute att = ds.createAttribute( attrName, str_type, att_space ); */
/*   att.write( str_type, attrValue ); */
/* } */


void Lgm_WriteDoubleAttr( hid_t DataSet, char *AttrNAme, double Val ) {

    hid_t   space, Attr;
    hsize_t Dims[1] = {1};
    herr_t  status;

    space = H5Screate_simple(1, Dims, NULL);
    Attr    = H5Acreate( DataSet, AttrNAme, H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT);
    status  = H5Awrite( Attr, H5T_NATIVE_DOUBLE, &Val );
    status  = H5Aclose( Attr );
    status  = H5Sclose( space );

    return;
}



/**
 *  \brief
 *      Gets an array of strings from an HDF5 dataset.
 *
 *  \details
 *      Given a file handle obtained with H5Fopen(), and a string indicating
 *      the dataset, this routine determines: 1) the size and number of
 *      strings, 2) allocates memory for them, and 3) reads the data from the
 *      hdf5 file into the string. Typical calling sequence would be;
 *
 *          \code
 *          hsize_t dims[2];
 *          char **IsoTimes;
 *          if ( (file = H5Fopen( InFile, H5F_ACC_RDONLY, H5P_DEFAULT )) > 0 ) {
 *              //Read ISO time
 *              IsoTimes = Get_StringDataset_1D( file, "/Epoch", dims);
 *              HDFclose( file );
 *          }
 *          ... Do stuff with IsoTimes array ...
 *          LGM_ARRAY_2D_FREE( IsoTimes );
 *          \endcode
 *
 *
 *      \param[in]      file       An hdf5 file handle (e.g. as returned by a call
 *                                 like 'file = H5Fopen( "MyDataFile.h5", H5F_ACC_RDONLY, H5P_DEFAULT );').
 *      \param[in]      DataSet    The dataset in the HDF5 file (e.g. "/IsoTime" ).
 *      \param[out]     Dims       The size of the array. (Dims[0] is the number of
 *                                 strings, Dims[1] is the fixed length of each string (including the
 *                                 null-term character).
 *
 *      \return         Returns a pointer to an array of null-terminated
 *                      strings. The user must free this with LGM_ARRAY_2D_FREE( ).
 *
 *      \note           A 1D array of strings in C is actually a 2D array of
 *                      chars. This is why this uses 2D memory allocation macros.
 *
 *      \author         Mike Henderson
 *      \date           2011
 *
 */
char **Get_StringDataset_1D( hid_t file, char *DataSet, hsize_t *Dims ) {

    size_t      size;
    H5T_cset_t  cset;
    H5T_str_t   strpad;
    hid_t       dataset, dataspace, datatype;
    herr_t      status;
    int         rank, status_n, i, j;
    char        **buf1;
    char        **buf2;

    dataset   = H5Dopen( file, DataSet, H5P_DEFAULT );
    dataspace = H5Dget_space( dataset );
    rank      = H5Sget_simple_extent_ndims( dataspace );
    if ( rank != 1 ) {
        printf("Dataset is not rank 1. Rank = %d\n", rank);
        exit(0);
    }
    status_n = H5Sget_simple_extent_dims( dataspace, Dims, NULL );
    datatype = H5Dget_type(dataset);
    size     = H5Tget_size(datatype);
    strpad   = H5Tget_strpad(datatype);
    cset     = H5Tget_cset(datatype);

    LGM_ARRAY_2D(  buf1, (int)Dims[0], (int)size, char );

    hid_t atype = H5Tcopy (H5T_C_S1);
    status  = H5Tset_size( atype, size);
    status  = H5Tset_strpad( atype, strpad );
    status  = H5Tset_cset( atype, cset );
    Dims[1] = size;

    status   = H5Dread( dataset, atype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &buf1[0][0] );

    if ( strpad != H5T_STR_NULLTERM ) {

        LGM_ARRAY_2D(  buf2, (int)Dims[0], (int)(size+1), char );

        for ( i=0; i<Dims[0]; i++ ) {
            for ( j=0; j<size; j++ ) {
                buf2[i][j] = buf1[i][j];
            }
            buf2[i][size] = '\0';
        }
        H5Dclose(dataset);
        H5Sclose(dataspace);
        LGM_ARRAY_2D_FREE( buf1 );
        return( buf2 );

    } else {

        H5Dclose(dataset);
        H5Sclose(dataspace);
        return( buf1 );

    }

}


/**
 *  \brief
 *      Gets a 1D array of doubles from an HDF5 dataset.
 *
 *  \details
 *      Given a file handle obtained with H5Fopen(), and a string indicating
 *      the dataset, this routine determines: 1) the size of the 1D array
 *      2) allocates memory for them, and 3) reads the data from the
 *      hdf5 file into the array. Typical calling sequence would be;
 *
 *
 *          \code
 *          hsize_t dims[1];
 *          double *MyData;
 *          if ( (file = H5Fopen( InFile, H5F_ACC_RDONLY, H5P_DEFAULT )) > 0 ) {
 *              //Read a double array
 *              MyData = Get_DoubleDataset_1D( file, "/Epoch", dims);
 *              HDFclose( file );
 *          }
 *          ... Do stuff with MyData array ...
 *          LGM_ARRAY_1D_FREE( MyData );
 *          \endcode
 *
 *
 *
 *      \param[in]      file   An hdf5 file handle (e.g. as returned by a call
 *                             like 'file = H5Fopen( "MyDataFile.h5", H5F_ACC_RDONLY, H5P_DEFAULT );').
 *      \param[in]      DataSet    The dataset in the HDF5 file (e.g. "/IsoTime" ).
 *      \param[out]     Dims   The size of the array. (Dims[0] is the number of
 *                             elements in the 1D array.)
 *
 *      \return         Returns a pointer to an array of doubles.
 *                      The user must free this with LGM_ARRAY_1D_FREE( ).
 *
 *      \author         Mike Henderson
 *      \date           2011
 *
 */
double *Get_DoubleDataset_1D( hid_t file, char *DataSet, hsize_t *Dims ) {

    hid_t     dataset, dataspace;
    herr_t    status;
    int       rank, status_n;
    double  *buf;

    dataset   = H5Dopen( file, DataSet, H5P_DEFAULT );
    dataspace = H5Dget_space( dataset );
    rank      = H5Sget_simple_extent_ndims( dataspace );
    if ( rank != 1 ) {
        printf("Dataset is not rank 1. Rank = %d\n", rank);
        exit(0);
    }
    status_n = H5Sget_simple_extent_dims( dataspace, Dims, NULL );

    LGM_ARRAY_1D(  buf, Dims[0], double );
    status   = H5Dread( dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &buf[0] );

    H5Dclose(dataset);
    H5Sclose(dataspace);

    return( buf );

}

/**
 *  \brief
 *      Gets a 2D array of doubles from an HDF5 dataset.
 *
 *  \details
 *      Given a file handle obtained with H5Fopen(), and a string indicating
 *      the dataset, this routine determines: 1) the size of the 2D array
 *      2) allocates memory for them, and 3) reads the data from the
 *      hdf5 file into the array. Typical calling sequence would be;
 *
 *
 *          \code
 *          hsize_t dims[2];
 *          double **MyData;
 *          if ( (file = H5Fopen( InFile, H5F_ACC_RDONLY, H5P_DEFAULT )) > 0 ) {
 *              //Read a double array
 *              MyData = Get_DoubleDataset_1D( file, "/Epoch", dims);
 *              HDFclose( file );
 *          }
 *          ... Do stuff with MyData array ...
 *          LGM_ARRAY_2D_FREE( MyData );
 *          \endcode
 *
 *
 *
 *      \param[in]      file   An hdf5 file handle (e.g. as returned by a call
 *                             like 'file = H5Fopen( "MyDataFile.h5", H5F_ACC_RDONLY, H5P_DEFAULT );').
 *      \param[in]      DataSet    The dataset in the HDF5 file (e.g. "/IsoTime" ).
 *      \param[out]     Dims   The size of the array. (Dims contains the size in each dimension.)
 *
 *      \return         Returns a pointer to a 2D array of doubles.
 *                      The user must free this with LGM_ARRAY_2D_FREE( ).
 *
 *      \author         Mike Henderson
 *      \date           2011
 *
 */
double **Get_DoubleDataset_2D( hid_t file, char *DataSet, hsize_t *Dims ) {

    hid_t     dataset, dataspace;
    herr_t    status;
    int       rank, status_n;
    double  **buf;

    dataset   = H5Dopen( file, DataSet, H5P_DEFAULT );
    dataspace = H5Dget_space( dataset );
    rank      = H5Sget_simple_extent_ndims( dataspace );
    if ( rank != 2 ) {
        printf("Dataset is not rank 2. Rank = %d\n", rank);
        exit(0);
    }
    status_n = H5Sget_simple_extent_dims( dataspace, Dims, NULL );

    LGM_ARRAY_2D(  buf, Dims[0], Dims[1], double );
    status   = H5Dread( dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &buf[0][0] );

    H5Dclose(dataset);
    H5Sclose(dataspace);

    return( buf );

}

double ***Get_DoubleDataset_3D( hid_t file, char *DataSet, hsize_t *Dims ) {

    hid_t     dataset, dataspace;
    herr_t    status;
    int       rank, status_n;
    double  ***buf;

    dataset   = H5Dopen( file, DataSet, H5P_DEFAULT );
    dataspace = H5Dget_space( dataset );
    rank      = H5Sget_simple_extent_ndims( dataspace );
    if ( rank != 3 ) {
        printf("Dataset is not rank 3. Rank = %d\n", rank);
        exit(0);
    }
    status_n = H5Sget_simple_extent_dims( dataspace, Dims, NULL );

    LGM_ARRAY_3D(  buf, Dims[0], Dims[1], Dims[2], double );
    status   = H5Dread( dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &buf[0][0][0] );

    H5Dclose(dataset);
    H5Sclose(dataspace);

    return( buf );

}

double ****Get_DoubleDataset_4D( hid_t file, char *DataSet, hsize_t *Dims ) {

    hid_t     dataset, dataspace;
    herr_t    status;
    int       rank, status_n;
    double  ****buf;

    dataset   = H5Dopen( file, DataSet, H5P_DEFAULT );
    dataspace = H5Dget_space( dataset );
    rank      = H5Sget_simple_extent_ndims( dataspace );
    if ( rank != 4 ) {
        printf("Dataset is not rank 4. Rank = %d\n", rank);
        exit(0);
    }
    status_n = H5Sget_simple_extent_dims( dataspace, Dims, NULL );

    LGM_ARRAY_4D(  buf, Dims[0], Dims[1], Dims[2], Dims[3], double );
    status   = H5Dread( dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &buf[0][0][0][0] );

    H5Dclose(dataset);
    H5Sclose(dataspace);

    return( buf );

}


/*
 *  Creates an extendible dataset (with zero "rows"). Returns handle to
 *  dataset. Also returns handle to the dataspace.  Both of these need to be
 *  closed by user when done. (E.g. H5Sclose( DataSpace ); and H5Dclose( DataSet ); )
 */
hid_t   CreateExtendableRank1DataSet( hid_t File, char *DataSetName, hid_t Type, hid_t *DataSpace ){

    int             Rank;
    hsize_t         Dims[4], MaxDims[4], ChunkDims[4];
    hid_t           cparms, status, DataSet;


    Rank         = 1;
    Dims[0]      = 0;
    MaxDims[0]   = H5S_UNLIMITED;
    ChunkDims[0] = 1;
    *DataSpace   = H5Screate_simple( Rank, Dims, MaxDims );

    cparms  = H5Pcreate( H5P_DATASET_CREATE );
    status  = H5Pset_chunk( cparms, Rank, ChunkDims );

    DataSet = H5Dcreate( File, DataSetName, Type, *DataSpace, H5P_DEFAULT, cparms, H5P_DEFAULT );


    return( DataSet );

}


hid_t   CreateExtendableRank2DataSet( hid_t File, char *DataSetName, int Cols, hid_t Type, hid_t *DataSpace ){

    int             Rank;
    hsize_t         Dims[4], MaxDims[4], ChunkDims[4];
    hid_t           cparms, status, DataSet;


    Rank         = 2;
    Dims[0]      = 0;               Dims[1]      = Cols;
    MaxDims[0]   = H5S_UNLIMITED;   MaxDims[1]   = Cols;
    ChunkDims[0] = 1;               ChunkDims[1] = Cols;
    *DataSpace   = H5Screate_simple( Rank, Dims, MaxDims );

    cparms  = H5Pcreate( H5P_DATASET_CREATE );
    status  = H5Pset_chunk( cparms, Rank, ChunkDims );

    DataSet = H5Dcreate( File, DataSetName, Type, *DataSpace, H5P_DEFAULT, cparms, H5P_DEFAULT );


    return( DataSet );

}



/*
 * Convenience rotuinem to create a string data type
 * close it with     H5Tclose( atype );
 */
hid_t CreateStrType( int StrLength ) {

    hid_t   status, atype;

    atype   = H5Tcopy( H5T_C_S1 );
    status  = H5Tset_size( atype, StrLength );
    status  = H5Tset_strpad( atype, H5T_STR_NULLTERM );
    status  = H5Tset_cset( atype, H5T_CSET_ASCII );

    return( atype );

}




hid_t   CreateSimpleRank1DataSet( hid_t File, char *DataSetName, int n, hid_t Type, hid_t *DataSpace ){

    int             Rank;
    hsize_t         Dims[4];
    hid_t           status, DataSet;


    Rank         = 1;
    Dims[0]      = n;
    *DataSpace   = H5Screate_simple( Rank, Dims, NULL );
    DataSet      = H5Dcreate( File, DataSetName, Type, *DataSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );


    return( DataSet );

}



hid_t   CreateSimpleRank2DataSet( hid_t File, char *DataSetName, int n, int m, hid_t Type, hid_t *DataSpace ){

    int             Rank;
    hsize_t         Dims[4];
    hid_t           status, DataSet;


    Rank         = 2;
    Dims[0]      = n;   Dims[1]      = m;
    *DataSpace   = H5Screate_simple( Rank, Dims, NULL );
    DataSet      = H5Dcreate( File, DataSetName, Type, *DataSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );


    return( DataSet );

}
