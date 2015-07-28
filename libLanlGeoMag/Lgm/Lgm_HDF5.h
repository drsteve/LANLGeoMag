#ifndef LGM_HDF5_H
#define LGM_HDF5_H

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <ctype.h>
#include "Lgm/Lgm_DynamicMemory.h"
#define H5_NO_DEPRECATED_SYMBOLS
#include <hdf5.h>



/*
 * Macros for extending Rank 2 HDF5 datasets by one "row".
 */
#define LGM_HDF5_EXTEND_RANK1_DATASET( File, DataSetName, nRow, Type, buf ) {\
{\
    hsize_t     Rank, NewSize[4], Offset[4], SlabSize[4];\
    hid_t       DataSet, DataSpace, MemSpace;\
    herr_t      status __attribute__((unused));\
    \
    Rank        = 1;\
    NewSize[0]  = nRow+1;\
    Offset[0]   = nRow;\
    SlabSize[0] = 1;\
    \
    DataSet   = H5Dopen( File, DataSetName, H5P_DEFAULT );\
    status    = H5Dextend( DataSet, &NewSize[0] );\
    \
    DataSpace = H5Dget_space( DataSet );\
    status    = H5Sselect_hyperslab( DataSpace, H5S_SELECT_SET, Offset, NULL, SlabSize, NULL );\
    \
    MemSpace  = H5Screate_simple( Rank, SlabSize, NULL );\
    status    = H5Dwrite( DataSet, Type, MemSpace, DataSpace, H5P_DEFAULT, buf );\
    \
    status = H5Sclose( MemSpace );\
    status = H5Dclose( DataSet );\
    status = H5Sclose( DataSpace );\
    \
}\
}

#define LGM_HDF5_EXTEND_RANK2_DATASET( File, DataSetName, nRow, nCol, Type, buf ) {\
{\
    hsize_t     Rank, NewSize[4], Offset[4], SlabSize[4];\
    hid_t       DataSet, DataSpace, MemSpace;\
    herr_t      status __attribute__((unused));\
    \
    Rank        = 2;\
    NewSize[0]  = nRow+1;  NewSize[1]  = nCol;\
    Offset[0]   = nRow;    Offset[1]   = 0;\
    SlabSize[0] = 1;       SlabSize[1] = nCol;\
    \
    DataSet   = H5Dopen( File, DataSetName, H5P_DEFAULT );\
    status    = H5Dextend( DataSet, &NewSize[0] );\
    \
    DataSpace = H5Dget_space( DataSet );\
    status    = H5Sselect_hyperslab( DataSpace, H5S_SELECT_SET, Offset, NULL, SlabSize, NULL );\
    \
    MemSpace  = H5Screate_simple( Rank, SlabSize, NULL );\
    status    = H5Dwrite( DataSet, Type, MemSpace, DataSpace, H5P_DEFAULT, buf );\
    \
    status = H5Sclose( MemSpace );\
    status = H5Dclose( DataSet );\
    status = H5Sclose( DataSpace );\
    \
}\
}


#define LGM_HDF5_EXTEND_RANK3_DATASET( File, DataSetName, nRow, nCol, nDepth, Type, buf ) {\
{\
    hsize_t     Rank, NewSize[4], Offset[4], SlabSize[4];\
    hid_t       DataSet, DataSpace, MemSpace;\
    herr_t      status __attribute__((unused));\
    \
    Rank        = 3;\
    NewSize[0]  = nRow+1;  NewSize[1]  = nCol;      NewSize[2]  = nDepth;\
    Offset[0]   = nRow;    Offset[1]   = 0;         Offset[2]   = 0;\
    SlabSize[0] = 1;       SlabSize[1] = nCol;      SlabSize[2] = nDepth;\
    \
    DataSet   = H5Dopen( File, DataSetName, H5P_DEFAULT );\
    status    = H5Dextend( DataSet, &NewSize[0] );\
    \
    DataSpace = H5Dget_space( DataSet );\
    status    = H5Sselect_hyperslab( DataSpace, H5S_SELECT_SET, Offset, NULL, SlabSize, NULL );\
    \
    MemSpace  = H5Screate_simple( Rank, SlabSize, NULL );\
    status    = H5Dwrite( DataSet, Type, MemSpace, DataSpace, H5P_DEFAULT, buf );\
    \
    status = H5Sclose( MemSpace );\
    status = H5Dclose( DataSet );\
    status = H5Sclose( DataSpace );\
    \
}\
}





/*
 * Macros for overwriting a "row" in a Rank 1 and 2  HDF5 dataset.
 */
#define LGM_HDF5_OVERWRITE_RANK1_DATASET( File, DataSetName, iRow, Type, buf ) {\
{\
    hsize_t     Rank, NewSize[4], Offset[4], SlabSize[4];\
    hid_t       DataSet, DataSpace, MemSpace;\
    herr_t      status;\
    \
    Rank        = 1;\
    Offset[0]   = iRow;\
    SlabSize[0] = 1;\
    \
    DataSet   = H5Dopen( File, DataSetName, H5P_DEFAULT );\
    DataSpace = H5Dget_space( DataSet );\
    status    = H5Sselect_hyperslab( DataSpace, H5S_SELECT_SET, Offset, NULL, SlabSize, NULL );\
    MemSpace  = H5Screate_simple( Rank, SlabSize, NULL );\
    status    = H5Dwrite( DataSet, Type, MemSpace, DataSpace, H5P_DEFAULT, buf );\
    \
    status = H5Sclose( MemSpace );\
    status = H5Dclose( DataSet );\
    status = H5Sclose( DataSpace );\
    \
}\
}

#define LGM_HDF5_OVERWRITE_RANK2_DATASET( File, DataSetName, iRow, nCol, Type, buf ) {\
{\
    hsize_t     Rank, NewSize[4], Offset[4], SlabSize[4];\
    hid_t       DataSet, DataSpace, MemSpace;\
    herr_t      status;\
    \
    Rank        = 2;\
    Offset[0]   = iRow;    Offset[1]   = 0;\
    SlabSize[0] = 1;       SlabSize[1] = nCol;\
    \
    DataSet   = H5Dopen( File, DataSetName, H5P_DEFAULT );\
    DataSpace = H5Dget_space( DataSet );\
    status    = H5Sselect_hyperslab( DataSpace, H5S_SELECT_SET, Offset, NULL, SlabSize, NULL );\
    MemSpace  = H5Screate_simple( Rank, SlabSize, NULL );\
    status    = H5Dwrite( DataSet, Type, MemSpace, DataSpace, H5P_DEFAULT, buf );\
    \
    status = H5Sclose( MemSpace );\
    status = H5Dclose( DataSet );\
    status = H5Sclose( DataSpace );\
    \
}\
}



/*
 * Macros for writing simple Rank 1 and 2  HDF5 dataset.
 */
#define LGM_HDF5_WRITE_SIMPLE_RANK1_DATASET( File, DataSetName, n, Type, buf ) {\
{\
    int         Rank;\
    hsize_t     Offset[4], SlabSize[4];\
    hid_t       DataSpace, DataSet, MemSpace;\
    \
    Rank = 1;\
    Offset[0]   = 0;\
    SlabSize[0] = n;\
\
    DataSet   = H5Dopen( file, DataSetName, H5P_DEFAULT );\
    DataSpace = H5Dget_space( DataSet );\
    status    = H5Sselect_hyperslab( DataSpace, H5S_SELECT_SET, Offset, NULL, SlabSize, NULL );\
\
    MemSpace  = H5Screate_simple( Rank, SlabSize, NULL );\
    status    = H5Dwrite( DataSet, Type, MemSpace, DataSpace, H5P_DEFAULT, buf );\
    status    = H5Sclose( MemSpace );\
    status    = H5Sclose( DataSpace );\
    status    = H5Dclose( DataSet );\
    \
}\
}

#define LGM_HDF5_WRITE_SIMPLE_RANK2_DATASET( File, DataSetName, n, m, Type, buf ) {\
{\
    int         Rank;\
    hsize_t     Offset[4], SlabSize[4];\
    hid_t       DataSpace, DataSet, MemSpace;\
    \
    Rank = 2;\
    Offset[0]   = 0;  Offset[1] = 0;\
    SlabSize[0] = n;  SlabSize[1] = m;\
\
    DataSet   = H5Dopen( file, DataSetName, H5P_DEFAULT );\
    DataSpace = H5Dget_space( DataSet );\
    status    = H5Sselect_hyperslab( DataSpace, H5S_SELECT_SET, Offset, NULL, SlabSize, NULL );\
\
    MemSpace  = H5Screate_simple( Rank, SlabSize, NULL );\
    status    = H5Dwrite( DataSet, Type, MemSpace, DataSpace, H5P_DEFAULT, buf );\
    status    = H5Sclose( MemSpace );\
    status    = H5Sclose( DataSpace );\
    status    = H5Dclose( DataSet );\
    \
}\
}


#define LGM_HDF5_WRITE_SIMPLE_RANK1_STRING_DATASET( File, DataSetName, n, Type, buf ) {\
{\
    int         i, Rank;\
    hsize_t     Offset[4], SlabSize[4];\
    hid_t       DataSpace, DataSet, MemSpace;\
    \
    Rank = 1;\
    SlabSize[0] = 1;\
\
    DataSet   = H5Dopen( file, DataSetName, H5P_DEFAULT );\
    DataSpace = H5Dget_space( DataSet );\
\
    MemSpace  = H5Screate_simple( Rank, SlabSize, NULL );\
    for (i=0;i<n;i++) {\
        Offset[0] = i;\
        status    = H5Sselect_hyperslab( DataSpace, H5S_SELECT_SET, Offset, NULL, SlabSize, NULL );\
        status    = H5Dwrite( DataSet, Type, MemSpace, DataSpace, H5P_DEFAULT, &buf[i][0] );\
    }\
    status    = H5Sclose( MemSpace );\
    status    = H5Sclose( DataSpace );\
    status    = H5Dclose( DataSet );\
    \
}\
}






char     **Get_StringDataset_1D( hid_t file, char *Str, hsize_t *Dims );
double    *Get_DoubleDataset_1D( hid_t file, char *Str, hsize_t *Dims );
double   **Get_DoubleDataset_2D( hid_t file, char *Str, hsize_t *Dims );
double  ***Get_DoubleDataset_3D( hid_t file, char *Str, hsize_t *Dims );
double ****Get_DoubleDataset_4D( hid_t file, char *Str, hsize_t *Dims );

void Lgm_WriteStringAttr( hid_t DataSet, char *AttrNAme, char *Str );
void Lgm_WriteDoubleAttr( hid_t DataSet, char *AttrNAme, double Val );

hid_t   CreateExtendableRank1DataSet( hid_t File, char *DataSetName, hid_t Type, hid_t *DataSpace );
hid_t   CreateExtendableRank2DataSet( hid_t File, char *DataSetName, int Cols, hid_t Type, hid_t *DataSpace );
hid_t   CreateExtendableRank3DataSet( hid_t File, char *DataSetName, int nCol, int nDepth, hid_t Type, hid_t *DataSpace );
hid_t   CreateStrType( int StrLength );


#endif
