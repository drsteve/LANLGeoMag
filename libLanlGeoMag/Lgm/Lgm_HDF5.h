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

char     **Get_StringDataset_1D( hid_t file, char *Str, hsize_t *Dims );
double    *Get_DoubleDataset_1D( hid_t file, char *Str, hsize_t *Dims );
double   **Get_DoubleDataset_2D( hid_t file, char *Str, hsize_t *Dims );
double  ***Get_DoubleDataset_3D( hid_t file, char *Str, hsize_t *Dims );
double ****Get_DoubleDataset_4D( hid_t file, char *Str, hsize_t *Dims );

void Lgm_WriteStringAttr( hid_t DataSet, char *AttrNAme, char *Str );
void Lgm_WriteDoubleAttr( hid_t DataSet, char *AttrNAme, double Val );

#endif
