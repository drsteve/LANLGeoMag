#ifndef CDF_UTILS_H
#define CDF_UTILS_H

#include <cdf.h>
#include <stdlib.h>
#include <stdio.h>


/*
 *   Function prototypes
 */
void StatusHandler( CDFstatus status );
void DynamicErrorHandle( CDFid id );
void DataTypeToStr( long dataType, char *str );
void CDF_Info( char *Filename, int verbosity );


#endif
