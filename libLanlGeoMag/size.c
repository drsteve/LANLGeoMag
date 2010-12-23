

#include <stddef.h>

/*******************************************************************************
 *  size.c  - simple uility routine needed by the ctypes pyton wrapper
 *  to know the sizes of varialbe types in C on the platform the library
 *  is installed.
 *
 *  Sized needed:
 *  size_t
 *
 *  Brian Larsen
 *  balarsen@lanl.gov
 *  Los Alamos National Lab
 *
 *  V1 23Dec2010 (BAL)
 ******************************************************************************/

int size_t_size(void){
    return sizeof(size_t);
}
