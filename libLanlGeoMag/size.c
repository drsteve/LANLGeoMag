/*******************************************************************************
 *  size.c  - simple uility routine needed by the ctypes pyton wrapper
 *  to know the sizes of variable types in C on the platform the library
 *  is installed.
 *
 *  Sized needed:
 *  size_t
 *  MagModelInfo
 *
 *  Brian Larsen
 *  balarsen@lanl.gov
 *  Los Alamos National Lab
 *
 *  V1 23Dec2010 (BAL)
 ******************************************************************************/

#include <stddef.h>

#include "Lgm/Lgm_MagModelInfo.h"
#include "Lgm/Lgm_CTrans.h"
#include "Lgm/Lgm_Octree.h"


#include <gsl/gsl_spline.h>


int size_t_size(void){
    return sizeof(size_t);
}

int size_MagModelInfo(void) {
    return sizeof(Lgm_MagModelInfo);
}

int size_CTrans(void) {
    return sizeof(Lgm_CTrans);
}

int size_Vector(void) {
    return sizeof(Lgm_Vector);
}

int size_DateTime(void) {
    return sizeof(Lgm_DateTime);
}

int size_gsl_interp_accel(void) {
    return sizeof(gsl_interp_accel);
}


int size_gsl_interp_type(void) {
    return sizeof(gsl_interp_type);
}


int size_gsl_interp(void) {
    return sizeof(gsl_interp);
}

int size_gsl_spline(void) {
    return sizeof(gsl_spline);
}

int size_Lgm_OctreeData(void) {
    return sizeof(Lgm_OctreeData);
}

int size_Lgm_OctreeCell(void) {
    return sizeof(Lgm_OctreeCell);
}

int size_pQueue(void) {
    return sizeof(pQueue);
}

int size_Lgm_LeapSeconds(void) {
    return sizeof(Lgm_LeapSeconds);
}



//int size_(void) {
//    return sizeof();
//}
//




/*
int size_DateAndTime(void) {
    return sizeof(Lgm_DateAndTime);
}*/



