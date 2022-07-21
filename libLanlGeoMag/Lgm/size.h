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

#include "Lgm_MagModelInfo.h"
#include "Lgm_CTrans.h"
#include "Lgm_Octree.h"

#include <gsl/gsl_spline.h>


int size_t_size(void);
int size_MagModelInfo(void);
int size_CTrans(void);
int size_Vector(void);
int size_DateTime(void);
int size_gsl_interp_accel(void);
int size_gsl_interp_type(void);
int size_gsl_interp(void);
int size_gsl_spline(void);
int size_Lgm_OctreeData(void);
int size_Lgm_OctreeCell(void);
int size_pQueue(void);
int size_Lgm_LeapSeconds(void); 
