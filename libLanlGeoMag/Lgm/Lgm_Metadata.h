
#ifndef LGM_JSONHEADER_H
#define LGM_JSONHEADER_H

#define Lgm_metadata_SUCCESS 1
#define Lgm_metadata_FAILURE 0

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Lgm/Lgm_DynamicMemory.h"
#include "Lgm/Lgm_Metadata.h"

typedef struct Lgm_metadata_attr {
  int array;  // boolean for if this is an array of metadata
  int string; // boolean for if this is a string attr
  char *name;
  char *value;
} Lgm_metadata_attr;

typedef struct Lgm_metadata_variable {
  long start_column;
  int dimension;
  long n_attrs;
  int data;
  char* name;
  Lgm_metadata_attr *attributes;
} Lgm_metadata_variable;


#define Lgm_metadata_MAXHEADERLEN 100000
#define Lgm_metadata_MAXSINGLEVARHEADER 10000
#define Lgm_metadata_TABCHAR "    "

void Lgm_metadata_initvar(Lgm_metadata_variable *var, int dimension, int data, char* name);

char *Lgm_metadata_JSONheader(int n_vars, ...);

char* Lgm_metadata_toJSON(Lgm_metadata_variable *var, short last, char comment);

char *Lgm_metadata_stringArrayToString(int len, char** instring);

char *Lgm_metadata_intArrayToString(int len, int* invals);

char *Lgm_metadata_doubleArrayToString(int len, double* invals);

int Lgm_metadata_addIntAttr(Lgm_metadata_variable *var, char *name, int len, int *invals);

int Lgm_metadata_addDoubleAttr(Lgm_metadata_variable *var, char *name, int len, double *invals);

int Lgm_metadata_addStringAttrArray(Lgm_metadata_variable *var, char *name, int len, char **invals );

int Lgm_metadata_addStringAttrArray2(Lgm_metadata_variable *var, char *name, int len, ... );

int Lgm_metadata_addStringAttr(Lgm_metadata_variable *var, char *name, char * value, int array);

int Lgm_metadata_AttrInVar(Lgm_metadata_variable *var, char *name);
  






#endif

