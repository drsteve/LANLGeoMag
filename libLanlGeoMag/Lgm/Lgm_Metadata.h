

#ifndef LGM_JSONHEADER_H
#define LGM_JSONHEADER_H

#include <stdarg.h>

#define Lgm_metadata_SUCCESS 1
#define Lgm_metadata_FAILURE 0


typedef struct Lgm_metadata_attr {
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

void Lgm_metadata_initvar(int dimension, int data, char* name, Lgm_metadata_variable *var);

int Lgm_metadata_addAttr(char *name, char *value, Lgm_metadata_variable *var);

char* Lgm_metadata_toJSON(Lgm_metadata_variable *var, short last, char comment);

char *Lgm_metadata_JSONheader(int n_vars, ...);


#endif

