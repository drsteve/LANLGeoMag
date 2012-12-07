

#ifndef LGM_JSONHEADER_H
#define LGM_JSONHEADER_H


typedef struct Lgm_metadata_attr {
  char *name;
  char *value;
} Lgm_metadata_attr;

typedef struct Lgm_metadata_variable {
  long start_column;
  int dimension;
  long n_attrs;
  char* name;
  Lgm_metadata_attr *attributes;
} Lgm_metadata_variable;

void Lgm_metadata_initvar(int start_column, int dumension, char* name, Lgm_metadata_variable *var);

void Lgm_metadata_addAttr(char *name, char *value, Lgm_metadata_variable *var);

char* Lgm_metadata_toJSON(Lgm_metadata_variable *var);


#endif

