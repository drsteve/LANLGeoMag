

#ifndef LGM_JSONHEADER_H
#define LGM_JSONHEADER_H

typedef struct Lgm_json_attr {
  char *name;
  char *value;
} Lgm_json_attr;

typedef struct Lgm_json_variable {
  long start_column;
  int dimension;
  char* name;
  Lgm_json_attr *attributes;
  long n_attrs;
} Lgm_json_variable;

void Lgm_json_initvar(long start_column, int dumension, char* name, Lgm_json_variable *var);

void Lgm_json_addAttr(char *name, char *value, Lgm_json_variable *var);

char* Lgm_json_toString(Lgm_json_variable *var);


#endif

