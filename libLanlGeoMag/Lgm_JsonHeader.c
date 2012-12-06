
#include <stdlib.h>
#include "Lgm/Lgm_JsonHeader.h"


/* typedef struct Lgm_json_attr { */
/*   char *name; */
/*   char *value; */
/* } json_attr; */

/* typedef struct Lgm_json_variable { */
/*   long start_column; */
/*   int dimension; */
/*   char* name; */
/*   Lgm_json_attr *attributes; */
/*   long n_attrs; */
/* } json_variable; */

/* void Lgm_json_initvar(long start_column, int dumension, char* name, Lgm_json_variable *var); */

/* void Lgm_json_addAttr(char *name, char *value, json_variable *var); */

/* void Lgm_json_toString(Lgm_json_variable *var, char* stringout); */



void Lgm_json_initvar(long start_column, int dimension, char* name, Lgm_json_variable *var) {
  var->dimension = dimension;
  var->start_column = start_column;
  var->name = name;
}

void Lgm_json_addAttr(char *name, char *value, Lgm_json_variable *var) {
  Lgm_json_attr atr;
  Lgm_json_variable newvar;

  atr.name = name;
  atr.value = value;

  var->n_attrs++;
  
  var->attributes = (Lgm_json_attr*) calloc(var->n_attrs, sizeof(Lgm_json_attr));
 
  

  /* Lgm_json_initvar(var->start_column, var->dimension, var->name,  newvar); */



}


