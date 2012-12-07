
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "Lgm/Lgm_Metadata.h"


/* #  "DOY":              { "DESCRIPTION": "Ordinal Day of Year.", */
/* #                               "NAME": "DOY", */
/* #                              "TITLE": "Day Of Year", */
/* #                              "LABEL": "DOY", */
/* #                       "START_COLUMN": 2, */
/* #                              "UNITS": "Days", */
/* #                          "VALID_MIN": 0, */
/* #                          "VALID_MAX": 366 }, */


#define MAXHEADERLEN 10000


char* Lgm_metadata_toJSON(Lgm_metadata_variable *var) {
  char Buffer[MAXHEADERLEN]; 
  char *outstr;
  int i, j, k;
  size_t L;
  char *p;

  //TODO make the # optional
  // start the block
  p = Buffer;
  p += sprintf(p, "# {\n");
  
  // add the name and start attrs
  p += sprintf(p, "#  \"%s\":\t{ ", var->name);
  // go through each attribute and add it to the string
  for (i=0;i<var->n_attrs;i++) {
    // add the name
    p += sprintf(p, " \"%s\": ", (var->attributes)[i].name);
    // add the value
    p += sprintf(p, "\"%s\"  ", (var->attributes)[i].value);
    if (i+1 < var->n_attrs)
      p += sprintf(p, ", \n# \t\t " );
    else
      p += sprintf(p, "}, \n");
  }

  outstr = malloc( (strlen(Buffer) + 1)*sizeof(char));
  L = strlen(Buffer);
  strncpy(outstr, Buffer, L);
  outstr[L] = '\0';
  return (outstr);
}


void Lgm_metadata_initvar(int start_column, int dimension, char* name, Lgm_metadata_variable *var) {
  char int_str[80];  // this limits to 999 columms or dimensions
  var->dimension = dimension;
  var->start_column = start_column;
  var->name = name;
  var->n_attrs = 0;
  var->attributes = NULL;

  sprintf(int_str, "%d", var->start_column);
  Lgm_metadata_addAttr("START_COLUMN", int_str, var);
  sprintf(int_str, "%d", var->dimension);
  Lgm_metadata_addAttr("DIMENSION", int_str, var);

}

void Lgm_metadata_addAttr(char *name, char *value, Lgm_metadata_variable *var) {
  Lgm_metadata_attr atr;
  Lgm_metadata_variable newvar;
  Lgm_metadata_attr* old_attr;

  atr.name = name;
  atr.value = value;

  old_attr = var->attributes;

  var->n_attrs++;
  var->attributes = (Lgm_metadata_attr*) calloc(var->n_attrs, sizeof(Lgm_metadata_attr));
  // put the ones back that were there
  {
    int i;
    for (i=0; i<var->n_attrs-1;i++) { // the -1 is since we have a ++ above
      var->attributes[i].name = old_attr[i].name;
      var->attributes[i].value = old_attr[i].value;
    }
  }
  var->attributes[var->n_attrs - 1].name = malloc((1+strlen(name))*sizeof(char));
  strcpy(var->attributes[var->n_attrs - 1].name, name);
  var->attributes[var->n_attrs - 1].value = malloc((1+strlen(value))*sizeof(char));
  strcpy(var->attributes[var->n_attrs - 1].value, value);

  printf("%s  n_attrs: %d  %d bytes\n", var->name, var->n_attrs, sizeof(var->attributes));
 printf("  added %s:%s\n", var->attributes[var->n_attrs - 1].name, var->attributes[var->n_attrs - 1].value);
  

  /* Lgm_metadata_initvar(var->start_column, var->dimension, var->name,  newvar); */

 free(old_attr);

}


