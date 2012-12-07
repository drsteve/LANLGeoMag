
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

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

/* char* Lgm_json_toString(Lgm_json_variable *var); */


    /* fprintf( fp, "#  \"Pfs_ED_MLT\":       { \"DESCRIPTION\": \"Magnetic Local Time of Southern Footpoint in Eccentric Dipole Coordinates.\",\n"); */
    /* fprintf( fp, "#                               \"NAME\": \"Pfs_ED_MLT\",\n"); */
    /* fprintf( fp, "#                              \"TITLE\": \"ED MLT of Southern Footpoint\",\n"); */
    /* fprintf( fp, "#                              \"LABEL\": \"%s ED MLT (Hours)\",\n", ExtModel ); */
    /* fprintf( fp, "#                              \"UNITS\": \"Hours\",\n"); */
    /* fprintf( fp, "#                       \"START_COLUMN\": %d,\n", nCol++); */
    /* fprintf( fp, "#                          \"VALID_MIN\": 0.0,\n"); */
    /* fprintf( fp, "#                          \"VALID_MAX\": 24.0,\n"); */
    /* fprintf( fp, "#                         \"FILL_VALUE\": -1e31 },\n"); */
    /* fprintf( fp, "#\n"); */


/* #  "DOY":              { "DESCRIPTION": "Ordinal Day of Year.", */
/* #                               "NAME": "DOY", */
/* #                              "TITLE": "Day Of Year", */
/* #                              "LABEL": "DOY", */
/* #                       "START_COLUMN": 2, */
/* #                              "UNITS": "Days", */
/* #                          "VALID_MIN": 0, */
/* #                          "VALID_MAX": 366 }, */


#define MAXHEADERLEN 10000
#define eos(s) ((s)+strlen(s))


char* Lgm_json_toString(Lgm_json_variable *var) {
  char Buffer[MAXHEADERLEN]; 
  char *outstr;
  int i, j, k;
  // start the block
  sprintf(eos(Buffer), "# {\n");
  // add the name and start attrs
  sprintf(eos(Buffer), "#  \"%s\":\t{ ", var->name);
  // go through each attribute and add it to the string
  for (i=0;i<var->n_attrs;i++) {
    // add the name
    sprintf(eos(Buffer), " \"%s\": ", (var->attributes)[i].name);
    // add the value
    sprintf(eos(Buffer), "\"%s\"  ", (var->attributes)[i].value);
    if (i+1 < var->n_attrs)
      sprintf(eos(Buffer), ", \n# \t\t " );
    else
      sprintf(eos(Buffer), "}, \n");
}

  outstr = malloc(strlen(Buffer + 1));
  strncpy(outstr, Buffer, strlen(Buffer));
  return (outstr);
}


void Lgm_json_initvar(long start_column, int dimension, char* name, Lgm_json_variable *var) {
  char *int_str = malloc(3+1);  // this limits to 999 columms or dimensions
  var->dimension = dimension;
  var->start_column = start_column;
  var->name = name;
  var->n_attrs = 0;
  var->attributes = NULL;

  sprintf(int_str, "%d", start_column);
  Lgm_json_addAttr("START_COLUMN", int_str, var);
  sprintf(int_str, "%d", dimension);
  Lgm_json_addAttr("DIMENSION", int_str, var);

}

void Lgm_json_addAttr(char *name, char *value, Lgm_json_variable *var) {
  Lgm_json_attr atr;
  Lgm_json_variable newvar;
  Lgm_json_attr* old_attr;

  atr.name = name;
  atr.value = value;

  old_attr = var->attributes;

  var->n_attrs++;
  var->attributes = (Lgm_json_attr*) calloc(var->n_attrs, sizeof(Lgm_json_attr));
  // put the ones back that were there
  {
    int i;
    for (i=0; i<var->n_attrs-1;i++) { // the -1 is since we have a ++ above
      var->attributes[i].name = old_attr[i].name;
      var->attributes[i].value = old_attr[i].value;
    }
  }

  var->attributes[var->n_attrs - 1].name = name;
  var->attributes[var->n_attrs - 1].value = value;

  printf("%s  n_attrs: %d  %d bytes\n", var->name, var->n_attrs, sizeof(var->attributes));
 printf("  added %s:%s\n", var->attributes[var->n_attrs - 1].name, var->attributes[var->n_attrs - 1].value);
  

  /* Lgm_json_initvar(var->start_column, var->dimension, var->name,  newvar); */



}


