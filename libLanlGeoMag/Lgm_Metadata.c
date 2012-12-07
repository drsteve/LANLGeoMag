
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "Lgm/Lgm_Metadata.h"


#define Lgm_metadata_MAXHEADERLEN 100000
#define Lgm_metadata_MAXSINGLEVARHEADER 10000
#define Lgm_metadata_TABCHAR "    "

char *Lgm_metadata_JSONheader(int n_vars, ...) {
  va_list hv;
  char Buffer[Lgm_metadata_MAXHEADERLEN];
  char *outstr;
  size_t L;

  va_start(hv, n_vars);
  
  // loop over n_vars vars making them all into JSON headers
  {
    int i;
    char *p;
    char *tmp_str;
    Lgm_metadata_variable* var_tmp;

    // start the block
    p = Buffer;
    p += sprintf(p, "%c {\n", '#');
    
    for (i=0; i<n_vars; i++) {
      var_tmp = va_arg(hv, Lgm_metadata_variable* );

      tmp_str = Lgm_metadata_toJSON(var_tmp, i == (n_vars-1), '#');
      p += sprintf(p, "%s", tmp_str);
      if (i < n_vars-1)
	p += sprintf(p, ",\n");
      else 
	p += sprintf(p, "\n");
	
    }
    p += sprintf(p, "%c } # ENDJSON\n", '#');
  }
  va_end(hv);
  L = strlen(Buffer);
  outstr = malloc( (L + 1)*sizeof(char));
  strncpy(outstr, Buffer, L);
  outstr[L] = '\0';
  return (outstr);
}

char* Lgm_metadata_toJSON(Lgm_metadata_variable *var, short last, char comment) {
  char Buffer[Lgm_metadata_MAXSINGLEVARHEADER]; 
  char *outstr;
  size_t L;
  char *p;

  p = Buffer;

  // add the name and start attrs
  p += sprintf(p, "%c  \"%s\":%s{ ", comment, var->name, Lgm_metadata_TABCHAR);
  // go through each attribute and add it to the string
  {
    int i;
    for (i=0;i<var->n_attrs;i++) {
      // add the name
      p += sprintf(p, " \"%s\": ", (var->attributes)[i].name);
      // add the value
      p += sprintf(p, "\"%s\"  ", (var->attributes)[i].value);
      // what ending do I need?
      if (i+1 < var->n_attrs)
	p += sprintf(p, ",\n%c %s%s", comment, Lgm_metadata_TABCHAR, Lgm_metadata_TABCHAR );
      else
	p += sprintf(p, "}");
      /* if (!last) */
      /* 	p += sprintf(p, ",\n"); */
      /* else */
      /* 	p += sprintf(p, "\n"); */
    }
  }

  L = strlen(Buffer);
  outstr = malloc( (L + 1)*sizeof(char));
  strncpy(outstr, Buffer, L);
  outstr[L] = '\0';
  return (outstr);
}


void Lgm_metadata_initvar(int start_column, int dimension, char* name, Lgm_metadata_variable *var) {
  char int_str[80]; 
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

 free(old_attr);

}


