
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "Lgm/Lgm_Metadata.h"
#include "Lgm/Lgm_DynamicMemory.h"

#define Lgm_metadata_MAXHEADERLEN 100000
#define Lgm_metadata_MAXSINGLEVARHEADER 10000
#define Lgm_metadata_TABCHAR "    "

/******************************************************************************/

void Lgm_metadata_initvar(Lgm_metadata_variable *var, int dimension, int data, char* name) {
  var->dimension = dimension;
  var->name = name;
  var->n_attrs = 0;
  var->attributes = NULL;
  var->data = data;
  Lgm_metadata_addIntAttr(var, "DIMENSION", 1, &dimension); // the & since it wants a pointer
}

/******************************************************************************/

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
    int current_col=0;
    char int_str[80]; 

    // start the block
    p = Buffer;
    p += sprintf(p, "%c {\n", '#');
    
    for (i=0; i<n_vars; i++) {
      var_tmp = va_arg(hv, Lgm_metadata_variable* );
      if (var_tmp->data) {
	Lgm_metadata_addIntAttr(var_tmp, "START_COLUMN", 1, &current_col); // the & since it wants a pointer
	current_col += var_tmp->dimension;
      }

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

/******************************************************************************/

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
    char sep;
    int i;
    for (i=0;i<var->n_attrs;i++) {
      // what sep are we using " or nothing
      if ((var->attributes)[i].value[0] == '[') // from Lgm_metadata_addIntAttr()
	sep = ' ';
      else
	sep = '"';
      // add the name
      p += sprintf(p, " %c%s%c: ", sep, (var->attributes)[i].name, sep);
      // add the value
      p += sprintf(p, "%c%s%c  ", sep, (var->attributes)[i].value, sep);
      // what ending do I need?
      if (i+1 < var->n_attrs)
	p += sprintf(p, ",\n%c %s%s", comment, Lgm_metadata_TABCHAR, Lgm_metadata_TABCHAR );
      else
	p += sprintf(p, "}");
    }
  }

  L = strlen(Buffer);
  outstr = malloc( (L + 1)*sizeof(char));
  strncpy(outstr, Buffer, L);
  outstr[L] = '\0';
  return (outstr);
}

/******************************************************************************/

char *Lgm_metadata_stringArraytoString(int len, char** instring) {
  char Buffer[Lgm_metadata_MAXHEADERLEN];
  char *outstr;
  char *p;
  size_t L;
  int i;

  p = Buffer;
  for (i=0; i<len-1; i++){
    p += sprintf(p, "%s, ", instring[i]);
  }
  p += sprintf(p, "%s, ", instring[i]);  // i should have the right value in it

  L = strlen(Buffer);
  outstr = malloc( (L + 1)*sizeof(char));
  strncpy(outstr, Buffer, L);
  outstr[L] = '\0';                                                                                     
  return (outstr);
}

/******************************************************************************/

char *Lgm_metadata_intArrayToString(int len, int* invals) {
  char Buffer[Lgm_metadata_MAXHEADERLEN];
  char *outstr;
  char *p;
  size_t L;
  int i;

  p = Buffer;
  for (i=0; i<len-1; i++){
    p += sprintf(p, "%d, ", invals[i]);
  }
  p += sprintf(p, "%d, ", invals[i]);  // i should have the right value in it

  L = strlen(Buffer);
  outstr = malloc( (L + 1)*sizeof(char));
  strncpy(outstr, Buffer, L);
  outstr[L] = '\0';                                                                                     
  return (outstr);
}

/******************************************************************************/

char *Lgm_metadata_doubleArrayToString(int len, int* invals) {
  char Buffer[Lgm_metadata_MAXHEADERLEN];
  char *outstr;
  char *p;
  size_t L;
  int i;

  p = Buffer;
  for (i=0; i<len-1; i++){
    p += sprintf(p, "%lf, ", invals[i]);
  }
  p += sprintf(p, "%lf, ", invals[i]);  // i should have the right value in it

  L = strlen(Buffer);
  outstr = malloc( (L + 1)*sizeof(char));
  strncpy(outstr, Buffer, L);
  outstr[L] = '\0';                                                                                     
  return (outstr);
}

/******************************************************************************/

int Lgm_metadata_addIntAttr(Lgm_metadata_variable *var, char *name, int len, int *invals) {
  char *outstr;
  // make sure the name is not already there, if it is ignore it
  if (Lgm_metadata_AttrInVar(var, name))
    return (Lgm_metadata_FAILURE); // it is there we are done

  outstr = Lgm_metadata_intArrayToString(len,invals);

  return (Lgm_metadata_addStringAttr(var, name, outstr, len) );
}

/******************************************************************************/

int Lgm_metadata_addStringAttrArray(Lgm_metadata_variable *var, char *name, int len, char **invals ) {
  char *outstr;
  // make sure the name is not already there, if it is ignore it
  if (Lgm_metadata_AttrInVar(var, name))
    return (Lgm_metadata_FAILURE); // it is there we are done

  outstr = Lgm_metadata_stringArrayToString(len, invals);

  return (Lgm_metadata_addStringAttr(var, name, outstr, len) );
}

/******************************************************************************/

int Lgm_metadata_addStringAttrArray2(Lgm_metadata_variable *var, char *name, int len, ... ) {
  char **outstrarr;
  va_list hv;
  int i;

  // make sure the name is not already there, if it is ignore it
  if (Lgm_metadata_AttrInVar(var, name))
    return (Lgm_metadata_FAILURE); // it is there we are done

  va_start(hv, len);
  LGM_ARRAY_2D( outstrarr,len, 80, char ); // this 80 is a limitiation here!!!
  for (i=0;i<len;i++) {
    outstrarr[i] = va_arg(hv, char*);
  }
  va_end(hv);
  
  return (Lgm_metadata_addStringAttrArray(var, name, len, outstrarr));
}

/******************************************************************************/

int Lgm_metadata_addStringAttr(Lgm_metadata_variable *var, char *name, char * value, int array) {
  Lgm_metadata_attr atr;
  Lgm_metadata_variable newvar;
  Lgm_metadata_attr* old_attr;

  // make sure the name is not already there, if it is ignore this call
  if (Lgm_metadata_AttrInVar(var, name))
    return (Lgm_metadata_FAILURE); // it is there we are done

  atr.name = name;
  atr.value = value;
  if (array)
    atr.array = 1;
  else
    atr.array = 0;

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

  /* Lgm_metadata_initvar(var->start_column, var->dimension, var->name,  newvar); */
 free(old_attr);
 return (Lgm_metadata_SUCCESS);
}


/******************************************************************************/

int Lgm_metadata_AttrInVar(Lgm_metadata_variable *var, char *name) {
  
  // make sure the name is not already there, if it is ignore it 
  int i;
  for (i=0; i<var->n_attrs; i++) {
    if (strcmp(var->attributes[i].name, name) == 0) { // there are equal
      return (1); // name is in the variable
    }
  }
  return (0);  // name is not in the variable
}



