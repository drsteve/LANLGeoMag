#include <stdio.h>

#include <Lgm_Metadata.h>

int main(void) {
  
  /* Create some data that we would like to write out a header for
   * - specify the metadata, dimensions, and if data exists then write it out
   */

  /* need a meatadat variable for each entry */
  Lgm_metadata_variable energy_meta;
  Lgm_metadata_variable counts_meta;
  Lgm_metadata_variable meta_meta;
  // and a string to put the answer
  char* json_header; 

  // create a 1 dimensional variable, with data, a name, and the variable
  Lgm_metadata_initvar(1, // Dimension
		       1, // has data (bool)
		       "ENERGY", // Name
		       &energy_meta // the variable to store it
		       );
  // add an attribute to the variable (all arguments must be strings)
  Lgm_metadata_addAttr("UNITS", "MeV", &energy_meta);
  Lgm_metadata_addAttr("LABEL", "Energy [MeV]", &energy_meta);
  Lgm_metadata_addAttr("COMMENT", "I AM A COMMENT ON THE ENERGY", &energy_meta);

  // create the JSON header
  //   this is a variadic function, number of vars and all the vars
  json_header = Lgm_metadata_JSONheader(1, &energy_meta);
  printf("\n\n**************************\n");
  printf("%s", json_header);

  printf("**************************\n\n");  

  Lgm_metadata_initvar(1, 1, "COUNTS", &counts_meta);
  Lgm_metadata_addAttr("LABEL", "Counts", &counts_meta);

  // make a string with both vars in it, order is specified in call
  json_header = Lgm_metadata_JSONheader(2, &energy_meta, &counts_meta);


  printf("\n\n**************************\n");
  printf("%s", json_header);

  printf("**************************\n\n");

  // create a variable with no data
  Lgm_metadata_initvar(1, 0, "INSTITUTION", &meta_meta);
  Lgm_metadata_addAttr("NAME", "Los Alamos National Laboratory", &meta_meta);
  Lgm_metadata_addAttr("GROUP", "ISR-1 Space Science and Applications", &meta_meta);

  // order is again set by the order in this call
  json_header = Lgm_metadata_JSONheader(3, &meta_meta, &energy_meta, &counts_meta);

  printf("\n\n**************************\n");
  printf("%s", json_header);

  printf("**************************\n\n");


  return (0);
}

