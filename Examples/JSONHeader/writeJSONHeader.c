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
  Lgm_metadata_initvar(&energy_meta, // the variable to store it
		       1,            // Dimension
		       1,            // has data (bool)
		       "ENERGY"      // Name 
		       );
  // add an attribute to the variable (all arguments must be strings)
  Lgm_metadata_addStringAttr(&energy_meta, // the variable to add it to
			     "UNITS",      // name
			     "MeV",        // value
			     0             // boolean it is not an array 
			     );
  Lgm_metadata_addStringAttr(&energy_meta, "LABEL", "Energy [MeV]",0);
  Lgm_metadata_addStringAttr(&energy_meta, "COMMENT", "I AM A COMMENT ON THE ENERGY", 0);

  // create the JSON header
  //   this is a variadic function, number of vars and all the vars
  json_header = Lgm_metadata_JSONheader(1, &energy_meta);
  printf("\n\n**************************\n");
  printf("%s", json_header);

  printf("**************************\n\n");

  Lgm_metadata_initvar(&counts_meta, 1, 1, "COUNTS");
  Lgm_metadata_addStringAttr(&counts_meta, "LABEL", "Counts", 0);

  // make a string with both vars in it, order is specified in call
  json_header = Lgm_metadata_JSONheader(2, &energy_meta, &counts_meta);


  printf("\n\n**************************\n");
  printf("%s", json_header);

  printf("**************************\n\n");

  // create a variable with no data
  Lgm_metadata_initvar(&meta_meta, 1, 0, "INSTITUTION");
  Lgm_metadata_addStringAttr(&meta_meta, "NAME", "Los Alamos National Laboratory", 0);
  Lgm_metadata_addStringAttr(&meta_meta, "GROUP", "ISR-1 Space Science and Applications", 0);

  // order is again set by the order in this call
  json_header = Lgm_metadata_JSONheader(3, &meta_meta, &energy_meta, &counts_meta);

  printf("\n\n**************************\n");
  printf("%s", json_header);

  printf("**************************\n\n");

  Lgm_metadata_addStringAttrArray2(&meta_meta, "SPACECRAFT", 2, "RBSPA", "RBSPB");

  json_header = Lgm_metadata_JSONheader(3, &meta_meta, &energy_meta, &counts_meta);
  printf("\n\n**************************\n");
  printf("%s", json_header);

  printf("**************************\n\n");

  {
    int invals[] = {1,2,3};
    double invals2[] = {1.1, 2.2, 3.3};
    Lgm_metadata_addIntAttr( &meta_meta, "THREE_INTS", 3, invals);
    Lgm_metadata_addDoubleAttr( &meta_meta, "THREE_DOUBLES", 3, invals2);
  }

  json_header = Lgm_metadata_JSONheader(3, &meta_meta, &energy_meta, &counts_meta);
  printf("\n\n**************************\n");
  printf("%s", json_header);

  printf("**************************\n\n");



  return (0);
}

