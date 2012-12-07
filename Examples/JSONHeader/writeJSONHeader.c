
#include <Lgm_Metadata.h>
#include <Lgm_Utils.h>
#include <stdio.h>

/* typedef struct Lgm_metadata_variable { */
/*   long start_column; */
/*   int dimension; */
/*   char* name; */
/*   json_attr *attributes; */
/*   long n_attrs; */
/* } json_variable; */


int main(void) {
  
  /* Create some data that we would like to write out a header for
   * - specify the metadata, dimensions, and start_column then write it out
   */

  double energy[10];
  double counts[10];
  Lgm_metadata_variable energy_meta;
  Lgm_metadata_variable counts_meta;
  char* json_header; 

  // populate the energy and counts
  Lgm_LinSpace(1, 100, 10, energy);
  Lgm_LogSpace(100, 1000, 10, counts);

  Lgm_metadata_initvar(0, 1, "ENERGY", &energy_meta);
  Lgm_metadata_addAttr("UNITS", "MeV", &energy_meta);
  Lgm_metadata_addAttr("LABEL", "Energy [MeV]", &energy_meta);
  Lgm_metadata_addAttr("COMMENT", "I AM A COMMENT ON THE ENERGY", &energy_meta);

  json_header = Lgm_metadata_JSONheader(1, &energy_meta);
  

  Lgm_metadata_initvar(0, 1, "COUNTS", &counts_meta);
  Lgm_metadata_addAttr("LABEL", "Counts", &counts_meta);

  json_header = Lgm_metadata_JSONheader(2, &energy_meta, &counts_meta);


  printf("\n\n**************************\n");
  printf("%s", json_header);

  printf("**************************\n\n");


  return (0);
}

