
#include <Lgm_JsonHeader.h>
#include <Lgm_Utils.h>
#include <stdio.h>

/* typedef struct Lgm_json_variable { */
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
  Lgm_json_variable energy_var;
  Lgm_json_variable counts_var;
  char* h_energy; 
  char* h_counts;

  // populate the energy and counts
  Lgm_LinSpace(1, 100, 10, energy);
  Lgm_LogSpace(100, 1000, 10, counts);

  Lgm_json_initvar(0, 1, "ENERGY", &energy_var);
  printf("%s\n", energy_var.name);
  Lgm_json_addAttr("UNITS", "MeV", &energy_var);
  Lgm_json_addAttr("LABEL", "Energy [MeV]", &energy_var);
  printf("%d\n", (int)energy_var.n_attrs);
  printf("%s\n", energy_var.attributes[0].name);
  printf("%s\n", energy_var.attributes[1].name);
  h_energy = Lgm_json_toString(&energy_var); 
  printf("\n\n%s\n\n", h_energy);


  Lgm_json_initvar(0, 1, "COUNTS", &counts_var);
  printf("%s\n", counts_var.name);
  Lgm_json_addAttr("LABEL", "Counts", &counts_var);
  printf("%d\n", (int)counts_var.n_attrs);
  printf("%s\n", counts_var.attributes[0].name);
  h_counts = Lgm_json_toString(&counts_var); 
  printf("\n\n%s\n\n", h_counts);




  return (0);
}

