
#include <stdio.h>
#include <Lgm_Utils.h>

int main(void) {
  double array[10];
  int i;
  
  Lgm_LogSpace(1, 100, 10, array);

  printf("Lgm_LogSpace:\n");
  for (i=0;i<10;i++)
    printf("%f\n", array[i]);
  printf("\n");


  return (0);
}


