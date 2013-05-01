
#include <stdio.h>
#include <Lgm_Utils.h>

int main(void) {
  double array[10];
  int i;
  
  Lgm_LogSpace(1, 100, 10, array);

  printf("Lgm_LogSpace:\n");
  for (i=0;i<10;i++)
    printf("%f\n", array[i]);
  {
    double min, max;
    Lgm_MinMax(array, 10, &min, &max);
    printf("Min: %lf\n", min);
    printf("Max: %lf\n", max);
  }
  printf("\n");



  Lgm_LinSpace(1, 100, 10, array);

  printf("Lgm_LinSpace:\n");
  for (i=0;i<10;i++)
    printf("%f\n", array[i]);
  {
    double min, max;
    Lgm_MinMax(array, 10, &min, &max);
    printf("Min: %lf\n", min);
    printf("Max: %lf\n", max);
  }
  printf("\n");



  return (0);
}


