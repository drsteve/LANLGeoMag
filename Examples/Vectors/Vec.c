
#include <stdio.h>

#include <Lgm_Vec.h>


int main(void) {
  Lgm_Vector *v;
  v = Lgm_CreateVector(1,2,3);
  printf("x: %lf y:%lf z:%lf mag:%lf\n", v->x, v->y, v->z, Lgm_Magnitude(v));

  Lgm_NormalizeVector(v);
  printf("x: %lf y:%lf z:%lf mag:%lf\n", v->x, v->y, v->z, Lgm_Magnitude(v));

  Lgm_ForceMagnitude(v, 2);
  printf("x: %lf y:%lf z:%lf mag:%lf\n", v->x, v->y, v->z, Lgm_Magnitude(v));

  Lgm_SetVecVal(v, 3);
  printf("x: %lf y:%lf z:%lf mag:%lf\n", v->x, v->y, v->z, Lgm_Magnitude(v));

  Lgm_SetVecElements(v, 1,2,3);
  printf("x: %lf y:%lf z:%lf mag:%lf\n", v->x, v->y, v->z, Lgm_Magnitude(v));

  Lgm_FreeVector(v);
  return(0);
}


