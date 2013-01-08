
#include <stdio.h>

#include <Lgm_Vec.h>


int main(void) {
  Lgm_Vector *v;
  double Lat, Lon, r;
  double Arr[3];

  v = Lgm_CreateVector(1,2,3);
  Lgm_PrintVector(v);

  Lgm_NormalizeVector(v);
  Lgm_PrintVector(v);

  Lgm_ForceMagnitude(v, 2);
  Lgm_PrintVector(v);

  Lgm_SetVecVal(v, 3);
  Lgm_PrintVector(v);

  Lgm_SetVecElements(v, 1,2,3);
  Lgm_PrintVector(v);

  Lgm_CartToSphCoords(v, &Lat, &Lon, &r);
  printf("Lat:%lf Lon:%lf r:%lf\n", Lat, Lon, r);

  Lgm_SphToCartCoords(Lat, Lon, r, v );
  Lgm_PrintVector(v);

  Lgm_VecToArr(v, Arr);
  printf("Arr[0]:%lf Arr[1]:%lf Arr[2]:%lf\n", Arr[0], Arr[1], Arr[2]);

  Arr[0] = -1;
  Lgm_ArrToVec(Arr, v);
  Lgm_PrintVector(v);


  Lgm_FreeVector(v);
  return(0);
}


