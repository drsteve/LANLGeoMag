
#include <stdio.h>

#include <Lgm_Vec.h>


int main(void) {
  Lgm_Vector *v, *v2, *v3, *v4;
  double Lat, Lon, r, ang;
  double Arr[3];

  printf("Make a Vector:\n");
  v = Lgm_CreateVector(1,2,3);
  Lgm_PrintVector(v);

  printf("Normalize it:\n");
  Lgm_NormalizeVector(v);
  Lgm_PrintVector(v);

  printf("Force it to have a magnitude:\n");
  Lgm_ForceMagnitude(v, 2);
  Lgm_PrintVector(v);

  printf("Set its value:\n");
  Lgm_SetVecVal(v, 3);
  Lgm_PrintVector(v);

  printf("Set its elements:\n");
  Lgm_SetVecElements(v, 1,2,3);
  Lgm_PrintVector(v);

  printf("Change it to spherical:\n");
  Lgm_CartToSphCoords(v, &Lat, &Lon, &r);
  printf("Lat:%lf Lon:%lf r:%lf\n", Lat, Lon, r);

  printf("And back to Cart:\n");
  Lgm_SphToCartCoords(Lat, Lon, r, v );
  Lgm_PrintVector(v);

  printf("Make it an array:\n");
  Lgm_VecToArr(v, Arr);
  printf("Arr[0]:%lf Arr[1]:%lf Arr[2]:%lf\n", Arr[0], Arr[1], Arr[2]);

  printf("Edit the array and make it a vector:\n");
  Arr[0] = -1;
  Lgm_ArrToVec(Arr, v);
  Lgm_PrintVector(v);

  printf("Divide it by a scalar:\n");
  Lgm_VecDivideScalar(v, 2);
  Lgm_PrintVector(v);

  printf("Multiply by a scalar:\n");
  Lgm_VecMultiplyScalar(v, 2);
  Lgm_PrintVector(v);

  printf("Make a new vector and get the angle between them:\n");
  v2 = Lgm_CreateVector(1,0,0);
  Lgm_SetVecElements(v, 0, 1, 0);
  ang = Lgm_VectorAngle(v, v2);
  printf("Angle:%lf\n", ang);

  printf("Make a new vector and get the angle between them:\n");
  Lgm_SetVecElements(v, 0, 1, 0);
  Lgm_SetVecElements(v2, 0, 1, 0);
  ang = Lgm_VectorAngle(v, v2);
  printf("Angle:%lf\n", ang);

  Lgm_SetVecElements(v, -2, 3, 1);
  Lgm_SetVecElements(v2, 0, 4, 0);
  v3 = Lgm_CreateVector(-1, 3, 3);

  r = Lgm_ScalarTripleProduct(v, v2, v3);
  printf("The ScalarTripleProduct of\n");
  Lgm_PrintVector(v);
  Lgm_PrintVector(v2);
  Lgm_PrintVector(v3);
  printf("  is %lf\n", r);

  Lgm_SetVecElements(v,  1, 2, 3);
  Lgm_SetVecElements(v2, 1, -3, 4);
  Lgm_SetVecElements(v3, 2, 4, -2);

  v4 = Lgm_CreateVector(0,0,0);

  Lgm_VectorTripleProduct(v, v2, v3, v4);
  printf("The VectorTripleProduct of\n");
  Lgm_PrintVector(v);
  Lgm_PrintVector(v2);
  Lgm_PrintVector(v3);
  printf("  is \n");
  Lgm_PrintVector(v4);



  Lgm_FreeVector(v);
  Lgm_FreeVector(v2);
  Lgm_FreeVector(v3);
  Lgm_FreeVector(v4);
  return(0);
}


