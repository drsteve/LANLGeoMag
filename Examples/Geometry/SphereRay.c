


#include <stdio.h>

#include <Lgm_Objects.h>
#include <Lgm_Vec.h>




int main(void) {

  int intersect=0;

  double tmin=0, tmax=1000, t;

  Lgm_Vector *vec1, *vec2, *vec3;
  
  SphereType sphere;
  RayType ray;
  // make a sphere at the origin
  vec1 = Lgm_CreateVector(0,0,0);
  sphere.Origin = *vec1;
  sphere.Radius = 2.0;
  // make a ray pointing at the origin
  vec2 = Lgm_CreateVector(10,0,0);
  vec3 = Lgm_CreateVector(-1,0,0);
  ray.Origin = *vec2;
  ray.Direction = *vec3;
   

  intersect = Lgm_SphereIntersect(&sphere, &ray, &tmin, &tmax, &t);

    
  Lgm_FreeVector(vec1);
  Lgm_FreeVector(vec2);
  Lgm_FreeVector(vec3);
  

  return(0);
}



// int Lgm_SphereIntersect( SphereType *Sphere, RayType *Ray, double *tmin, double *tmax, double *t );
