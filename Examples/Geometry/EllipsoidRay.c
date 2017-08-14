


#include <stdio.h>

#include <Lgm_Objects.h>
#include <Lgm_Vec.h>


int main(void) {

    int intersect = 0;

    double tmin = 0, tmax = 1000, t, volume;

    Lgm_Vector *vec1, *vec2, *vec3;

    EllipsoidType ellipsoid;
    RayType ray;
    // make a sphere at the origin
    vec1 = Lgm_CreateVector(0, 0, 0);
    ellipsoid.Origin = *vec1;
    ellipsoid.Radius_a = 2.0;
    ellipsoid.Radius_b = 2.0;
    ellipsoid.Radius2_a = ellipsoid.Radius_a * ellipsoid.Radius_a;
    ellipsoid.Radius2_b = ellipsoid.Radius_b * ellipsoid.Radius_b;
    // make a ray pointing at the origin
    vec2 = Lgm_CreateVector(10, 0, 0);
    vec3 = Lgm_CreateVector(-1, 0, 0);
    ray.Origin = *vec2;
    ray.Direction = *vec3;

    // print out some info on the geometry
    printf("Ray origin: %lf  %lf  %lf\n", ray.Origin.x, ray.Origin.y, ray.Origin.z);
    printf("Ray direction: %lf  %lf  %lf\n", ray.Direction.x, ray.Direction.y, ray.Direction.z);


    intersect = Lgm_EllipsoidIntersect(&ellipsoid, &ray, &tmin, &tmax, &t);

    if (intersect) {
        printf("The ray did intersect the ellipsoid, %lf  from the ray origin\n", t);
        printf("Closest intersection point, %lf  %lf  %lf\n",
               ray.Origin.x + ray.Direction.x * t,
               ray.Origin.y + ray.Direction.y * t,
               ray.Origin.z + ray.Direction.z * t);
    } else
        printf("The ray did not intersect the ellipsoid, %lf\n", t);


    volume = Lgm_EllipsoidVolume(&ellipsoid);

    printf("\n\nThe Ellipsoid has an area of %lf\n", volume);

    Lgm_FreeVector(vec1);
    Lgm_FreeVector(vec2);
    Lgm_FreeVector(vec3);


    return (0);
}



// int Lgm_SphereIntersect( SphereType *Sphere, RayType *Ray, double *tmin, double *tmax, double *t );
