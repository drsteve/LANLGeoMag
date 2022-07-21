#ifndef LGM_OBJECTS_H
#define LGM_OBJECTS_H
#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif
#include "Lgm_Vec.h"

/* 
 *  structure to hold a sphere's parameters 
 */
typedef struct SphereType {
        Lgm_Vector      Origin;             /* origin of sphere */
        double          Radius;             /* radius of sphere */
        double          Radius2;            /* radius^2 of sphere */
        int             Precision;          /* Sphere precision (controls how many facets to draw) */
        int             nVertex;            /* Number of vertices */
        Lgm_Vector      *Vertex;            /* vertices. Need to alloc as an array */
        Lgm_Vector      *Color;             /* color. Need to alloc as an array */
} SphereType;


/* 
 *  structure to hold a ellipsoid's parameters 
 *  Do we need all the extra stuff?
 */
typedef struct EllipsoidType {
        Lgm_Vector      Origin;             /* origin of ellipsoid */
        double          Radius_a;           /* equatorial radius of ellipsoid */
        double          Radius_b;           /* polar radius of ellipsoid */
        double          Radius2_a;          /* radius_a^2 of ellipsoid */
        double          Radius2_b;          /* radius_b^2 of ellipsoid */
        int             Precision;          /* Ellipsoid precision (controls how many facets to draw) */
        int             nVertex;            /* Number of vertices */
        Lgm_Vector      *Vertex;            /* vertices. Need to alloc as an array */
        Lgm_Vector      *Color;             /* color. Need to alloc as an array */
} EllipsoidType;


/* 
 *   structure to hold a plane's parameters 
 */
typedef struct PlaneType {
        Lgm_Vector          Normal;         /* normal to plane */
        double          d;              /* RHS of plane equation */
} PlaneType;




/* 
 *  structure to hold a box's parameters 
 */
typedef struct BoxType {
	Lgm_Vector		Bounds[2];		/* Corner points of box */
	double		xmin, xmax;
	double		ymin, ymax;
	double		zmin, zmax;
} BoxType;


/*
 *  Define Ray structure.
 *  In this representation, a point that is t units along the ray path
 *  (from the Origin) is given simply by;
 *             
 *              r.x = t * Ray.Direction.x
 *              r.y = t * Ray.Direction.y
 *              r.z = t * Ray.Direction.z
 *
 *  where r is a vector. 
 * 
 */
typedef struct RayType {


    Lgm_Vector 	Origin;			/*  GSM position (in Re) Where ray starts. 
					 *  I.e. the spacecraft location.           
					 */

    Lgm_Vector 	Direction;		/*  GSM look-direction of ray. 
					 *  Must be a unit vector.  
					 */

    double 	StartDistance;		/*  Distance (in Re) along ray to start of 
					 *  integration region. Always positive.
					 *  Since simulation volume will generally
					 *  be a sphere, this value will be obtained
					 *  via a ray/sphere intersection test.
					 */

    double 	EndDistance;		/*  Distance (in Re) along ray to end of 
					 *  integration region. Always positive.
					 *  Since simulation volume will generally
					 *  be a sphere, this value will be obtained
					 *  via a ray/sphere intersection test.
				         */

    Lgm_Vector	InvDirection;  /* 
					            *  InvDirection.x = 1.0/Direction.x;
					            *  InvDirection.y = 1.0/Direction.y;
					            *  InvDirection.z = 1.0/Direction.z;
					            */

    int		Sign[3];    /* 
					     *  Sign[0] = ( InvDirection.x < 0.0 );
					     *  Sign[1] = ( InvDirection.y < 0.0 );
					     *  Sign[2] = ( InvDirection.z < 0.0 );
					     */

    Lgm_Vector 	P[10];	/*
					     *  Convenient place to store intersection Points
					     */

    double	u[10];
    double	v[10];
    double 	t[10];		/*
					     *  t parameters for intersection points
					     *  u and v are barycentric coords of triangle at intersection point
					     */

    int 	n;			/*
					     *  number of valid intersection points stored
					     */

    double	Earth_t;    /* distance to earth intersection point */
    int flag;


} RayType;

//int Lgm_BoxIntersect( BoxType *Box, RayType *Ray, double t0, double t1 );
//int Lgm_SphereIntersect( SphereType *Sphere, RayType *Ray, double *tmin, double *tmax, double *t );
int Lgm_EllipsoidIntersect( EllipsoidType *Ellipsoid, RayType *Ray, double *tmin, double *tmax, double *t );
double Lgm_EllipsoidVolume( EllipsoidType *Ellipsoid);
EllipsoidType *Lgm_CreateEllipsoid( Lgm_Vector *Origin, double a, double b);
void Lgm_FreeEllipsoid( EllipsoidType *Ellipsoid );

#endif

