#ifndef LGM_OCTREE 
#define LGM_OCTREE

#include "Lgm_Vec.h"

#define     OCTREE_MAX_LEVELS           16
#define     OCTREE_ROOT_LEVEL           15
#define     OCTREE_MAX_VAL              (32768.0)
#define     OCTREE_MAX_DATA_PER_OCTANT  10

#define     TRUE    1
#define     FALSE   0


#define     OCTREE_KNN_SUCCESS              1
#define     OCTREE_KNN_TOO_FEW_NNS          0
#define     OCTREE_KNN_NOT_ENOUGH_DATA     -1

#define     OCTREE_IS_NULL                 -2


typedef struct _Lgm_OctreeData {

    Lgm_Vector          Position;
    Lgm_Vector          B;
    double              Dist2;

} Lgm_OctreeData;







typedef struct _Lgm_OctreeCell {

    unsigned int        xLocationCode;  // X Location Code
    unsigned int        yLocationCode;  // Y Location Code
    unsigned int        zLocationCode;  // Z Location Code

    unsigned int        Level;          // keeps track of what level we are on

    Lgm_Vector          Center;         // Center of cube
    double              h;              // half width of cube face

    struct _Lgm_OctreeCell  *Parent;        // points to parent cell
    struct _Lgm_OctreeCell  *Octant;        // points to block of 8 children cells


    unsigned long int   nDataBelow;     // Keeps track of data contained in all children below.
    unsigned long int   nData;          // Number of data items in this cell. (0 for non-leaf nodes)
    Lgm_OctreeData      *Data;          // Cell data. This will be NULL except for 
                                        // the leafs.

} Lgm_OctreeCell;


typedef struct _pQueue {

    Lgm_OctreeCell  *Obj;
    double          MinDist2;   // Minimum possible distance^2 between object and query point.
                                // If the object is a point, then its the actual distance^2.

    int             IsPoint;    // If this is TRUE, then a data point is stored in Obj.Data[j]
                                // Else its a non-leaf node.
    int             j;          // Index where data point is stored.

    struct _pQueue *Prev;       // for linked list
    struct _pQueue *Next;       // for linked list

} pQueue;


void            Binary( unsigned int n, char *Str );
void            Lgm_OctreeFreeBranch( Lgm_OctreeCell *Cell );
void            Lgm_FreeOctree( Lgm_OctreeCell *ot );
Lgm_OctreeCell  *Lgm_CreateOctreeRoot( );
Lgm_OctreeCell  *Lgm_OctreeTraverseToLocCode( Lgm_OctreeCell *Cell, unsigned int ChildLevel, unsigned int xLocationCode, unsigned int yLocationCode, unsigned int zLocationCode );
Lgm_OctreeCell  *Lgm_LocateNearestCell( Lgm_OctreeCell *Root, Lgm_Vector *q );
double          MinDist( Lgm_OctreeCell *Cell, Lgm_Vector *q );
double          InsertCell( Lgm_OctreeCell *Cell, Lgm_Vector *q, pQueue **PQ, double MaxDist2 );
void            InsertPoint( Lgm_OctreeCell *Cell, int j, Lgm_Vector *q, pQueue **PQ );
Lgm_OctreeCell  *DescendTowardClosestLeaf( Lgm_OctreeCell *Node, pQueue **PQ, Lgm_Vector *q, double MaxDist2 );
pQueue          *PopObj( pQueue **PQ );
int             Lgm_Octree_kNN( Lgm_Vector *q, Lgm_OctreeCell *Root, int K, int *Kgot, double MaxDist2, Lgm_OctreeData *kNN );
Lgm_OctreeCell  *CreateNewOctants( Lgm_OctreeCell *Parent );
void            SubDivideVolume( Lgm_OctreeCell *Vol );
Lgm_OctreeCell  *Lgm_InitOctree( Lgm_Vector *ObjectPoints, Lgm_Vector *ObjectData, unsigned long int N, double *Min, double *Max, double *Diff );



    

#endif

/*
 *    $Id: Lgm_Octree.h 46 2010-10-01 20:46:59Z mgh $
 */


