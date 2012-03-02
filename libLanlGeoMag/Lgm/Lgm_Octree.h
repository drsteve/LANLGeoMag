#ifndef LGM_OCTREE
#define LGM_OCTREE

#include "Lgm_Vec.h"

#define     OCTREE_MAX_LEVELS           16
#define     OCTREE_ROOT_LEVEL           15
#define     OCTREE_MAX_VAL              (32768.0)
#define     OCTREE_MAX_DATA_PER_OCTANT  4           // Too high slows things down. Fastest seems to be 4-5 range...

#define     TRUE    1
#define     FALSE   0


#define     OCTREE_KNN_SUCCESS              1
#define     OCTREE_KNN_TOO_FEW_NNS          0
#define     OCTREE_KNN_NOT_ENOUGH_DATA     -1

#define     OCTREE_IS_NULL                 -2


typedef struct _Lgm_OctreeData {

    unsigned long int   Id;
    Lgm_Vector          Position;
    Lgm_Vector          B;
    double              Dist2;

} Lgm_OctreeData;







/**
 *
 * This structure contains all the information needed to represent a node of
 * the Octree.
 *
 */
typedef struct _Lgm_OctreeCell {

    unsigned int        xLocationCode;  //<! X Location Code
    unsigned int        yLocationCode;  //<! Y Location Code
    unsigned int        zLocationCode;  //<! Z Location Code

    unsigned int        Level;          //<! keeps track of what level we are on

    Lgm_Vector          Center;         //<! Center of cube
    double              h;              //<! half width of cube face

    struct _Lgm_OctreeCell  *Parent;    //<! points to parent cell
    struct _Lgm_OctreeCell  *Octant;    //<! points to block of 8 children cells


    unsigned long int   nDataBelow;     //<! Keeps track of data contained in all children below.
    unsigned long int   nData;          //<! Number of data items in this cell. (0 for non-leaf nodes)
    Lgm_OctreeData      *Data;          //<! Cell data. This will be NULL except for the leafs.

} Lgm_OctreeCell;



/**
 *
 * An over-arching structure to contain both the Octree and information about
 * it. For example, when contrsucting the Octree, the positions are all scaled
 * to fit between -1 and 1. The Min and Max values are found and the x, y, z
 * components of the positions are scaled as follows;
 *
 *      \f[ x_new = (x-Min)/(Max-Min) \f]
 *
 * The scaling values are needed to convert back to un-scaled positions.
 *
 */
typedef struct _Lgm_Octree {

    unsigned long int n;            //<! Total number of points in the octree
    double            Min;          //<! Min value for scaling positions
    double            Max;          //<! Max value for scaling positions
    double            Diff;         //<! Max-Min
    long int          kNN_Lookups;  //<! Numbenr of kNN lookups performed

    Lgm_OctreeCell   *Root;         //<! Pointer to the Root node of the Octree

} Lgm_Octree;


/**
 *
 * This structure holds the "Priority Queue" information that is used as part
 * of the kNN (k Nearest Neighbor) search algorithm.
 *
 */
typedef struct _pQueue {

    Lgm_OctreeCell  *Obj;
    double          MinDist2;   //<! Minimum possible distance^2 between object and query point.
                                //<! If the object is a point, then its the actual distance^2.

    int             IsPoint;    //<! If this is TRUE, then a data point is stored in Obj.Data[j] Else its a non-leaf node.
    int             j;          //<! Index where data point is stored.

    struct _pQueue *Prev;       //<! Next node in linked list
    struct _pQueue *Next;       //<! Prev node in linked list

} pQueue;


void            Binary( unsigned int n, char *Str );
void            Lgm_OctreeFreeBranch( Lgm_OctreeCell *Cell );
void            Lgm_FreeOctree( Lgm_Octree *ot );
Lgm_OctreeCell  *Lgm_CreateOctreeRoot( );
Lgm_OctreeCell  *Lgm_OctreeTraverseToLocCode( Lgm_OctreeCell *Cell, unsigned int ChildLevel, unsigned int xLocationCode, unsigned int yLocationCode, unsigned int zLocationCode );
Lgm_OctreeCell  *Lgm_LocateNearestCell( Lgm_OctreeCell *Root, Lgm_Vector *q );
double          MinDist( Lgm_OctreeCell *Cell, Lgm_Vector *q );
double          InsertCell( Lgm_OctreeCell *Cell, Lgm_Vector *q, pQueue **PQ, double MaxDist2 );
void            InsertPoint( Lgm_OctreeCell *Cell, int j, Lgm_Vector *q, pQueue **PQ );
Lgm_OctreeCell  *DescendTowardClosestLeaf( Lgm_OctreeCell *Node, pQueue **PQ, Lgm_Vector *q, double MaxDist2 );
pQueue          *PopObj( pQueue **PQ );
int             Lgm_Octree_kNN( Lgm_Vector *q, Lgm_Octree *Octree, int K, int *Kgot, double MaxDist2, Lgm_OctreeData *kNN );
Lgm_OctreeCell  *CreateNewOctants( Lgm_OctreeCell *Parent );
void            SubDivideVolume( Lgm_OctreeCell *Vol );
Lgm_Octree      *Lgm_InitOctree( Lgm_Vector *ObjectPoints, Lgm_Vector *ObjectData, unsigned long int N );
void            Lgm_OctreeScalePosition( Lgm_Vector *u, Lgm_Vector *v, Lgm_Octree *Octree);
void            Lgm_OctreeUnScalePosition( Lgm_Vector *v, Lgm_Vector *u, Lgm_Octree *Octree);
void            Lgm_OctreeScaleDistance( double u, double *v, Lgm_Octree *Octree);
void            Lgm_OctreeUnScaleDistance( double v, double *u, Lgm_Octree *Octree);





#endif
