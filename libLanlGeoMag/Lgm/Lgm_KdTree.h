#ifndef LGM_KDTREE
#define LGM_KDTREE

#include "Lgm_PriorityQueue.h"

#define     KDTREE_MAX_LEVEL            1000
#define     KDTREE_ROOT_LEVEL           0
#define     KDTREE_MAX_DATA_PER_NODE    1           // Algorithm now depends on there being no more than 1 point per leaf....

#define     TRUE    1
#define     FALSE   0


#define     KDTREE_KNN_SUCCESS              1
#define     KDTREE_KNN_TOO_FEW_NNS          0
#define     KDTREE_KNN_NOT_ENOUGH_DATA     -1

#define     KDTREE_IS_NULL                 -2


#define     LGM_KDTREE_SPLIT_SEQUENTIAL     0
#define     LGM_KDTREE_SPLIT_RANDOM         1
#define     LGM_KDTREE_SPLIT_MAXRANGE       2



typedef struct _Lgm_KdTreeData {

    unsigned long int    Id;
    int                  D;         // The dimension of the space
    double              *Position;  // D-dimension position vector
    void                *Object;    // Pointer to Object associated with this point
    double              Dist2;      

} Lgm_KdTreeData;







/**
 *
 * This structure contains all the information needed to represent a node of
 * the KdTree.
 *
 */
typedef struct _Lgm_KdTreeNode {

    unsigned int         Level;         //<! keeps track of what level we are on
    unsigned int         D;             //<! keeps track of number of dimensions
    unsigned int         d;             //<! keeps track of the splitting dimension. (-1 for root node).
    double              *Min;           //<! Array of Minima - min value for each dimension for the data stored in or below this node.
    double              *Max;           //<! Array of Maxima - max value for each dimension for the data stored in or below this node.
    double              *Diff;          //<! Array of Ranges (Max-Min) for each dimension for the data stored in or below this node.

    struct _Lgm_KdTreeNode  *Parent;    //<! points to parent node
    struct _Lgm_KdTreeNode  *Left;      //<! points to left sub-tree node
    struct _Lgm_KdTreeNode  *Right;     //<! points to right sub-tree node


    unsigned long int   nDataBelow;     //<! Keeps track of data contained in all children below.
    unsigned long int   nData;          //<! Number of data items in this node. (0 for non-leaf nodes)
    Lgm_KdTreeData      *Data;          //<! Node data. This will be NULL except for the leafs.

} Lgm_KdTreeNode;



/**
 *
 * An over-arching structure to contain both the KdTree and information about
 * it. For example, when contrsucting the KdTree, the positions are all scaled
 * to fit between -1 and 1. The Min and Max values are found and the x, y, z
 * components of the positions are scaled as follows;
 *
 *      \f[ x_new = (x-Min)/(Max-Min) \f]
 *
 * The scaling values are needed to convert back to un-scaled positions.
 *
 */
typedef struct _Lgm_KdTree {

    unsigned long int n;             //<! Total number of points in the octree
    double            Min;           //<! Min value for scaling positions
    double            Max;           //<! Max value for scaling positions
    double            Diff;          //<! Max-Min
    long int          kNN_Lookups;   //<! Numbenr of kNN lookups performed
    int               SplitStrategy; //<! Strategy for doing dimension splitting. (Can be one of LGM_KDTREE_SPLIT_SEQUENTIAL, LGM_KDTREE_SPLIT_RANDOM, LGM_KDTREE_SPLIT_MAXRANGE)
    Lgm_pQueue       *PQ;            //<! Heap-based priority queue.

    Lgm_KdTreeNode   *Root;          //<! Pointer to the Root node of the KdTree

} Lgm_KdTree;


/**
 *
 * This structure holds the "Priority Queue" node information that is used as part
 * of the kNN (k Nearest Neighbor) search algorithm. 
 *
 */
typedef struct _Lgm_KdTree_pQueue_Node {

    Lgm_KdTreeNode  *Obj;
    double          MinDist2;           //<! Minimum possible distance^2 between object and query point.
                                        //<! If the object is a point, then its the actual distance^2.

    int             IsPoint;            //<! If this is TRUE, then a data point is stored in Obj.Data[j] Else its a non-leaf node.
    int             j;                  //<! Index where data point is stored.

//    struct _Lgm_KdTree_pQueue *Prev;    //<! Next node in linked list
//    struct _Lgm_KdTree_pQueue *Next;    //<! Prev node in linked list

} Lgm_KdTree_pQueue_Node;



void                Lgm_KdTree_SubDivideVolume( Lgm_KdTreeNode *t, Lgm_KdTree *kt );
Lgm_KdTree         *Lgm_KdTree_Init( double **Positions, void **Objects, unsigned long int N, int D ) ;
Lgm_KdTreeNode     *Lgm_CreateKdTreeRoot( int D );
int                 Lgm_KdTree_kNN( double *q_in, int D, Lgm_KdTree *KdTree, int K, int *Kgot, double MaxDist2, Lgm_KdTreeData *kNN );
double              Lgm_KdTree_MinDist( Lgm_KdTreeNode *Node, double *q );
double              Lgm_KdTree_InsertNode( Lgm_KdTreeNode *Node, double *q, Lgm_pQueue *PQ, double MaxDist2 );
void                Lgm_KdTree_InsertPoint( Lgm_KdTreeNode *Node, int j, double *q, Lgm_pQueue *PQ );
//Lgm_KdTree_pQueue  *Lgm_KdTree_PopObj( Lgm_KdTree_pQueue **PQ );
void                Lgm_KdTree_DescendTowardClosestLeaf( Lgm_KdTreeNode *Node, Lgm_pQueue *PQ, double *q, double MaxDist2 );
void                Lgm_KdTree_PrintPQ( Lgm_pQueue *PQ );


#endif
