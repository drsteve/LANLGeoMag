#ifndef LGM_PRIORITYQUEUE
#define LGM_PRIORITYQUEUE

#include <stdio.h>
#include <stdlib.h>


/**
 *
 * This structure defines a "node" in a "Priority Queue". 
 *
 */
typedef struct Lgm_pQueue_Node {

    double  key;    //<! A number defining the priority level. Lower means higher priority.
    void  *Data;

} Lgm_pQueue_Node;


/**
 *
 * This structure stores info needed to represent a heap-based "Priority Queue". 
 * The HeapArray is dynamically sized via both at creation via Lgm_pQueue_Create() and 
 * potentially when an insert needs a bigger array via Lgm_pQueue_Insert().
 *
 *  The HeapArray stores the heap tree in standard order -- I.e., index 0 is
 *  not used, the root node (highest priority node) is stored at index 1.
 *  Nodes are stored as follows;
 *  
 *      Parent(i) = i/2, LeftChild(i) = 2i, RightChild(i) = 2i+1
 *
 */
typedef struct Lgm_pQueue {

    long int            nHeapArray;     //<! number of elements allocated in the HeapArray
    Lgm_pQueue_Node    *HeapArray;      //<! point to an array of nodes.
    long int            HeapSize;       //<! Number of elements in the Heap

} Lgm_pQueue;


Lgm_pQueue *Lgm_pQueue_Create( long int n );
void        Lgm_pQueue_Destroy( Lgm_pQueue *p );
void        Lgm_pQueue_Insert( Lgm_pQueue_Node *X, Lgm_pQueue *p );
void        Lgm_pQueue_PercolateDown( long int i, Lgm_pQueue *p );
int         Lgm_pQueue_Pop( Lgm_pQueue_Node *X, Lgm_pQueue *p );



#endif
