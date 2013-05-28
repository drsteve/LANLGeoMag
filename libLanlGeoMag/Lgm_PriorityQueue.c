#pragma GCC push_options
#pragma GCC optimize ("O3")

#include "Lgm/Lgm_PriorityQueue.h"


/**
 *  Create a Heap-based priority queue.
 */
Lgm_pQueue *Lgm_pQueue_Create( long int n ) {

    Lgm_pQueue *p = (Lgm_pQueue *)calloc( 1, sizeof(Lgm_pQueue) );

    /*
     * Allocate memory for the heap array. Here we create an array of pointers
     * to void so that they can point to anything.
     */
    p->nHeapArray = n;
    p->HeapArray  = (Lgm_pQueue_Node *)calloc( n, sizeof(Lgm_pQueue_Node) );

    /*
     * Init heap size to 0.
     */
    p->HeapSize = 0;

    return( p );

}


/**
 *  Free a priority queue.
 */
void Lgm_pQueue_Destroy( Lgm_pQueue *p ) {
    free( p->HeapArray );
    free( p );
    return;
}


void Lgm_pQueue_Insert( Lgm_pQueue_Node *X, Lgm_pQueue *p ) {
    int i, n;

    if ( p->HeapSize >= p->nHeapArray-1 ) {
        n = 2*p->nHeapArray;
        p->HeapArray  = (Lgm_pQueue_Node *)realloc( p->HeapArray, n*sizeof(Lgm_pQueue_Node) );
        if ( !p->HeapArray ) {
            printf("Lgm_pQueue_Insert: Memory allocation failure. Trying to realloc %d elements\n");
            exit(1);
        } else {
            p->nHeapArray = n;
            //printf("Realloced HeapArray. Size = %d\n", p->nHeapArray);
        }
        
    }

    for ( i = ++(p->HeapSize); (i > 1)&&( X->key < p->HeapArray[i/2].key); i /= 2 ) {
        p->HeapArray[i] = p->HeapArray[i/2];
    }
    p->HeapArray[i] = *X;

    return;
}


void Lgm_pQueue_PercolateDown( long int i, Lgm_pQueue *p ) {
    int             C;
    Lgm_pQueue_Node t;

    t = p->HeapArray[i];
    for ( ; 2*i <= p->HeapSize; i = C ) {
        C = 2*i;
        if ( (C != p->HeapSize) && ( p->HeapArray[C+1].key  < p->HeapArray[C].key ) )  C++;
        if ( p->HeapArray[C].key < t.key ) {
            p->HeapArray[i] = p->HeapArray[C];
        } else {
            break;
        }
    }
    p->HeapArray[i] = t;
    return;
}

int Lgm_pQueue_Pop( Lgm_pQueue_Node *X, Lgm_pQueue *p ) {


    if ( p->HeapSize <= 0 ) return(0);

    // return the first element
    *X = p->HeapArray[1];

    // move last element to top
    p->HeapArray[1] = p->HeapArray[ (p->HeapSize)-- ];

    // reorder tree
    Lgm_pQueue_PercolateDown( 1, p );
    
    return(1);

}


int Lgm_pQueue_Peek( Lgm_pQueue_Node *X, Lgm_pQueue *p ) {


    if ( p->HeapSize <= 0 ) return(0);

    // return the first element
    *X = p->HeapArray[1];

    return(1);

}


#pragma GCC pop_options
