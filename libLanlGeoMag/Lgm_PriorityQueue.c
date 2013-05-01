#include <stdio.h>
#include <stdlib.h>

typedef struct Lgm_pQueue_Node {

    double  key;
    void  *Data;

} Lgm_pQueue_Node;


typedef struct Lgm_pQueue {

    long int            nHeapArray;
    Lgm_pQueue_Node    *HeapArray;

    long int            HeapSize;     // number of items currently on the heap.

} Lgm_pQueue;


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
    int i;
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

void Lgm_pQueue_Pop( Lgm_pQueue_Node *X, Lgm_pQueue *p ) {


    // return the first element
    *X = p->HeapArray[1];

    // move last element to top
    p->HeapArray[1] = p->HeapArray[ (p->HeapSize)-- ];

    // reorder tree
    Lgm_pQueue_PercolateDown( 1, p );
    
    return;

}




int main(){

    int              i;
    Lgm_pQueue      *p;
    Lgm_pQueue_Node  X, Y;

    p = Lgm_pQueue_Create( 1000 );

    X.key = 1.0;    Lgm_pQueue_Insert( &X, p );
    X.key = 25.0;   Lgm_pQueue_Insert( &X, p );
    X.key = 17.0;   Lgm_pQueue_Insert( &X, p );
    X.key = 3.0;    Lgm_pQueue_Insert( &X, p );
    X.key = 36.0;   Lgm_pQueue_Insert( &X, p );
    X.key = 100.0;  Lgm_pQueue_Insert( &X, p );
    X.key = 7.0;    Lgm_pQueue_Insert( &X, p );
    X.key = 0.5;    Lgm_pQueue_Insert( &X, p );
    X.key = 1000.0; Lgm_pQueue_Insert( &X, p );
    X.key = 2.0;    Lgm_pQueue_Insert( &X, p );
    X.key = 19.0;   Lgm_pQueue_Insert( &X, p );

for (i=0; i<=p->HeapSize; i++ ) printf("%g ", p->HeapArray[i].key ); printf("\n");
    

    Lgm_pQueue_Pop( &Y, p ); printf("Y.key = %g\n", Y.key); //for (i=0; i<=p->HeapSize; i++ ) printf("%g ", p->HeapArray[i].key ); printf("\n");
    Lgm_pQueue_Pop( &Y, p ); printf("Y.key = %g\n", Y.key); //for (i=0; i<=p->HeapSize; i++ ) printf("%g ", p->HeapArray[i].key ); printf("\n");
    Lgm_pQueue_Pop( &Y, p ); printf("Y.key = %g\n", Y.key); //for (i=0; i<=p->HeapSize; i++ ) printf("%g ", p->HeapArray[i].key ); printf("\n");
    Lgm_pQueue_Pop( &Y, p ); printf("Y.key = %g\n", Y.key); //for (i=0; i<=p->HeapSize; i++ ) printf("%g ", p->HeapArray[i].key ); printf("\n");
    Lgm_pQueue_Pop( &Y, p ); printf("Y.key = %g\n", Y.key); //for (i=0; i<=p->HeapSize; i++ ) printf("%g ", p->HeapArray[i].key ); printf("\n");
    Lgm_pQueue_Pop( &Y, p ); printf("Y.key = %g\n", Y.key); //for (i=0; i<=p->HeapSize; i++ ) printf("%g ", p->HeapArray[i].key ); printf("\n");
    Lgm_pQueue_Pop( &Y, p ); printf("Y.key = %g\n", Y.key); //for (i=0; i<=p->HeapSize; i++ ) printf("%g ", p->HeapArray[i].key ); printf("\n");
    Lgm_pQueue_Pop( &Y, p ); printf("Y.key = %g\n", Y.key); //for (i=0; i<=p->HeapSize; i++ ) printf("%g ", p->HeapArray[i].key ); printf("\n");
    Lgm_pQueue_Pop( &Y, p ); printf("Y.key = %g\n", Y.key); //for (i=0; i<=p->HeapSize; i++ ) printf("%g ", p->HeapArray[i].key ); printf("\n");
    Lgm_pQueue_Pop( &Y, p ); printf("Y.key = %g\n", Y.key); //for (i=0; i<=p->HeapSize; i++ ) printf("%g ", p->HeapArray[i].key ); printf("\n");
        



    Lgm_pQueue_Destroy( p );

    return(0);
}
