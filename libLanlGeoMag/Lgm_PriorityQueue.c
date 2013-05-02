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
            printf("Lgm_pQueue_Insert: Memory allocation failure. Tryingm to realloc %d elements\n");
            exit(1);
        } else {
            p->nHeapArray = n;
printf("realloced HeapArray. Size = %d\n", p->nHeapArray);
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




int main(){

    int              i;
    Lgm_pQueue      *p;
    Lgm_pQueue_Node  X, Y;
    double          *v, a[3], b[3], c[5], d[2];

    a[0] = 1.0;
    a[1] = 2.0;
    a[2] = 3.0;

    b[0] = 4.0;
    b[1] = 5.0;
    b[2] = 6.0;

    c[0] = 7.0;
    c[1] = 8.0;
    c[2] = 9.0;
    c[1] = 10.0;
    c[2] = 11.0;

    d[0] = 12.0;
    d[1] = 13.0;

    p = Lgm_pQueue_Create( 3 );

    X.key = 1.0;    X.Data = (void *) a; Lgm_pQueue_Insert( &X, p );
    X.key = 25.0;   X.Data = (void *) b; Lgm_pQueue_Insert( &X, p );
    X.key = 17.0;   X.Data = (void *) c; Lgm_pQueue_Insert( &X, p );
    X.key = 16.2;   X.Data = (void *) d; Lgm_pQueue_Insert( &X, p );


for (i=0; i<=p->HeapSize; i++ ) printf("%g ", p->HeapArray[i].key ); printf("\n");
    

for (i=0; i<2000; ++i){
    if ( Lgm_pQueue_Pop( &Y, p ) )  {
        v = (double *)Y.Data;
        printf("Y.key = %g       v[1] = %g\n", Y.key, v[1]); //for (i=0; i<=p->HeapSize; i++ ) printf("%g ", p->HeapArray[i].key ); printf("\n");
    }
}
        



    Lgm_pQueue_Destroy( p );

    return(0);
}
