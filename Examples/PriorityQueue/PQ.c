#include "Lgm/Lgm_PriorityQueue.h"
int main(){

    int              i;
    Lgm_pQueue      *p;
    Lgm_pQueue_Node  X, Y;
    double          *v, a[3], b[3], c[5], d[2];

    // define some, arbitrary data arrays
    a[0] = 1.0; a[1] = 2.0; a[2] = 3.0;
    b[0] = 4.0; b[1] = 5.0; b[2] = 6.0;
    c[0] = 7.0; c[1] = 8.0; c[2] = 9.0; c[1] = 10.0; c[2] = 11.0;
    d[0] = 12.0; d[1] = 13.0;

    // create a PQ with 3 elements in the heap initially
    p = Lgm_pQueue_Create( 3 );

    // Add 4 items. (An automatic realloc should happen).
    X.key = 1.0;    X.Data = (void *)a; Lgm_pQueue_Insert( &X, p );
    X.key = 25.0;   X.Data = (void *)b; Lgm_pQueue_Insert( &X, p );
    X.key = 17.0;   X.Data = (void *)c; Lgm_pQueue_Insert( &X, p );
    X.key = 16.2;   X.Data = (void *)d; Lgm_pQueue_Insert( &X, p );


    // print the contents of the HeapArray keys.
    for (i=0; i<=p->HeapSize; i++ ) printf("%g ", p->HeapArray[i].key ); printf("\n");


    // test popping things off (20 is alot more than we inserted -- the extras test proper returns when nothing is there...)
    for (i=0; i<20; ++i){
        if ( Lgm_pQueue_Pop( &Y, p ) )  {
            v = (double *)Y.Data;
            printf("Y.key = %g       v[1] = %g\n", Y.key, v[1]); //for (i=0; i<=p->HeapSize; i++ ) printf("%g ", p->HeapArray[i].key ); printf("\n");
        }
    }   
         
    
    // destroy PQ
    Lgm_pQueue_Destroy( p ); 
    
    return(0);

}
