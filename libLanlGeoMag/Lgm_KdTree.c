#include "Lgm/Lgm_KdTree.h"
#include "Lgm/quicksort.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/timeb.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>


/**
 *   Store given N-dimensional data into a D-dimensional KD-tree data structure.
 *
 *   Given arrays of positions and data, this routine recursively partitions
 *   the data into a kdtree data structure. 
 *
 *      \param[in]      Points      An array of position vectors in D-dimensional space. ObjectPoints[d][n] is the dth component of the nth point.
 *      \param[in]      Objects     An array of objects in D-dimensional space. Objects[n] is the nth pointer to an object.
 *      \param[in]      N           Number of points.
 *      \param[in]      D           Number of dimensions.
 *
 *      \returns        returns a pointer to the a KdTree structure. User is
 *                      responsible to freeing this with Lgm_FreeKdTree( )
 *
 *      \author         Mike Henderson
 *      \date           2013
 *
 */
Lgm_KdTree *Lgm_KdTree_Init( double **Positions, void **Objects, unsigned long int N, int D ) {


    int                 d;
    unsigned long int   j;
    Lgm_KdTreeNode      *t;
    Lgm_KdTree          *kt;

    /*
     * Initially dump all of the data into a single node.
     */
    t = Lgm_CreateKdTreeRoot( D );

    t->nData      = N;
    t->nDataBelow = t->nData;
    t->Data       = (Lgm_KdTreeData *) calloc( t->nData, sizeof( Lgm_KdTreeData ) );


    /*
     *  Loop over number of data points. Also keep track of ranges in each dimension.
     */
    for (j=0; j<t->nData; j++){

        t->Data[j].Position = (double *) calloc( D, sizeof(double) );

        t->Data[j].Id = j;
        for (d=0; d<D; d++) {
            if (Positions[d][j] < t->Min[d] ) t->Min[d] = Positions[d][j];
            if (Positions[d][j] > t->Max[d] ) t->Max[d] = Positions[d][j];
            t->Data[j].Position[d] = Positions[d][j];
        }
        //t->Data[j].Object = Objects[j];

    }
    for (d=0; d<D; d++) t->Diff[d] = t->Max[d] - t->Min[d];




    kt = (Lgm_KdTree *) calloc( 1, sizeof( Lgm_KdTree) );

    kt->kNN_Lookups   = 0;
    kt->SplitStrategy = LGM_KDTREE_SPLIT_MAXRANGE;
    
    Lgm_KdTree_SubDivideVolume( t, kt );
    kt->Root = t;

    return( kt );
    

}



/**
 *  Create a root-level node for an kdtree.
 *
 *
 *      \returns        void
 *
 *      \author         Mike Henderson
 *      \date           2013
 *
 */
Lgm_KdTreeNode *Lgm_CreateKdTreeRoot( int D ) {

    int              d;
    Lgm_KdTreeNode  *Node;

    /*
     * allocate
     */
    Node = (Lgm_KdTreeNode *) calloc( 1, sizeof( Lgm_KdTreeNode ) );
    Node->nData  = 0;
    Node->D      = D;
    Node->d      = -1; // no split dimension yet
    Node->Parent = NULL;
    Node->Level  = KDTREE_ROOT_LEVEL;

    Node->Min    = (double *) calloc( D, sizeof( double ) );
    Node->Max    = (double *) calloc( D, sizeof( double ) );
    Node->Diff   = (double *) calloc( D, sizeof( double ) );

    for ( d=0; d<D; d++ ) {
        Node->Min[d]  =  9e99;
        Node->Max[d]  = -9e99;
        Node->Diff[d] =  9e99;
    }

    return( Node );

}



/**
 *  Recursively subdivide volume. Stop when the number of data points per node
 *  is lower than the threshold given by KDTREE_MAX_DATA_PER_NODE.
 *
 *    \param[in]     Node in the kdtree to subdivide.
 *
 *    \returns       void
 *
 *    \author        Mike Henderson
 *    \date          2013
 *
 */
void Lgm_KdTree_SubDivideVolume( Lgm_KdTreeNode *t, Lgm_KdTree *kt ) {

    int                 Level, q, d, D, NewLevel;
    unsigned long int   j, nLeft, nRight;
    double              *Pos, V, c, Diff, MaxDiff;
    unsigned long int   *Idx;


    /*
     *  get dimension and set new Level
     */
    D        = t->D;
    Level    = t->Level;
    NewLevel = Level + 1;


    /*
     * Determine dimension to split on.
     * Could do;
     */
    if        ( kt->SplitStrategy == LGM_KDTREE_SPLIT_SEQUENTIAL ) {

        /*
         * Split sequential on each dimesion
         */
        q = Level%D;

    } else if ( kt->SplitStrategy == LGM_KDTREE_SPLIT_RANDOM ) {

        /*
         * Split randomly on each dimesion
         */
        q = (int)( rand()/(double)RAND_MAX*D );

    } else if ( kt->SplitStrategy == LGM_KDTREE_SPLIT_MAXRANGE ) {

        /*
         * Split on dimesion with largest range of values
         */
        MaxDiff = -1.0; q = 0;
        for (d=0; d<D; d++){
            if ( t->Diff[d] > MaxDiff ) {
                MaxDiff = t->Diff[d];
                q = d;
            }
        }
        //printf("Level = %d Split on %d dimension. MaxDiff = %g\n", Level, q, MaxDiff);

    } else {

        printf(" Unknown value of SplitStrategy, got %d (?)\n", kt->SplitStrategy );

    }



    /*
     * Sort the data in the q-dimension.
     */
    Pos = (double *) calloc( t->nData, sizeof( double ) );
    Idx = (unsigned long int *) calloc( t->nData, sizeof( unsigned long int ) );
    for ( j=0; j<t->nData; j++ ) {
        Idx[j] = j;
        Pos[j] = t->Data[j].Position[q];
    }
    quicksort2uli( t->nData, Pos-1, Idx-1 );


    /*
     *  Create new left and right nodes.
     */
    t->Left  = (Lgm_KdTreeNode *) calloc( 1, sizeof( Lgm_KdTreeNode ) );
    t->Right = (Lgm_KdTreeNode *) calloc( 1, sizeof( Lgm_KdTreeNode ) );



    /*
     * Assign first half of the sorted data to the left, and the rest to the right node.
     */
    nRight = t->nData/2;
    nLeft  = t->nData - nRight;
    //V = 0.5*(Pos[nLeft] + Pos[nRight]);
    //V = Pos[nLeft];
    //printf("Split at: %g %g %g\n", Pos[nLeft], Pos[nRight], V );

    t->Left->Data       = (Lgm_KdTreeData *) calloc( nLeft, sizeof( Lgm_KdTreeData ) );
    t->Left->Min        = (double *) calloc( D, sizeof( double ) );
    t->Left->Max        = (double *) calloc( D, sizeof( double ) );
    t->Left->Diff       = (double *) calloc( D, sizeof( double ) );
    t->Left->nData      = nLeft;
    t->Left->nDataBelow = nLeft;
    t->Left->Level      = NewLevel;
    t->Left->Left       = NULL;
    t->Left->Right      = NULL;
    t->Left->D          = D;

/*
    for ( d=0; d<D; d++ ){
        t->Left->Min[d]  = t->Min[d];
        t->Left->Max[d]  = t->Max[d];
        t->Left->Diff[d] = t->Diff[d];
    }
    t->Left->Min[q]  = Pos[0];
    t->Left->Max[q]  = Pos[nLeft-1];
    t->Left->Diff[q] = t->Left->Max[q] - t->Left->Min[q];
*/
    for (d=0; d<D; d++) {
        t->Left->Min[d] =  9e99;
        t->Left->Max[d] = -9e99;
    }
    for (j=0; j<nLeft; j++){
        t->Left->Data[j].Position = (double *) calloc( D, sizeof(double) );
        t->Left->Data[j].Id = t->Data[ Idx[j] ].Id;
        for (d=0; d<D; d++) {
            c = t->Data[ Idx[j] ].Position[d];
            if ( c < t->Left->Min[d] ) t->Left->Min[d] = c;
            if ( c > t->Left->Max[d] ) t->Left->Max[d] = c;
            t->Left->Data[j].Position[d] = c;
        }
        t->Left->Data[j].Object = t->Data[ Idx[j] ].Object;
    }
    for (d=0; d<D; d++) t->Left->Diff[d] = t->Left->Max[d] - t->Left->Min[d];









    t->Right->Data       = (Lgm_KdTreeData *) calloc( nRight, sizeof( Lgm_KdTreeData ) );
    t->Right->Min        = (double *) calloc( D, sizeof( double ) );
    t->Right->Max        = (double *) calloc( D, sizeof( double ) );
    t->Right->Diff       = (double *) calloc( D, sizeof( double ) );
    t->Right->nData      = nRight;
    t->Right->nDataBelow = nRight;
    t->Right->Level      = NewLevel;
    t->Right->Left       = NULL;
    t->Right->Right      = NULL;
    t->Right->D          = D;
/*
    for ( d=0; d<D; d++ ){
        t->Right->Min[d]  = t->Min[d];
        t->Right->Max[d]  = t->Max[d];
        t->Right->Diff[d] = t->Diff[d];
    }
    t->Right->Min[q]  = Pos[nLeft];
    t->Right->Max[q]  = Pos[t->nData-1];
    t->Right->Diff[q] = t->Right->Max[q] - t->Right->Min[q];
*/
    for (d=0; d<D; d++) {
        t->Right->Min[d] =  9e99;
        t->Right->Max[d] = -9e99;
    }
    for (j=0; j<nRight; j++){
        t->Right->Data[j].Position = (double *) calloc( D, sizeof(double) );
        t->Right->Data[j].Id = t->Data[ Idx[j+nLeft] ].Id;
        for (d=0; d<D; d++) {
            c = t->Data[ Idx[j+nLeft] ].Position[d];
            if ( c < t->Right->Min[d] ) t->Right->Min[d] = c;
            if ( c > t->Right->Max[d] ) t->Right->Max[d] = c;
            t->Right->Data[j].Position[d] = c;
        }
        t->Right->Data[j].Object = t->Data[ Idx[j+nLeft] ].Object;
    }
    for (d=0; d<D; d++) t->Right->Diff[d] = t->Right->Max[d] - t->Right->Min[d];
    
//for (d=0; d<D; d++)
//printf("Left:  npnts: %ld   q: %d  Min[%d] Max[%d] = %g %g\n", nLeft, q, d, d, t->Left->Min[d], t->Left->Max[d] );
//printf("\n");
    //printf("Left:  npnts: %ld   Min[%d] Max[%d] = %g %g\n", nLeft, q, q, t->Left->Min[q], t->Left->Max[q] );
    //if ( nLeft < 4 ) {
    //    for (j=0; j<nLeft; j++){
    //        printf(" ( %g, %g, %g )", t->Left->Data[j].Position[0], t->Left->Data[j].Position[1], t->Left->Data[j].Position[2]);
    //    }
    //    printf("\n");
    //}
//for (d=0; d<D; d++)
//printf("Right:  npnts: %ld   q: %d  Min[%d] Max[%d] = %g %g\n", nRight, q, d, d, t->Right->Min[d], t->Right->Max[d] );
//printf("\n");
    //printf("Right: npnts: %ld   Min[%d] Max[%d] = %g %g\n", nRight, q, q, t->Right->Min[q], t->Right->Max[q] );
    //if ( nRight < 4 ) {
    //    for (j=0; j<nRight; j++){
    //        printf(" ( %g, %g, %g )", t->Right->Data[j].Position[0], t->Right->Data[j].Position[1], t->Right->Data[j].Position[2]);
    //    }
     //   printf("\n");
    //}
    //printf("\n");


    free( Pos );
    free( Idx );


    /*
     * Zero the Data count in this node.
     * Free memory allocated to Parent data field -- its not a leaf anymore.
     */
    for (j=0; j<t->nData; j++) free( t->Data[j].Position );
    free( t->Data ); t->Data = NULL;
    t->nDataBelow = t->nData;
    t->nData = 0;


    /*
     *  Subdivide if there are too many objects in a node
     */
    if ( ( t->Left->Level  < KDTREE_MAX_LEVEL ) && ( t->Left->nData  > KDTREE_MAX_DATA_PER_NODE ) ) Lgm_KdTree_SubDivideVolume( t->Left, kt );
    if ( ( t->Right->Level < KDTREE_MAX_LEVEL ) && ( t->Right->nData > KDTREE_MAX_DATA_PER_NODE ) ) Lgm_KdTree_SubDivideVolume( t->Right, kt );


    return;

}


/**
 *  Finds the k Nearest Neighbors (kNN) of a query point q given that the set
 *  of data is stored as a KdTree. 
 *
 *    \param[in]     q          Query position (D-dimensional) . I.e. the point we want to find NNs for.
 *    \param[in]     Root       Root node of KdTree.
 *    \param[in]     K          Number of NNs to find.
 *    \param[in]     MaxDist2   Threshold distance^2 beyond which we give up on finding
 *                              NNs.  (i.e. could find them, but we arent interested
 *                              because they'd be too far from our query point to be
 *                              useful).
 *    \param[out]    Kgot       Number of NNs (within MaxDist2) that we actually found.
 *    \param[out]    kNN        List of kNN Data items. Sorted (closest first).
 *
 *    \returns       KDTREE_KNN_SUCCESS         Search succeeded.
 *                   KDTREE_KNN_TOO_FEW_NNS     Search terminated because we couldnt find K NNs that were close enough.
 *                   KDTREE_KNN_NOT_ENOUGH_DATA KdTree doesnt contain enough data points.
 *
 *    \author        Mike Henderson
 *    \date          2013
 *
 */
int Lgm_KdTree_kNN( double *q_in, int D, Lgm_KdTree *KdTree, int K, int *Kgot, double MaxDist2, Lgm_KdTreeData *kNN ) {

    int                  k, d, done;
    Lgm_KdTree_pQueue   *PQ, *p;
    Lgm_KdTreeNode      *Root;
    double              *q;
    //static count=0;


    Root = KdTree->Root;
    *Kgot = 0;


    /*
     * Check to see if there are enough points.
     * Bailout with error flag -1 if not.
     */
    if ( Root->nDataBelow < K ) return( KDTREE_KNN_NOT_ENOUGH_DATA );


    /*
     *  Add Root Node to the Priority Queue.
     */
    PQ = NULL;
    //Lgm_KdTreeScalePosition( q_in, &q, KdTree );
    q = (double *) calloc( D, sizeof(double) );
    for (d=0; d<D; d++) q[d] = q_in[d];
    Lgm_KdTree_InsertNode( Root, q, &PQ, MaxDist2 );


    /*
     *  Process items on the Priority Queue until we are done.
     */
    k    = 0;       // havent found any NNs yet.
    done = FALSE;
    while( !done ) {

//        Lgm_KdTree_PrintPQ( &PQ ); //only for debugging


        /*
         * Pop the highest priority item off of the Priority Queue.
         * This is the closest object to the query point. Note that
         * it could be a cell or it could be a point.
         */
        p = Lgm_KdTree_PopObj( &PQ );

        if ( p ) {

            if ( p->IsPoint ) {
                /*
                 * Since a point is now the closest object, it must be one of
                 * the kNN.  But, if its too far away, then we are done with
                 * the search -- we couldnt find K NNs close enough to the
                 * query point
                 */
                if ( p->MinDist2 > MaxDist2 ) {
                    while ( (p = Lgm_KdTree_PopObj( &PQ )) ) {free( p );}
                    free( q );
                    return( KDTREE_KNN_TOO_FEW_NNS );
                }

                /*
                 * Otherwise, add this point as the next NN and continue.
                 */
                kNN[ k   ]       = p->Obj->Data[p->j];
                kNN[ k++ ].Dist2 = p->MinDist2; // save the dist2 into the data struct
                *Kgot = k;

            } else {

                /*
                 * If the object is node, then descend one level closer to
                 * the leaf node that is closest to the query point. This also
                 * adds all nodes we encounter along the way to the PQ. Ignore
                 * any nodes that is farther than MaxDist2 away from query point.
                 */
                Lgm_KdTree_DescendTowardClosestLeaf( p->Obj, &PQ, q, MaxDist2 );

            }

        } else {

            // The priority queue is empy -- there is nothing more to search.
            done = TRUE;

        }

        /*
         * object is no longer needed so free up the space we allocated for it.
         */
        free( p );

        if ( k >= K ) done = TRUE;


        //Lgm_KdTree_PrintPQ( &PQ ); //only for debugging

    }

    /*
     * Free all remaining objects on the PQ
     */
    while ( (p = Lgm_KdTree_PopObj(&PQ)) ) {free( p );}


    //Lgm_KdTree_PrintPQ( &PQ ); //only for debugging

    ++(KdTree->kNN_Lookups);


    /*
     *  return success
     */
    free( q );
    return( KDTREE_KNN_SUCCESS );

}








/**
 *  Insert a Node object into the Priority Queue (PQ). This computes distance
 *  and puts it into the PQ at the right place to maintain a `mindist' order
 *  (the top object of the queue. This routine is used by Lgm_KdTree_kNN().
 *
 *      \param[in]      Node        Pointer to a node in the kdtree
 *      \param[in]      q           The D-dimensional query point.
 *      \param[in,out]  PQ          The "priority queue".
 *      \param[in]      MaxDist2    The maximum distance (squared) to care
 *                                  about. Square distances beyond this value are ignored.
 *
 *      \returns        The distance between the point and the node in
 *                      normalized coordinates.
 *
 *      \author         Mike Henderson
 *      \date           2013
 *
 */
double Lgm_KdTree_InsertNode( Lgm_KdTreeNode *Node, double *q, Lgm_KdTree_pQueue **PQ, double MaxDist2 ) {

    int                 done, InsertAtEnd;
    double              dist;
    Lgm_KdTree_pQueue  *new, *p;

    /*
     * If the cell's MinDist2 > MaxDist2, then we should just ignore
     * the entire cell becuase none of its contents can possibly have points
     * that are < MaxDist2
     */
    dist = Lgm_KdTree_MinDist( Node, q );
    if ( dist > MaxDist2 ) return(1e31);



    // Allocate an element for the PQ
    new = ( Lgm_KdTree_pQueue *) calloc( 1, sizeof( Lgm_KdTree_pQueue) );
//++total_callocs;
//printf("Lgm_KdTree_InsertNode = %d    MinDist2 = %g\n", total_callocs, dist);

    // Add Info
    new->Obj      = Node;
    new->MinDist2 = dist;
    //printf("Lgm_KdTree_InsertNode: MinDist2 = %g   Level = %d\n", dist, Node->Level);
    new->IsPoint  = FALSE;

    if ( *PQ == NULL ) {
        // nothing in the Priority Queue.
        *PQ = new;
        new->Prev = NULL;
        new->Next = NULL;
    //printf("InsertNode: 1.\n"); //Lgm_KdTree_PrintPQ( PQ ); //only for debugging
        return( dist );
    } else {

        // find where to add the new element
        p = *PQ;
        done = FALSE;
        InsertAtEnd = FALSE; // always insert new object before the object pointed at by p unless this is true
                             // then it goes at the end.
        while ( !done ) {
            if ( new->MinDist2 <= p->MinDist2 ) {
                done = TRUE;
            } else if ( p->Next == NULL ) {
                done = TRUE;
                InsertAtEnd = TRUE;
            } else {
                p = p->Next;
            }
            //printf("InsertNode: 2.\n"); //Lgm_KdTree_PrintPQ( PQ ); //only for debugging
        }



        //if ( p == NULL ) {
        if ( InsertAtEnd ) {
            // Insert at bottom (i.e. as last element)
            p->Next   = new;
            new->Next = NULL;
            new->Prev = p;
            //printf("InsertNode: 3.\n"); //Lgm_KdTree_PrintPQ( PQ ); //only for debugging
        } else {
            if ( p->Prev == NULL ) {
                // Insert at top (i.e. as first element)
                new->Next = p;
                new->Prev = NULL;
                p->Prev   = new;
                *PQ        = new;
            } else {
                // Insert before p
                new->Next     = p;
                new->Prev     = p->Prev;
                p->Prev->Next = new;
                p->Prev       = new;
            }
            //printf("InsertNode: 4.\n"); //Lgm_KdTree_PrintPQ( PQ ); //only for debugging
        }

        //printf("InsertNode: 5.\n"); //Lgm_KdTree_PrintPQ( PQ ); //only for debugging
    }

    return( dist );


}

/**
 *  Insert a point object into the priority queue. This computes distance and
 *  puts it into the priority queue at the right place.  This routine is used
 *  by Lgm_KdTree_kNN().
 *
 *      \param[in]      Node        Pointer to the node in the kdtree that contains the data point to add.
 *      \param[in]      j           array index of the data point (contained within the given cell) to add.
 *      \param[in]      q           The D-dimensional query point. (D is carried along in each node).
 *      \param[in,out]  PQ          The "priority queue".
 *
 *      \returns        void
 *
 *      \author         Mike Henderson
 *      \date           2013
 *
 */
void Lgm_KdTree_InsertPoint( Lgm_KdTreeNode *Node, int j, double *q, Lgm_KdTree_pQueue **PQ ) {

    int                 d, D;
    Lgm_KdTree_pQueue  *new, *p;
    double              delta, MinDist2;

    D = Node->D;

    // compute distance^2 between point and query point
    for ( MinDist2=0.0, d=0; d<D; d++ ) {
        delta     = q[d] - Node->Data[j].Position[d];
        MinDist2 += delta*delta;
    }


    // Allocate an element for the PQ
    new = (Lgm_KdTree_pQueue *) calloc( 1, sizeof(Lgm_KdTree_pQueue) );
//++total_callocs;
//printf("Lgm_KdTree_InsertPoint_callocs = %d\n", total_callocs);
//printf("MinDist2 = %g\n", MinDist2);

    // Add Info
    new->Obj      = Node;
    new->MinDist2 = MinDist2;
//printf("Lgm_KdTree_InsertPoint: MinDist2 = %g\n", MinDist2);
    new->IsPoint  = TRUE;
    new->j        = j;

    if ( *PQ == NULL ) {
        // nothing in the Priority Queue.
        *PQ = new;
        new->Next = NULL;
        new->Prev = NULL;
        //printf("1.\n"); Lgm_KdTree_PrintPQ( PQ ); //only for debugging
        return;
    } else {

        // find where to add the new element
        p = *PQ;
        while ( p ) {

            if ( (new->MinDist2 <= p->MinDist2) && (p->Prev != NULL) ) {
                // Insert before p
                new->Next     = p;
                new->Prev     = p->Prev;
                p->Prev->Next = new;
                p->Prev       = new;
        //printf("2.\n"); Lgm_KdTree_PrintPQ( PQ ); //only for debugging
                return;
            } else if ( (new->MinDist2 <= p->MinDist2) && (p->Prev == NULL) ) {
                // Insert at top (i.e. as first element)
                new->Next = p;
                new->Prev = NULL;
                p->Prev   = new;
                *PQ        = new;
        //printf("3.\n"); Lgm_KdTree_PrintPQ( PQ ); //only for debugging
                return;
            } else if ( (new->MinDist2 > p->MinDist2) && (p->Next == NULL) ) {
                // Insert at bottom (i.e. as last element)
                p->Next   = new;
                new->Next = NULL;
                new->Prev = p;
        //printf("4.\n"); Lgm_KdTree_PrintPQ( PQ ); //only for debugging
                return;
                return;
            }
            p = p->Next;
        }
    }


}









/**
 *  Pop an Object off the top of the PQ. This routine
 *  is used by Lgm_KdTree_kNN().
 *
 *      \param[in,out]  PQ          The "priority queue".
 *
 *      \returns        pointer to item popped off the priority queue.
 *
 *      \author         Mike Henderson
 *      \date           2013
 *
 */
Lgm_KdTree_pQueue *Lgm_KdTree_PopObj( Lgm_KdTree_pQueue **PQ ) {

    Lgm_KdTree_pQueue *p;

    p = *PQ;

    if ( p != NULL ) {
        if ( p->Next != NULL ) {
            *PQ = p->Next;
            (*PQ)->Prev = NULL;
        } else {
            *PQ = NULL;
        }
    }

    return( p );

}


/**
 *  This routine computes the minimum distance between a point and a KdTree
 *  node (each node maintains a region of D-dimensional space). Basically, we
 *  compute the closest distance between the query point and any point on the
 *  hyper-rectangle defined in the node.
 *
 *      \param[in]      Node    Pointer to a node in the kdtree
 *      \param[in]      q       The D-dimensional query point.
 *
 *      \returns        The minimum distance between the point and the node.
 *
 *      \author         Mike Henderson
 *      \date           2013
 *
 */
double  Lgm_KdTree_MinDist( Lgm_KdTreeNode *Node, double *q ) {

    int     d, D;
    double  distance2, delta;

    D = Node->D;

    /*
     * The region in each node is defined by the Min[] and Max[] arrays.
     * E.g., Min[0]/Max[0] is the bounds in the zeroeth dimension, etc.
     */

    for ( distance2=0.0, d=0; d<D; d++ ) {
        if      ( q[d] < Node->Min[d] ) { delta = Node->Min[d] - q[d]; distance2 += delta*delta; }
        else if ( q[d] > Node->Max[d] ) { delta = q[d] - Node->Max[d]; distance2 += delta*delta; }
    }

    //printf("Lgm_KdTree_MinDist(): distance2 = %g\n", distance2);
    return( distance2 );

}

/**
 *   Descend to the leaf node that is closest to the query point. This routine
 *   is used by Lgm_KdTree_kNN().
 *
 *      \param[in]      Node        Pointer to a node in the kdtree
 *      \param[in]      PQ          The "priority queue"3D query point.
 *      \param[in]      q           The 3D query point.
 *      \param[in]      MaxDist2    The maximum distance (squared) to care
 *                                  about. Square distances beyond this value are ignored.
 *
 *      \returns        pointer one level closer to leaf node that is closest to query point.
 *
 *      \author         Mike Henderson
 *      \date           2013
 *
 */
void Lgm_KdTree_DescendTowardClosestLeaf( Lgm_KdTreeNode *Node, Lgm_KdTree_pQueue **PQ, double *q, double MaxDist2 ) {

    unsigned long int   j;
    double              dL, dR;
    Lgm_KdTreeNode     *L, *R;


    L = Node->Left;
    R = Node->Right;



    if ( L && R ) {

        /*
         * This is not a leaf node -- should be no points contained in this
         * node. There are both left and right subtrees defined.  Throw both on
         * the queue for further processing. Make sure the closest one is
         * searched first.
         */

        dL = Lgm_KdTree_MinDist( L, q );
        dR = Lgm_KdTree_MinDist( R, q );
        if ( dL < dR ) {
            if (dL < MaxDist2) Lgm_KdTree_InsertNode( L, q, PQ, MaxDist2 );
            if (dR < MaxDist2) Lgm_KdTree_InsertNode( R, q, PQ, MaxDist2 );
        } else {
            if (dR < MaxDist2) Lgm_KdTree_InsertNode( R, q, PQ, MaxDist2 );
            if (dL < MaxDist2) Lgm_KdTree_InsertNode( L, q, PQ, MaxDist2 );
        }

    } else if ( L ) {

        /*
         * This is not a leaf node -- should be no points contained in this
         * node. Only left subtree defined.  Throw both on the queue for
         * further processing. 
         */
        dL = Lgm_KdTree_MinDist( L, q );
        if (dL < MaxDist2) Lgm_KdTree_InsertNode( L, q, PQ, MaxDist2 );

    } else if ( R ) {

        /*
         * This is not a leaf node -- should be no points contained in this
         * node. Only right subtree defined.  Throw both on the queue for
         * further processing. 
         */
        dR = Lgm_KdTree_MinDist( R, q );
        if (dR < MaxDist2) Lgm_KdTree_InsertNode( R, q, PQ, MaxDist2 );

    } else {

        /*
         * This is a leaf node. Add all points found here onto Priority Queue.
         */
        for (j=0; j<Node->nData; j++) Lgm_KdTree_InsertPoint( Node, j, q, PQ );

    }

    return;
}

/**
 *   Prints the contents of the priority queue. For dubugging only. This routine
 *   is used by Lgm_KdTree_kNN().
 *
 *      \param[in]  PQ          The "priority queue".
 *
 *      \returns        void
 *
 *      \author         Mike Henderson
 *      \date           2013
 *
 */
void Lgm_KdTree_PrintPQ( Lgm_KdTree_pQueue **PQ ) {

    long int    i=0;
    Lgm_KdTree_pQueue      *p;

    p = *PQ;

    printf("---- Contents of Priority Queue  ---------\n");
    while ( p != NULL ) {

    if ( p->IsPoint ) {
        printf("Item # %03ld  Point: MinDist2 = %g\n", i++, p->MinDist2);
    } else {
        printf("Item # %03ld   Node: MinDist2 = %g\n", i++, p->MinDist2);
    }
    p = p->Next;

    }
    printf("\n\n");
}


