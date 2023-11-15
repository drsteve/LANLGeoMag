/*! \file   Lgm_KdTree.c
 *  \brief  Set of routines for creating kdtrees (K-dimensional trees) and find k Nearest Neighbors.
 *  \author M.G. Henderson
 *  \date   2013
 */

#pragma GCC push_options
#pragma GCC optimize ("O3")

#include "Lgm/Lgm_KdTree.h"
#include "Lgm/quicksort.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/timeb.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

// We are missing the Lgm_KdTree_Free() routine???



/** 
 *   \brief
 *      Store given N-dimensional data into a D-dimensional KD-tree data structure.
 *
 *   \details
 *      Given arrays of positions and data, this routine recursively partitions
 *      the data into a kdtree data structure. 
 *
 *   \param[in]      Points      An array of position vectors in D-dimensional
 *                               space. ObjectPoints[d][n] is the dth component of the nth point.
 *
 *   \param[in]      Objects     An array of objects in D-dimensional space.
 *                               Objects[n] is the nth pointer to an object.  This is an array of 'void *'
 *                               pointers. This allows the user to use any object type here, so long as
 *                               they are properly typecast.
 *
 *   \param[in]      N           Number of points.
 *
 *   \param[in]      D           Number of dimensions.
 *
 *   \returns        returns a pointer to the a KdTree structure. User is
 *                   responsible to freeing this with Lgm_KdTree_Free( )
 *
 *   \author         Mike Henderson
 *   \date           2013
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
        t->Data[j].Object = Objects[j];
        //double *bbb;
        //bbb = (double *)Objects[j];
        //printf("bbb = %g %g %g\n", bbb[0], bbb[1], bbb[2]);

    }
    for (d=0; d<D; d++) t->Diff[d] = t->Max[d] - t->Min[d];




    kt = (Lgm_KdTree *) calloc( 1, sizeof( Lgm_KdTree) );

    kt->kNN_Lookups   = 0;
    kt->SplitStrategy = LGM_KDTREE_SPLIT_MAXRANGE;
    
    Lgm_KdTree_SubDivideVolume( t, kt );
    kt->Root = t;

    kt->PQN = Lgm_pQueue_Create( 5000 );
    kt->PQP = Lgm_pQueue_Create( 5000 );

    return( kt );
    

}

/*
 * Frees a Lgm_KdTree
NOT FINISHED!!!!!
NEED TO TRAVERSE TREE AND DEALLOCATE EVERYTHING.
 */
void Lgm_KdTree_Free( Lgm_KdTree *kt ) {


    if ( kt == NULL ) return;
    if ( kt->Root != NULL ) {
        Lgm_FreeKdTreeNode( kt->Root );
    }

    Lgm_pQueue_Destroy( kt->PQN );
    Lgm_pQueue_Destroy( kt->PQP );

    free( kt );

}




/*
 * Copy a KdTree structure. 
 *
 * This is NOT a full copy.  Here we do not copy the actual tree. Instead, we
 * copy the other items like PQN, PQP, etc...  The idea is that we want to be
 * able to use a tree for lookups in parallel threads. In this mode, the tree
 * and the data in the tree are not modified.
 *
 */
Lgm_KdTree *Lgm_KdTree_CopyLite( Lgm_KdTree *ks ) {

    Lgm_KdTree *kt;

    // alloc mem for target 
    kt = (Lgm_KdTree *) calloc( 1, sizeof( Lgm_KdTree) );

    kt->kNN_Lookups   = 0;
    kt->SplitStrategy = ks->SplitStrategy;
    
    /* copy pointer to the actual root nood of tree.
     * Do not free this.
     */
    kt->Root = ks->Root;

    /*
     * When done, these will need to be properly free'd, without freeing the tree
     * (the master copy should only free it once.)
     */
    kt->PQN = Lgm_pQueue_Create( 5000 );
    kt->PQP = Lgm_pQueue_Create( 5000 );

    return( kt );
    

}

/*
 * Frees just the extra bits alloced by Lgm_KdTree_CopyLite(). Leaves initial
 * tree unfreed.
 */
void Lgm_KdTree_FreeLite( Lgm_KdTree *kt ) {

    if ( kt == NULL ) return;

    if ( kt->PQN != NULL ) Lgm_pQueue_Destroy( kt->PQN );
    if ( kt->PQP != NULL ) Lgm_pQueue_Destroy( kt->PQP );

    free( kt );

}



/**
 *  \brief
 *      Create a root-level node for an kdtree.
 *
 *  \details
 *      This routine allocates memory for a root-level node of a KdTree with
 *      dimension D. The returned pointer to an Lgm_KdTreeNode must be
 *      destroyed with a call to Lgm_FreeKdTree().
 *
 *   \param[in]      D      Integer dimension of the KdTree.
 *
 *   \returns        pointer to Lgm_KdTreeNode
 *
 *   \author         Mike Henderson
 *   \date           2013
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
 *  \brief
 *      Recursively free the binary KdTree created int Lgm_KdTree_Init()
 *
 *  \details
 *      Recursive de-allocation of all memory used.
 *
 *   \param[in]      *Node      Pointer to Root Node.
 *
 *   \returns        void
 *
 *   \author         Mike Henderson
 *   \date           2017
 *
 */
void Lgm_FreeKdTreeNode( Lgm_KdTreeNode *Node ) {

    int j;

    if ( Node == NULL ) return;

    if ( Node->nData > 0 ) {

        /*
         * Its a leaf node that stores data -- free it. (Should not have Left
         * or Right set).
         */
        if ( Node->Data != NULL ) {
            for (j=0; j<Node->nData; j++) {
                if ( Node->Data[j].Position != NULL ) {
                    free( Node->Data[j].Position );
                }
            }
            free( Node->Data );
        }

    } else {

        if ( Node->Left  != NULL ) Lgm_FreeKdTreeNode( Node->Left  );
        if ( Node->Right != NULL ) Lgm_FreeKdTreeNode( Node->Right );

    }

    if ( Node->Min  != NULL ) free( Node->Min  );
    if ( Node->Max  != NULL ) free( Node->Max  );
    if ( Node->Diff != NULL ) free( Node->Diff );
    free( Node );

    return;

}




/**
 *  \brief
 *      Recursively subdivide volume. Stop when the number of data points per node
 *      is lower than the threshold given by KDTREE_MAX_DATA_PER_NODE.
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
    double             *Pos;
    double               V, c, MaxDiff;
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
    V      = Pos[nLeft];
    t->d      = q;   // splitting dimension
    t->CutVal = V;   // splitting value (values less than this go to the left, values greater than or equal go the right).
//printf("d cut = %d %g\n", d, t->CutVal);

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
//free( t->Min ); 
//free( t->Max ); 
//free( t->Diff ); 
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
 *  \brief
 *      Finds the k Nearest Neighbors (kNN) of a query point q given that the set
 *      of data is stored as a KdTree. 
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
int Lgm_KdTree_kNN( double *q, int D, Lgm_KdTree *KdTree, int K, int *Kgot, double MaxDist2, Lgm_KdTreeData *kNN ) {

    int                  k, index;
    Lgm_KdTreeNode      *Root, *p;
    double              dist, maxd2;
    Lgm_pQueue          *PQP;
    Lgm_pQueue_Node     FarthestPoint;

    maxd2 = MaxDist2;

    // reset priority queue heap for points
    PQP = KdTree->PQP;
    PQP->HeapSize = 0;

    Root = KdTree->Root;
    *Kgot = 0;

    /*
     * Check to see if there are enough points.
     * Bailout with error flag -1 if not.
     */
    if ( Root->nDataBelow < K ) return( KDTREE_KNN_NOT_ENOUGH_DATA );

    Lgm_KdTree_DepthFirstSearch( Root, PQP, K, q, &maxd2 );

// these are in reverse order -- fix...
    k = 0;
    while ( Lgm_pQueue_Pop( &FarthestPoint, PQP ) && (k <= K) ) { 

        dist  = -FarthestPoint.key;
        index = FarthestPoint.index;
        p     = (Lgm_KdTreeNode *)FarthestPoint.Data;

        if ( dist <= MaxDist2 ) {
            kNN[ k   ]       = p->Data[index];
            kNN[ k++ ].Dist2 = dist; // save the dist2 into the data struct
            *Kgot = k;
        }
        
    }

    //Lgm_KdTree_PrintPQ( &PQ ); //only for debugging
// not thread safe ?
//    ++(KdTree->kNN_Lookups);

    /*
     *  return success
     */
    if (k==K) {
        return( KDTREE_KNN_SUCCESS );
    } else {
        return( KDTREE_KNN_TOO_FEW_NNS );
    }

}


/**
 *  \brief
 *      This routine computes the minimum distance between a point and a KdTree
 *      node (each node maintains a region of D-dimensional space). 
 *
 *  \details
 *      Basically, we compute the closest distance between the query point and
 *      any point on the hyper-rectangle defined in the node.
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
 *  \brief
 *      Determine if a search of this subtree is waranted.
 *
 *  \details
 *      Rather than compute the entire minimum distance between the point and
 *      the hyper-rectangular boundaing box, this routine bails early if the
 *      min-dist^2 exceeds a maximum threshold.
 *
 *      \param[in]      Node    Pointer to a node in the kdtree
 *      \param[in]      q       The D-dimensional query point.
 *      \param[in]      md2     Max dist^2 to consider.
 *
 *      \returns        The minimum distance between the point and the node.
 *
 *      \author         Mike Henderson
 *      \date           2013
 *
 */
int  Lgm_KdTree_DoSearch( Lgm_KdTreeNode *Node, double *q, double md2 ) {

    int     d;
    double  distance2, delta, qd;


    /*
     * The region in each node is defined by the Min[] and Max[] arrays.
     * E.g., Min[0]/Max[0] is the bounds in the zeroeth dimension, etc.
     */

    for ( distance2=0.0, d=0; d<Node->D; d++ ) {
        qd = *(q + d);
        if      ( qd < Node->Min[d] ) { delta = Node->Min[d] - qd; distance2 += delta*delta; if ( distance2 > md2 ) return( FALSE ); }
        else if ( qd > Node->Max[d] ) { delta = qd - Node->Max[d]; distance2 += delta*delta; if ( distance2 > md2 ) return( FALSE ); }
    }
    if ( distance2 < md2 ) return( TRUE );

}

/**
 *  \brief
 *      Recursively descend to the leaf node that is closest to the query point
 *      floowed by the farthest node. This routine is used by Lgm_KdTree_kNN().
 *
 *      \param[in]      Node        Pointer to a node in the kdtree
 *      \param[in]      PQ          The "priority queue"
 *      \param[in]      q           The D-dim query point.
 *      \param[in]      MaxDist2    The maximum distance (squared) to care
 *                                  about. Square distances beyond this value are ignored.
 *
 *      \returns        pointer one level closer to leaf node that is closest to query point.
 *
 *      \author         Mike Henderson
 *      \date           2013
 *
 */
void Lgm_KdTree_DepthFirstSearch( Lgm_KdTreeNode *Node, Lgm_pQueue *PQP, int K, double *q, double *MaxDist2 ) {

    int                 i, d;
    double              *pos, dist, ddd, q_d, v_d, diff;
    Lgm_KdTreeNode     *L, *R;
    Lgm_pQueue_Node     New, FarthestPoint;


    
    L = Node->Left;
    R = Node->Right;


    New.index = -1;

// This could probably be sped up with recursion....?
// Note that we add both thew L and R at once. Recursion would wait to add the farthest choice which could have been pruned out by then....

    if ( L && R ) {

        /*
         * This is not a leaf node -- should be no points contained in this
         * node. There are both left and right subtrees defined.  Search the
         * best one first. Then search the farthest one.
         */
        q_d = q[ Node->d ]; // The component of the query vector in the cut-dimension.
        v_d = Node->CutVal; // The value of the cutting line (in dimension Node->d).

    
        if ( q_d < v_d ) {
            // follow the left subtree first -- its the closest.
            Lgm_KdTree_DepthFirstSearch( L, PQP, K, q, MaxDist2 );

            /*
             * Only search the other subtree if there is a chance a point can
             * be found there.  So far, we know that the point is on the left
             * side of the cut line. Therefore, the dist to the cutline is
             * greater than the current MaxDist2 value, the right subtree isa
             * too far away.  Note: this test is just based on the component
             * distances, not the true distance. We couyld try the true
             * distance though...
             */
            diff = q_d - v_d;
            if ( diff*diff < *MaxDist2 ) {
                 //if ( Lgm_KdTree_MinDist( R, q ) < *MaxDist2 ) Lgm_KdTree_DepthFirstSearch( R, PQP, K, q, MaxDist2 ); 
                 if ( Lgm_KdTree_DoSearch( R, q, *MaxDist2 ) ) Lgm_KdTree_DepthFirstSearch( R, PQP, K, q, MaxDist2 ); 
            }

        } else {
            // follow the right subtree first -- its the closest
            Lgm_KdTree_DepthFirstSearch( R, PQP, K, q, MaxDist2 );

            // Search left if there is a chance of finding a point.
            diff = q_d - v_d;
            if ( diff*diff < *MaxDist2 ) {
                //if ( Lgm_KdTree_MinDist( L, q ) < *MaxDist2 ) Lgm_KdTree_DepthFirstSearch( L, PQP, K, q, MaxDist2 ); 
                if ( Lgm_KdTree_DoSearch( L, q, *MaxDist2 ) ) Lgm_KdTree_DepthFirstSearch( L, PQP, K, q, MaxDist2 ); 
            }
        }

    } else if ( L )  {

        Lgm_KdTree_DepthFirstSearch( L, PQP, K, q, MaxDist2 );

    } else if ( R ) {

        Lgm_KdTree_DepthFirstSearch( R, PQP, K, q, MaxDist2 );

    } else {


        /*
         *  This is a leaf node that may contain multiple points. Add them all
         *  to the Points Priority Queue. Maintain in max order (by negating the distance).
         */
        for (i=0; i<Node->nData; i++ ) {

            pos = Node->Data[i].Position;

            if ( PQP->HeapSize < K ) {
                // add point no matter what if we dont have a full set yet.
                for ( dist=0.0, d=0; d<Node->D; d++) { diff = pos[d] - q[d]; dist += diff*diff; }
                New.key   = -dist;
                New.index = i;
                New.Data  = (void *)Node;
                Lgm_pQueue_Insert( &New, PQP );
                if ( PQP->HeapSize >= K ) *MaxDist2 = -PQP->HeapArray[1].key; // reset the MaxDist2 value -- all closer points must be closer than this.
            } else {
                // only add point if its closer than the farthest.
                ddd = -PQP->HeapArray[1].key;
                for ( dist=0.0, d=0; d<Node->D; d++) { 
                    diff = pos[d] - q[d]; 
                    dist += diff*diff; 
                    if ( dist >= ddd ) break; // its already bigger -- so stop expending effort.
                }

                if ( dist < ddd ) {
                    New.key   = -dist;
                    New.index = i;
                    New.Data  = (void *)Node;

                    Lgm_pQueue_Pop( &FarthestPoint, PQP );  // remove the farthest
                    Lgm_pQueue_Insert( &New, PQP );         // insert new (closer) point
                    *MaxDist2 = -PQP->HeapArray[1].key;     // reset the MaxDist2 value -- all closer points must be closer than this.
                }


            }

        }

    }

    return;
}

/**
 *  \brief
 *      Prints the contents of the priority queue. For dubugging only. This routine
 *      is used by Lgm_KdTree_kNN(). THIS ROUTINE IS NO LONGER VALID..
 *
 *
 *      \param[in]  PQ          The "priority queue".
 *
 *      \returns        void
 *
 *      \author         Mike Henderson
 *      \date           2013
 *
 */
void Lgm_KdTree_PrintPQ( Lgm_pQueue *PQ ) {

    long int         i=0;
    double           dist;
    Lgm_pQueue_Node  *p;
    Lgm_KdTreeNode   *N;



    printf("---- Contents of Priority Queue  ---------\n");
    for ( i=1; i<=PQ->HeapSize; i++ ) {

        p = &(PQ->HeapArray[i]);
        

        dist = p->key;
        N    = (Lgm_KdTreeNode *)p->Data;

        if ( N ) {
        if ( N->nData > 0 ) {
            printf("Point: MinDist2 = %g\n", dist );
        } else {
            printf("Node: MinDist2 = %g\n", dist );
        }
        }

    }

    printf("\n\n");
}


/**
 *  \brief
 *      Descend to the leaf node that is closest to the query point. This routine
 *      is used by Lgm_KdTree_kNN2().  Allows for non-recursive version of depth-first search.
 *
 *      \param[in]      Node        Pointer to a node in the kdtree
 *      \param[in]      PQ          The "priority queue"
 *      \param[in]      q           The D-dim query point.
 *      \param[in]      MaxDist2    The maximum distance (squared) to care
 *                                  about. Square distances beyond this value are ignored.
 *
 *      \returns        pointer one level closer to leaf node that is closest to query point.
 *
 *      \author         Mike Henderson
 *      \date           2013
 *
 */
void Lgm_KdTree_DescendTowardClosestLeaf2( Lgm_KdTreeNode *Node, Lgm_pQueue *PQN, Lgm_pQueue *PQP, int K, double *q, double *MaxDist2 ) {

    int                 i, d;
    double              dL, dR, *pos, dist, g, ddd;
    Lgm_KdTreeNode     *L, *R;
    Lgm_pQueue_Node     New, FarthestPoint;


    L = Node->Left;
    R = Node->Right;

    New.index = -1;

    if ( L && R ) {

        /*
         * This is not a leaf node -- should be no points contained in this
         * node. There are both left and right subtrees defined.  Search the
         * best one first. Then search the farthest one.
         */
        dL = Lgm_KdTree_MinDist( L, q );
        dR = Lgm_KdTree_MinDist( R, q );
        if ( dL < dR ) {
            if (dL < *MaxDist2) { New.key  = dL; New.Data = (void *)L; Lgm_pQueue_Insert( &New, PQN ); }
            if (dR < *MaxDist2) { New.key  = dR; New.Data = (void *)R; Lgm_pQueue_Insert( &New, PQN ); }
        } else {
            if (dR < *MaxDist2) { New.key  = dR; New.Data = (void *)R; Lgm_pQueue_Insert( &New, PQN ); }
            if (dL < *MaxDist2) { New.key  = dL; New.Data = (void *)L; Lgm_pQueue_Insert( &New, PQN ); }
        }

    } else if ( L ) {

        /*
         * This is not a leaf node -- should be no points contained in this
         * node. Only left subtree defined.  Throw both on the queue for
         * further processing. 
         */
        dL = Lgm_KdTree_MinDist( L, q );
        if (dL < *MaxDist2) { New.key  = dL; New.Data = (void *)L; Lgm_pQueue_Insert( &New, PQN ); }

    } else if ( R ) {

        /*
         * This is not a leaf node -- should be no points contained in this
         * node. Only right subtree defined.  Throw both on the queue for
         * further processing. 
         */
        dR = Lgm_KdTree_MinDist( R, q );
        if (dR < *MaxDist2) { New.key  = dR; New.Data = (void *)R; Lgm_pQueue_Insert( &New, PQN ); }

    } else {

        /*
         *  This is a leaf node that may contain multiple points. Add them all
         *  to the Points Priority Queue. Maintain in max order (by negating the distance).
         */
        for (i=0; i<Node->nData; i++ ) {

            pos = Node->Data[i].Position;
            for ( dist=0.0, d=0; d<Node->D; d++) { g = pos[d] - q[d]; dist += g*g; }

            if ( PQP->HeapSize < K ) {
                // add point no matter what if we dont have a full set yet.
                New.key   = -dist;
                New.index = i;
                New.Data  = (void *)Node;
                Lgm_pQueue_Insert( &New, PQP );
                if ( PQP->HeapSize >= K ) *MaxDist2 = -PQP->HeapArray[1].key; // reset the MaxDist2 value -- all closer points must be closer than this.
            } else {
                // only add point if its closer than the farthest.
                ddd = -PQP->HeapArray[1].key;
                if ( dist < ddd ) {
                    New.key   = -dist;
                    New.index = i;
                    New.Data  = (void *)Node;

                    Lgm_pQueue_Pop( &FarthestPoint, PQP );  // remove the farthest
                    Lgm_pQueue_Insert( &New, PQP );         // insert new (closer) point
                    *MaxDist2 = -PQP->HeapArray[1].key;     // reset the MaxDist2 value -- all closer points must be closer than this.
                }
            }

        }

    }

    return;
}




/**
 *  \brief
 *      Finds the k Nearest Neighbors (kNN) of a query point q given that the set
 *      of data is stored as a KdTree.  Non-recursive version.
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
int Lgm_KdTree_kNN2( double *q_in, int D, Lgm_KdTree *KdTree, int K, int *Kgot, double MaxDist2, Lgm_KdTreeData *kNN ) {

    int                  k, d, done, index;
    Lgm_KdTreeNode      *Root, *p;
    double              *q;
    double              dist, maxd2;
    double              FarthestPointDist2, ClosestNodeDist2;
    Lgm_pQueue          *PQN, *PQP;
    Lgm_pQueue_Node     New, ClosestNode, FarthestPoint;

    maxd2 = MaxDist2;

    // reset priority queue heap for nodes
    PQN = KdTree->PQN;
    PQN->HeapSize = 0;

    // reset priority queue heap for points
    PQP = KdTree->PQP;
    PQP->HeapSize = 0;



    Root = KdTree->Root;
    *Kgot = 0;


    /*
     * Check to see if there are enough points.
     * Bailout with error flag -1 if not.
     */
    if ( Root->nDataBelow < K ) return( KDTREE_KNN_NOT_ENOUGH_DATA );

    q = (double *) calloc( D, sizeof(double) );

    /*
     *  Add Root Node to the Priority Queue.
     */
    for (d=0; d<D; d++) q[d] = q_in[d];
    dist = Lgm_KdTree_MinDist( Root, q );
    if ( dist > maxd2 ) { free(q); return(0); }

    New.key   = dist;  // lowest dist is highest priority
    New.index = -1;
    New.Data  = (void *)Root;
    Lgm_pQueue_Insert( &New, PQN );




    /*
     *  Process items on the Priority Queue until we are done.
     */
    k    = 0;       // havent found any NNs yet.
    done = FALSE;
    while( !done ) {

        /*
         * Pop the highest priority item off of the Priority Queue.
         * This is the closest object to the query point. 
         */
        //if ( Lgm_pQueue_Peek( &ClosestNode, PQN ) ) {
        if ( PQN->HeapSize > 0 ) {

            //ClosestNodeDist2  = ClosestNode.key;
            ClosestNodeDist2  = PQN->HeapArray[1].key;

            if ( ClosestNodeDist2 < maxd2 ) {

                Lgm_pQueue_Pop( &ClosestNode, PQN );
                ClosestNodeDist2  = ClosestNode.key;
                p = (Lgm_KdTreeNode *)ClosestNode.Data;

                Lgm_KdTree_DescendTowardClosestLeaf2( p, PQN, PQP, K, q, &maxd2 );

            } else {

                done = TRUE;

            }
            if ( PQP->HeapSize >= K ) {
                //Lgm_pQueue_Peek( &FarthestPoint, PQP );
                //FarthestPointDist2  = -FarthestPoint.key;
                FarthestPointDist2  = -PQP->HeapArray[1].key;
                if ( FarthestPointDist2 < ClosestNodeDist2 ) done = TRUE;
            }

        } else {

            done = TRUE;

        }



    }

// these are in reverse order -- fix...
    k = 0;
    while ( Lgm_pQueue_Pop( &FarthestPoint, PQP ) && (k <= K) ) {

        dist  = -FarthestPoint.key;
        index = FarthestPoint.index;
        p     = (Lgm_KdTreeNode *)FarthestPoint.Data;

        kNN[ k   ]       = p->Data[index];
        kNN[ k++ ].Dist2 = dist; // save the dist2 into the data struct
        *Kgot = k;

    }




    //Lgm_KdTree_PrintPQ( &PQ ); //only for debugging

// Not thread safe?
    ++(KdTree->kNN_Lookups);


    /*
     *  return success
     */
    free( q );
    if (k==K) {
        return( KDTREE_KNN_SUCCESS );
    } else {
        return( KDTREE_KNN_TOO_FEW_NNS );
    }

}
#pragma GCC pop_options
