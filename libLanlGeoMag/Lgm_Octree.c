#include "Lgm/Lgm_Octree.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/timeb.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

struct timeb  StartTime;
double ElapsedTime2( struct timeb StartTime );

/**
 *   Store given 3D data into an octree data structure.
 *
 *   Given 1D arrays of vector positions and vector data (e.g. B-field), this
 *   routine recursively partitions the data into an octree data structure. The
 *   octree can be used to very efficiently find nearest neighbors (using
 *   Lgm_Octree_kNN() ).
 *
 *      \param[in]      ObjectPoints   An array of position vectors
 *      \param[in]      ObjectData     An array of data vectors (e.g. B-field)
 *      \param[in]      N              Number of ObjectPoints
 *
 *      \returns        returns a pointer to the an Octree structure. User is
 *                      responsible to freeing this with Lgm_FreeOctree( )
 *
 *      \author         Mike Henderson
 *      \date           2009-2012
 *
 */
Lgm_Octree *Lgm_InitOctree( Lgm_Vector *ObjectPoints, Lgm_Vector *ObjectData, unsigned long int N ) {

    double              Min, Max, Diff;
    unsigned long int   j;
    Lgm_OctreeCell      *t;
    Lgm_Octree          *ot;

    t = Lgm_CreateOctreeRoot();
    t->nData      = N;
    t->nDataBelow = t->nData;
    t->Data       = (Lgm_OctreeData *) calloc( t->nData, sizeof( Lgm_OctreeData ) );
    Max = -9e99;
    Min = 9e99;
    for (j=0; j<t->nData; j++){
        t->Data[j].Position.x = ObjectPoints[j].x;
        t->Data[j].Position.y = ObjectPoints[j].y;
        t->Data[j].Position.z = ObjectPoints[j].z;
        t->Data[j].B.x = ObjectData[j].x;
        t->Data[j].B.y = ObjectData[j].y;
        t->Data[j].B.z = ObjectData[j].z;
        /*
         * Find scaling for data
         */
        if ( ObjectPoints[j].x > Max ) Max = ObjectPoints[j].x;
        if ( ObjectPoints[j].y > Max ) Max = ObjectPoints[j].y;
        if ( ObjectPoints[j].z > Max ) Max = ObjectPoints[j].z;
        if ( ObjectPoints[j].x < Min ) Min = ObjectPoints[j].x;
        if ( ObjectPoints[j].y < Min ) Min = ObjectPoints[j].y;
        if ( ObjectPoints[j].z < Min ) Min = ObjectPoints[j].z;
    }

    /*
     * Scale Data (positions need to be in range of [0-1]
     */
    Diff = Max - Min;
    for (j=0; j<t->nData; j++){
        t->Data[j].Position.x = (t->Data[j].Position.x - Min)/Diff;
        t->Data[j].Position.y = (t->Data[j].Position.y - Min)/Diff;
        t->Data[j].Position.z = (t->Data[j].Position.z - Min)/Diff;
    }


    SubDivideVolume( t );

    ot = (Lgm_Octree *) calloc( 1, sizeof( Lgm_Octree) );
    ot->Min  = Min;
    ot->Max  = Max;
    ot->Diff = Diff;
    ot->Root = t;


    return( ot );

}




/**
 *   Creates a binary string representation of an unsigned int.
 *
 *      \param[in]      n       An unsigned int
 *      \param[out]     Str     Pointer to a string (i.e. an array of chars)
 *
 *      \returns        void
 *
 *      \author         Mike Henderson
 *      \date           2009-2012
 *
 */
void Binary( unsigned int n, char *Str ) {
    int             i, j;
    unsigned int    mask;
    for (j=0, i=15; i>=0; i--, j++){
        mask = 1<<i;
        Str[j] = ( n & mask ) ? '1' : '0';
        if (j==7) Str[++j] = ' ';
    }
    Str[j] = '\0';
}




/**
 *  Recursively free all nodes below the given node.
 *
 *      \param[in]      n       Pointer to a cell in an octree.
 *
 *      \returns        void
 *
 *      \author         Mike Henderson
 *      \date           2009-2012
 *
 */
void Lgm_OctreeFreeBranch( Lgm_OctreeCell *Cell ){

    int i;

    if ( Cell->Octant == NULL ) {

        /*
         *  Then this is a leaf. Free its contents then return.
         */
        if ( Cell->Data ) {
            free( Cell->Data );
            Cell->Data = NULL;
        }

    } else {

        for (i=0; i<8; ++i) Lgm_OctreeFreeBranch( &(Cell->Octant[i]) );
        free( Cell->Octant );
        Cell->Octant = NULL;

    }


    return;


}



/**
 *   Destroys an octree.
 *
 *   Given a pointer to an octree, this routine destroys the structure by
 *   freeing up all allocated memory.
 *
 *      \param[in]      ot   Pointer to an initialized octree
 *
 *      \returns        void
 *
 *      \author         Mike Henderson
 *      \date           2009-2012
 *
 */
void Lgm_FreeOctree( Lgm_Octree *ot ) {
    Lgm_OctreeFreeBranch( ot->Root );
    free( ot->Root );
    free( ot );
}


/**
 *  Create a root-level node for an octree.
 *
 *
 *      \returns        void
 *
 *      \author         Mike Henderson
 *      \date           2009-2012
 *
 */
Lgm_OctreeCell *Lgm_CreateOctreeRoot( ) {

    char            Str[20];
    Lgm_OctreeCell  *Cell;

    /*
     * allocate
     */
    Cell = (Lgm_OctreeCell *) calloc( 1, sizeof( Lgm_OctreeCell ) );
    Cell->nData = 0;
    Cell->Parent = NULL;
    Cell->Octant = NULL;
    Cell->xLocationCode = 0;
    Cell->yLocationCode = 0;
    Cell->zLocationCode = 0;
    Cell->Level = OCTREE_ROOT_LEVEL;
    Cell->h = 0.5;          // Root Cell has cube with face half-width = 0.5
    Cell->Center.x = 0.5;   // Center is at (0.5, 0.5, 0.5)
    Cell->Center.y = 0.5;
    Cell->Center.z = 0.5;

//    Binary( Cell->xLocationCode, Str ); printf("Root  xLocationCode = %s\n", Str);
//    Binary( Cell->yLocationCode, Str ); printf("Root  yLocationCode = %s\n", Str);
//    Binary( Cell->zLocationCode, Str ); printf("Root  zLocationCode = %s\n\n", Str);

    return( Cell );

}







/**
 *   Find the cell in the octree that is identified by the given x, y, and z Location Codes.
 *
 *      \param[in]      Cell            Pointer to a node in the octree
 *      \param[in]      ChildLevel      XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
 *      \param[in]      xLocationCode   unsigned int representation of the binary x location code to find.
 *      \param[in]      yLocationCode   unsigned int representation of the binary y location code to find.
 *      \param[in]      zLocationCode   unsigned int representation of the binary z location code to find.
 *
 *      \returns        Pointer to the cell that was found.
 *
 *      \author         Mike Henderson
 *      \date           2009-2012
 *
 */
Lgm_OctreeCell *Lgm_OctreeTraverseToLocCode( Lgm_OctreeCell *Cell, unsigned int ChildLevel, unsigned int xLocationCode, unsigned int yLocationCode, unsigned int zLocationCode ) {

    unsigned int    BranchBit, Index;

    while( Cell->Octant ) { // loop until we reach a leaf (i.e. until Octant is NULL)

        // set branch bit for the children of Cell
        BranchBit  = 1<<ChildLevel;

        // determine which octant to enter
        Index = (((xLocationCode & BranchBit) >> ChildLevel)
                        + 2*( ((yLocationCode & BranchBit) >> ChildLevel)
                            + 2*((zLocationCode & BranchBit) >> ChildLevel) ) );

        Cell = &(Cell->Octant[Index]);
        --ChildLevel;

    }

    return( Cell );

}



/**
 *   Returns a pointer to the cell that contains the query point
 *
 *      \param[in]      Root    Pointer to the root node of the octree
 *      \param[in]      q       The 3D query point.
 *
 *      \returns        Pointer to the cell that was found.
 *
 *      \author         Mike Henderson
 *      \date           2009-2012
 *
 */
Lgm_OctreeCell *Lgm_LocateNearestCell( Lgm_OctreeCell *Root, Lgm_Vector *q ){

    unsigned int    xLocationCode = (unsigned int)(q->x * OCTREE_MAX_VAL);
    unsigned int    yLocationCode = (unsigned int)(q->y * OCTREE_MAX_VAL);
    unsigned int    zLocationCode = (unsigned int)(q->z * OCTREE_MAX_VAL);
    unsigned int    ChildLevel;
    Lgm_OctreeCell  *Cell, *Result;

    Cell       = Root;
    ChildLevel = Cell->Level - 1;
    Result     = Lgm_OctreeTraverseToLocCode( Cell, ChildLevel, xLocationCode, yLocationCode, zLocationCode );

    return( Result );

}















/**
 *  This routine computes the minimum distance between a point and an octree
 *  cell (i.e. cube).  The dimensions and center of the cube are stored in the
 *  Lgm_OctreeCell structure.
 *
 *      \param[in]      Cell    Pointer to a cell (i.e. node) in the octree
 *      \param[in]      q       The 3D query point.
 *
 *      \returns        The minimum distance between the point and the cell in
 *                      normalized coordinates.
 *
 *      \author         Mike Henderson
 *      \date           2009-2012
 *
 */
double  MinDist( Lgm_OctreeCell *Cell, Lgm_Vector *q ) {

    double  px, py, pz, d, ph, mh, distance2 = 0.0;

    ph  = Cell->h;  // + half width of a cube face
    mh  = -ph;      // - half width of a cube face

    // assume cube is centered at origin, by shifting query point.
    // i.e. transform point into frame where cube is centered at origin
    px = q->x - Cell->Center.x;
    py = q->y - Cell->Center.y;
    pz = q->z - Cell->Center.z;
//printf("MinDist(): q = %g %g %g\n", q->x, q->y, q->z );

    if      ( px < mh ) { d = px+ph; distance2 += d*d; }
    else if ( px > ph ) { d = px-ph; distance2 += d*d; }

    if      ( py < mh ) { d = py+ph; distance2 += d*d; }
    else if ( py > ph ) { d = py-ph; distance2 += d*d; }

    if      ( pz < mh ) { d = pz+ph; distance2 += d*d; }
    else if ( pz > ph ) { d = pz-ph; distance2 += d*d; }

//printf("MinDist(): distance2 = %g\n", distance2);
    return( distance2 );

}




/**
 *  Insert a Cell object into the Priority Queue (PQ). This computes distance
 *  and puts it into the PQ at the right place to maintain a `mindist' order
 *  (the top object of the queue. This routine is used by Lgm_Octree_kNN().
 *
 *      \param[in]      Cell        Pointer to a cell (i.e. node) in the octree
 *      \param[in]      q           The 3D query point.
 *      \param[in,out]  PQ          The "priority queue"3D query point.
 *      \param[in]      MaxDist2    The maximum distance (squared) to care
 *                                  about. Square distances beyond this value are ignored.
 *
 *      \returns        The distance between the point and the cell in
 *                      normalized coordinates.
 *
 *      \author         Mike Henderson
 *      \date           2009-2012
 *
 */
double InsertCell( Lgm_OctreeCell *Cell, Lgm_Vector *q, pQueue **PQ, double MaxDist2 ) {

    int     done, InsertAtEnd;
    double  dist;
    pQueue  *new, *p;


    /*
     * If the cell's MinDist2 > MaxDist2, then we should just ignore
     * the entire cell becuase none of its contents can possibly have points
     * that are < MaxDist2
     */
    dist = MinDist( Cell, q );
    if ( dist > MaxDist2 ) return(9e99);



    // Allocate an element for the PQ
    new = (pQueue *) calloc( 1, sizeof(pQueue) );

    // Add Info
    new->Obj      = Cell;
    new->MinDist2 = dist;
//printf("InsertCell: MinDist2 = %g   Level = %d\n", dist, Cell->Level);
    new->IsPoint  = FALSE;

    if ( *PQ == NULL ) {
        // nothing in the Priority Queue.
        *PQ = new;
	    new->Prev = NULL;
	    new->Next = NULL;
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
        }


        //if ( p == NULL ) {
        if ( InsertAtEnd ) {
            // Insert at bottom (i.e. as last element)
            p->Next   = new;
            new->Next = NULL;
            new->Prev = p;
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
        }

    }

    return( dist );


}

/**
 *  Insert a point object into the priority queue. This computes distance and
 *  puts it into the priority queue at the right place.  This routine is used
 *  by Lgm_Octree_kNN().
 *
 *      \param[in]      Cell        Pointer to the cell (i.e. node) in the octree that contains the data point to add.
 *      \param[in]      j           array index of the data point (contained within the given cell) to add.
 *      \param[in]      q           The 3D query point.
 *      \param[in,out]  PQ          The "priority queue"3D query point.
 *
 *      \returns        void
 *
 *      \author         Mike Henderson
 *      \date           2009-2012
 *
 */
void InsertPoint( Lgm_OctreeCell *Cell, int j, Lgm_Vector *q, pQueue **PQ ) {

    pQueue  *new, *p;
    double  x, y, z, dx, dy, dz, MinDist2;

    MinDist2 = 0.0;

    // compute distance^2 between point and query point
    x = Cell->Data[j].Position.x; dx = q->x-x; MinDist2 += dx*dx;
    y = Cell->Data[j].Position.y; dy = q->y-y; MinDist2 += dy*dy;
    z = Cell->Data[j].Position.z; dz = q->z-z; MinDist2 += dz*dz;




    // Allocate an element for the PQ
    new = (pQueue *) calloc( 1, sizeof(pQueue) );

    // Add Info
    new->Obj      = Cell;
    new->MinDist2 = MinDist2;
//printf("InsertPoint: MinDist2 = %g\n", MinDist2);
    new->IsPoint  = TRUE;
    new->j        = j;

    if ( *PQ == NULL ) {
        // nothing in the Priority Queue.
        *PQ = new;
        new->Next = NULL;
        new->Prev = NULL;
        return;
    } else {
        // find where to add the new element
        p = *PQ;
        while ( p ) {
            if ( new->MinDist2 < p->MinDist2 ) {
                if ( p->Prev == NULL ) {
                    // Insert at top (i.e. as first element)
                    new->Next = p;
                    new->Prev = NULL;
                    p->Prev   = new;
                    *PQ        = new;
                    return;
                } else if ( p->Next == NULL ) {
                    // Insert at bottom (i.e. as last element)
                    p->Next   = new;
                    new->Next = NULL;
                    new->Prev = p;
                    return;
                } else {
                    // Insert before p
                    new->Next     = p;
                    new->Prev     = p->Prev;
                    p->Prev->Next = new;
                    p->Prev       = new;
                    return;
                }
            }
            p = p->Next;
        }
    }


}



/**
 *   Descend to the leaf node that is closest to the query point. This routine
 *   is used by Lgm_Octree_kNN().
 *
 *      \param[in]      Node        Pointer to a node (cell) in thge octree
 *      \param[in]      PQ          The "priority queue"3D query point.
 *      \param[in]      q           The 3D query point.
 *      \param[in]      MaxDist2    The maximum distance (squared) to care
 *                                  about. Square distances beyond this value are ignored.
 *
 *      \returns        pointer to closest leaf node.
 *
 *      \author         Mike Henderson
 *      \date           2009-2012
 *
 */
Lgm_OctreeCell *DescendTowardClosestLeaf( Lgm_OctreeCell *Node, pQueue **PQ, Lgm_Vector *q, double MaxDist2 ) {

    unsigned long int   j;
    double              dist, min_dist;
    int                 i, min_i;
    unsigned int        ChildLevel, BranchBit, Index;
    Lgm_OctreeCell     *Cell, *Result;

    Cell = Node;
    ChildLevel = Cell->Level - 1;


//    while( Cell->Octant ) { // loop until we reach a leaf (i.e. until Octant is NULL)
    if ( Cell->Octant ) { // do this if its a leaf cell

        min_dist = 9e99;
        min_i    = 0;
        for (i=0; i<8; i++){

            // Add Cell to Priority Queue
            dist = InsertCell( &(Cell->Octant[i]), q, PQ, MaxDist2 );

            // keep track of closest octant
            if (dist <= min_dist) {
                min_dist = dist;
//printf("HHHHHHH: min_dist = %g\n", min_dist);
                min_i    = i;
            }

        }

        Cell = &(Cell->Octant[min_i]);

    } else {

        // Add point(s) (if any) to Priority Queue
        for (j=0; j<Cell->nData; j++){
            InsertPoint( Cell, j, q, PQ );
        }

    }

    return( Cell );

}

/**
 *   Prints the contents of the priority queue. For dubugging only. This routine
 *   is used by Lgm_Octree_kNN().
 *
 *      \param[in]  PQ          The "priority queue"3D query point.
 *
 *      \returns        void
 *
 *      \author         Mike Henderson
 *      \date           2009-2012
 *
 */
void PrintPQ( pQueue **PQ ) {

    long int    i=0;
    pQueue      *p;

    p = *PQ;

    printf("---- Contents of Priority Queue  ---------\n");
    while ( p != NULL ) {

	if ( p->IsPoint ) {
	    printf("Item # %03ld  Point: MinDist2 = %g\n", i++, p->MinDist2);
	} else {
	    printf("Item # %03ld   Cell: MinDist2 = %g\n", i++, p->MinDist2);
	}
	p = p->Next;

    }
    printf("\n\n");
}

/**
 *  Pop an Object off the top of the PQ. This routine
 *  is used by Lgm_Octree_kNN().
 *
 *      \param[in,out]  PQ          The "priority queue"3D query point.
 *
 *      \returns        pointer to item popped off the priority queue.
 *
 *      \author         Mike Henderson
 *      \date           2009-2012
 *
 */
pQueue *PopObj( pQueue **PQ ) {

    pQueue *p;

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
 *  Finds the k Nearest Neighbors (kNN) of a query point q given that the set
 *  of data is stored as an Octree. All distances are in the normalized units
 *  (i.e. scaled from [0.0-1.0] ).
 *
 *    \param[in]     q           Query position. I.e. the point we want to find NNs for.
 *    \param[in]     Root       Root node of Octree.
 *    \param[in]     K          Number of NNs to find.
 *    \param[in]     MaxDist2   Threshold distance^2 beyond which we give up on finding
 *                              NNs.  (i.e. could find them, but we arent interested
 *                              because they'd be too far from our query point to be
 *                              useful).
 *    \param[out]    Kgot       Number of NNs (within MaxDist2) that we actually found.
 *    \param[out]    kNN        List of kNN Data items. Sorted (closest first).
 *
 *    \returns       OCTREE_KNN_SUCCESS         Search succeeded.
 *                   OCTREE_KNN_TOO_FEW_NNS     Search terminated because we couldnt find K NNs that were close enough.
 *                   OCTREE_KNN_NOT_ENOUGH_DATA Octree doesnt contain enough data points.
 *
 *    \author        Mike Henderson
 *    \date          2009-2012
 *
 */
int Lgm_Octree_kNN( Lgm_Vector *q_in, Lgm_Octree *Octree, int K, int *Kgot, double MaxDist2, Lgm_OctreeData *kNN ) {

    int             k, done;
    pQueue          *PQ, *p;
    Lgm_OctreeCell  *Root;
    Lgm_Vector      q;


    Root = Octree->Root;
    *Kgot = 0;


//MaxDist2 = .05;

    /*
     * Check to see if there are enough points.
     * Bailout with error flag -1 if not.
     */
    if (Root->nDataBelow < K ) return( OCTREE_KNN_NOT_ENOUGH_DATA );


    /*
     *  Add Root Node to the Priority Queue.
     */
    PQ = NULL;
    Lgm_OctreeScalePosition( q_in, &q, Octree );
    InsertCell( Root, &q, &PQ, MaxDist2 );

    /*
     *  Process items on the Priority Queue until we are done.
     */
    k    = 0;       // havent found any NNs yet.
    done = FALSE;
    while( !done ) {

	    //PrintPQ( &PQ ); //only for debugging

        /*
         * Pop the highest priority item off of the Priority Queue.
         * This is the closest object to the query point. Note that
         * it could be a cell or it could be a point.
         */
        p = PopObj( &PQ );

        if ( p ) {

            if ( p->IsPoint ) {
                /*
                 * Since a point is now the closest object, it must be one of
                 * the kNN.  But, if its too far away, then we are done with
                 * the search -- we couldnt find K NNs close enough to the
                 * query point
                 */
                if ( p->MinDist2 > MaxDist2 ) {
                    while ( (p = PopObj(&PQ)) ) free( p );
                    return( OCTREE_KNN_TOO_FEW_NNS );
                }

                /*
                 * Otherwise, add this point as the next NN and continue.
                 */
                kNN[ k   ]       = p->Obj->Data[p->j];
                kNN[ k++ ].Dist2 = p->MinDist2; // save the dist2 into the data struct
                *Kgot = k;

            } else {

                /*
                 * If the object is cell, then descend one level closer to
                 * the leaf node that is closest to the query point. This also
                 * adds all cells we encounter along the way to the PQ. Ignore
                 * Any cell that is farther than MaxDist2 away from query point.
                 */
                DescendTowardClosestLeaf( p->Obj, &PQ, &q, MaxDist2 );

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

    }

    /*
     * Free all remaining objects on the PQ
     */
    while ( (p = PopObj(&PQ)) ) free( p );


    /*
     *  return success
     */
    return( OCTREE_KNN_SUCCESS );

}







/**
 *  Allocate 8 new octants in a parent cell.
 *
 *    \param[in]     Parent     Cell in which new octants are to be created.
 *
 *    \returns       pointer to array of 8 newly allocated octants.
 *
 *    \author        Mike Henderson
 *    \date          2009-2012
 *
 */
Lgm_OctreeCell *CreateNewOctants( Lgm_OctreeCell *Parent ) {

    int             i;
    unsigned int    xMask, yMask, zMask, BranchBit;
    unsigned int    xParentCode, yParentCode, zParentCode;
    unsigned int    xRight, yRight, zRight, Level;
    char            Str[20];
    double          xMin, yMin, zMin;
    Lgm_OctreeCell  *Cell;

    Level     = Parent->Level - 1;
    BranchBit = 1<<Level;

    /*
     * allocate 8 cells at once
     */
    Cell = (Lgm_OctreeCell *) calloc( 8, sizeof( Lgm_OctreeCell ) );

    /*
     * For each cell, set nData = 0
     * Set Parent
     * Set Octant to NULL
     * And determine the x,y,z Location Codes.
     */
    xMask = 1; yMask = 1<<1; zMask = 1<<2;
    xParentCode = Parent->xLocationCode;
    yParentCode = Parent->yLocationCode;
    zParentCode = Parent->zLocationCode;
//    printf("Parent: Center, h = %g %g %g  %g\n", Parent->Center.x, Parent->Center.y, Parent->Center.z, Parent->h);
    for (i=0; i<8; i++){
        Cell[i].nData = 0;
        Cell[i].Parent = Parent;
        Cell[i].Octant = NULL;

        /*
         * determine if this octant is left or right (i.e. in each coord)
         */
        xRight = i&xMask;
        yRight = (i&yMask)>>1;
        zRight = (i&zMask)>>2;

        Cell[i].xLocationCode = ( xRight ) ? xParentCode + BranchBit : xParentCode;
        Cell[i].yLocationCode = ( yRight ) ? yParentCode + BranchBit : yParentCode;
        Cell[i].zLocationCode = ( zRight ) ? zParentCode + BranchBit : zParentCode;
        Cell[i].Level = Level;

        /*
         * Determine center of cell and width/2 (h=half face width)
         */
        Cell[i].h = 0.5*Parent->h;
        Cell[i].Center.x = (double)Cell[i].xLocationCode/(double)OCTREE_MAX_VAL + Cell[i].h;
        Cell[i].Center.y = (double)Cell[i].yLocationCode/(double)OCTREE_MAX_VAL + Cell[i].h;
        Cell[i].Center.z = (double)Cell[i].zLocationCode/(double)OCTREE_MAX_VAL + Cell[i].h;
//    	printf("Child: Center, h = %g %g %g  %g\n", Cell[i].Center.x, Cell[i].Center.y, Cell[i].Center.z, Cell[i].h);



//        printf("Cell[%d].Center = %g %g %g    h = %g\n", i, Cell[i].Center.x, Cell[i].Center.y, Cell[i].Center.z, Cell[i].h);
//        Binary( Cell[i].xLocationCode, Str ); printf("Index = %d  xLocationCode = %s\n", i, Str);
//        Binary( Cell[i].yLocationCode, Str ); printf("Index = %d  yLocationCode = %s\n", i, Str);
//        Binary( Cell[i].zLocationCode, Str ); printf("Index = %d  zLocationCode = %s\n\n", i, Str);

    }

    return( Cell );

}


/**
 *  Recursively subdivide volume. Stop when the number of data points ber cell
 *  is lower than the threshold given by OCTREE_MAX_DATA_PER_OCTANT.
 *
 *    \param[in]     Node in the octree to subdivide.
 *
 *    \returns       void
 *
 *    \author        Mike Henderson
 *    \date          2009-2012
 *
 */
void SubDivideVolume( Lgm_OctreeCell *Vol ) {

    int                 Level;
    unsigned int        i, j, BranchBit, xLocationCode, yLocationCode, zLocationCode;
    unsigned long int   ii;
    unsigned int        *Octant;
    char                Str[20];


    /*
     *  Create 8 new Octant cells
     *  CreateNewOctant() should make nData = 0 as initial value.
     */
    Vol->Octant  = CreateNewOctants( Vol );

    /*
     * Loop through the Data and figure out which
     * octant they belong to. Do this first so we know
     * how much memory we need to allocate to each child's
     * data field.
     */
    Octant = (unsigned int *) calloc( Vol->nData, sizeof( unsigned int ) );
    Level = Vol->Level-1;
    BranchBit = 1<<Level;
//printf("Level = %d\n", Level);
//printf("BranchBit = %d\n", BranchBit);
//Binary( BranchBit, Str );
//printf("BranchBit = %s\n", Str);
    for ( ii=0; ii<Vol->nData; ii++ ){

        // code to find what octant the object is in. Call it j
        xLocationCode = (unsigned int)(Vol->Data[ii].Position.x*OCTREE_MAX_VAL);
        yLocationCode = (unsigned int)(Vol->Data[ii].Position.y*OCTREE_MAX_VAL);
        zLocationCode = (unsigned int)(Vol->Data[ii].Position.z*OCTREE_MAX_VAL);
//printf("Vol->Data[ii].Position.x, Vol->Data[ii].Position.y, Vol->Data[ii].Position.z = %g %g %g\n", Vol->Data[ii].Position.x, Vol->Data[ii].Position.y, Vol->Data[ii].Position.z);
        j = (((xLocationCode & BranchBit) >> Level) + 2*( ((yLocationCode & BranchBit) >> Level) + 2*((zLocationCode & BranchBit) >> Level) ) );
        Octant[ii] = j;
        ++(Vol->Octant[j].nData);
//printf("%4d: j = %d\n", ii, j);

    }

    /*
     *  Alloc mem to hold the data in each child cell
     */
    for ( i=0; i<8; i++ ) {
//        printf("nData = %d\n", Vol->Octant[i].nData);
        Vol->Octant[i].Data = calloc( Vol->Octant[i].nData, sizeof( Lgm_OctreeData ) );
        Vol->Octant[i].nData = 0; // reset to zero (will recount below)
    }

    /*
     *  Copy the data into the children cells
     */
    for ( ii=0; ii<Vol->nData; ii++ ){
        j = Octant[ii];
        Vol->Octant[j].Data[ (Vol->Octant[j].nData)++ ] = Vol->Data[ii];

    }


    /*
     * Zero the Data count in this cell.
     * Free Octant -- not needed anymore.
     * Free memory allocated to Parent data field -- its not a leaf anymore.
     */
    Vol->nDataBelow = Vol->nData;
    Vol->nData = 0;
    free( Octant ); Octant = NULL;
    free( Vol->Data ); Vol->Data = NULL;



    /*
     *  Subdivide if there are too many objects in an octant cell
     */
    for ( i=0; i<8; i++ ) {
        if ( (Vol->Octant[i].Level > 0 ) && (Vol->Octant[i].nData > OCTREE_MAX_DATA_PER_OCTANT ) ) SubDivideVolume( &(Vol->Octant[i]) );
    }


    return;


}

/**
 *  Scale an input vector position. E.g., \f$ v = (u-Min)/ (Max-Min)\f$.
 *
 *    \param[in]     u       input vector position.
 *    \param[out]    v       output position in normalized coordinates.
 *    \param[in]     Octree  The octree structure (which contains the scaling values).
 *
 *    \returns       void
 *
 *    \author        Mike Henderson
 *    \date          2009-2012
 *
 */
void Lgm_OctreeScalePosition( Lgm_Vector *u, Lgm_Vector *v, Lgm_Octree *Octree ) {
    v->x = (u->x - Octree->Min)/Octree->Diff;
    v->y = (u->y - Octree->Min)/Octree->Diff;
    v->z = (u->z - Octree->Min)/Octree->Diff;
}

/**
 *  Un-scale an input vector position. E.g., \f$ u = v*Diff + Min\f$.
 *
 *    \param[in]     v       input position in normalized coordinates.
 *    \param[out]    u       output position in original coordinates.
 *    \param[in]     Octree  The octree structure (which contains the scaling values).
 *
 *    \returns       void
 *
 *    \author        Mike Henderson
 *    \date          2009-2012
 *
 */
void Lgm_OctreeUnScalePosition( Lgm_Vector *v, Lgm_Vector *u, Lgm_Octree *Octree ) {
    u->x = v->x*Octree->Diff + Octree->Min;
    u->y = v->y*Octree->Diff + Octree->Min;
    u->z = v->z*Octree->Diff + Octree->Min;
}


/**
 *  Scale a distance. E.g., \f$ u = v/Diff \f$.
 *
 *    \param[in]     u       input distance (in original coords).
 *    \param[out]    v       output position in scaled coordinates.
 *    \param[in]     Octree  The octree structure (which contains the scaling values).
 *
 *    \returns       void
 *
 *    \author        Mike Henderson
 *    \date          2009-2012
 *
 */
void Lgm_OctreeScaleDistance( double u, double *v, Lgm_Octree *Octree ) {
    *v = u/Octree->Diff;
}

/**
 *  Un-scale a distance. E.g., \f$ u = v*Diff \f$.
 *
 *    \param[in]     u       input distance (in scaled coords).
 *    \param[out]    v       output position in original coordinates.
 *    \param[in]     Octree  The octree structure (which contains the scaling values).
 *
 *    \returns       void
 *
 *    \author        Mike Henderson
 *    \date          2009-2012
 *
 */
void Lgm_OctreeUnScaleDistance( double v, double *u, Lgm_Octree *Octree ) {
    *u = v*Octree->Diff;
}



