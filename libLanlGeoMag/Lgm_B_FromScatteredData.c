/*! \file Lgm_B_FromScatteredData.c
 *  \brief Routines to compute B-field from a sattered set of data (for example, B defined on meshes).
 *
 *
 *
 */

#include "Lgm/Lgm_MagModelInfo.h"
#include "Lgm/Lgm_Octree.h"
#include "Lgm/Lgm_DFI_RBF.h"
#include "Lgm/qsort.h"

#define int_lt(a,b) ((*a)<(*b))



/**
 *      \param[in]         v  - array of position vectors
 *      \param[in]         B  - array of B-field vectors at the corresponding v's
 *      \param[in,out]  Info  - Pointer to Lgm_MagModelInfo structure.
 *
 *      \return Always returns 1. Fix this...?
 *
 *      \author  M. G. Henderson
 *      \date    January 25, 2012
 *
 *
 */
int Lgm_B_FromScatteredData( Lgm_Vector *v, Lgm_Vector *B, Lgm_MagModelInfo *Info ) {

    int               K, Kgot, n_data, i;
    double            eps, d;
    Lgm_OctreeData   *kNN;
    Lgm_DFI_RBF_Info *rbf;
    Lgm_Vector       *v_data, *B_data;


    /*
     *  Make sure Octree has been initialized
     */
//    if ( Info->Octree->Root == NULL ) return( OCTREE_IS_NULL );



    /*
     *  Allocate space for the K Nearest Neighbors.
     */
    K = Info->Octree_kNN_k;
    LGM_ARRAY_1D( kNN, K, Lgm_OctreeData );


    /*
     *  Find the K Nearest Neighbors.
     */
    Lgm_Octree_kNN( v, Info->Octree, K, &Kgot, 1.0*1.0, kNN );



    /*
     *   repack data into arrays.
     */
    n_data = Kgot;
    LGM_ARRAY_1D( v_data, n_data, Lgm_Vector );
    LGM_ARRAY_1D( B_data, n_data, Lgm_Vector );
    for ( i=0; i<n_data; i++){
        Lgm_OctreeUnScalePosition( &(kNN[i].Position), &v_data[i], Info->Octree );
        B_data[i] = kNN[i].B;
    }




    /*
     *  Initialize the  Divergence Free Interp. RBF stuff
     */
    //Lgm_OctreeUnScaleDistance( 1.0, &d, Info->Octree );
    //eps = 1.0/d;
    eps = 0.01;

//eps = 0.1;
    //rbf = Lgm_DFI_RBF_Init( v_data, B_data, n_data, eps, LGM_RBF_GAUSSIAN );
    //rbf = Lgm_DFI_RBF_Init( v_data, B_data, n_data, eps, LGM_RBF_MULTIQUADRIC );


    /*
     *  Evaluate Divergence Free Interpolation
     */
    Lgm_DFI_RBF_Eval( v, B, rbf );




    /*
     *  Cleanup. Free rbf, kNN, etc..
     */
    Lgm_DFI_RBF_Free( rbf );
    LGM_ARRAY_1D_FREE( kNN );
    LGM_ARRAY_1D_FREE( v_data );
    LGM_ARRAY_1D_FREE( B_data );



    return( 1 );

}

/*
 *  Setup the hash table used in Lgm_B_FromScatteredData().
 */
void Lgm_B_FromScatteredData_SetUp( Lgm_MagModelInfo *Info ) {

    if ( Info->rbf_ht_alloced ) Lgm_B_FromScatteredData_TearDown( Info );
    Info->rbf_ht         = NULL;
    Info->rbf_ht_alloced = FALSE;
    Info->RBF_nHashFinds = 0;
    Info->RBF_nHashAdds  = 0;
}


/*
 *  Iterates over all the entries in the hash table and 1) deletes them from
 *  the hash table, then 2) free the structure itself.
 */
void Lgm_B_FromScatteredData_TearDown( Lgm_MagModelInfo *Info ) {

    Lgm_DFI_RBF_Info *rbf, *rbf_tmp;

    HASH_ITER( hh, Info->rbf_ht, rbf, rbf_tmp ) {
        HASH_DELETE( hh, Info->rbf_ht, rbf );
        Lgm_DFI_RBF_Free( rbf );
    }

    if ( Info->Octree_kNN_Alloced > 0 ) {
        LGM_ARRAY_1D_FREE( Info->Octree_kNN );                                                                                                                                                               
    }

}


/** This routine is the same as Lgm_B_FromScatteredData(), except that here we
 *  assume the B's have had the dipole substracted already. The final result
 *  will have a dipole added back in.
 *
 *      \param[in]         v  - array of position vectors
 *      \param[in]         B  - array of B-field vectors at the corresponding v's
 *      \param[in,out]  Info  - Pointer to Lgm_MagModelInfo structure.
 *
 *      \return Always returns 1. Fix this...?
 *
 *      \author  M. G. Henderson
 *      \date    January 25, 2012
 *
 *
 */
int Lgm_B_FromScatteredData2( Lgm_Vector *v, Lgm_Vector *B, Lgm_MagModelInfo *Info ) {

    int                 K, Kgot, n_data, i;
    double              eps, d;
    Lgm_OctreeData     *kNN;
    Lgm_DFI_RBF_Info   *rbf;                      // single structure.
    Lgm_Vector         *v_data, *B_data, B1, B2;
    unsigned long int  *I_data;
    unsigned long int  *LookUpKey;                // key comprised of an array of unsigned long int Id's
    int                 KeyLength;                // length (in bytes) of LookUpKey


    /*
     *  Make sure Octree has been initialized
     */
//    if ( Info->Octree->Root == NULL ) return( OCTREE_IS_NULL );


if ( Lgm_Magnitude( v ) > 1.5 ) {

    /*
     *  Allocate space for the K Nearest Neighbors.
     *  This is probably a bit wasteful...
     *  Should cache this array.
     */
    K = Info->Octree_kNN_k;
    if (K > 2) {

        if ( Info->Octree_kNN_Alloced == 0 ) {

            LGM_ARRAY_1D( Info->Octree_kNN, K, Lgm_OctreeData );
            Info->Octree_kNN_Alloced = K;

        } else if ( K != Info->Octree_kNN_Alloced ) {

            /*
             * kNN is allocated but K has changed. Realloc.
             */
            LGM_ARRAY_1D_FREE( Info->Octree_kNN );
            LGM_ARRAY_1D( Info->Octree_kNN, K, Lgm_OctreeData );

        }
        Info->Octree_kNN_Alloced = K;

    } else {

        printf("Lgm_B_FromScatteredData2(): Error. Not enough nearest neighbors specified: Info->Octree_kNN_k = %d\n", Info->Octree_kNN_k);
        exit(-1);

    }






    /*
     *  Find the K Nearest Neighbors.
     */
    Lgm_Octree_kNN( v, Info->Octree, K, &Kgot, Info->Octree_kNN_MaxDist2, Info->Octree_kNN );



    /*
     *  From the K nearest neighbors, construct a key for the hash-table.
     *  Probably should sort them so that different permutations of the same k
     *  NN's will be identified as the same set. But lets worry about that
     *  later.
     *
     */
    LGM_ARRAY_1D( LookUpKey, Kgot, unsigned long int );
    for ( i=0; i<Kgot; i++ ) LookUpKey[i] = Info->Octree_kNN[i].Id;
    KeyLength = Kgot*sizeof( unsigned long int );
    QSORT( unsigned long int, LookUpKey, Kgot, int_lt );
    //quicksort_uli( (long int)Kgot, LookUpKey-1 );



    /*
     *  Look up the key in the hash-table to see if its already there.
     *  If it exists, we bypass refitting the RBF weights.
     */
//printf("Searching for: " );
//for(i=0;i<Kgot; i++) printf(" %ld ", LookUpKey[i] );
//printf("    (KeyLength = %d)\n", KeyLength);
//printf("Searching for: " );
//for(i=0;i<Kgot; i++) printf(" %ld ", LookUpKey[i] );
//printf("    (KeyLength = %d)\n\n", KeyLength);
    HASH_FIND( hh, Info->rbf_ht, LookUpKey, KeyLength, rbf );
    ++(Info->RBF_nHashFinds);




    /*
     * If key didnt exist in hash-table, we need to compute RBF weeights, package up info
     * into a structure and add it to the hash table.
     */
    if ( rbf == NULL ) {

        /*
         *   repack data into arrays. This is wasteful also.
         */
        n_data = Kgot;
        LGM_ARRAY_1D( I_data, n_data, unsigned long int );
        LGM_ARRAY_1D( v_data, n_data, Lgm_Vector );
        LGM_ARRAY_1D( B_data, n_data, Lgm_Vector );
        for ( i=0; i<n_data; i++){
            Lgm_OctreeUnScalePosition( &(Info->Octree_kNN[i].Position), &v_data[i], Info->Octree );
            B_data[i] = Info->Octree_kNN[i].B;
            I_data[i] = Info->Octree_kNN[i].Id;
            //printf("Position[%d] = %g %g %g    B_data[%d] = %g %g %g\n", i, v_data[i].x, v_data[i].y, v_data[i].z, i, B_data[i].x, B_data[i].y, B_data[i].z );
        }
        QSORT( unsigned long int, I_data, n_data, int_lt );
        //quicksort_uli( (long int)n_data, I_data-1 );


        /*
         *  Construct the rbf structure. We dont free these until the hash
         *  table is done with. Note that the hash table will be nthe only
         *  reference to the pointer.  To free, use
         *  Lgm_B_FromScatteredData_TearDown().
         */
        eps = 0.01;
        rbf = Lgm_DFI_RBF_Init( I_data, v_data, B_data, n_data, eps, LGM_RBF_GAUSSIAN );

        LGM_ARRAY_1D_FREE( I_data );
        LGM_ARRAY_1D_FREE( v_data );
        LGM_ARRAY_1D_FREE( B_data );

//        printf("Adding item to hash table\n");
        HASH_ADD_KEYPTR( hh, Info->rbf_ht, rbf->LookUpKey, KeyLength, rbf );
        ++(Info->RBF_nHashAdds);

    } else {

        //printf("got one\n");

    }




    /*
     *  Evaluate Divergence Free Interpolation
     */
    Lgm_DFI_RBF_Eval( v, &B1, rbf );


    /*
     *  Cleanup. Free rbf, kNN, etc..
     */
    LGM_ARRAY_1D_FREE( LookUpKey );

} else {
    B1.x = B1.y = B1.z = 0.0;
}



    switch ( Info->InternalModel ){

        case LGM_CDIP:
                        Lgm_B_cdip( v, &B2, Info );
                        break;
        case LGM_EDIP:
                        Lgm_B_edip( v, &B2, Info );
                        break;
        case LGM_IGRF:
                        Lgm_B_igrf( v, &B2, Info );
                        break;
        default:
                        fprintf(stderr, "Lgm_B_FromScatteredData2: Unknown internal model (%d)\n", Info->InternalModel );
                        break;

    }


    B->x = B1.x + B2.x;
    B->y = B1.y + B2.y;
    B->z = B1.z + B2.z;


    return( 1 );

}
