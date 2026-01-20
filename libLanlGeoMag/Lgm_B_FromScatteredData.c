/*! \file Lgm_B_FromScatteredData.c
 *  \brief Routines to compute B-field from a sattered set of data (for example, B defined on meshes).
 *
 *
 *
 */

#include "Lgm/Lgm_KdTree.h"
#include "Lgm/Lgm_MagModelInfo.h"
#include "Lgm/Lgm_Octree.h"
#include "Lgm/Lgm_RBF.h"
#include "Lgm/qsort.h"

#define int_lt(a,b) ((*a)<(*b))
void SimplifyPointCloud( Lgm_Vector *v_data, Lgm_Vector *B_data, Lgm_Vector *E_data, int n_data, Lgm_Vector *v_data2, Lgm_Vector *B_data2, Lgm_Vector *E_data2, int *n_data2, double s2min_min, double R, double nR );
double EstimateGridRes( Lgm_Vector *v );




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
 *  Set the radial basis func to use  (and its eps value)
 */
void Lgm_B_FromScatteredData_SetRbf( Lgm_MagModelInfo *Info, double eps, int RbfType ) {

    Info->RBF_Eps  = eps;
    Info->RBF_Type = RbfType;

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
        Info->Octree_kNN_Alloced = 0;
    }
    if ( Info->KdTree_kNN_Alloced > 0 ) {
        LGM_ARRAY_1D_FREE( Info->KdTree_kNN );
        Info->KdTree_kNN_Alloced = 0;
    }

}


/*
 *  Iterates over all the entries in the hash table and 1) deletes them from
 *  the hash table, then 2) free the structure itself.
 *
 *  Really should unify the rbf structures..
 */
void Lgm_B_FromScatteredData4_TearDown( Lgm_MagModelInfo *Info ) {

    Lgm_Vec_RBF_Info *rbf, *rbf_tmp;

    HASH_ITER( hh, Info->vec_rbf_ht, rbf, rbf_tmp ) {
        HASH_DELETE( hh, Info->vec_rbf_ht, rbf );
        Lgm_Vec_RBF_Free( rbf );
    }

    if ( Info->Octree_kNN_Alloced > 0 ) {
        LGM_ARRAY_1D_FREE( Info->Octree_kNN );
        Info->Octree_kNN_Alloced = 0;
    }
    if ( Info->KdTree_kNN_Alloced > 0 ) {
        LGM_ARRAY_1D_FREE( Info->KdTree_kNN );
        Info->KdTree_kNN_Alloced = 0;
    }

}



/*
 *  Setup the hash table used in Lgm_B_FromScatteredData().
 */
void Lgm_B_FromScatteredData5_SetUp( Lgm_MagModelInfo *Info ) {


Lgm_Vec_RBF_Info ***b;

    if ( Info->rbf_ht_alloced ) Lgm_B_FromScatteredData5_TearDown( Info );
    Info->vec_rbf_ht     = NULL;
    Info->vec_rbf_e_ht   = NULL;
    Info->rbf_ht_alloced = FALSE;
    Info->RBF_nHashFinds = 0;
    Info->RBF_nHashAdds  = 0;


    Info->RBF_CB.n           = 0;
    Info->RBF_CB.nEntries    =      0; // Initial entries
    Info->RBF_CB.N           = 200000; // Max entries
    Info->RBF_CB.Buf1        = (Lgm_Vec_RBF_Info **)calloc( Info->RBF_CB.N, sizeof( Lgm_Vec_RBF_Info *) );
    Info->RBF_CB.Buf2        = (Lgm_Vec_RBF_Info **)calloc( Info->RBF_CB.N, sizeof( Lgm_Vec_RBF_Info *) );

    Info->vec_rbf_ht_size    =    0.0; // Initial size in MB
    Info->vec_rbf_ht_maxsize = 2000.0; // Max size in MB
//    Info->vec_rbf_ht_maxsize = 200.0; // Max size in MB
    Info->RBF_CB.oldest_i    = 0;
    Info->RBF_CB.newest_i    = -1;
    

}
/*
 *  Iterates over all the entries in the hash table and 1) deletes them from
 *  the hash table, then 2) free the structure itself.
 *
 *  Really should unify the rbf structures..
 */
void Lgm_B_FromScatteredData5_TearDown( Lgm_MagModelInfo *Info ) {

    Lgm_Vec_RBF_Info *rbf, *rbf_tmp;

    // free the B-field RBFs
    HASH_ITER( hh, Info->vec_rbf_ht, rbf, rbf_tmp ) {
        HASH_DELETE( hh, Info->vec_rbf_ht, rbf );
        Lgm_Vec_RBF_Free( rbf );
    }

    // free the E-field RBFs
    HASH_ITER( hh, Info->vec_rbf_e_ht, rbf, rbf_tmp ) {
        HASH_DELETE( hh, Info->vec_rbf_e_ht, rbf );
        Lgm_Vec_RBF_Free( rbf );
    }

    if ( Info->Octree_kNN_Alloced > 0 ) {
        LGM_ARRAY_1D_FREE( Info->Octree_kNN );
        Info->Octree_kNN_Alloced = 0;
    }
    if ( Info->KdTree_kNN_Alloced > 0 ) {
        LGM_ARRAY_1D_FREE( Info->KdTree_kNN );
        Info->KdTree_kNN_Alloced = 0;
    }

    Info->RBF_CB.n           = 0;
    Info->RBF_CB.nEntries    =      0; // Initial entries
    Info->vec_rbf_ht_size    =    0.0; // Initial size in MB
    free( Info->RBF_CB.Buf1 );
    free( Info->RBF_CB.Buf2 );
    Info->RBF_CB.oldest_i    = 0;
    Info->RBF_CB.newest_i    = -1;

}




/** This routine is the same as Lgm_B_FromScatteredData(), except that here we
 * ssume the B's have had the dipole substracted already. The final result
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
         * If key didnt exist in hash-table, we need to compute RBF weights, package up info
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
             *  table is done with. Note that the hash table will be the only
             *  reference to the pointer.  To free, use
             *  Lgm_B_FromScatteredData_TearDown().
             */
    eps = 0.01;
    eps = 0.1;
            rbf = Lgm_DFI_RBF_Init( I_data, v_data, B_data, n_data, eps, LGM_RBF_GAUSSIAN );

            LGM_ARRAY_1D_FREE( I_data );
            LGM_ARRAY_1D_FREE( v_data );
            LGM_ARRAY_1D_FREE( B_data );

            //printf("Adding item to hash table\n");
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
         *  Evaluate derivatives of Divergence Free Interpolation 
         *  put in the Info structure.
         */
        Lgm_DFI_RBF_Derivs_Eval( v, &Info->RBF_dBdx, &Info->RBF_dBdy, &Info->RBF_dBdz, rbf );


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




/** This routine is the same as Lgm_B_FromScatteredData2(), except that here we
 *  use kdtree instead of octree.
 *
 *      \param[in]         v  - array of position vectors
 *      \param[in]         B  - array of B-field vectors at the corresponding v's
 *      \param[in,out]  Info  - Pointer to Lgm_MagModelInfo structure.
 *
 *      \return Always returns 1. Fix this...?
 *
 *      \author  M. G. Henderson
 *      \date    July 2, 2015
 *
 *
 */
int Lgm_B_FromScatteredData3( Lgm_Vector *v, Lgm_Vector *B, Lgm_MagModelInfo *Info ) {

    int                 K, Kgot, n_data, i;
    double             *eps, d;
    Lgm_KdTreeData     *kNN;
    Lgm_DFI_RBF_Info   *rbf;                      // single structure.
    Lgm_Vector         *v_data, *B_data, B1, B2;
    Lgm_Vector          b, Grad_B1, GradB_dipole; 
    unsigned long int  *I_data;
    unsigned long int  *LookUpKey;                // key comprised of an array of unsigned long int Id's
    int                 KeyLength;                // length (in bytes) of LookUpKey
    int                 (*Dipole)(Lgm_Vector*, Lgm_Vector*, Lgm_MagModelInfo*);         // tmp Pointer to Bfield function


    /*
     *  Make sure KdTree has been initialized
     */
    if ( Lgm_Magnitude( v ) > 1.5 ) {

        /*
         *  Allocate space for the K Nearest Neighbors.
         *  This is probably a bit wasteful...
         *  Should cache this array.
         */
        K = Info->KdTree_kNN_k;
        if (K > 2) {

            if ( Info->KdTree_kNN_Alloced == 0 ) {

                LGM_ARRAY_1D( Info->KdTree_kNN, K, Lgm_KdTreeData );
                Info->KdTree_kNN_Alloced = K;

            } else if ( K != Info->KdTree_kNN_Alloced ) {

                /*
                 * kNN is allocated but K has changed. Realloc.
                 */
                LGM_ARRAY_1D_FREE( Info->KdTree_kNN );
                LGM_ARRAY_1D( Info->KdTree_kNN, K, Lgm_KdTreeData );

            }
            Info->KdTree_kNN_Alloced = K;

        } else {

            printf("Lgm_B_FromScatteredData3(): Error. Not enough nearest neighbors specified: Info->KdTree_kNN_k = %d\n", Info->KdTree_kNN_k);
            exit(-1);

        }






        /*
         *  Find the K Nearest Neighbors.
         */
        double q[3];
        q[0] = v->x; q[1] = v->y; q[2] = v->z;
        Lgm_KdTree_kNN( q, 3, Info->KdTree, K, &Kgot, Info->KdTree_kNN_MaxDist2, Info->KdTree_kNN );
        //printf("K, Kgot = %d %d    q = %g %g %g  Info->KdTree_kNN_MaxDist2 = %g \n", K, Kgot, q[0], q[1], q[2], Info->KdTree_kNN_MaxDist2 );



        // not needed? Put under verbosity setting?
        for ( i=0; i<Kgot; i++ ) {
            if ( Info->KdTree_kNN[i].Dist2 > Info->KdTree_kNN_MaxDist2){
                printf("ERROR:Info->KdTree_kNN[i].Dist2 = %g\n", Info->KdTree_kNN[i].Dist2);
            }
        }


        /*
         *  From the K nearest neighbors, construct a key for the hash-table.
         *  Probably should sort them so that different permutations of the same k
         *  NN's will be identified as the same set. But lets worry about that
         *  later.
         *
         */
        LGM_ARRAY_1D( LookUpKey, Kgot, unsigned long int );
        for ( i=0; i<Kgot; i++ ) LookUpKey[i] = Info->KdTree_kNN[i].Id;
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
        rbf = NULL;
        HASH_FIND( hh, Info->rbf_ht, LookUpKey, KeyLength, rbf );
        ++(Info->RBF_nHashFinds);




        /*
         * If key didnt exist in hash-table, we need to compute RBF weights, package up info
         * into a structure and add it to the hash table.
         */
        if ( rbf == NULL ) {

            //printf("Did not find key - computing new rbf coeffs\n\n");

            /*
             *   repack data into arrays. This is wasteful also.
             */
            n_data = Kgot;
            LGM_ARRAY_1D( eps,    n_data, double );
            LGM_ARRAY_1D( I_data, n_data, unsigned long int );
            LGM_ARRAY_1D( v_data, n_data, Lgm_Vector );
            LGM_ARRAY_1D( B_data, n_data, Lgm_Vector );
            double *bbb;
            for ( i=0; i<n_data; i++){
                v_data[i].x = Info->KdTree_kNN[i].Position[0];
                v_data[i].y = Info->KdTree_kNN[i].Position[1];
                v_data[i].z = Info->KdTree_kNN[i].Position[2];

                bbb = (double *)Info->KdTree_kNN[i].Object;
                B_data[i].x = bbb[0];
                B_data[i].y = bbb[1];
                B_data[i].z = bbb[2];

                I_data[i] = Info->KdTree_kNN[i].Id;
                // eps not used yet for this version
                eps[i]    = Info->RBF_Eps;
            }
            QSORT( unsigned long int, I_data, n_data, int_lt );

/*
double dx, dy, dz, d2, d2min;
int j;
d2min = 1e6;
for ( i=0; i<n_data; i++){
  for ( j=0; i<n_data; i++){
    if (i!=j){

        dx = v_data[i].x - v_data[j].x;
        dy = v_data[i].x - v_data[j].x;
        dz = v_data[i].x - v_data[j].x;
        d2 = dx*dx + dy*dy + dz*dz;
        if ( (d2 > 0.0) && ( d2 < d2min) ) {
            d2min = d2;
        }
    
    }
  }
}
//printf("d2min = %g\n", d2min);
Info->RBF_Eps = 1.0/(d2min);
*/


            /*
             *  Construct the rbf structure. We dont free these until the hash
             *  table is done with. Note that the hash table will be the only
             *  reference to the pointer.  To free, use
             *  Lgm_B_FromScatteredData_TearDown().
             */
            rbf = Lgm_DFI_RBF_Init( I_data, v_data, B_data, n_data, Info->RBF_Eps, Info->RBF_Type ); 


            LGM_ARRAY_1D_FREE(    eps );
            LGM_ARRAY_1D_FREE( I_data );
            LGM_ARRAY_1D_FREE( v_data );
            LGM_ARRAY_1D_FREE( B_data );

            //printf("Adding item to hash table\n");
            HASH_ADD_KEYPTR( hh, Info->rbf_ht, rbf->LookUpKey, KeyLength, rbf );
            ++(Info->RBF_nHashAdds);

        } else {

            //printf("got one\n");
            //printf("        Found: " );
            //for(i=0;i<Kgot; i++) printf(" %ld ", rbf->LookUpKey[i] );
            //printf("\n\n");

        }




        /*
         *  Evaluate Divergence Free Interpolation
         */
        Lgm_DFI_RBF_Eval( v, &B1, rbf );
        //printf("Evaluating DFI_RBF with rbf = %p  at v = %g %g %g   B1 = %g %g %g\n\n\n", rbf, v->x, v->y, v->z, B1.x, B1.y, B1.z);

        /*
         *  Evaluate derivatives of Divergence Free Interpolation 
         *  put in the Info structure.
         */
        if ( Info->RBF_CompGradAndCurl ) {
            Lgm_DFI_RBF_Derivs_Eval( v, &Info->RBF_dBdx, &Info->RBF_dBdy, &Info->RBF_dBdz, rbf );
        }


        /*
         *  Cleanup. Free rbf, kNN, etc..
         */
        LGM_ARRAY_1D_FREE( LookUpKey );

    } else {

        B1.x = B1.y = B1.z = 0.0;

    }



    // Save Bfield, so we can temporaily swap it for an internal model (so we can compute GradB)
    switch ( Info->InternalModel ){

        case LGM_CDIP:
                        Lgm_B_cdip( v, &B2, Info );
                        Dipole = Lgm_B_cdip;
                        break;
        case LGM_EDIP:
                        Lgm_B_edip( v, &B2, Info );
                        Dipole = Lgm_B_edip;
                        break;
        case LGM_IGRF:
                        Lgm_B_igrf( v, &B2, Info );
                        Dipole = Lgm_B_igrf;
                        break;
        default:
                        fprintf(stderr, "Lgm_B_FromScatteredData3: Unknown internal model (%d)\n", Info->InternalModel );
                        break;

    }

    // Total B
    B->x = B1.x + B2.x;
    B->y = B1.y + B2.y;
    B->z = B1.z + B2.z;


    if ( Info->RBF_CompGradAndCurl ) {

        Lgm_Vector u, u0, Bvec;
        int DerivScheme, N;
        Lgm_Vector dBdx, dBdy, dBdz;
        double f1x[7], f2x[7], f3x[7];
        double f1y[7], f2y[7], f3y[7];
        double f1z[7], f2z[7], f3z[7];
        double H, Bmag;
        double h = 1e-3;

        u0 = *v;
        N = 1; DerivScheme = LGM_DERIV_TWO_POINT;
        
        for (i=-N; i<=N; ++i){
            if ( i != 0 ) {
                u = u0; H = (double)i*h; u.x += H;
                Dipole( &u, &Bvec, Info );
                f1x[i+N] = Bvec.x;
                f2x[i+N] = Bvec.y;
                f3x[i+N] = Bvec.z;
            }
        }
        for (i=-N; i<=N; ++i){
            if ( i != 0 ) {
                u = u0; H = (double)i*h; u.y += H;
                Dipole( &u, &Bvec, Info );
                f1y[i+N] = Bvec.x;
                f2y[i+N] = Bvec.y;
                f3y[i+N] = Bvec.z;
            }
        }
        for (i=-N; i<=N; ++i){
            if ( i != 0 ) {
                u = u0; H = (double)i*h; u.z += H;
                Dipole( &u, &Bvec, Info );
                Bmag = Lgm_Magnitude( &Bvec );
                f1z[i+N] = Bvec.x;
                f2z[i+N] = Bvec.y;
                f3z[i+N] = Bvec.z;
            }
        }
        if (DerivScheme == LGM_DERIV_SIX_POINT){

            dBdx.x = (f1x[6] - 9.0*f1x[5] + 45.0*f1x[4] - 45.0*f1x[2] + 9.0*f1x[1] - f1x[0])/(60.0*h);
            dBdx.y = (f2x[6] - 9.0*f2x[5] + 45.0*f2x[4] - 45.0*f2x[2] + 9.0*f2x[1] - f2x[0])/(60.0*h);
            dBdx.z = (f3x[6] - 9.0*f3x[5] + 45.0*f3x[4] - 45.0*f3x[2] + 9.0*f3x[1] - f3x[0])/(60.0*h);

            dBdy.x = (f1y[6] - 9.0*f1y[5] + 45.0*f1y[4] - 45.0*f1y[2] + 9.0*f1y[1] - f1y[0])/(60.0*h);
            dBdy.y = (f2y[6] - 9.0*f2y[5] + 45.0*f2y[4] - 45.0*f2y[2] + 9.0*f2y[1] - f2y[0])/(60.0*h);
            dBdy.z = (f3y[6] - 9.0*f3y[5] + 45.0*f3y[4] - 45.0*f3y[2] + 9.0*f3y[1] - f3y[0])/(60.0*h);

            dBdz.x = (f1z[6] - 9.0*f1z[5] + 45.0*f1z[4] - 45.0*f1z[2] + 9.0*f1z[1] - f1z[0])/(60.0*h);
            dBdz.y = (f2z[6] - 9.0*f2z[5] + 45.0*f2z[4] - 45.0*f2z[2] + 9.0*f2z[1] - f2z[0])/(60.0*h);
            dBdz.z = (f3z[6] - 9.0*f3z[5] + 45.0*f3z[4] - 45.0*f3z[2] + 9.0*f3z[1] - f3z[0])/(60.0*h);

        } else if (DerivScheme == LGM_DERIV_FOUR_POINT){

            dBdx.x = (-f1x[4] + 8.0*f1x[3] - 8.0*f1x[1] + f1x[0])/(12.0*h);
            dBdx.y = (-f2x[4] + 8.0*f2x[3] - 8.0*f2x[1] + f2x[0])/(12.0*h);
            dBdx.z = (-f3x[4] + 8.0*f3x[3] - 8.0*f3x[1] + f3x[0])/(12.0*h);

            dBdy.x = (-f1y[4] + 8.0*f1y[3] - 8.0*f1y[1] + f1y[0])/(12.0*h);
            dBdy.y = (-f2y[4] + 8.0*f2y[3] - 8.0*f2y[1] + f2y[0])/(12.0*h);
            dBdy.z = (-f3y[4] + 8.0*f3y[3] - 8.0*f3y[1] + f3y[0])/(12.0*h);

            dBdz.x = (-f1z[4] + 8.0*f1z[3] - 8.0*f1z[1] + f1z[0])/(12.0*h);
            dBdz.y = (-f2z[4] + 8.0*f2z[3] - 8.0*f2z[1] + f2z[0])/(12.0*h);
            dBdz.z = (-f3z[4] + 8.0*f3z[3] - 8.0*f3z[1] + f3z[0])/(12.0*h);

        } else if (DerivScheme == LGM_DERIV_TWO_POINT){

            dBdx.x = (f1x[2] - f1x[0])/(2.0*h);
            dBdx.y = (f2x[2] - f2x[0])/(2.0*h);
            dBdx.z = (f3x[2] - f3x[0])/(2.0*h);

            dBdy.x = (f1y[2] - f1y[0])/(2.0*h);
            dBdy.y = (f2y[2] - f2y[0])/(2.0*h);
            dBdy.z = (f3y[2] - f3y[0])/(2.0*h);

            dBdz.x = (f1z[2] - f1z[0])/(2.0*h);
            dBdz.y = (f2z[2] - f2z[0])/(2.0*h);
            dBdz.z = (f3z[2] - f3z[0])/(2.0*h);

        } else {
            printf("huh?\n");
            exit(0);
        }




        // Add in the contribution from external field.
        dBdx.x += Info->RBF_dBdx.x; 
        dBdx.y += Info->RBF_dBdx.y; 
        dBdx.z += Info->RBF_dBdx.z;

        dBdy.x += Info->RBF_dBdy.x; 
        dBdy.y += Info->RBF_dBdy.y; 
        dBdy.z += Info->RBF_dBdy.z;

        dBdz.x += Info->RBF_dBdz.x; 
        dBdz.y += Info->RBF_dBdz.y; 
        dBdz.z += Info->RBF_dBdz.z;
/*
dBdx.x = Info->RBF_dBdx.x; 
dBdx.y = Info->RBF_dBdx.y; 
dBdx.z = Info->RBF_dBdx.z;

dBdy.x = Info->RBF_dBdy.x; 
dBdy.y = Info->RBF_dBdy.y; 
dBdy.z = Info->RBF_dBdy.z;

dBdz.x = Info->RBF_dBdz.x; 
dBdz.y = Info->RBF_dBdz.y; 
dBdz.z = Info->RBF_dBdz.z;
*/

        // Compute final GradB
        b = *B; Lgm_NormalizeVector( &b );
        Info->RBF_Grad_B.x = Lgm_DotProduct( &b, &dBdx );
        Info->RBF_Grad_B.y = Lgm_DotProduct( &b, &dBdy );
        Info->RBF_Grad_B.z = Lgm_DotProduct( &b, &dBdz );

        // Compute Curl_B
        Info->RBF_Curl_B.x = Info->RBF_dBdy.z - Info->RBF_dBdz.y;
        Info->RBF_Curl_B.y = Info->RBF_dBdz.x - Info->RBF_dBdx.z;
        Info->RBF_Curl_B.z = Info->RBF_dBdx.y - Info->RBF_dBdy.x;
    }



    return( 1 );

}






/** This routine is the same as Lgm_B_FromScatteredData3(), except that here we
 *  use the non-DFI RBF funcs
 *
 *      \param[in]         v  - array of position vectors
 *      \param[in]         B  - array of B-field vectors at the corresponding v's
 *      \param[in,out]  Info  - Pointer to Lgm_MagModelInfo structure.
 *
 *      \return Always returns 1. Fix this...?
 *
 *      \author  M. G. Henderson
 *      \date    July 15, 2015
 *
 *
 */
int Lgm_B_FromScatteredData4( Lgm_Vector *v, Lgm_Vector *B, Lgm_MagModelInfo *Info ) {

    int                 K, Kgot, n_data, i;
    double              *eps, d;
    Lgm_KdTreeData     *kNN;
    Lgm_Vec_RBF_Info   *rbf;                      // single structure.
    Lgm_Vector         *v_data, *B_data, B1, B2;
    Lgm_Vector          b, Grad_B1, GradB_dipole; 
    unsigned long int  *I_data;
    unsigned long int  *LookUpKey;                // key comprised of an array of unsigned long int Id's
    int                 KeyLength;                // length (in bytes) of LookUpKey
    int                 (*Dipole)(Lgm_Vector*, Lgm_Vector*, Lgm_MagModelInfo*);         // tmp Pointer to Bfield function


    /*
     *  Make sure KdTree has been initialized
     */
    if ( Lgm_Magnitude( v ) > 1.5 ) {

        /*
         *  Allocate space for the K Nearest Neighbors.
         *  This is probably a bit wasteful...
         *  Should cache this array.
         */
        K = Info->KdTree_kNN_k;
        if (K > 2) {

            if ( Info->KdTree_kNN_Alloced == 0 ) {

                LGM_ARRAY_1D( Info->KdTree_kNN, K, Lgm_KdTreeData );
                Info->KdTree_kNN_Alloced = K;

            } else if ( K != Info->KdTree_kNN_Alloced ) {

                /*
                 * kNN is allocated but K has changed. Realloc.
                 */
                LGM_ARRAY_1D_FREE( Info->KdTree_kNN );
                LGM_ARRAY_1D( Info->KdTree_kNN, K, Lgm_KdTreeData );

            }
            Info->KdTree_kNN_Alloced = K;

        } else {

            printf("Lgm_B_FromScatteredData4(): Error. Not enough nearest neighbors specified: Info->KdTree_kNN_k = %d\n", Info->KdTree_kNN_k);
            exit(-1);

        }


        /*
         *  Find the K Nearest Neighbors.
         */
        double q[3];
        q[0] = v->x; q[1] = v->y; q[2] = v->z;
        Lgm_KdTree_kNN( q, 3, Info->KdTree, K, &Kgot, Info->KdTree_kNN_MaxDist2, Info->KdTree_kNN );


        // not needed? Put under verbosity setting?
        for ( i=0; i<Kgot; i++ ) {
            if ( Info->KdTree_kNN[i].Dist2 > Info->KdTree_kNN_MaxDist2){
                printf("Lgm_B_FromScatteredData4(): ERROR - Info->KdTree_kNN[i].Dist2 = %g\n", Info->KdTree_kNN[i].Dist2);
            }
        }


        /*
         *  From the K nearest neighbors, construct a key for the hash-table.
         *  Probably should sort them so that different permutations of the same k
         *  NN's will be identified as the same set. But lets worry about that
         *  later.
         *
         */
        LGM_ARRAY_1D( LookUpKey, Kgot, unsigned long int );
        for ( i=0; i<Kgot; i++ ) LookUpKey[i] = Info->KdTree_kNN[i].Id;
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
        rbf = NULL;
        HASH_FIND( hh, Info->vec_rbf_ht, LookUpKey, KeyLength, rbf );
        ++(Info->RBF_nHashFinds);


//if (0==1){
//}
//rbf = NULL;


        /*
         * If key didnt exist in hash-table, we need to compute RBF weights, package up info
         * into a structure and add it to the hash table.
         */
        if ( rbf == NULL ) {



            //printf("Did not find key - computing new rbf coeffs\n\n");

            /*
             *   repack data into arrays. This is wasteful also.
             */
            n_data = Kgot;
            LGM_ARRAY_1D( eps,    n_data, double );
            LGM_ARRAY_1D( I_data, n_data, unsigned long int );
            LGM_ARRAY_1D( v_data, n_data, Lgm_Vector );
            LGM_ARRAY_1D( B_data, n_data, Lgm_Vector );
            double *bbb;
            for ( i=0; i<n_data; i++){
                v_data[i].x = Info->KdTree_kNN[i].Position[0];
                v_data[i].y = Info->KdTree_kNN[i].Position[1];
                v_data[i].z = Info->KdTree_kNN[i].Position[2];

                bbb = (double *)Info->KdTree_kNN[i].Object;
                B_data[i].x = bbb[0];
                B_data[i].y = bbb[1];
                B_data[i].z = bbb[2];

                I_data[i] = Info->KdTree_kNN[i].Id;
                eps[i]    = Info->RBF_Eps;
            }
            QSORT( unsigned long int, I_data, n_data, int_lt );

/*
double dx, dy, dz, d2, d2min;
int j;
d2min = 1e6;
for ( i=0; i<n_data; i++){
  for ( j=0; i<n_data; i++){
    if (i!=j){

        dx = v_data[i].x - v_data[j].x;
        dy = v_data[i].x - v_data[j].x;
        dz = v_data[i].x - v_data[j].x;
        d2 = dx*dx + dy*dy + dz*dz;
        if ( (d2 > 0.0) && ( d2 < d2min) ) {
            d2min = d2;
        }
    
    }
  }
}
//printf("d2min = %g\n", d2min);
Info->RBF_Eps = 1.0/(d2min*10.0);
*/


            /*
             *  Construct the rbf structure. We dont free these until the hash
             *  table is done with. Note that the hash table will be the only
             *  reference to the pointer.  To free, use
             *  Lgm_B_FromScatteredData_TearDown().
             */
            rbf = Lgm_Vec_RBF_Init( I_data, v_data, B_data, eps, eps, eps, n_data, Info->RBF_DoPoly, Info->RBF_Type ); 


            LGM_ARRAY_1D_FREE(    eps );
            LGM_ARRAY_1D_FREE( I_data );
            LGM_ARRAY_1D_FREE( v_data );
            LGM_ARRAY_1D_FREE( B_data );

            //printf("Adding item to hash table\n");
            HASH_ADD_KEYPTR( hh, Info->vec_rbf_ht, rbf->LookUpKey, KeyLength, rbf );
            ++(Info->RBF_nHashAdds);

        } else {

            //printf("got one\n");
            //printf("        Found: " );
            //for(i=0;i<Kgot; i++) printf(" %ld ", rbf->LookUpKey[i] );
            //printf("\n\n");

        }




        /*
         *  Evaluate Divergence Free Interpolation
         */
        Lgm_Vec_RBF_Eval( v, &B1, rbf );
        //printf("Evaluating with rbf = %p  at v = %g %g %g   B1 = %g %g %g\n\n\n", rbf, v->x, v->y, v->z, B1.x, B1.y, B1.z);

        /*
         *  Evaluate derivatives of Divergence Free Interpolation 
         *  put in the Info structure.
         */
        if ( Info->RBF_CompGradAndCurl ) {
            Lgm_Vec_RBF_Derivs_Eval( v, &Info->RBF_dBdx, &Info->RBF_dBdy, &Info->RBF_dBdz, rbf );
        }


        /*
         *  Cleanup. Free rbf, kNN, etc..
         */
        LGM_ARRAY_1D_FREE( LookUpKey );

    } else {

        B1.x = B1.y = B1.z = 0.0;

    }



    // Save Bfield, so we can temporaily swap it for an internal model (so we can compute GradB)
    switch ( Info->InternalModel ){

        case LGM_CDIP:
                        Lgm_B_cdip( v, &B2, Info );
                        Dipole = Lgm_B_cdip;
                        break;
        case LGM_EDIP:
                        Lgm_B_edip( v, &B2, Info );
                        Dipole = Lgm_B_edip;
                        break;
        case LGM_IGRF:
                        Lgm_B_igrf( v, &B2, Info );
                        Dipole = Lgm_B_igrf;
                        break;
        default:
                        fprintf(stderr, "Lgm_B_FromScatteredData4(): Unknown internal model (%d)\n", Info->InternalModel );
                        break;

    }

    // Total B
    B->x = B1.x + B2.x;
    B->y = B1.y + B2.y;
    B->z = B1.z + B2.z;
//B->x = B1.x;
//B->y = B1.y;
//B->z = B1.z;


    //if ( Info->RBF_CompGradAndCurl ) {
{

        Lgm_Vector u, u0, Bvec;
        int DerivScheme, N;
        Lgm_Vector dBdx, dBdy, dBdz;
        double f1x[7], f2x[7], f3x[7];
        double f1y[7], f2y[7], f3y[7];
        double f1z[7], f2z[7], f3z[7];
        double H, Bmag;
        double h = 1e-3;

        u0 = *v;
        N = 1; DerivScheme = LGM_DERIV_TWO_POINT;
        
        for (i=-N; i<=N; ++i){
            if ( i != 0 ) {
                u = u0; H = (double)i*h; u.x += H;
                Dipole( &u, &Bvec, Info );
                f1x[i+N] = Bvec.x;
                f2x[i+N] = Bvec.y;
                f3x[i+N] = Bvec.z;
            }
        }
        for (i=-N; i<=N; ++i){
            if ( i != 0 ) {
                u = u0; H = (double)i*h; u.y += H;
                Dipole( &u, &Bvec, Info );
                f1y[i+N] = Bvec.x;
                f2y[i+N] = Bvec.y;
                f3y[i+N] = Bvec.z;
            }
        }
        for (i=-N; i<=N; ++i){
            if ( i != 0 ) {
                u = u0; H = (double)i*h; u.z += H;
                Dipole( &u, &Bvec, Info );
                Bmag = Lgm_Magnitude( &Bvec );
                f1z[i+N] = Bvec.x;
                f2z[i+N] = Bvec.y;
                f3z[i+N] = Bvec.z;
            }
        }
        if (DerivScheme == LGM_DERIV_SIX_POINT){

            dBdx.x = (f1x[6] - 9.0*f1x[5] + 45.0*f1x[4] - 45.0*f1x[2] + 9.0*f1x[1] - f1x[0])/(60.0*h);
            dBdx.y = (f2x[6] - 9.0*f2x[5] + 45.0*f2x[4] - 45.0*f2x[2] + 9.0*f2x[1] - f2x[0])/(60.0*h);
            dBdx.z = (f3x[6] - 9.0*f3x[5] + 45.0*f3x[4] - 45.0*f3x[2] + 9.0*f3x[1] - f3x[0])/(60.0*h);

            dBdy.x = (f1y[6] - 9.0*f1y[5] + 45.0*f1y[4] - 45.0*f1y[2] + 9.0*f1y[1] - f1y[0])/(60.0*h);
            dBdy.y = (f2y[6] - 9.0*f2y[5] + 45.0*f2y[4] - 45.0*f2y[2] + 9.0*f2y[1] - f2y[0])/(60.0*h);
            dBdy.z = (f3y[6] - 9.0*f3y[5] + 45.0*f3y[4] - 45.0*f3y[2] + 9.0*f3y[1] - f3y[0])/(60.0*h);

            dBdz.x = (f1z[6] - 9.0*f1z[5] + 45.0*f1z[4] - 45.0*f1z[2] + 9.0*f1z[1] - f1z[0])/(60.0*h);
            dBdz.y = (f2z[6] - 9.0*f2z[5] + 45.0*f2z[4] - 45.0*f2z[2] + 9.0*f2z[1] - f2z[0])/(60.0*h);
            dBdz.z = (f3z[6] - 9.0*f3z[5] + 45.0*f3z[4] - 45.0*f3z[2] + 9.0*f3z[1] - f3z[0])/(60.0*h);

        } else if (DerivScheme == LGM_DERIV_FOUR_POINT){

            dBdx.x = (-f1x[4] + 8.0*f1x[3] - 8.0*f1x[1] + f1x[0])/(12.0*h);
            dBdx.y = (-f2x[4] + 8.0*f2x[3] - 8.0*f2x[1] + f2x[0])/(12.0*h);
            dBdx.z = (-f3x[4] + 8.0*f3x[3] - 8.0*f3x[1] + f3x[0])/(12.0*h);

            dBdy.x = (-f1y[4] + 8.0*f1y[3] - 8.0*f1y[1] + f1y[0])/(12.0*h);
            dBdy.y = (-f2y[4] + 8.0*f2y[3] - 8.0*f2y[1] + f2y[0])/(12.0*h);
            dBdy.z = (-f3y[4] + 8.0*f3y[3] - 8.0*f3y[1] + f3y[0])/(12.0*h);

            dBdz.x = (-f1z[4] + 8.0*f1z[3] - 8.0*f1z[1] + f1z[0])/(12.0*h);
            dBdz.y = (-f2z[4] + 8.0*f2z[3] - 8.0*f2z[1] + f2z[0])/(12.0*h);
            dBdz.z = (-f3z[4] + 8.0*f3z[3] - 8.0*f3z[1] + f3z[0])/(12.0*h);

        } else if (DerivScheme == LGM_DERIV_TWO_POINT){

            dBdx.x = (f1x[2] - f1x[0])/(2.0*h);
            dBdx.y = (f2x[2] - f2x[0])/(2.0*h);
            dBdx.z = (f3x[2] - f3x[0])/(2.0*h);

            dBdy.x = (f1y[2] - f1y[0])/(2.0*h);
            dBdy.y = (f2y[2] - f2y[0])/(2.0*h);
            dBdy.z = (f3y[2] - f3y[0])/(2.0*h);

            dBdz.x = (f1z[2] - f1z[0])/(2.0*h);
            dBdz.y = (f2z[2] - f2z[0])/(2.0*h);
            dBdz.z = (f3z[2] - f3z[0])/(2.0*h);

        } else {
            printf("huh?\n");
            exit(0);
        }




        // Add in the contribution from external field.
        dBdx.x += Info->RBF_dBdx.x; 
        dBdx.y += Info->RBF_dBdx.y; 
        dBdx.z += Info->RBF_dBdx.z;

        dBdy.x += Info->RBF_dBdy.x; 
        dBdy.y += Info->RBF_dBdy.y; 
        dBdy.z += Info->RBF_dBdy.z;

        dBdz.x += Info->RBF_dBdz.x; 
        dBdz.y += Info->RBF_dBdz.y; 
        dBdz.z += Info->RBF_dBdz.z;
/*
*/
/*
dBdx.x = Info->RBF_dBdx.x; 
dBdx.y = Info->RBF_dBdx.y; 
dBdx.z = Info->RBF_dBdx.z;

dBdy.x = Info->RBF_dBdy.x; 
dBdy.y = Info->RBF_dBdy.y; 
dBdy.z = Info->RBF_dBdy.z;

dBdz.x = Info->RBF_dBdz.x; 
dBdz.y = Info->RBF_dBdz.y; 
dBdz.z = Info->RBF_dBdz.z;
*/

        // Compute final GradB
        b = *B; Lgm_NormalizeVector( &b );
        Info->RBF_Grad_B.x = Lgm_DotProduct( &b, &dBdx );
        Info->RBF_Grad_B.y = Lgm_DotProduct( &b, &dBdy );
        Info->RBF_Grad_B.z = Lgm_DotProduct( &b, &dBdz );

        // Compute Curl_B
        Info->RBF_Curl_B.x = Info->RBF_dBdy.z - Info->RBF_dBdz.y;
        Info->RBF_Curl_B.y = Info->RBF_dBdz.x - Info->RBF_dBdx.z;
        Info->RBF_Curl_B.z = Info->RBF_dBdx.y - Info->RBF_dBdy.x;
    }



    return( 1 );

}



/** This routine is the same as Lgm_B_FromScatteredData4(), except that here we
 *  additionally compute the electric field and store it in the Info structure.
 *
 *      \param[in]         v  - array of position vectors
 *      \param[in]         B  - array of B-field vectors at the corresponding v's
 *      \param[in,out]  Info  - Pointer to Lgm_MagModelInfo structure.
 *
 *      \return Always returns 1. Fix this...?
 *
 *      \author  M. G. Henderson
 *      \date    July 15, 2015
 *
 *
 */
int Lgm_B_FromScatteredData5( Lgm_Vector *v, Lgm_Vector *B, Lgm_MagModelInfo *Info ) {

    int                 K, Kgot, n_data, i, j, k;
    double              R, d, dx, dy, dz, d2x, d2y, d2z, d2, ux, uy, uz, vx, vy, vz;
    double              d2min_p[4][1000], d2min_m[4][1000], d2min_min[4], d2min_max[4], d2avg[4];
    double              d2xmin_p[1000], d2xmin_m[1000];
    double              d2ymin_p[1000], d2ymin_m[1000];
    double              d2zmin_p[1000], d2zmin_m[1000];
    int                 navg;
    double              *eps_x, *eps_y, *eps_z;
    Lgm_KdTreeData     *kNN;
    Lgm_Vec_RBF_Info   *rbf, *rbf_e;               // single structure.
    Lgm_Vector         *v_data, *B_data, *E_data, B1, B2;
    Lgm_Vector          b, Grad_B1, GradB_dipole; 
    unsigned long int  *I_data;
    unsigned long int  *LookUpKey;                // key comprised of an array of unsigned long int Id's
    int                 KeyLength;                // length (in bytes) of LookUpKey
    int                 (*Dipole)(Lgm_Vector*, Lgm_Vector*, Lgm_MagModelInfo*);              // tmp Pointer to Bfield function


//if (Info->KdTree_kNN_MaxDist2 < 1e8) printf("AHA Info->KdTree_kNN_MaxDist2 = %g\n", Info->KdTree_kNN_MaxDist2 );
    


    /*
     *  Make sure KdTree has been initialized
     */
    if ( !isfinite(v->x) || !isfinite(v->y) || !isfinite(v->z)  ) {
        printf("Lgm_B_FromScatteredData5(): Error. Input point is not finite!  v = %g %g %g\n", v->x, v->y, v->z );
        B->x = B->y = B->z = 0.0;
        return(-1);
    }

    R = Lgm_Magnitude( v );
    if ( R > 1e6 ) {
        printf("Lgm_B_FromScatteredData5(): Error. Input point is too large!  v = %g %g %g\n", v->x, v->y, v->z );
        B->x = B->y = B->z = 0.0;
        return(-1);
    }


    if ( R > 1.5 ) {

        /*
         *  Allocate space for the K Nearest Neighbors.
         *  This is probably a bit wasteful...
         *  Should cache this array.
         */
        K = Info->KdTree_kNN_k;
        if (K > 2) {

            if ( Info->KdTree_kNN_Alloced == 0 ) {

                LGM_ARRAY_1D( Info->KdTree_kNN, K, Lgm_KdTreeData );
                Info->KdTree_kNN_Alloced = K;

            } else if ( K != Info->KdTree_kNN_Alloced ) {

                /*
                 * kNN is allocated but K has changed. Realloc.
                 */
                LGM_ARRAY_1D_FREE( Info->KdTree_kNN );
                LGM_ARRAY_1D( Info->KdTree_kNN, K, Lgm_KdTreeData );

            }
            Info->KdTree_kNN_Alloced = K;

        } else {

            printf("Lgm_B_FromScatteredData5(): Error. Not enough nearest neighbors specified: Info->KdTree_kNN_k = %d\n", Info->KdTree_kNN_k);
            exit(-1);

        }





double GridRes = EstimateGridRes( v );
//Info->KdTree_kNN_MaxDist2 = GridRes*GridRes*2.25;





        /*
         *  Find the K Nearest Neighbors.
         */
        double q[3];
        q[0] = v->x; q[1] = v->y; q[2] = v->z;
double mike = Info->KdTree_kNN_MaxDist2;
        Lgm_KdTree_kNN( q, 3, Info->KdTree, K, &Kgot, Info->KdTree_kNN_MaxDist2, Info->KdTree_kNN );
        //printf("K, Kgot = %d %d    q = %g %g %g  Info->KdTree_kNN_MaxDist2 = %g GridRes = %g\n", K, Kgot, q[0], q[1], q[2], Info->KdTree_kNN_MaxDist2, GridRes );
if (Kgot < 1) printf("AHA   q = %g %g %g   Kgot = %d   mike = %g Info->KdTree_kNN_MaxDist2 = %g\n", v->x, v->y, v->z, Kgot, mike, Info->KdTree_kNN_MaxDist2 );





        // not needed? Put under verbosity setting?
        for ( i=0; i<Kgot; i++ ) {
            if ( Info->KdTree_kNN[i].Dist2 > Info->KdTree_kNN_MaxDist2){
                printf("Lgm_B_FromScatteredData5(): ERROR - Info->KdTree_kNN[i].Dist2 = %g\n", Info->KdTree_kNN[i].Dist2);
            }
        }


        /*
         *  From the K nearest neighbors, construct a key for the hash-table.
         *  Probably should sort them so that different permutations of the same k
         *  NN's will be identified as the same set. But lets worry about that
         *  later.
         *
         */
if (Kgot < 1) printf("AHA   Kgot = %d   Info->KdTree_kNN_MaxDist2 = %g\n", Kgot, Info->KdTree_kNN_MaxDist2 );
//exit(0);
//B->x = 0.0; B->y = 0.0; B->z = 0.0; return(-1);

        LGM_ARRAY_1D( LookUpKey, Kgot, unsigned long int );
        for ( i=0; i<Kgot; i++ ) LookUpKey[i] = Info->KdTree_kNN[i].Id;
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
        rbf   = NULL;
        HASH_FIND( hh, Info->vec_rbf_ht,   LookUpKey, KeyLength, rbf );
        ++(Info->RBF_nHashFinds);

        rbf_e = NULL;
        HASH_FIND( hh, Info->vec_rbf_e_ht, LookUpKey, KeyLength, rbf_e );
        ++(Info->RBF_nHashFinds);




        /*
         * If key didnt exist in hash-table, we need to compute RBF weights, package up info
         * into a structure and add it to the hash table.
         *
         *  We should only need to check the rbf for the B-field, since E and
         *  B are always done together, if B is there E should be too...
         */
        if ( rbf == NULL ) {

            //printf("Did not find key - computing new rbf coeffs. Kgot = %d\n\n", Kgot );

            n_data = Kgot;




            /*
             *   repack data into arrays. This is wasteful also.
             */
            LGM_ARRAY_1D( eps_x,  n_data, double );
            LGM_ARRAY_1D( eps_y,  n_data, double );
            LGM_ARRAY_1D( eps_z,  n_data, double );
            LGM_ARRAY_1D( I_data, n_data, unsigned long int );
            LGM_ARRAY_1D( v_data, n_data, Lgm_Vector );
            LGM_ARRAY_1D( B_data, n_data, Lgm_Vector );
            LGM_ARRAY_1D( E_data, n_data, Lgm_Vector );
            double *bbb;
            for ( i=0; i<n_data; i++){
                v_data[i].x = Info->KdTree_kNN[i].Position[0];
                v_data[i].y = Info->KdTree_kNN[i].Position[1];
                v_data[i].z = Info->KdTree_kNN[i].Position[2];

                bbb = (double *)Info->KdTree_kNN[i].Object;
//printf("i=%d\n", i);
//printf("bbb[0]=%g\n", bbb[0]);
                B_data[i].x = bbb[0];
                B_data[i].y = bbb[1];
                B_data[i].z = bbb[2];
                E_data[i].x = bbb[3];
                E_data[i].y = bbb[4];
                E_data[i].z = bbb[5];

                I_data[i] = Info->KdTree_kNN[i].Id;
                //eps[i]    = Info->RBF_Eps;
            }
            QSORT( unsigned long int, I_data, n_data, int_lt );








            /*
             * For each point, find its nerest neighbor in the +/- directions for each dimension.
             */
            double d2min[4][1000];
            for ( i=0; i<n_data; i++){
                d2min[3][i] = 9e99;
                d2min_p[0][i] = 9e99; d2min_p[1][i] = 9e99; d2min_p[2][i] = 9e99;
                d2min_m[0][i] = 9e99; d2min_m[1][i] = 9e99; d2min_m[2][i] = 9e99;
                d2xmin_p[i] = 9e99; d2ymin_p[i] = 9e99; d2zmin_p[i] = 9e99;
                d2xmin_m[i] = 9e99; d2ymin_m[i] = 9e99; d2zmin_m[i] = 9e99;

                ux = v_data[i].x;
                uy = v_data[i].y;
                uz = v_data[i].z;

                for ( j=0; j<n_data; j++){

                    if (i!=j){

                        vx = v_data[j].x;
                        vy = v_data[j].y;
                        vz = v_data[j].z;

                        dx = ux - vx; dy = uy - vy; dz = uz - vz;
                        d2x = dx*dx; d2y = dy*dy; d2z = dz*dz; d2 = d2x + d2y + d2z;


                        /*
                         * find closest neighbor in the +/- directions separately. An also keep track of the 
                         * corresponding distance in the dimension.
                         */
                        if ( d2 < d2min[3][i] ) d2min[3][i] = d2;

                        if ( (dx >  1e-3) && (d2 < d2min_p[0][i]) ) { d2min_p[0][i] = d2; d2xmin_p[i] = d2x; }
                        if ( (dx < -1e-3) && (d2 < d2min_m[0][i]) ) { d2min_m[0][i] = d2; d2xmin_m[i] = d2x; }

                        if ( (dy >  1e-3) && (d2 < d2min_p[1][i]) ) { d2min_p[1][i] = d2; d2ymin_p[i] = d2y; }
                        if ( (dy < -1e-3) && (d2 < d2min_m[1][i]) ) { d2min_m[1][i] = d2; d2ymin_m[i] = d2y; }

                        if ( (dz >  1e-3) && (d2 < d2min_p[2][i]) ) { d2min_p[2][i] = d2; d2zmin_p[i] = d2z; }
                        if ( (dz < -1e-3) && (d2 < d2min_m[2][i]) ) { d2min_m[2][i] = d2; d2zmin_m[i] = d2z; }
                            
                    }

                }

            }
                            

            for ( i=0; i<n_data; i++) {
                if ( d2xmin_p[i] > 1e99 ) { d2xmin_p[i] = d2xmin_m[i]; } else if ( d2xmin_m[i] > 1e99 ) { d2xmin_m[i] = d2xmin_p[i]; }
                if ( d2ymin_p[i] > 1e99 ) { d2ymin_p[i] = d2ymin_m[i]; } else if ( d2ymin_m[i] > 1e99 ) { d2ymin_m[i] = d2ymin_p[i]; }
                if ( d2zmin_p[i] > 1e99 ) { d2zmin_p[i] = d2zmin_m[i]; } else if ( d2zmin_m[i] > 1e99 ) { d2zmin_m[i] = d2zmin_p[i]; }
            }





            /*
             * Now assign epsilon vals to each RBF.
             */
            for ( i=0; i<n_data; i++) {
                d2 = (d2xmin_p[i] < d2xmin_m[i]) ? d2xmin_p[i] : d2xmin_m[i]; eps_x[i] = 1.0/(d2*8.0*8.0);
d2 = d2min[3][i];
eps_x[i] = 1.0/(d2*8.0*8.0);
eps_y[i] = 1.0/(d2*8.0*8.0);
eps_z[i] = 1.0/(d2*8.0*8.0);
//                d2 = (d2ymin_p[i] < d2ymin_m[i]) ? d2ymin_p[i] : d2ymin_m[i]; eps_y[i] = 1.0/(d2*8.0*8.0);
//                d2 = (d2zmin_p[i] < d2zmin_m[i]) ? d2zmin_p[i] : d2zmin_m[i]; eps_z[i] = 1.0/(d2*8.0*8.0);

            }

for ( i=0; i<n_data; i++) eps_x[i] = 1.0/(49.0*GridRes*GridRes);
for ( i=0; i<n_data; i++) eps_y[i] = 1.0/(49.0*GridRes*GridRes);
for ( i=0; i<n_data; i++) eps_z[i] = 1.0/(49.0*GridRes*GridRes);


int n_data2;
double ff, rat;
Lgm_Vector *v_data2, *B_data2, *E_data2;
LGM_ARRAY_1D( v_data2, n_data, Lgm_Vector );
LGM_ARRAY_1D( B_data2, n_data, Lgm_Vector );
LGM_ARRAY_1D( E_data2, n_data, Lgm_Vector );
rat = d2min_min[3]/d2min_max[3];
// caching may break with this...
if ( (0==1) && (( rat < 0.8) || (rat > 1.2)) ) {
    ff = 0.0;
    for (i=0; i<3; i++) { if ( d2min_max[i] > ff )  ff = d2min_max[i]; }
ff = d2min_max[3];
    SimplifyPointCloud( v_data, B_data, E_data, n_data, v_data2, B_data2, E_data2, &n_data2, 0.9*ff, 4.0, 20 );
            for (k=0; k<4; k++ ) { d2avg[k] = 0.0; navg = 0; }
            for ( i=0; i<n_data2; i++){

                for (k=0; k<4; k++ ) { d2min[k][i] = 1e6; }

                for ( j=0; j<n_data2; j++){
                    if (i!=j){

                        dx = v_data2[i].x - v_data2[j].x;
                        dy = v_data2[i].y - v_data2[j].y;
                        dz = v_data2[i].z - v_data2[j].z;
                        d2x = dx*dx; d2y = dy*dy; d2z = dz*dz; d2 = d2x + d2y + d2z;

                        d2avg[0] += d2x; d2avg[1] += d2y; d2avg[2] += d2z; d2avg[3] += d2;
                        ++navg;

                        if ( (d2x > 0.0) && (d2x < d2min[0][i]) ) d2min[0][i] = d2x; // x2
                        if ( (d2y > 0.0) && (d2y < d2min[1][i]) ) d2min[1][i] = d2y; // y2
                        if ( (d2z > 0.0) && (d2z < d2min[2][i]) ) d2min[2][i] = d2z; // z2
                        if ( (d2  > 0.0) && (d2  < d2min[3][i]) ) d2min[3][i] = d2;  // d2

                    }
                }

            }
            for (k=0; k<4; k++ ) d2avg[k] /= (double)navg;


            /*
             * From all of the d2min_i vals, find the minimum and maximum of these. I.e., its the
             * "smallest and largest spacing between two neighbors".
             */
            for ( k=0; k<4; k++ ) { d2min_min[k] = 1e31; d2min_max[k] = 0.0; }
            for ( i=0; i<n_data2; i++){
                for ( k=0; k<4; k++ ) {
                    if ( d2min[k][i] > d2min_max[k] ) d2min_max[k] = d2min[k][i];
                    if ( d2min[k][i] < d2min_min[k] ) d2min_min[k] = d2min[k][i];
                }
            }




            /*
             * Now assign epsilon vals to each RBF.
             */
//            rat = d2min_min[0]/d2min_max[0];
//            for ( i=0; i<n_data2; i++) eps_x[i] = 1.0/(d2avg[0]*8.0*8.0) * d2min_max[0]/d2min[0][i];
//            rat = d2min_min[1]/d2min_max[1];
//            for ( i=0; i<n_data2; i++) eps_y[i] = 1.0/(d2avg[1]*8.0*8.0) * d2min_max[1]/d2min[1][i];
//            rat = d2min_min[2]/d2min_max[2];
//            for ( i=0; i<n_data2; i++) eps_z[i] = 1.0/(d2avg[2]*8.0*8.0) * d2min_max[2]/d2min[2][i];
//
////for ( i=0; i<n_data; i++) eps_x[i] = 1.0/(d2avg[3]*5.0*5.0);
////for ( i=0; i<n_data; i++) eps_y[i] = 1.0/(d2avg[3]*5.0*5.0);
////for ( i=0; i<n_data; i++) eps_z[i] = 1.0/(d2avg[3]*5.0*5.0);
//////for ( i=0; i<n_data2; i++) eps_x[i] = 1.0/(d2avg[3]*5.0*5.0) * d2min_max[0]/d2min[0][i];
//////for ( i=0; i<n_data2; i++) eps_y[i] = 1.0/(d2avg[3]*5.0*5.0) * d2min_max[1]/d2min[1][i];
//////for ( i=0; i<n_data2; i++) eps_z[i] = 1.0/(d2avg[3]*5.0*5.0) * d2min_max[2]/d2min[2][i];
//for ( i=0; i<n_data2; i++) {
//    d2 = d2min[0][i]*8.0*8.0; if ( d2 < 64.0 ) d2 = 64.0; eps_x[i] = 1.0/d2;
//    d2 = d2min[1][i]*8.0*8.0; if ( d2 < 64.0 ) d2 = 64.0; eps_y[i] = 1.0/d2;
//    d2 = d2min[2][i]*8.0*8.0; if ( d2 < 64.0 ) d2 = 64.0; eps_z[i] = 1.0/d2;
//}

    rbf   = Lgm_Vec_RBF_Init( I_data, v_data2, B_data2, eps_x, eps_y, eps_z, n_data2, Info->RBF_DoPoly, Info->RBF_Type ); 
    rbf_e = Lgm_Vec_RBF_Init( I_data, v_data2, E_data2, eps_x, eps_y, eps_z, n_data2, Info->RBF_DoPoly, Info->RBF_Type ); 
} else {






            /*
             *  Construct the rbf structure. We dont free these until the hash
             *  table is done with. Note that the hash table will be the only
             *  reference to the pointer.  To free, use
             *  Lgm_B_FromScatteredData_TearDown().
             */
            rbf   = Lgm_Vec_RBF_Init( I_data, v_data, B_data, eps_x, eps_y, eps_z, n_data, Info->RBF_DoPoly, Info->RBF_Type ); 
            rbf_e = Lgm_Vec_RBF_Init( I_data, v_data, E_data, eps_x, eps_y, eps_z, n_data, Info->RBF_DoPoly, Info->RBF_Type ); 
}


            LGM_ARRAY_1D_FREE( eps_x  );
            LGM_ARRAY_1D_FREE( eps_y  );
            LGM_ARRAY_1D_FREE( eps_z  );
            LGM_ARRAY_1D_FREE( I_data );
            LGM_ARRAY_1D_FREE( v_data );
            LGM_ARRAY_1D_FREE( B_data );
            LGM_ARRAY_1D_FREE( E_data );
LGM_ARRAY_1D_FREE( v_data2 );
LGM_ARRAY_1D_FREE( B_data2 );
LGM_ARRAY_1D_FREE( E_data2 );




            //printf("Adding item to hash table\n");
            HASH_ADD_KEYPTR( hh, Info->vec_rbf_ht,   rbf->LookUpKey, KeyLength, rbf );
            HASH_ADD_KEYPTR( hh, Info->vec_rbf_e_ht, rbf->LookUpKey, KeyLength, rbf_e );
            ++(Info->RBF_nHashAdds);

            //printf("Info->vec_rbf_ht_maxsize, Info->vec_rbf_ht_size = %g %g    nEntries = %ld\n", Info->vec_rbf_ht_maxsize, Info->vec_rbf_ht_size, Info->RBF_CB.nEntries );

            int              oldest_i, newest_i;
            double           size, size_e; // memory size in MB
            Lgm_Vec_RBF_Info *oldest_rbf, *tmp_rbf;
            Lgm_Vec_RBF_Info *oldest_rbf_e, *tmp_rbf_e;

            /*
             * Remove rbfs if addition of new ones would exceed memory threshold. Dont try to remove if there arent any there.
             */
            while ( ( (rbf->size + Info->vec_rbf_ht_size) > Info->vec_rbf_ht_maxsize ) && (Info->RBF_CB.nEntries > 0) ) {

                // locate oldest entry
                oldest_i     = Info->RBF_CB.oldest_i;
                oldest_rbf   = Info->RBF_CB.Buf1[ oldest_i ];
                oldest_rbf_e = Info->RBF_CB.Buf2[ oldest_i ];
                if ( oldest_rbf != NULL ) {

                    // there is something there, so delete it from HT
                    // (there should always be something there)
                    HASH_DELETE( hh, Info->vec_rbf_ht,   oldest_rbf   );
                    HASH_DELETE( hh, Info->vec_rbf_e_ht, oldest_rbf_e );

                    // find its size, free it and update Info->vec_rbf_ht_size
                    size   = oldest_rbf->size;
                    size_e = oldest_rbf->size;
                    Lgm_Vec_RBF_Free( oldest_rbf );
                    Lgm_Vec_RBF_Free( oldest_rbf_e );
                    Info->RBF_CB.Buf1[ oldest_i ] = NULL;
                    Info->RBF_CB.Buf2[ oldest_i ] = NULL;
                    Info->vec_rbf_ht_size -= size;
                    Info->vec_rbf_ht_size -= size_e;
                    --Info->RBF_CB.nEntries;

                    // oldest will now be next element in the buffer. "n-1" is the max index slot defined so far.
                    ++oldest_i; if ( oldest_i > Info->RBF_CB.n-1 ) oldest_i = 0; // increment and wrap if necessary
                    Info->RBF_CB.oldest_i = oldest_i;

                }

            }


            /*
             * We now have enough space to add the new entry. Add it to the index just "above" the newest.
             * "N-1" is the max index slot of the whole buf.
             */
            newest_i = Info->RBF_CB.newest_i;
            ++newest_i; if ( newest_i > Info->RBF_CB.N-1 ) newest_i = 0; // increment and wrap if necessary
            tmp_rbf   = Info->RBF_CB.Buf1[ newest_i ];
            tmp_rbf_e = Info->RBF_CB.Buf2[ newest_i ];

            if ( tmp_rbf == NULL ) {
                //printf("appned...\n");
                // open slot -- just add it
                Info->RBF_CB.Buf1[ newest_i ] = rbf;
                Info->RBF_CB.Buf2[ newest_i ] = rbf_e;
                Info->RBF_CB.newest_i = newest_i;
                ++Info->RBF_CB.nEntries;
                if ( newest_i >= Info->RBF_CB.n ) Info->RBF_CB.n = newest_i+1;
                //printf("oldest_1, newest_i = %d %d\n", Info->RBF_CB.oldest_i, Info->RBF_CB.newest_i );
                //printf("\n\n");

            } else {
                // occupied slot -- this must the oldest delete what's there first and then add it
                //printf("replace...\n");
                size = tmp_rbf->size;
                HASH_DELETE( hh, Info->vec_rbf_ht, tmp_rbf );
                Info->vec_rbf_ht_size -= size;
                Lgm_Vec_RBF_Free( tmp_rbf );

                size = tmp_rbf->size;
                HASH_DELETE( hh, Info->vec_rbf_e_ht, tmp_rbf_e );
                Info->vec_rbf_ht_size -= size;
                Lgm_Vec_RBF_Free( tmp_rbf_e );

                Info->RBF_CB.Buf1[ newest_i ] = rbf;
                Info->RBF_CB.Buf2[ newest_i ] = rbf_e;
                Info->RBF_CB.newest_i = newest_i;

                oldest_i = newest_i+1; if ( oldest_i > Info->RBF_CB.N-1 ) oldest_i = 0;
                Info->RBF_CB.oldest_i = oldest_i;
                //printf("oldest_1, newest_i = %d %d\n", Info->RBF_CB.oldest_i, Info->RBF_CB.newest_i );
                //printf("\n\n");
            }

            Info->vec_rbf_ht_size += rbf->size; // increment size
            Info->vec_rbf_ht_size += rbf_e->size; // increment size
            
















//rat = d2min_min[3]/d2min_max[3];
//if ( ( rat < 0.8) || (rat > 1.2) ) {
//
//    B->x = 0.0;
//    B->y = 0.0;
//    B->z = 0.0;
//    return;
//}

        } else {

            //printf("got one\n");
            //printf("        Found: " );
            //for(i=0;i<Kgot; i++) printf(" %ld ", rbf->LookUpKey[i] );
            //printf("\n\n");

        }




        /*
         *  Evaluate Divergence Free Interpolation
         */
        Lgm_Vec_RBF_Eval( v, &B1, rbf );
        Lgm_Vec_RBF_Eval( v, &Info->RBF_E, rbf_e );
        //printf("Evaluating with rbf = %p rbf_e = %p  at v = %g %g %g   B1 = %g %g %g   E = %g %g %g\n\n\n", rbf, rbf_e, v->x, v->y, v->z, B1.x, B1.y, B1.z, Info->RBF_E.x, Info->RBF_E.y, Info->RBF_E.z );


        /*
         *  Evaluate derivatives of Divergence Free Interpolation 
         *  put in the Info structure.
         */
            Lgm_Vec_RBF_Derivs_Eval( v, &Info->RBF_dBdx, &Info->RBF_dBdy, &Info->RBF_dBdz, rbf   );
            Lgm_Vec_RBF_Derivs_Eval( v, &Info->RBF_dEdx, &Info->RBF_dEdy, &Info->RBF_dEdz, rbf_e );
        if ( Info->RBF_CompGradAndCurl ) {
        }


        /*
         *  Cleanup. Free rbf, kNN, etc..
         */
        LGM_ARRAY_1D_FREE( LookUpKey );

    } else {

        B1.x = B1.y = B1.z = 0.0;

    }



    // Save Bfield, so we can temporaily swap it for an internal model (so we can compute GradB)
    switch ( Info->InternalModel ){

        case LGM_CDIP:
                        Lgm_B_cdip( v, &B2, Info );
                        Dipole = Lgm_B_cdip;
                        break;
        case LGM_EDIP:
                        Lgm_B_edip( v, &B2, Info );
                        Dipole = Lgm_B_edip;
                        break;
        case LGM_IGRF:
                        Lgm_B_igrf( v, &B2, Info );
                        Dipole = Lgm_B_igrf;
                        break;
        default:
                        fprintf(stderr, "Lgm_B_FromScatteredData5(): Unknown internal model (%d)\n", Info->InternalModel );
                        break;

    }

    // Total B
    B->x = B1.x + B2.x;
    B->y = B1.y + B2.y;
    B->z = B1.z + B2.z;
B->x = B1.x;
B->y = B1.y;
B->z = B1.z;


    //if ( Info->RBF_CompGradAndCurl ) {
{

        Lgm_Vector u, u0, Bvec;
        int DerivScheme, N;
        Lgm_Vector dBdx, dBdy, dBdz;
        double f1x[7], f2x[7], f3x[7];
        double f1y[7], f2y[7], f3y[7];
        double f1z[7], f2z[7], f3z[7];
        double H, Bmag;
        double h = 1e-3;

        u0 = *v;
        N = 1; DerivScheme = LGM_DERIV_TWO_POINT;
        
        for (i=-N; i<=N; ++i){
            if ( i != 0 ) {
                u = u0; H = (double)i*h; u.x += H;
                Dipole( &u, &Bvec, Info );
                f1x[i+N] = Bvec.x;
                f2x[i+N] = Bvec.y;
                f3x[i+N] = Bvec.z;
            }
        }
        for (i=-N; i<=N; ++i){
            if ( i != 0 ) {
                u = u0; H = (double)i*h; u.y += H;
                Dipole( &u, &Bvec, Info );
                f1y[i+N] = Bvec.x;
                f2y[i+N] = Bvec.y;
                f3y[i+N] = Bvec.z;
            }
        }
        for (i=-N; i<=N; ++i){
            if ( i != 0 ) {
                u = u0; H = (double)i*h; u.z += H;
                Dipole( &u, &Bvec, Info );
                Bmag = Lgm_Magnitude( &Bvec );
                f1z[i+N] = Bvec.x;
                f2z[i+N] = Bvec.y;
                f3z[i+N] = Bvec.z;
            }
        }
        if (DerivScheme == LGM_DERIV_SIX_POINT){

            dBdx.x = (f1x[6] - 9.0*f1x[5] + 45.0*f1x[4] - 45.0*f1x[2] + 9.0*f1x[1] - f1x[0])/(60.0*h);
            dBdx.y = (f2x[6] - 9.0*f2x[5] + 45.0*f2x[4] - 45.0*f2x[2] + 9.0*f2x[1] - f2x[0])/(60.0*h);
            dBdx.z = (f3x[6] - 9.0*f3x[5] + 45.0*f3x[4] - 45.0*f3x[2] + 9.0*f3x[1] - f3x[0])/(60.0*h);

            dBdy.x = (f1y[6] - 9.0*f1y[5] + 45.0*f1y[4] - 45.0*f1y[2] + 9.0*f1y[1] - f1y[0])/(60.0*h);
            dBdy.y = (f2y[6] - 9.0*f2y[5] + 45.0*f2y[4] - 45.0*f2y[2] + 9.0*f2y[1] - f2y[0])/(60.0*h);
            dBdy.z = (f3y[6] - 9.0*f3y[5] + 45.0*f3y[4] - 45.0*f3y[2] + 9.0*f3y[1] - f3y[0])/(60.0*h);

            dBdz.x = (f1z[6] - 9.0*f1z[5] + 45.0*f1z[4] - 45.0*f1z[2] + 9.0*f1z[1] - f1z[0])/(60.0*h);
            dBdz.y = (f2z[6] - 9.0*f2z[5] + 45.0*f2z[4] - 45.0*f2z[2] + 9.0*f2z[1] - f2z[0])/(60.0*h);
            dBdz.z = (f3z[6] - 9.0*f3z[5] + 45.0*f3z[4] - 45.0*f3z[2] + 9.0*f3z[1] - f3z[0])/(60.0*h);

        } else if (DerivScheme == LGM_DERIV_FOUR_POINT){

            dBdx.x = (-f1x[4] + 8.0*f1x[3] - 8.0*f1x[1] + f1x[0])/(12.0*h);
            dBdx.y = (-f2x[4] + 8.0*f2x[3] - 8.0*f2x[1] + f2x[0])/(12.0*h);
            dBdx.z = (-f3x[4] + 8.0*f3x[3] - 8.0*f3x[1] + f3x[0])/(12.0*h);

            dBdy.x = (-f1y[4] + 8.0*f1y[3] - 8.0*f1y[1] + f1y[0])/(12.0*h);
            dBdy.y = (-f2y[4] + 8.0*f2y[3] - 8.0*f2y[1] + f2y[0])/(12.0*h);
            dBdy.z = (-f3y[4] + 8.0*f3y[3] - 8.0*f3y[1] + f3y[0])/(12.0*h);

            dBdz.x = (-f1z[4] + 8.0*f1z[3] - 8.0*f1z[1] + f1z[0])/(12.0*h);
            dBdz.y = (-f2z[4] + 8.0*f2z[3] - 8.0*f2z[1] + f2z[0])/(12.0*h);
            dBdz.z = (-f3z[4] + 8.0*f3z[3] - 8.0*f3z[1] + f3z[0])/(12.0*h);

        } else if (DerivScheme == LGM_DERIV_TWO_POINT){

            dBdx.x = (f1x[2] - f1x[0])/(2.0*h);
            dBdx.y = (f2x[2] - f2x[0])/(2.0*h);
            dBdx.z = (f3x[2] - f3x[0])/(2.0*h);

            dBdy.x = (f1y[2] - f1y[0])/(2.0*h);
            dBdy.y = (f2y[2] - f2y[0])/(2.0*h);
            dBdy.z = (f3y[2] - f3y[0])/(2.0*h);

            dBdz.x = (f1z[2] - f1z[0])/(2.0*h);
            dBdz.y = (f2z[2] - f2z[0])/(2.0*h);
            dBdz.z = (f3z[2] - f3z[0])/(2.0*h);

        } else {
            printf("huh?\n");
            exit(0);
        }




        // Add in the contribution from external field.
        dBdx.x += Info->RBF_dBdx.x; 
        dBdx.y += Info->RBF_dBdx.y; 
        dBdx.z += Info->RBF_dBdx.z;

        dBdy.x += Info->RBF_dBdy.x; 
        dBdy.y += Info->RBF_dBdy.y; 
        dBdy.z += Info->RBF_dBdy.z;

        dBdz.x += Info->RBF_dBdz.x; 
        dBdz.y += Info->RBF_dBdz.y; 
        dBdz.z += Info->RBF_dBdz.z;
dBdx.x = Info->RBF_dBdx.x; 
dBdx.y = Info->RBF_dBdx.y; 
dBdx.z = Info->RBF_dBdx.z;

dBdy.x = Info->RBF_dBdy.x; 
dBdy.y = Info->RBF_dBdy.y; 
dBdy.z = Info->RBF_dBdy.z;

dBdz.x = Info->RBF_dBdz.x; 
dBdz.y = Info->RBF_dBdz.y; 
dBdz.z = Info->RBF_dBdz.z;
/*
*/

        // Compute final GradB
        b = *B; Lgm_NormalizeVector( &b );
        Info->RBF_Grad_B.x = Lgm_DotProduct( &b, &dBdx );
        Info->RBF_Grad_B.y = Lgm_DotProduct( &b, &dBdy );
        Info->RBF_Grad_B.z = Lgm_DotProduct( &b, &dBdz );

        // Compute Curl_B
        Info->RBF_Curl_B.x = Info->RBF_dBdy.z - Info->RBF_dBdz.y;
        Info->RBF_Curl_B.y = Info->RBF_dBdz.x - Info->RBF_dBdx.z;
        Info->RBF_Curl_B.z = Info->RBF_dBdx.y - Info->RBF_dBdy.x;


        // Compute Curl_E
        Info->RBF_Curl_E.x = Info->RBF_dEdy.z - Info->RBF_dEdz.y;
        Info->RBF_Curl_E.y = Info->RBF_dEdz.x - Info->RBF_dEdx.z;
        Info->RBF_Curl_E.z = Info->RBF_dEdx.y - Info->RBF_dEdy.x;
    }



    return( 1 );

}







/*
 * Input a point cloud (With data). Output a simplified cloud.
 */
void SimplifyPointCloud( Lgm_Vector *v_data, Lgm_Vector *B_data, Lgm_Vector *E_data, int n_data,
                            Lgm_Vector *v_data2, Lgm_Vector *B_data2, Lgm_Vector *E_data2, int *n_data2, double s2min_min, double R, double nR ) {


    int         i, j, n2;
    double      xsum, ysum, zsum, x, y, z, x2, y2, z2, r2, min, max;
    double      u, v, w, u2, v2, w2, s2, s2min, ro, ri, ro2, ri2, rInc;
    Lgm_Vector  c0;

    /*
     * Compute centroid of point cloud.
     */
/*
    xsum = ysum = zsum = 0.0;
    for ( i=0; i<n_data; i++){
        xsum += v_data[i].x;
        ysum += v_data[i].y;
        zsum += v_data[i].z;
    }
    c0.x = xsum/(double)n_data;
    c0.y = ysum/(double)n_data;
    c0.z = zsum/(double)n_data;
*/


    /*
     * Compute center of bounding box.
     */
    max = -1e31; min = 1e31;
    for ( i=0; i<n_data; i++){
        if ( v_data[i].x < min ) min = v_data[i].x;
        if ( v_data[i].x > max ) max = v_data[i].x;
    }
    c0.x = (min+max)/2.0;

    max = -1e31; min = 1e31;
    for ( i=0; i<n_data; i++){
        if ( v_data[i].y < min ) min = v_data[i].y;
        if ( v_data[i].y > max ) max = v_data[i].y;
    }
    c0.y = (min+max)/2.0;

    max = -1e31; min = 1e31;
    for ( i=0; i<n_data; i++){
        if ( v_data[i].z < min ) min = v_data[i].z;
        if ( v_data[i].z > max ) max = v_data[i].z;
    }
    c0.z = (min+max)/2.0;



    rInc = R/(double)nR;


    /*
     * Test point cloud simplification in areas where mesh changes resolution.
     *
     * Loop over concentric speherical shells and only add points in that are
     * not too close to those that are already there. Start with largest shell
     * and work your way in towards the centroid.
     *
     */
    n2 = 0;
    //for ( ro = R; ro >= 0.0; ro -= rInc ){
    for ( ro = rInc; ro <= R; ro += rInc ){

        ro2 = ro*ro;

        ri = ro - rInc;
        ri2 = ri*ri;

        for ( i=0; i<n_data; i++){

            x = v_data[i].x - c0.x; x2 = x*x;
            y = v_data[i].y - c0.y; y2 = y*y;
            z = v_data[i].z - c0.z; z2 = z*z;
            r2 = x2 + y2 + z2; // dist squared between centroid and current point

            if ( (r2 <= ro2) && (r2 > ri2) ) {

                /*
                 * it is a candidate for inclusion because it is  in the shell,
                 * but check to see if its far enough away from all the other
                 * points (if any) we have so far.
                 */
                s2min = 1e31;
                for ( j=0; j<n2; j++ ) {
                    u = v_data[i].x - v_data2[j].x; u2 = u*u;
                    v = v_data[i].y - v_data2[j].y; v2 = v*v;
                    w = v_data[i].z - v_data2[j].z; w2 = w*w;
                    s2 = u2 + v2 + w2;
                    if ( s2 < s2min ) s2min = s2;
                }

                if ( s2min > s2min_min ) {
                    // add point to list.
                    v_data2[n2].x = v_data[i].x;
                    v_data2[n2].y = v_data[i].y;
                    v_data2[n2].z = v_data[i].z;
                    B_data2[n2].x = B_data[i].x;
                    B_data2[n2].y = B_data[i].y;
                    B_data2[n2].z = B_data[i].z;
                    E_data2[n2].x = E_data[i].x;
                    E_data2[n2].y = E_data[i].y;
                    E_data2[n2].z = E_data[i].z;
                    ++n2;
                }


            }


        }

    }
    *n_data2 = n2;

    printf("\n\n\n-------------------------------------------\n");
    for ( i=0; i<n_data; i++) printf("Original Point: %g %g %g\n", v_data[i].x, v_data[i].y, v_data[i].z);
    printf("\n");
    for ( i=0; i<*n_data2; i++) printf("Included Point: %g %g %g\n", v_data2[i].x, v_data2[i].y, v_data2[i].z);
    printf("Original number of points: %d  Final number of points: %d\n\n\n", n_data, *n_data2 );

}




double EstimateGridRes( Lgm_Vector *v ){

    double x, y, z;
    double xlo1, xhi1, ylo1, yhi1, zlo1, zhi1;
    double xlo2, xhi2, ylo2, yhi2, zlo2, zhi2;
    double xlo3, xhi3, ylo3, yhi3, zlo3, zhi3;
    double xlo4, xhi4, ylo4, yhi4, zlo4, zhi4;
    double xlo5, xhi5, ylo5, yhi5, zlo5, zhi5;
    double xlo6, xhi6, ylo6, yhi6, zlo6, zhi6;

    x = v->x; y = v->y; z = v->z;

    xlo1 =   -4.0; xhi1 =   4.0; ylo1 =  -4.0; yhi1 =  4.0; zlo1 = -4.0; zhi1 = 4.0;
    xlo2 =   -8.0; xhi2 =   8.0; ylo2 =  -8.0; yhi2 =  8.0; zlo2 = -8.0; zhi2 = 8.0;
    xlo3 =  -32.0; xhi3 =  12.0; ylo3 = -16.0; yhi3 = 16.0; zlo3 = -16.0; zhi3 = 16.0;
    xlo4 =  -40.0; xhi4 =  24.0; ylo4 = -24.0; yhi4 = 24.0; zlo4 = -24.0; zhi4 = 24.0;
    xlo5 =  -48.0; xhi5 =  32.0; ylo5 = -32.0; yhi5 = 24.0; zlo5 = -32.0; zhi5 = 32.0;
    xlo6 = -112.0; xhi6 = -48.0; ylo6 = -16.0; yhi6 = 16.0; zlo6 = -16.0; zhi6 = 16.0;

    if        ( ( x >= xlo1 ) && ( x <= xhi1 ) && ( y >= ylo1 ) && ( y <= yhi1 ) && ( z >= zlo1 ) && ( z<= zhi1 ) ){

        return( 0.0625 );

    } else if ( ( x >= xlo2 ) && ( x <= xhi2 ) && ( y >= ylo2 ) && ( y <= yhi2 ) && ( z >= zlo2 ) && ( z<= zhi2 ) ){

        return( 0.125 );

    } else if ( ( x >= xlo3 ) && ( x <= xhi3 ) && ( y >= ylo3 ) && ( y <= yhi3 ) && ( z >= zlo3 ) && ( z<= zhi3 ) ){

        return( 0.25 );

    } else if ( ( x >= xlo4 ) && ( x <= xhi4 ) && ( y >= ylo4 ) && ( y <= yhi4 ) && ( z >= zlo4 ) && ( z<= zhi4 ) ){

        return( 0.5 );

    } else if ( ( x >= xlo5 ) && ( x <= xhi5 ) && ( y >= ylo5 ) && ( y <= yhi5 ) && ( z >= zlo5 ) && ( z<= zhi5 ) ){

        return( 1.0 );

    } else if ( ( x >= xlo6 ) && ( x <= xhi6 ) && ( y >= ylo6 ) && ( y <= yhi6 ) && ( z >= zlo6 ) && ( z<= zhi6 ) ){


        return( 1.0 );

    } else {

        return( 2.0 );

    }


}




/** This routine is the same as Lgm_B_FromScatteredData4(), except that here we
 *  use the DFI RBF funcs. Experimental... Testing....
 *
 *      \param[in]         v  - array of position vectors
 *      \param[in]         B  - array of B-field vectors at the corresponding v's
 *      \param[in,out]  Info  - Pointer to Lgm_MagModelInfo structure.
 *
 *      \return Always returns 1. Fix this...?
 *
 *      \author  M. G. Henderson
 *      \date    March 8, 2019
 *
 *
 */
int Lgm_B_FromScatteredData6( Lgm_Vector *v, Lgm_Vector *B, Lgm_MagModelInfo *Info ) {

    int                 K, Kgot, n_data, i;
    double              *eps, d;
    Lgm_KdTreeData     *kNN;
    Lgm_DFI_RBF_Info   *rbf;                      // single structure.
    Lgm_Vector         *v_data, *B_data, B1, B2;
    Lgm_Vector          b, Grad_B1, GradB_dipole; 
    unsigned long int  *I_data;
    unsigned long int  *LookUpKey;                // key comprised of an array of unsigned long int Id's
    int                 KeyLength;                // length (in bytes) of LookUpKey
    int                 (*Dipole)(Lgm_Vector*, Lgm_Vector*, Lgm_MagModelInfo*);         // tmp Pointer to Bfield function


    /*
     *  Make sure KdTree has been initialized
     */
    if ( Lgm_Magnitude( v ) > 1.5 ) {

        /*
         *  Allocate space for the K Nearest Neighbors.
         *  This is probably a bit wasteful...
         *  Should cache this array.
         */
        K = Info->KdTree_kNN_k;
        if (K > 2) {

            if ( Info->KdTree_kNN_Alloced == 0 ) {

                LGM_ARRAY_1D( Info->KdTree_kNN, K, Lgm_KdTreeData );
                Info->KdTree_kNN_Alloced = K;

            } else if ( K != Info->KdTree_kNN_Alloced ) {

                /*
                 * kNN is allocated but K has changed. Realloc.
                 */
                LGM_ARRAY_1D_FREE( Info->KdTree_kNN );
                LGM_ARRAY_1D( Info->KdTree_kNN, K, Lgm_KdTreeData );

            }
            Info->KdTree_kNN_Alloced = K;

        } else {

            printf("Lgm_B_FromScatteredData4(): Error. Not enough nearest neighbors specified: Info->KdTree_kNN_k = %d\n", Info->KdTree_kNN_k);
            exit(-1);

        }


        /*
         *  Find the K Nearest Neighbors.
         */
        double q[3];
        q[0] = v->x; q[1] = v->y; q[2] = v->z;
        Lgm_KdTree_kNN( q, 3, Info->KdTree, K, &Kgot, Info->KdTree_kNN_MaxDist2, Info->KdTree_kNN );


        // not needed? Put under verbosity setting?
        for ( i=0; i<Kgot; i++ ) {
            if ( Info->KdTree_kNN[i].Dist2 > Info->KdTree_kNN_MaxDist2){
                printf("Lgm_B_FromScatteredData4(): ERROR - Info->KdTree_kNN[i].Dist2 = %g\n", Info->KdTree_kNN[i].Dist2);
            }
        }


        /*
         *  From the K nearest neighbors, construct a key for the hash-table.
         *  Probably should sort them so that different permutations of the same k
         *  NN's will be identified as the same set. But lets worry about that
         *  later.
         *
         */
        LGM_ARRAY_1D( LookUpKey, Kgot, unsigned long int );
        for ( i=0; i<Kgot; i++ ) LookUpKey[i] = Info->KdTree_kNN[i].Id;
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
        rbf = NULL;
        HASH_FIND( hh, Info->rbf_ht, LookUpKey, KeyLength, rbf );
        ++(Info->RBF_nHashFinds);


//if (0==1){
//}
//rbf = NULL;


        /*
         * If key didnt exist in hash-table, we need to compute RBF weights, package up info
         * into a structure and add it to the hash table.
         */
        if ( rbf == NULL ) {



            //printf("Did not find key - computing new rbf coeffs\n\n");

            /*
             *   repack data into arrays. This is wasteful also.
             */
            n_data = Kgot;
            LGM_ARRAY_1D( eps,    n_data, double );
            LGM_ARRAY_1D( I_data, n_data, unsigned long int );
            LGM_ARRAY_1D( v_data, n_data, Lgm_Vector );
            LGM_ARRAY_1D( B_data, n_data, Lgm_Vector );
            double *bbb;
            for ( i=0; i<n_data; i++){
                v_data[i].x = Info->KdTree_kNN[i].Position[0];
                v_data[i].y = Info->KdTree_kNN[i].Position[1];
                v_data[i].z = Info->KdTree_kNN[i].Position[2];

                bbb = (double *)Info->KdTree_kNN[i].Object;
                B_data[i].x = bbb[0];
                B_data[i].y = bbb[1];
                B_data[i].z = bbb[2];

                I_data[i] = Info->KdTree_kNN[i].Id;
                eps[i]    = Info->RBF_Eps;
            }
            QSORT( unsigned long int, I_data, n_data, int_lt );

/*
double dx, dy, dz, d2, d2min;
int j;
d2min = 1e6;
for ( i=0; i<n_data; i++){
  for ( j=0; i<n_data; i++){
    if (i!=j){

        dx = v_data[i].x - v_data[j].x;
        dy = v_data[i].x - v_data[j].x;
        dz = v_data[i].x - v_data[j].x;
        d2 = dx*dx + dy*dy + dz*dz;
        if ( (d2 > 0.0) && ( d2 < d2min) ) {
            d2min = d2;
        }
    
    }
  }
}
//printf("d2min = %g\n", d2min);
Info->RBF_Eps = 1.0/(d2min*10.0);
*/


            /*
             *  Construct the rbf structure. We dont free these until the hash
             *  table is done with. Note that the hash table will be the only
             *  reference to the pointer.  To free, use
             *  Lgm_B_FromScatteredData_TearDown().
             */
            rbf = Lgm_DFI_RBF_Init( I_data, v_data, B_data, n_data, Info->RBF_Eps, Info->RBF_Type ); 


            LGM_ARRAY_1D_FREE(    eps );
            LGM_ARRAY_1D_FREE( I_data );
            LGM_ARRAY_1D_FREE( v_data );
            LGM_ARRAY_1D_FREE( B_data );

            //printf("Adding item to hash table\n");
            HASH_ADD_KEYPTR( hh, Info->rbf_ht, rbf->LookUpKey, KeyLength, rbf );
            ++(Info->RBF_nHashAdds);

        } else {

            //printf("got one\n");
            //printf("        Found: " );
            //for(i=0;i<Kgot; i++) printf(" %ld ", rbf->LookUpKey[i] );
            //printf("\n\n");

        }




        /*
         *  Evaluate Divergence Free Interpolation
         */
        Lgm_DFI_RBF_Eval( v, &B1, rbf );
        //printf("Evaluating with rbf = %p  at v = %g %g %g   B1 = %g %g %g\n\n\n", rbf, v->x, v->y, v->z, B1.x, B1.y, B1.z);

        /*
         *  Evaluate derivatives of Divergence Free Interpolation 
         *  put in the Info structure.
         */
        if ( Info->RBF_CompGradAndCurl ) {
            Lgm_DFI_RBF_Derivs_Eval( v, &Info->RBF_dBdx, &Info->RBF_dBdy, &Info->RBF_dBdz, rbf );
        }


        /*
         *  Cleanup. Free rbf, kNN, etc..
         */
        LGM_ARRAY_1D_FREE( LookUpKey );

    } else {

        B1.x = B1.y = B1.z = 0.0;

    }



    // Save Bfield, so we can temporaily swap it for an internal model (so we can compute GradB)
    switch ( Info->InternalModel ){

        case LGM_CDIP:
                        Lgm_B_cdip( v, &B2, Info );
                        Dipole = Lgm_B_cdip;
                        break;
        case LGM_EDIP:
                        Lgm_B_edip( v, &B2, Info );
                        Dipole = Lgm_B_edip;
                        break;
        case LGM_IGRF:
                        Lgm_B_igrf( v, &B2, Info );
                        Dipole = Lgm_B_igrf;
                        break;
        default:
                        fprintf(stderr, "Lgm_B_FromScatteredData4(): Unknown internal model (%d)\n", Info->InternalModel );
                        break;

    }

    // Total B
    B->x = B1.x + B2.x;
    B->y = B1.y + B2.y;
    B->z = B1.z + B2.z;
//B->x = B1.x;
//B->y = B1.y;
//B->z = B1.z;


    //if ( Info->RBF_CompGradAndCurl ) {
{

        Lgm_Vector u, u0, Bvec;
        int DerivScheme, N;
        Lgm_Vector dBdx, dBdy, dBdz;
        double f1x[7], f2x[7], f3x[7];
        double f1y[7], f2y[7], f3y[7];
        double f1z[7], f2z[7], f3z[7];
        double H, Bmag;
        double h = 1e-3;

        u0 = *v;
        N = 1; DerivScheme = LGM_DERIV_TWO_POINT;
        
        for (i=-N; i<=N; ++i){
            if ( i != 0 ) {
                u = u0; H = (double)i*h; u.x += H;
                Dipole( &u, &Bvec, Info );
                f1x[i+N] = Bvec.x;
                f2x[i+N] = Bvec.y;
                f3x[i+N] = Bvec.z;
            }
        }
        for (i=-N; i<=N; ++i){
            if ( i != 0 ) {
                u = u0; H = (double)i*h; u.y += H;
                Dipole( &u, &Bvec, Info );
                f1y[i+N] = Bvec.x;
                f2y[i+N] = Bvec.y;
                f3y[i+N] = Bvec.z;
            }
        }
        for (i=-N; i<=N; ++i){
            if ( i != 0 ) {
                u = u0; H = (double)i*h; u.z += H;
                Dipole( &u, &Bvec, Info );
                Bmag = Lgm_Magnitude( &Bvec );
                f1z[i+N] = Bvec.x;
                f2z[i+N] = Bvec.y;
                f3z[i+N] = Bvec.z;
            }
        }
        if (DerivScheme == LGM_DERIV_SIX_POINT){

            dBdx.x = (f1x[6] - 9.0*f1x[5] + 45.0*f1x[4] - 45.0*f1x[2] + 9.0*f1x[1] - f1x[0])/(60.0*h);
            dBdx.y = (f2x[6] - 9.0*f2x[5] + 45.0*f2x[4] - 45.0*f2x[2] + 9.0*f2x[1] - f2x[0])/(60.0*h);
            dBdx.z = (f3x[6] - 9.0*f3x[5] + 45.0*f3x[4] - 45.0*f3x[2] + 9.0*f3x[1] - f3x[0])/(60.0*h);

            dBdy.x = (f1y[6] - 9.0*f1y[5] + 45.0*f1y[4] - 45.0*f1y[2] + 9.0*f1y[1] - f1y[0])/(60.0*h);
            dBdy.y = (f2y[6] - 9.0*f2y[5] + 45.0*f2y[4] - 45.0*f2y[2] + 9.0*f2y[1] - f2y[0])/(60.0*h);
            dBdy.z = (f3y[6] - 9.0*f3y[5] + 45.0*f3y[4] - 45.0*f3y[2] + 9.0*f3y[1] - f3y[0])/(60.0*h);

            dBdz.x = (f1z[6] - 9.0*f1z[5] + 45.0*f1z[4] - 45.0*f1z[2] + 9.0*f1z[1] - f1z[0])/(60.0*h);
            dBdz.y = (f2z[6] - 9.0*f2z[5] + 45.0*f2z[4] - 45.0*f2z[2] + 9.0*f2z[1] - f2z[0])/(60.0*h);
            dBdz.z = (f3z[6] - 9.0*f3z[5] + 45.0*f3z[4] - 45.0*f3z[2] + 9.0*f3z[1] - f3z[0])/(60.0*h);

        } else if (DerivScheme == LGM_DERIV_FOUR_POINT){

            dBdx.x = (-f1x[4] + 8.0*f1x[3] - 8.0*f1x[1] + f1x[0])/(12.0*h);
            dBdx.y = (-f2x[4] + 8.0*f2x[3] - 8.0*f2x[1] + f2x[0])/(12.0*h);
            dBdx.z = (-f3x[4] + 8.0*f3x[3] - 8.0*f3x[1] + f3x[0])/(12.0*h);

            dBdy.x = (-f1y[4] + 8.0*f1y[3] - 8.0*f1y[1] + f1y[0])/(12.0*h);
            dBdy.y = (-f2y[4] + 8.0*f2y[3] - 8.0*f2y[1] + f2y[0])/(12.0*h);
            dBdy.z = (-f3y[4] + 8.0*f3y[3] - 8.0*f3y[1] + f3y[0])/(12.0*h);

            dBdz.x = (-f1z[4] + 8.0*f1z[3] - 8.0*f1z[1] + f1z[0])/(12.0*h);
            dBdz.y = (-f2z[4] + 8.0*f2z[3] - 8.0*f2z[1] + f2z[0])/(12.0*h);
            dBdz.z = (-f3z[4] + 8.0*f3z[3] - 8.0*f3z[1] + f3z[0])/(12.0*h);

        } else if (DerivScheme == LGM_DERIV_TWO_POINT){

            dBdx.x = (f1x[2] - f1x[0])/(2.0*h);
            dBdx.y = (f2x[2] - f2x[0])/(2.0*h);
            dBdx.z = (f3x[2] - f3x[0])/(2.0*h);

            dBdy.x = (f1y[2] - f1y[0])/(2.0*h);
            dBdy.y = (f2y[2] - f2y[0])/(2.0*h);
            dBdy.z = (f3y[2] - f3y[0])/(2.0*h);

            dBdz.x = (f1z[2] - f1z[0])/(2.0*h);
            dBdz.y = (f2z[2] - f2z[0])/(2.0*h);
            dBdz.z = (f3z[2] - f3z[0])/(2.0*h);

        } else {
            printf("huh?\n");
            exit(0);
        }




        // Add in the contribution from external field.
        dBdx.x += Info->RBF_dBdx.x; 
        dBdx.y += Info->RBF_dBdx.y; 
        dBdx.z += Info->RBF_dBdx.z;

        dBdy.x += Info->RBF_dBdy.x; 
        dBdy.y += Info->RBF_dBdy.y; 
        dBdy.z += Info->RBF_dBdy.z;

        dBdz.x += Info->RBF_dBdz.x; 
        dBdz.y += Info->RBF_dBdz.y; 
        dBdz.z += Info->RBF_dBdz.z;
/*
*/
/*
dBdx.x = Info->RBF_dBdx.x; 
dBdx.y = Info->RBF_dBdx.y; 
dBdx.z = Info->RBF_dBdx.z;

dBdy.x = Info->RBF_dBdy.x; 
dBdy.y = Info->RBF_dBdy.y; 
dBdy.z = Info->RBF_dBdy.z;

dBdz.x = Info->RBF_dBdz.x; 
dBdz.y = Info->RBF_dBdz.y; 
dBdz.z = Info->RBF_dBdz.z;
*/

        // Compute final GradB
        b = *B; Lgm_NormalizeVector( &b );
        Info->RBF_Grad_B.x = Lgm_DotProduct( &b, &dBdx );
        Info->RBF_Grad_B.y = Lgm_DotProduct( &b, &dBdy );
        Info->RBF_Grad_B.z = Lgm_DotProduct( &b, &dBdz );

        // Compute Curl_B
        Info->RBF_Curl_B.x = Info->RBF_dBdy.z - Info->RBF_dBdz.y;
        Info->RBF_Curl_B.y = Info->RBF_dBdz.x - Info->RBF_dBdx.z;
        Info->RBF_Curl_B.z = Info->RBF_dBdx.y - Info->RBF_dBdy.x;
    }



    return( 1 );

}
/*
 *  Setup the hash table used in Lgm_B_FromScatteredData6().
 */
void Lgm_B_FromScatteredData6_SetUp( Lgm_MagModelInfo *Info ) {


Lgm_DFI_RBF_Info ***b;

    if ( Info->rbf_ht_alloced ) Lgm_B_FromScatteredData6_TearDown( Info );
    Info->rbf_ht     = NULL;
    Info->vec_rbf_e_ht   = NULL;
    Info->rbf_ht_alloced = FALSE;
    Info->RBF_nHashFinds = 0;
    Info->RBF_nHashAdds  = 0;


    Info->RBF_CB.n           = 0;
    Info->RBF_CB.nEntries    =      0; // Initial entries
    Info->RBF_CB.N           = 200000; // Max entries
    Info->RBF_CB.Buf1        = (Lgm_Vec_RBF_Info **)calloc( Info->RBF_CB.N, sizeof( Lgm_Vec_RBF_Info *) );
    Info->RBF_CB.Buf2        = (Lgm_Vec_RBF_Info **)calloc( Info->RBF_CB.N, sizeof( Lgm_Vec_RBF_Info *) );

    Info->dfi_rbf_ht_size    =    0.0; // Initial size in MB
    Info->dfi_rbf_ht_maxsize = 2000.0; // Max size in MB
    Info->vec_rbf_ht_size    =    0.0; // Initial size in MB
    Info->vec_rbf_ht_maxsize = 2000.0; // Max size in MB
//    Info->vec_rbf_ht_maxsize = 200.0; // Max size in MB
    Info->RBF_CB.oldest_i    = 0;
    Info->RBF_CB.newest_i    = -1;
    

}
/*
 *  Iterates over all the entries in the hash table and 1) deletes them from
 *  the hash table, then 2) free the structure itself.
 *
 *  Really should unify the rbf structures..
 */
void Lgm_B_FromScatteredData6_TearDown( Lgm_MagModelInfo *Info ) {

    Lgm_DFI_RBF_Info *rbf, *rbf_tmp;

    HASH_ITER( hh, Info->rbf_ht, rbf, rbf_tmp ) {
        HASH_DELETE( hh, Info->rbf_ht, rbf );
        Lgm_DFI_RBF_Free( rbf );
    }

    if ( Info->Octree_kNN_Alloced > 0 ) {
        LGM_ARRAY_1D_FREE( Info->Octree_kNN );
        Info->Octree_kNN_Alloced = 0;
    }
    if ( Info->KdTree_kNN_Alloced > 0 ) {
        LGM_ARRAY_1D_FREE( Info->KdTree_kNN );
        Info->KdTree_kNN_Alloced = 0;
    }

}
