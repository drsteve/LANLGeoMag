/*! \file Lgm_B_FromScatteredData.c
 *  \brief Routines to compute B-field from a sattered set of data (for example, B defined on meshes).
 *
 *
 *
 */

#include "Lgm/Lgm_MagModelInfo.h"
#include "Lgm/Lgm_Octree.h"
#include "Lgm/Lgm_KdTree.h"
#include "Lgm/Lgm_RBF.h"
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
    }
    if ( Info->KdTree_kNN_Alloced > 0 ) {
        LGM_ARRAY_1D_FREE( Info->KdTree_kNN );
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
    }
    if ( Info->KdTree_kNN_Alloced > 0 ) {
        LGM_ARRAY_1D_FREE( Info->KdTree_kNN );
    }

}



/*
 *  Setup the hash table used in Lgm_B_FromScatteredData().
 */
void Lgm_B_FromScatteredData5_SetUp( Lgm_MagModelInfo *Info ) {

    if ( Info->rbf_ht_alloced ) Lgm_B_FromScatteredData5_TearDown( Info );
    Info->vec_rbf_ht     = NULL;
    Info->rbf_ht_alloced = FALSE;
    Info->RBF_nHashFinds = 0;
    Info->RBF_nHashAdds  = 0;
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
    }
    if ( Info->KdTree_kNN_Alloced > 0 ) {
        LGM_ARRAY_1D_FREE( Info->KdTree_kNN );
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
    double              eps, d;
    Lgm_KdTreeData     *kNN;
    Lgm_DFI_RBF_Info   *rbf;                      // single structure.
    Lgm_Vector         *v_data, *B_data, B1, B2;
    Lgm_Vector          b, Grad_B1, GradB_dipole; 
    unsigned long int  *I_data;
    unsigned long int  *LookUpKey;                // key comprised of an array of unsigned long int Id's
    int                 KeyLength;                // length (in bytes) of LookUpKey
    int                 (*Dipole)();         // tmp Pointer to Bfield function


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
            Lgm_DFI_RBF_Derivs_Eval( v, &Info->RBF_dBdx, &Info->RBF_dBdy, &Info->RBF_dBdz, rbf );
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
                        fprintf(stderr, "Lgm_B_FromScatteredData3: Unknown internal model (%d)\n", Info->InternalModel );
                        break;

    }

    // Total B
    B->x = B1.x + B2.x;
    B->y = B1.y + B2.y;
    B->z = B1.z + B2.z;


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
*/
dBdx.x = Info->RBF_dBdx.x; 
dBdx.y = Info->RBF_dBdx.y; 
dBdx.z = Info->RBF_dBdx.z;

dBdy.x = Info->RBF_dBdy.x; 
dBdy.y = Info->RBF_dBdy.y; 
dBdy.z = Info->RBF_dBdy.z;

dBdz.x = Info->RBF_dBdz.x; 
dBdz.y = Info->RBF_dBdz.y; 
dBdz.z = Info->RBF_dBdz.z;

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
    double              eps, d;
    Lgm_KdTreeData     *kNN;
    Lgm_Vec_RBF_Info   *rbf;                      // single structure.
    Lgm_Vector         *v_data, *B_data, B1, B2;
    Lgm_Vector          b, Grad_B1, GradB_dipole; 
    unsigned long int  *I_data;
    unsigned long int  *LookUpKey;                // key comprised of an array of unsigned long int Id's
    int                 KeyLength;                // length (in bytes) of LookUpKey
    int                 (*Dipole)();         // tmp Pointer to Bfield function


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
        //printf("K, Kgot = %d %d    q = %g %g %g  Info->KdTree_kNN_MaxDist2 = %g \n", K, Kgot, q[0], q[1], q[2], Info->KdTree_kNN_MaxDist2 );


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
            }
            QSORT( unsigned long int, I_data, n_data, int_lt );

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
/*
*/


            /*
             *  Construct the rbf structure. We dont free these until the hash
             *  table is done with. Note that the hash table will be the only
             *  reference to the pointer.  To free, use
             *  Lgm_B_FromScatteredData_TearDown().
             */
            rbf = Lgm_Vec_RBF_Init( I_data, v_data, B_data, n_data, Info->RBF_Eps, Info->RBF_Type ); 


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
            Lgm_Vec_RBF_Derivs_Eval( v, &Info->RBF_dBdx, &Info->RBF_dBdy, &Info->RBF_dBdz, rbf );
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

    int                 K, Kgot, n_data, i;
    double              eps, d;
    Lgm_KdTreeData     *kNN;
    Lgm_Vec_RBF_Info   *rbf, *rbf_e;               // single structure.
    Lgm_Vector         *v_data, *B_data, *E_data, B1, B2;
    Lgm_Vector          b, Grad_B1, GradB_dipole; 
    unsigned long int  *I_data;
    unsigned long int  *LookUpKey;                // key comprised of an array of unsigned long int Id's
    int                 KeyLength;                // length (in bytes) of LookUpKey
    int                 (*Dipole)();              // tmp Pointer to Bfield function


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

            printf("Lgm_B_FromScatteredData5(): Error. Not enough nearest neighbors specified: Info->KdTree_kNN_k = %d\n", Info->KdTree_kNN_k);
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

            //printf("Did not find key - computing new rbf coeffs\n\n");

            /*
             *   repack data into arrays. This is wasteful also.
             */
            n_data = Kgot;
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
                B_data[i].x = bbb[0];
                B_data[i].y = bbb[1];
                B_data[i].z = bbb[2];
                E_data[i].x = bbb[3];
                E_data[i].y = bbb[4];
                E_data[i].z = bbb[5];

                I_data[i] = Info->KdTree_kNN[i].Id;
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
        if ( (d2 > .25*.25) && ( d2 < d2min) ) {
            d2min = d2;
        }
    
    }
  }
}
//printf("d2min = %g\n", d2min);
Info->RBF_Eps = 1.0/(d2min*4.0);
*/


            /*
             *  Construct the rbf structure. We dont free these until the hash
             *  table is done with. Note that the hash table will be the only
             *  reference to the pointer.  To free, use
             *  Lgm_B_FromScatteredData_TearDown().
             */
            rbf   = Lgm_Vec_RBF_Init( I_data, v_data, B_data, n_data, Info->RBF_Eps, Info->RBF_Type ); 
            rbf_e = Lgm_Vec_RBF_Init( I_data, v_data, E_data, n_data, Info->RBF_Eps, Info->RBF_Type ); 


            LGM_ARRAY_1D_FREE( I_data );
            LGM_ARRAY_1D_FREE( v_data );
            LGM_ARRAY_1D_FREE( B_data );
            LGM_ARRAY_1D_FREE( E_data );

            //printf("Adding item to hash table\n");
            HASH_ADD_KEYPTR( hh, Info->vec_rbf_ht,   rbf->LookUpKey, KeyLength, rbf );
            HASH_ADD_KEYPTR( hh, Info->vec_rbf_e_ht, rbf->LookUpKey, KeyLength, rbf_e );
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


        // Compute Curl_E
        Info->RBF_Curl_E.x = Info->RBF_dEdy.z - Info->RBF_dEdz.y;
        Info->RBF_Curl_E.y = Info->RBF_dEdz.x - Info->RBF_dEdx.z;
        Info->RBF_Curl_E.z = Info->RBF_dEdx.y - Info->RBF_dEdy.x;
    }



    return( 1 );

}


