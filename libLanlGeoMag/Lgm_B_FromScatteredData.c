/*! \file Lgm_B_FromScatteredData.c
 *  \brief Routines to compute B-field from a sattered set of data (for example, B defined on meshes).
 *
 *
 *
 */

#include "Lgm/Lgm_MagModelInfo.h"
#include "Lgm/Lgm_Octree.h"
#include "Lgm/Lgm_DFI_RBF.h"



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
    rbf = Lgm_DFI_RBF_Init( v_data, B_data, n_data, eps, LGM_RBF_GAUSSIAN );
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

    int               K, Kgot, n_data, i;
    double            eps, d;
    Lgm_OctreeData   *kNN;
    Lgm_DFI_RBF_Info *rbf;
    Lgm_Vector       *v_data, *B_data, B1, B2;


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
B1.x = B1.y = B1.z = 0.0;
    for ( i=0; i<n_data; i++){
        Lgm_OctreeUnScalePosition( &(kNN[i].Position), &v_data[i], Info->Octree );
        B_data[i] = kNN[i].B;
B1.x += B_data[i].x;
B1.y += B_data[i].y;
B1.z += B_data[i].z;
    }

B1.x /= (double)n_data;
B1.y /= (double)n_data;
B1.z /= (double)n_data;




if (0==1){

    /*
     *  Initialize the  Divergence Free Interp. RBF stuff
     */
    //Lgm_OctreeUnScaleDistance( 1.0, &d, Info->Octree );
    //eps = 1.0/d;
    eps = 0.001;

//eps = 0.1;
    rbf = Lgm_DFI_RBF_Init( v_data, B_data, n_data, eps, LGM_RBF_GAUSSIAN );
    //rbf = Lgm_DFI_RBF_Init( v_data, B_data, n_data, eps, LGM_RBF_MULTIQUADRIC );


    /*
     *  Evaluate Divergence Free Interpolation
     */
    Lgm_DFI_RBF_Eval( v, &B1, rbf );




    /*
     *  Cleanup. Free rbf, kNN, etc..
     */
    Lgm_DFI_RBF_Free( rbf );
}
    LGM_ARRAY_1D_FREE( kNN );
    LGM_ARRAY_1D_FREE( v_data );
    LGM_ARRAY_1D_FREE( B_data );




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
