/* Copyright (c) 2008 Michael G. Henderson <mghenderson@lanl.gov>
 *
 * Permission to use, copy, modify, distribute, and sell this software and its
 * documentation for any purpose is hereby granted without fee, provided that
 * the above copyright notice appear in all copies and that both that
 * copyright notice and this permission notice appear in supporting
 * documentation.  No representations are made about the suitability of this
 * software for any purpose.  It is provided "as is" without express or
 * implied warranty.
 *
 *
 *
 *  Computes B-field from a sattered set of data. E.g. from B defined on meshes
 *
 *
 *
 */

#include "Lgm/Lgm_MagModelInfo.h"
#include "Lgm/Lgm_Octree.h"
#include <stdlib.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
/*
 * 3D basis funcs (order = 1   i.e. linear)
 * There are only 4 basis funcs: 1, x, y, z
 */
void DFI_BasisFuncLin( double *b, Lgm_Vector *P ) {
    b[0] = 1.0; b[1] = P->x; b[2] = P->y; b[3] = P->z;
    return;
}
void DFI_FuncLin( Lgm_Vector *P, gsl_vector *C, Lgm_Vector *V ){
    int         i;
    double      b[4];
    DFI_BasisFuncLin( b, P );
    V->x = V->y = V->z = 0.0;
    for (i=0;  i<4; i++) V->x += b[i  ]*C->data[i];
    for (i=4;  i<8; i++) V->y += b[i-4]*C->data[i];
    for (i=8; i<12; i++) V->z += b[i-8]*C->data[i];
    return;
}
void DFI_DivBasisFuncLin( double *bdiv, Lgm_Vector *P ) {
    bdiv[0] = 0.0; bdiv[1] = 1.0; bdiv[2] = 0.0; bdiv[3] = 0.0; 
    bdiv[4] = 0.0; bdiv[5] = 0.0; bdiv[6] = 1.0; bdiv[7] = 0.0; 
    bdiv[8] = 0.0; bdiv[9] = 0.0; bdiv[10] = 0.0; bdiv[11] = 1.0; 
    return;
}




/*
 *  3D basis funcs (order = 2   i.e. quadratic)
 *  There are 10 basis funcs: 1, x, y, z, x^2, y^2, z^2, xy, yz, xz
 *  Using this one may require: (1) data that isnt noisy, (2) lots of extra
 *  points (e.g. more than just ~10).
 */
void DFI_BasisFunc( double *b, Lgm_Vector *P ) {
    b[0] = 1.0; b[1] = P->x; b[2] = P->y; b[3] = P->z;
    b[4] = P->x*P->x; b[5] = P->y*P->y; b[6] = P->z*P->z;
    b[7] = P->x*P->y; b[8] = P->y*P->z; b[9] = P->x*P->z;
    return;
}

void DFI_Func( Lgm_Vector *P, gsl_vector *C, Lgm_Vector *V ){
    int         i;
    double      b[10];
    DFI_BasisFunc( b, P );
    V->x = V->y = V->z = 0.0;
    for (i=0;  i<10; i++) V->x += b[i   ]*C->data[i];
    for (i=10; i<20; i++) V->y += b[i-10]*C->data[i];
    for (i=20; i<30; i++) V->z += b[i-20]*C->data[i];
    return;
}

/*
 * Components needed for divergence of basis functions.
 */
void DFI_DivBasisFunc( double *bdiv, Lgm_Vector *P ) {
    bdiv[0] = 0.0; bdiv[1] = 1.0; bdiv[2] = 0.0;
    bdiv[3] = 0.0; bdiv[4] = 2.0*P->x; bdiv[5] = 0.0;
    bdiv[6] = 0.0; bdiv[7] = P->y; bdiv[8] = 0.0; bdiv[9] = P->z;
    bdiv[10] = 0.0; bdiv[11] = 0.0; bdiv[12] = 1.0;
    bdiv[13] = 0.0; bdiv[14] = 0.0; bdiv[15] = 2.0*P->y;
    bdiv[16] = 0.0; bdiv[17] = P->x; bdiv[18] = P->z; bdiv[19] = 0.0;
    bdiv[20] = 0.0; bdiv[21] = 0.0; bdiv[22] = 0.0;
    bdiv[23] = 1.0; bdiv[24] = 0.0; bdiv[25] = 0.0;
    bdiv[26] = 2.0*P->z; bdiv[27] = 0.0; bdiv[28] = P->y; bdiv[29] = P->x;
    return;
}


int Lgm_DivFreeInterp( Lgm_Vector *q, Lgm_OctreeData *kNN, int K, Lgm_Vector *v ) {

    int             N, M, MB, rows, i, j, k;
    double          *B, *Phi, *W, b[10], bdiv[30], eps;
    gsl_matrix      *A;
    gsl_vector      *C, *Phi2, *tau;
double xx, yy, zz, Dist2;

//printf("K = %d\n", K);
//printf("q = %g %g %g\n", q->x, q->y, q->z);

    if ( K < 8 ) return( -1 ); // under-determined system



    N   = K;    // # of data points (each is a 3D vector)
    MB  = 10;   // # of basis functions 
    M   = 3*MB; // there is a set of basis funcs for each dimension
    eps = 0.0025; // param in weight function (W = 1/(d^2 + eps) )


    /*
     * We have a system of equations such that:
     *      B W c = W Phi
     * B is:
     *           [  b(p1)   0     0   ]
     *           [  0     b(p1)   0   ]
     *           [  0       0   b(p1) ]
     *       B = [          ...       ]  (4*N+1 (row) x M (col) matrix )
     *           [  b(pN)   0     0   ]
     *           [  0     b(pN)   0   ]
     *           [  0       0   b(pN) ]
     *
     *        W = Diag( w1,w2, ..., wN)  (4*Nx4*N diag matrix)
     *  
     *              [ ux( p1 ) ]
     *              [ uy( p1 ) ]
     *              [ uz( p1 ) ]
     *              [ ux( p2 ) ]
     *              [ uy( p2 ) ]   (M x 1 matrix)
     *       Phi =  [ uz( p2 ) ]
     *              [    .     ]
     *              [ ux( pN ) ]
     *              [ uy( pN ) ]
     *              [ uz( pN ) ]
     *  
     *  
     *              [ cx1 ]
     *              [ cx2 ]
     *              [  .  ]
     *              [ cxN ]
     *              [ cy1 ]
     *              [ cy2 ]
     *       c =    [  .  ]   (M x 1 matrix)
     *              [ cyN ]
     *              [ cz1 ]
     *              [ cz2 ]
     *              [  .  ]
     *              [ czN ]
     *  
     *  To solve for c, use normal equations: B^T W^2 B c = B^T W Phi
     *  And solve for c using QR method.
     *  
     */

    /*
     * Set up B matrix (I guess we are absorbing W into it here)
     */
    B = (double *)calloc( (4*N+1)*M, sizeof(double) );
    Phi = (double *) calloc( 4*N+1, sizeof( double ) );
    W   = (double *) calloc( N, sizeof( double ) );
    for ( rows=0,i=0; i<N; i++ ){
        DFI_BasisFunc( b, &kNN[i].Position );
xx = kNN[i].Position.x-q->x;
yy = kNN[i].Position.y-q->y;
zz = kNN[i].Position.z-q->z;
Dist2 = xx*xx +yy*yy + zz*zz;
        W[i] = 1.0/( Dist2*Dist2 + eps );
//printf("kNN: x, y, z, dist2, W =  %g %g %g %g %g %g %g %g\n", kNN[i].Position.x, kNN[i].Position.y, kNN[i].Position.z, kNN[i].Dist2, W[i], kNN[i].B.x, kNN[i].B.y, kNN[i].B.z);

        // X-comp
        for (j=0;j<10;j++) *(B + rows*M + j) = W[i]*b[j];
        for (j=10;j<30;j++) *(B + rows*M + j) = 0.0;
        Phi[rows++] = W[i]*kNN[i].B.x;

        // Y-comp
        for (j=0; j<10;j++) *(B + rows*M + j) = 0.0;
        for (j=10;j<20;j++) *(B + rows*M + j) = W[i]*b[j-10];
        for (j=20;j<30;j++) *(B + rows*M + j) = 0.0;
        Phi[rows++] = W[i]*kNN[i].B.y;

        // Z-comp
        for (j=0; j<20;j++) *(B + rows*M + j) = 0.0;
        for (j=20;j<30;j++) *(B + rows*M + j) = W[i]*b[j-20];
        Phi[rows++] = W[i]*kNN[i].B.z;

        // Div = 0 constraint
//        DFI_DivBasisFunc( bdiv, &kNN[i].Position );
//        for (j=0;j<30;j++) *(B + rows*M + j) = W[i]*bdiv[j];
//        Phi[rows++] = 0.0;
    }

    // Div = 0 constraint at the query point
//    DFI_DivBasisFunc( bdiv, q );
//    for (j=0;j<30;j++) *(B + rows*M + j) = bdiv[j]/eps;
//    Phi[rows++] = 0.0;


//for (i=0; i<rows; i++){
//  for (j=0; j<M; j++){
//    printf(" %10g", *(B+  i*M + j ) );
//  }
//  printf("\n");
//}





//printf("rows = %d\n", rows);
//printf("(4*N+1)*M = %d\n", (4*N+1)*M);
//printf("4*N+1 = %d\n", 4*N+1);
//
//for (k=0; k<rows; k++) {
//    printf("Phi = %g\n", Phi[k]);
//}

    /*
     * Compute B^T B  (=A) and its Cholesky decomposition
     */
    A = gsl_matrix_calloc( M, M );
    for (i=0; i<M; i++){
	for (j=0; j<M; j++){
	    A->data[ i*M + j ] = 0.0;
            for (k=0; k<rows; k++) A->data[ i*M + j ] += (*(B + k*M + i) * *(B + k*M + j));
        }
    }

//for (i=0; i<M; i++){
//  for (j=0; j<M; j++){
//    printf(" %10g", A->data[ i*M + j] );
//  }
//  printf("\n");
//}
//  printf("\n");

    tau = gsl_vector_alloc( M );
    gsl_linalg_QR_decomp( A, tau );

//    printf(" gsl_linalg_cholesky_decomp = %d\n", gsl_linalg_cholesky_decomp( A ));
//    printf(" GSL_EDO. = %d\n", GSL_EDOM);

//for (i=0; i<M; i++){
//  for (j=0; j<M; j++){
//    printf(" %10g", A->data[ i*M + j] );
//  }
//  printf("\n");
//}
//  printf("\n");

    /*
     * Compute B^T Phi
     */
    Phi2 = gsl_vector_alloc( M );
    for (j=0; j<M; j++){
        Phi2->data[ j ] = 0.0;
        for (k=0; k<rows; k++) Phi2->data[ j ] += *(B + k*M + j) * Phi[k];
//printf("Phi2->data[%d] = %g\n", j, Phi2->data[j]);
    }

    /*
     *  Solve for C
     */
    C = gsl_vector_alloc( M );
    gsl_linalg_QR_solve( A, tau, Phi2, C );
//    gsl_linalg_cholesky_solve( A, Phi2, C );


//for (i=0; i<M; i++){
//printf("C->data[i] = %g\n", C->data[i]);
//}
    /*
     * Compute vector value for the query point using the function defined by C
     */
    DFI_Func( q, C, v );


    free( B );
    gsl_vector_free( C );
    gsl_vector_free( Phi2 );
    gsl_vector_free( tau );
    gsl_matrix_free( A );
    free( Phi );
    free( W );

    return( 1 );
    
}



int Lgm_DivFreeInterp2( Lgm_Vector *q, Lgm_OctreeData *kNN, int K, Lgm_Vector *v ) {

    int             N, M, MB, rows, i, j, k;
    double          *B, *Phi, *W, b[4], bdiv[12], eps;
    gsl_matrix      *A;
    gsl_vector      *C, *Phi2, *tau, *resid;
double xx, yy, zz, Dist2;


    if ( K < 8 ) return( -1 ); // under-determined system



    N   = K;    // # of data points (each is a 3D vector)
    MB  = 4;   // # of basis functions 
    M   = 3*MB; // there is a set of basis funcs for each dimension
    eps = 0.00025; // param in weight function (W = 1/(d^2 + eps) )



    /*
     * Set up B matrix (I guess we are absorbing W into it here)
     */ 
    B = (double *)calloc( (4*N+1)*M, sizeof(double) );
    Phi = (double *) calloc( 4*N+1, sizeof( double ) );
    W   = (double *) calloc( N, sizeof( double ) );
    for ( rows=0,i=0; i<N; i++ ){
        DFI_BasisFuncLin( b, &kNN[i].Position );
xx = kNN[i].Position.x-q->x;
yy = kNN[i].Position.y-q->y;
zz = kNN[i].Position.z-q->z;
Dist2 = xx*xx +yy*yy + zz*zz;
        //W[i] = 1.0/( kNN[i].Dist2*kNN[i].Dist2 + eps );
        W[i] = 1.0/( Dist2*Dist2 + eps );
        //W[i] = 1.0/( kNN[i].Dist2 + eps );

        // X-comp
        for (j=0;j<4;j++) *(B + rows*M + j) = W[i]*b[j];
        for (j=4;j<12;j++) *(B + rows*M + j) = 0.0;
        Phi[rows++] = W[i]*kNN[i].B.x;

        // Y-comp
        for (j=0; j<4;j++) *(B + rows*M + j) = 0.0;
        for (j=4;j<8;j++) *(B + rows*M + j) = W[i]*b[j-4];
        for (j=8;j<12;j++) *(B + rows*M + j) = 0.0;
        Phi[rows++] = W[i]*kNN[i].B.y;

        // Z-comp
        for (j=0; j<8;j++) *(B + rows*M + j) = 0.0;
        for (j=8;j<12;j++) *(B + rows*M + j) = W[i]*b[j-8];
        Phi[rows++] = W[i]*kNN[i].B.z;

        // Div = 0 constraint
        DFI_DivBasisFuncLin( bdiv, &kNN[i].Position );
        for (j=0;j<12;j++) *(B + rows*M + j) = W[i]*bdiv[j];
        Phi[rows++] = 0.0;
    }

    // Div = 0 constraint at the query point
    DFI_DivBasisFuncLin( bdiv, q );
    for (j=0;j<12;j++) *(B + rows*M + j) = bdiv[j]/eps;
    Phi[rows++] = 0.0;

    /*
     * Compute B^T B  (=A) and its Cholesky decomposition
     */
    A = gsl_matrix_calloc( M, M );
    for (i=0; i<M; i++){
	for (j=0; j<M; j++){
	    A->data[ i*M + j ] = 0.0;
            //for (k=0; k<rows; k++) A->data[ i*M + j ] += (*(B + k*M + i) * *(B + k*M + j));
            for (k=0; k<rows; k++) A->data[ i*M + j ] += (*(B + k*M + i) * *(B + k*M + j));
        }
    }

    tau = gsl_vector_alloc( M );
    gsl_linalg_QR_decomp( A, tau );

//gsl_linalg_cholesky_decomp( A );
//    printf(" gsl_linalg_cholesky_decomp = %d\n", gsl_linalg_cholesky_decomp( A ));
//    printf(" GSL_EDO. = %d\n", GSL_EDOM);

    /*
     * Compute B^T Phi
     */
    Phi2 = gsl_vector_alloc( M );
    for (j=0; j<M; j++){
        Phi2->data[ j ] = 0.0;
        for (k=0; k<rows; k++) Phi2->data[ j ] += *(B + k*M + j) * Phi[k];
    }

    /*
     *  Solve for C
     */
    C = gsl_vector_alloc( M );
resid = gsl_vector_alloc( M );
    gsl_linalg_QR_lssolve( A, tau, Phi2, C, resid );
//  gsl_linalg_cholesky_solve( A, Phi2, C );


    /*
     * Compute vector value for the query point using the function defined by C
     */
    DFI_FuncLin( q, C, v );


    free( B );
    gsl_vector_free( C );
    gsl_vector_free( Phi2 );
    gsl_vector_free( tau );
gsl_vector_free(resid);
    gsl_matrix_free( A );
    free( Phi );
    free( W );

    return( 1 );
    
}





/*
 * This routine takes a point and returns the magnetic field there. The
 * magnetic field is computed by interpolating from a large set of scattered
 * pre-defined data points. In principle, this is a trivial task, but the
 * obvious brute force method ends up being phenomenally slow if the dataset is
 * large. The main problem boils down to efficiently finding the required
 * number of nearest neighbors that we need to do the interpolation.
 * Classically, this problem is referred to as the "k Nearest Neighbor Problem"
 * (or the "kNN Problem"). Here we use an octree-based algorithm to radically
 * accelerate the search for nearest neighbors.
 *
 * We also use a novel interpolation scheme that attempts to enforce divB=0.
 *
 *	There is an insidious problem with using this approach on certain types of meshes!
 *
 *	For example, when using an Euler potential mesh, what happens is that
 *	for query points that are far from the earth, the field lines diverge apart.
 *	This leads to a situation where all the NNs we get, come from the same (or
 *	nearly the same) meridian. This means our "interpolation" actually turns into
 *	a very bad extrapolation! 
 *
 *	Some possible fixes:
 *		1) Try to fill in more points as a pre-processing step.
 *		2) Try to generate a better distribution of points to
 *	           interpolate with. E.g., we could find some number of
 *		   NNs at the real query point. Then also some additional ones 
 *		   at other meridians.
 *	
 *      Lets try 2) here...
 *
 */
int Lgm_B_FromScatteredData( Lgm_Vector *v, Lgm_Vector *B, Lgm_MagModelInfo *Info ) {

    int             k, Flag;
    Lgm_Vector      P, P_m, P_p, v_m, v_p, vsm, v_m_gsm, v_p_gsm, Bcdip;
    double          r, CosTheta, SinTheta, Phi, Phi_m, Phi_p;
    double	    NormalizedMaxDist, MaxDist2, w, w1, w2;
    Lgm_OctreeData  *kNN;


    if ( Lgm_Magnitude(v) < 3.0 ) {
        Lgm_B_cdip( v, B, Info );
        return(1);
    }

    /*
     *  Make sure Octree has been initialized
     */
    if ( Info->OctreeRoot == NULL ) return( OCTREE_IS_NULL );



    /*
     *  Allocate some space for the kNNs. Actually we will generate 3*k because 
     *  we are going to try the +/- meridians to get a better cloud of points.
     */
    kNN = (Lgm_OctreeData *)calloc( 2+3*Info->Octree_kNN_k, sizeof( Lgm_OctreeData ) );

    /*
     *  determine the two additional points around which to find kNNs
     */
    Lgm_Convert_Coords( v, &vsm, GSM_TO_SM, Info->c );
    r = Lgm_Magnitude( &vsm );
    CosTheta = vsm.z/r;
    SinTheta = sqrt( 1.0 - CosTheta*CosTheta );
    Phi = atan2( vsm.y, vsm.x );

    Phi_m = Phi - 1.0*2.0*M_PI/61.0; // hard-wired based on Sorin's mesh files -- he has 61 beta values or "meridians"
    v_m.x = r*SinTheta*cos( Phi_m );
    v_m.y = r*SinTheta*sin( Phi_m );
    v_m.z = r*CosTheta;

    Phi_p = Phi + 1.0*2.0*M_PI/61.0; // hard-wired based on Sorin's mesh files -- he has 61 beta values or "meridians"
    v_p.x = r*SinTheta*cos( Phi_p );
    v_p.y = r*SinTheta*sin( Phi_p );
    v_p.z = r*CosTheta;

    




    /*
     *  The positions in the octree have components that are all scaled to be in the range
     *  [0.0, 1.0]. So we need to scale our point as well.
     */
//printf( "v = <%g %g %g>\n", v->x, v->y, v->z );
    P.x = (v->x - Info->OctreeScaleMin)/Info->OctreeScaleDiff;
    P.y = (v->y - Info->OctreeScaleMin)/Info->OctreeScaleDiff;
    P.z = (v->z - Info->OctreeScaleMin)/Info->OctreeScaleDiff;

    Lgm_Convert_Coords( &v_m, &v_m_gsm, SM_TO_GSM, Info->c );
    P_m.x = (v_m_gsm.x - Info->OctreeScaleMin)/Info->OctreeScaleDiff;
    P_m.y = (v_m_gsm.y - Info->OctreeScaleMin)/Info->OctreeScaleDiff;
    P_m.z = (v_m_gsm.z - Info->OctreeScaleMin)/Info->OctreeScaleDiff;

    Lgm_Convert_Coords( &v_p, &v_p_gsm, SM_TO_GSM, Info->c );
    P_p.x = (v_p_gsm.x - Info->OctreeScaleMin)/Info->OctreeScaleDiff;
    P_p.y = (v_p_gsm.y - Info->OctreeScaleMin)/Info->OctreeScaleDiff;
    P_p.z = (v_p_gsm.z - Info->OctreeScaleMin)/Info->OctreeScaleDiff;


    NormalizedMaxDist = Info->Octree_kNN_MaxDist/Info->OctreeScaleDiff;
    MaxDist2 = NormalizedMaxDist*NormalizedMaxDist;
    if (MaxDist2 > 1.0) MaxDist2 = 1.0;

MaxDist2 = 1.0;


    /*
     * Given these points and the octree, find the "Octree_kNN_k" nearest neighbors.
     */
    int kgot = 0;
    Flag = Lgm_Octree_kNN( &P, Info->OctreeRoot, Info->Octree_kNN_k, &k, MaxDist2, kNN );
//printf("P = %g, %g %g   Flag = %d\n", P.x, P.y, P.z, Flag);
//    if ( Flag != OCTREE_KNN_SUCCESS ) return( Flag );
    kgot += k;

    Flag = Lgm_Octree_kNN( &P_m, Info->OctreeRoot, Info->Octree_kNN_k, &k, MaxDist2, kNN+kgot );
//printf("P_m = %g, %g %g   Flag = %d\n", P_m.x, P_m.y, P_m.z, Flag);
//printf("Flag = %d\n", Flag);
//    if ( Flag != OCTREE_KNN_SUCCESS ) return( Flag );
    kgot += k;
    

    Flag = Lgm_Octree_kNN( &P_p, Info->OctreeRoot, Info->Octree_kNN_k, &k, MaxDist2, kNN+kgot );
//printf("P_p = %g, %g %g   Flag = %d\n", P_p.x, P_p.y, P_p.z, Flag);
//printf("Flag = %d\n", Flag);
//    if ( Flag != OCTREE_KNN_SUCCESS ) return( Flag );
    kgot += k;










int i;
/*
FILE *fpout;
fpout = fopen("NN.dat", "w");
sx = sy = sz = 0.0;
*/
for (i=0; i<kgot; i++){

/*
	v_p_gsm.x = kNN[i].Position.x*Info->OctreeScaleDiff + Info->OctreeScaleMin;
	v_p_gsm.y = kNN[i].Position.y*Info->OctreeScaleDiff + Info->OctreeScaleMin;
	v_p_gsm.z = kNN[i].Position.z*Info->OctreeScaleDiff + Info->OctreeScaleMin;
	fprintf(fpout, "%g %g\n", v_p_gsm.x, v_p_gsm.z);
xx = P_m.x-kNN[i].Position.x;
yy = P_m.y-kNN[i].Position.y;
zz = P_m.z-kNN[i].Position.z;
NewDist2 = xx*xx + yy*yy + zz*zz;
	printf("%d Dist2 = %g, kNN[%d].Position = %g %g %g   kNN[i].B[%d] = %g %g %g\n", i, kNN[i].Dist2, i, kNN[i].Position.x, kNN[i].Position.y, kNN[i].Position.z, i, kNN[i].B.x, kNN[i].B.y, kNN[i].B.z);
*/
/*
	printf("%d Dist2 = %g, NewDist2 = %g kNN[%d].Position = %g %g %g   v = %g %g %g,  kNN[i].B[%d] = %g %g %g", i, kNN[i].Dist2, NewDist2, i, kNN[i].Position.x, kNN[i].Position.y, kNN[i].Position.z, v_p_gsm.x, v_p_gsm.y, v_p_gsm.z, i, kNN[i].B.x, kNN[i].B.y, kNN[i].B.z);
	if (v_p_gsm.z > 0.0) printf(" ****\n");
	else printf("\n");

	sx += kNN[i].B.x;
	sy += kNN[i].B.y;
	sz += kNN[i].B.z;
*/
}
/*
sx /= (double)kgot;
sy /= (double)kgot;
sz /= (double)kgot;
printf("-----------------------------------------------\n");
printf("Avg = %g %g %g\n", sx, sy, sz);
fclose(fpout);

*/


double xx, yy, zz, NewDist2;
xx = v->x-kNN[0].Position.x;
yy = v->y-kNN[0].Position.y;
zz = v->z-kNN[0].Position.z;
NewDist2 = xx*xx+yy*yy+zz*zz;

    /*
     * Interpolate the points
     */
    if ( Info->Octree_kNN_InterpMethod == LINEAR_DFI ){

//printf("kgot = %d\n", kgot);
        Lgm_DivFreeInterp2( &P, kNN, kgot, B );

        //if (kgot == 3*Info->Octree_kNN_k) {
        //    Lgm_DivFreeInterp2( &P, kNN, kgot, B );
        //} else if (kgot > 2*Info->Octree_kNN_k ){
//        if ( NewDist2 < .25 ) {
//            Lgm_DivFreeInterp2( &P, kNN, kgot, B );
//        } else {//if (kgot > 2*Info->Octree_kNN_k ){
// add a boundary point
//kNN[kgot].Position.x = kNN[0].Position.x;
//kNN[kgot].Position.y = kNN[0].Position.y;
//kNN[kgot].Position.z = 10.0;
//kNN[kgot].B.x = 0.0;
//kNN[kgot].B.y = 0.0;
//kNN[kgot].B.z = -1.0;
//kgot++;
/*
kNN[kgot].Position.x = kNN[0].Position.x;
kNN[kgot].Position.y = kNN[0].Position.y;
kNN[kgot].Position.z = -10.0;
kNN[kgot].B.x = 0.0;
kNN[kgot].B.y = 0.0;
kNN[kgot].B.z = -1.0;
kgot++;
*/
//            Lgm_DivFreeInterp2( &P, kNN, kgot, B );
/*
            w  = (double)(3*Info->Octree_kNN_k);
            w1 = (double)kgot;
            w2 = (double)(3*Info->Octree_kNN_k - kgot);
            Lgm_B_cdip( v, &Bcdip, Info );
            B->x = (w1*B->x + w2*Bcdip.x)/w;
            B->y = (w1*B->y + w2*Bcdip.y)/w;
            B->z = (w1*B->z + w2*Bcdip.z)/w;
*/
//        } /*else {
//            Lgm_B_cdip( v, B, Info );
//        }*/
//printf("B = %g %g %g\n", B->x, B->y, B->z);
        

    } else if ( Info->Octree_kNN_InterpMethod == QUADRATIC_DFI ) { 

        Lgm_DivFreeInterp( &P, kNN, kgot, B );

    } else if ( Info->Octree_kNN_InterpMethod == NEWTON_INTERP ) {

        Lgm_Vector  f0, f1, f2, f3, A1, A2, A3;
        double	    dx, dy, dz, dx2, dy2, dz2, D10, D12, D20, D30, D31, D32;
        double	    dv0, dv1, dv2, dv3;
        
        f0.x = kNN[0].B.x; f0.y = kNN[0].B.y; f0.z = kNN[0].B.z;
        f1.x = kNN[1].B.x; f1.y = kNN[1].B.y; f1.z = kNN[1].B.z;
        f2.x = kNN[2].B.x; f2.y = kNN[2].B.y; f2.z = kNN[2].B.z;
        f3.x = kNN[3].B.x; f3.y = kNN[3].B.y; f3.z = kNN[3].B.z;

        dx = kNN[1].Position.x - kNN[0].Position.x; dx2 = dx*dx; 
        dy = kNN[1].Position.y - kNN[0].Position.y; dy2 = dy*dy; 
        dz = kNN[1].Position.z - kNN[0].Position.z; dz2 = dz*dz;
        D10 = sqrt( dx2 + dy2 +dz2 );

        dx = kNN[1].Position.x - kNN[2].Position.x; dx2 = dx*dx; 
        dy = kNN[1].Position.y - kNN[2].Position.y; dy2 = dy*dy; 
        dz = kNN[1].Position.z - kNN[2].Position.z; dz2 = dz*dz;
        D12 = sqrt( dx2 + dy2 +dz2 );
        
        dx = kNN[2].Position.x - kNN[0].Position.x; dx2 = dx*dx; 
        dy = kNN[2].Position.y - kNN[0].Position.y; dy2 = dy*dy; 
        dz = kNN[2].Position.z - kNN[0].Position.z; dz2 = dz*dz;
        D20 = sqrt( dx2 + dy2 +dz2 );
        
        dx = kNN[3].Position.x - kNN[0].Position.x; dx2 = dx*dx; 
        dy = kNN[3].Position.y - kNN[0].Position.y; dy2 = dy*dy; 
        dz = kNN[3].Position.z - kNN[0].Position.z; dz2 = dz*dz;
        D30 = sqrt( dx2 + dy2 +dz2 );
        
        dx = kNN[3].Position.x - kNN[1].Position.x; dx2 = dx*dx; 
        dy = kNN[3].Position.y - kNN[1].Position.y; dy2 = dy*dy; 
        dz = kNN[3].Position.z - kNN[1].Position.z; dz2 = dz*dz;
        D31 = sqrt( dx2 + dy2 +dz2 );
        
        dx = kNN[3].Position.x - kNN[2].Position.x; dx2 = dx*dx; 
        dy = kNN[3].Position.y - kNN[2].Position.y; dy2 = dy*dy; 
        dz = kNN[3].Position.z - kNN[2].Position.z; dz2 = dz*dz;
        D32 = sqrt( dx2 + dy2 +dz2 );

        A1.x = (f1.x-f0.x)/D10; 
        A1.y = (f1.y-f0.y)/D10; 
        A1.z = (f1.z-f0.z)/D10;

        A2.x = ((f2.x-f0.x)/D20 + (f0.x-f1.x)/D10)/D12; 
        A2.y = ((f2.y-f0.y)/D20 + (f0.y-f1.y)/D10)/D12; 
        A2.z = ((f2.z-f0.z)/D20 + (f0.z-f1.z)/D10)/D12;

        A3.x = (D12*((f3.x-f0.x)/D30+(f0.x-f1.x)/D10)+D31*((f1.x-f0.x)/D10+(f0.x-f2.x)/D20))/D12/D31/D32;
        A3.y = (D12*((f3.y-f0.y)/D30+(f0.y-f1.y)/D10)+D31*((f1.y-f0.y)/D10+(f0.y-f2.y)/D20))/D12/D31/D32;
        A3.z = (D12*((f3.z-f0.z)/D30+(f0.z-f1.z)/D10)+D31*((f1.z-f0.z)/D10+(f0.z-f2.z)/D20))/D12/D31/D32;

        dv0 = sqrt( kNN[0].Dist2 );
        dv1 = sqrt( kNN[1].Dist2 );
        dv2 = sqrt( kNN[2].Dist2 );
        dv3 = sqrt( kNN[3].Dist2 );

        B->x = f0.x + A1.x*dv0 + A2.x*dv0*dv1 + A3.x*dv0*dv1*dv2;
        B->y = f0.y + A1.y*dv0 + A2.y*dv0*dv1 + A3.y*dv0*dv1*dv2;
        B->z = f0.z + A1.z*dv0 + A2.z*dv0*dv1 + A3.z*dv0*dv1*dv2;

    }

    free( kNN );


    ++Info->nFunc;

    return( 1 );

}


/*
 *    $Id: B_FromScatteredData.c 45 2010-10-01 20:43:29Z mgh $
 */
