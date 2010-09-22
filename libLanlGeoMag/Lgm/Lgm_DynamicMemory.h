#ifndef LGM_DYNAMIC_MEMORY_H
#define LGM_DYNAMIC_MEMORY_H

/**
 *  Macros for dynamically defining multi-dimensional arrays. There are three
 *  sets of macros here:
 *          
 *      1) Macros that allocate the contiguous block of memory for you (via
 *         calloc - so memory is initialized to zeros). These are:
 *
 *              LGM_ARRAY_1D( A, n1, type );
 *              LGM_ARRAY_2D( A, n1, n2, type );
 *              LGM_ARRAY_3D( A, n1, n2, n3, type );
 *              LGM_ARRAY_4D( A, n1, n2, n3, n4, type );
 *              LGM_ARRAY_5D( A, n1, n2, n3, n4, n5, type );
 *
 *          A typical use would be:
 *
 *              double ***A;
 *              int       n1, n2, n3;
 *
 *              n1 = 5; n2 = 6; n3 = 523;
 *              LGM_ARRAY_3D( A, n1, n2, n3, double );
 *                  ... do stuff ...
 *              LGM_ARRAY_3D_FREE( A );
 *
 *          Note that because this is a macro rather than a function, ANY data
 *          type can be used including user-defined objects like structures.
 *          E.g., we could make a 2D array of Lgm_Vector's as follows;
 *
 *              Lgm_Vector **A;
 *              int          n1, n2;
 *
 *              n1 = 51; n2 = 62;
 *              LGM_ARRAY_2D( A, n1, n2, Lgm_Vector );
 *                  ... do stuff ...
 *              LGM_ARRAY_2D_FREE( A );
 *
 *      2) Macros to free the arrays. We saw examples of there use above. The Macros are:
 *
 *              LGM_ARRAY_1D_FREE( A );
 *              LGM_ARRAY_2D_FREE( A );
 *              LGM_ARRAY_3D_FREE( A );
 *              LGM_ARRAY_4D_FREE( A );
 *              LGM_ARRAY_5D_FREE( A );
 *
 *      3) Macros that define the array structure for you, but require the user
 *         to provide a pointer to an already allocated contiguous block of
 *         memory. (Use the same FREE macros).
 *
 *              LGM_ARRAY_FROM_DATA_1D( A, n1, type );
 *              LGM_ARRAY_FROM_DATA_2D( A, n1, n2, type );
 *              LGM_ARRAY_FROM_DATA_3D( A, n1, n2, n3, type );
 *              LGM_ARRAY_FROM_DATA_4D( A, n1, n2, n3, n4, type );
 *              LGM_ARRAY_FROM_DATA_5D( A, n1, n2, n3, n4, n5, type );
 *
 *
 *          A typical use would be:
 *
 *              double ***A, *pdata;
 *              int       n1, n2, n3;
 *
 *              n1 = 5; n2 = 6; n3 = 523;
 *              pdata = (double *)calloc( n1*n2*n3, sizeof( double ) );
 *              LGM_ARRAY_3D( A, pdata, n1, n2, n3, double );
 *                  ... do stuff ...
 *              LGM_ARRAY_3D_FREE( A );
 *
 *          Why would we need this? One use would be to make it easier to
 *          read/write dynamic arrays to binary files.  For example, if the
 *          arrays were defined statically, its trivial to do atomic writes and
 *          reads:
 *          
 *              double  A[4][5];
 *              fd = open( "file.dat", O_CREAT|O_WRONLY, 0664 );
 *              write( fd, A, 4*5*sizeof( double ) );
 *              close(fd);
 *
 *              fd = open( "file.dat", O_RDONLY );
 *              read( fd, A, 4*5*sizeof( double ) );
 *              close(fd);
 *
 *          But this cant be done when we define them dynamically. But using
 *          these macros we could do something like:
 *
 *              To write:
 *                  int       n1, n2;
 *                  double  **A;
 *                  n1 = 4; n2 = 5;
 *                  LGM_ARRAY_2D( A, n1, n2, double );
 *                  fd = open( "file.dat", O_CREAT|O_WRONLY, 0664 );
 *                  write( fd, &n1, sizeof( n1 ) );
 *                  write( fd, &n2, sizeof( n2 ) );
 *                  write( fd, A, n1*n2*sizeof( double ) );
 *                  close(fd);
 *
 *
 *              To read:
 *                  int       n1, n2;
 *                  double  **A, *pdata;
 *                  fd = open( "file.dat", O_RDONLY );
 *                  read( fd, &n1, sizeof( n1 ) );
 *                  read( fd, &n2, sizeof( n2 ) );
 *                  pdata = (double *)calloc( n1*n2, sizeof( double ) );
 *                  read( fd, pdata, n1*n2*sizeof( double ) );
 *                  close(fd);
 *                  LGM_ARRAY_FROM_DATA_2D( A, pdata, n1, n2, double );
 *              
 *
 */


#include <stdlib.h>
#include <stdio.h>
#if !defined(__APPLE__)
#include <malloc.h>
#endif

/*
 * Macros for dynamically allocating/freeing 1D arrays of any type
 */
#define LGM_ARRAY_1D( prow, col, type ) {\
    register type *pdata;\
    pdata = (type *)calloc( col, sizeof( type ) );\
    if ( pdata == (type *)NULL ) {\
        fprintf(stderr, "LGM_ARRAY_2D Macro: Could not allocate space for data\n");\
        exit(1);\
    }\
    prow = pdata;\
}\

#define LGM_ARRAY_1D_FREE( pdata ) {\
    free( pdata );\
}\


/*
 * Macros for dynamically allocating/freeing 2D arrays of any type
 */
#define LGM_ARRAY_2D( prow, row, col, type ) {\
    register type *pdata;\
    int      i;\
    pdata = (type *)calloc( row*col, sizeof( type ) );\
    if ( pdata == (type *)NULL ) {\
        fprintf(stderr, "LGM_ARRAY_2D Macro: Could not allocate space for data\n");\
        exit(1);\
    }\
    prow = (type **)calloc( row, sizeof( type * ));\
    if ( prow == (type **)NULL ) {\
        fprintf(stderr, "LGM_ARRAY_2D Macro: Could not allocate space for row pointers\n");\
        exit(1);\
    }\
    for (i=0; i<row; i++){\
        prow[i] = pdata;\
        pdata += col;\
    }\
}\

#define LGM_ARRAY_2D_FREE( prow ) {\
    free( *prow );\
    free( prow );\
}\




/*
 * Macros for dynamically allocating/freeing 3D arrays of any type
 */
#define LGM_ARRAY_3D( pgrid, grid, row, col, type ) {\
    register type **prow, *pdata;\
    int      i;\
    pdata = (type *)calloc( grid*row*col, sizeof( type ) );\
    if ( pdata == (type *)NULL ) {\
        fprintf(stderr, "LGM_ARRAY_3D Macro: Could not allocate space for data\n");\
        exit(1);\
    }\
    prow = (type **)calloc( grid*row, sizeof( type * ));\
    if ( prow == (type **)NULL ) {\
        fprintf(stderr, "LGM_ARRAY_3D Macro: Could not allocate space for row pointers\n");\
        exit(1);\
    }\
    pgrid = (type ***)calloc( grid, sizeof( type * ));\
    if ( pgrid == (type ***)NULL ) {\
        fprintf(stderr, "LGM_ARRAY_3D Macro: Could not allocate space for row pointers\n");\
        exit(1);\
    }\
    for (i=0; i<grid*row; i++){\
        prow[i] = pdata;\
        pdata += col;\
    }\
    for (i=0; i<grid; i++){\
        pgrid[i] = prow;\
        prow += row;\
    }\
}\
    
#define LGM_ARRAY_3D_FREE( pa ) {\
    free( **pa );\
    free( *pa );\
    free( pa );\
}\





/*
 * Macros for dynamically allocating/freeing 4D arrays of any type
 */
#define LGM_ARRAY_4D( pn4, n4, n3, n2, n1, type ) {\
    long int i;\
    register type ***pn3, **pn2, *pdata;\
    pdata = (type *)calloc( n4*n3*n2*n1, sizeof( type ) );\
    if ( pdata == (type *)NULL ) {\
        fprintf(stderr, "LGM_ARRAY_4D Macro: Could not allocate space for data\n");\
        exit(1);\
    }\
    pn2 = (type **)calloc( n4*n3*n2, sizeof( type * ));\
    if ( pn2 == (type **)NULL ) {\
        fprintf(stderr, "LGM_ARRAY_4D Macro: Could not allocate space for n2 pointers\n");\
        exit(1);\
    }\
    pn3 = (type ***)calloc( n4*n3, sizeof( type * ));\
    if ( pn3 == (type ***)NULL ) {\
        fprintf(stderr, "LGM_ARRAY_4D Macro: Could not allocate space for n3 pointers\n");\
        exit(1);\
    }\
    pn4 = (type ****)calloc( n4, sizeof( type * ));\
    if ( pn4 == (type ****)NULL ) {\
        fprintf(stderr, "LGM_ARRAY_4D Macro: Could not allocate space for n4 pointers\n");\
        exit(1);\
    }\
    for (i=0; i<n4*n3*n2; i++){\
        pn2[i] = pdata;\
        pdata += n1;\
    }\
    for (i=0; i<n4*n3; i++){\
        pn3[i] = pn2;\
        pn2 += n2;\
    }\
    for (i=0; i<n4; i++){\
        pn4[i] = pn3;\
        pn3 += n3;\
    }\
    pn4 = &pn4[0];\
}\
    
#define LGM_ARRAY_4D_FREE( pn4 ) {\
    free( ***pn4 );\
    free( **pn4 );\
    free( *pn4 );\
    free( pn4 );\
}\





/*
 * Macros for dynamically allocating/freeing 5D arrays of any type
 */
#define LGM_ARRAY_5D( pn5, n5, n4, n3, n2, n1, type ) {\
    long int i;\
    register type ****pn4, ***pn3, **pn2, *pdata;\
    pdata = (type *)calloc( n5*n4*n3*n2*n1, sizeof( type ) );\
    if ( pdata == (type *)NULL ) {\
        fprintf(stderr, "LGM_ARRAY_5D Macro: Could not allocate space for data\n");\
        exit(1);\
    }\
    pn2 = (type **)calloc( n5*n4*n3*n2, sizeof( type * ));\
    if ( pn2 == (type **)NULL ) {\
        fprintf(stderr, "LGM_ARRAY_5D Macro: Could not allocate space for n2 pointers\n");\
        exit(1);\
    }\
    pn3 = (type ***)calloc( n5*n4*n3, sizeof( type * ));\
    if ( pn3 == (type ***)NULL ) {\
        fprintf(stderr, "LGM_ARRAY_5D Macro: Could not allocate space for n3 pointers\n");\
        exit(1);\
    }\
    pn4 = (type ****)calloc( n5*n4, sizeof( type * ));\
    if ( pn4 == (type ****)NULL ) {\
        fprintf(stderr, "LGM_ARRAY_5D Macro: Could not allocate space for n4 pointers\n");\
        exit(1);\
    }\
    pn5 = (type *****)calloc( n5, sizeof( type * ));\
    if ( pn5 == (type *****)NULL ) {\
        fprintf(stderr, "LGM_ARRAY_5D Macro: Could not allocate space for n4 pointers\n");\
        exit(1);\
    }\
    for (i=0; i<n5*n4*n3*n2; i++){\
        pn2[i] = pdata;\
        pdata += n1;\
    }\
    for (i=0; i<n5*n4*n3; i++){\
        pn3[i] = pn2;\
        pn2 += n2;\
    }\
    for (i=0; i<n5*n4; i++){\
        pn4[i] = pn3;\
        pn3 += n3;\
    }\
    for (i=0; i<n5; i++){\
        pn5[i] = pn4;\
        pn4 += n4;\
    }\
}\


#define LGM_ARRAY_5D_FREE( pn5 ) {\
    free( ****pn5 );\
    free( ***pn5 );\
    free( **pn5 );\
    free( *pn5 );\
    free( pn5 );\
}\





/*
 *   More macros for dynamically allocating/freeing 1D arrays of any type.
 *   These versions require the user to provide a properly calloc'd block of
 *   data.  Why would we need this? One case would be to simplify unformated
 *   read/write of arrays to/from binary files. 
 */

#define LGM_ARRAY_FROM_DATA_1D( prow, pdata, col, type ) {\
    if ( pdata == (type *)NULL ) {\
        fprintf(stderr, "LGM_ARRAY_FROM_DATA_1D Macro: pdata is NULL\n");\
        exit(1);\
    }\
    prow = pdata;\
}\


/*
 * Macros for dynamically allocating/freeing 2D arrays of any type
 */
#define LGM_ARRAY_FROM_DATA_2D( prow, pdata, row, col, type ) {\
    int      i;\
    if ( pdata == (type *)NULL ) {\
        fprintf(stderr, "LGM_ARRAY_FROM_DATA_2D Macro: pdata is NULL\n");\
        exit(1);\
    }\
    prow = (type **)calloc( row, sizeof( type * ));\
    if ( prow == (type **)NULL ) {\
        fprintf(stderr, "LGM_ARRAY_FROM_DATA_2D Macro: pdata is NULL\n");\
        exit(1);\
    }\
    for (i=0; i<row; i++){\
        prow[i] = pdata;\
        pdata += col;\
    }\
}\


/*
 * Macros for dynamically allocating/freeing 3D arrays of any type
 */
#define LGM_ARRAY_FROM_DATA_3D( pgrid, pdata, grid, row, col, type ) {\
    register type **prow;\
    int      i;\
    if ( pdata == (type *)NULL ) {\
        fprintf(stderr, "LGM_ARRAY_FROM_DATA_3D Macro: Could not allocate space for data\n");\
        exit(1);\
    }\
    prow = (type **)calloc( grid*row, sizeof( type * ));\
    if ( prow == (type **)NULL ) {\
        fprintf(stderr, "LGM_ARRAY_FROM_DATA_3D Macro: Could not allocate space for row pointers\n");\
        exit(1);\
    }\
    pgrid = (type ***)calloc( grid, sizeof( type * ));\
    if ( pgrid == (type ***)NULL ) {\
        fprintf(stderr, "LGM_ARRAY_FROM_DATA_3D Macro: Could not allocate space for row pointers\n");\
        exit(1);\
    }\
    for (i=0; i<grid*row; i++){\
        prow[i] = pdata;\
        pdata += col;\
    }\
    for (i=0; i<grid; i++){\
        pgrid[i] = prow;\
        prow += row;\
    }\
}\
    



/*
 * Macros for dynamically allocating/freeing 4D arrays of any type
 */
#define LGM_ARRAY_FROM_DATA_4D( pn4, pdata, n4, n3, n2, n1, type ) {\
    long int i;\
    register type ***pn3, **pn2;\
    if ( pdata == (type *)NULL ) {\
        fprintf(stderr, "LGM_ARRAY_FROM_DATA_4D Macro: Could not allocate space for data\n");\
        exit(1);\
    }\
    pn2 = (type **)calloc( n4*n3*n2, sizeof( type * ));\
    if ( pn2 == (type **)NULL ) {\
        fprintf(stderr, "LGM_ARRAY_FROM_DATA_4D Macro: Could not allocate space for n2 pointers\n");\
        exit(1);\
    }\
    pn3 = (type ***)calloc( n4*n3, sizeof( type * ));\
    if ( pn3 == (type ***)NULL ) {\
        fprintf(stderr, "LGM_ARRAY_FROM_DATA_4D Macro: Could not allocate space for n3 pointers\n");\
        exit(1);\
    }\
    pn4 = (type ****)calloc( n4, sizeof( type * ));\
    if ( pn4 == (type ****)NULL ) {\
        fprintf(stderr, "LGM_ARRAY_FROM_DATA_4D Macro: Could not allocate space for n4 pointers\n");\
        exit(1);\
    }\
    for (i=0; i<n4*n3*n2; i++){\
        pn2[i] = pdata;\
        pdata += n1;\
    }\
    for (i=0; i<n4*n3; i++){\
        pn3[i] = pn2;\
        pn2 += n2;\
    }\
    for (i=0; i<n4; i++){\
        pn4[i] = pn3;\
        pn3 += n3;\
    }\
    pn4 = &pn4[0];\
}\
    



/*
 * Macro for dynamically allocating/freeing 5D arrays of any type User must
 * provide the data block. (e.g. pdata = (type *)calloc( n5*n4*n3*n2*n1, sizeof( type ) ); )
 */
#define LGM_ARRAY_FROM_DATA_5D( pn5, pdata, n5, n4, n3, n2, n1, type ) {\
    long int i;\
    register type ****pn4, ***pn3, **pn2;\
    if ( pdata == (type *)NULL ) {\
        fprintf(stderr, "LGM_ARRAY_FROM_DATA_5D Macro: Could not allocate space for data\n");\
        exit(1);\
    }\
    pn2 = (type **)calloc( n5*n4*n3*n2, sizeof( type * ));\
    if ( pn2 == (type **)NULL ) {\
        fprintf(stderr, "LGM_ARRAY_FROM_DATA_5D Macro: Could not allocate space for n2 pointers\n");\
        exit(1);\
    }\
    pn3 = (type ***)calloc( n5*n4*n3, sizeof( type * ));\
    if ( pn3 == (type ***)NULL ) {\
        fprintf(stderr, "LGM_ARRAY_FROM_DATA_5D Macro: Could not allocate space for n3 pointers\n");\
        exit(1);\
    }\
    pn4 = (type ****)calloc( n5*n4, sizeof( type * ));\
    if ( pn4 == (type ****)NULL ) {\
        fprintf(stderr, "LGM_ARRAY_FROM_DATA_5D Macro: Could not allocate space for n4 pointers\n");\
        exit(1);\
    }\
    pn5 = (type *****)calloc( n5, sizeof( type * ));\
    if ( pn5 == (type *****)NULL ) {\
        fprintf(stderr, "LGM_ARRAY_FROM_DATA_5D Macro: Could not allocate space for n4 pointers\n");\
        exit(1);\
    }\
    for (i=0; i<n5*n4*n3*n2; i++){\
        pn2[i] = pdata;\
        pdata += n1;\
    }\
    for (i=0; i<n5*n4*n3; i++){\
        pn3[i] = pn2;\
        pn2 += n2;\
    }\
    for (i=0; i<n5*n4; i++){\
        pn4[i] = pn3;\
        pn3 += n3;\
    }\
    for (i=0; i<n5; i++){\
        pn5[i] = pn4;\
        pn4 += n4;\
    }\
}\





#endif
