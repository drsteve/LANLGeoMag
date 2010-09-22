#include <stdlib.h>
#include <stdio.h>
/* #if !defined(__APPLE__) */
/* #include <malloc.h> */
/* #endif */

/*
 * Macros for dynamically allocating/freeing 1D arrays of any type
 */
#define ARRAY_1D( prow, col, type ) {\
\
    register type *pdata;\
\
    pdata = (type *) calloc( col, sizeof( type ) );\
    if ( pdata == (type *) NULL ) {\
        fprintf(stderr, "ARRAY_2D Macro: Could not allocate space for data\n");\
        exit(1);\
    }\
    prow = pdata;\
\
}\

#define ARRAY_1D_FREE( pdata ) {\
\
    free( pdata );\
\
}\


/*
 * Macros for dynamically allocating/freeing 2D arrays of any type
 */
#define ARRAY_2D( prow, row, col, type ) {\
\
    register type *pdata;\
    int      i;\
\
    pdata = (type *) calloc( row*col, sizeof( type ) );\
    if ( pdata == (type *) NULL ) {\
        fprintf(stderr, "ARRAY_2D Macro: Could not allocate space for data\n");\
        exit(1);\
    }\
\
    prow = (type **) calloc( row, sizeof( type * ));\
    if ( prow == (type **) NULL ) {\
        fprintf(stderr, "ARRAY_2D Macro: Could not allocate space for row pointers\n");\
        exit(1);\
    }\
\
    for (i=0; i<row; i++){\
        prow[i] = pdata;\
        pdata += col;\
    }\
\
}\

#define ARRAY_2D_FREE( prow ) {\
\
    free( *prow );\
    free( prow );\
\
}\




/*
 * Macros for dynamically allocating/freeing 3D arrays of any type
 */
#define ARRAY_3D( pgrid, grid, row, col, type ) {\
\
    register type **prow, *pdata;\
    int      i;\
\
    pdata = (type *) calloc( grid*row*col, sizeof( type ) );\
    if ( pdata == (type *) NULL ) {\
        fprintf(stderr, "ARRAY_3D Macro: Could not allocate space for data\n");\
        exit(1);\
    }\
\
    prow = (type **) calloc( grid*row, sizeof( type * ));\
    if ( prow == (type **) NULL ) {\
        fprintf(stderr, "ARRAY_3D Macro: Could not allocate space for row pointers\n");\
        exit(1);\
    }\
\
    pgrid = (type ***) calloc( grid, sizeof( type * ));\
    if ( pgrid == (type ***) NULL ) {\
        fprintf(stderr, "ARRAY_3D Macro: Could not allocate space for row pointers\n");\
        exit(1);\
    }\
\
    for (i=0; i<grid*row; i++){\
        prow[i] = pdata;\
        pdata += col;\
    }\
\
    for (i=0; i<grid; i++){\
        pgrid[i] = prow;\
        prow += row;\
    }\
\
}\
    
#define ARRAY_3D_FREE( pa ) {\
\
    free( **pa );\
    free( *pa );\
    free( pa );\
\
}\





/*
 * Macros for dynamically allocating/freeing 4D arrays of any type
 */
#define ARRAY_4D( pn4, n4, n3, n2, n1, type ) {\
\
    long int i;\
    register type ***pn3, **pn2, *pdata;\
\
    pdata = (type *) calloc( n4*n3*n2*n1, sizeof( type ) );\
    if ( pdata == (type *) NULL ) {\
        fprintf(stderr, "ARRAY_4D Macro: Could not allocate space for data\n");\
        exit(1);\
    }\
\
    pn2 = (type **) calloc( n4*n3*n2, sizeof( type * ));\
    if ( pn2 == (type **) NULL ) {\
        fprintf(stderr, "ARRAY_4D Macro: Could not allocate space for n2 pointers\n");\
        exit(1);\
    }\
\
    pn3 = (type ***) calloc( n4*n3, sizeof( type * ));\
    if ( pn3 == (type ***) NULL ) {\
        fprintf(stderr, "ARRAY_4D Macro: Could not allocate space for n3 pointers\n");\
        exit(1);\
    }\
\
    pn4 = (type ****) calloc( n4, sizeof( type * ));\
    if ( pn4 == (type ****) NULL ) {\
        fprintf(stderr, "ARRAY_4D Macro: Could not allocate space for n4 pointers\n");\
        exit(1);\
    }\
\
    for (i=0; i<n4*n3*n2; i++){\
        pn2[i] = pdata;\
        pdata += n1;\
    }\
\
    for (i=0; i<n4*n3; i++){\
        pn3[i] = pn2;\
        pn2 += n2;\
    }\
\
    for (i=0; i<n4; i++){\
        pn4[i] = pn3;\
        pn3 += n3;\
    }\
\
    pn4 = &pn4[0];\
\
}\
    
#define ARRAY_4D_FREE( pn4 ) {\
\
    free( ***pn4 );\
    free( **pn4 );\
    free( *pn4 );\
    free( pn4 );\
\
}\





/*
 * Macros for dynamically allocating/freeing 5D arrays of any type
 */
#define ARRAY_5D( pn5, n5, n4, n3, n2, n1, type ) {\
\
    long int i;\
    register type ****pn4, ***pn3, **pn2, *pdata;\
\
    pdata = (type *) calloc( n5*n4*n3*n2*n1, sizeof( type ) );\
    if ( pdata == (type *) NULL ) {\
        fprintf(stderr, "ARRAY_5D Macro: Could not allocate space for data\n");\
        exit(1);\
    }\
\
    pn2 = (type **) calloc( n5*n4*n3*n2, sizeof( type * ));\
    if ( pn2 == (type **) NULL ) {\
        fprintf(stderr, "ARRAY_5D Macro: Could not allocate space for n2 pointers\n");\
        exit(1);\
    }\
\
    pn3 = (type ***) calloc( n5*n4*n3, sizeof( type * ));\
    if ( pn3 == (type ***) NULL ) {\
        fprintf(stderr, "ARRAY_5D Macro: Could not allocate space for n3 pointers\n");\
        exit(1);\
    }\
\
    pn4 = (type ****) calloc( n5*n4, sizeof( type * ));\
    if ( pn4 == (type ****) NULL ) {\
        fprintf(stderr, "ARRAY_5D Macro: Could not allocate space for n4 pointers\n");\
        exit(1);\
    }\
\
    pn5 = (type *****) calloc( n5, sizeof( type * ));\
    if ( pn5 == (type *****) NULL ) {\
        fprintf(stderr, "ARRAY_5D Macro: Could not allocate space for n4 pointers\n");\
        exit(1);\
    }\
\
    for (i=0; i<n5*n4*n3*n2; i++){\
        pn2[i] = pdata;\
        pdata += n1;\
    }\
\
    for (i=0; i<n5*n4*n3; i++){\
        pn3[i] = pn2;\
        pn2 += n2;\
    }\
\
    for (i=0; i<n5*n4; i++){\
        pn4[i] = pn3;\
        pn3 += n3;\
    }\
\
    for (i=0; i<n5; i++){\
        pn5[i] = pn4;\
        pn4 += n4;\
    }\
\
}\


#define ARRAY_5D_FREE( pn5 ) {\
\
    free( ****pn5 );\
    free( ***pn5 );\
    free( **pn5 );\
    free( *pn5 );\
    free( pn5 );\
\
}\





