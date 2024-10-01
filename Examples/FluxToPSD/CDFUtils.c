#include <argp.h>
#include <cdf.h>
#include <gsl/gsl_spline.h>
#include <Lgm/Lgm_CTrans.h>
#include <Lgm/Lgm_DynamicMemory.h>
#include <Lgm/Lgm_Quat.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "CDFUtils.h"


/*
 * Status handler.
 */
void StatusHandler( CDFstatus status ) {

    char message[CDF_STATUSTEXT_LEN+1];

    if (status < CDF_WARN) {
      printf ("An error has occurred, halting...\n");
      CDFgetStatusText (status, message);
      printf ("%s\n", message);
      exit (1);
    } else if (status < CDF_OK) {
        printf ("Warning, function may not have completed as expected...\n");
        CDFgetStatusText (status, message);
        printf ("%s\n", message);
    } else if (status > CDF_OK) {
          printf ("Function completed successfully, but be advised that...\n");
          CDFgetStatusText (status, message);
          printf ("%s\n", message);
    }
    return;
}

/*
 * Dynamic allocation error handler.
 */
void DynamicErrorHandle( CDFid id ) {

    CDFstatus status;
    printf ("An error has occurred while doing malloc in CDF calls, halting...\n");
    status = CDFcloseCDF (id);
    if (status != CDF_OK) StatusHandler (status);

    return;

}

void DataTypeToStr( long dataType, char *str ) {

    switch( dataType ) {

        case 1:
                sprintf(str, "CDF_INT1");
                break;

        case 2:
                sprintf(str, "CDF_INT2");
                break;

        case 4:
                sprintf(str, "CDF_INT4");
                break;

        case 8:
                sprintf(str, "CDF_INT8");
                break;

        case 11:
                sprintf(str, "CDF_UINT1");
                break;

        case 12:
                sprintf(str, "CDF_UINT2");
                break;

        case 14:
                sprintf(str, "CDF_UINT4");
                break;

        case 21:
                sprintf(str, "CDF_REAL4");
                break;

        case 22:
                sprintf(str, "CDF_REAL8");
                break;

        case 31:
                sprintf(str, "CDF_EPOCH");
                break;

        case 32:
                sprintf(str, "CDF_EPOCH16");
                break;

        case 33:
                sprintf(str, "CDF_TIME_TT2000");
                break;

        case 41:
                sprintf(str, "CDF_BYTE");
                break;

        case 44:
                sprintf(str, "CDF_FLOAT");
                break;

        case 45:
                sprintf(str, "CDF_DOUBLE");
                break;

        case 51:
                sprintf(str, "CDF_CHAR");
                break;

        case 52:
                sprintf(str, "CDF_UCHAR");
                break;



        default:
                sprintf(str, "UNKNOWN");
                break;

    }


    return;
}


void CDF_Info( char *Filename, int verbosity ){

    int         ii, i, j, q, strsize, flag;
    CDFid       id;
    CDFstatus   status;
    char        Str[4096];
    char        attrName[CDF_ATTR_NAME_LEN256+1], subincrement, copyright[CDF_COPYRIGHT_LEN+1];
    long        version, release, increment, numDims, dimSizes[CDF_MAX_DIMS];
    long        datatype, numElements, numRecs, attrScope, maxgEntry, maxrEntry, maxzEntry;
    long        encoding, majority, maxrRec, numrVars, maxzRec, numzVars, numAttrs;
    char        varName[CDF_VAR_NAME_LEN256+1];
    long        dataType, numElems, recVary, dimVarys[CDF_MAX_DIMS];
    long        ilen;
    void        *buffer = NULL;


    /*
     * Get the current CDF library version number.
     */
    status = CDFgetLibraryVersion( &version, &release, &increment, &subincrement );
    if ( status != CDF_OK ) StatusHandler( status );

    status = CDFopenCDF( Filename, &id );
    if ( status != CDF_OK ) StatusHandler( status );

    status = CDFgetCopyright( id, copyright );
    if (status != CDF_OK) StatusHandler( status );

    /* Get the CDF version that was used to create this CDF file */
    status = CDFgetVersion( id, &version, &release, &increment );
    if (status != CDF_OK) StatusHandler( status );

    status = CDFinquireCDF (id,
                            &numDims,  dimSizes,        /* only good for rVars */
                            &encoding, &majority,
                            &maxrRec,  &numrVars,       /* only good for rVars */
                            &maxzRec,  &numzVars, &numAttrs);
    if (status != CDF_OK) StatusHandler( status );


    printf( "\t\tInformation for CDF File: %s\n", Filename );
    printf( "\t\t----------------------------------------------------------------------------------------------\n");
    printf( "\t\t                           CDF file name: %s\n", Filename );
    printf( "\t\t                     CDF library version: CDF %ld.%ld.%ld\n", version, release, increment );
    //printf( "\t\t                           CDF copyright:\n\t%s\n", copyright );
    printf( "\t\t                     Number of variables: %ld rVars, %ld zVars\n", numrVars, numzVars );
    printf( "\t\tNumber of attributes (global & variable): %ld\n", numAttrs );
    printf( "\t\t        Max record number for zVariables: %ld\n", maxzRec );



    /*
     * Read global attributes in the CDF. attrScope can be GLOBAL_SCOPE or VARIABLE_SCOPE.
     */
    printf( "\n\t\t                       Global attributes: \n " );
    for ( i=0; i<numAttrs; i++ ) {
        status = CDFinquireAttr (id, i, attrName, &attrScope, &maxgEntry, &maxrEntry, &maxzEntry);
        if ( status < CDF_OK ) StatusHandler( status );

        if ( attrScope == GLOBAL_SCOPE ) {
            if ( maxgEntry < 0 ) {
                printf( "\t\t                                         (empty attribute)\n" );
            } else {
                for ( j=0; j<=maxgEntry; j++ ) {

                    status = CDFinquireAttrgEntry( id, i, (long)j, &datatype, &numElements );
                    if ( status != CDF_OK ) StatusHandler( status );

                    status = CDFgetDataTypeSize( datatype, &ilen );
                    ilen = ilen * numElements;
                    if ( (datatype == CDF_CHAR) || (datatype == CDF_UCHAR) ) ++ilen;
                    buffer = (void *)malloc( ilen );
                    if ( buffer == NULL ) DynamicErrorHandle( id );

                    status = CDFgetAttrgEntry( id, i, (long)j, buffer );
                    if ( status != CDF_OK ) StatusHandler( status );
                    ((char *)buffer)[numElements] = '\0';

                    printf( "\t\t" );
                    if (j==0){
                        strsize = strlen(attrName);
                        for (q=0; q<45-strsize; q++) printf( " " );
                        printf( "%s: %s\n", attrName, (char*)buffer );
                    } else {
                        for (q=0; q<45; q++) printf( " " );
                        printf( "    %s\n", (char*)buffer );
                    }

                    free( buffer );
                }
            }
        }
    }



    /*
     * Read zVariable attributes in the CDF. attrScope can be GLOBAL_SCOPE or VARIABLE_SCOPE.
     */
if (0==1){
    printf( "\n\t\t                     Variable attributes: \n " );
    for ( i=0; i<numAttrs; i++ ) {
        status = CDFinquireAttr (id, i, attrName, &attrScope, &maxgEntry, &maxrEntry, &maxzEntry);
        if ( status < CDF_OK ) StatusHandler( status );

        if ( attrScope == VARIABLE_SCOPE ) {
            if ( maxzEntry < 0 ) {
                printf( "\t\t                                         (empty attribute)\n" );
            } else {
                for ( flag=0, j=0; j<=maxzEntry; j++ ) {


                    numElements = 0;
                    status = CDFinquireAttrzEntry( id, i, (long)j, &datatype, &numElements );
                    if ( status == CDF_OK ) {

                        status = CDFgetzVarName( id, (long)j, varName );
                        if ( status != CDF_OK ) StatusHandler( status );


                        status = CDFgetDataTypeSize( datatype, &ilen );
                        ilen = ilen * numElements;
                        if ( (datatype == CDF_CHAR) || (datatype == CDF_UCHAR) ) ++ilen;
                        if (ilen < 1 ) ilen = 1;
                        buffer = (void *)malloc( ilen );
                        if ( buffer == NULL ) DynamicErrorHandle( id );

                        status = CDFgetAttrzEntry( id, i, (long)j, buffer );
                        if ( status != CDF_OK ) StatusHandler( status );
                        if (numElements>0) {
                            ((char *)buffer)[numElements] = '\0';
                        } else {
                            ((char *)buffer)[0] = '\0';
                        }

                        if (flag==0){
                            printf( "\t\t" );
                            strsize = strlen(attrName);
                            for (q=0; q<45-strsize; q++) printf( " " );
                            printf( "%s:   %s => %s\n", attrName, varName, (char*)buffer );
                            flag = 1;
                        } else {
                            printf( "\t\t" );
                            for (q=0; q<45; q++) printf( " " );
                            printf( "    %s => %s\n", varName, (char*)buffer );
                        }

                        free( buffer );
                    }
                }
            }
        }
    }
}


    /*
     * Read r variable names in the CDF.
     */
    for ( i=0; i<numrVars; i++ ) {

        CDFvarInquire( id, i, varName, &dataType, &numElems, &recVary, &dimVarys );
        if ( status != CDF_OK ) StatusHandler( status );

        printf("r varName %d: %s\n", i, varName );


    }


    /*
     * Read z variable names in the CDF.
     */
    printf( "\n\t\t                There are %ld z variables:\n", numzVars);
    for ( i=0; i<numzVars; i++ ) {

        status = CDFinquirezVar( id, i, varName, &dataType, &numElems, &numDims, dimSizes, &recVary, dimVarys );
        if ( status != CDF_OK ) StatusHandler( status );

        status = CDFgetzVarNumRecsWritten( id, i, &numRecs );
        if (status != CDF_OK) StatusHandler( status );

        if ( verbosity > 1 ) {
            printf( "\n\n\t\t                              zVarName %d: %s\n", i, varName );
            printf( "\t\t                           ----------------------------------------\n");
            DataTypeToStr( dataType, Str );
            printf( "\t\t\t\t                               Data Type: %s\n", Str );
            printf( "\t\t\t\t                      Number of Elements: %ld\n", numElems );
            printf( "\t\t\t\t                    Number of Dimensions: %ld", numDims );
            if ( numDims > 0 ){
                printf("    ( dimSizes: " );
                for (j=0; j<numDims; j++) {
                    if (j!=numDims-1) printf("%ld, ", dimSizes[j] );
                    else              printf("%ld", dimSizes[j] );
                }
                printf(") \n" );
            } else {
                printf("\n" );
            }
            printf( "\t\t\t\t                         Record Variance: %ld\n", recVary );
            printf( "\t\t\t\t                       Number of Records: %ld\n", numRecs );
            printf( "\t\t\t\t                     Variable Attributes: \n" );



            /*
             * Find out which attributes belongm to this variable
             */
            for ( ii=0; ii<numAttrs; ii++ ) {
                status = CDFinquireAttr (id, ii, attrName, &attrScope, &maxgEntry, &maxrEntry, &maxzEntry);
                if ( status < CDF_OK ) StatusHandler( status );

                if ( attrScope == VARIABLE_SCOPE ) {

                    numElements = 0;
                    status = CDFinquireAttrzEntry( id, ii, (long)i, &datatype, &numElements );
                    if ( status == CDF_OK ) {

                        status = CDFgetzVarName( id, (long)i, varName );
                        if ( status != CDF_OK ) StatusHandler( status );


                        status = CDFgetDataTypeSize( datatype, &ilen );
                        ilen = ilen * numElements;
                        if ( (datatype == CDF_CHAR) || (datatype == CDF_UCHAR) ) ++ilen;
                        if (ilen < 1 ) ilen = 1;
                        buffer = (void *)malloc( ilen );
                        if ( buffer == NULL ) DynamicErrorHandle( id );

                        status = CDFgetAttrzEntry( id, ii, (long)i, buffer );
                        if ( status != CDF_OK ) StatusHandler( status );
                        if (numElements>0) {
                            ((char *)buffer)[numElements] = '\0';
                        } else {
                            ((char *)buffer)[0] = '\0';
                        }

                        printf( "\t\t\t\t" );
                        strsize = strlen(attrName);
                        for (q=0; q<45-strsize; q++) printf( " " );
                        printf( "%s: %s\n", attrName, (char*)buffer );

                        free( buffer );
                    }
                }
            }

        } else {
            printf( "\t\t                              zVarName %d: %s\n", i, varName );
        }
    }

}
