#include "Lgm/Lgm_Misc.h"

/*
 *  Routine to replace a substring with another string (of potentially
 *  different length) User provides OrigStr, SubStr, RepStr and these dont get
 *  modified in any way.
 *
 *  NewStr is allocated and returned. 
 *  User is responsible for freeing the memory.
 *
 */
void Lgm_ReplaceSubString2( char **OutStr, char *OrigStr, char *SubStr, char *RepStr ) {

    int     nOrigStr, nStr, nSubStr, nRepStr, n, nNewStr, nCopy1, nCopy2, done;
    char    *p, *Str, *Str2, *p_remaining, *NewStr;
    int     q;


    // get sizes of strings 
    nSubStr = strlen( SubStr );
    nRepStr = strlen( RepStr );



    // Make a copy of the original string
    nStr = strlen( OrigStr );
    Str  = (char *)calloc( nStr+1, sizeof(char) );
    strcpy( Str, OrigStr );
    printf("Str = %s n = %d\n", Str, nStr);


    // Allocate space for the NewStr -- initially we only need enough space for the null terminator
    q = 1;
    NewStr  = (char *)calloc( q, sizeof(char) );
    NewStr[q-1] = '\0'; // terminate it.
    nNewStr = strlen( NewStr );



    /*
     * Set pointer to look at the part of Str the we have left to
     * examine/replace initially it is just the start
     */
    p_remaining = Str;
    done = 0;
    while ( !done ) {

        /*
         *  Using the string pointed to by p_remaining, attempt to locate an
         *  occurrence of SubStr
         */
        if ( !(p = strstr( p_remaining, SubStr )) ) {

            /*
             *   We didn't find an occurrence of SubStr in the String pointed to
             *   by p_remaining. So we are done.
             *   But we do need to append p_remaining onto NewStr....
             */
            nCopy1 = strlen( p_remaining );
            q += nCopy1;
            NewStr = (char *)realloc( NewStr, q*sizeof(char) );

            // Append characters from start of p_remaining NewStr
            strncpy( &NewStr[nNewStr], p_remaining, nCopy1 ); // this copies over everything up to and including the char before where p is pointing.
            NewStr[nNewStr+nCopy1] = '\0'; // terminate it.
            nNewStr = strlen( NewStr );

            done = 1;

        } else {

            /*
             * We have found an occurrence of the substring. The pointer p
             * is pointing at it in the string p_remaining
             *
             * We need to:
             *
             *      1) copy into the NewStr, everything in the p_remaining
             *         from its start up until the char just before this point
             *
             *      2) copy into the NewStr, the new RepStr as a replacement
             *
             *      3) find the place in p_remaining that follows the SubStr
             *         we just replaced.
             *
             */
            //determine how many chars we need to copy over
            nCopy1 = p-p_remaining;
            nCopy2 = nRepStr;

            q += nCopy1;
            q += nCopy2;

            // realloc the size of the NewStr array to accomodate the part we need to copy + the size of the RepStr that will get copied
            NewStr = (char *)realloc( NewStr, q*sizeof(char) );

            // Append characters from start of p_remaining NewStr
            strncpy( &NewStr[nNewStr], p_remaining, nCopy1 ); // this copies over everything up to and including the char before where p is pointing.
            NewStr[nNewStr+nCopy1] = '\0'; // terminate it.
            nNewStr = strlen( NewStr );

            // Append characters from RepStr
            strncpy( &NewStr[nNewStr], RepStr, nCopy2 ); // this copies over contents of the Replacement string
            NewStr[nNewStr+nCopy2] = '\0'; // terminate it.
            nNewStr = strlen(NewStr);

            p_remaining = p + nSubStr;

        }

    }

    free( Str );


    *OutStr = NewStr;

    return;

}



/*
 *  Legacy version... For backwards compatibility.
 *  Routine to replace a substring with another string (of potentially
 *  different length) User provides OrigStr, SubStr, RepStr and these dont get
 *  modified in any way.
 *
 *  In this version, the new string is allocated by the user outside the routine.
 *  It must be big enough to hold the replaced substrings, so is not as safe as the other version.
 */
void Lgm_ReplaceSubString( char *OutStr, char *OrigStr, char *SubStr, char *RepStr ) {

    int     nOrigStr, nStr, nSubStr, nRepStr, n, nNewStr, nCopy1, nCopy2, done;
    char    *p, *Str, *Str2, *p_remaining, *NewStr;
    int     q;


    // get sizes of strings 
    nSubStr = strlen( SubStr );
    nRepStr = strlen( RepStr );



    // Make a copy of the original string
    nStr = strlen( OrigStr );
    Str  = (char *)calloc( nStr+1, sizeof(char) );
    strcpy( Str, OrigStr );
    printf("Str = %s n = %d\n", Str, nStr);


    // Allocate space for the NewStr -- initially we only need enough space for the null terminator
    q = 1;
    NewStr  = (char *)calloc( q, sizeof(char) );
    NewStr[q-1] = '\0'; // terminate it.
    nNewStr = strlen( NewStr );



    /*
     * Set pointer to look at the part of Str the we have left to
     * examine/replace initially it is just the start
     */
    p_remaining = Str;
    done = 0;
    while ( !done ) {

        /*
         *  Using the string pointed to by p_remaining, attempt to locate an
         *  occurrence of SubStr
         */
        if ( !(p = strstr( p_remaining, SubStr )) ) {

            /*
             *   We didn't find an occurrence of SubStr in the String pointed to
             *   by p_remaining. So we are done.
             *   But we do need to append p_remaining onto NewStr....
             */
            nCopy1 = strlen( p_remaining );
            q += nCopy1;
            NewStr = (char *)realloc( NewStr, q*sizeof(char) );

            // Append characters from start of p_remaining NewStr
            strncpy( &NewStr[nNewStr], p_remaining, nCopy1 ); // this copies over everything up to and including the char before where p is pointing.
            NewStr[nNewStr+nCopy1] = '\0'; // terminate it.
            nNewStr = strlen( NewStr );

            done = 1;

        } else {

            /*
             * We have found an occurrence of the substring. The pointer p
             * is pointing at it in the string p_remaining
             *
             * We need to:
             *
             *      1) copy into the NewStr, everything in the p_remaining
             *         from its start up until the char just before this point
             *
             *      2) copy into the NewStr, the new RepStr as a replacement
             *
             *      3) find the place in p_remaining that follows the SubStr
             *         we just replaced.
             *
             */
            //determine how many chars we need to copy over
            nCopy1 = p-p_remaining;
            nCopy2 = nRepStr;

            q += nCopy1;
            q += nCopy2;

            // realloc the size of the NewStr array to accomodate the part we need to copy + the size of the RepStr that will get copied
            NewStr = (char *)realloc( NewStr, q*sizeof(char) );

            // Append characters from start of p_remaining NewStr
            strncpy( &NewStr[nNewStr], p_remaining, nCopy1 ); // this copies over everything up to and including the char before where p is pointing.
            NewStr[nNewStr+nCopy1] = '\0'; // terminate it.
            nNewStr = strlen( NewStr );

            // Append characters from RepStr
            strncpy( &NewStr[nNewStr], RepStr, nCopy2 ); // this copies over contents of the Replacement string
            NewStr[nNewStr+nCopy2] = '\0'; // terminate it.
            nNewStr = strlen(NewStr);

            p_remaining = p + nSubStr;

        }

    }

    free( Str );


    strcpy( OutStr, NewStr );
    free( NewStr);

    return;

}



