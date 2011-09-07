#include "Lgm/Lgm_Misc.h"

void Lgm_ReplaceSubString( char *NewStr, char *OrigStr, char *SubStr, char *RepStr ) {

    int     nOrigStr, nStr, nSubStr, nRepStr, n, nNewStr, done;
    char    *p, *Str;


    // Make a copy of the original string
    nStr = strlen( OrigStr );
    Str  = (char *)calloc( nStr+1, sizeof(char) );
    strcpy( Str, OrigStr );
    strcpy( NewStr, OrigStr );

    done = 0;
    while ( !done ) {

        if ( !(p = strstr( Str, SubStr )) ) {

            done = 1;

        } else {

            // get sizes of strings
            nStr    = strlen( Str );
            nSubStr = strlen( SubStr );
            nRepStr = strlen( RepStr );


            // Copy characters from start of Str start to start of Orig
            strncpy( NewStr, Str, p-Str );
            NewStr[p-Str] = '\0';
            sprintf( NewStr + (p-Str), "%s%s", RepStr, p + nSubStr );
            nNewStr = strlen(NewStr);

            Str = realloc( Str, (nNewStr+1)*sizeof(char) );
            strcpy( Str, NewStr );
        }

    }

    free( Str );

    return;

}



