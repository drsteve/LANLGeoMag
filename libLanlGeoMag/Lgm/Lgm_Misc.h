#ifndef LGM_MISC_H
#define LGM_MISC_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


void Lgm_ReplaceSubString2( char **OutStr, char *OrigStr, char *SubStr, char *RepStr );
void Lgm_ReplaceSubString( char *OutStr, char *OrigStr, char *SubStr, char *RepStr );
char *Lgm_StrToLower( char *str, int nmax );
char *Lgm_StrToUpper( char *str, int nmax );

#endif
