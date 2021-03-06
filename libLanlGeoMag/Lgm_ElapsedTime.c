/*! \file Lgm_ElapsedTime.c
 *
 *  \brief  Utility routines to compute and print run-times and elapsed times. (Useful for timing long runs.)
 *
 */
#include <strings.h>
#include "Lgm/Lgm_ElapsedTime.h"


void Lgm_ElapsedTimeInit( Lgm_ElapsedTimeInfo *t, int red, int grn, int blu ) {

    int ri, gi, bi, index;

    t->RunStartTime = time( NULL );

    ri = 5*red/255;
    gi = 5*grn/255;
    bi = 5*blu/255;

    index = ri*36 + gi*6 + bi +16;
    t->ColorizeText = ( index != 0 ) ? 1 : 0;

    /*
     *  If we are colorizing, set the "Color code strings" appropriately
     */
    if ( t->ColorizeText ) {

        sprintf( t->ColorStart, "[38;5;%dm", index );
        sprintf( t->ColorEnd, "[0m" );

    } else {

        strcpy( t->ColorStart, "" );
        strcpy( t->ColorEnd, "" );

    }

    

}



double Lgm_SetElapsedTimeStr( Lgm_ElapsedTimeInfo *t ) {

    int     Days, Hours, Minutes;
    double  Seconds, s;
    time_t  CurrentTime;

    CurrentTime = time( NULL );
    s           = difftime( CurrentTime, t->RunStartTime );
    Seconds     = s;

    Days     = (int)(Seconds/86400.0);
    Seconds -= Days*86400.0;

    Hours    = (int)(Seconds/3600.0);
    Seconds -= Hours*3600.0;

    Minutes  = (int)(Seconds/60.0);
    Seconds -= Minutes*60.0;

    sprintf( t->ElapsedTimeStr, "%03d:%02d:%02d:%02d", Days, Hours, Minutes, (int)(Seconds+0.5)  );

    return( s );

}
void Lgm_PrintElapsedTime( Lgm_ElapsedTimeInfo *t ) {

    Lgm_SetElapsedTimeStr( t );
    printf( "%s*****  Elapsed Time (DD:HH:MM:SS): %s  *****%s\n", t->ColorStart, t->ElapsedTimeStr, t->ColorEnd  );

}



void Lgm_SetCurrentTimeStr( Lgm_ElapsedTimeInfo *t ) {

    struct tm   tp;
    time_t      Time;
    char        Buf[40], Str[80], *p;

    Time = time( NULL );
    localtime_r( &Time, &tp );
    sprintf( Str, "%s", asctime_r( &tp, Buf ) );
    p = index( Str, '\n' ); *p = '\0';
    sprintf( t->CurrentTimeStr, "%s", Str );
    sprintf( t->CurrentTimeStr2, "%4d-%02d-%02dT%02d:%02d:%02d", tp.tm_year+1900, tp.tm_mon+1, tp.tm_mday, tp.tm_hour, tp.tm_min, tp.tm_sec );

}
void Lgm_PrintCurrentTime( Lgm_ElapsedTimeInfo *t ) {

    Lgm_SetCurrentTimeStr( t );
    printf( "%s*****  Current Time: %s  *****%s\n", t->ColorStart, t->CurrentTimeStr2, t->ColorEnd  );

}
