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


double Lgm_PrintElapsedTime( Lgm_ElapsedTimeInfo *t ) {

    int     Days, Hours, Minutes;
    double  Seconds;
    time_t  CurrentTime;

    CurrentTime = time( NULL );
    Seconds     = difftime( CurrentTime, t->RunStartTime );

    Days     = (int)(Seconds/86400.0);
    Seconds -= Days*86400.0;

    Hours    = (int)(Seconds/3600.0);
    Seconds -= Hours*3600.0;

    Minutes  = (int)(Seconds/60.0);
    Seconds -= Minutes*60.0;

    printf( "%s*****  Elapsed Time (DD:HH:MM:SS): %03d:%02d:%02d:%02d  *****%s\n", t->ColorStart, Days, Hours, Minutes, (int)(Seconds+0.5), t->ColorEnd  );

    return( Seconds );

}


void Lgm_PrintCurrentTime( Lgm_ElapsedTimeInfo *t ) {

    struct tm   *tp;
    time_t      Time;
    char        Line[80], *p;

    Time = time( NULL );
    tp   = localtime( &Time );
    sprintf( Line, "%s", asctime( tp ) );
    p = index( Line, '\n' ); *p = '\0';
    printf( "%s*****  Current Time: %s  *****%s\n", t->ColorStart, Line, t->ColorEnd  );


}
