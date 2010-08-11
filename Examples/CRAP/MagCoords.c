#define MAIN
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <Lgm_CTrans.h>
#include "MagEphemInfo.h"
#include "DynamicMemory.h"





void WriteMagEphemData2( char *Filename, _MagEphemInfo *MagEphemInfo ){

    int fd;

    fd = open( Filename, O_CREAT|O_WRONLY|O_TRUNC, 0644 );
    write( fd, MagEphemInfo, sizeof( *MagEphemInfo ) );
    close( fd );

}





int main( int argc, char *argv[] ){

    char          HostName[128], UserName[128];
    double        UTC, Alpha[1000], a;
    long int      Date;
    int           nAlpha;
    char          Filename2[1024];
    Lgm_Vector    Psm, P;
    Lgm_CTrans    *c = Lgm_init_ctrans(0);
    _MagEphemInfo *MagEphemInfo   = (_MagEphemInfo *)calloc( 1, sizeof( *MagEphemInfo) );

    gethostname( HostName, 126 );
    getlogin_r( UserName, 126 );
    printf( "HostName = %s\n", HostName );
    printf( "UserName = %s\n", UserName );
 
    // Date and UTC
    Date       = 19800625;
    UTC        = 19.0;
    Lgm_Set_Coord_Transforms( Date, UTC, c );

    // Position in SM
    Psm.x = -3.0; Psm.y = 0.0; Psm.z = 0.0;
    Psm.x = -6.6; Psm.y = 0.0; Psm.z = 0.0;
    Psm.x = -1.5; Psm.y = 0.0; Psm.z = 0.0;
    Psm.x = -1.25; Psm.y = 0.0; Psm.z = 0.0;
    Psm.x = -1.05; Psm.y = 0.0; Psm.z = 0.0;
    Lgm_Convert_Coords( &Psm, &P, SM_TO_GSM, c );

    // Create array of Pitch Angles to compute
    for (nAlpha=0,a=1.0; a<=90.0; a+=0.1,++nAlpha) {
        Alpha[nAlpha] = a ;
        printf("Alpha[%d] = %g\n", nAlpha, Alpha[nAlpha]);
    }
//nAlpha = 1;
//Alpha[0] = 2.80;


    /*
     * Compute L*s, Is, Bms, Footprints, etc...
     * These quantities are stored in the MagEphemInfo Structure
     */
    ComputeFieldLineQuantities( Date, UTC, &P, nAlpha, Alpha, MagEphemInfo );

//    // Dump data
//    sprintf( Filename2, "test.dat" );
//    WriteMagEphemData2( Filename2, MagEphemInfo );

    free(c);

    return(0);

}





