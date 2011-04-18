#include <stdio.h>
main(){

    char    Filename[1024];
    char    Command[1024];
    int     n, j, N;
    FILE    *fp1;


    fp1 = fopen("List.txt", "r");
    system( "mkdir -p tmp" );
    system( "rm tmp/*.png" );

    n = 0;
    N = 5;
    while ( fscanf( fp1, "%s", Filename ) != EOF ){
        for (j=0; j<N; j++){
            sprintf(Command, "cp %s tmp/%03d.png", Filename, n);
            system(Command);
            ++n;
        }
    }

    //sprintf( Command, "ffmpeg -b 16192k -f image2 -i ./%%03d.png Movie.mp4" );
    sprintf( Command, "ffmpeg -b 32384k -f image2 -i tmp/%%03d.png Movie.mp4" );
    printf( "Command = %s\n", Command);
    system( Command );


    

}
