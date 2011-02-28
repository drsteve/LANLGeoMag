#include <Lgm_CTrans.h>
int main( ) {
    long int        sec;
    int             s, m, smax, i;
    Lgm_CTrans      *c = Lgm_init_ctrans( 0 );
    Lgm_DateTime    *TAI, *UTC;
    double          ns, u0, u1, t;





    TAI = Lgm_DateTime_Create( 1990, 1, 1, 12.0, LGM_TIME_SYS_TAI, c );
    UTC = Lgm_DateTime_Create( 1990, 1, 1, 12.0, LGM_TIME_SYS_UTC, c );

    printf("          TAI                             UTC\n");
    printf("-------------------------      --------------------------\n");

    for (ns = 2300; ns < 2700; ns += 25 ) {
        TAI->Time  = ns/360000.0;
        TAI->JD = Lgm_JD( TAI->Year, TAI->Month, TAI->Day, TAI->Time, LGM_TIME_SYS_TAI, c );
        Lgm_TAI_to_UTC( TAI, UTC, c );
        Lgm_Print_DateTime( *TAI, 4, 2 ); printf("\t"); Lgm_Print_DateTime( *UTC, 4, 2 ); printf("\n");
    }
    
    

    
    
    // free structures
    Lgm_free_ctrans( c ); 
    free( TAI );
    free( UTC );

    return(0);

}

