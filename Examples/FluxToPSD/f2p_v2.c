#include <stdio.h>
#include <math.h>
#include <Lgm_MagModelInfo.h>
#include <Lgm_FluxToPsd.h>
void DumpSvg( int,  char * );
int main( ) {

    Lgm_CTrans     *c = Lgm_init_ctrans( 0 );
    Lgm_DateTime    d;
    Lgm_Vector      u, v;
    //Lgm_FluxToPsd  *f2p = Lgm_F2P_CreateFluxToPsd(1);
    Lgm_FluxToPsd  *f2p;
    Lgm_PsdToFlux  *p2f;

    double          ***J, *E, *A, **J_DIFF;
    int             nE, nA, i, j, n, k;

    double          *Mu, *K;
    int             nMu, nK;

    double          f, p2c2, df, sa;

    double          S[100], a, b, r;

    char            Line[5000], ISO_TIME[200][80], Command[512];

    long int        pos;
    FILE            *fp;




    /*
     *  Make up some fake flux data... Lets take a two-component relativistic
     *  maxwellian as described by Cayton et al. [1989].
     *
     *  Maxwellian 1; n= 1e-4 cm^-3; T = 200 keV
     *  Maxwellian 2; n= 5e-3 cm^-3; T =  25 keV
     *
     */
    nE = 10;
    nA = 18;
    LGM_ARRAY_1D( E, nE, double );
    LGM_ARRAY_1D( A, nA, double );
    LGM_ARRAY_3D( J, 200, nE, nA, double );

    /*
     *  Choose a set of Mu's and K's 
     *  and compute Phase Space Density at these constant Mu's and K's
     */
    nMu = 36;
    nK  = 36;
    LGM_ARRAY_1D( Mu, nMu, double );
    LGM_ARRAY_1D( K, nK, double );
    for (i=0; i<nMu; i++) Mu[i] = 1.0 + 60.0*i;


    a = 0.01; b = 10.0;
    r = pow( b/a, 1.0/((double)(nK-1)));
    S[0] = a;
    for (j=1; j<nK; j++) S[j] = S[j-1]*r;
    //for (j=0; j<nK; j++) K[j] = S[nK-1-j];
    for (j=0; j<nK; j++) K[j] = S[j];

a = 1.0; b = 2000.0;
r = pow( b/a, 1.0/((double)(nMu-1)));
S[0] = a;
for (j=1; j<nMu; j++) S[j] = S[j-1]*r;
for (j=0; j<nMu; j++) Mu[j] = S[j];




    /*
     * Set date/time and position
     */
    //d = Lgm_DateTime_Create( 2001, 1, 1, 0.0, LGM_TIME_SYS_UTC, c );
    u.x = 6.6*cos(133.0*RadPerDeg); 
    u.y = 6.6*sin(133.0*RadPerDeg);
    u.z = 0.0;



    /*
     * Read in a sample GEO file
     */
    n = 0;
    fp = fopen("/data2/DREAM/data/LANL-97A/2002/20021211_LANL-97A_l3_data.txt", "r");
    //fp = fopen("puke.txt", "r");
    pos = ftell( fp );
    n = 0;
    while ( fgets( Line, 4096, fp ) != NULL ) {
        if ( Line[0] == '\n' ) {
            // do nothing
        } else if ( Line[0] != '#' ) {
            // backup so the line is available to read again.
            fseek( fp, pos, SEEK_SET );

            // read in data
            for (j=0; j<18; j++) {
                fscanf( fp, "%s", ISO_TIME[n] );
                if (j==0) printf("Read data for Time: %s\n", ISO_TIME[n] );
                fscanf( fp, "%lf", &A[j] );
                for (i=0; i<10; i++) {
                    fscanf( fp, "%lf", &J[n][i][j] );
                    //printf( "%lf ", J[n][i][j] );
                }
                //printf("\n");
            }
            ++n;
        } else if ( strstr( Line, "Energy Channels:" ) != NULL ) {
            // extract energy channel info from header.
            sscanf( Line, "%*[^:]:%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &E[0], &E[1], &E[2], &E[3], &E[4], &E[5], &E[6], &E[7], &E[8], &E[9] );
        }
        pos = ftell( fp );
    }
    fclose( fp );



    LGM_ARRAY_2D( J_DIFF, nE, nA, double );


    for (k=0; k<n; k++){

            /*
             * Go Forward
             */
            f2p = Lgm_F2P_CreateFluxToPsd(1);
            f2p->Extrapolate = 0;

            IsoTimeStringToDateTime( ISO_TIME[k], &d, c );
            printf("Date = %ld   UTC = %g\n", d.Date, d.Time );
            Lgm_Set_Coord_Transforms( d.Date, d.Time, c );

            Lgm_Convert_Coords( &u, &v, WGS84_TO_GSM, c );
            Lgm_F2P_SetDateTimeAndPos( &d, &v, f2p );

            // Add diff. flux data/info to f2p structure.
            Lgm_F2P_SetFlux( J[k], E, nE, A, nA, f2p ); 

printf("k = %d\n", k);
            //  Compute Phase Space Density at the constant Mu's and K's
            Lgm_F2P_GetPsdAtConstMusAndKs( Mu, nMu, K, nK, f2p );


            sprintf( Command, "convert -sample 400x400\! Lgm_FluxToPsd_FLUX_EA.gif Lgm_FluxToPsd_FLUX_EA_%03d.gif", k); system( Command );
            sprintf( Command, "convert -sample 400x400\! Lgm_FluxToPsd_PSD_EA.gif Lgm_FluxToPsd_PSD_EA_%03d.gif", k); system( Command );
            sprintf( Command, "convert -sample 400x400\! Lgm_FluxToPsd_PSD_MK.gif Lgm_FluxToPsd_PSD_MK_%03d.gif", k); system( Command );

            sprintf( Command, "cp Lgm_FluxToPsd_FLUX_EA.info Lgm_FluxToPsd_FLUX_EA_%03d.info", k); system( Command );
            sprintf( Command, "cp Lgm_FluxToPsd_PSD_EA.info Lgm_FluxToPsd_PSD_EA_%03d.info", k); system( Command );
            sprintf( Command, "cp Lgm_FluxToPsd_PSD_MK.info Lgm_FluxToPsd_PSD_MK_%03d.info", k); system( Command );

            DumpSvg( k, ISO_TIME[k] );
            sprintf( Command, "inkscape -e %03d.png mike.svg", k ); system( Command );

            




            /*
             * Go Backward
             */
            p2f = Lgm_P2F_CreatePsdToFlux(1);
            p2f->Extrapolate = 0;
            Lgm_P2F_SetDateTimeAndPos( &d, &v, p2f );
            // Add PSD data/info to p2f structure.
            Lgm_P2F_SetPsd( f2p->PSD_MK, f2p->Mu, f2p->nMu, f2p->K, f2p->nK, p2f ); 
            //  Compute Phase Space Density at the constant Mu's and K's
            Lgm_P2F_GetFluxAtConstEsAndAs( E, nE, A, nA, p2f );

            sprintf( Command, "convert -sample 400x400\! Lgm_PsdToFlux_FLUX_EA.gif Lgm_PsdToFlux_FLUX_EA_%03d.gif", k); system( Command );
            sprintf( Command, "convert -sample 400x400\! Lgm_PsdToFlux_PSD_EA.gif Lgm_PsdToFlux_PSD_EA_%03d.gif", k); system( Command );
            sprintf( Command, "convert -sample 400x400\! Lgm_PsdToFlux_PSD_MK.gif Lgm_PsdToFlux_PSD_MK_%03d.gif", k); system( Command );

            sprintf( Command, "cp Lgm_PsdToFlux_FLUX_EA.info Lgm_PsdToFlux_FLUX_EA_%03d.info", k); system( Command );
            sprintf( Command, "cp Lgm_PsdToFlux_PSD_EA.info Lgm_PsdToFlux_PSD_EA_%03d.info", k); system( Command );
            sprintf( Command, "cp Lgm_PsdToFlux_PSD_MK.info Lgm_PsdToFlux_PSD_MK_%03d.info", k); system( Command );



            /*
             * Difference
             */
            for ( i=0; i<nE; i++ ){
                for ( j=0; j<nA; j++ ){
                    if ( ( p2f->FLUX_EA[i][j] > 0.0 ) && ( J[k][i][j] > 0.0 ) ){
                        J_DIFF[i][j] = 100.0*fabs( (log10(J[k][i][j]) - log10(p2f->FLUX_EA[i][j]))/log10(p2f->FLUX_EA[i][j]) );
                    } else {
                        J_DIFF[i][j] = 0.0;
                    }
                }
            }
            DumpGif( "J_DIFF", nA, nE, J_DIFF );
            sprintf( Command, "convert -sample 400x400\! J_DIFF.gif J_DIFF_%03d.gif", k); system( Command );
            sprintf( Command, "cp J_DIFF.info J_DIFF_%03d.info", k); system( Command );



            Lgm_P2F_FreePsdToFlux( p2f );
            Lgm_F2P_FreeFluxToPsd( f2p );


    }
        

    LGM_ARRAY_1D_FREE( E );
    LGM_ARRAY_1D_FREE( A );
    LGM_ARRAY_3D_FREE( J );
    LGM_ARRAY_2D_FREE( J_DIFF );
    LGM_ARRAY_1D_FREE( Mu );
    LGM_ARRAY_1D_FREE( K );


    Lgm_free_ctrans( c );


    return(0);
}

