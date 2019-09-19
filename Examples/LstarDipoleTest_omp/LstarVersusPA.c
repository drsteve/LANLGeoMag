#define MAIN
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <Lgm_CTrans.h>
#include <Lgm_MagEphemInfo.h>
#include <Lgm_QinDenton.h>

#define KP_DEFAULT 0

int main( int argc, char *argv[] ){

    double           L, I, Bm, M;
    double           UTC, Alpha[1000], a, dist, JD;
    long int         Date;
    int              nAlpha, Kp, Colorize, i, j;
    char             Filename[1024];
    Lgm_Vector       Psm, P;
    Lgm_MagEphemInfo *MagEphemInfo;
    Lgm_QinDentonOne p;
    FILE             *fpout;


    // Create array of Pitch Angles to compute
    for (nAlpha=0,a=1.0; a<=90.0; a+=5.0, ++nAlpha) {
        Alpha[nAlpha] = a ;
    }
// Override with wahtever you want here.
//nAlpha = 1;
//Alpha[0] = 83.0;

    if ( nAlpha > 0 ){
        MagEphemInfo = Lgm_InitMagEphemInfo(0, nAlpha);
    } else {
        // doesnt seem to like allocating zero size...
        MagEphemInfo = Lgm_InitMagEphemInfo(0, 1);
    }



    // Date and UTC -- pick a time and date when tilt angle, psi ~ 0
    Date       = 20180524;
    UTC        = 23.562;
    JD = Lgm_Date_to_JD( Date, UTC, MagEphemInfo->LstarInfo->mInfo->c );
    Lgm_Set_Coord_Transforms( Date, UTC, MagEphemInfo->LstarInfo->mInfo->c );
    printf("Geo-Dipole Tile Angle: = %g Degrees\n", MagEphemInfo->LstarInfo->mInfo->c->psi*DegPerRad);

    // Position in SM
    Psm.x = -6.00; Psm.y = 0.0; Psm.z = 0.0;
    Lgm_Convert_Coords( &Psm, &P, SM_TO_GSM, MagEphemInfo->LstarInfo->mInfo->c );






    //USER INPUT STUFF
    Lgm_SetMagEphemLstarQuality( 3, 96, MagEphemInfo );
    MagEphemInfo->SaveShellLines = TRUE;
    MagEphemInfo->LstarInfo->VerbosityLevel = 2;
    MagEphemInfo->LstarInfo->mInfo->VerbosityLevel = 2;
    MagEphemInfo->LstarInfo->LstarMoment = LGM_LSTAR_MOMENT_CDIP;

    Kp = 1;
    MagEphemInfo->LstarInfo->mInfo->Kp = ( Kp >= 0 ) ? Kp : KP_DEFAULT;
    if ( MagEphemInfo->LstarInfo->mInfo->Kp > 5 ) MagEphemInfo->LstarInfo->mInfo->Kp = 5;
    MagEphemInfo->LstarInfo->mInfo->Kp = 0;
    
    // If you want T89 instead of dipole
    //Lgm_Set_MagModel( LGM_CDIP, LGM_EXTMODEL_T89, MagEphemInfo->LstarInfo->mInfo );

    // Test CDIP
    Lgm_Set_MagModel( LGM_CDIP, LGM_EXTMODEL_NULL, MagEphemInfo->LstarInfo->mInfo );


    Lgm_Set_LossConeHeight( MagEphemInfo->LstarInfo->mInfo, 3.0 );
    MagEphemInfo->LstarInfo->ISearchMethod = 2;
    Lgm_get_QinDenton_at_JD( JD, &p, 1, 1 );
    Lgm_set_QinDenton( &p, MagEphemInfo->LstarInfo->mInfo );

    /*
     * Compute L*s, Is, Bms, Footprints, etc...
     * These quantities are stored in the MagEphemInfo Structure
     */
    Colorize = TRUE;
    MagEphemInfo->LstarInfo->LSimpleMax = 15.0; // Extends threshold for doing the calucation
    MagEphemInfo->LstarInfo->VerbosityLevel = 1;
    MagEphemInfo->LstarInfo->mInfo->VerbosityLevel = 0;
    //MagEphemInfo->LstarInfo->mInfo->Lgm_MagStep_BS_atol = 1e-5;
    //MagEphemInfo->LstarInfo->mInfo->Lgm_MagStep_BS_rtol = 0.0;
    //MagEphemInfo->LstarInfo->ShabanskyHandling = LGM_SHABANSKY_IGNORE;
    MagEphemInfo->LstarInfo->ShabanskyHandling = LGM_SHABANSKY_HALVE_I;
    Lgm_ComputeLstarVersusPA( Date, UTC, &P, nAlpha, Alpha, Colorize, MagEphemInfo );


    /*
     * Dump results
     */
    sprintf( Filename, "DipoleTestResults_%d_Meth%d.dat", MagEphemInfo->LstarQuality, MagEphemInfo->LstarInfo->ISearchMethod);
    dist = Lgm_Magnitude(&Psm);
    fpout = fopen(Filename, "w");
    for (i=0; i<nAlpha; i++ ){
//        if ( MagEphemInfo->Lstar[i] > 0.0 ) {
            fprintf(fpout, "%.9lf %.9lf %.9lf %g\n", MagEphemInfo->Alpha[i], dist, MagEphemInfo->Lstar[i], dist-MagEphemInfo->Lstar[i] );
            printf("%.9lf %.9lf %.9lf %g %g\n", MagEphemInfo->Alpha[i], dist, MagEphemInfo->Lstar[i], dist-MagEphemInfo->Lstar[i], Re*(dist-MagEphemInfo->Lstar[i]) );
//        }
    }
    fclose(fpout);


    /*
     * Shows how to access the actual drift shells
     */
    double aa, B0_eq, Beq, Beq2, Sin2_AlphaEq, sa, sa2;
    int    k;
    sprintf( Filename, "DriftShells2.dat" );
    fpout = fopen(Filename, "w");
    for (i=0; i<nAlpha; i++ ){
        if ( MagEphemInfo->Lstar[i] > 0.0 ) {
            for (j=0; j<MagEphemInfo->nShellPoints[i]; j++ ){ //Loop over shell field lines

                // Find the (approximate) Bmin point on the FL segments that
                // define the shell -- note that this could be different from
                // the global Bmin if its a Shabansky orbit.
                Beq2 = 1e99;
                for (k=0; k<MagEphemInfo->nFieldPnts[i][j]; k++ ) {
                    if ( MagEphemInfo->Bmag[i][j][k] < Beq2 ) Beq2 = MagEphemInfo->Bmag[i][j][k];
                }

                //Keep track of what the equatorial pitch angle is at each longitud at each longitude.
                aa = MagEphemInfo->Alpha[i]; sa = sin( aa*RadPerDeg ); sa2 = sa*sa; // initial PA.
                B0_eq = MagEphemInfo->Bmin; // Initial Bmin
                Beq   = Lgm_Magnitude( &MagEphemInfo->Shell_Bmin[i][j] );
                Sin2_AlphaEq =  Beq2/B0_eq * sa2;

                fprintf(fpout, "%d %g %g %g %.9lf %.9lf %.9lf %.9lf\n", j, Beq, Beq2, MagEphemInfo->Alpha[i], DegPerRad*asin( sqrt( Sin2_AlphaEq ) ), MagEphemInfo->ShellSphericalFootprint_Pn[i][j].x, 
                                                                           MagEphemInfo->ShellSphericalFootprint_Pn[i][j].y, 
                                                                           MagEphemInfo->ShellSphericalFootprint_Pn[i][j].z );

            }
        }
    }
    fclose(fpout);
    


    /*
     * Dump the MagEphemInfo structure in order to visualize it with VieDriftShell App
     */
    WriteMagEphemInfoStruct( "LstarTest.dat", nAlpha, MagEphemInfo );

    Lgm_FreeMagEphemInfo( MagEphemInfo );

    return(0);

}





