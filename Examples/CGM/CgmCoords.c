#include <Lgm_CTrans.h>
#include <Lgm_MagModelInfo.h>

int main( ) {
    Lgm_MagModelInfo *m = Lgm_InitMagInfo();
    Lgm_Vector       Ugeo, Ugsm;
    long int         Date;
    double           UTC, CosLat, L;
    double           glat, glon, gr, galt;
    double           CgmLat, CgmLon, CgmRad, CgmAlt;
    FILE             *fp;

    Lgm_MagModelInfo_Set_MagModel( LGM_IGRF, LGM_EXTMODEL_NULL, m );
    m->Kp = 5;

    Date = 19800101; UTC  = 0.0;
    Lgm_Set_Coord_Transforms( Date, UTC, m->c );

    fp = fopen("input2.txt", "r");

    //printf("MagModel = LGM_IGRF + T89(kP5)\n");
    printf("MagModel = LGM_IGRF \n");
    printf("Re = %g\n", WGS84_A);
    printf("%10s %10s %10s %10s     %10s %10s %10s %10s     %10s\n", "glat", "glon", "galt", "gr", "CgmLat", "CgmLon", "CgmAlt", "CgmRad", "CgmL" );

    while( fscanf( fp, "%lf %lf %lf", &glat, &glon, &galt ) != EOF ) {

        gr   = 1.0 + galt/WGS84_A;
        Ugeo.x = gr*cos(RadPerDeg*glat)*cos(RadPerDeg*glon);
        Ugeo.y = gr*cos(RadPerDeg*glat)*sin(RadPerDeg*glon);
        Ugeo.z = gr*sin(RadPerDeg*glat);

        Lgm_Convert_Coords( &Ugeo, &Ugsm, GEO_TO_GSM, m->c );
        Lgm_CDMAG_TO_CGM( &Ugsm, &CgmLat, &CgmLon, &CgmRad, m );

        CgmAlt  = (CgmRad - 1.0)*WGS84_A;
        CosLat   = cos( CgmLat*RadPerDeg );

        L    = CgmRad/(CosLat*CosLat);
        printf("%10g %10g %10g %10g     %10g %10g %10g %10g     %10g\n", glat, glon, galt, gr, CgmLat, CgmLon, CgmAlt, CgmRad, L );

    }

    Lgm_FreeMagInfo( m );
    fclose(fp);

    return(0);
}

