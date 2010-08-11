#include "MagEphemInfo.h"

void WriteMagEphemData( char *Filename, long int Date, double UT, double Lat, double Lon, double Rad, int Mode ){

    int     i, j;
    char    Str[80];
    FILE    *fp;


    if ( Mode ) {

        /*
         * Write Header
         */
        fp = fopen(Filename, "wb");
        fprintf(fp, "nAlpha:  %d; ", MagEphemInfo->nAlpha );  for (i=0; i<MagEphemInfo->nAlpha; i++)  fprintf(fp, " %g", MagEphemInfo->Alpha[i]); fprintf(fp, "\n" );

        sprintf( Str, "Date" ); fprintf(fp, "%-10.10s", Str); 
        sprintf( Str, "UT" );   fprintf(fp, " %13.13s", Str);
        sprintf( Str, "Lat" );  fprintf(fp, " %8.8s", Str);
        sprintf( Str, "Lon" );  fprintf(fp, " %8.8s", Str);
        sprintf( Str, "Rad" );  fprintf(fp, " %8.8s", Str);
        fprintf(fp, "    ");
        for (i=0; i<MagEphemInfo->nAlpha; i++) { sprintf( Str, "L*%d", i); fprintf(fp, " %8.8s", Str); } fprintf(fp, "    ");
        for (i=0; i<MagEphemInfo->nAlpha; i++) { sprintf( Str, "Bm*%d", i); fprintf(fp, " %10.10s", Str); } fprintf(fp, "    ");
        for (i=0; i<MagEphemInfo->nAlpha; i++) { sprintf( Str, "I%d", i); fprintf(fp, " %8.8s", Str); } fprintf(fp, "    ");
        for (i=0; i<MagEphemInfo->nAlpha; i++) { sprintf( Str, "K*%d", i); fprintf(fp, " %8.8s", Str); } fprintf(fp, "    ");
        for (i=0; i<MagEphemInfo->nAlpha; i++) { sprintf( Str, "Sb%d", i); fprintf(fp, " %8.8s", Str); } fprintf(fp, "    ");
        for (i=0; i<MagEphemInfo->nAlpha; i++) { sprintf( Str, "Tb%d", i); fprintf(fp, " %8.8s", Str); } fprintf(fp, "    ");
        for (i=0; i<MagEphemInfo->nAlpha; i++) { sprintf( Str, "Lm%d", i); fprintf(fp, " %8.8s", Str); } fprintf(fp, "    ");
        fprintf(fp, "\n");
        fclose( fp );

    } 
    
    fp = fopen(Filename, "ab");
    fprintf(fp, "%-10ld", Date);
    fprintf(fp, " %13.8lf",  UT);
    fprintf(fp, " %8g",  Lat);
    fprintf(fp, " %8g",  Lon);
    fprintf(fp, " %8g",  Rad);
    fprintf(fp, "    ");
    for (i=0; i<MagEphemInfo->nAlpha; i++) { fprintf(fp, " %8g", MagEphemInfo->Lstar[i]); } fprintf(fp, "    ");
    for (i=0; i<MagEphemInfo->nAlpha; i++) { fprintf(fp, " %10.8g", MagEphemInfo->Bm[i]); } fprintf(fp, "    ");
    for (i=0; i<MagEphemInfo->nAlpha; i++) { fprintf(fp, " %8g", MagEphemInfo->I[i]); } fprintf(fp, "    ");
    for (i=0; i<MagEphemInfo->nAlpha; i++) { fprintf(fp, " %8g", MagEphemInfo->K[i]); } fprintf(fp, "    ");
    for (i=0; i<MagEphemInfo->nAlpha; i++) { fprintf(fp, " %8g", MagEphemInfo->Sb[i]); } fprintf(fp, "    ");
    for (i=0; i<MagEphemInfo->nAlpha; i++) { fprintf(fp, " %8g", MagEphemInfo->Tb[i]); } fprintf(fp, "    ");
    for (i=0; i<MagEphemInfo->nAlpha; i++) { fprintf(fp, " %8g", MagEphemInfo->LMcIlwain[i]); } 

    for (i=0; i<MagEphemInfo->nAlpha; i++) { fprintf(fp, " %8g %8g %8g", MagEphemInfo->Pmn_gsm[i].x, MagEphemInfo->Pmn_gsm[i].y, MagEphemInfo->Pmn_gsm[i].z); } fprintf(fp, "    ");
    for (i=0; i<MagEphemInfo->nAlpha; i++) { fprintf(fp, " %8g %8g %8g", MagEphemInfo->Pms_gsm[i].x, MagEphemInfo->Pms_gsm[i].y, MagEphemInfo->Pms_gsm[i].z); } fprintf(fp, "    ");

    for (i=0; i<MagEphemInfo->nAlpha; i++) { 
        fprintf(fp, " %8d", MagEphemInfo->nShellPoints[i]);
        for (j=0; j<MagEphemInfo->nShellPoints[i]; j++) { 
            fprintf(fp, " %8g %8g %8g", MagEphemInfo->ShellFootprint_Pn[i][j].x, MagEphemInfo->ShellFootprint_Pn[i][j].y, MagEphemInfo->ShellFootprint_Pn[i][j].z );
            fprintf(fp, " %8g %8g %8g", MagEphemInfo->ShellFootprint_Ps[i][j].x, MagEphemInfo->ShellFootprint_Ps[i][j].y, MagEphemInfo->ShellFootprint_Ps[i][j].z );
            fprintf(fp, " %8g %8g %8g", MagEphemInfo->ShellMirror_Pn[i][j].x, MagEphemInfo->ShellMirror_Pn[i][j].y, MagEphemInfo->ShellMirror_Pn[i][j].z );
            fprintf(fp, " %8g %8g %8g", MagEphemInfo->ShellMirror_Ps[i][j].x, MagEphemInfo->ShellMirror_Ps[i][j].y, MagEphemInfo->ShellMirror_Ps[i][j].z );
        }
    }

    fprintf(fp, "\n");
    fclose( fp );

}


void ReadMagEphemData( char *Filename, long int Date, double UT, double Lat, double Lon, double Rad, int Mode ){

    int     i, j;
    char    Str[80], *Header;
    FILE    *fp;


    /*
     *  Read off the header
     */
    Header = (char *)calloc( 30000, sizeof(char));
    fp = fopen(Filename, "rb");
    fgets( Header, 30000, fp );
    fgets( Header, 30000, fp );
    free( Header );

    fscanf(fp, "%ld", &MagEphemInfo->Date);
    fscanf(fp, " %lf",  &MagEphemInfo->UT);
//    fscanf(fp, " %lf",  &MagEphemInfo->Lat);
//    fscanf(fp, " %lf",  &MagEphemInfo->Lon);
//    fscanf(fp, " %lf",  &MagEphemInfo->Rad);
    fscanf(fp, " %lf",  &Lat);
    fscanf(fp, " %lf",  &Lon);
    fscanf(fp, " %lf",  &Rad);
//MagEphemInfo->nAlpha=9;
    for (i=0; i<MagEphemInfo->nAlpha; i++) { fscanf(fp, " %lf", &MagEphemInfo->Lstar[i]); }
    for (i=0; i<MagEphemInfo->nAlpha; i++) { fscanf(fp, " %lf", &MagEphemInfo->Bm[i]); }
    for (i=0; i<MagEphemInfo->nAlpha; i++) { fscanf(fp, " %lf", &MagEphemInfo->I[i]); }
    for (i=0; i<MagEphemInfo->nAlpha; i++) { fscanf(fp, " %lf", &MagEphemInfo->K[i]); }
    for (i=0; i<MagEphemInfo->nAlpha; i++) { fscanf(fp, " %lf", &MagEphemInfo->Sb[i]); }
    for (i=0; i<MagEphemInfo->nAlpha; i++) { fscanf(fp, " %lf", &MagEphemInfo->Tb[i]); }
    for (i=0; i<MagEphemInfo->nAlpha; i++) { fscanf(fp, " %lf", &MagEphemInfo->LMcIlwain[i]); } 

    for (i=0; i<MagEphemInfo->nAlpha; i++) { fscanf(fp, " %lf %lf %lf", &MagEphemInfo->Pmn_gsm[i].x, &MagEphemInfo->Pmn_gsm[i].y, &MagEphemInfo->Pmn_gsm[i].z); }
    for (i=0; i<MagEphemInfo->nAlpha; i++) { fscanf(fp, " %lf %lf %lf", &MagEphemInfo->Pms_gsm[i].x, &MagEphemInfo->Pms_gsm[i].y, &MagEphemInfo->Pms_gsm[i].z); }

    for (i=0; i<MagEphemInfo->nAlpha; i++) { 
        fscanf(fp, " %d", &MagEphemInfo->nShellPoints[i]);
printf("n = %d\n", MagEphemInfo->nShellPoints[i]);
        for (j=0; j<MagEphemInfo->nShellPoints[i]; j++) { 
            fscanf(fp, " %lf %lf %lf", &MagEphemInfo->ShellFootprint_Pn[i][j].x, &MagEphemInfo->ShellFootprint_Pn[i][j].y, &MagEphemInfo->ShellFootprint_Pn[i][j].z );
            fscanf(fp, " %lf %lf %lf", &MagEphemInfo->ShellFootprint_Ps[i][j].x, &MagEphemInfo->ShellFootprint_Ps[i][j].y, &MagEphemInfo->ShellFootprint_Ps[i][j].z );
            fscanf(fp, " %lf %lf %lf", &MagEphemInfo->ShellMirror_Pn[i][j].x, &MagEphemInfo->ShellMirror_Pn[i][j].y, &MagEphemInfo->ShellMirror_Pn[i][j].z );
            fscanf(fp, " %lf %lf %lf", &MagEphemInfo->ShellMirror_Ps[i][j].x, &MagEphemInfo->ShellMirror_Ps[i][j].y, &MagEphemInfo->ShellMirror_Ps[i][j].z );
        }
    }
MagEphemInfo->nAlpha=9;

    fclose( fp );

}
