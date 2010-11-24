#include "Lgm/Lgm_CTrans.h"
#include "Lgm/Lgm_Sgp.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

static char *MonStr1[] = {"January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December"};
//static char *MonStr2[] = {"Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"};

int LgmSgp_TleChecksum( char *Line ) {

    int i, sum, val;

    if ( strlen( Line ) < 69 ) return(-1);

    for (sum=0, i=0; i<68; i++){
        if ( Line[i] == '-' ) val = 1;
        else if ( (Line[i] >= '0') && (Line[i] <= '9') ) val = (int)Line[i]-'0';
        else val = 0;
        sum += val;
    }

    return(sum%10);

}

int LgmSgp_ReadTlesFromFile( char *Filename, int *nTLEs, _SgpTLE *TLEs, int Verbosity ) {

    char    *Line0, *Line1, *Line2, *ptr;
    int     nLines, done, TLEValid, CheckSumRead1, CheckSumComp1;
    int     CheckSumRead2, CheckSumComp2;
    FILE    *fp;

//    *nTLEs = 0;

    if ( (fp = fopen( Filename, "rb" )) == NULL ) {

        if (Verbosity > 0) printf("LgmSgp_ReadTwoLineElements: Couldn't open Two Line Element file: %s\n", Filename);
        return(0);

    } else {

        Line0 = (char *)calloc(256, sizeof(char));
        Line1 = (char *)calloc(256, sizeof(char));
        Line2 = (char *)calloc(256, sizeof(char));

        nLines    = 0;
        done      = FALSE;
        while( !done ) {

            TLEValid = TRUE;

            if ( fgets( Line0, 256, fp ) == NULL) {         // Try to read in Line 0
                done = TRUE;
            } else if ( fgets( Line1, 256, fp ) == NULL) {  // Try to read in Line 1
                done = TRUE;
            } else if ( fgets( Line2, 256, fp ) == NULL) {  // Try to read in Line 2
                done = TRUE;
            }

            if (!done){

                nLines += 3;

                /*
                 * fgets also reads in the newline. Also it appears that the
                 * TLE files you often get are dos format with the CRLF (i.e.
                 * '\r\n') line terminators instead of the more sane linux/unix
                 * LF (i.e. '\n' ) line terminator. So, to fix this lets just
                 * replace any occurances of either with a string terminator
                 * (i.e. a '\0').  Then the string will be properly terminated
                 * if its a dos or linux style txt file.
                 */
                if ( (ptr = strstr(Line0, "\n")) != NULL ) *ptr = '\0'; if ( (ptr = strstr(Line0, "\r")) != NULL ) *ptr = '\0';
                if ( (ptr = strstr(Line1, "\n")) != NULL ) *ptr = '\0'; if ( (ptr = strstr(Line1, "\r")) != NULL ) *ptr = '\0';
                if ( (ptr = strstr(Line2, "\n")) != NULL ) *ptr = '\0'; if ( (ptr = strstr(Line2, "\r")) != NULL ) *ptr = '\0';

                if (Verbosity > 3){
                    printf("Line0 = |%s|\n", Line0);
                    printf("Line1 = |%s|\n", Line1);
                    printf("Line2 = |%s|\n", Line2);
                }


                /*
                 * Now, validate the lines. Must be a valid TLE to get added to
                 * the list of TLEs.
                 */
                CheckSumRead1 = (int)(Line1[68]-'0');
                CheckSumComp1 = LgmSgp_TleChecksum( Line1 );
                
                CheckSumRead2 = (int)(Line2[68]-'0');
                CheckSumComp2 = LgmSgp_TleChecksum( Line2 );

                if ( Line1[0] != '1' ) {
                    if (Verbosity > 0) printf("LgmSgp_ReadTwoLineElements: Format Error in Two Line Element.\n");
                    if (Verbosity > 0) printf("    In file \"%s\", at line number %d, an assumed TLE Line#1 doesn't\n", Filename, nLines-1);
                    if (Verbosity > 0) printf("    start with a `1'.\n");
                    if (Verbosity > 0) printf("    Line1=\"%s\"\n\n", Line1);
                    TLEValid = FALSE;
                }
                if ( CheckSumRead1 != CheckSumComp1 ) {
                    if (Verbosity > 0) printf("LgmSgp_ReadTwoLineElements: Checksum Error in Two Line Element.\n");
                    if (Verbosity > 0) printf("    In file \"%s\", at line number %d, an assumed TLE Line#1 doesn't\n", Filename, nLines-1);
                    if (Verbosity > 0) printf("    have the correct Checksum. Read: %d, Computed: %d.\n", CheckSumRead1, CheckSumComp1);
                    if (Verbosity > 0) printf("    Line1=\"%s\"\n\n", Line1);
                    TLEValid = FALSE;
                }
                if ( Line2[0] != '2' ) {
                    if (Verbosity > 0) printf("LgmSgp_ReadTwoLineElements: Format Error in Two Line Element.\n");
                    if (Verbosity > 0) printf("    In file \"%s\", at line number %d, an assumed TLE Line#2 doesn't\n", Filename, nLines);
                    if (Verbosity > 0) printf("    start with a `2'.\n");
                    if (Verbosity > 0) printf("    Line2=\"%s\"\n\n", Line2);
                    TLEValid = FALSE;
                }
                if ( CheckSumRead2 != CheckSumComp2 ) {
                    if (Verbosity > 0) printf("LgmSgp_ReadTwoLineElements: Checksum Error in Two Line Element.\n");
                    if (Verbosity > 0) printf("    In file \"%s\", at line number %d, an assumed TLE Line#2 doesn't\n", Filename, nLines);
                    if (Verbosity > 0) printf("    have the correct Checksum. Read: %d, Computed: %d.\n", CheckSumRead2, CheckSumComp2);
                    if (Verbosity > 0) printf("    Line2=\"%s\"\n\n", Line2);
                    TLEValid = FALSE;
                }

//                if (TLEValid) {

                    if (Verbosity > 3) printf("\n\t\tAdding valid TLE to list. nTLEs = %d\n", *nTLEs);
                    if (Verbosity > 3) printf("\t\t--------------------------------------\n");

                    // Finally decode the lines and dump into structure -- must be valid TLE at this point
                    Lgm_SgpDecodeTle( Line0, Line1, Line2, &(TLEs[*nTLEs]), Verbosity );

                    ++(*nTLEs);
//                }

            }
            

        } 

        fclose(fp);
        free(Line0);
        free(Line1);
        free(Line2);

    }

    return(1);

}

int LgmSgp_ReadTlesFromStrings( char *Line0, char *Line1, char *Line2, int *nTLEs, _SgpTLE *TLEs, int Verbosity ) {

    char    *ptr;
    int     TLEValid, CheckSumRead1, CheckSumComp1;
    int     CheckSumRead2, CheckSumComp2;

//    *nTLEs = 0;

    /*
     * fgets also reads in the newline. Also it appears that the
     * TLE files you often get are dos format with the CRLF (i.e.
     * '\r\n') line terminators instead of the more sane linux/unix
     * LF (i.e. '\n' ) line terminator. So, to fix this lets just
     * replace any occurances of either with a string terminator
     * (i.e. a '\0').  Then the string will be properly terminated
     * if its a dos or linux style txt file.
     */
    if ( (ptr = strstr(Line0, "\n")) != NULL ) *ptr = '\0'; if ( (ptr = strstr(Line0, "\r")) != NULL ) *ptr = '\0';
    if ( (ptr = strstr(Line1, "\n")) != NULL ) *ptr = '\0'; if ( (ptr = strstr(Line1, "\r")) != NULL ) *ptr = '\0';
    if ( (ptr = strstr(Line2, "\n")) != NULL ) *ptr = '\0'; if ( (ptr = strstr(Line2, "\r")) != NULL ) *ptr = '\0';


    if (Verbosity > 3){
        printf("Line0 = |%s|\n", Line0);
        printf("Line1 = |%s|\n", Line1);
        printf("Line2 = |%s|\n", Line2);
    }


    /*
     * Now, validate the lines. Must be a valid TLE to get added to
     * the list of TLEs.
     */
    CheckSumRead1 = (int)(Line1[68]-'0');
    CheckSumComp1 = LgmSgp_TleChecksum( Line1 );
    
    CheckSumRead2 = (int)(Line2[68]-'0');
    CheckSumComp2 = LgmSgp_TleChecksum( Line2 );

    TLEValid = TRUE;
    if ( Line1[0] != '1' ) {
        if (Verbosity > 0) printf("LgmSgp_ReadTwoLineElements: Format Error in Two Line Element.\n");
        if (Verbosity > 0) printf("    TLE Line#1 doesn't\n" );
        if (Verbosity > 0) printf("    start with a `1'.\n");
        if (Verbosity > 0) printf("    Line1=\"%s\"\n\n", Line1);
        TLEValid = FALSE;
    }
    if ( CheckSumRead1 != CheckSumComp1 ) {
        if (Verbosity > 0) printf("LgmSgp_ReadTwoLineElements: Checksum Error in Two Line Element.\n");
        if (Verbosity > 0) printf("    TLE Line#1 doesn't\n");
        if (Verbosity > 0) printf("    have the correct Checksum. Read: %d, Computed: %d.\n", CheckSumRead1, CheckSumComp1);
        if (Verbosity > 0) printf("    Line1=\"%s\"\n\n", Line1);
//        TLEValid = FALSE;
    }
    if ( Line2[0] != '2' ) {
        if (Verbosity > 0) printf("LgmSgp_ReadTwoLineElements: Format Error in Two Line Element.\n");
        if (Verbosity > 0) printf("    TLE Line#2 doesn't\n");
        if (Verbosity > 0) printf("    start with a `2'.\n");
        if (Verbosity > 0) printf("    Line2=\"%s\"\n\n", Line2);
        TLEValid = FALSE;
    }
    if ( CheckSumRead2 != CheckSumComp2 ) {
        if (Verbosity > 0) printf("LgmSgp_ReadTwoLineElements: Checksum Error in Two Line Element.\n");
        if (Verbosity > 0) printf("    TLE Line#2 doesn't\n");
        if (Verbosity > 0) printf("    have the correct Checksum. Read: %d, Computed: %d.\n", CheckSumRead2, CheckSumComp2);
        if (Verbosity > 0) printf("    Line2=\"%s\"\n\n", Line2);
//        TLEValid = FALSE;
    }

    if (TLEValid) {

        if (Verbosity > 3) printf("\n\t\tAdding valid TLE to list. nTLEs = %d\n", *nTLEs);
        if (Verbosity > 3) printf("\t\t--------------------------------------\n");

        // Finally decode the lines and dump into structure -- must be valid TLE at this point
        Lgm_SgpDecodeTle( Line0, Line1, Line2, &(TLEs[*nTLEs]), Verbosity );

        ++(*nTLEs);
    }

    return(1);

}


/*
 * Given Lines 0-2, decode the vals and populate vars in the _SgpTLE structure
 */
void Lgm_SgpDecodeTle( char *Line0, char *Line1, char *Line2, _SgpTLE *TLE, int Verbosity ) {

    char    str[40], *p;
    int     d2, yy, yyyy;
    int     hh, mm, ss, i, ll;
    long    d1;
    Lgm_CTrans *c = Lgm_init_ctrans(0); // need this for JD calcs

    /*
     *  DECODE LINE 0
     */
    // Save the 3 lines
    strcpy( TLE->Line0, Line0 );
    strcpy( TLE->Line1, Line1 );
    strcpy( TLE->Line2, Line2 );

    // Object name
    strcpy( TLE->Name, Line0 );
    strcpy( TLE->ObjectType, "SAT"); 
    if ( (p = strstr( TLE->Name, " DEB" )) ) {
        if ( ( p[4] == ' ' ) || ( p[4] == '\0' ) ) strcpy( TLE->ObjectType, "DEB");
    } 
    if  ( strstr( TLE->Name, "R/B" ) ) strcpy( TLE->ObjectType, "R/B");
    if (Verbosity > 3) printf("\t\tObject Name     = %s (Object Type = %s)\n", TLE->Name, TLE->ObjectType);

    /*
     *  DECODE LINE 1
     */

    // Object ID number
    strncpy( str, Line1+2, 5 ); str[5] = '\0'; // be careful! strncpy doesnt NULL-terminate strings
    TLE->IdNumber = atoi(str);
    if (Verbosity > 3) printf("\t\tIdNumber        = %d\n", TLE->IdNumber);

    // Elset Classification
    TLE->ElsetClass = Line1[7];
    if (Verbosity > 3) printf("\t\tElsetClass      = %c\n", TLE->ElsetClass);

    //  International Designator
    strncpy( TLE->IntDesig, Line1+9, 8 ); TLE->IntDesig[8] = '\0';
    yy = (int)(TLE->IntDesig[0] -'0')*10 + (int)(TLE->IntDesig[1] -'0');
    yyyy = (yy<50) ? 2000+yy : 1900+yy;
    sprintf( TLE->IntDesig2, "%4d-%s", yyyy, &(TLE->IntDesig[2]) );
    // strip off trailing spaces if any..
    if ( (ll = strlen(TLE->IntDesig2)) > 0 ){
        for (i=ll-1; i>=0; i--){
            if (TLE->IntDesig2[i] == ' ') {
                TLE->IntDesig2[i] = '\0';
            } else {
                break;
            }
        }
    }
    if (Verbosity > 3) printf("\t\tIntDesig        = %s  (%s)\n", TLE->IntDesig, TLE->IntDesig2);

    // Element Set Epoch (UTC) (also decode into other date formats)
    strncpy( str, Line1+18, 14 ); str[14] = '\0';
    sscanf( str, "%lf", &(TLE->ElementSetEpoch));
    if (Verbosity > 3) printf("\t\tElementSetEpoch = %14.8lf\n", TLE->ElementSetEpoch);
    TLE->Date = (long int)TLE->ElementSetEpoch;
    if ( ((int)(TLE->Date/1000.0)) < 50 ) TLE->Date += 2000000;
    d1 = (int)TLE->ElementSetEpoch;
    TLE->UT = (TLE->ElementSetEpoch - (double)d1)*24.0;
    if (Verbosity > 3) printf("\t\t          Date  = %ld\n", TLE->Date);
    if (Verbosity > 3) printf("\t\t            UT  = %lf\n", TLE->UT);
    Lgm_Doy( TLE->Date, &TLE->Year, &TLE->Month, &TLE->Day, &TLE->Doy);
    if (Verbosity > 3) printf("\t\t          Year  = %d\n", TLE->Year);
    if (Verbosity > 3) printf("\t\t         Month  = %d (%s)\n", TLE->Month, MonStr1[TLE->Month-1]);
    if (Verbosity > 3) printf("\t\t           Day  = %d\n", TLE->Day);
    if (Verbosity > 3) printf("\t\t           Doy  = %d\n", TLE->Doy);
    Lgm_DayOfWeek( TLE->Year, TLE->Month, TLE->Day, TLE->Dow);
    if (Verbosity > 3) printf("\t\t           Dow  = %s\n", TLE->Dow);
    TLE->JD = Lgm_JD( TLE->Year, TLE->Month, TLE->Day, TLE->UT, LGM_TIME_SYS_UTC, c );
    if (Verbosity > 3) printf("\t\t    Julian Date = %.12lf\n", TLE->JD);
    TLE->YYYYDDDdFRAC = TLE->Year*1000.0 + (double)TLE->Doy + TLE->UT/24.0;
    if (Verbosity > 3) printf("\t\t   YYYYDDD.FRAC = %.12lf\n", TLE->YYYYDDDdFRAC);
    Lgm_UT_to_HMS( TLE->UT, &hh, &mm, &ss );
    sprintf(TLE->EpochStr, "%4d%02d%02d %02d:%02d:%02d", TLE->Year, TLE->Month, TLE->Day, hh, mm, ss );
    if (Verbosity > 3) printf("\t\t       EpochStr = %s\n", TLE->EpochStr);
        

    //  1st Derivative of the Mean Motion with respect to Time
    strncpy( str, Line1+33, 10 ); str[10] = '\0';
    sscanf( str, "%lf", &(TLE->dMMdT1));
    if (Verbosity > 3) printf("\t\tdMMdT1          = %14.8lf\n", TLE->dMMdT1);
    
    //  2nd Derivative of the Mean Motion with respect to Time
    strncpy( str, Line1+44, 6 ); str[6] = '\0'; d1 = atol( str ); // mantissa (*1e5 because theres no leading decimal yet)
    strncpy( str, Line1+50, 2 ); str[2] = '\0'; d2 = atoi( str ); // exponent
    TLE->dMMdT2 = (double)d1*pow(10.0, (double)d2-5.0 );
    if (Verbosity > 3) printf("\t\tdMMdT2          = %14.5e\n", TLE->dMMdT2);

    //  B* Drag Term
    strncpy( str, Line1+53, 6 ); str[6] = '\0'; d1 = atol( str ); // mantissa (*1e5 because theres no leading decimal yet)
    strncpy( str, Line1+59, 2 ); str[2] = '\0'; d2 = atoi( str ); // exponent
    TLE->BstarDrag = (double)d1*pow(10.0, (double)d2-5.0 );
    if (Verbosity > 3) printf("\t\tBstarDrag       = %14.5e\n", TLE->BstarDrag);

    // Element Set Type
    TLE->ElementSetType = (int)(Line1[62] - '0');
    if (Verbosity > 3) printf("\t\tElementSetType  = %d\n", TLE->ElementSetType);
    
    // Element Number
    strncpy( str, Line1+64, 4 ); str[4] = '\0';
    TLE->ElementSetNum = atoi(str);
    if (Verbosity > 3) printf("\t\tElementSetNum   = %d\n", TLE->ElementSetNum);
    
    // Line1 Checksum
    TLE->Line1CheckSum = (int)(Line1[68] - '0');
    if (Verbosity > 3) printf("\t\tLine1CheckSum   = %d (passed)\n", TLE->Line1CheckSum);

    /*
     *  DECODE LINE 2
     */

    //  Orbit Inclination (degrees)
    strncpy( str, Line2+8, 8 ); str[8] = '\0';
    sscanf( str, "%lf", &(TLE->Inclination));
    if (Verbosity > 3) printf("\t\tInclination     = %14.8lf (degrees)\n", TLE->Inclination);

    // Right Ascension of Ascending Node (degrees)
    strncpy( str, Line2+17, 8 ); str[8] = '\0';
    sscanf( str, "%lf", &(TLE->RAofAscNode));
    if (Verbosity > 3) printf("\t\tRAofAscNode     = %14.8lf (degrees)\n", TLE->RAofAscNode);
    
    // Eccentricity
    strncpy( str, Line2+26, 7 ); str[7] = '\0'; d1 = atol( str ); // mantissa (*1e7 because theres no leading decimal yet)
    TLE->Eccentricity = (double)d1*1e-7;
    if (Verbosity > 3) printf("\t\tEccentricity    = %14.8lf\n", TLE->Eccentricity);
    
    // Argument of Perigee (degrees)
    strncpy( str, Line2+34, 8 ); str[8] = '\0';
    sscanf( str, "%lf", &(TLE->ArgOfPerigee));
    if (Verbosity > 3) printf("\t\tArgOfPerigee    = %14.8lf (degrees)\n", TLE->ArgOfPerigee);
    
    // Mean Anomaly (degrees)
    strncpy( str, Line2+43, 8 ); str[8] = '\0';
    sscanf( str, "%lf", &(TLE->MeanAnomaly));
    if (Verbosity > 3) printf("\t\tMeanAnomaly     = %14.8lf (degrees)\n", TLE->MeanAnomaly);
    
    // Mean Motion (revolutions/day)
    strncpy( str, Line2+52, 11 ); str[11] = '\0';
    sscanf( str, "%lf", &(TLE->MeanMotion));
    TLE->Period = 1440.0/TLE->MeanMotion; // period in minutes
    if (Verbosity > 3) printf("\t\tMeanMotion      = %14.8lf (revs/day)  (Period = %g min)\n", TLE->MeanMotion, TLE->Period);
    
    // Revolution Number at Epoch
    strncpy( str, Line2+63, 5 ); str[5] = '\0';
    sscanf( str, "%d", &(TLE->RevNumAtEpoch));
    if (Verbosity > 3) printf("\t\tRevNumAtEpoch   = %d\n", TLE->RevNumAtEpoch);
    
    // Line2 Checksum
    TLE->Line2CheckSum = (int)(Line2[68] - '0');
    if (Verbosity > 3) printf("\t\tLine2CheckSum   = %d (passed)\n", TLE->Line2CheckSum);

    Lgm_free_ctrans( c );

}


/* -----------------------------------------------------------------------------
*
* function gstime
*
* this function finds the greenwich sidereal time.
*
* author : david vallado 719-573-2600 1 mar 2001
*
* inputs description range / units
* jdut1 - julian date in ut1 days from 4713 bc
*
* outputs :
* gstime - greenwich sidereal time 0 to 2pi rad
*
* locals :
* temp - temporary variable for doubles rad
* tut1 - julian centuries from the
* jan 1, 2000 12 h epoch (ut1)
*
* coupling :
* none
*
* references :
* vallado 2004, 191, eq 3-45
* --------------------------------------------------------------------------- */
double LgmSgp_gstime( double jdut1) {

    const double    deg2rad = M_PI / 180.0;
    double          temp, tut1;

    tut1 = (jdut1 - 2451545.0) / 36525.0;
    temp = -6.2e-6*tut1*tut1*tut1 + 0.093104*tut1*tut1 + (876600.0*3600 + 8640184.812866)*tut1 + 67310.54841; // sec
    temp = fmod(temp*deg2rad / 240.0, M_2PI); //360/86400 = 1/240, to deg, to rad

    // ------------------------ check quadrants ---------------------
    if (temp < 0.0) temp += M_2PI;

    return( temp );

} // end gstime


/* -----------------------------------------------------------------------------
*
* procedure dpper
*
* this procedure provides deep space long period periodic contributions
* to the mean elements. by design, these periodics are zero at epoch.
* this used to be dscom which included initialization, but it's really a
* recurring function.
*
* author : david vallado 719-573-2600 28 jun 2005
*
* inputs (used):
*     from structure, s: e3, ee2, peo, pgho, pho, pinco, plo, se2, se3, sgh2, sgh3, sgh4, 
*                     sh2, sh3, si2, si3, sl2, sl3, sl4, t, xh2, xh3, xi2, xi3, xl2, 
*                     xl3, xl4, zmol, zmos

*     inclo - inclination - needed for lyddane modification (pass by value -- not changed)
*     init  - init flag (pass by value -- not changed)
*     ep    - eccentricity 0.0 - 1.0
*     nodep - right ascension of ascending node
*     argpp - argument of perigee
*     mp    - mean anomaly
*
* outputs :
*     ep    - eccentricity 0.0 - 1.0
*     inclp - inclination
*     nodep - right ascension of ascending node
*     argpp - argument of perigee
*     mp    - mean anomaly
*
----------------------------------------------------------------------------*/
void LgmSgp_dpper( double inclo, char init, double *ep, double *inclp, double *nodep, double *argpp, double *mp, _SgpInfo *s ){


    /* --------------------- local variables ------------------------ */
    char        ildm;
    double      alfdp, betdp, cosip, cosop, dalf, dbet, dls;
    double      f2, f3, pe, pgh, ph, pinc, pl;
    double      sel, ses, sghl, sghs, shll, shs, sil;
    double      sinip, sinop, sinzf, sis, sll, sls, xls;
    double      xnoh, zf, zm, zel, zes, znl, zns, t;

    /* ---------------------- constants ----------------------------- */
    zns = 1.19459e-5;
    zes = 0.01675;
    znl = 1.5835218e-4;
    zel = 0.05490;

    /* --------------- calculate time varying periodics ----------- */
    t = s->t;
    zm = s->zmos + zns*t;

    // be sure that the initial call has time set to zero
    if (s->init == 'y') zm = s->zmos;
    zf     = zm + 2.0*zes*sin(zm);
    sinzf = sin(zf);
    f2    = 0.5*sinzf*sinzf - 0.25;
    f3    = -0.5*sinzf*cos(zf);
    ses   = s->se2* f2 + s->se3*f3;
    sis   = s->si2*f2 + s->si3*f3;
    sls   = s->sl2*f2 + s->sl3*f3 + s->sl4*sinzf;
    sghs  = s->sgh2*f2 + s->sgh3*f3 + s->sgh4*sinzf;
    shs   = s->sh2*f2 + s->sh3*f3;
    zm    = s->zmol + znl*t;
    if (s->init == 'y') zm = s->zmol;
    zf    = zm + 2.0*zel*sin(zm);
    sinzf = sin(zf);
    f2    = 0.5*sinzf*sinzf - 0.25;
    f3    = -0.5*sinzf*cos(zf);
    sel   = s->ee2*f2 + s->e3*f3;
    sil   = s->xi2*f2 + s->xi3*f3;
    sll   = s->xl2*f2 + s->xl3*f3 + s->xl4*sinzf;
    sghl  = s->xgh2*f2 + s->xgh3*f3 + s->xgh4*sinzf;
    shll  = s->xh2*f2 + s->xh3*f3;
    pe    = ses + sel;
    pinc  = sis + sil;
    pl    = sls + sll;
    pgh   = sghs + sghl;
    ph    = shs + shll;


    if (s->init == 'n') {

        // 0.2 rad = 11.45916 deg
        // sgp4fix for lyddane choice
        // add next three lines to set up use of original inclination per strn3 ver
        ildm = 'y';
        if (inclo >= 0.2) ildm = 'n';
        pe    -= s->peo;
        pinc  -= s->pinco;
        pl    -= s->plo;
        pgh   -= s->pgho;
        ph    -= s->pho;
        *inclp += pinc;
        *ep    += pe;
        sinip = sin(*inclp);
        cosip = cos(*inclp);
        /* ----------------- apply periodics directly ------------ */
        // sgp4fix for lyddane choice
        // strn3 used original inclination - this is technically feasible
        // gsfc used perturbed inclination - also technically feasible
        // probably best to readjust the 0.2 limit value and limit discontinuity
        // use next line for original strn3 approach and original inclination
        // if (inclo >= 0.2)
        // use next line for gsfc version and perturbed inclination
        if (*inclp >= 0.2) {

            ph    /= sinip;
            pgh   -= cosip*ph;
            *argpp += pgh;
            *nodep += ph;
            *mp    += pl;

        } else {

            /* ---- apply periodics with lyddane modification ---- */
            sinop = sin(*nodep);
            cosop = cos(*nodep);
            alfdp = sinip*sinop;
            betdp = sinip*cosop;
            dalf  =  ph*cosop + pinc*cosip*sinop;
            dbet  = -ph*sinop + pinc*cosip*cosop;

            alfdp += dalf;
            betdp += dbet;
            *nodep = fmod(*nodep, M_2PI);
            xls   = *mp + *argpp + cosip* (*nodep);
            dls   = pl + pgh - pinc*( *nodep)*sinip;
            xls   += dls;
            xnoh  = *nodep;
            *nodep = atan2(alfdp, betdp);
            if (fabs(xnoh - *nodep) > M_PI) {
                if ( *nodep < xnoh )
                    *nodep = *nodep + M_2PI;
                else
                    *nodep = *nodep - M_2PI;
            }
            *mp   += pl;
            *argpp = xls - *mp - cosip * *nodep;

        }

    } // if s->init == 'n'

    //#include "debug1.cpp"


} // end LgmSgp_dpper



/*-----------------------------------------------------------------------------
*
* procedure dscom
*
* this procedure provides deep space common items used by both the secular
* and periodics subroutines. input is provided as shown. this routine
* used to be called dpper, but the functions inside weren't well organized.
*
* author : david vallado 719-573-2600 28 jun 2005
*
* inputs :
* epoch -
* ep - eccentricity
* argpp - argument of perigee
* tc -
* inclp - inclination
* nodep - right ascension of ascending node
* np - mean motion
*
* outputs :
*     sinim , cosim , sinomm , cosomm , snodm , cnodm
*     day -
*     e3 -
*     ee2 -
*     em - eccentricity
*     emsq - eccentricity squared
*     gam -
*     peo -
*     pgho -
*     pho -
*     pinco -
*     plo -
*     rtemsq -
*     se2, se3 -
*     sgh2, sgh3, sgh4 -
*     sh2, sh3, si2, si3, sl2, sl3, sl4 -
*     s1, s2, s3, s4, s5, s6, s7 -
*     ss1, ss2, ss3, ss4, ss5, ss6, ss7, sz1, sz2, sz3 -
*     sz11, sz12, sz13, sz21, sz22, sz23, sz31, sz32, sz33 -
*     xgh2, xgh3, xgh4, xh2, xh3, xi2, xi3, xl2, xl3, xl4 -
*     nm - mean motion
*     z1, z2, z3, z11, z12, z13, z21, z22, z23, z31, z32, z33 -
*     zmol -
*     zmos -
*
* locals :
* a1, a2, a3, a4, a5, a6, a7, a8, a9, a10 -
* betasq -
* cc -
* ctem, stem -
* x1, x2, x3, x4, x5, x6, x7, x8 -
* xnodce -
* xnoi -
* zcosg , zsing , zcosgl , zsingl , zcosh , zsinh , zcoshl , zsinhl ,
* zcosi , zsini , zcosil , zsinil ,
* zx -
* zy -
*
* coupling :
* none.
*
* references :
* hoots, roehrich, norad spacetrack report #3 1980
* hoots, norad spacetrack report #6 1986
* hoots, schumacher and glover 2004
* vallado, crawford, hujsak, kelso 2006
----------------------------------------------------------------------------*/
void LgmSgp_dscom( double epoch, double ep, double argpp, double tc, double inclp, double nodep, double np,
                    double *snodm, double *cnodm, double *sinim, double *cosim, double *sinomm,
                    double *cosomm,double *day, double *e3, double *ee2, double *em,
                    double *emsq, double *gam, double *peo, double *pgho, double *pho,
                    double *pinco, double *plo, double *rtemsq, double *se2, double *se3,
                    double *sgh2, double *sgh3, double *sgh4, double *sh2, double *sh3,
                    double *si2, double *si3, double *sl2, double *sl3, double *sl4,
                    double *s1, double *s2, double *s3, double *s4, double *s5,
                    double *s6, double *s7, double *ss1, double *ss2, double *ss3,
                    double *ss4, double *ss5, double *ss6, double *ss7, double *sz1,
                    double *sz2, double *sz3, double *sz11, double *sz12, double *sz13,
                    double *sz21, double *sz22, double *sz23, double *sz31, double *sz32,
                    double *sz33, double *xgh2, double *xgh3, double *xgh4, double *xh2,
                    double *xh3, double *xi2, double *xi3, double *xl2, double *xl3,
                    double *xl4, double *nm, double *z1, double *z2, double *z3,
                    double *z11, double *z12, double *z13, double *z21, double *z22,
                    double *z23, double *z31, double *z32, double *z33, double *zmol, double *zmos) {


    /* -------------------------- constants ------------------------- */
    const double zes    = 0.01675;
    const double zel    = 0.05490;
    const double c1ss   = 2.9864797e-6;
    const double c1l    = 4.7968065e-7;
    const double zsinis = 0.39785416;
    const double zcosis = 0.91744867;
    const double zcosgs = 0.1945905;
    const double zsings = -0.98088458;

    /* --------------------- local variables ------------------------ */
    int     lsflg;
    double  a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, betasq, cc, ctem, stem;
    double  x1, x2, x3, x4, x5, x6, x7, x8, xnodce, xnoi, zcosg, zcosgl, zcosh, zcoshl;
    double  zcosi, zcosil, zsing, zsingl, zsinh, zsinhl, zsini, zsinil, zx, zy;


    *nm     = np;
    *em     = ep;
    *snodm  = sin(nodep);
    *cnodm  = cos(nodep);
    *sinomm = sin(argpp);
    *cosomm = cos(argpp);
    *sinim  = sin(inclp);
    *cosim  = cos(inclp);
    *emsq   = *em * *em;
    betasq = 1.0 - *emsq;
    *rtemsq = sqrt(betasq);

    /* ----------------- initialize lunar solar terms --------------- */
    *peo    = 0.0;
    *pinco  = 0.0;
    *plo    = 0.0;
    *pgho   = 0.0;
    *pho    = 0.0;
    *day    = epoch + 18261.5 + tc / 1440.0;
    xnodce = fmod(4.5236020 - 9.2422029e-4 * *day, M_2PI);
    stem   = sin(xnodce);
    ctem   = cos(xnodce);
    zcosil = 0.91375164 - 0.03568096 * ctem;
    zsinil = sqrt(1.0 - zcosil * zcosil);
    zsinhl = 0.089683511 * stem / zsinil;
    zcoshl = sqrt(1.0 - zsinhl * zsinhl);
    *gam    = 5.8351514 + 0.0019443680 * *day;
    zx     = 0.39785416 * stem / zsinil;
    zy     = zcoshl * ctem + 0.91744867 * zsinhl * stem;
    zx     = atan2(zx, zy);
    zx     = *gam + zx - xnodce;
    zcosgl = cos(zx);
    zsingl = sin(zx);


    /* ------------------------- do solar terms --------------------- */
    zcosg = zcosgs;
    zsing = zsings;
    zcosi = zcosis;
    zsini = zsinis;
    zcosh = *cnodm;
    zsinh = *snodm;
    cc = c1ss;
    xnoi = 1.0/(*nm);
    for (lsflg = 1; lsflg <= 2; lsflg++) {
        a1  = zcosg*zcosh + zsing*zcosi*zsinh;
        a3  = -zsing*zcosh + zcosg*zcosi*zsinh;
        a7  = -zcosg*zsinh + zsing*zcosi*zcosh;
        a8  = zsing*zsini;
        a9  = zsing*zsinh + zcosg*zcosi*zcosh;
        a10 = zcosg*zsini;
        a2  = *cosim*a7 + *sinim*a8;
        a4  = *cosim*a9 + *sinim*a10;
        a5  = -*sinim*a7 + *cosim*a8;
        a6  = -*sinim*a9 + *cosim*a10;
        x1  = a1 * *cosomm + a2 * *sinomm;
        x2  = a3 * *cosomm + a4 * *sinomm;
        x3  = -a1 * *sinomm + a2 * *cosomm;
        x4  = -a3 * *sinomm + a4 * *cosomm;
        x5  = a5 * *sinomm;
        x6  = a6 * *sinomm;
        x7  = a5 * *cosomm;
        x8  = a6 * *cosomm;
        *z31 = 12.0*x1*x1 - 3.0*x3*x3;
        *z32 = 24.0*x1*x2 - 6.0*x3*x4;
        *z33 = 12.0*x2*x2 - 3.0*x4*x4;
        *z1  = 3.0*(a1*a1 + a2*a2) + *z31 * *emsq;
        *z2  = 6.0*(a1*a3 + a2*a4) + *z32 * *emsq;
        *z3  = 3.0*(a3*a3 + a4*a4) + *z33 * *emsq;
        *z11 = -6.0*a1*a5 + *emsq*(-24.0*x1*x7-6.0*x3*x5);
        *z12 = -6.0*(a1*a6 + a3*a5) + *emsq * (-24.0*(x2*x7 + x1*x8) - 6.0*(x3*x6 + x4*x5));
        *z13 = -6.0*a3*a6 + *emsq*(-24.0*x2*x8 - 6.0*x4*x6);
        *z21 = 6.0*a2*a5 + *emsq*(24.0*x1*x5 - 6.0*x3*x7);
        *z22 = 6.0*(a4*a5 + a2*a6) + *emsq * (24.0*(x2*x5 + x1*x6) - 6.0*(x4*x7 + x3*x8));
        *z23 = 6.0*a4*a6 + *emsq*(24.0*x2*x6 - 6.0*x4*x8);
        *z1  = *z1 + *z1 + betasq * *z31;
        *z2  = *z2 + *z2 + betasq * *z32;
        *z3  = *z3 + *z3 + betasq * *z33;
        *s3  = cc*xnoi;
        *s2  = -0.5 * *s3 / *rtemsq;
        *s4  = *s3 * *rtemsq;
        *s1  = -15.0 * *em * *s4;
        *s5  = x1*x3 + x2*x4;
        *s6  = x2*x3 + x1*x4;
        *s7  = x2*x4 - x1*x3;


        /* ----------------------- do lunar terms ------------------- */
        if (lsflg == 1) {
            *ss1 = *s1; *ss2 = *s2; *ss3 = *s3; *ss4 = *s4; *ss5 = *s5; *ss6 = *s6; *ss7 = *s7;
            *sz1 = *z1; *sz2 = *z2; *sz3 = *z3;
            *sz11 = *z11; *sz12 = *z12; *sz13 = *z13; *sz21 = *z21; *sz22 = *z22; *sz23 = *z23; *sz31 = *z31; *sz32 = *z32; *sz33 = *z33;
            zcosg = zcosgl; zsing = zsingl; zcosi = zcosil; zsini = zsinil;
            zcosh = zcoshl*(*cnodm) + zsinhl*(*snodm);
            zsinh = *snodm*zcoshl - *cnodm*zsinhl;
            cc = c1l;
        }

    }

    *zmol = fmod(4.7199672 + 0.22997150 * *day - *gam, M_2PI);
    *zmos = fmod(6.2565837 + 0.017201977 * *day, M_2PI);

    /* ------------------------ do solar terms ---------------------- */
    *se2  = 2.0 * *ss1 * *ss6;
    *se3  = 2.0 * *ss1 * *ss7;
    *si2  = 2.0 * *ss2 * *sz12;
    *si3  = 2.0 * *ss2*(*sz13 - *sz11);
    *sl2  = -2.0 * *ss3 * *sz2;
    *sl3  = -2.0 * *ss3*(*sz3 - *sz1);
    *sl4  = -2.0 * *ss3*(-21.0 - 9.0 * *emsq)*zes;
    *sgh2 = 2.0 * *ss4 * *sz32;
    *sgh3 = 2.0 * *ss4*(*sz33 - *sz31);
    *sgh4 = -18.0 * *ss4*zes;
    *sh2  = -2.0 * *ss2 * *sz22;
    *sh3  = -2.0 * *ss2*(*sz23 - *sz21);

    /* ------------------------ do lunar terms ---------------------- */
    *ee2  = 2.0 * *s1 * *s6;
    *e3   = 2.0 * *s1 * *s7;
    *xi2  = 2.0 * *s2 * *z12;
    *xi3  = 2.0 * *s2*(*z13 - *z11);
    *xl2  = -2.0 * *s3 * *z2;
    *xl3  = -2.0 * *s3*(*z3 - *z1);
    *xl4  = -2.0 * *s3*(-21.0 - 9.0 * *emsq)*zel;
    *xgh2 = 2.0 * *s4 * *z32;
    *xgh3 = 2.0 * *s4*(*z33 - *z31);
    *xgh4 = -18.0 * *s4*zel;
    *xh2  = -2.0 * *s2 * *z22;
    *xh3  = -2.0 * *s2*(*z23 - *z21);

    //#include "debug2.cpp"

} // end dscom


/*-----------------------------------------------------------------------------
*
* procedure dsinit
*
* this procedure provides deep space contributions to mean motion dot due
* to geopotential resonance with half day and one day orbits.
*
* author : david vallado 719-573-2600 28 jun 2005
*
* inputs :
* cosim, sinim-
* emsq - eccentricity squared
* argpo - argument of perigee
* s1, s2, s3, s4, s5 -
* ss1, ss2, ss3, ss4, ss5 -
* sz1, sz3, sz11, sz13, sz21, sz23, sz31, sz33 -
* t - time
* tc -
* gsto - greenwich sidereal time rad
* mo - mean anomaly
* mdot - mean anomaly dot (rate)
* no - mean motion
* nodeo - right ascension of ascending node
* nodedot - right ascension of ascending node dot (rate)
* xpidot -
* z1, z3, z11, z13, z21, z23, z31, z33 -
* eccm - eccentricity
* argpm - argument of perigee
* inclm - inclination
* mm - mean anomaly
* xn - mean motion
* nodem - right ascension of ascending node
*
* outputs :
* em - eccentricity
* argpm - argument of perigee
* inclm - inclination
* mm - mean anomaly
* nm - mean motion
* nodem - right ascension of ascending node
* irez - flag for resonance 0-none, 1-one day, 2-half day
* atime -
* d2201, d2211, d3210, d3222, d4410, d4422, d5220, d5232, d5421, d5433 -
* dedt -
* didt -
* dmdt -
* dndt -
* dnodt -
* domdt -
* del1, del2, del3 -
* ses , sghl , sghs , sgs , shl , shs , sis , sls
* theta -
* xfact -
* xlamo -
* xli -
* xni
*
* locals :
* ainv2 -
* aonv -
* cosisq -
* eoc -
* f220, f221, f311, f321, f322, f330, f441, f442, f522, f523, f542, f543 -
* g200, g201, g211, g300, g310, g322, g410, g422, g520, g521, g532, g533 -
* sini2 -
* temp -
* temp1 -
* theta -
* xno2 -
*
* coupling :
* LgmSgp_GetGravConst
*
* references :
* hoots, roehrich, norad spacetrack report #3 1980
* hoots, norad spacetrack report #6 1986
* hoots, schumacher and glover 2004
* vallado, crawford, hujsak, kelso 2006
----------------------------------------------------------------------------*/


void LgmSgp_dsinit( int whichconst, double cosim, double emsq, double argpo,
                        double s1, double s2, double s3, double s4, double s5, double sinim, double ss1, 
                        double ss2, double ss3, double ss4, double ss5, double sz1, double sz3, double sz11, 
                        double sz13, double sz21, double sz23, double sz31, double sz33, double t, double tc, 
                        double gsto, double mo, double mdot, double no, double nodeo, double nodedot, 
                        double xpidot, double z1, double z3, double z11, double z13, double z21, double z23, 
                        double z31, double z33, double ecco, double eccsq, 

                        double *em, double *argpm, double *inclm, double *mm, double *nm, double *nodem, 
                        int *irez, double *atime, double *d2201, double *d2211, double *d3210, double *d3222, 
                        double *d4410, double *d4422, double *d5220, double *d5232, double *d5421, double *d5433, 
                        double *dedt, double *didt, double *dmdt, double *dndt, double *dnodt, double *domdt, 
                        double *del1, double *del2, double *del3, double *xfact, double *xlamo, double *xli, double *xni) {



    /* --------------------- local variables ------------------------ */
    double  ainv2, aonv=0.0, cosisq, eoc, f220, f221, f311;
    double  f321, f322, f330, f441, f442, f522, f523;
    double  f542, f543, g200, g201, g211, g300, g310;
    double  g322, g410, g422, g520, g521, g532, g533;
    double  ses, sgs, sghl, sghs, shs, shll, sis;
    double  sini2, sls, temp, temp1, theta, xno2, q22;
    double  q31, q33, root22, root44, root54, rptim, root32;
    double  root52, x2o3, xke, znl, emo, zns, emsqo;
    double  tumin, radiusearthkm, j2, j3, j4, j3oj2;



    q22    = 1.7891679e-6;
    q31    = 2.1460748e-6;
    q33    = 2.2123015e-7;
    root22 = 1.7891679e-6;
    root44 = 7.3636953e-9;
    root54 = 2.1765803e-9;
    rptim  = 4.37526908801129966e-3; // this equates to 7.29211514668855e-5 rad/sec
    root32 = 3.7393792e-7;
    root52 = 1.1428639e-7;
    x2o3   = 2.0 / 3.0;
    znl    = 1.5835218e-4;
    zns    = 1.19459e-5;



    // sgp4fix identify constants and allow alternate values
    LgmSgp_GetGravConst( whichconst, &tumin, &radiusearthkm, &xke, &j2, &j3, &j4, &j3oj2 );


    /* -------------------- deep space initialization ------------ */
    *irez = 0;
    if ((*nm < 0.0052359877) && (*nm > 0.0034906585))        *irez = 1;
    if ((*nm >= 8.26e-3) && (*nm <= 9.24e-3) && (*em >= 0.5)) *irez = 2;


    /* ------------------------ do solar terms ------------------- */
    ses  = ss1*zns*ss5;
    sis  = ss2*zns*(sz11 + sz13);
    sls  = -zns*ss3*(sz1 + sz3 - 14.0 - 6.0*emsq);
    sghs = ss4*zns*(sz31 + sz33 - 6.0);
    shs  = -zns*ss2*(sz21 + sz23);
    // sgp4fix for 180 deg incl
    if ((*inclm < 5.2359877e-2) || (*inclm > M_PI - 5.2359877e-2)) shs = 0.0;
    if (sinim != 0.0) shs = shs / sinim;
    sgs = sghs - cosim*shs;



    /* ------------------------- do lunar terms ------------------ */
    *dedt = ses + s1*znl*s5;
    *didt = sis + s2*znl*(z11 + z13);
    *dmdt = sls - znl*s3*(z1 + z3 - 14.0 - 6.0*emsq);
    sghl = s4*znl*(z31 + z33 - 6.0);
    shll = -znl*s2*(z21 + z23);
    // sgp4fix for 180 deg incl
    if ((*inclm < 5.2359877e-2) || (*inclm > M_PI - 5.2359877e-2)) shll = 0.0;
    *domdt = sgs + sghl;
    *dnodt = shs;
    if (sinim != 0.0) {
        *domdt = *domdt - cosim / sinim*shll;
        *dnodt = *dnodt + shll / sinim;
    }



    /* ----------- calculate deep space resonance effects -------- */
    *dndt  = 0.0;
    theta = fmod(gsto + tc*rptim, M_2PI);
    *em    += *dedt*t;
    *inclm += *didt*t;
    *argpm += *domdt*t;
    *nodem += *dnodt*t;
    *mm    += *dmdt*t;
    // sgp4fix for negative inclinations
    // the following if statement should be commented out
    //if (*inclm < 0.0)
    // {
    // *inclm = -*inclm;
    // *argpm = *argpm - pi;
    // *nodem = *nodem + pi;
    // }



    /* -------------- initialize the resonance terms ------------- */
    if (*irez != 0) {

        aonv = pow(*nm/xke, x2o3);

        /* ---------- geopotential resonance for 12 hour orbits ------ */
        if (*irez == 2) {
            cosisq = cosim*cosim;
            emo    = *em;
            *em     = ecco;
            emsqo  = emsq;
            emsq   = eccsq;
            eoc    = *em*emsq;
            g201   = -0.306 - (*em - 0.64)*0.440;
            if (*em <= 0.65) {
                g211 = 3.616 - 13.2470 * *em + 16.2900*emsq;
                g310 = -19.302 + 117.3900 * *em - 228.4190*emsq + 156.5910*eoc;
                g322 = -18.9068 + 109.7927 * *em - 214.6334*emsq + 146.5816*eoc;
                g410 = -41.122 + 242.6940 * *em - 471.0940*emsq + 313.9530*eoc;
                g422 = -146.407 + 841.8800 * *em - 1629.014*emsq + 1083.4350*eoc;
                g520 = -532.114 + 3017.977 * *em - 5740.032*emsq + 3708.2760*eoc;
            } else {
                g211 = -72.099 + 331.819 * *em - 508.738*emsq + 266.724*eoc;
                g310 = -346.844 + 1582.851 * *em - 2415.925*emsq + 1246.113*eoc;
                g322 = -342.585 + 1554.908 * *em - 2366.899*emsq + 1215.972*eoc;
                g410 = -1052.797 + 4758.686 * *em - 7193.992*emsq + 3651.957*eoc;
                g422 = -3581.690 + 16178.110 * *em - 24462.770*emsq + 12422.520*eoc;
                if (*em > 0.715) g520 =-5149.66 + 29936.92 * *em - 54087.36*emsq + 31324.56*eoc;
                else            g520 = 1464.74 - 4664.75 * *em + 3763.64*emsq;
            }


            if (*em < 0.7) {
                g533 = -919.22770 + 4988.6100 * *em - 9064.7700*emsq + 5542.21*eoc;
                g521 = -822.71072 + 4568.6173 * *em - 8491.4146*emsq + 5337.524*eoc;
                g532 = -853.66600 + 4690.2500 * *em - 8624.7700*emsq + 5341.4*eoc;
            } else {
                g533 =-37995.780 + 161616.52 * *em - 229838.20*emsq + 109377.94*eoc;
                g521 =-51752.104 + 218913.95 * *em - 309468.16*emsq + 146349.42*eoc;
                g532 =-40023.880 + 170470.89 * *em - 242699.48*emsq + 115605.82*eoc;
            }


            sini2= sinim*sinim;
            f220 = 0.75*(1.0 + 2.0*cosim+cosisq);
            f221 = 1.5*sini2;
            f321 = 1.875*sinim*(1.0 - 2.0*cosim - 3.0*cosisq);
            f322 = -1.875*sinim*(1.0 + 2.0*cosim - 3.0*cosisq);
            f441 = 35.0*sini2*f220;
            f442 = 39.3750*sini2*sini2;
            f522 = 9.84375*sinim*(sini2*(1.0 - 2.0*cosim- 5.0*cosisq) + 0.33333333*(-2.0 + 4.0*cosim + 6.0*cosisq) );
            f523 = sinim*(4.92187512*sini2*(-2.0 - 4.0*cosim + 10.0*cosisq) + 6.56250012*(1.0+2.0*cosim - 3.0*cosisq));
            f542 = 29.53125*sinim*(2.0 - 8.0*cosim+cosisq * (-12.0 + 8.0*cosim + 10.0*cosisq));
            f543 = 29.53125*sinim*(-2.0 - 8.0*cosim+cosisq * (12.0 + 8.0*cosim - 10.0*cosisq));
            xno2  = *nm * *nm;
            ainv2 = aonv*aonv;
            temp1 = 3.0*xno2*ainv2;
            temp  = temp1*root22;
            *d2201 = temp*f220*g201;
            *d2211 = temp*f221*g211;
            temp1 = temp1*aonv;
            temp  = temp1*root32;
            *d3210 = temp*f321*g310;
            *d3222 = temp*f322*g322;
            temp1 = temp1*aonv;
            temp  = 2.0*temp1*root44;
            *d4410 = temp*f441*g410;
            *d4422 = temp*f442*g422;
            temp1 = temp1*aonv;
            temp  = temp1*root52;
            *d5220 = temp*f522*g520;
            *d5232 = temp*f523*g532;
            temp  = 2.0*temp1*root54;
            *d5421 = temp*f542*g521;
            *d5433 = temp*f543*g533;
            *xlamo = fmod(mo + nodeo + nodeo-theta - theta, M_2PI);
            *xfact = mdot + *dmdt + 2.0*(nodedot + *dnodt - rptim) - no;
            *em    = emo;
            emsq  = emsqo;
        }



        /* ---------------- synchronous resonance terms -------------- */
        if (*irez == 1) {
            g200  = 1.0 + emsq*(-2.5 + 0.8125*emsq);
            g310  = 1.0 + 2.0*emsq;
            g300  = 1.0 + emsq*(-6.0 + 6.60937*emsq);
            f220  = 0.75*(1.0 + cosim)*(1.0 + cosim);
            f311  = 0.9375*sinim*sinim*(1.0 + 3.0*cosim) - 0.75*(1.0 + cosim);
            f330  = 1.0 + cosim;
            f330  = 1.875*f330*f330*f330;
            *del1  = 3.0 * *nm * *nm*aonv*aonv;
            *del2  = 2.0 * *del1*f220*g200*q22;
            *del3  = 3.0 * *del1*f330*g300*q33*aonv;
            *del1  = *del1*f311*g310*q31*aonv;
            *xlamo = fmod(mo + nodeo + argpo - theta, M_2PI);
            *xfact = mdot + xpidot - rptim + *dmdt + *domdt + *dnodt - no;
        }


        /* ------------ for sgp4, initialize the integrator ---------- */
        *xli   = *xlamo;
        *xni   = no;
        *atime = 0.0;
        *nm    = no + *dndt;

    }


//#include "debug3.cpp"


} // end dsinit




/*-----------------------------------------------------------------------------
*
* procedure dspace
*
* this procedure provides deep space contributions to mean elements for
* perturbing third body. these effects have been averaged over one
* revolution of the sun and moon. for earth resonance effects, the
* effects have been averaged over no revolutions of the satellite.
* (mean motion)
*
* author : david vallado 719-573-2600 28 jun 2005
*
* inputs :
* d2201, d2211, d3210, d3222, d4410, d4422, d5220, d5232, d5421, d5433.
* dedt,
* del1, del2, del3,
* didt,
* dmdt,
* dnodt,
* domdt,
* irez - flag for resonance 0-none, 1-one day, 2-half day
* argpo - argument of perigee
* argpdot - argument of perigee dot (rate)
* t - time
* tc -
* gsto - gst
* xfact -
* xlamo -
* no - mean motion
* atime -
* em - eccentricity
* ft -
* argpm - argument of perigee
* inclm - inclination
* xli -
* mm - mean anomaly
* xni - mean motion
* nodem - right ascension of ascending node
*
* outputs :
*     atime -
*     em - eccentricity
*     argpm - argument of perigee
*     inclm - inclination
*     xli -
*     mm - mean anomaly
*     xni -
*     nodem - right ascension of ascending node
*     dndt -
*     nm - mean motion
*
----------------------------------------------------------------------------*/

void LgmSgp_dspace( double tc, double *atime, double *em, double *argpm, 
        double *inclm, double *xli, double *mm, double *xni, double *nodem, double *dndt, double *nm, _SgpInfo *s ) {





    int     iretn, iret;
    double  delt, ft, theta, x2li, x2omi, xl, xldot, xnddt, xndt, xomi, t;
    double  g22, g32, g44, g52, g54, fasx2, fasx4, fasx6, rptim, step2, stepn, stepp;


    ft    = 0.0;
    fasx2 = 0.13130908;
    fasx4 = 2.8843198;
    fasx6 = 0.37448087;
    g22   = 5.7686396;
    g32   = 0.95240898;
    g44   = 1.8014998;
    g52   = 1.0508330;
    g54   = 4.4108898;
    rptim = 4.37526908801129966e-3; // this equates to 7.29211514668855e-5 rad/sec
    stepp = 720.0;
    stepn = -720.0;
    step2 = 259200.0;




    /* ----------- calculate deep space resonance effects ----------- */
    t = s->t;
    *dndt   = 0.0;
    theta  = fmod(s->gsto + tc * rptim, M_2PI);
    *em    += s->dedt*t;
    *inclm += s->didt*t;
    *argpm += s->domdt*t;
    *nodem += s->dnodt*t;
    *mm    += s->dmdt*t;


    // sgp4fix for negative inclinations
    // the following if statement should be commented out
    // if (*inclm < 0.0)
    // {
    // *inclm = -*inclm;
    // argpm = argpm - pi;
    // *nodem = *nodem + pi;
    // }
    /* - update resonances : numerical (euler-maclaurin) integration - */
    /* ------------------------- epoch restart ---------------------- */
    // sgp4fix for propagator problems
    // the following integration works for negative time steps and periods
    // the specific changes are unknown because the original code was so convoluted
    ft = 0.0;
    *atime = 0.0;
    if ( s->irez != 0 ) {

        if ((*atime == 0.0) || ((t >= 0.0) && (*atime < 0.0)) || ((t < 0.0) && (*atime >= 0.0))) {
            if (t >= 0.0) delt = stepp;
            else          delt = stepn;
            *atime = 0.0;
            s->xni   = s->no;
            s->xli   = s->xlamo;
        }

        iretn = 381; // added for do loop
        iret  = 0; // added for loop

        while ( iretn == 381 ) {

            if ((fabs(t) < fabs(*atime)) || (iret == 351)) {
                if (t >= 0.0) delt = stepn;
                else          delt = stepp;
                iret  = 351;
                iretn = 381;
            } else {
                if (t > 0.0) delt = stepp; // error if prev if has *atime:=0.0 and t:=0.0 (ge)
                else         delt = stepn;
                if (fabs(t - *atime) >= stepp) {
                    iret  = 0;
                    iretn = 381;
                } else {
                    ft    = t - *atime;
                    iretn = 0;
                }
            }



            /* ------------------- dot terms calculated ------------- */
            /* ----------- near - synchronous resonance terms ------- */
            if (s->irez != 2) {
                xndt  = s->del1*sin(s->xli - fasx2) + s->del2*sin(2.0*(s->xli - fasx4)) + s->del3*sin(3.0*(s->xli - fasx6));
                xldot = s->xni + s->xfact;
                xnddt = s->del1*cos(s->xli - fasx2) + 2.0*s->del2*cos(2.0*(s->xli - fasx4)) + 3.0*s->del3*cos(3.0*(s->xli - fasx6));
                xnddt = xnddt*xldot;
            } else {
                /* --------- near - half-day resonance terms -------- */
                xomi  = s->argpo + s->argpdot * *atime;
                x2omi = xomi + xomi;
                x2li  = s->xli + s->xli;
                xndt  = s->d2201*sin(x2omi + s->xli - g22) + s->d2211*sin(s->xli - g22) + s->d3210*sin(xomi + s->xli - g32) + s->d3222*sin(-xomi + s->xli - g32)
                        + s->d4410*sin(x2omi + x2li - g44)+ s->d4422*sin(x2li - g44) + s->d5220*sin(xomi + s->xli - g52) + s->d5232*sin(-xomi + s->xli - g52)
                        + s->d5421*sin(xomi + x2li - g54) + s->d5433*sin(-xomi + x2li - g54);
                xldot = s->xni + s->xfact;
                xnddt = s->d2201*cos(x2omi + s->xli - g22) + s->d2211*cos(s->xli - g22) 
                        + s->d3210*cos(xomi + s->xli - g32) + s->d3222*cos(-xomi + s->xli - g32) 
                        + s->d5220*cos(xomi + s->xli - g52) + s->d5232*cos(-xomi + s->xli - g52) 
                        + 2.0*(s->d4410*cos(x2omi + x2li - g44) 
                        + s->d4422*cos(x2li - g44) + s->d5421*cos(xomi + x2li - g54) 
                        + s->d5433*cos(-xomi + x2li - g54));
                xnddt = xnddt*xldot;
            }

            /* ----------------------- integrator ------------------- */
            if (iretn == 381) {
                s->xli   += xldot*delt + xndt*step2;
                s->xni   += xndt*delt + xnddt*step2;
                *atime += delt;
            }

        } // while iretn = 381


        *nm = s->xni + xndt*ft + xnddt*ft*ft*0.5;
        xl = s->xli + xldot*ft + xndt*ft*ft*0.5;

        if (s->irez != 1) {
            *mm = xl - 2.0*(*nodem) + 2.0*theta;
            *dndt = *nm - s->no;
        } else {
            *mm = xl - *nodem - *argpm + theta;
            *dndt = *nm - s->no;
        }
        *nm = s->no + *dndt;
    }


    //#include "debug4.cpp"


} // end dsspace



/*-----------------------------------------------------------------------------
*
* procedure initl
*
* this procedure initializes the spg4 propagator. all the initialization is
* consolidated here instead of having multiple loops inside other routines.
*
* author : david vallado 719-573-2600 28 jun 2005
*
* inputs :
*     ecco - eccentricity 0.0 - 1.0
*     epoch - epoch time in days from jan 0, 1950. 0 hr
*     inclo - inclination of satellite
*     no - mean motion of satellite
*     satn - satellite number
*
* outputs :
*     ainv - 1.0 / a
*     ao - semi major axis
*     con41 -
*     con42 - 1.0 - 5.0 cos(i)
*     cosio - cosine of inclination
*     cosio2 - cosio squared
*     eccsq - eccentricity squared
*     method - flag for deep space 'd', 'n'
*     omeosq - 1.0 - ecco * ecco
*     posq - semi-parameter squared
*     rp - radius of perigee
*     rteosq - square root of (1.0 - ecco*ecco)
*     sinio - sine of inclination
*     gsto - gst at time of observation rad
*     no - mean motion of satellite
*
*
* coupling :
* LgmSgp_GetGravConst
* gstime - find greenwich sidereal time from the julian date
*
* references :
* hoots, roehrich, norad spacetrack report #3 1980
* hoots, norad spacetrack report #6 1986
* hoots, schumacher and glover 2004
* vallado, crawford, hujsak, kelso 2006
----------------------------------------------------------------------------*/
void LgmSgp_initl( int satn, int whichconst, double ecco, double epoch, double inclo, 
                    double *no, char *method, double *ainv, double *ao, double *con41, double *con42, 
                    double *cosio, double *cosio2, double *eccsq, double *omeosq, double *posq, 
                    double *rp, double *rteosq, double *sinio, double *gsto) {

    /* --------------------- local variables ------------------------ */
    double  ak, d1, del, adel, po, x2o3, j2, xke, tumin, radiusearthkm, j3, j4, j3oj2;


    /* ----------------------- earth constants ---------------------- */
    // sgp4fix identify constants and allow alternate values
    LgmSgp_GetGravConst( whichconst, &tumin, &radiusearthkm, &xke, &j2, &j3, &j4, &j3oj2 );
    x2o3 = 2.0 / 3.0;


    /* ------------- calculate auxillary epoch quantities ---------- */
    *eccsq  = ecco*ecco;
    *omeosq = 1.0 - *eccsq;
    *rteosq = sqrt(*omeosq);
    *cosio  = cos(inclo);
    *cosio2 = (*cosio)*(*cosio);


    /* ------------------ un-kozai the mean motion ----------------- */
    ak      = pow(xke/(*no), x2o3);
    d1      = 0.75*j2*(3.0*(*cosio2) - 1.0)/((*rteosq)*(*omeosq));
    del     = d1/(ak*ak);
    adel    = ak*(1.0 - del*del - del*(1.0/3.0 + 134.0*del*del/81.0));
    del     = d1/(adel*adel);
    *no      = *no/(1.0 + del);
    *ao      = pow(xke/(*no), x2o3);
    *sinio   = sin(inclo);
    po      = *ao*(*omeosq);
    *con42   = 1.0 - 5.0*(*cosio2);
    *con41   = -(*con42)-(*cosio2)-(*cosio2);
    *ainv    = 1.0/(*ao);
    *posq    = po*po;
    *rp      = *ao*(1.0 - ecco);
    *method = 'n';
    *gsto    = LgmSgp_gstime(epoch + 2433281.5);

    //#include "debug5.cpp"


} // end LgmSgp_initl





/*-----------------------------------------------------------------------------
*
* procedure sgp4init
*
* this procedure initializes variables for sgp4.
*
* author : david vallado 719-573-2600 28 jun 2005
*
* inputs :
* satn - satellite number
* bstar - sgp4 type drag coefficient kg/m2er
* ecco - eccentricity
* epoch - epoch time in days from jan 0, 1950. 0 hr
* argpo - argument of perigee (output if ds)
* inclo - inclination
* mo - mean anomaly (output if ds)
* no - mean motion
* nodeo - right ascension of ascending node
*
* outputs :
* s->- common values for subsequent calls
* return code - non-zero on error.
* 1 - mean elements, ecc >= 1.0 or ecc < -0.001 or a < 0.95 er
* 2 - mean motion less than 0.0
* 3 - pert elements, ecc < 0.0 or ecc > 1.0
* 4 - semi-latus rectum < 0.0
* 5 - epoch elements are sub-orbital
* 6 - satellite has decayed
*
*
* coupling :
* LgmSgp_GetGravConst-
* initl -
* dscom -
* dpper -
* dsinit -
* sgp4 -
*
* references :
* hoots, roehrich, norad spacetrack report #3 1980
* hoots, norad spacetrack report #6 1986
* hoots, schumacher and glover 2004
* vallado, crawford, hujsak, kelso 2006
----------------------------------------------------------------------------*/


int LgmSgp_SGP4_Init( _SgpInfo *s, _SgpTLE *t ) {


    /* --------------------- local variables ------------------------ */
    double  ao, ainv, con42, cosio, sinio, cosio2, eccsq, omeosq, posq, rp, rteosq, cnodm;
    double  snodm, cosim, sinim, cosomm, sinomm, cc1sq, cc2, cc3, coef, coef1, cosio4, day;
    double  dndt, em, emsq, eeta, etasq, gam, argpm, nodem, inclm, mm, nm, perige, pinvsq;
    double  psisq, qzms24, rtemsq, s1, s2, s3, s4, s5, s6, s7, sfour, ss1, ss2, ss3, ss4;
    double  ss5, ss6, ss7, sz1, sz2, sz3, sz11, sz12, sz13, sz21, sz22, sz23, sz31, sz32;
    double  sz33, tc, temp, temp1, temp2, temp3, tsi, xpidot, xhdot1, z1, z2, z3, z11, z12;
    double  z13, z21, z22, z23, z31, z32, z33, qzms2t, ss, j2, j3oj2, j4, x2o3;//, r[3], v[3];
    double tumin, radiusearthkm, xke, j3, tmp, tmp2, tmp4;

    int     whichconst, satn;
    double  epoch, xbstar, xecco, xargpo, xinclo, xmo, xno, xnodeo;


    // mgh - transfer vars here.
    whichconst = s->GravConst;
whichconst = SGP_wgs72;
    satn   = t->IdNumber;
//    epoch  = t->ElementSetEpoch;
epoch  = t->JD-2433281.5;
    xbstar = t->BstarDrag;
    xecco  = t->Eccentricity;

    xargpo = t->ArgOfPerigee;
    xinclo = t->Inclination;
    xmo    = t->MeanAnomaly;
    xno    = t->MeanMotion;
    xnodeo = t->RAofAscNode;
   

    
double xpdotp = 1440.0/M_2PI;
xnodeo *= M_PI/180.0;
xargpo *= M_PI/180.0;
xmo    *= M_PI/180.0;
xinclo *= M_PI/180.0;
xno     = xno/xpdotp;
// WHERE ARE THESE SET IN THE NEW STUFF?
//    s->XNDT2O  =  s->XNDT2O*tmp;
//    s->XNDD6O  =  s->XNDD6O*tmp/SGP_XMNPDA;

    



    /* ------------------------ initialization --------------------- */
    // sgp4fix divisor for divide by zero check on inclination
    const double temp4 = 1.0 + cos(M_PI-1.0e-9);


    /* ----------- set all near earth variables to zero ------------ */
    s->isimp = 0;       s->method = 'n';    s->aycof = 0.0;
    s->con41 = 0.0;     s->cc1 = 0.0;       s->cc4 = 0.0;
    s->cc5 = 0.0;       s->d2 = 0.0;        s->d3 = 0.0;
    s->d4 = 0.0;        s->delmo = 0.0;     s->eta = 0.0;
    s->argpdot = 0.0;   s->omgcof = 0.0;    s->sinmao = 0.0;
    s->t = 0.0;         s->t2cof = 0.0;     s->t3cof = 0.0;
    s->t4cof = 0.0;     s->t5cof = 0.0;     s->x1mth2 = 0.0;
    s->x7thm1 = 0.0;    s->mdot = 0.0;      s->nodedot = 0.0;
    s->xlcof = 0.0;     s->xmcof = 0.0;     s->nodecf = 0.0;


    /* ----------- set all deep space variables to zero ------------ */
    s->irez = 0;        s->d2201 = 0.0;     s->d2211 = 0.0;
    s->d3210 = 0.0;     s->d3222 = 0.0;     s->d4410 = 0.0;
    s->d4422 = 0.0;     s->d5220 = 0.0;     s->d5232 = 0.0;
    s->d5421 = 0.0;     s->d5433 = 0.0;     s->dedt = 0.0;
    s->del1 = 0.0;      s->del2 = 0.0;      s->del3 = 0.0;
    s->didt = 0.0;      s->dmdt = 0.0;      s->dnodt = 0.0;
    s->domdt = 0.0;     s->e3 = 0.0;        s->ee2 = 0.0;
    s->peo = 0.0;       s->pgho = 0.0;      s->pho = 0.0;
    s->pinco = 0.0;     s->plo = 0.0;       s->se2 = 0.0;
    s->se3 = 0.0;       s->sgh2 = 0.0;      s->sgh3 = 0.0;
    s->sgh4 = 0.0;      s->sh2 = 0.0;       s->sh3 = 0.0;
    s->si2 = 0.0;       s->si3 = 0.0;       s->sl2 = 0.0;
    s->sl3 = 0.0;       s->sl4 = 0.0;       s->gsto = 0.0;
    s->xfact = 0.0;     s->xgh2 = 0.0;      s->xgh3 = 0.0;
    s->xgh4 = 0.0;      s->xh2 = 0.0;       s->xh3 = 0.0;
    s->xi2 = 0.0;       s->xi3 = 0.0;       s->xl2 = 0.0;
    s->xl3 = 0.0;       s->xl4 = 0.0;       s->xlamo = 0.0;
    s->zmol = 0.0;      s->zmos = 0.0;      s->atime = 0.0;
    s->xli = 0.0;       s->xni = 0.0;


    // sgp4fix - note the following variables are also passed directly via s->
    // it is possible to streamline the sgp4init call by deleting the "x"
    // variables, but the user would need to set the s->* values first. we
    // include the additional assignments in case twoline2rv is not used.
    s->bstar = xbstar;
    s->ecco  = xecco;
    s->argpo = xargpo;
    s->inclo = xinclo;
    s->mo    = xmo;
    s->no    = xno;
    s->nodeo = xnodeo;

    /* ------------------------ earth constants ----------------------- */
    // sgp4fix identify constants and allow alternate values
    LgmSgp_GetGravConst( whichconst, &tumin, &radiusearthkm, &xke, &j2, &j3, &j4, &j3oj2 );
    ss     = 78.0 / radiusearthkm + 1.0;
    qzms2t = pow(((120.0 - 78.0) / radiusearthkm), 4);
    x2o3   = 2.0 / 3.0;
    s->init = 'y';
    s->t = 0.0;


    LgmSgp_initl( satn, whichconst, s->ecco, epoch, s->inclo, 
            &(s->no), &(s->method), &ainv, &ao, &(s->con41), &con42, 
            &cosio, &cosio2, &eccsq, &omeosq, &posq, 
            &rp, &rteosq, &sinio, &(s->gsto));
    s->error = 0;
    if (rp < 1.0) {
        // printf("# *** satn%d epoch elts sub-orbital ***\n", satn);
        s->error = 5;
    }



    if ((omeosq >= 0.0 ) || ( s->no >= 0.0)) {

        s->isimp = (rp < (220.0/radiusearthkm + 1.0)) ? 1 : 0;
        sfour  = ss;
        qzms24 = qzms2t;
        perige = (rp - 1.0) * radiusearthkm;

        /* - for perigees below 156 km, s and qoms2t are altered - */
        if (perige < 156.0) {
            sfour  = perige - 78.0;
            if (perige < 98.0) sfour = 20.0;
            qzms24 = pow(((120.0 - sfour) / radiusearthkm), 4.0);
            sfour  = sfour / radiusearthkm + 1.0;
        }

        pinvsq = 1.0 / posq;
        tsi    = 1.0 / (ao - sfour);
        s->eta = ao * s->ecco * tsi;
        etasq = s->eta * s->eta;
        eeta  = s->ecco * s->eta;
        psisq = fabs(1.0 - etasq);
        tmp2  = tsi*tsi; tmp4=tmp2*tmp2;
        coef  = qzms24 * tmp4;
        coef1 = coef / pow(psisq, 3.5);
        cc2   = coef1*s->no*(ao*(1.0 + 1.5*etasq + eeta*(4.0 + etasq)) + 0.375*j2
                *tsi/psisq*s->con41*(8.0 + 3.0*etasq*(8.0 + etasq)));
        s->cc1 = s->bstar * cc2;
        cc3 = 0.0;
        if (s->ecco > 1.0e-4) cc3 = -2.0*coef*tsi*j3oj2*s->no*sinio/s->ecco;
        s->x1mth2 = 1.0 - cosio2;
        s->cc4 = 2.0*s->no*coef1*ao*omeosq
                        *(s->eta*(2.0 + 0.5*etasq) + s->ecco*(0.5 + 2.0*etasq) 
                            - j2*tsi/(ao*psisq)*(-3.0*s->con41*(1.0 - 2.0*eeta + etasq*(1.5 - 0.5*eeta)) 
                            + 0.75*s->x1mth2*(2.0*etasq - eeta*(1.0 + etasq))*cos(2.0*s->argpo)));

        s->cc5 = 2.0*coef1*ao*omeosq*(1.0 + 2.75*(etasq + eeta) + eeta*etasq);
        cosio4 = cosio2 * cosio2;
        temp1  = 1.5 * j2 * pinvsq * s->no;
        temp2  = 0.5 * temp1 * j2 * pinvsq;
        temp3  = -0.46875 * j4 * pinvsq * pinvsq * s->no;
        s->mdot    = s->no + 0.5*temp1*rteosq*s->con41 + 0.0625*temp2*rteosq*(13.0 - 78.0*cosio2 + 137.0*cosio4);
        s->argpdot = -0.5*temp1*con42 + 0.0625*temp2*(7.0 - 114.0*cosio2 + 395.0*cosio4) + temp3*(3.0 - 36.0*cosio2 + 49.0*cosio4);
        xhdot1     = -temp1*cosio;
        s->nodedot = xhdot1 + (0.5*temp2*(4.0 - 19.0*cosio2) + 2.0*temp3*(3.0 - 7.0*cosio2))*cosio;
        xpidot     = s->argpdot+ s->nodedot;
        s->omgcof  = s->bstar * cc3 * cos(s->argpo);
        s->xmcof   = 0.0;
        if (s->ecco > 1.0e-4) s->xmcof = -x2o3*coef*s->bstar/eeta;
        s->nodecf  = 3.5*omeosq*xhdot1*s->cc1;
        s->t2cof   = 1.5*s->cc1;

        // sgp4fix for divide by zero with xinco = 180 deg
        if (fabs(cosio+1.0) > 1.5e-12) s->xlcof = -0.25*j3oj2*sinio*(3.0 + 5.0*cosio)/(1.0 + cosio);
        else                           s->xlcof = -0.25*j3oj2*sinio*(3.0 + 5.0*cosio)/temp4;

        s->aycof  = -0.5*j3oj2*sinio;
        tmp       = 1.0 + s->eta*cos(s->mo);
        s->delmo  = tmp*tmp*tmp;
        s->sinmao = sin(s->mo);
        s->x7thm1 = 7.0*cosio2 - 1.0;

        /* --------------- deep space initialization ------------- */
        if ((2.0*M_PI / s->no) >= 225.0) {
            s->method = 'd';
            s->isimp  = 1;
            tc        = 0.0;
            inclm     = s->inclo;
            // mgh - this is the stupidest coding Ive seen in a long time...
            LgmSgp_dscom( epoch, s->ecco, s->argpo, tc, s->inclo, s->nodeo, s->no, 
                            &snodm, &cnodm, &sinim, &cosim, &sinomm, &cosomm, &day, &(s->e3), &(s->ee2), &em, &emsq, &gam, &(s->peo), 
                            &(s->pgho), &(s->pho), &(s->pinco), &(s->plo), &rtemsq, &(s->se2), &(s->se3), &(s->sgh2), &(s->sgh3), 
                            &(s->sgh4), &(s->sh2), &(s->sh3), &(s->si2), &(s->si3), &(s->sl2), &(s->sl3), &(s->sl4), &s1, &s2, &s3, 
                            &s4, &s5, &s6, &s7, &ss1, &ss2, &ss3, &ss4, &ss5, &ss6, &ss7, &sz1, &sz2, &sz3, &sz11, &sz12, &sz13, 
                            &sz21, &sz22, &sz23, &sz31, &sz32, &sz33, &(s->xgh2), &(s->xgh3), &(s->xgh4), &(s->xh2), &(s->xh3), 
                            &(s->xi2), &(s->xi3), &(s->xl2), &(s->xl3), &(s->xl4), &nm, &z1, &z2, &z3, &z11, &z12, &z13, &z21, 
                            &z22, &z23, &z31, &z32, &z33, &(s->zmol), &(s->zmos) );

            LgmSgp_dpper( inclm, s->init, &(s->ecco), &(s->inclo), &(s->nodeo), &(s->argpo), &(s->mo), s );

            argpm = 0.0;
            nodem = 0.0;
            mm    = 0.0;
            LgmSgp_dsinit( whichconst, cosim, emsq, s->argpo, s1, s2, s3, s4, s5, sinim, ss1, ss2, ss3, ss4, ss5, sz1, sz3, 
                            sz11, sz13, sz21, sz23, sz31, sz33, s->t, tc, s->gsto, s->mo, s->mdot, s->no, s->nodeo, s->nodedot, 
                            xpidot, z1, z3, z11, z13, z21, z23, z31, z33, s->ecco, eccsq, 
                            &em, &argpm, &inclm, &mm, &nm, &nodem, &s->irez, &s->atime, &s->d2201, &s->d2211, &s->d3210, &s->d3222, &s->d4410, 
                            &s->d4422, &s->d5220, &s->d5232, &s->d5421, &s->d5433, &s->dedt, &s->didt, &s->dmdt, &dndt, &s->dnodt, &s->domdt, 
                            &s->del1, &s->del2, &s->del3, &s->xfact, &s->xlamo, &s->xli, &s->xni);
        }




        /* ----------- set variables if not deep space ----------- */
        if ( s->isimp != 1 ) {
            cc1sq    = s->cc1*s->cc1;
            s->d2    = 4.0*ao*tsi*cc1sq;
            temp     = s->d2*tsi*s->cc1/3.0;
            s->d3    = (17.0*ao + sfour)*temp;
            s->d4    = 0.5*temp*ao*tsi*(221.0*ao + 31.0*sfour)*s->cc1;
            s->t3cof = s->d2 + 2.0*cc1sq;
            s->t4cof = 0.25*(3.0*s->d3 + s->cc1*(12.0*s->d2 + 10.0*cc1sq));
            s->t5cof = 0.2*(3.0*s->d4 + 12.0*s->cc1*s->d3 + 6.0*s->d2*s->d2 + 15.0*cc1sq*(2.0*s->d2 + cc1sq));
        }


    } // if omeosq = 0 ...





    /* finally propogate to zero epoch to initialise all others. */
//FIX ME
s->GravConst = whichconst;
if (s->error == 0) LgmSgp_SGP4( 0.0, s );
    s->init = 'n';

    //#include "debug6.cpp"
    return( s->error );



} // end sgp4init



/* -----------------------------------------------------------------------------
*
* function LgmSgp_GetGravConst
*
* this function gets constants for the propagator. note that mu is identified to
* facilitiate comparisons with newer models.
*
* author : david vallado 719-573-2600 21 jul 2006
*
* mods: Mike Henderson 20, feb 2009 -- changed arg passing
*
* inputs :
* whichconst - which set of constants to use 72, 84
*
* outputs :
* tumin - minutes in one time unit
* radiusearthkm - radius of the earth in km
* xke - reciprocal of tumin
* j2, j3, j4 - un-normalized zonal harmonic values
* j3oj2 - j3 divided by j2
*
* locals :
* mu - earth gravitational parameter
*
* coupling :
* none
*
* references :
* norad spacetrack report #3
* vallado, crawford, hujsak, kelso 2006
--------------------------------------------------------------------------- */
void LgmSgp_GetGravConst( int whichconst, double *tumin, double *radiusearthkm,
                    double *xke, double *j2, double *j3, double *j4, double *j3oj2 ) {

    double mu;

    switch ( whichconst ) {


            // -- wgs-72 low precision str#3 constants --
            case SGP_wgs72old:
                *radiusearthkm = 6378.135; // km
                *xke = 0.0743669161;
                *tumin = 1.0 / *xke;
                *j2 = 0.001082616;
                *j3 = -0.00000253881;
                *j4 = -0.00000165597;
                *j3oj2 = *j3 / *j2;
                break;
            // ------------ wgs-72 constants ------------
            case SGP_wgs72:
                mu = 398600.8; // in km3 / s2
                *radiusearthkm = 6378.135; // km
                *xke = 60.0 / sqrt( (*radiusearthkm) * (*radiusearthkm) * (*radiusearthkm)/mu);
                *tumin = 1.0 / *xke;
                *j2 = 0.001082616;
                *j3 = -0.00000253881;
                *j4 = -0.00000165597;
                *j3oj2 = *j3 / *j2;
                break;
            case SGP_wgs84:
                // ------------ wgs-84 constants ------------
                mu = 398600.5; // in km3 / s2
                *radiusearthkm = 6378.137; // km
                *xke = 60.0 / sqrt( (*radiusearthkm) * (*radiusearthkm) * (*radiusearthkm)/mu);
                *tumin = 1.0 / *xke;
                *j2 = 0.00108262998905;
                *j3 = -0.00000253215306;
                *j4 = -0.00000161098761;
                *j3oj2 = *j3 / *j2;
                break;
            default:
                fprintf(stderr,"unknown gravity option (%d)\n",whichconst);
                break;
    }

} // end getgravcons




/*-----------------------------------------------------------------------------
*
* procedure sgp4
*
* this procedure is the sgp4 prediction model from space command. this is an
* updated and combined version of sgp4 and sdp4, which were originally
* published separately in spacetrack report #3. this version follows the nasa
* release on the internet. there are a few fixes that are added to correct
* known errors in the existing implementations.
*
* author :  David Vallado (719-573-2600)            28 jun 2005
*               Original version in "Spacetrack Report #3 Revisited"
*
* Modifications:
*           Mike Henderson (mghenderson@lanl.gov)   20 Feb 2009
*
*             - Changed manner in which structure is passed.
*               Despite using a structure, the original version had an insane 
*               number of function args. Also strange way of passing...
*
*             - Changed an instance of pow(XXX, 3) to simple multiplies.
*
*             - Integrated code with my own TLE reader.
*
*             - Integrated code with my own TLE reader.
*
*           
*
* inputs :
*   s      - initialised structure from sgp4init() call.
*   tsince - time eince epoch (minutes)
*
* outputs (via structure):
*   r - position vector km
*   v - velocity km/sec
*
* return code - non-zero on error.
*   1 - mean elements, ecc >= 1.0 or ecc < -0.001 or a < 0.95 er
*   2 - mean motion less than 0.0
*   3 - pert elements, ecc < 0.0 or ecc > 1.0
*   4 - semi-latus rectum < 0.0
*   5 - epoch elements are sub-orbital
*   6 - satellite has decayed
*
*
* references :
*     hoots, roehrich, norad spacetrack report #3 1980
*     hoots, norad spacetrack report #6 1986
*     hoots, schumacher and glover 2004
*     vallado, crawford, hujsak, kelso 2006 (This is the Revisiting Spacetrack Report #3)
----------------------------------------------------------------------------*/
//int sgp4( int whichconst, elsetrec& s-> double tsince, double r[3], double v[3]) {
int LgmSgp_SGP4( double tsince, _SgpInfo *s ) {

    double  am, argpdf, argpm, argpp, axnl, aynl, betal, cnod, cos2u, coseo1, cosi, cosim;
    double  cosip, cosisq, cossu, cosu, delm, delomg, dndt, ecose, el2, em, emsq, eo1, ep;
    double  esine, inclm, j2, j3, j3oj2, j4, mm, mp, mrt=0.0, mvt, nm, nodedf, nodem;
    double  nodep, pl, radiusearthkm, rdotl, rl, rvdot, rvdotl, sin2u, sineo1, sini, sinim;
    double  sinip, sinsu, sinu, snod, su, t2, t3, t4, tc, tem5, temp, temp1, temp2, tempa;
    double  tempe, templ, tumin, u, ux, uy, uz, vkmpersec, vx, vy, vz, x2o3, xinc;
    double  xincp, xke, xl, xlm, xmdf, xmx, xmy, xnode, tmp;
    int     ktr;


    /* ------------------ set mathematical constants --------------- */
    // sgp4fix divisor for divide by zero check on inclination
    const double temp4 = 1.0 + cos(M_PI-1.0e-9);
    x2o3 = 2.0 / 3.0;
    // sgp4fix identify constants and allow alternate values
    LgmSgp_GetGravConst( s->GravConst, &tumin, &radiusearthkm, &xke, &j2, &j3, &j4, &j3oj2 );
    vkmpersec = radiusearthkm * xke/60.0;



    /* --------------------- clear sgp4 error flag ----------------- */
    s->t     = tsince;
    s->error = 0;



    /* ------- update for secular gravity and atmospheric drag ----- */
    xmdf   = s->mo + s->mdot*s->t;
    argpdf = s->argpo + s->argpdot*s->t;
    nodedf = s->nodeo + s->nodedot*s->t;
    argpm  = argpdf;
    mm     = xmdf;
    t2     = s->t*s->t;
    nodem  = nodedf + s->nodecf*t2;
    tempa  = 1.0 - s->cc1*s->t;
    tempe  = s->bstar*s->cc4*s->t;
    templ  = s->t2cof*t2;

    if ( s->isimp != 1 ) {
        delomg = s->omgcof*s->t;
        tmp    = 1.0 + s->eta*cos(xmdf);
        delm   = s->xmcof*(tmp*tmp*tmp - s->delmo);
        temp   = delomg + delm;
        mm     = xmdf + temp;
        argpm  = argpdf - temp;
        t3     = t2*s->t;
        t4     = t3*s->t;
        tempa  = tempa - s->d2*t2 - s->d3*t3 - s->d4*t4;
        tempe  = tempe + s->bstar*s->cc5*(sin(mm) - s->sinmao);
        templ  = templ + s->t3cof*t3 + t4*(s->t4cof + s->t*s->t5cof);
    }

    nm    = s->no;
    em    = s->ecco;
    inclm = s->inclo;

    if ( s->method == 'd' ) {
        tc = s->t;
        //LgmSgp_dspace( tc, &(s->atime), &em, &argpm, &inclm, &mm, &nodem, &dndt, &nm, s);
        LgmSgp_dspace( tc, &(s->atime), &em, &argpm, &inclm, &(s->xli), &mm, &(s->xni), &nodem, &dndt, &nm, s ); 

    } // if method = d

    if (nm <= 0.0) {
        // printf("# error nm %f\n", nm);
        s->error = 2;
    }

    am = pow((xke/nm), x2o3)*tempa*tempa;
    nm = xke/pow(am, 1.5);
    em = em - tempe;
    // fix tolerance for error recognition
    if ((em >= 1.0) || (em < -0.001) || (am < 0.95)) {
        // printf("# error em %f\n", em);
        s->error = 1;
    }

    if (em < 0.0) em = 1.0e-6;
    mm   += s->no * templ;
    xlm   = mm + argpm + nodem;
    emsq  = em * em;
    temp  = 1.0 - emsq;
    nodem = fmod(nodem, M_2PI);
    argpm = fmod(argpm, M_2PI);
    xlm   = fmod(xlm, M_2PI);
    mm    = fmod(xlm - argpm - nodem, M_2PI);


    /* ----------------- compute extra mean quantities ------------- */
    sinim = sin(inclm);
    cosim = cos(inclm);



    /* -------------------- add lunar-solar periodics -------------- */
    ep    = em;
    xincp = inclm;
    argpp = argpm;
    nodep = nodem;
    mp    = mm;
    sinip = sinim;
    cosip = cosim;

    if ( s->method == 'd' ) {

        LgmSgp_dpper( s->inclo, 'n', &ep, &xincp, &nodep, &argpp, &mp, s);

        if (xincp < 0.0) {
            xincp = -xincp;
            nodep = nodep + M_PI;
            argpp = argpp - M_PI;
        }
        if ( (ep < 0.0) || (ep > 1.0) ) {
            // printf("# error ep %f\n", ep);
            s->error = 3;
        }

    } // if method = d

    /* -------------------- long period periodics ------------------ */
    if ( s->method == 'd' ) {
        sinip = sin(xincp);
        cosip = cos(xincp);
        s->aycof = -0.5*j3oj2*sinip;
        // sgp4fix for divide by zero for xincp = 180 deg
        if (fabs(cosip+1.0) > 1.5e-12)  s->xlcof = -0.25*j3oj2*sinip*(3.0 + 5.0*cosip)/(1.0 + cosip);
        else                            s->xlcof = -0.25*j3oj2*sinip*(3.0 + 5.0*cosip)/temp4;
    }
    axnl = ep*cos(argpp);
    temp = 1.0/(am*(1.0 - ep*ep));
    aynl = ep*sin(argpp) + temp*s->aycof;
    xl   = mp + argpp + nodep + temp*s->xlcof*axnl;




    /* --------------------- solve kepler's equation --------------- */
    u    = fmod(xl - nodep, M_2PI);
    eo1  = u;
    tem5 = 9999.9;
    ktr  = 1;
    // sgp4fix for kepler iteration
    // the following iteration needs better limits on corrections
    while ( (fabs(tem5) >= 1.0e-12) && (ktr <= 10) ) {
        sineo1 = sin(eo1);
        coseo1 = cos(eo1);
        tem5   = 1.0 - coseo1*axnl - sineo1*aynl;
        tem5   = (u - aynl*coseo1 + axnl*sineo1 - eo1)/tem5;
        if(fabs(tem5) >= 0.95) tem5 = (tem5 > 0.0) ? 0.95 : -0.95;
        eo1 += tem5;
        ++ktr;
    }


    /* ------------- short period preliminary quantities ----------- */
    ecose = axnl*coseo1 + aynl*sineo1;
    esine = axnl*sineo1 - aynl*coseo1;
    el2   = axnl*axnl + aynl*aynl;
    pl    = am*(1.0-el2);

    if (pl < 0.0) {

        // printf("# error pl %f\n", pl);
        s->error = 4;

    } else {

        rl     = am*(1.0 - ecose);
        rdotl  = sqrt(am)*esine/rl;
        rvdotl = sqrt(pl)/rl;
        betal  = sqrt(1.0 - el2);
        temp   = esine/(1.0 + betal);
        sinu   = am/rl*(sineo1 - aynl - axnl*temp);
        cosu   = am/rl*(coseo1 - axnl + aynl*temp);
        su     = atan2(sinu, cosu);
        sin2u  = (cosu + cosu)*sinu;
        cos2u  = 1.0 - 2.0*sinu*sinu;
        temp   = 1.0/pl;
        temp1  = 0.5*j2*temp;
        temp2  = temp1*temp;



        /* -------------- update for short period periodics ------------ */
        if ( s->method == 'd' ) {
            cosisq    = cosip*cosip;
            s->con41  = 3.0*cosisq-1.0;
            s->x1mth2 = 1.0-cosisq;
            s->x7thm1 = 7.0*cosisq-1.0;
        }
        mrt   = rl*(1.0 - 1.5*temp2*betal* s->con41) + 0.5*temp1*s->x1mth2*cos2u;
        su    = su - 0.25*temp2*s->x7thm1*sin2u;
        xnode = nodep + 1.5*temp2*cosip*sin2u;
        xinc  = xincp + 1.5*temp2*cosip*sinip*cos2u;
        mvt   = rdotl - nm*temp1*s->x1mth2*sin2u/xke;
        rvdot = rvdotl + nm*temp1*(s->x1mth2*cos2u + 1.5*s->con41)/xke;

        /* --------------------- orientation vectors ------------------- */
        sinsu = sin(su);
        cossu = cos(su);

        snod  = sin(xnode);
        cnod  = cos(xnode);

        sini  = sin(xinc);
        cosi  = cos(xinc);

        xmx = -snod * cosi;
        xmy = cnod * cosi;

        ux = xmx * sinsu + cnod * cossu;
        uy = xmy * sinsu + snod * cossu;
        uz = sini * sinsu;

        vx = xmx * cossu - cnod * sinsu;
        vy = xmy * cossu - snod * sinsu;
        vz = sini * cossu;


        /* --------- position and velocity (in km and km/sec) ---------- */
        s->X  = (mrt * ux)* radiusearthkm;
        s->Y  = (mrt * uy)* radiusearthkm;
        s->Z  = (mrt * uz)* radiusearthkm;

        s->VX = (mvt * ux + rvdot * vx) * vkmpersec;
        s->VY = (mvt * uy + rvdot * vy) * vkmpersec;
        s->VZ = (mvt * uz + rvdot * vz) * vkmpersec;

    } // if pl > 0

    // sgp4fix for decaying satellites
    if (mrt < 1.0) {
        // printf("# decay condition %11.6f \n",mrt);
        s->error = 6;
    }

    //#include "debug7.cpp"
    return( s->error );


} // end sgp4



/*
 *   $Id$
 */

