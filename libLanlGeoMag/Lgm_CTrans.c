#include "Lgm/Lgm_CTrans.h"
#include "Lgm/Lgm_Quat.h"
#include "config.h"

#include <ctype.h>
#include <time.h>

#define USE_HIGH_ACCURACY_SUN 1

#ifndef LGM_EOP_DATA_DIR
#warning "hard-coding LGM_EOP_DATA_DIR because it was not in config.h"
#define LGM_EOP_DATA_DIR    /usr/local/share/LanlGeoMag/EopData
#endif



void Lgm_free_ctrans( Lgm_CTrans *c ) {



    free( c->l.LeapSecondDates );
    free( c->l.LeapSecondJDs );
    free( c->l.LeapSeconds );

    free(c);

}


void Lgm_ctransDefaults(Lgm_CTrans *c, int Verbose) {
    /*
     * Load Leap Seconds Stuff
     */
    Lgm_LoadLeapSeconds( c );

    c->Verbose     = Verbose;

    /*
     * Set Number of Nutation series terms to 106
     */
    c->nNutationTerms = 106;

    /*
     *  Initialize Earth Orientation Parameters to defaults They can be
     *  over-ridden later with actual values.  This allows users to use the
     *  less accurate defaults if EOP data is unavailable or not
     *  needed/desired.
     */
    c->DUT1  = 0.0; // seconds
    c->xp    = 0.0; // radians
    c->yp    = 0.0; // radians
    c->ddPsi = 0.0; // radians
    c->ddEps = 0.0; // radians

}


Lgm_CTrans *Lgm_init_ctrans( int Verbose ) {

    Lgm_CTrans  	*c;

    c = (Lgm_CTrans *) calloc (1, sizeof(*c));
    Lgm_ctransDefaults(c, Verbose);

    return c;

}

/**
 *  The Lgm_CTrans structure has pointers in it, so simple
 *  asignments (e.g. *t = *s) are dangerous. Here we make sure that
 *  the target gets an independent copy of the structure.
 */
Lgm_CTrans *Lgm_CopyCTrans( Lgm_CTrans *s ) {

    Lgm_CTrans *t;
    int         n;


    if ( s == NULL) {
        printf("Lgm_CopyCTrans: Error, source structure is NULL\n");
        return;
    }

    /*
     * Allocate memory for target
     */
    t = (Lgm_CTrans *)calloc( 1, sizeof(Lgm_CTrans) );

    /*
     * Do memcpy. Note that for things that are dynamically allocated in
     * Lgm_CTrans structure, this will copy pointers to to memory that belong
     * to the source. Currently, the LeapSecond stuff is dynamically allocated,
     * so we need to allocate our own memory for those and then do a memcpy on
     * the contents.
     */
    memcpy( t, s, sizeof(Lgm_CTrans) );

    /*
     *  Now, copy the LeapSeconds stuff properly (they were dyn. allocated).
     *  (memcpy's args are (dest, src, size))
     */
    n = t->l.nLeapSecondDates = s->l.nLeapSecondDates;
    t->l.LeapSecondDates  = (long int *) malloc( n*sizeof(long int) );
    t->l.LeapSecondJDs    = (double *)   malloc( n*sizeof(double) );
    t->l.LeapSeconds      = (double *)   malloc( n*sizeof(double) );
    memcpy( t->l.LeapSecondDates, s->l.LeapSecondDates, n*sizeof(long int) );
    memcpy( t->l.LeapSecondJDs,   s->l.LeapSecondJDs, n*sizeof(double) );
    memcpy( t->l.LeapSeconds,     s->l.LeapSeconds, n*sizeof(double) );

    return( t );

}



/**
 *  Converts RA and DEC to unit vector in cartesian coords
 */
void Lgm_Radec_to_Cart(double ra, double dec, Lgm_Vector *r) {

    double  CosDec;


    /*
     *  Convert ra/dec from degrees to radians
     */
    ra  *= RadPerDeg;
    dec *= RadPerDeg;

    CosDec = cos(dec);

    /*
     *  Compute cartesian coordinates (in GEI)
     */
    r->x = CosDec * cos(ra);
    r->y = CosDec * sin(ra);
    r->z = sin(dec);

}

double Lgm_angle2pi(double angle) {

    int n;
    double a = M_2PI;

    if (angle < 0.0){
	    n = (int)(angle/a) - 1;
	    return(angle-n*a);
    } else if (angle > a){
	    n = (int)(angle/a);
	    return(angle-n*a);
    } else{
	    return(angle);
    }

}

double Lgm_angle360(double angle) {

    int n;

    if (angle < 0.0){
	    n = (int)(angle/360.0) - 1;
	    return(angle-n*360.0);
    } else if (angle > 360.0){
	    n = (int)(angle/360.0);
	    return(angle-n*360.0);
    } else{
	    return(angle);
    }

}



/**
 *  Compute the Julian Day number for the given date.
 *  Julian Date is the number of days since noon of Jan 1 4713 B.C.
 */
long int Lgm_JD_to_Date(double jd, int *ny, int *nm, int *nd, double *UT) {

    long int    I, A, B, C, D, E, G, Date;
    double      F;

    jd += 0.5;

    I = (int)jd;
    F = jd - (double)I;

    if ( I > 2299160 ) {
        A = (int)( ((double)I - 1867216.25)/36524.25 );
        B = I+1+A-(int)( (double)A/4.0 );
    } else {
        B = I;
    }

    C = B+1524;

    D = (int)( ((double)C - 122.1)/365.25 );
    E = (int)( 365.25*D );

    G = (int)( ((double)C-(double)E)/30.6001 );

    /*
     *  Day of month
     */
    *nd = C - E - (int)( 30.6001*G );
    *UT = 24.0*F;

    if ( G <= 13 ) {
        *nm = G-1;
    } else if ( G > 13 ) {
        *nm = G-13;
    }

    if ( *nm > 2 ) {
        *ny = D-4716;
    } else if ( *nm <= 2 ) {
        *ny = D-4715;
    }

    Date = *ny*10000 + *nm*100 + *nd;

    return( Date );

}

long int Lgm_MJD_to_Date(double mjd, int *ny, int *nm, int *nd, double *UT) {
    return( Lgm_JD_to_Date( mjd+2400000.5, ny, nm, nd, UT ) );
}

double Lgm_hour24(double hour) {
    int n;

    if (hour < 0.0){
	    n = (int)(hour/24.0) - 1;
	    return(hour-n*24.0);
    } else if (hour > 24.0){
	    n = (int)(hour/24.0);
	    return(hour-n*24.0);
    } else{
	    return(hour);
    }
}

double Lgm_kepler(double M, double e) {
    int 	n=0;
    double 	E, Eold, eps = 1.0e-8;

    E = M + e*sin(M);
    do{
	    Eold = E;
	    E = Eold + (M-Eold+e*sin(Eold)) / (1.0-e*cos(Eold));
	    ++n;
    } while( (fabs(E-Eold) > eps) && (n < 100) );
    return(E);
}



/*
 *  This routine takes date and time and sets up all the transformation
 *  matrices to do coord transformations.
 *
 *   	Date    -- (long int) 	containing date ( e.g. YYYYMMDD or YYYYDDD )
 *	    UT	    -- (double)	    decimal UT in hours ( e.g. 12.0 is noon ).
 *      c	    -- (Lgm_CTrans)	structure containing releveant trans. data
 */
void Lgm_Set_Coord_Transforms(long int date, double UTC, Lgm_CTrans *c) {

    int    	    year, month, day, doy;
    double 	    TU, gmst, gast, sn, cs, Time;
    double 	    varep, varpi, spsi;
    double 	    eccen, epsilon;
    double 	    days, M, E, nu, lambnew;
    double 	    r0, earth_sun_distance;
    double 	    RA, DEC, RA_Moon, DEC_Moon;
    double 	    gclat, glon, psi;
    double 	    Lmoon_0, P0, N0, Imoon, Lmoon, Mmoon_m, Nmoon, Cmoon;
    double 	    Emoon_nu, Amoon_e, Amoon_3, Mmoon_mp, Emoon_c, amoon, emoon;
    double 	    Amoon_4, Lmoon_p, Vmoon, Lmoon_pp, Nmoon_p, LambdaMoon, BetaMoon;
    double	    g[13][13], h[13][13], Tmp[3][3];
    double  	l, r, b;
    double  	varep90, varpi90;
    double  	Zeta, Theta, Zee;
    double  	SinZee, CosZee, SinZeta, CosZeta, SinTheta, CosTheta;
    Lgm_Vector 	S, K, Y, Z, D, Dmod, Dgsm, u_mod, u_tod, u_gei, Zgeo, X;
    double      RA_tod, DEC_tod;
    double      RA_gei, DEC_gei;
    double	    sxp, cxp, syp, cyp;
    double      sdp, cdp, se, ce, set, cet;
    double      T_UT1, T2_UT1, T3_UT1, T_TT, T2_TT, T3_TT, T4_TT, T5_TT;
    double      cos_epsilon, sin_epsilon, sin_lambnew, sin_b, cos_b, sin_l, tmp;
    int 	    i, j, N;


    /*
     * Set UTC values
     */
    Lgm_Doy( date, &year, &month, &day, &doy);
    date          = year*10000+month*100+day;
    c->UTC.Date   = date;
    c->UTC.Time   = UTC;
    c->UTC.Year   = year;
    c->UTC.Month  = month;
    c->UTC.Day    = day;
    c->UTC.Doy    = doy;
    c->UTC.Dow    = Lgm_DayOfWeek( year, month, day, c->UTC.DowStr );
    c->UTC.JD     = Lgm_JD( c->UTC.Year, c->UTC.Month, c->UTC.Day, c->UTC.Time, LGM_TIME_SYS_UTC, c );
    c->UTC.T      = (c->UTC.JD - 2451545.0)/36525.0;
    c->UTC.fYear = (double)c->UTC.Year + ((double)c->UTC.Doy - 1.0 + c->UTC.Time/24.0)/(365.0 + (double)Lgm_LeapYear(c->UTC.Year));

    // set DAT
    c->DAT        = Lgm_GetLeapSeconds( c->UTC.JD, c );




    /*
     * Set UT1 values
     */
    c->UT1      = c->UTC; // initially make DateTime structures the same
    c->UT1.Time = c->UTC.Time + c->DUT1/3600.0;
    c->UT1.JD   = Lgm_JD( c->UT1.Year, c->UT1.Month, c->UT1.Day, c->UT1.Time, LGM_TIME_SYS_UT1, c );
    c->UT1.Date = Lgm_JD_to_Date( c->UT1.JD, &c->UT1.Year, &c->UT1.Month, &c->UT1.Day, &Time );
    c->UT1.Time = Lgm_hour24( c->UT1.Time ); // Keep the time we had rather than gettingm it back from Lgm_JD_to_Date(). This limits roundoff error
    c->UT1.T    = (c->UT1.JD - 2451545.0)/36525.0;
    T_UT1       = c->UT1.T;
    T2_UT1      = T_UT1*T_UT1;
    T3_UT1      = T_UT1*T2_UT1;


    /*
     * Set TAI values
     */
    c->TAI      = c->UTC; // initially make DateTime structures the same
    c->TAI.Time = c->UTC.Time + c->DAT/3600.0;
    c->TAI.JD   = Lgm_JD( c->TAI.Year, c->TAI.Month, c->TAI.Day, c->TAI.Time, LGM_TIME_SYS_TAI, c );
    c->TAI.Date = Lgm_JD_to_Date( c->TAI.JD, &c->TAI.Year, &c->TAI.Month, &c->TAI.Day, &Time );
    c->TAI.Time = Lgm_hour24( c->TAI.Time ); // Keep the time we had rather than gettingm it back from Lgm_JD_to_Date(). This limits roundoff error
    c->TAI.T    = (c->TAI.JD - 2451545.0)/36525.0;


    /*
     * Set TT values
     */
    c->TT       = c->TAI; // initially make DateTime structure the same as TAI
    c->TT.Time  = c->TAI.Time + 32.184/3600.0;
    c->TT.JD    = Lgm_JD( c->TT.Year, c->TT.Month, c->TT.Day, c->TT.Time, LGM_TIME_SYS_TT, c );
    c->TT.Date  = Lgm_JD_to_Date( c->TT.JD, &c->TT.Year, &c->TT.Month, &c->TT.Day, &Time );
    c->TT.Time  = Lgm_hour24( c->TT.Time ); // Keep the time we had rather than gettingm it back from Lgm_JD_to_Date(). This limits roundoff error
    c->TT.T     = (c->TT.JD - 2451545.0)/36525.0;
    T_TT        = c->TT.T;
    T2_TT       = T_TT*T_TT;
    T3_TT       = T_TT*T2_TT;
    T4_TT       = T2_TT*T2_TT;
    T5_TT       = T3_TT*T2_TT;



    Lgm_TT_to_TDB( &c->TT, &c->TDB, c );



    /*
     *  Compute Greenwich Mean Sidereal Time (gmst)
     *  T0 = number of Julian centuries since 2000 January 1.5 to 0h on the date we are interested in.
     *  T  = number of Julian centuries since 2000 January 1.5 to the UT on the date we are interested in.
     *  (See Astronomical Almanac or e.g. Vallado.)
     */
    c->gmst = fmod( (67310.54841 + (876600.0*3600.0 + 8640184.812866)*T_UT1  + 0.093104*T2_UT1 - 6.2e-6*T3_UT1)/3600.0, 24.0);
    gmst    = c->gmst*15.0*RadPerDeg; // convert to radians for ease later on








    /*
     *
     *   Construct Transformation Matrix from MOD to GSE  systems
     *
     *
     *   First compute:
     *          mean ecliptic longitude of sun at epoch TU (varep)
     *          elciptic longitude of perigee at epoch TU (varpi)
     *          eccentricity of orbit at epoch TU (eccen)
     *
     *   The TU here is the number of Julian centuries since
     *   1900 January 0.0 (= 2415020.0)
     */
//    TDT   = Lgm_UTC_to_TT( c->JD, UT );
//    c->TT = TDT;
//UPDATE THIS TO IAU-76/FK-5
    TU    = ( Lgm_JD( c->TT.Year,  c->TT.Month,  c->TT.Day, c->TT.Time, LGM_TIME_SYS_TT, c) - 2415020.0 )/36525.0;
    varep = (279.6966778 + 36000.76892*TU + 0.0003025*TU*TU)*RadPerDeg;
    varpi = (281.2208444 + 1.719175*TU + 0.000452778*TU*TU)*RadPerDeg;
    eccen = 0.01675104 - 0.0000418*TU - 0.000000126*TU*TU;
    c->eccentricity = eccen;



    /*
     *  Compute the Obliquity of the Ecliptic at epoch TU
     *  The TU in this formula is the number of Julian
     *  centuries since epoch 2000 January 1.5
     */
    //epsilon    = (84381.448 - 46.8150*T_TT - 0.00059*T2_TT + 0.001813*T3_TT)*RadPerDeg/3600.0; // IAU-76/FK5 value
    c->epsilon    = 84381.448 - 46.8150*T_TT - 0.00059*T2_TT + 0.001813*T3_TT; // IAU-76/FK5 value (arcsec)


    /*
     * Compute:
     *          Number of Days since epoch 1990.0 (days)
     *          The Mean Anomaly (M)
     *          The True Anomaly (nu)
     *	    The Eccentric Anomaly via Keplers equation (E)
     */
    days    = Lgm_JD( c->TT.Year,  c->TT.Month,  c->TT.Day, c->TT.Time, LGM_TIME_SYS_TT, c) - 2447891.5;
    varep90 = 279.403303*RadPerDeg;
    varpi90 = 282.768422*RadPerDeg;
    M       = Lgm_angle2pi(M_2PI/365.242191*days );
    M       = Lgm_angle2pi( M + varep90 - varpi90 );
    E       = Lgm_kepler(M, eccen);
    nu      = 2.0*atan( sqrt((1.0+eccen)/(1.0-eccen))*tan(E/2.0) );
    lambnew = Lgm_angle2pi(nu + varpi90);

    c->lambda_sun = lambnew;


    /*
     *  Compute distance from earth to the sun
     */
    r0 = 1.495985e8;  /* in km */
    earth_sun_distance  = r0*(1-eccen*eccen)/(1.0 + eccen*cos(nu))/Re;
    c->earth_sun_dist = earth_sun_distance;


    /*
     * Compute Right Ascension and Declination of the Sun
     * Low precision...
     */
    epsilon = c->epsilon*RadPerArcSec; // radians
    cos_epsilon = cos(epsilon); sin_epsilon = sin(epsilon); sin_lambnew = sin(lambnew);
    RA  = Lgm_angle360(atan2(sin_lambnew*cos_epsilon, cos(lambnew))*DegPerRad);
    DEC = asin(sin_epsilon*sin_lambnew)*DegPerRad;
    c->RA_sun  = RA;
    c->DEC_sun = DEC;


    if ( USE_HIGH_ACCURACY_SUN ) {
//        T = (Lgm_JD( year, month, day, c->TT, c ) - 2451545.0)/36525.0;
        Lgm_SunPosition( T_TT, &l, &r, &b );
        c->lambda_sun_ha = l;       // high accuracy values
        c->r_sun_ha      = r*AU/Re; // high accuracy values
        c->beta_sun_ha   = b;       // high accuracy values
        sin_b = sin(b); cos_b = cos(b); sin_l = sin(l);
        RA  = Lgm_angle360( atan2( cos_epsilon*cos_b*sin_l-sin_epsilon*sin_b, cos_b*cos(l) )*DegPerRad );
        DEC = asin( sin_epsilon*cos_b*sin_l + cos_epsilon*sin_b )*DegPerRad;
    	c->RA_sun_ha  = RA;  // high accuracy RA
    	c->DEC_sun_ha = DEC; // high accuracy DEC
    }




    /*
     * Compute Precession angles to go from J2000 Epoch (we call this GEI here) to Mean Of Date (MOD) system.
     * Mean of Date has corrections for precession, but not nutation.
     * (e.g. see "Reduction for Precession -- rigorous formulae" in Astronomical Almanac 2004 page B18)
     */
//    TU = (Lgm_JD( year, month, day, TDT, c ) - 2451545.0)/36525.0;
//    TU2 = TU*TU;
//    TU3 = TU*TU2;
//    Zeta  = (0.6406161*TU + 0.0000839*TU2 + 5.0e-6*TU3)*RadPerDeg;
//    Zee   = (0.6406161*TU + 0.0003041*TU2 + 5.1e-6*TU3)*RadPerDeg;
//    Theta = (0.5567530*TU - 0.0001185*TU2 - 1.16e-5*TU3)*RadPerDeg;
    c->Zeta  = 0.6406161*T_TT + 0.0000839*T2_TT + 5.0e-6*T3_TT;   // degrees
    c->Zee   = 0.6406161*T_TT + 0.0003041*T2_TT + 5.1e-6*T3_TT;   // degrees
    c->Theta = 0.5567530*T_TT - 0.0001185*T2_TT - 1.16e-5*T3_TT;  // degrees
    Zeta = c->Zeta*RadPerDeg; Theta = c->Theta*RadPerDeg; Zee = c->Zee*RadPerDeg;
    SinZee   = sin(Zee);   CosZee   = cos(Zee);
    SinZeta  = sin(Zeta);  CosZeta  = cos(Zeta);
    SinTheta = sin(Theta); CosTheta = cos(Theta);
    c->Agei_to_mod[0][0] = CosZeta*CosTheta*CosZee-SinZeta*SinZee; c->Agei_to_mod[1][0] = -SinZeta*CosTheta*CosZee-CosZeta*SinZee; c->Agei_to_mod[2][0] = -SinTheta*CosZee;
    c->Agei_to_mod[0][1] = CosZeta*CosTheta*SinZee+SinZeta*CosZee; c->Agei_to_mod[1][1] = -SinZeta*CosTheta*SinZee+CosZeta*CosZee; c->Agei_to_mod[2][1] = -SinTheta*SinZee;
    c->Agei_to_mod[0][2] = CosZeta*SinTheta;                       c->Agei_to_mod[1][2] = -SinZeta*SinTheta;                       c->Agei_to_mod[2][2] = CosTheta;
    for (i=0; i<3; i++){
        for (j=0; j<3; j++){
            c->Amod_to_gei[i][j] = c->Agei_to_mod[j][i];
        }
    }


    /*
     * Set up transformation to go from MOD->TOD
     * This requires calculating nutation quantities.
     * Here we do the full 106-term series. Could limit this if its too slow...
     */
//    c->OmegaMoon = (450160.398036 - 6962890.5431*T_TT + 7.4722*T2_TT + 0.007702*T3_TT - 0.00005939*T4_TT)*RadPerArcSec;
    c->OmegaMoon = fmod( 125.04455501 - (5.0*360.0    + 134.1361851)*T_TT + 0.0020756*T2_TT +  2.139e-6*T3_TT, 360.0); // degrees

    // Get dPSi and dEps in arcsec
    Lgm_Nutation( c->TT.T, c->nNutationTerms, &(c->dPsi), &(c->dEps) ); // does 106-term nutation series
    c->dPsi += c->ddPsi;    // apply EOP corrections (arcsec)
    c->dEps += c->ddEps;    // apply EOP corrections (arcsec)
    // Equation of the equinoxes EQ_Eq
    c->EQ_Eq = c->dPsi*cos(epsilon) + 0.00264*sin( c->OmegaMoon*RadPerDeg ) + 0.000063*sin( 2.0*c->OmegaMoon*RadPerDeg ); // arcsec
    // true obliquity of ecliptic
    c->epsilon_true = c->dEps + c->epsilon;  // arcsec
//printf("c->epsilon_true/3600.0 = %15.10lf\n", c->epsilon_true/3600.0);
//printf("c->dEps/3600.0 + c->epsilon/3600.0 = %15.10lf\n", c->dEps/3600.0 + c->epsilon/3600.0);
    tmp = c->dPsi*RadPerArcSec;           sdp = sin( tmp ); cdp = cos( tmp );
    tmp = c->epsilon*RadPerArcSec;        se  = sin( tmp ); ce  = cos( tmp );
    tmp = c->epsilon_true*RadPerArcSec;   set = sin( tmp ); cet = cos( tmp );
    c->Atod_to_mod[0][0] =  cdp;    c->Atod_to_mod[1][0] = sdp*cet;              c->Atod_to_mod[2][0] = set*sdp;
    c->Atod_to_mod[0][1] = -sdp*ce; c->Atod_to_mod[1][1] = cet*cdp*ce + set*se;  c->Atod_to_mod[2][1] = set*cdp*ce - se*cet;
    c->Atod_to_mod[0][2] = -sdp*se; c->Atod_to_mod[1][2] = cet*cdp*se - set*ce;  c->Atod_to_mod[2][2] = set*se*cdp + cet*ce;
    for (i=0; i<3; i++){
        for (j=0; j<3; j++){
            c->Amod_to_tod[i][j] = c->Atod_to_mod[j][i];
        }
    }






    /*
     * Set up transformation between PEF and TOD
     * Need to compute True (Aparent?) Sidereal Time GAST
     */
    c->gast = c->gmst*15.0 + c->EQ_Eq/3600.0;
    gast = c->gast*RadPerDeg; cs = cos( gast ); sn = sin( gast );
    c->Atod_to_pef[0][0] =  cs; c->Atod_to_pef[1][0] =  sn; c->Atod_to_pef[2][0] = 0.0;
    c->Atod_to_pef[0][1] = -sn; c->Atod_to_pef[1][1] =  cs; c->Atod_to_pef[2][1] = 0.0;
    c->Atod_to_pef[0][2] = 0.0; c->Atod_to_pef[1][2] = 0.0; c->Atod_to_pef[2][2] = 1.0;

    c->Apef_to_tod[0][0] =  cs; c->Apef_to_tod[1][0] = -sn; c->Apef_to_tod[2][0] = 0.0;
    c->Apef_to_tod[0][1] =  sn; c->Apef_to_tod[1][1] =  cs; c->Apef_to_tod[2][1] = 0.0;
    c->Apef_to_tod[0][2] = 0.0; c->Apef_to_tod[1][2] = 0.0; c->Apef_to_tod[2][2] = 1.0;






    /*
     * Set up transformation between TOD and TEME
     * Utod = Rz( dPsiCosEps ) Uteme
     *    CHECK THIS. Do we use Full equation of equinoxes???
     *    OR just dPsi*cos( eps )????
     */
    //sn = sin( c->dPsiCosEps ); cs = cos( c->dPsiCosEps );
//    double ang = c->dPsi * cos(c->epsilon_true);
//    sn = sin( ang ); cs = cos( ang );
//    c->Ateme_to_tod[0][0] =  cs; c->Ateme_to_tod[1][0] =  sn; c->Ateme_to_tod[2][0] = 0.0;
//    c->Ateme_to_tod[0][1] = -sn; c->Ateme_to_tod[1][1] =  cs; c->Ateme_to_tod[2][1] = 0.0;
//    c->Ateme_to_tod[0][2] = 0.0; c->Ateme_to_tod[1][2] = 0.0; c->Ateme_to_tod[2][2] = 1.0;
//
//    c->Atod_to_teme[0][0] =  cs; c->Atod_to_teme[1][0] = -sn; c->Atod_to_teme[2][0] = 0.0;
//    c->Atod_to_teme[0][1] =  sn; c->Atod_to_teme[1][1] =  cs; c->Atod_to_teme[2][1] = 0.0;
//    c->Atod_to_teme[0][2] = 0.0; c->Atod_to_teme[1][2] = 0.0; c->Atod_to_teme[2][2] = 1.0;

    /*
     * Lets set this up differently because Vallado and others just arent clear on how exactly to
     * go from TEME to TOD (too many approximations seem to be available). Instead, follow Vallado's
     * recommendation and go the PEF instead. This is also good because we have gmst already!
     * See Vallado's book and his article "Revisiting Spacetrack Report #3".
     *   Upef = Rz( gmst82 ) Uteme
     */
      sn = sin( gmst ); cs = cos( gmst );
      c->Ateme_to_pef[0][0] =  cs; c->Ateme_to_pef[1][0] =  sn; c->Ateme_to_pef[2][0] = 0.0;
      c->Ateme_to_pef[0][1] = -sn; c->Ateme_to_pef[1][1] =  cs; c->Ateme_to_pef[2][1] = 0.0;
      c->Ateme_to_pef[0][2] = 0.0; c->Ateme_to_pef[1][2] = 0.0; c->Ateme_to_pef[2][2] = 1.0;

      c->Apef_to_teme[0][0] =  cs; c->Apef_to_teme[1][0] = -sn; c->Apef_to_teme[2][0] = 0.0;
      c->Apef_to_teme[0][1] =  sn; c->Apef_to_teme[1][1] =  cs; c->Apef_to_teme[2][1] = 0.0;
      c->Apef_to_teme[0][2] = 0.0; c->Apef_to_teme[1][2] = 0.0; c->Apef_to_teme[2][2] = 1.0;






    /*
     * Set up transformation between WGS84 and PEF
     * Uwgs84 = Ry( -xp )Rz( -yp ) Upef
     */
    sxp = sin( c->xp*RadPerArcSec ); cxp = cos( c->xp*RadPerArcSec );
    syp = sin( c->yp*RadPerArcSec ); cyp = cos( c->yp*RadPerArcSec );
    c->Apef_to_wgs84[0][0] =  cxp;     c->Apef_to_wgs84[1][0] =  sxp*syp; c->Apef_to_wgs84[2][0] = sxp*cyp;
    c->Apef_to_wgs84[0][1] =  0.0;     c->Apef_to_wgs84[1][1] =      cyp; c->Apef_to_wgs84[2][1] = -syp;
    c->Apef_to_wgs84[0][2] = -sxp;     c->Apef_to_wgs84[1][2] =  cxp*syp; c->Apef_to_wgs84[2][2] = cxp*cyp;

    c->Awgs84_to_pef[0][0] =  cxp;     c->Awgs84_to_pef[1][0] =      0.0; c->Awgs84_to_pef[2][0] = -sxp;
    c->Awgs84_to_pef[0][1] =  sxp*syp; c->Awgs84_to_pef[1][1] =      cyp; c->Awgs84_to_pef[2][1] =  cxp*syp;
    c->Awgs84_to_pef[0][2] =  sxp*cyp; c->Awgs84_to_pef[1][2] =     -syp; c->Awgs84_to_pef[2][2] =  cxp*cyp;


//    c->Apef_to_wgs84[0][0] =    1.0; c->Apef_to_wgs84[1][0] =    0.0; c->Apef_to_wgs84[2][0] =  c->xp;
//    c->Apef_to_wgs84[0][1] =    0.0; c->Apef_to_wgs84[1][1] =    1.0; c->Apef_to_wgs84[2][1] = -c->yp;
//    c->Apef_to_wgs84[0][2] = -c->xp; c->Apef_to_wgs84[1][2] =  c->yp; c->Apef_to_wgs84[2][2] =    1.0;
//
//    c->Awgs84_to_pef[0][0] =    1.0; c->Awgs84_to_pef[1][0] =     0.0; c->Awgs84_to_pef[2][0] = -c->xp;
//    c->Awgs84_to_pef[0][1] =    0.0; c->Awgs84_to_pef[1][1] =     1.0; c->Awgs84_to_pef[2][1] =  c->yp;
//    c->Awgs84_to_pef[0][2] =  c->xp; c->Awgs84_to_pef[1][2] =  -c->yp; c->Awgs84_to_pef[2][2] =    1.0;






    /*
     * Compute Right Ascension and Declination of the Moon
     */
    Lmoon_0 = 318.351648;
    P0 = 36.340410;
    N0 = 318.510107;
    Imoon = 5.145396;
    days  = Lgm_JD( c->TT.Year,  c->TT.Month,  c->TT.Day, c->TT.Time, LGM_TIME_SYS_TT, c) - 2447891.5;
    Lmoon = 13.1763966*days + Lmoon_0;
    Mmoon_m = Lmoon - 0.1114041*days - P0;
    Nmoon = N0 - 0.0529539*days;
    Cmoon = Lmoon*RadPerDeg - lambnew;
    Emoon_nu = 1.2739*sin(2.0*Cmoon - Mmoon_m*RadPerDeg);
    Amoon_e = 0.1858*sin(M);
    Amoon_3 = 0.37*sin(M);
    Mmoon_mp = Mmoon_m + Emoon_nu - Amoon_e - Amoon_3;
    Emoon_c = 6.2886*sin(Mmoon_mp*RadPerDeg);
    Amoon_4 = 0.214*sin(2.0*Mmoon_mp*RadPerDeg);
    Lmoon_p = Lmoon + Emoon_nu + Emoon_c - Amoon_e + Amoon_4;
    Vmoon = 0.6583*sin(2.0*(Lmoon_p*RadPerDeg - lambnew));
    Lmoon_pp = Lmoon_p + Vmoon;
    Nmoon_p = Nmoon - 0.16*sin(M);
    LambdaMoon = atan2( sin((Lmoon_pp - Nmoon_p)*RadPerDeg)*cos(Imoon*RadPerDeg),
			cos((Lmoon_pp - Nmoon_p)*RadPerDeg) ) + Nmoon_p*RadPerDeg;
    BetaMoon = asin(sin((Lmoon_pp - Nmoon_p)*RadPerDeg)*sin(Imoon*RadPerDeg));
    RA_Moon  = Lgm_angle360(atan2(sin(LambdaMoon)*cos(epsilon)-tan(BetaMoon)*sin(epsilon),
		cos(LambdaMoon))*DegPerRad);
    DEC_Moon = asin( sin(BetaMoon)*cos(epsilon) + cos(BetaMoon)*sin(epsilon)*sin(LambdaMoon))*DegPerRad;
    c->RA_moon = RA_Moon;
    c->DEC_moon = DEC_Moon;


    /*
     * Compute Phase of the Moon
     */
    c->MoonPhase = 0.5*(1.0 - cos(Lmoon_pp*RadPerDeg - lambnew));


    /*
     * Compute Earth-Moon distance
     */
    amoon = 384401.0;
    emoon = 0.054900;
    c->EarthMoonDistance = amoon*(1.0 - emoon*emoon)/(1.0 + emoon*cos((Mmoon_mp + Emoon_c)*RadPerDeg));
    c->EarthMoonDistance /= Re;




    /*
     *  Define MOD <--> GSE Transformation
     */

    // Direction of the Sun in MOD coords
    S.x = cos(lambnew);
    S.y = cos(epsilon)*sin(lambnew);
    S.z = sin(epsilon)*sin(lambnew);
    c->Sun = S;

    // Compute the direction of the Ecliptic Pole in the MOD system
    K.x =  0.0;
    K.y = -sin(epsilon);
    K.z =  cos(epsilon);
    c->EcPole = K;

    //  Taking X in direction of sun and Z parrallel to K (ec. pole),
    //  compute the Y axis direction.
    Lgm_CrossProduct(&K, &S, &Y);

    //  Create transformation matrix from GEI to GSE (c->Amod_to_gse[cols][rows])
    c->Amod_to_gse[0][0] = S.x, c->Amod_to_gse[1][0] = S.y, c->Amod_to_gse[2][0] = S.z;
    c->Amod_to_gse[0][1] = Y.x, c->Amod_to_gse[1][1] = Y.y, c->Amod_to_gse[2][1] = Y.z;
    c->Amod_to_gse[0][2] = K.x, c->Amod_to_gse[1][2] = K.y, c->Amod_to_gse[2][2] = K.z;

    c->Agse_to_mod[0][0] = S.x, c->Agse_to_mod[1][0] = Y.x, c->Agse_to_mod[2][0] = K.x;
    c->Agse_to_mod[0][1] = S.y, c->Agse_to_mod[1][1] = Y.y, c->Agse_to_mod[2][1] = K.y;
    c->Agse_to_mod[0][2] = S.z, c->Agse_to_mod[1][2] = Y.z, c->Agse_to_mod[2][2] = K.z;









    /*
     *  Initialize IGRF to get (among other things) the dipole moment.
     *  We really should do this more carefully. Make sure the lat/lons
     *  are geocentric?
     */
    N = 10;
//    c->UTC.fYear = (double)c->year + ((double)c->doy + c->UTC/24.0)/(365.0 + (double)Lgm_LeapYear(c->year));
    c->Lgm_IGRF_FirstCall = TRUE;
    Lgm_InitIGRF( g, h, N, 1, c );
    gclat = c->CD_gcolat;
    glon  = c->CD_glon;






    /*
     *   Construct Transformation Matrix from MOD to GSM
     *
     *   First compute the location of the earths D axis based on IGRF
     *   coefficients. Note: we compute Dhat as thoughj it were sticking out of the north
     *   hemisphere. It really goes out the south, but its less confusing this way...
     *
     */
    c->CD_gcolat = gclat*DegPerRad;
    c->CD_glon   = glon*DegPerRad;
    D.x   = cos(glon)*sin(gclat);
    D.y   = sin(glon)*sin(gclat);
    D.z   = cos(gclat);

    //  Find D (i.e. dipole) axis in MOD system
    Dmod.x = cos(gmst)*D.x - sin(gmst)*D.y;
    Dmod.y = sin(gmst)*D.x + cos(gmst)*D.y;
    Dmod.z = D.z;

    // Compute Ygsm axis in MOD system
    Lgm_CrossProduct(&Dmod, &S, &Y);

    //  Normalize the Y vector
    Lgm_NormalizeVector(&Y);

    //  Compute Zgsm axis in MOD system
    Lgm_CrossProduct(&S, &Y, &Z);

    c->Amod_to_gsm[0][0] = S.x, c->Amod_to_gsm[1][0] = S.y, c->Amod_to_gsm[2][0] = S.z;
    c->Amod_to_gsm[0][1] = Y.x, c->Amod_to_gsm[1][1] = Y.y, c->Amod_to_gsm[2][1] = Y.z;
    c->Amod_to_gsm[0][2] = Z.x, c->Amod_to_gsm[1][2] = Z.y, c->Amod_to_gsm[2][2] = Z.z;

    c->Agsm_to_mod[0][0] = S.x, c->Agsm_to_mod[1][0] = Y.x, c->Agsm_to_mod[2][0] = Z.x;
    c->Agsm_to_mod[0][1] = S.y, c->Agsm_to_mod[1][1] = Y.y, c->Agsm_to_mod[2][1] = Z.y;
    c->Agsm_to_mod[0][2] = S.z, c->Agsm_to_mod[1][2] = Y.z, c->Agsm_to_mod[2][2] = Z.z;








    /*
     *  Compute the Dipole Tilt angle. First convert Dmod into GSM.
     *  Then psi is just the angle between Dgsm and Zgsm
     */
    Dgsm.x = c->Amod_to_gsm[0][0]*Dmod.x + c->Amod_to_gsm[1][0]*Dmod.y + c->Amod_to_gsm[2][0]*Dmod.z;
    Dgsm.y = c->Amod_to_gsm[0][1]*Dmod.x + c->Amod_to_gsm[1][1]*Dmod.y + c->Amod_to_gsm[2][1]*Dmod.z;
    Dgsm.z = c->Amod_to_gsm[0][2]*Dmod.x + c->Amod_to_gsm[1][2]*Dmod.y + c->Amod_to_gsm[2][2]*Dmod.z;

    spsi       = sqrt( Dgsm.x*Dgsm.x + Dgsm.y*Dgsm.y );
    psi        = asin( spsi );
    psi       *= ( Dgsm.x < 0.0 ) ? -1.0 : 1.0;
    c->psi     = psi;
    c->cos_psi = cos( psi );
    c->sin_psi = sin( psi );
    c->tan_psi = c->sin_psi/c->cos_psi;


    /*
     *  Construct trans. matrices for GSM <-> SM
     */
    c->Agsm_to_sm[0][0] = c->cos_psi, c->Agsm_to_sm[1][0] = 0.0, c->Agsm_to_sm[2][0] = -c->sin_psi;
    c->Agsm_to_sm[0][1] = 0.0,        c->Agsm_to_sm[1][1] = 1.0, c->Agsm_to_sm[2][1] =  0.0;
    c->Agsm_to_sm[0][2] = c->sin_psi, c->Agsm_to_sm[1][2] = 0.0, c->Agsm_to_sm[2][2] =  c->cos_psi;
    Lgm_Transpose( c->Agsm_to_sm, c->Asm_to_gsm );


    /*
     *  Now do WGS84 (i.e. GEO) <-> CDMAG
     *    Z: parallel to dipole axis.
     *    Y: perp to both Z and Rot axis.
     *    X: completes
     */
    Zgeo.x = 0.0; Zgeo.y = 0.0; Zgeo.z = 1.0;
    Lgm_CrossProduct(&Zgeo, &D, &Y); // Ycdmag-hat in geo coords.
    Lgm_NormalizeVector(&Y);
    Lgm_CrossProduct(&Y, &D, &X);    // Xcdmag-hat in geo coords.
    Lgm_NormalizeVector(&X);
    c->Awgs84_to_cdmag[0][0] = X.x, c->Awgs84_to_cdmag[1][0] = X.y, c->Awgs84_to_cdmag[2][0] =  X.z;
    c->Awgs84_to_cdmag[0][1] = Y.x, c->Awgs84_to_cdmag[1][1] = Y.y, c->Awgs84_to_cdmag[2][1] =  Y.z;
    c->Awgs84_to_cdmag[0][2] = D.x, c->Awgs84_to_cdmag[1][2] = D.y, c->Awgs84_to_cdmag[2][2] =  D.z;

    c->Acdmag_to_wgs84[0][0] = X.x, c->Acdmag_to_wgs84[1][0] = Y.x, c->Acdmag_to_wgs84[2][0] =  D.x;
    c->Acdmag_to_wgs84[0][1] = X.y, c->Acdmag_to_wgs84[1][1] = Y.y, c->Acdmag_to_wgs84[2][1] =  D.y;
    c->Acdmag_to_wgs84[0][2] = X.z, c->Acdmag_to_wgs84[1][2] = Y.z, c->Acdmag_to_wgs84[2][2] =  D.z;






    /*  Construct transformation matrix from GSM to GSE and vice versa ;
     * 	Agsm_to_gse = Amod_to_gse * Agsm_to_mod
     */
    Lgm_MatTimesMat( c->Amod_to_gse, c->Agsm_to_mod, c->Agsm_to_gse );
    Lgm_Transpose( c->Agsm_to_gse, c->Agse_to_gsm);
/*
    for (j=0; j<3; ++j){
        for (i=0; i<3; ++i){
            c->Agsm_to_gse[i][j] = c->Amod_to_gse[0][j]*c->Amod_to_gsm[0][i]
                                    + c->Amod_to_gse[1][j]*c->Amod_to_gsm[1][i]
                                    + c->Amod_to_gse[2][j]*c->Amod_to_gsm[2][i];
            c->Agse_to_gsm[i][j] = c->Amod_to_gsm[0][j]*c->Amod_to_gse[0][i]
                                    + c->Amod_to_gsm[1][j]*c->Amod_to_gse[1][i]
                                    + c->Amod_to_gsm[2][j]*c->Amod_to_gse[2][i];
        }
    }
*/


    /*  Construct transformation matrix from GSM to WGS84 and vice versa ;
     * 	Agsm_to_wgs84 = Amod_to_wgs84 * Agsm_to_mod
     * 	Awgs84_to_gsm = Transpose( Agsm_to_wgs84 )
     */
    Lgm_MatTimesMat( c->Amod_to_wgs84, c->Agsm_to_mod, c->Agsm_to_wgs84 );
    Lgm_Transpose( c->Agsm_to_wgs84, c->Awgs84_to_gsm );
/*
    for (j=0; j<3; ++j){
        for (i=0; i<3; ++i){
            c->Agsm_to_wgs84[i][j] = c->Amod_to_wgs84[0][j]*c->Agsm_to_mod[i][0]
                                    + c->Amod_to_wgs84[1][j]*c->Agsm_to_mod[i][1]
                                    + c->Amod_to_wgs84[2][j]*c->Agsm_to_mod[i][2];
            c->Awgs84_to_gsm[i][j] = c->Amod_to_gsm[0][j]*c->Awgs84_to_mod[i][0]
                                    + c->Amod_to_gsm[1][j]*c->Awgs84_to_mod[i][1]
                                    + c->Amod_to_gsm[2][j]*c->Awgs84_to_mod[i][2];
        }
    }
*/



    /*
     *  Construct Transformation Matricies between  WGS84 and MOD
     *  and between WGS84 and GEI
     */
    Lgm_MatTimesMat( c->Apef_to_tod, c->Awgs84_to_pef, Tmp );
    Lgm_MatTimesMat( c->Atod_to_mod, Tmp, c->Awgs84_to_mod );
    Lgm_MatTimesMat( c->Amod_to_gei, c->Awgs84_to_mod, c->Awgs84_to_gei );
    Lgm_Transpose( c->Awgs84_to_mod, c->Amod_to_wgs84 );
    Lgm_Transpose( c->Awgs84_to_gei, c->Agei_to_wgs84 );







    if ( c->Verbose ) {

	printf("Time Quantitites:\n");
        printf("    fYear (UTC)       = %lf\n", c->UTC.fYear);
//        printf("    Date              = %ld\n", c->Date);
        printf("    UTC               = %8ld  %.14lf ( ", c->UTC.Date, c->UTC.Time); Lgm_Print_HMSdp( c->UTC.Time, TRUE, 8 ); printf(" )\n");
        printf("    UT1               = %8ld  %.14lf ( ", c->UT1.Date, c->UT1.Time); Lgm_Print_HMSdp( c->UT1.Time, TRUE, 8 ); printf(" )\n");
        printf("    TAI               = %8ld  %.14lf ( ", c->TAI.Date, c->TAI.Time); Lgm_Print_HMSdp( c->TAI.Time, TRUE, 8 ); printf(" )\n");
        printf("    TT = TAI+32.184s  = %8ld  %.14lf ( ", c->TT.Date,  c->TT.Time);  Lgm_Print_HMSdp( c->TT.Time,  TRUE, 8 ); printf(" )\n");
        printf("    TDB               = %8ld  %.14lf ( ", c->TDB.Date, c->TDB.Time); Lgm_Print_HMSdp( c->TDB.Time, TRUE, 8 ); printf(" )\n");

        printf("    DUT1 = UT1-UTC    = %.7lf seconds\n", c->DUT1);
        printf("    DAT  = TAI-UTC    = %.7lf seconds\n", c->DAT);
        printf("    JD (UTC)          = %.12lf\n", c->UTC.JD);
        printf("    JD (UT1)          = %.12lf\n", c->UT1.JD);
        printf("    JD (TT)           = %.12lf\n", c->TT.JD);
        printf("    T  (UTC)          = %15.12lf Julian Centuries of UTC\n", c->UTC.T);
        printf("    T  (UT1)          = %15.12lf Julian Centuries of UT1\n", c->UT1.T);
        printf("    T  (TT)           = %15.12lf Julian Centuries of TT\n",  c->TT.T);

        printf("    Year   (UTC)      = %d\n", c->UTC.Year);
        printf("    Month  (UTC)      = %d\n", c->UTC.Month);
        printf("    Day    (UTC)      = %d\n", c->UTC.Day);
        printf("    Doy    (UTC)      = %d\n", c->UTC.Doy);
        printf("    Dow    (UTC)      = %d\n", c->UTC.Dow);
        printf("    Dowstr (UTC)      = %s\n", c->UTC.DowStr);

        printf("    gmst (hours)      = %.7lf ( ", c->gmst); Lgm_Print_HMSd( c->gmst ); printf(" )\n");
        printf("    gmst (degrees)    = %.7lf ( ", c->gmst*15.0); Lgm_Print_DMSd( c->gmst*15.0 ); printf(" )\n");
        printf("    gast (hours)      = %.7lf ( ", c->gast/15.0); Lgm_Print_HMSd( c->gast/15.0 ); printf(" )\n");
        printf("    gast (degrees)    = %.7lf ( ", c->gast); Lgm_Print_DMSd( c->gast ); printf(" )\n");

	printf("\n");
	printf("Eccentricity and Obliquity:\n");
        printf("    eccentricity                      = %.8lf\n", c->eccentricity);
        printf("    epsilon mean (obliq. of ecliptic) = %.8lf ( ", c->epsilon/3600.0);      Lgm_Print_DMSd( c->epsilon/3600.0 ); printf(" )\n");
        printf("    epsilon true (obliq. of ecliptic) = %.8lf ( ", c->epsilon_true/3600.0); Lgm_Print_DMSd( c->epsilon_true/3600.0 ); printf(" )\n");

	printf("\n");
	printf("Precession Quantities:\n");
        printf("    Zeta              = %g ( ", c->Zeta);  Lgm_Print_DMSd( c->Zeta );  printf(" )\n");
        printf("    Zee               = %g ( ", c->Zee);   Lgm_Print_DMSd( c->Zee );   printf(" )\n");
        printf("    Theta             = %g ( ", c->Theta); Lgm_Print_DMSd( c->Theta ); printf(" )\n");

	printf("\n");
	printf("Nutation Quantities:\n");
        printf("    dPsi (w.o. corrections)           = %.8lf ( ", (c->dPsi - c->ddPsi)/3600.0); Lgm_Print_DMSd( (c->dPsi - c->ddPsi)/3600.0 ); printf(" )\n");
        printf("    dEps (w.o. corrections)           = %.8lf ( ", (c->dEps - c->ddEps)/3600.0); Lgm_Print_DMSd( (c->dEps - c->ddEps)/3600.0 ); printf(" )\n");
        printf("    ddPsi (EOP correction)            = %.8lf ( ", c->ddPsi/3600.0);             Lgm_Print_DMSd( c->ddPsi/3600.0 ); printf(" )\n");
        printf("    ddEps (EOP correction)            = %.8lf ( ", c->ddEps/3600.0);             Lgm_Print_DMSd( c->ddEps/3600.0 ); printf(" )\n");
        printf("    dPsi (w. corrections)             = %.8lf ( ", c->dPsi/3600.0);              Lgm_Print_DMSd( c->dPsi/3600.0 ); printf(" )\n");
        printf("    dEps (w. corrections)             = %.8lf ( ", c->dEps/3600.0);              Lgm_Print_DMSd( c->dEps/3600.0 ); printf(" )\n");
        printf("    epsilon true (obliq. of ecliptic) = %.8lf ( ", c->epsilon_true/3600.0);      Lgm_Print_DMSd( c->epsilon_true/3600.0 ); printf(" )\n");
        printf("    Equation of the Equinox           = %.8lf ( ", c->EQ_Eq/3600.0);             Lgm_Print_DMSd( c->EQ_Eq/3600.0 ); printf(" )\n");

	printf("\n");
	printf("Low Accuracy Position of Sun:\n");
        printf("    lambda_sun      = %15lf  ( ",  c->lambda_sun*DegPerRad);	Lgm_Print_DMSd( c->lambda_sun*DegPerRad );  printf(" )\n");
        printf("    earth_sun_dist  = %15lf Re\n", c->earth_sun_dist);
        printf("    beta_sun        = %14g   ( ",  0.0);	        Lgm_Print_DMSd( 0.0 );            printf(" )\n");
        printf("    RA_sun  (MOD)   = %15lf  ( ",  c->RA_sun);          Lgm_Print_HMSd( c->RA_sun/15.0 ); printf(" )\n");
        printf("    DEC_sun (MOD)   = %15lf  ( ",  c->DEC_sun);         Lgm_Print_DMSd( c->DEC_sun );     printf(" )\n");
        Lgm_Radec_to_Cart( c->RA_sun, c->DEC_sun, &u_mod );
        Lgm_Convert_Coords( &u_mod, &u_tod, MOD_TO_TOD, c );
        RA_tod = atan2( u_tod.y, u_tod.x)*DegPerRad; if (RA_tod<0.0) RA_tod += 360.0;
        DEC_tod = asin( u_tod.z )*DegPerRad;
        printf("    RA_sun  (TOD)   = %15lf  ( ",  RA_tod);      Lgm_Print_HMSd( RA_tod/15.0 ); printf(" )\n");
        printf("    DEC_sun (TOD)   = %15lf  ( ",  DEC_tod);     Lgm_Print_DMSd( DEC_tod );     printf(" )\n");
        Lgm_Convert_Coords( &u_mod, &u_gei, MOD_TO_GEI2000, c );
        RA_gei = atan2( u_gei.y, u_gei.x)*DegPerRad; if (RA_gei<0.0) RA_gei += 360.0;
        DEC_gei = asin( u_gei.z )*DegPerRad;
        printf("    RA_sun  (J2000) = %15lf  ( ",  RA_gei);      Lgm_Print_HMSd( RA_gei/15.0 ); printf(" )\n");
        printf("    DEC_sun (J2000) = %15lf  ( ",  DEC_gei);     Lgm_Print_DMSd( DEC_gei );     printf(" )\n");

	printf("\n");
	printf("High Accuracy Position of Sun:\n");
        printf("    lambda_sun_ha   = %15lf  ( ",  c->lambda_sun_ha*DegPerRad ); Lgm_Print_DMSd( c->lambda_sun_ha*DegPerRad );  printf(" )\n");
        printf("    r_sun_ha        = %15lf Re\n", c->r_sun_ha );
        printf("    beta_sun_ha     = %14g   ( ",  c->beta_sun_ha*DegPerRad );   Lgm_Print_DMSd( c->beta_sun_ha*DegPerRad );    printf(" )\n");
        printf("    RA_sun  (MOD)   = %15lf  ( ",  c->RA_sun_ha);      Lgm_Print_HMSd( c->RA_sun_ha/15.0 ); printf(" )\n");
        printf("    DEC_sun (MOD)   = %15lf  ( ",  c->DEC_sun_ha);     Lgm_Print_DMSd( c->DEC_sun_ha );     printf(" )\n");
        Lgm_Radec_to_Cart( c->RA_sun_ha, c->DEC_sun_ha, &u_mod );
        Lgm_Convert_Coords( &u_mod, &u_tod, MOD_TO_TOD, c );
        RA_tod = atan2( u_tod.y, u_tod.x)*DegPerRad; if (RA_tod<0.0) RA_tod += 360.0;
	    DEC_tod = asin( u_tod.z )*DegPerRad;
        printf("    RA_sun  (TOD)   = %15lf  ( ",  RA_tod);      Lgm_Print_HMSd( RA_tod/15.0 ); printf(" )\n");
        printf("    DEC_sun (TOD)   = %15lf  ( ",  DEC_tod);     Lgm_Print_DMSd( DEC_tod );     printf(" )\n");
        Lgm_Convert_Coords( &u_mod, &u_gei, MOD_TO_GEI2000, c );
        RA_gei = atan2( u_gei.y, u_gei.x)*DegPerRad; if (RA_gei<0.0) RA_gei += 360.0;
        DEC_gei = asin( u_gei.z )*DegPerRad;
        printf("    RA_sun  (J2000) = %15lf  ( ",  RA_gei);      Lgm_Print_HMSd( RA_gei/15.0 ); printf(" )\n");
        printf("    DEC_sun (J2000) = %15lf  ( ",  DEC_gei);     Lgm_Print_DMSd( DEC_gei );     printf(" )\n");


	printf("\n");
	printf("Sun vector and Ecliptic Pole in GEI2000:\n");
        printf("    Sun               = (%lf, %lf, %lf)\n", c->Sun.x, c->Sun.y, c->Sun.z);
        printf("    EcPole            = (%lf, %lf, %lf)\n", c->EcPole.x, c->EcPole.y, c->EcPole.z);

	printf("\n");
	printf("Geo-dipole tilt angle:\n");
        printf("    psi                      = %lf  ( ", c->psi*DegPerRad);     Lgm_Print_DMSd( c->psi*DegPerRad );     printf(" )\n");
        printf("    sin_psi                  = %lf\n", c->sin_psi);
        printf("    cos_psi                  = %lf\n", c->cos_psi);
        printf("    tan_psi                  = %lf\n", c->tan_psi);

	printf("\n");
	printf("Position of Moon:\n");
        printf("   RA_moon                      = %lf  ( ", c->RA_moon);        Lgm_Print_HMSd( c->RA_moon/15.0 ); printf(" )\n");
        printf("   DEC_moon                     = %lf  ( ", c->DEC_moon);       Lgm_Print_DMSd( c->DEC_moon );     printf(" )\n");
        printf("   EarthMoonDistance            = %lf\n", c->EarthMoonDistance);
        printf("   MoonPhase                    = %lf\n", c->MoonPhase);


	printf("\n");
	printf("IGRF-derived quantities:\n");
        printf("    M_cd              = %lf\n", c->M_cd);
        printf("    M_cd_McIllwain    = %lf\n", c->M_cd_McIllwain);
        printf("    CD_gcolat         = %lf (deg.)  ( ", c->CD_gcolat); Lgm_Print_DMSd( c->CD_gcolat ); printf(" )\n");
        printf("    CD_glon           = %lf (deg.)  ( ", c->CD_glon);   Lgm_Print_DMSd( c->CD_glon ); printf(" )\n");
        printf("    ED_x0             = %lf  Re  (%lf km)\n", c->ED_x0, c->ED_x0*Re);
        printf("    ED_y0             = %lf  Re  (%lf km)\n", c->ED_y0, c->ED_y0*Re);
        printf("    ED_z0             = %lf  Re  (%lf km)\n", c->ED_z0, c->ED_z0*Re);


	printf("\n");
	printf("Transformation Matrices:\n");
        printf("                        [ %15.8lf  %15.8lf  %15.8lf ]\n",   c->Amod_to_gse[0][0], c->Amod_to_gse[1][0], c->Amod_to_gse[2][0]);
        printf("    Amod_to_gse       = [ %15.8lf  %15.8lf  %15.8lf ]\n",   c->Amod_to_gse[0][1], c->Amod_to_gse[1][1], c->Amod_to_gse[2][1]);
        printf("                        [ %15.8lf  %15.8lf  %15.8lf ]\n\n", c->Amod_to_gse[0][2], c->Amod_to_gse[1][2], c->Amod_to_gse[2][2]);

        printf("                        [ %15.8lf  %15.8lf  %15.8lf ]\n",   c->Amod_to_gsm[0][0], c->Amod_to_gsm[1][0], c->Amod_to_gsm[2][0]);
        printf("    Amod_to_gsm       = [ %15.8lf  %15.8lf  %15.8lf ]\n",   c->Amod_to_gsm[0][1], c->Amod_to_gsm[1][1], c->Amod_to_gsm[2][1]);
        printf("                        [ %15.8lf  %15.8lf  %15.8lf ]\n\n", c->Amod_to_gsm[0][2], c->Amod_to_gsm[1][2], c->Amod_to_gsm[2][2]);

        printf("                        [ %15.8lf  %15.8lf  %15.8lf ]\n",   c->Agei_to_wgs84[0][0], c->Agei_to_wgs84[1][0], c->Agei_to_wgs84[2][0]);
        printf("    Agei_to_wgs84     = [ %15.8lf  %15.8lf  %15.8lf ]\n",   c->Agei_to_wgs84[0][1], c->Agei_to_wgs84[1][1], c->Agei_to_wgs84[2][1]);
        printf("                        [ %15.8lf  %15.8lf  %15.8lf ]\n\n", c->Agei_to_wgs84[0][2], c->Agei_to_wgs84[1][2], c->Agei_to_wgs84[2][2]);

        printf("                        [ %15.8lf  %15.8lf  %15.8lf ]\n",   c->Agse_to_mod[0][0], c->Agse_to_mod[1][0], c->Agse_to_mod[2][0]);
        printf("    Agse_to_mod       = [ %15.8lf  %15.8lf  %15.8lf ]\n",   c->Agse_to_mod[0][1], c->Agse_to_mod[1][1], c->Agse_to_mod[2][1]);
        printf("                        [ %15.8lf  %15.8lf  %15.8lf ]\n\n", c->Agse_to_mod[0][2], c->Agse_to_mod[1][2], c->Agse_to_mod[2][2]);

        printf("                        [ %15.8lf  %15.8lf  %15.8lf ]\n",   c->Agse_to_gsm[0][0], c->Agse_to_gsm[1][0], c->Agse_to_gsm[2][0]);
        printf("    Agse_to_gsm       = [ %15.8lf  %15.8lf  %15.8lf ]\n",   c->Agse_to_gsm[0][1], c->Agse_to_gsm[1][1], c->Agse_to_gsm[2][1]);
        printf("                        [ %15.8lf  %15.8lf  %15.8lf ]\n\n", c->Agse_to_gsm[0][2], c->Agse_to_gsm[1][2], c->Agse_to_gsm[2][2]);

        printf("                        [ %15.8lf  %15.8lf  %15.8lf ]\n",   c->Awgs84_to_gei[0][0], c->Awgs84_to_gei[1][0], c->Awgs84_to_gei[2][0]);
        printf("    Awgs84_to_gei     = [ %15.8lf  %15.8lf  %15.8lf ]\n",   c->Awgs84_to_gei[0][1], c->Awgs84_to_gei[1][1], c->Awgs84_to_gei[2][1]);
        printf("                        [ %15.8lf  %15.8lf  %15.8lf ]\n\n", c->Awgs84_to_gei[0][2], c->Awgs84_to_gei[1][2], c->Awgs84_to_gei[2][2]);

        printf("                        [ %15.8lf  %15.8lf  %15.8lf ]\n",   c->Agsm_to_mod[0][0], c->Agsm_to_mod[1][0], c->Agsm_to_mod[2][0]);
        printf("    Agsm_to_mod       = [ %15.8lf  %15.8lf  %15.8lf ]\n",   c->Agsm_to_mod[0][1], c->Agsm_to_mod[1][1], c->Agsm_to_mod[2][1]);
        printf("                        [ %15.8lf  %15.8lf  %15.8lf ]\n\n", c->Agsm_to_mod[0][2], c->Agsm_to_mod[1][2], c->Agsm_to_mod[2][2]);

        printf("                        [ %15.8lf  %15.8lf  %15.8lf ]\n",   c->Agsm_to_sm[0][0], c->Agsm_to_sm[1][0], c->Agsm_to_sm[2][0]);
        printf("    Agsm_to_sm        = [ %15.8lf  %15.8lf  %15.8lf ]\n",   c->Agsm_to_sm[0][1], c->Agsm_to_sm[1][1], c->Agsm_to_sm[2][1]);
        printf("                        [ %15.8lf  %15.8lf  %15.8lf ]\n\n", c->Agsm_to_sm[0][2], c->Agsm_to_sm[1][2], c->Agsm_to_sm[2][2]);

        printf("                        [ %15.8lf  %15.8lf  %15.8lf ]\n",   c->Agsm_to_gse[0][0], c->Agsm_to_gse[1][0], c->Agsm_to_gse[2][0]);
        printf("    Agsm_to_gse       = [ %15.8lf  %15.8lf  %15.8lf ]\n",   c->Agsm_to_gse[0][1], c->Agsm_to_gse[1][1], c->Agsm_to_gse[2][1]);
        printf("                        [ %15.8lf  %15.8lf  %15.8lf ]\n\n", c->Agsm_to_gse[0][2], c->Agsm_to_gse[1][2], c->Agsm_to_gse[2][2]);

        printf("                        [ %15.8lf  %15.8lf  %15.8lf ]\n",   c->Asm_to_gsm[0][0], c->Asm_to_gsm[1][0], c->Asm_to_gsm[2][0]);
        printf("    Asm_to_gsm        = [ %15.8lf  %15.8lf  %15.8lf ]\n",   c->Asm_to_gsm[0][1], c->Asm_to_gsm[1][1], c->Asm_to_gsm[2][1]);
        printf("                        [ %15.8lf  %15.8lf  %15.8lf ]\n\n", c->Asm_to_gsm[0][2], c->Asm_to_gsm[1][2], c->Asm_to_gsm[2][2]);

        printf("                        [ %15.8lf  %15.8lf  %15.8lf ]\n",   c->Agei_to_mod[0][0], c->Agei_to_mod[1][0], c->Agei_to_mod[2][0]);
        printf("    Agei_to_mod       = [ %15.8lf  %15.8lf  %15.8lf ]\n",   c->Agei_to_mod[0][1], c->Agei_to_mod[1][1], c->Agei_to_mod[2][1]);
        printf("                        [ %15.8lf  %15.8lf  %15.8lf ]\n\n", c->Agei_to_mod[0][2], c->Agei_to_mod[1][2], c->Agei_to_mod[2][2]);

        printf("                        [ %15.8lf  %15.8lf  %15.8lf ]\n",   c->Amod_to_gei[0][0], c->Amod_to_gei[1][0], c->Amod_to_gei[2][0]);
        printf("    Amod_to_gei       = [ %15.8lf  %15.8lf  %15.8lf ]\n",   c->Amod_to_gei[0][1], c->Amod_to_gei[1][1], c->Amod_to_gei[2][1]);
        printf("                        [ %15.8lf  %15.8lf  %15.8lf ]\n\n", c->Amod_to_gei[0][2], c->Amod_to_gei[1][2], c->Amod_to_gei[2][2]);

        printf("                        [ %15.8lf  %15.8lf  %15.8lf ]\n",   c->Amod_to_tod[0][0], c->Amod_to_tod[1][0], c->Amod_to_tod[2][0]);
        printf("    Amod_to_tod       = [ %15.8lf  %15.8lf  %15.8lf ]\n",   c->Amod_to_tod[0][1], c->Amod_to_tod[1][1], c->Amod_to_tod[2][1]);
        printf("                        [ %15.8lf  %15.8lf  %15.8lf ]\n\n", c->Amod_to_tod[0][2], c->Agei_to_mod[1][2], c->Amod_to_tod[2][2]);

        printf("                        [ %15.8lf  %15.8lf  %15.8lf ]\n",   c->Atod_to_mod[0][0], c->Atod_to_mod[1][0], c->Atod_to_mod[2][0]);
        printf("    Atod_to_mod       = [ %15.8lf  %15.8lf  %15.8lf ]\n",   c->Atod_to_mod[0][1], c->Atod_to_mod[1][1], c->Atod_to_mod[2][1]);
        printf("                        [ %15.8lf  %15.8lf  %15.8lf ]\n\n", c->Atod_to_mod[0][2], c->Atod_to_mod[1][2], c->Atod_to_mod[2][2]);

        printf("                        [ %15.8lf  %15.8lf  %15.8lf ]\n",   c->Atod_to_pef[0][0], c->Atod_to_pef[1][0], c->Atod_to_pef[2][0]);
        printf("    Atod_to_pef       = [ %15.8lf  %15.8lf  %15.8lf ]\n",   c->Atod_to_pef[0][1], c->Atod_to_pef[1][1], c->Atod_to_pef[2][1]);
        printf("                        [ %15.8lf  %15.8lf  %15.8lf ]\n\n", c->Atod_to_pef[0][2], c->Atod_to_pef[1][2], c->Atod_to_pef[2][2]);

        printf("                        [ %15.8lf  %15.8lf  %15.8lf ]\n",   c->Apef_to_tod[0][0], c->Apef_to_tod[1][0], c->Apef_to_tod[2][0]);
        printf("    Apef_to_tod       = [ %15.8lf  %15.8lf  %15.8lf ]\n",   c->Apef_to_tod[0][1], c->Apef_to_tod[1][1], c->Apef_to_tod[2][1]);
        printf("                        [ %15.8lf  %15.8lf  %15.8lf ]\n\n", c->Apef_to_tod[0][2], c->Apef_to_tod[1][2], c->Apef_to_tod[2][2]);

        printf("                        [ %15.8e  %15.8e  %15.8e ]\n",   c->Ateme_to_pef[0][0], c->Ateme_to_pef[1][0], c->Ateme_to_pef[2][0]);
        printf("    Ateme_to_pef      = [ %15.8e  %15.8e  %15.8e ]\n",   c->Ateme_to_pef[0][1], c->Ateme_to_pef[1][1], c->Ateme_to_pef[2][1]);
        printf("                        [ %15.8e  %15.8e  %15.8e ]\n\n", c->Ateme_to_pef[0][2], c->Ateme_to_pef[1][2], c->Ateme_to_pef[2][2]);

        printf("                        [ %15.8e  %15.8e  %15.8e ]\n",   c->Apef_to_teme[0][0], c->Apef_to_teme[1][0], c->Apef_to_teme[2][0]);
        printf("    Apef_to_teme      = [ %15.8e  %15.8e  %15.8e ]\n",   c->Apef_to_teme[0][1], c->Apef_to_teme[1][1], c->Apef_to_teme[2][1]);
        printf("                        [ %15.8e  %15.8e  %15.8e ]\n\n", c->Apef_to_teme[0][2], c->Apef_to_teme[1][2], c->Apef_to_teme[2][2]);

        printf("                        [ %15.8e  %15.8e  %15.8e ]\n",   c->Awgs84_to_pef[0][0], c->Awgs84_to_pef[1][0], c->Awgs84_to_pef[2][0]);
        printf("    Awgs84_to_pef     = [ %15.8e  %15.8e  %15.8e ]\n",   c->Awgs84_to_pef[0][1], c->Awgs84_to_pef[1][1], c->Awgs84_to_pef[2][1]);
        printf("                        [ %15.8e  %15.8e  %15.8e ]\n\n", c->Awgs84_to_pef[0][2], c->Awgs84_to_pef[1][2], c->Awgs84_to_pef[2][2]);

        printf("                        [ %15.8e  %15.8e  %15.8e ]\n",   c->Apef_to_wgs84[0][0], c->Apef_to_wgs84[1][0], c->Apef_to_wgs84[2][0]);
        printf("    Apef_to_wgs84     = [ %15.8e  %15.8e  %15.8e ]\n",   c->Apef_to_wgs84[0][1], c->Apef_to_wgs84[1][1], c->Apef_to_wgs84[2][1]);
        printf("                        [ %15.8e  %15.8e  %15.8e ]\n\n", c->Apef_to_wgs84[0][2], c->Apef_to_wgs84[1][2], c->Apef_to_wgs84[2][2]);
    }

}


/*
 *  Routine converts coords from one system to another. All of the major
 *  systems aout to be defined -- see ctrans.h for flag definitions.
 *
 *     Lgm_Vector *u;      -- input vector
 *     Lgm_Vector *v;      -- output vector
 *     int     flag;       -- flag to specify the type of coord transformation
 *                            (see LgmCtrans.h for description of whats available)
 *     Lgm_CTrans *c;      -- structure that holds all the coord trans info
 *
 */
void Lgm_Convert_Coords(Lgm_Vector *u, Lgm_Vector *v, int flag, Lgm_CTrans *c) {

    Lgm_Vector 	w, r, z, e, d;
    int 	    inflag, outflag;

    /*
     *  Figure out what system the input coords are in
     *  and what system the output coords should be in...
     */
    inflag = (int)((float)flag/100.0);
    outflag = flag - inflag*100;

    /*
     * Catch case where in coords are same as out coords
     */
    if (inflag == outflag) {
        v->x = u->x;
        v->y = u->y;
        v->z = u->z;
        return;
    }

    /*
     *  Cant have trans matrices from everything to everything.
     *  So we pick a common system to transform in and out of.
     *  Lets choose MOD since its pretty common
     */


    /*
     *   convert coords to a common system: use MOD...
     */
    switch(inflag){
        case EME2000_COORDS: 	// in coords are in EME2000 (aka ICRF2000 or GEI2000)
             // aka ICRF2000_COORDS
             // aka GEI2000_COORDS
            Lgm_MatTimesVec(c->Agei_to_mod, u, &w);
            break;
        case MOD_COORDS: 	// in coords are in MOD
            w.x = u->x, w.y = u->y, w.z = u->z;
            break;
        case TOD_COORDS: 	// in coords are in TOD
            Lgm_MatTimesVec(c->Atod_to_mod, u, &w);
            break;
        case TEME_COORDS: 	// in coords are in TEME
            Lgm_MatTimesVec(c->Ateme_to_pef, u, &r);
            Lgm_MatTimesVec(c->Apef_to_tod, &r, &z);
            Lgm_MatTimesVec(c->Atod_to_mod, &z, &w);
            break;
        case PEF_COORDS: 	// in coords are in PEF
            Lgm_MatTimesVec(c->Apef_to_tod, u, &r);
            Lgm_MatTimesVec(c->Atod_to_mod, &r, &w);
            break;
        case WGS84_COORDS: 	// in coords are in WGS84
            Lgm_MatTimesVec(c->Awgs84_to_pef, u, &r);
            Lgm_MatTimesVec(c->Apef_to_tod, &r, &z);
            Lgm_MatTimesVec(c->Atod_to_mod, &z, &w);
            break;
        case GSE_COORDS: 	// in coords are in GSE
            Lgm_MatTimesVec(c->Agse_to_mod, u, &w);
            break;
        case GSM_COORDS: 	// in coords are in GSM
            Lgm_MatTimesVec(c->Agsm_to_mod, u, &w);
            break;
        case SM_COORDS: 	// in coords are in SM
            Lgm_MatTimesVec(c->Asm_to_gsm,  u, &r);
            Lgm_MatTimesVec(c->Agsm_to_mod, &r, &w);
            break;
        case CDMAG_COORDS: 	// in coords are in CDMAG
            Lgm_MatTimesVec(c->Acdmag_to_wgs84, u, &r);
            Lgm_MatTimesVec(c->Awgs84_to_mod, &r, &w);
            break;
        case EDMAG_COORDS: 	// in coords are in EDMAG
	        // Convert the  dipole offset to CDMAG coords first
	        d.x = c->ED_x0; d.y = c->ED_y0; d.z = c->ED_z0; // offset vector in wgs84 (i.e. in geo)
	        Lgm_MatTimesVec(c->Awgs84_to_cdmag, &d, &e);
            // now Rcd = d_offset_cd + Red
            z.x = u->x + e.x; z.y = u->y + e.y; z.z = u->z + e.z; // Rcd = R0 + Red  (z should be in CD now)
            Lgm_MatTimesVec(c->Acdmag_to_wgs84, &z, &r);
            Lgm_MatTimesVec(c->Awgs84_to_mod, &r, &w);
            break;
    }

    // w now holds the coords in MOD -- now transform to desired system

    switch(outflag){
        case EME2000_COORDS: 	// out coords are in EME2000 (aka ICRF2000 or GEI2000)
            Lgm_MatTimesVec(c->Amod_to_gei, &w, v);
            break;
        case MOD_COORDS: 	// out coords are in MOD
            v->x = w.x, v->y = w.y, v->z = w.z;
            break;
        case TOD_COORDS: 	// out coords are in TOD
            Lgm_MatTimesVec(c->Amod_to_tod, &w, v);
            break;
        case TEME_COORDS: 	// out coords are in TEME
            Lgm_MatTimesVec(c->Amod_to_tod, &w, &r);
            Lgm_MatTimesVec(c->Atod_to_pef, &r, &z);
            Lgm_MatTimesVec(c->Apef_to_teme, &z, v);
            break;
        case PEF_COORDS: 	// out coords are in PEF
            Lgm_MatTimesVec(c->Amod_to_tod, &w, &r);
            Lgm_MatTimesVec(c->Atod_to_pef, &r, v);
            break;
        case WGS84_COORDS: 	// out coords are in WGS84
            Lgm_MatTimesVec(c->Amod_to_tod, &w, &r);
            Lgm_MatTimesVec(c->Atod_to_pef, &r, &z);
            Lgm_MatTimesVec(c->Apef_to_wgs84, &z, v);
            break;
        case GSE_COORDS: 	// out coords are in GSE
            Lgm_MatTimesVec(c->Amod_to_gse, &w, v);
            break;
        case GSM_COORDS: 	// out coords are in GSM
            Lgm_MatTimesVec(c->Amod_to_gsm, &w, v);
            break;
        case SM_COORDS: 	// out coords are in SM
            Lgm_MatTimesVec(c->Amod_to_gsm, &w, &r);
            Lgm_MatTimesVec(c->Agsm_to_sm,  &r, v);
            break;
        case CDMAG_COORDS: 	// out coords are in CDMAG
            Lgm_MatTimesVec(c->Amod_to_wgs84, &w, &r);
            Lgm_MatTimesVec(c->Awgs84_to_cdmag, &r, v);
            break;
        case EDMAG_COORDS: 	// out coords are in ED
            Lgm_MatTimesVec(c->Amod_to_wgs84, &w, &r);
            Lgm_MatTimesVec(c->Awgs84_to_cdmag, &r, &z);
	        // Convert the  dipole offset to CDMAG coords
	        d.x = c->ED_x0; d.y = c->ED_y0; d.z = c->ED_z0; // offset vector in wgs84 (i.e. in geo)
	        Lgm_MatTimesVec(c->Awgs84_to_cdmag, &d, &e);
            // now Rcd - d_offset_cd = Red
            v->x = z.x - e.x; v->y = z.y - e.y; v->z = z.z - e.z; // Rcd - R0 = Red  (z should be in ED now)
            break;
        }


}


/*
 * Input:
 *          GeodLat:    in degrees
 *          GeodLong:   in degrees
 *          GeodHeight: in km
 *
 * Output:
 *          v: in Re
 */
void Lgm_GEOD_to_WGS84( double GeodLat, double GeodLong, double GeodHeight, Lgm_Vector *v ) {

    double      lam, phi, h, CosPhi, SinPhi, CosLam, SinLam, Chi;

    lam = GeodLat*RadPerDeg;
    phi = GeodLong*RadPerDeg;
    h   = GeodHeight;

    CosPhi = cos(phi); SinPhi = sin(phi); SinLam = sin(lam); CosLam = cos(lam);
    Chi = sqrt( 1.0 - WGS84_E2*SinLam*SinLam );

    // convert to GEO (in units of Re)
    v->x = (WGS84_A/Chi + h)*CosLam*CosPhi/WGS84_A;
    v->y = (WGS84_A/Chi + h)*CosLam*SinPhi/WGS84_A;
    v->z = (WGS84_A*(1.0-WGS84_E2)/Chi + h)*SinLam/WGS84_A;

}

void Lgm_WGS84_to_GEOD( Lgm_Vector *uin, double *GeodLat, double *GeodLong, double *GeodHeight ) {

    Lgm_Vector u;
    double  r, r2, z2, F, G, G2, G3, c, s, tt, tt2, P;
    double  Q, ro, U, V, zo;

    u = *uin;
    Lgm_ScaleVector( &u, WGS84_A ); // convert to km

    r2 = u.x*u.x + u.y*u.y;
    r = sqrt(r2);
    z2 = u.z*u.z;
    F = 54.0*WGS84_B2*z2;
    G = r2 + WGS84_1mE2*z2 - WGS84_E2*WGS84_A2mB2;
    G2 = G*G; G3 = G2*G;
    c = (WGS84_E4*F*r2)/G3;
    s = pow( 1.0 + c + sqrt(c*c + 2.0*c), M_OneThird );
    tt = s + 1.0/s + 1.0; tt2 = tt*tt;
    P = F/( 3.0*tt2*G2 );

    Q = sqrt( 1.0 + 2.0*WGS84_E4*P );
    ro = -(WGS84_E2*P*r)/(1.0+Q) + sqrt( (0.5*WGS84_A2)*(1.0+1.0/Q) - (WGS84_1mE2*P*z2)/(Q*(1.0+Q)) - 0.5*P*r2 );

    tt = (r - WGS84_E2*ro); tt2 = tt*tt;
    U   = sqrt( tt2 + z2 );
    V   = sqrt( tt2 + WGS84_1mE2*z2 );
    zo  = (WGS84_B2*u.z)/(WGS84_A*V);

    *GeodLat    = DegPerRad*atan( (u.z + WGS84_EP2*zo)/r );  // geodetic latitude
    *GeodLong   =  DegPerRad*atan2( u.y, u.x );              // geodetic longitude (same as GEO)
    *GeodHeight = U*( 1.0 - WGS84_B2/(WGS84_A*V));           // geodetic height (km)
}

// same as above, but just returns height (saves a few expensive atans)
void Lgm_WGS84_to_GeodHeight( Lgm_Vector *uin, double *GeodHeight ) {

    Lgm_Vector u;
    double  r, r2, z2, F, G, G2, G3, c, s, tt, tt2, P;
    double  Q, ro, U, V, zo;

    u = *uin;
    Lgm_ScaleVector( &u, WGS84_A ); // convert to km

    r2 = u.x*u.x + u.y*u.y;
    r = sqrt(r2);
    z2 = u.z*u.z;
    F = 54.0*WGS84_B2*z2;
    G = r2 + WGS84_1mE2*z2 - WGS84_E2*WGS84_A2mB2;
    G2 = G*G; G3 = G2*G;
    c = (WGS84_E4*F*r2)/G3;
    s = pow( 1.0 + c + sqrt(c*c + 2.0*c), M_OneThird );
    tt = s + 1.0/s + 1.0; tt2 = tt*tt;
    P = F/( 3.0*tt2*G2 );

    Q = sqrt( 1.0 + 2.0*WGS84_E4*P );
    ro = -(WGS84_E2*P*r)/(1.0+Q) + sqrt( (0.5*WGS84_A2)*(1.0+1.0/Q) - (WGS84_1mE2*P*z2)/(Q*(1.0+Q)) - 0.5*P*r2 );

    tt = (r - WGS84_E2*ro); tt2 = tt*tt;
    U   = sqrt( tt2 + z2 );
    V   = sqrt( tt2 + WGS84_1mE2*z2 );
    zo  = (WGS84_B2*u.z)/(WGS84_A*V);

    *GeodHeight = U*( 1.0 - WGS84_B2/(WGS84_A*V));           // geodetic height (km)
}



/*
 *  These rotuines are included here, just because its trivial
 *  to do it here.
 */
void Lgm_B_igrf_ctrans(Lgm_Vector *v, Lgm_Vector *B, Lgm_CTrans *c) {

    double  st, ct, sp, cp;
    double  B_r, B_theta, B_phi;
    double  r, theta, phi;
    Lgm_Vector  w, Bgeo;


    /*
     *  Convert Coords to GEO
     */
    Lgm_Convert_Coords( v, &w, GSM_TO_WGS84, c );


    /*
     *  compute GEO (WGS84) geocentric speherical coords. (r, theta, phi)  theta is colat
     */
    r     = sqrt(w.x*w.x + w.y*w.y + w.z*w.z );
    theta = acos( w.z / r );
    phi   = atan2(w.y, w.x);


    st    = sin( theta ); ct    = cos( theta );
    sp    = sin( phi );   cp    = cos( phi );


    /*
     *   Compute spherical components of B in GEO (WGS84)
     */
    w.x = r; w.y = theta; w.z = phi;
    Lgm_IGRF( &w, &Bgeo, c );
    B_r = Bgeo.x; B_theta = Bgeo.y; B_phi = Bgeo.z;



    /*
     *  Convert Bgeo to cartesian
     */
    Bgeo.x = B_r*st*cp + B_theta*ct*cp - B_phi*sp;
    Bgeo.y = B_r*st*sp + B_theta*ct*sp + B_phi*cp;
    Bgeo.z = B_r*ct    - B_theta*st;

    /*
     *  Convert Coords to GSM
     */
    Lgm_Convert_Coords( &Bgeo, B, WGS84_TO_GSM, c );


}

void Lgm_B_JensenCain1960_ctrans(Lgm_Vector *v, Lgm_Vector *B, Lgm_CTrans *c) {

    double  st, ct, sp, cp;
    double  B_r, B_theta, B_phi;
    double  r, theta, phi;
    Lgm_Vector  w, Bgeo;


    /*
     *  Convert Coords to GEO
     */
    Lgm_Convert_Coords( v, &w, GSM_TO_WGS84, c );


    /*
     *  compute GEO (WGS84) geocentric speherical coords. (r, theta, phi)  theta is colat
     */
    r     = sqrt(w.x*w.x + w.y*w.y + w.z*w.z );
    theta = acos( w.z / r );
    phi   = atan2(w.y, w.x);


    st    = sin( theta ); ct    = cos( theta );
    sp    = sin( phi );   cp    = cos( phi );


    /*
     *   Compute spherical components of B in GEO (WGS84)
     */
    w.x = r; w.y = theta; w.z = phi;
    Lgm_JensenCain1960( &w, &Bgeo, c );
    B_r = Bgeo.x; B_theta = Bgeo.y; B_phi = Bgeo.z;



    /*
     *  Convert Bgeo to cartesian
     */
    Bgeo.x = B_r*st*cp + B_theta*ct*cp - B_phi*sp;
    Bgeo.y = B_r*st*sp + B_theta*ct*sp + B_phi*cp;
    Bgeo.z = B_r*ct    - B_theta*st;

    /*
     *  Convert Coords to GSM
     */
    Lgm_Convert_Coords( &Bgeo, B, WGS84_TO_GSM, c );


}



void Lgm_B_cdip_ctrans(Lgm_Vector *v, Lgm_Vector *B, Lgm_CTrans *c) {

    double  x_sm, y_sm, z_sm;
    double  rho, theta, phi;
    Lgm_Vector  Bsm;
    double  B_rho, B_theta, B_phi;
    double  M, rho2, rho3;
    double  cp, sp, ct, st;

    M = c->M_cd;

    /*
     *  compute SM coords from GSM coords
     */
    x_sm = v->x*c->cos_psi - v->z*c->sin_psi;
    y_sm = v->y;
    z_sm = v->x*c->sin_psi + v->z*c->cos_psi;


    /*
     *  convert x_sm, y_sm, and, z_sm to spherical coords.
     */
    rho   = sqrt(x_sm*x_sm + y_sm*y_sm + z_sm*z_sm);
    phi   = atan2(y_sm, x_sm);
    theta = acos(z_sm / rho);


    /*
     *  compute centered dipole field in spherical coords
     *  i.e. (B_rho, B_theta, B_phi)
     */
    rho2     = rho*rho;
    rho3     = rho*rho2;
    cp       = cos( phi );
    sp       = sin( phi );
    ct       = cos( theta );
    st       = sin( theta );
    B_rho    = -2.0*M*ct/rho3;
    B_phi    =  0.0;
    B_theta  = -1.0*M*st/rho3;


    /*
     *   Transform (B_rho, B_theta, B_phi) -> (B_xsm, B_ysm, B_zsm)  (still SM)
     */
    Bsm.x = B_rho*st*cp + B_theta*ct*cp;
    Bsm.y = B_rho*st*sp + B_theta*ct*sp;
    Bsm.z = B_rho*ct    - B_theta*st;


    /*
     *  Transform (B_xsm, B_ysm, B_zsm) -> (B_x, B_y, B_z) i.e. trans. to GSM
     */
    B->x =  Bsm.x*c->cos_psi + Bsm.z*c->sin_psi;
    B->y =  Bsm.y;
    B->z = -Bsm.x*c->sin_psi + Bsm.z*c->cos_psi;


}


void Lgm_B_edip_ctrans(Lgm_Vector *v, Lgm_Vector *B, Lgm_CTrans *c) {

    double  x_sm, y_sm, z_sm;
    double  x_ed, y_ed, z_ed;
    double  rho, theta, phi;
    Lgm_Vector  Bsm;
    double  B_rho, B_theta, B_phi;
    double  M, rho2, rho3;
    double  cp, sp, ct, st;

    M = c->M_cd;

    /*
     *  compute SM coords from GSM coords
     */
    x_sm = v->x*c->cos_psi - v->z*c->sin_psi;
    y_sm = v->y;
    z_sm = v->x*c->sin_psi + v->z*c->cos_psi;


    /*
     *  compute ED coords from SM coords  (i.e. offset the dipole)
     */
    x_ed = x_sm - c->ED_x0;
    y_ed = y_sm - c->ED_y0;
    z_ed = z_sm - c->ED_z0;


    /*
     *  convert to spherical coords.
     */
    rho   = sqrt(x_ed*x_ed + y_ed*y_ed + z_ed*z_ed);
    phi   = atan2(y_ed, x_ed);
    theta = acos(z_ed / rho);


    /*
     *  compute centered dipole field in spherical coords
     *  i.e. (B_rho, B_theta, B_phi)
     */
    rho2     = rho*rho;
    rho3     = rho*rho2;
    cp       = cos( phi );
    sp       = sin( phi );
    ct       = cos( theta );
    st       = sin( theta );
    B_rho    = -2.0*M*ct/rho3;
    B_phi    =  0.0;
    B_theta  = -1.0*M*st/rho3;


    /*
     *   Transform (B_rho, B_theta, B_phi) -> (B_xsm, B_ysm, B_zsm)  (still SM)
     */
    Bsm.x = B_rho*st*cp + B_theta*ct*cp;
    Bsm.y = B_rho*st*sp + B_theta*ct*sp;
    Bsm.z = B_rho*ct    - B_theta*st;


    /*
     *  Transform (B_xsm, B_ysm, B_zsm) -> (B_x, B_y, B_z) i.e. trans. to GSM
     */
    B->x =  Bsm.x*c->cos_psi + Bsm.z*c->sin_psi;
    B->y =  Bsm.y;
    B->z = -Bsm.x*c->sin_psi + Bsm.z*c->cos_psi;


}


/*
 *  Compute the Julian Day number for the given date.
 *  Julian Date is the number of days since noon of Jan 1 4713 B.C.
 */
void Lgm_jd_to_ymdh ( double JD, long int *Date, int *year, int *month, int *day, double *UT ) {

    int     I, A, B, C, D, E, y, G;
    double  F, m, d;


    JD += 0.5;
    I = (int)JD;
    F = JD - (double)I;

    if ( I > 2299160 ) {
        A = (int)( ((double)I-1867216.25)/36524.25 );
        B = I + 1 + A - (int)( A/4 );
    } else {
        B = I;
    }

    C = B + 1524;
    D = (int)( ((double)C - 122.1)/365.25 );
    E = (int)( 365.25 * D );
    G = (int)( (double)((C-E))/30.6001 );
    d = (double)(C - E) + F - (double)( (int)( 30.6001*G ) );

    if ( (double)G < 13.5 ) {
        m = (double)G - 1.0;
    } else {
        m = (double)G - 13.0;
    }
    if ( m > 2.5 ) {
        y = D - 4716;
    } else {
        y = D - 4715;
    }

    *day = (int)d;
    *UT  = (d - *day)*24.0;
    *month = (int)m;
    *year = y;
    *Date = *year*10000 + *month*100 + *day;
}

//?
void Lgm_mjd_to_ymdh( double MJD, long int *Date, int *year, int *month, int *day, double *UT ) {
    Lgm_mjd_to_ymdh( MJD+2400000.5, Date, year, month, day, UT );
}



/*
 * milli-seconds are rounded up to nearest whole milli-second
 */
void Lgm_UT_to_hmsms( double UT, int *Hour, int *Min, int *Sec, int *MilliSec ) {

    *Hour = (int)UT; UT = (UT - *Hour)*60.0;
    *Min  = (int)UT; UT = (UT - *Min)*60.0;
    *Sec  = (int)UT; UT = (UT - *Sec)*1000.0;
    *MilliSec = (int)(UT+0.5);

    if (*MilliSec == 1000) {
        *MilliSec -= 1000;
        *Sec += 1;
    }

    if (*Sec == 60) {
        *Sec = *Sec-60;
        *Min += 1;
    }

    if ( *Min == 60 ) {
        *Min -= 60.0;
        *Hour += 1;
    }
}


/*
 * seconds are rounded up to nearest whole second
 */
void Lgm_UT_to_HMS( double UT, int *HH, int *MM, int *SS ) {

    int     Days, Hours, Minutes;
    double  Seconds;

    Seconds = UT*3600.0;

    Days     = (int)(Seconds/86400.0);
    Seconds -= Days*86400.0;

    Hours    = (int)(Seconds/3600.0);
    Seconds -= Hours*3600.0;

    Minutes  = (int)(Seconds/60.0);
    Seconds -= Minutes*60.0;

    *HH = Hours;
    *MM = Minutes;
    *SS = (int)(Seconds+0.5);

}

void Lgm_UT_to_HMSd( double UT, int *sgn, int *HH, int *MM, double *SS ) {

    int     Days, Hours, Minutes;
    double  Seconds;

    *sgn = (UT < 0.0) ? -1 : 1;
    UT = fabs(UT);

    Seconds = UT*3600.0;

    Days     = (int)(Seconds/86400.0);
    Seconds -= Days*86400.0;

    Hours    = (int)(Seconds/3600.0);
    Seconds -= Hours*3600.0;

    Minutes  = (int)(Seconds/60.0);
    Seconds -= Minutes*60.0;
//printf("Hours, Minutes, Seconds = %d %d %lf\n", Hours, Minutes, Seconds);

    if ( Seconds >= 60.0 ) {
        Seconds = 0.0;
        Minutes += 1;
    }
    if ( Minutes >= 60.0 ) {
        Minutes = 0;
        Hours += 1;
    }

    *HH = Hours;
    *MM = Minutes;
    *SS = Seconds;
//printf("Hours, Minutes, Seconds = %d %d %lf\n", Hours, Minutes, Seconds);

}

void Lgm_D_to_DMS( double D, int *DD, int *MM, int *SS ) {

    int     sgn;
    double  f;

    sgn = (D < 0.0) ? -1 : 1;
    D = fabs(D);

    *DD = (int)D; // degrees

    f = (D - (double)*DD)*60.0;  // decimal minutes
    *MM = (int)f;               // number of whole minutes

    f = (f - (double)*MM)*60.0; // decimal seconds
    *SS = (int)(f+0.5);

    if ( *SS == 60 ){
        *SS = 0;
        ++(*MM);
    }

    if ( *MM == 60 ){
        *MM = 0;
        ++(*DD);
    }




    *DD *= sgn;   // put sign back onto whole degrees part only

}


// return a decimal degree broken down int whole degrees, whole arc-minutes and
// decimal arc-seconds. Also return the sign of the original number.
void Lgm_D_to_DMSd( double D, int *sgn, int *DD, int *MM, double *SS ) {

    double  f;

    *sgn = (D < 0.0) ? -1 : 1;
    D = fabs(D);

    *DD = (int)D; // degrees

    f = (D - (double)*DD)*60.0;  // decimal minutes
    *MM = (int)f;               // number of whole minutes

    f = (f - (double)*MM)*60.0; // decimal seconds
    *SS = f;

}


double Lgm_GetCurrentJD( Lgm_CTrans *c ){

    time_t      CurrentTime;
    struct tm   tp;
    double      UTC, JD;


    CurrentTime = time( NULL );
    gmtime_r( &CurrentTime, &tp );
    UTC   = (double)(tp.tm_hour + tp.tm_min/60.0 +  tp.tm_sec/3600.0);
    JD   = Lgm_JD( tp.tm_year+1900, tp.tm_mon+1, tp.tm_mday, UTC, LGM_TIME_SYS_UTC, c );

    return( JD );

}

double Lgm_GetCurrentMJD( Lgm_CTrans *c ){
    return( Lgm_GetCurrentJD( c ) - 2400000.5 );
}


// types out hours/minutes/seconds with unicode characters.
void Lgm_Print_HMS( double d ){
    int HH, MM, SS;
    Lgm_UT_to_HMS( d, &HH, &MM, &SS );
    printf("%02d\u02b0 %02d\u1d50 %02d\u02e2", HH, MM, SS);
}
/*
 * This routine prints out time in different ways.
 *
 * int UnicodeHMS  -- Input
 *                    Controls the style of printing. If TRUE it
 *                    prints with h, m, s as superscripts.
 *                    Othwewise it prints with colons.
 * int p           -- Input
 *                    Controls how many decimals to print for the seconds.
 *                    Max of 20 (-- although prob more than machine prec.)
 */
void Lgm_Print_HMSdp( double d, int UnicodeHMS, int p ){
    int         HH, MM, SS, sgn;
    double      S, SFRAC;
    char        SecFracStr[30];


    Lgm_UT_to_HMSd( d, &sgn, &HH, &MM, &S );

//printf("HH, MM, S = %d %d %.20lf\n", HH, MM, S);
    S     = fabs(S)+pow(10.0, -1.0-p );
    SS    = (int)S;
    SFRAC = S-(double)SS;
//printf("HH, MM, S = %d %d %.20lf\n", HH, MM, S);


//printf("S = %lf SFRAC = %lf\n", S, SFRAC);
    if (( p <= 0 ) && ( SFRAC >= 0.5)) { SS += 1; }

//printf("S = %lf\n", S);

    sprintf( SecFracStr, "%.*lf", (p>=20)?20:p, SFRAC );
    if (SS==60)  { SS=0; ++MM; }
    if (MM==60)  { MM=0; ++HH; }


    // probably a cleaner way to do this(?)
    if ( p <= 0 ) {
        if ( UnicodeHMS ) {
            if (sgn<0) printf("-%02d\u02b0 %02d\u1d50 %02d\u02e2", HH, MM, SS );
            else       printf(" %02d\u02b0 %02d\u1d50 %02d\u02e2", HH, MM, SS );
        } else {
            if (sgn<0) printf("-%02d:%02d:%02d", HH, MM, SS );
            else       printf(" %02d:%02d:%02d", HH, MM, SS );
        }
    } else {
        if ( UnicodeHMS ) {
            if (sgn<0) printf("-%02d\u02b0 %02d\u1d50 %02d\u02e2.%s", HH, MM, SS, strstr(SecFracStr, ".")+1 );
            else       printf(" %02d\u02b0 %02d\u1d50 %02d\u02e2.%s", HH, MM, SS, strstr(SecFracStr, ".")+1 );
        } else {
            if (sgn<0) printf("-%02d:%02d:%02d.%s", HH, MM, SS, strstr(SecFracStr, ".")+1 );
            else       printf(" %02d:%02d:%02d.%s", HH, MM, SS, strstr(SecFracStr, ".")+1 );
        }
    }

}


// types out hours/minutes/seconds with unicode characters.
// does decimals on the seconds
void Lgm_Print_HMSd( double d ){
    int HH, MM, SS, MS, sgn;
    double S;
    Lgm_UT_to_HMSd( d, &sgn, &HH, &MM, &S );
    S = fabs(S);
    SS = (int)S;
    MS = (int)((S-(double)SS)*1000.0+0.5);
    if (MS==1000){ MS=0; ++SS; }
    if (SS==60)  { SS=0; ++MM; }
    if (MM==60)  { MM=0; ++HH; }
    if (sgn<0) printf("-%02d\u02b0 %02d\u1d50 %02d\u02e2.%03d", HH, MM, SS, MS);
    else       printf(" %02d\u02b0 %02d\u1d50 %02d\u02e2.%03d", HH, MM, SS, MS);

}

// types out degrees/arc-minutes/arc-seconds with unicode characters.
void Lgm_Print_DMS( double d ){
    int HH, MM, SS;
    Lgm_D_to_DMS( d, &HH, &MM, &SS );
    printf("%02d\u00b0 %02d\u2032 %02d\u2033", HH, MM, SS);
}

// types out degrees/arc-minutes/arc-seconds with unicode characters.
// does decimals on the arc-seconds
void Lgm_Print_DMSd( double d ){
    int DD, MM, SS, MS, sgn;
    double S;
    Lgm_D_to_DMSd( d, &sgn, &DD, &MM, &S );
    SS = (int)S;
    MS = (int)((S-(double)SS)*1000.0+0.5);
    if (MS==1000){ MS=0; ++SS; }
    if (SS==60)  { SS=0; ++MM; }
    if (MM==60)  { MM=0; ++DD; }

    if ( sgn < 0 ) printf("-%02d\u00b0 %02d\u2032 %02d\u2033.%03d", DD, MM, SS, MS);
    else           printf(" %02d\u00b0 %02d\u2032 %02d\u2033.%03d", DD, MM, SS, MS);
}

/*
 *  Given Cartesian CDMAG vector, convert to R, MLAT, MLON, MLT
 *  Note: MLT is defined to be the difference (in hours) between the magnetic
 *  longitude of the point in question and the magnetic longitude of the
 *  anti-solar point. Its coord system dependent.
 */
void Lgm_CDMAG_to_R_MLAT_MLON_MLT( Lgm_Vector *u, double *R, double *MLAT, double *MLON, double *MLT, Lgm_CTrans *c ) {

    double      SunMlon;
    Lgm_Vector  v, w;

    // Get R and MLAT
    w     = *u; // copy vector -- so we dont normalize original
    *R    = Lgm_NormalizeVector( &w );
    *MLAT = asin( u->z/(*R) )*DegPerRad;

    // Find Longitude of Sun vector in CDMAG coords.
    // And the Long thats 180deg. from it. (Because thats what MLT is reckoned from).
    Lgm_Convert_Coords( &(c->Sun), &v, MOD_TO_CDMAG, c );
    SunMlon = atan2( v.y, v.x )*DegPerRad; // in range -180 to 180
    SunMlon += 180.0; // in range 0 to 360

    // Find Longitude of input vector in CDMAG coords.
    *MLON = atan2( u->y, u->x )*DegPerRad; // in range -180 to 180
    if (*MLON < 0.0)  *MLON += 360.0;       // puts into range 0 to 360

    // Compute MLT (add 24 before the fmod so we dont get negative).
    *MLT = fmod( (*MLON-SunMlon)/15.0+24.0, 24.0 );

}


/*
 * returns cartesian vector in CDMAG coords.
 *  MLT in decimal hours.
 *  MLON in degrees.
 *  MLAT in degrees.
 *  R in Re
 *
 */
void Lgm_R_MLAT_MLT_to_CDMAG( double R, double MLAT, double MLT, Lgm_Vector *u, Lgm_CTrans *c ) {

    Lgm_Vector  v;
    double 	Mlon, SunMlon, phi, the, ct;

    // Find Longitude of Sun vector in CDMAG coords.
    // And the Long thats 180deg. from it. (Because thats what MLT is reckoned from).
    Lgm_Convert_Coords( &(c->Sun), &v, MOD_TO_CDMAG, c );
    SunMlon = atan2( v.y, v.x )*DegPerRad; // in range -180 to 180
    SunMlon += 180.0; // in range 0 to 360

    Mlon = (MLT*15.0 + SunMlon);
    phi = Mlon * RadPerDeg;
    the = MLAT*RadPerDeg;
    ct  = cos( the );


    u->x = R*cos( phi )*ct;
    u->y = R*sin( phi )*ct;
    u->z = R*sin( the );

}



/*
 *  Given Cartesian EDMAG vector, convert to R, MLAT, MLON, MLT
 */
void Lgm_EDMAG_to_R_MLAT_MLON_MLT( Lgm_Vector *u, double *R, double *MLAT, double *MLON, double *MLT, Lgm_CTrans *c ) {

    double      SunMlon;
    Lgm_Vector  v, w;

    // Get R and MLAT
    w     = *u; // copy vector -- so we dont normalize original
    *R    = Lgm_NormalizeVector( &w );
    *MLAT = asin( u->z/(*R) )*DegPerRad;

    // Find Longitude of Sun vector in EDMAG coords.
    // And the Long thats 180deg. from it. (Because thats what MLT is reckoned from).
    Lgm_Convert_Coords( &(c->Sun), &v, MOD_TO_EDMAG, c );
    SunMlon = atan2( v.y, v.x )*DegPerRad; // in range -180 to 180
    SunMlon += 180.0; // in range 0 to 360

    // Find Longitude of input vector in EDMAG coords.
    *MLON = atan2( u->y, u->x )*DegPerRad; // in range -180 to 180
    if (*MLON < 0.0)  *MLON += 360.0;       // puts into range 0 to 360

    // Compute MLT (add 24 before the fmod so we dont get negative).
    *MLT = fmod( (*MLON-SunMlon)/15.0+24.0, 24.0 );

}

/*
 * returns cartesian vector in EDMAG coords.
 *  	R in Re.
 *  	MLAT in degrees.
 *  	MLT in decimal hours.
 */
void Lgm_R_MLAT_MLT_to_EDMAG( double R, double MLAT, double MLT, Lgm_Vector *u, Lgm_CTrans *c ) {

    Lgm_Vector  v;
    double 	Mlon, SunMlon, phi, the, ct;

    // Find Longitude of Sun vector in EDMAG coords.
    // And the Long thats 180deg. from it. (Because thats what MLT is reckoned from).
    Lgm_Convert_Coords( &(c->Sun), &v, MOD_TO_EDMAG, c );
    SunMlon = atan2( v.y, v.x )*DegPerRad; // in range -180 to 180
    SunMlon += 180.0; // in range 0 to 360

    Mlon = (MLT*15.0 + SunMlon);
    phi = Mlon * RadPerDeg;
    the = MLAT*RadPerDeg;
    ct  = cos( the );


    u->x = R*cos( phi )*ct;
    u->y = R*sin( phi )*ct;
    u->z = R*sin( the );

}



/*
 *  Convert GLAT/GLON to CDMAG MLAT/MLON/MLT
 */
void Lgm_GLATLON_TO_CDMLATLONMLT( double GLAT, double GLON, double *MLAT, double *MLON, double *MLT, Lgm_CTrans *c ) {

    Lgm_Vector  v, w, z;
    double 	Mlon, SunMlon, phi, the, ct;

    // Find Longitude of Sun vector in CDMAG coords.
    // And the Long thats 180deg. from it. (Because thats what MLT is reckoned from).
    Lgm_Convert_Coords( &(c->Sun), &v, MOD_TO_CDMAG, c );
    SunMlon = atan2( v.y, v.x )*DegPerRad; // in range -180 to 180
    SunMlon += 180.0; // in range 0 to 360

    // Convert GLAT/GLON to cartesian (assume r==1)
    phi = GLON*RadPerDeg;
    the = GLAT*RadPerDeg;
    ct  = cos(the);
    w.x = cos( phi )*ct;
    w.y = sin( phi )*ct;
    w.z = sin( the );

    // Convert to CDMAG
    Lgm_Convert_Coords( &w, &z, WGS84_TO_CDMAG, c );
    *MLAT = asin( z.z )*DegPerRad;
    *MLON = atan2( z.y, z.x )*DegPerRad;
    if (*MLON < 0.0)  *MLON += 360.0;       // puts into range 0 to 360

    // Compute MLT (add 24 before the fmod so we dont get negative).
    *MLT = fmod( (*MLON-SunMlon)/15.0+24.0, 24.0 );

}

/*
 *  Convert GLAT/GLON to EDMAG MLAT/MLON/MLT
 */
void Lgm_GLATLON_TO_EDMLATLONMLT( double GLAT, double GLON, double *MLAT, double *MLON, double *MLT, Lgm_CTrans *c ) {

    Lgm_Vector  v, w, z;
    double 	Mlon, SunMlon, phi, the, ct;

    // Find Longitude of Sun vector in CDMAG coords.
    // And the Long thats 180deg. from it. (Because thats what MLT is reckoned from).
    Lgm_Convert_Coords( &(c->Sun), &v, MOD_TO_EDMAG, c );
    SunMlon = atan2( v.y, v.x )*DegPerRad; // in range -180 to 180
    SunMlon += 180.0; // in range 0 to 360

    // Convert GLAT/GLON to cartesian (assume r==1)
    phi = GLON*RadPerDeg;
    the = GLAT*RadPerDeg;
    ct  = cos(the);
    w.x = cos( phi )*ct;
    w.y = sin( phi )*ct;
    w.z = sin( the );

    // Convert to CDMAG
    Lgm_Convert_Coords( &w, &z, WGS84_TO_EDMAG, c );
    *MLAT = asin( z.z )*DegPerRad;
    *MLON = atan2( z.y, z.x )*DegPerRad;
    if (*MLON < 0.0)  *MLON += 360.0;       // puts into range 0 to 360

    // Compute MLT (add 24 before the fmod so we dont get negative).
    *MLT = fmod( (*MLON-SunMlon)/15.0+24.0, 24.0 );

}


int MonthStrToNum( char *str ) {
    char *p = Lgm_StrToLower( str, strlen(str) ); // force to lower case
    if      ( !strncmp( p, "jan", 3 ) ) { return(  1 ); }
    else if ( !strncmp( p, "feb", 3 ) ) { return(  2 ); }
    else if ( !strncmp( p, "mar", 3 ) ) { return(  3 ); }
    else if ( !strncmp( p, "apr", 3 ) ) { return(  4 ); }
    else if ( !strncmp( p, "may", 3 ) ) { return(  5 ); }
    else if ( !strncmp( p, "jun", 3 ) ) { return(  6 ); }
    else if ( !strncmp( p, "jul", 3 ) ) { return(  7 ); }
    else if ( !strncmp( p, "aug", 3 ) ) { return(  8 ); }
    else if ( !strncmp( p, "sep", 3 ) ) { return(  9 ); }
    else if ( !strncmp( p, "oct", 3 ) ) { return( 10 ); }
    else if ( !strncmp( p, "nov", 3 ) ) { return( 11 ); }
    else if ( !strncmp( p, "dec", 3 ) ) { return( 12 ); }
    else                                { return( LGM_ERROR ); }
}


char *Lgm_StrToLower( char *str, int nmax ) {
    int  n = 0;
    char *p = str;
    while ( (p != '\0' ) && (n<nmax) ){
	    *p = (unsigned char)tolower(*p);
	    ++p;
	    ++n;
	}
    return( str );
}

char *Lgm_StrToUpper( char *str, int nmax ) {
    int  n = 0;
    char *p = str;
    while ( (p != '\0' ) && (n<nmax) ){
	    *p = (unsigned char)toupper(*p);
	    ++p;
	    ++n;
	}
    return( str );
}
