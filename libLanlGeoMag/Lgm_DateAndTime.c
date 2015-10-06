/*! \file Lgm_DateAndTime.c
 *  
 *  \brief  Collection of routines for data and time calculations. (Also see rotuines in Lgm_CTrans.c).
 *
 */
#ifdef HAVE_CONFIG_H 
#include "config.h" 
#endif
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Lgm/Lgm_CTrans.h"





/*
 * returns 1 if the year is a leap year, 0 otherwise.
 */
int Lgm_LeapYear(int year) {

    if ((year%100 == 0)&&(year%400 != 0)) return(0);
    else if (year%4 == 0) return(1);
    else return(0);

}


/*!
 * Routines for dealing with leap seconds. Also time conversions that require
 * knowledge about leap seconds.
 *
 *
 * Lgm_GetLeapSeconds()
 * Lgm_IsLeapSecondDay()
 * Lgm_LoadLeapSeconds()
 *
 *
 * Leap seconds are added when necessary. First preference is given to
 * opportunities at the end of December and the end of June. However, secondary
 * preference is also given to opportunities at the end of March and September
 * if needed. Since leap seconds were introduced in 1972, so far only dates in
 * December and June have been used to include leap seconds.
 *
 *
 *  Needed to convert between TT or TAI, TDB and UTC. Some defs:
 *
 *      JD  -- Julian Date
 *      MJD -- Modified Julian Date (JD - 2400000.5)
 *      UT  -- Universal Time (before 1960 astro calcs were done with UT)
 *      ET  -- Ephemeris Time (then ET replaced it)
 *      TDT -- Terrestrial Dynamical Time (TDT replaced ET in 1981)
 *      TT  -- Terrestrial Time (in 1991 TDT was renamed to be TT)
 *      UTC -- Universal Time Coordinated
 *      TAI -- International Atomic Time
 *      TDB -- Terrestrial Barycentric Time
 *      UT1 -- Universal Time (UT1 is a corrected version of UT0)
 *      UT0 -- Universal Time (not used in this uncorrected form)
 *      dT  -- different between TT and UT1 (i.e. dT = TT-UT1)
 *             (this was 32.184s when TAI was introduced in 1958 -- hence the
 *              definitions below).
 *             
 *
 * Related by:
 *
 *      TAI = UTC + dAT
 *      TT = TAI + 32.184 (i.e.  TT = UTC + dAT + 32.184 )
 *
 * This routine simply determines what the dAT value should be for a given JD.
 * Note that after 1972 they are just leap seconds, but before they are
 * non-integral (leap seconds (werent invented yet?)
 *
 */
double Lgm_GetLeapSeconds( double JD, Lgm_CTrans *c ) {

    Lgm_LeapSeconds *l;
    int             i;
    double          MJD, TAI_Minus_UTC = 0.0;

    l = &(c->l);

    /*
     * For dates 1972 to present, leap seconds are stored in the
     * Lgm_LeapSeconds structure.
     * Parse list in reverse order.
     */
    for (i=l->nLeapSecondDates-1; i>=0; i--) {
        // do > and not >= here because the 86400th second is still before the leap second
        if ( JD >= l->LeapSecondJDs[i] ) {
            TAI_Minus_UTC = l->LeapSeconds[i];
            return( TAI_Minus_UTC );
        }
    }

    /*
     * Date is pre 1972.
     */
    // Different system (leap seconds introduced in 1972). Calculate MJD first.
    MJD = JD - 2400000.5;
    if        ( JD >= 2439887.5 ) TAI_Minus_UTC = 4.2131700 + (MJD-39126.0)*0.002592;  // 1968 FEB 1 =JD 2439887.5 TAI-UTC= 4.2131700 S + (MJD - 39126.) X 0.002592  S
    else if   ( JD >= 2439126.5 ) TAI_Minus_UTC = 4.3131700 + (MJD-39126.0)*0.002592;  // 1966 JAN 1 =JD 2439126.5 TAI-UTC= 4.3131700 S + (MJD - 39126.) X 0.002592  S
    else if   ( JD >= 2439004.5 ) TAI_Minus_UTC = 3.8401300 + (MJD-38761.0)*0.002592;  // 1965 SEP 1 =JD 2439004.5 TAI-UTC= 3.8401300 S + (MJD - 38761.) X 0.001296  S
    else if   ( JD >= 2438820.5 ) TAI_Minus_UTC = 3.6401300 + (MJD-38761.0)*0.001296;  // 1965 MAR 1 =JD 2438820.5 TAI-UTC= 3.6401300 S + (MJD - 38761.) X 0.001296  S
    else if   ( JD >= 2438942.5 ) TAI_Minus_UTC = 3.7401300 + (MJD-38761.0)*0.001296;  // 1965 JUL 1 =JD 2438942.5 TAI-UTC= 3.7401300 S + (MJD - 38761.) X 0.001296  S
    else if   ( JD >= 2438761.5 ) TAI_Minus_UTC = 3.5401300 + (MJD-38761.0)*0.001296;  // 1965 JAN 1 =JD 2438761.5 TAI-UTC= 3.5401300 S + (MJD - 38761.) X 0.001296  S
    else if   ( JD >= 2438639.5 ) TAI_Minus_UTC = 3.4401300 + (MJD-38761.0)*0.001296;  // 1964 SEP 1 =JD 2438639.5 TAI-UTC= 3.4401300 S + (MJD - 38761.) X 0.001296  S
    else if   ( JD >= 2438395.5 ) TAI_Minus_UTC = 3.2401300 + (MJD-38761.0)*0.001296;  // 1964 JAN 1 =JD 2438395.5 TAI-UTC= 3.2401300 S + (MJD - 38761.) X 0.001296  S
    else if   ( JD >= 2438486.5 ) TAI_Minus_UTC = 3.3401300 + (MJD-38761.0)*0.001296;  // 1964 APR 1 =JD 2438486.5 TAI-UTC= 3.3401300 S + (MJD - 38761.) X 0.001296  S
    else if   ( JD >= 2438334.5 ) TAI_Minus_UTC = 1.9458580 + (MJD-37665.0)*0.0011232; // 1963 NOV 1 =JD 2438334.5 TAI-UTC= 1.9458580 S + (MJD - 37665.) X 0.0011232 S
    else if   ( JD >= 2437665.5 ) TAI_Minus_UTC = 1.8458580 + (MJD-37665.0)*0.0011232; // 1962 JAN 1 =JD 2437665.5 TAI-UTC= 1.8458580 S + (MJD - 37665.) X 0.0011232 S
    else if   ( JD >= 2437300.5 ) TAI_Minus_UTC = 1.4228180 + (MJD-37300.0)*0.001296;  // 1961 JAN 1 =JD 2437300.5 TAI-UTC= 1.4228180 S + (MJD - 37300.) X 0.001296  S
    else if   ( JD >= 2437512.5 ) TAI_Minus_UTC = 1.3728180 + (MJD-37300.0)*0.001296;  // 1961 AUG 1 =JD 2437512.5 TAI-UTC= 1.3728180 S + (MJD - 37300.) X 0.001296  S

    return( TAI_Minus_UTC );

}




/*
 * Returns a 1 if the given date is a LeapSecond date, 0 otherwise. I.e. -- a
 * date on which a leap second was added.
 */
int Lgm_IsLeapSecondDay( long int Date, double *SecondsInDay, Lgm_CTrans *c ) {

    Lgm_LeapSeconds *l;
    int             i;

    l = &(c->l);

    /*
     * Parse list in reverse order.
     */
    *SecondsInDay = 86400.0;
    for (i=l->nLeapSecondDates-1; i>=0; i--) {
        if ( Date == l->LeapSecondDates[i] ) {
            if (i > 0) {
                if ( (l->LeapSeconds[i] - l->LeapSeconds[i-1]) > 0.0) {
                    *SecondsInDay = 86401.0;
                } else {
                    *SecondsInDay = 86399.0; // implies leap second deleted (this has never happened yet).
                }
            }
            return(1); 
        }
    }
    return(0);

}




/*
 * Reads in the file containing leap second info and packs the results into a
 * Lgm_LeapSeconds Structure contained in the Lgm_CTrans structure.
 */
int Lgm_LoadLeapSeconds( Lgm_CTrans  *c ) {

    Lgm_LeapSeconds *l;
    int              i, n=0, N=50, Year, Month, Day;
    char             Line[513], LeapSecondFile[512];
    double           JD, Time;
    FILE             *fp;

    l = &(c->l);

    sprintf( LeapSecondFile, "%s/%s", LGM_EOP_DATA_DIR, "/Lgm_LeapSecondDates.dat");
    //printf("File = %s\n", LeapSecondFile );
    
    if ( (fp = fopen( LeapSecondFile, "r" )) != NULL ){

        l->nLeapSecondDates = 0;
        l->LeapSecondDates  = (long int *) malloc( N*sizeof(long int) );
        l->LeapSecondJDs    = (double *)   malloc( N*sizeof(double) );
        l->LeapSeconds      = (double *)   malloc( N*sizeof(double) );
        
        // read off header
        for (i=0;i<17; i++) fgets( Line, 512, fp );

        // read off data
        while ( fscanf( fp, "%ld %lf %lf", &l->LeapSecondDates[n], &l->LeapSecondJDs[n], &l->LeapSeconds[n] ) != EOF) { 

            //printf("l->LeapSecondDates[%d] = %ld, l->LeapSecondJDs[%d] = %lf, l->LeapSeconds[%d] = %g\n", 
            //        n, l->LeapSecondDates[n], n, l->LeapSecondJDs[n], n, l->LeapSeconds[n] );
            ++n; 

            // increase number of elements if needed
            if ( n>(N-1) ){
                N = (int)(N*1.2);
                if ( (l->LeapSecondDates = (long int *)realloc( (void *)l->LeapSecondDates, N*sizeof(long int) )) == NULL ) {
                    printf("Lgm_LoadLeapSeconds: Memory allocation problem (in realloc)\n");
                    fclose(fp);
                    return(LGM_ERROR);
                }
                if ( (l->LeapSecondJDs = (double *)realloc( (void *)l->LeapSecondJDs, N*sizeof(double) )) == NULL ) {
                    printf("Lgm_LoadLeapSeconds: Memory allocation problem (in realloc)\n");
                    fclose(fp);
                    return(LGM_ERROR);
                }
                if ( (l->LeapSeconds = (double *)realloc( (void *)l->LeapSeconds, N*sizeof(double) )) == NULL ) {
                    printf("Lgm_LoadLeapSeconds: Memory allocation problem (in realloc)\n");
                    fclose(fp);
                    return(LGM_ERROR);
                }
            }
        }
        fclose(fp);

    } else {
        /*
         * Provides a fallback in case the Lgm_LeapSecondDates.dat was not present.
         */
        N = 27;
        l->nLeapSecondDates = N;
        l->LeapSecondDates  = (long int *) malloc( N*sizeof(long int) );
        l->LeapSecondJDs    = (double *)   malloc( N*sizeof(double) );
        l->LeapSeconds      = (double *)   malloc( N*sizeof(double) );
        l->LeapSecondDates[0]  = 19720101, l->LeapSecondJDs[0]  = 2441317.5, l->LeapSeconds[0]  = 10.0;
        l->LeapSecondDates[1]  = 19720701, l->LeapSecondJDs[1]  = 2441499.5, l->LeapSeconds[1]  = 11.0;
        l->LeapSecondDates[2]  = 19730101, l->LeapSecondJDs[2]  = 2441683.5, l->LeapSeconds[2]  = 12.0;
        l->LeapSecondDates[3]  = 19740101, l->LeapSecondJDs[3]  = 2442048.5, l->LeapSeconds[3]  = 13.0;
        l->LeapSecondDates[4]  = 19750101, l->LeapSecondJDs[4]  = 2442413.5, l->LeapSeconds[4]  = 14.0;
        l->LeapSecondDates[5]  = 19760101, l->LeapSecondJDs[5]  = 2442778.5, l->LeapSeconds[5]  = 15.0;
        l->LeapSecondDates[6]  = 19770101, l->LeapSecondJDs[6]  = 2443144.5, l->LeapSeconds[6]  = 16.0;
        l->LeapSecondDates[7]  = 19780101, l->LeapSecondJDs[7]  = 2443509.5, l->LeapSeconds[7]  = 17.0;
        l->LeapSecondDates[8]  = 19790101, l->LeapSecondJDs[8]  = 2443874.5, l->LeapSeconds[8]  = 18.0;
        l->LeapSecondDates[9]  = 19800101, l->LeapSecondJDs[9]  = 2444239.5, l->LeapSeconds[9]  = 19.0;
        l->LeapSecondDates[10] = 19810701, l->LeapSecondJDs[10] = 2444786.5, l->LeapSeconds[10] = 20.0;
        l->LeapSecondDates[11] = 19820701, l->LeapSecondJDs[11] = 2445151.5, l->LeapSeconds[11] = 21.0;
        l->LeapSecondDates[12] = 19830701, l->LeapSecondJDs[12] = 2445516.5, l->LeapSeconds[12] = 22.0;
        l->LeapSecondDates[13] = 19850701, l->LeapSecondJDs[13] = 2446247.5, l->LeapSeconds[13] = 23.0;
        l->LeapSecondDates[14] = 19880101, l->LeapSecondJDs[14] = 2447161.5, l->LeapSeconds[14] = 24.0;
        l->LeapSecondDates[15] = 19900101, l->LeapSecondJDs[15] = 2447892.5, l->LeapSeconds[15] = 25.0;
        l->LeapSecondDates[16] = 19910101, l->LeapSecondJDs[16] = 2448257.5, l->LeapSeconds[16] = 26.0;
        l->LeapSecondDates[17] = 19920701, l->LeapSecondJDs[17] = 2448804.5, l->LeapSeconds[17] = 27.0;
        l->LeapSecondDates[18] = 19930701, l->LeapSecondJDs[18] = 2449169.5, l->LeapSeconds[18] = 28.0;
        l->LeapSecondDates[19] = 19940701, l->LeapSecondJDs[19] = 2449534.5, l->LeapSeconds[19] = 29.0;
        l->LeapSecondDates[20] = 19960101, l->LeapSecondJDs[20] = 2450083.5, l->LeapSeconds[20] = 30.0;
        l->LeapSecondDates[21] = 19970701, l->LeapSecondJDs[21] = 2450630.5, l->LeapSeconds[21] = 31.0;
        l->LeapSecondDates[22] = 19990101, l->LeapSecondJDs[22] = 2451179.5, l->LeapSeconds[22] = 32.0;
        l->LeapSecondDates[23] = 20060101, l->LeapSecondJDs[23] = 2453736.5, l->LeapSeconds[23] = 33.0;
        l->LeapSecondDates[24] = 20090101, l->LeapSecondJDs[24] = 2454832.5, l->LeapSeconds[24] = 34.0;
        l->LeapSecondDates[25] = 20120701, l->LeapSecondJDs[25] = 2456109.5, l->LeapSeconds[25] = 35.0;
        l->LeapSecondDates[26] = 20150701, l->LeapSecondJDs[26] = 2457204.5, l->LeapSeconds[26] = 36.0;
        printf("Lgm_LoadLeapSeconds: Could not open Lgm_LeapSecondDates.dat file!\n");
        printf("                     Setting the leap second values that I know about\n");
        printf("                     (latest leap second I know about was introduced on\n");
        printf("                     %8ld (total leap seconds is %g))\n", l->LeapSecondDates[N-1], l->LeapSeconds[N-1]);
        //return(LGM_ERROR);
    }

    l->nLeapSecondDates = n;

    /*
     * The LeapSecondDates values are the dates that the given corrections came into effect. But the actual leap second
     * date is the day before -- thats when the extra second was tacked on. So subtract a day from each of the dates.
     */
    for (i=0; i<l->nLeapSecondDates; i++){
        JD = l->LeapSecondJDs[i]+0.5; // Noon on the given date
        Lgm_JD_to_Date( JD-1.0, &Year, &Month, &Day, &Time );
        l->LeapSecondDates[i] = Year*10000 + Month*100 + Day;
    }

    return(TRUE);

}




/*
 * GPS <-> TAI Conversions
 */
void Lgm_TAI_to_GPS( Lgm_DateTime *TAI, Lgm_DateTime *GPS, Lgm_CTrans *c ) {

    double  Time;

    *GPS = *TAI;
    GPS->Time -= 19.0/3600.0;
    GPS->JD   = Lgm_JD( GPS->Year, GPS->Month, GPS->Day, GPS->Time, LGM_TIME_SYS_GPS, c );
    GPS->Date = Lgm_JD_to_Date( GPS->JD, &GPS->Year, &GPS->Month, &GPS->Day, &Time );
    // Keep the time we had rather than getting it back from Lgm_JD_to_Date(). This limits roundoff error
    GPS->Time = Lgm_hour24( GPS->Time ); 
    GPS->T    = (GPS->JD - 2451545.0)/36525.0;
    GPS->TimeSystem = LGM_TIME_SYS_GPS;

}

void Lgm_GPS_to_TAI( Lgm_DateTime *GPS, Lgm_DateTime *TAI, Lgm_CTrans *c ) {

    double  Time;

    *TAI = *GPS;
    TAI->Time += 19.0/3600.0;
    TAI->JD   = Lgm_JD( TAI->Year, TAI->Month, TAI->Day, TAI->Time, LGM_TIME_SYS_TAI, c );
    TAI->Date = Lgm_JD_to_Date( TAI->JD, &TAI->Year, &TAI->Month, &TAI->Day, &Time );
    // Keep the time we had rather than getting it back from Lgm_JD_to_Date(). This limits roundoff error
    TAI->Time = Lgm_hour24( TAI->Time ); 
    TAI->T    = (TAI->JD - 2451545.0)/36525.0;
    TAI->TimeSystem = LGM_TIME_SYS_TAI;

}




/*
 * UTC <-> GPS conversions
 */
void Lgm_UTC_to_GPS( Lgm_DateTime *UTC, Lgm_DateTime *GPS, Lgm_CTrans *c ) {

    Lgm_DateTime    TAI;

    // First go UTC -> TAI
    Lgm_UTC_to_TAI( UTC, &TAI, c );

    // Then go TAI->GPS
    Lgm_TAI_to_GPS( &TAI, GPS, c );

    GPS->TimeSystem = LGM_TIME_SYS_GPS;

}

void Lgm_GPS_to_UTC( Lgm_DateTime *GPS, Lgm_DateTime *UTC, Lgm_CTrans *c ) {

    Lgm_DateTime    TAI;

    // First go GPS -> TAI
    Lgm_GPS_to_TAI( GPS, &TAI, c );

    // Then go TAI->UTC
    Lgm_TAI_to_UTC( &TAI, UTC, c );

    UTC->TimeSystem = LGM_TIME_SYS_UTC;

}





/*
 * GpsSeconds Routines
 * -------------------
 * Routines to compute number of SI seconds since 0h Jan 6, 1980 Some refer to
 * this as "GPS Time" but this is not a standard definition. (E.g., GPS Time is
 * often measured in whole weeks since 0h Jan 6, 1980 plus seconds into the
 * current week. But the week number will eventually roll over sometime soon?
 * The Weeks/Seconds thing is kind of useless for our purposes...)
 */
double  Lgm_GPS_to_GpsSeconds( Lgm_DateTime *GPS ) {
    double  JDN, Seconds;
    JDN  = Lgm_JDN( GPS->Year, GPS->Month, GPS->Day );
    Seconds = (JDN - LGM_JD_GPS0)*86400.0 + GPS->Time*3600.0;
    return( Seconds );
}

void Lgm_GpsSeconds_to_GPS( double GpsSeconds, Lgm_DateTime *GPS ) {

    long int dJDN;
    double   tmp, JDN;
    dJDN = (long int)(GpsSeconds/86400.0);
    GPS->Time = (GpsSeconds - dJDN*86400.0)/3600.0;
    JDN = LGM_JD_GPS0 + dJDN;
    Lgm_jd_to_ymdh( JDN, &GPS->Date, &GPS->Year, &GPS->Month, &GPS->Day, &tmp);
    
    GPS->JD   = JDN + GPS->Time/24.0;
    GPS->T    = (GPS->JD - 2451545.0)/36525.0;
    GPS->TimeSystem = LGM_TIME_SYS_GPS;

}

void Lgm_GpsSeconds_to_UTC( double GpsSeconds, Lgm_DateTime *UTC, Lgm_CTrans *c ) {
    Lgm_DateTime GPS;
    Lgm_GpsSeconds_to_GPS( GpsSeconds, &GPS );
    Lgm_GPS_to_UTC( &GPS, UTC, c );
    UTC->TimeSystem = LGM_TIME_SYS_UTC;
}

double Lgm_UTC_to_GpsSeconds( Lgm_DateTime *UTC, Lgm_CTrans *c ) {
    Lgm_DateTime GPS;
    Lgm_UTC_to_GPS( UTC, &GPS, c );
    return( Lgm_GPS_to_GpsSeconds( &GPS) );
}


/*
 * TaiSeconds Routines
 * -------------------
 * Routines to compute number of SI seconds since 0h Jan 1, 1958. Some refer to
 * this as "TAI", but this (IMO) is a distortion of what TAI means. TAI is
 * really defined in terms of how many seconds it differs from UTC. But some
 * have essentially re-named TAI to mean number of seconds since 0h Jan 1,
 * 1958.  That date is releveant, but it is actually just meant to be the date
 * when the difference between TAI and UTC was zero.
 *
 * Here, I call the "number of SI seconds since 0h Jan 1, 1958" TaiSeconds.
 */
double  Lgm_TAI_to_TaiSeconds( Lgm_DateTime *TAI ) {
    double  JDN, Seconds;
    JDN  = Lgm_JDN( TAI->Year, TAI->Month, TAI->Day );
    Seconds = (JDN - LGM_JD_TAI0)*86400.0 + TAI->Time*3600.0;
    return( Seconds );
}

void Lgm_TaiSeconds_to_TAI( double TaiSeconds, Lgm_DateTime *TAI ) {

    long int dJDN;
    double   tmp, JDN;
    dJDN = (long int)(TaiSeconds/86400.0);
    TAI->Time = (TaiSeconds - dJDN*86400.0)/3600.0;
    JDN = LGM_JD_TAI0 + dJDN;
    Lgm_jd_to_ymdh( JDN, &TAI->Date, &TAI->Year, &TAI->Month, &TAI->Day, &tmp);
    
    TAI->JD   = JDN + TAI->Time/24.0;
    TAI->T    = (TAI->JD - 2451545.0)/36525.0;
    TAI->TimeSystem = LGM_TIME_SYS_TAI;

}

void Lgm_TaiSeconds_to_UTC( double TaiSeconds, Lgm_DateTime *UTC, Lgm_CTrans *c ) {
    Lgm_DateTime TAI;
    Lgm_TaiSeconds_to_TAI( TaiSeconds, &TAI );
    Lgm_TAI_to_UTC( &TAI, UTC, c );
    UTC->TimeSystem = LGM_TIME_SYS_UTC;
}

double Lgm_UTC_to_TaiSeconds( Lgm_DateTime *UTC, Lgm_CTrans *c ) {
    Lgm_DateTime TAI;
    Lgm_UTC_to_TAI( UTC, &TAI, c );
    return( Lgm_TAI_to_TaiSeconds( &TAI) );
}





/*
 * UTC <-> TAI conversion rotuines
 *
 * Lgm_UTC_to_TAI() - converts UTC to TAI using appropriaste number of leapseconds
 * Lgm_TAI_to_UTC() - converts TAI to UTC using appropriaste number of leapseconds
 *
 */
void Lgm_UTC_to_TAI( Lgm_DateTime *UTC, Lgm_DateTime *TAI, Lgm_CTrans *c ) {

    double  Time, DAT, JD;
    int     Year, Month, Day, sgn;

    // Set relavent fields in TAI structure.
    *TAI = *UTC;
    TAI->TimeSystem = LGM_TIME_SYS_TAI;
    TAI->DaySeconds = 86400.0;

    // Get leap seconds
    JD = Lgm_JD( UTC->Year, UTC->Month, UTC->Day, UTC->Time, LGM_TIME_SYS_UTC, c );
    DAT = Lgm_GetLeapSeconds( JD, c );

    // Add leap seconds to UTC time to get TAI time
    TAI->Time = UTC->Time + DAT/3600.0;

    // Convert Year,Month,Day,Time TAI into a Julian Date.(Note that the Time
    // argument to Lgm_JD() doesnt have to be between 0 and 24. This allows for
    // date roll-overs one way or the other.)
    TAI->JD  = Lgm_JD( TAI->Year, TAI->Month, TAI->Day, TAI->Time, LGM_TIME_SYS_TAI, c );

    // From the new JD, compute the Gregorian Date that the TAI Time occurred on.
    // This routin returns Time back out, but we will ignore it.
    TAI->Date = Lgm_JD_to_Date( TAI->JD, &TAI->Year, &TAI->Month, &TAI->Day, &Time );

    // Instead, just remap the time we have to be between 0 and 24 (limits roundoff error)
    TAI->Time = Lgm_hour24( TAI->Time );

    // Compute Julian Centuries since J2000.
    TAI->T = (TAI->JD - 2451545.0)/36525.0;

    // Set Day Of Week and Day Of Week String
    TAI->Dow = Lgm_DayOfWeek( TAI->Year, TAI->Month, TAI->Day, TAI->DowStr );

    // Set Hour, Minute, Seconds
    Lgm_UT_to_HMSd( TAI->Time, &sgn, &TAI->Hour, &TAI->Minute, &TAI->Second );

    // Set Day Of Year
    Lgm_Doy( TAI->Date, &Year, &Month, &Day, &TAI->Doy);

    TAI->TimeSystem = LGM_TIME_SYS_TAI;

    return;
}

/*
 *  TAI -> UTC 
 *  This rotuine is not too complicated, but its not trivial. Here is an
 *  outline of the problem:
 *
 *  Let,
 *
 *      u = UTC time
 *      t = TAI time
 *      f = Number of leap seconds. Its a function of u.
 *
 *  Then to convert u->t, we simply compute;
 *
 *      t = u + f(u)
 *
 *  This is straight forward and gives us no problems. However, if we need to go
 *  t->u, the we have to solve,
 *
 *      u = t - f(u)
 *
 *  The problem is that we need f(u) before we know what u is. Most of the time
 *  this is no problem because we can just assume the dates are the same. But
 *  if t->u maps a TAI time back across a leapsecond (from one date to a
 *  previous date), there will be an inconsistency in f -- we will have used a
 *  value for f that was too large.  We need to detect these cases and adjust
 *  as necessary...
 *
 *
 */
void Lgm_TAI_to_UTC( Lgm_DateTime *TAI, Lgm_DateTime *UTC, Lgm_CTrans *c ) {

    double          t_target, f, d, tim;
    Lgm_DateTime    t;

    /*
     * Keep track of the TAI time we want to convert.
     */
    t_target = TAI->Time;

    /*
     * estimate of JD for UTC
     */
    d = TAI->JD;
//printf("d = %lf\n", d);
//printf("%d %d %d\n", TAI->Year, TAI->Month, TAI->Day);
//printf("JD = %lf\n", Lgm_JD( TAI->Year, TAI->Month, TAI->Day, TAI->Time, LGM_TIME_SYS_TAI, c ));

    /* 
     * Estimate UTC using the JD of the TAI time. This is a naive first guess
     * at the answer.
     */
    f = Lgm_GetLeapSeconds( d, c );
    UTC->Time  = TAI->Time - f/3600.0; 
    UTC->Year  = TAI->Year; UTC->Month = TAI->Month; UTC->Day = TAI->Day;

    /*
     * convert UTC back to TAI in order verify that we have the right answer. If we
     * dont, we'll have to do some more work.
     */
    UTC->JD   = Lgm_JD( UTC->Year, UTC->Month, UTC->Day, UTC->Time, LGM_TIME_SYS_UTC, c );
    UTC->Date = Lgm_JD_to_Date( UTC->JD, &UTC->Year, &UTC->Month, &UTC->Day, &tim );
    Lgm_IsLeapSecondDay( UTC->Date, &UTC->DaySeconds, c );
    UTC->Time = Lgm_RemapTime( UTC->Time, UTC->DaySeconds );
    Lgm_UTC_to_TAI( UTC, &t, c );


    /* 
     *  If subtracting the leap seconds took us back over a leapsecond
     *  (positive or minus), then the time will be off. We need to adjust d so
     *  that the calculation uses the correct number of leap seconds.
     */
    if ( fabs( t_target - t.Time ) > 0.1 ) {

        if ( (t_target - t.Time) > 0.1 ) {
            // occurs because the date we used to get leapseconds has more leap
            // seconds in it that we should have used.
            --d;
        } else {
            // occurs because the date we used to get leapseconds has fewer
            // leap seconds in it that we should have used. (deletions of leap
            // seconds has not occured yet, but theortetically could one
            // year...)
            ++d;
        }

        f = Lgm_GetLeapSeconds( d, c );
        UTC->Time  = TAI->Time - f/3600.0; 
        UTC->Year  = TAI->Year; UTC->Month = TAI->Month; UTC->Day = TAI->Day;
        UTC->JD = Lgm_JD( UTC->Year, UTC->Month, UTC->Day, UTC->Time, LGM_TIME_SYS_UTC, c );
        UTC->Date = Lgm_JD_to_Date( UTC->JD, &UTC->Year, &UTC->Month, &UTC->Day, &tim );
        Lgm_IsLeapSecondDay( UTC->Date, &UTC->DaySeconds, c );
        UTC->Time = Lgm_RemapTime( UTC->Time, UTC->DaySeconds );
    }
    

    // Set relavent fields in UTC structure.
    UTC->TimeSystem = LGM_TIME_SYS_UTC;
    UTC->T          = (UTC->JD - 2451545.0)/36525.0;
    UTC->TimeSystem = LGM_TIME_SYS_UTC;
    return;
}



/*
 *  TT -> TAI
 *  (Dont really need leap seconds here, but pass CTrans structure so all the args are the same.)
 */
void Lgm_TT_to_TAI( Lgm_DateTime *TT, Lgm_DateTime *TAI, Lgm_CTrans *c ) {

    double  Time;

    *TAI      = *TT; // initially make DateTime structures the same
    TAI->Time = TT->Time - 32.184/3600.0;
    TAI->JD   = Lgm_JD( TAI->Year, TAI->Month, TAI->Day, TAI->Time, LGM_TIME_SYS_TAI, c );
    TAI->Date = Lgm_JD_to_Date( TAI->JD, &TAI->Year, &TAI->Month, &TAI->Day, &Time );
    // Keep the time we had rather than getting it back from Lgm_JD_to_Date(). This limits roundoff error
    TAI->Time = Lgm_hour24( TAI->Time ); 
    TAI->T    = (TAI->JD - 2451545.0)/36525.0;
    TAI->TimeSystem = LGM_TIME_SYS_TAI;
    return;

}

// Does Lgm_JD need leapseconds for TAI ??? check that its working ok

/*
 *  TAI -> TT
 *  (Dont really need leap seconds here, but pass CTrans structure so all the args are the same.)
 */
void Lgm_TAI_to_TT( Lgm_DateTime *TAI, Lgm_DateTime *TT, Lgm_CTrans *c ) {

    double  Time;

    *TT      = *TAI; // initially make DateTime structures the same
    TT->Time = TAI->Time + 32.184/3600.0;
    TT->JD   = Lgm_JD( TT->Year, TT->Month, TT->Day, TT->Time, LGM_TIME_SYS_TT, c );
    TT->Date = Lgm_JD_to_Date( TT->JD, &TT->Year, &TT->Month, &TT->Day, &Time );
    // Keep the time we had rather than getting it back from Lgm_JD_to_Date(). This limits roundoff error
    TT->Time = Lgm_hour24( TT->Time ); 
    TT->T    = (TT->JD - 2451545.0)/36525.0;
    TT->TimeSystem = LGM_TIME_SYS_TT;
    return;

}










/*
 *  TT -> TDB ( IAU 1976 Resolution Version )
 */
void Lgm_TT_to_TDB( Lgm_DateTime *TT, Lgm_DateTime *TDB, Lgm_CTrans *c ) {

    double  Mearth, Time;

    *TDB     = *TT; // initially make DateTime structures the same

    Mearth    = (357.53 + 35999.050*TT->T)*RadPerDeg;
    TDB->Time = TT->Time + ( 0.001658*sin(Mearth) + 0.000014*sin(2.0*Mearth) )/3600.0;
    TDB->JD   = Lgm_JD( TDB->Year, TDB->Month, TDB->Day, TDB->Time, LGM_TIME_SYS_TDB, c );
    TDB->Date = Lgm_JD_to_Date( TDB->JD, &TDB->Year, &TDB->Month, &TDB->Day, &Time );
    // Keep the time we had rather than getting it back from Lgm_JD_to_Date(). This limits roundoff error
    TDB->Time = Lgm_hour24( TDB->Time );
    TDB->T    = (TDB->JD - 2451545.0)/36525.0;
    TDB->TimeSystem = LGM_TIME_SYS_TDB;

    return;
}



/*
 *  TT -> TDB ( IAU 2006 Resolution Version )
 */
void Lgm_TT_to_TDB_IAU2006( Lgm_DateTime *TT, Lgm_DateTime *TDB, Lgm_CTrans *c ) {

    double  Time;

    *TDB     = *TT; // initially make DateTime structures the same
    TDB->Time = TT->Time + (0.001658*sin( 628.3076*TT->T + 6.2401 )
                         + 0.000022*sin( 575.3385*TT->T + 4.2970 ) + 0.000014*sin( 1256.6152*TT->T + 6.1969 )
                         + 0.000005*sin( 606.9777*TT->T + 4.0212 ) + 0.000005*sin( 52.9691*TT->T + 0.4444 )
                         + 0.000002*sin( 21.3299*TT->T + 5.5431 ) + 0.000010*sin( 628.3076*TT->T + 4.2490 ))/3600.0;

    TDB->JD   = Lgm_JD( TDB->Year, TDB->Month, TDB->Day, TDB->Time, LGM_TIME_SYS_TDB, c );
    TDB->Date = Lgm_JD_to_Date( TDB->JD, &TDB->Year, &TDB->Month, &TDB->Day, &Time );
    // Keep the time we had rather than getting it back from Lgm_JD_to_Date(). This limits roundoff error
    TDB->Time = Lgm_hour24( TDB->Time );
    TDB->T    = (TDB->JD - 2451545.0)/36525.0;
    TDB->TimeSystem = LGM_TIME_SYS_TDB;
    return;
}



/*
 * TDB -> TT  
 * Unfortunately, this is a transendental equation, have to solve numerically.
 * The two times should be very close -- to within a few ms of each other...
 */
void Lgm_TDB_to_TT( Lgm_DateTime *TDB, Lgm_DateTime *TT, Lgm_CTrans *c ) {

    int             done;
    double          max, tdb_target, Ta, Tb, Tc, Ba, Bb, Bc, Fa, Fb, Fc, Time;
    Lgm_DateTime    *T, *B;

    // start out with approximation that TT ~= TDB. This is our "initial guess".
    *TT = *TDB;



    tdb_target = TDB->Time;


    T = (Lgm_DateTime *) calloc( 1, sizeof(Lgm_DateTime) );
    B = (Lgm_DateTime *) calloc( 1, sizeof(Lgm_DateTime) );

    // set T to our initial guess
    *T  = *TDB;


    /*
     * set up bracket
     */
    max = (0.001658 + 0.000014)/3600.0; // maximum possible diff...

    //(Ta, Ba)
    Ta = tdb_target - max;
    T->Time = Ta; Lgm_TT_to_TDB( T, B, c ); Ba = B->Time;
    // Lgm_TT_to_TDB() will recast the times if they cross date boundaries.(undo this.)
    if ( (Ba - Ta) > 12.0 ) Ba -= 24.0;
    if ( (Ta - Ba) > 12.0 ) Ba += 24.0;

    //(Tc, Bc)
    Tc = tdb_target + max;
    T->Time = Tc; Lgm_TT_to_TDB( T, B, c ); Bc = B->Time;
    // Lgm_TT_to_TDB() will recast the times if they cross date boundaries.(undo this.)
    if ( (Bc - Tc) > 12.0 ) Bc -= 24.0;
    if ( (Tc - Bc) > 12.0 ) Bc += 24.0;



    /*
       F is what we want to make zero
       F = TDB - TDB(TT)
       if we converge to a soln, it should be zero.
    */
    Fa = tdb_target - Ba;
    Fc = tdb_target - Bc;


    done = FALSE;
    while ( !done ) {


        // bisect and evaluate result at new point
        Tb = 0.5*(Ta+Tc);
        T->Time = Tb; Lgm_TT_to_TDB( T, B, c ); Bb = B->Time;
        if ( (Bb - Tb) > 12.0 ) Bb -= 24.0;
        if ( (Tb - Bb) > 12.0 ) Bb += 24.0;
        Fb = tdb_target - Bb;
            
        if ( Fb*Fa > 0.0 ) { 
            // we are on the `a' side
            Ta = Tb;
            Fa = Fb;
        } else {
            // we are on the `c' side
            Tc = Tb;
            Fc = Fb;
        } 

        if ( fabs(Ta-Tc)*86400.0 < 1e-8 ){
            done = TRUE;
        }


    }
    free( T ); 
    free( B );


    // Fill in some of the other fields correctly.
    TT->Time = 0.5*(Ta+Tc);
    TT->JD   = Lgm_JD( TT->Year, TT->Month, TT->Day, TT->Time, LGM_TIME_SYS_TT, c );
    TT->Date = Lgm_JD_to_Date( TT->JD, &TT->Year, &TT->Month, &TT->Day, &Time );
    // Keep the time we had rather than getting it back from Lgm_JD_to_Date(). This limits roundoff error
    TT->Time = Lgm_hour24( TT->Time );
    TT->T    = (TT->JD - 2451545.0)/36525.0;
    TT->TimeSystem = LGM_TIME_SYS_TT;

    return;
}


void Lgm_UTC_to_TT( Lgm_DateTime *UTC, Lgm_DateTime *TT, Lgm_CTrans *c ) {

    Lgm_DateTime TAI;
    double       Time;

//printf("\tLgm_UTC_to_TT: UTC =  "); Lgm_Print_DateTime( *UTC, 4, 8 ); printf("\n");


    TAI = *UTC;
    Lgm_UTC_to_TAI( UTC, &TAI, c );
//printf("\tLgm_UTC_to_TT: TAI =  "); Lgm_Print_DateTime( TAI, 4, 8 ); printf("\n");
    *TT = TAI;
    TT->Time = TAI.Time + 32.184/3600.0;
    TT->JD   = Lgm_JD( TT->Year, TT->Month, TT->Day, TT->Time, LGM_TIME_SYS_TT, c );
    TT->Date = Lgm_JD_to_Date( TT->JD, &TT->Year, &TT->Month, &TT->Day, &Time );
    // Keep the time we had rather than getting it back from Lgm_JD_to_Date(). This limits roundoff error
    TT->Time = Lgm_hour24( TT->Time );
    TT->T    = (TT->JD - 2451545.0)/36525.0;
//printf("\tLgm_UTC_to_TT: TT =  "); Lgm_Print_DateTime( *TT, 4, 8 ); printf("\n");
    TT->TimeSystem = LGM_TIME_SYS_TT;

    return;
}




void Lgm_TT_to_UTC( Lgm_DateTime *TT, Lgm_DateTime *UTC, Lgm_CTrans *c ) {

    Lgm_DateTime TAI;
    double       Time;

//printf("\tLgm_TT_to_UTC: TT =  "); Lgm_Print_DateTime( *TT, 4, 8 ); printf("\n");

    // Go TT -> TAI
    TAI = *TT;
    Lgm_TT_to_TAI( TT, &TAI, c );
//printf("\tLgm_TT_to_UTC: TAI =  "); Lgm_Print_DateTime( TAI, 4, 8 ); printf("\n");

    // Go TAI -> UTC
    Lgm_TAI_to_UTC( &TAI, UTC, c );

    // Fill in some of the other fields correctly
    UTC->JD   = Lgm_JD( UTC->Year, UTC->Month, UTC->Day, UTC->Time, LGM_TIME_SYS_UTC, c );
    UTC->Date = Lgm_JD_to_Date( UTC->JD, &UTC->Year, &UTC->Month, &UTC->Day, &Time );
    // Keep the time we had rather than getting it back from Lgm_JD_to_Date(). This limits roundoff error
    UTC->Time = Lgm_hour24( UTC->Time );
    UTC->T    = (UTC->JD - 2451545.0)/36525.0;
//printf("\tLgm_TT_to_UTC: UTC =  "); Lgm_Print_DateTime( *UTC, 4, 8 ); printf("\n");
    UTC->TimeSystem = LGM_TIME_SYS_UTC;

    return;
}

void Lgm_UTC_to_TDB( Lgm_DateTime *UTC, Lgm_DateTime *TDB, Lgm_CTrans *c ) {
    Lgm_DateTime TT;
//printf("\tLgm_UTC_to_TDB: UTC: "); Lgm_Print_DateTime( *UTC, 4, 8 ); printf("\n");
    Lgm_UTC_to_TT( UTC, &TT, c );
//printf("\tLgm_UTC_to_TDB: TT:  "); Lgm_Print_DateTime( TT, 4, 8 ); printf("\n");
    Lgm_TT_to_TDB( &TT, TDB, c );
//printf("\tLgm_UTC_to_TDB: TDB: "); Lgm_Print_DateTime( *TDB, 4, 8 ); printf("\n");
    TDB->TimeSystem = LGM_TIME_SYS_TDB;
    return;
}

void Lgm_TDB_to_UTC( Lgm_DateTime *TDB, Lgm_DateTime *UTC, Lgm_CTrans *c ) {
    Lgm_DateTime TT;
//printf("\tLgm_TDB_to_UTC: TDB: "); Lgm_Print_DateTime( *TDB, 4, 8 ); printf("\n");
    Lgm_TDB_to_TT( TDB, &TT, c );
//printf("\tLgm_TDB_to_UTC: TT:  "); Lgm_Print_DateTime( TT, 4, 8 ); printf("\n");
    Lgm_TT_to_UTC( &TT, UTC, c );
//printf("\tLgm_TDB_to_UTC: UTC: "); Lgm_Print_DateTime( *UTC, 4, 8 ); printf("\n");
    UTC->TimeSystem = LGM_TIME_SYS_UTC;
    return;
}










Lgm_DateTime *Lgm_DateTime_Create( int Year, int Month, int Day, double Time, int TimeSystem, Lgm_CTrans *c ) {

    double          t;
    Lgm_DateTime    *d;

    d = calloc( 1, sizeof(Lgm_DateTime) );

    d->Date   = Year*10000 + Month*100 + Day;
    d->Time   = Time;
    Lgm_Doy( d->Date, &d->Year, &d->Month, &d->Day, &d->Doy );
    d->JD     = Lgm_JD( d->Year, d->Month, d->Day, d->Time, TimeSystem, c );
    d->Date   = Lgm_JD_to_Date( d->JD, &d->Year, &d->Month, &d->Day, &t );
    // we could just set UTC->Time to the t from the previous call, but t will have
    // suffered round off errors. The following reduces roundoff errors.
    if      ( (d->Time - t) > 12.0 ) { d->Time -= 24.0; }
    else if ( (t - d->Time) > 12.0 ) { d->Time += 24.0; }

    if ( TimeSystem == LGM_TIME_SYS_UTC ) {
        Lgm_IsLeapSecondDay( d->Date, &d->DaySeconds, c );
    } else {
        d->DaySeconds = 86400.0;
    }

    d->Dow    = Lgm_DayOfWeek( d->Year, d->Month, d->Day, d->DowStr );
    d->Time   = Lgm_RemapTime( d->Time, d->DaySeconds );
    d->T      = (d->JD - 2451545.0)/36525.0;
    d->fYear  = (double)d->Year + ((double)d->Doy + d->Time/24.0)/(365.0 + (double)Lgm_LeapYear(d->Year));

    d->TimeSystem = TimeSystem;

    return(d);

}
void Lgm_DateTime_Destroy( Lgm_DateTime *d ) {
    free( d );
}

    

/*
 *  Fills up a UTC DateTime Structure.
 */
int Lgm_Make_UTC( long int Date, double Time, Lgm_DateTime *UTC, Lgm_CTrans *c ) {

    double  t;
    int     year, month, day, doy;

    Lgm_Doy( Date, &year, &month, &day, &doy );
    UTC->Date   = year*10000 + month*100 + day;
    UTC->Year   = year;
    UTC->Month  = month;
    UTC->Day    = day;
    UTC->Doy    = doy;
    UTC->Time   = Time;
//printf("Lgm_Make_UTC: Date, Year, Month, Day, Doy, Time = %ld, %d %d %d %d %lf\n", UTC->Date, UTC->Year, UTC->Month, UTC->Day, UTC->Doy, UTC->Time);

    // convert to JD.
    UTC->JD     = Lgm_JD( UTC->Year, UTC->Month, UTC->Day, UTC->Time, LGM_TIME_SYS_UTC, c );
//printf("Lgm_Make_UTC: Date, Year, Month, Day, Doy, Time = %ld, %d %d %d %d %lf\n", UTC->Date, UTC->Year, UTC->Month, UTC->Day, UTC->Doy, UTC->Time);

    UTC->Date   = Lgm_JD_to_Date( UTC->JD, &UTC->Year, &UTC->Month, &UTC->Day, &t );
    
    // redo the Doy calcs in case we got a new date (e.g. due to Time beingm > 24.0)
    Lgm_Doy( UTC->Date, &UTC->Year, &UTC->Month, &UTC->Day, &UTC->Doy );
//printf("Lgm_Make_UTC: Date, Year, Month, Day, Doy, Time, t = %ld, %d %d %d %d %lf %lf\n", UTC->Date, UTC->Year, UTC->Month, UTC->Day, UTC->Doy, UTC->Time, t);

    // we could just set UTC->Time to the t from the previous call, but t will have
    // suffered round off errors. The following reduces roundoff errors.
    if      ( (UTC->Time - t) > 12.0 ) { UTC->Time -= 24.0; }
    else if ( (t - UTC->Time) > 12.0 ) { UTC->Time += 24.0; }
//printf("Lgm_Make_UTC: Date, Year, Month, Day, Doy, Time = %ld, %d %d %d %d %lf\n", UTC->Date, UTC->Year, UTC->Month, UTC->Day, UTC->Doy, UTC->Time);
    Lgm_IsLeapSecondDay( UTC->Date, &UTC->DaySeconds, c );
    UTC->Time   = Lgm_RemapTime( UTC->Time, UTC->DaySeconds );
//printf("Lgm_Make_UTC: Date, Year, Month, Day, Doy, Time = %ld, %d %d %d %d %lf\n", UTC->Date, UTC->Year, UTC->Month, UTC->Day, UTC->Doy, UTC->Time);
    UTC->T      = (UTC->JD - 2451545.0)/36525.0;
    UTC->fYear  = (double)UTC->Year + ((double)UTC->Doy + UTC->Time/24.0)/(365.0 + (double)Lgm_LeapYear(UTC->Year));
    UTC->TimeSystem = LGM_TIME_SYS_UTC;

    t = UTC->Time;
    UTC->Hour   = (int)t;    t = (t - UTC->Hour)*60.0;
    UTC->Minute = (int)t;    t = (t - UTC->Minute)*60.0;
    UTC->Second = t;

    UTC->Dow    = Lgm_DayOfWeek( UTC->Year, UTC->Month, UTC->Day, UTC->DowStr );
    UTC->Week   = Lgm_ISO_WeekNumber( UTC->Year, UTC->Month, UTC->Day, &UTC->wYear );

    return( 1 ); // eventually return FALSE if invalid date.

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
void Lgm_Print_DateTime( Lgm_DateTime *DT, int Style, int p ){
    char *Str= calloc( 128, sizeof(char) );
    Lgm_DateTimeToString( Str, DT, Style, p );
    printf("%s", Str);
    free(Str);
}
void Lgm_DateTimeToString( char *Str, Lgm_DateTime *DT, int Style, int p ){
    int         HH, MM, SS, sgn, leapsec;
    int         Year, Month, d;
    double      S, SFRAC, Seconds, TotalSeconds;
    char        SecFracStr[30];

    leapsec = (int)(DT->DaySeconds - 86400.0);

    // Total Number of Seconds - leapsecs
    Seconds = DT->Time*3600.0;
    TotalSeconds = Seconds;
    
    // Hours
    HH = (int)(Seconds/3600.0);
    Seconds -= HH*3600.0;

    // Minutes
    MM = (int)(Seconds/60.0);
    Seconds -= MM*60.0;

//    Lgm_UT_to_HMSd( DT->Time, &sgn, &HH, &MM, &S );
    sgn = 1;


    S     = fabs(Seconds)+pow(10.0, -1.0-p );
    SS    = (int)S;
    SFRAC = S-(double)SS;


    if (( p <= 0 ) && ( SFRAC >= 0.5)) { SS += 1; }


    sprintf( SecFracStr, "%.*lf", (p>=20)?20:p, SFRAC );
    if ( leapsec && (TotalSeconds>=86400.0) && (TotalSeconds<86401.0) ){
        HH = 23; MM = 59; SS = 60;
    } else {
        if (SS==(60))  { SS=0; ++MM; }
        if (MM==60)  { MM=0; ++HH; }
    }

     d = DT->Date; Year = d/10000;
     d = d - Year*10000; Month = d/100;
     d = d - Month*100;

    // probably a cleaner way to do this(?)
    // Styles 0, 1, and 2 are ISO compliant...(check?)
    if ( p <= 0 ) {
        if ( Style == 0 ) {
            if ( DT->TimeSystem == LGM_TIME_SYS_UTC )
                sprintf( Str, "%02d-%02d-%02dT%02d:%02d:%02dZ", Year, Month, d, HH, MM, SS );
            else 
                sprintf( Str, "%02d-%02d-%02dT%02d:%02d:%02d", Year, Month, d, HH, MM, SS );
        } else if ( Style == 1 ) {
            if ( DT->TimeSystem == LGM_TIME_SYS_UTC )
                sprintf( Str, "%8ldT%02d%02d%02dZ", DT->Date, HH, MM, SS );
            else 
                sprintf( Str, "%8ldT%02d%02d%02d", DT->Date, HH, MM, SS );
        } else if ( Style == 2 ) {
            if ( DT->TimeSystem == LGM_TIME_SYS_UTC )
                sprintf( Str, "%8ld %02d:%02d:%02dZ", DT->Date, HH, MM, SS );
            else 
                sprintf( Str, "%8ld %02d:%02d:%02d", DT->Date, HH, MM, SS );
        } else if ( Style == 3 ) {
            if (sgn<0) sprintf( Str, "%8ld -%02d:%02d:%02d", DT->Date, HH, MM, SS );
            else       sprintf( Str, "%8ld  %02d:%02d:%02d", DT->Date, HH, MM, SS );
        } else {
            if (sgn<0) sprintf( Str, "%8ld -%02d\u02b0 %02d\u1d50 %02d\u02e2", DT->Date, HH, MM, SS );
            else       sprintf( Str, "%8ld  %02d\u02b0 %02d\u1d50 %02d\u02e2", DT->Date, HH, MM, SS );
        }
    } else {
        if ( Style == 0 ) {
            if ( DT->TimeSystem == LGM_TIME_SYS_UTC )
                sprintf( Str, "%02d-%02d-%02dT%02d:%02d:%02d.%sZ", Year, Month, d, HH, MM, SS, strstr(SecFracStr, ".")+1 );
            else 
                sprintf( Str, "%02d-%02d-%02dT%02d:%02d:%02d.%s", Year, Month, d, HH, MM, SS, strstr(SecFracStr, ".")+1 );
        } else if ( Style == 1 ) {
            if ( DT->TimeSystem == LGM_TIME_SYS_UTC )
                sprintf( Str, "%8ldT%02d%02d%02d.%sZ", DT->Date, HH, MM, SS, strstr(SecFracStr, ".")+1 );
            else 
                sprintf( Str, "%8ldT%02d%02d%02d.%s", DT->Date, HH, MM, SS, strstr(SecFracStr, ".")+1 );
        } else if ( Style == 2 ) {
            if ( DT->TimeSystem == LGM_TIME_SYS_UTC )
                sprintf( Str, "%8ld %02d:%02d:%02d.%sZ", DT->Date, HH, MM, SS, strstr(SecFracStr, ".")+1 );
            else 
                sprintf( Str, "%8ld %02d:%02d:%02d.%s", DT->Date, HH, MM, SS, strstr(SecFracStr, ".")+1 );
        } else if ( Style == 3 ) {
            if (sgn<0) sprintf( Str, "%8ld -%02d:%02d:%02d.%s", DT->Date, HH, MM, SS, strstr(SecFracStr, ".")+1 );
            else       sprintf( Str, "%8ld  %02d:%02d:%02d.%s", DT->Date, HH, MM, SS, strstr(SecFracStr, ".")+1 );
        } else {
            if (sgn<0) sprintf( Str, "%8ld -%02d\u02b0 %02d\u1d50 %02d\u02e2.%s", DT->Date, HH, MM, SS, strstr(SecFracStr, ".")+1 );
            else       sprintf( Str, "%8ld  %02d\u02b0 %02d\u1d50 %02d\u02e2.%s", DT->Date, HH, MM, SS, strstr(SecFracStr, ".")+1 );
        }
    }

}


void Lgm_Print_SimpleTime( Lgm_DateTime *DT, int p, char *Str ){
    int         HH, MM, SS, sgn, leapsec;
    int         Year, Month, d;
    double      S, SFRAC, Seconds, TotalSeconds;
    char        SecFracStr[30];

    leapsec = (int)(DT->DaySeconds - 86400.0);

    // Total Number of Seconds - leapsecs
    Seconds = DT->Time*3600.0;
    TotalSeconds = Seconds;
    
    // Hours
    HH = (int)(Seconds/3600.0);
    Seconds -= HH*3600.0;

    // Minutes
    MM = (int)(Seconds/60.0);
    Seconds -= MM*60.0;

//    Lgm_UT_to_HMSd( DT.Time, &sgn, &HH, &MM, &S );
    sgn = 1;


    S     = fabs(Seconds)+pow(10.0, -1.0-p );
    SS    = (int)S;
    SFRAC = S-(double)SS;


    if (( p <= 0 ) && ( SFRAC >= 0.5)) { SS += 1; }


    sprintf( SecFracStr, "%.*lf", (p>=20)?20:p, SFRAC );
    if ( leapsec && (TotalSeconds>=86400.0) && (TotalSeconds<86401.0) ){
        HH = 23; MM = 59; SS = 60;
   } else {
        if (SS==(60))  { SS=0; ++MM; }
        if (MM==60)  { MM=0; ++HH; }
    }

     d = DT->Date; Year = d/10000;
     d = d - Year*10000; Month = d/100;
     d = d - Month*100;

    // probably a cleaner way to do this(?)
    if ( p <= 0 ) {
        sprintf(Str, "%02d:%02d:%02d", HH, MM, SS );
    } else {
        sprintf(Str, "%02d:%02d:%02d.%s", HH, MM, SS, strstr(SecFracStr, ".")+1 );
    }

}


/*

 *  Compute the Julian Date for the given Gregorian Date and Time.
 *  Julian Date is the number of days since noon of Jan 1 4713 B.C.
 *  Only need to worry about leap seconds if TimeSystem is UTC.
 *  
 *  Computes the Julian Date for a given time. Recall that JDN is an integer
 *  corresponding to noon on the given Gregorian date.  Thus, 
 *
 *                  JDN - 0.5 
 *
 *  is the Julian Date at the beginning of the day. And,
 *
 *                 JDN - 0.5 + FractionOfDay 
 *
 *  is the Julian date corresponding to the given Time. Note that the Julian
 *  Date can be constructed from times in different time systems (e.g. UTC,
 *  TAI, TT, TDB, etc.). The reason this matters, is that leap seconds are
 *  occasionally introduced in UTC but not in the other time sytems. Thus, the
 *  number of seconds in a day can change in UTC.  There will be 86400 seconds
 *  on most days, but on leap second days there will be 86401 seconds. This
 *  routine first computes the number of secons in a day, then computes the
 *  fraction of a day. This fraction is then added on to the Julian Date at the
 *  beginning of the day (i.e. JDN-0.5).
 *  
 */
double Lgm_JD( int Year, int Month, int Day, double Time, int TimeSystem, Lgm_CTrans *c ) {

    long int    Date, JDN;
    double      ns, FracDay, JD;

    // Get JDN
    JDN = Lgm_JDN( Year, Month, Day );
//printf("JDN, Year, Month, Day = %ld %d %d %d\n", JDN, Year, Month, Day);

    // Form Date
    Date = Year*10000 + Month*100 + Day;

    // Set number of seconds in this day -- only applies to UTC system
    if ( TimeSystem == LGM_TIME_SYS_UTC ) {
        Lgm_IsLeapSecondDay( Date, &ns, c );
        // Compute fractional part of current day
        FracDay = Time*3600.0/ns;
    } else {
        ns = 86400.0;
        // Compute fractional part of current day
        FracDay = Time/24.0;
    }
//printf("ns = %lf\n", ns);


    // Shift back to midnight and then add in the Fractional part
    JD = JDN - 0.5 + FracDay;
//printf("FracDay = %lf\n", FracDay);

    return( JD );

}



/*
 *  Compute the Julian Day Number for the given Gregorian Date.
 *  Julian Date is the number of days since noon of Jan 1 4713 B.C.
 *  See Expl. Suppl. to Astron. Almanac, p. 604  (2006 edition).
 *  The long integer number returned here is the Julian Day Number
 *  for Noon on the given calendar date.
 */
long int Lgm_JDN( int Year, int Month, int Day ) {
    int q = (Month < 3 ) ? -1 : 0;
    // Compute integer part of julian date
    return(  (1461*(Year + 4800 + q))/4 + (367*(Month - 2 - 12*q))/12
            - (3*((Year + 4900 + q)/100))/4 + Day - 32075 );
}



double Lgm_MJD(int ny, int nm, int nd, double UT, int TimeSystem, Lgm_CTrans *c ) {
    return( Lgm_JD( ny, nm, nd, UT, TimeSystem, c ) - 2400000.5 );
}


double Lgm_Date_to_JD( long int Date, double UT, Lgm_CTrans *c) {

    int year, month, day, doy;

    Lgm_Doy( Date, &year, &month, &day, &doy);
    return( Lgm_JD( year, month, day, UT, LGM_TIME_SYS_UTC, c) );

}

int Lgm_DayOfYear(int year, int month, int day, Lgm_CTrans *c) {
    return((int)(Lgm_JD(year, month, day, 0.0, LGM_TIME_SYS_UTC, c) - Lgm_JD(year, 1, 0, 0.0, LGM_TIME_SYS_UTC, c)));
}




/*
 * Here the number goes from 0-6 Sun-Sat.
 * But ISO 8601 goes from: 1 is Monday to 7 is Sunday
 * In the Lgm_DateTime struct we will use ISO defs thats why the last bit that changes 0 to 7 if we get a zero....
 */
int Lgm_DayOfWeek( int Year, int Month, int Day, char *dowstr ) {

    int q, m, K, J, Y, h, d;

    q = Day; 
    m = Month; 
    Y = Year;
    if (m<3){
        --Y; 
        m+=12;
    }
    K = Y%100;
    J = Year/100;                                             

    h = ( q + (((m+1)*26)/10) + Y + (Y/4) + 6*(Y/100) + (Y/400))%7;          
    d = ((h+5)%7)+1;          

    if      ( d == 1 ) {  strcpy(dowstr, "Mon");  }
    else if ( d == 2 ) {  strcpy(dowstr, "Tue");  }
    else if ( d == 3 ) {  strcpy(dowstr, "Wed");  }
    else if ( d == 4 ) {  strcpy(dowstr, "Thu");  }
    else if ( d == 5 ) {  strcpy(dowstr, "Fri");  }
    else if ( d == 6 ) {  strcpy(dowstr, "Sat");  }
    else if ( d == 7 ) {  strcpy(dowstr, "Sun");  }


    return(d);

}

/*
 * Compute the JDN of the start of week 1 for the given year
 * This relies on the fact that Jan 4 is always in the first week of its year.
 */
long int Lgm_JDNofWeek1( int Year ) {
    char     str[10];
    int      Dow_JAN4 = Lgm_DayOfWeek( Year, 1, 4, str );
    long int JDN_JAN4 = Lgm_JDN( Year, 1, 4 );
    return( JDN_JAN4 + 1 - Dow_JAN4 );
}

/*
 * Compute ISO week number
 */
int Lgm_ISO_WeekNumber( int Year, int Month, int Day, int *ISO_WeekYear ) {

    int      ISO_Year  = Year;
    long int JDN_DEC29 = Lgm_JDN( Year, 12, 29 );
    long int JDN       = Lgm_JDN( Year, Month, Day );
    long int JDN_Week1;

    if ( JDN >= JDN_DEC29 ) {
        JDN_Week1 = Lgm_JDNofWeek1( ISO_Year+1 );
        if ( JDN < JDN_Week1 ){
            JDN_Week1 = Lgm_JDNofWeek1( ISO_Year );
        } else {
            ++ISO_Year;
        }
    } else {
        JDN_Week1 = Lgm_JDNofWeek1( ISO_Year );
        if ( JDN < JDN_Week1 ) JDN_Week1 = Lgm_JDNofWeek1( --ISO_Year );
    }

    *ISO_WeekYear = ISO_Year;
    return( (JDN-JDN_Week1)/7 + 1 );

}


/*
 * Compute the last week of the year. It can be 52 or 53.
 * This relies on the fact that Dec 28 is always in the last week of its year.
 */
int Lgm_MaxWeekNumber( int Year ) {
    int tmp;
    return( Lgm_ISO_WeekNumber( Year, 12, 28, &tmp ) );
}

/*
 * Convert ISO_WeekYear/Week/Dow back to Date
 */
void Lgm_ISO_YearWeekDow_to_Date( int ISO_WeekYear, int Week, int Dow, long int *Date, int *Year, int *Month, int *Day ) {
    long int    JDN_Week1 = Lgm_JDNofWeek1( ISO_WeekYear );
    long int    JDN       = JDN_Week1 + 7*(Week-1) + Dow - 1;
    double      tmp;
    Lgm_jd_to_ymdh( JDN, Date, Year, Month, Day, &tmp);
}



/*
 *  Routine converts between day number and mmdd formats. Input date can be in any of the following
 *  formats;
 *      yyddd   (where yy is assumed to be a year in the 20th century)
 *      yymmdd  (where yy is assumed to be a year in the 20th century)
 *      yyyyddd   (where yyyy is assumed to be between 1000 A.D. and 9999 A.D.)
 *      yyyymmdd  (where yyyy is assumed to be between 1000 A.D. and 9999 A.D.)
 *
 *  Should be good enough for space physics! Need separate routines or add a flag or something
 *  if you want it to be completely general.... Its not completely general now but it is
 *  convenient.
 *
 */
int Lgm_Doy( long int date, int *YY, int *MM, int *DD, int *DOY) {

    int ddd, n, i, m;
    double f;
    static double days1[] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    static double days2[] = {31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};

    if ((date < 1e5)||((date < 1e7)&&(date > 1e6))){    /* then the date is in yyddd or yyyyddd format */
    *YY = (int)(date/1000.0);
    *DOY = date - *YY * 1000;
    if (date < 1e5) *YY += 1900;

    ddd = 0, n = 0;
    if (Lgm_LeapYear(*YY)){
        while ( (ddd+days2[n]) < *DOY) {
            ddd += days2[n];
            ++n;
        }
    }
    else{
        while ( (ddd+days1[n]) < *DOY) {
            ddd += days1[n];
            ++n;
        }
    }
    *MM = n+1;
    *DD = *DOY - ddd;
    }
    else{           /* then the date is in yymmdd or yyyymmdd format */
    f = date/10000.0;
    *YY = (int)f;

    f = date - *YY * 10000;
    *MM = (int)(f/100.0);
    m = *MM - 1;

    f = f - *MM * 100;
    *DD = (int)f;

    if (Lgm_LeapYear(  (date < 1e7) ? 1900 + *YY : *YY  )){
        for (i=0, ddd=0; i<m; ++i) ddd += days2[i];
    }
    else{
        for (i=0, ddd=0; i<m; ++i) ddd += days1[i];
    }

    *DOY = ddd + *DD;

    if (date < 1e7) *YY += 1900;
    }

    return( *DOY );

}


int Lgm_IsValidDate( long int Date ) {

    int yy, mm, dd, doy;
    static int mdays[]  = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    static int mdays2[] = {31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};

    if (Date < 0) return(0);

    Lgm_Doy(Date, &yy, &mm, &dd, &doy);

    if ((mm < 1)||(mm > 12)) return(0);

    if (dd < 1) return(0);

    if (Lgm_LeapYear(yy)){
        if (dd > mdays2[mm-1]) return(0);
    } else {
        if (dd > mdays[mm-1]) return(0);
    }

    /*
     *   Must be a valid date (at least in terms of month and day)
     */
    return(1);

}

double Lgm_RemapTime( double Time, double SecondsInADay ) {

    double  Seconds;

    Seconds = Time*3600.0;

    if ( Seconds < 0.0 )            Seconds += SecondsInADay;
    if ( Seconds >= SecondsInADay ) Seconds -= SecondsInADay;
    
    Time = Seconds/3600.0;

    return(Time);

}




/*
 * We really should examine the TimeSystem field in the Lgm_DateTime structure
 * and convert appropriately . But for now lets just assume input is in UTc.
 */
double Lgm_UTC_to_TDBSeconds( Lgm_DateTime *UTC, Lgm_CTrans *c ){
    return( Lgm_TDBSecSinceJ2000( UTC, c ) );
}


/*
 * TDBSecSinceJ2000 is like SPICE's 'ephemeris time'. Probably differs in 7th
 * decimal place or so. Convenience rotuine to get TdbSecSinceJ2000 driectly
 * from UTC rather than from TDB
 */
double Lgm_TDBSecSinceJ2000( Lgm_DateTime *UTC, Lgm_CTrans *c ){

    long int        JDN;
    double          DaySeconds, Seconds;
    Lgm_DateTime    TDB;

    // Convet to TDB
    Lgm_UTC_to_TDB( UTC, &TDB, c );

    // Get JDN
    JDN = Lgm_JDN( TDB.Year, TDB.Month, TDB.Day );

    // Get seconds past start of day (i.e. midnight)
    DaySeconds = TDB.Time*3600.0;

    // Compute Seconds since J2000
    Seconds = ((JDN-LGM_JD_J2000)*86400 - 43200) + DaySeconds;

    return( Seconds );

}



double  Lgm_TDB_to_TdbSecSinceJ2000( Lgm_DateTime *TDB ) {
    double  JDN, Seconds;
    JDN  = Lgm_JDN( TDB->Year, TDB->Month, TDB->Day );
    Seconds = ((JDN - LGM_JD_J2000)*86400.0 - 43200.0) + TDB->Time*3600.0;
    return( Seconds );
}

void Lgm_TdbSecSinceJ2000_to_TDB( double TdbSeconds, Lgm_DateTime *TDB ) {

    long int dJDN;
    double   tmp, JDN;
    dJDN = (long int)((TdbSeconds + 43200.0)/86400.0);
    TDB->Time = (TdbSeconds - dJDN*86400.0 + 43200.0)/3600.0;
    JDN = LGM_JD_J2000 + dJDN;

    Lgm_jd_to_ymdh( JDN, &TDB->Date, &TDB->Year, &TDB->Month, &TDB->Day, &tmp);
    
    TDB->JD   = JDN + TDB->Time/24.0;
    TDB->T    = (TDB->JD - 2451545.0)/36525.0;
    TDB->TimeSystem = LGM_TIME_SYS_TDB;

}

void Lgm_TdbSecSinceJ2000_to_UTC( double TdbSeconds, Lgm_DateTime *UTC, Lgm_CTrans *c ) {
    Lgm_DateTime TDB;
    Lgm_TdbSecSinceJ2000_to_TDB( TdbSeconds, &TDB );
    Lgm_TDB_to_UTC( &TDB, UTC, c );
    UTC->TimeSystem = LGM_TIME_SYS_UTC;
}

double Lgm_UTC_to_TdbSecSinceJ2000( Lgm_DateTime *UTC, Lgm_CTrans *c ) {
    Lgm_DateTime TDB;
    Lgm_UTC_to_TDB( UTC, &TDB, c );
    return( Lgm_TDB_to_TdbSecSinceJ2000( &TDB) );
}


/*
 *  Convert JD directly to a UTC DateTime structure.
 */
void Lgm_JD_to_DateTime( double JD, Lgm_DateTime *UTC, Lgm_CTrans *c ){

    long int Date;
    int      Year, Month, Day;
    double  Time;

    Date = Lgm_JD_to_Date( JD, &Year, &Month, &Day, &Time );
    Lgm_Make_UTC( Date, Time, UTC, c );

}

/*
 *  Convert MJD directly to a UTC DateTime structure
 */
void Lgm_MJD_to_DateTime( double MJD, Lgm_DateTime *UTC, Lgm_CTrans *c ) {
    double   JD;

    JD = MJD + 2400000.5;
    Lgm_JD_to_DateTime( JD, UTC, c );
}
