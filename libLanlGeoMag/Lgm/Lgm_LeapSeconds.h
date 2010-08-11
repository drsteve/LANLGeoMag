#ifndef LGM_LEAPSECONDS_H
#define LGM_LEAPSECONDS_H

#define STRINGIFY(x) #x
#define EXPAND(x) STRINGIFY(x)


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define FALSE           0
#define TRUE            1


typedef struct Lgm_LeapSeconds {

    int         nLeapSecondDates;   /*
                                     * Number of leap second dates.
                                     */

    long int    *LeapSecondDates;   /*
                                     * Array for holdin the Dates on which leap
                                     * seconds were added
                                     */

    double      *LeapSecondJDs;     /*
                                     * Array for holdin the Julian Dates on
                                     * which leap seconds were added
                                     */



    double      *LeapSeconds;       /* 
                                     * The actual number of leap seconds that
                                     * went into effect on the given date
                                     */


} Lgm_LeapSeconds;


int         Lgm_LoadLeapSeconds( Lgm_LeapSeconds  *l );
double      Lgm_GetLeapSeconds( double JD, Lgm_LeapSeconds *l );
int         Lgm_IsLeapSecondDay( long int Date, Lgm_LeapSeconds *l );
double      Lgm_UTC_to_TAI( double JD, double UTC, Lgm_LeapSeconds *l );
double      Lgm_TAI_to_UTC( double JD, double TAI, Lgm_LeapSeconds *l );
double      Lgm_UTC_to_TT( double JD, double UTC, Lgm_LeapSeconds *l );
double      Lgm_TT_to_UTC( double JD, double TT, Lgm_LeapSeconds *l );
double      Lgm_TT_to_TDB( double JD, double TT, Lgm_LeapSeconds *l );
double      Lgm_TDB_to_TT( double JD, double TDB, Lgm_LeapSeconds *l );
double      Lgm_UTC_to_TDB( double JD, double UTC, Lgm_LeapSeconds *l );
double      Lgm_TDB_to_UTC( double JD, double TDB, Lgm_LeapSeconds *l );



#endif
