#ifndef LGM_CTRANS_H
#define LGM_CTRANS_H

#ifndef STRINGIFY
#define STRINGIFY(x) #x
#endif

#ifndef EXPAND
#define EXPAND(x) STRINGIFY(x)
#endif


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Lgm_Vec.h"
#include "Lgm_WGS84.h"
#include "Lgm_JPLeph.h"

#define DegPerRad       57.295779513082320876798154814105
#define RadPerDeg        0.017453292519943295769236907568
#define RadPerArcSec     4.848136811095359935899141023579e-6    /* radians per arc-second */

#ifndef M_SQRTPI
#define M_SQRTPI         1.772453850905516027298167483341       /*  sqrt(PI)      */
#endif

#define M_SQRTPI_2       0.886226925452758013649083741671       /*  sqrt(PI)/2.0  */
#define M_1_SQRTPI       0.564189583547756286948079451561       /*  1/sqrt(PI)    */
#ifndef M_2PI
#define M_2PI            6.283185307179586476925286766559       /*  2PI           */
#endif
#define M_OneThird       0.333333333333333333333333333333       /* 1.0/3.0        */
#define M_SQRT_3        1.7320508075688772935274463415058       /* sqrt(3)        */

#ifndef FALSE
#define FALSE           0
#endif

#ifndef TRUE
#define TRUE            1
#endif



#define Re              6378.137 	    // This is consistent with equatorial
                                        // radius in the WGS84 model [km]

#define AU              149.5978700e6   // Astronomical Unit in km

#define SOLAR_RADIUS    696342.0        // km (from SOHO meas.)
#define LUNAR_RADIUS    1738.14         // km equatorial radius


#define LGM_NO_ECLIPSE          0
#define LGM_PENUMBRAL_ECLIPSE   1
#define LGM_UMBRAL_ECLIPSE      2


#define LGM_GOLD                1.61803398874989484820  // Phi -- Golden ratio
#define LGM_1_OVER_GOLD         0.61803398874989484820  // 1/Phi
#define LGM_1_MINUS_1_OVER_GOLD 0.38196601125010515180  // 1-1/Phi

#define LGM_ERROR             -1  // generic error flag
#define LGM_FILL_VALUE    (-1e31) // Value to flag undefined or bad data
#define LGM_LARGE_NEG_DBL (-9e99) // A large negative double value
#define LGM_LARGE_POS_DBL (9e99)  // A large positive double value

#define  LGM_JD_J2000   2451545.0 // Julian Date of J2000 ( i.e. 12:00:00 TT on Jan 1 2000 (or 2000 Jan 1.5 TT ) )
#define  LGM_JD_GPS0    2444245.0 // Julian Date of introduction of GPS time. (0h, Jan 6, 1980)
#define  LGM_JD_TAI0    2436205.0 // Julian Date of introduction of TAI time. (0h, Jan 1, 1958)

#define LGM_EPH_LOW_ACCURACY    0
#define LGM_EPH_HIGH_ACCURACY   1
#define LGM_EPH_DE              2
#define LGM_PN_IAU76            12
#define LGM_PN_IAU06            12


/*
 * Define Time Systems
 */
#define LGM_TIME_SYS_UTC    0
#define LGM_TIME_SYS_TAI    1
#define LGM_TIME_SYS_GPS    2
#define LGM_TIME_SYS_TT     3
#define LGM_TIME_SYS_TDB    4
#define LGM_TIME_SYS_UT1    5



/*
 * Definition of Coordinate Systems.
 * We follow the International Astronomical Union 1976/FK-5 conventions
 * (IAU-76/FK5). E.g. see Vallado, "Fundamentals of Astrodynamics and
 * Applications" Third Edition, 2007.  The newer IAU-2000 theory has similar
 * (but different in details) formulae from these -- be careful not to mix
 * them.
 */

#define EME2000_COORDS  1   //!<  Earth Mean Equator and Equinox of Epoch J2000 (EME2000 system)
#define ICRF2000_COORDS 1   //!<  aka International Celestial Reference Frame (ICRF)
#define GEI2000_COORDS  1   //!<  aka GEI (Geocentric Equatorial Inertial) at Epoch J2000
                            //!<
                            //!<      Z-axis - parallel to mean rot axis of Earth
                            //!<               at fixed epoch J2000
                            //!<      X-axis - points to direction of mean vernal
                            //!<               equinox at the fixed epoch J2000
                            //!<      Y-axis - completes right handed system
                            //!<

#define MOD_COORDS      2   //!< Mean Of Date (MOD) system.
                            //!<
                            //!< This is same as EME2000_COORDS, except that the
                            //!< mean rotation axis and mean equinox of date are
                            //!< used to define the system.  Transforming between
                            //!< J2000 and MOD involves roatations using the
                            //!< so-called Precession angles (Zeta, Zee, Theta).
                            //!< The sun position (and other quantities) that we
                            //!< compute natively comes out in MOD coords.
                            //!< Transformation between MOD and EME2000
                            //!<      Umod = Rz(-Zee)Ry(Theta)Rz(-Zeta) Ueme2000

#define TOD_COORDS      3   //!< True Of Date (TOD) system.
                            //!<
                            //!< This is same as MOD, except that the true
                            //!< rotation axis and true equinox of date are used
                            //!< to define the system.  Transforming between MOD
                            //!< and TOD involves roatations using the Nutation
                            //!< corrections delta-Psi and delta-Eps.
                            //!< Transformation between TOD and MOD
                            //!<      Utod = Rx( -(Eps+dEps) )Rz(-dPsi)Rx(Eps) Umod

#define TEME_COORDS     4   //!< True Equator, Mean Equinox (of Date).
                            //!<
                            //!< This is a hybrid system using true equat., but
                            //!< mean equinox.  This is the system that the
                            //!< output of the SGP4 orbit propagator uses
                            //!< implicitly. The system comes about implicitly
                            //!< when mean sidereal time is used in place of true
                            //!< sidereal time when transforming between an
                            //!< earth-fixed system and an interial one.  Its not
                            //!< a system we would necessarily want to work in,
                            //!< but its what SGP4 (implicitly) outputs due to
                            //!< the nature of how the TLEs and orbit
                            //!< calculations are done.
                            //!< Transformation between TEME and TOD
                            //!<      Utod = Rz( dPsi cos(Eps) ) Uteme


#define PEF_COORDS      5   //!< Pseudo Earth Fixed (PEF) system
                            //!<
                            //!< This system has Z-axis aligned with
                            //!< instantaneous rotation axis of the Earth and
                            //!< X-axis to Greenwich meridian. Its closely
                            //!< related to TOD system. PEF to TOD are related by
                            //!< a single rotation around z-axis by angle
                            //!< correspondingm to true sidereal time. This is
                            //!< called Pseudo-Fixed because in reality the
                            //!< rotation axis moves slightly relative to the
                            //!< ground (pole wander). So its not a truely
                            //!< Earth-fixed system.
                            //!< Transformation:
                            //!<      Upef = Rz( Theta ) Utod
                            //!<      (Theta is true sidereal time; Theta = Mean Theta + Eqn or Equinoxes)


#define WGS84_COORDS    6   //!< World Geodetic System 1984 (WGS84)
#define ITRF_COORDS     6   //!< ( aka International Terrestrial Reference Frame)
#define GEO_COORDS      6   //!< ( aka "GEO" coords -- my term )
                            //!< (WGS84 and ITRF are not *exactly* the same, but
                            //!< supposedly agree at the cm level and within
                            //!< uncertainties of WGS84 are often considered to
                            //!< be equivelant)
                            //!<
                            //!< This is a truely Earth-fixed system. It's Z-axis
                            //!< is the true International Reference Pole and
                            //!< Meridian. And these are fixed in the crust of
                            //!< the Earth.  To go from PEF to WGS84, you need to
                            //!< know the angular offsets of the true pole from
                            //!< the PEF pole. These angle,  (xp and yp) are part
                            //!< of the so-called Earth-Observing-Parameters
                            //!< (EOP). They cant be predicted ahead of time, so
                            //!< we often just assume they are zero.  Ignoring
                            //!< them can introduce offsets on the ground of
                            //!< about 10m or so. The WGS84 system is the sytem
                            //!< that most GPS equipment uses(check?). Note also
                            //!< that geodetic latitude is defined differently
                            //!< from geocenteric latitude due to the flattening
                            //!< of the Earth in the z-direction. If the EOP values
                            //!< are not applied then WGS84 is same as PEF.
                            //!< Transformation:
                            //!<      Uwgs84 = Ry(-xp)Rz(-yp) Upef


#define GSE_COORDS      7   //!< Geocentric Solar Ecliptic (mean of date quantities?)

#define GSM_COORDS      8   //!< Geocentric Solar Magnetospheric (X to sun, mag dipole in X-Z plane)

#define SM_COORDS       9   //!< Solar Magnetic (Mag dipole in Z, sun in X-Z plane)

#define EDMAG_COORDS   10   //!< Eccentric Dipole Coords (true dipole is offset from center)

#define CDMAG_COORDS   11   //!< Centered Dipole Coords

#define GSE2000_COORDS 12   //!< GSE2000 coordinates (for MMS mission)
                            //!< X is J2000 Sun vector, Y is cross-product of X and ecliptic pole at J2000 epoch
                            //!< Z is cross-product of X and Y



/*
 * These are defined to indicate what coordinate conversion you want If START
 * is the starting coordinate system and END is the desired ending coordinate
 * system the START_TO_END encodes the desired transformation.
 */
#define EME2000_TO_EME2000      101
#define EME2000_TO_ICRF2000     101
#define EME2000_TO_GEI2000      101
#define EME2000_TO_MOD          102
#define EME2000_TO_TOD          103
#define EME2000_TO_TEME         104
#define EME2000_TO_PEF          105
#define EME2000_TO_WGS84        106
#define EME2000_TO_ITRF         106
#define EME2000_TO_GEO          106
#define EME2000_TO_GSE          107
#define EME2000_TO_GSM          108
#define EME2000_TO_SM           109
#define EME2000_TO_EDMAG        110
#define EME2000_TO_CDMAG        111
#define EME2000_TO_GSE2000      112

#define ICRF2000_TO_EME2000     101
#define ICRF2000_TO_ICRF2000    101
#define ICRF2000_TO_GEI2000     101
#define ICRF2000_TO_MOD         102
#define ICRF2000_TO_TOD         103
#define ICRF2000_TO_TEME        104
#define ICRF2000_TO_PEF         105
#define ICRF2000_TO_WGS84       106
#define ICRF2000_TO_ITRF        106
#define ICRF2000_TO_GEO         106
#define ICRF2000_TO_GSE         107
#define ICRF2000_TO_GSM         108
#define ICRF2000_TO_SM          109
#define ICRF2000_TO_EDMAG       110
#define ICRF2000_TO_CDMAG       111
#define ICRF2000_TO_GSE2000     112

#define GEI2000_TO_EME2000      101
#define GEI2000_TO_ICRF2000     101
#define GEI2000_TO_GEI2000      101
#define GEI2000_TO_MOD          102
#define GEI2000_TO_TOD          103
#define GEI2000_TO_TEME         104
#define GEI2000_TO_PEF          105
#define GEI2000_TO_WGS84        106
#define GEI2000_TO_ITRF         106
#define GEI2000_TO_GEO          106
#define GEI2000_TO_GSE          107
#define GEI2000_TO_GSM          108
#define GEI2000_TO_SM           109
#define GEI2000_TO_EDMAG        110
#define GEI2000_TO_CDMAG        111
#define GEI2000_TO_GSE2000      112

#define MOD_TO_EME2000          201
#define MOD_TO_ICRF2000         201
#define MOD_TO_GEI2000          201
#define MOD_TO_MOD              202
#define MOD_TO_TOD              203
#define MOD_TO_TEME             204
#define MOD_TO_PEF              205
#define MOD_TO_WGS84            206
#define MOD_TO_ITRF             206
#define MOD_TO_GEO          	206
#define MOD_TO_GSE              207
#define MOD_TO_GSM              208
#define MOD_TO_SM               209
#define MOD_TO_EDMAG            210
#define MOD_TO_CDMAG            211
#define MOD_TO_GSE2000          212

#define TOD_TO_EME2000          301
#define TOD_TO_ICRF2000         301
#define TOD_TO_GEI2000          301
#define TOD_TO_MOD              302
#define TOD_TO_TOD              303
#define TOD_TO_TEME             304
#define TOD_TO_PEF              305
#define TOD_TO_WGS84            306
#define TOD_TO_ITRF             306
#define TOD_TO_GEO              306
#define TOD_TO_GSE              307
#define TOD_TO_GSM              308
#define TOD_TO_SM               309
#define TOD_TO_EDMAG            310
#define TOD_TO_CDMAG            311
#define TOD_TO_GSE2000          312

#define TEME_TO_EME2000         401
#define TEME_TO_ICRF2000        401
#define TEME_TO_GEI2000         401
#define TEME_TO_MOD             402
#define TEME_TO_TOD             403
#define TEME_TO_TEME            404
#define TEME_TO_PEF             405
#define TEME_TO_WGS84           406
#define TEME_TO_ITRF            406
#define TEME_TO_GEO             406
#define TEME_TO_GSE             407
#define TEME_TO_GSM             408
#define TEME_TO_SM              409
#define TEME_TO_EDMAG           410
#define TEME_TO_CDMAG           411
#define TEME_TO_GSE2000         412

#define PEF_TO_EME2000          501
#define PEF_TO_ICRF2000         501
#define PEF_TO_GEI2000          501
#define PEF_TO_MOD              502
#define PEF_TO_TOD              503
#define PEF_TO_TEME             504
#define PEF_TO_PEF              505
#define PEF_TO_WGS84            506
#define PEF_TO_ITRF             506
#define PEF_TO_GEO              506
#define PEF_TO_GSE              507
#define PEF_TO_GSM              508
#define PEF_TO_SM               509
#define PEF_TO_EDMAG            510
#define PEF_TO_CDMAG            511
#define PEF_TO_GSE2000          512

#define WGS84_TO_EME2000        601
#define WGS84_TO_ICRF2000       601
#define WGS84_TO_GEI2000        601
#define WGS84_TO_MOD            602
#define WGS84_TO_TOD            603
#define WGS84_TO_TEME           604
#define WGS84_TO_PEF            605
#define WGS84_TO_WGS84          606
#define WGS84_TO_ITRF           606
#define WGS84_TO_GEO            606
#define WGS84_TO_GSE            607
#define WGS84_TO_GSM            608
#define WGS84_TO_SM             609
#define WGS84_TO_EDMAG          610
#define WGS84_TO_CDMAG          611
#define WGS84_TO_GSE2000        612

#define ITRF_TO_EME2000         601
#define ITRF_TO_ICRF2000        601
#define ITRF_TO_GEI2000         601
#define ITRF_TO_MOD             602
#define ITRF_TO_TOD             603
#define ITRF_TO_TEME            604
#define ITRF_TO_PEF             605
#define ITRF_TO_WGS84           606
#define ITRF_TO_ITRF            606
#define ITRF_TO_GEO             606
#define ITRF_TO_GSE             607
#define ITRF_TO_GSM             608
#define ITRF_TO_SM              609
#define ITRF_TO_EDMAG           610
#define ITRF_TO_CDMAG           611
#define ITRF_TO_GSE2000         612

#define GEO_TO_EME2000          601
#define GEO_TO_ICRF2000         601
#define GEO_TO_GEI2000          601
#define GEO_TO_MOD              602
#define GEO_TO_TOD              603
#define GEO_TO_TEME             604
#define GEO_TO_PEF              605
#define GEO_TO_WGS84            606
#define GEO_TO_ITRF             606
#define GEO_TO_GEO              606
#define GEO_TO_GSE              607
#define GEO_TO_GSM              608
#define GEO_TO_SM               609
#define GEO_TO_EDMAG            610
#define GEO_TO_CDMAG            611
#define GEO_TO_GSE2000          612

#define GSE_TO_EME2000          701
#define GSE_TO_ICRF2000         701
#define GSE_TO_GEI2000          701
#define GSE_TO_MOD              702
#define GSE_TO_TOD              703
#define GSE_TO_TEME             704
#define GSE_TO_PEF              705
#define GSE_TO_WGS84            706
#define GSE_TO_ITRF             706
#define GSE_TO_GEO              706
#define GSE_TO_GSE              707
#define GSE_TO_GSM              708
#define GSE_TO_SM               709
#define GSE_TO_EDMAG            710
#define GSE_TO_CDMAG            711
#define GSE_TO_GSE2000          712

#define GSM_TO_EME2000          801
#define GSM_TO_ICRF2000         801
#define GSM_TO_GEI2000          801
#define GSM_TO_MOD              802
#define GSM_TO_TOD              803
#define GSM_TO_TEME             804
#define GSM_TO_PEF              805
#define GSM_TO_WGS84            806
#define GSM_TO_ITRF             806
#define GSM_TO_GEO              806
#define GSM_TO_GSE              807
#define GSM_TO_GSM              808
#define GSM_TO_SM               809
#define GSM_TO_EDMAG            810
#define GSM_TO_CDMAG            811
#define GSM_TO_GSE2000          812

#define SM_TO_EME2000           901
#define SM_TO_ICRF2000          901
#define SM_TO_GEI2000           901
#define SM_TO_MOD               902
#define SM_TO_TOD               903
#define SM_TO_TEME              904
#define SM_TO_PEF               905
#define SM_TO_WGS84             906
#define SM_TO_ITRF              906
#define SM_TO_GEO               906
#define SM_TO_GSE               907
#define SM_TO_GSM               908
#define SM_TO_SM                909
#define SM_TO_EDMAG             910
#define SM_TO_CDMAG             911
#define SM_TO_GSE2000           912

#define EDMAG_TO_EME2000       1001
#define EDMAG_TO_ICRF2000      1001
#define EDMAG_TO_GEI2000       1001
#define EDMAG_TO_MOD           1002
#define EDMAG_TO_TOD           1003
#define EDMAG_TO_TEME          1004
#define EDMAG_TO_PEF           1005
#define EDMAG_TO_WGS84         1006
#define EDMAG_TO_ITRF          1006
#define EDMAG_TO_GEO           1006
#define EDMAG_TO_GSE           1007
#define EDMAG_TO_GSM           1008
#define EDMAG_TO_SM            1009
#define EDMAG_TO_EDMAG         1010
#define EDMAG_TO_CDMAG         1011
#define EDMAG_TO_GSE2000       1012

#define CDMAG_TO_EME2000       1101
#define CDMAG_TO_ICRF2000      1101
#define CDMAG_TO_GEI2000       1101
#define CDMAG_TO_MOD           1102
#define CDMAG_TO_TOD           1103
#define CDMAG_TO_TEME          1104
#define CDMAG_TO_PEF           1105
#define CDMAG_TO_WGS84         1106
#define CDMAG_TO_ITRF          1106
#define CDMAG_TO_GEO           1106
#define CDMAG_TO_GSE           1107
#define CDMAG_TO_GSM           1108
#define CDMAG_TO_SM            1109
#define CDMAG_TO_EDMAG         1110
#define CDMAG_TO_CDMAG         1111
#define CDMAG_TO_GSE2000       1112

#define GSE2000_TO_EME2000     1201
#define GSE2000_TO_ICRF2000    1201
#define GSE2000_TO_GEI2000     1201
#define GSE2000_TO_MOD         1202
#define GSE2000_TO_TOD         1203
#define GSE2000_TO_TEME        1204
#define GSE2000_TO_PEF         1205
#define GSE2000_TO_WGS84       1206
#define GSE2000_TO_ITRF        1206
#define GSE2000_TO_GEO         1206
#define GSE2000_TO_GSE         1207
#define GSE2000_TO_GSM         1208
#define GSE2000_TO_SM          1209
#define GSE2000_TO_EDMAG       1210
#define GSE2000_TO_CDMAG       1211
#define GSE2000_TO_GSE2000     1212

typedef struct Lgm_LeapSeconds {

    int         nLeapSecondDates;   //!< Number of leap second dates.

    long int    *LeapSecondDates;   //!< Array for holdin the Dates on which leap seconds were added

    double      *LeapSecondJDs;     //!< Array for holdin the Julian Dates on which leap seconds were added

    double      *LeapSeconds;       //!< The actual number of leap seconds that  went into effect on the given date


} Lgm_LeapSeconds;


typedef struct Lgm_DateTime {

    long int    Date;       //!< In basic ISO format (YYYYMMDD or YYYYDDD) Represented as a single long int

    int         Year;       //!< 4-digit year

    int         Month;      //!< [1-12]

    int         Day;        //!< Day Of Month [1-31]

    int         Doy;        //!< Day Of Year [1-31]

    double      Time;       //!< Decimal value of time in hours

    int         Hour;       //!< Hours [0-23]

    int         Minute;     //!< Minutes [0-59]

    double      Second;     //!< Seconds [0-60] (the 60 accommodates leap seconds)

    int         Week;       //!< ISO Week number [1-53]

    int         wYear;      //!< ISO Year associated with the ISO Week Number (can be different from Year)

    int         Dow;        //!< ISO Day Of Week number [1-7]

    char        DowStr[10]; //!< ISO Day Of Week number [1-7]

    double      fYear;      //!< Decimal year (e.g. 2004.2345)

    double      JD;         //!< Julian Date

    double      T;          //!< Julian Centuries since J2000
                            //!< for this time system

    double      DaySeconds; //!< Number of seconds in the day.

    int         TZD_sgn;    //!< Sign of Time zone offset
    int         TZD_hh;     //!< Time zone offset hours
    int         TZD_mm;     //!< Time zone offset minutes

    int         TimeSystem; //!< e.g. LGM_UTC, LGM_UT1, LGM_TAI, LGM_GPS, LGM_TT, LGM_TDB, LGM_TCG, etc..

} Lgm_DateTime;



typedef struct Lgm_CTrans {

    int               Verbose;

    Lgm_LeapSeconds   l;            //!< Structure containing Leap Second Info

    Lgm_JPLephemInfo *jpl;          //!< Structure containing JPL Ephem info
    int              jpl_initialized;  //!< Flag to indicate jpl has been initialized

    Lgm_DateTime     UT1;           //!< A corrected version of UT0.
                                    /**< A corrected version of UT0.
                                     *   UT is the mean solar time at Greenwich.
                                     *   UT0 is a version of UT that uses data
                                     *   from many different ground stations.
                                     *   UT1 is a version of UT0 in which
                                     *   corrections for polar motion have been
                                     *   made so that time is independant of
                                     *   observing location. There is also a UT2,
                                     *   but we wont use UT0 or UT2 here.
                                     *   Units: Decimal hours
                                     */

    Lgm_DateTime     UTC;           //!< Universal Time Coordinated.
                                    /**< Universal Time Coordinated.
                                     *   Most commonly used time system. Derived
                                     *   from atomic time. It is maintained to be
                                     *   within +/- 0.9s of UT1 (via addition or
                                     *   subtraction(?) of leap seconds).
                                     *   Units: Decimal hours
                                     */

    double           DUT1;          //!< Difference between UT1 and UTC ( DUT1 = UT1 - UTC.)
                                    /**<
                                     *   Difference between UT1 and UTC.
                                     *          DUT1 = UT1 - UTC.
                                     *   This is monitored and reported as part
                                     *   of the Earth Orientation Parameters
                                     *   (EOP). Can be predicted a short time
                                     *   into the future, but definitive values
                                     *   only available retrospectively.  We set
                                     *   this value to 0.0 by default. Thus in
                                     *   the absence of EOP data, we assume its
                                     *   initial value.
                                     *   Units: Decimal seconds
                                     */

    double          LOD;            //!< Length Of Day (LOD). Seconds in a day - 86400.
                                    /**<
                                     *   Length Of Day (LOD). Its the amount of extra
                                     *   time in seconds that the current day has. Not
                                     *   predictable. Part of EOP values.
                                     */

    Lgm_DateTime    TAI;            //!< International Atomic Time. (TAI = UTC + DAT).
                                    /**<
                                      *  International Atomic Time.
                                      *     TAI = UTC + DAT
                                      */

    Lgm_DateTime    GPS;            //!< GPS time (GPS = TAI - 19s).
                                    /**<
                                      *  GPS time
                                      *     GPS = TAI - 19s
                                      */

    double          DAT;            //!< Difference between TAI and UTC. (DAT = TAI - UTC.)
                                    /**<
                                      *  Difference between TAI and UTC.
                                      *      DAT = TAI - UTC
                                      *  DAT is essentially the number of leap seconds
                                      *  and are an integral number of whole seconds.
                                      *  Units: Decimal seconds.
                                      */

    Lgm_DateTime    TT;             //!< Terestrial Time (TT).
                                    /**<
                                      *  Terestrial Time (TT).
                                      *  Essentially the same thing as
                                      *  "Terrestrial Dynamical Time (TDT) or
                                      *  Ephmeris Time (ET). Its defined to be,
                                      *       TT = TAI + 32.184s
                                      */

    Lgm_DateTime    TDB;            //!< Barycentric Dynamical Time.
                                    /**<
                                      *  Barycentric Dynamical Time.
                                      *  Not used here.
                                      */

    Lgm_DateTime    TCG;            //!< Geocentric Coordinate Time.
                                    /**<
                                      *  Geocentric Coordinate Time.
                                      *  Not used here.
                                      */


    double      gmst;               //!< Greenwich Mean Sidereal Time.
                                    /**<
                                      *  Greenwich Mean Sidereal Time
                                      *  Units: in radians
                                      */

    double      gast;               //!< Greenwich Apparent Sidereal Time.
                                    /**<
                                      *  Greenwich Apparent Sidereal Time
                                      *  Units: in radians
                                      */


    double      xp, yp;             //!< EOP Pole wander parameters.
                                    /**<
                                      *  Pole wander parameters.
                                      *  part of EOP data.
                                      *  Units: radians
                                      */

    double      epsilon;            //!< Mean Obliquity of the Ecliptic.
                                    /**<
                                     *   Mean Obliquity of the Ecliptic
                                     *   (in radians)
                                     */

    double      epsilon_true;       /**<
                                     *  True Obliquity of the Ecliptic
                                     *  \f$\epsilon_{true} = \epsilon + dEps\f$
                                     * (in radians)
                                     */

    double      eccentricity;       /**< Eccentricity of Earth-Sun orbit */
    double      mean_anomaly;       /**< Mean anomaly of Earth-Sun orbit */
    double      true_anomaly;       /**< Mean anomaly of Earth-Sun orbit */

    double      lambda_sun;         /**< Ecliptic Long. of Sun (in radians) */
    double      beta_sun;           /**< Ecliptic elevation of Sun (in radians) */
    double      earth_sun_dist;     /**< Earth-Sun distance (in units of Earth radii) */
    double      RA_sun;             /**< Right Ascention of Sun (in degrees) */
    double      DEC_sun;            /**< Declination of Sun (in degrees) */

    Lgm_Vector  Sun;                /**< direction of Sun in MOD system (unit vector) */
    Lgm_Vector  SunJ2000;           /**< direction of Sun in GEI2000 system (unit vector) */
    Lgm_Vector  EcPole;             /**< direction of Ecliptic Pole in MOD system (unit vector) */
    double      psi;                /**< Geodipole tilt angle, \f$\psi\f$ (in radians) */
    double      sin_psi;            /**< \f$\sin(\psi)\f$ */
    double      cos_psi;            /**< \f$\cos(\psi)\f$ */
    double      tan_psi;            /**< \f$\tan(\psi)\f$ */
    Lgm_Vector  MoonJ2000;           /**< direction of Moon in GEI2000 system (unit vector) */
    double      RA_moon;            /**< Right Ascension of Moon (in degrees) */
    double      DEC_moon;           /**< Declination of Moon (in degrees) */
    double      MoonPhase;          /**< The Phase of the Moon (in days) */
    double      EarthMoonDistance;  /**< Distance between the Earth and Moon (in earth-radii) */

    int         ephModel;           /**< Model to use for Sun and Moon positions */
    int         pnModel;            /**< Precession-nutation model */

    /*
     *  The following are various important parameters derived from
     *  the IGRF field. Note that these are the basis for defining
     *  Mag coord systems. That's why they are here and not somewhere else...
     */
    double      M_cd;               /**< centered  dipole Magnetic moment of epoch. (nT Re^3)*/
    double      M_cd_McIlwain;      /**<  magnetic dipole moment used by McIlwain to compute L. Sometimes want to use this for consistency? */
    double      M_cd_2010;          /**<  magnetic dipole moment for epoch 2010-01-01T00:00:00Z. Used by us to compute L*. */
    double      CD_gcolat;          /**< Geographic colat of centered dipole axis (deg.)  */
    double      CD_glon;            /**< Geographic long. of centered dipole axis (deg.)  */
    double      ED_x0;              /**<  x-comp of dipole displacement from center. Used in eccentric dipole field. */
    double      ED_y0;              /**<  y-comp of dipole displacement from center. Used in eccentric dipole field. */
    double      ED_z0;              /**<  z-comp of dipole displacement from center. Used in eccentric dipole field. */


    /*
     * Precession Angles
     */
    double      Zeta;               /**< Precession angle, \f$\zeta\f$ */
    double      Theta;              /**< Precession angle, \f$\theta\f$ */
    double      Zee;                /**< Precession angle, \f$z\f$ */

    /*
     * Some things for nutation reduction
     */
    int         nNutationTerms; /**< number of terms to usek in the dPsi/dEps Nutation series. */
    double      dPsi;
    double      dEps;
    double      dPsiCosEps;
    double      dPsiSinEps;
    double      ddPsi;      /**< radians additional corrections to dPsi -- part of EOP data */
    double      ddEps;      /**< radians additional corrections to dEps -- part of EOP data */
    double      EQ_Eq;      /**< Equation of the equinoxes. */
    double      OmegaMoon;  /**< Ascending node of Moon. */
    double      dX;         /**< for IUA-2000A reduction (not used yet) */
    double      dY;         /**< for IUA-2000A reduction (not used yet) */



    /*
     * Transformation matrices between various coord systems
     */
    double      Agei_to_mod[3][3];
    double      Amod_to_gei[3][3];
    double      Amod_to_tod[3][3];
    double      Atod_to_mod[3][3];
    double      Ateme_to_pef[3][3];
    double      Apef_to_teme[3][3];
    double      Apef_to_tod[3][3];
    double      Atod_to_pef[3][3];
    double      Awgs84_to_pef[3][3];
    double      Apef_to_wgs84[3][3];
    double      Agse_to_mod[3][3];
    double      Amod_to_gse[3][3];
    double      Agse2000_to_gei[3][3];
    double      Agei_to_gse2000[3][3];
    double      Asm_to_gsm[3][3];
    double      Agsm_to_sm[3][3];
    double      Agsm_to_mod[3][3];
    double      Amod_to_gsm[3][3];
    double      Agsm_to_gse[3][3];
    double      Agse_to_gsm[3][3];
    double      Awgs84_to_mod[3][3];
    double      Amod_to_wgs84[3][3];
    double      Awgs84_to_gei[3][3];
    double      Agei_to_wgs84[3][3];
    double      Agsm_to_wgs84[3][3];
    double      Awgs84_to_gsm[3][3];
    double      Awgs84_to_cdmag[3][3];
    double      Acdmag_to_wgs84[3][3];



    /*
     *  These variables are needed to make IGRF Calls reentrant/thread-safe.
     */
    int         Lgm_IGRF_FirstCall;
    double      Lgm_IGRF_OldYear;
    double      Lgm_IGRF_g[14][14];
    double      Lgm_IGRF_h[14][14];
    double      Lgm_IGRF_R[14][14];
    double      Lgm_IGRF_K[14][14];
    double      Lgm_IGRF_S[14][14];
    double      Lgm_IGRF_TwoNm1_Over_NmM[14][14];
    double      Lgm_IGRF_NpMm1_Over_NmM[14][14];
    double      Lgm_IGRF_SqrtNM1[14][14];
    double      Lgm_IGRF_SqrtNM2[14][14];



} Lgm_CTrans;


void        Lgm_free_ctrans( Lgm_CTrans *c );
void        Lgm_free_ctrans_children( Lgm_CTrans *c );
Lgm_CTrans  *Lgm_init_ctrans( int );
Lgm_CTrans  *Lgm_CopyCTrans( Lgm_CTrans *s );
void        Lgm_ctransDefaults(Lgm_CTrans *, int);


void        Lgm_Radec_to_Cart( double, double, Lgm_Vector * );
double      Lgm_angle2pi( double );
double      Lgm_angle360( double );

int         Lgm_LeapYear( int );

long int    Lgm_JDN( int Year, int Month, int Day );

double      Lgm_JD( int Year, int Month, int Day, double Time, int TimeSystem, Lgm_CTrans *c );
long int    Lgm_JD_to_Date(double jd, int *ny, int *nm, int *nd, double *UT);
double      Lgm_Date_to_JD( long int Date, double UT, Lgm_CTrans *c );
void        Lgm_jd_to_ymdh ( double JD, long int *Date, int *year, int *month, int *day, double *UT );
int         Lgm_DayOfYear(int year, int month, int day, Lgm_CTrans *c);

double      Lgm_MJD( int Year, int Month, int Day, double Time, int TimeSystem, Lgm_CTrans *c );
long int    Lgm_MJD_to_Date(double mjd, int *ny, int *nm, int *nd, double *UT);
void        Lgm_mjd_to_ymdh ( double MJD, long int *Date, int *year, int *month, int *day, double *UT );

double      Lgm_hour24( double );
double      Lgm_kepler( double, double );
double      Lgm_Dipole_Tilt(long int date, double UTC);
void        Lgm_Set_CTrans_Options( int ephModel, int pnModel, Lgm_CTrans *c );
void        Lgm_Set_Coord_Transforms( long int, double, Lgm_CTrans * );
void        Lgm_ComputeSun( Lgm_CTrans *c );
void        Lgm_ComputeMoon( Lgm_CTrans *c );
void        Lgm_Convert_Coords(Lgm_Vector *, Lgm_Vector *, int, Lgm_CTrans * );
int         Lgm_IsValidDate( long int );
int         Lgm_Doy( long int, int *, int *, int *, int * );
void        Lgm_UT_to_hmsms( double UT, int *HH, int *MM, int *SS, int *MilliSec );
void        Lgm_UT_to_HMS( double UT, int *HH, int *MM, int *SS );
void        Lgm_UT_to_HMSd( double UT, int *sgn, int *HH, int *MM, double *SS );
void        Lgm_D_to_DMS( double D, int *DD, int *MM, int *SS );
void        Lgm_D_to_DMSd( double D, int *sgn, int *DD, int *MM, double *SS );
void        Lgm_Print_HMS( double d );
void        Lgm_Print_HMSd( double d );
void        Lgm_Print_HMSdp( double d, int UnicodeHMS, int p );
void        Lgm_Print_DMS( double d );
void        Lgm_Print_DMSd( double d );
double      Lgm_GetCurrentJD( Lgm_CTrans *c );
double      Lgm_GetCurrentMJD( Lgm_CTrans *c );
void        Lgm_SunPosition( double T, double *l, double *r, double *b );
void        Lgm_GLATLON_TO_CDMLATLONMLT( double GLAT, double GLON, double *MLAT, double *MLON, double *MLT, Lgm_CTrans *c );
void        Lgm_GLATLON_TO_EDMLATLONMLT( double GLAT, double GLON, double *MLAT, double *MLON, double *MLT, Lgm_CTrans *c );
void        Lgm_CDMAG_to_R_MLAT_MLON_MLT( Lgm_Vector *u, double *R, double *MLAT, double *MLON, double *MLT, Lgm_CTrans *c );
void        Lgm_R_MLAT_MLT_to_CDMAG( double R, double MLAT, double MLT, Lgm_Vector *u, Lgm_CTrans *c );
void        Lgm_EDMAG_to_R_MLAT_MLON_MLT( Lgm_Vector *u, double *R, double *MLAT, double *MLON, double *MLT, Lgm_CTrans *c );
void        Lgm_R_MLAT_MLT_to_EDMAG( double R, double MLAT, double MLT, Lgm_Vector *u, Lgm_CTrans *c );
void        Lgm_WGS84_to_GEOD( Lgm_Vector *uin, double *GeodLat, double *GeodLong, double *GeodHieght );
void        Lgm_WGS84_to_GeodHeight( Lgm_Vector *uin, double *GeodHieght );
void        Lgm_GEOD_to_WGS84( double GeodLat, double GeodLong, double GeodHieght, Lgm_Vector *v );
void        Lgm_Nutation( double T_TT, double nTerms, double *dPSi, double *dEps );
int         IsoTimeStringToDateTime( char *TimeString, Lgm_DateTime *d, Lgm_CTrans *c );
int         MonthStrToNum( char *str );
char       *Lgm_StrToLower( char *str, int nmax );
char       *Lgm_StrToUpper( char *str, int nmax );

/*
 * Leap Second Related  Routines
 */
Lgm_DateTime *Lgm_DateTime_Create( int Year, int Month, int Day, double Time, int TimeSystem, Lgm_CTrans *c );
void          Lgm_DateTime_Destroy( Lgm_DateTime *d );
int           Lgm_Make_UTC( long int Date, double Time, Lgm_DateTime *UTC, Lgm_CTrans *c );
int           Lgm_LoadLeapSeconds( Lgm_CTrans *c );
double        Lgm_GetLeapSeconds( double JD, Lgm_CTrans *c );
int           Lgm_IsLeapSecondDay( long int Date, double *SecondsInDay, Lgm_CTrans *c );
void          Lgm_UTC_to_TAI( Lgm_DateTime *UTC, Lgm_DateTime *TAI, Lgm_CTrans *c );
void          Lgm_TAI_to_UTC( Lgm_DateTime *TAI, Lgm_DateTime *UTC, Lgm_CTrans *c );
void          Lgm_TT_to_TAI( Lgm_DateTime *TT, Lgm_DateTime *TAI, Lgm_CTrans *c );
void          Lgm_TAI_to_TT( Lgm_DateTime *TAI, Lgm_DateTime *TT, Lgm_CTrans *c );
void          Lgm_TT_to_TDB( Lgm_DateTime *TT, Lgm_DateTime *TDB, Lgm_CTrans *c );
void          Lgm_TDB_to_TT( Lgm_DateTime *TDB, Lgm_DateTime *TT, Lgm_CTrans *c );
void          Lgm_UTC_to_TT( Lgm_DateTime *UTC, Lgm_DateTime *TT, Lgm_CTrans *c );
void          Lgm_TT_to_UTC( Lgm_DateTime *TT, Lgm_DateTime *UTC, Lgm_CTrans *c );
void          Lgm_TAI_to_GPS( Lgm_DateTime *TAI, Lgm_DateTime *GPS, Lgm_CTrans *c );
void          Lgm_GPS_to_TAI( Lgm_DateTime *GPS, Lgm_DateTime *TAI, Lgm_CTrans *c );
void          Lgm_UTC_to_GPS( Lgm_DateTime *UTC, Lgm_DateTime *GPS, Lgm_CTrans *c );
void          Lgm_GPS_to_UTC( Lgm_DateTime *GPS, Lgm_DateTime *UTC, Lgm_CTrans *c );
void          Lgm_Print_DateTime( Lgm_DateTime *DT, int Style, int p );
void          Lgm_DateTimeToString( char *Str, Lgm_DateTime *DT, int Style, int p );
void          Lgm_Print_SimpleTime( Lgm_DateTime *DT, int p, char * );
//int         Lgm_DayofWeek( int, int, int, char *, Lgm_CTrans *c );
int           Lgm_DayOfWeek( int Year, int Month, int Day, char *dowstr );
long int      Lgm_JDNofWeek1( int Year );
int           Lgm_MaxWeekNumber( int Year );
int           Lgm_ISO_WeekNumber( int Year, int Month, int Day, int *ISO_WeekYear );
void          Lgm_ISO_YearWeekDow_to_Date( int ISO_WeekYear, int Week, int Dow, long int *Date, int *Year, int *Month, int *Day );
double        Lgm_RemapTime( double Time, double SecondsInADay );

double        Lgm_GPS_to_GpsSeconds( Lgm_DateTime *GPS );
void          Lgm_GpsSeconds_to_GPS( double GpsSeconds, Lgm_DateTime *GPS );
void          Lgm_GpsSeconds_to_UTC( double GpsSeconds, Lgm_DateTime *UTC, Lgm_CTrans *c );
double        Lgm_UTC_to_GpsSeconds( Lgm_DateTime *UTC, Lgm_CTrans *c );

double        Lgm_TAI_to_TaiSeconds( Lgm_DateTime *TAI );
void          Lgm_TaiSeconds_to_GPS( double TaiSeconds, Lgm_DateTime *TAI );
void          Lgm_TaiSeconds_to_UTC( double TaiSeconds, Lgm_DateTime *UTC, Lgm_CTrans *c );
double        Lgm_UTC_to_TaiSeconds( Lgm_DateTime *UTC, Lgm_CTrans *c );


double        Lgm_TDB_to_TdbSecSinceJ2000( Lgm_DateTime *TDB );
void          Lgm_TdbSecSinceJ2000_to_TDB( double TdbSeconds, Lgm_DateTime *TDB );
void          Lgm_TdbSecSinceJ2000_to_UTC( double TdbSeconds, Lgm_DateTime *UTC, Lgm_CTrans *c );
double        Lgm_UTC_to_TdbSecSinceJ2000( Lgm_DateTime *UTC, Lgm_CTrans *c );


double        TAISecondsSinceJ2000( double UTCDaysSinceJ2000, Lgm_CTrans *c );
double        UTCDaysSinceJ2000( double TAI, Lgm_CTrans *c );

double        Lgm_TDBSecSinceJ2000( Lgm_DateTime *UTC, Lgm_CTrans *c );

void          Lgm_JD_to_DateTime( double JD, Lgm_DateTime *UTC, Lgm_CTrans *c );
void          Lgm_MJD_to_DateTime( double MJD, Lgm_DateTime *UTC, Lgm_CTrans *c );


double        Lgm_CdipMirrorLat( double Alpha0 );

int           Lgm_EarthEclipse( Lgm_Vector *Psc, Lgm_CTrans *c );
int           Lgm_MoonEclipse( Lgm_Vector *Psc, Lgm_CTrans *c );



double  Lgm_UTC_to_TTSeconds( Lgm_DateTime *UTC, Lgm_CTrans *c );
double  Lgm_TTSecSinceJ2000( Lgm_DateTime *UTC, Lgm_CTrans *c );
double  Lgm_TT_to_TTSecSinceJ2000( Lgm_DateTime *TT ) ;
void    Lgm_TTSecSinceJ2000_to_TT( double TTSeconds, Lgm_DateTime *TT ) ;
void    Lgm_TTSecSinceJ2000_to_UTC( double TTSeconds, Lgm_DateTime *UTC, Lgm_CTrans *c ) ;
double  Lgm_UTC_to_TTSecSinceJ2000( Lgm_DateTime *UTC, Lgm_CTrans *c ) ;



//double      Lgm_UTC_to_TT( double
//double      Lgm_TT_to_UTC( double JD, double TT, Lgm_LeapSeconds *l );
//double      Lgm_TT_to_TDB( double JD, double TT, Lgm_LeapSeconds *l );
//double      Lgm_TDB_to_TT( double JD, double TDB, Lgm_LeapSeconds *l );
//double      Lgm_UTC_to_TDB( double JD, double UTC, Lgm_LeapSeconds *l );
//double      Lgm_TDB_to_UTC( double JD, double TDB, Lgm_LeapSeconds *l );



/*
 *  Internal B-field prototypes
 */
void    Lgm_B_igrf_ctrans(Lgm_Vector *, Lgm_Vector *, Lgm_CTrans *);
void    Lgm_B_cdip_ctrans(Lgm_Vector *, Lgm_Vector *, Lgm_CTrans *);
void    Lgm_B_edip_ctrans(Lgm_Vector *v, Lgm_Vector *B, Lgm_CTrans *c);
void    Lgm_B_JensenCain1960_ctrans(Lgm_Vector *v, Lgm_Vector *B, Lgm_CTrans *c);


/*
 *  IGRF prototypes
 */
double  Lgm_Factorial( int );
void    Lgm_InitIGRF( double g[14][14], double h[14][14], int N, int Flag, Lgm_CTrans *c );
void    Lgm_InitPnm( double ct, double st, double R[14][14], double P[14][14], double dP[14][14], int N, Lgm_CTrans *c );
void    Lgm_InitTrigmp( double , double , double *, double *, int );
void    Lgm_PolFunInt( double *, double *, int, double, double *, double * );
void    Lgm_RatFunInt( double *, double *, int, double, double *, double * );
void    Lgm_IGRF( Lgm_Vector *, Lgm_Vector *, Lgm_CTrans * );
void    _Lgm_IGRF( Lgm_Vector *, Lgm_Vector *, Lgm_CTrans * );
void    _Lgm_IGRF2( Lgm_Vector *, Lgm_Vector *, Lgm_CTrans * );
void    _Lgm_IGRF3( Lgm_Vector *, Lgm_Vector *, Lgm_CTrans * );
void    _Lgm_IGRF4( Lgm_Vector *, Lgm_Vector *, Lgm_CTrans * );

void   Lgm_InitdPnm( double P[14][14], double dP[14][14], int N, Lgm_CTrans *c );
void   Lgm_InitSqrtFuncs( double SqrtNM1[14][14], double SqrtNM2[14][14], int N );
void   Lgm_InitK( double K[14][14], int N );
void   Lgm_InitS( double S[14][14], int N );


/*
 *  Additional Sph Harmonic models
 */
void    Lgm_InitSphHarm( int Model, double g[14][14], double h[14][14], int N, int Flag, Lgm_CTrans *c );
void    Lgm_JensenCain1960( Lgm_Vector *v, Lgm_Vector *B, Lgm_CTrans *c );
void    _Lgm_JensenCain1960( Lgm_Vector *v, Lgm_Vector *B, Lgm_CTrans *c );


void Lgm_Terminator( double GLON, double *GLAT_N, int *nRoots, double alpha, Lgm_CTrans *c );


#endif
