#include <EXTERN.h>
#include <perl.h>
#define ISO_YYYYMMDDTHHMMSS 1
#define ISO_YYYYWwwDTHHMMSS 2
#define ISO_YYYYMMDDTHHMM   3
#define ISO_YYYYWwwDTHHMM   4
#define ISO_YYYYMMDD        5
#define ISO_YYYYWwwD        6
#define ISO_YYYYMM          7
#define ISO_YYYYWww         8

static PerlInterpreter *my_perl;

ParseTimeString( char *TimeString ) {

    int     Year, wYear, Month, Day, Hours, Minutes, TZD_sgn, TZD_hh, TZD_mm, Week, DayOfWeek;
    int     ISOFormat;
    double  Seconds;
    STRLEN  n_a;
    char    *embedding[] = { "", "-e", "0" };
    char    Str[6000];


    my_perl = perl_alloc();
    perl_construct( my_perl );

    perl_parse(my_perl, NULL, 3, embedding, NULL);
    perl_run(my_perl);

    /* 
     * ISO 8601 Time Formats are strictly formated, but have a number of variations
     * as follows;
     *
     *     YYYY
     *     YYYY-MM
     *     YYYY-MM-DD (hypens optional)
     *     YYYY-MM-DDThh:mmTZD (hypens and colons optional)
     *     YYYY-MM-DDThh:mm:ssTZD (hypens and colons optional)
     *     YYYY-Www (hypens can be omitted)
     *     YYYY-Www-D (hypens can be omitted)
     *
     * where,
     * 
     *     YYYY = four-digit year
     *     MM   = two-digit month (01=January, etc.)
     *     DD   = two-digit day of month (01 through 31)
     *     hh   = two digits of hour (00 through 23) (am/pm NOT allowed)
     *     mm   = two digits of minute (00 through 59)
     *     ss   = two digits of second (00 through 59)
     *     s    = one or more digits representing a decimal fraction of a second
     *     TZD  = time zone designator (Z or +hh:mm or -hh:mm) (colons optional)
     *     ww   = ww is the week of the year
     *     D    = ww is day of week 
     * 
     * The hypens and colons are optional above except in the YYYY-MM format.
     *
     * YYYY can be positive or negative but only if users agree on it (so by
     * default its not in the standard). Here we assume it to be positive.
     *
     * Also, its not clear that if you omit one hyphen you have to omit them
     * all.  Or that omitting one colon means you have to omit them all. Check
     * the standard on this, but here I allow you to omit any number of them.
     * E.g., 2009-06-21T05:12:34.01Z and 2009-0621T0512:34.01Z should be
     * treated the same here.
     *
     * Combined Date and Time formats must have the T separating them, but you
     * can also have separate date and time formats and we recognize both E.g.,
     * 2009-06-21T05:12:34.01Z versus 2009-06-21 05:12:34.01Z
     */
    sprintf( Str, " $str = '%s';\n\
                    $str =~ m/^\\s*(\\S*\\s*\\S*)\\s*$/; $str = $1; # Strip leading and trailing whitespace\n\
                    $str =~ s/\\s{2,}/ /g; # Compress out extraneous whitespace\n\
                    print \"str = $str\\n\";\n\
                    my $ISOFormat = 0; my $Year = -9999; my $Month = -9999; my $Day = -9999; my $Hour = -9999; my $Minute = -9999; my $Second = -9999; my $TZD = -9999;\n\
                    if ( $str =~ m/^(\\d{4})-?(\\d{2})-?(\\d{2})[T ](\\d{2}):?(\\d{2}):?(\\d{2}[\\.]?\\d*)([+-]\\d{2}:?\\d{2}|Z)/ ) { #  YYYY[-]MM[-]DDThh[:]mm[:]ss[.ssssss](TZD or Z)\n\
                        $Year = $1; $Month = $2; $Day = $3; $Hour = $4; $Minute = $5; $Second = $6; $TZD = $7; $ISOFormat = 1;\n\
                    } elsif ( $str =~ m/^(\\d{4})-?W(\\d{2})-?(\\d{1})[T ](\\d{2}):?(\\d{2}):?(\\d{2}[\\.]?\\d*)([+-]\\d{2}:?\\d{2}|Z)/ ) { #  YYYY[-]Www[-]DThh[:]mm[:]ss[.ssssss](TZD or Z)\n\
                        $Year = $1; $Week = $2; $DayOfWeek = $3; $Hour = $4; $Minute = $5; $Second = $6; $TZD = $7; $ISOFormat = 2;\n\
                    } elsif ( $str =~ m/^(\\d{4})-?(\\d{2})-?(\\d{2})[T ](\\d{2}):?(\\d{2})([+-]\\d{2}:?\\d{2}|Z)$/ ) { #  YYYY[-]MM[-]DDThh[:]mm(TZD or Z)\n\
                        $Year = $1; $Month = $2; $Day = $3; $Hour = $4; $Minute = $5; $TZD = $6; $ISOFormat = 3;\n\
                    } elsif ( $str =~ m/^(\\d{4})-?W(\\d{2})-?(\\d{1})[T ](\\d{2}):?(\\d{2})([+-]\\d{2}:?\\d{2}|Z)$/ ) { #  YYYY[-]Www[-]DThh[:]mm(TZD or Z)\n\
                        $Year = $1; $Week = $2; $DayOfWeek = $3; $Hour = $4; $Minute = $5; $TZD = $6; $ISOFormat = 4;\n\
                    } elsif ( $str =~ m/^(\\d{4})-?(\\d{2})-?(\\d{2})$/ ) { # YYYY[-]MM[-]DD\n\
                        $Year = $1; $Month = $2; $Day = $3;  $ISOFormat = 5;\n\
                    } elsif ( $str =~ m/^(\\d{4})-?W(\\d{2})-?(\\d{1})$/ ) { # YYYY[-]Www[-]D\n\
                        $Year = $1; $Week = $2; $DayOfWeek = $3;  $ISOFormat = 6;\n\
                    } elsif ( $str =~ m/^(\\d{4})-(\\d{2})$/ ) { # YYYY-MM\n\
                        $Year = $1; $Month = $2;  $ISOFormat = 7;\n\
                    } elsif ( $str =~ m/^(\\d{4})-W(\\d{2})$/ ) { # YYYY-Www\n\
                        $Year = $1; $Week = $2;  $ISOFormat = 8;\n\
                    } elsif ( $str =~ m/^(\\d{4})$/ ) { # YYYY\n\
                        $Year = $1; $ISOFormat = 9;\n\
                    }\n\
                    if ( $TZD =~ m/^Z$/ ) {\n\
                        $TZD_sgn = +1; $TZD_hh = 0; $TZD_mm = 0;\n\
                    } elsif ( $TZD =~ m/^([+-])(\\d{2}):?(\\d{2})$/ ) {\n\
                        $TZD_sgn = ($1 == '+') ? +1 : -1; $TZD_hh = $2; $TZD_mm = $3;\n\
                    } elsif ( $TZD =~ m/^([+-])(\\d{2})$/ ) {\n\
                        $TZD_sgn = ($1 == '+') ? +1 : -1; $TZD_hh = $2; $TZD_mm = 0;\n\
                    } else { \n\
                        $TZD_sgn = +1; $TZD_hh = 0; $TZD_mm = 0;\n\
                    }\n\
                    $oYear = $Year; $oMonth = $Month; $oDay = $Day; $oHour = $Hour; $oMinute = $Minute; $oSecond = $Second; $oTZD = $TZD; $oTZD_sgn = $TZD_sgn; $oTZD_hh = $TZD_hh; $oTZD_mm = $TZD_mm; $oWeek = $Week; $oDayOfWeek = $DayOfWeek; $oISOFormat = $ISOFormat;\n", TimeString);
    

    eval_pv( Str, TRUE );
    printf( "Year      = %d\n",  Year = SvIV(get_sv("oYear", FALSE)) );
    printf( "Month     = %d\n",  Month = SvIV(get_sv("oMonth", FALSE)) );
    printf( "Day       = %d\n",  Day = SvIV(get_sv("oDay", FALSE)) );
    printf( "Hours     = %d\n",  Hours = SvIV(get_sv("oHour", FALSE)) );
    printf( "Minutes   = %d\n",  Minutes = SvIV(get_sv("oMinute", FALSE)) );
    printf( "Seconds   = %lf\n", Seconds = atof(SvPV(get_sv("oSecond", FALSE), n_a)));
    printf( "TZD_sgn   = %d\n",  TZD_sgn = SvIV(get_sv("oTZD_sgn", FALSE)) );
    printf( "TZD_hh    = %d\n",  TZD_hh = SvIV(get_sv("oTZD_hh", FALSE)) );
    printf( "TZD_mm    = %d\n",  TZD_mm = SvIV(get_sv("oTZD_mm", FALSE)) );
    printf( "Week      = %d\n",  Week = SvIV(get_sv("oWeek", FALSE)) );
    printf( "DayOfWeek = %d\n",  DayOfWeek = SvIV(get_sv("oDayOfWeek", FALSE)) );
    printf( "ISOFormat = %d\n",  ISOFormat = SvIV(get_sv("oISOFormat", FALSE)) );

    perl_destruct(my_perl);
    perl_free(my_perl);


    /*
     * Test for sane numbers
     */
    if        ( (ISOFormat == ISO_YYYYMMDDTHHMMSS) || (ISOFormat == ISO_YYYYMMDDTHHMM) ) {

        if (ISOFormat == ISO_YYYYMMDDTHHMM) Seconds = 0.0; // assume seconds are zero

        if ( ( Year < 0 ) || ( Year > 9999 ) ) {
            printf( "Year Out of Range!\n", Year );
            return( -1 );
        }
        if ( ( Month < 1 ) || ( Month > 12 ) ) {
            printf( "Month Out of Range!\n", Month );
            return( -1 );
        }
        if ( ( Day < 1 ) || ( Day > 31 ) ) {
            printf( "Day Out of Range!\n", Day );
            return( -1 );
        }
        if ( ( Hours < 0 ) || ( Hours > 23 ) ) {
            printf( "Hours Out of Range!\n", Hours );
            return( -1 );
        }
        if ( ( Minutes < 0 ) || ( Minutes > 59 ) ) {
            printf( "Minutes Out of Range!\n", Minutes );
            return( -1 );
        }
        if ( ( Seconds < 0 ) || ( Seconds > 60 ) ) { // Yes seconds can be 60 to accomodate leap seconds
            printf( "Seconds Out of Range!\n", Seconds ); 
            return( -1 );
        }

    } else if ( (ISOFormat = ISO_YYYYWwwDTHHMM) || (ISOFormat = ISO_YYYYWwwDTHHMMSS) ) {

        if (ISOFormat = ISO_YYYYWwwDTHHMMSS) Seconds = 0.0; // assume seconds are zero
        wYear = Year;

        if ( ( wYear < 0 ) || ( wYear > 9999 ) ) {
            printf( "wYear Out of Range!\n", wYear );
            return( -1 );
        }
        if ( ( Week < 1 ) || ( Week > 53 ) ) {
            printf( "Week Out of Range!\n", Week );
            return( -1 );
        }
        if ( ( DayOfWeek < 1 ) || ( DayOfWeek > 7 ) ) {
            printf( "DayOfWeek Out of Range!\n", DayOfWeek );
            return( -1 );
        }
        if ( ( Hours < 0 ) || ( Hours > 23 ) ) {
            printf( "Hours Out of Range!\n", Hours );
            return( -1 );
        }
        if ( ( Minutes < 0 ) || ( Minutes > 59 ) ) {
            printf( "Minutes Out of Range!\n", Minutes );
            return( -1 );
        }
        if ( ( Seconds < 0 ) || ( Seconds > 60 ) ) { // Yes seconds can be 60 to accomodate leap seconds
            printf( "Seconds Out of Range!\n", Seconds ); 
            return( -1 );
        }

    } else if ( (ISOFormat = ISO_YYYYMM) || (ISOFormat = ISO_YYYYMMDD) ) {
        
        if ( ISOFormat = ISO_YYYYMM ) Day = 1; //assume day is 1st day on Month

        if ( ( Year < 0 ) || ( Year > 9999 ) ) {
            printf( "Year Out of Range!\n", Year );
            return( -1 );
        }
        if ( ( Month < 1 ) || ( Month > 12 ) ) {
            printf( "Month Out of Range!\n", Month );
            return( -1 );
        }
        if ( ( Day < 1 ) || ( Day > 31 ) ) {
            printf( "Day Out of Range!\n", Day );
            return( -1 );
        }

    } else if ( (ISOFormat = ISO_YYYYWww) || (ISOFormat = ISO_YYYYWwwD) ) {

        if ( ISOFormat = ISO_YYYYWww ) DayOfWeek = 1; //assume day is 1st day on Week
        wYear = Year;

        if ( ( wYear < 0 ) || ( wYear > 9999 ) ) {
            printf( "Year Out of Range!\n", wYear );
            return( -1 );
        }
        if ( ( Week < 1 ) || ( Week > 53 ) ) {
            printf( "Week Out of Range!\n", Week );
            return( -1 );
        }
        if ( ( DayOfWeek < 1 ) || ( DayOfWeek > 7 ) ) {
            printf( "DayOfWeek Out of Range!\n", DayOfWeek );
            return( -1 );
        }

    }




}



main (int argc, char **argv, char **env) {

    char    TimeString[128];

    sprintf(TimeString, "2009-06-21T05:12:34.123456789Z");
    sprintf(TimeString, "2009-06-21T05:1234.01Z");
    sprintf(TimeString, "2009-06-21                           0512:34.01Z");
    sprintf(TimeString, "2009-0621T99:12:34.123456789-03:34");
    sprintf(TimeString, "2009-W52-7T2312:34.01Z");
    printf("TimeString = %s\n", TimeString );


    ParseTimeString( TimeString );



}
