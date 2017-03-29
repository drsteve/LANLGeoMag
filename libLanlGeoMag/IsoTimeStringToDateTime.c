/*! \file IsoTimeStringToDateTime.c
 *
 *  \brief A Perl-based ISO 8601 time-parser. Results stored in a Lgm_DateTime structure.
 *
 *
 *
 *  \author M.G. Henderson
 *  \date   2010-2011
 *
 *
 *
 */

#if ENABLE_PERL
#include <EXTERN.h>
#include <perl.h>
#endif /* ENABLE_PERL */

#include "Lgm/Lgm_CTrans.h"
#include "time.h"


/*
 *  Allowed formats for PERL parser
 */
#define PERL_ERROR           0       
#define PERL_YYYYMMDDTHHMMSS 1       // 4-digit years
#define PERL_YYMMDDTHHMMSS   2       // 2-digit years
#define PERL_YYYYDDDTHHMMSS  3       // 4-digit years
#define PERL_YYDDDTHHMMSS    4       // 2-digit years
#define PERL_YYYYWwwDTHHMMSS 5       // 4-digit years
#define PERL_YYWwwDTHHMMSS   6       // 2-digit years
#define PERL_YYYYMMDDTHHMM   7       // 4-digit years
#define PERL_YYMMDDTHHMM     8       // 2-digit years
#define PERL_YYYYWwwDTHHMM   9       // 4-digit years
#define PERL_YYWwwDTHHMM     10      // 2-digit years
#define PERL_YYYYMMDD        11      // 4-digit years (assumes UTC == 0)
#define PERL_YYMMDD          12      // 2-digit years (assumes UTC == 0)
#define PERL_YYYYDDD         13      // 4-digit years (assumes UTC == 0)
#define PERL_YYDDD           14      // 2-digit years (assumes UTC == 0)
#define PERL_YYYYWwwD        15      // 4-digit years (assumes UTC == 0)
#define PERL_YYWwwD          16      // 2-digit years (assumes UTC == 0)
#define PERL_YYYYMM          17      // 4-digit years (YYYY-MM only YYYYMM not allowed) (assumes UTC == 0)
#define PERL_YYYYWww         19      // 4-digit years (assumes UTC == 0)
#define PERL_YYWww           20      // 2-digit years (assumes UTC == 0)
#define PERL_YYYY            21      // 4-digit year 
#define PERL_YY              22      // 2-digit year 

/*
 *  Allowed formats for C parser
 */
#define C_ERROR           -9999  
#define C_YYYYMMDDTHHMMSS 0       // YYYY-MM-DDThh:mm:ss
#define C_YYYYDDDTHHMMSS  1       // YYYY-DDDThh:mm:ss
#define C_YYYYWwwDTHHMMSS 2       // YYYY-Www-DThh:mm:ss
#define C_YYYYMMDD        3       // YYYY-MM-DDThh:mm:ss
#define C_YYYYDDD         4       // YYYY-DDDThh:mm:ss
#define C_YYYYWwwD        5       // YYYY-Www-DThh:mm:ss

/*
 * Notes.
 *
 *  In ISO 8601, you can delete the least significant fields all the way to the
 *  2-digit century. For example,
 *
 *      YYYY-MM-DDTHH:MM:SS.SSSSSS+HH:MM
 *      YYYY-MM-DDTHH:MM:SS.SSSSSS+HH
 *      YYYY-MM-DDTHH:MM:SS.SSSSSS
 *      YYYY-MM-DDTHH:MM:SS
 *      YYYY-MM-DDTHH:MM
 *      YYYY-MM-DDTHH
 *      YYYY-MM-DD
 *      YYYY-MM
 *      YYYY
 *      YY
 *
 *  The 2 Y's in the last case are not the usual 2-digit years, they are the
 *  century (e.g. if YYYY is 1996, then the last case above would be YY=19).
 *
 *  The separators can be omitted for all of these except for the YYYY-MM case.
 *  That one always has to have a separator or it will get confused with other
 *  formats.
 *
 *  The original ISO 8601 allowed "truncated representations" meaning that
 *  2-digit years could be used. But apparently the latest version od ISO 8601
 *  got rid of that.
 *  
 *  The C-parser is more restrictive than the PERL parser. The C parser
 *  requires all punctuation ('-' and ':'). It also does not allow two-digit
 *  years (YY) or incomplete time designations (e.g. no seconds, or no months)
 *
 */

#if ENABLE_PERL
static PerlInterpreter *my_perl;
#endif /* ENABLE_PERL */

int IsoTimeStringToDateTime( char *TimeString, Lgm_DateTime *d, Lgm_CTrans *c ) {

#if ENABLE_PERL

    long int    Date;
    double      Offset, Time;
    int         tyear, tday, tmonth, sgn;
    int         Year, wYear, Month, Day, Hours, Minutes, TZD_sgn, TZD_hh, TZD_mm, Week, DayOfWeek, DayOfYear;
    int         ISOFormat, TZDError=FALSE, MaxWeek, InvalidDate=FALSE, IsLeapSecondDay;
    double      Seconds;
    STRLEN      n_a;
    char        *embedding[] = { "", "-e", "0" };
    char        Str[6000];


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

    int AllowTruncatedReps = 1;

    sprintf( Str, " $str = '%s';\n"
                  " $str =~ s/^\\s+|\\s+$//g;\n"
                  " $str =~ s/\\s{2,}/ /g; # Compress out extraneous whitespace\n"
                  " #print \"str = [$str]\\n\";\n"
                  " my $ISOFormat = 0; my $Year = -9999; my $DayOfYear = -9999; my $Month = -9999; my $Day = -9999; my $Hour = -9999; my $Minute = -9999; my $Second = -9999; my $TZD = 'Z';\n", TimeString );
            

                /*
                 *  PERL_YYYYMMDDTHHMMSS  
                 *
                 *  Match ISO strings of the form:
                 *
                 *           YYYY[-]MM[-]DDThh[:]mm[:]ss[.ssssss](TZD or Z or nothing)
                 *
                 *  Examples: 
                 *        2010-10-31T12:34:56.789-06:00
                 *        20101031T123456.789-0600
                 */
                strcat( Str, 
                    " if ( $str =~ m/^(\\d{4})-?(\\d{2})-?(\\d{2})[T ](\\d{2}):?(\\d{2}):?(\\d{2}[\\.]?\\d*)([+-].*|Z?)/ ) { \n\
                         $Year = $1; $Month = $2; $Day = $3; $Hour = $4; $Minute = $5; $Second = $6; $TZD = $7; $ISOFormat = 1;\n");


                /*
                 *  PERL_YYYYDDDTHHMMSS
                 *
                 *  Match ISO strings of the form:
                 *
                 *           YYYY[-]DDDThh[:]mm[:]ss[.ssssss](TZD or Z or nothing)
                 *
                 *  Examples: 
                 *        2010-297T12:34:56.789-06:00
                 *        1996-297T01:54:12.4534789+05:00
                 *        2010297T123456.789-06:00
                 *        1996297T015412.4534789+05:00
                 */
                strcat( Str, 
                    " } elsif ( $str =~ m/^(\\d{4})-?(\\d{3})[T ](\\d{2}):?(\\d{2}):?(\\d{2}[\\.]?\\d*)([+-].*|Z?)/ ) { \n\
                            $Year = $1; $DayOfYear = $2; $Hour = $3; $Minute = $4; $Second = $5; $TZD = $6; $ISOFormat = 3;\n"
                );


                  /*
                   *  PERL_YYYYWwwDTHHMMSS
                   *
                   *  Match ISO strings of the form:
                   *
                   *        YYYY[-]Www[-]DThh[:]mm[:]ss[.ssssss](TZD or Z or nothing)
                   *
                   *  Examples:  
                   *       1986-W13-2T09:45:32+00
                   *       1986W132T094532Z
                   */
                strcat( Str, 
                  " } elsif ( $str =~ m/^(\\d{4})-?W(\\d{2})-?(\\d{1})[T ](\\d{2}):?(\\d{2}):?(\\d{2}[\\.]?\\d*)([+-].*|Z?)/ ) { \n\
                        $Year = $1; $Week = $2; $DayOfWeek = $3; $Hour = $4; $Minute = $5; $Second = $6; $TZD = $7; $ISOFormat = 5;\n"
                );


                  /*
                   *  PERL_YYYYMMDDTHHMM
                   *
                   *  Match ISO strings of the form:
                   *
                   *        YYYY[-]MM[-]DDThh[:]mm(TZD or Z or nothing)
                   *
                   *  Examples: 
                   *        2010-10-31T12:34-06:00
                   *        20101031T1234-0600
                   */
                strcat( Str, 
                  " } elsif ( $str =~ m/^(\\d{4})-?(\\d{2})-?(\\d{2})[T ](\\d{2}):?(\\d{2})([+-].*|Z?)$/ ) { \n\
                        $Year = $1; $Month = $2; $Day = $3; $Hour = $4; $Minute = $5; $Second = 0.0; $TZD = $6; $ISOFormat = 7;\n"
                );


                  /*
                   *  PERL_YYYYWwwDTHHMM
                   *
                   *  Match ISO strings of the form:
                   *
                   *        YYYY[-]Www[-]DThh[:]mm(TZD or Z or nothing)
                   *
                   *  Examples: 
                   *       1986-W13-2T09:45+00
                   *       1986W132T0945Z
                   */
                strcat( Str, 
                  " } elsif ( $str =~ m/^(\\d{4})-?W(\\d{2})-?(\\d{1})[T ](\\d{2}):?(\\d{2})([+-].*|Z?)$/ ) { \n\
                        $Year = $1; $Week = $2; $DayOfWeek = $3; $Hour = $4; $Minute = $5; $Second = 0.0; $TZD = $6; $ISOFormat = 9;\n"
                );


                  /*
                   *  PERL_YYYYMMDD
                   *
                   *  Match ISO strings of the form:
                   *
                   *        YYYY[-]MM[-]DD
                   *
                   *  Examples: 
                   *        1976-04-23
                   *        19760423
                   *        2006-04-23
                   *        20060423
                   */
                strcat( Str, 
                  " } elsif ( $str =~ m/^(\\d{4})-?(\\d{2})-?(\\d{2})$/ ) {\n\
                        $Year = $1; $Month = $2; $Day = $3;  $Hour = $Minute = 0; $Second = 0.0; $ISOFormat = 11;\n"
                );


                  /*
                   *  PERL_YYYYMM
                   *
                   *  Match ISO strings of the form:
                   *
                   *        YYYY-MM
                   *
                   *  Examples: 
                   *        1976-12
                   *        197612 <- this is not allowed in ISO 8601
                   * (Keep this one before PERL_YYMMDD)
                   */
                strcat( Str, 
                  " } elsif ( $str =~ m/^(\\d{4})-(\\d{2})$/ ) { \n\
                        $Year = $1; $Month = $2;  $Day = 1; $Hour = $Minute = 0; $Second = 0.0;$ISOFormat = 17;\n"
                );


                  /*
                   *  PERL_YYYYDDD
                   *
                   *  Match ISO strings of the form:
                   *
                   *        YYYY[-]DDD
                   *
                   *  Examples: 
                   *        1976-142
                   *        1976142
                   *        2006-142
                   *        2006142
                   */
                strcat( Str, 
                  " } elsif ( $str =~ m/^(\\d{4})-?(\\d{3})$/ ) {\n\
                        $Year = $1; $DayOfYear = $2; $Hour = $Minute = 0; $Second = 0.0;$ISOFormat = 13;\n"
                );

                  /*
                   *  PERL_YYYYWwwD
                   *
                   *  Match ISO strings of the form:
                   *
                   *        YYYY[-]Www[-]D
                   *
                   *  Examples: 
                   *        1976-W12-6
                   *        1976W126
                   */
                strcat( Str, 
                  " } elsif ( $str =~ m/^(\\d{4})-?W(\\d{2})-?(\\d{1})$/ ) {\n\
                        $Year = $1; $Week = $2; $DayOfWeek = $3;  $Hour = $Minute = 0; $Second = 0.0; $ISOFormat = 15;\n"
                );


                  /*
                   *  PERL_YYYYWww
                   *  
                   *  Match ISO strings of the form:
                   *  
                   *        YYYY-Www
                   *
                   *  Examples: 
                   *        1976-W32
                   */
                strcat( Str, 
                  " } elsif ( $str =~ m/^(\\d{4})-W(\\d{2})$/ ) { \n\
                        $Year = $1; $Week = $2; $DayOfWeek = 1; $Hour = $Minute = 0; $Second = 0.0; $ISOFormat = 19;\n"
                );








                /*
                 *  PERL_YYYY
                 *  
                 *  Match ISO strings of the form:
                 *  
                 *        YYYY
                 *
                 *  Examples: 
                 *        1976
                 */
                strcat( Str, 
                  " } elsif ( $str =~ m/^(\\d{4})$/ ) { \n\
                        $Year = $1; $Month = 1; $Day = 1; $Hour = $Minute = 0; $Second = 0.0; $ISOFormat = 21;\n"
                );


                /*
                 *  PERL_YY
                 *  
                 *  Match ISO strings of the form:
                 *  
                 *        YY
                 * 
                 *  Note: this is the century not a two-digit year.
                 *
                 *  Examples: 
                 *        19
                 */
                strcat( Str, 
                        " } elsif ( $str =~ m/^(\\d{2})$/ ) { \n\
                                $Year = $1*100; $Month = 1; $Day = 1; $Hour = $Minute = 0; $Second = 0.0; $ISOFormat = 22;\n"
                );





                if ( AllowTruncatedReps ) {

                    /*
                     *  PERL_YYMMDDTHHMMSS
                     *
                     *  Match ISO strings of the form:
                     *
                     *           YY[-]MM[-]DDThh[:]mm[:]ss[.ssssss](TZD or Z or nothing)
                     *
                     *  Examples: 
                     *        10-10-31T12:34:56.789-06:00
                     *        96-03-21T01:54:12.4534789+05:00
                     */
                    strcat( Str,    
                        " } elsif ( $str =~ m/^(\\d{2})-?(\\d{2})-?(\\d{2})[T ](\\d{2}):?(\\d{2}):?(\\d{2}[\\.]?\\d*)([+-].*|Z?)/ ) { \n\
                                $Year = ($1 > 50) ? $1+1900:$1+2000; $Month = $2; $Day = $3; $Hour = $4; $Minute = $5; $Second = $6; $TZD = $7; $ISOFormat = 2;\n"
                    );

                    /*
                     *  PERL_YYDDDTHHMMSS
                     *
                     *  Match ISO strings of the form:
                     *
                     *           YY[-]DDDThh[:]mm[:]ss[.ssssss](TZD or Z or nothing)
                     *
                     *  Examples: 
                     *        10-297T12:34:56.789-06:00
                     *        96-297T01:54:12.4534789+05:00
                     *        10297T123456.789-06:00
                     *        96297T015412.4534789+05:00
                     */
                    strcat( Str, 
                        " } elsif ( $str =~ m/^(\\d{2})-?(\\d{3})[T ](\\d{2}):?(\\d{2}):?(\\d{2}[\\.]?\\d*)([+-].*|Z?)/ ) { \n\
                                $Year = ($1 > 50) ? $1+1900:$1+2000; $DayOfYear = $2; $Hour = $3; $Minute = $4; $Second = $5; $TZD = $6; $ISOFormat = 4;\n"
                    );

                    /*
                     *  PERL_YYWwwDTHHMMSS
                     *
                     *  Match ISO strings of the form:
                     *
                     *    YY[-]Www[-]DThh[:]mm[:]ss[.ssssss](TZD or Z or nothing)
                     *
                     *  Examples:  
                     *       86-W13-2T09:45:32
                     *       86W132T094532
                     */
                    strcat( Str, 
                        " } elsif ( $str =~ m/^(\\d{2})-?W(\\d{2})-?(\\d{1})[T ](\\d{2}):?(\\d{2}):?(\\d{2}[\\.]?\\d*)([+-].*|Z?)/ ) { \n\
                            $Year = ($1 > 50) ? $1+1900:$1+2000; $Week = $2; $DayOfWeek = $3; $Hour = $4; $Minute = $5; $Second = $6; $TZD = $7; $ISOFormat = 6;\n"
                    );

                    /*
                     *  PERL_YYMMDDTHHMM
                     *
                     *  Match ISO strings of the form:
                     *
                     *        YY[-]MM[-]DDThh[:]mm(TZD or Z or nothing)
                     *
                     *  Examples: 
                     *        10-10-31T12:34-06:00
                     *        101031T1234-0600
                     *        86-10-31T12:34-06:00
                     *        661031T1234-0600
                     */
                    strcat( Str, 
                        " } elsif ( $str =~ m/^(\\d{2})-?(\\d{2})-?(\\d{2})[T ](\\d{2}):?(\\d{2})([+-].*|Z?)$/ ) { \n\
                            $Year = ($1 > 50) ? $1+1900:$1+2000; $Month = $2; $Day = $3; $Hour = $4; $Minute = $5; $Second = 0.0; $TZD = $6; $ISOFormat = 8;\n"
                    );

                    /*
                     *  PERL_YYWwwDTHHMM
                     *
                     *  Match ISO strings of the form:
                     *
                     *        YY[-]Www[-]DThh[:]mm(TZD or Z or nothing)
                     *
                     *  Examples: 
                     *       86-W13-2T09:45+00
                     *       86W132T0945Z
                     *       09-W13-2T09:45+00
                     *       09W132T0945Z
                     */
                    strcat( Str, 
                        " } elsif ( $str =~ m/^(\\d{2})-?W(\\d{2})-?(\\d{1})[T ](\\d{2}):?(\\d{2})([+-].*|Z?)$/ ) { \n\
                                    $Year = ($1 > 50) ? $1+1900:$1+2000; $Week = $2; $DayOfWeek = $3; $Hour = $4; $Minute = $5; $Second = 0.0; $TZD = $6; $ISOFormat = 10;\n"
                    );
                    /*
                     *  PERL_YYMMDD
                     *
                     *  Match ISO strings of the form:
                     *
                     *        YY[-]MM[-]DD
                     *
                     *  Examples: 
                     *        76-04-23
                     *        760423
                     *        06-04-23
                     *        060423
                     */
                    strcat( Str, 
                        " } elsif ( $str =~ m/^(\\d{2})-?(\\d{2})-?(\\d{2})$/ ) {\n\
                                    $Year = ($1 > 50) ? $1+1900:$1+2000; $Month = $2; $Day = $3; $Hour = $Minute = 0; $Second = 0.0; $ISOFormat = 12;\n"
                    );

                    /*
                     *  PERL_YYDDD
                     *
                     *  Match ISO strings of the form:
                     *
                     *        YY[-]DDD
                     *
                     *  Examples: 
                     *        76-142
                     *        76142
                     *        06-142
                     *        06142
                     */
                    strcat( Str, 
                        " } elsif ( $str =~ m/^(\\d{2})-?(\\d{3})$/ ) {\n\
                                    $Year = ($1 > 50) ? $1+1900:$1+2000; $DayOfYear = $2; $Hour = $Minute = 0; $Second = 0.0;$ISOFormat = 14;\n"
                    );

                    /*
                     *  PERL_YYWwwD
                     *
                     *  Match ISO strings of the form:
                     *
                     *        YY[-]Www[-]D
                     *
                     *  Examples: 
                     *        76-W12-6
                     *        76W126
                     */
                    strcat( Str, 
                        " } elsif ( $str =~ m/^(\\d{2})-?W(\\d{2})-?(\\d{1})$/ ) {\n\
                                    $Year = ($1 > 50) ? $1+1900:$1+2000; $Week = $2; $DayOfWeek = $3; $Hour = $Minute = 0; $Second = 0.0; $ISOFormat = 16;\n"
                    );


                    /*
                     *  PERL_YYMM
                     *
                     *  Match ISO strings of the form:
                     *
                     *    YY-MM
                     *  
                     *  The hyphen is mandatory here.
                     *  Examples: 
                     *        76-12
                     */
                    strcat( Str, 
                        " } elsif ( $str =~ m/^(\\d{2})-(\\d{2})$/ ) { \n\
                                    $Year = ($1 > 50) ? $1+1900:$1+2000; $Day = 1; $Month = $2;  $Hour = $Minute = 0; $Second = 0.0;$ISOFormat = 18;\n"
                    );

                    /*
                     *  PERL_YYWww
                     *  
                     *  Match ISO strings of the form:
                     *  
                     *        YY-Www
                     *
                     *  Examples: 
                     *        76-W32
                     */
                    strcat( Str, 
                        " } elsif ( $str =~ m/^(\\d{2})-W(\\d{2})$/ ) { \n\
                                    $Year = ($1 > 50) ? $1+1900:$1+2000; $Week = $2; $DayOfWeek = 1; $Hour = $Minute = 0; $Second = 0.0; $ISOFormat = 20;\n"
                    );

                }





                strcat( Str, 
                  " }\n\
                    $TZDError = 0;\n\
                    if ( $TZD =~ m/^$/ ) {\n\
                        $TZD_sgn = +1; $TZD_hh = 0; $TZD_mm = 0;\n\
                    } elsif ( $TZD =~ m/^Z$/ ) {\n\
                        $TZD_sgn = +1; $TZD_hh = 0; $TZD_mm = 0;\n\
                    } elsif ( $TZD =~ m/^([+-])(\\d{2}):?(\\d{2})$/ ) {\n\
                        $TZD_sgn = ($1 eq '+') ? +1 : -1; $TZD_hh = $2; $TZD_mm = $3;\n\
                    } elsif ( $TZD =~ m/^([+-])(\\d{2})$/ ) {\n\
                        $TZD_sgn = ($1 eq '+') ? +1 : -1; $TZD_hh = $2; $TZD_mm = 0;\n\
                    } else { \n\
                        $TZDError = 1;\n\
                    }\n\
                    $oYear = $Year; $oMonth = $Month; $oDay = $Day; $oHour = $Hour; $oMinute = $Minute; $oSecond = $Second; $oTZD = $TZD; $oTZD_sgn = $TZD_sgn; $oTZD_hh = $TZD_hh; $oTZD_mm = $TZD_mm; $oWeek = $Week; $oDayOfWeek = $DayOfWeek; $oDayOfYear = $DayOfYear; $oTZDError = $TZDError; $oISOFormat = $ISOFormat;\n");
    

    eval_pv( Str, TRUE );

    Year      = SvIV(get_sv("oYear", FALSE));
    DayOfYear = SvIV(get_sv("oDayOfYear", FALSE));
    Month     = SvIV(get_sv("oMonth", FALSE));
    Day       = SvIV(get_sv("oDay", FALSE));
    Hours     = SvIV(get_sv("oHour", FALSE));
    Minutes   = SvIV(get_sv("oMinute", FALSE));
    Seconds   = atof(SvPV(get_sv("oSecond", FALSE), n_a));
    TZD_sgn   = SvIV(get_sv("oTZD_sgn", FALSE));
    TZD_hh    = SvIV(get_sv("oTZD_hh", FALSE));
    TZD_mm    = SvIV(get_sv("oTZD_mm", FALSE));
    Week      = SvIV(get_sv("oWeek", FALSE));
    DayOfWeek = SvIV(get_sv("oDayOfWeek", FALSE));
    TZDError  = SvIV(get_sv("oTZDError", FALSE));
    ISOFormat = SvIV(get_sv("oISOFormat", FALSE));

    perl_destruct(my_perl);
    perl_free(my_perl);


    // Assume time is UTC for now...
    if ( (ISOFormat == PERL_YYYYWwwDTHHMMSS) || (ISOFormat == PERL_YYWwwDTHHMMSS)
            || (ISOFormat == PERL_YYYYWwwDTHHMM) || (ISOFormat == PERL_YYWwwDTHHMM)
            || (ISOFormat == PERL_YYYYWwwD) || (ISOFormat == PERL_YYWwwD)
            || (ISOFormat == PERL_YYYYWww) || (ISOFormat == PERL_YYWww) ) {

        d->wYear  = Year;
        d->Week   = Week;
        d->Dow    = DayOfWeek;
        Lgm_ISO_YearWeekDow_to_Date( d->wYear, Week, DayOfWeek, &Date, &Year, &Month, &Day );
        Lgm_Doy( Date, &Year, &Month, &Day, &DayOfYear );
        Lgm_DayOfWeek( Year, Month, Day, d->DowStr );
        MaxWeek = Lgm_MaxWeekNumber( d->wYear );
        if ( Week > MaxWeek ) {
            printf("IsoTimeStringToDateTime: Invalid Date. Week number %d does not exist in year %d\n", d->Week, d->wYear );
            InvalidDate = TRUE;
        } 

    } else {

        d->Dow        = Lgm_DayOfWeek( Year, Month, Day, d->DowStr );
        d->Week       = Lgm_ISO_WeekNumber( Year, Month, Day, &d->wYear );

        if ( DayOfYear > 0 ){
            Date = Year*1000 + DayOfYear;
        } else {
            Date = Year*10000 + Month*100 + Day;
        }
        if ( Lgm_IsValidDate( Date ) ) {
            Lgm_Doy( Date, &Year, &Month, &Day, &DayOfYear );
        } else {
            printf( "Date is invalid\n" );
            InvalidDate = TRUE;
        }

    }

#else /* ENABLE_PERL */

  // Use a simplified C ISO date/time parser if not using PERL

  struct tm timeinfo;
  char *ret_val = NULL;
  const char *format[] = {
  "%Y-%m-%dT%H:%M:%S%Z",
  "%Y-%jT%H:%M:%S%Z",
  "%Y-W%U-%uT%H:%M:%S%Z",
  "%Y-%m-%d",
  "%Y-%j",                       
  "%Y-W%U-%u"};

  int i;
  int N = 6;                  // Number of date formats supported

  long int Date;              // 8-digit date in the format YYYYMMDD
  double Offset;              // Time zone designator offset from UTC
  double Time;                // Decimal hours during day [0 to 23.999...]
  int tyear, tmonth, tday;    // Temporary year, month, and day; only used for Lgm_DOY return
  int TZDError=FALSE;         // Time zone designator parsing error flag
  int MaxWeek;                // Maximum number of weeks in the specified year
  int InvalidDate=FALSE;      // Invalid date flag
  int IsLeapSecondDay;        // Flag to indicate if the day has a leap second included

  int ISOFormat = -9999;      // ISO format type (defined above)
  int Year = -9999;           // 4-digit year
  int Month = -9999;          // 2-digit month [1-12]
  int Day = -9999;            // 2-digit day of month [1-31]
  int DayOfYear = -9999;      // 3-digit day of year [1-366]
  int DayOfWeek = -9999;      // 1-digit day of week [1-7]
  int Week = -9999;           // 2-digit week number [1-53]
  int Hours = -9999;          // 2-digit hour [0-23]
  int Minutes = -9999;        // 2-digit minute [0-59]
  double Seconds = -9999;     // Decimal seconds
  char *TZD_sgn_str;          // Time zone designator sign (+/-/Z)
  int TZD_sgn = 0;            // Time zone designator integer (1, -1)
  int TZD_hh = 0;             // Time zone offset hours
  int TZD_mm = 0;             // Time zone offset minutes

  int DateOnlyFlag = 0;       // 1 if format is date-only; 0 otherwise.

  int scan_val;
  Lgm_DateTime *d1;
  d1 = calloc(1, sizeof(Lgm_DateTime));

  // Use strptime only to determine if the string matches the expected format
  // strptime() is not used for collecting the time structure because
  // it does not handle fractional seconds and has problems with time zones
  // strptime() has the advantage that it will accept formats with or without
  // timezone information
  
  // Loop over formats and find the proper format
  for(i=0; i<N; i++) {
    
    // Check for format
    ret_val = strptime(TimeString, format[i], &timeinfo);
    if (ret_val != NULL) {

      // Define the ISO format flag
      ISOFormat = i;
      break;
    }
    
  }

  // Check for unsupported partial formats
  // If the ISO format includes a T and any subsequent
  // time information but matches the date only
  // ISO formats then set the ISOFormat to the error code
  if(ISOFormat > 2) {
    if(strstr(TimeString, "T") != NULL) {
      ISOFormat = -9999;
    }
  }

  // Use sscanf() is to collect the appropriate time information
  // based on the format
  switch(ISOFormat) {

  case -9999:  // No acceptable format found
    printf("Unable to parse ISO time string.\n");
    printf("No matching ISO format found\n");
    break;

    /*
     *  Match ISO strings of the form:
     *
     *           YYYY-MM-DDThh:mm:ss[.ssssss](TZD or Z or nothing)
     *
     *  Examples: 
     *        2010-10-31T12:34:56.789-06:00
     *        2010-10-31T12:34:56.789Z
     */
  case 0:
    scan_val = sscanf(TimeString, "%4d-%2d-%2dT%2d:%2d:%lf%m[Z+-]%2d:%2d",
		      &Year, &Month, &Day, &Hours, &Minutes, &Seconds,
		      &TZD_sgn_str, &TZD_hh, &TZD_mm);
    break;

    /*
     *  Match ISO strings of the form:
     *
     *           YYYY-DDDThh:mm:ss[.ssssss](TZD or Z or nothing)
     *
     *  Examples: 
     *        2010-297T12:34:56.789-06:00
     *        1996-297T01:54:12.4534789+05:00
     */
  case 1:
    scan_val = sscanf(TimeString, "%4d-%3dT%2d:%2d:%lf%m[Z+-]%2d:%2d",
		      &Year, &DayOfYear, &Hours, &Minutes, &Seconds,
		      &TZD_sgn_str, &TZD_hh, &TZD_mm);
    break;

    /*
     *  Match ISO strings of the form:
     *
     *        YYYY-Www-DThh:mm:ss[.ssssss](TZD or Z or nothing)
     *
     *  Examples:  
     *       1986-W13-2T09:45:32+00
     */
  case 2:
    scan_val = sscanf(TimeString, "%4d-W%2d-%1dT%2d:%2d:%lf%m[Z+-]%2d:%2d",
		      &Year, &Week, &DayOfWeek, &Hours, &Minutes, &Seconds,
		      &TZD_sgn_str, &TZD_hh, &TZD_mm);
    break;

    /*
     *  Match ISO strings of the form:
     *
     *        YYYY-MM-DD
     *
     *  Examples: 
     *        1976-04-23
     *        19760423
     *        2006-04-23
     *        20060423
     */
  case 3:
    DateOnlyFlag = 1;
    TZD_sgn_str = NULL;
    scan_val = sscanf(TimeString, "%4d-%2d-%2d",
		      &Year, &Month, &Day);
    break;
    
    /*
     *  Match ISO strings of the form:
     *
     *        YYYY-DDD
     *
     *  Examples: 
     *        1976-142
     *        2006-142
     */
  case 4:
    DateOnlyFlag = 1;
    TZD_sgn_str = NULL;
    scan_val = sscanf(TimeString, "%4d-%3d",
		      &Year, &DayOfYear);
    break;

    /*
     *  Match ISO strings of the form:
     *
     *        YYYY-Www-D
     *
     *  Examples: 
     *        1976-W12-6
     */
  case 5:
    DateOnlyFlag = 1;
    TZD_sgn_str = NULL;
    scan_val = sscanf(TimeString, "%4d-W%2d-%1d",
		      &Year, &Week, &DayOfWeek);
    break;
    
  }

    // Set time and timezone if only date is given
    if(DateOnlyFlag) {

      // Assume 00:00:00 UTC
      Hours = 0;
      Minutes = 0;
      Seconds = 0.0;

      // Set UTC timezone
      TZD_sgn = 0;
      TZD_hh = 0;
      TZD_mm = 0;

    }

    // Fill the timezone information it none is given
    if(TZD_sgn_str == NULL) {
      
      // Assume UTC
      TZD_sgn = 1;
      TZD_hh = 0;
      TZD_mm = 0;

    } else {

      // Convert string to integer
      if (!strcmp(TZD_sgn_str, "+")) {
	TZD_sgn = 1;
      } else if (!strcmp(TZD_sgn_str, "-")) {
	TZD_sgn = -1;
      } else if (!strcmp(TZD_sgn_str, "Z")) {
	TZD_sgn = 1;
      } else {
	TZDError = 1;	
      }

    }

    // Assume time is UTC for now...
    if ((ISOFormat == C_YYYYWwwDTHHMMSS) || (ISOFormat == C_YYYYWwwD)) {

      d->wYear = Year;
      d->Week = Week;
      d->Dow = DayOfWeek;
      Lgm_ISO_YearWeekDow_to_Date( d->wYear, Week, DayOfWeek, &Date, &Year, &Month, &Day );
      Lgm_Doy( Date, &Year, &Month, &Day, &DayOfYear );
      Lgm_DayOfWeek( Year, Month, Day, d->DowStr );
      MaxWeek = Lgm_MaxWeekNumber( d->wYear );
      if ( Week > MaxWeek ) {
	printf("IsoTimeStringToDateTime: Invalid Date.\n");
	printf("Week number %d does not exist in year %d\n", d->Week, d->wYear);
	InvalidDate = TRUE;
      } 

    } else {

      d->Dow = Lgm_DayOfWeek( Year, Month, Day, d->DowStr );
      d->Week = Lgm_ISO_WeekNumber( Year, Month, Day, &d->wYear );

      if ( DayOfYear > 0 ){
	Date = Year*1000 + DayOfYear;
      } else {
	Date = Year*10000 + Month*100 + Day;
      }
      if ( Lgm_IsValidDate( Date ) ) {
	Lgm_Doy( Date, &Year, &Month, &Day, &DayOfYear );
      } else {
	printf( "Date is invalid\n" );
	InvalidDate = TRUE;
      }

    }

#endif /* ENABLE_PERL */

    // determine time zone offset
    Offset = (double)TZD_hh + (double)TZD_mm/60.0;

    // set UTC. At this point this may not be in the range [0,24] -- this is what we want.
    Time = Hours + Minutes/60.0 + Seconds/3600.0 - TZD_sgn*Offset;

    // Stuff date and time into JD. (JD takes Time as it is even if its not in range [0,24])
    d->JD     = Lgm_JD( Year, Month, Day, Time, LGM_TIME_SYS_UTC, c );

    // back out the date and time. Note date could be different from what we
    // started with due to time not in range [0,24]
    d->Date   = Lgm_JD_to_Date( d->JD, &d->Year, &d->Month, &d->Day, &Time );

    // Redo the Week stuff, as dates may have changed...
    d->Dow    = Lgm_DayOfWeek( d->Year, d->Month, d->Day, d->DowStr );
    d->Week   = Lgm_ISO_WeekNumber( d->Year, d->Month, d->Day, &d->wYear );


    // set hours, minutes, seconds based on originally parsed vals. (Dont
    // decode the "Time" value returned from Lgm_JD_to_Date() -- its not as
    // precise.)
    d->Hour   = Hours-TZD_sgn*TZD_hh; 
    d->Minute = Minutes-TZD_sgn*TZD_mm;
    d->Second = Seconds;


    // correct minutes for roll-overs
    if ( d->Minute < 0 ) {
        d->Minute += 60; 
        --d->Hour;
    } else if ( d->Minute > 59 ) {
        d->Minute -= 60;
        ++d->Hour;
    }

    // correct hours for roll-overs
    if ( d->Hour <  0 ) {
        d->Hour += 24; 
    } else if ( d->Hour > 23 ) {
        d->Hour -= 24;
    }

    // reconstruct Time
    d->Time       = (double)d->Hour + (double)d->Minute/60.0 + d->Second/3600.0;

    // Set Doy -p- ignore other things returned
    Lgm_Doy( d->Date, &tyear, &tmonth, &tday, &d->Doy);

    // The time is UTC now
    d->TimeSystem = LGM_TIME_SYS_UTC;

    // save TZ info
    d->TZD_sgn    = TZD_sgn;
    d->TZD_hh     = TZD_hh;
    d->TZD_mm     = TZD_mm;

    // save # secs in a day
    IsLeapSecondDay = Lgm_IsLeapSecondDay( d->Date, &d->DaySeconds, c );

    // Julian centuries
    d->T          = (d->JD - 2451545.0)/36525.0;

    // decimal year
    d->fYear      = (double)d->Year + ((double)d->Doy - 1.0 + d->Time/24.0)/(365.0 + (double)Lgm_LeapYear(d->Year));
    
    /*
     * Test for sane numbers
     */
    if ( ( d->Year < 0 ) || ( d->Year > 9999 ) ) {
        printf( "Year Out of Range! (Year = %d)\n", d->Year );
        InvalidDate = TRUE;
    }
    if ( ( d->Doy < 1 ) || ( d->Doy > 365+Lgm_LeapYear(d->Year) ) ) {
        printf( "DayOfYear Out of Range! (DayOfYear = %d)\n", d->Doy );
        InvalidDate = TRUE;
    }
    if ( ( d->Month < 1 ) || ( d->Month > 12 ) ) {
        printf( "Month Out of Range! (Month = %d)\n", d->Month );
        InvalidDate = TRUE;
    }
    if ( ( d->Day < 1 ) || ( d->Day > 31 ) ) {
        printf( "Day Out of Range! (Day = %d)\n", d->Day );
        InvalidDate = TRUE;
    }
    if ( ( d->Hour < 0 ) || ( d->Hour > 23 ) ) {
        printf( "Hours Out of Range! (Hours = %d)\n", d->Hour );
        InvalidDate = TRUE;
    }
    if ( ( d->Minute < 0 ) || ( d->Minute > 59 ) ) {
        printf( "Minutes Out of Range! (Minutes = %d)\n", d->Minute );
        InvalidDate = TRUE;
    }
    if ( d->Second < 0.0 ) { 
        printf( "Seconds Out of Range! (Seconds = %lf)\n", d->Second ); 
        InvalidDate = TRUE;
    }
    if ( d->Second >= 60.0 ) { // Yes seconds can be 60 to accomodate leap seconds
        if ( IsLeapSecondDay && (d->Minute == 59) ) {
            if ( d->Second >= 61.0 ){
                // not sure we ever get here, because the date will have rolled over above the call to Lgm_IsLeapSecondDay()
                printf( "Seconds Out of Range! (Seconds = %lf) on a leap second day (must be less than 61.0)\n", d->Second ); 
                InvalidDate = TRUE;
            }
        } else {
            printf( "Seconds Out of Range! (Seconds = %lf)\n", d->Second ); 
            InvalidDate = TRUE;
        }
    }
    if ( ( d->wYear < 0 ) || ( d->wYear > 9999 ) ) {
        printf( "Year (%d) Out of Range!\n", d->wYear );
        InvalidDate = TRUE;
    }
    if ( ( d->Week < 1 ) || ( d->Week > 53 ) ) {
        printf( "Week (%d) Out of Range!\n", d->Week );
        InvalidDate = TRUE;
    }
    if ( ( d->Dow < 1 ) || ( d->Dow > 7 ) ) {
        printf( "DayOfWeek (%d) Out of Range!\n", d->Dow );
        InvalidDate = TRUE;
    }

    if (TZDError ) {
        printf( "Unrecognized time zone field!\n" );
        return(-1);
    }

    if ( InvalidDate ) {
        d->Date = -1;
        d->Time = 0.0;
        printf( "Invalid Date!\n" );
        return(-1);
    }


    return( 1 );

}


// For testing purposes

/* void test_iso_parse(char *TimeString) { */

/*   int     Flag; */
/*   Lgm_CTrans *c = Lgm_init_ctrans(0); */
/*   Lgm_DateTime *d; */
/*   d = calloc(1, sizeof(Lgm_DateTime)); */

/*   printf("TimeString = %s\n", TimeString); */
/*   Flag = IsoTimeStringToDateTime(TimeString, d, c); */
/*   printf("Flag = %d\n\n", Flag); */

/* } */

/* int main (int argc, char **argv, char **env) { */

/*    int     Flag; */
/*    char    TimeString[128]; */
/*    Lgm_CTrans *c = Lgm_init_ctrans(0); */
/*    Lgm_DateTime *d; */
/*    d = calloc(1, sizeof(Lgm_DateTime)); */

/*    test_iso_parse("2009-06-21T05:12:34.123456789Z"); */
/*    test_iso_parse("2009-06-21T05:12:34.123456789"); */
/*    test_iso_parse("2010-11-03T16:15:02-06:00"); */
/*    test_iso_parse("2010-11-03T16:15:02-06"); */
/*    test_iso_parse("   2010-11-03T16:15:02   "); */
/*    test_iso_parse("2010-11-03"); */
/*    test_iso_parse("2016-179"); */
/*    test_iso_parse("2016-179T13:01:03.234"); */
/*    test_iso_parse("2016-W03-6"); */
/*    test_iso_parse("2016-W03-1T03:59:59.353"); */

/*    // Additional PERL only ISO formats */
/* #if 0 //ENABLE_PERL */
/*    test_iso_parse("20090621T0545Z"); */
/*    test_iso_parse("19850412T101530"); */
/*    test_iso_parse("850412T101530"); */
/*    test_iso_parse("2007-04-05T14:30"); */
/*    test_iso_parse("2007-04-05T14:30Z"); */
/*    test_iso_parse("2007-04-05T14:30B"); */
/*    test_iso_parse("    2007-04-05T14:30Z      "); */
/* #endif */

/*    return(0); */

/* } */
