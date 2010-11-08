#include <EXTERN.h>
#include <perl.h>
#include <Lgm_CTrans.h>

/*
 *  Dates based on Months and Days
 */
#define ISO_ERROR           0       

#define ISO_YYYYMMDDTHHMMSS 1       // 4-digit years
#define ISO_YYMMDDTHHMMSS   2       // 2-digit years
#define ISO_YYYYDDDTHHMMSS  3       // 4-digit years
#define ISO_YYDDDTHHMMSS    4       // 2-digit years

#define ISO_YYYYWwwDTHHMMSS 5       // 4-digit years
#define ISO_YYWwwDTHHMMSS   6       // 2-digit years

#define ISO_YYYYMMDDTHHMM   7       // 4-digit years
#define ISO_YYMMDDTHHMM     8       // 2-digit years

#define ISO_YYYYWwwDTHHMM   9       // 4-digit years
#define ISO_YYWwwDTHHMM     10      // 2-digit years


#define ISO_YYYYMMDD        11      // 4-digit years (assumes UTC == 0)
#define ISO_YYMMDD          12      // 2-digit years (assumes UTC == 0)
#define ISO_YYYYDDD         13      // 4-digit years (assumes UTC == 0)
#define ISO_YYDDD           14      // 2-digit years (assumes UTC == 0)

#define ISO_YYYYWwwD        15      // 4-digit years (assumes UTC == 0)
#define ISO_YYWwwD          16      // 2-digit years (assumes UTC == 0)

#define ISO_YYYYMM          17      // 4-digit years (YYYY-MM only YYYYMM not allowed) (assumes UTC == 0)
//ALLOWED? #define ISO_YYMM            18      // 2-digit years (is this allowed?) (assumes UTC == 0)

#define ISO_YYYYWww         19      // 4-digit years (assumes UTC == 0)
#define ISO_YYWww           20      // 2-digit years (assumes UTC == 0)

#define ISO_YYYY            21      // 4-digit year 
#define ISO_YY              22      // 2-digit year 

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
 *
 *
 *
 *
 */


static PerlInterpreter *my_perl;

int ParseTimeString( char *TimeString ) {

    long int    Date;
    int         Year, wYear, Month, Day, Hours, Minutes, TZD_sgn, TZD_hh, TZD_mm, Week, DayOfWeek, DayOfYear;
    int         ISOFormat, TZDError;
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
                  " print \"str = [$str]\\n\";\n"
                  " my $ISOFormat = 0; my $Year = -9999; my $DayOfYear = -9999; my $Month = -9999; my $Day = -9999; my $Hour = -9999; my $Minute = -9999; my $Second = -9999; my $TZD = 'Z';\n", TimeString );
            

                /*
                 *  ISO_YYYYMMDDTHHMMSS  
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
                 *  ISO_YYMMDDTHHMMSS
                 *
                 *  Match ISO strings of the form:
                 *
                 *           YY[-]MM[-]DDThh[:]mm[:]ss[.ssssss](TZD or Z or nothing)
                 *
                 *  Examples: 
                 *        10-10-31T12:34:56.789-06:00
                 *        96-03-21T01:54:12.4534789+05:00
                 */
                if ( AllowTruncatedReps ) {
                    strcat( Str,    
                        " } elsif ( $str =~ m/^(\\d{2})-?(\\d{2})-?(\\d{2})[T ](\\d{2}):?(\\d{2}):?(\\d{2}[\\.]?\\d*)([+-].*|Z?)/ ) { \n\
                                $Year = ($1 > 50) ? $1+1900:$1+2000; $Month = $2; $Day = $3; $Hour = $4; $Minute = $5; $Second = $6; $TZD = $7; $ISOFormat = 2;\n"
                    );
                }

                /*
                 *  ISO_YYYYDDDTHHMMSS
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
                 *  ISO_YYDDDTHHMMSS
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
                if ( AllowTruncatedReps ) {
                    strcat( Str, 
                        " } elsif ( $str =~ m/^(\\d{2})-?(\\d{3})[T ](\\d{2}):?(\\d{2}):?(\\d{2}[\\.]?\\d*)([+-].*|Z?)/ ) { \n\
                                $Year = ($1 > 50) ? $1+1900:$1+2000; $DayOfYear = $2; $Hour = $3; $Minute = $4; $Second = $5; $TZD = $6; $ISOFormat = 4;\n"
                    );
                }









                  /*
                   *  ISO_YYYYWwwDTHHMMSS
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
                   *  ISO_YYWwwDTHHMMSS
                   *
                   *  Match ISO strings of the form:
                   *
                   *    YY[-]Www[-]DThh[:]mm[:]ss[.ssssss](TZD or Z or nothing)
                   *
                   *  Examples:  
                   *       86-W13-2T09:45:32
                   *       86W132T094532
                   */
                if ( AllowTruncatedReps ) {
                    strcat( Str, 
                    " } elsif ( $str =~ m/^(\\d{2})-?W(\\d{2})-?(\\d{1})[T ](\\d{2}):?(\\d{2}):?(\\d{2}[\\.]?\\d*)([+-].*|Z?)/ ) { \n\
                            $Year = ($1 > 50) ? $1+1900:$1+2000; $Week = $2; $DayOfWeek = $3; $Hour = $4; $Minute = $5; $Second = $6; $TZD = $7; $ISOFormat = 6;\n"
                    );
                }











                  /*
                   *  ISO_YYYYMMDDTHHMM
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
                 *  ISO_YYMMDDTHHMM
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
                if ( AllowTruncatedReps ) {
                    strcat( Str, 
                      " } elsif ( $str =~ m/^(\\d{2})-?(\\d{2})-?(\\d{2})[T ](\\d{2}):?(\\d{2})([+-].*|Z?)$/ ) { \n\
                            $Year = ($1 > 50) ? $1+1900:$1+2000; $Month = $2; $Day = $3; $Hour = $4; $Minute = $5; $Second = 0.0; $TZD = $6; $ISOFormat = 8;\n"
                    );
                }











                  /*
                   *  ISO_YYYYWwwDTHHMM
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
                 *  ISO_YYWwwDTHHMM
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
                if ( AllowTruncatedReps ) {
                    strcat( Str, 
                            " } elsif ( $str =~ m/^(\\d{2})-?W(\\d{2})-?(\\d{1})[T ](\\d{2}):?(\\d{2})([+-].*|Z?)$/ ) { \n\
                                    $Year = ($1 > 50) ? $1+1900:$1+2000; $Week = $2; $DayOfWeek = $3; $Hour = $4; $Minute = $5; $Second = 0.0; $TZD = $6; $ISOFormat = 10;\n"
                    );
                }










                  /*
                   *  ISO_YYYYMMDD
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
                   *  ISO_YYMMDD
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
                if ( AllowTruncatedReps ) {
                strcat( Str, 
                  " } elsif ( $str =~ m/^(\\d{2})-?(\\d{2})-?(\\d{2})$/ ) {\n\
                        $Year = ($1 > 50) ? $1+1900:$1+2000; $Month = $2; $Day = $3; $Hour = $Minute = 0; $Second = 0.0; $ISOFormat = 12;\n"
                );
                }


                  /*
                   *  ISO_YYYYDDD
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
                   *  ISO_YYDDD
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
                if ( AllowTruncatedReps ) {
                strcat( Str, 
                  " } elsif ( $str =~ m/^(\\d{2})-?(\\d{3})$/ ) {\n\
                        $Year = ($1 > 50) ? $1+1900:$1+2000; $DayOfYear = $2; $Hour = $Minute = 0; $Second = 0.0;$ISOFormat = 14;\n"
                );
                }













                  /*
                   *  ISO_YYYYWwwD
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
                   *  ISO_YYWwwD
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
                   *  ISO_YYYYMM
                   *
                   *  Match ISO strings of the form:
                   *
                   *        YYYY[-]MM
                   *
                   *  Examples: 
                   *        1976-12
                   *        197612
                   */
                strcat( Str, 
                  " } elsif ( $str =~ m/^(\\d{4})-(\\d{2})$/ ) { \n\
                        $Year = $1; $Month = $2;  $Day = 1; $Hour = $Minute = 0; $Second = 0.0;$ISOFormat = 17;\n"
                );

                  /*
                   *  ISO_YYMM
                   *
                   *  Match ISO strings of the form:
                   *
                   *    YY-MM
                   *  
                   *  The hyphen is mandatory here.
                   *  Examples: 
                   *        76-12
                   */
                if ( AllowTruncatedReps ) {
                strcat( Str, 
                  " } elsif ( $str =~ m/^(\\d{2})-(\\d{2})$/ ) { \n\
                        $Year = ($1 > 50) ? $1+1900:$1+2000; $Day = 1; $Month = $2;  $Hour = $Minute = 0; $Second = 0.0;$ISOFormat = 18;\n"
                );
                }









                  /*
                   *  ISO_YYYYWww
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
                   *  ISO_YYWww
                   *  
                   *  Match ISO strings of the form:
                   *  
                   *        YY-Www
                   *
                   *  Examples: 
                   *        76-W32
                   */
                if ( AllowTruncatedReps ) {
                strcat( Str, 
                  " } elsif ( $str =~ m/^(\\d{2})-W(\\d{2})$/ ) { \n\
                        $Year = ($1 > 50) ? $1+1900:$1+2000; $Week = $2; $DayOfWeek = 1; $Hour = $Minute = 0; $Second = 0.0; $ISOFormat = 20;\n"
                );
                }







                /*
                 *  ISO_YYYY
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
                 *  ISO_YY
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


    if ( (ISOFormat == ISO_YYYYWwwDTHHMMSS) || (ISOFormat == ISO_YYWwwDTHHMMSS)
            || (ISOFormat == ISO_YYYYWwwDTHHMM) || (ISOFormat == ISO_YYWwwDTHHMM)
            || (ISOFormat == ISO_YYYYWwwD) || (ISOFormat == ISO_YYWwwD)
            || (ISOFormat == ISO_YYYYWww) || (ISOFormat == ISO_YYWww) ) {

        printf( "Warning: need to add code to conveert weeks\n");

    } else {

        if ( DayOfYear > 0 ){
            Date = Year*1000 + DayOfYear;
        } else {
            Date = Year*10000 + Month*100 + Day;
        }
        if ( Lgm_IsValidDate( Date ) ) {
            Lgm_Doy( Date, &Year, &Month, &Day, &DayOfYear );
        } else {
            printf( "Date is invalid\n" );
            return(-1);
        }

    }


    printf( "Year      = %d\n",  Year );
    printf( "DayOfYear = %d\n",  DayOfYear );
    printf( "Month     = %d\n",  Month );
    printf( "Day       = %d\n",  Day );
    printf( "Hours     = %d\n",  Hours );
    printf( "Minutes   = %d\n",  Minutes );
    printf( "Seconds   = %lf\n", Seconds );
    printf( "TZD_sgn   = %d\n",  TZD_sgn );
    printf( "TZD_hh    = %d\n",  TZD_hh );
    printf( "TZD_mm    = %d\n",  TZD_mm );
    printf( "Week      = %d\n",  Week );
    printf( "DayOfWeek = %d\n",  DayOfWeek );
    printf( "TZDError  = %d\n",  TZDError );
    printf( "ISOFormat = %d\n",  ISOFormat );



    /*
     * Test for sane numbers
     */
//    if ( (ISOFormat == ISO_YYYYMMDDTHHMMSS) 

        if ( ( Year < 0 ) || ( Year > 9999 ) ) {
            printf( "Year Out of Range! (Year = %d)\n", Year );
            return( -1 );
        }
        if ( ( DayOfYear < 1 ) || ( DayOfYear > 365+Lgm_LeapYear(Year) ) ) {
            printf( "DayOfYear Out of Range! (DayOfYear = %d)\n", DayOfYear );
            return( -1 );
        }
        if ( ( Month < 1 ) || ( Month > 12 ) ) {
            printf( "Month Out of Range! (Month = %d)\n", Month );
            return( -1 );
        }
        if ( ( Day < 1 ) || ( Day > 31 ) ) {
            printf( "Day Out of Range! (Day = %d)\n", Day );
            return( -1 );
        }
        if ( ( Hours < 0 ) || ( Hours > 23 ) ) {
            printf( "Hours Out of Range! (Hours = %d)\n", Hours );
            return( -1 );
        }
        if ( ( Minutes < 0 ) || ( Minutes > 59 ) ) {
            printf( "Minutes Out of Range! (Minutes = %d)\n", Minutes );
            return( -1 );
        }
        if ( ( Seconds < 0 ) || ( Seconds > 60 ) ) { // Yes seconds can be 60 to accomodate leap seconds
            printf( "Seconds Out of Range! (Seconds = %lf)\n", Seconds ); 
            return( -1 );
        }
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

/*
    } else if ( (ISOFormat == ISO_YYYYDDDTHHMMSS) || (ISOFormat == ISO_YYDDDTHHMMSS) 
        || (ISOFormat == ISO_YYYYDDD) || (ISOFormat == ISO_YYDDD) ) {

        if ( ( Year < 0 ) || ( Year > 9999 ) ) {
            printf( "Year Out of Range! (Year = %d)\n", Year );
            return( -1 );
        }
        if ( ( DayOfYear < 1 ) || ( DayOfYear > 365+Lgm_LeapYear(Year) ) ) {
            printf( "DayOfYear Out of Range! (DayOfYear = %d)\n", DayOfYear );
            return( -1 );
        }
        if ( ( Hours < 0 ) || ( Hours > 23 ) ) {
            printf( "Hours Out of Range! (Hours = %d)\n", Hours );
            return( -1 );
        }
        if ( ( Minutes < 0 ) || ( Minutes > 59 ) ) {
            printf( "Minutes Out of Range! (Minutes = %d)\n", Minutes );
            return( -1 );
        }
        if ( ( Seconds < 0 ) || ( Seconds > 60 ) ) { // Yes seconds can be 60 to accomodate leap seconds
            printf( "Seconds Out of Range! (Seconds = %lf)\n", Seconds ); 
            return( -1 );
        }


    } else if ( (ISOFormat = ISO_YYYYWwwDTHHMMSS) || (ISOFormat = ISO_YYYYWwwDTHHMM) 
            || (ISOFormat = ISO_YYWwwDTHHMMSS) || (ISOFormat = ISO_YYWwwDTHHMM) ) {

        if (ISOFormat == ISO_YYYYWwwDTHHMMSS) Seconds = 0.0; // assume seconds are zero
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
        
        if ( ISOFormat == ISO_YYYYMM ) Day = 1; //assume day is 1st day of Month

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

        if ( ISOFormat == ISO_YYYYWww ) DayOfWeek = 1; //assume day is 1st day of Week
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
*/

    if (TZDError ) {
        printf( "Unrecognized time zone field!\n" );
        return(-1);
    }


    return( 1 );


}



int main (int argc, char **argv, char **env) {

    int     Flag;
    char    TimeString[128];

    sprintf(TimeString, "2009-06-21T05:12:34.123456789Z");
    sprintf(TimeString, "2009-06-21                           0512:34.01Z");
    sprintf(TimeString, "2009-W52-7T2312:34.01Z");
    sprintf(TimeString, "2009-W52-1T2312:34.01Z");
    sprintf(TimeString, "2009-06-21T05:1234.01Z");
    sprintf(TimeString, "2009-0621T19:12:34.123456789-03:34");

//    sprintf(TimeString, "2009-06-21T05:12:34.123456789Z"); printf("TimeString = %s\n", TimeString ); Flag = ParseTimeString( TimeString ); printf( "Flag = %d\n", Flag ); 
//    sprintf(TimeString, "2009-06-21T05:12:34.123456789"); printf("TimeString = %s\n", TimeString ); Flag = ParseTimeString( TimeString ); printf( "Flag = %d\n", Flag ); 
//    sprintf(TimeString, "20090621T0545Z"); printf("TimeString = %s\n", TimeString ); Flag = ParseTimeString( TimeString ); printf( "Flag = %d\n", Flag ); 
//    sprintf(TimeString, "2010-11-03T16:15:02-06:00"); printf("TimeString = %s\n", TimeString ); Flag = ParseTimeString( TimeString ); printf( "Flag = %d\n\n", Flag ); 
//    sprintf(TimeString, "2010-11-03T16:15:02-06"); printf("TimeString = %s\n", TimeString ); Flag = ParseTimeString( TimeString ); printf( "Flag = %d\n\n", Flag ); 
//    sprintf(TimeString, "2010-11-03T16:15:02-6"); printf("TimeString = %s\n", TimeString ); Flag = ParseTimeString( TimeString ); printf( "Flag = %d\n\n", Flag ); 
//    sprintf(TimeString, "2010-11-03T16:15:02"); printf("TimeString = %s\n", TimeString ); Flag = ParseTimeString( TimeString ); printf( "Flag = %d\n\n", Flag ); 
//    sprintf(TimeString, "19850412T101530"); printf("TimeString = %s\n", TimeString ); Flag = ParseTimeString( TimeString ); printf( "Flag = %d\n\n", Flag ); 
//    sprintf(TimeString, "850412T101530"); printf("TimeString = %s\n", TimeString ); Flag = ParseTimeString( TimeString ); printf( "Flag = %d\n\n", Flag ); 
    sprintf(TimeString, argv[1]); printf("TimeString = %s\n", TimeString ); Flag = ParseTimeString( TimeString ); printf( "Flag = %d\n\n", Flag ); 

//    sprintf(TimeString, "2007-04-05T14:30"); printf("TimeString = %s\n", TimeString ); Flag = ParseTimeString( TimeString ); printf( "Flag = %d\n", Flag ); 
//    sprintf(TimeString, "    2007-04-05T14:30Z      "); printf("TimeString = %s\n", TimeString ); Flag = ParseTimeString( TimeString ); printf( "Flag = %d\n", Flag ); 
//    sprintf(TimeString, "2007-04-05T14:30Z"); printf("TimeString = %s\n", TimeString ); Flag = ParseTimeString( TimeString ); printf( "Flag = %d\n", Flag ); 
//    sprintf(TimeString, "2007-04-05T14:30B"); printf("TimeString = %s\n", TimeString ); Flag = ParseTimeString( TimeString ); printf( "Flag = %d\n", Flag ); 

    return(0);

}
