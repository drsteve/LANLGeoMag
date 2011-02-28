#!/usr/bin/python

import re


#/*
# *  Dates based on Months and Days
# */
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
#//ALLOWED? #define ISO_YYMM            18      // 2-digit years (is this allowed?) (assumes UTC == 0)

#define ISO_YYYYWww         19      // 4-digit years (assumes UTC == 0)
#define ISO_YYWww           20      // 2-digit years (assumes UTC == 0)

#define ISO_YYYY            21      // 4-digit year 
#define ISO_YY              22      // 2-digit year 

# 
#  Notes.
# 
#   In ISO 8601, you can delete the least significant fields all the way to the
#   2-digit century. For example,
# 
#       YYYY-MM-DDTHH:MM:SS.SSSSSS+HH:MM
#       YYYY-MM-DDTHH:MM:SS.SSSSSS+HH
#       YYYY-MM-DDTHH:MM:SS.SSSSSS
#       YYYY-MM-DDTHH:MM:SS
#       YYYY-MM-DDTHH:MM
#       YYYY-MM-DDTHH
#       YYYY-MM-DD
#       YYYY-MM
#       YYYY
#       YY
# 
#   The 2 Y's in the last case are not the usual 2-digit years, they are the
#   century (e.g. if YYYY is 1996, then the last case above would be YY=19).
# 
#   The separators can be omitted for all of these except for the YYYY-MM case.
#   That one always has to have a separator or it will get confused with other
#   formats.
# 
#   The original ISO 8601 allowed "truncated representations" meaning that
#   2-digit years could be used. But apparently the latest version od ISO 8601
#   got rid of that.
#   
# 
# 
# 
# 
#


#      
#      ISO 8601 Time Formats are strictly formated, but have a number of variations
#      as follows;
#     
#          YYYY
#          YYYY-MM
#          YYYY-MM-DD (hypens optional)
#          YYYY-MM-DDThh:mmTZD (hypens and colons optional)
#          YYYY-MM-DDThh:mm:ssTZD (hypens and colons optional)
#          YYYY-Www (hypens can be omitted)
#          YYYY-Www-D (hypens can be omitted)
#     
#      where,
#      
#          YYYY = four-digit year
#          MM   = two-digit month (01=January, etc.)
#          DD   = two-digit day of month (01 through 31)
#          hh   = two digits of hour (00 through 23) (am/pm NOT allowed)
#          mm   = two digits of minute (00 through 59)
#          ss   = two digits of second (00 through 59)
#          s    = one or more digits representing a decimal fraction of a second
#          TZD  = time zone designator (Z or +hh:mm or -hh:mm) (colons optional)
#          ww   = ww is the week of the year
#          D    = ww is day of week 
#      
#      The hypens and colons are optional above except in the YYYY-MM format.
#     
#      YYYY can be positive or negative but only if users agree on it (so by
#      default its not in the standard). Here we assume it to be positive.
#     
#      Also, its not clear that if you omit one hyphen you have to omit them
#      all.  Or that omitting one colon means you have to omit them all. Check
#      the standard on this, but here I allow you to omit any number of them.
#      E.g., 2009-06-21T05:12:34.01Z and 2009-0621T0512:34.01Z should be
#      treated the same here.
#     
#      Combined Date and Time formats must have the T separating them, but you
#      can also have separate date and time formats and we recognize both E.g.,
#      2009-06-21T05:12:34.01Z versus 2009-06-21 05:12:34.01Z
#

AllowTruncatedReps = 1;

str = '2010-11-08T123456.345'

ISO_YYYYMMDDTHHMMSS = re.compile( r"^(\d{4})-?(\d{2})-?(\d{2})[T ](\d{2}):?(\d{2}):?(\d{2}[\.]?\d*)([+-].*|Z?)$" )
ISO_YYYYMMDDTHHMM   = re.compile( r"^(\d{4})-?(\d{2})-?(\d{2})[T ](\d{2}):?(\d{2})([+-].*|Z?)$" )
ISO_YYYYMMDDTHH     = re.compile( r"^(\d{4})-?(\d{2})-?(\d{2})[T ](\d{2})([+-].*|Z?)$" )





if (ISO_YYYYMMDDTHHMMSS.findall( str )):
    parts = ISO_YYYYMMDDTHHMMSS.findall( str )
elif (parts = ISO_YYYYMMDDTHHMM.findall( str )):
    print parts

    













