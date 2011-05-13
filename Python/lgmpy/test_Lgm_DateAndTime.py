#!/usr/bin/env python

"""
Test suite for the Lgm_DateAndTime file
These are automaticlly wrapped, tests added as regesrrion on the C

`Note` very partial at this time


@author: Brian Larsen
@organization: LANL
@contact: balarsen@lanl.gov

@version: V1: 13-May-2010 (BAL)
"""

from __future__ import division

import ctypes
import unittest
import datetime

import numpy

import Lgm_CTrans
import Lgm_Vector
from Lgm_Wrap import Lgm_LeapYear, Lgm_GetLeapSeconds, Lgm_IsLeapSecondDay,  Lgm_LoadLeapSeconds, Lgm_GPS_to_GpsSeconds, Lgm_DateTime, Lgm_GpsSeconds_to_GPS

class Lgm_DateAndTime_tests(unittest.TestCase):
    def setUp(self):
        super(Lgm_DateAndTime_tests, self).setUp()
        self.c = Lgm_CTrans.Lgm_CTrans()

    def tearDown(self):
        super(Lgm_DateAndTime_tests, self).tearDown()
        del self.c

    def test_Lgm_LeapYear(self):
        """leap year calculator should work"""
        years = range(1999, 2011)
        ans = [0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0]
        for i, val in enumerate(years):
            self.assertEqual(ans[i], Lgm_LeapYear(val))

    def test_Lgm_GetLeapSeconds(self):
        """Lgm_GetLeapSeconds should have known output"""
        # 2009/6/1 is 34s
        self.assertEqual(34, Lgm_GetLeapSeconds(2454984.0, ctypes.pointer(self.c)))

    def test_Lgm_IsLeapSecondDay(self):
        """Lgm_IsLeapSecondDay should have known output"""
        sec_in_day = ctypes.c_double()
        # 1 means it it is a leap second day
        self.assertEqual(1, Lgm_IsLeapSecondDay( 19920630, ctypes.pointer(sec_in_day), ctypes.pointer(self.c) ) )
        self.assertEqual(86401, sec_in_day.value)
        self.assertEqual(0, Lgm_IsLeapSecondDay( 19920629, ctypes.pointer(sec_in_day), ctypes.pointer(self.c) ) )
        self.assertEqual(86400, sec_in_day.value)

    def test_Lgm_LoadLeapSeconds(self):
        """Lgm_LoadLeapSeconds should work and not change"""
        self.assertTrue(Lgm_LoadLeapSeconds(ctypes.pointer(self.c)))
        # The LeapSecondDates values are the dates that the given corrections came into effect. But the actual leap second
        # date is the day before -- thats when the extra second was tacked on. So subtract a day from each of the dates.
        LeapSecondDates = [19720101, 19720701, 19730101, 19740101, 19750101,
                           19760101, 19770101, 19780101, 19790101, 19800101,
                           19810701, 19820701, 19830701, 19850701, 19880101,
                           19900101, 19910101, 19920701, 19930701, 19940701,
                           19960101, 19970701, 19990101, 20060101, 20090101,]
        LeapSecondJDs = [2441317.5, 2441499.5, 2441683.5, 2442048.5, 2442413.5,
                         2442778.5, 2443144.5, 2443509.5, 2443874.5, 2444239.5,
                         2444786.5, 2445151.5, 2445516.5, 2446247.5, 2447161.5,
                         2447892.5, 2448257.5, 2448804.5, 2449169.5, 2449534.5,
                         2450083.5, 2450630.5, 2451179.5, 2453736.5, 2454832.5,]
        LeapSeconds = [10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0,
                       19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0,
                       28.0, 29.0, 30.0, 31.0, 32.0, 33.0, 34.0, ]
        LeapSecondDates = Lgm_CTrans.dateLongToDate(LeapSecondDates)
        LeapSecondDates = [val - datetime.timedelta(days=1) for val in LeapSecondDates]
        LeapSecondDates = Lgm_CTrans.dateToDateLong(LeapSecondDates)

        for i, (v1, v2, v3) in enumerate(zip(LeapSecondDates, LeapSecondJDs, LeapSeconds)):
            self.assertEqual(v1, self.c.l.LeapSecondDates[i])
            self.assertEqual(v2, self.c.l.LeapSecondJDs[i])
            self.assertEqual(v3, self.c.l.LeapSeconds[i])

    def test_Lgm_GPS_to_GpsSeconds(self):
        """Lgm_GPS_to_GpsSeconds should give known values"""
        date = Lgm_DateTime(20001223, 2000, 12, 23, 0, 13+34/60+40/60/60, 13, 34, 40)
        self.assertEqual(661613680.0, Lgm_GPS_to_GpsSeconds(ctypes.pointer(date)))

    def test_Lgm_GpsSeconds_to_GPS(self):
        """Lgm_GpsSeconds_to_GPS should give known answer"""
        date = Lgm_DateTime()
        Lgm_GpsSeconds_to_GPS(661613680.0, ctypes.pointer(date))
        self.assertEqual(20001223, date.Date)
        self.assertAlmostEqual(13.577777777777778, date.Time)


if __name__ == '__main__':
    unittest.main()
