#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Test suite for the IsoTimeStringToDateTime

@author: Brian Larsen
"""

import unittest
import ctypes
import numpy as np

from lgmpy.Lgm_Wrap import IsoTimeStringToDateTime, Lgm_DateTime
from lgmpy import Lgm_CTrans

all = ['IsoTimeStringToDateTime_Tests']

class IsoTimeStringToDateTime_Tests(unittest.TestCase):
    def test_Lgm_IsoTimeStringToDateTime(self):
        """regression on IsoTimeStringToDateTime"""
        lgmct = Lgm_CTrans.Lgm_CTrans()
        lgmdt = Lgm_DateTime()
        str1 = np.asarray('2001')
        IsoTimeStringToDateTime(str1.ctypes.data_as(ctypes.POINTER(ctypes.c_char)), 
                                ctypes.pointer(lgmdt), 
                                ctypes.pointer(lgmct))
        self.assertEqual(lgmdt.Date, 20010101)
        self.assertEqual(lgmdt.Day, 1)
        self.assertEqual(lgmdt.DaySeconds, 86400.0)
        self.assertEqual(lgmdt.Dow, 1)
        self.assertEqual(lgmdt.DowStr, 'Mon')
        self.assertEqual(lgmdt.Doy, 1)
        self.assertEqual(lgmdt.Hour, 0)
        self.assertEqual(lgmdt.Minute, 0)
        self.assertEqual(lgmdt.Second, 0)
        self.assertAlmostEqual(lgmdt.T, 0.010006844626967831)
        self.assertEqual(lgmdt.TZD_hh, 0)
        self.assertEqual(lgmdt.TZD_mm, 0)
        self.assertEqual(lgmdt.TZD_sgn, 1)
        self.assertEqual(lgmdt.Time, 0.0)
        self.assertEqual(lgmdt.TimeSystem, 0)
        self.assertEqual(lgmdt.Week, 1)
        self.assertEqual(lgmdt.Year, 2001)
        self.assertEqual(lgmdt.fYear, 2001.0)
        self.assertEqual(lgmdt.wYear, 2001)


if __name__ == '__main__':
    unittest.main()

