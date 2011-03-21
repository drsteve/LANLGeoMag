#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Test suite for the LstarVersusPA

@author: Brian Larsen
@organization: LANL
@contact: balarsen@lanl.gov

@version: V1: 21-Mar-2011 (BAL)
"""

import unittest
import datetime

import numpy

import LstarVersusPA
from Lgm_Wrap import Lgm_LeapSeconds, size_Lgm_LeapSeconds, size_CTrans
from Lgm_Wrap import Lgm_DateTime, size_DateTime, Lgm_Set_Coord_Transforms
from Lgm_Wrap import Lgm_Convert_Coords, GSM_TO_SM, SM_TO_GSM

class LstarVersusPA_Tests(unittest.TestCase):
    def setUp(self):
        super(LstarVersusPA_Tests, self).setUp()
        self.date = datetime.datetime(2010, 10, 12)

    def tearDown(self):
        super(LstarVersusPA_Tests, self).tearDown()

    def testLgm_B_T89_1(self):
        """This is a regression functional test for LstarVersusPA (regression)"""
        ans = LstarVersusPA.LstarVersusPA([-4.2, 1, 1], self.date, alpha = 90, Kp = 4, coord_system='SM', Bfield = 'Lgm_B_T89', LstarQuality = 1)
        self.assertAlmostEqual(5.02624236376, ans[90]['LHilton'])
        self.assertAlmostEqual(5.02636496111, ans[90]['LMcIlwain'])
        self.assertAlmostEqual(4.44658097, ans[90]['Lstar'][0])
        self.assertAlmostEqual(5.80944453, ans[90]['Lsimple'][0])

    def testLgm_B_T89_2(self):
        """This is a regression functional test for LstarVersusPA (regression)"""
        ans = LstarVersusPA.LstarVersusPA([-4.2, 1, 1], self.date, alpha = 90, Kp = 5, coord_system='SM', Bfield = 'Lgm_B_T89', LstarQuality = 1)
        self.assertAlmostEqual(5.17724847477, ans[90]['LHilton'], places=6)
        self.assertAlmostEqual(5.17742158538, ans[90]['LMcIlwain'], places=6)
        self.assertAlmostEqual(4.35534221, ans[90]['Lstar'][0])
        self.assertAlmostEqual(5.46755074, ans[90]['Lsimple'][0])

    def testLgm_B_OP77_1(self):
        """This is a regression functional test for LstarVersusPA (regression)"""
        ans = LstarVersusPA.LstarVersusPA([-4.2, 1, 1], self.date, alpha = 90, Kp = 4, coord_system='SM', Bfield = 'Lgm_B_OP77', LstarQuality = 1)
        self.assertAlmostEqual(4.85874360104, ans[90]['LHilton'])
        self.assertAlmostEqual(4.8588254004, ans[90]['LMcIlwain'])
        self.assertAlmostEqual(4.55655149, ans[90]['Lstar'][0])
        self.assertAlmostEqual(6.29265777, ans[90]['Lsimple'][0])

    def testLgm_B_OP77_2(self):
        """This is a regression functional test for LstarVersusPA (regression)"""
        # should be no change with Kp
        ans = LstarVersusPA.LstarVersusPA([-4.2, 1, 1], self.date, alpha = 90, Kp = 5, coord_system='SM', Bfield = 'Lgm_B_OP77', LstarQuality = 1)
        self.assertAlmostEqual(4.85874360104, ans[90]['LHilton'])
        self.assertAlmostEqual(4.8588254004, ans[90]['LMcIlwain'])
        self.assertAlmostEqual(4.55655149, ans[90]['Lstar'][0])
        self.assertAlmostEqual(6.29265777, ans[90]['Lsimple'][0])

class Lstar_Data_Tests(unittest.TestCase):
    def setUp(self):
        super(Lstar_Data_Tests, self).setUp()

    def tearDown(self):
        super(Lstar_Data_Tests, self).tearDown()

    def test_init(self):
        """Lstar_Data can be created"""
        a = LstarVersusPA.Lstar_Data()
        self.assertTrue(hasattr(a, 'attrs'))
        self.assertTrue(a['position'] == {})



if __name__ == '__main__':
    unittest.main()
