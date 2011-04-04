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

import Lstar, magcoords

class Lstar_Tests(unittest.TestCase):
    def setUp(self):
        super(Lstar_Tests, self).setUp()
        self.date = datetime.datetime(2010, 10, 12)

    def tearDown(self):
        super(Lstar_Tests, self).tearDown()

    def testLgm_B_T89_1(self):
        """This is a regression functional test for LstarVersusPA (regression)"""
        ans = Lstar.get_Lstar([-4.2, 1, 1], self.date, alpha = 90, Kp = 4, coord_system='SM', Bfield = 'Lgm_B_T89', LstarQuality = 1)
        self.assertAlmostEqual(5.02320512, ans[90]['LHilton'])
        self.assertAlmostEqual(5.02332627, ans[90]['LMcIlwain'])
        self.assertAlmostEqual(4.44428332, ans[90]['Lstar'][0])
        self.assertAlmostEqual(5.80944453, ans[90]['Lsimple'][0])

    def testLgm_B_T89_2(self):
        """This is a regression functional test for LstarVersusPA (regression)"""
        ans = Lstar.get_Lstar([-4.2, 1, 1], self.date, alpha = 90, Kp = 5, coord_system='SM', Bfield = 'Lgm_B_T89', LstarQuality = 1)
        self.assertAlmostEqual(5.17317697, ans[90]['LHilton'], places=6)
        self.assertAlmostEqual(5.17334804, ans[90]['LMcIlwain'], places=6)
        self.assertAlmostEqual(4.35283590, ans[90]['Lstar'][0], places=6)
        self.assertAlmostEqual(5.46755074, ans[90]['Lsimple'][0], places=6)

    def testLgm_B_OP77_1(self):
        """This is a regression functional test for LstarVersusPA (regression)"""
        ans = Lstar.get_Lstar([-4.2, 1, 1], self.date, alpha = 90, Kp = 4, coord_system='SM', Bfield = 'Lgm_B_OP77', LstarQuality = 1)
        self.assertAlmostEqual(4.85646725, ans[90]['LHilton'])
        self.assertAlmostEqual(4.85654807, ans[90]['LMcIlwain'])
        self.assertAlmostEqual(4.55863046, ans[90]['Lstar'][0])
        self.assertAlmostEqual(6.29265777, ans[90]['Lsimple'][0])

    def testLgm_B_OP77_2(self):
        """This is a regression functional test for LstarVersusPA (regression)"""
        # should be no change with Kp
        ans = Lstar.get_Lstar([-4.2, 1, 1], self.date, alpha = 90, Kp = 5, coord_system='SM', Bfield = 'Lgm_B_OP77', LstarQuality = 1)
        self.assertAlmostEqual(4.85646725, ans[90]['LHilton'])
        self.assertAlmostEqual(4.85654807, ans[90]['LMcIlwain'])
        self.assertAlmostEqual(4.55863046, ans[90]['Lstar'][0])
        self.assertAlmostEqual(6.29265777, ans[90]['Lsimple'][0])

class Lstar_Data_Tests(unittest.TestCase):
    def setUp(self):
        super(Lstar_Data_Tests, self).setUp()

    def tearDown(self):
        super(Lstar_Data_Tests, self).tearDown()

    def test_init(self):
        """Lstar_Data can be created"""
        a = Lstar.Lstar_Data()
        self.assertTrue(hasattr(a, 'attrs'))
        self.assertTrue(a['position'] == {})

class Lvalue_Tests(unittest.TestCase):
    def setUp(self):
        super(Lvalue_Tests, self).setUp()
        self.date = datetime.datetime(2010, 10, 12)

    def tearDown(self):
        super(Lvalue_Tests, self).tearDown()

    def testL_B_T89_(self):
        """Test that Lvalue code returns same value of L returned by Lstar"""
        ans = Lstar.get_Lstar([-4.2, 1, 1], self.date, alpha = 90, Kp = 4, coord_system='GSM', Bfield = 'Lgm_B_T89', LstarQuality = 3)
        ans2 = magcoords.Lvalue([-4.2, 1, 1], self.date, alpha=90, Kp=4, Bfield = 'Lgm_B_T89', method='Hilton')
        ans3 = magcoords.Lvalue([-4.2, 1, 1], self.date, alpha=90, Kp=4, Bfield = 'Lgm_B_T89', method='McIlwain')
        self.assertAlmostEqual(ans2['L'], ans[90]['LHilton'], places=5)
        self.assertAlmostEqual(ans3['L'], ans[90]['LMcIlwain'], places=5)

    def test_L_cdip(self):
        LH = magcoords.Lvalue([-6, 0, 0], self.date, alpha=90, Bfield='Lgm_B_cdip', method='Hilton', coord_system='SM')
        LM = magcoords.Lvalue([-6, 0, 0], self.date, alpha=90, Bfield='Lgm_B_cdip', method='McIlwain', coord_system='SM')
        self.assertAlmostEqual(LH['L'], 6.00000, places=6)
        self.assertAlmostEqual(LM['L'], 6.00000, places=6)

if __name__ == '__main__':
    unittest.main()
