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

    def test_Lstar_input(self):
        """Lstar does some inpout checking"""
        self.assertRaises(TypeError, Lstar.get_Lstar, [1,2,3], 'bad')
        self.assertRaises(TypeError, Lstar.get_Lstar, 'bad', coord_system = 'GSM')
        self.assertRaises(TypeError, Lstar.get_Lstar, 'bad', coord_system = 'SM')
        self.assertRaises(NotImplementedError, Lstar.get_Lstar, [1,2,3], self.date, coord_system = 'bad')

    def testLgm_B_T89_1(self):
        """This is a regression functional test for LstarVersusPA (regression)"""
        ans = Lstar.get_Lstar([-4.2, 1, 1], self.date, alpha = 90, Kp = 4, coord_system='SM', Bfield = 'Lgm_B_T89', LstarQuality = 1)
        self.assertAlmostEqual(5.02319883, ans[90]['LHilton'])
        self.assertAlmostEqual(5.02331998, ans[90]['LMcIlwain'])
        self.assertAlmostEqual(4.44427733, ans[90]['Lstar'][0])
        self.assertAlmostEqual(5.809442661, ans[90]['Lsimple'][0])

    def testLgm_B_T89_2(self):
        """This is a regression functional test for LstarVersusPA (regression)"""
        ans = Lstar.get_Lstar([-4.2, 1, 1], self.date, alpha = 90, Kp = 5, coord_system='SM', Bfield = 'Lgm_B_T89', LstarQuality = 1)
        self.assertAlmostEqual(5.17315141, ans[90]['LHilton'], places=6)
        self.assertAlmostEqual(5.17332246, ans[90]['LMcIlwain'], places=6)
        self.assertAlmostEqual(4.35281882, ans[90]['Lstar'][0], places=6)
        self.assertAlmostEqual(5.467548617, ans[90]['Lsimple'][0], places=6)

    def testLgm_B_OP77_1(self):
        """This is a regression functional test for LstarVersusPA (regression)"""
        ans = Lstar.get_Lstar([-4.2, 1, 1], self.date, alpha = 90, Kp = 4, coord_system='SM', Bfield = 'Lgm_B_OP77', LstarQuality = 1)
        self.assertAlmostEqual(4.856172936788191, ans[90]['LHilton'], places=6)
        self.assertAlmostEqual(4.856253625894013, ans[90]['LMcIlwain'], places=6)
        self.assertAlmostEqual(4.558376476155467, ans[90]['Lstar'][0], places=6)
        self.assertAlmostEqual(6.292656284529973, ans[90]['Lsimple'][0], places=6)

    def testLgm_B_OP77_2(self):
        """This is a regression functional test for LstarVersusPA (regression)"""
        # should be no change with Kp
        ans = Lstar.get_Lstar([-4.2, 1, 1], self.date, alpha = 90, Kp = 5, coord_system='SM', Bfield = 'Lgm_B_OP77', LstarQuality = 1)
        self.assertAlmostEqual(4.856172936788191, ans[90]['LHilton'], places=6)
        self.assertAlmostEqual(4.856253625894013, ans[90]['LMcIlwain'], places=6)
        self.assertAlmostEqual(4.558376476155467, ans[90]['Lstar'][0], places=6)
        self.assertAlmostEqual(6.292656284529973, ans[90]['Lsimple'][0], places=6)

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
        self.assertAlmostEqual(ans2['L'], ans[90]['LHilton'], places=4)
        self.assertAlmostEqual(ans3['L'], ans[90]['LMcIlwain'], places=4)

    def test_L_cdip(self):
        LH = magcoords.Lvalue([-6, 0, 0], self.date, alpha=90, Bfield='Lgm_B_cdip', method='Hilton', coord_system='SM')
        LM = magcoords.Lvalue([-6, 0, 0], self.date, alpha=90, Bfield='Lgm_B_cdip', method='McIlwain', coord_system='SM')
        self.assertAlmostEqual(LH['L'], 6.00000, places=6)
        self.assertAlmostEqual(LM['L'], 6.00000, places=6)

if __name__ == '__main__':
    unittest.main()
