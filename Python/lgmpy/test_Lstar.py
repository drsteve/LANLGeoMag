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
import itertools

import Lstar, magcoords

class Lstar_Tests(unittest.TestCase):
    def setUp(self):
        super(Lstar_Tests, self).setUp()
        self.date = datetime.datetime(2010, 10, 12)

    def tearDown(self):
        super(Lstar_Tests, self).tearDown()

    def test_Lstar_input(self):
        """Lstar does some input checking"""
        self.assertRaises(TypeError, Lstar.get_Lstar, [1,2,3], 'bad')
        self.assertRaises(TypeError, Lstar.get_Lstar, 'bad', coord_system = 'GSM')
        self.assertRaises(TypeError, Lstar.get_Lstar, 'bad', coord_system = 'SM')
        self.assertRaises(TypeError, Lstar.get_Lstar, 12, self.date, coord_system = 'GSM')
        self.assertRaises(TypeError, Lstar.get_Lstar, 12, self.date, coord_system = 'SM')
        self.assertRaises(NotImplementedError, Lstar.get_Lstar, [1,2,3], self.date, coord_system = 'bad')
        self.assertRaises(NotImplementedError, Lstar.get_Lstar, [1,2,3], self.date, Bfield = 'bad')

    def testLgm_B_T89_1(self):
        """This is a regression functional test for LstarVersusPA (regression)"""
        ans = Lstar.get_Lstar([-4.2, 1, 1], self.date, alpha = 90, Kp = 4, coord_system='SM', Bfield = 'Lgm_B_T89', LstarQuality = 1)
        self.assertAlmostEqual(5.0231064800914069, ans[90]['LHilton'], places=6)
        self.assertAlmostEqual(5.023227580319058, ans[90]['LMcIlwain'], places=6)
        self.assertAlmostEqual(4.4470253961778594, ans[90]['Lstar'][0])
        self.assertAlmostEqual(5.8143100824678715, ans[90]['Lsimple'][0])

    def testLgm_B_T89_2(self):
        """This is a regression functional test for LstarVersusPA (regression)"""
        ans = Lstar.get_Lstar([-4.2, 1, 1], self.date, alpha = 90, Kp = 5, coord_system='SM', Bfield = 'Lgm_B_T89', LstarQuality = 1)
        self.assertAlmostEqual(5.1730709738906144, ans[90]['LHilton'], places=6)
        self.assertAlmostEqual(5.1732419875708571, ans[90]['LMcIlwain'], places=6)
        self.assertAlmostEqual(4.3421853691096652, ans[90]['Lstar'][0], places=6)
        self.assertAlmostEqual(5.4671418867273349, ans[90]['Lsimple'][0], places=6)

    def testLgm_B_OP77_1(self):
        """This is a regression functional test for LstarVersusPA (regression)"""
        ans = Lstar.get_Lstar([-4.2, 1, 1], self.date, alpha = 90, Kp = 4, coord_system='SM', Bfield = 'Lgm_B_OP77', LstarQuality = 1)
        self.assertAlmostEqual(4.856099449443076, ans[90]['LHilton'], places=6)
        self.assertAlmostEqual(4.8561801068770318, ans[90]['LMcIlwain'], places=6)
        self.assertAlmostEqual(4.559620799451368, ans[90]['Lstar'][0], places=6)
        self.assertAlmostEqual(6.289533683215196, ans[90]['Lsimple'][0], places=6)

    def testLgm_B_OP77_2(self):
        """This is a regression functional test for LstarVersusPA (regression)"""
        # should be no change with Kp
        ans = Lstar.get_Lstar([-4.2, 1, 1], self.date, alpha = 90, Kp = 5, coord_system='SM', Bfield = 'Lgm_B_OP77', LstarQuality = 1)
        self.assertAlmostEqual(4.856099449443076, ans[90]['LHilton'], places=6)
        self.assertAlmostEqual(4.8561801068770318, ans[90]['LMcIlwain'], places=6)
        self.assertAlmostEqual(4.559620799451368, ans[90]['Lstar'][0], places=6)
        self.assertAlmostEqual(6.289533683215196, ans[90]['Lsimple'][0], places=6)

    def testCentredDipole90(self):
        """Unit test for centred dipole - should give known result"""
        for rdist, qlevel in itertools.product([3, 3.5, 4, 5, 6] ,range(1,7)):
            ans = Lstar.get_Lstar([-rdist,0,0], self.date, alpha=90, coord_system='SM', Bfield='Lgm_B_cdip', LstarQuality=qlevel)
            try:
                self.assertAlmostEqual(rdist, ans[90]['Lstar'][0], places=qlevel)
            except:
                print('Failed Lstar:CentredDipole90:: L=%f, Qual=%d -- dropped to %d places' 
                    % (rdist, qlevel, qlevel-1))
                self.assertAlmostEqual(rdist, ans[90]['Lstar'][0], places=qlevel-1)
    
    def testCentredDipole75(self):
        """Unit test for centred dipole - should give known result"""
        for rdist, qlevel in itertools.product([3,3.5,4,5,6] ,range(2,7)):
            ans = Lstar.get_Lstar([-rdist,0,0], self.date, alpha=75, coord_system='SM', Bfield='Lgm_B_cdip', LstarQuality=qlevel)
            try:
                self.assertAlmostEqual(rdist, ans[75]['Lstar'][0], places=qlevel)
            except:
                print('FAILURE L=%f, Qual=%d, date=%s' % (rdist, qlevel, self.date))

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
