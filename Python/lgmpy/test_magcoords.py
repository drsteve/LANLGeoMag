#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Test suite for the LstarVersusPA

@author: Brian Larsen
@organization: LANL
@contact: balarsen@lanl.gov

@version: V1: 4-Apr-2011 (BAL)
"""

import unittest
import datetime

import numpy

import magcoords
import Lgm_Vector

class magcoords_Tests(unittest.TestCase):
    def setUp(self):
        super(magcoords_Tests, self).setUp()

    def tearDown(self):
        super(magcoords_Tests, self).tearDown()

    def test_end2end(self):
        """quick functional tests (regression)"""
        ans = [-3.608026916281572, 4.163336342344337e-17, -1.7268878861662331]
        numpy.testing.assert_array_almost_equal(ans,
                            magcoords.coordTrans([-4,0,0], datetime.datetime(2009,1,1),'SM','GSM'), decimal=6)
        ans = [-3.9999999999999991, 4.0592529337857286e-16, 8.8817841970012523e-16]
        numpy.testing.assert_array_almost_equal(ans,
                            magcoords.coordTrans([-3.608026916281573, 2.5673907444456745e-16,
            -1.7268878861662329], datetime.datetime(2009,1,1),'GSM','SM'), decimal=6)
        ans = [-0.8908435824705201, 0.09891125760552615, -0.4434120822553017]
        numpy.testing.assert_array_almost_equal(ans,
            magcoords.coordTrans([-4, 0, 1], datetime.datetime(2009,1,1), 'WGS84', 'GSM'), decimal=6)
        ans = [36.7677624344546, -2.329852593118382, 4.12310562561766]
        numpy.testing.assert_array_almost_equal(ans,
            magcoords.coordTrans([-4, 0, 1], datetime.datetime(2009,1,1), 'GSM', 'WGS84'), decimal=6)

    def test_Lvalue(self):
        """Lvalue should have known output"""
        ans = {'I': 10.899260245224534, 'L': 7.9509964756205083}
        vals = magcoords.Lvalue([-4, 0, 1], datetime.datetime(2009,1,1))
        for key in vals:
            self.assertAlmostEqual(vals[key], ans[key], places=4)
        ans = {'I': 0.69125074722140079, 'L': 4.5100680346100681}
        vals = magcoords.Lvalue([-4, 0, 1], datetime.datetime(2009,1,1), coord_system='SM')
        for key in vals:
            self.assertAlmostEqual(vals[key], ans[key], places=4)

    def test_Lvalue_extended_out(self):
        """Lvalue has an extended out (regression)"""
        ans = {'Blocal': 630.0509017184838,
         'Bmin': 22.033891107554684,
         'Bmirr': 630.0509017184838,
         'I': 10.899260245224534,
         'L': 7.9509964756205083,
         'M': 29966.895576135077}
        vals = magcoords.Lvalue([-4, 0, 1], datetime.datetime(2009,1,1), extended_out=True)
        for key in vals:
            self.assertAlmostEqual(vals[key], ans[key], places=4)

    def test_Lvalue_raises(self):
        """Test the exceptions in Lvalue"""
        self.assertRaises(NotImplementedError, magcoords.Lvalue,
                          [-4, 0, 1],
                          datetime.datetime(2009,1,1),
                          Bfield='bad Bfield')
        self.assertRaises(TypeError, magcoords.Lvalue,
                          [-4, 0, 1],
                          'bad datetime')
        self.assertRaises(NotImplementedError, magcoords.Lvalue,
                          [-4, 0, 1],
                          datetime.datetime(2009,1,1),
                          coord_system='bad coords')


    def test_input(self):
        """there is some input checking"""
        self.assertRaises(NotImplementedError, magcoords.coordTrans, [-3.608026916281573, 2.5673907444456745e-16,
            -1.7268878861662329], datetime.datetime(2009,1,1),'GSM','BAD')

    def test_Pin_pos_in(self):
        """I don't understand why this is here but here is a regression test"""
        in_pos = Lgm_Vector.Lgm_Vector(-4,0,0)
        ans = [-3.608026916281572, 4.163336342344337e-17, -1.7268878861662331]
        numpy.testing.assert_array_almost_equal(ans,
                            magcoords.coordTrans(in_pos, datetime.datetime(2009,1,1),'SM','GSM'), decimal=6)

    def test_raises(self):
        """coordTrans has various exceptions checked"""
        self.assertRaises(TypeError, magcoords.coordTrans, [-3.608026916281573, 2.5673907444456745e-16,
            -1.7268878861662329], 'bad date','SM','GSM')


if __name__ == '__main__':
    unittest.main()
