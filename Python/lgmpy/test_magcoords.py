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

from lgmpy import magcoords
from lgmpy import Lgm_Vector

class magcoords_Tests(unittest.TestCase):
    def setUp(self):
        super(magcoords_Tests, self).setUp()

    def tearDown(self):
        super(magcoords_Tests, self).tearDown()

    def test_end2end(self):
        """quick functional tests (regression)"""
        ans = [-3.608026916281572, 4.163336342344337e-17, -1.7268878861662331]
        numpy.testing.assert_allclose(ans,
                            magcoords.coordTrans([-4,0,0], datetime.datetime(2009,1,1),
                            'SM','GSM'), atol=1e-8)
        ans = [-3.9999999999999991, 1.87350135e-16, 0.00000000e+00]
        numpy.testing.assert_allclose(ans,
                            magcoords.coordTrans([-3.608026916281573, 2.5673907444456745e-16,
            -1.7268878861662329], datetime.datetime(2009,1,1),'GSM','SM'), atol=1e-8)
        ans = [-0.8908435824705201, 0.09891125760552615, -0.4434120822553017]
        numpy.testing.assert_allclose(ans,
            magcoords.coordTrans([-4, 0, 1], datetime.datetime(2009,1,1),
            'WGS84', 'GSM'), rtol=1e-6)
        ans = [36.7677624344546, -2.329852593118382, 4.12310562561766]
        numpy.testing.assert_allclose(ans,
            magcoords.coordTrans([-4, 0, 1], datetime.datetime(2009,1,1),
            'GSM', 'WGS84'), rtol=1e-6)

    def test_Lvalue(self):
        """Lvalue should have known output"""
        ans = {'I': 10.9298583451352, 'L': 7.966548339}
        vals = magcoords.Lvalue([-4, 0, 1], datetime.datetime(2009,1,1))
        for key in vals:
            self.assertAlmostEqual(vals[key], ans[key], places=4)
        ans = {'I': 0.6915739190021, 'L': 4.51510833989}
        vals = magcoords.Lvalue([-4, 0, 1], datetime.datetime(2009,1,1), coord_system='SM')
        for key in vals:
            self.assertAlmostEqual(vals[key], ans[key], places=4)

    def test_Lvalue_extended_out(self):
        """Lvalue has an extended out (regression)"""
        ans = {'Blocal': 628.0753565024,
         'Bmin': 21.686322939,
         'Bmirr': 628.0753565024,
         'I': 10.929858345135273,
         'L': 7.9665483393055,
         'M': 29966.895576135077,
         'MLon': 180.0,
         'MLT': 0.0,
         }
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
        self.assertRaises(KeyError, magcoords.Lvalue,
                          [-4, 0, 1],
                          datetime.datetime(2009,1,1),
                          coord_system='bad coords')


    def test_input(self):
        """there is some input checking"""
        self.assertRaises(KeyError, magcoords.coordTrans, [-3.608026916281573, 2.5673907444456745e-16,
            -1.7268878861662329], datetime.datetime(2009,1,1),'GSM','BAD')

    def test_Pin_pos_in(self):
        """I don't understand why this is here but here is a regression test"""
        in_pos = Lgm_Vector.Lgm_Vector(-4,0,0)
        ans = [-3.608026916281572, 4.163336342344337e-17, -1.7268878861662331]
        numpy.testing.assert_allclose(ans,
                            magcoords.coordTrans(in_pos, datetime.datetime(2009,1,1),'SM','GSM'), atol=1e-8)

    def test_raises(self):
        """coordTrans has various exceptions checked"""
        self.assertRaises(TypeError, magcoords.coordTrans, [-3.608026916281573, 2.5673907444456745e-16,
            -1.7268878861662329], 'bad date','SM','GSM')

    def test_GEO_GSE(self):
        """regression test on GEO_TO_GSM"""
        expected = [-1.4268697964755936, -0.6973756761170802, 3.3878768794432226]
        numpy.testing.assert_allclose(expected,
                            magcoords.coordTrans([1,2,3], datetime.datetime(2012, 3, 4), 'GEO', 'GSE'))
        numpy.testing.assert_allclose([1,2,3],
                            magcoords.coordTrans(expected, datetime.datetime(2012, 3, 4), 'GSE', 'GEO'))



if __name__ == '__main__':
    unittest.main()
