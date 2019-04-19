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
        ans = [-3.607999564981572, 1.3877763423443373e-17, -1.726945030766233]
        numpy.testing.assert_allclose(ans,
                            magcoords.coordTrans([-4,0,0], datetime.datetime(2009,1,1),
                            'SM','GSM'), atol=1e-8)
        ans = [-4., 0., 0.]
        numpy.testing.assert_allclose(ans,
                            magcoords.coordTrans([-3.607999564981572, 1.3877763423443373e-17,
            -1.726945030766233], datetime.datetime(2009,1,1),'GSM','SM'), atol=1e-8)
        #Use SphToCart/CartToSph to convert between the C values (always
        #cartesian even in WGS84) and the Python ones (WGS84 is spherical)
        ans = [-0.8908435824709999, 0.09890080939498999, -0.443414412792]
        numpy.testing.assert_allclose(ans,
            magcoords.coordTrans([-4, 0, 1], datetime.datetime(2009,1,1),
            'WGS84', 'GSM'), rtol=1e-6)
        ans = [36.76783195260621, -2.3294531575340423, 4.1231056256205525]
        numpy.testing.assert_allclose(ans,
            magcoords.coordTrans([-4, 0, 1], datetime.datetime(2009,1,1),
            'GSM', 'WGS84'), rtol=1e-6)

    def test_end2end_jpl_eph(self):
        """quick functional tests with JPL DE421"""
        #Straight from the C
        cases = ["2015-03-08T19:18:31.786778 479114378970778386 GEI2000 GSE -72417.672749 -51921.225510 -29798.643397 -58063.063795 -73566.510184 -6688.896414",
                 "2015-03-11T11:58:40.170326 479347187354325578 GSE2000 SM -20164.747998 -52414.140552 -20241.497772 -20870.283837 -46078.983969 -31697.569507",
                 "2015-03-12T04:25:43.605928 479406410789928239 SM GSM 75645.819112 -40422.353445 -35182.256254 81628.221639 -40422.353445 -17231.208208"]
        for c in cases:
            dt, tt, fromsys, tosys, xin, yin, zin, xout, yout, zout = c.split()
            dt = datetime.datetime.strptime(dt, '%Y-%m-%dT%H:%M:%S.%f')
            inputs = [float(xin), float(yin), float(zin)]
            expected = [float(xout), float(yout), float(zout)]
            actual = magcoords.coordTrans(
                inputs, dt, fromsys, tosys, de_eph=True)
            numpy.testing.assert_allclose(expected, actual, atol=1e-8)

    def test_Lvalue(self):
        """Lvalue should have known output"""
        #Default Kp is 2 in both the C and the Python
        ans = {'I': 10.9312487946278, 'L': 7.96706218379283}
        vals = magcoords.Lvalue([-4, 0, 1], datetime.datetime(2009,1,1))
        for key in vals:
            self.assertAlmostEqual(vals[key], ans[key], places=4)
        ans = {'I': 0.69161497288978, 'L': 4.51513147096729}
        vals = magcoords.Lvalue([-4, 0, 1], datetime.datetime(2009,1,1), coord_system='SM')
        for key in vals:
            self.assertAlmostEqual(vals[key], ans[key], places=4)

    def test_Lvalue_extended_out(self):
        """Lvalue has an extended out (regression)"""
        ans = {
         'Blocal': 628.07992515195,
         'Bmin': 21.678823206256,
         'Bmirr': 628.07992515195,
         'I': 10.9312487946278,
         'L': 7.96706218379283,
         'M': 29966.8851195941,
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
        ans = [-3.607999564981572, 1.3877763423443373e-17, -1.726945030766233]
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

    def test_Lvalues_lstar(self):
        """Test Hilton/McIlwain for the values used in L* tests"""
        #testLgm_B_T89_1
        L = magcoords.Lvalue([-4.2, 1, 1], datetime.datetime(2010, 10, 12),
                             alpha=90, Kp=4, method='McIlwain',
                             coord_system='SM', Bfield='Lgm_B_T89')['L']
        self.assertEqual(5.0308451893, L)


if __name__ == '__main__':
    unittest.main()
