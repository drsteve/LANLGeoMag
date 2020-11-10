#!/usr/bin/env python

"""
Test suite for the Lgm_T89 file

@author: Brian Larsen
@organization: LANL
@contact: balarsen@lanl.gov

@version: V1: 20-Dec-2010 (BAL)
"""

import unittest
import datetime
import itertools

import numpy

from lgmpy import Lgm_T89
from lgmpy import Lgm_Vector

class Lgm_T89_T89_tests(unittest.TestCase):
    """
    Tests related to Lgm_T89.T89 wrapper
    @author: Brian Larsen
    @organization: LANL
    @contact: balarsen@lanl.gov

    @version: V1: 11-Jan-2011 (BAL)
    """
    def setUp(self):
        super(Lgm_T89_T89_tests, self).setUp()
        self.pos = [-6.6, 0, 0]
        self.kp = 0
        self.dt = datetime.datetime(2005, 8, 31, 9, 0, 0)
    def tearDown(self):
        super(Lgm_T89_T89_tests, self).tearDown()

    def test_T89(self):
        """the T89 simple static wrapper should work (regression)"""
        ans = [-18.926081,  -1.857515,  80.052116]
        numpy.testing.assert_allclose( Lgm_T89.T89(self.pos, self.dt, self.kp),
            [-18.926081,  -1.857515,  80.052116], rtol=1e-3, atol=0)
        ans = numpy.array([ans*2])
        numpy.testing.assert_allclose(Lgm_T89.T89([self.pos]*2,
                                    [self.dt]*2,
                                    [self.kp]*2) , ans.reshape(2,3), rtol=1e-3, atol=0)

class Lgm_T89Tests(unittest.TestCase):
    """
    Tests related to Lgm_T89
    @author: Brian Larsen
    @organization: LANL
    @contact: balarsen@lanl.gov

    @version: V1: 04-Jan-2011 (BAL)
    """
    def setUp(self):
        super(Lgm_T89Tests, self).setUp()
        self.pos = [-6.6, 0, 0]
        self.kp = 0
        self.dt = datetime.datetime(2005, 8, 31, 9, 0, 0)
        self.ans = [[-18.927103, -1.857488, 80.052164],
                    [-20.78416, -1.857488, 74.281026],
                    [-22.513859, -1.857488, 70.183396],
                    [-26.367951, -1.857488, 64.303916],
                    ]

    def tearDown(self):
        super(Lgm_T89Tests, self).tearDown()

    def test_T89_1(self):
        """First simple in/out tests of T89 (compared to C)"""
        for i, kp in enumerate(range(4)):
            B = Lgm_T89.Lgm_T89(self.pos, self.dt, kp)
            out = [B['B'].x, B['B'].y, B['B'].z]
            numpy.testing.assert_allclose(self.ans[i], out, atol=1e-5)

    def test_kp_checking(self):
        """for T89 Kp is between 0 and 5 inclusive"""
        self.assertRaises(ValueError, Lgm_T89.Lgm_T89, self.pos, self.dt, 10)
        self.assertRaises(ValueError, Lgm_T89.Lgm_T89, self.pos, self.dt, 6)
        self.assertRaises(ValueError, Lgm_T89.Lgm_T89, self.pos, self.dt, -1)
        self.assertRaises(ValueError, Lgm_T89.Lgm_T89, self.pos, self.dt, [-1, 0])
        self.assertRaises(ValueError, Lgm_T89.Lgm_T89, self.pos, self.dt, [7, 0])

    def test_list_in(self):
        """Make sure that list inputs work correctly (regression)"""
        for i, kp in enumerate(range(4)):
            a = Lgm_T89.Lgm_T89([self.pos]*2, [self.dt]*2, [kp]*2)
            B = a.calc_B()
            B = [val.tolist() for val in B]
            Bv = list(itertools.chain.from_iterable(B))
            numpy.testing.assert_allclose(numpy.tile(self.ans[i], 2), Bv, atol=1e-5)
        self.assertEqual(len(B), 2)

    def test_pos_checking(self):
        """Lgm_T89 pos argument has checking"""
        self.assertRaises(TypeError, Lgm_T89.Lgm_T89, 'bad', self.dt, self.kp)

    def test_time_checking(self):
        """Lgm_T89 time argument has checking"""
        self.assertRaises(TypeError, Lgm_T89.Lgm_T89, self.pos, 'bad', self.kp)

    def test_internal_model(self):
        """Lgm_T89 internal_model argument has checking"""
        self.assertRaises(ValueError, Lgm_T89.Lgm_T89, self.pos, self.dt,
                          self.kp, INTERNAL_MODEL=4)
        self.assertRaises(ValueError, Lgm_T89.Lgm_T89, self.pos, self.dt,
                          self.kp, INTERNAL_MODEL='bla')
        a = Lgm_T89.Lgm_T89(self.pos, self.dt, self.kp, INTERNAL_MODEL=1)
        self.assertEqual(a.attrs['internal_model'], 1)
        a = Lgm_T89.Lgm_T89(self.pos, self.dt, self.kp, INTERNAL_MODEL='LGM_CDIP')
        self.assertEqual(a.attrs['internal_model'], 0)
        a = Lgm_T89.Lgm_T89(self.pos, self.dt, self.kp, INTERNAL_MODEL='LGM_EDIP')
        self.assertEqual(a.attrs['internal_model'], 1)
        a = Lgm_T89.Lgm_T89(self.pos, self.dt, self.kp, INTERNAL_MODEL='LGM_IGRF')
        self.assertEqual(a.attrs['internal_model'], 2)

    def test_internal_model_values(self):
        """Lgm_T89 internal_model agrument changes the ouput values"""
        # values form C test run
        #        Date = 20090101; UTC  = 0.0; mInfo->Kp = 2; 
        #        u.x = 3.0; u.y = 4.0; u.z = 3.0
        ans = {}
        ans['IGRF'] = [-123.994193, -73.755223, 77.986257]
        ans['CDIP'] = [-123.603627, -75.262655, 77.740804]
        ans['EDIP'] = [-122.127541, -75.432165, 79.463172]
        pos = [3.0, 4.0, 3.0]
        dt = datetime.datetime(2009, 1, 1, 0)
        kp = 2
        a = Lgm_T89.T89(pos, dt, kp, INTERNAL_MODEL='LGM_IGRF')
        numpy.testing.assert_allclose(ans['IGRF'], a)
        a = Lgm_T89.T89(pos, dt, kp, INTERNAL_MODEL='LGM_CDIP')
        numpy.testing.assert_allclose(ans['CDIP'], a)
        a = Lgm_T89.T89(pos, dt, kp, INTERNAL_MODEL='LGM_EDIP')
        numpy.testing.assert_allclose(ans['EDIP'], a)

    def test_coord_system(self):
        """Lgm_T89 only inpout GSM for now (regression)"""
        self.assertRaises(NotImplementedError, Lgm_T89.Lgm_T89, self.pos,
                          self.dt, self.kp, coord_system='SM')

    def test_lengths(self):
        """all scalars or all lists"""
        self.assertRaises(ValueError, Lgm_T89.Lgm_T89,
                          [self.pos, self.pos], self.dt, self.kp)
        self.assertRaises(ValueError, Lgm_T89.Lgm_T89,
                          [self.pos]*2, [self.dt]*2, [self.kp]*3)



if __name__ == '__main__':
    unittest.main()
