#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Test suite for the Lgm_OP77 file

@author: Brian Larsen
@organization: LANL
@contact: balarsen@lanl.gov

@version: V1: 20-Dec-2010 (BAL)
"""

import unittest
import datetime
import itertools

import numpy as np

from lgmpy import Lgm_OP77
from lgmpy import Lgm_Vector

class Lgm_OP77_OP77_tests(unittest.TestCase):
    """
    Tests related to Lgm_OP77.OP77 wrapper
    @author: Brian Larsen
    @organization: LANL
    @contact: balarsen@lanl.gov

    @version: V1: 21-Mar-2011 (BAL)
    """
    def setUp(self):
        super(Lgm_OP77_OP77_tests, self).setUp()
        self.pos = [-6.6, 0, 0]
        self.dt = datetime.datetime(2005, 8, 31, 9, 0, 0)
    def tearDown(self):
        super(Lgm_OP77_OP77_tests, self).tearDown()

    def test_OP77(self):
        """the OP77 simple static wrapper should work (regression)"""
        np.testing.assert_allclose( Lgm_OP77.OP77(self.pos, self.dt, ),
            [-18.34951, -1.857515, 85.610048], rtol=1e-3, atol=0)
        ans = np.array([[-18.34951, -1.857515, 85.610048]*2])
        np.testing.assert_allclose(Lgm_OP77.OP77([self.pos]*2,
                                    [self.dt]*2), ans.reshape(2,3), rtol=1e-3, atol=0)


class Lgm_OP77Tests(unittest.TestCase):
    """
    Tests related to Lgm_OP77
    @author: Brian Larsen
    @organization: LANL
    @contact: balarsen@lanl.gov

    @version: V1: 21-Mar-2011 (BAL)
    """
    def setUp(self):
        super(Lgm_OP77Tests, self).setUp()
        self.pos = [-6.6, 0, 0]
        self.dt = datetime.datetime(2005, 8, 31, 9, 0, 0)
        self.ans = np.array([-18.350502, -1.857488, 85.610065])

    def tearDown(self):
        super(Lgm_OP77Tests, self).tearDown()

    def test_OP77_1(self):
        """First simple in/out tests of OP77 (regression)"""
        B = Lgm_OP77.Lgm_OP77(self.pos, self.dt)
        np.testing.assert_allclose(self.ans, np.array(B['B'].tolist()),
                                   rtol=1e-8, atol=1e-5)

    def test_list_in(self):
        """Make sure that list inputs work correctly (regression)"""
        a = Lgm_OP77.Lgm_OP77([self.pos]*2, [self.dt]*2)
        B = a.calc_B()
        B = [val.tolist() for val in B]
        Bv = list(itertools.chain.from_iterable(B))
        np.testing.assert_allclose(np.tile(self.ans, 2), Bv, atol=1e-5)
        self.assertEqual(len(B), 2)

    def test_pos_checking(self):
        """Lgm_OP77 pos agrument has checking"""
        self.assertRaises(TypeError, Lgm_OP77.Lgm_OP77, 'bad', self.dt, )

    def test_time_checking(self):
        """Lgm_OP77 time agrument has checking"""
        self.assertRaises(TypeError, Lgm_OP77.Lgm_OP77, self.pos, 'bad', )

    def test_internal_model(self):
        """Lgm_OP77 internal_model agrument has checking"""
        self.assertRaises(ValueError, Lgm_OP77.Lgm_OP77, self.pos, self.dt,
                           INTERNAL_MODEL=4)
        self.assertRaises(ValueError, Lgm_OP77.Lgm_OP77, self.pos, self.dt,
                           INTERNAL_MODEL='bla')
        a = Lgm_OP77.Lgm_OP77(self.pos, self.dt,  INTERNAL_MODEL=1)
        self.assertEqual(a.attrs['internal_model'], 1)
        a = Lgm_OP77.Lgm_OP77(self.pos, self.dt,  INTERNAL_MODEL='LGM_CDIP')
        self.assertEqual(a.attrs['internal_model'], 0)
        a = Lgm_OP77.Lgm_OP77(self.pos, self.dt,  INTERNAL_MODEL='LGM_EDIP')
        self.assertEqual(a.attrs['internal_model'], 1)
        a = Lgm_OP77.Lgm_OP77(self.pos, self.dt,  INTERNAL_MODEL='LGM_IGRF')
        self.assertEqual(a.attrs['internal_model'], 2)

    def test_coord_system(self):
        """Lgm_OP77 only inpout GSM for now (regression)"""
        self.assertRaises(NotImplementedError, Lgm_OP77.Lgm_OP77, self.pos,
                          self.dt,  coord_system='SM')

    def test_lengths(self):
        """all scalars or all lists"""
        self.assertRaises(ValueError, Lgm_OP77.Lgm_OP77,
                          [self.pos, self.pos], self.dt, )

    def test_exceptions(self):
        """the OP77 __init__ does some input checking"""
        self.assertRaises(TypeError, Lgm_OP77.Lgm_OP77, 'bad pos', self.dt)
        self.assertRaises(ValueError, Lgm_OP77.Lgm_OP77, [self.pos]*2, [self.dt])


if __name__ == '__main__':
    unittest.main()
