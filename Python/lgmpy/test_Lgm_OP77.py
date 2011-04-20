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

import numpy

import Lgm_OP77
import Lgm_Vector

class Lgm_OP77_OP77(unittest.TestCase):
    """
    Tests related to Lgm_OP77.OP77 wrapper
    @author: Brian Larsen
    @organization: LANL
    @contact: balarsen@lanl.gov

    @version: V1: 21-Mar-2011 (BAL)
    """
    def setUp(self):
        super(Lgm_OP77_OP77, self).setUp()
        self.pos = [-6.6, 0, 0]
        self.dt = datetime.datetime(2005, 8, 31, 9, 0, 0)
    def tearDown(self):
        super(Lgm_OP77_OP77, self).tearDown()

    def test_OP77(self):
        """the OP77 simple static wrapper should work (regression)"""
        numpy.testing.assert_array_almost_equal( Lgm_OP77.OP77(self.pos, self.dt, ),
            [-18.399869809240773, -1.8653978151945658, 85.951037824012488])
        ans = numpy.array([[-18.399869809240773, -1.8653978151945658, 85.951037824012488]*2])
        numpy.testing.assert_array_almost_equal(Lgm_OP77.OP77([self.pos]*2,
                                    [self.dt]*2), ans.reshape(2,3))

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
    def tearDown(self):
        super(Lgm_OP77Tests, self).tearDown()

    def test_pos2Lgm_Vector(self):
        """pos2Lgm_Vector should have known output"""
        a = Lgm_OP77.Lgm_OP77(self.pos, self.dt, )
        self.assertEqual(list(a['Position']), a._Vpos.tolist())
        self.assertTrue(isinstance(a._Vpos, Lgm_Vector.Lgm_Vector))
        b = Lgm_OP77.Lgm_OP77(a._Vpos, self.dt, )
        self.assertEqual(a._Vpos, b._Vpos)
        # above tested thought __init__ below is raw
        self.assertEqual(a._pos2Lgm_Vector([1, 2, 3]),
                         Lgm_Vector.Lgm_Vector(1, 2, 3) )

    def test_OP77_1(self):
        """First simple in/out tests of OP77 (regression)"""
        ans = numpy.array([-18.399869809240773, -1.8653978151945658, 85.951037824012488])
        B = Lgm_OP77.Lgm_OP77(self.pos, self.dt)
        numpy.testing.assert_array_almost_equal(ans, numpy.array(B['B'].tolist()))

    def test_list_in(self):
        """Make sure that list inputs work correctly (regression)"""
        ans = [-18.399869809240773, -1.8653978151945658, 85.95103782401249]
        a = Lgm_OP77.Lgm_OP77([self.pos]*2, [self.dt]*2)
        B = a.calc_B()
        B = [val.tolist() for val in B]
        Bv = list(itertools.chain.from_iterable(B))
        Av = list(itertools.chain.from_iterable([ans]*2))
        for i, val in enumerate(Bv):
            self.assertAlmostEqual(Av[i], val, 5)
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




if __name__ == '__main__':
    unittest.main()
