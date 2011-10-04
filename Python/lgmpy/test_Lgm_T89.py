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

import Lgm_T89
import Lgm_Vector

class Lgm_T89_T89(unittest.TestCase):
    """
    Tests related to Lgm_T89.T89 wrapper
    @author: Brian Larsen
    @organization: LANL
    @contact: balarsen@lanl.gov

    @version: V1: 11-Jan-2011 (BAL)
    """
    def setUp(self):
        super(Lgm_T89_T89, self).setUp()
        self.pos = [-6.6, 0, 0]
        self.kp = 0
        self.dt = datetime.datetime(2005, 8, 31, 9, 0, 0)
    def tearDown(self):
        super(Lgm_T89_T89, self).tearDown()

    def test_T89(self):
        """the T89 simple static wrapper should work (regression)"""
        numpy.testing.assert_allclose( Lgm_T89.T89(self.pos, self.dt, self.kp),
            [-18.976439213243122, -1.8653978086481349, 80.39310505873112])
        ans = numpy.array([[-18.976439213243122, -1.8653978086481349, 80.39310505873112]*2])
        numpy.testing.assert_allclose(Lgm_T89.T89([self.pos]*2,
                                    [self.dt]*2,
                                    [self.kp]*2) , ans.reshape(2,3))

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
    def tearDown(self):
        super(Lgm_T89Tests, self).tearDown()

    def test_pos2Lgm_Vector(self):
        """pos2Lgm_Vector should have known output"""
        a = Lgm_T89.Lgm_T89(self.pos, self.dt, self.kp)
        self.assertEqual(list(a['Position']), a._Vpos.tolist())
        self.assertTrue(isinstance(a._Vpos, Lgm_Vector.Lgm_Vector))
        b = Lgm_T89.Lgm_T89(a._Vpos, self.dt, self.kp)
        self.assertEqual(a._Vpos, b._Vpos)
        # above tested thought __init__ below is raw
        self.assertEqual(a._pos2Lgm_Vector([1, 2, 3]),
                         Lgm_Vector.Lgm_Vector(1, 2, 3) )
        self.assertRaises(NotImplementedError, a._pos2Lgm_Vector,
                          numpy.array([1, 2, 3]))

    def test_T89_1(self):
        """First simple in/out tests of T89 (regression)"""
        ans = [[-18.976439213243122, -1.8653978086481349, 80.39310505873112],
            [-20.833382716404937, -1.8653978086481349, 74.62194986821649 ],
            [-22.562973770829664, -1.8653978086481349, 70.52430794391046],
            [-26.416839227182663, -1.8653978086481349, 64.64479976507458], ]
        for i, kp in enumerate(range(4)):
            B = Lgm_T89.Lgm_T89(self.pos, self.dt, kp)
            self.assertAlmostEqual(ans[i][0], B['B'].x, 5)
            self.assertAlmostEqual(ans[i][1], B['B'].y, 5)
            self.assertAlmostEqual(ans[i][2], B['B'].z, 5)

    def test_kp_checking(self):
        """for T89 Kp is between 0 and 5 inclusive"""
        self.assertRaises(ValueError, Lgm_T89.Lgm_T89, self.pos, self.dt, 10)
        self.assertRaises(ValueError, Lgm_T89.Lgm_T89, self.pos, self.dt, 6)
        self.assertRaises(ValueError, Lgm_T89.Lgm_T89, self.pos, self.dt, -1)
        self.assertRaises(ValueError, Lgm_T89.Lgm_T89, self.pos, self.dt, [-1, 0])
        self.assertRaises(ValueError, Lgm_T89.Lgm_T89, self.pos, self.dt, [7, 0])

    def test_list_in(self):
        """Make sure that list inputs work correctly (regression)"""
        ans = [[-18.976439213243122, -1.8653978086481349, 80.39310505873112],
            [-20.833382716404937, -1.8653978086481349, 74.62194986821649 ],
            [-22.562973770829664, -1.8653978086481349, 70.52430794391046],
            [-26.416839227182663, -1.8653978086481349, 64.64479976507458], ]
        for i, kp in enumerate(range(4)):
            a = Lgm_T89.Lgm_T89([self.pos]*2, [self.dt]*2, [kp]*2)
            B = a.calc_B()
            B = [val.tolist() for val in B]
            Bv = list(itertools.chain.from_iterable(B))
            Av = list(itertools.chain.from_iterable([ans[i]]*2))
            for i, val in enumerate(Bv):
                self.assertAlmostEqual(Av[i], val, 5)
        self.assertEqual(len(B), 2)

    def test_pos_checking(self):
        """Lgm_T89 pos agrument has checking"""
        self.assertRaises(TypeError, Lgm_T89.Lgm_T89, 'bad', self.dt, self.kp)

    def test_time_checking(self):
        """Lgm_T89 time agrument has checking"""
        self.assertRaises(TypeError, Lgm_T89.Lgm_T89, self.pos, 'bad', self.kp)

    def test_internal_model(self):
        """Lgm_T89 internal_model agrument has checking"""
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
