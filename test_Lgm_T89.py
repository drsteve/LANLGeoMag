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

import Lgm
from _Lgm import lib
import Lgm_Vector
import Lgm_T89


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

    def tearDown(self):
        super(Lgm_T89Tests, self).tearDown()

    def test_T89_1(self):
        """First simple in/out tests of T89 (regression)"""
        dt = datetime.datetime(2005, 8, 31, 9, 0, 0)
        pos = [-6.6, 0, 0]
        ans = [[-18.97193562594128, -1.8611995170538265, 80.3933831714847],
            [-20.828853435515278, -1.8611995170538265, 74.62222419125123 ],
            [-22.558420054640656, -1.8611995170538265, 70.52457956754591],
            [-26.412234393652312, -1.8611995170538265, 64.64506509634843],
            [-32.16112573540407, -1.8611995170538265, 60.078415300152216],
            [-45.379156657247805, -1.8611995170538265, 49.36315537639906] ]
        for i, kp in enumerate(range(6)):
            a = Lgm_T89.Lgm_T89(pos, dt, kp)
            B = a.calc_B()
            self.assertAlmostEqual(ans[i][0], B.x)
            self.assertAlmostEqual(ans[i][1], B.y)
            self.assertAlmostEqual(ans[i][2], B.z)

    def test_kp_checking(self):
        """for T89 Kp is between 0 and 5 inclusive"""
        dt = datetime.datetime(2005, 8, 31, 9, 0, 0)
        pos = [-6.6, 0, 0]
        self.assertRaises(ValueError, Lgm_T89.Lgm_T89, pos, dt, 10)
        self.assertRaises(ValueError, Lgm_T89.Lgm_T89, pos, dt, 6)
        self.assertRaises(ValueError, Lgm_T89.Lgm_T89, pos, dt, -1)
        self.assertRaises(ValueError, Lgm_T89.Lgm_T89, pos, dt, [-1, 0])
        self.assertRaises(ValueError, Lgm_T89.Lgm_T89, pos, dt, [7, 0])

    def test_list_in(self):
        """Make sure that list inputs work correctly"""
        dt = datetime.datetime(2005, 8, 31, 9, 0, 0)
        pos = [-6.6, 0, 0]
        ans = [[-18.97193562594128, -1.8611995170538265, 80.3933831714847],
            [-20.828853435515278, -1.8611995170538265, 74.62222419125123 ],
            [-22.558420054640656, -1.8611995170538265, 70.52457956754591],
            [-26.412234393652312, -1.8611995170538265, 64.64506509634843],
            [-32.16112573540407, -1.8611995170538265, 60.078415300152216],
            [-45.379156657247805, -1.8611995170538265, 49.36315537639906] ]
        for i, kp in enumerate(range(6)):
            a = Lgm_T89.Lgm_T89([pos]*2, [dt]*2, [kp]*2)
            B = a.calc_B()
            B = [val.tolist() for val in B]
            Bv = list(itertools.chain.from_iterable(B))
            Av = list(itertools.chain.from_iterable([ans[i]]*2))
            for i, val in enumerate(Bv):
                self.assertAlmostEqual(Av[i], val)



    def test_pos_checking(self):
        """Lgm_T89 pos agrument has checking"""
        pos = numpy.array([-6.6, 0, 0])
        kp = 0
        dt = datetime.datetime(2005, 8, 31, 9, 0, 0)
        self.assertRaises(TypeError, Lgm_T89.Lgm_T89, pos, dt, kp)

    def test_time_checking(self):
        """Lgm_T89 time agrument has checking"""
        pos = [-6.6, 0, 0]
        kp = 0
        dt = datetime.datetime(2005, 8, 31, 9, 0, 0)
        self.assertRaises(TypeError, Lgm_T89.Lgm_T89, pos, 'bad', kp)

    def test_intermal_model(self):
        """Lgm_T89 internal_model agrument has checking"""
        pos = [-6.6, 0, 0]
        kp = 0
        dt = datetime.datetime(2005, 8, 31, 9, 0, 0)
        self.assertRaises(ValueError, Lgm_T89.Lgm_T89, pos, dt, kp, INTERNAL_MODEL=4)
        self.assertRaises(ValueError, Lgm_T89.Lgm_T89, pos, dt, kp, INTERNAL_MODEL='bla')
        a = Lgm_T89.Lgm_T89(pos, dt, kp, INTERNAL_MODEL=1)
        self.assertEqual(a.INTERNAL_MODEL, 1)
        a = Lgm_T89.Lgm_T89(pos, dt, kp, INTERNAL_MODEL='LGM_CDIP')
        self.assertEqual(a.INTERNAL_MODEL, 0)
        a = Lgm_T89.Lgm_T89(pos, dt, kp, INTERNAL_MODEL='LGM_EDIP')
        self.assertEqual(a.INTERNAL_MODEL, 1)
        a = Lgm_T89.Lgm_T89(pos, dt, kp, INTERNAL_MODEL='LGM_IGRF')
        self.assertEqual(a.INTERNAL_MODEL, 2)

    def test_coord_system(self):
        """Lgm_T89 only inpout GSM for now (regression)"""
        pos = [-6.6, 0, 0]
        kp = 0
        dt = datetime.datetime(2005, 8, 31, 9, 0, 0)
        self.assertRaises(NotImplementedError, Lgm_T89.Lgm_T89, pos, dt, kp, coord_system='SM')

    def test_lengths(self):
        """all scalars or all lists"""
        pos = [-6.6, 0, 0]
        kp = 0
        dt = datetime.datetime(2005, 8, 31, 9, 0, 0)
        self.assertRaises(ValueError, Lgm_T89.Lgm_T89, [pos, pos], dt, kp)
        self.assertRaises(ValueError, Lgm_T89.Lgm_T89, [pos]*2, [dt]*2, [kp]*3)



if __name__ == '__main__':
    unittest.main()
