#!/usr/bin/env python

"""
Test suite for the Lgm_CTrans file <<This is an important one>>

@author: Brian Larsen
@organization: LANL
@contact: balarsen@lanl.gov

@version: V1: 20-Dec-2010 (BAL)
"""

import ctypes
import unittest
import datetime

import numpy

import Lgm_CTrans
import Lgm_Vector
from Lgm_Wrap import Lgm_LeapSeconds, size_Lgm_LeapSeconds, size_CTrans
from Lgm_Wrap import Lgm_DateTime, size_DateTime, Lgm_Set_Coord_Transforms
from Lgm_Wrap import Lgm_Convert_Coords, GSM_TO_SM, SM_TO_GSM

class Lgm_CoordsTests(unittest.TestCase):
    """
    Tests related to Lgm_CTrans.Lgm_Coords
    @author: Brian Larsen
    @organization: LANL
    @contact: balarsen@lanl.gov

    @version: V1: 20-Dec-2010 (BAL)
    """
    def setUp(self):
        super(Lgm_CoordsTests, self).setUp()

    def tearDown(self):
        super(Lgm_CoordsTests, self).tearDown()

    def test_init(self):
        """Lgm_Coords does some input checking"""
        self.assertRaises(NotImplementedError, Lgm_CTrans.Lgm_Coords, [1, 2, 3], units='km')
        self.assertRaises(NotImplementedError, Lgm_CTrans.Lgm_Coords, [1, 2, 3], system='bad')
        self.assertRaises(ValueError, Lgm_CTrans.Lgm_Coords, [0.1, 0.5, 0.6])
        coords = Lgm_CTrans.Lgm_Coords([1, 2, 3])
        self.assertEqual(coords[:], [1,2,3])
        self.assertEqual(coords.system, 'GSM') # with will catch a default change

class Lgm_CTransTests(unittest.TestCase):
    """
    Tests related to Lgm_CTrans
    @author: Brian Larsen
    @organization: LANL
    @contact: balarsen@lanl.gov

    @version: V1: 20-Dec-2010 (BAL)
    """
    def setUp(self):
        super(Lgm_CTransTests, self).setUp()

    def tearDown(self):
        super(Lgm_CTransTests, self).tearDown()

    def test_sizeLgm_CTrans(self):
        """Lgm_CTrans c and python must have same size"""
        #print("c: %d" % (size_CTrans()) )
        self.assertEqual(ctypes.sizeof(Lgm_CTrans.Lgm_CTrans), size_CTrans())

    def test_sizeLgm_LeapSeconds(self):
        """Lgm_LeapSeconds c/python must be same size"""
        self.assertEqual(ctypes.sizeof(Lgm_LeapSeconds),
                         size_Lgm_LeapSeconds())

    def test_sizeDateTime(self):
        """DateTime c and python should match size"""
        self.assertEqual(ctypes.sizeof(Lgm_DateTime), size_DateTime())

    def test_Lgm_Convert_CoordsBack(self):
        """Lgm_Convert_Coords should give known output (and back)"""
        # taken from CoordQuickStart.c
        Date = 20040812
        UTC = 12.34567
        Ugsm = Lgm_Vector.Lgm_Vector(x = -6.6, y = 3.4, z = -2.3)
        Usm = Lgm_Vector.Lgm_Vector()
        c = Lgm_CTrans.Lgm_CTrans()
        Lgm_Set_Coord_Transforms(Date, UTC, ctypes.pointer(c))
        Lgm_Convert_Coords(ctypes.pointer(Ugsm), ctypes.pointer(Usm),
                           GSM_TO_SM, ctypes.pointer(c))
        self.assertAlmostEqual(-5.5352494753370127, Usm.x, places=5)
        self.assertAlmostEqual( 3.3999999999999995, Usm.y, places=5)
        self.assertAlmostEqual(-4.2674363786448328, Usm.z, places=5)
        Lgm_Convert_Coords(ctypes.pointer(Usm), ctypes.pointer(Ugsm),
                           SM_TO_GSM, ctypes.pointer(c))
        self.assertAlmostEqual(-6.6, Ugsm.x)
        self.assertAlmostEqual( 3.4, Ugsm.y)
        self.assertAlmostEqual(-2.3, Ugsm.z)

    def test_dateToDateLong(self):
        """dateToDateLong should give known results for known input"""
        d1 = datetime.datetime(2000, 12, 12)
        self.assertEqual(20001212L, Lgm_CTrans.dateToDateLong(d1))
        self.assertEqual([20001212L, 20001212L], Lgm_CTrans.dateToDateLong([d1, d1]))
        self.assertEqual([20001212L, 20001212L],
            Lgm_CTrans.dateToDateLong(numpy.array([d1, d1])).tolist())
        self.assertEqual(20001212L, Lgm_CTrans.dateToDateLong(d1.date()))

    def test_dateToFPHours(self):
        """dateToFPHours should give known output for known input"""
        d1 = datetime.datetime(2000, 12, 12)
        self.assertEqual(0.0, Lgm_CTrans.dateToFPHours(d1))
        for val in range(24):
            self.assertEqual(val,
                Lgm_CTrans.dateToFPHours(datetime.datetime(2000, 12, 12, val)))
        self.assertEqual(12.5,
            Lgm_CTrans.dateToFPHours(datetime.datetime(2000, 12, 12, 12, 30)))
        self.assertEqual(12.25,
            Lgm_CTrans.dateToFPHours(datetime.datetime(2000, 12, 12, 12, 15)))
        self.assertEqual(12.75,
            Lgm_CTrans.dateToFPHours(datetime.datetime(2000, 12, 12, 12, 45)))
        d1 = datetime.datetime(2000, 12, 12, 12, 30)
        self.assertEqual([12.5, 12.5], Lgm_CTrans.dateToFPHours([d1, d1]))
        self.assertEqual([12.5, 12.5],
            Lgm_CTrans.dateToFPHours(numpy.array([d1, d1])).tolist())

    def test_init_vals(self):
        """__init__ sets upa few valiues (regression)"""
        tst = Lgm_CTrans.Lgm_CTrans()
        self.assertEqual(tst.nNutationTerms, 106)
        self.assertEqual(tst.DUT1, 0.0)
        self.assertEqual(tst.xp, 0.0)
        self.assertEqual(tst.yp, 0.0)
        self.assertEqual(tst.ddPsi, 0.0)
        self.assertEqual(tst.ddEps, 0.0)


if __name__ == '__main__':
    unittest.main()
