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

import Lgm
import Lgm_CTrans
import Lgm_Vector
from _Lgm import lib

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
        #print("c: %d" % (lib.size_CTrans()) )
        self.assertEqual(ctypes.sizeof(Lgm_CTrans.Lgm_CTrans), lib.size_CTrans())

    def test_sizeLgm_LeapSeconds(self):
        """Lgm_LeapSeconds c/python must be same size"""
        self.assertEqual(ctypes.sizeof(Lgm_CTrans.Lgm_LeapSeconds),
                         lib.size_Lgm_LeapSeconds())

    def test_sizeDateTime(self):
        """DateTime c and python should match size"""
        self.assertEqual(ctypes.sizeof(Lgm_CTrans.Lgm_DateTime), lib.size_DateTime())

    def test_Lgm_DateTime(self):
        """Lgm_DateTime should create and set defaults"""
        a = Lgm_CTrans.Lgm_DateTime()
        self.assertEqual(a.nNutationTerms, 106)
        self.assertEqual(a.DUT1, 0.0)
        self.assertEqual(a.ddPsi, 0.0)

    def test_Lgm_Convert_Coords(self):
        """Lgm_Convert_Coords should give known output (regression)"""
        # taken from CoordQuickStart.c
        Date = 20040812
        UTC = 12.34567
        Ugsm = Lgm_Vector.Lgm_Vector(x = -6.6, y = 3.4, z = -2.3)
        Usm = Lgm_Vector.Lgm_Vector()
        c = Lgm_CTrans.Lgm_CTrans()
        lib.Lgm_Set_Coord_Transforms(Date, UTC, c)
        lib.Lgm_Convert_Coords(Ugsm, Usm, Lgm_CTrans.GSM_TO_SM, c)
        self.assertAlmostEqual(-5.5352494753370127, Usm.x)
        self.assertAlmostEqual( 3.3999999999999995, Usm.y)
        self.assertAlmostEqual(-4.2674363786448328, Usm.z)

    def test_Lgm_Convert_CoordsBack(self):
        """Lgm_Convert_Coords should give known output (and back)"""
        # taken from CoordQuickStart.c
        Date = 20040812
        UTC = 12.34567
        Ugsm = Lgm_Vector.Lgm_Vector(x = -6.6, y = 3.4, z = -2.3)
        Usm = Lgm_Vector.Lgm_Vector()
        c = Lgm_CTrans.Lgm_CTrans()
        lib.Lgm_Set_Coord_Transforms(Date, UTC, c)
        lib.Lgm_Convert_Coords(Ugsm, Usm, Lgm_CTrans.GSM_TO_SM, c)
        self.assertAlmostEqual(-5.5352494753370127, Usm.x)
        self.assertAlmostEqual( 3.3999999999999995, Usm.y)
        self.assertAlmostEqual(-4.2674363786448328, Usm.z)
        lib.Lgm_Convert_Coords(Usm, Ugsm, Lgm_CTrans.SM_TO_GSM, c)
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


if __name__ == '__main__':
    unittest.main()
