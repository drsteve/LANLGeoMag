#!/usr/bin/env python

"""
Test suite for the Lgm_CTrans file <<This is an important one>>

@author: Brian Larsen
@organization: LANL
@contact: balarsen@lanl.gov

@version: V1: 20-Dec-2010 (BAL)
"""

import unittest
import _Lgm_CTrans
import _Lgm_Vector
import _Lgm
import math

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
        self.lgm = _Lgm._Lgm()

    def tearDown(self):
        super(Lgm_CTransTests, self).tearDown()

    def test_Lgm_Convert_Coords(self):
        """Lgm_Convert_Coords should give known output (regression)"""
        # taken from CoordQuickStart.c
        Date = 20040812
        UTC = 12.34567
        Ugsm = _Lgm_Vector.Lgm_Vector(x = -6.6, y = 3.4, z = -2.3)
        Usm = _Lgm_Vector.Lgm_Vector()
        c = _Lgm_CTrans.Lgm_CTrans()
        self.lgm.lib.Lgm_Set_Coord_Transforms(Date, UTC, c)
        self.lgm.lib.Lgm_Convert_Coords(Ugsm, Usm, _Lgm_CTrans.GSM_TO_SM, c)
        self.assertAlmostEqual(-5.5352494753370127, Usm.x)
        self.assertAlmostEqual( 3.3999999999999995, Usm.y)
        self.assertAlmostEqual(-4.2674363786448328, Usm.z)

    def test_Lgm_Convert_CoordsBack(self):
        """Lgm_Convert_Coords should give known output (and back)"""
        # taken from CoordQuickStart.c
        Date = 20040812
        UTC = 12.34567
        Ugsm = _Lgm_Vector.Lgm_Vector(x = -6.6, y = 3.4, z = -2.3)
        Usm = _Lgm_Vector.Lgm_Vector()
        c = _Lgm_CTrans.Lgm_CTrans()
        self.lgm.lib.Lgm_Set_Coord_Transforms(Date, UTC, c)
        self.lgm.lib.Lgm_Convert_Coords(Ugsm, Usm, _Lgm_CTrans.GSM_TO_SM, c)
        self.assertAlmostEqual(-5.5352494753370127, Usm.x)
        self.assertAlmostEqual( 3.3999999999999995, Usm.y)
        self.assertAlmostEqual(-4.2674363786448328, Usm.z)
        self.lgm.lib.Lgm_Convert_Coords(Usm, Ugsm, _Lgm_CTrans.SM_TO_GSM, c)
        self.assertAlmostEqual(-6.6, Ugsm.x)
        self.assertAlmostEqual( 3.4, Ugsm.y)
        self.assertAlmostEqual(-2.3, Ugsm.z)


if __name__ == '__main__':
    unittest.main()
