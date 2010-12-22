#!/usr/bin/env python

"""
Test suite for the Lgm_Vector file

@author: Brian Larsen
@organization: LANL
@contact: balarsen@lanl.gov

@version: V1: 20-Dec-2010 (BAL)
"""


import unittest
import Lgm_CTrans
import Lgm_Vector
import _Lgm
import math


class Lgm_VecTests(unittest.TestCase):
    """
    Tests related to Lgm_Vector
    @author: Brian Larsen
    @organization: LANL
    @contact: balarsen@lanl.gov

    @version: V1: 20-Dec-2010 (BAL)
    """

    def setUp(self):
        super(Lgm_VecTests, self).setUp()
        self.lgm = _Lgm._Lgm()

    def tearDown(self):
        super(Lgm_VecTests, self).tearDown()

    def test_Lgm_CrossProduct(self):
        """Lgm_CrossProduct should give known output"""
        vec1 = Lgm_Vector.Lgm_Vector(1,2,3)
        vec2 = Lgm_Vector.Lgm_Vector(3,1,0)
        vecans = Lgm_Vector.Lgm_Vector()
        self.lgm.lib.Lgm_CrossProduct(vec1, vec2, vecans)
        self.assertEqual(-3.0, vecans.x)
        self.assertEqual( 9.0, vecans.y)
        self.assertEqual(-5.0, vecans.z)

    def test_Lgm_DotProduct(self):
        """Lgm_DotProduct should give known output"""
        vec1 = Lgm_Vector.Lgm_Vector(1,2,3)
        vec2 = Lgm_Vector.Lgm_Vector(3,1,0)
        self.assertEqual(5, self.lgm.lib.Lgm_DotProduct(vec1, vec2))
        vec1 = Lgm_Vector.Lgm_Vector(1,2,3)
        vec2 = Lgm_Vector.Lgm_Vector(1,2,3)
        self.assertEqual(14, self.lgm.lib.Lgm_DotProduct(vec1, vec2))

    def test_Lgm_NormalizeVector(self):
        """Lgm_NormalizeVector should give known output"""
        vec1 = Lgm_Vector.Lgm_Vector(3,1,0)
        magans = self.lgm.lib.Lgm_NormalizeVector(vec1)
        self.assertAlmostEqual(3.1622776601683795, magans)
        self.assertAlmostEqual(0.94868329805051377, vec1.x)
        self.assertAlmostEqual(0.31622776601683794, vec1.y)
        self.assertAlmostEqual(0.0, vec1.z)

    def test_Lgm_ScaleVector(self):
        """Lgm_ScaleVector should give known output"""
        vec1 = Lgm_Vector.Lgm_Vector(3,1,0)
        self.lgm.lib.Lgm_ScaleVector(vec1, 10.5)
        self.assertAlmostEqual(31.5, vec1.x)
        self.assertAlmostEqual(10.5, vec1.y)
        self.assertAlmostEqual(0.0, vec1.z)

    def test_Lgm_Magnitude(self):
        """Lgm_Magnitude should give known output"""
        vec1 = Lgm_Vector.Lgm_Vector(3,1,0)
        magans = self.lgm.lib.Lgm_Magnitude(vec1)
        self.assertAlmostEqual(3.1622776601683795, magans)

    def test_Lgm_Vector(self):
        """Lgm_Vector has a x, y, z"""
        self.assertTrue(hasattr(Lgm_Vector.Lgm_Vector, 'x'))
        self.assertTrue(hasattr(Lgm_Vector.Lgm_Vector, 'y'))
        self.assertTrue(hasattr(Lgm_Vector.Lgm_Vector, 'z'))

    def test_Lgm_Vector_Type(self):
        """Lgm_Vector should be of type LgmDouble"""
        vec1 = Lgm_Vector.Lgm_Vector(3,1,0)
        self.assertTrue(isinstance(vec1.x, float))
        self.assertTrue(isinstance(vec1.y, float))
        self.assertTrue(isinstance(vec1.z, float))
        self.assertRaises(TypeError, Lgm_Vector.Lgm_Vector, 'bad', 0, 1)
        self.assertRaises(TypeError, Lgm_Vector.Lgm_Vector, 1, 'bad', 0)
        self.assertRaises(TypeError, Lgm_Vector.Lgm_Vector, 5, 3, 'bad')


# Tests not written for these yet
#'Lgm_VecSub': [None, Lgm_VectorP, Lgm_VectorP, Lgm_VectorP],
#'Lgm_VecAdd': [None, Lgm_VectorP, Lgm_VectorP, Lgm_VectorP],
#'Lgm_VecDiffMag': [LgmDouble, Lgm_VectorP, Lgm_VectorP],
#'Lgm_ForceMagnitude': [None, Lgm_VectorP, LgmDouble],
#'Lgm_MatTimesVec': [None, LgmDouble * 3 * 3, Lgm_VectorP, Lgm_VectorP],
#'Lgm_Transpose' : [None, LgmDouble * 3 * 3, LgmDouble * 3 * 3],
#'Lgm_MatTimesMat' : [None, LgmDouble * 3 * 3, LgmDouble * 3 * 3, LgmDouble * 3 * 3],


if __name__ == '__main__':
    unittest.main()
