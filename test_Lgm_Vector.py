#!/usr/bin/env python

import unittest
import ctypes

import numpy

import Lgm_Vector

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

    def tearDown(self):
        super(Lgm_VecTests, self).tearDown()

    def testLgm_Vector(self):
        """Lgm_Vector has a x, y, z"""
        self.assertTrue(hasattr(Lgm_Vector.Lgm_Vector, 'x'))
        self.assertTrue(hasattr(Lgm_Vector.Lgm_Vector, 'y'))
        self.assertTrue(hasattr(Lgm_Vector.Lgm_Vector, 'z'))

    def testLgm_Vector_Type(self):
        """Lgm_Vector should be of type LgmDouble"""
        vec1 = Lgm_Vector.Lgm_Vector(3,1,0)
        self.assertTrue(isinstance(vec1.x, float))
        self.assertTrue(isinstance(vec1.y, float))
        self.assertTrue(isinstance(vec1.z, float))
        self.assertRaises(TypeError, Lgm_Vector.Lgm_Vector, 'bad', 0, 1)
        self.assertRaises(TypeError, Lgm_Vector.Lgm_Vector, 1, 'bad', 0)
        self.assertRaises(TypeError, Lgm_Vector.Lgm_Vector, 5, 3, 'bad')

    def test_Vectorsize(self):
        """Make sure that the Lgm c and python are the same size"""
        lib = __import__('Lgm_Wrap', fromlist=['size_Vector'])
        self.assertEqual(ctypes.sizeof(Lgm_Vector.Lgm_Vector), lib.size_Vector())
        del lib

class Lgm_VectorTestsWrap(unittest.TestCase):
    """
    Tests related to Lgm_Vector
    @author: Brian Larsen
    @organization: LANL
    @contact: balarsen@lanl.gov

    @version: V1: 22-Dec-2010 (BAL)
    """

    def setUp(self):
        super(Lgm_VectorTestsWrap, self).setUp()
        self.vec1 = Lgm_Vector.Lgm_Vector(1, 2, 3)
        self.vec2 = Lgm_Vector.Lgm_Vector(3, 4, 5)

    def tearDown(self):
        super(Lgm_VectorTestsWrap, self).tearDown()

    def test_eq(self):
        """__eq__ has known output"""
        self.assertTrue(self.vec1 == self.vec1)
        self.assertFalse(self.vec1 == self.vec2)
        self.vec2.x = self.vec1.x
        self.assertFalse(self.vec1 == self.vec2)
        self.vec2.y = self.vec1.y
        self.assertFalse(self.vec1 == self.vec2)
        self.assertTrue(self.vec1 == [1, 2, 3])
        self.assertFalse(self.vec1 == [3, 2, 3])
        self.assertRaises(TypeError, self.vec1.__eq__, ['bad'])
        self.assertRaises(TypeError, self.vec1.__eq__, 'bad')

    def test_gt(self):
        """__gt__ has known output"""
        self.assertTrue(self.vec2 > self.vec1)
        self.assertFalse(self.vec1 > self.vec2)
        self.assertFalse(self.vec1 > self.vec1)

    def test_ge(self):
        """__ge__ has known output"""
        self.assertTrue(self.vec2 >= self.vec1)
        self.assertFalse(self.vec1 >= self.vec2)
        self.assertTrue(self.vec1 >= self.vec1)

    def test_lt(self):
        """__lt__ has known output"""
        self.assertTrue(self.vec1 < self.vec2)
        self.assertFalse(self.vec2 < self.vec1)
        self.assertFalse(self.vec1 < self.vec1)

    def test_le(self):
        """__le__ has known output"""
        self.assertTrue(self.vec1 <= self.vec2)
        self.assertFalse(self.vec2 <= self.vec1)
        self.assertTrue(self.vec1 <= self.vec1)

    def test_str(self):
        """str has known output"""
        self.vec1 = Lgm_Vector.Lgm_Vector()
        self.assertEqual(str(self.vec1), '[0.0, 0.0, 0.0]')

    def test_repr_(self):
        """__repr__ has a known output"""
        val = self.vec1.__repr__()
        self.assertEqual('[1.0, 2.0, 3.0]', val)

    def test_getattr(self):
        """__getattr__ has 1 element and 1 exception"""
        val = self.vec1.mag
        self.assertAlmostEqual(3.7416573867739413, self.vec1.mag)

    def test_getattrException(self):
        """only mag is implemented"""
        self.assertRaises(AttributeError, self.vec1.__getattribute__, 'bad')
        self.assertRaises(AttributeError, self.vec1.__getattr__, 'bad')

    def test_add(self):
        """add gives known output"""
        vec3 = self.vec1 + self.vec2
        self.assertEqual(4, vec3.x)
        self.assertEqual(6, vec3.y)
        self.assertEqual(8, vec3.z)

        vec3 = self.vec1 + 1.5
        self.assertEqual(2.5, vec3.x)
        self.assertEqual(3.5, vec3.y)
        self.assertEqual(4.5, vec3.z)

        self.assertRaises(ArithmeticError, self.vec1.__add__, 'bad' )

    def test_sub(self):
        """sub gives known output"""
        vec3 = self.vec1 - self.vec2
        self.assertEqual(-2, vec3.x)
        self.assertEqual(-2, vec3.y)
        self.assertEqual(-2, vec3.z)

        vec3 = self.vec1 - 1.5
        self.assertEqual(-0.5, vec3.x)
        self.assertEqual(0.5, vec3.y)
        self.assertEqual(1.5, vec3.z)

        self.assertRaises(ArithmeticError, self.vec1.__sub__, 'bad' )

    def test_mul(self):
        """mul gives known output"""
        ans = Lgm_Vector.Lgm_Vector()
        ans = self.vec1 * self.vec2
        self.assertEqual(-2.0, ans.x)
        self.assertEqual( 4.0, ans.y)
        self.assertEqual(-2.0, ans.z)

        vec3 = self.vec1 * 3
        self.assertEqual(3, vec3.x)
        self.assertEqual(6, vec3.y)
        self.assertEqual(9, vec3.z)

        self.assertRaises(ArithmeticError, self.vec1.__mul__, 'bad' )

    def test_div(self):
        """div gives known output"""
        ans = self.vec1 / 10.
        self.assertEqual(0.1, ans.x)
        self.assertEqual(0.2, ans.y)
        self.assertEqual(0.3, ans.z)

        self.assertRaises(ArithmeticError, self.vec1.__div__, 'bad' )

    def test_crossProduct(self):
        """crossProduct gives known output"""
        ans = self.vec1.crossProduct(self.vec2)
        self.assertEqual(-2.0, ans.x)
        self.assertEqual( 4.0, ans.y)
        self.assertEqual(-2.0, ans.z)
        ans = self.vec1.crossProduct(self.vec1)
        self.assertEqual(0.0, ans.x)
        self.assertEqual(0.0, ans.y)
        self.assertEqual(0.0, ans.z)

    def test_dotProduct(self):
        """dotProduct gives known output"""
        self.assertEqual(26.0, self.vec1.dotProduct(self.vec2))
        self.assertEqual(14.0, self.vec1.dotProduct(self.vec1))
        vec3 = Lgm_Vector.Lgm_Vector(1, 0, 0)
        vec4 = Lgm_Vector.Lgm_Vector(0, 1, 0)
        self.assertEqual(0.0, vec3.dotProduct(vec4))

    def test_magnitude(self):
        """magnitude gives known output"""
        self.assertAlmostEqual(3.7416573867739413, self.vec1.magnitude())

    def test_normalize(self):
        """Normalize should have known behaviour"""
        self.vec2.normalize()
        ans = [0.4242640687119285, 0.565685424949238, 0.7071067811865475]
        self.assertEqual(ans[0], self.vec2.x)
        self.assertEqual(ans[1], self.vec2.y)
        self.assertEqual(ans[2], self.vec2.z)

    def test_scale(self):
        """scale should have known output"""
        self.vec2.scale(2)
        ans = [6, 8, 10]
        self.assertEqual(ans[0], self.vec2.x)
        self.assertEqual(ans[1], self.vec2.y)
        self.assertEqual(ans[2], self.vec2.z)

    def test_diffMag(self):
        """diffMag should have known output (regression)"""
        ans = 3.4641016151377544
        self.assertAlmostEqual(ans, self.vec1.diffMag(self.vec2))

    def test_forceMagnitude(self):
        """forceMagnitude should have known output (regression)"""
        self.vec2.forceMagnitude(10)
        ans = [4.242640687119285, 5.65685424949238, 7.071067811865475]
        self.assertEqual(ans[0], self.vec2.x)
        self.assertEqual(ans[1], self.vec2.y)
        self.assertEqual(ans[2], self.vec2.z)

    def test_SphToCart(self):
        """SphToCart should known known output (regression)"""
        invals = [ [0, 0, 5],
            [90, 0, 5],
            [90, 90, 5],
            [90, 180, 5] ]
        ans = [ [5.0, 0.0, 0.0],
            [3.061616997868383e-16, 0.0, 5.0],
            [1.8746997283273223e-32, 3.061616997868383e-16, 5.0],
            [-3.061616997868383e-16, 3.7493994566546446e-32, 5.0] ]
        for i, val in enumerate(invals):
            vec3 = Lgm_Vector.SphToCart(*val)
            for j, val2 in enumerate(vec3.tolist()):
                self.assertAlmostEqual(ans[i][j], val2)
        # test an input check
        self.assertRaises(ValueError, Lgm_Vector.SphToCart, [1]*2, [2]*3, [3]*2)
        # test putting in lists
        ans_tst = Lgm_Vector.SphToCart(zip(*invals)[0], zip(*invals)[1], zip(*invals)[2] )
        for i, v1 in enumerate(ans_tst):
            for j, v2 in enumerate(v1.tolist()):
                self.assertAlmostEqual(ans[i][j], ans_tst[i].tolist()[j])

if __name__ == '__main__':
    unittest.main()
