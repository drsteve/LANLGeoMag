#!/usr/bin/env python2.6

import unittest
import Lgm_Vector


class Lgm_VecTests(unittest.TestCase):
    """Tests related to Lgm_Vector"""

    def setUp(self):
        super(Lgm_VecTests, self).setUp()
        self.lgm = Lgm_Vec.Lgm_Vec()
        self.vec1 = Lgm_Vec.Lgm_Vector_cls(1,2,3)
        self.vec2 = Lgm_Vec.Lgm_Vector_cls(3,1,0)
        self.vec3 = Lgm_Vec.Lgm_Vector_cls(7,8,6)
        self.vec4 = Lgm_Vec.Lgm_Vector_cls(8,3,9)
        self.vec5 = Lgm_Vec.Lgm_Vector_cls(2,5,1.4)


    def tearDown(self):
        super(Lgm_VecTests, self).tearDown()

    def test_Lgm_CrossProduct(self):
        """Lgm_CrossProduct should give known output"""
        vecans = Lgm_Vec.Lgm_Vector_cls()
        self.lgm.lib.Lgm_CrossProduct(self.vec1, self.vec2, vecans)
        self.assertEqual(-3.0, vecans.x)
        self.assertEqual( 9.0, vecans.y)
        self.assertEqual(-5.0, vecans.z)





if __name__ == '__main__':
    unittest.main()

