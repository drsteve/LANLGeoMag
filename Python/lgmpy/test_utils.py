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

from lgmpy import utils
from lgmpy import Lgm_Vector

class Utils_Tests(unittest.TestCase):
    """
    Tests related to Lgm.utils wrapper
    """
    def test_pos2Lgm_Vector(self):
        """pos2Lgm_Vector should have known output"""
        a = [1.2, 3.4, 4.3]
        self.assertEqual(Lgm_Vector.Lgm_Vector(*a), utils.pos2Lgm_Vector(a))
        self.assertEqual(Lgm_Vector.Lgm_Vector(*a), utils.pos2Lgm_Vector(Lgm_Vector.Lgm_Vector(*a)))
        self.assertEqual(Lgm_Vector.Lgm_Vector(*a), utils.pos2Lgm_Vector(np.asarray(a)))
        self.assertEqual([Lgm_Vector.Lgm_Vector(*a)]*2, utils.pos2Lgm_Vector([a]*2))
        
                
        


if __name__ == '__main__':
    unittest.main()

