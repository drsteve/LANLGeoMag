#!/usr/bin/env python

"""
Test suite for the Lgm_LeapSeconds file

@author: Brian Larsen
@organization: LANL
@contact: balarsen@lanl.gov

@version: V1: 10-Jan-2011 (BAL)
"""

import ctypes
import unittest

import Lgm
import Lgm_LeapSeconds
from _Lgm import lib

class Lgm_LeapSecondsTests(unittest.TestCase):
    """
    Tests related to Lgm_Octree
    @author: Brian Larsen
    @organization: LANL
    @contact: balarsen@lanl.gov

    @version: V1: 10-Jan-2011 (BAL)
    """

    def setUp(self):
        super(Lgm_LeapSecondsTests, self).setUp()

    def tearDown(self):
        super(Lgm_LeapSecondsTests, self).tearDown()

    def test_Lgm_LeapSeconds_size(self):
        """Lgm_LeapSeconds c/python must be the same size"""
        self.assertEqual(lib.size_Lgm_LeapSeconds(),
                         ctypes.sizeof(Lgm_LeapSeconds.Lgm_LeapSeconds) ) 





if __name__ == '__main__':
    unittest.main()
