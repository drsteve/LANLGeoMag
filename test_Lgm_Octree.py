#!/usr/bin/env python

"""
Test suite for the Lgm_Octree file

@author: Brian Larsen
@organization: LANL
@contact: balarsen@lanl.gov

@version: V1: 20-Dec-2010 (BAL)
"""

import ctypes
import unittest

import Lgm
import _Lgm_Octree
from _Lgm import lib

class Lgm_OctreeTests(unittest.TestCase):
    """
    Tests related to Lgm_Octree
    @author: Brian Larsen
    @organization: LANL
    @contact: balarsen@lanl.gov

    @version: V1: 10-Jan-2011 (BAL)
    """

    def setUp(self):
        super(Lgm_OctreeTests, self).setUp()

    def tearDown(self):
        super(Lgm_OctreeTests, self).tearDown()

    def test_Lgm_OctreeDataSize(self):
        """the Lgm_OctreeData struct must have same C/Python size"""
        self.assertEqual(lib.size_Lgm_OctreeData(),
                         ctypes.sizeof(_Lgm_Octree.Lgm_OctreeData) )

    def test_Lgm_OctreeCellSize(self):
        """the Lgm_OctreeCell struct must have same C/Python size"""
        self.assertEqual(lib.size_Lgm_OctreeCell(),
                         ctypes.sizeof(_Lgm_Octree.Lgm_OctreeCell) )

    def test_pQueueSize(self):
        """the pQueue struct must have same C/Python size"""
        self.assertEqual(lib.size_pQueue(),
                         ctypes.sizeof(_Lgm_Octree.pQueue) )



if __name__ == '__main__':
    unittest.main()
