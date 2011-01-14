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

import Lgm_MagModelInfo
from Lgm_Wrap import size_MagModelInfo

class Lgm_MagModelInfoTests(unittest.TestCase):
    """
    Tests related to Lgm_MagModelInfo
    @author: Brian Larsen
    @organization: LANL
    @contact: balarsen@lanl.gov

    @version: V1: 10-Jan-2011 (BAL)
    """

    def setUp(self):
        super(Lgm_MagModelInfoTests, self).setUp()

    def tearDown(self):
        super(Lgm_MagModelInfoTests, self).tearDown()

    def test_size_Lgm_MagModelInfo(self):
        """for Lgm_MagModelInfo the c and python sizes should match"""
        self.assertEqual(size_MagModelInfo(),
                         ctypes.sizeof(Lgm_MagModelInfo.Lgm_MagModelInfo) )




if __name__ == '__main__':
    unittest.main()
