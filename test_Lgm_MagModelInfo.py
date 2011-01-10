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
import math

import Lgm
import Lgm_MagModelInfo
from _Lgm import lib

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



if __name__ == '__main__':
    unittest.main()
