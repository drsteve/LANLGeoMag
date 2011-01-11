#!/usr/bin/env python

"""
Test suite for the _Lgm

@author: Brian Larsen
@organization: LANL
@contact: balarsen@lanl.gov

@version: V1: 11-Jan-2011 (BAL)
"""


import unittest
import ctypes
import sys

import Lgm
import Lgm_Vector

class Lgm_Tests(unittest.TestCase):
    def setUp(self):
        super(Lgm_Tests, self).setUp()
    def tearDown(self):
        super(Lgm_Tests, self).tearDown()



if __name__ == '__main__':
    unittest.main()
