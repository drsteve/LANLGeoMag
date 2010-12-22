#!/usr/bin/env python

"""
Test suite for the Lgm_LeapSeconds file

@author: Brian Larsen
@organization: LANL
@contact: balarsen@lanl.gov

@version: V1: 22-Dec-2010 (BAL)
"""


import unittest
import Lgm_CTrans
import Lgm_Vector
import Lgm_LeapSeconds
import _Lgm
import math


class Lgm_LeapSecondsTests(unittest.TestCase):
    """
    Tests related to Lgm_LeapSeconds
    @author: Brian Larsen
    @organization: LANL
    @contact: balarsen@lanl.gov

    @version: V1: 20-Dec-2010 (BAL)
    """

    def setUp(self):
        super(Lgm_LeapSecondsTests, self).setUp()
        self.lgm = _Lgm._Lgm()

    def tearDown(self):
        super(Lgm_LeapSecondsTests, self).tearDown()








if __name__ == '__main__':
    unittest.main()
