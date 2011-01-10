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

    def test_gslStructSizes(self):
        """The c and pyhtion sizes of the used GSL structs much be the same"""
        self.assertEqual(lib.size_gsl_interp_accel(), ctypes.sizeof(Lgm_MagModelInfo.gsl_interp_accel))
        self.assertEqual(lib.size_gsl_interp_type(), ctypes.sizeof(Lgm_MagModelInfo.gsl_interp_type))
        self.assertEqual(lib.size_gsl_interp(), ctypes.sizeof(Lgm_MagModelInfo.gsl_interp))
        self.assertEqual(lib.size_gsl_spline(), ctypes.sizeof(Lgm_MagModelInfo.gsl_spline))

    def test_size_Lgm_MagModelInfo(self):
        """for Lgm_MagModelInfo the c and python sizes should match"""
        print("c: %d" % (lib.size_MagModelInfo()) )
        self.assertEqual(lib.size_MagModelInfo(),
                         ctypes.sizeof(Lgm_MagModelInfo.Lgm_MagModelInfo) )


if __name__ == '__main__':
    unittest.main()
