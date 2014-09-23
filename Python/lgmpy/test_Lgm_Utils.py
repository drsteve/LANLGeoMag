#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Test suite for the Lgm_Utils

@author: Brian Larsen
"""

import unittest
import ctypes
import numpy as np

from lgmpy.Lgm_Wrap import Lgm_LogSpace, Lgm_Bisect

class Lgm_Utils_Tests(unittest.TestCase):
    def test_Lgm_LogSpace(self):
        """regression on Lgm_LogSpace"""
        a = np.empty([10], dtype=ctypes.c_double)
        start = 1.
        stop = 100.
        Lgm_LogSpace(start, stop, len(a), a.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))
        expected = [   1.        ,    1.66810054,    2.7825594 ,    4.64158883,
                       7.74263683,   12.91549665,   21.5443469 ,   35.93813664,
                       59.94842503,  100.        ]
        np.testing.assert_array_almost_equal(expected, a)

    def test_Lgm_Bisect(self):
        """regression on Lgm_Bisect"""
        a = np.linspace(0, 100, 34)
        fn = lambda v: Lgm_Bisect(a.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), v, len(a))
        self.assertEqual( fn(15.), 5  )
        self.assertEqual( fn(200.), 34  )
        self.assertEqual( fn(-5.), 0  )
        self.assertEqual( fn(50.), 17  )
        self.assertEqual( fn(34.), 12  )
        


if __name__ == '__main__':
    unittest.main()

