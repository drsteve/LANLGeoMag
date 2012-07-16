#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Test suite for the Lgm_Utils

@author: Brian Larsen
"""

import unittest
import ctypes
import numpy as np

from Lgm_Wrap import Lgm_LogSpace, Lgm_Bisect

class Lgm_Utils_Tests(unittest.TestCase):
    def test_Lgm_LogSpace(self):
        """regression on Lgm_LogSpace"""
        a = np.empty([10], dtype=ctypes.c_double)
        start = 1.
        stop = 100.
        Lgm_LogSpace(start, stop, len(a), a.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))
        expected = [  1.        ,   1.58489319,   2.51188643,   3.98107171,
         6.30957344,  10.        ,  15.84893192,  25.11886432,
        39.81071706,  63.09573445]
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

