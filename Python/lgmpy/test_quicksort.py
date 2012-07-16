#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Test suite for the quicksort

@author: Brian Larsen
"""

import unittest
import ctypes
import numpy as np

from Lgm_Wrap import quicksort, bubbleSort, quicksort_uli

class quicksort_Tests(unittest.TestCase):
    def test_quicksort(self):
        """regression on quicksort"""
        a = np.asarray([-999, 6,4,8,2], dtype=ctypes.c_double) # routine starts at index 1
        expected = a.copy()
        expected[1:] = np.sort(expected[1:])
        quicksort(len(a)-1, a.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))
        np.testing.assert_array_almost_equal(expected, a)
        
    def test_quicksort_uli(self):
        """regression on quicksort_uli"""
        a = np.asarray([999, 6,4,8,2], dtype=ctypes.c_ulong) # routine starts at index 1
        expected = a.copy()
        expected[1:] = np.sort(expected[1:])
        quicksort_uli(len(a)-1, a.ctypes.data_as(ctypes.POINTER(ctypes.c_ulong)))
        np.testing.assert_array_almost_equal(expected, a)

    def test_bubbleSort(self):
        """regression on bubbleSort"""
        a = np.asarray([-999, 6,4,8,2], dtype=ctypes.c_double) # routine starts at index 1
        expected = a.copy()
        expected[1:] = np.sort(expected[1:])
        bubbleSort(len(a), a.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))
        np.testing.assert_array_almost_equal(expected, a)
             

if __name__ == '__main__':
    unittest.main()

