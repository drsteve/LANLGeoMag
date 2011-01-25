#!/usr/bin/env python

"""
Test suite for the DataArray file

@author: Brian Larsen
@organization: LANL
@contact: balarsen@lanl.gov

@version: V1: 21-Jan-2010 (BAL)
"""

import unittest
import datetime
import itertools
import ctypes

import numpy

import Lgm_Data

class DataArrayTests(unittest.TestCase):
    """
    Tests related to DataArray class
    @author: Brian Larsen
    @organization: LANL
    @contact: balarsen@lanl.gov

    @version: V1: 21-Jan-2010 (BAL)
    """
    def setUp(self):
        super(DataArrayTests, self).setUp()
    def tearDown(self):
        super(DataArrayTests, self).tearDown()

    def test_Lgm_Data_repr__(self):
        """__repr__ is abstract and should raise NotImplementedError"""
        a = Lgm_Data.Lgm_Data([])
        self.assertRaises(NotImplementedError, a.__repr__)

    def test_Lgm_Array(self):
        """Lgm_Data should behave"""
        a = Lgm_Data.DataArray([1,2,3], attrs={'bla':'bla'})
        self.assertTrue(hasattr(a, 'attrs'))
        self.assertTrue(a.attrs['bla'], 'bla')


if __name__ == '__main__':
    unittest.main()
