#!/usr/bin/env python

import unittest
import ctypes
import numpy as np

import _Bfield_dict

class Bfield_dictTests(unittest.TestCase):
    """
    Regression catching for _Bfield_dict
    @author: Brian Larsen
    @organization: LANL
    @contact: balarsen@lanl.gov

    @version: V1: 21-Apr-2011 (BAL)
    """
    def setUp(self):
        super(Bfield_dictTests, self).setUp()

    def tearDown(self):
        super(Bfield_dictTests, self).tearDown()

    def test_length(self):
        """the _Bfield_dict dist should have known length"""
        self.assertEqual(len(_Bfield_dict.Bfield_dict), 4)

    def test_keys(self):
        """the _Bfield_dict dist should have known keys"""
        keys = ['Lgm_B_OP77', 'Lgm_B_edip', 'Lgm_B_T89', 'Lgm_B_cdip']
        for val in keys:
            self.assertTrue(val in _Bfield_dict.Bfield_dict)


if __name__ == '__main__':
    unittest.main()
