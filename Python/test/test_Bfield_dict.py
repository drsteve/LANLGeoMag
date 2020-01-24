#!/usr/bin/env python

import unittest
import ctypes
import numpy as np

from lgmpy import _Bfield_dict

class Bfield_dictTests(unittest.TestCase):
    """
    Regression catching for _Bfield_dict
    """
    def setUp(self):
        super(Bfield_dictTests, self).setUp()

    def tearDown(self):
        super(Bfield_dictTests, self).tearDown()

    def test_length(self):
        """the _Bfield_dict dist should have known length"""
        self.assertEqual(len(_Bfield_dict.Bfield_dict), 8)

    def test_keys(self):
        """the _Bfield_dict dist should have known keys"""
        keys = ['Lgm_B_OP77', 'Lgm_B_edip', 'Lgm_B_T89', 'Lgm_B_T89c', 'Lgm_B_cdip', 'Lgm_B_T96', 'Lgm_B_Dungey', 'Lgm_B_T01S']
        self.assertEqual(
            sorted(keys), sorted(_Bfield_dict.Bfield_dict.keys()))


if __name__ == '__main__':
    unittest.main()
