#!/usr/bin/env python

import unittest
import ctypes

import numpy

import Lgm_Wrap

class Lgm_WrapTests(unittest.TestCase):
    """
    Tests related to generate file Lgm_Wrap
    @author: Brian Larsen
    @organization: LANL
    @contact: balarsen@lanl.gov

    @version: V1: 28-Feb-2011 (BAL)
    """
    def setUp(self):
        super(Lgm_WrapTests, self).setUp()

    def tearDown(self):
        super(Lgm_WrapTests, self).tearDown()

    def test_MutableString(self):
        """Tests of the builtin MutableString"""
        a = Lgm_Wrap.MutableString('asd')
        self.assertEqual(a, 'asd')
        a[0] = 'x'
        self.assertEqual(a, 'xsd')
        self.assertRaises(TypeError, hash, a)



if __name__ == '__main__':
    unittest.main()
