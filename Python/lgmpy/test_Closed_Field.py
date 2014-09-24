#!/usr/bin/env python
# -*- coding: utf-8 -*-

import unittest
import datetime

import numpy

from lgmpy import Closed_Field

class Closed_FieldTests(unittest.TestCase):
    """
    Tests related to Closed_FieldTests
    """
    def setUp(self):
        super(Closed_FieldTests, self).setUp()
        self.date = datetime.datetime(2010, 12, 12)

    def tearDown(self):
        super(Closed_FieldTests, self).tearDown()

    def test_input_checking(self):
        """Closed_Field does input params count checking"""
        self.assertRaises(RuntimeError, Closed_Field, '')
        self.assertRaises(RuntimeError, Closed_Field, '', '', '')

    def test_field_model(self):
        """test whether field not implemented raises exception"""
        self.assertRaises(NotImplementedError, Closed_Field,
                          [1,2,3], self.date, bfield='Lgm_B_TS04')

    def test_position(self):
        """position must be an iterable"""
        self.assertRaises(TypeError, Closed_Field,
                          1, self.date)

    def test_Closed_Field_Output(self):
        """Closed_Field should give known results (regression)"""
        self.assertEqual(Closed_Field([1,2,2], self.date), 'LGM_CLOSED')
        self.assertEqual(Closed_Field([.1,.1,.1], self.date), 'LGM_INSIDE_EARTH')
        self.assertEqual(Closed_Field([1,1,10], self.date), 'LGM_OPEN_N_LOBE')
        self.assertEqual(Closed_Field([1,1,9], self.date), 'LGM_OPEN_N_LOBE')
        self.assertEqual(Closed_Field([12,1,9], self.date), 'LGM_OPEN_IMF')
        self.assertEqual(Closed_Field([4,1,-9], self.date), 'LGM_OPEN_S_LOBE')
        self.assertTrue(Closed_Field([.1,1,0], self.date) in ['LGM_BAD_TRACE', 'LGM_INSIDE_EARTH'])
        # still haven't testing this one
        #'LGM_TARGET_HEIGHT_UNREACHABLE'

    def test_extended_out(self):
        """Closed_Field extended_out flag should have known behaviour (regression)"""
        data = Closed_Field([1,2,2], self.date, extended_out = True)
        self.assertEqual(data[0], 'LGM_CLOSED')
        numpy.testing.assert_allclose(data[1],
                [-0.062054,  0.394555,  0.930685], rtol=1e-4)
        numpy.testing.assert_allclose(data[2],
                [ 0.76665,  0.33534, -0.57132], rtol=1e-4)
        numpy.testing.assert_allclose(data[3],
                [ 2.27943481,  2.80962393,  1.12454726], rtol=1e-4)
        self.assertAlmostEqual(data[4], 3.737398946245622, places=5)


if __name__ == '__main__':
    unittest.main()
