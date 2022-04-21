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
        ans, north_fp, south_fp, minB, L \
            = Closed_Field([1,2,2], self.date, extended_out=True)
        self.assertEqual(ans, 'LGM_CLOSED')
        numpy.testing.assert_allclose(north_fp,
                [-0.0620421299, 0.39455411901, 0.9306480813], rtol=1e-5)
        numpy.testing.assert_allclose(south_fp,
                [0.7666220541, 0.3352533669, -0.571429965], rtol=1e-5)
        numpy.testing.assert_allclose(minB,
                [2.279625676, 2.809834544, 1.1244815844], rtol=1e-5)
        #This value isn't tested in the C
        #So this is essentially a regression test on _simpleL
        self.assertAlmostEqual(L, 3.7371417478198037, places=5)


if __name__ == '__main__':
    unittest.main()
