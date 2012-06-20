#!/usr/bin/env python
# -*- coding: utf-8 -*-

import unittest
import datetime

import numpy

import Closed_Field

class Closed_FieldTests(unittest.TestCase):
    """
    Tests related to Closed_FieldTests
    @author: Brian Larsen
    @organization: LANL
    @contact: balarsen@lanl.gov

    @version: V1: 1-Mar-2011 (BAL)
    """
    def setUp(self):
        super(Closed_FieldTests, self).setUp()
        self.date = datetime.datetime(2010, 12, 12)

    def tearDown(self):
        super(Closed_FieldTests, self).tearDown()

    def test_input_checking(self):
        """Closed_Field does input params count checking"""
        self.assertRaises(RuntimeError, Closed_Field.Closed_Field, '')
        self.assertRaises(RuntimeError, Closed_Field.Closed_Field, '', '', '')

    def test_coord_system(self):
        """Only GSM is implemented thus far"""
        self.assertRaises(NotImplementedError, Closed_Field.Closed_Field,
                          [1,2,3], self.date, coord_system='SM')

    def test_field_model(self):
        """only Lgm_B_OP77 and Lgm_B_T89 is implemented now"""
        self.assertRaises(NotImplementedError, Closed_Field.Closed_Field,
                          [1,2,3], self.date, bfield='Lgm_B_T96')

    def test_position(self):
        """position must be an iterable"""
        self.assertRaises(TypeError, Closed_Field.Closed_Field,
                          1, self.date)

    def test_Closed_Field_Output(self):
        """Closed_Field should give known results (regression)"""
        self.assertEqual(Closed_Field.Closed_Field([1,2,2], self.date), 'LGM_CLOSED')
        self.assertEqual(Closed_Field.Closed_Field([.1,.1,.1], self.date), 'LGM_INSIDE_EARTH')
        self.assertEqual(Closed_Field.Closed_Field([1,1,10], self.date), 'LGM_OPEN_N_LOBE')
        self.assertEqual(Closed_Field.Closed_Field([1,1,9], self.date), 'LGM_OPEN_N_LOBE')
        self.assertEqual(Closed_Field.Closed_Field([12,1,9], self.date), 'LGM_OPEN_IMF')
        self.assertEqual(Closed_Field.Closed_Field([4,1,-9], self.date), 'LGM_OPEN_S_LOBE')
        self.assertTrue(Closed_Field.Closed_Field([.1,1,0], self.date) in ['LGM_BAD_TRACE', 'LGM_INSIDE_EARTH'])
        # still haven't testing this one
        #'LGM_TARGET_HEIGHT_UNREACHABLE'

    def test_extended_out(self):
        """Closed_Field extended_out flag should have known behaviour (regression)"""
        data = Closed_Field.Closed_Field([1,2,2], self.date, extended_out = True)
        self.assertEqual(data[0], 'LGM_CLOSED')
        numpy.testing.assert_allclose(data[1],
                [-0.0621066,  0.39457349,  0.93069306], rtol=1e-4)
        numpy.testing.assert_allclose(data[2],
                [ 0.76665539,  0.33534993, -0.57132811], rtol=1e-4)
        numpy.testing.assert_allclose(data[3],
                [ 2.27943481,  2.80962393,  1.12454726], rtol=1e-4)
        self.assertAlmostEqual(data[4], 3.7377517203450235, places=5)


if __name__ == '__main__':
    unittest.main()
