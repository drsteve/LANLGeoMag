#!/usr/bin/env python

import unittest
import datetime

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

    def test_coord_system(self):
        """Only GSM is implemented thus far"""
        self.assertRaises(NotImplementedError, Closed_Field.Closed_Field,
                          [1,2,3], self.date, coord_system='SM')

    def test_field_model(self):
        """only Lgm_B_OP77 is implemented now"""
        self.assertRaises(NotImplementedError, Closed_Field.Closed_Field,
                          [1,2,3], self.date, bfield='Lgm_B_T89')

    def test_position(self):
        """position must be an iterable"""
        self.assertRaises(TypeError, Closed_Field.Closed_Field,
                          1, self.date)

    def test_Closed_Field_Ouptut(self):
        """Closed_Field should give known results (regression)"""
        self.assertEqual(Closed_Field.Closed_Field([1,2,2], self.date), 'LGM_CLOSED')
        self.assertEqual(Closed_Field.Closed_Field([1,0,0], self.date), 'LGM_INSIDE_EARTH')
        self.assertEqual(Closed_Field.Closed_Field([1,1,10], self.date), 'LGM_OPEN_N_LOBE')
        self.assertEqual(Closed_Field.Closed_Field([1,1,9], self.date), 'LGM_CLOSED')
        # need to find some tests that hit the rest of the options
        #'LGM_OPEN_IMF'
        #'LGM_CLOSED'
        #'LGM_OPEN_N_LOBE'
        #'LGM_OPEN_S_LOBE'
        #'LGM_INSIDE_EARTH'
        #'LGM_TARGET_HEIGHT_UNREACHABLE'


if __name__ == '__main__':
    unittest.main()


