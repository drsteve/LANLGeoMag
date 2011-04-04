#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Test suite for the LstarVersusPA

@author: Brian Larsen
@organization: LANL
@contact: balarsen@lanl.gov

@version: V1: 4-Apr-2011 (BAL)
"""

import unittest
import datetime

import numpy

import magcoords

class magcoords_Tests(unittest.TestCase):
    def setUp(self):
        super(magcoords_Tests, self).setUp()

    def tearDown(self):
        super(magcoords_Tests, self).tearDown()

    def test_end2end(self):
        """quick functional tests (regression)"""
        ans = [-3.608026916281572, 4.163336342344337e-17, -1.7268878861662331]
        numpy.testing.assert_array_almost_equal(ans,
                            magcoords.coordTrans([-4,0,0], datetime.datetime(2009,1,1),'SM','GSM'), decimal=6)
        ans = [-3.9999999999999991, 4.0592529337857286e-16, 8.8817841970012523e-16]
        numpy.testing.assert_array_almost_equal(ans,
                            magcoords.coordTrans([-3.608026916281573, 2.5673907444456745e-16,
            -1.7268878861662329], datetime.datetime(2009,1,1),'GSM','SM'), decimal=6)

    def test_input(self):
        """there is some input checking"""
        self.assertRaises(NotImplementedError, magcoords.coordTrans, [-3.608026916281573, 2.5673907444456745e-16,
            -1.7268878861662329], datetime.datetime(2009,1,1),'GSM','BAD')



if __name__ == '__main__':
    unittest.main()
