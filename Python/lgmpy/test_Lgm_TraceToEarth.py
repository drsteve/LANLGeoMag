#!/usr/bin/env python

"""
Test suite for the Lgm_T89 file

@author: Brian Larsen
@organization: LANL
@contact: balarsen@lanl.gov

@version: V1: 20-Dec-2010 (BAL)
"""

import unittest
import datetime
import itertools

import numpy

from lgmpy import Lgm_TraceToEarth
from lgmpy import Lgm_Vector


class Lgm_TraceToEarthTests(unittest.TestCase):
    """
    Tests related to Lgm_TraceToEarth
    @author: Brian Larsen
    @organization: LANL
    @contact: balarsen@lanl.gov

    """
    def setUp(self):
        super(Lgm_TraceToEarthTests, self).setUp()
        self.pos = [-6.6, 0, 0]
        self.kp = 2
        self.dt = datetime.datetime(2005, 8, 31, 9, 0, 0)
    def tearDown(self):
        super(Lgm_TraceToEarthTests, self).tearDown()

    def test_TraceToEarth(self):
        """TraceToEarth should given known result"""
        ans = ([-0.389267, 0.004424, 0.952456],
                [-0.533227, 0.097428, -0.874284])
        numpy.testing.assert_allclose(ans,
                                      Lgm_TraceToEarth.TraceToEarth(self.pos,
                                                                    self.dt,
                                                                    MAGMODEL_ARGS={'Kp':self.kp},
                                                                    TargetHeight=200), rtol=1e-04, atol=0)
        


if __name__ == '__main__':
    unittest.main()
