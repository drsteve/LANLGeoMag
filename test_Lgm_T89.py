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

import Lgm
from _Lgm import lib
import Lgm_Vector
import Lgm_T89


class Lgm_T89Tests(unittest.TestCase):
    """
    Tests related to Lgm_T89
    @author: Brian Larsen
    @organization: LANL
    @contact: balarsen@lanl.gov

    @version: V1: 04-Jan-2011 (BAL)
    """

    def setUp(self):
        super(Lgm_T89Tests, self).setUp()

    def tearDown(self):
        super(Lgm_T89Tests, self).tearDown()

    def test_T89_1(self):
        """First simple in/out tests of T89 (regression)"""
        dt = datetime.datetime(2005, 8, 31, 9, 0, 0)
        pos = [-6.6, 0, 0]
        ans = [[-18.97193562594128, -1.8611995170538265, 80.3933831714847],
            [-20.828853435515278, -1.8611995170538265, 74.62222419125123 ],
            [-22.558420054640656, -1.8611995170538265, 70.52457956754591],
            [-26.412234393652312, -1.8611995170538265, 64.64506509634843],
            [-32.16112573540407, -1.8611995170538265, 60.078415300152216],
            [-45.379156657247805, -1.8611995170538265, 49.36315537639906] ]
        for i, kp in enumerate(range(6)):
            a = Lgm_T89.Lgm_T89(pos, dt, kp)
            B = a.calc_B()
            self.assertAlmostEqual(ans[i][0], B.x)
            self.assertAlmostEqual(ans[i][1], B.y)
            self.assertAlmostEqual(ans[i][2], B.z)


if __name__ == '__main__':
    unittest.main()
