#!/usr/bin/env python

"""
Test suite for the Lgm_LeapSeconds file

@author: Brian Larsen
@organization: LANL
@contact: balarsen@lanl.gov

@version: V1: 22-Dec-2010 (BAL)
"""


import unittest
import Lgm_LeapSeconds
import _Lgm
from Lgm_Types import *


class Lgm_LeapSecondsTests(unittest.TestCase):
    """
    Tests related to Lgm_LeapSeconds
    @author: Brian Larsen
    @organization: LANL
    @contact: balarsen@lanl.gov

    @version: V1: 20-Dec-2010 (BAL)
    """

    def setUp(self):
        super(Lgm_LeapSecondsTests, self).setUp()
        self.lgm = _Lgm._Lgm()

    def tearDown(self):
        super(Lgm_LeapSecondsTests, self).tearDown()


    def test_Lgm_LeapSeconds(self):
        """Lgm_LeapSeconds has a nLeapSecondDates, LeapSecondDates, LeapSecondJDs, LeapSeconds"""
        self.assertTrue(hasattr(Lgm_LeapSeconds.Lgm_LeapSeconds, 'nLeapSecondDates'))
        self.assertTrue(hasattr(Lgm_LeapSeconds.Lgm_LeapSeconds, 'LeapSecondDates'))
        self.assertTrue(hasattr(Lgm_LeapSeconds.Lgm_LeapSeconds, 'LeapSecondJDs'))
        self.assertTrue(hasattr(Lgm_LeapSeconds.Lgm_LeapSeconds, 'LeapSeconds'))


    def test_Lgm_LeapSeconds_Type(self):
        """Lgm_LeapSeconds input should be of correct types """
        b = LgmLongP()
        c = LgmDoubleP()
        d = LgmDoubleP()
        ans = Lgm_LeapSeconds.Lgm_LeapSeconds(5, b, c, d)
        self.assertTrue(isinstance(ans.nLeapSecondDates, int))
        self.assertTrue(isinstance(ans.LeapSecondDates, LgmLongP))
        self.assertTrue(isinstance(ans.LeapSecondJDs, LgmDoubleP))
        self.assertTrue(isinstance(ans.LeapSeconds, LgmDoubleP))









if __name__ == '__main__':
    unittest.main()
