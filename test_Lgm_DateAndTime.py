#!/usr/bin/env python

"""
Test suite for the Lgm_DateAndTime file

@author: Brian Larsen
@organization: LANL
@contact: balarsen@lanl.gov

@version: V1: 22-Dec-2010 (BAL)
"""


import unittest
import Lgm_DateAndTime
import _Lgm
from Lgm_Types import *


class Lgm_DateAndTimeTests(unittest.TestCase):
    """
    Tests related to Lgm_DateAndTime
    @author: Brian Larsen
    @organization: LANL
    @contact: balarsen@lanl.gov

    @version: V1: 20-Dec-2010 (BAL)
    """

    def setUp(self):
        super(Lgm_DateAndTimeTests, self).setUp()
        self.lgm = _Lgm._Lgm()

    def tearDown(self):
        super(Lgm_DateAndTimeTests, self).tearDown()


    def test_Lgm_DateAndTime(self):
        """Lgm_DateAndTime has a nLeapSecondDates, LeapSecondDates, LeapSecondJDs, LeapSeconds"""
        self.assertTrue(hasattr(Lgm_DateAndTime.Lgm_DateAndTime, 'nLeapSecondDates'))
        self.assertTrue(hasattr(Lgm_DateAndTime.Lgm_DateAndTime, 'LeapSecondDates'))
        self.assertTrue(hasattr(Lgm_DateAndTime.Lgm_DateAndTime, 'LeapSecondJDs'))
        self.assertTrue(hasattr(Lgm_DateAndTime.Lgm_DateAndTime, 'LeapSeconds'))


    def test_Lgm_DateAndTime_Type(self):
        """Lgm_DateAndTime input should be of correct types """
        b = LgmLongP()
        c = LgmDoubleP()
        d = LgmDoubleP()
        ans = Lgm_DateAndTime.Lgm_DateAndTime(5, b, c, d)
        self.assertTrue(isinstance(ans.nLeapSecondDates, int))
        self.assertTrue(isinstance(ans.LeapSecondDates, LgmLongP))
        self.assertTrue(isinstance(ans.LeapSecondJDs, LgmDoubleP))
        self.assertTrue(isinstance(ans.LeapSeconds, LgmDoubleP))









if __name__ == '__main__':
    unittest.main()
