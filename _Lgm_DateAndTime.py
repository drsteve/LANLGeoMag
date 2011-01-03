"""
Lgm_DateAndTime module, this contains the necessary code Leap Seconds


@author: Brian Larsen
@organization: LANL
@contact: balarsen@lanl.gov

@version: V1: 20-Dec-2010 (BAL)
"""


import ctypes
from Lgm_Types import LgmLong, LgmDouble, LgmInt



class Lgm_DateAndTime(ctypes.Structure):
    @classmethod
    def assign_fields(cls):
        cls._fields_ = [ ("nLeapSecondDates", LgmInt), # Number of leap second dates.
            ("LeapSecondDates", LgmLong), # Array for holdin the Dates on which leap
                                      # seconds were added
            ("LeapSecondJDs", LgmDouble), # Array for holdin the Julian Dates on
                                       # which leap seconds were added
            ("LeapSeconds", LgmDouble) ] # The actual number of leap seconds that
                                         # went into effect on the given date
