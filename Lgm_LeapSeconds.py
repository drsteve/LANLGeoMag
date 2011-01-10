
import ctypes

import Lgm
from _Lgm import lib
from Lgm_Types import LgmInt, LgmLongP, LgmDoubleP


class Lgm_LeapSeconds(ctypes.Structure):
    """
    LanlGeoMag class for holding leap second info

    @author: Brian Larsen
    @organization: LANL
    @contact: balarsen@lanl.gov

    @version: V1: 10-Jan-2011 (BAL)
    """
    @classmethod
    def assign_fields(cls):
        cls._fields_ = [ \
        ("nLeapSecondDates", LgmInt),
        ("LeapSecondDates", LgmLongP),
        ("LeapSecondJDs", LgmDoubleP),
        ("LeapSeconds", LgmDoubleP) ]
