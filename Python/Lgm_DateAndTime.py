"""
Lgm_DateAndTime module, this contains the necessary code Leap Seconds


@author: Brian Larsen
@organization: LANL
@contact: balarsen@lanl.gov

@version: V1: 20-Dec-2010 (BAL)
"""


import ctypes
import datetime

import numpy

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

def dateToDateLong(inval):
    """
    convert a python date or datetime object to a Date (long) object that
    LanlGeoMag Likes to use

    @author: Brian Larsen
    @organization: LANL
    @contact: balarsen@lanl.gov

    @version: V1: 04-Jan-2011 (BAL)
    """
    try:
        if len(inval) > 1:
            if isinstance(inval, numpy.ndarray):
                return numpy.array([long(val.strftime('%Y%m%d')) for val in inval])
            else:
                return [long(val.strftime('%Y%m%d')) for val in inval]
    except:
        return long(inval.strftime('%Y%m%d'))

