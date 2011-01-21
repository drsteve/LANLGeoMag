"""
Lgm_CTrans module, this contains the necessary code for coordinate
transformations in Lgm

@author: Brian Larsen
@organization: LANL
@contact: balarsen@lanl.gov

@version: V1: 20-Dec-2010 (BAL)
"""

from __future__ import division

import math

import numpy
import spacepy.toolbox as tb

import Lgm_Vector
from Lgm_Wrap import Lgm_CTrans


class Lgm_Coords(list):
    """
    Base class for coordinate transforms

    @author: Brian Larsen
    @organization: LANL
    @contact: balarsen@lanl.gov

    @version: V1: 18-Jan-2011 (BAL)
    """
    def __init__(self, inval, system='GSM', units='Re'):
        if len(inval) != 3 and not isinstance(inval, list):
            raise(NotImplementedError('So far only one 3 element list is supported as input' ) )
        if units != 'Re':
            raise(NotImplementedError('Only Re units supoorted so far' ) )
        if tb.hypot(*inval) < 1.0:
            raise(ValueError('Invalid position'))
        if system != 'GSM':
            raise(NotImplementedError('Only GSM coordinated supoorted so far' ) )

        self.system = system
        self[:] = inval


class Lgm_CTrans(Lgm_CTrans):
    def __init__(self):
        self.nNutationTerms = 106;

         #
         #  Initialize Earth Orientation Parameters to defaults They can be
         #  over-ridden later with actual values.  This allows users to use the
         #  less accurate defaults if EOP data is unavailable or not
         #  needed/desired.
         #
        self.DUT1  = 0.0# ; // seconds
        self.xp    = 0.0#; // radians
        self.yp    = 0.0#; // radians
        self.ddPsi = 0.0#; // radians
        self.ddEps = 0.0#; // radians

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

def dateToFPHours(inval):
    """
    convert a python datetime object to a Floating point hours (double) object that
    LanlGeoMag Likes to use

    @author: Brian Larsen
    @organization: LANL
    @contact: balarsen@lanl.gov

    @version: V1: 06-Jan-2011 (BAL)
    """
    try:
        if len(inval) > 1:
            lst = [val.hour + val.minute/60 +
                                    val.second/60/60 +
                                    val.microsecond/60/60/1000000 for val in inval]
            if isinstance(inval, numpy.ndarray):
                return numpy.array(lst)
            else:
                return lst
    except:
        return inval.hour + inval.minute/60 + \
                                    inval.second/60/60 + \
                                    inval.microsecond/60/60/1000000
