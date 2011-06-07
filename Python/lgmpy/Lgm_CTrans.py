"""
Overview
--------
Lgm_CTrans module, this contains the necessary code for coordinate
transformations in Lgm

    Authors
    -------
    Brian Larsen - LANL
"""

from __future__ import division

from ctypes import pointer
import datetime

import numpy
import spacepy.toolbox as tb

from Lgm_Wrap import Lgm_CTrans, Lgm_ctransDefaults


class Lgm_Coords(list):
    """
    Base class for coordinate transforms

    Parameters
    ----------
    inval : list
        3-element list containing a position
    system : str, optional
        coordinate system of the input, default='GSM'
    units : str, optional
        units of the position, default=Re

    Returns
    -------
    out : Lgm_Coords class
    """
    def __init__(self, inval, system='GSM', units='Re'):
        if len(inval) != 3 and not isinstance(inval, list):
            raise(NotImplementedError('So far only one 3 element list is supported as input' ) )
        if units != 'Re':
            raise(NotImplementedError('Only Re units supported so far' ) )
        if inval[0] < 1 and inval[1] < 1 and inval[2] < 1:
            raise(ValueError('Invalid position'))
        if system != 'GSM':
            raise(NotImplementedError('Only GSM coordinated supported so far' ) )

        self.system = system
        self[:] = inval

class Lgm_CTrans(Lgm_CTrans):
    """
    Class to wrap the C structure Lgm_CTrans.  Default values are filled by the
    library

    Parameters
    ----------
    verbose : bool, optional
        switch to the C library on verbose

    Returns
    -------
    out : Lgm_CTrans class
    """
    def __init__(self, Verbose=False):
        # initialize to the values set in c so we don't have to maintain two places
        Lgm_ctransDefaults(pointer(self), Verbose)

def dateToDateLong(inval):
    """
    convert a python date or datetime object to a Date (long) object that
    LanlGeoMag Likes to use

    Parameters
    ----------
    inval : datetime, ndarray
        datetime object to change to a long (or an array of datetimes to change)

    Returns
    -------
    out : long, list
        long or list of longs corresponding to the datetime

    Examples
    --------
    >>> from lgmpy import Lgm_CTrans
    >>> import datetime
    >>> Lgm_CTrans.dateToDateLong(datetime.datetime(2000, 12, 13))
    20001213L
    """
    try:
        if len(inval) > 1:
            if isinstance(inval, numpy.ndarray):
                return numpy.array([long(val.strftime('%Y%m%d')) for val in inval])
            else:
                return [long(val.strftime('%Y%m%d')) for val in inval]
    except:
        return long(inval.strftime('%Y%m%d'))

def dateLongToDate(inval):
    """
    convert a python date or datetime object to a Date (long) object that
    LanlGeoMag Likes to use

    Parameters
    ----------
    inval : dateLong, int
        int to change to a datetime object

    Returns
    -------
    out : datetime, list
        datetime or list of datetime corresponding to the dateLongs
    """
    try:
        if len(inval) > 1:
            if isinstance(inval, numpy.ndarray):
                return numpy.array([datetime.datetime.strptime(str(val), '%Y%m%d') for val in inval])
            else:
                return [datetime.datetime.strptime(str(val), '%Y%m%d') for val in inval]
    except:
        return datetime.datetime.strptime(str(inval), '%Y%m%d')

def dateToFPHours(inval):
    """
    convert a python datetime object to a Floating point hours (double) object that
    LanlGeoMag Likes to use

    Parameters
    ----------
    inval : datetime, ndarray
        datetime object to change to floating point hours
        (or an array of datetimes to change)

    Returns
    -------
    out : float, list
        floating point hours or a list of fp hours

    Examples
    --------
    >>> from lgmpy import Lgm_CTrans
    >>> import datetime
    >>> Lgm_CTrans.dateToFPHours(datetime.datetime(2000, 12, 13))
    0.0
    >>> Lgm_CTrans.dateToFPHours(datetime.datetime(2000, 12, 13, 12))
    12.0
    >>> Lgm_CTrans.dateToFPHours(datetime.datetime(2000, 12, 13, 12, 30))
    12.5
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
    except TypeError:
        return inval.hour + inval.minute/60 + \
                                    inval.second/60/60 + \
                                    inval.microsecond/60/60/1000000
