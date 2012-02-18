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

from ctypes import pointer, c_double
import datetime
import itertools

import spacepy.datamodel as dm
import numpy

import Lgm_MagModelInfo
import Lgm_Vector
from Lgm_Wrap import Lgm_CTrans, Lgm_ctransDefaults, Lgm_free_ctrans_children, Lgm_Convert_Coords, \
                GSM_TO_WGS84, WGS84_TO_EDMAG, Lgm_Set_Coord_Transforms, Lgm_EDMAG_to_R_MLAT_MLON_MLT


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
        
    def __del__(self):
        Lgm_free_ctrans_children(pointer(self))
    
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
        if isinstance(inval, numpy.ndarray):
            inval = inval.item()
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
        if isinstance(inval, numpy.ndarray):
            inval = inval.item()
        return inval.hour + inval.minute/60 + \
                                    inval.second/60/60 + \
                                    inval.microsecond/60/60/1000000


def GSMtoMLT(gsm, dt):
    """
    convert GSM values to MLT in the lgm way
    
    Parameters
    ----------
    gsm : array_like
        Nx3 array_like of the GSM position
    dt : array_like
        N elementarray_like of datetime objects
    
    Returns
    -------
    out : numpy.array
        N element array of the MLT values
    """
    def doConv(gsm, dt):
        Pgsm = Lgm_Vector.Lgm_Vector(*gsm)
        Pwgs = Lgm_Vector.Lgm_Vector()
        Pmlt = Lgm_Vector.Lgm_Vector()
        # can use a smaller structure here        
        mmi = Lgm_MagModelInfo.Lgm_MagModelInfo()
        Lgm_Set_Coord_Transforms( dateToDateLong(dt), dateToFPHours(dt), mmi.c) # dont need pointer as it is one

        Lgm_Convert_Coords( pointer(Pgsm), pointer(Pwgs), GSM_TO_WGS84, mmi.c )
        Lgm_Convert_Coords( pointer(Pwgs), pointer(Pmlt), WGS84_TO_EDMAG, mmi.c )
        R, MLat, MLon, MLT = c_double(), c_double(), c_double(), c_double(),
        Lgm_EDMAG_to_R_MLAT_MLON_MLT( pointer(Pmlt),  pointer(R), pointer(MLat), pointer(MLon),
            pointer(MLT), mmi.c)
        return MLT.value
        
    gsm_ = numpy.asanyarray(gsm)
    dt_ = numpy.asanyarray(dt)
    if gsm_.ndim == 2:
        if gsm_.shape[1] != 3:
            raise(ValueError("Invalid vector shape"))
        if gsm_.shape[0] != dt_.size:
            if dt_.size == 1:
                dt_ = dm.dmarray([dt_]*gsm_.shape[0])
            else:
                raise(ValueError("Array size mismatch"))
        ans = dm.dmarray(numpy.empty(len(dt_)), dtype=numpy.double, attrs={'coord_system': 'EDMAG'})
        for ii, (gsm_val, dt_val) in enumerate(itertools.izip(gsm_, dt_)):
            ans[ii] = doConv(gsm_val, dt_val)
    else:
        ans = dm.dmarray(doConv(gsm_, dt_), attrs={'coord_system': 'EDMAG'})
    return ans

