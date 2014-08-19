"""
Overview
--------
Lgm_CTrans module, this contains the necessary code for coordinate
transformations in Lgm

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
                GSM_TO_WGS84, WGS84_TO_EDMAG, Lgm_Set_Coord_Transforms, Lgm_EDMAG_to_R_MLAT_MLON_MLT, \
                Lgm_Dipole_Tilt


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
    convert a python date or datetime object to a LanlGeoMag's (long int) date format

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
        if isinstance(inval, numpy.ndarray):
            retval = numpy.array([long(val.strftime('%Y%m%d')) for val in inval])
        else:
            retval = [long(val.strftime('%Y%m%d')) for val in inval]
        if len(retval)==1: retval=retval[0]
    except: #not iterable
        retval = long(inval.strftime('%Y%m%d'))
    return retval

def dateLongToDate(inval):
    """
    convert a python date or datetime object from a Date (long) object

    Parameters
    ----------
    inval : dateLong, int
        int/long to change to a datetime object

    Returns
    -------
    out : datetime, list
        datetime or list of datetime corresponding to the dateLongs
    """
    try:
        if isinstance(inval, numpy.ndarray):
            retval = numpy.array([datetime.datetime.strptime('{0}'.format(val), '%Y%m%d') for val in inval])
        else:
            retval = [datetime.datetime.strptime('{0}'.format(val), '%Y%m%d') for val in inval]
        if len(retval)==1: retval=retval[0]
    except: #not iterable
        retval = datetime.datetime.strptime('{0}'.format(inval), '%Y%m%d')
    return retval

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
        lst = [val.hour + val.minute/60 +
                                val.second/60/60 +
                                val.microsecond/60/60/1000000 for val in inval]
        if isinstance(inval, numpy.ndarray):
            retval = numpy.array(lst)
        else:
            retval = lst
        if len(retval)==1: retval=retval[0]
    except TypeError:
        retval = inval.hour + inval.minute/60 + inval.second/60/60 + inval.microsecond/60/60/1000000
    return retval


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
        cT = pointer(Lgm_CTrans())
        Lgm_Set_Coord_Transforms( dateToDateLong(dt), dateToFPHours(dt), cT)

        Lgm_Convert_Coords( pointer(Pgsm), pointer(Pwgs), GSM_TO_WGS84, cT )
        Lgm_Convert_Coords( pointer(Pwgs), pointer(Pmlt), WGS84_TO_EDMAG, cT )
        R, MLat, MLon, MLT = c_double(), c_double(), c_double(), c_double(),
        Lgm_EDMAG_to_R_MLAT_MLON_MLT( pointer(Pmlt),  pointer(R), pointer(MLat), pointer(MLon),
            pointer(MLT), cT)
        return MLT.value

    gsm_ = numpy.asanyarray(gsm)
    if isinstance(dt, datetime.datetime):
        dt_ = numpy.asanyarray([dt])
    else:
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
        if dt_.size==1:
            ans = dm.dmarray([doConv(gsm_, dt_)], attrs={'coord_system': 'EDMAG'})
        else:
            ans = dm.dmarray(doConv(gsm_, dt_), attrs={'coord_system': 'EDMAG'})
    return ans


def getDipoleTilt(date, degrees=False):
    """
    returns the dipole tile angle in radians

    Parameters
    ==========
    date : datetime or iterable of datetime
        the date for the computation

    Returns
    =======
    out : array or float
        array of dipole tile angles

    Examples
    ========
    from lgmpy import Lgm_CTrans
    import datetime
    Lgm_CTrans.getDipoleTilt(datetime.datetime(2000, 1, 1))
    # -0.45114989442541137
    """
    if hasattr(date, 'year'):
        date = [date]
    ans = numpy.empty(len(date), dtype=float)
    for i, d in enumerate(date):
        datelong = dateToDateLong(d)
        utc = dateToFPHours(d)
        ans[i] = Lgm_Dipole_Tilt(datelong, utc)
    if len(ans) == 1:
        ans = ans[0]
    if degrees: ans = numpy.rad2deg(ans)
    return ans
