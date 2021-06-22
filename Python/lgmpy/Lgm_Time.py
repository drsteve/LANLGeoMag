#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
Module with Python wrappers for LanlGeoMag time conversion routines

This should be either managed within classes or used as the backend to Spacepy's time module
'''

import datetime as dt
from ctypes import pointer, c_int, c_double
from . import Lgm_CTrans, Lgm_Vector, magcoords
from .Lgm_Wrap import Lgm_Date_to_JD, LGM_TIME_SYS_UTC, Lgm_init_ctrans, \
    Lgm_Convert_Coords, Lgm_Set_Coord_Transforms
try:
    import spacepy.time as spt
except ImportError:
    _spacepy = False
else:
    _spacepy = True

def datetime2Lgm(*args):
    '''Convenience function to convert Python datetime to Lgm date and utc format
    
    Parameters
    ----------
    time_in: datetime.datetime
        the input time to convert to a different system
    
    Outputs
    -------
    datelong: long
        LanlGeoMag long integer date representation (YYYYMMDD)
    utc: float
        LanlGeoMag floating point time-of-day (UTC) representation
    '''
    if len(args) == 1:
        time_in = args[0]
    else:
        raise ValueError('Invalid arguments supplied in function call') 
    
    if _spacepy:
        if isinstance(time_in, spt.Ticktock):
            raise NotImplementedError('SpacePy compatibility not yet enabled')
    try:
        datelong = Lgm_CTrans.dateToDateLong(time_in)
        utc = Lgm_CTrans.dateToFPHours(time_in)
    except AttributeError:
        raise TypeError("Date must be a datetime object")
    
    return datelong, utc
    
def Lgm2jd(datelong, utc, cTrans=None):
    '''
    convert from Lgm float UTC time representation to julian day
    '''
    if not cTrans:
        c = pointer(Lgm_CTrans.Lgm_CTrans())
    else:
        c = cTrans
    JD = Lgm_Date_to_JD(datelong, utc, c)
    return JD