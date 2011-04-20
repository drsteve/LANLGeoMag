# -*- coding: utf-8 -*-
"""
this is a cheezy scrpt that needs to be wrapped up with the rest
but I need it now and dont feel like it today

@author: Brian Larsen
@organization: LANL
@contact: balarsen@lanl.gov

@version: V1: 11-Jan-2011 (BAL)
"""

import itertools
from ctypes import pointer, c_double

import numpy
from spacepy import datamodel

try:
    #this block makes sure the tests in working directory run on local files not installed versions
    #TODO: I don't like this -- so how can this be done better????
    from .Lgm_Wrap import Lgm_Set_Coord_Transforms, Lgm_Convert_Coords, Lgm_McIlwain_L
    from .Lgm_Wrap import TEME_TO_WGS84, WGS84_TO_GSM, GSM_TO_WGS84, WGS84_TO_GSE, GSE_TO_WGS84, SM_TO_GSM, GSM_TO_SM, GSM_TO_GSE, GSE_TO_GSM
    from . import Lgm_Vector, Lgm_CTrans, Lgm_MagModelInfo
    from .Lstar import Lstar_Data
except:
    from lgmpy.Lgm_Wrap import Lgm_Set_Coord_Transforms, Lgm_Convert_Coords, Lgm_McIlwain_L
    from Lgm_Wrap import TEME_TO_WGS84, WGS84_TO_GSM, GSM_TO_WGS84, WGS84_TO_GSE, GSE_TO_WGS84, SM_TO_GSM, GSM_TO_SM, GSM_TO_GSE, GSE_TO_GSM
    from lgmpy import Lgm_Vector, Lgm_CTrans, Lgm_MagModelInfo
    from lgmpy.Lstar import Lstar_Data
    
from _Bfield_dict import Bfield_dict

conv_dict = {'SM_GSM': SM_TO_GSM,
                 'GSM_SM': GSM_TO_SM,
                 'WGS84_GSM': WGS84_TO_GSM,
                 'GSM_WGS84': GSM_TO_WGS84,
                 'WGS84_GSE': WGS84_TO_GSE,
                 'GSE_WGS84': GSE_TO_WGS84,
                 'GSM_GSE': GSM_TO_GSE,
                 'GSE_GSM': GSE_TO_GSM,
                 'TEME_WGS84': TEME_TO_WGS84}

def coordTrans(*args):
    ''' Convert coordinates between almost any system using LanlGeoMag
    
    Input arguments:
    ----------------
    position - [list] a three element vector of positions in input coord system
    time -  [datetime] a datimetime object representing the time at the desired conversion
    system_in - [str] a string giving the acronym for the input coordinate system
    system_out - [str] a string giving the acronym for the desired output coordinate system
    
    Example:
    --------
    >>> import datetime
    >>> coordTrans([-4,0,0], datetime.datetime(2009,1,1),'SM','GSM')
    [-3.608026916281573, 2.5673907444456745e-16, -1.7268878861662329]
    >>> coordTrans([-3.608026916281573, 2.5673907444456745e-16, 
        -1.7268878861662329], datetime.datetime(2009,1,1),'GSM','SM')
    [-3.9999999999999991, 4.0592529337857286e-16, 8.8817841970012523e-16]
    
    
    TODO: extend interface to get necessary args from a MagModel or cTrans structure
    '''
    
    def doConversion(pos_in, to_str, cdict, cTrans):
        if to_str not in cdict: raise NotImplementedError('Coordinate transform not implemented')
        try:
            Pin = Lgm_Vector.Lgm_Vector(*pos_in)
        except:
            Pin = pos_in
        Pout = Lgm_Vector.Lgm_Vector()
        Lgm_Convert_Coords( pointer(Pin), pointer(Pout), cdict[to_str], cTrans )      
        return Pout
    
    #change args to named vars to make easier to process different input types
    pos_in = args[0]
    time_in = args[1]
    in_sys = args[2]
    out_sys = args[3]
    
    # change datetime to Lgm Datelong and UTC
    mInfo = Lgm_MagModelInfo.Lgm_MagModelInfo()
    try:
        datelong = Lgm_CTrans.dateToDateLong(time_in)
        utc = Lgm_CTrans.dateToFPHours(time_in)
        Lgm_Set_Coord_Transforms( datelong, utc, mInfo.c) # dont need pointer as it is one
    except AttributeError:
        raise(TypeError("Date must be a datetime object"))
    
    to_str = in_sys+'_'+out_sys

    ## do this as WGS uses cartesian but needs to be converted from desired spherical input
    if 'WGS' in in_sys:
        XYZ = Lgm_Vector.SphToCart(*pos_in)
        SPH = Lgm_Vector.Lgm_Vector(XYZ.x,XYZ.y, XYZ.z)
        Pout = doConversion(SPH, to_str, cdict=conv_dict, cTrans=mInfo.c)
    else:
        Pout = doConversion(pos_in, to_str, conv_dict, mInfo.c)
    
    if 'WGS' in out_sys:
        nlat, nlon, nrad = Lgm_Vector.CartToSph(*Pout.tolist())
        Pout = Lgm_Vector.Lgm_Vector(nlat, nlon, nrad)
    
    return Pout.tolist()

def Lvalue(*args, **kwargs):
    '''
    Calling syntax: Lvalue(pos, datetime, )
    
    TODO: Add docstring
    '''
    defaults = {'alpha': 90,
                'Bfield': 'Lgm_B_T89',
                'method': 'Hilton',
                'Kp': 2,
                'coord_system': 'GSM',
                'extended_out': False}

    #replace missing kwargs with defaults
    for dkey in defaults:
        if dkey not in kwargs:
            kwargs[dkey] = defaults[dkey]
    
    method_dict = {'Hilton': 1, 'McIlwain': 0} 
    
    # change datetime to Lgm Datelong and UTC
    mInfo = Lgm_MagModelInfo.Lgm_MagModelInfo()
    mInfo.Kp = kwargs['Kp']
    try:
        Bfield_dict[kwargs['Bfield']](pointer(mInfo))
    except KeyError:
        raise(NotImplementedError("Only Bfield=%s currently supported" % Bfield_dict.keys()))
    try:
        datelong = Lgm_CTrans.dateToDateLong(args[1])
        utc = Lgm_CTrans.dateToFPHours(args[1])
        Lgm_Set_Coord_Transforms( datelong, utc, mInfo.c) # dont need pointer as it is one
    except AttributeError:
        raise(TypeError("Date must be a datetime object"))
    #else:
        #ans['Epoch'] = datamodel.dmarray([args[1]])
    
    if kwargs['coord_system'] != 'GSM':
        Pin = Lgm_Vector.Lgm_Vector(*args[0])
        Pgsm = Lgm_Vector.Lgm_Vector()
        conv_str = kwargs['coord_system']+'_GSM'
        try:
            #TODO: Fix this so it uses the (generalized) coordTrans function above
            Lgm_Convert_Coords( pointer(Pin), pointer(Pgsm), conv_dict[conv_str], mInfo.c )
        except:
            raise NotImplementedError('Coordinate system not implemented')
    else:
        Pgsm = Lgm_Vector.Lgm_Vector(*args[0])
    
    Iout = c_double()
    Bm = c_double()
    M = c_double()
    ans = Lgm_McIlwain_L(datelong, utc, pointer(Pgsm), kwargs['alpha'], 
                            method_dict[kwargs['method']],
                            pointer(Iout), pointer(Bm), pointer(M),
                            pointer(mInfo))

    #TODO: decide on format for output -- perhaps use datamodel and have method as attribute?
    #maybe only return I for extended_out flag=True
    if kwargs['extended_out']:
        return {'L': ans, 'I': Iout.value, 'Bmin': mInfo.Bmin, 'Blocal': mInfo.Blocal, 
                'Bmirr': mInfo.Bm, 'M': M.value}
    else:
        return {'L': ans, 'I': Iout.value}
    