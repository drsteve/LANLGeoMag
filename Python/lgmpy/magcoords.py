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

from lgmpy.Lgm_Wrap import Lgm_Set_Coord_Transforms, Lgm_Convert_Coords, WGS84_TO_GSM, SM_TO_GSM, Lgm_McIlwain_L, Lgm_Set_Lgm_B_OP77, Lgm_Set_Lgm_B_T89
from lgmpy import Lgm_Vector, Lgm_CTrans, Lgm_MagModelInfo
from _Bfield_dict import Bfield_dict

def LatLon2GSM(lat, lon, rad, time):
    XYZ = Lgm_Vector.SphToCart(lat, lon, rad)
    c = Lgm_CTrans.Lgm_CTrans()
    Ugsm = Lgm_Vector.Lgm_Vector(0,0,0)
    Ull = Lgm_Vector.Lgm_Vector(XYZ.x,XYZ.y, XYZ.z)
    Lgm_Set_Coord_Transforms(Lgm_CTrans.dateToDateLong(time),
                    Lgm_CTrans.dateToFPHours(time), pointer(c))
    Lgm_Convert_Coords(pointer(Ull), pointer(Ugsm), WGS84_TO_GSM, pointer(c))  # maybe this is just GEO
    return Ugsm.tolist()



def Lvalue(*args, **kwargs):
    '''
    Calling syntax: Lvalue(pos, datetime, )
    
    \param[in]        Date        Date in format (e.g. 20101231). 
    \param[in]        UTC         Universal Time (Coordinated) in decimal hours (e.g. 23.5).
    \param[in]        u           Position to compute L-shell.
    \param[in]        Alpha       Pitch angle to compute L for. In degrees.
    \param[in]        Type        Flag to indicate which alogorithm to use (0=original McIlwain; else use Hilton's formula).
    \param[out]       I           The integral invariant, I that was computed along the way.
    \param[out]       Bm          The mirror magnetic field value, Bm that was computed along the way.
    \param[out]       M           The dipole magnetic moment used to compute L = f(I, Bm, M)
    \param[in,out]    mInfo       Properly initialized Lgm_MagModelInfo structure. (A number of otherm usefull things will have been set in mInfo).
    
    \return           L           McIlwain L-shell parameter (a dimensioless number).
    
    '''
    defaults = {'alpha': 90,
                'Bfield': 'Lgm_B_T89',
                'method': 'Hilton',
                'Kp': 2,
                'coord_system': 'GSM'}

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
        raise(NotImplementedError("Only Bfield='OP77, T89' currently supported"))
    try:
        datelong = Lgm_CTrans.dateToDateLong(args[1])
        utc = Lgm_CTrans.dateToFPHours(args[1])
        Lgm_Set_Coord_Transforms( datelong, utc, mInfo.c) # dont need pointer as it is one
    except AttributeError:
        raise(TypeError("Date must be a datetime object"))
    #else:
        #ans['Epoch'] = datamodel.dmarray([args[1]])
    
    if kwargs['coord_system'] == 'SM':
        Psm = Lgm_Vector.Lgm_Vector(*args[0])
        Pgsm = Lgm_Vector.Lgm_Vector()
        Lgm_Convert_Coords( pointer(Psm), pointer(Pgsm), SM_TO_GSM, mInfo.c )
    else:
        Pgsm = Lgm_Vector.Lgm_Vector(*args[0])
    
    Iout = c_double()
    Bm = c_double()
    M = c_double()
     
    ans = Lgm_McIlwain_L(datelong, utc, pointer(Pgsm), kwargs['alpha'], 
                            method_dict[kwargs['method']],
                            pointer(Iout), pointer(Bm), pointer(M),
                            pointer(mInfo))
                            
    return {'L': ans, 'I': Iout.value}
    