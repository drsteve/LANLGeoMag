# -*- coding: utf-8 -*-
"""
Overview
--------
this is a coordinate transformation class based on LangGeoMag

It is still just a a partial solution


"""
__author__ = 'Steve Morley, Brian Larsen (Python), Mike Henderson (C) - LANL'

import itertools
from ctypes import pointer, c_double

import numpy as np
from spacepy import datamodel


from lgmpy.Lgm_Wrap import Lgm_Set_Coord_Transforms, Lgm_Convert_Coords, Lgm_McIlwain_L
from Lgm_Wrap import TEME_TO_WGS84, WGS84_TO_GSM, GSM_TO_WGS84, WGS84_TO_GSE, GSE_TO_WGS84, SM_TO_GSM, GSM_TO_SM, GSM_TO_GSE, GSE_TO_GSM, MOD_TO_GSM, GSE_TO_SM, SM_TO_GSE
from lgmpy import Lgm_Vector, Lgm_CTrans, Lgm_MagModelInfo
from lgmpy.Lstar import Lstar_Data

from _Bfield_dict import Bfield_dict

conv_dict = {'SM_GSM': SM_TO_GSM,
                 'GSM_SM': GSM_TO_SM,
                 'SM_GSE': SM_TO_GSE,
                 'GSE_SM': GSE_TO_SM,
                 'WGS84_GSM': WGS84_TO_GSM,
                 'GSM_WGS84': GSM_TO_WGS84,
                 'WGS84_GSE': WGS84_TO_GSE,
                 'GSE_WGS84': GSE_TO_WGS84,
                 'GSM_GSE': GSM_TO_GSE,
                 'GSE_GSM': GSE_TO_GSM,
                 'TEME_WGS84': TEME_TO_WGS84,
                 'MOD_GSM' : MOD_TO_GSM,
                 'SM_GSE' : SM_TO_GSE}

def coordTrans(*args):
    '''
    Convert coordinates between almost any system using LanlGeoMag

    Parameters
    ----------
    position : list
        a three element vector of positions in input coord system
    time : datetime
        a datimetime object representing the time at the desired conversion
    system_in : str
        a string giving the acronym for the input coordinate system
    system_out : str
        a string giving the acronym for the desired output coordinate system

    Returns
    -------
    out : list
        3-element list of the converted coordinate

    Examples
    --------
    >>> from lgmpy import magcoords
    >>> import datetime
    >>> magcoords.coordTrans([-4,0,0], datetime.datetime(2009,1,1),'SM','GSM')
    [-3.60802691..., 2.5673907444...e-16, -1.72688788616...]
    >>> magcoords.coordTrans([-3.608026916281573, 2.5673907444456745e-16, -1.7268878861662329], datetime.datetime(2009,1,1),'GSM','SM')
    [-3.99999999..., 4.0592529337...e-16, 8.8817841970...3e-16]

    TODO
    ----
    extend interface to get necessary args from a MagModel or cTrans structure
    '''

    def doConversion(pos_in, to_str, cdict, cTrans):
        """
        Function that does the conversion

        .. warning::

            This function is not called by the user

        Parameters
        ==========
        pos_in : Lgm_Vector
            the input position to convert
        to_str : str
            a string giving the acronym for the output coordinate system
        cdict : dict
            the coordinate transform dictionary
        cTrans : Lgm_CTrans
            the Lgm_CTrans object used in the conversion

        Returns
        =======
        out : Lgm_Vector
            vector of the new position
        """

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
        Lgm_Set_Coord_Transforms( datelong, utc, mInfo.c) # don't need pointer as it is one
    except AttributeError:
        raise(TypeError("Date must be a datetime object"))

    to_str = in_sys+'_'+out_sys

    ## do this as WGS uses Cartesian but needs to be converted from desired spherical input
    if 'WGS84' in in_sys:
        XYZ = Lgm_Vector.SphToCart(*pos_in)
        SPH = Lgm_Vector.Lgm_Vector(XYZ.x,XYZ.y, XYZ.z)
        Pout = doConversion(SPH, to_str, cdict=conv_dict, cTrans=mInfo.c)
    else:
        Pout = doConversion(pos_in, to_str, conv_dict, mInfo.c)

    if 'WGS84' in out_sys:
        nlat, nlon, nrad = Lgm_Vector.CartToSph(*Pout.tolist())
        Pout = Lgm_Vector.Lgm_Vector(nlat, nlon, nrad)

    return Pout.tolist()

def Lvalue(*args, **kwargs):
    '''
    Function to return the L-value of a position using either McIlwain or Hilton
    approximation

    Parameters
    ==========
    pos : list
        3-element position int he specified coord_system
    time : datetime
        date and time for the calculation
    alpha : float, optional
        the pitch angle for the L calculation, default=90
    Bfield : str, optional
        the magnetic field model to use, default=Lgm_B_T89
    method : str, optional
        the L-value formula to use, McIlwain or Hilton, default=Hilton
    Kp : int
        Kp index value for the calculation
    coord_system : str
        the input coordinate system, default=GSM
    extended_out : bool
        keyword specifying short or extended output, default=False

    Returns
    =======
    out : dict
        dictionary of the values, see examples for documentation of dictionary

    Examples
    ========
    >>> from lgmpy import magcoords
    >>> import datetime
    >>> magcoords.Lvalue([3, 0, 1], datetime.datetime(2000, 1, 1), extended_out=False)
    {'I': 0.2434969602..., 'L': 3.195481841...}

    The ``extended_out=False`` output is:
        - I : I value at the given point
        - L : L value at the given point

    >>> from lgmpy import magcoords
    >>> import datetime
    >>> magcoords.Lvalue([3, 0, 1], datetime.datetime(2000, 1, 1), extended_out=True)
    {'Blocal': 1024.1142193703838,
    'Bmin': 921.8869150...,
    'Bmirr': 1024.1142193...,
    'I': 0.24349696021...,
    'L': 3.1954818410...,
    'M': 30119.614287...}

    The ``extended_out=True`` output is:
        - I : I value at the given point
        - L : L value at the given point
        - Bmin : minimum B at the input position (nT)
        - Bmirror : mirror B at the input position (nT)
        - M : TODO what exactly and I, magnetic moment?
    '''
    defaults = {'alpha': 90.,
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
        sunPos_GSM = coordTrans(mInfo.c.contents.Sun, args[1], 'MOD', kwargs['coord_system'])
        SunMlon = np.rad2deg(np.arctan2( sunPos_GSM[1], sunPos_GSM[0]))
        SunMlon += 180.0 # flip to midnight

        MLON = np.rad2deg(np.arctan2( Pgsm.y, Pgsm.x))
        if (MLON < 0.0):
            MLON += 360.0

        MLT = np.mod( (MLON-SunMlon)/15.0+24.0, 24.0 )
        return {'L': ans, 'I': Iout.value, 'Bmin': mInfo.Bmin, 'Blocal': mInfo.Blocal,
                'Bmirr': mInfo.Bm, 'M': M.value, 'MLon':MLON, 'MLT':MLT}
    else:
        return {'L': ans, 'I': Iout.value}
