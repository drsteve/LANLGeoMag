# -*- coding: utf-8 -*-
from __future__ import division

"""
Overview
--------
perform tracing to see if a file line is closed

    Authors
    -------
    Steve Morley, Brian Larsen (python)
"""
from ctypes import pointer
import math, numpy

from Lgm_Wrap import Lgm_Trace, LGM_OPEN_IMF, LGM_CLOSED, LGM_OPEN_N_LOBE, LGM_BAD_TRACE
from Lgm_Wrap import LGM_OPEN_S_LOBE, LGM_INSIDE_EARTH, LGM_TARGET_HEIGHT_UNREACHABLE
from Lgm_Wrap import Lgm_Set_Lgm_B_OP77, Lgm_Set_Coord_Transforms, Lgm_Set_Lgm_B_T89
from Lgm_Wrap import Lgm_Convert_Coords, GSM_TO_SM, WGS84_A
import Lgm_Vector
import Lgm_MagModelInfo
import Lgm_CTrans
from _Bfield_dict import Bfield_dict

def _simpleL(position, MagModelInfo):
    """
    code to return a simple measure of L
    Inputs:
    - GSM position of northern most point (from Closed_Field)
    """
    # Get a simple measure of how big L is
    position_sm = Lgm_Vector.Lgm_Vector()
    Lgm_Convert_Coords( pointer(position), pointer(position_sm), GSM_TO_SM,
                       MagModelInfo.c );
    Lam = math.asin(position_sm.z/position_sm.magnitude())
    CosLam = math.cos(Lam)
    LSimple = (1.0+120.0/WGS84_A)/( CosLam*CosLam )
    return LSimple


# BAL: I believe the position is in GSM, check with Mike
def Closed_Field(*args, **kwargs):
    """
    Function to see if a field line is closed

    Either MagEphem or pos and date must be specified

    Parameters
    ----------
    MagEphem : Lgm_MagEphemInfo, optional
        If a populated Lgm_MagEphemInfo class is passed in the data is pulled
        from it
    pos : list, optional
        3-element list of the position in system coord_system
    date : datetime, optional
        date and time of the calculation
    height : float, optional
        height above the earth to consider a particle lost [km], default=100
    tol1 : float, optional
        TODO what do I set?  default=0.01
    tol2 : float, optional
        TODO what do I set?  default=1e-7
    bfield : str, optional
        The magnetic field model to use, default=Lgm_B_T89
    Kp : int, optional
        Kp index for the calculation, default=2
    coord_system : str
        the coordinate system of the input position, default=GSM
    extended_out : bool
        switch on extended vs regular output, see examples for details,
        default=False

    Returns
    -------
    out : str
        a string with the open of closed value for the input
            - LGM_OPEN_IMF
            - LGM_CLOSED
            - LGM_OPEN_N_LOBE
            - LGM_OPEN_S_LOBE
            - LGM_INSIDE_EARTH
            - LGM_TARGET_HEIGHT_UNREACHABLE

    Examples
    --------
    >>> from lgmpy import Closed_Field
    >>> import datetime
    >>> Closed_Field([3,1,0], datetime.datetime(2000, 12, 3))
    'LGM_CLOSED'
    >>> Closed_Field([6,1,12], datetime.datetime(2000, 12, 3))
    'LGM_OPEN_IMF'
    >>> Closed_Field([-16,1,5], datetime.datetime(2000, 12, 3))
    'LGM_OPEN_N_LOBE'
    """

    defaults = {'height': 100,
                'tol1': 0.01,
                'tol2': 1e-7,
                'bfield': 'Lgm_B_T89',
                'Kp': 2,
                'coord_system': 'GSM',
                'extended_out': False}

    #replace missing kwargs with defaults
    for dkey in defaults:
        if dkey not in kwargs:
            kwargs[dkey] = defaults[dkey]

    northern = Lgm_Vector.Lgm_Vector()
    southern = Lgm_Vector.Lgm_Vector()
    minB = Lgm_Vector.Lgm_Vector()
    #check for call w/ Lgm_MagModelInfo
    if len(args) == 1:
        MagEphemInfo = args[0]
        try:
            mmi = MagEphemInfo.LstarInfo.contents.mInfo.contents
        except AttributeError:
            raise(RuntimeError('Incorrect arguments specified'))
        dum = [MagEphemInfo.P.x, MagEphemInfo.P.y, MagEphemInfo.P.z]
        position = Lgm_Vector.Lgm_Vector(*dum)
    elif len(args) == 2:
        # input checking
        if kwargs['coord_system'] != 'GSM':
            raise(NotImplementedError('Different coord systems are not yet ready to use') )
        # could consider a Lgm_MagModelInfo param to use an existing one
        mmi = Lgm_MagModelInfo.Lgm_MagModelInfo()
        mmi.Kp = kwargs['Kp']
        try:
            Bfield_dict[kwargs['bfield']](pointer(mmi))
        except KeyError:
            raise(NotImplementedError("Only Bfield=%s currently supported" % Bfield_dict.keys()))

        datelong = Lgm_CTrans.dateToDateLong(args[1])
        utc = Lgm_CTrans.dateToFPHours(args[1])
        Lgm_Set_Coord_Transforms( datelong, utc, mmi.c) # dont need pointer as it is one

        try:
            position = Lgm_Vector.Lgm_Vector(*args[0])
        except TypeError:
            raise(TypeError('position must be an iterable') )
    else:
        raise(RuntimeError('Incorrect number of arguments specified'))

    ans = Lgm_Trace(pointer(position),
            pointer(southern),
            pointer(northern),
            pointer(minB), kwargs['height'], kwargs['tol1'], kwargs['tol2'], pointer(mmi) )

    L = numpy.nan #default to this is field not closed
    if ans == LGM_OPEN_IMF:
        retstr = 'LGM_OPEN_IMF'
    elif ans == LGM_CLOSED:
        retstr = 'LGM_CLOSED'
        L = _simpleL(northern, mmi)
    elif ans == LGM_OPEN_N_LOBE:
        retstr = 'LGM_OPEN_N_LOBE'
    elif ans == LGM_OPEN_S_LOBE:
        retstr = 'LGM_OPEN_S_LOBE'
    elif ans == LGM_INSIDE_EARTH:
        retstr = 'LGM_INSIDE_EARTH'
    elif ans == LGM_TARGET_HEIGHT_UNREACHABLE:
        retstr = 'LGM_TARGET_HEIGHT_UNREACHABLE'
    elif ans == LGM_BAD_TRACE:
        retstr = 'LGM_BAD_TRACE'
    if kwargs['extended_out']:
        return retstr, northern.tolist(), southern.tolist(), minB.tolist(), L
    else:
        return retstr
