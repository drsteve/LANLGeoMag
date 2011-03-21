# -*- coding: utf-8 -*-
from __future__ import division

"""
perform tracing to see if a file line is closed
"""
from ctypes import pointer
import math, numpy

from Lgm_Wrap import Lgm_Trace, LGM_OPEN_IMF, LGM_CLOSED, LGM_OPEN_N_LOBE
from Lgm_Wrap import LGM_OPEN_S_LOBE, LGM_INSIDE_EARTH, LGM_TARGET_HEIGHT_UNREACHABLE
from Lgm_Wrap import Lgm_B_OP77, Lgm_Set_Lgm_B_OP77, Lgm_Set_Coord_Transforms, Lgm_Set_Lgm_B_T89
from Lgm_Wrap import Lgm_Convert_Coords, GSM_TO_SM, WGS84_A
import Lgm_Vector
import Lgm_MagModelInfo
import Lgm_CTrans

def _simpleL(position, MagModelInfo):
    """
    code to return a simple measure of L
    Inputs:
    - GSM position of nothern most point (from Closed_Field)
    """
    # Get a simple measure of how big L is
    position_sm = Lgm_Vector.Lgm_Vector()
    Lgm_Convert_Coords( pointer(position), pointer(position_sm), GSM_TO_SM,
                       MagModelInfo.c );
    Lam = math.asin(position_sm.z/position_sm.magnitude())
    CosLam = math.cos(Lam)
    LSimple = (1.0+120.0/WGS84_A)/( CosLam*CosLam )
    return LSimple


# BAL: I belive the position is in GSM, check with Mike
def Closed_Field(*args, **kwargs):
    '''TODO: Need docstring'''
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
        mmi = MagEphemInfo.LstarInfo.contents.mInfo.contents
        dum = [MagEphemInfo.P.x, MagEphemInfo.P.y, MagEphemInfo.P.z]
        position = Lgm_Vector.Lgm_Vector(*dum)
    elif len(args) == 2:
        # input checking
        if kwargs['coord_system'] != 'GSM':
            raise(NotImplementedError('Different coord systems are not yet ready to use') )
        # could consider a Lgm_MagModelInfo param to use an existing one
        mmi = Lgm_MagModelInfo.Lgm_MagModelInfo()
        if kwargs['bfield'] == 'Lgm_B_OP77':
            Lgm_Set_Lgm_B_OP77( pointer(mmi) )
        elif kwargs['bfield'] == 'Lgm_B_T89':
            mmi.Kp = kwargs['Kp']
            Lgm_Set_Lgm_B_T89( pointer(mmi) )
        else:
            raise(NotImplementedError('Only Lgm_B_OP77 and Lgm_B_T89 implented so far'))

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
    if kwargs['extended_out']:
        return retstr, northern.tolist(), southern.tolist(), minB.tolist(), L
    else:
        return retstr
