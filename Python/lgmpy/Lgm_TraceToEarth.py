"""
Overview
--------
Python implementation of the LanlGeoMag TraceToEarth routine

"""
__author__ = 'Brian Larsen, Mike Henderson - LANL'

import datetime
import ctypes
from ctypes import pointer
import warnings

import numpy as np

import MagData
from Lgm_Wrap import LGM_CDIP, LGM_EDIP, LGM_IGRF, Lgm_Set_Coord_Transforms, \
    Lgm_B_T89, Lgm_Set_Lgm_B_cdip_InternalModel, Lgm_Set_Lgm_B_edip_InternalModel, Lgm_Set_Lgm_B_IGRF_InternalModel, LGM_EXTMODEL_NULL, LGM_EXTMODEL_T87, LGM_EXTMODEL_T89, LGM_EXTMODEL_T89c, LGM_EXTMODEL_T96, LGM_EXTMODEL_T01S, LGM_EXTMODEL_T02, LGM_EXTMODEL_TS04, LGM_EXTMODEL_TS07, LGM_EXTMODEL_OP77, Lgm_Set_Lgm_B_T89, Lgm_Set_Lgm_B_TS04, Lgm_Set_Lgm_B_OP77
from Lgm_Wrap import LGM_OPEN_IMF, LGM_CLOSED, LGM_OPEN_N_LOBE, LGM_OPEN_S_LOBE, LGM_INSIDE_EARTH, LGM_TARGET_HEIGHT_UNREACHABLE, LGM_BAD_TRACE
from Lgm_Wrap import Lgm_TraceToEarth


import Lgm_Vector
import Lgm_CTrans
import Lgm_MagModelInfo
from utils import pos2Lgm_Vector


# int Lgm_TraceToEarth( Lgm_Vector *u, Lgm_Vector *v, double TargetHeight, double sgn, double tol, Lgm_MagModelInfo *Info ) {
class Lgm_TraceToEarth_py(MagData.MagData):
    """
    Python implementation of the LanlGeoMag TraceToEarth routine.  This is
    the full wrapper, most users will not use this.

    Parameters
    ----------
    pos : list
        3-element list of the position to do the calculation
    time : datetime
        date and time of when to do the calculation
    direction : str, optional
        the direction to trace, default='North'
    coord_system : str, optional
        the coordinate system of the position, default=GSM
    INTERNAL_MODEL : str, optional
        the internal magnetic field model to use, default=LGM_IGRF
    EXTERNAL_MODEL : str, optional
        the external magnetic field model to use, default=LGM_T89
    MAGMODEL_ARGS : dict, optional
        kwarguments to the magnetic field model
    TargetHeight : float, optional
        target height for the tracing, default=120km
    tol : float, optional
        tolerance for tracing, default=1e-7

    Returns
    -------
    out : np.ndarray
        array of the footpoint in 

    See Also
    --------
    Lgm_TraceToEarth.TraceToEarth
    """    
    def __init__(self, pos, time, direction='NORTH', coord_system = 'GSM', INTERNAL_MODEL='LGM_IGRF', EXTERNAL_MODEL='LGM_EXTMODEL_T89', MAGMODEL_ARGS=None, TargetHeight=120, tol=1e-7):
        super(Lgm_TraceToEarth_py, self).__init__(Position=pos, Epoch=time, direction=direction, coord_system = coord_system, INTERNAL_MODEL=INTERNAL_MODEL, EXTERNAL_MODEL=EXTERNAL_MODEL,
                                               TargetHeight=TargetHeight,
                                               tol=tol)
        if not isinstance(pos, Lgm_Vector.Lgm_Vector) and \
            not isinstance(pos, list):
            raise(TypeError('pos must be a Lgm_Vector or list') )
        self._Vpos = pos2Lgm_Vector(pos)
        # time must be a datetime
        if not isinstance(time, datetime.datetime) and \
            not isinstance(time, list):
            raise(TypeError('time must be a datetime or list of datetime') )

        if direction[0].upper() == 'N':
            direction=1.0
        elif direction[0].upper() == 'S':
            direction=-1.0
        else:
            raise(ValueError("Did not understand direction, must be N or S"))

                    
        if INTERNAL_MODEL not in (LGM_CDIP,
                                  LGM_EDIP,
                                  LGM_IGRF) and \
            INTERNAL_MODEL not in ('LGM_CDIP',
                                  'LGM_EDIP',
                                  'LGM_IGRF'):
            raise(ValueError('INTERNAL_MODEL must be LGM_CDIP, LGM_EDIP, or LGM_IGRF') )
        if isinstance(INTERNAL_MODEL, str):
            INTERNAL_MODEL = eval(INTERNAL_MODEL)
        self.attrs['internal_model'] = INTERNAL_MODEL

        if coord_system != 'GSM':
            raise(NotImplementedError('Different coord systems are not yet ready to use') )

        self._mmi = Lgm_MagModelInfo.Lgm_MagModelInfo()
        
        # and actually set the internal model in Lgm
        if self.attrs['internal_model'] == LGM_CDIP:
            Lgm_Set_Lgm_B_cdip_InternalModel(pointer(self._mmi))
        elif self.attrs['internal_model'] == LGM_EDIP:
            Lgm_Set_Lgm_B_edip_InternalModel(pointer(self._mmi))
        elif self.attrs['internal_model'] == LGM_IGRF:
            Lgm_Set_Lgm_B_IGRF_InternalModel(pointer(self._mmi))

        # set the external field model
        if EXTERNAL_MODEL not in (LGM_EXTMODEL_NULL,
                                  LGM_EXTMODEL_T87,
                                  LGM_EXTMODEL_T89,
                                  LGM_EXTMODEL_T89c,
                                  LGM_EXTMODEL_T96,
                                  LGM_EXTMODEL_T01S,
                                  LGM_EXTMODEL_T02,
                                  LGM_EXTMODEL_TS04,
                                  LGM_EXTMODEL_TS07,
                                  LGM_EXTMODEL_OP77) and \
          EXTERNAL_MODEL not in ('LGM_EXTMODEL_NULL',
                                  'LGM_EXTMODEL_T87',
                                  'LGM_EXTMODEL_T89',
                                  'LGM_EXTMODEL_T89c',
                                  'LGM_EXTMODEL_T96',
                                  'LGM_EXTMODEL_T01S',
                                  'LGM_EXTMODEL_T02',
                                  'LGM_EXTMODEL_TS04',
                                  'LGM_EXTMODEL_TS07',
                                  'LGM_EXTMODEL_OP77'):
            raise(ValueError('INTERNAL_MODEL must be {0}'.format('LGM_EXTMODEL_NULL',
                                  'LGM_EXTMODEL_T87',
                                  'LGM_EXTMODEL_T89',
                                  'LGM_EXTMODEL_T89c',
                                  'LGM_EXTMODEL_T96',
                                  'LGM_EXTMODEL_T01S',
                                  'LGM_EXTMODEL_T02',
                                  'LGM_EXTMODEL_TS04',
                                  'LGM_EXTMODEL_TS07',
                                  'LGM_EXTMODEL_OP77')) )
        if isinstance(EXTERNAL_MODEL, str):
            EXTERNAL_MODEL = eval(EXTERNAL_MODEL)
        self.attrs['external_model'] = EXTERNAL_MODEL
          
        # and actually set the internal model in Lgm
        if self.attrs['external_model'] == LGM_EXTMODEL_T89:
            Lgm_Set_Lgm_B_T89(pointer(self._mmi))
            self._mmi.Kp = MAGMODEL_ARGS['Kp']
        elif self.attrs['external_model'] == LGM_EXTMODEL_OP77:
            Lgm_Set_Lgm_B_OP77(pointer(self._mmi))
        else:
            raise(NotImplementedError("The external model has not yet been implemented: {0}".format(self.attrs['external_model'])))

                        
        # either they are all one element or they are compatible lists no 1/2 way
        try:
            if len(self._Vpos) != len(self['Kp']) or \
                len(self._Vpos) != len(self['Epoch']) or \
                len(self['Epoch']) != len(self['Kp']):
                raise(ValueError('Inputs must be the same length, scalars or lists'))
        except TypeError:
            if isinstance(self._Vpos, list) and not isinstance(self['Kp'], list) \
                and not isinstance(self['Epoch'], list):
                    raise(ValueError('Inputs must be the same length, scalars or lists'))

        date = Lgm_CTrans.dateToDateLong(self['Epoch'])
        utc = Lgm_CTrans.dateToFPHours(self['Epoch'])
        Lgm_Set_Coord_Transforms( date, utc, self._mmi.c) # dont need pointer as it is one

                    
        ans = Lgm_Vector.Lgm_Vector(-1, -1, -1)

        retval = Lgm_TraceToEarth(pointer(self._Vpos),
                                ctypes.pointer(ans),
                                ctypes.c_double(TargetHeight),
                                ctypes.c_double(direction),
                                ctypes.c_double(tol),
                                pointer(self._mmi) )

        if retval == LGM_CLOSED:
            self['footpoint'] = ans
            self.attrs['retcode'] = retval
        elif retval == LGM_OPEN_IMF:
            self['footpoint'] = np.nan
            warnings.warn("LGM_OPEN_IMF")
            self.attrs['retcode'] = retval
        elif retval == LGM_OPEN_N_LOBE:
            self['footpoint'] = np.nan
            warnings.warn("LGM_OPEN_N_LOBE")
            self.attrs['retcode'] = retval
        elif retval == LGM_OPEN_S_LOBE:
            self['footpoint'] = np.nan
            warnings.warn("LGM_OPEN_S_LOBE")
            self.attrs['retcode'] = retval
        elif retval == LGM_INSIDE_EARTH:
            self['footpoint'] = np.nan
            warnings.warn("LGM_INSIDE_EARTH")
            self.attrs['retcode'] = retval
        elif retval == LGM_TARGET_HEIGHT_UNREACHABLE:
            self['footpoint'] = np.nan
            warnings.warn("LGM_TARGET_HEIGHT_UNREACHABLE")
            self.attrs['retcode'] = retval
        elif retval == LGM_BAD_TRACE:
            self['footpoint'] = np.nan
            warnings.warn("LGM_BAD_TRACE")
            self.attrs['retcode'] = retval
         


def TraceToEarth(pos, time,
                 MAGMODEL_ARGS=None,
                 coord_system = 'GSM',
                 INTERNAL_MODEL='LGM_IGRF',
                 EXTERNAL_MODEL='LGM_EXTMODEL_T89',
                 TargetHeight=120.,
                 tol=1e-7
                 ):
    """
    Easy wrapper to Lgm_TraceToEarth_py that returns both footpoints


    Parameters
    ----------
    pos : list
        3-element list of the position to do the calculation
    time : datetime
        date and time of when to do the calculation
    coord_system : str, optional
        the coordinate system of the position, default=GSM
    INTERNAL_MODEL : str, optional
        the internal magnetic field model to use, default=LGM_IGRF
    EXTERNAL_MODEL : str, optional
        the external magnetic field model to use, default=LGM_T89
    MAGMODEL_ARGS : dict, optional
        kwarguments to the magnetic field model, example for T89 {'Kp':5}
    TargetHeight : float, optional
        target height for the tracing, default=120km
    tol : float, optional
        tolerance for tracing, default=1e-7

    Returns
    -------
    out : tuple
        north and south footpoint 
    

    """
    a = Lgm_TraceToEarth_py(pos, time,
                direction='NORTH',
                coord_system = 'GSM',
                INTERNAL_MODEL=INTERNAL_MODEL,
                EXTERNAL_MODEL=EXTERNAL_MODEL,
                MAGMODEL_ARGS=MAGMODEL_ARGS,
                TargetHeight=TargetHeight, tol=tol)
    
    b = Lgm_TraceToEarth_py(pos, time,
                direction='SOUTH',
                coord_system = coord_system,
                INTERNAL_MODEL=INTERNAL_MODEL,
                EXTERNAL_MODEL=EXTERNAL_MODEL,
                MAGMODEL_ARGS=MAGMODEL_ARGS,
                TargetHeight=TargetHeight, tol=tol)
    return a['footpoint'].tolist(), b['footpoint'].tolist()
