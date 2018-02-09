"""
Overview
--------
Python implementation of the LanlGeoMag T89 Magnetic field model

"""
__author__ = 'Brian Larsen, Mike Henderson - LANL'

import datetime
from ctypes import pointer

import numpy as np

import MagData
from Lgm_Wrap import LGM_CDIP, LGM_EDIP, LGM_IGRF, Lgm_Set_Coord_Transforms, \
    Lgm_B_T89c, Lgm_Set_Lgm_B_cdip_InternalModel, Lgm_Set_Lgm_B_edip_InternalModel, Lgm_Set_Lgm_B_IGRF_InternalModel
import Lgm_Vector
import Lgm_CTrans
import Lgm_MagModelInfo
from utils import pos2Lgm_Vector

class Lgm_T89c(MagData.MagData):
    """
    Python implementation of the LanlGeoMag T89c Magnetic field model.  This is
    the full wrapper, most users will not use this.

    Parameters
    ----------
    pos : list
        3-element list of the position to do the calculation
    time : datetime
        date and time of when to do the calculation
    Kp : int
        Kp value for the calculation
    coord_system : str, optional
        the coordinate system of the position, default=GSM
    INTERNAL_MODEL : str, optional
        the internal magnetic field model to use, default=LGM_IGRF

    Returns
    -------
    out : Lgm_T89c object
        Lgm_T89c object will the magnetic field value and other information

    See Also
    --------
    Lgm_T89c.T89c

    """
    def __init__(self, pos, time, Kp, coord_system = 'GSM', INTERNAL_MODEL='LGM_IGRF',):
        # pos must be a Lgm_Vector or list of Lgm_Vectors
        super(Lgm_T89c, self).__init__(Position=pos, Epoch=time, Kp=Kp, coord_system = coord_system, INTERNAL_MODEL=INTERNAL_MODEL,)

        if not isinstance(pos, Lgm_Vector.Lgm_Vector) and \
            not isinstance(pos, list):
            raise(TypeError('pos must be a Lgm_Vector or list of Lgm_vectors') )
        self._Vpos = pos2Lgm_Vector(pos)

        # time must be a datetime
        if not isinstance(time, datetime.datetime) and \
            not isinstance(time, list):
            raise(TypeError('time must be a datetime or list of datetime') )

        try:
            for val in Kp:
                if val < 0 or val > 5:
                    raise(ValueError('T89c is only defined for integer Kp from 0 to 5') )
        except TypeError:
            if Kp < 0 or Kp > 5:
                raise(ValueError('T89c is only defined for integer Kp from 0 to 5') )

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

        #self.data = T89_Data(pos, time, Kp, coord_system, INTERNAL_MODEL)
        self['B'] = self.calc_B()
        #self['B'].attrs['Units'] = 'nT'
        #self['B'].attrs['System'] = 'GSM'

    def calc_B(self):
        try:
            ans = []
            for v1, v2, v3 in zip(self._Vpos, self['Epoch'], self['Kp']):
                date = Lgm_CTrans.dateToDateLong(v2)
                utc = Lgm_CTrans.dateToFPHours(v2)
                Lgm_Set_Coord_Transforms( date, utc, self._mmi.c) # dont need pointer as it is one
                B = Lgm_Vector.Lgm_Vector()
                self._mmi.Kp = v3
                retval = Lgm_B_T89c(pointer(v1), pointer(B), pointer(self._mmi))
                if retval != 1:
                    raise(RuntimeWarning('Odd return from Lgm_T89c') )
                ans.append(B)
            return ans
        except TypeError:
            date = Lgm_CTrans.dateToDateLong(self['Epoch'])
            utc = Lgm_CTrans.dateToFPHours(self['Epoch'])
            Lgm_Set_Coord_Transforms( date, utc, self._mmi.c) # dont need pointer as it is one
            B = Lgm_Vector.Lgm_Vector()
            self._mmi.Kp = self['Kp']
            retval = Lgm_B_T89c(pointer(self._Vpos), pointer(B), pointer(self._mmi) )
            if retval != 1:
                raise(RuntimeWarning('Odd return from Lgm_T89c') )
            return B


def T89c(pos, time, Kp, coord_system = 'GSM', INTERNAL_MODEL='LGM_IGRF',):
    """
    Easy wrapper to just return values without having to create an instance of
    Lgm_T89c.

    All input parameters can be either their type or a list of the that type, all
    inputs must be the same length

    Parameters
    ----------
    pos : list
        3-element list of the position to do the calculation
    time : datetime
        date and time of when to do the calculation
    Kp : int
        Kp value for the calculation
    coord_system : str, optional
        the coordinate system of the position, default=GSM
    INTERNAL_MODEL : str, optional
        the internal magnetic field model to use, default=LGM_IGRF

    Returns
    -------
    out : list
        Magnetic field vector in GSM coordinates in nT

    Examples
    --------
    >>> from lgmpy import Lgm_T89c
    >>> import datetime
    >>> Lgm_T89c.T89c([1,2,3], datetime.datetime(1999, 1, 16, 12, 34, 12), 3)
    [-509.297442..., -619.8167352..., -426.193633024...]


    See Also
    --------
    Lgm_T89c.Lgm_T89c

    """
    a = Lgm_T89c(pos, time, Kp, coord_system = coord_system,
                        INTERNAL_MODEL=INTERNAL_MODEL)
    try:
        ans = [val.tolist() for val in a['B']]
        return ans
    except TypeError:
        return a['B'].tolist()
