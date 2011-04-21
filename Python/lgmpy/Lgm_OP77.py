# -*- coding: utf-8 -*-
"""
Python implementation of the LanlGeoMag OP77 Magnetic field model


@author: Brian Larsen
@organization: LANL
@contact: balarsen@lanl.gov

@version: V1: 21-Mar-2011 (BAL)
"""
import datetime
from ctypes import pointer

import numpy as np

from Lgm_Wrap import LGM_CDIP, LGM_EDIP, LGM_IGRF, Lgm_Set_Coord_Transforms, \
    Lgm_B_OP77

import MagData
import Lgm_Vector
import Lgm_CTrans
import Lgm_MagModelInfo

class Lgm_OP77(MagData.MagData):
    """
    Python implementation of the LanlGeoMag OP77 Magnetic field model

    @author: Brian Larsen
    @organization: LANL
    @contact: balarsen@lanl.gov

    @version: V1: 21-Mar-2011 (BAL)
    """
    def __init__(self, pos, time, coord_system = 'GSM', INTERNAL_MODEL='LGM_IGRF',):
        super(Lgm_OP77, self).__init__(Position=pos, Epoch=time, coord_system = coord_system, INTERNAL_MODEL=INTERNAL_MODEL,)

        # pos must be an Lgm_Vector or list or sensible ndarray
        try:
            self._Vpos = self._pos2Lgm_Vector(pos)
            assert self._Vpos
        except:
            raise(TypeError('pos must be a Lgm_Vector or list of Lgm_vectors') )

        # time must be a datetime
        if not isinstance(time, datetime.datetime) and \
            not isinstance(time, list):
            raise(TypeError('time must be a datetime or list of datetime') )

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

        # either they are all one elemet or they are compatible lists no 1/2 way
        try:
            if len(self._Vpos) != len(self['Epoch']):
                raise(ValueError('Inputs must be the same length, scalars or lists'))
        except TypeError:
            if isinstance(self._Vpos, list) and not isinstance(self['Epoch'], list):
                raise(ValueError('Inputs must be the same length, scalars or lists'))

        self['B'] = self.calc_B()

    def calc_B(self):
        try:
            ans = []
            for v1, v2 in zip(self._Vpos, self['Epoch'], ):
                date = Lgm_CTrans.dateToDateLong(v2)
                utc = Lgm_CTrans.dateToFPHours(v2)
                Lgm_Set_Coord_Transforms( date, utc, self._mmi.c) # dont need pointer as it is one
                B = Lgm_Vector.Lgm_Vector()
                retval = Lgm_B_OP77(pointer(v1), pointer(B), pointer(self._mmi))
                if retval != 1:
                    raise(RuntimeWarning('Odd return from OP77.c') )
                ans.append(B)
            return ans
        except TypeError:
            date = Lgm_CTrans.dateToDateLong(self['Epoch'])
            utc = Lgm_CTrans.dateToFPHours(self['Epoch'])
            Lgm_Set_Coord_Transforms( date, utc, self._mmi.c) # dont need pointer as it is one
            B = Lgm_Vector.Lgm_Vector()
            retval = Lgm_B_OP77(pointer(self._Vpos), pointer(B), pointer(self._mmi) )
            if retval != 1:
                raise(RuntimeWarning('Odd return from OP77.c') )
            return B

    def _pos2Lgm_Vector(self, pos):
        if isinstance(pos, Lgm_Vector.Lgm_Vector):
            return pos
        if isinstance(pos, np.ndarray):
            pos = pos.tolist()
        if isinstance(pos, list):
            Vpos = []
            for val in pos:
                if not isinstance(val, list):
                    return Lgm_Vector.Lgm_Vector(pos[0], pos[1], pos[2])
                Vpos.append(Lgm_Vector.Lgm_Vector(val[0], val[1], val[2]))
            return Vpos

def OP77(pos, time, coord_system = 'GSM', INTERNAL_MODEL='LGM_IGRF',):
    """
    Easy wrapper to just return values without having to create an instance of
    Lgm_OP77

    @param pos: a list of 3 element lists of positions in coord_system system
    @type pos: list
    @param time: a datetime or list of datetime objects
    @type time: (list, datetime)

    @keyword coord_system: the name of the coord system to use (or Lgm number)
    @type coord_system: (str, int)
    @keyword INTERNAL_MODEL: the intermal magnetic field model to use (or Lgm number)
    @type INTERNAL_MODEL: (str, int)

    @retval: list of the x, y, z, compents of B in nT
    @rtype: list

    @author: Brian Larsen
    @organization: LANL
    @contact: balarsen@lanl.gov

    @version: V1: 21-Mar-2011 (BAL)
    """
    a = Lgm_OP77(pos, time, coord_system = coord_system,
                        INTERNAL_MODEL=INTERNAL_MODEL)
    try:
        ans = [val.tolist() for val in a['B']]
        return ans
    except TypeError:
        return a['B'].tolist()
