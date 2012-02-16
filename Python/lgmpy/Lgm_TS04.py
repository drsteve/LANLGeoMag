"""
Overview
--------
Python implementation of the LanlGeoMag T89 Magnetic field model

"""
__author__ = 'Brian Larsen, Mike Henderson - LANL'

import datetime
from ctypes import pointer
import ctypes

import numpy as np

import MagData
from Lgm_Wrap import LGM_CDIP, LGM_EDIP, LGM_IGRF, Lgm_Set_Coord_Transforms, Lgm_B_TS04, \
                    Lgm_Set_Lgm_B_cdip_InternalModel, Lgm_Set_Lgm_B_edip_InternalModel, Lgm_Set_Lgm_B_IGRF_InternalModel, \
                    Lgm_get_QinDenton_at_JD, Lgm_Date_to_JD, Lgm_B_TS04_opt

import lgmpy
import Lgm_Vector
import Lgm_CTrans
import Lgm_MagModelInfo


class Lgm_TS04_QD(MagData.MagData):
    #int Lgm_B_TS04_opt( Lgm_Vector *v, Lgm_Vector *B, Lgm_MagModelInfo *Info ) {
    def __init__(self, pos, time, coord_system = 'GSM', INTERNAL_MODEL='LGM_IGRF', verbose=False):
        super(Lgm_TS04_QD, self).__init__(Position=pos, Epoch=time)

        self.verbose = verbose

        if not isinstance(pos, (Lgm_Vector.Lgm_Vector, list, np.ndarray)):
            raise(TypeError('pos must be a Lgm_Vector or list of Lgm_vectors') )
        self._Vpos = self._pos2Lgm_Vector(pos)

        # time must be a datetime
        if not isinstance(time, (datetime.datetime, list)):
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

        # and actually set the internal model in Lgm
        if self.attrs['internal_model'] == LGM_CDIP:
            Lgm_Set_Lgm_B_cdip_InternalModel(pointer(self._mmi))
        elif self.attrs['internal_model'] == LGM_EDIP:
            Lgm_Set_Lgm_B_edip_InternalModel(pointer(self._mmi))
        elif self.attrs['internal_model'] == LGM_IGRF:
            Lgm_Set_Lgm_B_IGRF_InternalModel(pointer(self._mmi))

        #self.data = T89_Data(pos, time, Kp, coord_system, INTERNAL_MODEL)
        self['B'] = self.calc_B()

    def calc_B(self):
        ans = []
#            for v1, v2, v3 in zip(self._Vpos, self['Epoch'], self['Kp']):
        date = Lgm_CTrans.dateToDateLong(self['Epoch'])
        utc = Lgm_CTrans.dateToFPHours(self['Epoch'])
        Lgm_Set_Coord_Transforms( date, utc, self._mmi.c) # dont need pointer as it is one
        B = Lgm_Vector.Lgm_Vector()

        # Grab the QinDenton data
        # Lgm_get_QinDenton_at_JD( JD, &p, 1 );
        # JD = Lgm_Date_to_JD( Date, UTC, mInfo->c );
        JD = Lgm_Date_to_JD(date, utc, pointer(self._mmi.c))
        qd_one = lgmpy.Lgm_Wrap.Lgm_QinDentonOne()
        Lgm_get_QinDenton_at_JD( JD, pointer(qd_one), self.verbose)
        lgmpy.Lgm_Wrap.Lgm_set_QinDenton(pointer(qd_one), pointer(self._mmi.c))

        #int Lgm_B_TS04_opt( Lgm_Vector *v, Lgm_Vector *B, Lgm_MagModelInfo *Info ) {
        retval = Lgm_B_TS04(pointer(self._Vpos), pointer(B), pointer(self._mmi))
#        retval = Lgm_B_TS04_opt(pointer(self._Vpos), pointer(B), pointer(self._mmi))
        if retval != 1:
            raise(RuntimeWarning('Odd return from Lgm_TS04') )
        ans.append(B)
        return ans

    def _pos2Lgm_Vector(self, pos):
        if isinstance(pos, Lgm_Vector.Lgm_Vector):
            return pos
        if isinstance(pos, list):
            Vpos = []
            for val in pos:
                if not isinstance(val, list):
                    return Lgm_Vector.Lgm_Vector(pos[0], pos[1], pos[2])
                Vpos.append(Lgm_Vector.Lgm_Vector(val[0], val[1], val[2]))
            return Vpos
        if isinstance(pos, np.ndarray):
            return(self._pos2Lgm_Vector(pos.tolist())) 


    def _pos2Lgm_Vector(self, pos):
        if isinstance(pos, Lgm_Vector.Lgm_Vector):
            return pos
        if isinstance(pos, list):
            Vpos = []
            for val in pos:
                if not isinstance(val, list):
                    return Lgm_Vector.Lgm_Vector(pos[0], pos[1], pos[2])
                Vpos.append(Lgm_Vector.Lgm_Vector(val[0], val[1], val[2]))
            return Vpos
        if isinstance(pos, np.ndarray):
            return(self._pos2Lgm_Vector(pos.tolist())) 
                        
            

class Lgm_TS04(MagData.MagData):
    """
    Python implementation of the LanlGeoMag TS04 Magnetic field model.  This is
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
    out : Lgm_T89 object
        Lgm_T89 object will the magnetic field value and other information

    See Also
    --------
    Lgm_T89.T89

    """

    #int Lgm_B_TS04_opt( Lgm_Vector *v, Lgm_Vector *B, Lgm_MagModelInfo *Info ) {
    def __init__(self, pos, time, P, Dst, By, Bz, W, coord_system = 'GSM', INTERNAL_MODEL='LGM_IGRF',):
        #            parmod[1]  = Info->P; 	// Pressure in nPa
        #            parmod[2]  = Info->Dst;     // Dst in nPa
        #            parmod[3]  = Info->By; 	// IMF By in nT
        #            parmod[4]  = Info->Bz; 	// IMF Bz in nT
        #            parmod[5]  = Info->W[0];   // W1
        #            parmod[6]  = Info->W[1];   // W2
        #            parmod[7]  = Info->W[2];   // W3
        #            parmod[8]  = Info->W[3];   // W4
        #            parmod[9]  = Info->W[4];   // W5
        #            parmod[10] = Info->W[5];   // W6
        # pos must be a Lgm_Vector or list of Lgm_Vectors
        super(Lgm_TS04, self).__init__(Position=pos, Epoch=time, P=P, Dst=Dst, By=By, Bz=Bz, W=W, coord_system = coord_system, INTERNAL_MODEL=INTERNAL_MODEL,)

        if not isinstance(pos, Lgm_Vector.Lgm_Vector) and \
            not isinstance(pos, list):
            raise(TypeError('pos must be a Lgm_Vector or list of Lgm_vectors') )
        self._Vpos = self._pos2Lgm_Vector(pos)

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

        # and actually set the internal model in Lgm
        if self.attrs['internal_model'] == LGM_CDIP:
            Lgm_Set_Lgm_B_cdip_InternalModel(pointer(self._mmi))
        elif self.attrs['internal_model'] == LGM_EDIP:
            Lgm_Set_Lgm_B_edip_InternalModel(pointer(self._mmi))
        elif self.attrs['internal_model'] == LGM_IGRF:
            Lgm_Set_Lgm_B_IGRF_InternalModel(pointer(self._mmi))

        #self.data = T89_Data(pos, time, Kp, coord_system, INTERNAL_MODEL)
        self['B'] = self.calc_B()

    def calc_B(self):
        ans = []
#            for v1, v2, v3 in zip(self._Vpos, self['Epoch'], self['Kp']):
        date = Lgm_CTrans.dateToDateLong(self['Epoch'])
        utc = Lgm_CTrans.dateToFPHours(self['Epoch'])
        Lgm_Set_Coord_Transforms( date, utc, self._mmi.c) # dont need pointer as it is one
        B = Lgm_Vector.Lgm_Vector()
        self._mmi.Dst = self['Dst']
        self._mmi.By = self['By']
        self._mmi.Bz = self['Bz']
        arr_type = ctypes.c_double*6
        self._mmi.W = arr_type(*self['W'])
        
        #int Lgm_B_TS04_opt( Lgm_Vector *v, Lgm_Vector *B, Lgm_MagModelInfo *Info ) {
        retval = Lgm_B_TS04(pointer(self._Vpos), pointer(B), pointer(self._mmi))
        if retval != 1:
            raise(RuntimeWarning('Odd return from Lgm_T89') )
        ans.append(B)
        return ans

    def _pos2Lgm_Vector(self, pos):
        if isinstance(pos, Lgm_Vector.Lgm_Vector):
            return pos
        if isinstance(pos, list):
            Vpos = []
            for val in pos:
                if not isinstance(val, list):
                    return Lgm_Vector.Lgm_Vector(pos[0], pos[1], pos[2])
                Vpos.append(Lgm_Vector.Lgm_Vector(val[0], val[1], val[2]))
            return Vpos
        if isinstance(pos, np.ndarray):
            raise(NotImplementedError('Only lists can be input for position now') )

def TS04(pos, time, P, Dst, By, Bz, W, coord_system = 'GSM', INTERNAL_MODEL='LGM_IGRF',):
    """
    Easy wrapper to just return values without having to create an instance of
    Lgm_TS04.

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
    >>> from lgmpy import Lgm_T89
    >>> import datetime
    >>> Lgm_T89.T89([1,2,3], datetime.datetime(1999, 1, 16, 12, 34, 12), 3)
    [-509.297442..., -619.8167352..., -426.193633024...]


    See Also
    --------
    Lgm_T89.Lgm_T89

    """
    a = Lgm_TS04(pos, time, P, Dst, By, Bz, W, coord_system = coord_system,
                        INTERNAL_MODEL=INTERNAL_MODEL)
    try:
        ans = [val.tolist() for val in a['B']]
        return ans
    except TypeError:
        return a['B'].tolist()
