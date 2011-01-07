"""
Python implementation of the LanlGeoMag T89 Magnetic field model


@author: Brian Larsen
@organization: LANL
@contact: balarsen@lanl.gov

@version: V1: 23-Dec-2010 (BAL)
"""

import ctypes
import itertools
import datetime

import numpy as np
import spacepy.toolbox as tb
from pylab import griddata, pcolor, pcolormesh, gca, draw

import Lgm
from Lgm_Types import LgmInt
from _Lgm import lib
import Lgm_Vector
import Lgm_CTrans
import Lgm_MagModelInfo
import Lgm_DateAndTime

class Lgm_T89(object):
    """
    Python implementation of the LanlGeoMag T89 Magnetic field model

    @author: Brian Larsen
    @organization: LANL
    @contact: balarsen@lanl.gov

    @version: V1: 23-Dec-2010 (BAL)
    """
    def __init__(self, pos, time, Kp, coord_system = 'GSM', INTERNAL_MODEL='LGM_IGRF',):
        # pos must be a Lgm_Vector or list of Lgm_Vectors
        if not isinstance(pos, Lgm_Vector.Lgm_Vector) and \
            not isinstance(pos, list) and \
            not isinstance(pos, np.ndarray):
            raise(TypeError('pos must be a Lgm_Vector or (list, numpy.ndarray) of Lgm_vectors') )
        self.position = pos
        self._Vpos = self._pos2Lgm_Vector(pos)

        # time must be a datetime
        if not isinstance(time, datetime.datetime) and \
            not isinstance(time, list) and \
            not isinstance(time, np.ndarray):
            raise(TypeError('time must be a datetime or (list, numpy.ndarray) of datetime') )
        self.time = time

        #TODO implement this and remove this bit
        try:
            assert(len(pos))
            assert(len(time))
            raise(NotImplementedError('Multi element input not yet written') )
        except:
            pass

        if Kp < 0 or Kp > 5:
            raise(ValueError('T89 is only defined for integer Kp from 0 to 5') )
        self.Kp = Kp

        if INTERNAL_MODEL not in (Lgm_MagModelInfo.LGM_CDIP,
                                  Lgm_MagModelInfo.LGM_EDIP,
                                  Lgm_MagModelInfo.LGM_IGRF) and \
            INTERNAL_MODEL not in ('LGM_CDIP',
                                  'LGM_EDIP',
                                  'LGM_IGRF'):
            raise(ValueError('INTERNAL_MODEL must be LGM_CDIP, LGM_EDIP, or LGM_IGRF') )
        if isinstance(INTERNAL_MODEL, str):
            self.INTERNAL_MODEL = Lgm_MagModelInfo.__getattribute__(INTERNAL_MODEL)
        else:
            self.INTERNAL_MODEL = INTERNAL_MODEL

        if coord_system != 'GSM':
            raise(NotImplementedError('Different coord systems are not yet ready to use') )
        self.coord_system = coord_system

        self._mmi = Lgm_MagModelInfo.Lgm_MagModelInfo()

    def calc_B(self):
        date = Lgm_DateAndTime.dateToDateLong(self.time)
        utc = Lgm_DateAndTime.dateToFPHours(self.time)
        lib.Lgm_Set_Coord_Transforms( date, utc, self._mmi.c);
        self.B = Lgm_Vector.Lgm_Vector()
        self._mmi.Kp = self.Kp
        retval = lib.Lgm_B_T89(self._Vpos, self.B, self._mmi)
        if retval != 1:
            raise(RuntimeWarning('Odd return from Lgm_T89') )
        return self.B

    def _pos2Lgm_Vector(self, pos):
        if isinstance(pos, Lgm_Vector.Lgm_Vector):
            return pos
        if isinstance(pos, list):
            Vpos = []
            for val in pos:
                if not isinstance(val, list):
                    return Lgm_Vector.Lgm_Vector(pos[0], pos[1], pos[2])
                Vpos.append(Lgm_Vector.Lgm_Vector(val[0], val[1], val[2]))
        if isinstance(pos, np.ndarray):
            raise(NotImplementedError('Only lists can be imput for position now') )


##################################################
#
#Kp      Ux (Re)      Uy (Re)      Uz (Re)      Bx (nT)      By (nT)      Bz (nT)   Bmag (nT)
# 0         -6.6            0            0     -18.9764      -1.8654      80.3931      82.6235
# 1         -6.6            0            0     -20.8334      -1.8654      74.6219       77.498
# 2         -6.6            0            0      -22.563      -1.8654      70.5243      74.0692
# 3         -6.6            0            0     -26.4168      -1.8654      64.6448       69.859
# 4         -6.6            0            0     -32.1658      -1.8654      60.0782      68.1726
# 5         -6.6            0            0      -45.384      -1.8654      49.3629      67.0812
# 3           -1            0            0     -6244.66     -792.844      31094.6      31725.3
# 3         -1.5            0            0     -1630.83     -386.285      9019.71      9174.09
# 3           -2            0            0     -641.802     -157.771      3731.95      3790.02
# 3         -2.5            0            0     -317.069     -72.7738      1865.09      1893.25
# 3           -3            0            0     -181.087     -37.5704      1045.65      1061.88
# 3         -3.5            0            0     -114.734      -21.193      632.635      643.304
# 3           -4            0            0     -78.8141     -12.8114       403.97      411.785
# 3         -4.5            0            0     -57.8633     -8.18183      268.621      274.905
# 3           -5            0            0     -44.9333     -5.46261      184.424      189.897
# 3         -5.5            0            0     -36.5814     -3.78302       129.97      135.073
# 3           -6            0            0      -30.969     -2.70126      93.6087      98.6355
# 3         -6.5            0            0      -27.055     -1.97955      68.6516      73.8169
# 3           -7            0            0     -24.2235     -1.48334      51.1005      56.5706
##################################################
