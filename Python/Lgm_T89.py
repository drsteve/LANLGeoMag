"""
Python implementation of the LanlGeoMag T89 Magnetic field model


@author: Brian Larsen
@organization: LANL
@contact: balarsen@lanl.gov

@version: V1: 23-Dec-2010 (BAL)
"""

import Lgm
from _Lgm import lib

class Lgm_T89(object):
    """
    Python implementation of the LanlGeoMag T89 Magnetic field model

    @author: Brian Larsen
    @organization: LANL
    @contact: balarsen@lanl.gov

    @version: V1: 23-Dec-2010 (BAL)
    """
    pass



import Lgm
from _Lgm import lib
import ctypes
import Lgm_Vector
import _Lgm_MagModelInfo

pos = Lgm_Vector.Lgm_Vector(-6.6, 0, 0)
B = Lgm_Vector.Lgm_Vector()
mmi = _Lgm_MagModelInfo.Lgm_MagModelInfo()

print(lib.Lgm_B_T89(pos, B, mmi))
print(B.x, B.y, B.z)
