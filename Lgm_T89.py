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
import numpy as np
import spacepy.toolbox as tb
from pylab import griddata, pcolor, pcolormesh

#pos = Lgm_Vector.Lgm_Vector(-6.6, 1, 1)
pos = Lgm_Vector.Lgm_Vector(-6.6, 0.1, 0)
#pos = Lgm_Vector.Lgm_Vector(-6.6, 1, 1)

B = Lgm_Vector.Lgm_Vector()
mmi = _Lgm_MagModelInfo.Lgm_MagModelInfo()

lib.Lgm_B_T89(pos, B, mmi)
print(B.x, B.y, B.z)

ans = []
for x in np.linspace(-10, 10, 101):
    print('x:%f' % (x) )
    for y in np.linspace(-10, 10, 101):
        for z in np.linspace(-10, 10, 101):
            pos = Lgm_Vector.Lgm_Vector(x, y, z)
            mmi = _Lgm_MagModelInfo.Lgm_MagModelInfo()
            lib.Lgm_B_T89(pos, B, mmi)
            ans.append( [x, y, z, B.x, B.y, B.z] )
bx = np.array(zip(*ans)[3])
by = np.array(zip(*ans)[4])
bz = np.array(zip(*ans)[5])
x = np.array(zip(*ans)[0])
y = np.array(zip(*ans)[1])
z = np.array(zip(*ans)[2])

indx = [ not np.isnan(val) for val in bx]
indy = [ not np.isnan(val) for val in by]
indz = [ not np.isnan(val) for val in bz]
ind = np.array([v1 and v2 and v2 for v1, v2, v3 in zip(indx, indy, indz)])

bx = bx[ind]
by = by[ind]
bz = bz[ind]
x = x[ind]
y = y[ind]
z = z[ind]

gd = griddata(x, y, np.sqrt(bx**2+ by**2 + bz**2), np.linspace(-10, 10, 101), np.linspace(-10, 10, 101))
pcolormesh(gd)
