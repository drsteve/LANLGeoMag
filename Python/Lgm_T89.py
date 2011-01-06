"""
Python implementation of the LanlGeoMag T89 Magnetic field model


@author: Brian Larsen
@organization: LANL
@contact: balarsen@lanl.gov

@version: V1: 23-Dec-2010 (BAL)
"""

class Lgm_T89(object):
    """
    Python implementation of the LanlGeoMag T89 Magnetic field model

    @author: Brian Larsen
    @organization: LANL
    @contact: balarsen@lanl.gov

    @version: V1: 23-Dec-2010 (BAL)
    """
    pass



import ctypes
import itertools

import numpy as np
import spacepy.toolbox as tb
from pylab import griddata, pcolor, pcolormesh, gca, draw

import Lgm
from Lgm_Types import LgmInt
from _Lgm import lib
import Lgm_Vector
import Lgm_CTrans
import _Lgm_MagModelInfo


from pylab import *

pos = Lgm_Vector.Lgm_Vector(-6.6, 0.1, 0)

B = Lgm_Vector.Lgm_Vector()
mmi = _Lgm_MagModelInfo.Lgm_MagModelInfo()

## up good

#c = Lgm_CTrans.Lgm_CTrans()
#mmi.c = ctypes.pointer(c)


Date = 20050831
UTC  = 9.0


lib.Lgm_Set_Coord_Transforms( Date, UTC, mmi.c);

#mmi.Bfield = ctypes.pointer((ctypes.CFUNCTYPE(LgmInt))(lib.Lgm_B_T89))

pos = Lgm_Vector.Lgm_Vector(-6.6, 0.0, 0.0)

for val in [0, 1, 2, 3, 4, 5]:
    mmi.Kp = val
    lib.Lgm_B_T89(pos, B, mmi)
    print(mmi.Kp, pos.x, pos.y, pos.z, B.x, B.y, B.z, B.magnitude())



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


ans = np.zeros(shape=[1001*1001, 6])
z = 0.
count = 0
for x in np.linspace(-10, 10, 1001):
    print('x:%f' % (x) )
    for y in np.linspace(-10, 10, 1001):
        #for z in np.linspace(-10, 10, 101):
        pos = Lgm_Vector.Lgm_Vector(x, y, z)
            #mmi = _Lgm_MagModelInfo.Lgm_MagModelInfo()
        bla = lib.Lgm_B_T89(pos, B, mmi)
        ans[count, :] = [x, y, z, B.x, B.y, B.z]
        count = count + 1

bx = ans[:, 3]
by = ans[:, 4]
bz = ans[:, 5]
x = ans[:, 0]
y = ans[:, 1]
z = ans[:, 2]

ind = np.negative(np.isnan(bx) + np.isnan(by) + np.isnan(bz))

bx = bx[ind]
by = by[ind]
bz = bz[ind]
x = x[ind]
y = y[ind]
z = z[ind]

ind = np.sqrt(x**2 + y**2 + z**2) > 1
bx = bx[ind]
by = by[ind]
bz = bz[ind]
x = x[ind]
y = y[ind]
z = z[ind]

ind = ( np.sqrt(bx**2+ by**2 + bz**2) < 1000 )

bx = bx[ind]
by = by[ind]
bz = bz[ind]
x = x[ind]
y = y[ind]
z = z[ind]




z0ind = np.where(z==0)[0]

xb = np.array(tb.bin_center_to_edges(np.linspace(-10, 10, 101)))
yb = np.array(tb.bin_center_to_edges(np.linspace(-10, 10, 101)))

gd = griddata(x[z0ind], y[z0ind], np.log10(np.sqrt(bx[z0ind]**2+ by[z0ind]**2 + bz[z0ind]**2)), \
              np.linspace(-10, 10, 101), np.linspace(-10, 10, 101))
pcolormesh(xb, yb, gd)
ax = gca()
ax.set_xlabel('X (Re)')
ax.set_ylabel('Y (Re)')
cb = colorbar()
cb.set_label('Bmag (nT)')
draw()
