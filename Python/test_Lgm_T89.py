#!/usr/bin/env python

"""
Test suite for the Lgm_T89 file

@author: Brian Larsen
@organization: LANL
@contact: balarsen@lanl.gov

@version: V1: 20-Dec-2010 (BAL)
"""

import unittest
import Lgm
import Lgm_CTrans
import Lgm_Vector
import Lgm_T89
from _Lgm import lib

class Lgm_T89Tests(unittest.TestCase):
    """
    Tests related to Lgm_T89
    @author: Brian Larsen
    @organization: LANL
    @contact: balarsen@lanl.gov

    @version: V1: 04-Jan-2011 (BAL)
    """

    def setUp(self):
        super(Lgm_CTransTests, self).setUp()

    def tearDown(self):
        super(Lgm_CTransTests, self).tearDown()



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
#
#


if __name__ == '__main__':
    unittest.main()
