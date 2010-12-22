"""
Vector class for Lgm

@author: Brian Larsen
@organization: LANL
@contact: balarsen@lanl.gov

@version: V1: 20-Dec-2010 (BAL)
"""


import ctypes
from _Lgm_Types import *


class Lgm_Vector(ctypes.Structure):
    _fields_ = [ ( "x", LgmDouble ),
        ("y", LgmDouble),
        ("z", LgmDouble) ]

Lgm_VectorP = ctypes.POINTER(Lgm_Vector)
