"""
Lgm_Eop module, this contains the necessary code

@author: Brian Larsen
@organization: LANL
@contact: balarsen@lanl.gov

@version: V1: 21-Dec-2010 (BAL)
"""

import ctypes
from _Lgm_Types import *

class Lgm_NgaEopp(ctypes.Structure):
    @classmethod
    def assign_fields(cls):
            cls._fields_ = [("ta", LgmDouble),
            ("A", LgmDouble),
            ("B", LgmDouble),
            ("C1", LgmDouble),
            ("C2", LgmDouble),
            ("D1", LgmDouble),
            ("D2", LgmDouble),
            ("P1", LgmDouble),
            ("P2", LgmDouble),
            ("E", LgmDouble),
            ("F", LgmDouble),
            ("G1", LgmDouble),
            ("G2", LgmDouble),
            ("H1", LgmDouble),
            ("H2", LgmDouble),
            ("Q1", LgmDouble),
            ("Q1", LgmDouble),
            ("Q2", LgmDouble),
            ("tb", LgmDouble),
            ("I", LgmDouble),
            ("J", LgmDouble),
            ("K1", LgmDouble),
            ("K2", LgmDouble),
            ("K3", LgmDouble),
            ("K4", LgmDouble),
            ("L1", LgmDouble),
            ("L2", LgmDouble),
            ("L3", LgmDouble),
            ("L4", LgmDouble),
            ("R1", LgmDouble),
            ("R2", LgmDouble),
            ("R3", LgmDouble),
            ("R4", LgmDouble),
            ("dat", LgmInt),
            ("EOPPWk", LgmInt),
            ("teff", LgmInt) ]

class Lgm_Eop(ctypes.Structure):
    @classmethod
    def assign_fields(cls):
            cls._fields_ = [("Size", LgmLong),
            ("nEopVals", LgmLong),
            ("Verbosity", LgmInt),
            ("Date", LgmLong),
            ("MJD", LgmDoubleP),
            ("xp", LgmDoubleP),
            ("yp", LgmDoubleP),
            ("DUT1", LgmDoubleP),
            ("LOD", LgmDoubleP),
            ("dPsi", LgmDoubleP),
            ("dEps", LgmDoubleP),
            ("dX", LgmDoubleP),
            ("dY", LgmDoubleP),
            ("DAT", LgmDoubleP) ]



class Lgm_EopOne(ctypes.Structure):
    @classmethod
    def assign_fields(cls):
        cls._fields_ = [("Date", LgmLong),
        ("JD", LgmDouble),
        ("MJD", LgmDouble),
        ("UTC", LgmDouble),
        ("xp", LgmDouble),
        ("yp", LgmDouble),
        ("DUT1", LgmDouble),
        ("LOD", LgmDouble),
        ("dPsi", LgmDouble),
        ("eDps", LgmDouble),
        ("dX", LgmDouble),
        ("dY", LgmDouble),
        ("DAT", LgmDouble) ]
