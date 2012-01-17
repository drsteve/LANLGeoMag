"""
Overview
--------
Lgm_MagModelInfo module, this contains the necessary code Mag Model Info

    Authors
    -------
    Brian Larsen - LANL
"""
from ctypes import pointer

from Lgm_Wrap import Lgm_MagModelInfo, Lgm_InitMagInfoDefaults, Lgm_FreeMagInfo_children


class Lgm_MagModelInfo(Lgm_MagModelInfo):
    """
    LanlGeoMag class for holding Mag Model info needed by all models, low level
    wrapper for C.  Calls the library for the default values

    Returns
    =======
    out : Lgm_MagModelInfo class
    """
    def __init__(self):
        """
        Set the default values from Lgm
        """
        Lgm_InitMagInfoDefaults(pointer(self))

    def __del__(self):
        """
        clean up memory that is allocated in C
        """
        Lgm_FreeMagInfo_children(pointer(self))