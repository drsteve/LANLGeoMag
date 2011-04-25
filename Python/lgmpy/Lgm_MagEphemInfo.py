"""
Overview
--------
Lgm_MagModelInfo module, this contains the necessary code Mag Model Info

    Authors
    -------
    Brian Larsen - LANL
"""
from ctypes import pointer

from Lgm_Wrap import Lgm_MagEphemInfo, Lgm_InitMagEphemInfoDefaults


class Lgm_MagEphemInfo(Lgm_MagEphemInfo):
    """
    LanlGeoMag class for holding Mag Ephem info needed by all models, low level
    wrapper for C.  Calls the library for the default values

    Returns
    =======
    out : Lgm_MagEphemInfo class
    """
    def __init__(self, MaxPitchAngles, Verbosity=0):
        """
        Set the default values from Lgm
        """
        Lgm_InitMagEphemInfoDefaults(pointer(self), MaxPitchAngles, Verbosity)
