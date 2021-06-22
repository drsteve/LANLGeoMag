"""
Overview
--------
Lgm_QinDenton module, this contains the necessary code QinDenton wrapper code

Authors
-------
Brian Larsen - LANL
"""
from ctypes import pointer

from .Lgm_Wrap import Lgm_QinDenton, Lgm_init_QinDentonDefaults, Lgm_destroy_QinDenton_children


class Lgm_QinDenton(Lgm_QinDenton):
    """
    LanlGeoMag class for holding QinDenton info needed by models, low level
    wrapper for C.  Calls the library for the default values and memory allocation

    Returns
    =======
    out : Lgm_QinDenton class
    """
    def __init__(self, verbose=False):
        """
        Set the default values from Lgm
        """
        Lgm_init_QinDentonDefaults(pointer(self), verbose)

    def __del__(self):
        """
        clean up memory that is allocated in C
        """
        Lgm_destroy_QinDenton_children(pointer(self))