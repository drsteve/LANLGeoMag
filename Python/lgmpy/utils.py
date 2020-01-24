from __future__ import division

"""
collection of utility routines uses around lgmpy
"""
import numpy as np

from . import Lgm_Vector


all = ['pos2Lgm_Vector']

def pos2Lgm_Vector(pos):
    if isinstance(pos, Lgm_Vector.Lgm_Vector):
        return pos
    if isinstance(pos, np.ndarray):
        pos = pos.tolist()
    if isinstance(pos, list):
        Vpos = []
        for val in pos:
            if not isinstance(val, list):
                return Lgm_Vector.Lgm_Vector(pos[0], pos[1], pos[2])
            Vpos.append(Lgm_Vector.Lgm_Vector(val[0], val[1], val[2]))
        return Vpos



