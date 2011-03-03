from __future__ import division

"""
perform tracing to see if a file line is closed
"""
from ctypes import pointer
import math

from Lgm_Wrap import Lgm_Trace, LGM_OPEN_IMF, LGM_CLOSED, LGM_OPEN_N_LOBE
from Lgm_Wrap import LGM_OPEN_S_LOBE, LGM_INSIDE_EARTH, LGM_TARGET_HEIGHT_UNREACHABLE
from Lgm_Wrap import Lgm_B_OP77, Lgm_Set_Lgm_B_OP77, Lgm_Set_Coord_Transforms
from Lgm_Wrap import Lgm_Convert_Coords, GSM_TO_SM, WGS84_A
import Lgm_Vector
import Lgm_MagModelInfo
import Lgm_CTrans

def _simpleL(position, MagModelInfo):
    """
    code to return a simple measure of L
    Inputs:
    - GSM position of nothern most point (from Closed_Field)
    """
    # Get a simple measure of how big L is
    position_sm = Lgm_Vector.Lgm_Vector()
    Lgm_Convert_Coords( pointer(position), pointer(position_sm), GSM_TO_SM,
                       MagModelInfo.c );
    Lam = math.asin(position_sm.z/position_sm.magnitude())
    CosLam = math.cos(Lam)
    LSimple = (1.0+120.0/WGS84_A)/( CosLam*CosLam )
    return LSimple


# I belive the position is in GSM, check with Mike
def Closed_Field(position, date,
                 height=120,
                 tol1=0.01, tol2=1e-7,
                 bfield='Lgm_B_OP77', coord_system='GSM',
                 extended_out = False,
                 **params):
    # input checking
    if coord_system != 'GSM':
        raise(NotImplementedError('Different coord systems are not yet ready to use') )
    # could consider a Lgm_MagModelInfo param to use an existing one
    mmi = Lgm_MagModelInfo.Lgm_MagModelInfo()
    if bfield == 'Lgm_B_OP77':
        Lgm_Set_Lgm_B_OP77( pointer(mmi) )
    else:
        raise(NotImplementedError('Only Lgm_B_OP77 implented so far'))

    datelong = Lgm_CTrans.dateToDateLong(date)
    utc = Lgm_CTrans.dateToFPHours(date)
    Lgm_Set_Coord_Transforms( datelong, utc, mmi.c) # dont need pointer as it is one

    try:
        position = Lgm_Vector.Lgm_Vector(*position)
    except TypeError:
        raise(TypeError('position must be an iterable') )
    northern = Lgm_Vector.Lgm_Vector()
    southern = Lgm_Vector.Lgm_Vector()
    minB = Lgm_Vector.Lgm_Vector()

    ans = Lgm_Trace(pointer(position),
                    pointer(northern),
                    pointer(southern),
                    pointer(minB), height, tol1, tol2, pointer(mmi) )
    if ans == LGM_OPEN_IMF:
        retstr = 'LGM_OPEN_IMF'
    elif ans == LGM_CLOSED:
        retstr = 'LGM_CLOSED'
    elif ans == LGM_OPEN_N_LOBE:
        retstr = 'LGM_OPEN_N_LOBE'
    elif ans == LGM_OPEN_S_LOBE:
        retstr = 'LGM_OPEN_S_LOBE'
    elif ans == LGM_INSIDE_EARTH:
        retstr = 'LGM_INSIDE_EARTH'
    elif ans == LGM_TARGET_HEIGHT_UNREACHABLE:
        retstr = 'LGM_TARGET_HEIGHT_UNREACHABLE'
    if extended_out:
        L = _simpleL(northern, mmi)
        return retstr, northern.tolist(), southern.tolist(), minB.tolist(), L
    else:
        return retstr
