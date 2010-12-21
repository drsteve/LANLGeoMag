
import ctypes
from Lgm_Types import *
# from Lgm_CTrans import *



class Lgm_Vector(ctypes.Structure):
    _fields_ = [ ( "x", LgmDouble ),
        ("y", LgmDouble),
        ("z", LgmDouble) ]
Lgm_VectorP = ctypes.POINTER(Lgm_Vector)

