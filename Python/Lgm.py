
import ctypes
from Lgm_Structures import *

class Lgm(object):
    """
    Python wrapper for LanlGeoMag
    TODO incomplete wrapper, do more

    @author: Brian Larsen
    @organization: LANL
    @contact: balarsen@lanl.gov

    @version: V1: 06-Dec-2010 (BAL)
    """
    def __init__(self):
        self.lib = ctypes.CDLL('/usr/local/lib/libLanlGeoMag.dylib')

        self.call_dict = {
                'Lgm_CrossProduct': [None, Lgm_Vector_p, Lgm_Vector_p, Lgm_Vector_p],
                'Lgm_DotProduct': [LgmDouble, Lgm_Vector_p, Lgm_Vector_p],
                'Lgm_NormalizeVector': [LgmDouble, Lgm_Vector_p],
                'Lgm_ScaleVector': [None, Lgm_Vector_p, LgmDouble],
                'Lgm_Magnitude': [LgmDouble, Lgm_Vector_p],
                'Lgm_VecSub': [None, Lgm_Vector_p, Lgm_Vector_p, Lgm_Vector_p],
                'Lgm_VecAdd': [None, Lgm_Vector_p, Lgm_Vector_p, Lgm_Vector_p],
                'Lgm_VecDiffMag': [LgmDouble, Lgm_Vector_p, Lgm_Vector_p],
                'Lgm_ForceMagnitude': [None, Lgm_Vector_p, LgmDouble],
                'Lgm_LeapYear': [LgmInt, LgmInt],
                'Lgm_GetLeapSeconds': [LgmDouble, LgmDouble, Lgm_CTrans_p]
#                'Lgm_B_T89': [int, Lgm_Vector_p, Lgm_Vector_p, Lgm_MagModelInfo_p]
            }
        self.load_call_dict()

    def load_call_dict(self):
        """Use call dictionary to set argument and return types for functions

        Uses ctypes methods to set the return types and arguments for Lgm
        functions, as defined in L{call_dict}. Called on setup, but should
        be called again after updating L{call_dict}.

        @author: Jon Niehof
        @organization: LANL
        @contact: jniehof@lanl.gov

        @version: V1: 06-Dec-2010 (JN)
        """
        for funcname in self.call_dict:
            func = getattr(self.lib, funcname)
            args = self.call_dict[funcname]
            func.restype = args[0]
            if len(args) <= 1:
                func.argtypes = None
            else:
                func.argtypes = args[1:]
