
import ctypes
from Lgm_Types import *
from Lgm_Vector import *
from Lgm_CTrans import *


class _Lgm(object):
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
            # Lgm_Vector
                'Lgm_CrossProduct': [None, Lgm_VectorP, Lgm_VectorP, Lgm_VectorP],
                'Lgm_DotProduct': [LgmDouble, Lgm_VectorP, Lgm_VectorP],
                'Lgm_NormalizeVector': [LgmDouble, Lgm_VectorP],
                'Lgm_ScaleVector': [None, Lgm_VectorP, LgmDouble],
                'Lgm_Magnitude': [LgmDouble, Lgm_VectorP],
                'Lgm_VecSub': [None, Lgm_VectorP, Lgm_VectorP, Lgm_VectorP],
                'Lgm_VecAdd': [None, Lgm_VectorP, Lgm_VectorP, Lgm_VectorP],
                'Lgm_VecDiffMag': [LgmDouble, Lgm_VectorP, Lgm_VectorP],
                'Lgm_ForceMagnitude': [None, Lgm_VectorP, LgmDouble],
                'Lgm_MatTimesVec': [None, LgmDouble * 3 * 3, Lgm_VectorP, Lgm_VectorP],
                'Lgm_Transpose' : [None, LgmDouble * 3 * 3, LgmDouble * 3 * 3],
                'Lgm_MatTimesMat' : [None, LgmDouble * 3 * 3, LgmDouble * 3 * 3, LgmDouble * 3 * 3],
            # Lgm_LeapSeconds
                'Lgm_LeapYear': [LgmInt, LgmInt],
                'Lgm_GetLeapSeconds': [LgmDouble, LgmDouble, Lgm_CTransP],
                'Lgm_LoadLeapSeconds' : [LgmInt, Lgm_LeapSecondsP],
                'Lgm_GetLeapSeconds' : [LgmDouble, LgmDouble, Lgm_LeapSecondsP],
                'Lgm_IsLeapSecondDay' : [LgmInt, LgmInt, Lgm_LeapSecondsP],
                'Lgm_UTC_to_TAI' : [None, Lgm_DateTimeP, Lgm_DateTimeP, Lgm_CTransP],
                'Lgm_TAI_to_UTC' : [None, Lgm_DateTimeP, Lgm_DateTimeP, Lgm_CTransP],
                'Lgm_UTC_to_TT' : [None, Lgm_DateTimeP, Lgm_DateTimeP, Lgm_CTransP],
                'Lgm_TT_to_UTC' : [None, Lgm_DateTimeP, Lgm_DateTimeP, Lgm_CTransP],
                'Lgm_TT_to_TDB' : [None, Lgm_DateTimeP, Lgm_DateTimeP, Lgm_CTransP],
                'Lgm_TDB_to_TT' : [None, Lgm_DateTimeP, Lgm_DateTimeP, Lgm_CTransP],
                'Lgm_UTC_to_TDB' : [None, Lgm_DateTimeP, Lgm_DateTimeP, Lgm_CTransP],
                'Lgm_TDB_to_UTC' : [None, Lgm_DateTimeP, Lgm_DateTimeP, Lgm_CTransP],
                'Lgm_JD' : [LgmDouble, LgmInt, LgmInt, LgmInt, LgmDouble, LgmInt, Lgm_CTransP  ],
                #'' : [],
                #'' : [],
                #'' : [],
                #'' : [],
                #'' : [],
                #'' : [],
                #'' : [],

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
