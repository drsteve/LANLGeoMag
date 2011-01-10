"""
Main import to use the LanlGeoMag

import Lgm
then import what you want next

@author: Brian Larsen
@organization: LANL
@contact: balarsen@lanl.gov

@version: V1: 21-Dec-2010 (BAL)
"""

import ctypes
from _Lgm import lib
import Lgm_Vector
Lgm_Vector.Lgm_Vector.assign_fields()

from _Lgm_Eop import Lgm_NgaEopp, Lgm_Eop, Lgm_EopOne
from Lgm_DateAndTime import Lgm_DateAndTime

from Lgm_Types import LgmDouble, LgmInt, LgmLong, LgmDoubleP, LgmIntP, LgmLongP, LgmCharP
from Lgm_MagModelInfo import Lgm_MagModelInfo, Lgm_MagModelInfoP

from Lgm_CTrans import Lgm_DateTime, Lgm_CTrans




Lgm_VectorP = ctypes.POINTER(Lgm_Vector.Lgm_Vector)
Lgm_DateAndTimeP = ctypes.POINTER(Lgm_DateAndTime)
Lgm_DateTimeP = ctypes.POINTER(Lgm_DateTime)
Lgm_CTransP = ctypes.POINTER(Lgm_CTrans)
Lgm_NgaEoppP = ctypes.POINTER(Lgm_NgaEopp)
Lgm_EopP = ctypes.POINTER(Lgm_Eop)
Lgm_EopOneP = ctypes.POINTER(Lgm_EopOne)



# the order here matters, if you get a _fields_ is final then play with order
# it becomes final when pointers are created as Pyhton has to figure how big
# the object is (is the current theory)

Lgm_DateTime.assign_fields()
Lgm_DateAndTime.assign_fields()
Lgm_CTrans.assign_fields()
Lgm_Eop.assign_fields()
Lgm_EopOne.assign_fields()
Lgm_MagModelInfo.assign_fields()



# format is [Out type, In type1, In type2, In type3, ... ]
call_dict = {
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
        # Lgm_DateAndTime
        'Lgm_GetLeapSeconds': [LgmDouble, LgmDouble, Lgm_CTransP],
        'Lgm_LoadLeapSeconds' : [LgmInt, Lgm_CTransP],
        'Lgm_GetLeapSeconds' : [LgmDouble, LgmDouble, Lgm_CTransP],
        'Lgm_IsLeapSecondDay' : [LgmInt, LgmLong, LgmDoubleP, Lgm_CTransP],
        'Lgm_LeapYear': [LgmInt, LgmInt],
        'Lgm_JD' : [LgmDouble, LgmInt, LgmInt, LgmInt, LgmDouble, LgmInt, Lgm_CTransP  ],
        'Lgm_Date_to_JD' : [LgmDouble, LgmLong, LgmDouble, Lgm_CTransP],
        'Lgm_read_eop' : [None, Lgm_EopP],
        'Lgm_get_eop_at_JD' : [None, LgmDouble, Lgm_EopOneP, Lgm_EopP],
        'Lgm_set_eop' : [None, Lgm_EopOneP, Lgm_CTransP],
        'Lgm_JD_to_Date' : [LgmLong, LgmDouble, LgmIntP, LgmIntP, LgmIntP, LgmDoubleP],
        'Lgm_TAI_to_GPS' : [None, Lgm_CTransP, Lgm_CTransP, Lgm_CTransP],
        'Lgm_GPS_to_TAI' : [None, Lgm_CTransP, Lgm_CTransP, Lgm_CTransP],
        'Lgm_UTC_to_GPS' : [None, Lgm_CTransP, Lgm_CTransP, Lgm_CTransP],
        'Lgm_GPS_to_UTC' : [None, Lgm_CTransP, Lgm_CTransP, Lgm_CTransP],
        'Lgm_GPS_to_GpsSeconds' : [LgmDouble, Lgm_CTransP],
        'Lgm_GpsSeconds_to_GPS' : [None, LgmDouble, Lgm_CTransP],
        'Lgm_GpsSeconds_to_UTC' : [None, LgmDouble, Lgm_CTransP, Lgm_CTransP],
        'Lgm_UTC_to_GpsSeconds' : [LgmDouble, Lgm_CTransP, Lgm_CTransP],
        'Lgm_TAI_to_TaiSeconds' : [LgmDouble, Lgm_CTransP],
        'Lgm_TaiSeconds_to_TAI' : [None, LgmDouble, Lgm_CTransP],
        'Lgm_TaiSeconds_to_UTC' : [None, LgmDouble, Lgm_CTransP, Lgm_CTransP],
        'Lgm_UTC_to_TaiSeconds' : [LgmDouble, Lgm_CTransP, Lgm_CTransP],
        'Lgm_UTC_to_TAI' : [None, Lgm_CTransP, Lgm_CTransP, Lgm_CTransP],
        'Lgm_TAI_to_UTC' : [None, Lgm_CTransP, Lgm_CTransP, Lgm_CTransP],
        'Lgm_TT_to_TAI' : [None, Lgm_CTransP, Lgm_CTransP, Lgm_CTransP],
        'Lgm_TAI_to_TT' : [None, Lgm_CTransP, Lgm_CTransP, Lgm_CTransP],
        'Lgm_TT_to_TDB' : [None, Lgm_CTransP, Lgm_CTransP, Lgm_CTransP],
        'Lgm_TT_to_TDB_IAU2006' : [None, Lgm_CTransP, Lgm_CTransP, Lgm_CTransP],
        'Lgm_TDB_to_TT' : [None, Lgm_CTransP, Lgm_CTransP, Lgm_CTransP],
        'Lgm_UTC_to_TT' : [None, Lgm_CTransP, Lgm_CTransP, Lgm_CTransP],
        'Lgm_TT_to_UTC' : [None, Lgm_CTransP, Lgm_CTransP, Lgm_CTransP],
        'Lgm_UTC_to_TDB' : [None, Lgm_CTransP, Lgm_CTransP, Lgm_CTransP],
        'Lgm_TDB_to_UTC' : [None, Lgm_CTransP, Lgm_CTransP, Lgm_CTransP],
        'Lgm_DateTime_Create' : [Lgm_DateTime, LgmInt, LgmInt, LgmInt, LgmDouble, LgmInt, Lgm_CTransP],
        'Lgm_Make_UTC' : [LgmInt, LgmLong, LgmDouble, Lgm_DateAndTimeP, Lgm_CTransP],
        'Lgm_Print_DateTime' : [None, Lgm_DateTime, LgmInt, LgmInt],
        'Lgm_DateTimeToString' : [None, LgmCharP, Lgm_DateTime, LgmInt, LgmInt], # todo 80 was made up use a pointer
        #'Lgm_Print_SimpleTime' : [None, Lgm_CTransP, LgmInt, ctypes.create_string_buffer(80)],# todo 80 was made up use a pointer
        'Lgm_JDN' : [LgmLong, LgmInt, LgmInt, LgmInt],
        'Lgm_MJD' : [LgmDouble, LgmInt, LgmInt, LgmInt, LgmDouble, LgmInt, Lgm_CTransP ],
        'Lgm_Date_to_JD' : [LgmDouble, LgmLong, LgmDouble, Lgm_CTransP],
        'Lgm_DayOfYear' : [LgmInt, LgmInt, LgmInt, LgmInt, Lgm_CTransP],
        #'Lgm_DayOfWeek' : [LgmInt, LgmInt, LgmInt, LgmInt, ctypes.create_string_buffer(80)],# todo 80 was made up use a pointer
        'Lgm_JDNofWeek1' : [LgmLong, LgmInt],
        'Lgm_ISO_WeekNumber' : [LgmInt, LgmInt, LgmInt, LgmInt, LgmIntP],
        'Lgm_MaxWeekNumber' : [LgmInt, LgmInt],
        'Lgm_ISO_YearWeekDow_to_Date' : [None, LgmInt, LgmInt, LgmInt, LgmLongP, LgmIntP, LgmIntP, LgmIntP],
        'Lgm_Doy' : [LgmInt, LgmLong, LgmIntP, LgmIntP, LgmIntP, LgmIntP],
        'Lgm_IsValidDate' : [LgmInt, LgmLong],
        'Lgm_RemapTime' : [LgmDouble, LgmDouble, LgmDouble],
        'Lgm_UTC_to_TDBSeconds' : [LgmDouble, Lgm_CTransP, Lgm_CTransP],
        'Lgm_TDBSecSinceJ2000' : [LgmDouble, Lgm_CTransP, Lgm_CTransP],
        'Lgm_UTC_to_TAI' : [None, Lgm_CTransP, Lgm_CTransP, Lgm_CTransP],
        'Lgm_TAI_to_UTC' : [None, Lgm_CTransP, Lgm_CTransP, Lgm_CTransP],
        'Lgm_UTC_to_TT' : [None, Lgm_CTransP, Lgm_CTransP, Lgm_CTransP],
        'Lgm_TT_to_UTC' : [None, Lgm_CTransP, Lgm_CTransP, Lgm_CTransP],
        'Lgm_TT_to_TDB' : [None, Lgm_CTransP, Lgm_CTransP, Lgm_CTransP],
        'Lgm_TDB_to_TT' : [None, Lgm_CTransP, Lgm_CTransP, Lgm_CTransP],
        'Lgm_UTC_to_TDB' : [None, Lgm_CTransP, Lgm_CTransP, Lgm_CTransP],
        'Lgm_TDB_to_UTC' : [None, Lgm_CTransP, Lgm_CTransP, Lgm_CTransP],
        # coord transforms
        'Lgm_Set_Coord_Transforms' : [None, LgmLong, LgmDouble, Lgm_CTransP ],
        'Lgm_Convert_Coords' : [None, Lgm_VectorP, Lgm_VectorP, LgmInt, Lgm_CTransP],
        'size_t_size' : [LgmInt],
        'size_MagModelInfo' : [LgmInt],
        'size_CTrans' : [LgmInt],
        'size_Vector' : [LgmInt],
        'size_DateTime' : [LgmInt],
        'size_gsl_interp_accel' : [LgmInt],
        'size_gsl_interp_type' : [LgmInt],
        'size_gsl_interp' : [LgmInt],
        'size_gsl_spline' : [LgmInt],




        'Lgm_Set_Open_Limits' : [None, Lgm_MagModelInfoP, LgmDouble, LgmDouble, LgmDouble, LgmDouble, LgmDouble, LgmDouble ],
        'Lgm_Set_LossConeHeight' : [None, Lgm_MagModelInfoP, LgmDouble],
        'Lgm_Set_Octree_kNN_MaxDist' : [Lgm_MagModelInfoP, LgmDouble],
        'Lgm_Set_Octree_kNN_k' : [Lgm_MagModelInfoP, LgmInt],
        'Lgm_Set_Octree_kNN_InterpMethod' : [None, Lgm_MagModelInfoP, LgmInt],
        'Lgm_MagModelInfo_Set_Psw' : [None, LgmDouble, Lgm_MagModelInfoP],
        'Lgm_MagModelInfo_Set_Kp' : [LgmDouble, Lgm_MagModelInfoP],
        'Lgm_B_T89' : [LgmInt, Lgm_VectorP, Lgm_VectorP, Lgm_MagModelInfoP ] }
        #'' : [],
        #'' : [],
        #'' : [],
        #'' : [],
        #'' : [],


def load_call_dict(lib, call_dict):
    """Use call dictionary to set argument and return types for functions

    Uses ctypes methods to set the return types and arguments for Lgm
    functions, as defined in L{call_dict}. Called on setup, but should
    be called again after updating L{call_dict}.

    @author: Jon Niehof
    @organization: LANL
    @contact: jniehof@lanl.gov

    @version: V1: 06-Dec-2010 (JN)
    """
    for funcname in call_dict:
        func = getattr(lib, funcname)
        args = call_dict[funcname]
        func.restype = args[0]
        if len(args) <= 1:
            func.argtypes = None
        else:
            func.argtypes = args[1:]

load_call_dict(lib, call_dict)
