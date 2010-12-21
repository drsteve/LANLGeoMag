"""
Main Lgm Library file, contains all the fucntion definitions


@author: Brian Larsen
@organization: LANL
@contact: balarsen@lanl.gov

@version: V1: 20-Dec-2010 (BAL)
"""

import ctypes
from Lgm_Types import *
from Lgm_CTrans import *
from Lgm_Vector import *
from Lgm_Eop import *
import sys

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
        if sys.platform == 'win32':
            raise(NotImplementedError("Sorry windows isn't up and running"))
        elif sys.platform == 'darwin':
            libname = 'libLanlGeoMag.dylib'
        else:
            libname = 'libLanlGeoMag.so'

        self.lib = ctypes.CDLL('/usr/local/lib/' + libname )

        # format is [Out type, In type1, In type2, In type3, ... ]
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
                # Lgm_DateAndTime
                'Lgm_GetLeapSeconds': [LgmDouble, LgmDouble, Lgm_CTransP],
                'Lgm_LoadLeapSeconds' : [LgmInt, Lgm_LeapSecondsP],
                'Lgm_GetLeapSeconds' : [LgmDouble, LgmDouble, Lgm_LeapSecondsP],
                'Lgm_IsLeapSecondDay' : [LgmInt, LgmInt, Lgm_LeapSecondsP],
                'Lgm_LeapYear': [LgmInt, LgmInt],
                'Lgm_JD' : [LgmDouble, LgmInt, LgmInt, LgmInt, LgmDouble, LgmInt, Lgm_CTransP  ],
                'Lgm_Date_to_JD' : [LgmDouble, LgmLong, LgmDouble, Lgm_CTransP],
                'Lgm_read_eop' : [None, Lgm_EopP],
                'Lgm_get_eop_at_JD' : [None, LgmDouble, Lgm_EopOneP, Lgm_EopP],
                'Lgm_set_eop' : [None, Lgm_EopOneP, Lgm_CTransP],
                'Lgm_JD_to_Date' : [LgmLong, LgmDouble, LgmIntP, LgmIntP, LgmIntP, LgmDoubleP],
                'Lgm_TAI_to_GPS' : [None, Lgm_DateTimeP, Lgm_DateTimeP, Lgm_CTransP],
                'Lgm_GPS_to_TAI' : [None, Lgm_DateTimeP, Lgm_DateTimeP, Lgm_CTransP],
                'Lgm_UTC_to_GPS' : [None, Lgm_DateTimeP, Lgm_DateTimeP, Lgm_CTransP],
                'Lgm_GPS_to_UTC' : [None, Lgm_DateTimeP, Lgm_DateTimeP, Lgm_CTransP],
                'Lgm_GPS_to_GpsSeconds' : [LgmDouble, Lgm_DateTimeP],
                'Lgm_GpsSeconds_to_GPS' : [None, LgmDouble, Lgm_DateTimeP],
                'Lgm_GpsSeconds_to_UTC' : [None, LgmDouble, Lgm_DateTimeP, Lgm_CTransP],
                'Lgm_UTC_to_GpsSeconds' : [LgmDouble, Lgm_DateTimeP, Lgm_CTransP],
                'Lgm_TAI_to_TaiSeconds' : [LgmDouble, Lgm_DateTimeP],
                'Lgm_TaiSeconds_to_TAI' : [None, LgmDouble, Lgm_DateTimeP],
                'Lgm_TaiSeconds_to_UTC' : [None, LgmDouble, Lgm_DateTimeP, Lgm_CTransP],
                'Lgm_UTC_to_TaiSeconds' : [LgmDouble, Lgm_DateTimeP, Lgm_CTransP],
                'Lgm_UTC_to_TAI' : [None, Lgm_DateTimeP, Lgm_DateTimeP, Lgm_CTransP],
                'Lgm_TAI_to_UTC' : [None, Lgm_DateTimeP, Lgm_DateTimeP, Lgm_CTransP],
                'Lgm_TT_to_TAI' : [None, Lgm_DateTimeP, Lgm_DateTimeP, Lgm_CTransP],
                'Lgm_TAI_to_TT' : [None, Lgm_DateTimeP, Lgm_DateTimeP, Lgm_CTransP],
                'Lgm_TT_to_TDB' : [None, Lgm_DateTimeP, Lgm_DateTimeP, Lgm_CTransP],
                'Lgm_TT_to_TDB_IAU2006' : [None, Lgm_DateTimeP, Lgm_DateTimeP, Lgm_CTransP],
                'Lgm_TDB_to_TT' : [None, Lgm_DateTimeP, Lgm_DateTimeP, Lgm_CTransP],
                'Lgm_UTC_to_TT' : [None, Lgm_DateTimeP, Lgm_DateTimeP, Lgm_CTransP],
                'Lgm_TT_to_UTC' : [None, Lgm_DateTimeP, Lgm_DateTimeP, Lgm_CTransP],
                'Lgm_UTC_to_TDB' : [None, Lgm_DateTimeP, Lgm_DateTimeP, Lgm_CTransP],
                'Lgm_TDB_to_UTC' : [None, Lgm_DateTimeP, Lgm_DateTimeP, Lgm_CTransP],
                'Lgm_DateTime_Create' : [Lgm_DateTime, LgmInt, LgmInt, LgmInt, LgmDouble, LgmInt, Lgm_CTransP],
                'Lgm_Make_UTC' : [LgmInt, LgmLong, LgmDouble, Lgm_DateTimeP, Lgm_CTransP],
                'Lgm_Print_DateTime' : [None, Lgm_DateTime, LgmInt, LgmInt],
                #'Lgm_DateTimeToString' : [None, ctypes.create_string_buffer(80), Lgm_DateTime, LgmInt, LgmInt], # todo 80 was made up use a pointer
                #'Lgm_Print_SimpleTime' : [None, Lgm_DateTimeP, LgmInt, ctypes.create_string_buffer(80)],# todo 80 was made up use a pointer
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
                'Lgm_UTC_to_TDBSeconds' : [LgmDouble, Lgm_DateTimeP, Lgm_CTransP],
                'Lgm_TDBSecSinceJ2000' : [LgmDouble, Lgm_DateTimeP, Lgm_CTransP],
                'Lgm_UTC_to_TAI' : [None, Lgm_DateTimeP, Lgm_DateTimeP, Lgm_CTransP],
                'Lgm_TAI_to_UTC' : [None, Lgm_DateTimeP, Lgm_DateTimeP, Lgm_CTransP],
                'Lgm_UTC_to_TT' : [None, Lgm_DateTimeP, Lgm_DateTimeP, Lgm_CTransP],
                'Lgm_TT_to_UTC' : [None, Lgm_DateTimeP, Lgm_DateTimeP, Lgm_CTransP],
                'Lgm_TT_to_TDB' : [None, Lgm_DateTimeP, Lgm_DateTimeP, Lgm_CTransP],
                'Lgm_TDB_to_TT' : [None, Lgm_DateTimeP, Lgm_DateTimeP, Lgm_CTransP],
                'Lgm_UTC_to_TDB' : [None, Lgm_DateTimeP, Lgm_DateTimeP, Lgm_CTransP],
                'Lgm_TDB_to_UTC' : [None, Lgm_DateTimeP, Lgm_DateTimeP, Lgm_CTransP],
                # coord transforms
                'Lgm_Set_Coord_Transforms' : [None, LgmInt, LgmDouble, Lgm_CTransP ],
                'Lgm_Convert_Coords' : [None, Lgm_VectorP, Lgm_VectorP, LgmInt, Lgm_CTransP] }
                #'' : [],
                #'' : [],
                #'' : [],
                #'' : [],
                #'' : [],
                #'' : [],
                #'' : [],
                #'' : [],
                #'' : [],
                #'' : [],
                #'' : [],
                #'' : [],
                #'' : [],
                #'' : [],


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
