"""
Simple palce to keep all the types used by the Lgm wrapper

@author: Brian Larsen
@organization: LANL
@contact: balarsen@lanl.gov

@version: V1: 20-Dec-2010 (BAL)
"""

import ctypes


#set up types
LgmInt = ctypes.c_int
LgmIntP = ctypes.POINTER(LgmInt)

LgmLong = ctypes.c_long
LgmLongP = ctypes.POINTER(LgmLong)

LgmBoolean = ctypes.c_int
LgmBooleanP = ctypes.POINTER(LgmBoolean)

ConstLgmBoolean = LgmBoolean

LgmChar = ctypes.c_char
ConstLgmChar = LgmChar
LgmCharP = ctypes.c_char_p
ConstLgmCharP = ctypes.c_char_p

LgmDouble = ctypes.c_double
LgmDoubleP = ctypes.POINTER(LgmDouble)

ConstLgmDouble = LgmDouble

ConstLgmInt = LgmInt

c_types = [ctypes.c_byte, ctypes.c_short, ctypes.c_int,
           ctypes.c_long, ctypes.c_longlong]
c_sizes = [ctypes.sizeof(i) for i in c_types]

LgmFALSE = 0
LgmTRUE = 1
