"""
Simple palce to keep all the types used by the Lgm wrapper

@author: Brian Larsen
@organization: LANL
@contact: balarsen@lanl.gov

@version: V1: 20-Dec-2010 (BAL)
"""




#set up types
LgmLong = ctypes.c_long
LgmBoolean = ctypes.c_int
ConstLgmBoolean = LgmBoolean
LgmFALSE = 0
LgmTRUE = 1
LgmChar = ctypes.c_char
ConstLgmChar = LgmChar
LgmCharP = ctypes.c_char_p
ConstLgmCharP = ctypes.c_char_p
LgmDouble = ctypes.c_double
ConstLgmDouble = LgmDouble
LgmDoubleP = ctypes.POINTER(LgmDouble)
#idx = c_sizes.index(ctypes.sizeof(LgmDouble) / 2)
ConstLgmInt = LgmInt
LgmIntP = ctypes.POINTER(LgmInt)
Lgm_VectorP = ctypes.POINTER(Lgm_Vector)
LgmInt = ctypes.c_int

c_types = [ctypes.c_byte, ctypes.c_short, ctypes.c_int,
           ctypes.c_long, ctypes.c_longlong]
c_sizes = [ctypes.sizeof(i) for i in c_types]

