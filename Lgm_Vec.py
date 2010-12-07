
import ctypes

class Lgm_Vector_cls(ctypes.Structure):
    _fields_ = [ ( "x", ctypes.c_double ),
        ("y", ctypes.c_double),
        ("z", ctypes.c_double) ]
    #_argtypes_ = [ctypes.c_double, ctypes.c_double, ctypes.c_double]



#set up types
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
c_types = [ctypes.c_byte, ctypes.c_short, ctypes.c_int,
           ctypes.c_long, ctypes.c_longlong]
c_sizes = [ctypes.sizeof(i) for i in c_types]
idx = c_sizes.index(ctypes.sizeof(LgmDouble) / 2)
LgmInt = c_types[idx]
ConstLgmInt = LgmInt
LgmIntP = ctypes.POINTER(LgmInt)
Lgm_Vector = ctypes.POINTER(Lgm_Vector_cls)

class Lgm_Vec(object):
    """
    Python wrapper for the vector function in Lgm_vec.c
    TODO wrap the matrix operations too

    @author: Brian Larsen
    @organization: LANL
    @contact: balarsen@lanl.gov

    @version: V1: 06-Dec-2010 (BAL)
    """
    def __init__(self):
        self.lib = ctypes.CDLL('/usr/local/lib/libLanlGeoMag.dylib')

#        lgm.Lgm_Magnitude.argtypes=[ctypes.POINTER(Lgm_Vector)]
#lgm.Lgm_Magnitude.restype = ctypes.c_double
#vec = Lgm_Vector(1., 2., 3.)
#ans = lgm.Lgm_Magnitude(vec)
#
#if ans != 3.74165738677:
#    print("Answers not equal")
#    print("Got %f" % (ans))


        self.call_dict = {
            'Lgm_CrossProduct': [None, Lgm_Vector, Lgm_Vector, Lgm_Vector],
            'Lgm_DotProduct': [LgmDouble, Lgm_Vector, Lgm_Vector],
            'Lgm_NormalizeVector': [LgmDouble, Lgm_Vector],
            'Lgm_ScaleVector': [None, Lgm_Vector, LgmDouble],
            'Lgm_Magnitude': [LgmDouble, Lgm_Vector],
            'Lgm_VecSub': [None, Lgm_Vector, Lgm_Vector, Lgm_Vector],
            'Lgm_VecAdd': [None, Lgm_Vector, Lgm_Vector, Lgm_Vector],
            'Lgm_VecDiffMag': [LgmDouble, Lgm_Vector, Lgm_Vector],
            'Lgm_ForceMagnitude': [None, Lgm_Vector, LgmDouble],
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
