
import ctypes

class Lgm_Vector(ctypes.Structure):
    _fields_ = [ ( "x", ctypes.c_double ),
        ("y", ctypes.c_double),
        ("z", ctypes.c_double) ]


lgm = ctypes.CDLL('/usr/local/lib/libLanlGeoMag.dylib')
lgm.Lgm_Magnitude.argtypes=[ctypes.POINTER(Lgm_Vector)]
lgm.Lgm_Magnitude.restype = ctypes.c_double
vec = Lgm_Vector(1., 2., 3.)
ans = lgm.Lgm_Magnitude(vec)

if ans != 3.74165738677:
    print("Answers not equal")
    print("Got %f" % (ans))
