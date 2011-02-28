from ctypes import pointer, CFUNCTYPE, c_int, c_void_p, cast

from Lgm_Wrap import Lgm_Vector, Lgm_MagModelInfo, Lgm_Set_Coord_Transforms, Lgm_TraceToEarth, Lgm_B_OP77
from Lgm_Wrap import Lgm_Set_Lgm_B_igrf, Lgm_Set_Lgm_B_T01S, Lgm_Set_gm_B_TS04, Lgm_Set_Lgm_B_T89, Lgm_Set_Lgm_B_OP77


import Lgm_MagModelInfo
import Lgm_Vector



Date = 20100203
UTC  = 12.34567

MagInfo = Lgm_MagModelInfo.Lgm_MagModelInfo( )
Lgm_Set_Coord_Transforms( Date, UTC, MagInfo.c )

#MagInfo.Bfield = cast((CFUNCTYPE(c_int))(Lgm_B_OP77), c_void_p)
Lgm_Set_Lgm_B_OP77(pointer(MagInfo))

MagInfo.Kp = 5

MagInfo.VerbosityLevel = 5;

u = Lgm_Vector.Lgm_Vector(1., 2., 3.)
v = Lgm_Vector.Lgm_Vector(0., 0., 0.)

u.x = -1.6
u.y =  0.0
u.z =  0.0 # Re
u.x = -1.601
u.y =  0.0
u.z =  0.0 # Re
u.x = -6.6
u.y =  0.0
u.z =  0.0 # Re
Lgm_TraceToEarth( pointer(u), pointer(v), 120.0, -11.0, 1e-7, pointer(MagInfo) )
print(u , 'maps to ', v)

#
#printf( "u = %g %g %g Re\n", u.x, u.y, u.z );
#printf( "v = %g %g %g Re\n", v.x, v.y, v.z );
#printf( "Manitude(v) = %g\n", Lgm_Magnitude( &v ) );
#
#
#u.x = -2.6; u.y =  0.0;  u.z =  0.0; // Re
#Lgm_TraceToEarth( &u, &v, 120.0, -11.0, 1e-7, MagInfo );
#printf( "u = %g %g %g Re\n", u.x, u.y, u.z );
#printf( "v = %g %g %g Re\n", v.x, v.y, v.z );
#printf( "Manitude(v) = %g\n", Lgm_Magnitude( &v ) );
