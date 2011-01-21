import math
from ctypes import pointer, c_double, c_int

import numpy

from Lgm_Wrap import Lgm_Set_Coord_Transforms, SM_TO_GSM, Lgm_Convert_Coords, TRUE, Lgm_Set_Lgm_B_OP77, Lgm_LstarInfo, SetLstarTolerances, Lgm_Trace, GSM_TO_SM, WGS84_A, RadPerDeg, NewTimeLstarInfo, Lstar
import Lgm_Vector
import Lgm_CTrans
import Lgm_MagEphemInfo

ans = {}


Psm = Lgm_Vector.Lgm_Vector(0, 0, 0)
P = Lgm_Vector.Lgm_Vector(0, 0, 0)
c = Lgm_CTrans.Lgm_CTrans()

# this is the pitch angles to calculate
Alpha = range(1, 90, 20)  # 1...89

MagEphemInfo = Lgm_MagEphemInfo.Lgm_MagEphemInfo(len(Alpha), 0)

#// Date and UTC
Date       = 19800625;
UTC        = 19.0;
Lgm_Set_Coord_Transforms( Date, UTC, pointer(c) );

ans['Date'] = Date
ans['UTC'] = UTC

Psm.x = -6.6
Psm.y = 0.0
Psm.z = 0.0

Lgm_Convert_Coords( pointer(Psm), pointer(P), SM_TO_GSM, pointer(c) );
ans['PosSM'] = Psm.tolist()
ans['PosGSM'] = P.tolist()


MagEphemInfo.LstarQuality   = 3;
MagEphemInfo.SaveShellLines = TRUE;
MagEphemInfo.LstarInfo.contents.VerbosityLevel = 0;
MagEphemInfo.LstarInfo.contents.mInfo.contents.VerbosityLevel = 0;

Kp = 1;
#MagEphemInfo->LstarInfo->mInfo->Bfield        = Lgm_B_T89;
#MagEphemInfo->LstarInfo->mInfo->Bfield        = Lgm_B_cdip;
#MagEphemInfo->LstarInfo->mInfo->Bfield        = Lgm_B_OP77;
#MagEphemInfo->LstarInfo->mInfo->InternalModel = LGM_CDIP;
Lgm_Set_Lgm_B_OP77( MagEphemInfo.LstarInfo.contents.mInfo )

ans['Field'] = {}
ans['Field']['model'] = 'Lgm_B_OP77'
ans['Field']['Kp'] = Kp

MagEphemInfo.LstarInfo.contents.mInfo.contents.Kp = Kp

LstarInfo = Lgm_LstarInfo()

LstarInfo = MagEphemInfo.LstarInfo;


# Save Date, UTC to MagEphemInfo structure
MagEphemInfo.Date   = Date
MagEphemInfo.UTC    = UTC

# Save nAlpha, and Alpha array to MagEphemInfo structure
MagEphemInfo.nAlpha = len(Alpha)
MagEphemInfo.Alpha = (c_double*len(Alpha))(*Alpha)

# Set Tolerances
SetLstarTolerances(MagEphemInfo.LstarQuality, MagEphemInfo.LstarInfo )

# set coord transformation
Lgm_Set_Coord_Transforms(Date, UTC, MagEphemInfo.LstarInfo.contents.mInfo.contents.c)

# *  Blocal at sat location.
MagEphemInfo.P = P

Bvec = Lgm_Vector.Lgm_Vector(0,0,0)
# Get B at the point in question
MagEphemInfo.LstarInfo.contents.mInfo.contents.Bfield(pointer(P), pointer(Bvec),
                                                      MagEphemInfo.LstarInfo.contents.mInfo)


# save its magnitude in the structure
MagEphemInfo.B = Bvec.magnitude()

u=P
v1 = Lgm_Vector.Lgm_Vector(0,0,0)
v2 = Lgm_Vector.Lgm_Vector(0,0,0)
v3 = Lgm_Vector.Lgm_Vector(0,0,0)
vv1 = Lgm_Vector.Lgm_Vector(0,0,0)
TRACE_TOL = 1e-7

#  Compute Field-related quantities for each Pitch Angle.

if Lgm_Trace(pointer(u), pointer(v1), pointer(v2), pointer(v3),
          120, 0.01,
          TRACE_TOL, MagEphemInfo.LstarInfo.contents.mInfo) == 1:

    MagEphemInfo.Pmin = v3

    MagEphemInfo.Bmin = MagEphemInfo.LstarInfo.contents.mInfo.contents.Bmin
    # Get a simple measure of how big L is
    Lgm_Convert_Coords( pointer(v1), pointer(vv1), GSM_TO_SM,
                       MagEphemInfo.LstarInfo.contents.mInfo.contents.c );

    Lam = math.asin(vv1.z/vv1.magnitude())
    CosLam = math.cos(Lam)
    LSimple = (1.0+120.0/WGS84_A)/( CosLam*CosLam )
    ans['LSimple'] = LSimple


    # LOOP OVER PITCH ANGLES
    for i, pa in enumerate(Alpha):
        ans[str(pa)] = {}

        print(MagEphemInfo.LstarInfo.contents.PreStr)
        print(MagEphemInfo.LstarInfo.contents.PostStr)
        PreStr = MagEphemInfo.LstarInfo.contents.PreStr
        PostStr = MagEphemInfo.LstarInfo.contents.PostStr

    # Set Pitch Angle, sin, sin^2, and Bmirror
        sa = math.sin( pa*RadPerDeg )
        sa2 = sa*sa;

        print("%sComputing L* for Pitch Angle: Alpha[%d] = %g Date: %d   UTC: %g   Lsimple = %g%s\n" %
              ( PreStr, i, MagEphemInfo.Alpha[i], Date, UTC, LSimple, PostStr ))

        MagEphemInfo.LstarInfo.contents.mInfo.contents.Bm = MagEphemInfo.B/sa2
        ans[str(pa)]['Bm'] = MagEphemInfo.LstarInfo.contents.mInfo.contents.Bm

        Lgm_Set_Coord_Transforms( Date, UTC, MagEphemInfo.LstarInfo.contents.mInfo.contents.c )
        MagEphemInfo.LstarInfo.contents.PitchAngle = pa

        MagEphemInfo.Bm[i] = MagEphemInfo.LstarInfo.contents.mInfo.contents.Bm

        # Compute L*

        # USER SHOULD DECIDE THRESHOLD HERE
        if LSimple < 10.0:
            if MagEphemInfo.LstarInfo.contents.VerbosityLevel >= 2 :
                print("\n\n\t%sComputing L* for: UTC = %g PA = %d  (%g)%s\n" %
                      ( PreStr, UTC, i, MagEphemInfo.Alpha[i], PostStr ))
                print("    \t%s                  I   = %g PA = %d  (%g)%s\n" %
                    (PreStr, MagEphemInfo.I[i], i, MagEphemInfo.Alpha[i], PostStr ))
            LS_Flag = Lstar( pointer(v3), MagEphemInfo.LstarInfo)

            if MagEphemInfo.LstarInfo.contents.VerbosityLevel >= 2:
                print("\t%sUTC, L*          = %g %g%s\n" %
                       (PreStr, UTC, MagEphemInfo.LstarInfo.contents.LS, PostStr ))
                print("\t%sUTC, L*_McIlwain = %g %g%s\n" %
                       (PreStr, UTC, MagEphemInfo.LstarInfo.contents.LS_McIlwain_M, PostStr ))
                print("\t%sUTC, LSimple     = %g %g%s\n\n\n" %
                       (PreStr, UTC, LSimple, PostStr ))

            MagEphemInfo.Lstar[i] = MagEphemInfo.LstarInfo.contents.LS
            ans[str(pa)]['Lstar'] = MagEphemInfo.LstarInfo.contents.LS

            # Save results to the MagEphemInfo structure.
            MagEphemInfo.nShellPoints[i] = MagEphemInfo.LstarInfo.contents.nPnts
            # ShellI = numpy.ctypeslib.ndarray
            ## pull all this good extra info into numpy arrays
            ans[str(pa)]['ShellI'] = numpy.ctypeslib.ndarray([len(MagEphemInfo.LstarInfo.contents.I)],
                dtype=c_double, buffer=MagEphemInfo.LstarInfo.contents.I)
            ans[str(pa)]['ShellEllipsoidFootprint_Pn'] = \
                numpy.ctypeslib.ndarray(len(MagEphemInfo.LstarInfo.contents.Ellipsoid_Footprint_Pn),
                                        dtype=c_double,
                                        buffer=MagEphemInfo.LstarInfo.contents.Ellipsoid_Footprint_Pn)
            ans[str(pa)]['ShellEllipsoidFootprint_Ps'] = \
                numpy.ctypeslib.ndarray(len(MagEphemInfo.LstarInfo.contents.Ellipsoid_Footprint_Ps),
                                        dtype=c_double,
                                        buffer=MagEphemInfo.LstarInfo.contents.Ellipsoid_Footprint_Ps)

            ans[str(pa)]['ShellMirror_Pn'] = \
                numpy.ctypeslib.ndarray(len(MagEphemInfo.LstarInfo.contents.Mirror_Pn),
                                        dtype=c_double,
                                        buffer=MagEphemInfo.LstarInfo.contents.Mirror_Pn)
            ans[str(pa)]['ShellMirror_Ps'] = \
                numpy.ctypeslib.ndarray(len(MagEphemInfo.LstarInfo.contents.Mirror_Ps),
                                        dtype=c_double,
                                        buffer=MagEphemInfo.LstarInfo.contents.Mirror_Ps)
            ans[str(pa)]['ShellMirror_Ss'] = MagEphemInfo.LstarInfo.contents.mInfo.contents.Sm_South
            ans[str(pa)]['ShellMirror_Sn'] = MagEphemInfo.LstarInfo.contents.mInfo.contents.Sm_North
            ans[str(pa)]['nFieldPnts'] = \
                numpy.ctypeslib.ndarray(len(MagEphemInfo.LstarInfo.contents.nFieldPnts),
                                        dtype=c_int,
                                        buffer=MagEphemInfo.LstarInfo.contents.nFieldPnts)

            ans[str(pa)]['s_gsm'] = \
                numpy.ctypeslib.ndarray([len(MagEphemInfo.LstarInfo.contents.s_gsm),
                                         len(MagEphemInfo.LstarInfo.contents.s_gsm[0])],
                                        dtype=c_double,
                                        buffer=MagEphemInfo.LstarInfo.contents.s_gsm)
            ans[str(pa)]['Bmag'] = \
                numpy.ctypeslib.ndarray([len(MagEphemInfo.LstarInfo.contents.Bmag),
                                         len(MagEphemInfo.LstarInfo.contents.Bmag[0])],
                                        dtype=c_double,
                                        buffer=MagEphemInfo.LstarInfo.contents.Bmag)
            ans[str(pa)]['x_gsm'] = \
                numpy.ctypeslib.ndarray([len(MagEphemInfo.LstarInfo.contents.x_gsm),
                                         len(MagEphemInfo.LstarInfo.contents.x_gsm[0])],
                                        dtype=c_double,
                                        buffer=MagEphemInfo.LstarInfo.contents.x_gsm)
            ans[str(pa)]['y_gsm'] = \
                numpy.ctypeslib.ndarray([len(MagEphemInfo.LstarInfo.contents.y_gsm),
                                         len(MagEphemInfo.LstarInfo.contents.y_gsm[0])],
                                        dtype=c_double,
                                        buffer=MagEphemInfo.LstarInfo.contents.y_gsm)
            ans[str(pa)]['z_gsm'] = \
                numpy.ctypeslib.ndarray([len(MagEphemInfo.LstarInfo.contents.z_gsm),
                                         len(MagEphemInfo.LstarInfo.contents.z_gsm[0])],
                                        dtype=c_double,
                                        buffer=MagEphemInfo.LstarInfo.contents.z_gsm)

for i in range(len(Alpha)):
    print(Alpha[i], MagEphemInfo.Lstar[i])
