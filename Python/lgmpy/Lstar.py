# -*- coding: utf-8 -*-

"""
Class to compute Lstar from position, date, and pitch angle
for a given Kp and Bfield model

see: Lstar.get_Lstar
"""

from __future__ import division

import math
from ctypes import pointer, c_double, c_int
import datetime

import numpy

from spacepy import datamodel
import spacepy.toolbox as tb

from Lgm_Wrap import Lgm_Set_Coord_Transforms, SM_TO_GSM, Lgm_Convert_Coords, \
    Lgm_Set_Lgm_B_OP77, SetLstarTolerances, \
    RadPerDeg, Lgm_Set_Lgm_B_T89, \
    LFromIBmM_Hilton, LFromIBmM_McIlwain
from Lgm_Wrap import Lstar as Lgm_Lstar
import Lgm_Vector
import Lgm_CTrans
import Lgm_MagEphemInfo
import Closed_Field



class Lstar_Data(datamodel.SpaceData):
    """
    Class to store data and attributes

    @author: Brian Larsen
    @organization: LANL
    @contact: balarsen@lanl.gov

    @version: V1: 03-Mar-2011 (BAL)
    """
    def __init__(self, *args, **kwargs):
        super(Lstar_Data, self).__init__(*args, **kwargs)
        self['position'] = {}

    def __repr__(self):
        tb.dictree(self, verbose=True, attrs=True)
        return ''

def get_Lstar(pos, date, alpha = 90,
                  Kp = 2, coord_system='GSM',
                  Bfield = 'Lgm_B_OP77',
                  LstarThres = 10.0,  # beyond this Lsimple don't compute Lstar
                  extended_out = False,
                  LstarQuality = 3):

    # setup a datamodel object to hold the answer
    ans = Lstar_Data()

    # change datetime to Lgm Datelong and UTC
    try:
        datelong = Lgm_CTrans.dateToDateLong(date)
        utc = Lgm_CTrans.dateToFPHours(date)
    except AttributeError:
        raise(TypeError("Date must be a datetime object"))
    else:
        ans['Epoch'] = datamodel.dmarray([date])

    # pitch angles to calculate
    try:
        Alpha = list(alpha)
    except TypeError:
        Alpha = [alpha]

    # required setup
    MagEphemInfo = Lgm_MagEphemInfo.Lgm_MagEphemInfo(len(Alpha), 0)

    # setup a shortcut to MagModelInfo
    mmi = MagEphemInfo.LstarInfo.contents.mInfo.contents
    Lgm_Set_Coord_Transforms( datelong, utc, mmi.c) # dont need pointer as it is one

    # convert to **GSM**
    if coord_system == 'GSM':
        try:
            Pgsm = Lgm_Vector.Lgm_Vector(*pos)
        except TypeError:
            raise(TypeError("Position must be listlike" ) )
        ans['position']['GSM'] = datamodel.dmarray(pos, attrs={'units':'Re'})
    elif coord_system == 'SM':
        try:
            Psm = Lgm_Vector.Lgm_Vector(*pos)
        except TypeError:
            raise(TypeError("Position must be listlike" ) )
        Pgsm = Lgm_Vector.Lgm_Vector()
        Lgm_Convert_Coords( pointer(Psm), pointer(Pgsm), SM_TO_GSM, mmi.c )
        ans['position']['SM'] = datamodel.dmarray(pos, attrs={'units':'Re'})
        ans['position']['GSM'] = datamodel.dmarray(Pgsm.tolist(), attrs={'units':'Re'})
    else:
        raise(NotImplementedError("Only GSM or SM input currently supported"))

    # decide which field model to use, this is a keyword
    Bfield_dict = {'Lgm_B_OP77': Lgm_Set_Lgm_B_OP77,
                   'Lgm_B_T89': Lgm_Set_Lgm_B_T89}

    # save Kp
    # TODO maybe add some Kp checking
    ans['Kp'] = datamodel.dmarray([Kp])

    # Set the LstarQuality, TODO add a comment here on what each does
    MagEphemInfo.LstarQuality = LstarQuality;

    # L* in one place is L* in lots of places (for GPS set to False)
    MagEphemInfo.SaveShellLines = extended_out
    # TODO maybe not hardcoded, but for now its fine
    MagEphemInfo.LstarInfo.contents.VerbosityLevel = 0;
    MagEphemInfo.LstarInfo.contents.mInfo.contents.VerbosityLevel = 0;

    #MagEphemInfo->LstarInfo->mInfo->Bfield        = Lgm_B_T89;
    #MagEphemInfo->LstarInfo->mInfo->Bfield        = Lgm_B_cdip;
    #MagEphemInfo->LstarInfo->mInfo->Bfield        = Lgm_B_OP77;
    #MagEphemInfo->LstarInfo->mInfo->InternalModel = LGM_CDIP;

    ## decide which field model to use, this is a keyword
    try:
        Bfield_dict[Bfield](MagEphemInfo.LstarInfo.contents.mInfo)
    except KeyError:
        raise(NotImplementedError("Only Bfield='OP77, T89' currently supported"))

    MagEphemInfo.LstarInfo.contents.mInfo.contents.Kp = Kp
    # Save Date, UTC to MagEphemInfo structure ** is this needed?
    MagEphemInfo.Date   = datelong
    MagEphemInfo.UTC    = utc

    # Save nAlpha, and Alpha array to MagEphemInfo structure
    MagEphemInfo.nAlpha = len(Alpha)
    MagEphemInfo.Alpha = (c_double*len(Alpha))(*Alpha)

    # Set Tolerances
    SetLstarTolerances(LstarQuality, MagEphemInfo.LstarInfo )
    # *  Blocal at sat location
    MagEphemInfo.P = Pgsm

    Bvec = Lgm_Vector.Lgm_Vector()
    # Get B at the point in question
    MagEphemInfo.LstarInfo.contents.mInfo.contents.Bfield(pointer(Pgsm), pointer(Bvec),
                                                          MagEphemInfo.LstarInfo.contents.mInfo)
    ans['Bcalc'] = datamodel.dmarray(Bvec.tolist(), attrs={'units':'nT'})
    ans['Bcalc'].attrs['model'] = Bfield
    ans['Bcalc'].attrs['Kp'] = Kp
    ans['Bcalc'].attrs['coord_system'] = 'GSM'

    # save its magnitude in the structure
    MagEphemInfo.B = Bvec.magnitude()

    # check and see if the field line is closed before doing much work
    trace, northern, southern, minB, Lsimple = Closed_Field.Closed_Field(MagEphemInfo, extended_out=True)

    # presetup the ans[Angle] so that it can be filled correctly
    for pa in Alpha:
        ans[pa] = datamodel.SpaceData()
        ans[pa]['Lsimple'] = datamodel.dmarray([Lsimple])

    if trace != 'LGM_CLOSED':
        for pa in Alpha:
            ans[pa]['Lstar'] = datamodel.dmarray(numpy.nan, attrs={'info':trace})
        return ans
        # if this is not LGM_CLOSED then don't both with any pitch angle?  true?

    #  Save field-related quantities for each Pitch Angle.
    MagEphemInfo.Pmin = Lgm_Vector.Lgm_Vector(*minB)
    MagEphemInfo.Bmin = MagEphemInfo.LstarInfo.contents.mInfo.contents.Bmin
    ans[pa]['Bmin'] = datamodel.dmarray(minB, attrs={'units':'nT'})
    ans[pa]['Bmin'].attrs['coord_system'] = 'GSM'

    # LOOP OVER PITCH ANGLES
    for i, pa in enumerate(Alpha):

        tnow = datetime.datetime.now()
        PreStr = MagEphemInfo.LstarInfo.contents.PreStr
        PostStr = MagEphemInfo.LstarInfo.contents.PostStr

        # ***********************************************
        # *** not sure I fully understand this chunk, B at the mirror point*****
        # Set Pitch Angle, sin, sin^2, and Bmirror
        sa = math.sin( pa*RadPerDeg )
        sa2 = sa*sa;

        # print("{0}Computing L* for Pitch Angle: Alpha[{1}] = {2} Date: {3}   UTC: {4}   Lsimple = {5:4.3}{6}\n").format(PreStr, i, MagEphemInfo.Alpha[i], date, utc, Lsimple, PostStr )

        MagEphemInfo.LstarInfo.contents.mInfo.contents.Bm = MagEphemInfo.B/sa2
        ans[pa]['Bmirror'] = datamodel.dmarray(MagEphemInfo.LstarInfo.contents.mInfo.contents.Bm, attrs={'units':'nT'})
        ans[pa]['Bmirror'].attrs['coord_system'] = 'GSM'
        # ***********************************************

        # I tink this is already done
        # Lgm_Set_Coord_Transforms( Date, UTC, MagEphemInfo.LstarInfo.contents.mInfo.contents.c )

        MagEphemInfo.LstarInfo.contents.PitchAngle = pa
        MagEphemInfo.Bm[i] = MagEphemInfo.LstarInfo.contents.mInfo.contents.Bm
        # Compute L*
        if Lsimple < LstarThres:
            Ls_vec = Lgm_Vector.Lgm_Vector(*minB)

            LS_Flag = Lgm_Lstar( pointer(Ls_vec), MagEphemInfo.LstarInfo)

            lstarinf = MagEphemInfo.LstarInfo.contents #shortcut
            MagEphemInfo.LHilton.contents.value = LFromIBmM_Hilton(c_double(lstarinf.I[0]),
                                                c_double(MagEphemInfo.Bm[i]),
                                                c_double(lstarinf.mInfo.contents.c.contents.M_cd))
            ans[pa]['LHilton'] = MagEphemInfo.LHilton.contents.value
            MagEphemInfo.LMcIlwain.contents.value = LFromIBmM_McIlwain(c_double(lstarinf.I[0]),
                                                c_double(MagEphemInfo.Bm[i]),
                                                c_double(lstarinf.mInfo.contents.c.contents.M_cd))
            ans[pa]['LMcIlwain'] = MagEphemInfo.LMcIlwain.contents.value
            if LS_Flag == -2: # mirror below southern hemisphere mirror alt
                ans[pa]['Lstar'] = datamodel.dmarray([numpy.nan], attrs={'info':'S_LOSS'})
            elif LS_Flag == -1: # mirror below nothern hemisphere mirror alt
                ans[pa]['Lstar'] = datamodel.dmarray([numpy.nan], attrs={'info':'N_LOSS'})
            elif LS_Flag == 0: # valid calc
                ans[pa]['Lstar'] = datamodel.dmarray([lstarinf.LS], attrs={'info':'GOOD'}) # want better word?

            MagEphemInfo.Lstar[i] = lstarinf.LS

            # Save results to the MagEphemInfo structure.
            MagEphemInfo.nShellPoints[i] = lstarinf.nPnts
            ## pull all this good extra info into numpy arrays
            ans[pa]['I'] = datamodel.dmarray(lstarinf.I[0])

            if extended_out:
                ans[pa]['ShellI'] = \
                    datamodel.dmarray(numpy.ctypeslib.ndarray([len(lstarinf.I)],
                                            dtype=c_double, buffer=lstarinf.I) )
                ans[pa]['ShellEllipsoidFootprint_Pn'] = \
                    numpy.ctypeslib.ndarray(len(lstarinf.Ellipsoid_Footprint_Pn),
                                            dtype=c_double,
                                            buffer=lstarinf.Ellipsoid_Footprint_Pn)
                ans[pa]['ShellEllipsoidFootprint_Ps'] = \
                    numpy.ctypeslib.ndarray(len(lstarinf.Ellipsoid_Footprint_Ps),
                                            dtype=c_double,
                                            buffer=lstarinf.Ellipsoid_Footprint_Ps)

                ans[pa]['ShellMirror_Pn'] = \
                    numpy.ctypeslib.ndarray(len(lstarinf.Mirror_Pn),
                                            dtype=c_double,
                                            buffer=lstarinf.Mirror_Pn)
                ans[pa]['ShellMirror_Ps'] = \
                    numpy.ctypeslib.ndarray(len(lstarinf.Mirror_Ps),
                                            dtype=c_double,
                                            buffer=lstarinf.Mirror_Ps)
                ans[pa]['ShellMirror_Ss'] = lstarinf.mInfo.contents.Sm_South
                ans[pa]['ShellMirror_Sn'] = lstarinf.mInfo.contents.Sm_North
                ans[pa]['nFieldPnts'] = \
                    numpy.ctypeslib.ndarray(len(lstarinf.nFieldPnts),
                                            dtype=c_int,
                                            buffer=lstarinf.nFieldPnts)

                ans[pa]['s_gsm'] = \
                    numpy.ctypeslib.ndarray([len(lstarinf.s_gsm),
                                             len(lstarinf.s_gsm[0])],
                                            dtype=c_double,
                                            buffer=lstarinf.s_gsm)
                ans[pa]['Bmag'] = \
                    numpy.ctypeslib.ndarray([len(lstarinf.Bmag),
                                             len(lstarinf.Bmag[0])],
                                            dtype=c_double,
                                            buffer=lstarinf.Bmag)
                ans[pa]['x_gsm'] = \
                    numpy.ctypeslib.ndarray([len(lstarinf.x_gsm),
                                             len(lstarinf.x_gsm[0])],
                                            dtype=c_double,
                                            buffer=lstarinf.x_gsm)
                ans[pa]['y_gsm'] = \
                    numpy.ctypeslib.ndarray([len(lstarinf.y_gsm),
                                             len(lstarinf.y_gsm[0])],
                                            dtype=c_double,
                                            buffer=lstarinf.y_gsm)
                ans[pa]['z_gsm'] = \
                    numpy.ctypeslib.ndarray([len(lstarinf.z_gsm),
                                             len(lstarinf.z_gsm[0])],
                                            dtype=c_double,
                                            buffer=lstarinf.z_gsm)
                delT = datetime.datetime.now() - tnow
                ans[pa].attrs['Calc_Time'] = delT.seconds + delT.microseconds/1e6
    return ans

if __name__ == '__main__':
    date = datetime.datetime(2010, 10, 12)
    ans = get_Lstar([-4.2, 1, 1], date, alpha = 90, Kp = 4, coord_system='SM', Bfield = 'Lgm_B_T89', LstarQuality = 1, extended_out=True)
    print('Lgm_B_T89 Kp=4')
    print ans[90]['LHilton']
    print ans[90]['LMcIlwain']
    print ans[90]['Lstar']
    print ans[90]['Lsimple']

    ans = get_Lstar([-4.2, 1, 1], date, alpha = 90, Kp = 5, coord_system='SM', Bfield = 'Lgm_B_T89', LstarQuality = 1, extended_out=True)
    print('Lgm_B_T89 Kp=5')
    print ans[90]['LHilton']
    print ans[90]['LMcIlwain']
    print ans[90]['Lstar']
    print ans[90]['Lsimple']

    ans = get_Lstar([-4.2, 1, 1], date, alpha = 90, Kp = 4, coord_system='SM', Bfield = 'Lgm_B_OP77', LstarQuality = 1, extended_out=True)
    print('Lgm_B_OP77 Kp=4')
    print ans[90]['LHilton']
    print ans[90]['LMcIlwain']
    print ans[90]['Lstar']
    print ans[90]['Lsimple']

    ans = get_Lstar([-4.2, 1, 1], date, alpha = 90, Kp = 5, coord_system='SM', Bfield = 'Lgm_B_OP77', LstarQuality = 1, extended_out=True)
    print('Lgm_B_OP77 Kp=6')
    print ans[90]['LHilton']
    print ans[90]['LMcIlwain']
    print ans[90]['Lstar']
    print ans[90]['Lsimple']
