# -*- coding: utf-8 -*-

"""
Overview
--------
Class to compute Lstar using LanlGeoMag

This class computes Lstar, I, and other quantities for simple magnetic field models
  - eventually this will work with more complicated models as well, the C already does

"""

from __future__ import division

import math
from ctypes import pointer, c_double, c_int
import datetime

import numpy

from spacepy import datamodel
import spacepy.toolbox as tb

from Lgm_Wrap import Lgm_Set_Coord_Transforms, SM_TO_GSM, Lgm_Convert_Coords, \
    SetLstarTolerances, RadPerDeg, GSM_TO_WGS84, WGS84_TO_EDMAG,\
    LFromIBmM_Hilton, LFromIBmM_McIlwain, Lgm_EDMAG_to_R_MLAT_MLON_MLT, Lgm_FreeMagEphemInfo_Children
from Lgm_Wrap import Lstar as Lgm_Lstar
import Lgm_Vector
import Lgm_CTrans
import Lgm_MagEphemInfo
import Closed_Field
from _Bfield_dict import Bfield_dict

__author__ = 'Brian Larsen, Steve Morley (Python), Mike Henderson (C) - LANL'

class Lstar_Data(datamodel.SpaceData):
    """
    Class to store Lstar data and associated attributes
        - this class is a dictionary of numpy arrays that have attributes

    Examples
    ========
    >>> import spacepy.datamodel as dm
    >>> from lgmpy import Lstar
    >>> dat = Lstar.Lstar_Data()
    >>> dat.attrs['Author'] = 'John Doe'
    >>> dat['position'] = dm.dmarray([1,2,3], attrs={'system':'GSM'})
    >>> print(dat)
    +
    :|____Author (str [8])
    |____position (spacepy.datamodel.dmarray (3,))
        :|____system (str [3])
    ...

    See Also
    ========
      - spacepy.datamodel.Spacedata
      - spacepy.datamodel.dmarray

    """
    def __init__(self, *args, **kwargs):
        super(Lstar_Data, self).__init__(*args, **kwargs)
        self['position'] = {}

    def __repr__(self):
        tb.dictree(self, verbose=True, attrs=True)
        return ''

def get_Lstar(pos, date, alpha = 90.,
                  Kp = 2, coord_system='GSM',
                  Bfield = 'Lgm_B_OP77',
                  LstarThresh = 10.0,  # beyond this Lsimple don't compute Lstar
                  extended_out = False,
                  LstarQuality = 3):
    """
    This method does all the work for the calculation of Lstar.

    There are many options to set and a lot of output.  We have tried to describe
    it well here but certainly have missed something, if you can't understand something
    please contact the authors.

    Parameters
    ==========
    pos : list or array
        The position that Lstar should be calculated for
    date : datetime
        The date and time that Lstar should be calculated for
    alpha : float, optional
        Local pitch angle at ``pos`` to calculate Lstar, default (90)
    Kp : int, optional
        The Kp index to pass to the magnetic field model, used in T89, ignored
        in OP77, default (2)
    coord_system : str, optional
        The coordinate system of the input position, default (GSM)
    Bfield : str, optional
        The Magnetic field model to use for the calculation, default (Lgm_B_T89)
    LstarThresh : float, optional
        The calculation computes a simple L value and does not calculation Lstar
        if simple L is beyond this, default (10)
    extended_out : bool, optional
        Keyword that enables the output of significantly more information into
        the Lstar_Data output, default (False).  The extra information allows
        the points along the Lstar trace to be used since Lstar is the same at
        each point along the trace.
    LstarQuality : int
        The quality flag for the integrators in the calculation, this can have
        serious impact on run speed and Lstar value

    Returns
    =======
    out : Lstar_Data
        Returns a datamodel object that contains all the information for the run.
        See examples for the contents of this object for extended_out=False and True

    Examples
    ========
    >>> from lgmpy import Lstar
    >>> import datetime
    >>> dat = Lstar.get_Lstar([1,2,0], datetime.datetime(2000, 1, 2), extended_out=False)
    >>> print(dat)
    +
    |____90.0 (spacepy.datamodel.SpaceData [8])
         |____Bmin (spacepy.datamodel.dmarray ())
             :|____units (str [2])
         |____Bmirror (spacepy.datamodel.dmarray ())
             :|____coord_system (str [3])
             :|____units (str [2])
         |____I (spacepy.datamodel.dmarray ())
         |____LHilton (float)
         |____LMcIlwain (float)
         |____Lsimple (spacepy.datamodel.dmarray (1,))
         |____Lstar (spacepy.datamodel.dmarray (1,))
             :|____info (str [4])
         |____Pmin (spacepy.datamodel.dmarray (3,))
             :|____coord_system (str [3])
             :|____units (str [3])
    |____Bcalc (spacepy.datamodel.dmarray (3,))
        :|____Kp (int)
        :|____coord_system (str [3])
        :|____model (str [10])
        :|____units (str [2])
    |____Epoch (spacepy.datamodel.dmarray (1,))
    |____Kp (spacepy.datamodel.dmarray (1,))
    |____MLT (spacepy.datamodel.dmarray ())
        :|____coord_system (str [5])
    |____position (dict [1])
         |____GSM (spacepy.datamodel.dmarray (3,))
             :|____units (str [2])
    ...

    The data object is a dictionary with keys:
        - 90.0 : this is the pitch angle specified (multiple can be specified)
            - Bmin : the minimum B for that pitch angle and position
            - BMirror : the position of the mirror point for that position
            - I : the I value for that position
            - LHilton : L value calculated with the Hilton approximation
            - LMcIlwain: L values calculated with the McIlwain formula
            - Lsimple : a simple L value
            - Lstar : the value of Lstar for that position and pitch angle
            - Pmin : TODO what is Pmin?
        - Bcalc : information about the magnetic field model
        - Epoch : the time of the calculation
        - Kp : Kp used for the calculation
        - MLT : MLT value for the position and coord_system
        - position : the position of the calculation
            - GSM : the value in this system, if conversions are done they all appear here

    >>> from lgmpy import Lstar
    >>> import datetime
    >>> dat = Lstar.get_Lstar([1,2,0], datetime.datetime(2000, 1, 2), extended_out=True)
    >>> print(dat)
    +
    |____90.0 (spacepy.datamodel.SpaceData [21])
        :|____Calc_Time (float)
         |____Bmag (numpy.ndarray (100, 1000))
         |____Bmin (spacepy.datamodel.dmarray ())
             :|____units (str [2])
         |____Bmirror (spacepy.datamodel.dmarray ())
             :|____coord_system (str [3])
             :|____units (str [2])
         |____I (spacepy.datamodel.dmarray ())
         |____LHilton (float)
         |____LMcIlwain (float)
         |____Lsimple (spacepy.datamodel.dmarray (1,))
         |____Lstar (spacepy.datamodel.dmarray (1,))
             :|____info (str [4])
         |____Pmin (spacepy.datamodel.dmarray (3,))
             :|____coord_system (str [3])
             :|____units (str [3])
         |____ShellEllipsoidFootprint_Pn (numpy.ndarray (100,))
         |____ShellEllipsoidFootprint_Ps (numpy.ndarray (100,))
         |____ShellI (spacepy.datamodel.dmarray (100,))
         |____ShellMirror_Pn (numpy.ndarray (100,))
         |____ShellMirror_Ps (numpy.ndarray (100,))
         |____ShellMirror_Sn (float)
         |____ShellMirror_Ss (float)
         |____nFieldPnts (numpy.ndarray (100,))
         |____s_gsm (numpy.ndarray (100, 1000))
         |____x_gsm (numpy.ndarray (100, 1000))
         |____y_gsm (numpy.ndarray (100, 1000))
         |____z_gsm (numpy.ndarray (100, 1000))
    |____Bcalc (spacepy.datamodel.dmarray (3,))
        :|____Kp (int)
        :|____coord_system (str [3])
        :|____model (str [10])
        :|____units (str [2])
    |____Epoch (spacepy.datamodel.dmarray (1,))
    |____Kp (spacepy.datamodel.dmarray (1,))
    |____MLT (spacepy.datamodel.dmarray ())
        :|____coord_system (str [5])
    |____position (dict [1])
         |____GSM (spacepy.datamodel.dmarray (3,))
             :|____units (str [2])
    ...

    The data object is a dictionary with keys:
        - 90.0 : this is the pitch angle specified (multiple can be specified)
            - Calc_Time : the clock time that the calculation took to run
            - Bmag : the magnitude of B along the Lstar trace
            - Bmin : the minimum B for that pitch angle and position
            - BMirror : the position of the mirror point for that position
            - I : the I value for that position
            - LHilton : L value calculated with the Hilton approximation
            - LMcIlwain: L values calculated with the McIlwain formula
            - Lsimple : a simple L value
            - Lstar : the value of Lstar for that position and pitch angle
            - Pmin : TODO what am I?
            - ShellEllipsoidFootprint_Pn : TODO what am I?
            - ShellEllipsoidFootprint_Ps : TODO what am I?
            - ShellI : TODO what am I?
            - ShellMirror_Pn : TODO what am I?
            - ShellMirror_Ps : TODO what am I?
            - ShellMirror_Sn : TODO what am I?
            - ShellMirror_Ss : TODO what am I?
            - nFieldPnts : TODO what am I?
            - s_gsm : TODO what am I?
            - x_gsm : TODO what am I?
            - y_gsm : TODO what am I?
            - z_gsm : TODO what am I?
        - Bcalc : information about the magnetic field model
        - Epoch : the time of the calculation
        - Kp : Kp used for the calculation
        - MLT : MLT value for the position and coord_system
        - position : the position of the calculation
            - GSM : the value in this system, if conversions are done they all appear here


    """

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

    Pwgs = Lgm_Vector.Lgm_Vector()
    Pmlt = Lgm_Vector.Lgm_Vector()
    Lgm_Convert_Coords( pointer(Pgsm), pointer(Pwgs), GSM_TO_WGS84, mmi.c )
    Lgm_Convert_Coords( pointer(Pwgs), pointer(Pmlt), WGS84_TO_EDMAG, mmi.c )
    R, MLat, MLon, MLT = c_double(), c_double(), c_double(), c_double(),
    Lgm_EDMAG_to_R_MLAT_MLON_MLT( pointer(Pmlt),  pointer(R), pointer(MLat), pointer(MLon),
        pointer(MLT), mmi.c)
    ans['MLT'] = datamodel.dmarray(MLT.value, attrs={'coord_system': 'EDMAG'})

    # save Kp
    # TODO maybe add some Kp checking
    ans['Kp'] = datamodel.dmarray([Kp])

    # Set the LstarQuality, TODO add a comment here on what each does
    MagEphemInfo.LstarQuality = LstarQuality

    # L* in one place is L* in lots of places (for GPS set to False)
    MagEphemInfo.SaveShellLines = extended_out
    # TODO maybe not hardcoded, but for now its fine
    MagEphemInfo.LstarInfo.contents.VerbosityLevel = 0
    MagEphemInfo.LstarInfo.contents.mInfo.contents.VerbosityLevel = 0

    #MagEphemInfo->LstarInfo->mInfo->Bfield        = Lgm_B_T89;
    #MagEphemInfo->LstarInfo->mInfo->Bfield        = Lgm_B_cdip;
    #MagEphemInfo->LstarInfo->mInfo->Bfield        = Lgm_B_OP77;
    #MagEphemInfo->LstarInfo->mInfo->InternalModel = LGM_CDIP;

    ## decide which field model to use, this is a keyword
    try:
        Bfield_dict[Bfield](MagEphemInfo.LstarInfo.contents.mInfo)
    except KeyError:
        raise(NotImplementedError("Only Bfield=%s currently supported" % Bfield_dict.keys()))

    MagEphemInfo.LstarInfo.contents.mInfo.contents.Kp = Kp
    # Save Date, UTC to MagEphemInfo structure ** is this needed?
    MagEphemInfo.Date   = datelong
    MagEphemInfo.UTC    = utc

    # Save nAlpha, and Alpha array to MagEphemInfo structure
    MagEphemInfo.nAlpha = len(Alpha)
    for i in range(len(Alpha)):
        MagEphemInfo.Alpha[i] = Alpha[i]

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
        #sets up nans in case of Lstar failure
        ans[pa]['I'] = datamodel.dmarray(numpy.nan)
        ans[pa]['Lstar'] = datamodel.dmarray(numpy.nan, attrs={'info':trace})
        ans[pa]['LMcIlwain'] = datamodel.dmarray(numpy.nan)
        ans[pa]['LHilton'] = datamodel.dmarray(numpy.nan)
        ans[pa]['Bmin'] = datamodel.dmarray(numpy.nan, attrs={'units':'nT'})
        ans[pa]['Bmirror'] = datamodel.dmarray(numpy.nan, attrs={'units':'nT'})

    if trace != 'LGM_CLOSED':
        return ans
        # if this is not LGM_CLOSED then don't both with any pitch angle?  true?

    #  Save field-related quantities for each Pitch Angle.
    MagEphemInfo.Pmin = Lgm_Vector.Lgm_Vector(*minB)
    MagEphemInfo.Bmin = MagEphemInfo.LstarInfo.contents.mInfo.contents.Bmin
    ans[pa]['Pmin'] = datamodel.dmarray(minB, attrs={'units':'R_E'})
    ans[pa]['Pmin'].attrs['coord_system'] = 'GSM'
    ans[pa]['Bmin'] = datamodel.dmarray(MagEphemInfo.Bmin, attrs={'units':'nT'})

    # LOOP OVER PITCH ANGLES
    for i, pa in enumerate(Alpha):

        tnow = datetime.datetime.now()
        PreStr = MagEphemInfo.LstarInfo.contents.PreStr
        PostStr = MagEphemInfo.LstarInfo.contents.PostStr

        # ***********************************************
        # *** not sure I fully understand this chunk, B at the mirror point*****
        # Set Pitch Angle, sin, sin^2, and Bmirror
        sa = math.sin( pa*RadPerDeg )
        sa2 = sa*sa

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
        if Lsimple < LstarThresh:
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
            elif LS_Flag == -1: # mirror below northern hemisphere mirror alt
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

    Lgm_FreeMagEphemInfo_Children(pointer(MagEphemInfo))



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
