#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import glob
import sys
import itertools
import datetime as dt
from optparse import OptionParser, IndentedHelpFormatter
import optparse
import textwrap
from multiprocessing import Pool, cpu_count
import numpy as np
from ctypes import pointer, c_int, c_double
from lgmpy import Lgm_CTrans, Lgm_Vector, magcoords
from lgmpy.Lgm_Wrap import _SgpInfo, _SgpTLE, LgmSgp_ReadTlesFromStrings, Lgm_JD, LGM_TIME_SYS_UTC, Lgm_init_ctrans, LgmSgp_SGP4_Init, Lgm_Convert_Coords, Lgm_Set_Coord_Transforms, LgmSgp_SGP4, WGS84_A, TEME_TO_GSE, TEME_TO_GEO

import dateutil.parser as dup
import spacepy.time as spt
import spacepy.toolbox as tb
import FindTLEforGivenTime as fTLE


def findNamed(path, name):
    """ Finds the path for the 'name'.

    Parameters
    ==========
    path : str
        path of the file
    name : str
        name of the file or directory being searched for

    Returns
    =======
    path : str
        The path 'name' is found in (excludes 'name' in the string)
    """
    pp = os.path.split(path)
    if pp[-1] == '':
        return None
    if pp[-1] != name:
        path = findNamed(pp[0], name)
    return path


basepath = findNamed(os.path.dirname(os.path.abspath(__file__)), 'dream')
sys.path.insert(-1, os.path.join(basepath, 'bin'))
from dream_option_parser import defaults

defaults = {
    'Delta': 5, 
    'outpath': os.path.join(basepath, 'Spacecraft'),
    'TLEpath': os.path.join(os.path.sep, 'n', 'space_data', 'TLE_DATABASE')}

WGS72_A = 6378.135

# the TLEs used here are assumed to have the (3-line) format:
#Line0 = 'NAVSTAR 49 (USA 154)'
#Line1 = '1 26605U 00071A   07067.92696543 -.00000083  00000-0  10000-3 0  5386'
#Line2 = '2 26605 056.6140 012.5765 0031970 244.1062 115.5621 02.00556447 46349'

# set start date and end date here
# fdata = gps.GPSlevel1([dt.datetime(2007,12,1), dt.datetime(2007,12,8)],
#                      bird='ns41', dataDir='/home/smorley/projects/gps/level_1/')
#use_times = fdata['Epoch']


def parserSetup():
    """ Sets up the option parser.

    Returns
    =======
    parser : OptionParser
        A configured OptionParser object.
    """
    # Define a command-line option parser and add the options we need
    parser = OptionParser(usage="%prog [options]",
                          formatter=IndentedHelpFormatter(),
                          version="%prog Version 1.1 (June 6, 2013)")

    parser.add_option("-s", "--Start", dest="Start_ISO",
                            help="Start date in ISO format [required]",
                            metavar="YYYY-MM-DD")

    parser.add_option("-e", "--End", dest="End_ISO",
                            help="End date in ISO format [required]",
                            metavar="YYYY-MM-DD")

    parser.add_option("-b", "--Bird", dest="Bird",
                            help="Satellite Name [required for header]")

    parser.add_option("-n", "--SatNum", dest="SatNum",
                            help="NORAD tracking number of object [required]")

    parser.add_option(
        "-d",
        "--Delta",
        dest="Delta",
        help="Time increment [minutes] for ephemeris data [default 5 minutes]")

    parser.add_option(
        "-t",
        "--TLEpath",
        dest="TLEpath",
        help="Path to TLEs for given Bird [default /n/space_data/TLE_DATABASE/XXXXX_XXXXX/SatNum]")

    parser.add_option(
        "-o",
        "--outpath",
        dest="outpath",
        help="Path for file output (/Ephem/YEAR is auto-appended) [default .]")

    parser.add_option(
        "-m",
        "--Mission",
        dest="Mission",
        help="Mission Name (e.g. ns41's mission is GPS) [optional, defaults to same as Bird]")

    return parser


def getLatLonRadfromTLE(epochs, TLEpath, options):
    """ Reads Latitude and Longitude from TLE.

    Parameters
    ==========
    epochs : 
        List of Ticktock objects.
    TLEpath : 
        Path to TLE file.
    options : optparse.Values
        Organized options from the command line.

    Returns
    =======
    testlat, testlon, testrad : list, list, list
        Latitude, Longitude, and Radius, respectively, each being
        a list of float values.
    """
    # now do Mike's setup for getting coords from TLE using SGP4
    pos_in = [0, 0, 0]
    s = _SgpInfo()
    TLEs = _SgpTLE()

    # loop over all times
    testlat = np.asarray(epochs).copy()
    testlat.fill(0)
    testlon = testlat.copy()
    testrad = testlat.copy()
    testtdiff = testlat.copy()
    print('Fetching TLEs & converting for range {0} to {1}'.format(
        epochs[0].isoformat(), epochs[-1].isoformat()))
    for idx, c_date in enumerate(epochs):
        #print('Doing {0}'.format(c_date))
        # put into JD as SGP4 needs serial time
        c = Lgm_init_ctrans(0)

        # now do Mike's setup for getting coords from TLE using SGP4
        dstr = int(c_date.strftime('%Y%j')) + c_date.hour / 24.0 + \
            c_date.minute / 1440.0 + c_date.second / 86400.0
        globstat = os.path.join(TLEpath, '*.txt')
        TLEfiles = glob.glob(globstat)
        if not TLEfiles:
            raise IOError(
                'No TLE files found in {0}. Aborting...'.format(TLEpath))
        Line0, Line1, Line2 = fTLE.findTLEinfiles(
                TLEfiles,
                ParseMethod='UseSatelliteNumber',
                TargetEpoch=dstr,
                SatelliteNumber=options.SatNum,
                Verbose=False,
                PurgeDuplicates=True)

        # print("{0}\n{1}\n{2}\n\n".format(Line0,Line1,Line2))
        nTLEs = c_int(0)
        LgmSgp_ReadTlesFromStrings(
            Line0,
            Line1,
            Line2,
            pointer(nTLEs),
            pointer(TLEs),
            1)
        LgmSgp_SGP4_Init(pointer(s), pointer(TLEs))
        date = Lgm_CTrans.dateToDateLong(c_date)
        utc = Lgm_CTrans.dateToFPHours(c_date)
        JD = Lgm_JD(
            c_date.year,
            c_date.month,
            c_date.day,
            utc,
            LGM_TIME_SYS_UTC,
            c)

        # Set up the trans matrices
        Lgm_Set_Coord_Transforms(date, utc, c)
        # get SGP4 output, needs minutes-since-TLE-epoch
        tsince = (JD - TLEs.JD) * 1440.0
        LgmSgp_SGP4(tsince, pointer(s))

        pos_in[0] = s.X
        pos_in[1] = s.Y
        pos_in[2] = s.Z

        Pin = Lgm_Vector.Lgm_Vector(*pos_in)
        Pout = Lgm_Vector.Lgm_Vector()
        Lgm_Convert_Coords(pointer(Pin), pointer(Pout), TEME_TO_GEO, c)
        PoutPy = Pout.tolist()
        PoutPy[0] /= WGS84_A
        PoutPy[1] /= WGS84_A
        PoutPy[2] /= WGS84_A
        nlat, nlon, nrad = Lgm_Vector.CartToSph(*PoutPy)
        testlat[idx] = nlat
        testlon[idx] = nlon
        testrad[idx] = nrad
        testtdiff[idx] = tsince / 1440.0
    return testlat, testlon, testrad


def writeEphemOut(epochs, lat, lon, rad, options):
    """ Creates a neatly formatted Ephem text file.

    Parameters
    ==========
    epochs : list or spacepy.time.Ticktock
        List of Ticktock objects (or just one) with times.
    lat : list
        List of Latitudes for each epoch.
    lon : list
        List of Longitudes for each epoch.
    rad : list
        List of radii for each epoch (in Re).
    options : optparse.Values
        Organized options from the command line.

    Returns
    =======
    None 
    """
    bird = options.Bird
    # first chunk into days...
    if isinstance(epochs, np.ndarray):
        epochs = epochs.tolist()

    def dayincrement(startval):
        """ Generator that increments a datetime object by one day.

        Parameters
        ==========
        startval : datetime
            Starting date.

        Returns
        =======
        None 
        """
        # badly written generator
        d1 = startval
        while True:
            yield d1
            d1 += dt.timedelta(days=1)

    dictofdays, dictofcoords = {}, {}
    validday, daynum = True, 1
    daycc = dayincrement(epochs[0])
    while validday:
        dayval = daycc.next()
        dayinds = [i for (i, t) in enumerate(epochs)
                   if t.date() == dayval.date()]
        if dayinds:
            dictofdays[daynum - 1] = epochs[dayinds[0]:dayinds[-1] + 1]
            dictofcoords[daynum - 1] = {}
            dictofcoords[daynum - 1]['lat'] = lat[dayinds[0]:dayinds[-1] + 1]
            dictofcoords[daynum - 1]['rad'] = rad[dayinds[0]:dayinds[-1] + 1]
            dictofcoords[daynum - 1]['lon'] = lon[dayinds[0]:dayinds[-1] + 1]
            daynum += 1
        else:
            validday = False

    # then write in format corresponding to:
    # 2002-09-01T00:07:30.000Z    0.77103341    8.31399822    6.61799765
    # i.e. 0 strftime(%Y-%m-%dT%h:%M:%S.%3fZ) 1: 1.8f  2: 1.8f  3: 1.8f
    listofdays = dictofdays.keys()
    for day in listofdays:
        dayepochs = dictofdays[day]
        #print('Working on output for {0}'.format(dayepochs[0].strftime('%Y%m%d')))
        fname = '{0}_{1}_eph.txt'.format(dayepochs[0].strftime('%Y%m%d'),
                                         bird)
        eph_path = os.path.join(
            options.outpath, 'Ephem', str(
                dayepochs[0].year))
        if not os.path.isdir(eph_path):
            os.makedirs(eph_path)
            print('Creating path {0}'.format(eph_path))
        fname = os.path.join(eph_path, fname)
        print(fname)
        fh = open(fname, 'w')
        # first write header
        cadence = (dayepochs[1] - dayepochs[0]).seconds // 60
        fh.write(
            ''.join(
                makeHeader(
                    bird,
                    options.Mission,
                    cadence,
                    dayepochs[0])))
        # now write stuff here
        for ind, linetime in enumerate(dayepochs):
            linelat = dictofcoords[day]['lat'][ind]
            linelon = dictofcoords[day]['lon'][ind]
            linerad = dictofcoords[day]['rad'][ind]

            fdict = {'time': linetime.strftime('%Y-%m-%dT%H:%M:%S.%fZ'),
                     'lat': linelat, 'lon': linelon, 'rad': linerad}
            lineout = '{time}    {lat: 1.8f}    {lon: 1.8f}    {rad: 1.8f}\n'.format(**fdict)
            fh.write(lineout)
        fh.close()
        #print('Did output for {0}'.format(dayepochs[0].strftime('%Y%m%d')))

    return None


def splitMonths(inputiter):
    """ Takes input iterable of datetimes and splits into a list 
        of datetimes by month.

    Parameters
    ==========
    inputiter : list
        List of spacepy.time.Ticktock objects.

    Returns
    =======
    outlist : list
        A two dimensional list -- list containing monthly-lists of
        Ticktock objects.
    """
    outlist, tmp = [], []
    monthcnt = 0
    curr_mn = inputiter[0].month
    for tt in inputiter:
        if tt.month == curr_mn:
            tmp.append(tt)
        else:
            curr_mn = tt.month
            outlist.append(tmp)
            tmp = []
            tmp.append(tt)
    outlist.append(tmp)
    return(outlist)


def makeHeader(bird, mission, cadence, date):
    """ Given bird (e.g. 'ns41'), mission (GPS) and cadence 
        (in minutes) make header for Ephm txt files.

    Parameters
    ==========
    bird : str
        The bird/satellite.
    mission : 
        Mission name for the INFO field.
    cadence : 
        Cadence/frequency of information; placed at the end of the INFO field.
    date : 
        Date of the information.

    Returns
    =======
    hline : list
        Header for Ephem txt files.
    """
    hline = [
        '# INFO: {0} {1} min ephemeris\n'.format(
            mission,
            cadence),
        '# SAT: {0}\n'.format(bird),
        '# DATE: {0}\n'.format(
            date.strftime('%Y%m%d')),
        '#                    UTC     Lat (deg)     Lon (deg)      Rad (Re)\n']
    #print('Making header for {0}'.format(date.strftime('%Y%m%d')))

    return hline


def doGLLR_WEO(inputs):
    """ Wrapper for multiprocessing.

    Parameters
    ==========
    inputs : list
        Contains parameters for this function (for multi-processing). 
        Should contain epochs, TLEpath, and options, in respective order.

    Returns
    =======
    None
    """
    epochs, TLEpath, options = inputs[0], inputs[1], inputs[2]
    lat, lon, rad = getLatLonRadfromTLE(epochs, TLEpath, options)
    writeEphemOut(epochs, lat, lon, rad, options)


if __name__ == '__main__':

    parser = parserSetup()
    # Parse the args that were (potentially) given to us on the command line
    (options, in_args) = parser.parse_args()

    # roll defaults into options dict and do any checking
    valid_opt = 0
    for key in options.__dict__.keys():
        if options.__dict__[key] is not None:
            valid_opt += 1
        try:
            if options.__dict__[key] is None:
                options.__dict__[key] = defaults[key]
        except KeyError:
            pass

    # check for any input arguments, if not, print help
    if not valid_opt:
        parser.print_help()
        exit()

    if not (options.Bird and options.SatNum):
        print('Must specify valid bird and NORAD number')
        parser.print_help()
        exit()

    if options.outpath == defaults['outpath']:
        # if using the default output path, now add the satellite name to it...
        options.outpath = os.path.join(options.outpath, options.Bird)

    if not options.Mission:
        options.Mission = options.Bird

    if options.Start_ISO and options.End_ISO:
        start_day = dup.parse(options.Start_ISO)
        end_day = dup.parse(options.End_ISO)
        # test whether start==end: if so, one day is being requested
        if start_day == end_day:
            end_day = end_day + dt.timedelta(days=1, microseconds=-1)
        if end_day.hour == 0 and end_day.minute == 0:  # if time is 00:00 then go to end of day
            end_day = end_day + dt.timedelta(days=1, microseconds=-1)
    else:
        print('Start and End times must be valid ISO')
        parser.print_help()
        exit()

    # TLEpath = '/n/space_data/TLE_DATABASE/26000_26999/26605'  #ns41
    partDir = '{0}000_{0}999'.format(int(options.SatNum) // 1000)
    if os.path.isdir(os.path.join(options.TLEpath, partDir, options.SatNum)):
        TLEpath = os.path.join(options.TLEpath, partDir, options.SatNum)
    elif os.path.isdir(os.path.join(options.TLEpath, options.SatNum)):
        TLEpath = os.path.join(options.TLEpath, options.SatNum)
    else:
        print('Specified TLE path not valid')
        parser.print_help()
        exit()

    epochs = spt.tickrange(
        start_day,
        end_day,
        dt.timedelta(
            minutes=float(
                options.Delta)),
        dtype='UTC')
    epochs = epochs.UTC
    print('Getting times between: {0} and {1}'.format(epochs[0], epochs[-1]))

    # Set up number of cpus required for multiprocessing of different months
    ncpus = cpu_count() - 1
    nmonths = 1 + (epochs[-1].year - epochs[0].year) * \
        12 + epochs[-1].month - epochs[0].month
    if nmonths <= ncpus:
        ncpus = nmonths
    pool = Pool(ncpus)
    if ncpus > 1:
        print('Multiprocessing EphemFromTLE request')
        # list of (lists of days per month to be processed)
        mdateslist = splitMonths(epochs)
        inputs = list(itertools.product(mdateslist, [TLEpath], [options]))
        result = pool.map(doGLLR_WEO, inputs)
    else:
        lat, lon, rad = getLatLonRadfromTLE(epochs, TLEpath, options)
        writeEphemOut(epochs, lat, lon, rad, options)
