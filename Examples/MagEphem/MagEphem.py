#!/usr/bin/python

import sys
from optparse import OptionParser
import optparse
import time
import datetime
import os
import textwrap
import re
import math

#
# Kludgey way to create an enum -- which python doesnt natively support.
#
class Enum(set):
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError


class IndentedHelpFormatterWithNL(optparse.IndentedHelpFormatter):
#class IndentedHelpFormatterWithNL():
  def format_description(self, description):
    if not description: return ""
    desc_width = self.width - self.current_indent
    indent = " "*self.current_indent
# the above is still the same
    bits = description.split('\n')
    formatted_bits = [
      textwrap.fill(bit,
        desc_width,
        initial_indent=indent,
        subsequent_indent=indent)
      for bit in bits]
    result = "\n".join(formatted_bits) + "\n"
    return result

  def format_option(self, option):
    # The help for each option consists of two parts:
    #   * the opt strings and metavars
    #   eg. ("-x", or "-fFILENAME, --file=FILENAME")
    #   * the user-supplied help string
    #   eg. ("turn on expert mode", "read data from FILENAME")
    #
    # If possible, we write both of these on the same line:
    #   -x    turn on expert mode
    #
    # But if the opt string list is too long, we put the help
    # string on a second line, indented to the same column it would
    # start in if it fit on the first line.
    #   -fFILENAME, --file=FILENAME
    #       read data from FILENAME
    result = []
    opts = self.option_strings[option]
    opt_width = self.help_position - self.current_indent - 2
    if len(opts) > opt_width:
      opts = "%*s%s\n" % (self.current_indent, "", opts)
      indent_first = self.help_position
    else: # start help on same line as opts
      opts = "%*s%-*s  " % (self.current_indent, "", opt_width, opts)
      indent_first = 0
    result.append(opts)
    if option.help:
      help_text = self.expand_default(option)
# Everything is the same up through here
      help_lines = []
      for para in help_text.split("\n"):
        help_lines.extend(textwrap.wrap(para, self.help_width))
# Everything is the same after here
      result.append("%*s%s\n" % (
        indent_first, "", help_lines[0]))
      result.extend(["%*s%s\n" % (self.help_position, "", line)
        for line in help_lines[1:]])
    elif opts[-1] != "\n":
      result.append("\n")
    return "".join(result) 




#
# enum parse methods
#
ParseMethods = Enum(["UseSatelliteNumber", "UseIntDesignator", "None"])


#
# Set default values before we parse the command-line
#
CommonName      = ''
SatelliteNumber = ''
IntDesignator   = ''
ParseMethod     = ParseMethods.None
OutputFile      = ''
DumpToFile      = False
PurgeDuplicates = False
IntFieldModel   = 'IGRF'
ExtFieldModel   = 'T89'
Kp              = 5
Dst             = -20



#
# Define a command-line option parser and add the options we need
#
parser = OptionParser(  usage="%prog [options]",\
                        formatter=IndentedHelpFormatterWithNL(),\
                        version="%prog Version 1.01 (Movember 12, 2010)"  )

parser.add_option("-n", "--SatNumber",    dest="SatNumber",  
                        help="Satellite (or object) number ",
                        metavar="SAT_NUMBER")


parser.add_option("-s", "--Start",      dest="Start_ISO",    
                        help="Start date/time in ISO format",
                        metavar="YYYY-MM-DDTHH:MM:SS")

parser.add_option("-e", "--End",      dest="End_ISO",    
                        help="End date/time in ISO format",
                        metavar="YYYY-MM-DDTHH:MM:SS")

parser.add_option("-d", "--Delta",  dest="Delta_ISO",    
                        help="Time increment in ISO format",
                        metavar="HH:MM:SS")

parser.add_option("-q", "--Quality",    dest="Quality",  
                        help="Quality parameter for L* calculations [0-8]. Default is 3.",
                        metavar="QUALITY")

parser.add_option("-c", "--Colorize",   action="store_true", dest="Colorize", 
                        help="Print messages out in different colors when running multi-threaded.")

parser.add_option("-p", "--PitchAngles",    dest="PitchAngles",  
                        help="Range of pitch angles to compute. Default is \"5,90,5\" which gives 18 pitch angles from 5 to 90 degrees in 5 degree increments.",
                        metavar="PA_START, PA_END, PA_INC")

parser.add_option("-P", "--NoPitchAngles",   action="store_true", dest="NoPitchAngles", 
                        help="Do not compute pitch-angle dependent quantities.")

parser.add_option("-H", "--FootpointHeight",  dest="FootpointHeight",    
                        help="Geodetic altitude to use for footpoints. This is the height (in km) above the WGS84 geoid.",
                        metavar="HEIGHT_IN_KM")

parser.add_option("-a", "--Append",   action="store_true", dest="AppendMode", 
                        help="Append results to file (specified with the -o option).")


parser.add_option("-i", "--InputFile",   dest="InputFile", 
                        help="If specified, the times and positions to compute the magnetic ephemeris for will be read from INFILE. "\
                             "INFILE should have 4 columns: ISO DateTime String, and GEI_J2000 cartesian components in units of Re.",    
                        metavar="INFILE")

parser.add_option("-g", "--InputGeoFile",   dest="InputGeoFile", 
                        help="If specified, the times and positions to compute the magnetic ephemeris for will be read from INFILE. "\
                             "INFILE should have 4 columns: ISO DateTime String, and GEO Lat (Deg.), Lon (Deg.), Radius (Re).",    
                        metavar="INFILE")


parser.add_option("-o", "--OutputFile",   dest="OutputFile", 
                        help="Filename to dump output to. Stdout is assumed if "\
                             "no file given.",    
                        metavar="OUTFILE")

parser.add_option("-I", "--InternalFieldModel",   dest="IntFieldModel", 
                        help="Internal field model to use. Valid choices are:\n"\
                             "   CDIP - Centered Dipole (parameters taken from\n"\
                             "          first 3 IGRF coeffs.)\n"
                             "   EDIP - Eccentric dipole (parameters taken from\n"\
                             "          first 8 IGRF coeffs.)\n"
                             "   IGRF - International Geophysical Reference\n"\
                             "          Field\n"\
                             "Default is IGRF." )

parser.add_option("-E", "--ExternalFieldModel",   dest="ExtFieldModel", 
                        help="Internal field model to use. Valid choices are:\n"\
                             "   T87 - Tsyganenko 1987 model.\n"\
                             "   T89 - Tsyganenko 1989 model.\n"\
                             "   CDIP - Centered Dipole (parameters taken from\n"\
                             "          first 3 IGRF coeffs.)\n"
                             "   EDIP - Eccentric dipole (parameters taken from\n"\
                             "          first 8 IGRF coeffs.)\n"
                             "   IGRF - International Geophysical Reference\n"\
                             "          Field\n"\
                             "Default is T89." )



#
# Parse the args that were (potenitally) given to us on the command line
#
(options, args) = parser.parse_args()
print options
print args


if options.Colorize:
    Colorize = True
else:
    Colorize = False

if options.NoPitchAngles:
    DoPitchAngles = False
else:
    DoPitchAngles = True



if DoPitchAngles:
    if options.PitchAngles:
    #    pat = re.compile( r"\s*(\d+\.?\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*" )

        pat = re.compile( r"\s*([+-]?\d*\.?\d*)\s*,\s*([+-]?\d*\.?\d*)\s*,\s*([+-]?\d*\.?\d*)\s*" )



        m = re.match( pat, options.PitchAngles )
        if ( m != None ):
            l = m.groups()
            if ( len(l) != 3 ): 
                exit()
            else:
                ps = float(l[0])
                pe = float(l[1])
                pi = float(l[2])
                if (ps < 0.0)or(ps>90.0)or(pe < 0.0)or(pe>90.0):
                    print 'pitch angles must be in range [0,90]\n'
                    parser.print_help()
                    exit()
                if (math.fabs(pi) < .1):
                    print 'pitch angle increment too small\n'
                    parser.print_help()
                    exit()
            
        else:
            print 'Must provide 3 values for pitch angle range.\n'
            parser.print_help()
            exit()
    else:
        ps = 5.0
        pe = 90.0
        pi = 5.0
        
    if (pe >= ps):
        p = ps
        PitchAngles = []
        while ( p <= pe ):
            PitchAngles.append(p)
            p += pi
    else:
        p = ps
        PitchAngles = []
        while ( p >= pe ):
            PitchAngles.append(p)
            p -= pi

    nPitchAngles = len(PitchAngles)

    StrList = []
    StrList.append('{0:d}'.format(nPitchAngles))
    for p in PitchAngles:
        StrList.append(' {0:g}'.format(p))
    PA_Str = ''.join(StrList)
else:
    PA_Str = '0'


print PA_Str


#
# Set variables based on what we got on the command line
#

if options.InputFile != None:
    ReadFromFile = True
    InputFile = options.InputFile
else:
    ReadFromFile = False


if options.InputGeoFile != None:
    ReadFromGeoFile = True
    InputFile = options.InputGeoFile
else:
    ReadFromGeoFile = False




if options.OutputFile == None:
    print 'Must provide an output file via -o option\n'
    parser.print_help()
    exit()
else:
    OutputFile = options.OutputFile

    
if options.FootpointHeight != None:
    FootpointHeight = float(options.FootpointHeight)
else:
    FootpointHeight = 100.0

print FootpointHeight


# Quality setting
if options.Quality != None:
    Quality = int(options.Quality)
    if (Quality < 0): Quality = 0
    if (Quality > 8): Quality = 8
else:
    Quality = 3


if options.AppendMode:
    Append = True
else:
    Append = False

if options.IntFieldModel != None:
    IntFieldModel = options.IntFieldModel

if options.ExtFieldModel != None:
    ExtFieldModel = options.ExtFieldModel

PathName = os.path.dirname(sys.argv[0])
print PathName + '/puke'


if (ReadFromFile == True) or (ReadFromGeoFile == True):
        # 
        # Find correct Kp  and Dst
        # Kp as a floating point number is defined like: 4-, 5o, 5+ corresponds to (4.7, 5.0, 5.3)
        # We are letting MagEphemFromFile do the entire file which means that there will be a problem
        # with Kp/Dst stuff. Maybe we still want to do one line at a time.
        #


        # 
        # Add other information to the "input.txt" file. This is a slightly
        # different format to the one used for tle's
        #
        with open("input.txt", "w") as f:
            f.write('InputFile:'+InputFile+'\n')    # Input file
            f.write('OutputFile:'+OutputFile+'\n')    # Output file
            f.write('Internal Field Model:'+IntFieldModel+'\n') # The internal field model to use
            f.write('External Field Model:'+ExtFieldModel+'\n') # The external field model to use
            f.write('PitchAngles:'+PA_Str+'\n')
            f.write('Footpoint Height:'+str(FootpointHeight)+'\n')
            f.write('Colorize:'+str(Colorize)+'\n')
            f.write('Quality:'+str(Quality)+'\n')
            f.close()

            #exit()

        # 
        # Add start and end times to input.txt
        #
        if (ReadFromFile == True):
            if Append:
                os.system(PathName + '/MagEphemFromFile -a')
            else:
                os.system(PathName + '/MagEphemFromFile')
                Append = True
        else:
            if Append:
                os.system(PathName + '/MagEphemFromLatLonRad -a')
            else:
                os.system(PathName + '/MagEphemFromLatLonRad')
                Append = True

else:

    if options.SatNumber != None:
        ParseMethod     = ParseMethods.UseSatelliteNumber
        SatelliteNumber = options.SatNumber
        SatNum          = int( SatelliteNumber )
    else:
        print 'Must specify object catalog number via -n option.\n'
        parser.print_help()
        exit()



    # Start Date/Time
    if options.Start_ISO == None:
        print 'Must provide a start date/time via -s option\n'
        parser.print_help()
        exit()
    else:
        sdt = datetime.datetime.strptime(options.Start_ISO, "%Y-%m-%dT%H:%M:%S")
        s = int(time.mktime(sdt.timetuple()))
        #print s


    # End Date/Time
    if options.End_ISO == None:
        print 'Must provide a end date/time via -e option\n'
        parser.print_help()
        exit()
    else:
        edt = datetime.datetime.strptime(options.End_ISO, "%Y-%m-%dT%H:%M:%S")
        e = int(time.mktime(edt.timetuple()))
        #print e


    # Cadence
    if options.Delta_ISO == None:
        print 'Must provide a cadence or time increment via -d option\n'
        parser.print_help()
        exit()
    else:
        cad = datetime.datetime.strptime(options.Delta_ISO, "%H:%M:%S")
        delta = int(cad.hour*3600 + cad.minute*60 + cad.second)
        #print delta


    a = int(SatNum/1000)


    start = time.time()
    for t in range(s,e,delta):
        dt    = datetime.datetime.fromtimestamp(t)
        iso   = dt.isoformat()
        print "************************************************ Date: " + iso



        # 
        # Find correct TLE and dump to input3.txt
        #
        b     = dt.year
        c     = b-1
        fname1 = str.format('/home/mgh/TLE_DATABASE/{0}000_{0}999/{1}/{2}.txt', a, SatNum, b )
        fname2 = str.format('/home/mgh/TLE_DATABASE/{0}000_{0}999/{1}/{2}.txt', a, SatNum, c )
        command = str.format(PathName + '/FindTLEforGivenTime.py -n {0} -d {1} -o input3.txt {2} {3}', SatelliteNumber, iso, fname1, fname2)
        print command
        os.system(command)


        # 
        # Add other information to the "input.txt file
        #
        with open("input3.txt", "a") as f:
            f.write('OutputFile:'+OutputFile+'\n')    # Output file
            f.write('ISO Start Date/Time:'+iso+'\n')           # Start Date/Time
            f.write('ISO End   Date/Time:'+iso+'\n')           # End   Date/Time
            f.write('Cadence in seconds:'+str(delta)+'\n')          # Cadence in seconds
            f.write('Internal Field Model:'+IntFieldModel+'\n') # The internal field model to use
            f.write('External Field Model:'+ExtFieldModel+'\n') # The external field model to use
            f.write('PitchAngles:'+PA_Str+'\n')
            f.write('Footpoint Height:'+str(FootpointHeight)+'\n')
            f.write('Colorize:'+str(Colorize)+'\n')
            f.write('Quality:'+str(Quality)+'\n')
            f.close()

            #exit()

        # 
        # Add start and end times to input.txt
        #
        if Append:
            os.system(PathName + '/MagEphemFromTLE -a')
            
        else:
            os.system(PathName + '/MagEphemFromTLE')
            Append = True

        end = time.time()
        elapsed= end - start
        print "Elapsed Time: ", elapsed, "seconds"
            
        



