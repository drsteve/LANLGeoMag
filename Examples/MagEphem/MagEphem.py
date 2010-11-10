#!/usr/bin/python

import sys
from optparse import OptionParser
import time
import datetime
import os

#
# Kludgey way to create an enum -- which python doesnt natively support.
#
class Enum(set):
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError






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



#
# Define a command-line option parser and add the options we need
#
parser = OptionParser(  usage="%prog [options]",\
                        version="%prog Version 1.00 (Movember 10, 2010)"  )

parser.add_option("-s", "--Start",      dest="Start_ISO",    
                        help="Start date/time in ISO format",
                        metavar="YYYY-MM-DDTHH:MM:SS")

parser.add_option("-e", "--End",      dest="End_ISO",    
                        help="End date/time in ISO format",
                        metavar="YYYY-MM-DDTHH:MM:SS")

parser.add_option("-c", "--Cadence",  dest="Cadence_ISO",    
                        help="Time increment in ISO format",
                        metavar="HH:MM:SS")

parser.add_option("-n", "--SatNumber",    dest="SatNumber",  
                        help="Satellite (or object) number (including trailing U "\
                             "or S - for example 29155U).", 
                        metavar="SAT_NUMBER")

parser.add_option("-a", "--Append",   action="store_true", dest="AppendMode", 
                        help="Append results to file (specified with the -o option).")


parser.add_option("-o", "--OutputFile",   dest="OutputFile", 
                        help="Filename to dump output to. Stdout is assumed if "\
                             "no file given.",    
                        metavar="OUTFILE")



#
# Parse the args that were (potenitally) given to us on the command line
#
(options, args) = parser.parse_args()
print options
print args


#
# Set variables based on what we got on the command line
#
if options.SatNumber != None:
    ParseMethod     = ParseMethods.UseSatelliteNumber
    SatelliteNumber = options.SatNumber
    SatNum          = int( SatelliteNumber )
#elif options.ID != None:
#    ParseMethod   = ParseMethods.UseIntDesignator
#    IntDesignator = options.ID
else:
    print 'Must specify object catalog number via -n option.\n'
    parser.print_help()
    exit()


if options.OutputFile == None:
    print 'Must provide an output file via -o option\n'
    parser.print_help()
    exit()
else:
    OutputFile = options.OutputFile

if options.AppendMode:
    Append = True
else:
    Append = False




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
if options.Cadence_ISO == None:
    print 'Must provide a cadence or time increment via -c option\n'
    parser.print_help()
    exit()
else:
    cad = datetime.datetime.strptime(options.Cadence_ISO, "%H:%M:%S")
    delta = int(cad.hour*3600 + cad.minute*60 + cad.second)
    #print delta


a = int(SatNum/1000)


for t in range(s,e,delta):
    dt    = datetime.datetime.fromtimestamp(t)
    iso   = dt.isoformat()

    # 
    # Find correct TLE and dump to input.txt
    #
    b     = dt.year
    c     = b-1
    fname1 = str.format('/home/mgh/TLE_DATABASE/{0}000_{0}999/{1}/{2}.txt', a, SatNum, b )
    fname2 = str.format('/home/mgh/TLE_DATABASE/{0}000_{0}999/{1}/{2}.txt', a, SatNum, c )
    command = str.format('./FindTLEforGivenTime.py -n {0} -d {1} -o input.txt {2} {3}', SatelliteNumber, iso, fname1, fname2)
    print command
    os.system(command)

    # 
    # Add start and end times to input.txt
    #
    with open("input.txt", "a") as f:
        f.write(OutputFile+'\n')
        f.write(iso+'\n')
        f.write(iso+'\n')
        f.close()


    # 
    # Add start and end times to input.txt
    #
    if Append:
        os.system('./MagEphemFromTLE -a')
    else:
        os.system('./MagEphemFromTLE')
        Append = True
        
    



