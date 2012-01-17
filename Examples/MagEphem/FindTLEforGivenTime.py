#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
from optparse import OptionParser
import datetime

# dict of parse methods
ParseMethods = {'UseCommonName': "UseCommonName", 'UseSatelliteNumber': "UseSatelliteNumber", 
                    'UseIntDesignator': "UseIntDesignator", 'None': None}

def findTLEinfiles(args, **kwargs):
    defaults = dict(CommonName='', SatelliteNumber='', IntDesignator='',
                    ParseMethod=ParseMethods['None'], OutputFile='',
                    DumpToFile=False, PurgeDuplicates=False)
    #replace missing kwargs with defaults
    for dkey in defaults:
        if dkey not in kwargs:
            kwargs[dkey] = defaults[dkey]
    
    ParseMethod = kwargs['ParseMethod']
    
    # make a list of all the lines that we want to scan through This supports
    # globbing (i.e. if user enters multiple files or *.txt or some such thing, we
    # will read all the files before processing them).
    lines = []
    if args.__len__() > 0:
        # If there are args provided, assume they are the names of files we want to
        # scan through. The shell will already have expanded wildcards to multiple
        # entries here -- we just have to process them all.
        for filename in args:
            try:
                f = open(filename, 'r')
                lines += f.readlines()
                f.close()
            except:
                print 'Unable to read file: ', filename
    else:
        # If there are no args provided, assume input is coming from the stdin
        # (e.g. via a pipe like: cat *.txt | <program> )
        lines = sys.stdin.readlines()

    # Select output stream to dump to
    # Default is stdout, but if an OutputFile is specified, try to open it for
    # writing.
    if kwargs['DumpToFile']:
        try:
            f = open( kwargs['OutputFile'], "w" )
            sys.stdout = f
        except:
            raise IOError('Unable to open output file: %s' % OutputFile)

    # Scan off sets of 3 lines and test to see if the TLE is one we are looking
    # for.
    TLEs = []
    e = enumerate(lines)
    for i, line in e:
        Line0 = line            # Line 0 of a "three-line" element set
        Line1 = e.next()[1]     # Line 1 of a "three-line" element set
        Line2 = e.next()[1]     # Line 2 of a "three-line" element set

        # Add it to the list of TLEs if its one we are interested in For the
        # purposes of sorting later, we will place the lines slightly out of order.
        # Instead of putting Line0 as the first element, lets put it in as the
        # last.
        if ( ((ParseMethod == ParseMethods['UseCommonName']) and (CommonName.upper() in Line0.upper())) \
            or ((ParseMethod == ParseMethods['UseSatelliteNumber']) and (kwargs['SatelliteNumber'] in Line1[2:8])) \
            or ((ParseMethod == ParseMethods['UseIntDesignator']) and (IntDesignator in Line1[9:17]) )):
            # In some cases, there appears to be a formatting problem in the epoch field. For example,
            # one of the TLEs for GOES13 has "07 95.94735520" for the epoch when it should be
            # "07095.94735520". I.e. it seems that the ddd fields may have been printed without leading
            # zeros. We need to fix this now or later sorting may get screwed up. 
            Line1 = Line1[0:18] + Line1[18:32].replace(' ', '0') + Line1[32:69]
            TLE  = [Line1.strip(), Line2.strip(), Line0.strip()]
            TLEs.append( TLE )

    # Remove duplicates entries. If the list elements were hashable, we could just
    # do TLEs = list( set(TLEs) ), but the entries are themselves list which
    # apparently arent hashable. So we do it in a more brute force way. Also, the
    # sort order and the order we scan through (end to start) should ensure that of
    # duplicates, the ones that are retained have the most info in Line0
    # E.g., if we had:
    #
    #   GOES 13 (MORE INFO)
    #   1 29155U 06018A   10297.41394480 -.00000241  00000-0  10000-3 0  3276
    #   2 29155 000.1844 075.0961 0003327 026.7415 005.1518 01.00278452 16218
    #   GOES 13
    #   1 29155U 06018A   10297.41394480 -.00000241  00000-0  10000-3 0  3276
    #   2 29155 000.1844 075.0961 0003327 026.7415 005.1518 01.00278452 16218
    #
    # the entry with "MORE INFO" would be the one that gets retained.
    # Since all the TLEs should be for the same object at this point, sorting on Line1
    # kk
    if kwargs['PurgeDuplicates']:
        if TLEs:
            TLEs.sort()
            t = TLEs[-1]
            for i in range( len(TLEs)-2, -1, -1):
                if t[0:2] == TLEs[i][0:2]:
                    # remove if Line1 and Line2 are the same (i.e. Line0 can differ)
                    # recall that here the lines are stored as [line1, line2, line0]
                    del TLEs[i]
                else:
                    t = TLEs[i]

    # Now, scan through list and find the nearest match. Set t to the first one and
    # use that as a default in the event that all of the TLEs have epoch times
    # greater than our target time.
    t = TLEs[0] # default TLE to return  is the first one in the sorted list
    for i in range( len(TLEs)-1, -1, -1):
        t1 = TLEs[i][0]        # Line1 (its in the zeroeth slot)
        y = int( t1[18:20] )   # 2-digit epoch year
        e = float( t1[18:32] ) # epoch in yyddd.ffffff format

        # convert to 4-digit year (assumes any 2-digit year > 50 is in the 1900s
        # otherwise in the 2000s)
        if y > 50: Epoch = e+1900000.0
        else: Epoch = e+2000000.0

        if Epoch <= kwargs['TargetEpoch']:
            t = TLEs[i]
            break
    
    if kwargs['Verbose']:
        print t[2]
        print t[0]
        print t[1]

    # close OutputFile if we open one.
    if kwargs['DumpToFile']:
        f.close()
        
    return t[2], t[0], t[1]

if __name__ == '__main__':
    
    kwargs = {}
    
    # Define a command-line option parser and add the options we need
    parser = OptionParser(  usage="%prog [options] [filenames to scan]\n"\
                        "       cat [filenames to scan] | %prog [options]", 
                            version="%prog Version 1.00 (October 27, 2010)"  )

    parser.add_option("-d", "--ISO_DateTime",      dest="ISO_DateTime",    
                            help="Date the user wants a TLE for in ISO 8601 format.", 
                            metavar="YYYY-MM-DDTHH:MM:SS")

    parser.add_option("-c", "--CommonName",   dest="CommonName", 
                            help="Common name of object to filter for.",  
                            metavar="COMMON_NAME")

    parser.add_option("-n", "--SatNumber",    dest="SatNumber",  
                            help="Satellite (or object) number (including trailing U "\
                                "or S - for example 29155U).", 
                            metavar="SAT_NUMBER")

    parser.add_option("-i", "--ID",           dest="ID",         
                            help="International Designator (for example 06018A).",  
                            metavar="ID")

    parser.add_option("-o", "--OutputFile",   dest="OutputFile", 
                            help="Filename to dump output to. Stdout is assumed if "\
                                "no file given.",    
                            metavar="OUTFILE")
    
    parser.add_option("-v", "--Verbose",   dest="Verbose", 
                            help="Give verbose output [default=True] ",    
                            metavar="OUTFILE")

    # Parse the args that were (potenitally) given to us on the command line
    (options, args) = parser.parse_args()

    # Set variables based on what we got on the command line
    if options.CommonName != None:
        kwargs['ParseMethod'] = ParseMethods['UseCommonName']
        kwargs['CommonName']  = options.CommonName
    elif options.SatNumber != None:
        kwargs['ParseMethod']    = ParseMethods['UseSatelliteNumber']
        kwargs['SatelliteNumber'] = options.SatNumber
    elif options.ID != None:
        kwargs['ParseMethod']   = ParseMethods['UseIntDesignator']
        kwargs['IntDesignator'] = options.ID

    if options.OutputFile != None:
        kwargs['OutputFile'] = options.OutputFile
        kwargs['DumpToFile'] = True
    if options.Verbose != None:
        kwargs['Verbose'] = options.Verbose
    else:
        kwargs['Verbose'] = True

    kwargs['PurgeDuplicates'] = True

    if options.ISO_DateTime == None:
        print 'Must provide a date/time via -d option\n'
        parser.print_help()
        exit()
    else:
        dt = datetime.datetime.strptime(options.ISO_DateTime, "%Y-%m-%dT%H:%M:%S")
        kwargs['TargetEpoch'] = int(dt.strftime('%Y%j')) + dt.hour/24.0 + dt.minute/1440.0 + dt.second/86400.0
        
    # If ParseMethod is still None, then there was a probable getting command-line
    # args -- bail.
    if not kwargs['ParseMethod']:
        parser.print_help()
        raise ValueError('Parse options invalid') 
    
    Line0, Line1, Line2 = findTLEinfiles(args, **kwargs)
