"""
Main Lgm Library file, contains all the fucntion definitions


@author: Brian Larsen
@organization: LANL
@contact: balarsen@lanl.gov

@version: V1: 20-Dec-2010 (BAL)
"""

# import everything here then in each file just import as needed
import ctypes
import sys

if sys.platform == 'win32':
    raise(NotImplementedError("Sorry windows isn't up and running"))
elif sys.platform == 'darwin':
    libname = 'libLanlGeoMag.dylib'
else:
    libname = 'libLanlGeoMag.so'

lib = ctypes.CDLL('/usr/local/lib/' + libname )
