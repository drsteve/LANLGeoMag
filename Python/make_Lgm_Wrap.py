#!/usr/bin/env python

import subprocess
import sys

if sys.platform == 'darwin':
    subprocess.call(['ctypesgen.py',
                     '-l/usr/local/lib/libLanlGeoMag.dylib',
                     '/usr/local/include/Lgm/*.h',
                     '--no-macro-warnings',
                     '-o', 'OSX_Lgm_Wrap.py'])
else:
    subprocess.call(['ctypesgen.py',
                     '-L/usr/local/lib/',
                     '/usr/local/include/Lgm/*.h',
                     '--no-macro-warnings',
                     '-llibLanlGeoMag.so',
                     '-o Linux_Lgm_Wrap.py'])
