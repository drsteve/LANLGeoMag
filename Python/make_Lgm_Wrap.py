#!/usr/bin/env python

import subprocess
import sys
import glob
import shutil

headers = glob.glob('/usr/local/include/Lgm/*.h')

try:
    subprocess.call(['ctypesgen.py',
                     '-L/usr/local/lib/',
                     '-llibLanlGeoMag.dylib',
                     '--no-macro-warnings',
                     '-o', 'OSX_Lgm_Wrap.py.bak'] +
                     headers)
    subprocess.call(['ctypesgen.py',
                     '-L/usr/local/lib/',
                     '-llibLanlGeoMag.so',
                     '--no-macro-warnings',
                     '-o', 'LIN_Lgm_Wrap.py.bak'] +
                     headers)
except OSError:
    print("ctyoesgen not installed, rename one of the *_Lgm_Wrap.py.bak files as Lgm_Wrap")
finally:
    if sys.platform == 'darwin':
        shutil.copy('OSX_Lgm_Wrap.py.bak', 'lgmpy/Lgm_Wrap.py')
    elif 'linux' in sys.platform:
        shutil.copy('LIN_Lgm_Wrap.py.bak', 'lgmpy/Lgm_Wrap.py')
    else:
        raise(OSError("You are using an OS that I don't understand, sorry you are on your own"))
