#!/usr/bin/env python

import glob
import os.path
import shutil
import subprocess
import sys

try:
    p = subprocess.Popen(['pkg-config', '--cflags-only-I', 'lgm'],
                         stdout=subprocess.PIPE)
    output, errors = p.communicate()
    assert p.returncode == 0
    includes = output.split()
    p = subprocess.Popen(['pkg-config', '--libs-only-L', 'lgm'],
                         stdout=subprocess.PIPE)
    output, errors = p.communicate()
    assert p.returncode == 0
    libpaths = output.split()
except:
    includes = ['-I/usr/local/include/Lgm']
    libpaths = ['-L/usr/local/lib']
headerpath = next((i[2:] for i in includes if i[-3:] == 'Lgm'))
headers = glob.glob(os.path.join(headerpath, '*.h'))
#Some includes in Lgm have the Lgm/, so include level above Lgm
includes.extend((i[:-3] for i in includes if i[-3:] == 'Lgm' ))

try:
    subprocess.call(['ctypesgen.py',
                     '-llibLanlGeoMag.dylib',
                     '--no-macro-warnings',
                     '-o', 'OSX_Lgm_Wrap.py.bak'] +
                    includes + libpaths + 
                    headers)
    subprocess.call(['ctypesgen.py',
                     '-llibLanlGeoMag.so',
                     '--no-macro-warnings',
                     '-o', 'LIN_Lgm_Wrap.py.bak'] +
                    includes + libpaths + 
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
