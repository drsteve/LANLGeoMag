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
# go in and find the library and figure out its extension, either dylib or so
# TODO are there more coices that need to be handled?

tmp_name = []
for pth in libpaths:
    tmp_name.append(glob.glob(pth[2:] + '/' + 'libLanlGeoMag.*'))

# we assume it is only in one place (it better be)
while len(tmp_name) > 0:
    libname = tmp_name.pop()
    if len(libname) == 0:
        continue

while len(libname) > 0:
    name = libname.pop()
    if name.split(os.path.extsep)[-2] == '0':
        continue
    if name.split(os.path.extsep)[-1] == 'so':
        library_name = name
    elif name.split(os.path.extsep)[-1] == 'dylib':
        library_name = name

try:
    subprocess.call(['ctypesgen.py',
                     '-l' + library_name,
                     '--no-macro-warnings',
                     '-o', 'Lgm_Wrap.py.bak'] +
                    includes + libpaths +
                    headers)
except OSError:
    print("ctypesgen not installed, rename one of the *_Lgm_Wrap.py.bak files as Lgm_Wrap")
#finally:
#    print('copying ' + 'Lgm_Wrap.py.bak' + ' to ' + 'lgmpy/Lgm_Wrap.py')
#    shutil.copy('Lgm_Wrap.py.bak', 'lgmpy/Lgm_Wrap.py')
