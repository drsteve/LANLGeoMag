#!/usr/bin/env python

import glob
import os
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

ctypesgen_worked=False
DEVNULL = os.open(os.devnull, os.O_RDWR)
try:
    subprocess.check_call(['ctypesgen.py',
                           '-l' + library_name,
                           '--no-macro-warnings',
                           '-o', 'Lgm_Wrap.py.bak'] +
                          includes + libpaths +
                          headers,
                          stdout=DEVNULL, stderr=DEVNULL)
    ctypesgen_worked = True
except:
    #This is kinda ridiculous, because it tries to find something
    #not on the binary PATH
    interp = os.path.basename(sys.executable)
    for path in sys.path:
        ctypesgen = os.path.normpath(
            os.path.join(path, '..', '..', '..', 'bin', 'ctypesgen.py'))
        if not os.path.isfile(ctypesgen):
            continue
        try:
            subprocess.check_call([interp, ctypesgen,
                                   '-l' + library_name,
                                   '--no-macro-warnings',
                                   '-o', 'Lgm_Wrap.py.bak'] +
                                  includes + libpaths +
                                  headers,
                                  stdout=DEVNULL, stderr=DEVNULL)
        except:
            continue
        ctypesgen_worked = True
        break
os.close(DEVNULL)

outfile = os.path.join('lgmpy','Lgm_Wrap.py')
if os.path.exists('Lgm_Wrap.py.bak'):
    if not ctypesgen_worked:
        print('ctypesgen not found; using wrapper Lgm_Wrap.py.bak')
    else:
        print('Using ctypesgen-created wrapper Lgm_Wrap.py.bak')
    shutil.copy('Lgm_Wrap.py.bak', outfile)
else:
    if sys.platform == 'linux2':
        print('ctypesgen not found; using wrapper LIN_Lgm_Wrap.py.bak')
        shutil.copy('LIN_Lgm_Wrap.py.bak', outfile)
    elif sys.platform == 'darwin':
        print('ctypesgen not found; using wrapper OSX__Lgm_Wrap.py.bak')
        shutil.copy('OSX_Lgm_Wrap.py.bak', outfile)
    else:
        print("ctypesgen not found and no pre-made wrappers available.\n"
              "Please contact the development team.")

