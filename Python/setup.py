#!/usr/bin/env python

from distutils.core import setup
from distutils.command.build import build as _build
import distutils.dep_util
import glob
import os
import os.path
import subprocess
import sys


class build(_build):
    """Extends base build to make the ctypesgen wrapper"""

    def make_wrapper(self):
        """Generate the wrapper using ctypesgen"""
        if sys.version_info[0] == 3:
            print('Python 3 is not expected to work, trying anyway...')
        if not 'LIBDIR' in os.environ:
            raise RuntimeError('LIBDIR environment variable not set. '
                               'Use Makefile to build the python wrapper.')
        outlog = os.path.join(self.build_temp, 'ctypesgen.log')
        if not os.path.isdir(self.build_temp):
            os.makedirs(self.build_temp)
        outdir = os.path.join(self.build_lib, 'lgmpy')
        if not os.path.isdir(outdir):
            os.makedirs(outdir)
        outfile = os.path.join(outdir, 'Lgm_Wrap.py')
        libsrcdir = os.path.join('..', 'libLanlGeoMag')
        if not distutils.dep_util.newer_group(
            #Makefile will change if LIBDIR changes
            [os.path.join(libsrcdir, 'libLanlGeoMag.la'), 'Makefile'],
            outfile):
            return #wrapper up to date
        cmd = ['ctypesgen.py',
               '-lLanlGeoMag',
               '--no-macro-warnings',
               '-o', outfile,
               '--compile-libdir=' + os.path.join(libsrcdir, '.libs'),
               '--runtime-libdir=' + os.environ['LIBDIR'],
               '--includedir=' + os.path.join(libsrcdir, 'Lgm'),
               '--includedir=' + libsrcdir,
               ] + \
               list(glob.glob(os.path.join(libsrcdir, 'Lgm', '*.h')))
        print('running ctypesgen to create wrapper')
        FDLOG = os.open(outlog, os.O_WRONLY | os.O_CREAT)
        try:
            subprocess.check_call(cmd, stdout=FDLOG, stderr=FDLOG)
        except:
            raise
            raise RuntimeError(
                'Wrapper generation failed, see ' + outlog + ' for details.')
        finally:
            os.close(FDLOG)
        #TODO: Technically this should compile the wrapper .py to .pyc
        #if that's been asked for (look at build_py in distutils)
        
    def run(self):
        """Perform the build"""
        _build.run(self)
        self.make_wrapper()


setup(name='lgmpy',
      version='0.0',
      description='Python wrapper for Lgm',
      author='Brian Larsen',
      author_email='balarsen@lanl.gov',
      packages=['lgmpy'],
      cmdclass={'build': build,
                },
     )
