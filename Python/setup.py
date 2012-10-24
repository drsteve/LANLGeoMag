#!/usr/bin/env python

from distutils.core import setup
from distutils.command.build import build as _build
from distutils.command.install import install as _install
import os 


class build(_build):
    """Extends base build to make the ctypesgen wrapper"""

    def make_wrapper(self):
        """Generate the wrapper using ctypesgen"""
        #this should dump the output directly in the build dir, NOT including the pyc
        #or put it in a temp dir since install needs to move it over??
        #should check the generated wrapper against the age of the Lgm library to avoid thrash
        pass

    def run(self):
        """Perform the build"""
        _build.run(self)
        self.make_wrapper()


class install(_install):
    """Rewrite the ctypesgen wrapper and install it"""

    def rewrite_wrapper(self):
        """Point the wrapper at the final install dir, put it in build (compiled)"""
        pass

    def run(self):
        """Perform the install"""
        self.rewrite_wrapper() #must go first b/c it puts things in build
        _install.run(self)


setup(name='lgmpy',
      version='0.0',
      description='Python wrapper for Lgm',
      author='Brian Larsen',
      author_email='balarsen@lanl.gov',
      packages=['lgmpy'],
      cmdclass={'build': build,
                'install': install,
                },
     )
