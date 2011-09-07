#!/usr/bin/env python

from distutils.core import setup

# create the Lgm_Wrap file so that it will be installed
# TODO can this be done smarter?
#import make_Lgm_Wrap
#import os 
#os.system('python make_Lgm_Wrap.py')

# this is where extra stuff goes for the install
# like ctypesgen and the quilt patch if I want to go that way

setup(name='lgmpy',
      version='0.0',
      description='Python wrapper for Lgm',
      author='Brian Larsen',
      author_email='balarsen@lanl.gov',
      packages=['lgmpy']
     )
