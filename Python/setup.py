#!/usr/bin/env python

from distutils.core import setup

# this is where extra stuff goes for the install
# like ctypesgen and the quilt patch if I want to go that way

setup(name='lgmpy',
      version='0.0',
      description='Python wrapper for Lgm',
      author='Brian Larsen',
      author_email='balarsen@lanl.gov',
      packages=['.']
     )
