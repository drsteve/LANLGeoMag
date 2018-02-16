"""
LanlGeoMag python wrappers
See Mike Henderson about LanlGeoMag and Brian Larsen about the python wrappers


@author: Brian Larsen
@organization: LANL
@contact: balarsen@lanl.gov

@version: V1: 06-Dec-2010 (BAL)
"""
# put modules here that you want to be accessible through 'from spacepy import *'
__all__ = ['Lstar', 'Closed_Field', 'Lgm_T89', 'Lgm_T89c', 'Lgm_OP77']

# Expose definitions from modules in this package.
from lgmpy.Closed_Field import Closed_Field
from Lgm_OP77 import OP77
from Lgm_T89 import T89
from lgmpy.Lstar import get_Lstar
from lgmpy.magcoords import coordTrans

#package info
__version__ = '0.1dev'
__author__ = ['Brian Larsen', 'Steve Morley']
