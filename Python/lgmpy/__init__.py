"""
LanlGeoMag python wrappers
See Mike Henderon about LanlGeoMag and Brian Larsen about the python wrappers


@author: Brian Larsen
@organization: LANL
@contact: balarsen@lanl.gov

@version: V1: 06-Dec-2010 (BAL)
"""
# put modules here that you want to be accessible through 'from spacepy import *'
__all__ = ['Lstar', 'Closed_Field', 'Lgm_T89', 'Lgm_OP77']

# Expose definitions from modules in this package.
from lgmpy import Closed_Field, Lgm_OP77, MagData, Lgm_T89, Lgm_CTrans, Lgm_MagEphemInfo
from lgmpy import Lgm_Wrap, Lgm_MagModelInfo,  Lstar

#package info
__version__ = '0.1dev'
__author__ = ['Brian Larsen', 'Steve Morley']

