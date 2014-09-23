#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Master test suite for the package

@author: Brian Larsen
@organization: LANL
@contact: balarsen@lanl.gov

@version: V1: 20-Dec-2010 (BAL)
"""

import sys

try:
    import unittest_pretty as ut
except ImportError:
    import unittest as ut

from lgmpy.test_Lgm_Vector import *
from lgmpy.test_Lgm_CTrans import *
from lgmpy.test_Lgm_MagModelInfo import *
from lgmpy.test_Lgm_T89 import *
from lgmpy.test_Lgm_Wrap import *
from lgmpy.test_Closed_Field import *
from lgmpy.test_Lstar import *
from lgmpy.test_Lgm_OP77 import *
from lgmpy.test_magcoords import *
from lgmpy.test_Bfield_dict import *
from lgmpy.test_Lgm_DateAndTime import *
from lgmpy.test_quicksort import *
from lgmpy.test_Lgm_Utils import *
from lgmpy.test_IsoTimeStringToDateTime import *
# add others here as they exist

if __name__ == '__main__':
    ut.main()
