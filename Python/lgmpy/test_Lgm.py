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

from test_Lgm_Vector import *
from test_Lgm_CTrans import *
from test_Lgm_MagModelInfo import *
from test_Lgm_T89 import *
from test_Lgm_Data import *
from test_Lgm_Wrap import *
from test_Closed_Field import *
from test_Lstar import *
from test_Lgm_OP77 import *
# add others here as they exist

if __name__ == '__main__':
    ut.main()
