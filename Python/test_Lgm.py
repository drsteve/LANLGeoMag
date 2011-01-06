#!/usr/bin/env python

"""
Master test suite for the package

@author: Brian Larsen
@organization: LANL
@contact: balarsen@lanl.gov

@version: V1: 20-Dec-2010 (BAL)
"""

import unittest_pretty

from test_Lgm_Vector import *
from test_Lgm_CTrans import *
from test_Lgm_DateAndTime import *
#from test_Lgm_T89 import *
# add others here as they exist



if __name__ == '__main__':
    unittest_pretty.main()
