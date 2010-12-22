#!/usr/bin/env python

"""
Test suite for the Lgm_Vector file

@author: Brian Larsen
@organization: LANL
@contact: balarsen@lanl.gov

@version: V1: 20-Dec-2010 (BAL)
"""


import unittest
import _Lgm_CTrans
import _Lgm_Vector
import _Lgm
import math



# Tests not written for these yet
#'Lgm_VecSub': [None, Lgm_VectorP, Lgm_VectorP, Lgm_VectorP],
#'Lgm_VecAdd': [None, Lgm_VectorP, Lgm_VectorP, Lgm_VectorP],
#'Lgm_VecDiffMag': [LgmDouble, Lgm_VectorP, Lgm_VectorP],
#'Lgm_ForceMagnitude': [None, Lgm_VectorP, LgmDouble],
#'Lgm_MatTimesVec': [None, LgmDouble * 3 * 3, Lgm_VectorP, Lgm_VectorP],
#'Lgm_Transpose' : [None, LgmDouble * 3 * 3, LgmDouble * 3 * 3],
#'Lgm_MatTimesMat' : [None, LgmDouble * 3 * 3, LgmDouble * 3 * 3, LgmDouble * 3 * 3],


if __name__ == '__main__':
    unittest.main()
