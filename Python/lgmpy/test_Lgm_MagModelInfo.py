#!/usr/bin/env python

"""
Test suite for the Lgm_CTrans file <<This is an important one>>

"""

import ctypes
import unittest

from lgmpy import Lgm_MagModelInfo
from lgmpy.Lgm_Wrap import size_MagModelInfo, TRUE, LGM_ABSOLUTE_JUMP_METHOD, DQAGS, \
    DQAGP


class Lgm_MagModelInfoTests(unittest.TestCase):
    """
    Tests related to Lgm_MagModelInfo
    """

    def setUp(self):
        super(Lgm_MagModelInfoTests, self).setUp()

    def tearDown(self):
        super(Lgm_MagModelInfoTests, self).tearDown()

    def test_size_Lgm_MagModelInfo(self):
        """for Lgm_MagModelInfo the c and python sizes should match"""
        self.assertEqual(size_MagModelInfo(),
                         ctypes.sizeof(Lgm_MagModelInfo.Lgm_MagModelInfo) )

    def test_init_vals(self):
        """test that the defaults in init don't change (regression)"""
        mmi = Lgm_MagModelInfo.Lgm_MagModelInfo()
        self.assertEqual(mmi.Kp, 2)
        self.assertEqual(mmi.P, 2.1)
        self.assertEqual(mmi.nFunc, 0)
        self.assertEqual(mmi.B0, 1.00)
        self.assertEqual(mmi.B1, 0.8100)
        self.assertEqual(mmi.B2, 0.4065)
        self.assertEqual(mmi.SavePoints, 0)
        self.assertEqual(mmi.Hmax, 1.0)

        self.assertEqual(mmi.UseInterpRoutines, TRUE)
        self.assertEqual(mmi.Lgm_I_integrand_JumpMethod, LGM_ABSOLUTE_JUMP_METHOD)
        #
        #  Some inits for MagStep
        #
        # self.assertEqual(mmi.Lgm_MagStep_FirstTimeThrough, TRUE) # this went away in 3e24b5ae65a11ef9376fd4246f5eee3a68485d39
        # self.assertEqual(mmi.Lgm_MagStep_eps_old, -1.0) # this went away in 3e24b5ae65a11ef9376fd4246f5eee3a68485d39
        #
        # Set some default tolerances
        #
        self.assertEqual(mmi.Lgm_MagFlux_Integrator_epsabs, 0.0)
        self.assertEqual(mmi.Lgm_MagFlux_Integrator_epsrel, 1e-5)

        self.assertEqual(mmi.Lgm_LambdaIntegral_Integrator_epsabs, 0.0)
        self.assertEqual(mmi.Lgm_LambdaIntegral_Integrator_epsrel, 1e-3)

        self.assertEqual(mmi.Lgm_I_Integrator_epsrel, 0.0)
        self.assertEqual(mmi.Lgm_I_Integrator_epsabs, 1e-3)
        self.assertEqual(mmi.Lgm_I_Integrator, DQAGS)

        self.assertEqual(mmi.Lgm_Sb_Integrator_epsrel, 0.0)
        self.assertEqual(mmi.Lgm_Sb_Integrator_epsabs, 1e-4)
        self.assertEqual(mmi.Lgm_Sb_Integrator, DQAGP)

        self.assertEqual(mmi.Lgm_FindBmRadius_Tol, 1e-10)
        self.assertEqual(mmi.Lgm_FindShellLine_I_Tol, 1e-3)
        self.assertEqual(mmi.Lgm_TraceToMirrorPoint_Tol, 1e-7)




if __name__ == '__main__':
    unittest.main()
