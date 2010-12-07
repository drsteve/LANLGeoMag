
import unittest
import Lgm_Vec


class TimeTests(unittest.TestCase):
    """Tests related to converting times"""

    def setUp(self):
        """Load the leap second kernel"""
        super(TimeTests, self).setUp()
        naif.spice.load_kern(lsk_kern)
        naif.spice.load_kern(sclk_kern)

    def tearDown(self):
        """Unload the leap second kernel"""
        super(TimeTests, self).tearDown()
        naif.spice.load_kern(sclk_kern)
        naif.spice.unload_kern(lsk_kern)

    def testZeroEpoch(self):
        """Test the UT for ET=0"""
        et = naif.spice.time_to_ET(
            datetime.datetime(2000, 1, 1, 11, 58, 55, 816000))
        self.assertAlmostEqual(0.0, et, places=3)



if __name__ == '__main__':
    unittest.main()



#        lgm.Lgm_Magnitude.argtypes=[ctypes.POINTER(Lgm_Vector)]
#lgm.Lgm_Magnitude.restype = ctypes.c_double
#vec = Lgm_Vector(1., 2., 3.)
#ans = lgm.Lgm_Magnitude(vec)
#
#if ans != 3.74165738677:
#    print("Answers not equal")
#    print("Got %f" % (ans))
