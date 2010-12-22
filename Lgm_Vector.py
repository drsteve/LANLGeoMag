"""
Wrapper for the vector utilites in LanlGeoMag


@author: Brian Larsen
@organization: LANL
@contact: balarsen@lanl.gov

@version: V1: 22-Dec-2010 (BAL)
"""

import _Lgm


class Lgm_Vector(object):
    """
    Lgm_Vector class to wrap up the LanlGeoMag routines and add methods

    @todo: partilly finished as an example

    @ivar x: x-component of the vector
    @type a: double
    @ivar y: y-component of the vector
    @type a: double
    @ivar z: z-component of the vector
    @type a: double
    @ivar _lib: Lgm library access
    @type _lib: _Lgm


    @author: Brian Larsen
    @organization: LANL
    @contact: balarsen@lanl.gov

    @version: V1: 22-Dec-2010 (BAL)
    """
    def __init__(self, x=0.0, y=0.0, z=0.0):
        # get the lib
        self._lib = _Lgm._Lgm().lib
        self.x = x
        self.y = y
        self.z = z

    def magnitude(self):
        """
        return the magnitude of the vector

        @author: Brian Larsen
        @organization: LANL
        @contact: balarsen@lanl.gov

        @version: V1: 22-Dec-2010 (BAL)
        """
