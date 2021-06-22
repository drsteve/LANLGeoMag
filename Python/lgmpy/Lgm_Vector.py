# -*- coding: utf-8 -*-
"""
Overview
--------
Vector classes used in LanlGeoMag

This is a wrapper class for the C-structure and associated functions that
operate on the vectors.  There is full coverage for <, >, ==, +, -. /.

All the operations are done in the underlying C library (even though that is
silly for much of this).

"""
from __future__ import division
import sys
import copy
from ctypes import pointer, c_double

import numpy as np

from .Lgm_Wrap import Lgm_Vector, Lgm_VecSub, Lgm_ScaleVector, Lgm_NormalizeVector, \
    Lgm_CrossProduct, Lgm_Magnitude, Lgm_ForceMagnitude, Lgm_DotProduct, \
    Lgm_VecDiffMag, Lgm_VecAdd, Lgm_SphToCartCoords, Lgm_CartToSphCoords, Lgm_Slerp, \
    Lgm_InitSlerp, Lgm_SlerpInfo

__author__ = 'Brian Larsen (Python), Mike Henderson (C) - LANL'

class Lgm_Vector(Lgm_Vector):
    if sys.version_info.major<3:
        _typetup = (int, float, long)
    else:
        _typetup = (int, float)
    def __eq__(self, other):
        """
        if the components of a Vector are equal the vectors are equal
        """
        if isinstance(other, Lgm_Vector):
            if other.x != self.x:
                return False
            elif other.y != self.y:
                return False
            elif other.z != self.z:
                return False
            else:
                return True
        elif isinstance(other, list):
            try:
                other = Lgm_Vector(other[0], other[1], other[2])
            except:
                raise TypeError('Bad type: %s in __eq__ comparison' % (type(other)) )
            else:
                return self == other
        else:
            raise TypeError('Bad type: %s in __eq__ comparison' % (type(other)) )

    def __gt__(self, other):
        """
        if the magnitude is greater the vector is greater
        """
        # done this way because diffMag is not signed
        m1 = self.magnitude()
        m2 = other.magnitude()
        if m1 > m2:
            return True
        else:
            return False

    def __lt__(self, other):
        """
        if the magnitude is greater the vector is greater
        """
        # done this way because diffMag is not signed
        m1 = self.magnitude()
        m2 = other.magnitude()
        if m1 < m2:
            return True
        else:
            return False

    def __le__(self, other):
        """
        if the magnitude is greater the vector is greater
        """
        # done this way because diffMag is not signed
        m1 = self.magnitude()
        m2 = other.magnitude()
        if m1 <= m2:
            return True
        else:
            return False

    def __ge__(self, other):
        """
        if the magnitude is greater the vector is greater
        """
        # done this way because diffMag is not signed
        m1 = self.magnitude()
        m2 = other.magnitude()
        if m1 >= m2:
            return True
        else:
            return False

    def __str__(self):
        """
        print out the Lgm_Vector as a list

        @return: List representation of the components
        @rtype: str

        @author: Brian Larsen
        @organization: LANL
        @contact: balarsen@lanl.gov

        @version: V1: 22-Dec-2010 (BAL)
        """
        return str([self.x, self.y, self.z])

    def __repr__(self):
        """
        print out the Lgm_Vector as a list

        @return: Name of the vector object
        @rtype: str

        @author: Brian Larsen
        @organization: LANL
        @contact: balarsen@lanl.gov

        @version: V1: 22-Dec-2010 (BAL)
        """
        return str([self.x, self.y, self.z])

    def __getattr__(self, other):
        """
        getattr used for easy access to some methods (mag)

        @param other: attr to get from the object (mag)

        @author: Brian Larsen
        @organization: LANL
        @contact: balarsen@lanl.gov

        @version: V1: 22-Dec-2010 (BAL)
        """
        if other == 'mag':
            self.mag = self.magnitude()
        else:
            raise AttributeError(other)

    def __add__(self, other):
        """
        add two vectors together or add a scalar to each component

        @param other: other vector or scalar to add to vector
        @type other: (Lgm_Vector, int, long, float)

        @return: new vector
        @rtype: Lgm_Vector

        @author: Brian Larsen
        @organization: LANL
        @contact: balarsen@lanl.gov

        @version: V1: 22-Dec-2010 (BAL)
        """
        if isinstance(other, Lgm_Vector): # an other vector
            o_vec = Lgm_Vector(0, 0, 0)
            Lgm_VecAdd(pointer(o_vec), pointer(self), pointer(other))
            return o_vec
        elif isinstance(other, self._typetup):
            x = self.x + other
            y = self.y + other
            z = self.z + other
            return  Lgm_Vector(x, y, z)
        else:
            raise ArithmeticError("Cannot add type %s to a Lgm_Vector" % (type(other)))

    def __sub__(self, other):
        """
        subtract two vectors or subtract a scalar from each component

        @param other: other vector or scalar to subtract from vector
        @type other: (Lgm_Vector, int, long, float)

        @return: new vector
        @rtype: Lgm_Vector

        @author: Brian Larsen
        @organization: LANL
        @contact: balarsen@lanl.gov

        @version: V1: 22-Dec-2010 (BAL)
        """
        if isinstance(other, Lgm_Vector): # an other vector
            o_vec = Lgm_Vector(0, 0, 0)
            Lgm_VecSub(pointer(o_vec), pointer(self), pointer(other))
            return o_vec
        elif isinstance(other, self._typetup):
            x = self.x - other
            y = self.y - other
            z = self.z - other
            return  Lgm_Vector(x, y, z)
        else:
            raise ArithmeticError("Cannot subtract %s from a Lgm_Vector" % (type(other)))

    def __mul__(self, other):
        """
        compute the cross product of two vectors or multi a scalar to each component

        @param other: other vector or scalar to cross product or multiply to vector
        @type other: (Lgm_Vector, int, long, float)

        @return: new vector
        @rtype: Lgm_Vector

        @author: Brian Larsen
        @organization: LANL
        @contact: balarsen@lanl.gov

        @version: V1: 22-Dec-2010 (BAL)
        """
        if isinstance(other, Lgm_Vector): # another vector
            return self.crossProduct(other)
        elif isinstance(other, self._typetup):
            x = self.x * other
            y = self.y * other
            z = self.z * other
            return Lgm_Vector(x, y, z)
        else:
            raise ArithmeticError("Cannot subtract %s from a Lgm_Vector" % (type(other)))

    def __div__(self, other):
        """
        divide a scalar into each component

        @param other: scalar to divide into vector
        @type other: (int, long, float)

        @return: new vector
        @rtype: Lgm_Vector

        @author: Brian Larsen
        @organization: LANL
        @contact: balarsen@lanl.gov

        @version: V1: 22-Dec-2010 (BAL)
        """
        if isinstance(other, self._typetup):
            x = self.x / other
            y = self.y / other
            z = self.z / other
            return Lgm_Vector(x, y, z)
        else:
            raise ArithmeticError("Cannot divide an Lgm_Vector by type %s" % (type(other)))

    __truediv__ = __div__

    def __floordiv__(self, other):
        """
        floor divide a scalar into each component

        @param other: scalar to divide into vector
        @type other: (int, long, float)

        @return: new vector
        @rtype: Lgm_Vector

        @author: Steve Morley
        @organization: LANL
        """
        if isinstance(other, self._typetup):
            x = self.x // other
            y = self.y // other
            z = self.z // other
            return Lgm_Vector(x, y, z)
        else:
            raise ArithmeticError("Cannot divide an Lgm_Vector by type %s" % (type(other)))

    def tolist(self):
        """
        change a Lgm_Vector to a list

        Returns
        -------
        out : list
            the contents of the Lgm_Vector changed into a list

        Examples
        --------
        >>> from lgmpy import Lgm_Vector
        >>> dat = Lgm_Vector.Lgm_Vector(1,2,3)
        >>> dat.tolist()
        [1.0, 2.0, 3.0]
        """
        return [self.x, self.y, self.z]

    def toarray(self):
        """
        change a Lgm_Vector to a numpy array

        Returns
        -------
        out : list
            the contents of the Lgm_Vector changed into a list

        Examples
        --------
        >>> from lgmpy import Lgm_Vector
        >>> dat = Lgm_Vector.Lgm_Vector(1,2,3)
        >>> dat.tolist()
        [1.0, 2.0, 3.0]
        """
        return np.asarray([self.x, self.y, self.z])

    def crossProduct(self, other):
        """
        compute the cross product of two vectors

        Parameters
        ----------
        other : other Lgm_Vector to cross prod (on the right)

        Returns
        -------
        out : the cross product of the 2 vectors (Lgm_Vector)

        Examples
        --------
        >>> from lgmpy import Lgm_Vector
        >>> dat = Lgm_Vector.Lgm_Vector(1,2,3)
        >>> dat.crossProduct(dat)
        [0.0, 0.0, 0.0]
        >>> dat2 = Lgm_Vector.Lgm_Vector(3,2,1)
        >>> dat.crossProduct(dat2)
        [-4.0, 8.0, -4.0]
        """
        o_vec = Lgm_Vector(0, 0, 0)
        Lgm_CrossProduct( pointer(self), pointer(other), pointer(o_vec) )
        return o_vec

    def magnitude(self):
        """
        Find the magnitude of a vector

        Returns
        -------
        out : double
            A double value that is the magnitude of the vector

        Examples
        --------
        >>> from lgmpy import Lgm_Vector
        >>> dat = Lgm_Vector.Lgm_Vector(1,2,3)
        >>> dat.magnitude()
        3.74165738...
        """
        return Lgm_Magnitude( pointer(self) )

    def normalize(self):
        """
        Normalize the vector in place. Replaces the contents of self.

        Examples
        --------
        >>> from lgmpy import Lgm_Vector
        >>> import numpy as np
        >>> dat = Lgm_Vector.Lgm_Vector(1,2,3)
        >>> dat.normalize()
        >>> print(dat)
        [0.26726124..., 0.53452248..., 0.801783725...]
        """
        Lgm_NormalizeVector( pointer(self) )

    def dotProduct(self, other):
        """
        compute the dot product of two vectors

        Parameters
        ----------
        other : other Lgm_Vector to dot product

        Returns
        -------
        out : the dot product of the 2 vectors (double)

        Examples
        --------
        >>> from lgmpy import Lgm_Vector
        >>> dat = Lgm_Vector.Lgm_Vector(1,2,3)
        >>> dat.dotProduct(dat)
        14.0
        >>> dat2 = Lgm_Vector.Lgm_Vector(3,2,1)
        >>> dat.dotProduct(dat2)
        10.0
        """
        return Lgm_DotProduct(pointer(self), pointer(other) )

    def scale(self, val):
        """
        Scale a vector by a scalar and replace self

        Parameters
        ----------
        val : the value to scale the vector by (int, long, float)

        Examples
        --------
        >>> from lgmpy import Lgm_Vector
        >>> dat = Lgm_Vector.Lgm_Vector(1,2,3)
        >>> dat.scale(10)
        >>> print(dat)
        [10.0, 20.0, 30.0]
        """
        Lgm_ScaleVector( pointer(self), val)

    def diffMag(self, other):
        """
        Find Magnitude of difference between to vectors

        Parameters
        ----------
        other : another Lgm_Vector

        Returns
        -------
        out : the difference in magnitude (double)

        Examples
        --------
        >>> from lgmpy import Lgm_Vector
        >>> dat = Lgm_Vector.Lgm_Vector(1,2,3)
        >>> dat2 = Lgm_Vector.Lgm_Vector(3,2,1)
        >>> dat.diffMag(dat2)
        2.8284271247461903
        """
        return Lgm_VecDiffMag(pointer(self), pointer(other) )

    def forceMagnitude(self, val):
        """
        Force the vector to have a given magnitude in place

        Parameters
        ----------
        val : the magnitude of the vector

        Examples
        --------
        >>> from lgmpy import Lgm_Vector
        >>> dat = Lgm_Vector.Lgm_Vector(1,2,3)
        >>> dat.forceMagnitude(10)
        >>> print(dat)
        [2.672612..., 5.3452248..., 8.0178372...]
        """
        Lgm_ForceMagnitude(pointer(self), val)

    def slerp(self, other, frac):
        """
        In computer graphics, Slerp is shorthand for spherical linear interpolation, introduced by Ken Shoemake in the
        context of quaternion interpolation for the purpose of animating 3D rotation. It refers to constant-speed
        motion along a unit-radius great circle arc, given the ends and an interpolation parameter between 0 and 1.

        Parameters
        ----------
        other : another Lgm_Vector
        frac  : float
            Fractional distance between the two unit vectors to interpolate
            
        Returns
        -------
        out : Lgm_Vector
            New vector between self and other

        Examples
        --------
        >>> from lgmpy import Lgm_Vector
        >>> v1 = Lgm_Vector.Lgm_Vector(1,2,3)
        >>> v1.normalize() 
        >>> v2 = Lgm_Vector.Lgm_Vector(3,7,1)
        >>> v2.normalize()
        >>> print(v1.slerp(v2, 0)
        []
        >>> print(v1.slerp(v2, 1)
        []
        >>> print(v1.slerp(v2, 0.3)
        []
        """
        frac = float(frac)
        out = Lgm_Vector(0,0,0)
        si = Lgm_SlerpInfo()
        Lgm_InitSlerp(pointer(self), pointer(other), pointer(si))
        Lgm_Slerp( pointer(self), pointer(other), pointer(out), c_double(frac), pointer(si))
        return out
        
def SphToCart(lat, lon, rad):
    """
    takes an input Lat, Lon, Rad and returns x, y, z

    Parameters
    ----------
    lat : latitude (deg)
    lon : longitude (deg)
    rad : radius (Distance units)

    Returns
    -------
    out : Lgm_Vector
        A Lgm_Vector in Cartesian coords same units as the input

    Examples
    --------
    >>> from lgmpy import Lgm_Vector
    >>> Lgm_Vector.SphToCart(0,0,1)
    [1.0, 0.0, 0.0]

    >>> Lgm_Vector.SphToCart(45., 45., 1.)
    [0.5000000..., 0.5, 0.70710678...]
    """

    vec1 = Lgm_Vector(0, 0, 0)
    try:
        assert type(lat) == type(lon) == type(rad)
        llat, llon, lrad = len(lat), len(lon), len(rad)
        assert llat == llon == lrad
    except TypeError:
        Lgm_SphToCartCoords(lat, lon, rad, pointer(vec1))
        return vec1
    except AssertionError:
        raise ValueError('All input must be the same length and type')
    else:
        ans = []
        for v1, v2, v3 in zip(lat, lon, rad):
            Lgm_SphToCartCoords(v1, v2, v3, pointer(vec1))
            ans.append(copy.copy(vec1))
        return ans

def CartToSph(x, y, z):
    """
    takes an input x, y, z and returns Lat, Lon, Rad

    Parameters
    ----------
    x : x coordinate (Distance unit)
    y : y coordinate (Distance unit)
    z : z coordinate (Distance unit

    Returns
    -------
    lat : output latitude (deg)
    lon : output longitude (deg)
    rad : output radius (Distance unit)

    Examples
    --------
    >>> from lgmpy import Lgm_Vector
    >>> Lgm_Vector.CartToSph(0,0,1)
    (90.0, 0.0, 1.0)

    >>> import math
    >>> Lgm_Vector.CartToSph(math.sqrt(2),0.,math.sqrt(2))
    (45.0, 0.0, 2.0)
    """

    lat, lon, rad = c_double(), c_double(), c_double()
    try:
        assert type(x) == type(y) == type(z)
        lx, ly, lz = len(x), len(y), len(z)
        assert lx == ly == lz
    except TypeError:
        vec1 = Lgm_Vector(x, y, z)
        Lgm_CartToSphCoords(pointer(vec1), pointer(lat), pointer(lon), pointer(rad))
        return lat.value, lon.value, rad.value
    except AssertionError:
        raise ValueError('All input must be the same length and type')
    else:
        ans = []
        for v1, v2, v3 in zip(x, y, z):
            vec1 = Lgm_Vector(v1, v2, v3)
            Lgm_CartToSphCoords(pointer(vec1), pointer(lat), pointer(lon), pointer(rad))
            ans.append(copy.copy([lat.value, lon.value, rad.value]))
        return ans
