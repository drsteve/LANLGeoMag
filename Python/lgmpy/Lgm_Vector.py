# -*- coding: utf-8 -*-
"""
Vector class for Lgm

@author: Brian Larsen
@organization: LANL
@contact: balarsen@lanl.gov

@version: V1: 20-Dec-2010 (BAL)
"""
import itertools
import copy
from ctypes import pointer, c_double

from Lgm_Wrap import Lgm_Vector, Lgm_VecSub, Lgm_ScaleVector, Lgm_NormalizeVector, \
    Lgm_CrossProduct, Lgm_Magnitude, Lgm_ForceMagnitude, Lgm_DotProduct, \
    Lgm_VecDiffMag, Lgm_VecAdd, Lgm_SphToCartCoords, Lgm_CartToSphCoords

class Lgm_Vector(Lgm_Vector):
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
                raise(TypeError('Bad type: %s in __eq__ comparison' % (type(other)) ))
            else:
                return self == other
        else:
            raise(TypeError('Bad type: %s in __eq__ comparison' % (type(other)) ))

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

        @return: List repersentation of the components
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
        add two vectors together or add a scaler to each component

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
        elif isinstance(other, (int, float, long)):
            x = self.x + other
            y = self.y + other
            z = self.z + other
            return  Lgm_Vector(x, y, z)
        else:
            raise(ArithmeticError("Cannot add type %s to a Lgm_Vector" % (type(other))))

    def __sub__(self, other):
        """
        subtract two vectors or subtract a scaler from each component

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
        elif isinstance(other, (int, float, long)):
            x = self.x - other
            y = self.y - other
            z = self.z - other
            return  Lgm_Vector(x, y, z)
        else:
            raise(ArithmeticError("Cannot subtract %s from a Lgm_Vector" % (type(other))))

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
        elif isinstance(other, (int, float, long)):
            x = self.x * other
            y = self.y * other
            z = self.z * other
            return Lgm_Vector(x, y, z)
        else:
            raise(ArithmeticError("Cannot subtract %s from a Lgm_Vector" % (type(other))))

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
        if isinstance(other, (int, float, long)):
            x = self.x / other
            y = self.y / other
            z = self.z / other
            return Lgm_Vector(x, y, z)
        else:
            raise(ArithmeticError("Cannot subtract %s from a Lgm_Vector" % (type(other))))

    def tolist(self):
        """
        change an Lgm_Vector to a list
        """
        return [self.x, self.y, self.z]

    def crossProduct(self, other):
        """
        compute the cross product of two vectors

        @param other: other vector to cross prod (on the right)
        @type other: Lgm_Vector
        @return: the cross product of the 2 vectors
        @rtype: Lgm_Vector

        @author: Brian Larsen
        @organization: LANL
        @contact: balarsen@lanl.gov

        @version: V1: 22-Dec-2010 (BAL)
        """
        o_vec = Lgm_Vector(0, 0, 0)
        Lgm_CrossProduct( pointer(self), pointer(other), pointer(o_vec) )
        return o_vec

    def magnitude(self):
        """
        return the magnitude of the vector

        @return: Magnitude of the vector
        @rtype: double

        @author: Brian Larsen
        @organization: LANL
        @contact: balarsen@lanl.gov

        @version: V1: 22-Dec-2010 (BAL)
        """
        return Lgm_Magnitude( pointer(self) )

    def normalize(self):
        """
        Normalize the vector in place

        @author: Brian Larsen
        @organization: LANL
        @contact: balarsen@lanl.gov

        @version: V1: 22-Dec-2010 (BAL)
        """
        Lgm_NormalizeVector( pointer(self) )

    def dotProduct(self, other):
        """
        compute the dot product of two vectors

        @param other: other vector to cross prod (on the right)
        @type other: Lgm_Vector
        @return: the dot product of the 2 vectors
        @rtype: Lgm_Vector

        @author: Brian Larsen
        @organization: LANL
        @contact: balarsen@lanl.gov

        @version: V1: 22-Dec-2010 (BAL)
        """
        return Lgm_DotProduct(pointer(self), pointer(other) )

    def scale(self, val):
        """
        Scale a vector by a scalar

        @param val: the value to scale the vector by
        @type val: (int, long, float)

        @author: Brian Larsen
        @organization: LANL
        @contact: balarsen@lanl.gov

        @version: V1: 22-Dec-2010 (BAL)
        """
        Lgm_ScaleVector( pointer(self), val)

    def diffMag(self, other):
        """
        Find Magnitude of difference between to vectors

        @param other: other vector to compare
        @type other: Lgm_Vector
        @return: difference in magnitudes
        @rtype: float

        @author: Brian Larsen
        @organization: LANL
        @contact: balarsen@lanl.gov

        @version: V1: 22-Dec-2010 (BAL)
        """
        return Lgm_VecDiffMag(pointer(self), pointer(other) )

    def forceMagnitude(self, val):
        """
        Force the vector to have a given magnitude inplace

        @param val: the value for the new magnitude
        @type val: (int, long, float)

        @author: Brian Larsen
        @organization: LANL
        @contact: balarsen@lanl.gov

        @version: V1: 22-Dec-2010 (BAL)
        """
        Lgm_ForceMagnitude(pointer(self), val)

def SphToCart(lat, lon, rad):
    """takes an input Lat, Lon, Rad and returns x, y, z"""
    vec1 = Lgm_Vector(0, 0, 0)
    try:
        assert type(lat) == type(lon) == type(rad)
        llat, llon, lrad = len(lat), len(lon), len(rad)
        assert llat == llon == lrad
    except TypeError:
        Lgm_SphToCartCoords(lat, lon, rad, pointer(vec1))
        return vec1
    except AssertionError:
        raise(ValueError('All input must be the same length and type'))
    else:
        ans = []
        for v1, v2, v3 in itertools.izip(lat, lon, rad):
            Lgm_SphToCartCoords(v1, v2, v3, pointer(vec1))
            ans.append(copy.copy(vec1))
        return ans
        
def CartToSph(x, y, z):
    """takes an input x, y, z and returns Lat, Lon, Rad"""
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
        raise(ValueError('All input must be the same length and type'))
    else:
        ans = []
        for v1, v2, v3 in itertools.izip(x, y, z):
            vec1 = Lgm_Vector(v1, v2, v3)
            Lgm_CartToSphCoords(pointer(vec1), pointer(lat), pointer(lon), pointer(rad))
            ans.append(copy.copy([lat.value, lon.value, rad.value]))
        return ans
