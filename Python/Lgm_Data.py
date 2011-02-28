import numpy

class DataArray(numpy.ndarray):
    """
    Container for data within a object (from GPSbase)

    @author: Brian Larsen
    @organization: Los Alamos National Lab
    @contact: balarsen@lanl.gov

    @version: V1: 19-Nov-2010
    @version: V2: 24-Nov-2010 (S.Morley) changed to subclass ndarray
    """
    def __new__(cls, input_array, attrs={}):
        # Input array is an already formed ndarray instance
        # We first cast to be our class type
        obj = numpy.asarray(input_array).view(cls)
        # add the new attribute to the created instance
        obj.attrs = attrs
        # Finally, return the newly created object:
        return obj

    def __array_finalize__(self, obj):
        # see InfoArray.__array_finalize__ for comments
        if obj is None:
            return
        # TODO check to see if this line is needed
        self.attrs = getattr(obj, 'attrs', {})



class Lgm_Data(dict):
    """
    Base Lgm data storage (data model) class
    """
    def __repr__(self):
        """
        Abstract method, reimplement
        """
        raise(NotImplementedError("Abstract method called, reimplement __repr__"))
