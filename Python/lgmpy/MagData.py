"""
Overview
========
Data model for magnetic field data.

The SpaceData class is a derived dictionary used to hold information in a comman
format
"""
__author__ = 'Brian Larsen'

import numpy as np
import spacepy.toolbox as tb

from spacepy import datamodel

class MagData(datamodel.SpaceData):
    """
    Class to store data and attributes

    Parameters
    ==========
    args : arguments, optional
        passed directly to datamodel.SpaceData
    kwargs : keyword arguments, optional
        passed directly to datamodel.SpaceData

    See Also
    ========
    spacepy.datamodel.SpaceData
    """
    def __init__(self, *args, **kwargs):
        super(MagData, self).__init__(*args, **kwargs)

    def __repr__(self):
        tb.dictree(self, verbose=True, attrs=True)
        return ''
