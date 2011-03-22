"""
Data model for magnetic field data


@author: Brian Larsen
@organization: LANL
@contact: balarsen@lanl.gov

@version: V1: 21-Mar-2011 (BAL)
"""
import numpy as np
import spacepy.toolbox as tb

from spacepy import datamodel

class MagData(datamodel.SpaceData):
    """
    Class to store data and attributes

    @author: Brian Larsen
    @organization: LANL
    @contact: balarsen@lanl.gov

    @version: V1: 21-Mar-2011 (BAL)
    """
    def __init__(self, *args, **kwargs):
        super(MagData, self).__init__(*args, **kwargs)

    def __repr__(self):
        tb.dictree(self, verbose=True, attrs=True)
        return ''
