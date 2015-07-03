# nbody.py
# Jacob Hummel
"""
This module contains classes for reading Gadget2 nbody (dark matter)
particle data.
"""
import hdf5

class PartTypeDM(hdf5.PartTypeX):
    """
    Class for Dark Matter particles.
    Available to extend class PartTypeX for Dark Matter specific applications.
    """
    def __init__(self, file_id, units, **kwargs):
        super(PartTypeDM,self).__init__(file_id,1, units, **kwargs)

