# snapshot.py
# Jacob Hummel
"""
This module contains classes for reading Gadget2 HDF5 snapshot files.
"""
import os
import string
import threading
import Queue
import h5py
import numpy

from hdf5 import Header
from nbody import PartTypeNbody
from sph import PartTypeSPH

class File(object):
    """
    Class for Gadget2 HDF5 snapshot files.
    """
    def __init__(self, sim, filename, **kwargs):
        self.sim = sim
        self.filename = filename
        f = filename.replace('.hdf5','')
        self.number = int(f.split('_')[-1])
        self.file_id = h5py.File(filename, 'r')
        self.header = Header(self.file_id)
        kwargs['refine'] = kwargs.pop('refine_nbody', False)
        self.dm = PartTypeNbody(self.file_id, 1, sim, **kwargs)
        kwargs['refine'] = kwargs.pop('refine_gas', False)
        self.gas = PartTypeSPH(self.file_id, 0, sim, **kwargs)

    def __getstate__(self):
        result = self.__dict__.copy()
        del result['file_id']
        return result

#    def __setstate__(self, in_dict):
#        self.__dict__ = in_dict
#        self.file_id = h5py.File(self.filename, 'r')
#        self.header = hdf5.Header(self.file_id)
        
    def keys(self):
        for key in self.file_id.keys():
            print key
        
    def close(self):
        self.file_id.close()

