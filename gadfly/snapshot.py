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
        self.define_ptype('dm', 1, PartTypeNbody, **kwargs)
        kwargs['refine'] = kwargs.pop('refine_gas', False)
        self.define_ptype('gas', 0, PartTypeSPH, **kwargs)

    def __getstate__(self):
        result = self.__dict__.copy()
        self.close()
        del result['file_id']
        return result

    def define_ptype(self, groupname, ptype, ptype_class, **kwargs):
        vars(self)[groupname] = ptype_class(self.file_id, ptype, self.sim, **kwargs)
        
    def keys(self):
        for key in self.file_id.keys():
            print key
        
    def close(self):
        try:
            self.file_id.close()
        except(ValueError):
            print 'File already closed.'
            pass
