# snapshot.py
# Jacob Hummel
"""
This module contains classes for reading Gadget2 HDF5 snapshot files.
"""

import h5py
import hdf5
import nbody
import sph

class load:
    """
    Class for Gadget2 HDF5 snapshot files.
    """
    def __init__(self, filename):
        self.filename = filename
        self.number = int(filename[-8:-5])
        self.file_id = h5py.File(filename, 'r')
        self.header = hdf5.Header(self.file_id)
        self.dm = nbody.PartTypeDM(self.file_id)
        self.gas = sph.PartTypeSPH(self.file_id)
        
    def keys(self):
        for key in self.file_id.keys():
            print key
        
    def close(self):
        self.file_id.close()
