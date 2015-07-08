# hdf5.py
# Jacob Hummel
"""
This module contains classes for reading Gadget2 HDF5 snapshot data.
"""
import numpy
from pandas import Series, DataFrame

import units
import constants
import coordinates
import analyze
import visualize
        
class Header(object):
    """
    Class for header information from Gadget2 HDF5 snapshots.
    """
    def __init__(self, file_id):
        group = file_id['Header']
        for key in group.attrs.items():
            vars(self)[key[0]] = key[1]
        self.ScaleFactor = self.Time

    def keys(self):
        """
        Print keys for all items in the header.
        """
        for key in vars(self):
            print key

    def get(self,key):
        """
        Return the value for 'key' in the header.
        """
        return vars(self)[key]

class PartType(DataFrame):
    """
    Class for generic particle info.
    """
    def __init__(self, file_id, ptype, unit, **kwargs):
        group = file_id['PartType'+str(ptype)]
        for item in group.items():
            key = '_'+item[0].replace(' ', '_')
            vars(self)[key] = item[1]
        super(PartType, self).__init__(index=self._ParticleIDs.value)
        self._header = Header(file_id)
        self.units = unit
        self.__init_load_dict__()

    def __getstate__(self):
        result = self.__dict__.copy()
        del result['_Coordinates']
        del result['_ParticleIDs']
        del result['_Velocities']
        del result['_Masses']
        return result

    def __setstate__(self, in_dict):
        self.__dict__ = in_dict
        self.__init_load_dict__()

    def refine_dataset(self, criterion):
        self._drop_ids = self[criterion].index
        self.drop(self._drop_ids, inplace=True)
        print self.index.size, 'particles selected.'

    def load_PIDs(self):
        """
        Load Particle ID numbers
        """
        particleIDs = self._ParticleIDs.value
        if self._drop_ids is not None:
            particleIDs = Series(particleIDs, index=self._ParticleIDs.value)
            self['particleIDs'] = particleIDs.drop(self._drop_ids)
        else:
            self['particleIDs'] = particleIDs

    def get_PIDs(self):
        """
        Return Particle ID numbers
        """
        try:
            return self.particleIDs
        except AttributeError:
            self.load_PIDs()
            return self.particleIDs



    def load_quantity(self, *keys):
        """
        Load key(s) from the list of loadable/derivable keys.
        keys: arbitrary number of keys from the list.
        """
        for key in keys:
            load_func = self._load_dict[key]
            load_func()

    def load_all(self):
        """
        Loads all defined keys in self.loadable_keys
        """
        self.load_quantity(*self.loadable_keys)

    def cleanup(self, *exclude):
        """
        Clean up loaded data to save memory.
        exclude: properties to leave loaded.
        """
        if 'coordinates' in exclude:
            exclude += ('x', 'y', 'z')
        if 'velocities' in exclude:
            exclude += ('u', 'v', 'w')
        to_drop = [key for key in self.columns if key not in exclude]
        self.drop(to_drop, axis=1, inplace=True)

    def load_data(self, *properties, **kwargs):
        """
        Load a selection of keys and refine to highest resolution particles.
        Cleans up unrequested properties loaded for calculation of other 
        requested properties when finished.

        properties: arbitrary number of keys from the list.
        refine (default True): refine to highest resolution particles only.
        stride: If set, take every stride'th particle.
        """
        #load primary quantities first.
        print 'Loading data...'
        for p in properties:
            if p not in self._calculated:
                self.load_quantity(p)
        # load calculated quantities second.
        for p in properties:
            if p in self._calculated:
                self.load_quantity(p)

        stride = kwargs.pop('stride', None)
        if stride:
            for prop in properties:
                vars(self)[prop] = vars(self)[prop][::stride]

        # Cleanup to save memory
        self.cleanup(*properties)
