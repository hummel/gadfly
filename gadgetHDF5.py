# gadgetHDF5.py
# Jacob Hummel
"""
This module contains classes for reading Gadget2 HDF5 snapshot data.
"""

import h5py
import units

class Snapshot:
    """
    Class for Gadget2 HDF5 snapshot files
    """
    def __init__(self, filename):
        self.filename = filename
        self.file_id = h5py.File(filename, 'r')
        self.header = Header(self.file_id)
        self.gas = PartTypeX(self.file_id, 0)
        self.dm = PartTypeX(self.file_id, 1)
        
    def keys(self):
        for key in self.file_id.keys():
            print key
        
    def close(self):
        self.file_id.close()

class HDF5Group:
    """
    Base Class for HDF5 groups
    """
    def keys(self):
        """
        Print keys for all arrays in the hdf5 group.
        """
        for key in vars(self):
            print key

    def get(self,key):
        """
        Return the raw values of the array from the HDF5 group.
        """
        return vars(self)[key].value
        
class Header(HDF5Group):
    """
    Class for header information from Gadget2 HDF5 snapshots.
    """
    def __init__(self, file_id):
        group = file_id['Header']
        for key in group.attrs.items():
            vars(self)[key[0]] = key[1]
            
class PartTypeX(HDF5Group):
    """
    Class for generic particle info.
    """
    def __init__(self, file_id, ptype):
        group = file_id['PartType'+str(ptype)]
        for item in group.items():
            key = item[0].replace(' ', '_')
            vars(self)[key] = item[1]

    def getMass(self, conv=units.Mass_sun):
        """
        Default: Return Particle Masses in units of M_sun/h
        conv: use the given unit conversion from code units
        """
        return self.Masses.value*conv

    def getCoords(self, conv=units.Length_kpc):
        """
        Default: Return Particle Coordinates in units of kpc/h
        conv: use the given unit conversion from code units
        """
        return self.Coordinates.value*conv

    def getVelocities(self, conv=units.Velocity_km_s):
        """
        Default: Return Particle Velocities in units of km/s
        conv: use the given unit conversion from code units
        """
        return self.Velocities.value*conv

    def getPIDs(self):
        """
        Return Particle ID numbers
        """
        return self.ParticleIDs.value
