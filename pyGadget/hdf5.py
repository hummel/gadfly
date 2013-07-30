# hdf5.py
# Jacob Hummel
"""
This module contains classes for reading Gadget2 HDF5 snapshot data.
"""
import units
import constants

class HDF5Group(object):
    """
    Base Class for HDF5 groups
    """
    def __init__(self):
        super(HDF5Group,self).__init__()
        
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
        self.ScaleFactor = self.Time
            
class PartTypeX(HDF5Group):
    """
    Class for generic particle info.
    """
    def __init__(self, file_id, ptype):
        super(PartTypeX,self).__init__()
        self._header = Header(file_id)
        group = file_id['PartType'+str(ptype)]
        for item in group.items():
            key = '_'+item[0].replace(' ', '_')
            vars(self)[key] = item[1]

    def load_masses(self, conv=units.Mass_sun, no_h=True):
        """
        Load Particle Masses in units of M_sun (default)
        conv: unit conversion from code units
        no_h: if True, remove dependence on hubble parameter.
        """
        if no_h:
            h = self._header.HubbleParam
            self.masses = self._Masses.value*conv/h
        else:
            self.masses = self._Masses.value*conv

    def get_masses(self, conv=units.Mass_sun, no_h=True):
        """
        Return Particle Masses in units of M_sun (default)
        conv: unit conversion from code units
        no_h: if True, remove dependence on hubble parameter.
        """
        try:
            return self.masses
        except AttributeError:
            self.load_masses(conv,no_h)
            return self.masses

    def load_coords(self, conv=units.Length_kpc, no_h=True, comoving=False):
        """
        Load Particle Coordinates in units of kpc (default)
        conv: unit conversion from code units
        no_h: if True, remove dependence on hubble parameter.
        comoving: if False, remove dependence on scale factor.
        """
        self.coordinates = self._Coordinates.value*conv
        if no_h:
            h = self._header.HubbleParam
            self.coordinates = self.coordinates/h
        if not comoving:
            a = self._header.ScaleFactor
            self.coordinates = self.coordinates*a

    def get_coords(self, conv=units.Length_kpc, no_h=True, comoving=False):
        """
        Return Particle Coordinates in units of kpc (default)
        conv: unit conversion from code units
        no_h: if True, remove dependence on hubble parameter.
        comoving: if False, remove dependence on scale factor.
        """
        try:
            return self.coordinates
        except AttributeError:
            self.load_coords(conv, no_h, comoving)
            return self.coordinates

    def load_velocities(self, conv=units.Velocity_km_s):
        """
        Load Particle Velocities in units of km/s (default)
        conv: unit conversion from code units
        """
        self.velocities = self._Velocities.value*conv

    def get_velocities(self, conv=units.Velocity_km_s):
        """
        Return Particle Velocities in units of km/s (default)
        conv: unit conversion from code units
        """
        try:
            return self.velocities
        except AttributeError:
            self.load_velocities(conv)
            return self.velocities

    def load_PIDs(self):
        """
        Load Particle ID numbers
        """
        self.particleIDs = self._ParticleIDs.value

    def get_PIDs(self):
        """
        Load Particle ID numbers
        """
        try:
            return self.particleIDs
        except AttributeError:
            self.load_PIDs()
            return self.particleIDs

    def load_all(self, no_h=True, comoving=False):
        self.load_masses(no_h=no_h)
        self.load_coords(no_h=no_h, comoving=comoving)
        self.load_velocities()
        self.load_PIDs()
