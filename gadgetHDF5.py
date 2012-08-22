# gadgetHDF5.py
# Jacob Hummel
"""
This module contains classes for reading Gadget2 HDF5 snapshot data.
"""

import h5py
import units
import constants

class Snapshot:
    """
    Class for Gadget2 HDF5 snapshot files
    """
    def __init__(self, filename):
        self.filename = filename
        self.file_id = h5py.File(filename, 'r')
        self.header = Header(self.file_id)
        self.dm = PartTypeDM(self.file_id)
        self.gas = PartTypeSPH(self.file_id)
        
    def keys(self):
        for key in self.file_id.keys():
            print key
        
    def close(self):
        self.file_id.close()

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
            self.load_coords(conv)
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
        self.load_coordinates(no_h=no_h, comoving=comoving)
        self.load_velocities()
        self.load_PIDs()

class PartTypeDM(PartTypeX):
    """
    Class for Dark Matter particles.
    Available to extend class PartTypeX for Dark Matter specific applications.
    """
    def __init__(self, file_id):
        super(PartTypeDM,self).__init__(file_id,1)

class PartTypeSPH(PartTypeX):
    """
    Class for SPH particles.
    Extends class PartTypeX to include gas physics stuff.
    """
    def __init__(self, file_id):
        super(PartTypeSPH,self).__init__(file_id,0)

    def load_density(self, conv=units.Density_cgs, no_h=True, comoving=False):
        """
        Load Particle Densities in cgs units (default)
        conv: unit conversion from code units
        no_h: if True, remove dependence on hubble parameter.
        comoving: if False, remove dependence on scale factor.
        """
        try:
            del self.ndensity
        except AttributeError:
            pass
        self.density = self._Density.value*conv
        if no_h:
            h = self._header.HubbleParam
            self.density = self.density * h*h
        if not comoving:
            ainv = self._header.Redshift + 1 # 1/(scale factor)
            self.density = self.density * ainv**3

    def get_density(self, conv=units.Density_cgs, no_h=True, comoving=False):
        """
        Return Particle Densities in cgs units (default)
        conv: unit conversion from code units
        no_h: if True, remove dependence on hubble parameter.
        comoving: if False, remove dependence on scale factor.
        """
        try:
            return self.density
        except AttributeError:
            self.load_density(conv, no_h, comoving)
            return self.density

    def load_number_density(self, conv=units.Density_cgs, no_h=True, comoving=False):
        """
        Load Particle Number Densities in cgs units (default)
        conv: unit conversion from code units
        no_h: if True, remove dependence on hubble parameter.
        comoving: if False, remove dependence on scale factor.
        """
        try:
            del self.density
        except AttributeError:
            pass
        self.ndensity = self._Density.value*conv
        self.ndensity = self.ndensity * constants.X_h / constants.m_H
        if no_h:
            h = self._header.HubbleParam
            self.ndensity = self.density * h*h
        if not comoving:
            ainv = self._header.Redshift + 1 # 1/(scale factor)
            self.ndensity = self.ndensity * ainv**3

    def get_number_density(self, conv=units.Density_cgs, no_h=True, comoving=False):
        """
        Return Particle Densities in cgs units (default)
        conv: unit conversion from code units
        no_h: if True, remove dependence on hubble parameter.
        comoving: if False, remove dependence on scale factor.
        """
        try:
            return self.ndensity
        except AttributeError:
            self.load_number_density(conv, no_h, comoving)
            return self.ndensity

    def load_internal_energy(self, conv=units.Energy_cgs/units.Mass_g):
        """
        Load internal particle energies per unit mass in cgs units (default)
        conv: unit conversion from code units
        """
        self.energy = self._InternalEnergy.value*conv


    def get_internal_energy(self, conv=units.Density_cgs/units.Mass_g):
        """
        Return Particle Densities in cgs units (default)
        conv: unit conversion from code units
        """
        try:
            return self.energy
        except AttributeError:
            self.load_internal_energy(conv)
            return self.energy
