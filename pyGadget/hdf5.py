# hdf5.py
# Jacob Hummel
"""
This module contains classes for reading Gadget2 HDF5 snapshot data.
"""
import numpy
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
    def __init__(self, file_id, ptype, unit):
        super(PartTypeX,self).__init__()
        self._header = Header(file_id)
        self.units = unit
        group = file_id['PartType'+str(ptype)]
        for item in group.items():
            key = '_'+item[0].replace(' ', '_')
            vars(self)[key] = item[1]

        self._load_dict = {'masses':self.get_masses,
                           'coordinates':self.get_coords,
                           'velocities':self.get_velocities,
                           'particleIDs':self.get_PIDs}
        self.loadable_keys = self._load_dict.keys()
        self._calculated = []

    def load_masses(self, unit=None):
        """
        Load Particle Masses in units of M_sun (default set in units class)
        unit: unit conversion from code units
        """
        if unit:
            self.units.set_mass(unit)
        self.masses = self._Masses.value * self.units.mass_conv
        if self.units.remove_h:
            h = self._header.HubbleParam
            self.masses /= h

    def get_masses(self, unit=None):
        """
        Return Particle Masses in units of M_sun (default set in units class)
        unit: unit conversion from code units
        """
        if unit:
            if unit != self.units.mass_unit:
                self.load_masses(unit)
        try:
            return self.masses
        except AttributeError:
            self.load_masses(unit)
            return self.masses

    def load_coords(self, unit=None):
        """
        Load Particle Coordinates in units of kpc (default set in units class)
        unit: unit conversion from code units
        """
        if unit:
            self.units._set_coord_length(unit)
        self.coordinates = self._Coordinates.value * self.units.length_conv
        if self.units.remove_h:
            h = self._header.HubbleParam
            self.coordinates /= h
        if self.units.coordinate_system == 'physical':
            a = self._header.ScaleFactor
            self.coordinates *= a

    def get_coords(self, unit=None):
        """
        Return Particle Coordinates in units of kpc (default set in units class)
        unit: unit conversion from code units
        """
        if unit:
            if unit != self.units._coord_unit:
                self.load_coords(unit)
        try:
            return self.coordinates
        except AttributeError:
            self.load_coords(unit)
            return self.coordinates

    def load_velocities(self, unit=None):
        """
        Load Particle Velocities in units of km/s (default set in units class)
        unit: unit conversion from code units
        """
        if unit:
            self.units.set_velocity(unit)
        self.velocities = self._Velocities.value * self.units.velocity_conv

    def get_velocities(self, unit=None):
        """
        Return Particle Velocities in units of km/s (default set in units class)
        unit: unit conversion from code units
        """
        if unit:
            if unit != self.units.velocity_unit:
                self.load_velocities(unit)
        try:
            return self.velocities
        except AttributeError:
            self.load_velocities(unit)
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
        keys = vars(self).keys()
        for key in keys:
            if key in self.loadable_keys:
                if key not in exclude:
                    del vars(self)[key]

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

        # Refine if desired
        rf = kwargs.pop('refine', True)
        if rf:
            mass = self.get_masses()
            sinks = self.get_sinks()
            minimum = numpy.amin(mass)
            refined = numpy.where((mass <= minimum) | (sinks != 0.))[0]
            for prop in properties:
                # refine only if not already refined.
                if vars(self)[prop].size == mass.size:
                    vars(self)[prop] = vars(self)[prop][refined]
        stride = kwargs.pop('stride', None)
        if stride:
            for prop in properties:
                vars(self)[prop] = vars(self)[prop][::stride]

        # Cleanup to save memory
        self.cleanup(*properties)
