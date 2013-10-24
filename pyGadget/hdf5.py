# hdf5.py
# Jacob Hummel
"""
This module contains classes for reading Gadget2 HDF5 snapshot data.
"""
import numpy
import units
import constants
import analyze

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
        if self._refined is not None:
            self.masses = self.masses[self._refined]

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
        if self._refined is not None:
            self.coordinates = self.coordinates[self._refined]

    def calculate_spherical_coords(self, unit=None, **kwargs):
        """
        Load particle positions in spherical coordinates with default units of 
        kpc (default set in units class)
        unit: unit conversion from code units
        """
        if unit:
            if unit != self.units._coord_unit:
                self.load_coords(unit)
        try:
            xyz = self.coordinates
        except AttributeError:
            self.load_coords(unit)
            xyz = self.coordinates

        center = kwargs.get('center', None)
        if center:
            xyz =  analyze.center_box(xyz, center)
        else:
            centering = kwargs.get('centering', None)
            if centering in ['avg','max']:
                try:
                    dens = self.get_number_density()
                except AttributeError:
                    raise KeyError("Cannot density-center dark matter!")
                xyz = analyze.center_box(xyz, density=dens, **kwargs)
            else:
                xyz = analyze.center_box(xyz, **kwargs)

        r,theta,phi = analyze.transform_cart2sph(xyz[:,0],xyz[:,1],xyz[:,2])
        self.spherical_coords = numpy.column_stack((r,theta,phi))

    def get_coords(self, unit=None, **kwargs):
        """
        Return Particle Coordinates in units of kpc (default set in units class)
        unit: unit conversion from code units
        system: coordinate system to use.  Options are cartesian or spherical.
                Note that spherical coordinates require a cartesian coordinate
                'center' point to be calculated by analyze.find_center()
        """
        system = kwargs.pop('system','cartesian')
        if unit:
            if unit != self.units._coord_unit:
                self.load_coords(unit)
                if system == 'spherical':
                    self.calculate_spherical_coordinates(unit, **kwargs)

        if system == 'cartesian':
            try:
                return self.coordinates
            except AttributeError:
                self.load_coords(unit)
                return self.coordinates
        elif system == 'spherical':
            try:
                return self.spherical_coords
            except AttributeError:
                self.calculate_spherical_coords(unit, **kwargs)
                return self.spherical_coords
        else:
            raise KeyError("Coordinate system options are 'cartesian' or 'spherical'")

    def load_velocities(self, unit=None):
        """
        Load Particle Velocities in units of km/s (default set in units class)
        unit: unit conversion from code units
        """
        if unit:
            self.units.set_velocity(unit)
        self.velocities = self._Velocities.value * self.units.velocity_conv
        if self._refined is not None:
            self.velocities = self.velocities[self._refined]

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
        if self._refined is not None:
            self.particleIDs = self.particleIDs[self._refined]

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

        stride = kwargs.pop('stride', None)
        if stride:
            for prop in properties:
                vars(self)[prop] = vars(self)[prop][::stride]

        # Cleanup to save memory
        self.cleanup(*properties)
