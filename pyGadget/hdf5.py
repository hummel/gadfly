# hdf5.py
# Jacob Hummel
"""
This module contains classes for reading Gadget2 HDF5 snapshot data.
"""
import numpy
from pandas import DataFrame

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

            
class PartTypeX(DataFrame):
    """
    Class for generic particle info.
    """
    def __init__(self, file_id, ptype, unit):
        group = file_id['PartType'+str(ptype)]
        for item in group.items():
            key = '_'+item[0].replace(' ', '_')
            vars(self)[key] = item[1]
        super(PartTypeX,self).__init__()
        self._header = Header(file_id)
        self.units = unit
        self.__init_load_dict__()

    def __getstate__(self):
        result = self.__dict__.copy()
        del result['_Coordinates']
        del result['_ParticleIDs']
        del result['_Velocities']
        del result['_Masses']
        del result['_load_dict']
        del result['loadable_keys']
        del result['_calculated']
        return result

    def __setstate__(self, in_dict):
        self.__dict__ = in_dict
        self.__init_load_dict__()

    def __init_load_dict__(self):
        self._load_dict = {'masses':self.load_masses,
                           'coordinates':self.load_coords,
                           'velocities':self.load_velocities,
                           'particleIDs':self.load_PIDs}
        derived = {'spherical_coords':self.calculate_spherical_coords,
                   'cylindrical_coords':self.calculate_cylindrical_coords}
        self._load_dict.update(derived)
        self.loadable_keys = self._load_dict.keys()
        self._calculated = derived.keys()

    def load_PIDs(self):
        """
        Load Particle ID numbers
        """
        self['particleIDs'] = self._ParticleIDs.value
        if self._refined is not None:
            self.particleIDs = self.particleIDs[self._refined]

    def get_PIDs(self):
        """
        Return Particle ID numbers
        """
        try:
            return self.particleIDs
        except AttributeError:
            self.load_PIDs()
            return self.particleIDs

    def load_masses(self, unit=None):
        """
        Load Particle Masses in units of M_sun (default set in units class)
        unit: unit conversion from code units
        """
        if unit:
            self.units.set_mass(unit)
        self['masses'] = self._Masses.value * self.units.mass_conv
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
        xyz = self._Coordinates.value * self.units.length_conv
        if self.units.remove_h:
            h = self._header.HubbleParam
            xyz /= h
        if self.units.coordinate_system == 'physical':
            a = self._header.ScaleFactor
            xyz *= a
        if self._refined is not None:
            xyz = xyz[self._refined]
        self['coord_x'] = xyz[:,0]
        self['coord_y'] = xyz[:,1]
        self['coord_z'] = xyz[:,2]

    def load_velocities(self, unit=None):
        """
        Load Particle Velocities in units of km/s (default set in units class)
        unit: unit conversion from code units
        """
        if unit:
            self.units.set_velocity(unit)
        vxyz = self._Velocities.value * self.units.velocity_conv
        if self.units.coordinate_system == 'physical':
            a = self._header.ScaleFactor
            vxyz *= numpy.sqrt(a)
        if self._refined is not None:
            vxyz = vxyz[self._refined]
        self['velocity_x'] = vxyz[:,0]
        self['velocity_y'] = vxyz[:,1]
        self['velocity_z'] = vxyz[:,2]

    def orient_box(self, **kwargs):
        """
        Center and rotate box coordinates AND velocities according to 
        received kwargs 'center' and 'view'. If 'center' unspecified,
        looks for 'centering' kwargs and attempts to auto-center the box.
        """
        try:
            xyz = self.coordinates
        except AttributeError:
            self.load_coords()
            xyz = self.coordinates
        try:
            vxyz = self.velocities
        except AttributeError:
            self.load_velocities()
            vxyz = self.velocities

        center = kwargs.get('center', None)
        vcenter = kwargs.get('vcenter', None)
        centering = kwargs.get('centering', None)
        view = kwargs.get('view', None)
        dlim = kwargs.pop('dens_lim', 1e9)

        if center:
            xyz =  analyze.center_box(xyz, center)
        else:
            if centering in ['avg','max']:
                try:
                    dens = self.get_number_density()
                except AttributeError:
                    raise KeyError("Cannot density-center dark matter!")
                xyz = analyze.center_box(xyz, density=dens, **kwargs)
            elif centering == 'box':
                xyz = analyze.center_box(xyz, **kwargs)

        if vcenter:
            vxyz =  analyze.center_box(vxyz, vcenter)
        else:
            if centering in ['avg','max']:
                try:
                    dens = self.get_number_density()
                except AttributeError:
                    raise KeyError("Cannot density-center dark matter!")
                vxyz = analyze.center_box(vxyz, density=dens, **kwargs)
            else:
                print "Warning: Not Re-centering particle Velocities!"

        if view:
            if view == 'face':
                try:
                    dens = self.get_number_density()
                except AttributeError:
                    raise KeyError("Cannot density-center dark matter!")
                xyz, vxyz = visualize.set_view(view, xyz, velocity=vxyz,
                                               density=dens, dens_lim=dlim)
            else:
                xyz, vxyz = visualize.set_view(view, xyz, velocity=vxyz)

        self.coordinates = xyz
        self.velocities = vxyz

    def calculate_spherical_coords(self, c_unit=None, v_unit=None, **kwargs):
        """
        Load particle positions in spherical coordinates with default units of
        kpc (default set in units class)
        unit: unit conversion from code units
        """
        if (c_unit or v_unit):
            if (c_unit != self.units._coord_unit
                or v_unit != self.units.velocity_unit):
                self.load_coords(c_unit)
                self.load_velocities(v_unit)
        self.orient_box(**kwargs)

        xyz = self.coordinates
        r,theta,phi = coordinates.cartesian_to_spherical(xyz[:,0],
                                                         xyz[:,1],
                                                         xyz[:,2])
        self.spherical_coords = numpy.column_stack((r,theta,phi))

        vxyz = self.velocities
        vr,vtheta,vphi = coordinates.cartesian_to_spherical_velocities(xyz,vxyz)
        self.spherical_velocities = numpy.column_stack((vr,vtheta,vphi))

    def calculate_cylindrical_coords(self, c_unit=None, v_unit=None, **kwargs):
        """
        Load particle positions in cylindrical coordinates with default units of
        kpc (default set in units class)
        unit: unit conversion from code units
        """
        if (c_unit or v_unit):
            if (c_unit != self.units._coord_unit
                or v_unit != self.units.velocity_unit):
                self.load_coords(c_unit)
                self.load_velocities(v_unit)
        self.orient_box(**kwargs)
        xyz = self.coordinates
        r,theta,z = coordinates.cartesian_to_cylindrical(xyz[:,0],
                                                         xyz[:,1],
                                                         xyz[:,2])
        self.cylindrical_coords = numpy.column_stack((r,theta,z))

        vel = self.velocities

    def get_coords(self, unit=None, **kwargs):
        """
        Return Particle Coordinates in units of kpc (default set in units class)
        unit: unit conversion from code units
        system: coordinate system to use.  Options are cartesian or spherical.
                Note that spherical coordinates require a cartesian coordinate
                'center' point to be calculated by analyze.find_center()
        """
        system = kwargs.pop('system','cartesian')
        if (unit or len(kwargs) > 0):
            if (unit != self.units._coord_unit or len(kwargs) > 0):
                if system == 'cartesian':
                    self.load_coords(unit)
                    self.orient_box(**kwargs)
                elif system == 'spherical':
                    self.calculate_spherical_coords(c_unit=unit, **kwargs)
                elif system == 'cylindrical':
                    self.calculate_cylindrical_coords(c_unit=unit, **kwargs)

        if system == 'cartesian':
            try:
                return self.coordinates
            except AttributeError:
                self.load_coords(unit)
                self.orient_box(**kwargs)
                return self.coordinates
        elif system == 'spherical':
            try:
                return self.spherical_coords
            except AttributeError:
                self.calculate_spherical_coords(c_unit=unit, **kwargs)
                return self.spherical_coords
        elif system == 'cylindrical':
            try:
                return self.cylindrical_coords
            except AttributeError:
                self.calculate_cylindrical_coords(c_unit=unit, **kwargs)
                return self.cylindrical_coords
        else:
            raise KeyError("Coordinate system options: 'cartesian' "\
                           "'spherical' 'cylindrical'")

    def get_velocities(self, unit=None, **kwargs):
        """
        Return Particle Velocities in units of km/s (default set in units class)
        unit: unit conversion from code units
        system: coordinate system to use.
                Options are cartesian, spherical or cylindrical.
                Note that calculating a spherical/cylindrical coordinates
                require a center point.
        """
        system = kwargs.pop('system','cartesian')
        if (unit or len(kwargs) > 0):
            if (unit != self.units.velocity_unit or len(kwargs) > 0):
                if system == 'cartesian':
                    self.load_velocities(unit)
                    self.orient_box(**kwargs)
                elif system == 'spherical':
                    self.calculate_spherical_coords(v_unit=unit, **kwargs)
                elif system == 'cylindrical':
                    self.calculate_cylindrical_coords(v_unit=unit, **kwargs)

        if system == 'cartesian':
            try:
                return self.velocities
            except AttributeError:
                self.load_velocities(unit)
                self.orient_box(**kwargs)
                return self.velocities
        elif system == 'spherical':
            try:
                return self.spherical_velocities
            except AttributeError:
                self.calculate_spherical_coords(v_unit=unit, **kwargs)
                return self.spherical_velocities
        elif system == 'cylindrical':
            try:
                return self.cylindrical_velocities
            except AttributeError:
                self.calculate_cylindrical_coords(v_unit=unit, **kwargs)
                return self.cylindrical_velocities
        else:
            raise KeyError("Coordinate system options: 'cartesian' "\
                           "'spherical' 'cylindrical'")

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
        to_drop = [key for key in self.columns if key not in exclude]
        print to_drop
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
