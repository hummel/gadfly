# nbody.py
# Jacob Hummel
"""
This module contains classes for reading Gadget2 N-body (dark matter)
particle data.
"""
import numpy
from pandas import Series, DataFrame

from hdf5 import HDF5group
import units
import constants
import coordinates

class PartTypeNbody(HDF5group):
    """
    Class for N-body particles. 
    Extends: hdf5.HDF5group
    """
    def __init__(self, file_id, ptype, units, **kwargs):
        super(PartTypeNbody,self).__init__(file_id, ptype, units, **kwargs)
        self.__init_load_dict__()
        self._drop_ids = None
        self.refined = kwargs.pop('refine', False)
        if self.refined:
            print 'Turning on N-body particle refinement.'
            self.refine_dataset()

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

    def refine_dataset(self, *keys, **kwargs):
        if len(keys) < 1:
            keys = ['masses']
            self.load_data(*keys)
        criterion = kwargs.pop('criterion', None)
        if criterion is None:
	        criterion = (self.masses > self.masses.min())
	        super(PartTypeNbody, self).refine_dataset(criterion)

    def load_masses(self, unit=None):
        """
        Load Particle Masses in units of M_sun (default set in units class)
        unit: unit conversion from code units
        """
        if unit:
            self.units.set_mass(unit)
        masses = self._Masses.value * self.units.mass_conv
        if self.units.remove_h:
            h = self._header.HubbleParam
            masses /= h
        if self._drop_ids is not None:
            masses = Series(masses, index=self._ParticleIDs.value)
            self['masses'] = masses.drop(self._drop_ids)
        else:
            self['masses'] = masses

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
        xyz = DataFrame(xyz, index=self._ParticleIDs.value, columns=['x', 'y', 'z'])
        if self._drop_ids is not None:
            xyz.drop(self._drop_ids, inplace=True)
        self[['x', 'y', 'z']] = xyz

    def load_velocities(self, unit=None):
        """
        Load Particle Velocities in units of km/s (default set in units class)
        unit: unit conversion from code units
        """
        if unit:
            self.units.set_velocity(unit)
        uvw = self._Velocities.value * self.units.velocity_conv
        if self.units.coordinate_system == 'physical':
            a = self._header.ScaleFactor
            uvw *= numpy.sqrt(a)
        uvw = DataFrame(uvw, index=self._ParticleIDs.value, columns=['u', 'v', 'w'])
        if self._drop_ids is not None:
            uvw.drop(self._drop_ids, inplace=True)
        self[['u', 'v', 'w']] = uvw

    def orient_box(self, **kwargs):
        """
        Center and rotate box coordinates AND velocities according to 
        received kwargs 'center' and 'view'. If 'center' unspecified,
        looks for 'centering' kwargs and attempts to auto-center the box.
        """
        center = kwargs.get('center', None)
        vcenter = kwargs.get('vcenter', None)
        centering = kwargs.get('centering', None)
        view = kwargs.get('view', None)
        dlim = kwargs.pop('dens_lim', 1e9)
        xyz = ['x', 'y', 'z']
        uvw = ['u', 'v', 'w']
        try:
            pos_vel = self[xyz + uvw]
        except KeyError:
            self.load_coords()
            self.load_velocities()
            pos_vel = self[xyz + uvw]

        if center:
            if vcenter is None:
                print "Warning: Not Re-centering particle Velocities!"
            pos_vel = analyze.center_box(pos_vel, center, vcenter)
        else:
            if centering in ['avg','max']:
                try:
                    dens = self.get_number_density()
                except AttributeError:
                    raise KeyError("Cannot density-center dark matter!")
                pos_vel = analyze.center_box(pos_vel, density=dens, **kwargs)
            elif centering == 'box':
                pos_vel = analyze.center_box(pos_vel, **kwargs)
        self[pos_vel.keys()] = pos_vel

        if view:
            print 'Rotating Box...'
            if view == 'face':
                try:
                    dens = self.get_number_density()
                except AttributeError:
                    raise KeyError("Cannot density-center dark matter!")
                xyz, uvw = visualize.set_view(view, self[xyz], velocity=self[uvw],
                                               density=dens, dens_lim=dlim)
            else:
                xyz, uvw = visualize.set_view(view, self[xyz], velocity=self[uvw])
            print 'Rotation complete.'
            self[['x', 'y', 'z']] = xyz
            self[['u', 'v', 'w']] = uvw

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
        print 'Converting to spherical coordinates...'
        coordinates.cartesian_to_spherical(self)

        #uvw = self[['u', 'v', 'w']]
        #print 'Converting to spherical coordinate velocities...'
        #vr,vtheta,vphi = coordinates.cartesian_to_spherical_velocities(xyz,uvw)
        #self.spherical_velocities = numpy.column_stack((vr,vtheta,vphi))
        #print 'Done.'

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
        print 'Converting to cylindrical coordinates...'
        coordinates.cartesian_to_cylindrical(self)

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
                return self[['x', 'y', 'z']]
            except KeyError:
                self.load_coords(unit)
                self.orient_box(**kwargs)
                return self[['x', 'y', 'z']]
        elif system == 'spherical':
            try:
                return self.spherical_coords
            except KeyError:
                self.calculate_spherical_coords(c_unit=unit, **kwargs)
                return self.spherical_coords
        elif system == 'cylindrical':
            try:
                return self.cylindrical_coords
            except KeyError:
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
                return self[['u', 'v', 'w']]
            except KeyError:
                self.load_velocities(unit)
                self.orient_box(**kwargs)
                return self[['u', 'v', 'w']]
        elif system == 'spherical':
            try:
                return self.spherical_velocities
            except KeyError:
                self.calculate_spherical_coords(v_unit=unit, **kwargs)
                return self.spherical_velocities
        elif system == 'cylindrical':
            try:
                return self.cylindrical_velocities
            except KeyError:
                self.calculate_cylindrical_coords(v_unit=unit, **kwargs)
                return self.cylindrical_velocities
        else:
            raise KeyError("Coordinate system options: 'cartesian' "\
                           "'spherical' 'cylindrical'")
