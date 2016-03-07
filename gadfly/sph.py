# sph.py
# Jacob Hummel
"""
This module contains classes for reading Gadget2 SPH particle data.
"""
import numpy
from pandas import Series, DataFrame

from nbody import PartTypeNbody
import units

class PartTypeSPH(PartTypeNbody):
    """
    Class for SPH particles.
    Extends: nbody.PartTypeNbody
    """
    def __init__(self, file_id, ptype, sim, **kwargs):
        kwargs.pop('refine_nbody', None)
        super(PartTypeNbody,self).__init__(file_id, ptype, sim, **kwargs)
        self.__init_load_dict__()

        self.refined = kwargs.pop('refine', False)
        if self.refined:
            print 'Turning on gas particle refinement.'
            self.refine_dataset()

    def __getstate__(self):
        result = self.__dict__.copy()
        del result['_coordinates']
        del result['_particleIDs']
        del result['_velocities']
        del result['_masses']
        del result['_density']
        del result['_smoothing_length']
        del result['_internal_energy']
        del result['_load_dict']
        del result['loadable_keys']
        return result

    def __setstate__(self, in_dict):
        self.__dict__ = in_dict
        self.__init_load_dict__()

    def __init_load_dict__(self):
        super(PartTypeSPH,self).__init_load_dict__()
        sph_loaders = {'density':self.get_density,
                       'internal_energy':self.get_internal_energy,
                       'smoothing_length':self.get_smoothing_length
                       }
        self._load_dict.update(sph_loaders)
        self.loadable_keys = self._load_dict.keys()

    def refine_dataset(self, *keys, **kwargs):
        if len(keys) < 1:
            keys = ['masses', 'sink_value']
            self.load_data(*keys)
        criterion = kwargs.pop('criterion', None)
        if criterion is None:
            criterion = (self.masses > self.masses.min()) & (self.sink_value == 0.)
        super(PartTypeNbody, self).refine_dataset(criterion)

    def load_density(self, unit=None):
        """
        Load Particle Densities in cgs units (default set in units class)
        unit: unit conversion from code units
        """
        if unit:
            self.units.set_density(unit)
        density = self._density.value * self.units.density_conv
        if self.units.remove_h:
            h = self._header.HubbleParam
            density *=  h**2
        if self.units.coordinate_system == 'physical':
            ainv = self._header.Redshift + 1 # 1/(scale factor)
            density *= ainv**3
        if self._drop_ids is not None:
            density = Series(density, index=self._particleIDs.value)
            self['density'] = density.drop(self._drop_ids)
        else:
            self['density'] = density

    def get_density(self, unit=None):
        """
        Return Particle Densities in cgs units (default set in units class)
        unit: unit conversion from code units
        """
        if unit:
            if unit != self.units.density_unit:
                self.load_density(unit)
        try:
            return self.density
        except AttributeError:
            self.load_density(unit)
            return self.density

    def load_internal_energy(self, unit=None):
        """
        Load internal particle energies per unit mass in cgs units 
        (default set in units class)
        unit: unit conversion from code units
        """
        if unit:
            self.units.set_energy(unit)
        energy = self._internal_energy.value * self.units.energy_conv
        if self._drop_ids is not None:
            energy = Series(energy, index=self._particleIDs.value)
            self['internal_energy'] = energy.drop(self._drop_ids)
        else:
            self['internal_energy'] = energy

    def get_internal_energy(self, unit=None):
        """
        Return Particle Densities in cgs units (default set in units class)
        unit: unit conversion from code units
        """
        if unit:
            if unit != self.units.energy_unit:
                self.load_internal_energy(unit)
        try:
            return self.internal_energy
        except AttributeError:
            self.load_internal_energy(unit)
            return self.internal_energy

    def load_smoothing_length(self, unit=None):
        """
        Load Particle Smoothing Lengths in units of kpc.
        (default set in units class)
        unit: unit conversion from code units
        """
        if unit:
            self.units._set_smoothing_length(unit)
        hsml = self._smoothing_length.value * self.units.length_conv
        if self.units.remove_h:
            h = self._header.HubbleParam
            hsml /= h
        if self.units.coordinate_system == 'physical':
            a = self._header.ScaleFactor
            hsml *= a
        if self._drop_ids is not None:
            hsml = Series(hsml, index=self._particleIDs.value)
            self['smoothing_length'] = hsml.drop(self._drop_ids)
        else:
            self['smoothing_length'] = hsml

    def get_smoothing_length(self, unit=None):
        """
        Return Particle Smoothing Lengths in units of kpc.
        (default set in units class)
        unit: unit conversion from code units
        """
        if unit:
            if unit != self.units._smoothing_unit:
                self.load_smoothing_length(unit)
        try:
            return self.smoothing_length
        except AttributeError:
            self.load_smoothing_length(unit)
            return self.smoothing_length