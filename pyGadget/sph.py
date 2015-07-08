# sph.py
# Jacob Hummel
"""
This module contains classes for reading Gadget2 SPH particle data.
"""
import numpy
from pandas import Series, DataFrame

from nbody import PartTypeNbody
import units
import constants
import sink

class PartTypeSPH(PartTypeNbody):
    """
    Class for SPH particles.
    Extends: nbody.PartTypeNbody
    """
    def __init__(self, file_id, sim, **kwargs):
        kwargs.pop('refine_nbody', None)
        super(PartTypeNbody,self).__init__(file_id,0, sim, **kwargs)
        self.__init_load_dict__()

        self._drop_ids = None
        self.refined = kwargs.pop('refine', False)
        if self.refined:
            print 'Turning on gas particle refinement.'
            self.refine_dataset()

        self.sink_tracking = kwargs.pop('track_sinks', False)
        if self.sink_tracking:
            print 'Tracking sinks.'
            self.locate_sink_particles()

    def __getstate__(self):
        result = self.__dict__.copy()
        del result['_coordinates']
        del result['_particleIDs']
        del result['_velocities']
        del result['_masses']
        del result['_adiabatic_index']
        del result['_abundances']
        del result['_density']
        del result['_sink_value']
        del result['_smoothing_length']
        del result['_internal_energy']
        del result['_load_dict']
        del result['loadable_keys']
        del result['_calculated']
        return result

    def __setstate__(self, in_dict):
        self.__dict__ = in_dict
        self.__init_load_dict__()

    def __init_load_dict__(self):
        super(PartTypeSPH,self).__init_load_dict__()
        sph_loaders = {'density':self.get_density,
                       'ndensity':self.get_number_density,
                       'internal_energy':self.get_internal_energy,
                       'adiabatic_index':self.get_adiabatic_index,
                       'abundances':self.get_abundances,
                       'sink_value':self.get_sinks,
                       'smoothing_length':self.get_smoothing_length}
        sph_derived = {'h2frac':self.get_H2_fraction,
                       'HDfrac':self.get_HD_fraction,
                       'electron_frac':self.get_electron_fraction,
                       'temp':self.get_temperature,
                       'c_s':self.get_sound_speed,
                       't_ff':self.get_freefall_time,
                       'jeans_length':self.get_jeans_length}
        self._load_dict.update(sph_loaders)
        self._load_dict.update(sph_derived)
        self.loadable_keys = self._load_dict.keys()
        self._calculated.append(sph_derived.keys())

    def refine_dataset(self, *keys, **kwargs):
        if len(keys) < 1:
            keys = ['masses', 'sink_value']
            self.load_data(*keys)
        criterion = kwargs.pop('criterion', None)
        if criterion is None:
            criterion = (self.masses > self.masses.min()) & (self.sink_value == 0.)
        super(PartTypeNbody, self).refine_dataset(criterion)

    def locate_sink_particles(self):
        sinks = self.get_sinks()
        self._sink_indices = numpy.where(sinks != 0.)[0]
        print self._sink_indices.size, 'sinks found.'

    def get_sink_properties(self):
        mass = self.get_masses()
        pos = self.get_coords()
        vel = self.get_velocities()
        pid = self.get_PIDs()
        sinks = []
        for index in self._sink_indices:
            sinks.append(sink.Sink(m = mass[index], pid = pid[index],
                                   pos = pos[index], vel=vel[index],
                                   index=index))
        return sinks

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

    def load_number_density(self, unit=None):
        """
        Load Particle Number Densities in cgs units (default set in units class)
        unit: unit conversion from code units
        """
        if unit:
            self.units.set_density(unit)
        ndensity = self._density.value * self.units.density_conv \
                   * constants.X_h / constants.m_H
        if self.units.remove_h:
            h = self._header.HubbleParam
            ndensity *=  h**2
        if self.units.coordinate_system == 'physical':
            ainv = self._header.Redshift + 1 # 1/(scale factor)
            ndensity *= ainv**3
        if self._drop_ids is not None:
            ndensity = Series(ndensity, index=self._particleIDs.value)
            self['ndensity'] = ndensity.drop(self._drop_ids)
        else:
            self['ndensity'] = ndensity

    def get_number_density(self, unit=None):
        """
        Return Particle Densities in cgs units (default set in units class)
        unit: unit conversion from code units
        """
        if unit:
            if unit != self.units.density_unit:
                self.load_number_density(unit)
        try:
            return self.ndensity
        except AttributeError:
            self.load_number_density(unit)
            return self.ndensity

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

    def load_adiabatic_index(self):
        """
        Load particle adiabatic index.
        """
        gamma = self._adiabatic_index.value
        if self._drop_ids is not None:
            gamma = Series(gamma, index=self._particleIDs.value)
            self['adiabatic_index'] = gamma.drop(self._drop_ids)
        else:
            self['adiabatic_index'] = gamma


    def get_adiabatic_index(self):
        """
        Return particle adiabatic index.
        """
        try:
            return self.adiabatic_index
        except AttributeError:
            self.load_adiabatic_index()
            return self.adiabatic_index

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

    def load_sinks(self):
        """
        Load particle by particle sink flag values.
        """
        sinks = self._sink_value.value
        if self._drop_ids is not None:
            sinks = Series(sinks, index=self._particleIDs.value)
            self['sink_value'] = sinks.drop(self._drop_ids)
        else:
            self['sink_value'] = sinks

    def get_sinks(self):
        """
        Return particle by particle sink flag values.
        """
        try:
            return self.sink_value
        except AttributeError:
            self.load_sinks()
            return self.sink_value

    def load_abundances(self, tracked_species=None):
        """
        Load chemical abundances array.

        There are six abundances tracked for each particle.
        0:H2 1:HII 2:DII 3:HD 4:HeII 5:HeIII
        """
        default_species = ['H2', 'HII', 'DII', 'HD', 'HeII', 'HeIII']
        if tracked_species is None:
            tracked_species = default_species
        abundances = self._abundances.value
        abundances = DataFrame(abundances, index=self._particleIDs.value,
                               columns=tracked_species)
        if self._drop_ids is not None:
            abundances.drop(self._drop_ids)
        self[tracked_species] = abundances

    def get_abundances(self, *species, **kwargs):
        """
        Return requested species from chemical abundances array.

        There are six abundances tracked for each particle.
        0:H2 1:HII 2:DII 3:HD 4:HeII 5:HeIII
        """
        default_species = ['H2', 'HII', 'DII', 'HD', 'HeII', 'HeIII']
        all_species = kwargs.get('tracked_species', default_species)
        if len(species) < 1:
            try:
                return self[all_species]
            except(KeyError):
                self.load_abundances(**kwargs)
                return self[all_species]
        else:
            try:
                return self[list(species)]
            except(KeyError):
                self.load_abundances(**kwargs)
                return self[list(species)]

    def calculate_electron_fraction(self):
        """
        Calculate the free electron fraction.
        """
        # Chemical Abundances--> 0:H2I 1:HII 2:DII 3:HDI 4:HeII 5:HeIII
        abundances = self.get_abundances('HII', 'DII', 'HeII', 'HeIII')
        self['electron_fraction'] = abundances.HII + abundances.DII \
                                    + abundances.HeII + 2*abundances.HeIII

    def get_electron_fraction(self):
        """
        Return the free electron fraction.
        """
        try:
            return self.electron_fraction
        except AttributeError:
            self.calculate_electron_fraction()
            return self.electron_fraction

    def calculate_H2_fraction(self):
        """
        Calculate the molecular hydrogen fraction.
        """
        self['h2frac'] = 2 * self.get_abundances('H2')

    def get_H2_fraction(self):
        """
        Return the molecular hydrogen fraction.
        """
        try:
            return self.h2frac
        except AttributeError:
            self.calculate_H2_fraction()
            return self.h2frac

    def calculate_HD_fraction(self):
        """
        Calculate the molecular HD fraction.
        """
        self['HDfrac'] = 2 * self.get_abundances('HD')

    def get_HD_fraction(self):
        """
        Return the molecular HD fraction.
        """
        try:
            return self.HDfrac
        except AttributeError:
            self.calculate_HD_fraction()
            return self.HDfrac

    def calculate_temperature(self):
        """
        Calculate Particle Temperatures in degrees Kelvin.
        """
        gamma = self.get_adiabatic_index()
        energy = self.get_internal_energy('specific cgs')
        h2frac = self.get_H2_fraction()
        mu = (0.24/4.0) + ((1.0-h2frac)*0.76) + (h2frac*.76/2.0)
        mu = 1 / mu # mean molecular weight
        self['temperature'] = (mu * constants.m_H / constants.k_B) * energy * (gamma-1)

    def get_temperature(self):
        """
        Return Particle Temperatures in degrees Kelvin.
        """
        try:
            return self.temperature
        except AttributeError:
            self.calculate_temperature()
            return self.temperature

    def calculate_sound_speed(self):
        """
        Calculate the sound speed of the gas in cm/s.
        """
        temp = self.get_temperature()
        h2frac = self.get_H2_fraction()
        mu = (0.24/4.0) + ((1.0-h2frac)*0.76) + (h2frac*.76/2.0)
        mu = 1 / mu # mean molecular weight
        self.c_s = numpy.sqrt(constants.k_B*temp/(mu*constants.m_H))

    def get_sound_speed(self):
        """
        Return the sound speed of the gas in cm/s.
        """
        try:
            return self.c_s
        except AttributeError:
            self.calculate_sound_speed()
            return self.c_s

    def calculate_freefall_time(self):
        """
        Calculate the freefall time of the gas in s.
        """
        rho = self.get_density() # NOT number density!
        denominator = 32 * numpy.pi * constants.GRAVITY* rho
        self.t_ff = numpy.sqrt(3/denominator)

    def get_freefall_time(self):
        """
        Return the freefall time of the gas in s.
        """
        try:
            return self.t_ff
        except AttributeError:
            self.calculate_freefall_time()
            return self.t_ff

    def calculate_jeans_length(self):
        """
        Calculate Jeans Length for gas in cm.
        """
        c_s = self.get_sound_speed()
        t_ff = self.get_freefall_time()
        self.jeans_length = c_s * t_ff

    def get_jeans_length(self):
        """
        Return Jeans Length for gas in cm.
        """
        try:
            return self.jeans_length
        except AttributeError:
            self.calculate_jeans_length()
            return self.jeans_length
