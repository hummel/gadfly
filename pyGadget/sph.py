# sph.py
# Jacob Hummel
"""
This module contains classes for reading Gadget2 SPH particle data.
"""
import numpy
import units
import constants
import hdf5
import sink

class PartTypeSPH(hdf5.PartTypeX):
    """
    Class for SPH particles.
    Extends class PartTypeX to include gas physics stuff.
    """
    def __init__(self, file_id, units, **kwargs):
        super(PartTypeSPH,self).__init__(file_id,0, units)
        self.__init_load_dict__()
        self._drop_ids = None
        self.refine = kwargs.pop('refine_gas', True)
        if self.refine:
            print 'Turning on gas particle refinement.'
            self.choose_subset()

        self.sink_tracking = kwargs.pop('track_sinks', False)
        if self.sink_tracking:
            print 'Tracking sinks.'
            self.locate_sink_particles()

    def __getstate__(self):
        result = self.__dict__.copy()
        del result['_Coordinates']
        del result['_ParticleIDs']
        del result['_Velocities']
        del result['_Masses']
        del result['_Adiabatic_index']
        del result['_ChemicalAbundances']
        del result['_Density']
        del result['_SinkValue']
        del result['_SmoothingLength']
        del result['_InternalEnergy']
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
                       'internalenergy':self.get_internal_energy,
                       'adiabatic_index':self.get_adiabatic_index,
                       'abundances':self.get_abundances,
                       'sink_value':self.get_sinks,
                       'smoothing_length':self.get_smoothing_length,
                       'electron_frac':self.get_electron_fraction,
                       'h2frac':self.get_H2_fraction,
                       'HDfrac':self.get_HD_fraction}
        sph_derived = {'temp':self.get_temperature,
                       'c_s':self.get_sound_speed,
                       't_ff':self.get_freefall_time,
                       'jeans_length':self.get_jeans_length}
        self._load_dict.update(sph_loaders)
        self._load_dict.update(sph_derived)
        self.loadable_keys = self._load_dict.keys()
        self._calculated.append(sph_derived.keys())

    def choose_subset(self, keys=None, criterion=None):
        if keys is None:
            keys = ['masses', 'sink_value']
        self.load_data(*keys)
        if criterion is None:
            criterion = (self.masses > self.masses.min()) & (self.sink_value == 0.)
        self._drop_ids = self[criterion].index
        self.drop(self._drop_ids, inplace=True)
        print self.index.size, 'particles selected.'

    def locate_refined_particles(self):
        self.load_masses()
        self.load_sinks()
        m_min = self.masses.min()
        drop_ids = self[(self.masses > m_min) & (self.sink_value == 0.)].index
        self._drop_ids = drop_ids
        self.drop(drop_ids, inplace=True)
        print 'There are', self.index.size, 'highest resolution particles.'

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
        try:
            del self.ndensity
        except AttributeError:
            pass
        if unit:
            self.units.set_density(unit)
        self.density = self._Density.value * self.units.density_conv
        if self.units.remove_h:
            h = self._header.HubbleParam
            self.density = self.density * h*h
        if self.units.coordinate_system == 'physical':
            ainv = self._header.Redshift + 1 # 1/(scale factor)
            self.density = self.density * ainv**3
        if self._drop_ids is not None:
            self.density = self.density[self._drop_ids]

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
        try:
            del self.density
        except AttributeError:
            pass
        if unit:
            self.units.set_density(unit)
        self['ndensity'] = self._Density.value * self.units.density_conv
        self.ndensity = self.ndensity * constants.X_h / constants.m_H
        if self.units.remove_h:
            h = self._header.HubbleParam
            self.ndensity = self.ndensity * h*h
        if self.units.coordinate_system == 'physical':
            ainv = self._header.Redshift + 1 # 1/(scale factor)
            self.ndensity = self.ndensity * ainv**3
        if self._drop_ids is not None:
            self.ndensity = self.ndensity[self._drop_ids]

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
        self['internalenergy'] = self._InternalEnergy.value * self.units.energy_conv
        if self._drop_ids is not None:
            self.internalenergy = self.internalenergy[self._drop_ids]


    def get_internal_energy(self, unit=None):
        """
        Return Particle Densities in cgs units (default set in units class)
        unit: unit conversion from code units
        """
        if unit:
            if unit != self.units.energy_unit:
                self.load_internal_energy(unit)
        try:
            return self.internalenergy
        except AttributeError:
            self.load_internal_energy(unit)
            return self.internalenergy

    def load_adiabatic_index(self):
        """
        Load particle adiabatic index.
        """
        self['adiabatic_index'] = self._Adiabatic_index.value
        if self._drop_ids is not None:
            self.adiabatic_index = self.adiabatic_index[self._drop_ids]


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
        self['smoothing_length'] \
            = self._SmoothingLength.value * self.units.length_conv
        if self.units.remove_h:
            h = self._header.HubbleParam
            self.smoothing_length /= h
        if self.units.coordinate_system == 'physical':
            a = self._header.ScaleFactor
            self.smoothing_length *= a
        if self._drop_ids is not None:
            self.smoothing_length = self.smoothing_length[self._drop_ids]

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

    def load_abundances(self):
        """
        Load chemical abundances array.

        There are six abundances tracked for each particle.
        0:H2I 1:HII 2:DII 3:HDI 4:HeII 5:HeIII
        """
        self.abundances = self._ChemicalAbundances.value
        if self._drop_ids is not None:
            self.abundances = self.abundances[self._drop_ids]


    def get_abundances(self, species=None):
        """
        Return chemical abundances array.

        There are six abundances tracked for each particle.
        0:H2I 1:HII 2:DII 3:HDI 4:HeII 5:HeIII
        """
        abundance_dict = {'H2':0, 'HII':1, 'DII':2, 'HD':3, 'HeII':4, 'HeIII':5}
        try:
            abundances = self.abundances
        except AttributeError:
            self.load_abundances()
            abundances = self.abundances
        if species is None:
            return abundances
        elif isinstance(species, basestring):
            return abundances[:,abundance_dict[species]]
        elif isinstance(species, list):
            abunds = []
            for sp in species:
                abunds.append(abundances[:,abundance_dict[sp]])
            return abunds

    def load_electron_fraction(self):
        """
        Load the free electron fraction.
        """
        # Chemical Abundances--> 0:H2I 1:HII 2:DII 3:HDI 4:HeII 5:HeIII
        abundances = self.get_abundances()
        self.electron_frac = abundances[:,1] + abundances[:,2]
        self.electron_frac += abundances[:,4] + 2*abundances[:,5]

    def get_electron_fraction(self):
        """
        Return the free electron fraction.
        """
        try:
            return self.electron_fraction
        except AttributeError:
            self.load_electron_fraction()
            return self.electron_fraction

    def load_sinks(self):
        """
        Load particle by particle sink flag values.
        """
        sinks = self._SinkValue.value
        if self._drop_ids is not None:
            sinks= Series(sinks, index=self._ParticleIDs.value)
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

    def load_H2_fraction(self):
        """
        Load the molecular hydrogen fraction.
        """
        # Chemical Abundances--> 0:H2I 1:HII 2:DII 3:HDI 4:HeII 5:HeIII
        abundances = self.get_abundances()
        self['h2frac'] = 2*abundances[:,0]

    def get_H2_fraction(self):
        """
        Return the molecular hydrogen fraction.
        """
        try:
            return self.h2frac
        except AttributeError:
            self.load_H2_fraction()
            return self.h2frac

    def load_HD_fraction(self):
        """
        Load the molecular HD fraction.
        """
        # Chemical Abundances--> 0:H2I 1:HII 2:DII 3:HDI 4:HeII 5:HeIII
        abundances = self.get_abundances()
        self['HDfrac'] = 2*abundances[:,3]

    def get_HD_fraction(self):
        """
        Return the molecular HD fraction.
        """
        try:
            return self.HDfrac
        except AttributeError:
            self.load_HD_fraction()
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
