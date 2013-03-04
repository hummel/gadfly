# sph.py
# Jacob Hummel
"""
This module contains classes for reading Gadget2 SPH particle data.
"""
import numpy
import units
import constants
import hdf5

class PartTypeSPH(hdf5.PartTypeX):
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

    def load_number_density(self, conv=units.Density_cgs, no_h=True, 
                            comoving=False):
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
            self.ndensity = self.ndensity * h*h
        if not comoving:
            ainv = self._header.Redshift + 1 # 1/(scale factor)
            self.ndensity = self.ndensity * ainv**3

    def get_number_density(self, conv=units.Density_cgs, no_h=True, 
                           comoving=False):
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


    def get_internal_energy(self, conv=units.Energy_cgs/units.Mass_g):
        """
        Return Particle Densities in cgs units (default)
        conv: unit conversion from code units
        """
        try:
            return self.energy
        except AttributeError:
            self.load_internal_energy(conv)
            return self.energy

    def load_gamma(self):
        """
        Load particle adiabatic index.
        """
        self.gamma = self._Adiabatic_index.value


    def get_gamma(self):
        """
        Return particle adiabatic index.
        """
        try:
            return self.gamma
        except AttributeError:
            self.load_gamma()
            return self.gamma

    def load_abundances(self):
        """
        Load chemical abundances array.

        There are six abundances tracked for each particle.
        0:H2I 1:HII 2:DII 3:HDI 4:HeII 5:HeIII
        """
        self.abundances = self._ChemicalAbundances.value


    def get_abundances(self):
        """
        Return chemical abundances array.

        There are six abundances tracked for each particle.
        0:H2I 1:HII 2:DII 3:HDI 4:HeII 5:HeIII
        """
        try:
            return self.abundances
        except AttributeError:
            self.load_abundances()
            return self.abundances

    def load_sinks(self):
        """
        Load particle sink values.
        """
        self.sink_value = self._SinkValue.value


    def get_sinks(self):
        """
        Return particle sink values.
        """
        try:
            return self.sink_value
        except AttributeError:
            self.load_sinks()
            return self.sink_value

    def load_smoothing_length(self, conv=units.Length_kpc, 
                              no_h=True, comoving=False):
        """
        Load Particle Smoothing Lengths in units of kpc (default)
        conv: unit conversion from code units
        no_h: if True, remove dependence on hubble parameter.
        comoving: if False, remove dependence on scale factor.
        """
        self.smoothing_length = self._SmoothingLength.value*conv
        if no_h:
            h = self._header.HubbleParam
            self.smoothing_length = self.smoothing_length/h
        if not comoving:
            a = self._header.ScaleFactor
            self.smoothing_length = self.smoothing_length*a

    def get_smoothing_length(self, conv=units.Length_kpc, 
                             no_h=True, comoving=False):
        """
        Return Particle Smoothing Lengths in units of kpc (default)
        conv: unit conversion from code units
        no_h: if True, remove dependence on hubble parameter.
        comoving: if False, remove dependence on scale factor.
        """
        try:
            return self.smoothing_length
        except AttributeError:
            self.load_smoothing_length(conv, no_h, comoving)
            return self.smoothing_length

    def load_electron_fraction(self):
        """
        Load the free electron fraction.
        """
        # Chemical Abundances--> 0:H2I 1:HII 2:DII 3:HDI 4:HeII 5:HeIII
        abundances = self.get_abundances()
        self.electron_frac = abundances[:,1] + abundances[:,4] + abundances[:,5]

    def get_electron_fraction(self):
        """
        Return the free electron fraction.
        """
        try:
            return self.electron_frac
        except AttributeError:
            self.load_electron_fraction()
            return self.electron_frac

    def load_H2_fraction(self):
        """
        Load the molecular hydrogen fraction.
        """
        # Chemical Abundances--> 0:H2I 1:HII 2:DII 3:HDI 4:HeII 5:HeIII
        abundances = self.get_abundances()
        self.h2frac = 2*abundances[:,0]

    def get_H2_fraction(self):
        """
        Return the molecular hydrogen fraction.
        """
        try:
            return self.h2frac
        except AttributeError:
            self.load_H2_fraction()
            return self.h2frac

    def calculate_temperature(self):
        """
        Calculate Particle Temperatures in degrees Kelvin.
        """
        gamma = self.get_gamma()
        energy = self.get_internal_energy()
        h2frac = self.get_H2_fraction()
        mu = (0.24/4.0) + ((1.0-h2frac)*0.76) + (h2frac*.76/2.0)
        mu = 1 / mu # mean molecular weight
        self.temp = (mu * constants.m_H / constants.k_B) * energy * (gamma-1)

    def get_temperature(self):
        """
        Return Particle Temperatures in degrees Kelvin.
        """
        try:
            return self.temp
        except AttributeError:
            self.calculate_temperature()
            return self.temp

    def calculate_sound_speed(self):
        """
        Calculate the sound speed of the gas in cm/s.
        """
        temp = self.get_temperature()
        self.c_s = numpy.sqrt(constants.k_B*temp/constants.m_H)

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
        denominator = 32 * numpy.pi * constants.G * rho
        self.t_ff = numpy.sqrt(3/denominator)

    def get_freefall_time(self):
        """
        Return Jeans Length for gas in cm.
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

    def calculate_optical_depth(self,sigma):
        """
        Estimate optical depth of the gas based on Jeans Length.
        sigma:: cross-section of the species of interest.
        """
        n = self.get_number_density()
        L = self.get_jeans_length()
        self.tau = n*L*sigma

    def get_optical_depth(self,sigma):
        """
        Return estimate of the optical depth based on Jeans Length.
        sigma:: cross-section of the species of interest.
        """
        try:
            return self.tau
        except AttributeError:
            self.calculate_optical_depth(self,sigma)
            return self.tau
