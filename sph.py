# sph.py
# Jacob Hummel
"""
This module contains classes for reading Gadget2 SPH particle data.
"""

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
