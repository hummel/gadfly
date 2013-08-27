# units.py
# Jacob Hummel
"""
Physical cgs unit conversions for analyzing my Gadget2 HDF5 snapshot data.
"""
### Code units:
Mass_g = 1.989e43
Length_cm = 3.085678e21
Velocity_cgs= 1.0e5

### cgs Conversions
Time_s= Length_cm / Velocity_cgs
Density_cgs= Mass_g / Length_cm**3
Pressure_cgs= Mass_g / Length_cm/ Time_s**2
Energy_cgs= Mass_g * Length_cm**2 / Time_s**2

### Additional Mass Conversions
Mass_sun = Mass_g / 1.989e33

### Additional Distance Conversions
Length_AU = Length_cm / 1.49598e13
Length_pc = Length_cm / 3.085678e18
Length_kpc = Length_pc / 1e3

### Additional Velocity Conversions
Velocity_kms = Velocity_cgs / 1e5

### Additional Time Conversions
Time_yr = Time_s / 3.15569e7
Time_myr = Time_yr / 1e6
Time_gyr = Time_yr / 1e9

### Unit Selection Dictionaries

class Units(object):
    def __init__(self, **unitargs):
        super(Units,self).__init__()
        self.remove_h = not unitargs.pop('units_over_h', False)
        self.set_coordinate_system(unitargs.pop('coordinates', 'physical'))
        self.set_length(unitargs.pop('length', 'kpc'))
        self.set_mass(unitargs.pop('mass', 'solar'))
        self.set_time(unitargs.pop('time', 'yr'))
        self.set_velocity(unitargs.pop('velocity', 'kms'))
        self.set_density(unitargs.pop('density', 'cgs'))
        self.set_pressure(unitargs.pop('pressure', 'cgs'))
        self.set_energy(unitargs.pop('energy', 'cgs'))
        self._coord_unit = self.length_unit
        self._smoothing_unit = self.length_unit

    def set_coordinate_system(self,coordinates):
        if coordinates not in ['physical', 'comoving']:
            raise KeyError
        self.coordinate_system = coordinates

    def set_length(self, unit):
        Lengths = {'cm':Length_cm, 'AU':Length_AU,
                   'pc':Length_pc, 'kpc':Length_kpc}
        self.length_unit = unit
        self.length_conv = Lengths[self.length_unit]

    def _set_coord_length(self, unit):
        self.set_length(unit)
        self._coord_unit = unit

    def _set_smoothing_length(self, unit):
        self.set_length(unit)
        self._smoothing_unit = unit

    def set_mass(self, unit):
        Masses = {'g':Mass_g, 'solar':Mass_sun}
        self.mass_unit = unit
        self.mass_conv = Masses[self.mass_unit]

    def set_time(self, unit):
        Times = {'s':Time_s, 'yr':Time_yr, 'myr':Time_myr, 'gyr':Time_gyr}
        self.time_unit = unit
        self.time_conv = Times[self.time_unit]

    def set_velocity(self, unit):
        Velocities = {'cgs':Velocity_cgs, 'kms':Velocity_kms}
        self.velocity_unit = unit
        self.velocity_conv = Velocities[self.velocity_unit]

    def set_density(self, unit):
        Densities = {'cgs':Density_cgs}
        self.density_unit = unit
        self.density_conv = Densities[self.density_unit]

    def set_pressure(self, unit):
        Pressures = {'cgs':Pressure_cgs}
        self.pressure_unit = unit
        self.pressure_conv = Pressures[self.pressure_unit]

    def set_energy(self, unit):
        Energies = {'cgs':Energy_cgs, 'specific cgs':Energy_cgs/Mass_g}
        self.energy_unit = unit
        self.energy_conv = Energies[self.energy_unit]
