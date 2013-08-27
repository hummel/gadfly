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
Lengths = {'cm':Length_cm, 'AU':Length_AU, 'pc':Length_pc, 'kpc':Length_kpc}
Masses = {'g':Mass_g, 'solar':Mass_sun}
Times = {'s'Time_s, 'yr':Time_yr, 'myr':Time_myr, 'gyr':Time_gyr}
Velocities = {'cgs':Velocit_cgs, 'kms':Velocity_kms}
Densities = {'cgs':Density_cgs}
Pressures = {'cgs':Pressure_cgs}
Energies = {'cgs':Energy_cgs}

class Unit(object):
    def __init__(self, **unitargs):
        super(Units,self).__init__()
        self.length_unit = unitargs.pop('length', 'kpc')
        self.length_conv = Lengths[self.length_unit]
        self.mass_unit = unitargs.pop('mass', 'kpc')
        self.mass_conv = Masses[self.mass_unit]
        self.time_unit = unitargs.pop('time', 'kpc')
        self.time_conv = Times[self.time_unit]
        self.velocity_unit = unitargs.pop('velocity', 'kms')
        self.velocity_conv = Velocities[self.velocity_unit]
        self.density_unit = unitargs.pop('density', 'kpc')
        self.density_conv = Densities[self.density_unit]
        self.pressure_unit = unitargs.pop('pressure', 'kpc')
        self.pressure_conv = Pressures[self.pressure_unit]
        self.energy_unit = unitargs.pop('energy', 'kpc')
        self.energy_conv = Energies[self.energy_unit]
