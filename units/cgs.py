# cgs.py
# Jacob Hummel
"""
Physical cgs unit conversions for analyzing my Gadget2 HDF5 snapshot data.
"""
Mass_g = 1.989e43
Length_cm = 3.085678e21
Velocity_cgs= 1.0e5

Time_s= Length_cm / Velocity_cgs
Density_cgs= Mass_g / Length_cm**3
Pressure_cgs= Mass_g / Length_cm/ Time_s**2
Energy_cgs= Mass_g * Length_cm**2 / Time_s**2

Time_yr = Time_s/3.15569e7
Time_myr = Time_yr/1e6
Time_gyr = Time_yr/1e9

