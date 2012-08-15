# constants.py
# Jacob Hummel
"""
This module contains the fundamental constants needed for my Gadget2 HDF5 snapshot data.
All values in cgs units.
"""
from units import cgs
### Fundamental Constants
k_B = 1.3806e-16 # erg/K
m_H = 1.6726e-24 # g
GRAVITY = 6.6726e-8 # dyne * cm**2 / g**2
G = GRAVITY / cgs.Length_cm**3 * cgs.Mass_g * cgs.Time_s**2
X_h = 0.76 # Hydrogen Mass Fraction
