# analyze.py
# Jacob Hummel
import numpy
from . import units
from . import constants

#===============================================================================
def halo_properties(snapshot, 
                    r_start=3.08568e18, r_multiplier=1.01, verbose=True):
    h = snapshot.header.HubbleParam
    a = snapshot.header.ScaleFactor
    redshift = snapshot.header.Redshift

    length_unit = units.Length_cm
    dm_mass = snapshot.dm.get_masses()
    dm_pos = snapshot.dm.get_coords(length_unit)
    gas_mass = snapshot.gas.get_masses()
    gas_pos = snapshot.gas.get_coords(length_unit)
    dens = snapshot.gas.get_number_density()

    gasx = gas_pos[:,0]
    gasy = gas_pos[:,1]
    gasz = gas_pos[:,2]
    dmx = dm_pos[:,0]
    dmy = dm_pos[:,1]
    dmz = dm_pos[:,2]

    mass = numpy.concatenate((gas_mass, dm_mass))
    x = numpy.concatenate((gasx,dmx))
    y = numpy.concatenate((gasy,dmy))
    z = numpy.concatenate((gasz,dmz))

    center = gas_pos[dens.argmax()]
    if verbose: print 'Center:', center
    x = x - center[0]
    y = y - center[1]
    z = z - center[2]
    r = numpy.sqrt(numpy.square(x) + numpy.square(y) + numpy.square(z))
    
    halo_properties = []
    n = old_n = density = 0
    background_density = .27 * 9.31e-30 * (1+redshift)**3 #Omega_m * rho_crit(z)
    rmax = r_start
    while (density > 180 * background_density or n < 100):
        inR = numpy.where(r <= rmax)[0]
        n = inR.size
        if n > old_n:
            properties = []
            rpc = rmax/3.08568e18
            if verbose: print 'R = %.2e pc' %rpc,
            Mtot = mass[inR].sum()
            solar_masses = Mtot/1.989e33
            if verbose: print 'Mass enclosed: %.2e' %solar_masses,
            density = 3 * Mtot / (4*numpy.pi * rmax**3)
            delta = density/background_density
            if verbose: print 'delta: %.3f' %delta,
            if verbose: print 'n:', n
            energy = constants.G * Mtot**2 / rmax

            properties.append(rpc)
            properties.append(delta)
            properties.append(solar_masses)
            properties.append(density)
            properties.append(energy)
            properties.append(n)
            halo_properties.append(properties)
            old_n = n
        rmax *= r_multiplier
    
    return numpy.asarray(halo_properties)

