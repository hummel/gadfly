# radial_props.py
# Jacob Hummel
# script for reducing simulation data by radial properties of halo.
import sys
import numpy
import pandas as pd
import pyGadget as pyg

#===============================================================================
def radial_properties(sim, snap, **kwargs):
    r_start = kwargs.pop('r_start', 3.08568e14)
    r_multiplier = kwargs.pop('multiplier', 1.2)
    verbose = kwargs.pop('verbose', True)
    n_min = kwargs.pop('n_min', 32)

    length_unit = 'cm'
    mass_unit = 'g'

    snapshot = sim.load_snapshot(snap)
    redshift = snapshot.header.Redshift
    dens = snapshot.gas.get_density('cgs')
    dm_mass = snapshot.dm.get_masses(mass_unit)
    gas_mass = snapshot.gas.get_masses(mass_unit)
    temp = snapshot.gas.get_temperature()

    xyz = snapshot.gas.get_coords(length_unit)
    uvw = snapshot.gas.get_velocities('cgs')
    center, vcenter = pyg.analyze.find_center_vcenter(xyz[:,0],xyz[:,1],xyz[:,2],
                                              uvw[:,0],uvw[:,1],uvw[:,2],
                                              dens, centering='avg')
    gas_pos = snapshot.gas.get_coords(length_unit, system='spherical',
                                      center=center, vcenter=vcenter, view='face')
    vel = snapshot.gas.get_velocities(system='spherical')
    dm_pos = snapshot.dm.get_coords(length_unit, system='spherical',
                                    center=center, vcenter=vcenter)
    snapshot.close()

    gasr = gas_pos[:,0]
    dmr = dm_pos[:,0]
    r = numpy.concatenate((gasr,dmr))
    mass = numpy.concatenate((gas_mass, dm_mass))

    print 'Data loaded.  Analyzing...'
    halo_properties = analyze_halo(redshift, r, gasr, mass, gas_mass, temp, dens, vel,
                                   r_start, r_multiplier, n_min, verbose)
    print 'snapshot', snapshot.number, 'analyzed.'
    colnames = ('redshift', 'radius', 'delta', 'Mshell', 'gMshell', 'energy',
                'density', 'gdensity', 'rhoShell', 'grhoShell', 'dens_avg',
                'tshell', 'vr', 'vtheta', 'vphi', 'vr_sigma', 'vtheta_sigma',
                'vphi_sigma', 'npart')
    halo_properties = pd.DataFrame(halo_properties, columns=colnames)
    halo_properties['snapshot'] = snap
    return halo_properties

#@autojit
def analyze_halo(redshift, r, gasr, mass, gmass, temp, dens, vel,
                 r_start, r_multiplier, n_min, verbose):
    k_B = 1.3806e-16 # erg/K
    m_H = 1.6726e-24 # g
    GRAVITY = 6.6726e-8 # dyne * cm**2 / g**2
    halo_properties = []
    n = 0
    old_n = 0
    old_r = 0
    density = 0
    energy = 0
    # background density:: Omega_m * rho_crit(z)
    background_density = .27 * 9.31e-30 * (1+redshift)**3 
    rmax = r_start
    while (density > 178 * background_density or n < 100):
        inR = numpy.where(r <= rmax)[0]
        gasinR = numpy.where(gasr <= rmax)[0]
        n = inR.size
        if n > old_n + n_min:
            inShell = numpy.where((r > old_r) & (r <= rmax))[0]
            gasinShell = numpy.where((gasr > old_r) & (gasr <= rmax))[0]
            Mtot = mass[inR].sum()
            Mshell = mass[inShell].sum()
            energy += GRAVITY * Mtot * Mshell / rmax
            gMtot = gmass[gasinR].sum()
            gMshell = gmass[gasinShell].sum()
            density = 3 * Mtot / (4*numpy.pi * rmax**3)
            gdensity = 3 * gMtot / (4*numpy.pi * rmax**3)
            delta = density/background_density
            shell_vol = 4/3 * numpy.pi * (rmax**3 - old_r**3)
            rhoShell = Mshell / shell_vol
            grhoShell = gMshell / shell_vol
            dens_avg = dens[gasinShell].mean()
            tshell = pyg.analyze.reject_outliers(temp[gasinShell]).mean()
            vr = vel[gasinShell, 0].mean()
            vtheta = vel[gasinShell,1].mean()
            vphi = vel[gasinShell,2].mean()
            vr_sigma = vel[gasinShell, 0].std()
            vtheta_sigma = vel[gasinShell,1].std()
            vphi_sigma = vel[gasinShell,2].std()
            if verbose:
                rpc = rmax/3.08568e18
                Msun = Mtot/1.989e33
                print 'R = %.2e pc' %rpc,
                print 'Mass enclosed: %.2e' %Msun,
                print 'Energy: %.3e' %energy,
                print 'delta: %.3f' %delta
            halo_properties.append((redshift,rmax,delta,Mshell,gMshell,-energy,
                                    density,gdensity,rhoShell,grhoShell,dens_avg,
                                    tshell,vr,vtheta,vphi,
                                    vr_sigma,vtheta_sigma,vphi_sigma,n))
            old_n = n
            old_r = rmax
        rmax *= r_multiplier
    return halo_properties

def chunks(l, n):
    """ Yield successive n-sized chunks from l.
    """
    for i in xrange(0, len(l), n):
        yield l[i:i+n]

#===============================================================================
if __name__ == '__main__':
    sim = pyg.sim.Simulation(sys.argv[1])
    sim.units.set_length('cm')
    sim.units.set_mass('g')
    sim.units.set_velocity('cgs')
    sim.units.set_density('cgs')

    keys = sim.snapfiles.keys()
    keys.sort()
    all_frames = []
    for chunk in chunks(keys,25):
        frames = [radial_properties(sim,snap) for snap in chunk]
        all_frames += frames
        hprops = pd.concat(all_frames)

        store = pd.HDFStore(sim.plotpath +'/'+ sim.name+'/radial_properties.hdf5')
        store[sim.name] = hprops
        store.close()

