# analyze.py
# Jacob Hummel
import glob
import numpy
import statistics
from . import units
from . import constants

#===============================================================================
def find_center(x, y, z, dens, centering='avg', dens_limit=1e11, nparticles=100,
		verbose=True):
	if centering == 'avg':
		hidens = numpy.where(dens >= dens_limit)[0]
		while hidens.size < nparticles:
		    dens_limit /= 2
		    hidens = numpy.where(dens >= dens_limit)[0]
		if verbose:
		    print ('Center averaged over all particles with density '\
			       'greater than %.2e particles/cc' %dens_limit)
		#Center on highest density clump, rejecting outliers:
		x = x - numpy.average(statistics.reject_outliers(x[hidens]))
		y = y - numpy.average(statistics.reject_outliers(y[hidens]))
		z = z - numpy.average(statistics.reject_outliers(z[hidens]))
	elif centering == 'max':
		center = dens.argmax()
		x = x - x[center]
		y = y - y[center]
		z = z - z[center]
        return x,y,z

#===============================================================================
def compile_halos(directory):
    print directory+'haloz/????.npy'
    files = glob.glob(directory+'haloz/????.npy')
    files.sort()
    print files
    data = []
    for f in files:
        print f
        halo = numpy.load(f)
        print halo
        print halo.size
        if halo.size > 0:
            data.append(halo)

    array_lengths = [x.shape[0] for x in data]
    maxL = max(array_lengths)
    total = len(data)
    for i in range(total):
        data[i].resize([maxL,13], refcheck=False)
    datarray = numpy.concatenate([x for x in data])
    datarray = datarray.reshape(total,maxL,7)
    return datarray

#===============================================================================
def halo_properties(snapshot, 
                    r_start=3.08568e17, r_multiplier=1.01, verbose=True):
    h = snapshot.header.HubbleParam
    a = snapshot.header.ScaleFactor
    redshift = snapshot.header.Redshift

    length_unit = units.Length_cm
    dm_mass = snapshot.dm.get_masses()
    dm_pos = snapshot.dm.get_coords(length_unit)
    gas_mass = snapshot.gas.get_masses()
    gas_pos = snapshot.gas.get_coords(length_unit)
    dens = snapshot.gas.get_density()
    temp = snapshot.gas.get_temperature()

    gasx = gas_pos[:,0]
    gasy = gas_pos[:,1]
    gasz = gas_pos[:,2]
    dmx = dm_pos[:,0]
    dmy = dm_pos[:,1]
    dmz = dm_pos[:,2]
    del dm_pos

    mass = numpy.concatenate((gas_mass, dm_mass))
    del gas_mass
    del dm_mass
    x = numpy.concatenate((gasx,dmx))
    y = numpy.concatenate((gasy,dmy))
    z = numpy.concatenate((gasz,dmz))
    del dmx,dmy,dmz

    x,y,z = find_center(x,y,z,dens,'max',verbose=verbose)
    gasx,gasy,gasz = find_center(gasx,gasy,gasz,dens,'max',verbose=verbose)
    del dens
    del gas_pos
    r = numpy.sqrt(numpy.square(x) + numpy.square(y) + numpy.square(z))
    gasr = numpy.sqrt(numpy.square(gasx)
		      + numpy.square(gasy)
		      + numpy.square(gasz))
    del x,y,z,gasx,gasy,gasz
    
    halo_properties = []
    n = old_n = old_r = density = energy = 0
    background_density = .27 * 9.31e-30 * (1+redshift)**3 #Omega_m * rho_crit(z)
    rmax = r_start
    while (density > 180 * background_density or n < 50):
        inR = numpy.where(r <= rmax)[0]
        gasinR = numpy.where(gasr <= rmax)[0]
        n = inR.size
        if n > old_n:
            inShell = numpy.where(r[inR] > old_r)[0]
            gasinShell = numpy.where(r[gasinR] > old_r)[0]
            rpc = rmax/3.08568e18
            Mtot = mass[inR].sum()
	    Mshell = mass[inShell].sum()
            solar_masses = Mtot/1.989e33
            density = 3 * Mtot / (4*numpy.pi * rmax**3)
            delta = density/background_density
	    tshell = temp[gasinShell].mean()
	    tavg = temp[gasinR].mean()
	    tff = numpy.sqrt(3*numpy.pi/32/constants.GRAVITY/density)
	    cs = numpy.sqrt(constants.k_B * tavg / constants.m_H)
	    Lj = cs*tff
	    Mj = density * (4*numpy.pi/3) * Lj**3 / 1.989e33
	    energy += constants.GRAVITY * Mtot * Mshell / rmax
            if verbose: 
                print 'R = %.2e pc' %rpc,
		print 'Mass enclosed: %.2e' %solar_masses,
                print 'Energy: %.3e' %energy,
                print 'delta: %.3f' %delta
            if delta >= 178.0:
                halo_properties.append((redshift,rpc,delta,solar_masses,density,
                                        tavg,tshell,tff,cs,Lj,Mj,-energy,n))
            old_n = n
	    old_r = rmax
        rmax *= r_multiplier
    
    del r
    print 'snapshot', snapshot.number, 'analyzed.'
    return numpy.asarray(halo_properties)

