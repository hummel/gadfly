# halo.py
# Jacob Hummel
import os
import numpy
import sqlite3
import analyze
import constants

class Halo(object):
    def __init__(self, sim, **kwargs):
        super(Halo, self).__init__()
        self.sim = sim
        default = sim.plotpath +'/'+ sim.name
        if not os.path.exists(default):
            os.makedirs(default)
        dbfile = kwargs.pop('dbfile', default+'/halo.db')
        self.db = sqlite3.connect(dbfile)
        self.c = self.db.cursor()

    def load(self, snap, *hprops):
        fields = ''
        if hprops:
            for hprop in hprops[:-1]:
                fields += hprop + ', '
            fields += hprops[-1]
        else: 
            fields = '*'
        table = 'snapshot{:0>4}'.format(snap)
        try:
            self.c.execute("SELECT " + fields + " FROM " + table)
        except sqlite3.OperationalError:
            print "Warning: Halo data for this snapshot does not exist."
            print "Analyzing..."
            self.populate(snap, verbose=False)
            self.c.execute("SELECT " + fields + " FROM " + table)
            
        self.data = numpy.asarray(self.c.fetchall())

    def populate(self, snap, **kwargs):
        snapshot = self.sim.load_snapshot(snap)
        haloprops = radial_properties(snapshot, **kwargs)
        snapshot.close()
        table = 'snapshot{:0>4}'.format(snap)
        create = ("CREATE TABLE " + table +
                  "(redshift real, radius real, delta real, mass real,"\
                      " density real, Tavg real, Tshell real, tff real,"\
                      " cs real, Lj real, Mj real, gpe real, npart integer)")
        try:
            self.c.execute(create)
        except sqlite3.OperationalError:
            self.c.execute("DROP TABLE " + table)
            self.c.execute(create)
        insert = ("INSERT INTO " + table + "(redshift, radius, "\
                      "delta, mass, density, Tavg, Tshell, tff, cs, Lj, "\
                      "Mj, gpe, npart) VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?)")
        self.c.executemany(insert, haloprops)
        self.db.commit()

#===============================================================================
def radial_properties(snapshot, **kwargs):
    r_start = kwargs.pop('r_start', 3.08568e17)
    r_multiplier = kwargs.pop('multiplier', 1.01)
    verbose = kwargs.pop('verbose', True)
    n_min = kwargs.pop('n_min', 50)

    length_unit = 'cm'
    mass_unit = 'g'
    redshift = snapshot.header.Redshift
    dm_mass = snapshot.dm.get_masses(mass_unit)
    dm_pos = snapshot.dm.get_coords(length_unit)
    gas_mass = snapshot.gas.get_masses(mass_unit)
    gas_pos = snapshot.gas.get_coords(length_unit)
    dens = snapshot.gas.get_density('cgs')
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

    x,y,z = analyze.center_box(x,y,z, density=dens, centering='max',
                               verbose=verbose)
    gasx,gasy,gasz = analyze.center_box(gasx,gasy,gasz, density=dens,
                                        centering='max', verbose=verbose)
    del dens
    del gas_pos
    r = numpy.sqrt(numpy.square(x) + numpy.square(y) + numpy.square(z))
    gasr = numpy.sqrt(numpy.square(gasx)
		      + numpy.square(gasy)
		      + numpy.square(gasz))
    del x,y,z,gasx,gasy,gasz
    
    halo_properties = []
    n = old_n = old_r = density = energy = 0
    # background density:: Omega_m * rho_crit(z)
    background_density = .27 * 9.31e-30 * (1+redshift)**3 
    rmax = r_start
    while (density > 180 * background_density or n < n_min):
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
    return halo_properties


