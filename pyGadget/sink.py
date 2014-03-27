# sinks.py
# Jacob Hummel
"""
Classes and routines for analyzing sink data output by gadget.
"""
import os
import sys
import warnings
import numpy
from astropy.io import ascii
import sqlite3
#from numba import autojit

import units
import analyze
import constants
#===============================================================================

class Sink(object):
    def __init__(self,**properties):
        super(Sink,self).__init__()
        self.mass = properties.pop('m', None)
        self.pos = properties.pop('pos', None)
        self.velocity = properties.pop('vel', None)
        self.radius = properties.pop('r', None)
        self.energy = properties.pop('e', None)
        self.pressure = properties.pop('p', None)
        self.npart_acc = properties.pop('n', None)
        self.pid = properties.pop('pid', None)
        self.index = properties.pop('index', None)

        self.x = properties.pop('x', self.pos[0])
        self.y = properties.pop('y', self.pos[1])
        self.z = properties.pop('z', self.pos[2])
        self.vx = properties.pop('vx', self.velocity[0])
        self.vy = properties.pop('vy', self.velocity[1])
        self.vz = properties.pop('vz', self.velocity[2])

    def update_coordinates(self, x,y,z):
        self.x = x[self.index]
        self.y = y[self.index]
        self.z = z[self.index]

    def update_velocities(self, x,y,z):
        self.vx = x[self.index]
        self.vy = y[self.index]
        self.vz = z[self.index]

    def update_index(self, index):
        self.index = index

class SinkData(object):
    def __init__(self,path):
        super(SinkData,self).__init__()
        ### Read in the data
        try:
            sinkdata = ascii.read(path+'/sinkdat')
        except IOError:
            raise IOError("Specified sinkmasses file not found!")
        try:
            sinkmasses = ascii.read(path+'/sinkmasses')
        except IOError:
            raise IOError("Specified sinkmasses file not found!")

        self.time = sinkdata['col1']
        self.npart_acc = sinkdata['col2']
        self.radius = sinkdata['col3']
        self.part_internal_energy = sinkdata['col4']
        self.entropy = sinkdata['col5']
        self.part_id = sinkdata['col6']
        self.sink_id  = sinkdata['col7']
        self.pressure = sinkdata['col8']
        self.a = self.time # Scale Facor

        # Restrict to real sinks
        IDs = numpy.unique(sinkmasses['col2'])
        real = numpy.in1d(self.sink_id, IDs)
        for key in vars(self).keys():
            vars(self)[key] = vars(self)[key][real]

        h = 0.7 #Hubble Parameter
        h2 = h*h
        a3 = self.a**3
        ### Convert units
        self.time = self.time*units.Time_yr
        # npart_acc is a simple integer (no units)
        self.pressure = self.pressure*units.Pressure_cgs*h2/(a3**1.4)
        self.radius = self.radius*units.Length_AU*h

        good = numpy.where(self.radius > 10)[0]
        for key in vars(self).keys():
            vars(self)[key] = vars(self)[key][good]

class SinkHistory(SinkData):
    '''
    Select sink data for a single sink particle.

    possible options:
    nform: select the n(th) sink to form.
    ID: select sink by ID.
    '''
    def __init__(self, path, nform=None, id_=None):
        super(SinkHistory,self).__init__(path)
        unique = numpy.unique(self.sink_id)

        if((nform is None) and (id_ is None)):
            print "No sink specified: Selecting first sink to form..."
            nform = 1
        if nform:
            print "Key set: nform =", nform
            if id_: warnings.warn("nform key overrides id_")
            # Select n(th) sink to form
            new = []
            i = 0
            while len(new) < nform:
                if self.sink_id[i] not in new:
                    new.append(self.sink_id[i])
                i += 1
            id_ = new[-1]
        elif id_ is not None:
            print "Key set: id_ =", id_
        else:
            raise RuntimeError("Execution should not have reached this point!")
        print "Using sink ID", id_
        
        # Restrict to a single sink
        lines = numpy.where(self.sink_id == id_)
        for key in vars(self).keys():
            vars(self)[key] = vars(self)[key][lines]

        # Select final output for each timestep
        tsteps = numpy.unique(self.time)
        selection = []
        for t in tsteps:
            times = numpy.where(self.time == t)[0]
            selection.append(times[-1])
        for key in vars(self).keys():
                vars(self)[key] = vars(self)[key][selection]
        self.sink_id = id_
        self.nform = nform

        # Calculate sink mass at each timestep
        self.mass = numpy.zeros_like(self.time)
        for i in xrange(self.mass.size):
            self.mass[i] = 0.015*self.npart_acc[:i+1].sum()

        # Finally, record total number of sinks found.
        self.all_ids = unique

#===============================================================================
class AccretionDisk(object):
    def __init__(self, sim, sink, **kwargs):
        super(AccretionDisk, self).__init__()
        self.sim = sim
        self.sink = sink
        default = sim.plotpath +'/'+ sim.name
        if not os.path.exists(default):
            os.makedirs(default)
        dbfile = kwargs.pop('dbfile', default+'/disk{}.db'.format(sink.nform))
        self.db = sqlite3.connect(dbfile)
        self.c = self.db.cursor()

    def load(self, snap, *dprops, **kwargs):
        verbose = kwargs.pop('verbose', False)
        kwargs['verbose'] = verbose
        fields = ''
        if dprops:
            for dprop in dprops[:-1]:
                fields += dprop + ', '
            fields += dprops[-1]
        else:
            fields = '*'
        table = 'snapshot{:0>4}'.format(snap)
        command = "SELECT " + fields + " FROM " + table
        print '"'+command+'"'
        try:
            self.c.execute(command)
        except sqlite3.OperationalError:
            print "Warning: Error loading requested accretion disk data!"
            print "Recalculating..."
            self.populate(snap, **kwargs)
            self.c.execute(command)

        self.data = numpy.asarray(self.c.fetchall())
        if self.data.size < 1:
            print "Warning: Requested accretion disk data does not exist!"
            print "Calculating..."
            self.populate(snap, **kwargs)
            self.c.execute(command)
            self.data = numpy.asarray(self.c.fetchall())

    def populate(self, snap, **kwargs):
        self.sim.units.set_length('cm')
        self.sim.units.set_mass('g')
        self.sim.units.set_velocity('cgs')
        snapshot = self.sim.load_snapshot(snap, track_sinks=True)
        diskprops = disk_properties(snapshot, self.sink.sink_id, **kwargs)
        snapshot.gas.cleanup()
        snapshot.close()
        print "Saving table containing", len(diskprops), "entries to database."
        table = 'snapshot{:0>4}'.format(snap)
        create = ("CREATE TABLE " + table +
                  "(redshift real, "\
                  "radius real, "\
                  "density real, "\
                  "total_mass real, "\
                  "mass real, "\
                  "vrot real, "\
                  "vkepler real, "\
                  "tff real, "\
                  "T real, "\
                  "Tavg real, "\
                  "cs real, "\
                  "H real, "\
                  "Lj real, "\
                  "Mj real, "\
                  "npart integer)")
        insert = ("INSERT INTO " + table +
                  "(redshift, "\
                  "radius, "\
                  "density, "\
                  "total_mass, "\
                  "mass, "\
                  "vrot, "\
                  "vkepler, "\
                  "tff, "\
                  "T, "\
                  "Tavg, "\
                  "cs, "\
                  "H, "\
                  "Lj, "\
                  "Mj, "\
                  "npart) "\
                  "VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)")
        try:
            self.c.execute(create)
        except sqlite3.OperationalError:
            self.c.execute("DROP TABLE " + table)
            self.c.execute(create)
        self.c.executemany(insert, diskprops)
        self.db.commit()

#===============================================================================
def disk_properties(snapshot, sink_id, **kwargs):
    r_start = kwargs.pop('r_start', 50*1.49597871e13)
    r_step = kwargs.pop('r_step', 1.49597871e14)
    r_multiplier = kwargs.pop('multiplier', 1.2)
    verbose = kwargs.pop('verbose', True)
    n_min = kwargs.pop('n_min', 32)
    dens_lim = kwargs.pop('dens_lim', 1e6)
    orb_crit = kwargs.pop('orb_crit', 1.)
    disky_crit = kwargs.pop('disky_crit', .9)
    kepler_crit = kwargs.pop('kepler_crit', .8)

    redshift = snapshot.header.Redshift
    length_unit = 'cm'
    mass_unit = 'g'
    velocity_unit = 'cgs'
    snapshot.gas.units.set_velocity(velocity_unit)
    xyz = snapshot.gas.get_coords(length_unit)
    uvw = snapshot.gas.get_velocities()
    snapshot.update_sink_frame_ofR(xyz, uvw)

    dens = snapshot.gas.get_number_density('cgs')
    mass = snapshot.gas.get_masses(mass_unit)
    temp = snapshot.gas.get_temperature()
    csound = snapshot.gas.get_sound_speed()

    i = 0
    try:
        while snapshot.sinks[i].pid != sink_id:
            i += 1
    except IndexError:
        raise IndexError("Sink {} has not yet formed!".format(sink_id))
    else:
        sink = snapshot.sinks[i]
        sinkpos = (sink.x, sink.y, sink.z)
        sinkvel = (sink.vx, sink.vy, sink.vz)
        pos = snapshot.gas.get_coords(system='spherical',view='face',
                                      center=sinkpos, vcenter=sinkvel)
    vel = snapshot.gas.get_velocities(system='spherical')
    xyz = snapshot.gas.get_coords()
    uvw = snapshot.gas.get_velocities()
    snapshot.update_sink_frame_ofR(xyz, uvw)
    pos = numpy.nan_to_num(pos)
    vel = numpy.nan_to_num(vel)

    L9 = analyze.total_angular_momentum(*analyze.data_slice(dens > 1e9,
                                                            xyz,uvw,mass))
    uL9 = L9 / numpy.linalg.norm(L9)

    if dens_lim:
        arrs = [dens,pos,vel,xyz,uvw,mass,temp,csound]
        dslice = dens > dens_lim
        dens,pos,vel,xyz,uvw,mass,temp,csound = analyze.data_slice(dslice,
                                                                    *arrs)
    L = analyze.angular_momentum(xyz,uvw,mass)
    uL = L / numpy.linalg.norm(L, axis=1)[:, numpy.newaxis]
    disky = numpy.where((numpy.abs(vel[:,2])/numpy.abs(vel[:,0]) > orb_crit)
                    & (uL.dot(uL9) > disky_crit))[0]
    print disky.size, '"disky" particles',
    print '({}) percent'.format(float(disky.size)/dens.size * 100)

    inf = numpy.where(numpy.abs(vel[:,0])
                      / numpy.linalg.norm(vel, axis=1) > .5)[0]
    print inf.size, 'infalling particles',
    print '({}) percent'.format(float(inf.size)/dens.size * 100)

    orb = numpy.where((numpy.abs(vel[:,2])/numpy.abs(vel[:,0]) > orb_crit)
                   & (uL.dot(uL9) < disky_crit))[0]
    print orb.size, '"orbiting" particles',
    print '({}) percent'.format(float(orb.size)/dens.size * 100)

    notorb = numpy.where(numpy.abs(vel[:,2]) < numpy.abs(vel[:,0]))[0]
    print notorb.size, '"not orbiting" particles',
    print '({}) percent'.format(float(notorb.size)/dens.size * 100)

    print 'Data loaded.  Analyzing...'
    disk_properties = []
    vrot = vk = 1.
    n = old_n = r0 = 0
    r1 = r_start
    print "Starting at {:.2e} AU".format(r1/1.49e13)
    r2d = numpy.sqrt(xyz[:,0]**2 + xyz[:,1]**2)
    data = []
    while vrot/vk > kepler_crit  and n < disky.size:
        inR = numpy.where(pos[disky,0] <= r1)[0]
        n = inR.size
        if n > old_n and n > n_min:
            annulus = numpy.where((r2d[disky] > r0) & (r2d[disky] <= r1))[0]
            annulus = disky[annulus]
            radii = r2d[annulus]
            radius = radii.mean()
            vrot = vel[annulus, 2].mean()
            T = analyze.reject_outliers(temp[annulus]).mean()
            cs = analyze.reject_outliers(csound[annulus]).mean()
            H = cs * radius / vrot
            Mcyl = mass[annulus].sum()
            if annulus.size > 0:
                zmax = numpy.abs(xyz[annulus,2]).max()
            else:
                zmax = 0.0
            density = dens[annulus].mean()
            if numpy.isnan(density):
                density = dens.max()
            mdensity = density * constants.m_H / constants.X_h
            tff = numpy.sqrt(3 * numpy.pi / 32 / constants.GRAVITY / mdensity)

            tavg = analyze.reject_outliers(temp[inR]).mean()
            Mtot = mass[numpy.where(pos[:,0] <= radius)[0]].sum()
            vk = numpy.sqrt(6.6726e-8 * Mtot / radius) / 1e5

            Lj = cs*tff #Jeans Length
            Mj = mdensity * (4*numpy.pi/3) * Lj**3 / 1.989e33 #Jeans Mass
            H = cs * radius / vrot #Disk Scale Height

            rau = radius/1.49597871e13
            Msun = Mtot/1.989e33
            vrot = vrot / 1e5
            disk_properties.append((redshift,rau,density,Msun,Mcyl,
                                    vrot,vk,tff,T,tavg,cs,H,Lj,Mj,n))
            if verbose:
                print 'R = %.2e AU' %rau,
                print 'Mass enclosed: %.2e' %Msun,
                print 'density: %.3e' %density,
                print 'vrot: %.3e' %vrot,
                print 'npart: {}'.format(n)
        old_n = n
        r0 = r1
        r1 *= r_multiplier
    print 'snapshot', snapshot.number, 'analyzed.'
    return disk_properties
