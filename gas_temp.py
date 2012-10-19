# gas_temp.py
# Jacob Hummel

import os
import sys
import numpy
from matplotlib import pyplot
import pyGadget
#===============================================================================

def plot_temp(path, write_dir):
    snapshot = pyGadget.snapshot.load(path)
    h = snapshot.header.HubbleParam # H = 100*h
    redshift = snapshot.header.Redshift
    particle_mass = snapshot.gas.get_masses()
    dens = snapshot.gas.get_number_density()
    h2frac = snapshot.gas.get_H2_fraction()
    temp = snapshot.gas.get_temperature()
    snapshot.close()

    # Refine
    print 'Refining...'
    minimum = numpy.amin(particle_mass)
    stride = 100
    refined = numpy.where(particle_mass <= minimum)[0]#[0:-1:stride]
    dens = dens[refined]
    h2frac = h2frac[refined]
    temp = temp[refined]

    # Plot!
    print 'Plotting...'
    fig = pyplot.figure(1,(15,10))
    fig.clf()
    ax = fig.add_subplot(111)
    ax.scatter(dens, temp, s=.75, c='k', linewidths=0.0)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(2e-3, 1e13)
    ax.set_ylim(10, 2e4)

    pyplot.title('Redshift: %.2f' %(redshift,))
    ax.set_xlabel('n [cm^-3]')
    ax.set_ylabel('Temperature [K]')
    snap = path[-8:-5]
    pyplot.savefig(write_dir+snap+'-temp.png',
                   bbox_inches='tight')
#===============================================================================

if __name__ == '__main__':
    #pyplot.ioff()
    if len(sys.argv) < 4:
        print 'Usage: python gas_temp.py [simulation name] '\
            '[beginning snapshot] [final snapshot]'
        sys.exit()

    simulation = sys.argv[1]
    path = os.getenv('HOME')+'/sim/'+simulation+'/snapshot_'
    write_dir = os.getenv('HOME')+'/data/simplots/'+simulation+'/'

    start = int(sys.argv[2])
    stop = int(sys.argv[3])+1
    for snap in xrange(start,stop):
        fname = path + '{:0>3}'.format(snap)+'.hdf5'
        print 'loading', fname
        plot_temp(fname, write_dir)
