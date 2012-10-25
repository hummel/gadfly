# gas_temp.py
# Jacob Hummel

import os
import sys
import glob
import numpy
from matplotlib import pyplot
import pyGadget
#===============================================================================

def load_snapshot(path):
    snapshot = pyGadget.snapshot.load(path)
    snapshot.gas.load_masses()
    snapshot.gas.load_number_density()
    snapshot.gas.load_internal_energy()
    snapshot.gas.load_gamma()
    snapshot.gas.load_H2_fraction()
    snapshot.close()
    return snapshot

def plot_temp(snapshot):
    h = snapshot.header.HubbleParam
    redshift = snapshot.header.Redshift
    particle_mass = snapshot.gas.get_masses()
    dens = snapshot.gas.get_number_density()
    temp = snapshot.gas.get_temperature()

    # Refine
    print 'Refining...'
    minimum = numpy.amin(particle_mass)
    stride = 100
    refined = numpy.where(particle_mass <= minimum)[0]#[0:-1:stride]
    dens = dens[refined]
    temp = temp[refined]

    # Plot!
    print 'Plotting...'
    fig = pyplot.figure(figsize=(15,10))
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
    return fig

def save_plot(fig,path):
    pyplot.savefig(path+'-temp.png',
                   bbox_inches='tight')
#===============================================================================
def multitask(path,write_dir,start,stop):
    for snap in xrange(start,stop):
        fname = path + '{:0>3}'.format(snap)+'.hdf5'
        print 'loading', fname
        snapshot = load_snapshot(fname)
        fig = plot_temp(snapshot)
        save_plot(fig,write_dir+str(snapshot.number))

#===============================================================================

if __name__ == '__main__':
    pyplot.ioff()
    if ((len(sys.argv) not in [2,3,4]) or (sys.argv[1] == '-h')):
        print 'Usage::'
        print '   Option 1: python gas_temp.py [simulation name] '\
            '[beginning snapshot] [final snapshot]'
        print '   Option 2: python gas_temp.py [simulation name] '\
            '[single snapshot]'
        print '   Option 3: python gas_temp.py [simulation name] '\
            '(creates a plot for every snapshot)'
        sys.exit()

    simulation = sys.argv[1]
    path = os.getenv('HOME')+'/sim/'+simulation+'/snapshot_'
    write_dir = os.getenv('HOME')+'/data/simplots/'+simulation+'/'

    if len(sys.argv) == 3:
        snap = sys.argv[2]
        fname = path + '{:0>3}'.format(snap)+'.hdf5'
        print 'loading', fname
        snapshot = load_snapshot(fname)
        fig = plot_temp(snapshot)
        save_plot(fig,write_dir+snap)
        
    elif len(sys.argv) == 4:
        start = int(sys.argv[2])
        stop = int(sys.argv[3])+1
        multitask(path,write_dir,start,stop)

    else:
        files = glob.glob(path+'???.hdf5')
        files.sort()
        start = int(files[0][-8:-5])
        stop = start + len(files)
        multitask(path,write_dir,start,stop)

