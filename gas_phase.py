#!/usr/bin/env python
# gas_phase.py
# Jacob Hummel

import os
import sys
import glob
import numpy
import Queue
import subprocess
import multiprocessing as mp
from matplotlib import pyplot
import pyGadget

#===============================================================================
def load_snapshot(path):
    snapshot = pyGadget.snapshot.File(path)
    masses = snapshot.gas.get_masses()
    snapshot.gas.load_number_density()
    snapshot.gas.load_electron_fraction()
    snapshot.gas.calculate_temperature()
    snapshot.close()
    # Refine
    minimum = numpy.amin(masses)
    refined = numpy.where(masses <= minimum)[0]
    snapshot.gas.ndensity = snapshot.gas.ndensity[refined]
    snapshot.gas.temp = snapshot.gas.temp[refined]
    snapshot.gas.gamma = snapshot.gas.gamma[refined]
    snapshot.gas.h2frac = snapshot.gas.h2frac[refined]
    snapshot.gas.electron_frac = snapshot.gas.electron_frac[refined]
    
    # Cleanup to save memory
    del snapshot.gas.masses
    del snapshot.gas.abundances
    del snapshot.gas.energy

    return snapshot

def plot_gprops(snapshot,wpath):
    redshift = snapshot.header.Redshift
    dens = snapshot.gas.get_number_density()
    temp = snapshot.gas.get_temperature()
    electronfrac = snapshot.gas.get_electron_fraction()
    h2frac = snapshot.gas.get_H2_fraction()
    gamma = snapshot.gas.get_gamma()

    ### Create Plot!
    ### All plots vs density
    fig = pyplot.figure(1,(12,10))
    fig.clf()
    ax0 = fig.add_subplot(221)
    ax1 = fig.add_subplot(222)
    ax2 = fig.add_subplot(223)
    ax3 = fig.add_subplot(224)
    density_labels = (1e-2,1e0,1e2,1e4,1e6,1e8,1e10,1e12)

    # Temperature
    ax0.hexbin(dens,temp,gridsize=500,bins='log',xscale='log',yscale='log',
               mincnt=1)
    ax0.set_xscale('log')
    ax0.set_yscale('log')
    ax0.set_xlim(1e-2, 1e12)
    ax0.set_ylim(7, 1e4)
    ax0.set_xticks(density_labels)
    ax0.set_xlabel('n [cm$^{-3}$]')
    ax0.set_ylabel('Temperature [K]')

    # Free electron fraction
    ax1.hexbin(dens,electronfrac,gridsize=500,bins='log',xscale='log',
               yscale='log',mincnt=1)
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlim(1e-2, 1e12)
    ax1.set_ylim(1e-12, 1e-2)
    ax1.set_xticks(density_labels)
    ax1.set_xlabel('n [cm$^{-3}$]')
    ax1.set_ylabel('f$_{e^-}$')

    # Molecular Hydrogen Fraction
    ax2.hexbin(dens, h2frac,gridsize=500,bins='log',xscale='log',yscale='log',
               mincnt=1)
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.set_xlim(1e-2, 1e12)
    ax2.set_ylim(1e-7,2)
    ax2.set_xticks(density_labels)
    ax2.set_xlabel('n [cm$^{-3}$]')
    ax2.set_ylabel('f$_{H_2}$')

    # Adiabatic exponent
    ax3.hexbin(dens, gamma,gridsize=500,bins='log',xscale='log',yscale='log',
               mincnt=1)
    ax3.set_xscale('log')
    ax3.set_xlim(1e-2, 1e12)
    ax3.set_ylim(1,2)
    ax3.set_xticks(density_labels)
    ax3.set_xlabel('n [cm$^{-3}$]')
    ax3.set_ylabel('$\gamma_{ad}$')

    title = fig.suptitle('Redshift: %.3f' %(redshift,))
    fig.subplots_adjust(top=0.94, left=0.085, right=.915)
    pyplot.savefig(wpath+'-gprops.png', 
                   bbox_extra_artists=(title,),
                   bbox_inches='tight')

#===============================================================================
def multitask(path,write_dir,start,stop):
    maxprocs = mp.cpu_count()
    # Hack to take advantage of larger memory on r900 machines
    if 'r900' not in subprocess.check_output(['uname','-n']):
        maxprocs /= 2
    file_queue = Queue.Queue()
    data_queue = Queue.Queue(maxprocs/4)
    pyGadget.snapshot.Loader(load_snapshot, file_queue, data_queue).start()
    for snap in xrange(start,stop+1):
        fname = path + '{:0>3}'.format(snap)+'.hdf5'
        file_queue.put((fname,))
    file_queue.put(None)

    done = False
    while not done:
        jobs = []
        for i in xrange(maxprocs):
            snapshot = data_queue.get()
            if snapshot is None:
                done = True
                break
            else:
                wd = write_dir+'{:0>4}'.format(snapshot.number)
                p = mp.Process(target=plot_gprops, args=(snapshot, wd))
                jobs.append(p)
                p.start()
        print '\nClearing Queue!\n'
        for process in jobs:
            process.join()
        print '\nQueue Cleared!\n'
    for process in jobs:
        process.join()

#===============================================================================
if __name__ == '__main__':
    pyplot.ioff()
    if ((len(sys.argv) not in [2,3,4]) or (sys.argv[1] == '-h')):
        print 'Usage::'
        print '   Option 1: python gas_props.py [simulation name] '\
            '[beginning snapshot] [final snapshot]'
        print '   Option 2: python gas_props.py [simulation name] '\
            '[single snapshot]'
        print '   Option 3: python gas_props.py [simulation name] '\
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
        plot_gprops(snapshot,write_dir+'{:0>4}'.format(snap))
        
    elif len(sys.argv) == 4:
        start = int(sys.argv[2])
        stop = int(sys.argv[3])
        multitask(path,write_dir,start,stop)

    else:
        files0 = glob.glob(path+'???.hdf5')
        files1 = glob.glob(path+'1???.hdf5')
        files0.sort()
        files1.sort()
        files = files0 + files1
        start = int(files[0][-8:-5])
        stop = int(files[-1][-8:-5])
        multitask(path,write_dir,start,stop)

