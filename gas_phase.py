#!/usr/bin/env python
# gas_temp.py
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
def load_snapshot(path,key):
    if key not in ['temp', 'frac']:
        raise KeyError
    snapshot = pyGadget.snapshot.File(path)
    masses = snapshot.gas.get_masses()
    snapshot.gas.load_number_density()
    snapshot.gas.calculate_temperature()
    if key == 'frac':
        snapshot.gas.load_electron_fraction()
    snapshot.close()

    # Refine
    minimum = numpy.amin(masses)
    refined = numpy.where(masses <= minimum)[0]
    snapshot.gas.ndensity = snapshot.gas.ndensity[refined]
    snapshot.gas.temp = snapshot.gas.temp[refined]
    if key == 'frac':
            snapshot.gas.gamma = snapshot.gas.gamma[refined]
            snapshot.gas.h2frac = snapshot.gas.h2frac[refined]
            snapshot.gas.electron_frac = snapshot.gas.electron_frac[refined]

    # Cleanup to save memory
    del snapshot.gas.masses
    del snapshot.gas.abundances
    if key == 'temp':
        del snapshot.gas.energy
        del snapshot.gas.h2frac
        del snapshot.gas.gamma

    return snapshot

def prep_plot():
    fig = pyplot.figure(figsize=(12,8))
    fig.clf()
    return fig

def save_plot(snapshot, fig, wpath):
    redshift = snapshot.header.Redshift
    fig.suptitle('Redshift: %.2f' %(redshift,))
    fig.savefig(wpath,bbox_inches='tight')

def plot_temp(snapshot,wpath):
    dens = snapshot.gas.get_number_density()
    temp = snapshot.gas.get_temperature()
    fig = prep_plot()
    ax = fig.add_subplot(111)
    ax.hexbin(dens,temp,gridsize=500,bins='log',xscale='log',yscale='log',
              mincnt=1)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(2e-3, 1e12)
    ax.set_ylim(10, 2e4)
    ax.set_xlabel('n [cm^-3]')
    ax.set_ylabel('Temperature [K]')
    save_plot(snapshot, fig, wpath+'-temp.png')

def plot_gas_fraction(snapshot,wpath):
    dens = snapshot.gas.get_number_density()
    temp = snapshot.gas.get_temperature()
    electronfrac = snapshot.gas.get_electron_fraction()
    h2frac = snapshot.gas.get_H2_fraction()
    gamma = snapshot.gas.get_gamma()
    ### Create Plot!
    ### All plots vs density
    fig = prep_plot()
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(223)
    ax4 = fig.add_subplot(224)
    density_labels = (1e-2,1e0,1e2,1e4,1e6,1e8,1e10,1e12)
    # Temperature
    ax1.hexbin(dens,temp,gridsize=500,bins='log',xscale='log',yscale='log',
               mincnt=1)
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlim(1e-2, 1e12)
    ax1.set_ylim(7, 1e4)
    ax1.set_xticks(density_labels)
    ax1.set_xlabel('n [cm$^{-3}$]')
    ax1.set_ylabel('Temperature [K]')
    # Free electron fraction
    ax2.hexbin(dens,electronfrac,gridsize=500,bins='log',xscale='log',
               yscale='log',mincnt=1)
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.set_xlim(1e-2, 1e12)
    ax2.set_ylim(1e-12, 1e-2)
    ax2.set_xticks(density_labels)
    ax2.set_xlabel('n [cm$^{-3}$]')
    ax2.set_ylabel('f$_{e^-}$')
    # Molecular Hydrogen Fraction
    ax3.hexbin(dens, h2frac,gridsize=500,bins='log',xscale='log',yscale='log',
               mincnt=1)
    ax3.set_xscale('log')
    ax3.set_yscale('log')
    ax3.set_xlim(1e-2, 1e12)
    ax3.set_ylim(1e-7,2)
    ax3.set_xticks(density_labels)
    ax3.set_xlabel('n [cm$^{-3}$]')
    ax3.set_ylabel('f$_{H_2}$')

    # Adiabatic exponent
    ax4.hexbin(dens, gamma,gridsize=500,bins='log',xscale='log',yscale='log',
               mincnt=1)
    ax4.set_xscale('log')
    ax4.set_xlim(1e-2, 1e12)
    ax4.set_ylim(1,2)
    ax4.set_xticks(density_labels)
    ax4.set_xlabel('n [cm$^{-3}$]')
    ax4.set_ylabel('$\gamma_{ad}$')

    fig.subplots_adjust(top=0.94, left=0.085, right=.915)
    save_plot(snapshot, fig, wpath+'-gprops.png')

#===============================================================================
def multitask(key,path,write_dir,start,stop):
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
    if key == 'temp':
        plot_func = plot_temp
    elif key == 'frac':
        plot_func = plot_gas_fraction

    done = False
    while not done:
        jobs = []
        for i in xrange(maxprocs):
            snapshot = data_queue.get()
            if snapshot is None:
                done = True
                break
            else:
                wp = write_dir+'{:0>4}'.format(snapshot.number)
                p = mp.Process(target=plot_func, args=(snapshot, wp))
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
    if ((len(sys.argv) not in [3,4,5]) or (sys.argv[1] == '-h')):
        print 'Usage::'
        print '   Option 1: python gas_temp.py [key] [simulation name] '\
            '[beginning snapshot] [final snapshot]'
        print '   Option 2: python gas_temp.py [key] [simulation name] '\
            '[single snapshot]'
        print '   Option 3: python gas_temp.py [key] [simulation name] '\
            '(creates a plot for every snapshot)'
        sys.exit()

    key = sys.argv[1]
    if key not in ['temp','frac']:
        raise KeyError
    simulation = sys.argv[2]
    path = os.getenv('HOME')+'/sim/'+simulation+'/snapshot_'
    write_dir = os.getenv('HOME')+'/data/simplots/'+simulation+'/gas/'+key+'/'
    if not os.path.exists(write_dir):
        os.makedirs(write_dir)
                 
    if len(sys.argv) == 4:
        snap = sys.argv[3]
        fname = path + '{:0>3}'.format(snap)+'.hdf5'
        print 'loading', fname
        snapshot = load_snapshot(fname,key)
        if key == 'temp':
            plot_temp(snapshot,write_dir+'{:0>4}'.format(snap))
        elif key == 'frac':
            plot_gas_fraction(snapshot,write_dir+'{:0>4}'.format(snap))
        
    elif len(sys.argv) == 5:
        start = int(sys.argv[3])
        stop = int(sys.argv[4])
        multitask(key,path,write_dir,start,stop)

    else:
        files0 = glob.glob(path+'???.hdf5')
        files1 = glob.glob(path+'1???.hdf5')
        files0.sort()
        files1.sort()
        files = files0 + files1
        start = int(files[0][-8:-5])
        stop = int(files[-1][-8:-5])
        multitask(path,write_dir,start,stop)

