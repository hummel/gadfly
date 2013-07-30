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
        snapshot.gas.load_HD_fraction()
    snapshot.close()

    # Refine
    minimum = numpy.amin(masses)
    refined = numpy.where(masses <= minimum)[0]
    snapshot.gas.ndensity = snapshot.gas.ndensity[refined]
    snapshot.gas.temp = snapshot.gas.temp[refined]
    if key == 'frac':
            snapshot.gas.h2frac = snapshot.gas.h2frac[refined]
            snapshot.gas.HDfrac = snapshot.gas.HDfrac[refined]
            snapshot.gas.electron_frac = snapshot.gas.electron_frac[refined]

    # Cleanup to save memory
    del snapshot.gas.masses
    del snapshot.gas.abundances
    del snapshot.gas.gamma
    if key == 'temp':
        del snapshot.gas.energy
        del snapshot.gas.h2frac

    return snapshot

def prep_figure():
    fig = pyplot.figure(figsize=(12,8))
    fig.clf()
    return fig

def save_figure(fig, snapshot, wpath):
    redshift = snapshot.header.Redshift
    title = fig.suptitle('Redshift: %.2f' %(redshift,))
    fig.savefig(wpath,bbox_extra_artists=(title,),bbox_inches='tight')

def plot_temp(ax,snapshot):
    dens = snapshot.gas.get_number_density()
    temp = snapshot.gas.get_temperature()
    ax.hexbin(dens,temp,gridsize=500,bins='log',xscale='log',yscale='log',
              mincnt=1)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(2e-3, 1e12)
    ax.set_ylim(10, 2e4)
    ax.set_xlabel('n [cm^-3]')
    ax.set_ylabel('Temperature [K]')
    return ax

def plot_electron_frac(ax, snapshot):
    dens = snapshot.gas.get_number_density()
    efrac = snapshot.gas.get_electron_fraction()
    ax.hexbin(dens,efrac,gridsize=500,bins='log',xscale='log',
               yscale='log',mincnt=1)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(1e-2, 1e12)
    ax.set_ylim(1e-12, 1e-2)
    ax.set_xlabel('n [cm$^{-3}$]')
    ax.set_ylabel('f$_{e^-}$')
    return ax

def plot_h2frac(ax, snapshot):
    dens = snapshot.gas.get_number_density()
    h2frac = snapshot.gas.get_H2_fraction()
    ax.hexbin(dens, h2frac,gridsize=500,bins='log',xscale='log',yscale='log',
               mincnt=1)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(1e-2, 1e12)
    ax.set_ylim(1e-7,2)
    ax.set_xlabel('n [cm$^{-3}$]')
    ax.set_ylabel('f$_{H_2}$')
    return ax

def plot_HDfrac(ax, snapshot):
    dens = snapshot.gas.get_number_density()
    HDfrac = snapshot.gas.get_HD_fraction()
    ax.hexbin(dens,HDfrac,gridsize=500,bins='log',xscale='log',yscale='log',
               mincnt=1)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(1e-2, 1e12)
    ax.set_ylim(1e-11,1e-4)
    ax.set_xlabel('n [cm$^{-3}$]')
    ax.set_ylabel('f$_{HD}$')

def figure_temp(snapshot,wpath):
    fig = prep_figure()
    ax = fig.add_subplot(111)
    ax = plot_temp(ax, snapshot)
    save_figure(fig, snapshot, wpath+'-temp.png')

def figure_gas_fraction(snapshot,wpath):
    ### Create Plot!
    ### All plots vs density
    fig = prep_figure()
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(223)
    ax4 = fig.add_subplot(224)
    axes = [ax1,ax2,ax3,ax4]

    # Temperature
    ax1 = plot_temp(ax1,snapshot)
    # Free electron fraction
    ax2 = plot_electron_frac(ax2,snapshot)
    # Molecular Hydrogen Fraction
    ax3 = plot_h2frac(ax3,snapshot)
    # HD fraction
    ax4 = plot_HDfrac(ax4,snapshot)
    
    density_labels = (1e-2,1e0,1e2,1e4,1e6,1e8,1e10,1e12)
    for ax in axes:
        ax.set_xticks(density_labels)
    fig.subplots_adjust(top=0.94, left=0.085, right=.915)
    save_figure(fig, snapshot, wpath+'-frac.png')

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
        file_queue.put((fname,key))
    file_queue.put(None)
    if key == 'temp':
        plot_func = figure_temp
    elif key == 'frac':
        plot_func = figure_gas_fraction

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
            figure_temp(snapshot,write_dir+'{:0>4}'.format(snap))
        elif key == 'frac':
            figure_gas_fraction(snapshot,write_dir+'{:0>4}'.format(snap))
        
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
        multitask(key,path,write_dir,start,stop)

