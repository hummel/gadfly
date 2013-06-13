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
def load_snapshot(path):
    snapshot = pyGadget.snapshot.File(path)
    masses = snapshot.gas.get_masses()
    snapshot.gas.load_number_density()
    snapshot.gas.calculate_temperature()
    snapshot.close()
    # Refine
    minimum = numpy.amin(masses)
    refined = numpy.where(masses <= minimum)[0]
    snapshot.gas.ndensity = snapshot.gas.ndensity[refined]
    snapshot.gas.temp = snapshot.gas.temp[refined]
    # Cleanup to save memory
    del snapshot.gas.masses
    del snapshot.gas.abundances
    del snapshot.gas.energy
    del snapshot.gas.h2frac
    del snapshot.gas.gamma

    return snapshot

def plot_temp(snapshot,wpath):
    h = snapshot.header.HubbleParam
    redshift = snapshot.header.Redshift
    dens = snapshot.gas.get_number_density()
    temp = snapshot.gas.get_temperature()

    # Plot!
    fig = pyplot.figure(figsize=(15,10))
    fig.clf()
    ax = fig.add_subplot(111)
    ax.hexbin(dens,temp,gridsize=500,bins='log',xscale='log',yscale='log',
              mincnt=1)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(2e-3, 1e12)
    ax.set_ylim(10, 2e4)

    pyplot.title('Redshift: %.2f' %(redshift,))
    ax.set_xlabel('n [cm^-3]')
    ax.set_ylabel('Temperature [K]')
    pyplot.savefig(wpath+'-temp.png',
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
                p = mp.Process(target=plot_temp, args=(snapshot, wd))
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
        plot_temp(snapshot,write_dir+'{:0>4}'.format(snap))
        
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

