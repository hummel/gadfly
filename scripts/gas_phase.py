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
    snapshot = pyGadget.snapshot.File(path)
    if key == 'temp':
        snapshot.gas.load_data('ndensity','temp')
    elif key == 'frac':
        snapshot.gas.load_data('ndensity','temp','electron_frac',
                               'h2frac','HDfrac')
    else:
        raise KeyError
    snapshot.close()
    return snapshot

def figure_temp(snapshot,wpath):
    fig = pyGadget.plotting.Phase(snapshot)
    fig.plot('temp')
    fig.save(wpath+'-temp.png')

def figure_gas_fraction(snapshot,wpath):
    fig = pyGadget.plotting.Quad(snapshot)
    fig.plot('temp','electron_frac','h2frac','HDfrac')
    fig.save(wpath+'-frac.png')

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

