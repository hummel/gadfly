#!/usr/bin/env python
# density_vis.py
# Jacob Hummel

import os
import sys
import glob
import numpy
import Queue
import multiprocessing as mp

import pyGadget

#===============================================================================
def load_data(fname,length_unit,mass_unit):
    snapshot = pyGadget.snapshot.File(fname)
    snapshot.dm.load_masses(mass_unit)
    snapshot.dm.load_coords(length_unit)
    snapshot.gas.load_masses(mass_unit)
    snapshot.gas.load_coords(length_unit)
    snapshot.gas.load_number_density()
    return snapshot

#===============================================================================
#===============================================================================
def multitask(path,write_dir,start,stop,length_unit,mass_unit):
    maxprocs = mp.cpu_count()
    file_queue = Queue.Queue()
    data_in = Queue.Queue(2)
    data_out = Queue.Queue()
    pyGadget.snapshot.Loader(load_data, file_queue, data_in).start()
    for snap in xrange(start,stop+1):
        fname = path + '{:0>3}'.format(snap)+'.hdf5'
        file_queue.put((fname, length_unit, mass_unit))
    file_queue.put(None)
    procs = []
    done = False
    while not done:
        snapshot = data_in.get()
        if snapshot is None:
            done = True
            for process in procs:
                process.join()
        else:
            p = mp.Process(target=pyGadget.analyze.halo_properties,
                           args=(snapshot,))
                                              
                                              
            procs.append(p)
            while True:
                running_procs = 0
                for proc in procs:
                    if proc.is_alive():
                        running_procs +=1
                if running_procs < maxprocs:
                    print 'Plotting', snapshot.filename
                    p.start()
                    break

#===============================================================================
if __name__ == '__main__':
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
sinkpath = os.getenv('HOME')+'/data/sinks/'+simulation+'/'
write_dir = os.getenv('HOME')+'/data/simplots/'+simulation+'/'

mass_unit = pyGadget.units.Mass_g
length_unit = pyGadget.units.Length_cm

if len(sys.argv) == 3:
    snap = sys.argv[2]
    fname = path + '{:0>3}'.format(snap)+'.hdf5'
    print 'loading', fname
    snapshot = load_data(fname,length_unit, mass_unit)
    hprops = pyGadget.analyze.halo_properties(snapshot)

        
elif len(sys.argv) == 4:
    start = int(sys.argv[2])
    stop = int(sys.argv[3])+1
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



