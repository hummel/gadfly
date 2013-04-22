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
    snapshot.close()
    # Refine
    minimum = numpy.amin(snapshot.dm.masses)
    refined = numpy.where(snapshot.dm.masses <= minimum)[0]
    snapshot.dm.masses = snapshot.dm.masses[refined]
    snapshot.dm.coordinates = snapshot.dm.coordinates[refined]
    minimum = numpy.amin(snapshot.gas.masses)
    refined = numpy.where(snapshot.gas.masses <= minimum)[0]
    snapshot.gas.masses = snapshot.gas.masses[refined]
    snapshot.gas.coordinates = snapshot.gas.coordinates[refined]
    snapshot.gas.ndensity = snapshot.gas.ndensity[refined]
    return snapshot

#===============================================================================
def analyze_halo(snapshot, write_dir):
    halo = pyGadget.analyze.halo_properties(snapshot,verbose=False)
    numpy.save(write_dir+'/haloz/'+str(snapshot.number)+'.npy',halo)
    del snapshot.dm.masses
    del snapshot.dm.coordinates
    del snapshot.gas.masses
    del snapshot.gas.coordinates
    del snapshot.gas.ndensity

#===============================================================================
def compile_halos(directory):
    files0 = glob.glob(directory+'???.npy')
    files1 = glob.glob(directory+'1???.npy')
    files0.sort()
    files1.sort()
    files = files0 + files1
    data = []
    for f in files:
        halo = numpy.load(f)
        if halo.size > 0:
            data.append(halo)

    array_lengths = [x.shape[0] for x in data]
    maxL = max(array_lengths)
    total = len(data)
    for i in range(total):
        data[i].resize([maxL,7], refcheck=False)
    datarray = numpy.concatenate([x for x in data])
    datarray = datarray.reshape(total,maxL,7)
    return datarray

#===============================================================================
def multitask_serial(path,write_dir,start,stop,length_unit,mass_unit):
    file_queue = Queue.Queue()
    data_queue = Queue.Queue(3)
    pyGadget.snapshot.Loader(load_data, file_queue, data_queue).start()
    for snap in xrange(start,stop+1):
        fname = path + '{:0>3}'.format(snap)+'.hdf5'
        file_queue.put((fname, length_unit, mass_unit))
    file_queue.put(None)

    jobs = []
    halos = []
    done = False
    while not done:
        snapshot = data_queue.get()
        if snapshot is None:
            done = True
        else:
            halo = pyGadget.analyze.halo_properties(snapshot,verbose=False)
            halos.append(halo)
    return halos

#===============================================================================
def multitask_parallel(path,write_dir,start,stop,length_unit,mass_unit):
    maxjobs = mp.cpu_count()
    file_queue = Queue.Queue()
    data_queue = Queue.Queue(3)
    pyGadget.snapshot.Loader(load_data, file_queue, data_queue).start()
    for snap in xrange(start,stop+1):
        fname = path + '{:0>3}'.format(snap)+'.hdf5'
        file_queue.put((fname, length_unit, mass_unit))
    file_queue.put(None)

    jobs = []
    done = False
    while not done:
        snapshot = data_queue.get()
        if snapshot is None:
            done = True
        else:
            p = mp.Process(target=analyze_halo, args=(snapshot,write_dir))
            jobs.append(p)
            p.start()
    for process in jobs:
        process.join()

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
    stop = int(sys.argv[3])
    multitask_parallel(path,write_dir,start,stop,length_unit,mass_unit)
    data = compile_halos(write_dir)
    numpy.save(write_dir+'halo_properties.npy',data)
    
else:
    files0 = glob.glob(path+'???.hdf5')
    files1 = glob.glob(path+'1???.hdf5')
    files0.sort()
    files1.sort()
    files = files0 + files1
    start = int(files[0][-8:-5])
    stop = int(files[-1][-8:-5])
    multitask_parallel(path,write_dir,start,stop,length_unit,mass_unit)
    data = compile_halos(write_dir)
    numpy.save(write_dir+'halo_properties.npy',data)
