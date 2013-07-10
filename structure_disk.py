#!/usr/bin/env python
# structure_disk.py
# Jacob Hummel

import os
import sys
import glob
import numpy
import Queue
from matplotlib import pyplot

import pyGadget
global t0
#===============================================================================
def load_dens(fname,length_unit):
    global t0
    snapshot = pyGadget.snapshot.File(fname)
    snapshot.gas.load_masses()
    snapshot.gas.load_number_density()
    snapshot.gas.load_coords(length_unit)
    smL = snapshot.gas.get_smoothing_length(length_unit)
    sinkval = snapshot.gas.get_sinks()

    ### Sinks!
    ids = snapshot.sink_ids = numpy.where(sinkval != 0)[0] 
    if snapshot.sink_ids.size > 0:
        print snapshot.sink_ids.size,'sinks found.'
        if t0 is None:
            t0 = snapshot.header.Time*pyGadget.units.Time_yr
        # If sink smooting length is too inflated, artificially shrink it.
        snapshot.gas.smoothing_length[ids] *= 1e-6

    return snapshot

def project(snap, write_dir, boxsize, length_unit, *args, **kwargs):
    global t0
    folder = 'disk/'
    for view in ['xy', 'xz', 'yz']:
        wpath = write_dir + folder + view + '/'
        if not os.path.exists(wpath):
            os.makedirs(wpath)
        wpath += '{:0>4}-disk-'.format(snap.number) + view + '.png'
        x,y,z = pyGadget.visualize.density_projection(snap, boxsize, 1., 
                                                      length_unit, view, 'halo',
                                                      *args, **kwargs)
                                                      
        z = numpy.log10(z)
        zmin,zmax = (1e8,1e12)
        
        print 'Plotting '+view+'...'
        fig = pyplot.figure(1,(16,12))
        fig.clf()
        pyplot.imshow(z, extent=[x.min(),x.max(),y.min(),y.max()],
                      cmap=pyGadget.colormap.get_cmap('smooth','./colormap'))
        pyplot.clim(numpy.log10(zmin),numpy.log10(zmax))
        cb = pyplot.colorbar()
        ax = pyplot.gca()
        for sink in snap.sinks:
            ax.plot(sink[1], -sink[0], 'k+', ms=7, mew=1.5) #90-degree rotation
            ax.text(sink[1]+10, -sink[0]+5, '%.1f' %sink[3])
        ax.set_xlim(x.min(),x.max())
        ax.set_ylim(y.min(),y.max())
        ax.set_xlabel('AU')
        ax.set_ylabel('AU')
        cb.set_label('Log Number Density [cm$^{-3}$]')
        density_labels = (8,9,10,11,12)
        cb.set_ticks(density_labels)
        ax.text(-boxsize*.475,boxsize*.4625,'z: %.2f' %snap.header.Redshift,
                color='white',fontsize=18)
        if t0:
            t_acc = snap.header.Time*pyGadget.units.Time_yr - t0
            ax.text(boxsize*.275,boxsize*.4625,'t$_{form}$: %.1f yr' %t_acc,
                    color='white',fontsize=18)
        pyplot.draw()
        pyplot.savefig(wpath, 
                       bbox_inches='tight')
    snap.close()

def multitask(path, write_dir, start, stop, boxsize,length_unit,**kwargs):
    global t0
    file_queue = Queue.Queue()
    data_queue = Queue.Queue(3)
    # Start the file-load thread
    pyGadget.snapshot.Loader(load_dens, file_queue, data_queue).start()
    # Add filenames to the queue
    for i in xrange(start,stop+1):
        fname = path + '{:0>3}'.format(i) + '.hdf5'
        file_queue.put((fname,length_unit))
    file_queue.put(None)
    # Process images
    while 1:
        snapshot = data_queue.get()
        if snapshot is None:
            break # reached end of queue!
        project(snapshot,write_dir,boxsize,length_unit,**kwargs)
    print 'Done.'
    
#===============================================================================
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

length_unit = pyGadget.units.Length_AU
boxsize = 5e3
### Optional arguments (If you want to override defaults.)
pps = 500  # 'pixels' per side
sm = 1.7   # smoothing factor
dlim = 1e11 # density limit for finding halo center
np = 100  # number of particles to require for finding halo center
kwargs = {'pps':pps, 'sm':sm, 'dens_limit':dlim, 'nparticles':np}

try:
    sink = pyGadget.sinks.SingleSink(sinkpath)
    t0 = sink.time[0]
except IOError:
    t0 = None

if len(sys.argv) == 3:
    snap = sys.argv[2]
    fname = path + '{:0>3}'.format(snap)+'.hdf5'
    print 'loading', fname
    snapshot = load_dens(fname,length_unit)
    project(snapshot,write_dir,boxsize,length_unit,**kwargs)

elif len(sys.argv) == 4:
    start = int(sys.argv[2])
    stop = int(sys.argv[3])
    multitask(path,write_dir,start,stop,boxsize,length_unit,**kwargs)

else:
    files0 = glob.glob(path+'???.hdf5')
    files1 = glob.glob(path+'1???.hdf5')
    files0.sort()
    files1.sort()
    start = int(files0[0][-8:-5])
    stop = int(files0[-1][-8:-5])
    if files1:
        stop = int(files1[-1][-9:-5])
        
    multitask(path,write_dir,start,stop,boxsize,length_unit,**kwargs)
