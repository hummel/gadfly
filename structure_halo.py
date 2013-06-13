#!/usr/bin/env python
# density_vis.py
# Jacob Hummel

import os
import sys
import glob
import numpy
import Queue
from matplotlib import pyplot
from matplotlib import cm

import pyGadget
global t0
#===============================================================================
def load_dens(fname,length_unit):
    global t0
    snapshot = pyGadget.snapshot.File(fname)
    snapshot.gas.load_masses()
    snapshot.gas.load_number_density()
    snapshot.gas.load_coords(length_unit,comoving=False)
    smL = snapshot.gas.get_smoothing_length(length_unit,comoving=False)
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

def project(snap, write_dir, scale, cscale, *args, **kwargs):
    global t0
    boxsize = int("".join(ch if ch.isdigit() else "" for ch in scale))
    unit = "".join(ch if not ch.isdigit() else "" for ch in scale)
    length_unit = pyGadget.units.Lengths[unit]

    #(Scaling boxsize to physical size at redshift 25.)
    boxsize = 25 * boxsize / (1 + snap.header.Redshift)
    for suffix in ['-halo_'+scale+'-xy.png']:
        wpath = write_dir + '{:0>4}'.format(snap.number) + suffix
        view = suffix[-6:-4]
        x,y,z = pyGadget.visualize.density_projection(snap, boxsize, 1., 
                                                      length_unit, view, 'halo',
                                                      *args,**kwargs)
        z = numpy.log10(z)
        print 'Plotting '+view+'...'
        fig = pyplot.figure(1,(16,12))
        fig.clf()
        pyplot.imshow(z, extent=[x.min(),x.max(),y.min(),y.max()], cmap=cm.jet)
        pyplot.clim(cscale[0],cscale[1])
        pyplot.colorbar()
        ax = pyplot.gca()
        ax.set_xlim(x.min(),x.max())
        ax.set_ylim(y.min(),y.max())
        ax.set_xlabel('Physical Distance ['+unit+']')
        ax.set_ylabel('Physical Distance ['+unit+']')
        ax.text(-950,925,'z: %.2f' %snap.header.Redshift,
                color='white',fontsize=18)
        if t0:
            t_acc = snap.header.Time*pyGadget.units.Time_yr - t0
            ax.text(550,925,'t$_{form}$: %.1f yr' %t_acc,
                    color='white',fontsize=18)
        pyplot.title('z: %.2f' %snap.header.Redshift)
        pyplot.draw()
        pyplot.savefig(wpath, 
                       bbox_inches='tight')
    snap.close()

def multitask(path, start, stop, *args,**kwargs):
    global t0
    unit = "".join(ch if not ch.isdigit() else "" for ch in args[1])
    length_unit = pyGadget.units.Lengths[unit]

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
        project(snapshot,*args,**kwargs)
    print 'Done.'
    
#===============================================================================
if ((len(sys.argv) not in [3,4,5]) or (sys.argv[1] == '-h')):
    print 'Usage::'
    print '   Option 1: python gas_temp.py [simulation name] [length scale] '\
        '[beginning snapshot] [final snapshot]'
    print '   Option 1: python gas_temp.py [simulation name] [length scale] '\
        '[single snapshot]'
    print '   Option 1: python gas_temp.py [simulation name] [length scale] '\
        '(creates a plot for every snapshot)'
    sys.exit()

simulation = sys.argv[1]
path = os.getenv('HOME')+'/sim/'+simulation+'/snapshot_'
sinkpath = os.getenv('HOME')+'/data/sinks/'+simulation+'/'
write_dir = os.getenv('HOME')+'/data/simplots/'+simulation+'/'
colors = {'1kpc': (-2.5,3), '1000pc':(-2.5,3), '500pc':(-2.5,3), 
          '100pc':(-1,3.5), '10pc':(1,6), '1pc':(3.5,9)}
scaling = sys.argv[2]
cscaling = colors[scaling]

### Optional arguments (If you want to override defaults.)
pps = 500  # 'pixels' per side
sm = 1.7   # smoothing factor
dlim = 1e6 # density limit for finding halo center
np = 1000  # number of particles to require for finding halo center
kwargs = {'pps':pps, 'sm':sm, 'dens_limit':dlim, 'nparticles':np}

try:
    sink = pyGadget.sinks.SingleSink(sinkpath)
    t0 = sink.time[0]
except IOError:
    t0 = None

if len(sys.argv) == 4:
    snap = sys.argv[3]
    fname = path + '{:0>3}'.format(snap)+'.hdf5'
    print 'loading', fname
    unit = "".join(ch if not ch.isdigit() else "" for ch in scaling)
    length_unit = pyGadget.units.Lengths[unit]
    snapshot = load_dens(fname,length_unit)
    project(snapshot,write_dir,scaling,cscaling,**kwargs)

elif len(sys.argv) == 5:
    start = int(sys.argv[3])
    stop = int(sys.argv[4])
    multitask(path,start,stop,write_dir,scaling,cscaling,**kwargs)

else:
    files0 = glob.glob(path+'???.hdf5')
    files1 = glob.glob(path+'1???.hdf5')
    files0.sort()
    files1.sort()
    start = int(files0[0][-8:-5])
    stop = int(files0[-1][-8:-5])
    if files1:
        stop = int(files1[-1][-9:-5])
    multitask(path,start,stop,write_dir,scaling,cscaling,**kwargs)
